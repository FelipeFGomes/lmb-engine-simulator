use crate::base::constants::{_P_REF, _T_REF};
use crate::connector::valve::Valve;
use crate::core::traits::ZeroDim;
use crate::reaction::combustion;
use crate::reaction::combustion::{Combustion, WiebeFunction};
use crate::reaction::gas::Gas;
use crate::zero_dim::cylinder::Cylinder;
use super::json_reader::{JsonEngine, JsonFuel};
use crate::FlowRatio;
use ndarray::*;
use serde_json;
use std::f64::consts::PI;
use std::io::Write;

pub struct Engine {
    speed: f64,
    displacement: f64,
    bore: f64,
    eccentricity: f64,
    firing_order: String,
    combustion: Box<dyn Combustion>,
    injector: Option<Injector>,
    cylinders: Vec<Cylinder>,
    valves: Vec<Valve>,
    operat_param: OperationalParameters,
    sec_to_rad: f64,
}

type EngineOutput = Result<Engine, String>;

impl Engine {
    pub fn new(file_name: &str, gas: &Gas) -> EngineOutput {
        let json_engine = match Engine::reading_json(file_name) {
            Ok(eng) => eng,
            Err(err) => {
                let msg = format!(
                    "Error at Engine::new \n unable to parse {} file \n{}",
                    file_name, err
                );
                return Err(msg);
            }
        };

        let firing_order: Vec<f64> = json_engine
            .firing_order
            .split('-')
            .map(|s| -> f64 {
                match s.parse() {
                    Ok(s) => s,
                    Err(err) => {
                        eprintln!("Error at Engine:new \n unable to parse {} \n{}", s, err);
                        std::process::exit(1);
                    }
                }
            })
            .collect();

        if firing_order.len() != json_engine.cylinders.len() {
            let msg = format!(
                "Error at Engine::new 
            `firing_order` and `cylinders` must have the same length.
            `firing_order`: {}
            `cylinder`: {}",
                firing_order.len(),
                json_engine.cylinders.len()
            );
            return Err(msg);
        }

        let num_divisions = firing_order.iter().cloned().fold(0. / 0., f64::max); // gets highest value
        let division = 720.0 / num_divisions;

        // checking condition
        if json_engine.combustion.is_some() && json_engine.injector.is_none() {
            let msg = format!(
                "
                Combustion model can only be created if an \"Injector\" has been added.
                Add the \"Injector\" field on file: \"{}\"
                ",
                file_name
            );
            return Err(msg);
        }

        // instantiating combustion
        let mut air_gas = gas.clone();
        air_gas.X("O2:0.21, N2:0.79");
        let mut injector: Option<Injector> = None;
        let combustion: Box<dyn Combustion>;
        if let Some(inj_json) = &json_engine.injector {
            // creating Fuel obj
            let fuel = Fuel::new(&inj_json.fuel, gas);
            // creating Injector obj
            injector = Some(Injector::new(
                inj_json.inj_type.clone(),
                inj_json.air_fuel_ratio,
                &fuel,
                &air_gas,
            ));
            if let Some(comb) = &json_engine.combustion {
                if !gas.contains_specie(&inj_json.fuel.name) {
                    return Err(format!(
                        "the gas added to the cylinder does not contain specie \"{}\"",
                        inj_json.fuel.name
                    ));
                }
                // Combustion models
                if comb.model == "Two-zone model" {
                    let wiebe =
                        WiebeFunction::new(comb.wiebe.a, comb.wiebe.m, comb.wiebe.comb_duration);
                    combustion = Box::new(combustion::TwoZoneCombustion::new(
                        comb.comb_ini,
                        inj_json.air_fuel_ratio,
                        wiebe,
                        gas,
                        &fuel,
                        &air_gas,
                    )?);
                } else {
                    println!("WARNING! Model `{}` not found", comb.model);
                    println!("WARNING! No combustion will be used instead");
                    combustion = Box::new(combustion::NoCombustion::new());
                }
            } else {
                combustion = Box::new(combustion::NoCombustion::new());
            }
        } else {
            combustion = Box::new(combustion::NoCombustion::new());
        }

        // instantianting cylinders
        let mut cyl_names: Vec<String> = Vec::new();
        let mut cylinders: Vec<Cylinder> = Vec::new();
        for (cylinder, order) in json_engine.cylinders.iter().zip(firing_order.iter()) {
            cyl_names.push(cylinder.name.clone());
            cylinders.push(Cylinder::new(
                cylinder.name.clone(),
                180.0 + (order - 1.0) * division,
                &json_engine,
                &cylinder,
                &cylinder.intake_valves,
                &cylinder.exhaust_valves,
                combustion.clone(),
                injector.clone(),
                gas,
            )?);
        }

        // instantianting valves
        let mut valve_names: Vec<String> = Vec::new();
        let mut valves: Vec<Valve> = Vec::new();
        for cylinder in json_engine.cylinders.iter() {
            // intake valves
            for valve in cylinder.intake_valves.iter() {
                valve_names.push(valve.name.clone());
                let cyl = match cylinders.iter().find(|c| c.name() == cylinder.name) {
                    Some(c) => c,
                    None => {
                        let msg = format!("cylinder `{}` has not been added", cylinder.name);
                        return Err(msg);
                    }
                };
                valves.push(Valve::new(
                    valve.name.clone(),
                    valve.opening_angle,
                    valve.closing_angle,
                    valve.diameter * 1e-3,
                    valve.max_lift * 1e-3,
                    cyl,
                )?);
            }

            // exhaust valves
            for valve in cylinder.exhaust_valves.iter() {
                valve_names.push(valve.name.clone());
                let cyl = match cylinders.iter().find(|c| c.name() == cylinder.name) {
                    Some(c) => c,
                    None => {
                        let msg = format!("cylinder `{}` has not been added", cylinder.name);
                        return Err(msg);
                    }
                };
                valves.push(Valve::new(
                    valve.name.clone(),
                    valve.opening_angle,
                    valve.closing_angle,
                    valve.diameter * 1e-3,
                    valve.max_lift * 1e-3,
                    cyl,
                )?);
            }
        }

        let engine = Engine {
            speed: json_engine.speed,
            displacement: json_engine.displacement,
            bore: json_engine.bore,
            eccentricity: json_engine.eccentricity,
            firing_order: json_engine.firing_order.clone(),
            sec_to_rad: 2.0 * PI * json_engine.speed / 60.0,
            cylinders,
            valves,
            combustion,
            injector: injector,
            operat_param: OperationalParameters::new(),
        };
        Ok(engine)
    }

    pub fn advance(&mut self, dt: f64) {
        self.cylinders.iter_mut().for_each(|cyl| cyl.advance(dt));
    }

    pub fn calc_operational_param(
        &mut self,
        press: Vec<ArrayView1<f64>>, // bar
        vol: Vec<ArrayView1<f64>>,   // cm^3
    ) {
        let mut power = 0.0f64;
        let mut total_work = 0.0f64;
        for (p, v) in press.iter().zip(vol) {
            let mut dv: Array1<f64> = Array1::zeros(v.len() - 1);
            let mut p_mean: Array1<f64> = Array1::zeros(p.len() - 1);
            Zip::from(&mut dv)
                .and(v.slice(s![1..]))
                .and(v.slice(s![0..v.len() - 1]))
                .apply(|dv, v2, v1| *dv = (v2 - v1) * 1e-6);
            Zip::from(&mut p_mean)
                .and(p.slice(s![1..]))
                .and(p.slice(s![0..p.len() - 1]))
                .apply(|pm, p2, p1| *pm = 0.5 * (p2 + p1) * 1e5);
            let work = dv.dot(&(p_mean - _P_REF));
            total_work += work;
            power += work * self.speed / 120.0; // W
        }

        let mut total_fuel_mass = 0.0f64;
        self.cylinders()
            .iter()
            .for_each(|c| total_fuel_mass += c.fuel_mass());
        // println!("fuel mass: {:.4} [mg]", total_fuel_mass*1e6);

        let mut mass: f64 = 0.0;
        self.cylinders()
            .iter()
            .for_each(|c| mass += c.closed_phase_mass());
        let total_disp = self.displacement * (self.cylinders().len() as f64) * 1e-6;
        let r_ref = 287.0;
        let vol_effic = 100.0 * mass / (_P_REF * total_disp / (r_ref * _T_REF));

        let mut residual_mass: f64 = 0.0;
        self.cylinders()
            .iter()
            .for_each(|c| residual_mass += c.residual_mass_frac());
        residual_mass = 100.0 * residual_mass / (self.cylinders().len() as f64);

        let torque = power / (self.speed * PI / 30.0);
        let imep = total_work / self.displacement * 10.0; // displ is in cm³
        let thermal_effic: f64;
        if let Some(injector) = &self.injector {
            thermal_effic = 100.0 * total_work / (total_fuel_mass * injector.fuel().lhv());
        } else {
            thermal_effic = 0.0
        }

        self.operat_param.speed.push(self.speed);
        self.operat_param.power.push(power);
        self.operat_param.torque.push(torque);
        self.operat_param.imep.push(imep);
        self.operat_param.thermal_effic.push(thermal_effic);
        self.operat_param.vol_effic.push(vol_effic);
        self.operat_param.residual_mass.push(residual_mass);
    }

    pub fn write_performance_to(&self, file_name: &str) {
        let op = &self.operat_param;
        let mut data: Vec<String> = Vec::new();
        for i in 0..self.operat_param.speed.len() {
            data.push(format!(
                "{:.1}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\n",
                op.speed[i],
                op.power[i],
                op.torque[i],
                op.imep[i],
                op.thermal_effic[i],
                op.vol_effic[i],
                op.residual_mass[i],
            ));
        }
        let header = "Speed [RPM]\tPower [W]\tTorque [Nm]\tIMEP [bar]\tEfficiency [%]\tVolumetric effic [%]\tResidual mass [%]".to_string();
        let mut file = std::fs::File::create(file_name).expect("Error opening writing file");
        write!(file, "{}\n", header).expect("Unable to write data");
        write!(file, "{}", data.join("")).expect("Unable to write data");
    }

    pub fn update_cylinders_flow_ratio(&mut self, flow_ratio: Vec<Vec<(&str, &FlowRatio)>>) {
        for (cylinder, flow) in self.cylinders.iter_mut().zip(flow_ratio) {
            cylinder.update_flow_ratio(flow);
        }
    }
    /// Set engine speed, input in RPM
    pub fn set_speed(&mut self, speed: f64) {
        self.speed = speed;
        self.sec_to_rad = 2.0 * PI * speed / 60.0;
        self.cylinders.iter_mut().for_each(|c| c.set_speed(speed));
    }

    /// Set engine displacement, input in cm³
    pub fn set_displacement(&mut self, disp: f64) {
        self.displacement = disp;
        self.cylinders
            .iter_mut()
            .for_each(|c| c.set_displacement(disp * 1e-6));
    }

    /// Set cylinder `cyl` compression ratio
    pub fn set_compression_ratio_of(&mut self, cyl: &str, comp_ratio: f64) {
        if let Some(cylinder) = self.cylinders.iter_mut().find(|c| c.name() == cyl) {
            cylinder.set_compression_ratio(comp_ratio);
        } else {
            println!("Error at Engine::set_compression_ratio_of()");
            println!(" \"{}\" was not found", cyl);
            std::process::exit(1);
        }
    }

    /// Set combustion model for all cylinders
    pub fn set_combustion_model(&mut self, comb: &Box<dyn Combustion>) {
        if let Some(_) = &self.injector {
            self.combustion = comb.clone();
            self.cylinders
                .iter_mut()
                .for_each(|c| c.set_combustion_model(comb.clone()));
        } else {
            println!(
                "Error at Engine::set_combustion_model(): Engine does not contain an Injector"
            );
            println!(" A combustion model cannot be added to an engine without an Injector");
            std::process::exit(1);
        }
        self.combustion = comb.clone();
    }

    /// Set injectors relative air-fuel ratio, input between 0 and 1
    pub fn set_air_fuel_ratio(&mut self, afr: f64) {
        if let Some(inj) = &mut self.injector {
            inj.set_air_fuel_ratio(afr);
        } else {
            println!("Error at Engine:: set_air_fuel_ratio(): Injector does not exist");
            std::process::exit(1);
        }
    }

    /// Set if the species inside the cylinder should be storable
    pub fn store_species(&mut self, state: bool) {
        self.cylinders
            .iter_mut()
            .for_each(|c| c.set_store_species(state));
    }

    /// Returns a reference to the cylinders in the engine
    pub fn cylinders<'a>(&'a self) -> &'a Vec<Cylinder> {
        &self.cylinders
    }

    /// Returns a reference to the valves in the engine
    pub fn valves<'a>(&'a self) -> &'a Vec<Valve> {
        &self.valves
    }

    /// Returns cylinders bore in mm
    pub fn bore(&self) -> f64 {
        self.bore
    }

    /// Returns the firing order of the cylinders
    pub fn firing_order(&self) -> String {
        self.firing_order.clone()
    }

    /// Returns cylinders eccentricity in mm
    pub fn eccentricity(&self) -> f64 {
        self.eccentricity
    }

    /// Returns the combustion model name
    pub fn combustion_model(&self) -> String {
        self.combustion.model_name().to_string()
    }

    pub fn injector(&self) -> Option<&Injector> {
        if let Some(inj) = &self.injector {
            Some(inj)
        } else {
            None
        }
    }

    /// Returns the constant multiplier to change from seconds to crank-angle degrees
    pub fn sec_to_rad(&self) -> f64 {
        self.sec_to_rad
    }

    pub fn operat_param(&self) -> &OperationalParameters {
        &self.operat_param
    }

    fn reading_json(file_name: &str) -> serde_json::Result<JsonEngine> {
        let json_file = std::fs::read_to_string(file_name).expect("Unable to read file");
        let data: JsonEngine = serde_json::from_str(&json_file)?;
        Ok(data)
    }
}

pub struct OperationalParameters {
    speed: Vec<f64>,
    power: Vec<f64>,
    torque: Vec<f64>,
    imep: Vec<f64>,
    thermal_effic: Vec<f64>,
    vol_effic: Vec<f64>,
    residual_mass: Vec<f64>,
}

impl OperationalParameters {
    fn new() -> OperationalParameters {
        OperationalParameters {
            speed: Vec::new(),
            power: Vec::new(),
            torque: Vec::new(),
            imep: Vec::new(),
            thermal_effic: Vec::new(),
            vol_effic: Vec::new(),
            residual_mass: Vec::new(),
        }
    }
}

impl std::fmt::Display for OperationalParameters {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "
            Speed [RPM]: {:.1?}
            Power [W]:\t {:.2?}
            Torque [Nm]: {:.2?}
            IMEP [bar]:  {:.2?}
            effic [%]:\t {:.2?}
            vol_effic [%]: {:.2?}
            residual_mass [%]: {:.2?}
            ",
            self.speed,
            self.power,
            self.torque,
            self.imep,
            self.thermal_effic,
            self.vol_effic,
            self.residual_mass,
        )
    }
}
#[derive(Debug, Clone)]
pub struct Fuel {
    name: String,
    mole_weight: f64,
    air_moles: f64,
    lhv: f64,
    heat_vap: f64,
    comp: String,
    carbon: f64,
    hydrogen: f64,
    oxigen: f64,
    nitrogen: f64,
}

impl Fuel {
    fn new(fuel: &JsonFuel, gas: &Gas) -> Fuel {
        let mole_weight = gas.mole_weight_of(&fuel.name);
        let atoms = gas.atoms_of(&fuel.name);
        let lhv = match fuel.lhv {
            Some(lhv) => lhv,
            None => Fuel::calc_low_heat_value(&fuel.name),
        };

        let heat_vap = match fuel.heat_vap {
            Some(h_vap) => h_vap,
            None => 0.0,
        };

        let c = match atoms.get("C") {
            Some(x) => *x,
            None => 0.0,
        };
        let h = match atoms.get("H") {
            Some(x) => *x,
            None => 0.0,
        };
        let o = match atoms.get("O") {
            Some(x) => *x,
            None => 0.0,
        };

        let n = match atoms.get("N") {
            Some(x) => *x,
            None => 0.0,
        };

        let air_moles = c + 0.25 * h - 0.5 * o;
        let comp = format!("{}:1.0", fuel.name);

        Fuel {
            name: fuel.name.to_string(),
            mole_weight,
            air_moles,
            lhv,
            heat_vap,
            comp,
            carbon: c,
            hydrogen: h,
            oxigen: o,
            nitrogen: n,
        }
    }

    fn calc_low_heat_value(name: &str) -> f64 {
        if name == "CH4" {
            50e06
        } else if name == "C2H5OH" {
            25.858e06
        } else if name == "C8H18" {
            44.651e06
        } else {
            println!("Error! Fuel `{}` has no low heat value stored", name);
            std::process::exit(1);
        }
    }

    pub fn name<'a>(&'a self) -> &'a str {
        &self.name
    }
    pub fn mole_weight(&self) -> f64 {
        self.mole_weight
    }
    pub fn air_moles(&self) -> f64 {
        self.air_moles
    }
    pub fn lhv(&self) -> f64 {
        self.lhv
    }
    pub fn heat_vap(&self) -> f64 {
        self.heat_vap
    }
    pub fn composition<'a>(&'a self) -> &'a str {
        &self.comp
    }
    pub fn carbon(&self) -> f64 {
        self.carbon
    }
    pub fn hydrogen(&self) -> f64 {
        self.hydrogen
    }
    pub fn oxigen(&self) -> f64 {
        self.oxigen
    }
    pub fn nitrogen(&self) -> f64 {
        self.nitrogen
    }
}

#[derive(Debug, Clone)]
pub struct Injector {
    inj_type: String,
    air_fuel_ratio: f64,
    fuel: Fuel,
    air_comp: Gas,
    fuel_mass_frac: f64,
    injected_fuel: f64,
    afr_afr_stoich: f64,
}

impl Injector {
    pub fn new(inj_type: String, afr: f64, fuel: &Fuel, air_comp: &Gas) -> Injector {
        let air_molar_mass =
            (air_comp.mole_frac() / air_comp.mole_frac_of("O2") * air_comp.mole_weight()).sum();
        let total_mass = fuel.mole_weight() + afr * fuel.air_moles() * air_molar_mass;
        let fuel_mass_frac = fuel.mole_weight() / total_mass;
        let afr_stoich = fuel.air_moles() * air_molar_mass / fuel.mole_weight();
        let afr_afr_stoich = afr * afr_stoich;

        Injector {
            inj_type,
            air_fuel_ratio: afr,
            fuel: fuel.clone(),
            air_comp: air_comp.clone(),
            fuel_mass_frac,
            injected_fuel: 0.0,
            afr_afr_stoich,
        }
    }
    pub fn calc_port_injected_fuel(&self, intake_mass: f64) -> f64 {
        self.fuel_mass_frac * intake_mass
    }
    pub fn calc_direct_injected_fuel(&self, intake_mass: f64) -> f64 {
        intake_mass / self.afr_afr_stoich
    }
    pub fn set_injected_fuel(&mut self, fuel_mass: f64) {
        self.injected_fuel = fuel_mass;
    }
    pub fn inj_type<'a>(&'a self) -> &'a str {
        &self.inj_type
    }
    pub fn air_fuel_ratio(&self) -> f64 {
        self.air_fuel_ratio
    }
    pub fn fuel<'a>(&'a self) -> &'a Fuel {
        &self.fuel
    }
    pub fn injected_fuel(&self) -> f64 {
        self.injected_fuel
    }
    pub fn set_air_fuel_ratio(&mut self, afr: f64) {
        self.air_fuel_ratio = afr;
        // updating other properties
        let air_molar_mass = (self.air_comp.mole_frac() / self.air_comp.mole_frac_of("O2")
            * self.air_comp.mole_weight())
        .sum();
        let total_mass = self.fuel.mole_weight() + afr * self.fuel.air_moles() * air_molar_mass;
        let fuel_mass_frac = self.fuel.mole_weight() / total_mass;
        let afr_stoich = self.fuel.air_moles() * air_molar_mass / self.fuel.mole_weight();
        let afr_afr_stoich = afr * afr_stoich;
        self.fuel_mass_frac = fuel_mass_frac;
        self.afr_afr_stoich = afr_afr_stoich;
    }
}
