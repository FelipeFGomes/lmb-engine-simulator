use crate::numerics::ode_solvers as ode;
use crate::reaction::gas::Gas;
use crate::reaction::combustion::{Combustion};
use crate::engine::json_reader::{JsonEngine, JsonCylinder, JsonValve};
use crate::engine::engine::Injector;
use crate::core::traits::{ZeroDim, SaveData, ZeroD};
use crate::{BasicProperties, FlowRatio};
use ansi_term::Style;
use ndarray::*;
use std::f64::consts::PI;
use std::io::Write;

pub struct Cylinder {
    name: String,
    gas: Gas,
    mass: f64,       // [kg]
    volume: f64,     // [m³] - instant volume as function of crank angle
    angle: f64,      // [CA deg] - instant crank angle
    speed: f64,      // [RPS]
    fuel_mass: f64,  // [kg]
    sec_to_rad: f64, // constant: 2*PI*speed
    geometry: Geometry,
    crankshaft: Crankshaft,
    piston: Piston,
    head: Head,
    heat_transfer: HeatTransfer,
    injector: Option<Injector>,
    combustion: Box<dyn Combustion>,
    int_valves: ValvesInfo,
    exh_valves: ValvesInfo,
    store_species: bool,
    total_injected_fuel: f64,
    total_fresh_charge: f64,
    closed_phase_mass: f64,
    residual_mass_frac: f64,

    // blow_by: bool,
    // crevice: bool,
}

impl Cylinder {
    /// Creates a cylinder object. Inputs units must be: `mm`, `RPM` and `CA deg`.
    pub fn new(
        name: String,
        ini_angle: f64,
        engine_info: &JsonEngine,
        cylinder_info: &JsonCylinder,
        int_valves_info: &Vec<JsonValve>,
        exh_valves_info: &Vec<JsonValve>,
        combustion: Box<dyn Combustion>,
        injector: Option<Injector>,
        gas: &Gas,
    ) -> Result<Cylinder, String> {
        if engine_info.bore < 0.0 {
            return Err(format!("diameter cannot be lower than zero"));
        } else if engine_info.displacement < 0.0 {
            return Err(format!("displacement cannot be lower than zero"));
        } else if cylinder_info.compression_ratio < 0.0 {
            return Err(format!("compration ratio cannot be lower than zero"));
        } else if engine_info.speed < 0.0 {
            return Err(format!("engine speed cannot be lower than zero"));
        }

        // In SI uinits
        let speed = engine_info.speed;
        let conrod = engine_info.conrod * 1e-3; // [m]
        let eccentricity = engine_info.eccentricity * 1e-3; // [m]
        let diam = engine_info.bore * 1e-3; // [m]
        let displ = engine_info.displacement;
        let comp_ratio = cylinder_info.compression_ratio;
        let wall_temp = cylinder_info.wall_temperature; // [K] - Common value
        let geometry = Geometry::new(diam, displ * 1e-6, comp_ratio, wall_temp);
        let crank = 0.5 * geometry.stroke; //standard value
        let angle_tdc = (eccentricity / (conrod + crank)).asin(); // crank angle [deg] of the top-dead-center
        let crankshaft = Crankshaft::new(conrod, crank, eccentricity, angle_tdc);
        let piston = Piston::new(geometry.stroke, diam, wall_temp, speed);
        let head = Head {
            temperature: wall_temp,
            area: geometry.transverse_area,
        };
        let (volume, _) = Cylinder::calc_volume(&geometry, &crankshaft, ini_angle.to_radians());

        // instantiating valves
        let intake_gas_comp = "O2:0.21, N2:0.79".to_string();
        let exhaust_gas_comp = "N2:0.662586, H2O:0.202449, CO2:0.134965".to_string();
        let mut intake_valves_info = Vec::new();
        let mut exhaust_valves_info = Vec::new();
        for v in int_valves_info {
            intake_valves_info.push(ValveBasicInfo::new(v.name.clone(), v.opening_angle, v.closing_angle))
        }
        for v in exh_valves_info {
            exhaust_valves_info.push(ValveBasicInfo::new(v.name.clone(), v.opening_angle, v.closing_angle))
        }

        let int_valves = ValvesInfo::new(intake_valves_info, intake_gas_comp.clone());
        let exh_valves = ValvesInfo::new(exhaust_valves_info, exhaust_gas_comp.clone());

        let store_species = match cylinder_info.store_species {
            Some(s) => s,
            None => false,
        };

        Ok(Cylinder {
            name,
            gas: gas.clone(),
            mass: gas.P() * volume / (gas.R() * gas.T()),
            volume,
            angle: (ini_angle.to_radians() - crankshaft.angle_tdc),
            speed: speed / 60.0,
            fuel_mass: 0.0,
            sec_to_rad: 2.0 * PI * speed / 60.0,
            geometry,
            crankshaft,
            piston,
            head,
            heat_transfer: HeatTransfer {},
            combustion,
            injector,
            int_valves,
            exh_valves,
            store_species,
            total_injected_fuel: 0.0,
            total_fresh_charge: 0.0,
            closed_phase_mass: 0.0,
            residual_mass_frac: 0.0,
        })
    }
    /// Returns the instant volume and volume's derivative with crank angle radian, respectively.
    /// Input `angle` must be in radian.
    /// Outputs are `m³` and `m³/CA-rad`
    fn calc_volume(geometry: &Geometry, crankshaft: &Crankshaft, angle: f64) -> (f64, f64) {
        let angle = angle - crankshaft.angle_tdc;
        let sin_gama =
            (crankshaft.crank * angle.sin() - crankshaft.eccentricity) / crankshaft.conrod;
        let cos_gama = (1.0 - sin_gama * sin_gama).sqrt();
        let piston_position = crankshaft.mech_total_length
            - crankshaft.crank * angle.cos()
            - crankshaft.conrod * cos_gama;
        let volume = geometry.transverse_area * piston_position + geometry.clearance;
        let tmp = (crankshaft.crank / crankshaft.conrod) * (angle.cos() / cos_gama);
        let dpiston_position_dangle =
            crankshaft.crank * angle.sin() * (1.0 + tmp) - crankshaft.eccentricity * tmp;
        let d_volume = geometry.transverse_area * dpiston_position_dangle;
        (volume, d_volume)
    }

    /// `d_angle` in crank angle radians
    fn closed_phase(&mut self, d_angle: f64) -> (f64, f64, f64, f64, Array1<f64>) {
        // Closed Phase -----------------------------------------------------------
        if let Some(inj) = &mut self.injector {
            if inj.inj_type() == "port" {
                if self.fuel_mass == 0.0 {
                    self.fuel_mass = inj.injected_fuel();
                    self.total_injected_fuel = self.fuel_mass;
                    self.closed_phase_mass = self.mass;
                    self.residual_mass_frac = 1.0 - (self.total_fresh_charge/self.mass);
                    self.total_fresh_charge = 0.0;
                } else {
                    inj.set_injected_fuel(0.0);
                }
            } else if inj.inj_type() == "direct" {
                if self.fuel_mass == 0.0 {
                    self.fuel_mass = inj.calc_direct_injected_fuel(0.0);
                    self.total_injected_fuel = self.fuel_mass;
                    self.closed_phase_mass = self.mass;
                    self.residual_mass_frac = 1.0 - self.total_fresh_charge/self.mass;
                    self.total_fresh_charge = 0.0;
                } else {
                    inj.set_injected_fuel(0.0);
                }
            }
        } else {
            self.closed_phase_mass = self.mass;
            self.total_fresh_charge = 0.0;
        }
        let heat_combustion = self.combustion.get_heat_release_rate(&self.gas, self.fuel_mass, self.angle);
        let const_1 = 1.0 / (self.mass * self.gas.cv());
        let closed_phase_equations = |angle: &f64, x: &Array1<f64>, _: &Vec<f64>| -> Array1<f64> {
            // x[0] = P
            let (vol, d_vol) = Cylinder::calc_volume(&self.geometry, &self.crankshaft, *angle);    
            let temp = x[0]*vol/(self.mass*self.gas.R());
            let heat_transfer = self.heat_transfer.calculate(vol, temp, x[0], &self); // [J/s]
            let heat_transfer = heat_transfer / self.sec_to_rad; // [J/CA radian]
            let d_temp = const_1 * (heat_combustion + heat_transfer - x[0] * d_vol); // [K/CA radian]
            let d_press = x[0] * (d_temp / temp - d_vol / vol); // [Pa/CA radian]
            array![d_press]
        };
        let ini = array![self.gas.P()];

        let closed_phase_integrated = ode::rk4_step(
            closed_phase_equations,
            &ini,
            &Vec::new(),
            &self.angle,
            d_angle,
        );

        // let temp = closed_phase_integrated[0];  
        let press = closed_phase_integrated[0];       
        let (vol, _) = Cylinder::calc_volume(&self.geometry, &self.crankshaft, self.angle + d_angle);
        let temp = press*vol/(self.mass*self.gas.R());

        // Estimating final compositions:
        let mole_frac = self.combustion.update_composition(&mut self.gas, self.mass, self.angle + d_angle, press, vol);
        ( temp, press, self.mass, vol, mole_frac )
    }

    /// `d_angle` in crank angle radians
    fn open_phase(&mut self, d_angle: f64) -> (f64, f64, f64, f64, Array1<f64>) {
        // Open Phase -----------------------------------------------------------
        let cv = self.gas.cv();
        let cv_inv = 1.0 / cv;
        let mass_flow = (self.int_valves.flow_info.mass_flow + self.exh_valves.flow_info.mass_flow) / self.sec_to_rad; // [kg/CA radian]
        let enthalpy_flow = (self.int_valves.flow_info.enthalpy_flow + self.exh_valves.flow_info.enthalpy_flow) / self.sec_to_rad; // [J/kg/CA radian]
        let open_phase_equations = |angle: &f64, x: &Array1<f64>, _: &Vec<f64>| -> Array1<f64> {
            //x[0] = temperature, x[1] = mass
            let (vol, d_vol) = Cylinder::calc_volume(&self.geometry, &self.crankshaft, *angle);
            let press = x[1] * self.gas.R() * x[0] / vol;
            let heat_transfer = self.heat_transfer.calculate(vol, x[0], press, &self); // [J/s]
            let heat_transfer = heat_transfer / self.sec_to_rad; // [J/CA radian]
            let d_mass = mass_flow; // [kg/CA radian]
            let d_temp = cv_inv / x[1]
                * (heat_transfer - press * d_vol + enthalpy_flow - cv * x[0] * d_mass); // [K/CA radian]
            array![d_temp, d_mass]
        };

        // println!("mass flow: {}\t enthalpy flow: {}", mass_flow, enthalpy_flow);
        let ini_condition = array![self.gas.T(), self.mass];
        let open_phase_integrated = ode::rk4_step(
            open_phase_equations,
            &ini_condition,
            &Vec::new(),
            &self.angle,
            d_angle,
        );
        let temp = open_phase_integrated[0];
        let mass = open_phase_integrated[1];
        let (vol, _) = Cylinder::calc_volume(&self.geometry, &self.crankshaft, self.angle + d_angle);
        let press = mass * self.gas.R() * temp / vol;

        // Estimating final compositions:
        let dt = d_angle/self.sec_to_rad;

        // summing the fresh charge of all intake valves
        let mut fresh_charge_mass = 0.0; // kg
        self.int_valves.basic_info.iter_mut().for_each(|v| {
            fresh_charge_mass += v.get_charge(dt);
        });
        self.total_fresh_charge += fresh_charge_mass;

        // injecting fuel
        self.fuel_mass = 0.0;
        let fuel_mass: f64;
        let additional_mass: Vec<(f64, &str)>;
        if let Some(inj) = &mut self.injector {
            if inj.inj_type() == "port" {
                fuel_mass = inj.calc_port_injected_fuel(fresh_charge_mass);
            } else if inj.inj_type() == "direct" {
                fuel_mass = 0.0;
            } else {panic!("Unknown injector type!")}
            inj.set_injected_fuel(fuel_mass + inj.injected_fuel());
            additional_mass = vec![(fresh_charge_mass - fuel_mass, &self.int_valves.gas_comp), (fuel_mass, inj.fuel().composition())];
        } else {
            additional_mass = vec![(fresh_charge_mass, &self.int_valves.gas_comp)];
        }
        let new_mole_frac = self.gas.if_mixed_with(self.mass, additional_mass); 

        ( temp, press, mass, vol, new_mole_frac )
    }

    pub fn fuel_mass(&self) -> f64 {
        if let Some(_) = &self.injector {
            self.total_injected_fuel
        } else {
            0.0
        }
    }

    pub fn closed_phase_mass(&self) -> f64 {self.closed_phase_mass}

    pub fn residual_mass_frac(&self) -> f64 {self.residual_mass_frac}

    pub fn set_speed(&mut self, speed: f64) {
        if speed.is_sign_negative() {
            println!("Error at Cylinder::set_speed()");
            println!(" speed must be a positive value! {}", speed);
            std::process::exit(1);
        }
        self.speed = speed;
        self.sec_to_rad = 2.0 * PI * speed / 60.0;
        self.piston.mean_velocity = 2.0 * self.geometry.stroke * speed / 60.0;
    }

    pub fn set_displacement(&mut self, disp: f64) {
        if disp.is_sign_negative() {
            println!("Error at Cylinder::set_displacement()");
            println!(" displacement must be a positive value! {}", disp);
            std::process::exit(1);
        }
        self.geometry.displacement = disp;
        self.geometry.stroke = disp / self.geometry.transverse_area;
        self.geometry.clearance = disp / (self.geometry.compression_ratio - 1.0);
        self.geometry.total_volume = disp + self.geometry.clearance;
        self.piston.mean_velocity = 2.0 * self.geometry.stroke * self.speed;
    }

    pub fn set_compression_ratio(&mut self, comp_ratio: f64) {
        if comp_ratio.is_sign_negative() {
            println!("Error at Cylinder::set_compression_ratio()");
            println!(" compression ratio must be a positive value! {}", comp_ratio);
            std::process::exit(1);
        }
        self.geometry.compression_ratio = comp_ratio;
        self.geometry.clearance = self.geometry.displacement / (comp_ratio - 1.0);
        self.geometry.total_volume = self.geometry.displacement + self.geometry.clearance;
    }

    pub fn set_store_species(&mut self, state: bool) {
        self.store_species = state;
    }

    pub fn set_combustion_model(&mut self, comb: Box<dyn Combustion>) {
        self.combustion = comb;
    }

    /// Test the compression phase of a cylinder.
    pub fn _test_closed_phase(&mut self) {
        // setting up cylinder at bottle-dead-center
        self.angle = PI;    // BDC
        let (vol, _) = Cylinder::calc_volume(&self.geometry, &self.crankshaft, self.angle);
        self.volume = vol;

        // setting up the run 
        let d_angle:f64 = 0.01; // [CA deg]
        let limit = (180.0 / d_angle) as i32;
        let d_angle = d_angle.to_radians();
        let pvk = self.gas.P()*self.volume.powf(self.gas.k());
        let k = self.gas.k();
        println!("initial PV^K = {}", pvk);
        for _ in 0..limit {
            let new_prop = self.closed_phase(d_angle);
            let temp = new_prop.0;
            let press = new_prop.1;
            let mass = new_prop.2;
            let vol = new_prop.3;

            // update: T, P, V, angle, mass and composition
            self.angle = if self.angle + d_angle >= 4.0 * PI {
                self.angle + d_angle - 4.0 * PI
            } else {
                self.angle + d_angle
            };
            self.gas.TP(temp, press);
            self.mass = mass;
            self.volume = vol;
        }
        println!("max P = {}", pvk/self.volume.powf(k));
    }

    /// Estimates the cylinder volume and its derivative.
    pub fn _test_volume(&mut self) {
        println!("Running '{}' `_test_volume()` ", self.name());
        let d_angle = 0.1f64;
        let mut angle = 0.0f64;
        let mut result: Vec<String> = Vec::new();
        result.push(format!("crank-angle [deg]\tvolume [m³]\tdv [m³/CA-rad]\n"));
        for _ in 0..7200 {
            let (v, dv) =
                Cylinder::calc_volume(&self.geometry, &self.crankshaft, angle.to_radians());
            result.push(format!("{:.2}\t{}\t{}\n", angle, v, dv));
            angle +=d_angle
        }
        let mut file = std::fs::File::create("v_dv_test.txt").expect("Error opening writing file");
        write!(file, "{}", result.join("")).expect("Unable to write data");
        println!("Test fineshed!");
    }
}

impl ZeroDim for Cylinder {
    fn name<'a>(&'a self) -> &'a str {
        &self.name
    }
    fn get_state(&self) -> BasicProperties {
        BasicProperties {
            name: self.name(),
            pressure: self.gas.P(),
            temperature: self.gas.T(),
            cp: self.gas.cp(),
            cv: self.gas.cv(),
            cp_cv: self.gas.k(),
            gas_const: self.gas.R(),
            crank_angle: Some(self.angle),
        }
    }
    fn advance(&mut self, dt: f64) {
        let d_angle = dt * self.sec_to_rad; // [CA radian]
        let new_prop: (f64, f64, f64, f64, Array1<f64>);
        if self.int_valves.flow_info.mass_flow == 0.0 && self.exh_valves.flow_info.mass_flow == 0.0 {
            new_prop = self.closed_phase(d_angle);
        } else {
            new_prop = self.open_phase(d_angle);
        }
        let temp = new_prop.0;
        let press = new_prop.1;
        let mass = new_prop.2;
        let vol = new_prop.3;
        let mole_frac = new_prop.4;

        // update: T, P, V, angle, mass and composition
        self.angle = if self.angle + d_angle >= 4.0 * PI {
            self.angle + d_angle - 4.0 * PI
        } else {
            self.angle + d_angle
        };
        self.gas.TPX_array(temp, press, &mole_frac);
        self.mass = mass;
        self.volume = vol;
    }

    fn update_flow_ratio(&mut self, total_flow_ratio: Vec<(&str, &FlowRatio)>) {

        // find intake valves
        let mut intake_flow_ratio = FlowRatio::new();
        for valve in self.int_valves.basic_info.iter_mut() {
            match total_flow_ratio.iter().find(|(name, _)| **name == valve.name) {
                Some((_,flow_ratio)) => {
                    valve.mass_flow = flow_ratio.mass_flow;
                    intake_flow_ratio = &intake_flow_ratio + flow_ratio;
                },
                None => {
                    println!("Error at Cylinder::update_flow_ratio()");
                    println!(" object {} is not connected to one of the objects:", valve.name);
                    for (name, _) in total_flow_ratio.iter() {
                        print!(" '{}'", name);
                    }
                    std::process::exit(1)
                }
            }
        }

        // find exhaust valves
        let mut exhaust_flow_ratio = FlowRatio::new();
        for valve in self.exh_valves.basic_info.iter_mut() {
            match total_flow_ratio.iter().find(|(name, _)| **name == valve.name) {
                Some((_,flow_ratio)) => {
                    valve.mass_flow = flow_ratio.mass_flow;
                    exhaust_flow_ratio = &exhaust_flow_ratio + flow_ratio;
                },
                None => {
                    println!("Error at Cylinder::update_flow_ratio()");
                    println!(" object {} is not connected to one of the objects:", valve.name);
                    for (name, _) in total_flow_ratio.iter() {
                        print!(" '{}'", name);
                    }
                    std::process::exit(1)
                }
            }
        }
        self.int_valves.flow_info = intake_flow_ratio;
        self.exh_valves.flow_info = exhaust_flow_ratio;
    }
}

impl SaveData for Cylinder {
    fn get_headers(&self) -> String {
        if self.store_species {
            let hearder = "crank-angle [deg]\tpressure [bar]\ttemperature [K]\tvolume [cm³]\tmass [mg]".to_string();
            let species = self.gas.species().join("\t");
            format!("{}\t{}", hearder, species)
        } else {
            "crank-angle [deg]\tpressure [bar]\ttemperature [K]\tvolume [cm³]\tmass [mg]".to_string()
        }
        
    }
    fn num_storable_variables(&self) -> usize {
        let num_prop: usize = 5;
        if self.store_species {
            num_prop + self.gas.species().len()
        } else {num_prop}
    }
    fn get_storable_data(&self) -> Array1<f64> {
        if self.store_species {
            stack![Axis(0), 
            array![self.angle.to_degrees(), self.gas.P()/1e5, self.gas.T(), self.volume*1e6, self.mass*1e6],
            self.gas.mole_frac().clone()
            ]
        } else {
            array![self.angle.to_degrees(), self.gas.P()/1e5, self.gas.T(), self.volume*1e6, self.mass*1e6]
        }
    }
}

impl std::fmt::Display for Cylinder {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}: 
        gas: `{}` 
        angle: {:.2} [CA deg] (ref at TDC firing)
        temperature: {:.2} [K]
        pressure: {:.2} [Pa]
        volume: {:.1} [cm³]
        mass: {:.1} [mg]
        {} \t\t\t {} \t\t\t {}
        diameter: {:.1} [mm] \t\t\t mean velocity: {:.1} [m/s] \t\t temperature: {:.1} [K]
        stroke: {:.1} [mm] \t\t\t temperature: {:.1} [K] \t\t area: {:.1} [cm²]
        displacement: {:.1} [cm³] \t\t area: {:.1} [cm²] \t\t
        compression_ratio: {:.1}
        wall temperature: {:.1} [K]",
            Style::new().bold().paint(&self.name),
            self.gas.name(),
            self.angle.to_degrees(),
            self.gas.T(),
            self.gas.P(),
            self.volume * 1e6,
            self.mass * 1e6,
            Style::new().underline().paint("     Geometry     "),
            Style::new().underline().paint("      Piston      "),
            Style::new().underline().paint("       Head       "),
            self.geometry.diameter * 1e3,
            self.piston.mean_velocity,
            self.head.temperature,
            self.geometry.stroke * 1e3,
            self.piston.temperature,
            self.head.area * 1e4,
            self.geometry.displacement * 1e6,
            self.piston.area * 1e4,
            self.geometry.compression_ratio,
            self.geometry.wall_temp,
        )
    }
}

impl ZeroD for Cylinder {}

#[derive(Debug, Clone)]
struct Geometry {
    compression_ratio: f64, //[-]
    diameter: f64,          //[m]
    stroke: f64,            //[m]
    transverse_area: f64,   //[m²]
    displacement: f64,      //[m³]
    clearance: f64,         //[m³]
    total_volume: f64,      //[m³]
    wall_temp: f64,         //[K]
}

impl Geometry {
    /// Creates a `Geometry` object. Inputs must be in SI units
    fn new(diam: f64, displ: f64, comp_ratio: f64, wall_temp: f64) -> Geometry {
        let transverse_area = 0.25 * PI * diam * diam;
        let stroke = displ / transverse_area;
        let clearance = displ / (comp_ratio - 1.0);
        let total_volume = displ + clearance;
        Geometry {
            compression_ratio: comp_ratio,
            diameter: diam,
            stroke,
            transverse_area,
            displacement: displ,
            clearance,
            total_volume,
            wall_temp,
        }
    }
}

#[derive(Debug, Clone)]
struct Crankshaft {
    conrod: f64,            // [m]
    crank: f64,             // [m]
    eccentricity: f64,      // [m]
    mech_total_length: f64, // [m]
    angle_tdc: f64,         // [CA radian] - should be zero if engine's eccentricity equals to zero
}

impl Crankshaft {
    fn new(conrod: f64, crank: f64, eccentricity: f64, angle_tdc: f64) -> Crankshaft {
        Crankshaft {
            conrod,
            crank,
            eccentricity,
            mech_total_length: (conrod + crank) * angle_tdc.cos(),
            angle_tdc: angle_tdc.to_radians(),
        }
    }
}

#[derive(Debug, Clone)]
struct Piston {
    mean_velocity: f64, // [m/s]
    temperature: f64,   // [K]
    area: f64,          // [m^2]
}

impl Piston {
    fn new(stroke: f64, diameter: f64, temperature: f64, speed: f64) -> Piston {
        let mean_velocity = 2.0 * stroke * speed / 60.0;
        let area = 0.25 * PI * diameter * diameter;
        Piston {
            mean_velocity,
            temperature,
            area,
        }
    }
}

#[derive(Debug, Clone)]
struct Head {
    temperature: f64, // [K]
    area: f64,        // [m^2]
}

#[derive(Debug, Clone)]
struct HeatTransfer {}

impl HeatTransfer {
    fn calculate(&self, vol: f64, temp: f64, press: f64, cyl: &Cylinder) -> f64 {
        let heat_trans_coeff = 130.0
            * vol.powf(-0.06)
            * (press * 1e-5).powf(0.8)
            * temp.powf(-0.4)
            * (cyl.piston.mean_velocity + 1.4).powf(0.8);
        let area_sup = PI * cyl.geometry.diameter * (vol / cyl.geometry.transverse_area);
        let q_piston = heat_trans_coeff * cyl.piston.area * (cyl.piston.temperature - temp);
        let q_head = heat_trans_coeff * cyl.head.area * (cyl.head.temperature - temp);
        let q_cylinder = heat_trans_coeff * area_sup * (cyl.geometry.wall_temp - temp);
        q_piston + q_head + q_cylinder
    }
}

#[derive(Debug, Clone)]
struct ValvesInfo {
    basic_info: Vec<ValveBasicInfo>,
    flow_info: FlowRatio,
    gas_comp: String,
    opening: f64,
    closing: f64,
}

impl ValvesInfo {
    fn new(basic_info: Vec<ValveBasicInfo>, gas_comp: String) -> ValvesInfo {
        let mut opens_at = std::f64::INFINITY;
        let mut closes_at = std::f64::NEG_INFINITY;
        for v in basic_info.iter() {
            if v.opening_angle < opens_at {
                opens_at = v.opening_angle;
            }
            if v.closing_angle > closes_at {
                closes_at = v.closing_angle;
            }
        }
        ValvesInfo {
            basic_info,
            flow_info: FlowRatio::new(),
            gas_comp,
            opening: opens_at.to_radians(),
            closing: closes_at.to_radians(),
        }
    }
}

#[derive(Debug, Clone)]
struct ValveBasicInfo {
    name: String,
    opening_angle: f64,
    closing_angle: f64,
    mass_flow: f64,
    back_flow: f64,
}

impl ValveBasicInfo {
    fn new(name: String, opening_angle: f64, closing_angle: f64) -> ValveBasicInfo {
        ValveBasicInfo {
            name,
            opening_angle,
            closing_angle,
            mass_flow: 0.0,
            back_flow: 0.0,
        }
    }
    /// Returns the fresh charge through the valve, always positive sign.
    /// If flow is leaving the cylinder, it is added to the backflow counter. 
    fn get_charge(&mut self, dt: f64) -> f64 {
        let intake_mass = self.mass_flow*dt;
        if intake_mass >= 0.0 {
            let fresh_charge = intake_mass - self.back_flow;
            if fresh_charge >= 0.0 {
                self.back_flow = 0.0;
                return fresh_charge;
            } else {
                self.back_flow = -fresh_charge;
                return 0.0;
            }
        } else {
            self.back_flow -= intake_mass;
            0.0
        }
    }
}

