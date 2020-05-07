use crate::numerics::ode_solvers as ode;
use crate::reaction::gas::Gas;
use crate::zero_dim::zero_core::ZeroDim;
use crate::{BasicProperties, FlowRatio};
use ansi_term::Style;
use ndarray::*;
use std::f64::consts::PI;
use std::io::Write;

#[derive(Debug)]
pub struct Cylinder {
    name: String,
    gas: Gas,
    mass: f64,       // [kg]
    volume: f64,     // [m³] - instant volume as function of crank angle
    angle: f64,      // [CA deg] - instant crank angle
    speed: f64,      // [RPS]
    sec_to_rad: f64, // constant: 2*PI*speed
    geometry: Geometry,
    crankshaft: Crankshaft,
    piston: Piston,
    head: Head,
    heat_transfer: HeatTransfer,
    combustion: Combustion,
    flow_ratio: FlowRatio,
    // blow_by: bool,
    // crevice: bool,
}

impl Cylinder {
    /// Creates a cylinder object. Inputs units must be: `mm`, `RPM` and `CA deg`.
    pub fn new(
        name: String,
        diam: f64,
        displ: f64,
        comp_ratio: f64,
        conrod: f64,
        speed: f64,
        ini_angle: f64,
        eccentricity: f64,
        gas: &Gas,
    ) -> Result<Cylinder, String> {
        if diam < 0.0 {
            return Err(format!("diameter cannot be lower than zero"));
        } else if displ < 0.0 {
            return Err(format!("displacement cannot be lower than zero"));
        } else if comp_ratio < 0.0 {
            return Err(format!("compration ratio cannot be lower than zero"));
        } else if speed < 0.0 {
            return Err(format!("engine speed cannot be lower than zero"));
        }

        // In SI uinits
        let conrod = conrod * 1e-3; // [m]
        let eccentricity = eccentricity * 1e-3; // [m]
        let diam = diam * 1e-3; // [m]
        let wall_temp = 520.0; // [K] - Common value
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
        Ok(Cylinder {
            name,
            gas: gas.clone(),
            mass: gas.P() * volume / (gas.R() * gas.T()),
            volume,
            angle: (ini_angle.to_radians() - crankshaft.angle_tdc),
            speed: speed / 60.0,
            sec_to_rad: 2.0 * PI * speed / 60.0,
            geometry,
            crankshaft,
            piston,
            head,
            heat_transfer: HeatTransfer {},
            combustion: Combustion {},
            flow_ratio: FlowRatio::new(),
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
    fn closed_phase(&mut self, d_angle: f64) -> (f64, f64, f64, f64 ) {
        // Closed Phase -----------------------------------------------------------
        let heat_combustion = self.combustion.calculate();
        let const_1 = (self.gas.k() - 1.0) / (self.mass * self.gas.R());
        
        let closed_phase_equations = |angle: &f64, x: &Array1<f64>, _: &Vec<f64>| -> Array1<f64> {
            // x[0] = P, x[1] = T
            let (vol, d_vol) = Cylinder::calc_volume(&self.geometry, &self.crankshaft, *angle);    
            let heat_transfer = self.heat_transfer.calculate(vol, x[1], x[0], &self); // [J/s]
            let heat_transfer = heat_transfer / self.sec_to_rad; // [J/CA radian]
            let d_temp = const_1 * (heat_combustion + heat_transfer - x[0] * d_vol); // [K/CA radian]
            let d_press = x[0] * (d_temp / x[1] - d_vol / vol); // [Pa/CA radian]
            array![d_press, d_temp]
        };

        let ini = array![self.gas.P(), self.gas.T()];
        let closed_phase_integrated = ode::rk4_step(
            closed_phase_equations,
            &ini,
            &Vec::new(),
            &self.angle,
            d_angle,
        );

        let press = closed_phase_integrated[0];
        let temp = closed_phase_integrated[1];        
        let (vol, _) = Cylinder::calc_volume(&self.geometry, &self.crankshaft, self.angle + d_angle);
        ( temp, press, self.mass, vol )
    }

    /// `d_angle` in crank angle radians
    fn open_phase(&mut self, d_angle: f64) -> (f64, f64, f64, f64 ) {
        // Open Phase -----------------------------------------------------------
        let cv = self.gas.R() / (self.gas.k() - 1.0);
        let cv_inv = 1.0 / cv;
        let open_phase_equations = |angle: &f64, x: &Array1<f64>, _: &Vec<f64>| -> Array1<f64> {
            //x[0] = temperature, x[1] = mass
            let (vol, d_vol) = Cylinder::calc_volume(&self.geometry, &self.crankshaft, *angle);
            let press = x[1] * self.gas.R() * x[0] / vol;
            let heat_transfer = self.heat_transfer.calculate(vol, x[0], press, &self); // [J/s]
            let heat_transfer = heat_transfer / self.sec_to_rad; // [J/CA radian]
            let d_mass = self.flow_ratio.mass_flow / self.sec_to_rad; // [kg/CA radian]
            let d_temp = cv_inv / x[1]
                * (heat_transfer - press * d_vol + self.flow_ratio.enthalpy_flow / self.sec_to_rad
                    - cv * x[0] * d_mass); // [K/CA radian]
            array![d_temp, d_mass]
        };

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
        ( temp, press, mass, vol )
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
        let mut result: Vec<String> = Vec::new();
        let d_angle = d_angle.to_radians();
        let pvk = self.gas.P()*self.volume.powf(self.gas.k());
        let k = self.gas.k();
        println!("initial PV^K = {}", pvk);
        result.push(format!("crank-angle\ttemperature\tpressure\tvolume\tmass\n"));
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

            result.push(self._get_main_properties());
        }
        let mut file = std::fs::File::create(self.name()).expect("Error opening writing file");
        write!(file, "{}", result.join("")).expect("Unable to write data");
        println!("max P = {}", pvk/self.volume.powf(k));
    }

    pub fn _test_volume(&mut self) {
        let range: Array1<f64> = Array1::linspace(0.0, 720.0, 7200);
        let mut result: Vec<String> = Vec::new();
        result.push(format!("crank-angle [deg]\tvolume [m³]\tdv [m³/CA-rad]\n"));
        for angle in range.iter() {
            let (v, dv) =
                Cylinder::calc_volume(&self.geometry, &self.crankshaft, angle.to_radians());
            result.push(format!("{}\t{}\t{}\n", angle, v, dv));
        }
        let mut file = std::fs::File::create("v_dv").expect("Error opening writing file");
        write!(file, "{}", result.join("")).expect("Unable to write data");
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
        let new_prop: (f64, f64, f64, f64);
        if self.flow_ratio.mass_flow == 0.0 {
            new_prop = self.closed_phase(d_angle);
        } else {
            new_prop = self.open_phase(d_angle);
        }
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

    fn update_flow_ratio(&mut self, total_flow_ratio: FlowRatio) {
        self.flow_ratio = total_flow_ratio;
    }

    fn _get_main_properties(&self) -> String {
        format!(
            "{:.2}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\n",
            self.angle.to_degrees(),
            self.gas.T(),
            self.gas.P(),
            self.volume*1e6,
            self.mass,
        )
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

#[derive(Debug)]
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

#[derive(Debug)]
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

#[derive(Debug)]
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

#[derive(Debug)]
struct Head {
    temperature: f64, // [K]
    area: f64,        // [m^2]
}

#[derive(Debug)]
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

#[derive(Debug)]
struct Combustion {}

impl Combustion {
    fn calculate(&self) -> f64 {
        0.0
    }
}
