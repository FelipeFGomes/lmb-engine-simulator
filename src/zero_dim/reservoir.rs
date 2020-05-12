use crate::reaction::gas::Gas;
use crate::zero_dim::zero_core::{ZeroDim};
use crate::{BasicProperties, FlowRatio, StoreData};
use crate::numerics::ode_solvers as ode;
use ndarray::*;

pub struct Reservoir {
    name: String,
    gas: Gas,
    volume: f64,
    mass: f64,
    flow_ratio: FlowRatio,
    store_data: StoreData,
}

impl Reservoir {
    /// `volume` in m³
    pub fn new(name: String, gas: &Gas, volume: f64) -> Result<Reservoir, &'static str> {
        if volume <= 0.0 {
            return Err("`volume` must be greater than zero");
        }
        let header = "pressure [bar]\ttemperature [K]\tmass [kg]";
        Ok( Reservoir {
            name,
            gas: gas.clone(),
            volume,
            mass: gas.P()*volume/(gas.R()*gas.T()),
            flow_ratio: FlowRatio::new(),
            store_data: StoreData::new(header, 3),
        } )
    }
}


impl ZeroDim for Reservoir { 
    fn name<'a>(&'a self) -> &'a str {&self.name}
    fn get_state(&self) -> BasicProperties {
        BasicProperties {
            name: self.name(),
            pressure: self.gas.P(),
            temperature: self.gas.T(),
            cp: self.gas.cp(),
            cv: self.gas.cv(),
            cp_cv: self.gas.k(),
            gas_const: self.gas.R(),
            crank_angle: None,
        }
    }
    fn advance(&mut self, dt: f64) {
        let cv = self.gas.R() / (self.gas.k() - 1.0);
        let cv_inv = 1.0 / cv;
        let system_equations = |_: &f64, x: &Array1<f64>, _: &Vec<f64>| -> Array1<f64> {
            // x[0] = temperature, x[1] = mass
            let d_mass = self.flow_ratio.mass_flow;
            let d_temp = cv_inv / x[1] * (self.flow_ratio.enthalpy_flow - cv * x[0] * d_mass); // [K/CA radian]
            array![d_temp, d_mass]
        };

        // Runge-Kutta 4th order solution
        let ini_condition = array![self.gas.T(), self.mass];
        let integrated = ode::rk4_step(
            system_equations,
            &ini_condition,
            &Vec::new(),
            &0.0,
            dt,
        );
        let temp = integrated[0];
        let mass = integrated[1];
        
        // calculating new pressure
        let press = mass*self.gas.R()*temp/self.volume;
        self.gas.TP(temp, press);
        self.mass = mass;
        
        // storing data
        self.store_data.add_data(array![self.gas.P()/1e5, self.gas.T(), self.mass]);
    }
    fn update_flow_ratio(&mut self, total_flow_ratio: FlowRatio) {
        self.flow_ratio = total_flow_ratio;
    }
    fn write_to_file(
        &self,
        file_name: &str,
        index_range: (usize, usize),
        additional_data: Option<(String, ArrayView2<f64>)>,
    ) 
    {
        self.store_data.write_to_file(file_name, index_range, additional_data);
    }
}