use crate::reaction::gas::Gas;
use crate::zero_dim::zero_core::{ZeroDim};
use crate::{BasicProperties, FlowRatio};

pub struct Reservoir {
    name: String,
    gas: Gas,
    volume: f64,
    mass: f64,
    flow_ratio: FlowRatio,
}

impl Reservoir {
    pub fn new(name: String, gas: &Gas, volume: f64) -> Result<Reservoir, &'static str> {
        if volume <= 0.0 {
            return Err("`volume` must be greater than zero");
        }

        Ok( Reservoir {
            name,
            gas: gas.clone(),
            volume,
            mass: gas.P()/(volume*gas.R()*gas.T()),
            flow_ratio: FlowRatio::new(),
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
        let h_dot = self.flow_ratio.enthalpy_flow;
        let m_dot = self.flow_ratio.mass_flow;
        
        // Solution: Euler method
        let dtemp_dt = 1.0/(self.mass*self.gas.cv())*(h_dot - self.gas.e()*m_dot);
        let mass = self.mass + dt*m_dot;
        let temp = self.gas.T() + dt*dtemp_dt;
        
        // calculating new pressure
        let press = mass*self.gas.R()*temp*self.volume;
        self.gas.TP(temp, press);
        self.mass = mass;
        // println!("Advanced '{}'\n new state:\n {:?}", self.name(), self.get_state());
    }
    fn update_flow_ratio(&mut self, total_flow_ratio: FlowRatio) {
        self.flow_ratio = total_flow_ratio;
    }

    fn _get_main_properties(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\n",
            self.gas.T(),
            self.gas.P(),
            self.volume,
            self.mass,
        )
    }
    
}