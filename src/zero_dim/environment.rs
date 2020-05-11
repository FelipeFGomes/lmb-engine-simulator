use crate::reaction::gas::Gas;
use crate::zero_dim::zero_core::{ZeroDim};
use crate::{BasicProperties, FlowRatio};

pub struct Environment {
    name: String,
    gas: Gas,
    _mass: f64,
}

impl Environment {
    pub fn new(name: String, gas: &Gas) -> Result<Environment, &'static str> {
        Ok(Environment {
            name,
            gas: gas.clone(),
            _mass: std::f64::INFINITY,
        })
    } 
}

impl ZeroDim for Environment { 
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
    fn advance(&mut self, _: f64) {}
    fn update_flow_ratio(&mut self, _: FlowRatio) {}
}