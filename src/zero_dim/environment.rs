use crate::core::traits::{SaveData, ZeroD, ZeroDim};
use crate::reaction::gas::Gas;
use crate::{BasicProperties, FlowRatio};
use ndarray::*;

/// Zero-Dimensional struct with constant pressure, temperature and composition.
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
            crank_angle: None,
        }
    }
    fn advance(&mut self, _: f64) {}
    fn update_flow_ratio(&mut self, _: Vec<(&str, &FlowRatio)>) {}
}

impl SaveData for Environment {
    fn get_headers(&self) -> String {
        "pressure [bar]\ttemperature [K]".to_string()
    }
    fn num_storable_variables(&self) -> usize {
        2
    }
    fn get_storable_data(&self) -> Array1<f64> {
        array![self.gas.P() / 1e5, self.gas.T()]
    }
}

impl ZeroD for Environment {}
