//! # lmb_engine_simulator
//!
//! The `lmb_engine_simulator` crate provides an easy way to simulate engines and 1D gas dynamics.

use std::ops::Add;

mod core;
mod reaction;
mod base;
mod engine;
mod zero_dim;
mod connector;
mod one_dim;
mod numerics;

// Re-exporting
pub use crate::reaction::gas::Gas;
pub use crate::core::system_builder::SystemBuilder;
pub use crate::numerics::ode_solvers;

// Object's type. 'String' meant to store the name of the object
#[derive(Debug, PartialEq)]
pub enum ObjectType {
    ZeroDim,
    OneDim,
}

#[derive(Debug)]
pub struct BasicProperties<'a> {
    pub name: &'a str,
    pub pressure: f64,      // Pa
    pub temperature: f64,   // K
    pub cp: f64,            // J/(kg.K)
    pub cv: f64,            // J/(kg.K)
    pub cp_cv: f64,
    pub gas_const: f64,     // J/(kg.K)
    pub crank_angle: Option<f64>,   // CA rad
}

impl<'a> std::fmt::Display for BasicProperties<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}: 
        pressure: {} [Pa] 
        temperature: {} [K]
        cp: {} [J/(kg.K)]
        cv: {} [J/(kg.K)]
        cp/cv: {}
        R: {} [J/(kg.K)]", 
        self.name, 
        self.pressure,
        self.temperature, 
        self.cp, 
        self.cv, 
        self.cp_cv, 
        self.gas_const)
    }
}

#[derive(Debug)]
pub struct FlowRatio {
    pub mass_flow: f64,     // kg/s
    pub enthalpy_flow: f64, // J/s 
}

impl FlowRatio {
    pub fn new() -> FlowRatio {
        FlowRatio {
            mass_flow: 0.0,
            enthalpy_flow: 0.0,
        }
    }
}

impl Add for FlowRatio {
    type Output = FlowRatio;
    fn add(self, other: FlowRatio) -> FlowRatio {
        FlowRatio {
            mass_flow: self.mass_flow + other.mass_flow,
            enthalpy_flow: self.enthalpy_flow + other.enthalpy_flow,
        }
    }
}

impl<'a, 'b> Add<&'b FlowRatio> for &'a FlowRatio {
    type Output = FlowRatio;
    fn add(self, other: &'b FlowRatio) -> FlowRatio {
        FlowRatio {
            mass_flow: self.mass_flow + other.mass_flow,
            enthalpy_flow: self.enthalpy_flow + other.enthalpy_flow,
        }
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
