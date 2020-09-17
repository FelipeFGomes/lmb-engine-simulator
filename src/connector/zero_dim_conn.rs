#![allow(non_snake_case)]

use crate::core::traits::{Connector, SaveData, Conn};
use crate::{BasicProperties, FlowRatio};
use ndarray::*;

pub struct ZeroDimConn {
    name: String,
    area: f64,
    discharge_coeff: f64,
    connecting: Vec<String>,
    flow_ratio: Vec<FlowRatio>,
}

impl ZeroDimConn {
    pub fn new(name: &str, diam: f64, discharge_coeff: f64, connecting: Vec<String>) -> Result<ZeroDimConn, String> {
        if connecting.len() != 2 {
            return Err(format!("Connector `ZeroDimConn` must connect only two elements, connecting: {}", connecting.len()));
        }
        if discharge_coeff > 1.0 || discharge_coeff <= 0.0 {
            return Err(format!("`discharge_coeff` must be between 0.0 and 1.0: {}", discharge_coeff));
        }
        let flow_ratio = vec![FlowRatio::new(), FlowRatio::new()];
        Ok(ZeroDimConn {
            name: name.to_string(),
            area: 0.25*std::f64::consts::PI*diam*diam,
            discharge_coeff,
            connecting,
            flow_ratio,
        })
    }
}

impl Connector for ZeroDimConn {
    fn name<'a>(&'a self) -> &'a str {&self.name}
    fn connecting<'a>(&'a self) -> &'a Vec<String> {&self.connecting}
    fn connect_to(&mut self, elem_name: &str) -> Result<(),String> {
        self.connecting.push(elem_name.to_string());
        self.flow_ratio.push( FlowRatio::new() );
        if self.connecting.len() != 2 {
            return Err(
                format!("Wrong the number of connections. ZeroDimConn should connect only two elements"),
            );
        }
        Ok(())
    }
    fn update_flow_ratio(&mut self, prop: Vec<BasicProperties>, _step: f64) {
        // checking flow diretion
        let i_up: usize;
        let i_down: usize;
        let press_ratio = prop[0].pressure / prop[1].pressure;
        if press_ratio > 1.00000001 {
            i_up = 0;
            i_down = 1;
        } else if press_ratio < 0.99999991 {
            i_up = 1;
            i_down = 0;
        } else {
            self.flow_ratio.iter_mut().for_each(|f| *f = FlowRatio::new());
            return;
        }

        let P_up = prop[i_up].pressure;
        let T_up = prop[i_up].temperature;
        let k = prop[i_up].cp_cv;
        let R = prop[i_up].gas_const;
        let P_down = prop[i_down].pressure;

        let kp = k + 1.0;
        let km = k - 1.0;
        let P_du = P_down / P_up;
        let m_dot: f64;
        // estimating mass flow
        if P_du > (2.0 / kp).powf(k / km) {
            m_dot = self.discharge_coeff * self.area * P_up / (R * T_up).sqrt()
                * (2.0 * k / km * (P_du.powf(2.0 / k) - P_du.powf(kp / k))).sqrt();
        } else {
            // chocked flow: independent of downstream pressure
            m_dot = self.discharge_coeff * self.area * P_up / (R * T_up).sqrt() * (k * (2.0 / kp).powf(kp / km)).sqrt();
        }

        // updating `flow_ratio` for upstream objects
        let i = match self.connecting.iter().position(|name| *name == prop[i_up].name) {
            Some(i) => i,
            None => {
                println!("Error at ZeroDimConn::update_flow_ratio:");
                println!(" Objects:");
                for obj in prop.iter() {
                    print!(" '{}'", obj.name);
                }
                println!("are not connected to '{}'", self.name());
                std::process::exit(1)
            }
        };
        let ii = i; // store the position to use in the downstream
        self.flow_ratio[i].mass_flow = -m_dot;
        self.flow_ratio[i].enthalpy_flow = -m_dot*k*R/km*prop[i_up].temperature;

        // updating `flow_ratio` for downstream objects
        let i = match self.connecting.iter().position(|name| *name == prop[i_down].name) {
            Some(i) => i,
            None => {
                println!("Error at ZeroDimConn::update_flow_ratio():");
                println!(" Objects:");
                for obj in prop.iter() {
                    print!(" '{}'", obj.name);
                }
                println!("are not connected to '{}'", self.name());
                std::process::exit(1)
            }
        };
        self.flow_ratio[i].mass_flow = -self.flow_ratio[ii].mass_flow;
        self.flow_ratio[i].enthalpy_flow = -self.flow_ratio[ii].enthalpy_flow;
    }
    fn get_flow_ratio<'a>(&'a self, elem_name: &str) -> Result<&'a FlowRatio, String> {
        match self.connecting.iter().position(|name| name == elem_name) {
            Some(i) => Ok(&self.flow_ratio[i]),
            None => Err(format!("object '{}' was not found in '{}'", elem_name, self.name())),
        }
    }
}

impl SaveData for ZeroDimConn {
    fn get_headers(&self) -> String {"mass flow [kg/s]\tenthalpy flow [J/s]".to_string()}
    fn num_storable_variables(&self) -> usize {2}
    fn get_storable_data(&self) -> Array1<f64> {
        array![self.flow_ratio[0].mass_flow, self.flow_ratio[0].enthalpy_flow]
    }
}

impl Conn for ZeroDimConn {}