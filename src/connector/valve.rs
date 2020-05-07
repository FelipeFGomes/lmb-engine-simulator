#![allow(non_snake_case)]

use super::conn_core::Connector;
use crate::{BasicProperties, FlowRatio};
use crate::zero_dim::cylinder::Cylinder;
use crate::zero_dim::zero_core::ZeroDim;
use std::io::Write;
// use crate::ObjectType;

#[derive(Debug)]
pub struct Valve {
    name: String,
    opening_angle: f64,
    closing_angle: f64,
    delta_angle: f64,
    diameter: f64,
    max_lift: f64,
    valve_lift: ValveLift,
    flow_ratio: Vec<(String, FlowRatio)>,
    connecting: Vec<String>,
}

impl Valve {
    pub fn new(
        name: String,
        opening_angle: f64,
        closing_angle: f64,
        diameter: f64,
        max_lift: f64,
        cylinder: &Cylinder,
    ) -> Result<Valve, &'static str> {
        let delta_angle = closing_angle - opening_angle;
        let max_lift_diam_ratio = max_lift / diameter;
        let time_opened = if delta_angle < 0.0 {
            delta_angle + 720.0
        } else {
            delta_angle
        };
        let mut flow_ratio:Vec<(String, FlowRatio)> = Vec::new();
        let mut connecting: Vec<String> = Vec::new();
        connecting.push( cylinder.name().to_string() );
        flow_ratio.push( (cylinder.name().to_string(), FlowRatio::new()) );
        let valve_lift = ValveLift::new(max_lift_diam_ratio, time_opened);
        
        Ok(Valve {
            name,
            opening_angle,
            closing_angle,
            delta_angle,
            diameter,
            max_lift,
            valve_lift,
            flow_ratio,
            connecting,
        })
    }
    
    fn calc_throat_area(&self, angle: f64) -> f64 {
        let lift_diam = self.valve_lift.calc_lift(angle);
        let lift = lift_diam * self.diameter;
        if lift_diam <= 0.125 {
            2.22144146908 * lift * (self.diameter + 0.5 * lift)
        } else if lift_diam <= 0.274049 {
            3.33794219444
                * self.diameter
                * (lift * lift - 0.125 * lift * self.diameter
                    + 0.0078125 * self.diameter * self.diameter)
                    .sqrt()
        } else {
            0.73631077818 * self.diameter * self.diameter
        }
    }

    fn is_open(&self, chank_angle: f64) -> bool {
        // checking if valve is open
        if self.delta_angle < 0.0 {
            if chank_angle < self.opening_angle && chank_angle > self.closing_angle {
                false
            } else {
                true
            }
        } else {
            if chank_angle < self.opening_angle || chank_angle > self.closing_angle {
                false
            }
            else {
                true
            }
        }
    }

    fn set_flow_to_zero(&mut self) {
        self.flow_ratio.iter_mut().for_each(|(_, flow)| *flow = FlowRatio::new() );
    }

    pub fn _test_valve_area_and_lift(&self) {
        let file_name = format!("valve_test_{}", self.name());
        let mut file = std::fs::File::create(file_name).expect("Error opening writing file");
        let mut result: Vec<String> = Vec::new();
        result.push(format!("crank-angle [deg]\tarea[mm²]\tlift[mm]\n"));

        for angle in 0..720 {
            let area: f64;
            let lift: f64;
            if self.is_open(angle as f64) {
                let delta = angle as f64 - self.opening_angle;
                let angle_cam = if delta >= 0.0 { delta } else { delta + 720.0 };
                if angle == 340 {
                    let _a=0;
                }
                area = self.calc_throat_area(angle_cam);
                lift = self.valve_lift.calc_lift(angle_cam)*self.diameter;
            } else {
                area = 0.0;
                lift = 0.0;
            }
            result.push(format!("{}\t{}\t{}\n", angle, area*1e6, lift*1e3));
        } 
        write!(file, "{}", result.join("")).expect("Unable to write data");
    }
}

impl Connector for Valve {
    fn name<'a>(&'a self) -> &'a str {&self.name}
    fn connecting<'a>(&'a self) -> &'a Vec<String> {&self.connecting}
    fn connect_to(&mut self, elem_name: &str) -> Result<(),String> {
        self.connecting.push(elem_name.to_string());
        self.flow_ratio.push( (elem_name.to_string(), FlowRatio::new()) );
        if self.connecting.len() != 2 {
            return Err(
                format!("Wrong the number of connections. Valve should connect only two elements"),
            );
        }
        Ok(())
    }

    fn update_flow_ratio(&mut self, prop: Vec<BasicProperties>, _: f64) {

        if prop.len() != 2 {
            println!("Error at Valve::update_flow_ratio:");
            println!(" need BasicProperties from two objects: getting {}", prop.len());
            println!(" objects are: ");
            for obj in prop.iter() {
                    print!(" '{}'", obj.name);
            }
            std::process::exit(1)
        }

        let mut crank_angle = std::f64::NAN; 
        for info in prop.iter() {
            match info.crank_angle {
                Some(angle) => { crank_angle = angle },
                None => {}
            }
        }
        if crank_angle.is_nan() {
            println!("'{}' did not get crank-angle", self.name());
            std::process::exit(1)
        }

        let crank_angle = crank_angle.to_degrees();
        let thoat_area: f64; 
        // checking if valve is open
        if self.is_open(crank_angle) {
            let delta = crank_angle - self.opening_angle;
            let angle = if delta >= 0.0 { delta } else { delta + 720.0 };
            thoat_area = self.calc_throat_area(angle);
        } else {
            self.set_flow_to_zero();
            return
        }
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
            self.set_flow_to_zero();
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
            m_dot = thoat_area * P_up / (R * T_up).sqrt()
                * (2.0 * k / km * (P_du.powf(2.0 / k) - P_du.powf(kp / k))).sqrt();
        } else {
            // chocked flow: independent of downstream pressure
            m_dot = thoat_area * P_up / (R * T_up).sqrt() * (k * (2.0 / kp).powf(kp / km)).sqrt();
        }

        // updating `flow_ratio` for upstream objects
        let i = match self.flow_ratio.iter().position(|(name,_)| *name == prop[i_up].name) {
            Some(i) => i,
            None => {
                println!("Error at Valve::update_flow_ratio:");
                println!(" Objects:");
                for obj in prop.iter() {
                    print!(" '{}'", obj.name);
                }
                println!("are not connected to '{}'", self.name());
                std::process::exit(1)
            }
        };
        let ii = i; // store the position to use in the downstream
        self.flow_ratio[i].1.mass_flow = -m_dot;
        self.flow_ratio[i].1.enthalpy_flow = -m_dot*k*R/km*prop[i_up].temperature;

        // updating `flow_ratio` for downstream objects
        let i = match self.flow_ratio.iter().position(|(name,_)| *name == prop[i_down].name) {
            Some(i) => i,
            None => {
                println!("Error at Valve::update_flow_ratio:");
                println!(" Objects:");
                for obj in prop.iter() {
                    print!(" '{}'", obj.name);
                }
                println!("are not connected to '{}'", self.name());
                std::process::exit(1)
            }
        };
        self.flow_ratio[i].1.mass_flow = -self.flow_ratio[ii].1.mass_flow;
        self.flow_ratio[i].1.enthalpy_flow = -self.flow_ratio[ii].1.enthalpy_flow;

        // let msg = format!("`{}`:\n {:?}", self.name(), self.flow_ratio);
        // println!("{}", msg);
        // println!("P_up = {}, P_down = {}", P_up, P_down);
        
    }

    fn get_flow_ratio<'a>(&'a self, elem_name: &str) -> Result<&'a FlowRatio, String> {
        match self.flow_ratio.iter().find(|(obj_name, _)| obj_name == elem_name) {
            Some((_, data)) => Ok(data),
            None => Err(format!("object '{}' was not found in '{}'", elem_name, self.name())),
        }
    }
}

#[derive(Debug)]
struct ValveLift {
    time_opened: f64, //in crank angle degree
    acceler_ratio: f64,
    max_lift_diam_ratio: f64,
    consts: Vec<f64>,
    intervals: Vec<f64>,
}

impl ValveLift {
    fn new(max_lift_diam_ratio: f64, time_opened: f64) -> ValveLift {
        let acceler_ratio = -2.0;
        let a1 = 4.0 * max_lift_diam_ratio * (1.0 - acceler_ratio) / time_opened.powi(2);
        let b1 = a1 / acceler_ratio;
        let b2 = -b1 * time_opened;
        let b3 = -b2 * time_opened * (1.0 / (4.0 * (1.0 - acceler_ratio)));
        let c1 = a1;
        let c2 = -2.0 * a1 * time_opened;
        let c3 = a1 * time_opened * time_opened;
        let n = 2.0 * (1.0 - acceler_ratio);
        let intervals = vec![time_opened / n, (n - 1.0) * time_opened / n];

        ValveLift {
            time_opened,
            acceler_ratio,
            max_lift_diam_ratio,
            consts: vec![a1, b1, b2, b3, c1, c2, c3],
            intervals,
        }
    }

    /// Calculate valve lift over diameter, `L/D`. `angle` must be in the relative angle in degrees: `theta - opening_angle`
    fn calc_lift(&self, angle: f64) -> f64 {
        if angle <= self.intervals[0] {
            self.consts[0] * angle * angle
        } else if angle >= self.intervals[1] {
            self.consts[4] * angle * angle + self.consts[5] * angle + self.consts[6]
        } else {
            self.consts[1] * angle * angle + self.consts[2] * angle + self.consts[3]
        }
    }
}


