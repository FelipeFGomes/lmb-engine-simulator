use crate::reaction::gas::Gas;
use crate::zero_dim::cylinder::Cylinder;
use crate::zero_dim::zero_core::ZeroDim;
use crate::connector::valve::Valve;
use std::f64::consts::PI;
use serde::{Deserialize, Serialize};
use serde_json;

#[derive(Debug)]
pub struct Engine {
    speed: f64,
    eccentricity: f64,
    firing_order: String,
    sec_to_rad: f64,
}

type EngineOutput = Result<(Engine, Vec<Cylinder>, Vec<Valve>), String>;

impl Engine {
    pub fn new(file_name: &str, gas: &Gas) -> EngineOutput {
        let json_engine = match Engine::reading_json(file_name) {
            Ok(i) => i,
            Err(err) => {
                let msg = format!("Error at Engine::new \n unable to parse {} file \n{}", file_name, err);
                return Err(msg);
            }
        };

        let firing_order: Vec<f64> = json_engine.firing_order.split('-').map(|s| -> f64 {
                match s.parse() {
                    Ok(s) => s,
                    Err(err) => {
                        eprintln!("Error at Engine:new \n unable to parse {} \n{}", s, err);
                        std::process::exit(1);
                    }
                }
            }).collect();

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

        // instantianting cylinders
        let mut cylinders: Vec<Cylinder> = Vec::new();
        for (cylinder, order) in json_engine.cylinders.iter().zip(firing_order.iter()) {
            cylinders.push(Cylinder::new(
                cylinder.name.clone(),
                json_engine.cylinders_diameter,
                json_engine.displacement,
                cylinder.compression_ratio,
                json_engine.conrod,
                json_engine.speed,
                180.0 + (order - 1.0) * division,
                json_engine.eccentricity,
                gas,
            )?);
        }

        // instantianting valves
        let mut valves: Vec<Valve> = Vec::new();
        for valve in json_engine.valves.iter() {
            let cyl = match cylinders.iter().find(|c| c.name() == valve.connecting) {
                Some(c) => c,
                None => { 
                    let msg = format!("{} is not connected to any cylinder", valve.name);
                    return Err(msg); 
                },
            };
            valves.push(Valve::new(
                valve.name.clone(),
                valve.opening_angle,
                valve.closing_angle,
                valve.diameter*1e-3,
                valve.max_lift*1e-3,
                cyl,
            )?);
        }
        let engine = Engine {
            speed: json_engine.speed,
            eccentricity: json_engine.eccentricity,
            firing_order: json_engine.firing_order.clone(),
            sec_to_rad: 2.0 * PI * json_engine.speed / 60.0,
        };
        Ok( (engine, cylinders, valves) )
    }

    pub fn sec_to_rad(&self) -> f64 {
        self.sec_to_rad
    } 

    fn reading_json(file_name: &str) -> serde_json::Result<JsonEngine> {
        let json_file = std::fs::read_to_string(file_name).expect("Unable to read file");
        let data: JsonEngine = serde_json::from_str(&json_file)?;
        Ok(data)
    }
}

#[derive(Serialize, Deserialize, Debug)]
struct JsonEngine {
    speed: f64,
    eccentricity: f64,
    conrod: f64,             // [mm]
    displacement: f64,       // [cm³]
    cylinders_diameter: f64, // [mm]
    firing_order: String,
    cylinders: Vec<JsonCylinder>,
    valves: Vec<JsonValve>,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct JsonCylinder {
    name: String,
    compression_ratio: f64, // [-]
}

#[derive(Serialize, Deserialize, Debug)]
pub struct JsonValve {
    name: String,
    valve_type: String,
    opening_angle: f64,
    closing_angle: f64,
    diameter: f64,  // [mm]
    max_lift: f64,  // [mm]
    connecting: String,
}
