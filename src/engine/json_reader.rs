//! # json_reader
//! 
//! Set of structs designed to read an engine system from an .json file
//! 
//! **Attention when entering the variables in crank-angle degree!** 
//! The reference, where crank-angle is zero, is at top-dead-center (TDC) of compression phase and it only accepts positive numbers.
//! Therefore, the full cycle starts in 0 CA-deg and finishes at 720 CA-deg. 
 

use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Debug)]
/// Struct used to read the engine data from .json file.
pub struct JsonEngine {
    /// [RPM]
    pub speed: f64,  
    /// [mm]      
    pub eccentricity: f64, 
    /// [mm]
    pub conrod: f64,     
    /// [cm³]  
    pub displacement: f64, 
    /// [mm]
    pub bore: f64,      
    /// i.e "1-3-2"   
    pub firing_order: String,
    pub combustion: Option<JsonCombustion>,
    pub injector: Option<JsonInjector>,
    pub cylinders: Vec<JsonCylinder>,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct JsonCylinder {
    pub name: String,
    pub compression_ratio: f64,
    /// [K]
    pub wall_temperature: f64,
    pub store_species: Option<bool>,
    pub intake_valves: Vec<JsonValve>,
    pub exhaust_valves: Vec<JsonValve>,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct JsonValve {
    pub name: String,
    /// Crank-angle degree [CA-deg]
    pub opening_angle: f64,
    /// Crank-angle degree [CA-deg]
    pub closing_angle: f64,
    /// [mm]
    pub diameter: f64,
    /// [mm]
    pub max_lift: f64,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct JsonCombustion {
    pub model: String,
    /// Crank-angle degree [CA-deg]
    pub comb_ini: f64,
    pub wiebe: JsonWiebe,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct JsonInjector {
    pub inj_type: String,
    pub air_fuel_ratio: f64,
    pub fuel: JsonFuel,
}
#[derive(Serialize, Deserialize, Debug)]
pub struct JsonFuel {
    pub name: String,
    pub state: String,
    pub lhv: Option<f64>,
    pub heat_vap: Option<f64>,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct JsonWiebe {
    pub a: f64,
    pub m: f64,
    pub comb_duration: f64,
}
