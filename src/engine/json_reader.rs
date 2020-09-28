use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Debug)]
pub struct JsonEngine {
    pub speed: f64,
    pub eccentricity: f64,
    pub conrod: f64,       // [mm]
    pub displacement: f64, // [cm³]
    pub bore: f64,         // [mm]
    pub firing_order: String,
    pub combustion: Option<JsonCombustion>,
    pub injector: Option<JsonInjector>,
    pub cylinders: Vec<JsonCylinder>,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct JsonCylinder {
    pub name: String,
    pub compression_ratio: f64, // [-]
    pub wall_temperature: f64,
    pub store_species: Option<bool>,
    pub intake_valves: Vec<JsonValve>,
    pub exhaust_valves: Vec<JsonValve>,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct JsonValve {
    pub name: String,
    pub opening_angle: f64,
    pub closing_angle: f64,
    pub diameter: f64, // [mm]
    pub max_lift: f64, // [mm]
}

#[derive(Serialize, Deserialize, Debug)]
pub struct JsonCombustion {
    pub model: String,
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
