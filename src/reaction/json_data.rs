// File to read and treat the data in the .json files

#![allow(unused_variables)]
#![allow(non_snake_case)]

use serde::{Deserialize, Serialize};
use ndarray::prelude::*;

#[derive(Serialize, Deserialize, Debug)]
struct IdealGas {
    phase: Phase,
    species_data: Vec<SpeciesData>,
}

#[derive(Serialize, Deserialize, Debug)]
struct Phase {
    id: String,
    speciesArray: String,
    state: State,
}

#[derive(Serialize, Deserialize, Debug)]
struct SpeciesData {
    name: String,
    molecular_weight: f64,
    thermo: Option<Vec<PolynomInterp>>,
}

#[derive(Serialize, Deserialize, Debug)]
struct State {
    temperature: f64,
    pressure: f64,
    moleFractions: String,
}

#[derive(Serialize, Deserialize, Debug)]
#[allow(non_snake_case)]
pub struct PolynomInterp {
    Tmin: f64,
    Tmax: f64,
    len: usize,
    coeffs: Vec<f64>,
}
impl PolynomInterp {
    pub fn Tmin(&self) -> f64 {
        self.Tmin
    }
    pub fn Tmax(&self) -> f64 {
        self.Tmax
    }
    pub fn len(&self) -> usize {
        self.len
    }
    pub fn coeffs(&self) -> &Vec<f64> {
        &self.coeffs
    }
}

#[derive(Debug)]
pub struct OutputJson {
    pub name: String,
    pub species: Vec<String>,
    pub ini_temp: f64,
    pub ini_press: f64,
    pub mol_frac: Array1<f64>,
    pub species_molar_weight: Array1<f64>,
    pub thermo_interp: Vec<Vec<PolynomInterp>>,
}


pub fn read_and_treat_json(file_name: &str) -> OutputJson {
    // Reading .json file
    let json_file = std::fs::read_to_string(file_name).expect("Unable to read file");
    let gas: IdealGas = serde_json::from_str(&json_file).unwrap();

    let name = gas.phase.id.clone();
    let species = get_species(&gas);
    let (ini_temp, ini_press) = get_ini_state(&gas);
    let mol_frac = get_mol_frac(&gas, &species);
    let (species_molar_weight, thermo_interp) = get_thermo(gas, &species);

    OutputJson {
        name,
        species,
        ini_temp,
        ini_press,
        mol_frac,
        species_molar_weight,
        thermo_interp,
    }
}

fn get_species(gas: &IdealGas) -> Vec<String> {
    let name = gas.phase.id.clone();
    let species: Vec<String> = gas.phase.speciesArray.split_whitespace().map(|s| s.to_string()).collect();
    
    // Checking if file is appropriate
    if species.len() > gas.species_data.len() {
        panic!("not enough data for the species in 'speciesArray'");
    }
    species
}

fn get_ini_state(gas: &IdealGas) -> (f64, f64) {
    let ini_temp = gas.phase.state.temperature;
    let ini_press = gas.phase.state.pressure;
    (ini_temp, ini_press)
}

fn get_mol_frac(gas: &IdealGas, species: &Vec<String>) -> Array1<f64> {
    let mut mol_frac = Array::from_elem(species.len(), 0.);
    let strings = gas.phase.state.moleFractions.clone();
    let strings: Vec<String> = strings.replace(&[',', '\"'][..], "").split_whitespace().map(|s| s.to_string()).collect();

    for word in strings.iter() {
        let specie: Vec<&str> = word.split(":").collect();  // specie should be like ["O2", "0.21"]
        if species.iter().all( |n| species.contains(&specie[0].to_string()) ) {
            let (i, _) = species.iter().enumerate().find(|(i, s)| **s == *specie[0]).unwrap();
            mol_frac[i] = specie[1].parse().unwrap();
        } else {
            panic!("{} not found in 'speciesArray'", specie[0]);
        }
    }
    if mol_frac.sum() != 1.0 {
        panic!("mol_fraction must sum 1.0: mol_frac = {}", mol_frac.sum());
    }
    mol_frac
}

fn get_thermo(mut gas: IdealGas, species: &Vec<String>) -> (Array1<f64>, Vec<Vec<PolynomInterp>>) {
    let mut molecular_weight = Array::from_elem(species.len(), 0.);
    let mut thermo_interp: Vec<Vec<PolynomInterp>> = Vec::new();
    for (i, specie) in species.iter().enumerate() {
        for data in gas.species_data.iter_mut() {
            if *data.name == *specie {
                molecular_weight[i] = data.molecular_weight;
                thermo_interp.push( data.thermo.take().unwrap() );
            }
        }
    }
    (molecular_weight, thermo_interp)
}



