// File to read and treat the data in the .json files

#![allow(unused_variables)]
#![allow(non_snake_case)]

use serde::{Deserialize, Serialize};
use ndarray::prelude::*;
use crate::reaction::thermo::{ThermoInterp};

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
struct PolynomInterp {
    Tmin: f64,
    Tmax: f64,
    len: usize,
    coeffs: Option<Vec<f64>>,
}

impl PolynomInterp {
    fn validade(poly: &Option<Vec<PolynomInterp>>, specie_name: &str) {
        match poly {
            Some(t) => {
                if t.len() > 2 {
                    panic!("Specie {} must have only two temperature ranges", specie_name);
                }
                if t[0].len != t[1].len {
                    panic!("Specie {} must have coefficiets with the same size", specie_name);
                }
            },
            None => panic!("For {}, no thermo data detected in the file", specie_name)
        }
    }
    fn poly_interp_to_therm_interp(mut poly: Vec<PolynomInterp>, specie_name: &str) -> ThermoInterp {
        if poly[0].Tmax != poly[1].Tmin {
            panic!("For specie {}, discontinuity temperature range in the polynomial", specie_name);
        } else {
            let Tmid = poly[0].Tmax;
            let coeffs_low = Array::from( poly[0].coeffs.take().unwrap() );
            let coeffs_high = Array::from( poly[1].coeffs.take().unwrap() );
            let thermo = ThermoInterp::new(specie_name.to_string(),
                                            Tmid,
                                            coeffs_low,
                                            coeffs_high,);
            thermo.validate();
            thermo
        }
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
    pub thermo_interp: Vec<ThermoInterp>,
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

fn get_thermo(mut gas: IdealGas, species: &Vec<String>) -> (Array1<f64>, Vec<ThermoInterp>) {
    let mut molecular_weight = Array::from_elem(species.len(), 0.);
    let mut thermo_interp_array: Vec<ThermoInterp> = Vec::new();
    for (index, specie_name) in species.iter().enumerate() {
        for data in gas.species_data.iter_mut() {
            if *data.name == *specie_name {
                molecular_weight[index] = data.molecular_weight;
                PolynomInterp::validade(&data.thermo, &species[index]);
                let thermo_interp = 
                    PolynomInterp::poly_interp_to_therm_interp(data.thermo.take().unwrap(), &species[index]);
                thermo_interp.validate();
                thermo_interp_array.push(thermo_interp);
            }
        }
    }
    (molecular_weight, thermo_interp_array)
}

