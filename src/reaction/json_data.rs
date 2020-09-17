// File to read and treat the data in the .json files

#![allow(unused_variables)]
#![allow(non_snake_case)]

use std::collections::HashMap;
use crate::reaction::thermo::ThermoInterp;
use ndarray::prelude::*;
use serde::{Deserialize, Serialize};

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
    atoms: String,
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
                    panic!(
                        "Specie {} must have only two temperature ranges",
                        specie_name
                    );
                }
                if t[0].len != t[1].len {
                    panic!(
                        "Specie {} must have coefficiets with the same size",
                        specie_name
                    );
                }
            }
            None => panic!("For {}, no thermo data detected in the file", specie_name),
        }
    }
    fn poly_interp_to_therm_interp(
        mut poly: Vec<PolynomInterp>,
        specie_name: &str,
    ) -> ThermoInterp {
        if poly[0].Tmax != poly[1].Tmin {
            panic!(
                "For specie {}, discontinuity temperature range in the polynomial",
                specie_name
            );
        } else {
            let Tmid = poly[0].Tmax;
            let coeffs_low = Array::from(poly[0].coeffs.take().unwrap());
            let coeffs_high = Array::from(poly[1].coeffs.take().unwrap());
            let thermo = ThermoInterp::new(specie_name.to_string(), Tmid, coeffs_low, coeffs_high);
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
    pub species_atoms: Vec<HashMap<String, f64>>,
    pub thermo_interp: Vec<ThermoInterp>,
}

pub fn read_and_treat_json(file_name: &str) -> OutputJson {
    // Reading .json file
    let json_file = std::fs::read_to_string(file_name).expect("Unable to read file");
    let gas: IdealGas = serde_json::from_str(&json_file).unwrap();

    let name = gas.phase.id.clone();
    let species = get_species(&gas);
    let species_atoms = get_species_atoms(&gas, &species);
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
        species_atoms,
        thermo_interp,
    }
}

fn get_species(gas: &IdealGas) -> Vec<String> {
    let name = gas.phase.id.clone();
    let species: Vec<String> = gas
        .phase
        .speciesArray
        .replace(&[',', '\"'][..], "")
        .split_whitespace()
        .map(|s| s.to_string())
        .collect();

    // Checking if there is `species_data` for all species
    for specie in species.iter() {
        if !gas.species_data.iter().any(|g| g.name == *specie) {
            println!("Specie `{}` was not found in `species_data` in file `{}`", specie, gas.phase.id);
            std::process::exit(1);
        }
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
    let strings: Vec<String> = strings
        .replace(&[',', '\"'][..], "")
        .split_whitespace()
        .map(|s| s.to_string())
        .collect();

    for word in strings.iter() {
        let specie: Vec<&str> = word.split(":").collect(); // `specie` should be a Vec like ["O2", "0.21"]
        if species.contains(&specie[0].to_string()) {
            let i = species.iter().position(|s| **s == *specie[0]).unwrap();
            mol_frac[i] = match specie[1].parse() {
                Ok(s) => s,
                Err(err) => {  
                    println!("{}", err);
                    std::process::exit(1);
                }
            }
        } else {
            println!("Error at `get_mol_frac`.{} not found in 'speciesArray'", specie[0]);
            std::process::exit(1);
        }
    }
    if mol_frac.sum() != 1.0 {
        println!("mol_fraction must sum 1.0: mol_frac = {}", mol_frac.sum());
        std::process::exit(1);
    }
    mol_frac
}

fn get_species_atoms(gas: &IdealGas, all_species: &Vec<String>) -> Vec<HashMap<String, f64>> {
    let mut species_atoms: Vec<HashMap<String, f64>> = vec![HashMap::new(); all_species.len()];

    for (i, specie_name) in all_species.iter().enumerate() {
        let specie_data = gas.species_data.iter().find(|s| &s.name == specie_name ).unwrap();
        let string: Vec<String> = specie_data.atoms
            .replace(&[',', '\"'][..], "")
            .split_whitespace()
            .map(|s| s.to_string())
            .collect();
        let mut map: HashMap<String, f64> = HashMap::new();
        for word in string {
            let s: Vec<&str> = word.split(":").collect(); // `s` should be a Vec like ["O", "2"]
            map.insert(s[0].to_string(), s[1].parse().unwrap());
        }
        species_atoms[i] = map;
    }
    species_atoms
}

fn get_thermo(mut gas: IdealGas, species: &Vec<String>) -> (Array1<f64>, Vec<ThermoInterp>) {
    let mut molecular_weight = Array::from_elem(species.len(), 0.);
    let mut thermo_interp_array: Vec<ThermoInterp> = Vec::new();
    for (index, specie_name) in species.iter().enumerate() {
        for data in gas.species_data.iter_mut() {
            if *data.name == *specie_name {
                molecular_weight[index] = data.molecular_weight;
                PolynomInterp::validade(&data.thermo, &species[index]);
                let thermo_interp = PolynomInterp::poly_interp_to_therm_interp(
                    data.thermo.take().unwrap(),
                    &species[index],
                );
                thermo_interp.validate();
                thermo_interp_array.push(thermo_interp);
            }
        }
    }
    (molecular_weight, thermo_interp_array)
}
