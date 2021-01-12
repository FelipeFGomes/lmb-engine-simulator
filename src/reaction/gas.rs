//! # Gas Struct
//! The `Gas` struct is built from a .json file, use `air.json` as reference. 
//! In this file, there are two main structs: `"phase"` and `"species_data"`.
//! `"speciesArray"` determines which species can exist in the gas. 
//! `"state"` determines the initial state of the gas upon instantiation.
//! `"species_data"` is a vector containing information about all species in `"speciesArray"`.
//! Therefore, every specie declared in `"speciesArray"` must have your basic information declared in
//! `"species_data"`. The order in which the data is added is not importante. There can be more
//! `"specie_data"` than `"speciesArray"` but never the opposite.

#![allow(non_snake_case)]
use crate::base::constants::{R, _P_REF};
use crate::reaction::{
    json_data,
    thermo::{ThermoInterp, ThermoProp},
};
use ndarray::prelude::*;
use std::collections::HashMap;

/// Struct of an ideal gas
#[derive(Debug, Clone)]
pub struct Gas {
    name: String,
    species: Vec<String>,
    species_molar_weight: Array1<f64>,
    mole_frac: Array1<f64>,
    species_atoms: Vec<HashMap<String, f64>>,
    thermo_interp: Vec<ThermoInterp>,
    thermo_prop: ThermoProp,
    num_species: usize,
}

impl Gas {
    /// Creates a `Gas` object from a file .json file
    pub fn new(gas_file: &str) -> Gas {
        let json_output = json_data::read_and_treat_json(gas_file);
        let num_species = json_output.species.len();

        let mut gas = Gas {
            name: json_output.name,
            species: json_output.species,
            species_molar_weight: json_output.species_molar_weight,
            mole_frac: json_output.mol_frac,
            species_atoms: json_output.species_atoms,
            thermo_interp: json_output.thermo_interp,
            thermo_prop: ThermoProp::new(),
            num_species,
        };
        gas.TP(json_output.ini_temp, json_output.ini_press);
        gas
    }

    /// Set temperature and pressure. Thermo properties are recalculated  
    /// # Examples
    /// ```
    /// let temp = 350;
    /// let press = 2e5;
    /// gas.TP(temp, press);
    /// assert_eq!(350, gas.T());
    /// assert_eq!(2e5, gas.P());
    /// ```
    pub fn TP<'a>(&'a mut self, temp: f64, press: f64) -> &'a mut Self {
        self.thermo_prop.T = temp;
        self.thermo_prop.P = press;
        self.update_prop();
        self
    }

    /// Set temperature, pressure, and mole fraction of species. Thermo properties are recalculated  
    /// # Examples
    /// ```
    /// let temp = 350;
    /// let press = 2e5;
    /// let mole_frac = "O2:0.21, N2:0.79";
    /// gas.TPX(temp, press, mol_frac);
    /// assert_eq!(350, gas.T());
    /// assert_eq!(2e5, gas.P());
    /// assert_eq!(0.21, gas.mol_frac_of("O2"));
    /// ```
    pub fn TPX<'a>(&'a mut self, temp: f64, press: f64, mol_frac: &str) -> &'a mut Self {
        self.thermo_prop.T = temp;
        self.thermo_prop.P = press;
        self.X(mol_frac);
        self
    }

    /// Set mole fraction of species from `&str`. Thermo properties are recalculated  
    /// # Examples
    /// ```
    /// let mole_frac = "O2:0.21, N2:0.79";
    /// gas.X(mol_frac);
    /// assert_eq!(0.21, gas.mol_frac_of("O2"));
    /// assert_eq!(0.79, gas.mol_frac_of("N2"));
    /// ```
    pub fn X<'a>(&'a mut self, mole_frac: &str) -> &'a mut Self {
        let X = self.break_str_into_X_array(mole_frac);
        self.mole_frac = X;
        self.update_prop();
        self
    }

    /// Set mole fraction of species from ndarray::Array1<f64>, must be the same size.
    /// Thermo properties are recalculated
    pub fn X_array<'a>(&'a mut self, mol_frac: &Array1<f64>) -> &'a mut Self {
        if (mol_frac.sum() - 1.0).abs() > 1e-8 {
            println!(
                "Error!\n mol_fraction must sum 1.0: mol_frac = {}",
                mol_frac.sum()
            );
            std::process::exit(1);
        }
        self.mole_frac.assign(&mol_frac);
        self.update_prop();
        self
    }

    pub fn TPX_array<'a>(
        &'a mut self,
        temp: f64,
        press: f64,
        mol_frac: &Array1<f64>,
    ) -> &'a mut Self {
        self.thermo_prop.T = temp;
        self.thermo_prop.P = press;
        self.X_array(mol_frac);
        self
    }

    /// Returns the equivalent mole fraction, if a `mass`, in kg, of `self`
    /// were mixed with other gases in `add_gas<(mass, composition)>`
    /// # Examples
    /// ```
    /// gas.X("O2:1.0");
    /// let additional_gas = vec![(0.329337487, "N2:1.0")];
    /// let mix_mole_frac: Array1<f64> = gas.if_mixed_with(0.100, additional_gas);
    /// //mix_mole_frac has: O2:0.21, N2:0.79
    /// ```
    pub fn if_mixed_with(&self, mass: f64, add_gas: Vec<(f64, &str)>) -> Array1<f64> {
        if add_gas.iter().all(|(m, _)| *m == 0.0 as f64) {
            return self.mole_frac().clone();
        }
        let mut added_moles = Array::from_elem(self.num_species, 0.0);
        for (m, composition) in add_gas.iter() {
            let mole_frac = self.break_str_into_X_array(composition);
            let molar_weight = mole_frac.dot(&self.species_molar_weight);
            added_moles = added_moles + mole_frac * (*m) / molar_weight;
        }

        let new_mole_frac: Array1<f64>;
        let current_moles = self.mole_frac() * mass / self.M();
        let total_moles = current_moles.sum() + added_moles.sum();
        new_mole_frac = &(&current_moles + &added_moles) / total_moles;
        new_mole_frac
    }

    fn update_prop(&mut self) {
        let mut cp_array = Array::from_elem(self.species_molar_weight.len(), 0.);
        let mut h_array = Array::from_elem(self.species_molar_weight.len(), 0.);
        let mut s_array = Array::from_elem(self.species_molar_weight.len(), 0.);
        let mut cp_h_s: (f64, f64, f64);

        for (i, thermo_interp) in self.thermo_interp.iter().enumerate() {
            if self.T() < thermo_interp.Tmid() {
                cp_h_s = ThermoInterp::calc_thermo_properties(thermo_interp.coeffs_low(), self.T());
            } else {
                cp_h_s =
                    ThermoInterp::calc_thermo_properties(thermo_interp.coeffs_high(), self.T());
            }
            cp_array[i] = cp_h_s.0;
            h_array[i] = cp_h_s.1;
            s_array[i] = cp_h_s.2;
        }

        cp_array = R * cp_array; // [J/kmol/K]
        h_array = R * self.T() * h_array; // [J/kmol]
                                          // tmp: Array1<f64> = ln(mole_frac*P/P_ref)
        let tmp = self.mole_frac().mapv(|x| -> f64 {
            if x == 0.0 {
                0.0
            } else {
                (x * self.P() / _P_REF).ln()
            }
        });
        s_array = R * (s_array - tmp); // [J/kmol/K]
        let cv_array = &cp_array - R;
        // All properties in mass basis
        self.thermo_prop.M = self.mole_frac.dot(&self.species_molar_weight);
        self.thermo_prop.cp = self.mole_frac.dot(&(cp_array / self.M()));
        self.thermo_prop.cv = self.mole_frac.dot(&(cv_array / self.M()));
        self.thermo_prop.h = self.mole_frac.dot(&(h_array / self.M()));
        self.thermo_prop.s = self.mole_frac.dot(&(s_array / self.M()));
        self.thermo_prop.R = R / self.M();
        self.thermo_prop.k = self.cp() / self.cv();
        self.thermo_prop.rho = self.P() / (self.R() * self.T());
        self.thermo_prop.e = self.h() - self.P() / self.rho();
        self.thermo_prop.a = (self.k() * self.R() * self.T()).sqrt();
        self.thermo_prop.mu =
            (1.458e-6) * (self.T() * self.T() * self.T() / (self.T() + 110.4)).sqrt();
    }
    /// gas name
    pub fn name(&self) -> String {
        self.name.clone()
    }
    /// Vec with all species
    pub fn species(&self) -> &Vec<String> {
        &self.species
    }
    /// array of species mole fraction
    pub fn mole_frac(&self) -> &Array1<f64> {
        &self.mole_frac
    }
    /// array of all species molar weight [kg/kmol]
    pub fn mole_weight(&self) -> &Array1<f64> {
        &self.species_molar_weight
    }
    /// temperature [K]
    pub fn T(&self) -> f64 {
        self.thermo_prop.T
    }
    /// pressure [Pa]
    pub fn P(&self) -> f64 {
        self.thermo_prop.P
    }
    /// density [kg/m³]
    pub fn rho(&self) -> f64 {
        self.thermo_prop.rho
    }
    /// specific heat at constant pressure [J/kg.K]
    pub fn cp(&self) -> f64 {
        self.thermo_prop.cp
    }
    /// specific heat at constant volume [J/kg.K]
    pub fn cv(&self) -> f64 {
        self.thermo_prop.cv
    }
    /// ideal gas constant [J/kg.K]
    pub fn R(&self) -> f64 {
        self.thermo_prop.R
    }
    /// specific heat ratio (cp/cv)
    pub fn k(&self) -> f64 {
        self.thermo_prop.k
    }
    /// molar mass [kg/kmol]
    pub fn M(&self) -> f64 {
        self.thermo_prop.M
    }
    /// internal energy [J]
    pub fn e(&self) -> f64 {
        self.thermo_prop.e
    }
    /// enthalpy [J/kg]
    pub fn h(&self) -> f64 {
        self.thermo_prop.h
    }
    /// entropy [J/kg.K]
    pub fn s(&self) -> f64 {
        self.thermo_prop.s
    }
    /// sound speed [m/s]
    pub fn a(&self) -> f64 {
        self.thermo_prop.a
    }
    /// viscosity [m²/s]
    pub fn mu(&self) -> f64 {
        self.thermo_prop.mu
    }
    /// return the number os species
    pub fn num_species(&self) -> usize {
        self.num_species
    }
    /// check if `specie` contains in `species` vector
    pub fn contains_specie(&self, specie: &str) -> bool {
        self.species.contains(&specie.to_string())
    }
    /// get the index of `specie`
    pub fn get_specie_index(&self, specie: &str) -> usize {
        match self.species.iter().position(|s| s == specie) {
            Some(i) => i,
            None => {
                println!(
                    "Error at `get_specie_index()`. Specie \"{}\" not found",
                    specie
                );
                std::process::exit(1);
            }
        }
    }
    /// return the mole fraction of `specie`
    pub fn mole_frac_of(&self, specie: &str) -> f64 {
        let i = match self.species.iter().position(|s| s == specie) {
            Some(i) => i,
            None => {
                println!("Error at `mol_frac_of`. Specie \"{}\" not found", specie);
                std::process::exit(1);
            }
        };
        self.mole_frac[i]
    }
    /// return the mole weight of `specie`
    pub fn mole_weight_of(&self, specie: &str) -> f64 {
        let i = match self.species.iter().position(|s| s == specie) {
            Some(i) => i,
            None => {
                println!("Error at `mol_frac_of`. Specie {} not found", specie);
                std::process::exit(1);
            }
        };
        self.species_molar_weight[i]
    }
    /// return a `&HashMap` containing the `atoms_name` and `num_of_atoms` of `specie`
    pub fn atoms_of(&self, specie: &str) -> &HashMap<String, f64> {
        let i = match self.species.iter().position(|s| s == specie) {
            Some(i) => i,
            None => {
                println!("Error at `mol_frac_of`. Specie {} not found", specie);
                std::process::exit(1);
            }
        };
        &self.species_atoms[i]
    }

    fn break_str_into_X_array(&self, mole_frac: &str) -> Array1<f64> {
        let strings: Vec<String> = mole_frac
            .replace(&[',', '\"'][..], "")
            .split_whitespace()
            .map(|s| s.to_string())
            .collect();
        let mut X = Array::from_elem(self.species.len(), 0.);

        for word in strings.iter() {
            let specie: Vec<&str> = word.split(":").collect(); // specie should be like ["O2", "0.21"]
            if specie.len() != 2 {
                println!("Error setting mole fraction! \"{}\" is not valid", word);
                println!("Valid example: \"O2:0.21, N2:0.79\"");
                std::process::exit(1);
            }
            if self.species.contains(&specie[0].to_string()) {
                let i = self.species.iter().position(|s| s == specie[0]).unwrap();
                X[i] = specie[1].parse().unwrap();
            } else {
                println!(
                    "Error!\n Specie `{}` was not found in `species_data` in file `{}`",
                    specie[0],
                    self.name()
                );
                std::process::exit(1);
            }
        }
        if ((X.sum() - 1.0) as f64).abs() > 1e-8 {
            println!("Error!\n mol_fraction must sum 1.0: mol_frac = {}", X.sum());
            std::process::exit(1);
        }
        X
    }
}
