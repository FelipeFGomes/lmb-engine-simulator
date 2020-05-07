#![allow(non_snake_case)]

use crate::base::constants::R;
use crate::reaction::{
    json_data,
    thermo::{ThermoInterp, ThermoProp},
};
use ndarray::prelude::*;

#[derive(Debug, Clone)]
pub struct Gas {
    name: String,
    species: Vec<String>,
    species_molar_weight: Array1<f64>,
    mol_frac: Array1<f64>,
    thermo_interp: Vec<ThermoInterp>,
    thermo_prop: ThermoProp,
}

impl Gas {
    /// Creates a `Gas` object from a file
    pub fn new(gas_file: &str) -> Gas {
        let json_output = json_data::read_and_treat_json(gas_file);

        let mut gas = Gas {
            name: json_output.name,
            species: json_output.species,
            species_molar_weight: json_output.species_molar_weight,
            mol_frac: json_output.mol_frac * 1e-3,
            thermo_interp: json_output.thermo_interp,
            thermo_prop: ThermoProp::new(),
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
    /// let mol_frac = "O2:0.21, N2:0.79";
    /// gas.TPX(temp, press, mol_frac);
    /// assert_eq!(350, gas.T());
    /// assert_eq!(2e5, gas.P());
    /// assert_eq!(0.21, gas.X("O2"));
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
    /// let mol_frac = "O2:0.21, N2:0.79";
    /// gas.TP(temp, press);
    /// assert_eq!(0.21, gas.X("O2"));
    /// assert_eq!(0.79, gas.X("N2"));
    /// ```
    pub fn X<'a>(&'a mut self, mol_frac: &str) -> &'a mut Self {
        let strings: Vec<String> = mol_frac
            .replace(&[',', '\"'][..], "")
            .split_whitespace()
            .map(|s| s.to_string())
            .collect();
        let mut X = Array::from_elem(self.species.len(), 0.);

        for word in strings.iter() {
            let specie: Vec<&str> = word.split(":").collect(); // specie should be like ["O2", "0.21"]
            if self
                .species
                .iter()
                .all(|_| self.species.contains(&specie[0].to_string()))
            {
                let (i, _) = self
                    .species
                    .iter()
                    .enumerate()
                    .find(|(_, s)| **s == *specie[0])
                    .unwrap();
                X[i] = specie[1].parse().unwrap();
            } else {
                panic!("{} not found in 'speciesArray'", specie[0]);
            }
        }
        if X.sum() != 1.0 {
            panic!("mol_fraction must sum 1.0: mol_frac = {}", X.sum());
        }
        self.mol_frac = X;
        self.update_prop();
        self
    }

    /// Set mole fraction of species from ndarray::Array1<f64>, must be the same size.
    /// Thermo properties are recalculated
    pub fn X_array<'a>(&'a mut self, mol_frac: Array1<f64>) -> &'a mut Self {
        self.mol_frac.assign(&mol_frac);
        self
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
        s_array = R * s_array; // [J/kmol/K]
        let cv_array = &cp_array - R;
        // All properties in mass basis
        self.thermo_prop.M = self.mol_frac.dot(&self.species_molar_weight);
        self.thermo_prop.cp = self.mol_frac.dot(&(cp_array / self.M()));
        self.thermo_prop.cv = self.mol_frac.dot(&(cv_array / self.M()));
        self.thermo_prop.h = self.mol_frac.dot(&(h_array / self.M()));
        self.thermo_prop.s = self.mol_frac.dot(&(s_array / self.M()));
        self.thermo_prop.R = R / self.M() / 1000.0;
        self.thermo_prop.k = self.cp() / self.cv();
        self.thermo_prop.rho = self.P() / (self.R() * self.T());
        self.thermo_prop.e = self.h() - self.P() / self.rho();
        self.thermo_prop.a = (self.k() * self.R() * self.T()).sqrt();
        self.thermo_prop.mu =
            (1.458e-6) * (self.T() * self.T() * self.T() / (self.T() + 110.4)).sqrt();
    }

    pub fn name(&self) -> String {
        self.name.clone()
    }

    pub fn species(&self) -> &Vec<String> {
        &self.species
    }

    pub fn mol_frac(&self) -> ArrayView1<f64> {
        self.mol_frac.view()
    }

    pub fn T(&self) -> f64 {
        self.thermo_prop.T
    }

    pub fn P(&self) -> f64 {
        self.thermo_prop.P
    }

    pub fn rho(&self) -> f64 {
        self.thermo_prop.rho
    }

    pub fn cp(&self) -> f64 {
        self.thermo_prop.cp
    }

    pub fn cv(&self) -> f64 {
        self.thermo_prop.cv
    }

    pub fn R(&self) -> f64 {
        self.thermo_prop.R
    }

    pub fn k(&self) -> f64 {
        self.thermo_prop.k
    }

    pub fn M(&self) -> f64 {
        self.thermo_prop.M
    }

    pub fn e(&self) -> f64 {
        self.thermo_prop.e
    }

    pub fn h(&self) -> f64 {
        self.thermo_prop.h
    }

    pub fn s(&self) -> f64 {
        self.thermo_prop.s
    }

    pub fn a(&self) -> f64 {
        self.thermo_prop.a
    }

    pub fn mu(&self) -> f64 {
        self.thermo_prop.mu
    }
}
