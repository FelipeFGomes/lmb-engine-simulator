#![allow(non_snake_case)]

use ndarray::prelude::*;
use crate::reaction::json_data;

#[derive(Debug)]
pub struct Gas {
    name: String,
    species: Vec<String>,
    species_molar_weight: Array1<f64>,
    mol_frac: Array1<f64>, 
    thermo_interp: Vec<Vec<json_data::PolynomInterp>>,
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
            mol_frac: json_output.mol_frac*1e-3, 
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
    pub fn TP(&mut self, temp: f64, press: f64) {
        self.thermo_prop.T = temp;
        self.thermo_prop.P = press;
        self.calc_prop();
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
    pub fn TPX(&mut self, temp: f64, press: f64, mol_frac: String) {
        self.thermo_prop.T = temp;
        self.thermo_prop.P = press;
        self.X(mol_frac);
    }

    /// Set mole fraction of species. Thermo properties are recalculated  
    /// # Examples
    /// ```
    /// let mol_frac = "O2:0.21, N2:0.79";
    /// gas.TP(temp, press);
    /// assert_eq!(0.21, gas.X("O2"));
    /// assert_eq!(0.79, gas.X("N2"));
    /// ```
    pub fn X(&mut self, mol_frac: String) {
        let strings: Vec<String> = mol_frac.replace(&[',', '\"'][..], "")
            .split_whitespace().map(|s| s.to_string()).collect();
        let mut X = Array::from_elem(self.species.len(), 0.);

        for word in strings.iter() {
            let specie: Vec<&str> = word.split(":").collect();  // specie should be like ["O2", "0.21"]
            if self.species.iter().all( |_| self.species.contains(&specie[0].to_string()) ) {
                let (i, _) = self.species.iter().enumerate().find(|(_, s)| **s == *specie[0]).unwrap();
                X[i] = specie[1].parse().unwrap();
            } else {
                panic!("{} not found in 'speciesArray'", specie[0]);
            }
        }
        if X.sum() != 1.0 {
            panic!("mol_fraction must sum 1.0: mol_frac = {}", X.sum());
        }
        self.mol_frac = X;
        self.calc_prop();
    }

    fn calc_prop(&mut self) {
        let R = 8.3143;
        let T = self.thermo_prop.T;
        let mut cp_species = Array::from_elem(self.species_molar_weight.len(), 0.);
        for (specie_index, thermo_interp) in self.thermo_interp.iter().enumerate() {
            let polyinterp = thermo_interp.iter().find(|polyinterp| -> bool {
                polyinterp.Tmin() <= T && T <= polyinterp.Tmax() 
            });
            if polyinterp.is_some() {
                let poly = polyinterp.unwrap();
                let T = T/100.0; 
                for (power, i) in (0..poly.len()).rev().zip(0..poly.len()) {
                    cp_species[specie_index] += poly.coeffs()[i]*T.powi(power as i32);
                }
            } else {
                panic!("Unable to estimate properties: temperature outside polynomial range");
            }
        }
        let cv_species = &cp_species - R;
        self.thermo_prop.M = self.mol_frac.dot( &self.species_molar_weight.t() );
        self.thermo_prop.cp =  self.mol_frac.dot( &(cp_species/self.thermo_prop.M) );
        self.thermo_prop.cv = self.mol_frac.dot( &(cv_species/self.thermo_prop.M) );
        self.thermo_prop.R = R/self.thermo_prop.M;
        self.thermo_prop.k = self.thermo_prop.cp/self.thermo_prop.cv;
        self.thermo_prop.rho = self.thermo_prop.P/(self.thermo_prop.R*self.thermo_prop.T);
        self.thermo_prop.e = self.thermo_prop.cv*self.thermo_prop.T;
        self.thermo_prop.h = self.thermo_prop.cp*self.thermo_prop.T;
        self.thermo_prop.s = self.thermo_prop.cp*(self.thermo_prop.T/298.15).ln() 
            - self.thermo_prop.R*(self.thermo_prop.P/101325.0).ln();
        self.thermo_prop.a = (self.thermo_prop.k*self.thermo_prop.R*self.thermo_prop.T).sqrt();
        self.thermo_prop.mu = (1.458e-6)*(self.thermo_prop.T*self.thermo_prop.T*self.thermo_prop.T/(self.thermo_prop.T+110.4)).sqrt();
    }

    pub fn name(&self) -> String {
        self.name.clone()
    }

    pub fn species(&self) -> &Vec<String> {
        &self.species
    }

    pub fn mol_frac(&self) -> &Array1<f64> {
        &self.mol_frac
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



#[derive(Debug)]
struct ThermoProp {
    pub P: f64,   // pressure
    pub T: f64,   // temperature
    pub rho: f64, // density
    pub cp: f64,  // specific heat capacity - cp
    pub cv: f64,  // specific heat capacity - cv
    pub R: f64,   // ideal gas constant [J/(kg.K)]
    pub k: f64,   // cp/cv
    pub M: f64,   // molecular weight
    pub e: f64,   // internal energy
    pub h: f64,   // enthalpy
    pub s: f64,   // entropy
    pub a: f64,   // sound speed
    pub mu: f64,  // viscosity
}

impl ThermoProp {
    fn new() -> ThermoProp {
        ThermoProp {
            P: 0.0,   // pressure
            T: 0.0,   // temperature
            rho: 0.0, // density
            cp: 0.0,  // specific heat capacity - cp
            cv: 0.0,  // specific heat capacity - cv
            R: 0.0,   // ideal gas constant [J/(kg.K)]
            k: 0.0,   // cp/cv
            M: 0.0,   // molecular weight
            e: 0.0,   // internal energy
            h: 0.0,   // enthalpy
            s: 0.0,   // entropy
            a: 0.0,   // sound speed
            mu: 0.0,  // viscosity
        }
    }
}
