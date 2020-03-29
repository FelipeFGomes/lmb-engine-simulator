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
    
    pub fn TP(&mut self, temp: f64, press: f64) {
        self.thermo_prop.T = temp;
        self.thermo_prop.P = press;
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
        println!("{:?}, {:?}", self.mol_frac.shape(), self.species_molar_weight.shape()); 
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
