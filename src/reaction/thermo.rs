#![allow(non_snake_case)]

use ndarray::prelude::*;

#[derive(Debug, Clone)]
pub struct ThermoInterp {
    specie: String,
    Tmid: f64,
    coeffs_low: Array1<f64>,
    coeffs_high: Array1<f64>,
}


#[derive(Debug, Clone)]
pub struct ThermoProp {
    pub P: f64,   // pressure [Pa]
    pub T: f64,   // temperature [K]
    pub rho: f64, // density [kg/m^3]
    pub cp: f64,  // specific heat capacity - cp [J/(kg.K)]
    pub cv: f64,  // specific heat capacity - cv [J/(kg.K)]
    pub R: f64,   // ideal gas constant [J/(kg.K)]
    pub k: f64,   // cp/cv
    pub M: f64,   // molecular weight 
    pub e: f64,   // internal energy
    pub h: f64,   // enthalpy
    pub s: f64,   // entropy
    pub a: f64,   // sound speed
    pub mu: f64,  // viscosity
}

impl ThermoInterp {
    pub fn new(specie: String, Tmid: f64, coeffs_low: Array1<f64>, coeffs_high: Array1<f64>) -> ThermoInterp {
        ThermoInterp {
            specie,
            Tmid,
            coeffs_low,
            coeffs_high,
        }
    }
    
    pub fn validate(&self) {
        let (cp_low, h_low, s_low) = ThermoInterp::calc_thermo_properties(&self.coeffs_low, self.Tmid);
        let (cp_high, h_high, s_high) = ThermoInterp::calc_thermo_properties(&self.coeffs_high, self.Tmid);
        
        //cp
        let delta = cp_low - cp_high;
        if (delta/(cp_low.abs()+1.0E-4)).abs() > 0.01 {
            panic!("ThermoInterp.validate(),
                \nFor species {}, discontinuity in cp/R detected at Tmid = {}
                Value computed using low-temperature polynomial:  {}
                Value computed using high-temperature polynomial: {}",
                self.specie, self.Tmid, cp_low, cp_high);
        }
        //enthalpy
        let delta = h_low - h_high;
        if delta.abs()/cp_low.abs() > 0.001 {
            panic!("ThermoInterp.validate(),
                \nFor species {}, discontinuity in h/RT detected at Tmid = {}
                Value computed using low-temperature polynomial:  {}
                Value computed using high-temperature polynomial: {}",
                self.specie, self.Tmid, cp_low, cp_high);
        }
        //entropy
        let delta = s_low - s_high;
        if (delta/(s_low.abs()+cp_low)).abs() > 0.001 {
            panic!("ThermoInterp.validate(),
                \nFor species {}, discontinuity in s/R detected at Tmid = {}
                Value computed using low-temperature polynomial:  {}
                Value computed using high-temperature polynomial: {}",
                self.specie, self.Tmid, cp_low, cp_high);
        }
    }

    /// Calculate non-dimensional cp, enthalpy and entropy for a given temperature using 4th order 
    /// NASA polinomial
    pub fn calc_thermo_properties(coeff: &Array1<f64>, temp: f64) -> (f64, f64, f64) {
        let cT0 = coeff[0];
        let cT1 = coeff[1]*temp;
        let cT2 = coeff[2]*temp.powi(2);
        let cT3 = coeff[3]*temp.powi(3);
        let cT4 = coeff[4]*temp.powi(4);
        let cT5 = coeff[5]/temp;
        let cT6 = coeff[0]*temp.ln();

        let cp_R = cT0 + cT1 + cT2 + cT3 + cT4;
        let h_RT = cT0 + 0.5*cT1 + 1.0/3.0*cT2 + 0.25*cT3 + 0.20*cT4 + cT5;
        let s_R = cT6 + cT1 + 0.5*cT2 + 1.0/3.0*cT3 + 0.25*cT4 + coeff[6];
        (cp_R, h_RT, s_R)
    }

    pub fn Tmid(&self) -> f64 {
        self.Tmid
    }

    pub fn _specie(&self) -> &str {
        &self.specie
    }

    pub fn coeffs_low(&self) -> &Array1<f64> {
        &self.coeffs_low
    }

    pub fn coeffs_high(&self) -> &Array1<f64> {
        &self.coeffs_high
    }
}

impl ThermoProp {
    pub fn new() -> ThermoProp {
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
