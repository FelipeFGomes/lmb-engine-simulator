use super::gas::Gas;
use crate::engine::engine::{Fuel};
use ndarray::prelude::*;
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;
use dyn_clone::DynClone;

pub trait Combustion: DynClone {
    fn model_name<'a>(&'a self) -> &str;
    /// Returns the HRR in `[J/CA-radian]`
    fn get_heat_release_rate(&mut self, gas: &Gas, fuel_mass: f64, angle: f64) -> f64;
    /// Returns an `Array1<f64>` containing an updated mole fraction. The species index are the same as `gas`.
    fn update_composition(&mut self, gas: &Gas, mass: f64, time: f64, press: f64, vol: f64) -> Array1<f64>;
    /// Returns the initial of combustion phase in crank-angle radians
    fn ini_combustion(&self) -> f64;
}

dyn_clone::clone_trait_object!(Combustion);

#[derive(Debug, Clone)]
pub struct TwoZoneCombustion {
    model_name: String,
    unburned_zone: Gas,
    burned_zone: Gas,
    fuel: Fuel,
    ini_combustion: f64,
    end_combustion: f64,
    delta_angle: f64,
    comb_eficiency: f64,
    air_fuel_ratio: f64,
    is_comb_ready: bool,
    wiebe_function: WiebeFunction,
}
impl TwoZoneCombustion {
    pub fn new(ign_angle: f64, afr: f64, wiebe: &WiebeFunction, gas: &Gas, fuel: &Fuel, air_comp: &Gas) -> TwoZoneCombustion {
        let ign_angle = ign_angle.to_radians();
        let wiebe_function = WiebeFunction::new(wiebe.a, wiebe.m, wiebe.comb_duration.to_radians());
        let mut end_combustion = ign_angle + wiebe_function.comb_duration;
        if end_combustion > 4.0*PI {
            end_combustion -= 4.0*PI;
        }
        let fuel = fuel.clone();
        let air_o2_frac = air_comp.mole_frac_of("O2");
        
        // burned-zone mole fraction - Fixed composition: assuming complete combustion
        let mut new_mole_frac = Array::from_elem(gas.species().len(), 0.);
        if afr >= 1.0 {
            let mole_co2 = fuel.carbon() + air_comp.mole_frac_of("CO2")/air_o2_frac;
            let mole_h2o = 0.5*fuel.hydrogen() + air_comp.mole_frac_of("H2O")/air_o2_frac;
            let mole_n2 = 0.5*fuel.nitrogen() + afr*fuel.air_moles()*(air_comp.mole_frac_of("N2")/air_o2_frac);
            let mole_o2 = fuel.air_moles()*(afr - 1.0);
            let total_moles = mole_co2 + mole_h2o + mole_n2 + mole_o2;
            let i_co2 = gas.species().iter().position(|s| s == "CO2").unwrap();
            let i_h2o = gas.species().iter().position(|s| s == "H2O").unwrap();
            let i_n2 = gas.species().iter().position(|s| s == "N2").unwrap();
            let i_o2 = gas.species().iter().position(|s| s == "O2").unwrap();

            new_mole_frac[i_co2] = mole_co2/total_moles;
            new_mole_frac[i_h2o] = mole_h2o/total_moles;
            new_mole_frac[i_n2] = mole_n2/total_moles;
            new_mole_frac[i_o2] = mole_o2/total_moles;
        } else {
            let mole_co2 = fuel.carbon();
            let mole_h2o = 0.5*fuel.hydrogen();
            let mole_n2 = 0.5*fuel.nitrogen() + afr*fuel.air_moles()*79.0/21.0;
            let total_moles = mole_co2 + mole_h2o + mole_n2;
            let i_co2 = gas.species().iter().position(|s| s == "CO2").unwrap();
            let i_h2o = gas.species().iter().position(|s| s == "H2O").unwrap();
            let i_n2 = gas.species().iter().position(|s| s == "N2").unwrap();

            new_mole_frac[i_co2] = mole_co2/total_moles;
            new_mole_frac[i_h2o] = mole_h2o/total_moles;
            new_mole_frac[i_n2] = mole_n2/total_moles;
        }

        let mut burned_zone = gas.clone();
        burned_zone.X_array(&new_mole_frac);

        TwoZoneCombustion {
            model_name: "Two-zone model".to_string(),
            unburned_zone: gas.clone(),
            burned_zone,
            fuel,
            ini_combustion: ign_angle,
            end_combustion,
            delta_angle: end_combustion - ign_angle,
            comb_eficiency: 0.99*(-1.602+4.6509*afr-2.0746*(afr*afr)),
            air_fuel_ratio: afr,
            is_comb_ready: false,
            wiebe_function,
        }

    }

    fn has_started(&self, angle: f64) -> bool {
        // checking if combustion has started
        if self.delta_angle < 0.0 {
            if angle <= self.ini_combustion && angle >= self.end_combustion {
                false
            } else {
                true
            }
        } else {
            if angle <= self.ini_combustion || angle >= self.end_combustion {
                false
            }
            else {
                true
            }
        }
    }
}

impl Combustion for TwoZoneCombustion {
    fn model_name<'a>(&'a self) -> &str {&self.model_name}
    fn get_heat_release_rate(&mut self, gas: &Gas, fuel_mass: f64, angle: f64) -> f64 {
        if self.has_started(angle) {
            if !self.is_comb_ready {
                // setting new unburned gas
                self.unburned_zone = gas.clone();
                self.is_comb_ready = true;
                // println!("fuel mass: {:.5} [mg]\tmass: {:.5}[mg]", self.fuel.mass()*1e6, mass*1e6);
            }
            self.comb_eficiency*fuel_mass*self.fuel.lhv()*self.wiebe_function.derivative_burned_mass_frac(angle, self.ini_combustion)
        } else {
            self.is_comb_ready = false;
            0.0
        }
        
    }
    fn update_composition(&mut self, gas: &Gas, mass: f64, angle: f64, press: f64, vol: f64) -> Array1<f64>{
        // updating unburned-zone: Adiabatic model
        let rho_un = self.unburned_zone.rho()*(press/self.unburned_zone.P()).powf(1.0/self.unburned_zone.k());
        let temp_un = press/(rho_un*self.unburned_zone.R());
        
        // updating burned-zone     
        let burned_mass_frac:f64;
        if self.has_started(angle) && self.is_comb_ready {
            burned_mass_frac = self.wiebe_function.burned_mass_frac(angle, self.ini_combustion);
            if burned_mass_frac > 0.993 {
                return self.burned_zone.mole_frac().clone();
            }
        } else {
            return gas.mole_frac().clone();
        }

        // println!("angle: {:.2}\t Xb: {}", angle.to_degrees(), burned_mass_frac);
        let mass_bur = burned_mass_frac*mass;
        let mass_un = mass - mass_bur;
        let vol_un = mass_un/rho_un;
        let vol_bur = vol - vol_un;
        let rho_bur = mass_bur/vol_bur;
        let temp_bur = press/(rho_bur*self.burned_zone.R());

        // Updating compositions: Fixed composition
        self.unburned_zone.TP(temp_un, press);
        self.burned_zone.TP(temp_bur, press);

        // returning new gas mixture
        let unburned_moles = self.unburned_zone.mole_frac()/self.unburned_zone.M()*mass_un;
        let burned_moles = self.burned_zone.mole_frac()/self.burned_zone.M()*mass_bur;
        let total_moles = unburned_moles.sum() + burned_moles.sum();
        let gas_moles = unburned_moles + burned_moles;
        let new_mole_frac = &gas_moles/total_moles;
        new_mole_frac
    }
    fn ini_combustion(&self) -> f64 {self.ini_combustion}
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct WiebeFunction {
    m: f64,
    a: f64,
    comb_duration: f64,   // [crank angle degree]
}
impl WiebeFunction {
    fn new(a:f64, m:f64, comb_duration:f64) -> WiebeFunction {
        WiebeFunction {
            m,
            a,
            comb_duration: comb_duration,
        }
    }

    // angle and ini_combustion must be in crank-angle degree
    fn burned_mass_frac(&self, angle: f64, ini_combustion: f64) -> f64 {
        let a = self.a;
        let m = self.m;
        let comb_duration = self.comb_duration;
        let d_angle = if (angle - ini_combustion) >= 0.0 {
            angle - ini_combustion
        } else {
            angle - ini_combustion + 4.0*PI
        };
        1.0 - (-a*(d_angle/comb_duration).powf(m+1.0)).exp()
    }

    // angle and ini_combustion must be in crank-angle degree
    fn derivative_burned_mass_frac(&self, angle: f64, ini_combustion: f64) -> f64 {
        let a = self.a;
        let m = self.m;
        let comb_duration = self.comb_duration;
        let d_angle = if (angle - ini_combustion) >= 0.0 {
            angle - ini_combustion
        } else {
            angle - ini_combustion + 4.0*PI
        };
        let tmp = d_angle/comb_duration;
        a*(m+1.0)/comb_duration*tmp.powf(m)*(-a*tmp.powf(m+1.0)).exp()
    }

}

#[derive(Debug, Clone)]
pub struct NoCombustion {model_name: String}
impl NoCombustion {
    pub fn new() -> NoCombustion {NoCombustion{model_name: "no combustion model".to_string()}}
}
impl Combustion for NoCombustion {
    fn model_name<'a>(&'a self) -> &'a str {&self.model_name}
    fn get_heat_release_rate(&mut self, _: &Gas, _: f64, _: f64) -> f64 {0.0}
    fn update_composition(&mut self, gas: &Gas, _: f64, _: f64, _: f64, _: f64) -> Array1<f64> {gas.mole_frac().clone()}
    fn ini_combustion(&self) -> f64 {0.0}
}