use lmb_engine_simulator as lmb;
use lmb::Gas;
use std::time::{Instant};

fn main() {
    
    let mut gas1 = Gas::new("air.json");

    gas1.TPX(450.0, 101325.0, "O2:0.19, N2:0.77, AR:0.04");
    println!("Gas1:");
    println!("enthalpy = {}", gas1.h());
    println!("internal energy = {}", gas1.e());
    println!("entropy = {}", gas1.s());
    println!("cp = {}", gas1.cp());
    println!("mean mol. weight = {}", gas1.M());

    // let mut gas2 = Gas::new("air_lmb.json"); 
    // gas2.TP(1000.0, 100000.0);
    // println!("Gas2:");
    // println!("enthalpy = {}", gas2.h());
    // println!("entropy = {}", gas2.s());
    // println!("cp = {}", gas2.cp());
    // println!("k = {}", gas2.cp()/gas2.cv());

    // println!("Difference:");
    // println!("enthalpy = {}", gas1.h() - gas2.h());
    // println!("entropy = {}", gas1.s() - gas2.s());

}