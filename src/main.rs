#![allow(unused_variables)]

use lmb::Gas;
use lmb_engine_simulator as lmb;
use std::f64::consts::PI;
use std::io::Write;

// use std::time::Instant;

fn main() {
    let gas = Gas::new("air.json");
    // gas.TP(300.0, 100000.0);
    let mut builder = lmb::SystemBuilder::new();
    builder.add_environment("env_1", &gas);
    builder.add_environment("env_2", &gas);
    builder.add_engine("engine.json", &gas);
    builder.connect_from_to("valve_intake", "env_1");
    builder.connect_from_to("valve_exhaust", "env_2");
    
    let mut system = builder.build_system();

    let speed = 3000.0/60.0; // [RPS]
    let sec_to_rad = 2.0*PI*speed;
    let d_angle = 0.01; // [CA deg]
    let dt = (d_angle*PI/180.0)/sec_to_rad;
    let num_cycles = 10.0;
    let limit  = (num_cycles*(720.0/d_angle)).round();

    for _ in 0..limit as usize {
        system.advance(dt);
    }

    let first = system.stored_data.len() - (720.0/d_angle) as usize;
    let final_cycle = &system.stored_data[first..];
    let mut file = std::fs::File::create("result").expect("Error opening writing file");
    write!(file, "{}", final_cycle.join("")).expect("Unable to write data");

}
