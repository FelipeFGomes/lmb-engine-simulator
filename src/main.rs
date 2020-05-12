#![allow(unused_variables, unused_imports)]

use lmb::Gas;
use lmb_engine_simulator as lmb;
use std::f64::consts::PI;
use std::time::Instant;

fn main() {
    let gas = Gas::new("air.json");
    let mut builder = lmb::SystemBuilder::new();
    builder
        .add_environment("intake_env", &gas)
        .add_environment("exhaust_env", &gas)
        .add_engine("engine.json", &gas)
        .connect_from_to("valve_intake_1", "intake_env")
        .connect_from_to("valve_exhaust_1", "exhaust_env")
        .connect_from_to("valve_intake_2", "intake_env")
        .connect_from_to("valve_exhaust_2", "exhaust_env")
        .connect_from_to("valve_intake_3", "intake_env")
        .connect_from_to("valve_exhaust_3", "exhaust_env");
    let mut system = builder.build_system();

    // Calculating
    let now = Instant::now(); //measuring time
    system.advance_to_steady_state();
    println!("time taken: {:?}", now.elapsed());

    // writing to file
    system.write_to_file("cyl_1.txt", "cyl_1");
    system.write_to_file("cyl_2.txt", "cyl_2");
    system.write_to_file("cyl_3.txt", "cyl_3");
    system.write_to_file("valve_intake.txt", "valve_intake_1");
    system.write_to_file("valve_exhaust.txt", "valve_exhaust_1");

}
