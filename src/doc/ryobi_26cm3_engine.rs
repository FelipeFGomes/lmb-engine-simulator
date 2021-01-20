//! # Example of a 26cm³ single-cylinder engine
//!
//! ## The engine.json file
//!
//! ```
//! {
//!     "speed": 7500.0,
//!     "eccentricity": 0.0,
//!     "conrod": 58.2,
//!     "displacement": 26.2,
//!     "bore": 32.0,
//!     "firing_order": "1",
//!     "combustion": {
//!         "model": "Two-zone model",        
//!         "comb_ini": 710.0,
//!         "wiebe": {
//!             "m": 1.535,
//!             "a": 6.60,
//!             "comb_duration": 70.0
//!         }
//!     },
//!     "injector": {
//!         "inj_type": "port",
//!         "air_fuel_ratio": 0.896,
//!         "fuel": {
//!             "name": "C8H18",
//!             "state": "liquid",
//!             "lhv": 42.7e06,
//!             "heat_vap": 420.0e3
//!         }
//!     },
//!     "cylinders": [
//!         {
//!             "name": "cyl_1",
//!             "compression_ratio": 8.5,
//!             "wall_temperature": 470.0,
//!             "intake_valves": [
//!                 {
//!                     "name": "valve_int",
//!                     "opening_angle": 285.0,
//!                     "closing_angle": 630.0,
//!                     "diameter": 9.3,
//!                     "max_lift": 3.4
//!                 }
//!             ],
//!             "exhaust_valves": [
//!                 {   
//!                     "name": "valve_exh",
//!                     "opening_angle": 90.0,
//!                     "closing_angle": 435.0,
//!                     "diameter": 9.3,
//!                     "max_lift": 3.4
//!                 }
//!             ]
//!         }
//!     ]
//! }
//! ```
//!
//! ```
//! use lmb::Gas;
//! use lmb_engine_simulator as lmb;
//!
//! fn main() {
//!
//!     let gas = Gas::new("air.json");
//!     let gas_intake = Gas::new("air.json");
//!     let mut gas_exhaust = Gas::new("air.json");
//!     gas_exhaust.TPX(500.0, 101325.0, "N2:0.662586, H2O:0.202449, CO2:0.134965");
//!     
//!     let mut builder = lmb::SystemBuilder::new();
//!     builder.add_engine("engine_Ryobi.json", &gas)
//!         .add_environment("ambient", &gas)
//!         .add_reservoir("int_plenum", 250.0, &gas_intake)
//!         .add_reservoir("int_port", 2.04, &gas_intake)
//!         .add_reservoir("exh_port", 2.04, &gas_exhaust)
//!         .add_reservoir("exh_plenum_1", 95.0, &gas_exhaust)
//!         .add_reservoir("exh_plenum_2", 158.0, &gas_exhaust)
//!         .add_orifice("int_plenum -> amb 1", 9.5, 0.92, vec!["ambient", "int_plenum"])
//!         .add_orifice("int_plenum -> amb 2", 9.5, 0.92, vec!["ambient", "int_plenum"])
//!         .add_orifice("int_plenum -> int_port", 6.5, 0.78, vec!["int_plenum", "int_port"])
//!         .add_orifice("exh_plenum_2 -> amb", 8.9, 0.78, vec!["exh_plenum_2", "ambient"])
//!         .add_orifice("exh_plenum_1 -> exh_plenum_2", 12.7, 0.78, vec!["exh_plenum_1", "exh_plenum_2"])
//!         .add_orifice("exh_port -> exh_plenum_1", 10.8, 0.78, vec!["exh_port", "exh_plenum_1"])
//!         .connect_from_to("valve_int", "int_port")
//!         .connect_from_to("valve_exh", "exh_port");
//!     
//!     let mut system = builder.build_system();
//!
//!     // Calculating
//!     for speed in vec![5000.0, 5500.0, 6000.0, 6500.0, 7000.0, 7500.0, 8000.0, 8500.0, 9000.0] {
//!         system.engine_mut().unwrap().set_speed(speed);
//!         system.advance_to_steady_state();
//!
//!         // Writting data
//!         let folder_name = format!("./Ryobi_26_results/{:.0}_", speed);
//!         system.write_to_file( &(folder_name.clone() + "cylinder.txt"), "cyl_1", None);
//!         system.write_to_file(&(folder_name.clone() + "int_plenum.txt"), "int_plenum", None);
//!         system.write_to_file(&(folder_name.clone() + "int_port.txt"), "int_port", None);
//!         system.write_to_file(&(folder_name.clone() + "exh_port.txt"), "exh_port", None);
//!         system.write_to_file(&(folder_name.clone() + "exh_plenum_1.txt"), "exh_plenum_1", None);
//!         system.write_to_file(&(folder_name.clone() + "exh_plenum_2.txt"), "exh_plenum_2", None);
//!         system.write_to_file(&(folder_name.clone() + "int_valve.txt"), "valve_int", None);
//!         system.write_to_file(&(folder_name.clone() + "exh_valve.txt"), "valve_exh", None);
//!     }
//!
//!     system.engine().unwrap().write_performance_to("./Ryobi_26_results/engine_performance.txt");
//!
//! }
//! ```
//!
