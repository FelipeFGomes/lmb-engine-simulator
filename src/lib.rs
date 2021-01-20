//! # lmb_engine_simulator
//!
//! The `lmb_engine_simulator` crate provides an easy way to simulate internal combustion engines.
//!
//! This library employed the **Builder Pattern** so the user feels as she/he is acctually building
//! an engine system. To construct (build) the system, the struct [`SystemBuilder`](core/system_builder/struct.SystemBuilder.html)
//! is used to add the desired components. After finishing building, the method `build_system()` can
//! be used to return an object [`System`](core/system/struct.System.html) which is used to solve the components numerically.
//!
//! ## Current stage
//!
//! The library allows you to build [Gas](reaction/gas/index.html), [Zero Dimensional](zero_dim/index.html) and [Connector](connector/index.html) objects.
//! The basic logic of the program is to use the `SystemBuilder` to build all the objects and then
//! indicate how the objects are connected with each other. Basicly, all of the dimensional objects
//! require a `Gas` to be created and a `Gas` object is created from a .json file. Currently, the
//! only available one is `air.json`. However, they were made in a way to be easily contructed.
//! Once all dimensional objects and connections were added, the `SystemBuilder` object can use the
//! method `build_system()` to create an object of type `System`, which will carry all the
//! [Zero Dimensional](zero_dim/index.html) and [Connector](connector/index.html) objects. The `System`
//! object is responsible to manage the interaction between the objects, the simulation process, the storable variables.
//!
//! In order to simulate, two methods can be used: [`advance(dt)`](core/system/struct.System.html#method.advance) 
//! which advance the state of the objects by `dt` and [`advance_to_steady_state()`](core/system/struct.System.html#method.advance_to_steady_state)
//! which advances until the system reaches steady state. **OBS: So far, steady state conditions are not checked. 
//! The method runs long enough that almost every system will have reached steady state by then.**
//! For engine simulation, most commonly, it is used [`advance_to_steady_state()`](core/system/struct.System.html#method.advance_to_steady_state).
//! After the simulation is finished, all stored variables can only be accessed by writing them into a 
//! file via system method [`write_to_file()`](core/system/struct.System.html#method.write_to_file) 
//!
//! ### Example
//! A simple system with a [Reservoir](zero_dim/reservoir/struct.Reservoir.html) and [Environment](zero_dim/environment/struct.Environment.html)
//! connected by an [Orifice](connector/orifice/struct.Orifice.html) is created and simulated until steady state. After, the stored data is 
//! written into two files.
//! ```
//! use lmb::Gas;
//! use lmb_engine_simulator as lmb;
//! 
//! fn main() {
//!     let mut gas_ambient = Gas::new("air.json");
//!     gas_ambient.TPX(293.0, 2.0*101325.0, "O2:0.21, N2:0.79");
//!     let mut gas_chamber = Gas::new("air.json");
//!     gas_chamber.TPX(293.0, 101325.0, "O2:0.21, N2:0.79");
//!     let mut builder = lmb::SystemBuilder::new();
//!     builder
//!         .add_environment("ambient", &gas_ambient)
//!         .add_reservoir("chamber", 500.0, &gas_chamber)
//!         .add_orifice("orifice", 50.0, 0.9, vec!["ambient", "chamber"]);
//!
//!     let mut system = builder.build_system();
//!    
//!     // Calculating
//!     system.advance_to_steady_state();
//!    
//!     // Writting data
//!     system.write_to_file("chamber.txt", "chamber", None);
//!     system.write_to_file("orifice.txt", "orifice", None);
//! }

//! ```
//! 
//! ## Simulating engines
//! 
//! To add an engine to the system, the method [`add_engine("engine_file.json", "gas")`](core/system_builder/struct.SystemBuilder.html#method.add_engine)
//! of `SystemBuilder` must be used. The `engine_file.json` is read into the struct [`json_engine`](engine/json_reader/struct.JsonEngine.html). 
//! This file **must** have at least the following attributes:
//! * "speed" in RPM,
//! * "eccentricity" in mm,
//! * "conrod" in mm,
//! * "displacement" in cm³,
//! * "bore" in mm,
//! * "firing_order" as a string - i.e "1-3-2",
//! * "cylinders" as a vector of structs [`json_cylinder`](engine/json_reader/struct.JsonCylinder.html)
//! 
//! All possible attributes of the `engine_file.json` file can be found at the full documentation at [`Json Reader`](engine/json_reader/index.html).
//! **Attention when entering the variables in crank-angle degree!** The reference, where crank-angle is zero, 
//! is at top-dead-center (TDC) of compression phase and it only accepts positive numbers. 
//! Therefore, the full cycle starts in 0 CA-deg and finishes at 720 CA-deg. 
//! 
//! ### Example
//! Let's simulate a simple engine system with intake and exhaust manifolds as [environments](zero_dim/environment/struct.Environment.html)
//! connected to a single cylinder through only two valves (intake and exhaust). The `engine.json` file will be
//! ```
//! {
//!     "speed": 3000.0,
//!     "eccentricity": 0.0,
//!     "conrod": 145.6,
//!     "displacement": 400.0,
//!     "bore": 80.0,
//!     "firing_order": "1",
//!     "combustion": {
//!         "model": "Two-zone model",        
//!         "comb_ini": 690.0,
//!         "wiebe": {
//!             "m": 2.0,
//!             "a": 6.908,
//!             "comb_duration": 40.0
//!         }
//!     },
//!     "injector": {
//!         "inj_type": "port",
//!         "air_fuel_ratio": 1.0,
//!         "fuel": {
//!             "name": "C2H5OH",
//!             "state": "liquid",
//!             "lhv": 25.858e06,
//!             "heat_vap": 900.0e3
//!         }
//!     },
//!     "cylinders": [
//!         {
//!             "name": "cyl_1",
//!             "compression_ratio": 12.50,
//!             "wall_temperature": 520.0,
//!             "store_species": true,
//!             "intake_valves": [
//!                 {
//!                     "name": "valve_int",
//!                     "opening_angle": 340.0,
//!                     "closing_angle": 570.0,
//!                     "diameter": 30.93,
//!                     "max_lift": 9.30
//!                 }
//!             ],
//!             "exhaust_valves": [
//!                 {
//!                     "name": "valve_exh",
//!                     "opening_angle": 130.0,
//!                     "closing_angle": 375.0,
//!                     "diameter": 28.27,
//!                     "max_lift": 8.48
//!                 }
//!             ]
//!         }
//!     ]
//! }
//! ```
//! 
//! In this example, we added an [`injector`](engine/json_reader/struct.JsonInjector.html), with relative air fuel ratio equal 1.0 
//! and ethanol (C2H5OH) as fuel, and a [`combustion`](engine/json_reader/struct.JsonCombustion.html) model. 
//! If they are not added, the engine will run as a motoring. 
//! Right now, the **only combustion model implemented** is the "Two-zone model". Notice that the 
//! [`cylinder`](engine/json_reader/struct.JsonCylinder.html) requires both intake and exhaust [`valves`](engine/json_reader/struct.JsonValve.html) 
//! connected to it. In the `main.rs`, we will need to connect these valves to their ports with 
//! [`connect_from_to("valve_name", "object_name")`](core/system_builder/struct.SystemBuilder.html#method.connect_from_to) method.
//! 
//! The `main.rs` can be written as: 
//! 
//! ```
//! use lmb::Gas;
//! use lmb_engine_simulator as lmb;
//!
//! fn main() {
//!     let gas_intake = Gas::new("air.json");
//!     let mut gas_exhaust = Gas::new("air.json");
//!     gas_exhaust.TPX(500.0, 101325.0, "N2:0.662586, H2O:0.202449, CO2:0.134965");
//!     let mut builder = lmb::SystemBuilder::new();
//!     builder
//!         .add_engine("engine.json", &gas_intake)
//!         .add_environment("intake_port", &gas_intake)
//!         .add_environment("exhaust_port", &gas_exhaust)
//!         .connect_from_to("valve_int", "intake_port")
//!         .connect_from_to("valve_exh", "exhaust_port");
//!
//!     let mut system = builder.build_system();
//!    
//!     // Calculating
//!     system.advance_to_steady_state();
//!    
//!     // Writting data
//!     system.write_to_file("cylinder.txt", "cyl_1", None);
//!     system.write_to_file("intake_valve.txt", "valve_int", None);
//!     system.write_to_file("exhaust_valve.txt", "valve_exh", None);
//! }
//! ```
//! 
//! For a real-life engine simulation, see [Engine Examples](doc/Ryobi_26cm3_engine/index.html)


use crate::base::constants::MAX_ARRAY_LEN;
use ndarray::*;
use std::io::Write;
use std::ops::Add;

pub mod base;
pub mod connector;
pub mod core;
pub mod engine;
pub mod numerics;
pub mod one_dim;
pub mod reaction;
pub mod zero_dim;
mod doc;

// Re-exporting
pub use crate::core::system_builder::SystemBuilder;
pub use crate::engine::engine::Engine;
pub use crate::reaction::combustion;
pub use crate::reaction::gas::Gas;

// Object's type
#[derive(Debug, Clone)]
pub enum ObjectType {
    ZeroDim,
    OneDim,
    Connector,
    Cylinder,
}

#[derive(Debug, Clone)]
pub struct BasicProperties<'a> {
    pub name: &'a str,
    pub pressure: f64,    // Pa
    pub temperature: f64, // K
    pub cp: f64,          // J/(kg.K)
    pub cv: f64,          // J/(kg.K)
    pub cp_cv: f64,
    pub gas_const: f64,           // J/(kg.K)
    pub crank_angle: Option<f64>, // CA rad
}

impl<'a> std::fmt::Display for BasicProperties<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}: 
        pressure: {} [Pa] 
        temperature: {} [K]
        cp: {} [J/(kg.K)]
        cv: {} [J/(kg.K)]
        cp/cv: {}
        R: {} [J/(kg.K)]",
            self.name,
            self.pressure,
            self.temperature,
            self.cp,
            self.cv,
            self.cp_cv,
            self.gas_const
        )
    }
}

#[derive(Debug, Clone)]
pub struct FlowRatio {
    pub mass_flow: f64,     // kg/s
    pub enthalpy_flow: f64, // J/s
}

impl FlowRatio {
    pub fn new() -> FlowRatio {
        FlowRatio {
            mass_flow: 0.0,
            enthalpy_flow: 0.0,
        }
    }
}

impl Add for FlowRatio {
    type Output = FlowRatio;
    fn add(self, other: FlowRatio) -> FlowRatio {
        FlowRatio {
            mass_flow: self.mass_flow + other.mass_flow,
            enthalpy_flow: self.enthalpy_flow + other.enthalpy_flow,
        }
    }
}

impl<'a, 'b> Add<&'b FlowRatio> for &'a FlowRatio {
    type Output = FlowRatio;
    fn add(self, other: &'b FlowRatio) -> FlowRatio {
        FlowRatio {
            mass_flow: self.mass_flow + other.mass_flow,
            enthalpy_flow: self.enthalpy_flow + other.enthalpy_flow,
        }
    }
}

#[derive(Debug, Clone)]
pub struct StoreData {
    header: String,
    data: Array2<f64>,
    last_index: usize,
}

impl StoreData {
    fn new(header: &str, num_variables: usize) -> StoreData {
        StoreData {
            header: header.to_string(),
            data: Array::from_elem((MAX_ARRAY_LEN, num_variables), 0.),
            last_index: 0,
        }
    }

    fn add_data(&mut self, data: Array1<f64>) {
        if self.last_index == MAX_ARRAY_LEN - 1 {
            println!("Error! Maximum allow array length exceeded!");
            std::process::exit(1);
        }
        self.data.row_mut(self.last_index).assign(&data);
        self.last_index += 1;
    }

    fn get_data(&self, rows_range: (usize, usize), columns: Vec<usize>) -> Array2<f64> {
        let mut data: Array2<f64> = Array2::zeros((rows_range.1 - rows_range.0, columns.len()));
        for (i, col) in columns.iter().enumerate() {
            let filtered_data = self.data.slice(s![rows_range.0..rows_range.1, *col]);
            data.slice_mut(s![.., i]).assign(&filtered_data);
        }
        data
    }

    fn reset_data(&mut self) {
        let num_variables = self.data.ncols();
        self.data = Array::from_elem((MAX_ARRAY_LEN, num_variables), 0.);
        self.last_index = 0;
    }

    /// Write the stored data in `data` limited by the index `range`to a file,
    ///  the first line is the content in `header`.
    fn write_to_file(
        &self,
        file_name: &str,
        range: (usize, usize),
        additional_data: Option<(String, ArrayView2<f64>)>,
    ) {
        let data: ArrayView2<f64>;
        let tmp: Array2<f64>;
        let mut additional_header = String::from("");
        let filtered_data = self.data.slice(s![range.0..range.1, ..]);
        if let Some((header, add)) = additional_data {
            if filtered_data.nrows() != add.nrows() {
                println!(
                    "`additional_data` must have the same number of rows as the writable data"
                );
                println!(
                    " `additional_data`: {}, writable data: {}",
                    add.len(),
                    filtered_data.len()
                );
                std::process::exit(1);
            }
            tmp = stack![Axis(1), add, filtered_data];
            data = tmp.view();
            additional_header = header.to_string();
        } else {
            data = filtered_data;
        }
        let num_cols = data.ncols() - 1;
        let data: Vec<String> = data
            .indexed_iter()
            .map(|((_, j), d)| -> String {
                if j < num_cols {
                    format!("{:.6}", d) + "\t"
                } else {
                    format!("{:.6}", d) + "\n"
                }
            })
            .collect();
        let mut file = std::fs::File::create(file_name).expect("Error opening writing file");
        write!(file, "{}{}\n", additional_header, self.header).expect("Unable to write data");
        write!(file, "{}", data.join("")).expect("Unable to write data");
    }
}

#[derive(Debug, Clone)]
pub struct ObjectInfo {
    pub name: String,
    pub obj_type: ObjectType,
    pub index: usize,
    pub stored_data: StoreData,
}

impl ObjectInfo {
    pub fn new(
        name: String,
        obj_type: ObjectType,
        index: usize,
        stored_data: StoreData,
    ) -> ObjectInfo {
        ObjectInfo {
            name,
            obj_type,
            index,
            stored_data,
        }
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
