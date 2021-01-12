//! # lmb_engine_simulator
//!
//! The `lmb_engine_simulator` crate provides an easy way to simulate engines and 1D gas dynamics.
//!
//! This library employed the **Builder Pattern** so the user feels as she/he is acctually building
//! an engine system. To construct (build) the system, the struct `SystemBuilder`
//! is used to add the desired components. After finishing building, the method `build_system()` can
//! be used to return an object `System` which is used to solve the components numerically.
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
//! ## Examples
//! ```
//! let mut gas_ambient = Gas::new("air.json");
//! gas_ambient.TPX(293.0, 2.0*101325.0, "O2:0.21, N2:0.79");
//! let mut gas_chamber = Gas::new("air.json");
//! gas_chamber.TPX(293.0, 101325.0, "O2:0.21, N2:0.79");
//! let mut builder = lmb::SystemBuilder::new();
//! builder
//!     .add_environment("ambient", &gas_ambient)
//!     .add_reservoir("chamber", 500.0, &gas_chamber)
//!     .add_orifice("orifice", 50.0, 0.9, vec!["ambient", "chamber"]);
//!
//! let mut system = builder.build_system();
//!    
//! // Calculating
//! system.advance_to_steady_state();
//!    
//! // Writting data
//! system.write_to_file("chamber.txt", "chamber", None);
//! system.write_to_file("orifice.txt", "orifice", None);
//! ```
//! 
//! ## Simulation with engines
//! Simulations with engines are a bit more complex and are treated in details at [Engine Examples](doc/Ryobi_26cm3_engine/index.html)

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
