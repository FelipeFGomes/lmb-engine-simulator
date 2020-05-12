//! # lmb_engine_simulator
//!
//! The `lmb_engine_simulator` crate provides an easy way to simulate engines and 1D gas dynamics.
//! 
//! This library employed the **Builder Pattern** so the user feels as she/he is acctually building
//! an engine system. To construct (build) the system, the struct `SystemBuilder` 
//! is used to add the desired components. After finishing building, the method `build_system()` can
//! be used to return an object `System` which is used to solve the components numerically. 
//! 
//! # Example
//! ```
//! let gas = Gas::new("air.json");
//! let mut builder = lmb::SystemBuilder::new();
//! builder
//!     .add_reservoir("res_1", &gas, 0.5)      // object name, gas inside, volume
//!     .add_environment("exhaust_env", &gas);  // object name, gas inside
//! let system = builder.build_system();
//! ```

use std::ops::Add;
use std::io::Write;
use ndarray::*;

pub mod base;
pub mod connector;
pub mod core;
pub mod engine;
pub mod numerics;
pub mod one_dim;
pub mod reaction;
pub mod zero_dim;

// Re-exporting
pub use crate::core::system_builder::SystemBuilder;
pub use crate::reaction::gas::Gas;

// Object's type. 'String' meant to store the name of the object
#[derive(Debug, PartialEq)]
pub enum ObjectType {
    ZeroDim,
    OneDim,
    Connector,
}

#[derive(Debug)]
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

#[derive(Debug)]
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

#[derive(Debug)]
pub struct StoreData {
    header: String,
    data: Array2<f64>,
    last_index: usize,
}

impl StoreData {
    pub fn new(header: &str, num_variables: usize) -> StoreData {
        StoreData {
            header: header.to_string(),
            data: Array::from_elem( (500000, num_variables), 0.),
            last_index: 0,
        }
    }

    pub fn add_data(&mut self, data: Array1<f64>) {
        self.data.row_mut(self.last_index).assign(&data);
        self.last_index += 1;
    }

    pub fn reset_index<'a>(&'a mut self) -> &'a mut Self {
        self.last_index = 0;
        self
    }

    /// Write the stored data in `data` limited by the index `range`to a file,
    ///  the first line is the content in `header`. 
    pub fn write_to_file( &self,
        file_name: &str,
        range: (usize, usize),
        additional_data: Option<(String, ArrayView2<f64>)>) 
        {
            let data: ArrayView2<f64>;
            let tmp: Array2<f64>;
            let mut additional_header = String::from("");
            let filtered_data = self.data.slice(s![range.0..range.1, ..]);
            if let Some((header, add)) = additional_data {
                if filtered_data.nrows() != add.nrows() {
                    println!("`additional_data` must have the same number of rows as the writable data");
                    println!(" `additional_data`: {}, writable data: {}", add.len(), filtered_data.len());
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

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
