use ndarray::*;
use crate::{BasicProperties, FlowRatio};

// Super Traits
pub trait ZeroD: ZeroDim + SaveData {}
pub trait OneD: OneDim + SaveData {}
pub trait Conn: Connector + SaveData {}

pub trait ZeroDim {
    fn name<'a>(&'a self) -> &'a str;
    fn get_state(&self) -> BasicProperties;
    fn advance(&mut self, dt: f64);
    fn update_flow_ratio(&mut self, total_flow_ratio: Vec<(&str, &FlowRatio)>);
}

pub trait OneDim { 
    fn name<'a>(&'a self) -> &'a str;
    fn advance(&mut self, dt: f64) {
        println!("object {} is advancing its state with dt={}", self.name(), dt);
    }
    fn write_to_file(
        &self,
        _file_name: &str,
        _range: (usize, usize),
        _extra_data: Option<(String, ArrayView2<f64>)>,
    ) {
        println!(
            "`write_to_file` method has not been implemented to `{}`",
            self.name()
        );
    }
}

pub trait Connector {
    fn name<'a>(&'a self) -> &'a str;
    fn connecting<'a>(&'a self) -> &'a Vec<String>;
    fn connect_to(&mut self, elem_name: &str) -> Result<(),String>;
    fn update_flow_ratio(&mut self, info: Vec<BasicProperties>, _step: f64) {
        println!("updating {} with information {:?}", self.name(), info);
    }
    fn get_flow_ratio<'a>(&'a self, elem_name: &str) -> Result<&'a FlowRatio, String>;
}

pub trait SaveData {
    fn get_headers(&self) -> String;
    fn num_storable_variables(&self) -> usize;
    fn get_storable_data(&self) -> Array1<f64>;
}