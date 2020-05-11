use crate::{BasicProperties, FlowRatio};
use ndarray::*;

pub trait Connector {
    fn name<'a>(&'a self) -> &'a str;
    fn connecting<'a>(&'a self) -> &'a Vec<String>;
    fn connect_to(&mut self, elem_name: &str) -> Result<(),String>; 
    fn update_flow_ratio(&mut self, info: Vec<BasicProperties>, _step: f64) {
        println!("updating {} with information {:?}", self.name(), info);
    }
    fn get_flow_ratio<'a>(&'a self, elem_name: &str) -> Result<&'a FlowRatio, String>;
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