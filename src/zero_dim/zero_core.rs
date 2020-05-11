use crate::{BasicProperties, FlowRatio};
use ndarray::*;

pub trait ZeroDim {
    fn name<'a>(&'a self) -> &'a str;
    fn get_state(&self) -> BasicProperties;
    fn advance(&mut self, dt: f64);
    fn update_flow_ratio(&mut self, total_flow_ratio: FlowRatio);
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
