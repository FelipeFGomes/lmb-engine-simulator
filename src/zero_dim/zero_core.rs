use crate::{BasicProperties, FlowRatio};

pub trait ZeroDim {
    fn name<'a>(&'a self) -> &'a str;
    fn get_state(&self) -> BasicProperties;
    fn advance(&mut self, dt: f64);
    fn update_flow_ratio(&mut self, total_flow_ratio: FlowRatio);
    fn _get_main_properties(&self) -> String;
}
