use ndarray::*;

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