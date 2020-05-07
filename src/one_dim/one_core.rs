pub trait OneDim { 
    fn name<'a>(&'a self) -> &'a str;
    fn advance(&mut self, dt: f64) {
        println!("object {} is advancing its state with dt={}", self.name(), dt);
    }
}