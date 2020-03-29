pub mod reaction;

// Re-exporting
pub use crate::reaction::gas::Gas;

pub fn hi() {
    println!("hi!");
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
