//! # lmb_engine_simulator
//!
//! The `lmb_engine_simulator` crate provides an easy way to simulate engines and 1D gas dynamics.

pub mod reaction;

// Re-exporting
pub use crate::reaction::gas::Gas;

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
