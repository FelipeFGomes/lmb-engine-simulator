use lmb_engine_simulator as lmb;
use lmb::Gas;

fn main() {
    let gas = Gas::new("air_lmb.json");
    println!("{:?}", gas);
}