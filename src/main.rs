use lmb::Gas;
use lmb_engine_simulator as lmb;

fn main() {
    let mut gas_ambient = Gas::new("air.json");
    gas_ambient.TPX(293.0, 2.0*101325.0, "O2:0.21, N2:0.79");
    let mut gas_chamber = Gas::new("air.json");
    gas_chamber.TPX(293.0, 101325.0, "O2:0.21, N2:0.79");
    let mut builder = lmb::SystemBuilder::new();
    builder
        .add_environment("ambient", &gas_ambient)
        .add_reservoir("chamber", 500.0, &gas_chamber)
        .add_orifice("orifice", 50.0, 0.9, vec!["ambient", "chamber"]);

    let mut system = builder.build_system();
    
    // Calculating
    system.advance_to_steady_state();
    
    // Writting data
    system.write_to_file("chamber.txt", "chamber", None);
    system.write_to_file("orifice.txt", "orifice", None);
}
