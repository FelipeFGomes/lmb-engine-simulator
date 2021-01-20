use lmb::Gas;
use lmb_engine_simulator as lmb;

fn main() {
    let gas_intake = Gas::new("air.json");
    let mut gas_exhaust = Gas::new("air.json");
    gas_exhaust.TPX(500.0, 101325.0, "N2:0.662586, H2O:0.202449, CO2:0.134965");
    let mut builder = lmb::SystemBuilder::new();
    builder
        .add_engine("engine.json", &gas_intake)
        .add_environment("intake_port", &gas_intake)
        .add_environment("exhaust_port", &gas_exhaust)
        .connect_from_to("valve_int", "intake_port")
        .connect_from_to("valve_exh", "exhaust_port");

    let mut system = builder.build_system();
    
    // Calculating
    system.advance_to_steady_state();
    
    // Writting data
    system.write_to_file("cylinder.txt", "cyl_1", None);
    system.write_to_file("intake_valve.txt", "valve_int", None);
    system.write_to_file("exhaust_valve.txt", "valve_exh", None);
}
