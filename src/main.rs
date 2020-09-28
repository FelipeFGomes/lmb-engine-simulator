use lmb::Gas;
use lmb_engine_simulator as lmb;

fn main() {

    let mut gas = Gas::new("air.json");
    gas.TP(293.0, 101325.0);
    let mut gas_intake = Gas::new("air.json");
    gas_intake.TP(298.0, 101325.0);
    let mut gas_exhaust = Gas::new("air.json");
    gas_exhaust.TPX(700.0, 101325.0, "N2:0.662586, H2O:0.202449, CO2:0.134965");
    
    let mut builder = lmb::SystemBuilder::new();
    builder.add_engine("engine.json", &gas)
        .add_environment("int_port", &gas_intake)
        .add_environment("exh_port", &gas_exhaust)
        .connect_from_to("valve_int", "int_port")
        .connect_from_to("valve_exh", "exh_port");
    
    let mut system = builder.build_system();  
    
    // Calculating
    system.advance_to_steady_state();

    // Writting data
    system.write_to_file("cylinder.txt", "cyl_1", None);
    system.write_to_file("valve_int.txt", "valve_int", None);
    system.write_to_file("valve_exh.txt", "valve_exh", None);
}
