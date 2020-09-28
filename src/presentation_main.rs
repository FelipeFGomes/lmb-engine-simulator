use lmb::Gas;
use lmb_engine_simulator as lmb;
use lmb::combustion::{Combustion, TwoZoneCombustion};

fn main() {

    let mut gas = Gas::new("air.json");
    gas.TP(293.0, 101325.0);
    let mut gas_intake = Gas::new("air.json");
    gas_intake.TP(298.0, 101325.0);
    let mut gas_exhaust = Gas::new("air.json");
    gas_exhaust.TPX(700.0, 101325.0, "N2:0.662586, H2O:0.202449, CO2:0.134965");
    
    let mut builder = lmb::SystemBuilder::new();
    builder.add_engine("engine.json", &gas)
        .add_environment("ambient", &gas)
        .add_reservoir("int_plenum", 250.0, &gas_intake)
        .add_reservoir("int_port", 2.04, &gas_intake)
        .add_reservoir("exh_port", 2.04, &gas_exhaust)
        .add_reservoir("exh_plenum_1", 95.0, &gas_exhaust)
        .add_reservoir("exh_plenum_2", 158.0, &gas_exhaust)
        .add_connection_between_0D("amb 1 -> int_plenum", 9.5, 0.92, vec!["ambient", "int_plenum"])
        .add_connection_between_0D("amb 2 -> int_plenum", 9.5, 0.92, vec!["ambient", "int_plenum"])
        .add_connection_between_0D("int_plenum -> int_port", 6.5, 0.80, vec!["int_plenum", "int_port"])
        .add_connection_between_0D("amb -> exh_plenum_2", 8.9, 0.80, vec!["ambient", "exh_plenum_2"])
        .add_connection_between_0D("exh_plenum_2 -> exh_plenum_1", 12.7, 0.80, vec!["exh_plenum_2", "exh_plenum_1"])
        .add_connection_between_0D("exh_plenum_1 -> exh_port", 10.8, 0.80, vec!["exh_plenum_1", "exh_port"])
        .connect_from_to("valve_int", "int_port")
        .connect_from_to("valve_exh", "exh_port");
    
    let mut system = builder.build_system();  
    
    
    // Calculating
    let speeds = vec![5000.0, 5500.0, 6000.0, 6500.0, 7000.0, 7500.0, 8000.0, 8500.0, 9000.0];
    let folder_name = "./Ryobi_26_results/";
    for speed in speeds {
        system.engine_mut().unwrap().set_speed(speed);
        system.advance_to_steady_state();

        // Writting data
        let fp = format!("{}{:.0}_", folder_name, speed);
        system.write_to_file( &(fp.clone() + "cylinder.txt"), "cyl_1", None);
        system.write_to_file( &(fp.clone() + "valve_int.txt"), "valve_int", None);
        system.write_to_file( &(fp.clone() + "valve_exh.txt"), "valve_exh", None);
        system.write_to_file( &(fp.clone() + "int_port.txt"), "int_port", None);
        system.write_to_file( &(fp.clone() + "exh_port.txt"), "exh_port", None);
    
    }
    let fp = format!("{}engine_performance.txt", folder_name);
    system.engine().unwrap().write_performance_to(&fp);

    /*
    // Variable parameters
    let afr = 0.90;
    let fuel = system.engine().unwrap().injector().unwrap().fuel().clone();
    let air_comp = gas_intake.clone();
    let speed = vec![1200.0, 2400.0, 3600.0, 4800.0, 6000.0];
    let comb_ini = vec![712.0, 708.0, 710.0, 709.0, 706.0];
    let a = vec![6.49, 10.32, 8.2, 7.62, 6.02];
    let m = vec![1.72, 1.93, 1.82, 1.69, 1.64];
    let comb_duration = vec![32.0, 41.0, 39.0, 44.0, 44.0];

    for i in 0..speed.len() {
        // creating a Combustion model
        let wiebe = lmb::combustion::WiebeFunction::new(a[i], m[i], comb_duration[i]);
        let comb: Box<dyn Combustion> = Box::new( TwoZoneCombustion::new(comb_ini[i], afr, wiebe, &gas, &fuel, &air_comp) );
        system.engine_mut().unwrap().set_combustion_model(&comb);
        system.engine_mut().unwrap().set_air_fuel_ratio(afr);
        system.engine_mut().unwrap().set_speed(speed[i]);
        system.advance_to_steady_state();

        // Writting data
        let folder_name = format!("./results/{:.0}_", speed[i]);
        system.write_to_file( &(folder_name.clone() + "cylinder.txt"), "cyl_1", None);
        system.write_to_file( &(folder_name.clone() + "valve_int.txt"), "valve_int", None);
        system.write_to_file( &(folder_name.clone() + "valve_exh.txt"), "valve_exh", None);
        system.write_to_file( &(folder_name.clone() + "int_port.txt"), "int_port", None);
        system.write_to_file( &(folder_name.clone() + "exh_port.txt"), "exh_port", None);
    }
    
    system.engine().unwrap().write_performance_to("./results/engine_performance.txt");
    */
}
