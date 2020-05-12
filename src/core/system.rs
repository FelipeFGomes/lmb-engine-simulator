use crate::engine::engine::Engine;
use crate::ObjectType;
use crate::{BasicProperties, FlowRatio};
use ndarray::*;
use std::f64::consts::PI;

// Core Traits
use crate::connector::conn_core::Connector;
use crate::one_dim::one_core::OneDim;
use crate::zero_dim::zero_core::ZeroDim;

pub struct System {
    _objs_info: Vec<(String, ObjectType, usize)>,
    engine: Option<Engine>,
    zero_dim: Vec<Box<dyn ZeroDim>>,
    one_dim: Vec<Box<dyn OneDim>>,
    connector: Vec<Box<dyn Connector>>,
    zero_dim_connectors_index: Vec<Vec<usize>>,
    one_dim_connectors_index: Vec<Vec<usize>>,
    connector_objects_index: Vec<Vec<(ObjectType, usize)>>,
    cycle_start: usize,
    iterations_counter: usize,
    time: Array2<f64>,
}

impl System {
    pub fn new(
        _objs_info: Vec<(String, ObjectType, usize)>,
        engine: Option<Engine>,
        zero_dim: Vec<Box<dyn ZeroDim>>,
        one_dim: Vec<Box<dyn OneDim>>,
        connector: Vec<Box<dyn Connector>>,
    ) -> Result<System, String> {
        let mut system = System {
            _objs_info,
            engine,
            zero_dim,
            one_dim,
            connector,
            zero_dim_connectors_index: Vec::new(),
            one_dim_connectors_index: Vec::new(),
            connector_objects_index: Vec::new(),
            cycle_start: 0,
            iterations_counter: 0,
            time: Array::from_elem((500000, 1), 0.),
        };

        system = system.setup_indexes()?;

        Ok(system)
    }
    pub fn advance<'a>(&'a mut self, dt: f64) -> &'a mut Self {
        // Advancing ZeroDim objects
        self.zero_dim.iter_mut().for_each(|zd| zd.advance(dt));

        // Update connectors `flow_ratio`: require `basic_properties` of DimElements
        for (connector, obj_i) in self
            .connector
            .iter_mut()
            .zip(self.connector_objects_index.iter())
        {
            let mut basic_properties: Vec<BasicProperties> =
                Vec::with_capacity(connector.connecting().len());
            for (obj_type, obj_index) in obj_i.iter() {
                match obj_type {
                    ObjectType::ZeroDim => {
                        basic_properties.push(self.zero_dim[*obj_index].get_state())
                    }
                    _ => panic!("Error at System.advance method\n Unknown type!"),
                }
            }
            connector.update_flow_ratio(basic_properties, dt);
        }

        // Update ZeroDim objects: require `flow_ratio` from connectors
        for (zero_dim, conn_index_list) in self
            .zero_dim
            .iter_mut()
            .zip(self.zero_dim_connectors_index.iter())
        {
            let mut sum_flow_ratio = FlowRatio::new();
            for conn_index in conn_index_list.iter() {
                let flow_ratio = match self.connector[*conn_index].get_flow_ratio(zero_dim.name()) {
                    Ok(flow) => flow,
                    Err(err) => {
                        println!("Error at 'advance':\n {}", err);
                        std::process::exit(1)
                    }
                };
                sum_flow_ratio = &sum_flow_ratio + flow_ratio;
            }
            zero_dim.update_flow_ratio(sum_flow_ratio);
        }

        // Updating OneDim objects

        self
    }

    pub fn advance_to_steady_state<'a>(&'a mut self) -> &'a mut Self {
        let mut time = 0.0;
        let mut cycle = 0;
        let mut total_angle = 0.0;
        let max_cycles = 10;
        let max_time = 5.0; // sec
        let step: f64;
        let step_calculator: Box<dyn Fn(&System) -> f64>;

        // if no `one_dim` elements exist step is constant (estimated by engine or assumed arbitrarily),
        // otherwise the step is determined by `one_dim` objects with CFL condition
        if self.one_dim.is_empty() {
            // estimate by engine
            match &self.engine {
                Some(engine) => {
                    let d_angle = 0.1;
                    step = (d_angle * PI / 180.0) / engine.sec_to_rad();
                    step_calculator = Box::new(|_: &System| -> f64 { step })
                }
                None => {
                    println!("WARNING: At `advance_to_steady_state`!");
                    println!("WARNING: No one-dimensional nor engine objects were found!");
                    println!("WARNING: Time step assumed as 0.1ms");
                    step = 1e-4;
                    step_calculator = Box::new(|_: &System| -> f64 { step })
                }
            }
        } else {
            step_calculator = Box::new(System::get_time_step);
        }

        // stop conditions are: cycle > max_cycles or time > max_time
        loop {
            let step = step_calculator(self);
            match &self.engine {
                Some(e) => {
                    total_angle += step * e.sec_to_rad();
                    if total_angle > 4.0 * PI {
                        total_angle -= 4.0 * PI;
                        cycle += 1;
                        if cycle > max_cycles {
                            break;
                        }
                        self.cycle_start = self.iterations_counter;
                    }
                }
                None => {}
            }
            if time > max_time {
                break;
            }
            self.advance(step);
            time += step;
            self.time[[self.iterations_counter, 0]] = time;
            self.iterations_counter += 1;
        }
        self
    }

    pub fn write_to_file(&self, file_name: &str, obj_name: &str) {
        if let Some((_, obj_type, index)) =
            self._objs_info.iter().find(|(name, _, _)| name == obj_name)
        {
            let time = self
                .time
                .slice(s![self.cycle_start..self.iterations_counter, ..]);
            match obj_type {
                ObjectType::ZeroDim => self.zero_dim[*index].write_to_file(
                    file_name,
                    (self.cycle_start, self.iterations_counter),
                    Some(("time [s]\t".to_string(), time)),
                ),
                ObjectType::OneDim => self.one_dim[*index].write_to_file(
                    file_name,
                    (self.cycle_start, self.iterations_counter),
                    None,
                ),
                ObjectType::Connector => self.connector[*index].write_to_file(
                    file_name,
                    (self.cycle_start, self.iterations_counter),
                    Some(("time [s]\t".to_string(), time)),
                ),
            }
        } else {
            println!("Error at 'write_to_file'");
            println!(" {} was not found", obj_name);
            std::process::exit(1)
        }
    }

    pub fn _run_obj_tests(&self, obj_name: &str) {
        if !self._objs_info.iter().any(|(name,_,_)| name == obj_name) {
            println!(
                "Error at: _run_obj_tests:\n Object `{}` not found!",
                obj_name
            );
            std::process::exit(1);
        }
    }

    fn setup_indexes(mut self) -> Result<Self, String> {
        let mut zero_dim_connectors_index: Vec<Vec<usize>> = Vec::new();
        for zero in self.zero_dim.iter() {
            let indexes = match self.get_connectors_index(zero.name().to_string()) {
                Ok(x) => x,
                Err(err) => {
                    return Err(format!("{}", err));
                }
            };
            zero_dim_connectors_index.push(indexes);
        }

        let mut one_dim_connectors_index: Vec<Vec<usize>> = Vec::new();
        for one in self.one_dim.iter() {
            let indexes = match self.get_connectors_index(one.name().to_string()) {
                Ok(x) => x,
                Err(err) => {
                    return Err(format!("{}", err));
                }
            };
            one_dim_connectors_index.push(indexes);
        }

        let mut connector_objects_index: Vec<Vec<(ObjectType, usize)>> = Vec::new();
        for conn in self.connector.iter() {
            let indexes = match self.get_dim_obj_index(conn.connecting()) {
                Ok(x) => x,
                Err(err) => {
                    return Err(format!("{}", err));
                }
            };
            connector_objects_index.push(indexes);
        }
        self.zero_dim_connectors_index = zero_dim_connectors_index;
        self.one_dim_connectors_index = one_dim_connectors_index;
        self.connector_objects_index = connector_objects_index;
        Ok(self)
    }

    fn get_connectors_index(&self, obj_name: String) -> Result<Vec<usize>, String> {
        let mut conn_indexes: Vec<usize> = Vec::new();
        for (index, conn) in self.connector.iter().enumerate() {
            if conn.connecting().iter().any(|name| *name == obj_name) {
                conn_indexes.push(index);
            }
        }
        if conn_indexes.is_empty() {
            let msg = format!("Object '{}' is not connected to anything.", obj_name);
            return Err(msg);
        }
        Ok(conn_indexes)
    }

    fn get_dim_obj_index(
        &self,
        conn_list: &Vec<String>,
    ) -> Result<Vec<(ObjectType, usize)>, String> {
        let mut objs_indexes: Vec<(ObjectType, usize)> = Vec::new();
        for obj in conn_list.iter() {
            if self.zero_dim.iter().any(|z| z.name() == obj) {
                let index = self.zero_dim.iter().position(|z| z.name() == obj).unwrap();
                objs_indexes.push((ObjectType::ZeroDim, index));
            } else if self.one_dim.iter().any(|z| z.name() == obj) {
                let index = self.zero_dim.iter().position(|z| z.name() == obj).unwrap();
                objs_indexes.push((ObjectType::OneDim, index));
            } else {
                let msg = format!("Object '{}' was not found in the arrays", obj);
                return Err(msg);
            }
        }
        Ok(objs_indexes)
    }

    fn get_time_step(_system: &System) -> f64 {
        println!("using get_time_step function!");
        0.0
    }

    pub fn print_state(&self) {
        for obj in self.zero_dim.iter() {
            println!("{}", obj.get_state());
        }

        for (zero_dim, conn_index_list) in self
            .zero_dim
            .iter()
            .zip(self.zero_dim_connectors_index.iter())
        {
            let mut sum_flow_ratio = FlowRatio::new();
            for conn_index in conn_index_list.iter() {
                let flow_ratio = match self.connector[*conn_index].get_flow_ratio(zero_dim.name()) {
                    Ok(flow) => flow,
                    Err(err) => {
                        println!("Error at 'advance':\n {}", err);
                        std::process::exit(1)
                    }
                };
                sum_flow_ratio = &sum_flow_ratio + flow_ratio;
            }
            println!("`{}`, {:?}", zero_dim.name(), sum_flow_ratio);
        }
    }
}
