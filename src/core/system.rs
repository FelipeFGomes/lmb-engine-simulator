use crate::engine::engine::Engine;
use crate::{BasicProperties, FlowRatio};
use crate::{ObjectInfo, ObjectType};
use ndarray::*;
use std::f64::consts::PI;
use std::time::Instant;

// Core Traits
// use crate::core::traits::Connector;
// use crate::core::traits::OneDim;
use crate::core::traits::SaveData;
use crate::core::traits::ZeroDim;

// Types
type IndexOutput = Result<Vec<(ObjectType, usize)>, String>;

// Super Traits
use crate::core::traits::Conn;
use crate::core::traits::OneD;
use crate::core::traits::ZeroD;

pub struct System {
    objs_info: Vec<ObjectInfo>,
    engine: Option<Engine>,
    zero_dim: Vec<Box<dyn ZeroD>>,
    one_dim: Vec<Box<dyn OneD>>,
    connector: Vec<Box<dyn Conn>>,
    engine_connectors_index: Vec<Vec<(ObjectType, usize)>>,
    zero_dim_connectors_index: Vec<Vec<(ObjectType, usize)>>,
    one_dim_connectors_index: Vec<Vec<(ObjectType, usize)>>,
    connector_objects_index: Vec<Vec<(ObjectType, usize)>>,
    cycle_start: usize,
    iterations_counter: usize,
    time: Array2<f64>,
}

impl System {
    pub fn new(
        objs_info: Vec<ObjectInfo>,
        engine: Option<Engine>,
        zero_dim: Vec<Box<dyn ZeroD>>,
        one_dim: Vec<Box<dyn OneD>>,
        connector: Vec<Box<dyn Conn>>,
    ) -> Result<System, String> {
        let mut system = System {
            objs_info,
            engine,
            zero_dim,
            one_dim,
            connector,
            engine_connectors_index: Vec::new(),
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

    /// Advance all objects in System by `dt`. The state of the objects is stored.
    pub fn advance<'a>(&'a mut self, dt: f64) -> &'a mut Self {
        // Advancing ZeroDim objects
        self.zero_dim.iter_mut().for_each(|zd| zd.advance(dt));

        // Advancing OneDim objects

        // Advacing Engine
        if let Some(eng) = &mut self.engine {
            eng.advance(dt);
        }

        // check if masses less than zero

        // Update connectors `flow_ratio`: require `basic_properties` of DimElements
        for (connector, obj_index_list) in self
            .connector
            .iter_mut()
            .zip(self.connector_objects_index.iter())
        {
            let mut basic_properties: Vec<BasicProperties> =
                Vec::with_capacity(connector.connecting().len());
            for (obj_type, i) in obj_index_list.iter() {
                match obj_type {
                    ObjectType::ZeroDim => basic_properties.push(self.zero_dim[*i].get_state()),
                    ObjectType::Cylinder => basic_properties
                        .push(self.engine.as_ref().unwrap().cylinders()[*i].get_state()),
                    _ => panic!("Error at `System::advance()`\n Object of unknown type!"),
                }
            }
            connector.update_flow_ratio(basic_properties, dt);
        }

        // Update ZeroDim objects: require `flow_ratio` and `name `from connectors
        for (zero_dim, conn_index_list) in self
            .zero_dim
            .iter_mut()
            .zip(self.zero_dim_connectors_index.iter())
        {
            let mut total_flow_ratio: Vec<(&str, &FlowRatio)> =
                Vec::with_capacity(conn_index_list.len());
            for (_, i) in conn_index_list.iter() {
                let (flow_ratio, name) = match self.connector[*i].get_flow_ratio(zero_dim.name()) {
                    Ok(flow) => (flow, self.connector[*i].name()),
                    Err(err) => {
                        println!("Error at 'System::advance()': \n {}", err);
                        std::process::exit(1);
                    }
                };

                total_flow_ratio.push((name, flow_ratio));
            }
            zero_dim.update_flow_ratio(total_flow_ratio);
        }

        // Updating OneDim objects

        // Updating Engine objects
        if let Some(engine) = &mut self.engine {
            let mut cylinders_total_flow: Vec<Vec<(&str, &FlowRatio)>> =
                Vec::with_capacity(engine.cylinders().len());
            for (cyl, conn_index_list) in engine
                .cylinders()
                .iter()
                .zip(self.engine_connectors_index.iter())
            {
                let mut total_flow_ratio: Vec<(&str, &FlowRatio)> =
                    Vec::with_capacity(conn_index_list.len());
                for (_, i) in conn_index_list.iter() {
                    let (flow_ratio, name) = match self.connector[*i].get_flow_ratio(cyl.name()) {
                        Ok(flow) => (flow, self.connector[*i].name()),
                        Err(err) => {
                            println!("Error at 'System::advance()': \n {}", err);
                            std::process::exit(1);
                        }
                    };
                    total_flow_ratio.push((name, flow_ratio));
                }
                cylinders_total_flow.push(total_flow_ratio);
            }
            engine.update_cylinders_flow_ratio(cylinders_total_flow)
        }

        for data in self.objs_info.iter_mut() {
            match data.obj_type {
                ObjectType::ZeroDim => {
                    data.stored_data
                        .add_data(self.zero_dim[data.index].get_storable_data());
                }
                ObjectType::OneDim => {
                    data.stored_data
                        .add_data(self.one_dim[data.index].get_storable_data());
                }
                ObjectType::Connector => {
                    data.stored_data
                        .add_data(self.connector[data.index].get_storable_data());
                }
                ObjectType::Cylinder => {
                    if let Some(eng) = &self.engine {
                        data.stored_data
                            .add_data(eng.cylinders()[data.index].get_storable_data());
                    }
                }
            }
        }

        self
    }

    /// Advance all objects in System until steady state or until a max number of iterations/cycles is reached.
    /// All the stored data of objects are reseted
    pub fn advance_to_steady_state<'a>(&'a mut self) -> &'a mut Self {
        let now = Instant::now(); // measuring time

        // reseting StoredData from all objects:
        self.objs_info
            .iter_mut()
            .for_each(|obj| obj.stored_data.reset_data());
        self.cycle_start = 0;
        self.iterations_counter = 0;
        self.time = Array::from_elem((500000, 1), 0.);

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
                    let d_angle = 0.10;
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

        println!("\n\t\tSystem successfully advanced to steady state!");
        println!("number of iterations: {}", self.iterations_counter);
        println!("time taken: {:?}", now.elapsed());

        if let Some(engine) = &mut self.engine {
            let mut pressure: Vec<ArrayView1<f64>> = Vec::new();
            let mut volume: Vec<ArrayView1<f64>> = Vec::new();
            let mut data: Vec<Array2<f64>> = Vec::new();
            let range = (self.cycle_start, self.iterations_counter);
            let p_index: usize = 1;
            let v_index: usize = 3;
            for info in &self.objs_info {
                match info.obj_type {
                    ObjectType::Cylinder => {
                        data.push(info.stored_data.get_data(range, vec![p_index, v_index]));
                    }
                    _ => {}
                }
            }
            data.iter_mut().for_each(|d| {
                pressure.push( d.column(0) );
                volume.push(d.column(1))
            });

            engine.calc_operational_param(
                pressure,
                volume,
            );

            println!("Engine performance:{}", engine.operat_param());
        }
        self
    }

    /// Write the stored data from a object, `obj_name`, into a file `file_name`.
    /// If `Engine` has been added, only the data stored in the last cycle will be written as default.
    /// `_range` can be used to set the first and last index of the writable data.
    pub fn write_to_file(&self, file_name: &str, obj_name: &str, _range: Option<(usize, usize)>) {
        let range: (usize, usize);
        if let Some(r) = _range {
            range = r;
        } else {
            range = (self.cycle_start, self.iterations_counter);
        }
        let time = self.time.slice(s![range.0..range.1, ..]);

        if let Some(obj_info) = self.objs_info.iter().find(|info| info.name == obj_name) {
            obj_info.stored_data.write_to_file(
                file_name,
                range,
                Some(("time [s]\t".to_string(), time)),
            )
        } else {
            println!("Error at 'write_to_file'");
            println!(" {} was not found", obj_name);
            std::process::exit(1)
        }
    }

    pub fn engine<'a>(&'a self) -> Option<&'a Engine> {
        match &self.engine {
            Some(eng) => Some(eng),
            None => None,
        }
    }

    pub fn engine_mut<'a>(&'a mut self) -> Option<&'a mut Engine> {
        match &mut self.engine {
            Some(eng) => Some(eng),
            None => None,
        }
    }

    fn setup_indexes(mut self) -> Result<Self, String> {
        // finding which connectors are connected 0D objects
        let mut zero_dim_connectors_index: Vec<Vec<(ObjectType, usize)>> = Vec::new();
        for zero in self.zero_dim.iter() {
            let indexes = match self.get_connectors_index(zero.name().to_string()) {
                Ok(x) => x,
                Err(err) => {
                    return Err(format!("{}", err));
                }
            };
            zero_dim_connectors_index.push(indexes);
        }

        // finding which connectors are connected 1D objects
        let mut one_dim_connectors_index: Vec<Vec<(ObjectType, usize)>> = Vec::new();
        for one in self.one_dim.iter() {
            let indexes = match self.get_connectors_index(one.name().to_string()) {
                Ok(x) => x,
                Err(err) => {
                    return Err(format!("{}", err));
                }
            };
            one_dim_connectors_index.push(indexes);
        }

        // finding which objects are connected to the connectors
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

        // finding which connectors are connected the engine's cylinders
        if self.engine.is_some() {
            let mut engine_connectors_index: Vec<Vec<(ObjectType, usize)>> = Vec::new();
            for cyl in self.engine.as_ref().unwrap().cylinders() {
                let indexes = match self.get_connectors_index(cyl.name().to_string()) {
                    Ok(x) => x,
                    Err(err) => {
                        return Err(format!("{}", err));
                    }
                };
                engine_connectors_index.push(indexes);
            }
            self.engine_connectors_index = engine_connectors_index;
        }

        self.zero_dim_connectors_index = zero_dim_connectors_index;
        self.one_dim_connectors_index = one_dim_connectors_index;
        self.connector_objects_index = connector_objects_index;
        Ok(self)
    }

    fn get_connectors_index(&self, obj_name: String) -> IndexOutput {
        let mut conn_indexes: Vec<(ObjectType, usize)> = Vec::new();
        // searching connector vector
        for (index, conn) in self.connector.iter().enumerate() {
            if conn.connecting().iter().any(|name| *name == obj_name) {
                conn_indexes.push((ObjectType::Connector, index));
            }
        }
        if conn_indexes.is_empty() {
            let msg = format!("Object '{}' is not connected to anything.", obj_name);
            return Err(msg);
        }
        Ok(conn_indexes)
    }

    fn get_dim_obj_index(&self, conn_list: &Vec<String>) -> IndexOutput {
        let mut objs_indexes: Vec<(ObjectType, usize)> = Vec::new();
        for obj in conn_list.iter() {
            if self.zero_dim.iter().any(|z| z.name() == obj) {
                let index = self.zero_dim.iter().position(|z| z.name() == obj).unwrap();
                objs_indexes.push((ObjectType::ZeroDim, index));
            } else if self.one_dim.iter().any(|z| z.name() == obj) {
                let index = self.zero_dim.iter().position(|z| z.name() == obj).unwrap();
                objs_indexes.push((ObjectType::OneDim, index));
            } else if self.engine.is_some() {
                let index = self
                    .engine
                    .as_ref()
                    .unwrap()
                    .cylinders()
                    .iter()
                    .position(|z| z.name() == obj)
                    .unwrap();
                objs_indexes.push((ObjectType::Cylinder, index));
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

    pub fn _store_composition_of(&mut self, _obj_name: &str) {
        // match self._objs_info.iter().find(|(name,_,_)| name == obj_name) {
        //     Some((name, obj_type, i)) => {
        //         match obj_type {
        //         ObjectType::ZeroDim => self.zero_dim[*index].write_to_file(
        //             file_name,
        //             (self.cycle_start, self.iterations_counter),
        //             Some(("time [s]\t".to_string(), time)),
        //         ),
        //         _ =>
        //     }
        //     },
        //     None => {},
        // }
    }
}
