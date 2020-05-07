use crate::ObjectType;
use crate::{BasicProperties, FlowRatio};
// Core Traits
use crate::connector::conn_core::Connector;
use crate::one_dim::one_core::OneDim;
use crate::zero_dim::zero_core::ZeroDim;

pub struct System {
    _objs_name: Vec<String>,
    zero_dim: Vec<Box<dyn ZeroDim>>,
    one_dim: Vec<Box<dyn OneDim>>,
    connector: Vec<Box<dyn Connector>>,
    zero_dim_connectors_index: Vec<Vec<usize>>,
    one_dim_connectors_index: Vec<Vec<usize>>,
    connector_objects_index: Vec<Vec<(ObjectType, usize)>>,
    pub stored_data: Vec<String>,
}

impl System {
    pub fn new(
        _objs_name: Vec<String>,
        zero_dim: Vec<Box<dyn ZeroDim>>,
        one_dim: Vec<Box<dyn OneDim>>,
        connector: Vec<Box<dyn Connector>>,
    ) -> Result<System, String> {
        let mut system = System {
            _objs_name,
            zero_dim,
            one_dim,
            connector,
            zero_dim_connectors_index: Vec::new(),
            one_dim_connectors_index: Vec::new(),
            connector_objects_index: Vec::new(),
            stored_data: Vec::new(),
        };
        system = system.setup_indexes()?;
        Ok(system)
    }
    pub fn advance<'a>(&'a mut self, dt: f64) -> &'a mut Self {
        // Advancing ZeroDim objects
        self.zero_dim.iter_mut().for_each(|zd| zd.advance(dt));

        // Update connectors `flow_ratio`: require `basic_properties` of DimElements
        for (connector, obj_i) in self.connector.iter_mut().zip(self.connector_objects_index.iter()) 
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
        for (zero_dim, conn_index_list) in self.zero_dim.iter_mut().zip(self.zero_dim_connectors_index.iter()) 
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

        // Storing the data
        let cyl = match self.zero_dim.iter().find(|zd| zd.name() == "cyl_1") {
            Some(c) => c,
            None => {
                println!("'cyl_1' not found");
                std::process::exit(1)},
        };
        self.stored_data.push(cyl._get_main_properties());

        // Updating OneDim objects

        self
    }

    pub fn _run_obj_tests(&self, obj_name: &str) {
        if !self._objs_name.iter().any(|name| name == obj_name) {
            println!("Error at: _run_obj_tests:\n Object `{}` not found!", obj_name);
            std::process::exit(1);
        }
    }

    fn setup_indexes(mut self) -> Result<Self, String> {
        let mut zero_dim_connectors_index: Vec<Vec<usize>> = Vec::new();
        for zero in self.zero_dim.iter() {
            let indexes = match self.get_connectors_index(zero.name().to_string()) {
                Ok(x) => x,
                Err(err) => { return Err(format!("{}", err)); }
            };
            zero_dim_connectors_index.push(indexes);
        }

        let mut one_dim_connectors_index: Vec<Vec<usize>> = Vec::new();
        for one in self.one_dim.iter() {
            let indexes = match self.get_connectors_index(one.name().to_string()) {
                Ok(x) => x,
                Err(err) => { return Err(format!("{}", err)); }
            };
            one_dim_connectors_index.push(indexes);
        }

        let mut connector_objects_index: Vec<Vec<(ObjectType, usize)>> = Vec::new();
        for conn in self.connector.iter() {
            let indexes = match self.get_dim_obj_index(conn.connecting()) {
                Ok(x) => x,
                Err(err) => { return Err(format!("{}", err)); }
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
                objs_indexes.push( (ObjectType::ZeroDim, index) );
            } else if self.one_dim.iter().any(|z| z.name() == obj) {
                let index = self.zero_dim.iter().position(|z| z.name() == obj).unwrap();
                objs_indexes.push( (ObjectType::OneDim, index) );
            }
            else {
                let msg = format!("Object '{}' was not found in the arrays", obj);
                return Err(msg);
            }
        }
        Ok(objs_indexes)
    }

    pub fn print_state(&self) {
        for obj in self.zero_dim.iter() {
            println!("{}", obj.get_state());
        }

        for (zero_dim, conn_index_list) in self.zero_dim.iter().zip(self.zero_dim_connectors_index.iter()) 
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
