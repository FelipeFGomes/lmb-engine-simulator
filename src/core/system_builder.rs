use super::system::System;
use crate::connector;
use crate::one_dim;
use crate::reaction::gas::Gas;
use crate::zero_dim;
use crate::engine::engine::Engine;
// use crate::ObjectType;

// Core Traits
use connector::conn_core::Connector;
use one_dim::one_core::OneDim;
use zero_dim::zero_core::ZeroDim;

pub struct SystemBuilder {
    objs_name: Vec<String>,
    zero_dim: Vec<Box<dyn ZeroDim>>,
    one_dim: Vec<Box<dyn OneDim>>,
    connector: Vec<Box<dyn Connector>>,
}

impl SystemBuilder {
    pub fn new() -> SystemBuilder {
        SystemBuilder {
            objs_name: Vec::new(),
            zero_dim: Vec::new(),
            one_dim: Vec::new(),
            connector: Vec::new(),
        }
    }

    pub fn build_system(self) -> System {
        match System::new(self.objs_name, self.zero_dim, self.one_dim, self.connector) {
            Ok(s) => s,
            Err(err) => {
                println!("Error at 'SystemBuilder.build_system()':\n {}", err);
                std::process::exit(1)
            }
        }
    }

    pub fn add_engine<'a>(&'a mut self, file_name: &str, gas: &Gas) -> &'a mut Self {
        let engine = match Engine::new(file_name, gas) {
            Ok(e) => e,
            Err(err) => {
                println!("Error at 'add_engine':\n {}", err);
                std::process::exit(1)
            }
        };

        // Pushing cylinders to `zero_dim` vector
        for cylinder in engine.cylinders {
            if self.does_it_exist(cylinder.name()) {
                println!("Error at 'add_engine':");
                println!("Object with the same name already exists: {}", cylinder.name());
                std::process::exit(1)
            }
            self.objs_name.push(cylinder.name().to_string());
            self.zero_dim.push( Box::new(cylinder) );
        }

        // Pushing valves to `connector` vector
        for valve in engine.valves {
            if self.does_it_exist(valve.name()) {
                println!("Error at `add_engine`:");
                println!(" Object with the same name already exists: {}", valve.name());
                std::process::exit(1);
            }
            self.objs_name.push( valve.name().to_string() );
            self.connector.push( Box::new(valve) );
        }
        
        self
        
    }

    pub fn add_environment<'a>(&'a mut self, elem_name: &str, gas: &Gas) -> &'a mut Self {
        // checking if 'elem_name' already exists
        if self.does_it_exist(elem_name) {
            println!("Error at 'add_environment':");
            println!("Object with the same name already exists");
            std::process::exit(1)
        }
        // pushing environment
        let env = match zero_dim::environment::Environment::new(elem_name.to_string(), gas) {
            Ok(v) => v,
            Err(err) => {
                println!("Error at 'add_environment':\n {}", err);
                std::process::exit(1)
            }
        };
        self.zero_dim.push(Box::new(env));

        // adding to list of objects
        self.objs_name.push(elem_name.to_string());
        self
    }

    pub fn add_reservoir<'a>(
        &'a mut self,
        elem_name: &str,
        gas: &Gas,
        volume: f64,
    ) -> &'a mut Self {
        // checking if 'elem_name' already exists
        if self.does_it_exist(elem_name) {
            println!("Error at 'add_reservoir':");
            println!("Object with the same name already exists");
            std::process::exit(1)
        }
        // pushing reservoir
        let res = match zero_dim::reservoir::Reservoir::new(elem_name.to_string(), gas, volume) {
            Ok(v) => v,
            Err(err) => {
                println!("Error at 'add_reservoir':\n {}", err);
                std::process::exit(1)
            }
        };
        self.zero_dim.push(Box::new(res));

        // adding to list of objects
        self.objs_name.push(elem_name.to_string());
        self
    }

    pub fn connect_from_to<'a>(&'a mut self, connector: &str, elem_name: &str) -> &'a mut Self {
        //checking if 'elem_name' already exists
        if !self.does_it_exist(elem_name) {
            println!("Error at 'connect_from_to':");
            println!(" Object with name '{}' does not exist", elem_name);
            std::process::exit(1);
        }

        //finding `connector` in self.connector
        match self.connector.iter_mut().find(|c| c.name() == connector) {
            Some(c) => {
                match c.connect_to(elem_name) {
                    Ok(()) => {},
                    Err(err) => {
                        println!("Error at 'connect_from_to':");
                        println!(" {}", err);
                        std::process::exit(1);
                    }
                }
            },
            None => {
                println!("Error at 'connect_from_to':");
                println!(" Object with name '{}' does not exist", connector);
                std::process::exit(1);
            },
        };
        self
    }

    fn does_it_exist(&self, obj_name: &str) -> bool {
        self.objs_name.iter().any(|x| *x == obj_name)
    }

    fn _rs_in_path(elem: &str, path: &str) -> bool {
        let elem_file = format!("{}.rs", elem);
        let path_buf = std::path::Path::new(path).join(&elem_file);
        path_buf.exists()
    }
}
