#![allow(non_snake_case)]

use super::system::System;
use crate::zero_dim;
use crate::connector;
// use crate::one_dim;
use crate::engine::engine::Engine;
use crate::reaction::gas::Gas;
use crate::{ObjectInfo, ObjectType, StoreData};

// Core Traits
use crate::core::traits::Connector;
// use crate::core::traits::OneDim;
use crate::core::traits::SaveData;
use crate::core::traits::ZeroDim;

// Super Traits
use crate::core::traits::Conn;
use crate::core::traits::OneD;
use crate::core::traits::ZeroD;

pub struct SystemBuilder {
    objs_info: Vec<ObjectInfo>,
    engine: Option<Engine>,
    zero_dim: Vec<Box<dyn ZeroD>>,
    one_dim: Vec<Box<dyn OneD>>,
    connector: Vec<Box<dyn Conn>>,
}

impl SystemBuilder {

    /// Creates a `SystemBuilder`. This object is used to construct the desired system to simulate. 
    /// The construction is made by the object methods exclusively. 
    /// Once the building is finished, the system can be built using `build_system()` method.
    pub fn new() -> SystemBuilder {
        SystemBuilder {
            objs_info: Vec::new(),
            engine: None,
            zero_dim: Vec::new(),
            one_dim: Vec::new(),
            connector: Vec::new(),
        }
    }

    /// Build a `System`. `SystemBuilder` objects is consumed in the process.
    pub fn build_system(self) -> System {
        match System::new(
            self.objs_info,
            self.engine,
            self.zero_dim,
            self.one_dim,
            self.connector,
        ) {
            Ok(s) => s,
            Err(err) => {
                println!("Error at 'SystemBuilder::build_system()':\n {}", err);
                std::process::exit(1)
            }
        }
    }

    /// Add a `Engine` and its components from a `.json` file and `Gas` object. The mandatory components are `cylinders` and `valves`.
    /// The variable `gas` is clone into the `zero_dim::Cylinder` objects.
    pub fn add_engine<'a>(&'a mut self, file_name: &str, gas: &Gas) -> &'a mut Self {
        let engine = match Engine::new(file_name, gas) {
            Ok(eng) => eng,
            Err(err) => {
                println!("Error at 'add_engine':\n {}", err);
                std::process::exit(1)
            }
        };

        for (i, cyl) in engine.cylinders().iter().enumerate() {
            if self.does_it_exist(cyl.name()) {
                println!("Error at `add_engine`:");
                println!(" Object with the same name already exists: {}", cyl.name());
                std::process::exit(1);
            }
            self.objs_info.push(ObjectInfo::new(
                cyl.name().to_string(),
                ObjectType::Cylinder,
                i,
                StoreData::new(&cyl.get_headers(), cyl.num_storable_variables()),
            ));
        }

        for val in engine.valves() {
            if self.does_it_exist(val.name()) {
                println!("Error at `add_engine`:");
                println!(" Object with the same name already exists: {}", val.name());
                std::process::exit(1);
            }
            let i = self.connector.len();
            self.objs_info.push(ObjectInfo::new(
                val.name().to_string(),
                ObjectType::Connector,
                i,
                StoreData::new(&val.get_headers(), val.num_storable_variables()),
            ));
            self.connector.push(Box::new(val.clone()));
        }

        // Pushing ´engine´
        self.engine = Some(engine);
        self
    }

    /// Add a `zero_dim::Enrivonment`. It has constant temperature and pressure and infinite mass
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

        // adding to list of objects
        let i = self.zero_dim.len();
        self.objs_info.push(ObjectInfo::new(
                elem_name.to_string(),
                ObjectType::ZeroDim,
                i,
                StoreData::new(&env.get_headers(), env.num_storable_variables()),
            ));

        self.zero_dim.push(Box::new(env));
        self
    }

    /// Add a `zero_dim::Reservoir`. It corresponds to a plenum or a constant volume chamber. The input `volume`
    /// must be in cubic centimeters [cm³].
    pub fn add_reservoir<'a>(
        &'a mut self,
        elem_name: &str,
        volume: f64,
        gas: &Gas,
    ) -> &'a mut Self {
        // checking if 'elem_name' already exists
        if self.does_it_exist(elem_name) {
            println!("Error at 'add_reservoir':");
            println!("Object with the same name already exists");
            std::process::exit(1)
        }

        // pushing reservoir
        let res = match zero_dim::reservoir::Reservoir::new(elem_name.to_string(), gas, volume*1e-6) {
            Ok(v) => v,
            Err(err) => {
                println!("Error at 'add_reservoir':\n {}", err);
                std::process::exit(1)
            }
        };

        // adding to list of objects
        let i = self.zero_dim.len();
        self.objs_info.push(ObjectInfo::new(
            elem_name.to_string(),
            ObjectType::ZeroDim,
            i,
            StoreData::new(&res.get_headers(), res.num_storable_variables()),
        ));

        self.zero_dim.push(Box::new(res));

        self
    }

    /// Add a `connector::Orifice` connector. It connects two `ZeroDim` through a hole of diameter `diam` in mm.
    /// The discharge coefficient must be between 0 and 1.
    pub fn add_orifice<'a>(&'a mut self, elem_name:&str, diam: f64, discharge_coeff: f64, conn: Vec<&str>) -> &'a mut Self {
        // checking if 'elem_name' already exists
        if self.does_it_exist(elem_name) {
            println!("Error at 'add_connection_between_0D':");
            println!("Object with the same name already exists: `{}`", elem_name);
            std::process::exit(1)
        }

        // pushing connector
        let mut connecting: Vec<String> = Vec::new();
        for c in conn {connecting.push(c.to_string());}
        let zd_conn = match connector::orifice::Orifice::new(elem_name, diam*1e-3, discharge_coeff, connecting) {
            Ok(v) => v,
            Err(err) => {
                println!("Error at 'add_connection_between_0D':\n {}", err);
                std::process::exit(1)
            }
        };

        // adding to list of objects
        let i = self.connector.len();
        self.objs_info.push(ObjectInfo::new(
            zd_conn.name().to_string(),
            ObjectType::Connector,
            i,
            StoreData::new(&zd_conn.get_headers(), zd_conn.num_storable_variables()),
        ));
        self.connector.push(Box::new(zd_conn));
        self
    }

    /// Connect a `connector` object to an element object. The inputs must be the name of the connector and element as `&str`
    pub fn connect_from_to<'a>(&'a mut self, connector: &str, elem_name: &str) -> &'a mut Self {
        //checking if 'elem_name' already exists
        if !self.does_it_exist(elem_name) {
            println!("Error at 'connect_from_to':");
            println!(" Object with name '{}' does not exist", elem_name);
            std::process::exit(1);
        }

        let obj_info = self
            .objs_info
            .iter()
            .find(|info| info.name == connector)
            .unwrap();

        let i = obj_info.index;
        match obj_info.obj_type {
            ObjectType::Connector => {
                if let Err(err) = self.connector[i].connect_to(elem_name) {
                    println!("Error at 'SystemBuilder::connect_from_to()':");
                    println!(" {}", err);
                    std::process::exit(1);
                }
            }
            _ => {
                println!("Error at 'SystemBuilder::connect_from_to()':");
                println!("Unknown object type \"{:?}\" Cannot be used", obj_info.obj_type);
                std::process::exit(1);
            }
        }
        self
    }

    fn does_it_exist(&self, obj_name: &str) -> bool {
        self.objs_info.iter().any(|info| info.name == obj_name)
    }

    fn _rs_in_path(elem: &str, path: &str) -> bool {
        let elem_file = format!("{}.rs", elem);
        let path_buf = std::path::Path::new(path).join(&elem_file);
        path_buf.exists()
    }
}
