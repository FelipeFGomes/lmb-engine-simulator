//! Solver for system of ordinary differential equations (ODE) 

use ndarray::*;

/// Integrates a system of ODEs over a single time step using 4th order Runge-Kutta
/// 
/// `x` is the initial condition store in a 1D-ndarray form, `c` is a Vec<f64> with constant values of the problem 
/// 
/// # Examples 
/// 
/// Let's solve the famous Lorentz system of equations: 
/// ```
/// let lorentz_eqs = |t: &f, x: &Array1<f64>, c: &Vec<f64>| -> Array1<f64> {
///     array![ c[0] * (x[1] - x[0]),
///             x[0] * (c[1] - x[2]) - x[1],
///             x[0] * x[1] - c[2] * x[2] ]
///     };
/// let ini_state = array![0.1, 0.1, 0.1];
/// let consts = vec![10.0, 28.0, 8.0/3.0];
/// let solution = rk4_step( lorentz_eqs, &ini_state, &consts, &0.0, 1e-4 );
/// ```
pub fn rk4_step<F>(f: F, x: &Array1<f64>, c: &Vec<f64>, t: &f64, step: f64) -> Array1<f64> 
    where
    F: Fn(&f64, &Array1<f64>, &Vec<f64>) -> Array1<f64>,
{
    let tmp = step/2.0;
    let tmp_2 = tmp+t;
    let k1 = f( &t, x, c );
    let k2 = f( &tmp_2, &(x + &(&k1*tmp)), c );
    let k3 = f( &tmp_2, &(x + &(&k2*tmp)), c );
    let k4 = f( &t, &(x + &(&k3*step)), c );
    let f_out = x + &((step/6.0)*( k1 + 2.0*k2 + 2.0*k3 + k4 ));
    f_out
}

/// Integrates a system of ODEs over a single time step using explicit Euler method
pub fn euler_step<F>(f: F, x: &Array1<f64>, c: &Vec<f64>, t: &f64, step: f64) -> Array1<f64> 
    where
    F: Fn(&f64, &Array1<f64>, &Vec<f64>) -> Array1<f64>,
{
    let k1 = f( &t, x, c );
    let f_out = x + &(step*k1) ;
    f_out
}