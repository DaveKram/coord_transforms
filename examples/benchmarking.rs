extern crate coord_transforms;
extern crate nalgebra as na;
extern crate chrono;
extern crate rayon;
use coord_transforms::*;
use na::Vector3;
use na::Vector2;
use chrono::prelude::*;
use rayon::prelude::*;

fn main() {
	let num_iters = 1000000;
	test_lla2ecef(num_iters);
	test_cartesian2spherical(num_iters);
	test_ll2utm(num_iters);
}

#[allow(unused_variables)]
fn test_cartesian2spherical(num_iters: u32) {
	//Test normal loop
    let mut d1: DateTime<Local> = Local::now();
	for _ in 0..num_iters {
		let cart_vec: Vector3<f64> = Vector3::new(3.0, 4.0, 5.0);
		d3::cartesian2spherical(&cart_vec);
	}
	let mut d2: DateTime<Local> = Local::now();
	let mut dur = d2.signed_duration_since(d1);
	println!("(cartesian2spherical): Basic loop - num_iters: {} took {:?} microseconds", num_iters, dur.num_microseconds().unwrap());

	//Test rayon loop
	d1 = Local::now();
	(0..num_iters).into_par_iter().map(|i| { 
		let cart_vec: Vector3<f64> = Vector3::new(3.0, 4.0, 5.0);
		d3::cartesian2spherical(&cart_vec) 
	}).count();
	d2 = Local::now();
	dur = d2.signed_duration_since(d1);
	println!("(cartesian2spherical): Rayon loop - num_iters: {} took {:?} microseconds", num_iters, dur.num_microseconds().unwrap());
}

#[allow(unused_variables)]
fn test_lla2ecef(num_iters: u32) {
    //Test normal loop
	let mut d1: DateTime<Local> = Local::now();
	for _ in 0..num_iters {
		let ellipsoid = structs::geo_ellipsoid::geo_ellipsoid::new(structs::geo_ellipsoid::WGS84_SEMI_MAJOR_AXIS_METERS,
    										structs::geo_ellipsoid::WGS84_FLATTENING);
    	let lla_vec: Vector3<f64> = Vector3::new(0.8527087756759584, 0.04105401863784606, 1000.000000000);
		geo::lla2ecef(&lla_vec, &ellipsoid);
	}
	let mut d2: DateTime<Local> = Local::now();
	let mut dur = d2.signed_duration_since(d1);
	println!("(lla2ecef): Basic loop - num_iters: {} took {:?} microseconds", num_iters, dur.num_microseconds().unwrap());

	//Test rayon loop
	d1 = Local::now();
	(0..num_iters).into_par_iter().map(|i| { 
		let ellipsoid = structs::geo_ellipsoid::geo_ellipsoid::new(structs::geo_ellipsoid::WGS84_SEMI_MAJOR_AXIS_METERS,
    										structs::geo_ellipsoid::WGS84_FLATTENING);
    	let lla_vec: Vector3<f64> = Vector3::new(0.8527087756759584, 0.04105401863784606, 1000.000000000);
		geo::lla2ecef(&lla_vec, &ellipsoid) 
	}).count();
	d2 = Local::now();
	dur = d2.signed_duration_since(d1);
	println!("(lla2ecef): Rayon loop - num_iters: {} took {:?} microseconds", num_iters, dur.num_microseconds().unwrap());
}

#[allow(unused_variables)]
fn test_ll2utm(num_iters: u32) {
    //Test normal loop
	let mut d1: DateTime<Local> = Local::now();
	for _ in 0..num_iters {
		let ellipsoid = structs::geo_ellipsoid::geo_ellipsoid::new(structs::geo_ellipsoid::WGS84_SEMI_MAJOR_AXIS_METERS,
    										structs::geo_ellipsoid::WGS84_FLATTENING);
    	let ll_vec: Vector2<f64> = Vector2::new(1.3804121468, 0.20555336013);
		geo::ll2utm(&ll_vec, &ellipsoid);
	}
	let mut d2: DateTime<Local> = Local::now();
	let mut dur = d2.signed_duration_since(d1);
	println!("(ll2utm): Basic loop - num_iters: {} took {:?} microseconds", num_iters, dur.num_microseconds().unwrap());

	//Test rayon loop
	d1 = Local::now();
	(0..num_iters).into_par_iter().map(|i| { 
		let ellipsoid = structs::geo_ellipsoid::geo_ellipsoid::new(structs::geo_ellipsoid::WGS84_SEMI_MAJOR_AXIS_METERS,
    										structs::geo_ellipsoid::WGS84_FLATTENING);
    	let ll_vec: Vector2<f64> = Vector2::new(1.3804121468, 0.20555336013);
		geo::ll2utm(&ll_vec, &ellipsoid);
	}).count();
	d2 = Local::now();
	dur = d2.signed_duration_since(d1);
	println!("(ll2utm): Rayon loop - num_iters: {} took {:?} microseconds", num_iters, dur.num_microseconds().unwrap());
}