/*
    Note: This example "cascades" from top to bottom
    It attempts to use the previous converted version of the LLA vec to show
    that there are no errors converting between coordinate systems
*/

extern crate coord_transforms;
use coord_transforms::prelude::*;

fn main() {
    let lat_deg: f64 = 57.77348022;
    let lon_deg: f64 = 157.37338121;
    let alt_m: f64 = 1000.0;

    let lla_vec = Vector3::new(lat_deg.to_radians(), lon_deg.to_radians(), alt_m);
    println!("-----LLA (Base)-----");
    println!("Lat: {}", lla_vec.x.to_degrees());
    println!("Lon: {}", lla_vec.y.to_degrees());
    println!("Alt: {}", lla_vec.z);

    //Define ellipsoid for the Earth
    let ellipsoid = geo_ellipsoid::geo_ellipsoid::new(geo_ellipsoid::WGS84_SEMI_MAJOR_AXIS_METERS,
    										geo_ellipsoid::WGS84_FLATTENING);

    //Convert to Earth-Centered Earth-Fixed (ECEF)
    let ecef_vec = geo::lla2ecef(&lla_vec, &ellipsoid);
    println!("-----ECEF-----");
    println!("x: {}", ecef_vec.x);
    println!("y: {}", ecef_vec.y);
    println!("z: {}", ecef_vec.z);

    //Convert it back to LLA
    let lla_again_vec = geo::ecef2lla(&ecef_vec, &ellipsoid);
    println!("-----LLA Again-----");
    println!("Lat: {}", lla_again_vec.x.to_degrees());
    println!("Lon: {}", lla_again_vec.y.to_degrees());
    println!("Alt: {}", lla_again_vec.z);

    //Convert it to UTM
    let ll_vec = Vector2::new(lla_again_vec.x, lla_again_vec.y);
    let utm = geo::ll2utm(&ll_vec, &ellipsoid);
    println!("-----UTM-----");
    println!("Zone: {}", utm.get_zone());
    println!("Hemisphere: {:?}", utm.get_hem());
    println!("Easting: {}", utm.get_easting());
    println!("Northing: {}", utm.get_northing());
    println!("Convergence: {}", utm.get_convergence());
    println!("Scale: {}", utm.get_scale());

    //Convert UTM to LL
    let llfromutm_vec = geo::utm2ll(&utm, &ellipsoid);
    println!("-----LL from UTM-----");
    println!("Lat: {}", llfromutm_vec.x.to_degrees());
    println!("Lon: {}", llfromutm_vec.y.to_degrees());


    //Convert it to NED and back
    let lla_orig_vec = Vector3::new((lat_deg + 1.0).to_radians(), (lon_deg + 1.0).to_radians(), alt_m);
    let ned_vec = geo::lla2ned(&lla_orig_vec, &lla_again_vec, &ellipsoid);
    println!("-----NED-----");
    println!("LLA Origin: {} - {} - {}", lla_orig_vec.x.to_degrees(), lla_orig_vec.y.to_degrees(), lla_orig_vec.z);
    println!("N: {}", ned_vec.x);
    println!("E: {}", ned_vec.y);
    println!("D: {}", ned_vec.z);

    //Simple NED-ENU
    let simple_enu = geo::ned2enu(&ned_vec);
    println!("-----ENU-----");
    println!("LLA Origin: {} - {} - {}", lla_orig_vec.x.to_degrees(), lla_orig_vec.y.to_degrees(), lla_orig_vec.z);
    println!("E: {}", simple_enu.x);
    println!("N: {}", simple_enu.y);
    println!("U: {}", simple_enu.z);

    //Convert NED back into LLA
    let lla_revert_ned = geo::ned2lla(&lla_orig_vec, &ned_vec, &ellipsoid);
    println!("-----LLA (from NED)-----");
    println!("LLA Origin: {} - {} - {}", lla_orig_vec.x.to_degrees(), lla_orig_vec.y.to_degrees(), lla_orig_vec.z);
    println!("Lat: {}", lla_revert_ned.x.to_degrees());
    println!("Lon: {}", lla_revert_ned.y.to_degrees());
    println!("Alt: {}", lla_revert_ned.z);

    //Convert ENU back into LLA
    let lla_revert_enu = geo::enu2lla(&lla_orig_vec, &simple_enu, &ellipsoid);
    println!("-----LLA (from ENU)-----");
    println!("LLA Origin: {} - {} - {}", lla_orig_vec.x.to_degrees(), lla_orig_vec.y.to_degrees(), lla_orig_vec.z);
    println!("Lat: {}", lla_revert_enu.x.to_degrees());
    println!("Lon: {}", lla_revert_enu.y.to_degrees());
    println!("Alt: {}", lla_revert_enu.z);

}