use na::Vector3;
use std::f64::consts;
use structs::geo_ellipsoid;
use structs::utm_grid;

/// Converts 3-d ENU coordinates to 3-d NED coordinates
/// 
/// # Arguments
/// 
/// * `enu_vec` - Vector3 reference to the ENU vector (x, y, z)
/// 
/// # Return Value
/// 
/// * `nalgebra::Vector3<f64>` - x, y, z
/// 
/// # Formula
/// 
/// * x = y
/// * y = x
/// * z = -z
pub fn enu2ned(enu_vec: &Vector3<f64>) -> Vector3<f64> {
    let mut ret_vec: Vector3<f64> = Vector3::new(0.0, 0.0, 0.0);
    ret_vec.x = enu_vec.y;
    ret_vec.y = enu_vec.x;
    ret_vec.z = -enu_vec.z;
    ret_vec
}

/// Converts 3-d NED coordinates to 3-d ENU coordinates
/// 
/// # Arguments
/// 
/// * `ned_vec` - Vector3 reference to the NED vector (x, y, z)
/// 
/// # Return Value
/// 
/// * `nalgebra::Vector3<f64>` - x, y, z
/// 
/// # Formula
/// 
/// * x = y
/// * y = x
/// * z = -z
pub fn ned2enu(ned_vec: &Vector3<f64>) -> Vector3<f64> {
    let mut ret_vec: Vector3<f64> = Vector3::new(0.0, 0.0, 0.0);
    ret_vec.x = ned_vec.y;
    ret_vec.y = ned_vec.x;
    ret_vec.z = -ned_vec.z;
    ret_vec
}


/// Converts 3-d LLA coordinates to 3-d ECEF coordinates
/// 
/// # Arguments
/// 
/// * `lla_vec` - Vector3 reference to the LLA vector (latitude, longitude, altitude) (radians, radians, meters)
/// 
/// # Return Value
/// 
/// * `nalgebra::Vector3<f64>` - x, y, z
/// 
/// # Formula
/// 
/// * x = (N + h) * cos(lat) * cos(lon)
/// * y = (N + h) * cos(lat) * sin(lon)
/// * z = (( b^2 / a^2 ) * N + h) * sin(lat)
pub fn lla2ecef(lla_vec: &Vector3<f64>, ellipsoid: &geo_ellipsoid::geo_ellipsoid) -> Vector3<f64> {
	let mut ret_vec: Vector3<f64> = Vector3::new(0.0, 0.0, 0.0);
	let n = ellipsoid.get_semi_major_axis() / (1.0 - ellipsoid.get_first_ecc().powi(2) * lla_vec.x.sin().powi(2)).sqrt();
	ret_vec.x = (n + lla_vec.z) * lla_vec.x.cos() * lla_vec.y.cos();
	ret_vec.y = (n + lla_vec.z) * lla_vec.x.cos() * lla_vec.y.sin();
	ret_vec.z = ((ellipsoid.get_semi_minor_axis().powi(2) / ellipsoid.get_semi_major_axis().powi(2)) * n + lla_vec.z) * lla_vec.x.sin();
	ret_vec
}


/// Converts 3-d ECEF coordinates to 3-d LLA coordinates
/// 
/// # Arguments
/// 
/// * `ecef_vec` - Vector3 reference to the ECEF vector (x, y, z)
/// 
/// # Return Value
/// 
/// * `nalgebra::Vector3<f64>` - lat, long, alt (radians, radians, meters)
/// 
/// # Formula
/// 
/// * x = arctan((z + e'^2 * b * sin^3 (theta)) / (p - e^2 * a * cos^3 (theta)))
/// * y = arctan(y / x)
/// * z = (p  / cos(lat)) - N
pub fn ecef2lla(ecef_vec: &Vector3<f64>, ellipsoid: &geo_ellipsoid::geo_ellipsoid) -> Vector3<f64> {
    let mut ret_vec: Vector3<f64> = Vector3::new(0.0, 0.0, 0.0);
    let p = (ecef_vec.x.powi(2) + ecef_vec.y.powi(2)).sqrt();
    let theta = (ecef_vec.z * ellipsoid.get_semi_major_axis()).atan2(p * ellipsoid.get_semi_minor_axis());
    let x_top = ecef_vec.z + ellipsoid.get_second_ecc().powi(2) * ellipsoid.get_semi_minor_axis() * theta.sin().powi(3);
    let x_bot = p - ellipsoid.get_first_ecc().powi(2) * ellipsoid.get_semi_major_axis() * theta.cos().powi(3);
    ret_vec.x = x_top.atan2(x_bot);
    ret_vec.y = ecef_vec.y.atan2(ecef_vec.x);
    let n = ellipsoid.get_semi_major_axis() / (1.0 - ellipsoid.get_first_ecc().powi(2) * (ret_vec.x.sin() * ret_vec.x.sin())).sqrt();
    ret_vec.z = (p / ret_vec.x.cos()) - n;
    ret_vec
}


//TODO: Comment (include Karney method)
//TODO: Clean this up a lot. Check the math - currently outputs the correct data -> minus the Norway/Svalbard case
//TODO: Reference: http://www.movable-type.co.uk/scripts/latlong-utm-mgrs.html (c) Chris Veness 2014-2017 MIT Licence
pub fn lla2utm(lla_vec: &Vector3<f64>, ellipsoid: &geo_ellipsoid::geo_ellipsoid) -> utm_grid::utm_grid {
    let mut ret_utm = utm_grid::utm_grid::new();
    let zone = ((lla_vec.y.to_degrees() + 180.0) / 6.0).floor() + 1.0; //Longitudinal zone
    ret_utm.set_zone(zone as u32);
    let lon_cent_mer = ((zone - 1.0) * 6.0 - 180.0 + 3.0).to_radians();

    //TODO: Handle Norway/Svalbard exceptions
    
    //Determine easting/northing
    let phi = lla_vec.x;                    //latitude +- from equator
    let lambda = lla_vec.y - lon_cent_mer;  //longitude +- from central merdian
    let a = ellipsoid.get_semi_major_axis();
    let f = ellipsoid.get_flattening();
    let k0 = utm_grid::SCALE_FACTOR_CENTERAL_MERIDIAN;
    let e = (f * (2.0 - f)).sqrt();
    let n = f / (2.0 - f);
    let n2 = n * n;
    let n3 = n2 * n;
    let n4 = n3 * n;
    let n5 = n4 * n;
    let n6 = n5 * n;
    let cLambda = lambda.cos();
    let sLambda = lambda.sin();
    let tLambda = lambda.tan();
    let tau = phi.tan();
    let sigma = (e * ((e * tau) / (1.0 + tau.powi(2)).sqrt()).atanh()).sinh();
    let taup = tau * (1.0 + sigma.powi(2)).sqrt() - sigma * (1.0 + tau.powi(2)).sqrt();
    let xip = taup.atan2(cLambda);
    let etap = (sLambda / (taup.powi(2) + cLambda.powi(2)).sqrt()).asinh();
    let A = a / (1.0 + n) * (1.0 + (1.0/4.0) * n2 + (1.0 / 64.0) * n4 + (1.0/256.0) * n6);
    let alpha = 
        [0.0,
        1.0/2.0*n - 2.0/3.0*n2 + 5.0/16.0*n3 + 41.0/180.0*n4 - 127.0/288.0*n5 + 7891.0/37800.0*n6,
        13.0/48.0*n2 - 3.0/5.0*n3 + 557.0/1440.0*n4 + 281.0/630.0*n5 - 1983433.0/1935360.0*n6,
        61.0/240.0*n3 - 103.0/140.0*n4 + 15061.0/26880.0*n5 + 167603.0/181440.0*n6,
        4956.01/161280.0*n4 - 179.0/168.0*n5 + 6601661.0/7257600.0*n6,
        4729.0/80640.0*n5 - 3418889.0/1995840.0*n6,
        212378941.0/319334400.0*n6];
    let mut xi = xip;
    for j in 1..7 {
        xi += alpha[j] * (2.0 * (j as f64) * xip).sin() * (2.0 * (j as f64) * etap).cosh();
    }
    let mut eta = etap;
    for j in 1..7 {
        eta += alpha[j] * (2.0 * (j as f64) * xip).cos() * (2.0 * (j as f64) * etap).sinh();
    }
    let mut easting = k0 * A * eta;
    let mut northing = k0 * A * xi;

    let mut pp = 1.0;
    for j in 1..7 {
        pp += (2.0 * (j as f64) * alpha[j]) * (2.0 * (j as f64) * xip).cos() * (2.0 * (j as f64) * etap).cosh();
    }
    let mut qp = 0.0;
    for j in 1..7 {
        qp += (2.0 * (j as f64) * alpha[j]) * (2.0 * (j as f64) * xip).sin() * (2.0 * (j as f64) * etap).sinh();
    }

    let yp = (taup / (1.0 + taup.powi(2)).sqrt() * tLambda).atan();
    let ypp = qp.atan2(pp);
    let mut y = yp + ypp;

    let sphi = phi.sin();
    let kp = (1.0 - e.powi(2) * sphi.powi(2)).sqrt() * (1.0 + tau.powi(2)).sqrt() / (taup.powi(2) + cLambda.powi(2)).sqrt();
    let kpp = A / a * (pp.powi(2) + qp.powi(2)).sqrt();
    let mut k = k0 * kp * kpp;
    easting += utm_grid::FALSE_EASTING;
    if northing < 0.0 {
        northing += utm_grid::FALSE_NORTHING;
    }

    if lla_vec.x >= 0.0 {
        ret_utm.set_hem(utm_grid::hemisphere::NORTH);
    }else{
        ret_utm.set_hem(utm_grid::hemisphere::SOUTH);
    }

    easting = (easting * 1000000.0).round() / 1000000.0;    //nm precision
    northing = (northing * 1000000.0).round() / 1000000.0;  //nm precision
    y = (y * 1000000000.0).round() / 1000000000.0;
    k = (k * 1000000000000.0).round() / 1000000000000.0;

    ret_utm.set_easting(easting);
    ret_utm.set_northing(northing);
    ret_utm.set_convergence(y);
    ret_utm.set_scale(k);

    ret_utm
}

//Unit tests
#[cfg(test)]
mod tests {
	use super::*;
    use float_cmp::ApproxEqUlps;
    use float_cmp::ApproxEqRatio;
    #[test]
    fn test_enu2ned() {
        let enu_vec: Vector3<f64> = Vector3::new(3.0, 4.0, 5.0);
        let ned_vec = enu2ned(&enu_vec);

        let test_x = 4.0;
        let test_y = 3.0;
        let test_z = -5.0;
        assert!(ned_vec.x.approx_eq_ulps(&test_x, 2));
        assert!(ned_vec.y.approx_eq_ulps(&test_y, 2));
        assert!(ned_vec.z.approx_eq_ulps(&test_z, 2));
    }
    #[test]
    fn test_ned2enu() {
        let ned_vec: Vector3<f64> = Vector3::new(3.0, 4.0, 5.0);
        let enu_vec = ned2enu(&ned_vec);

        let test_x = 4.0;
        let test_y = 3.0;
        let test_z = -5.0;
        assert!(enu_vec.x.approx_eq_ulps(&test_x, 2));
        assert!(enu_vec.y.approx_eq_ulps(&test_y, 2));
        assert!(enu_vec.z.approx_eq_ulps(&test_z, 2));
    }
    #[test]
    fn test_lla2ecef() {
    	let ellipsoid = geo_ellipsoid::geo_ellipsoid::new(geo_ellipsoid::WGS84_SEMI_MAJOR_AXIS_METERS,
    										geo_ellipsoid::WGS84_FLATTENING);
        let lla_vec: Vector3<f64> = Vector3::new(0.8527087756759584, 0.04105401863784606, 1000.000000000);
        let ecef_vec = lla2ecef(&lla_vec, &ellipsoid);

        let test_x = 4201570.9492264455;
        let test_y = 172588.3449531975;
        let test_z = 4780835.4317144295;
        assert!(ecef_vec.x.approx_eq_ratio(&test_x, 0.0000000001));
        assert!(ecef_vec.y.approx_eq_ratio(&test_y, 0.0000000001));
        assert!(ecef_vec.z.approx_eq_ratio(&test_z, 0.0000000001));
    }
    #[test]
    fn test_ecef2lla() {
        let ellipsoid = geo_ellipsoid::geo_ellipsoid::new(geo_ellipsoid::WGS84_SEMI_MAJOR_AXIS_METERS,
                                            geo_ellipsoid::WGS84_FLATTENING);
        let ecef_vec: Vector3<f64> = Vector3::new(4201570.9492264455, 172588.3449531975, 4780835.4317144295);
        let lla_vec = ecef2lla(&ecef_vec, &ellipsoid);

        let test_x = 0.8527087756759584;
        let test_y = 0.04105401863784606;
        let test_z = 1000.000000000;
        assert!(lla_vec.x.approx_eq_ratio(&test_x, 0.0000000001));
        assert!(lla_vec.y.approx_eq_ratio(&test_y, 0.0000000001));
        assert!(lla_vec.z.approx_eq_ratio(&test_z, 0.0000000001));
    }
    #[test]
    fn test_lla2utm() {
        let ellipsoid = geo_ellipsoid::geo_ellipsoid::new(geo_ellipsoid::WGS84_SEMI_MAJOR_AXIS_METERS,
                                            geo_ellipsoid::WGS84_FLATTENING);
        let lla_vec: Vector3<f64> = Vector3::new(0.8527087756759584, 0.04105401863784606, 1000.000000000);
        let utm = lla2utm(&lla_vec, &ellipsoid);
        //TODO: Write asserts
    }
}
