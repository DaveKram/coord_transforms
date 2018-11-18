use na::Vector3;
use na::Vector2;
use na::Matrix3;
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
/// * `ellipsoid` - geo_ellipsoid reference to the ellipsoid
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
/// * `ellipsoid` - geo_ellipsoid reference to the ellipsoid
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

/// Converts 2-d LL coordinates to UTM Grid (using Karney method) - accruacy to within a few nanometers within 3900km of the central meridian
/// 
/// # Arguments
/// 
/// * `ll_vec` - Vector2 reference to the LL vector (latitude, longitude) (radians, radians)
/// * `ellipsoid` - geo_ellipsoid reference to the ellipsoid
/// 
/// # Return Value
/// 
/// * `utm_grid::utm_grid` - UTM Grid the latitude and longitude belong in
/// 
/// # Notes
/// 
/// * Based on code from here: http://www.movable-type.co.uk/scripts/latlong-utm-mgrs.html 
/// * Based on white paper from here: https://arxiv.org/abs/1002.1417
/// * (c) Chris Veness 2014-2017 MIT Licence
pub fn ll2utm(ll_vec: &Vector2<f64>, ellipsoid: &geo_ellipsoid::geo_ellipsoid) -> utm_grid::utm_grid {
    let mut ret_utm = utm_grid::utm_grid::new(0, utm_grid::hemisphere::NORTH, 0.0, 0.0, 0.0, 0.0);
    let mut zone = ((ll_vec.y.to_degrees() + 180.0) / 6.0).floor() + 1.0;
    let mut lon_cent_mer = ((zone - 1.0) * 6.0 - 180.0 + 3.0).to_radians();

    //Handle Norway/Svalbard exceptions
    let mgrs_lat_bands = ['C', 'D', 'E', 'F', 'G', 'H', 'J', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'X'];
    let lat_band = mgrs_lat_bands[(ll_vec.x.to_degrees() / 8.0 + 10.0).floor() as usize];
    let loncentmer_shift_rads = 0.10472;
    if zone == 31.0 && lat_band == 'V' && ll_vec.y >= 0.0523599 {
        zone += 1.0; lon_cent_mer += loncentmer_shift_rads;
    }
    if zone == 32.0 && lat_band == 'X' && ll_vec.y < 0.15708 {
        zone -= 1.0; lon_cent_mer -= loncentmer_shift_rads;
    }
    if zone == 32.0 && lat_band == 'X' && ll_vec.y >= 0.15708 {
        zone += 1.0; lon_cent_mer += loncentmer_shift_rads;
    }
    if zone == 34.0 && lat_band == 'X' && ll_vec.y < 0.366519 {
        zone -= 1.0; lon_cent_mer -= loncentmer_shift_rads;
    }
    if zone == 34.0 && lat_band == 'X' && ll_vec.y >= 0.366519 {
        zone += 1.0; lon_cent_mer += loncentmer_shift_rads;
    }
    if zone == 36.0 && lat_band == 'X' && ll_vec.y < 0.575959 {
        zone -= 1.0; lon_cent_mer -= loncentmer_shift_rads;
    }
    if zone == 36.0 && lat_band == 'X' && ll_vec.y >= 0.575959 {
        zone += 1.0; lon_cent_mer += loncentmer_shift_rads;
    }

    //Determine easting/northing
    let lambda = ll_vec.y - lon_cent_mer;
    let c_lambda = lambda.cos();
    let s_lambda = lambda.sin();
    let t_lambda = lambda.tan();
    let e = (ellipsoid.get_flattening() * (2.0 - ellipsoid.get_flattening())).sqrt();
    let n = ellipsoid.get_flattening() / (2.0 - ellipsoid.get_flattening());
    let n2 = n * n;
    let n3 = n * n * n;
    let n4 = n * n * n * n;
    let n5 = n * n * n * n * n;
    let n6 = n * n * n * n * n * n;
    let sigma = (e * ((e * ll_vec.x.tan()) / (1.0 + ll_vec.x.tan().powi(2)).sqrt()).atanh()).sinh();
    let taup = ll_vec.x.tan() * (1.0 + sigma.powi(2)).sqrt() - sigma * (1.0 + ll_vec.x.tan().powi(2)).sqrt();
    let xip = taup.atan2(c_lambda);
    let etap = (s_lambda / (taup.powi(2) + c_lambda.powi(2)).sqrt()).asinh();
    let a_maj = ellipsoid.get_semi_major_axis() / (1.0 + n) * (1.0 + (1.0 / 4.0) * n2 + (1.0 / 64.0) * n4 + (1.0 / 256.0) * n6);

    let alpha = 
        [0.0,
        1.0 / 2.0 * n - 2.0 / 3.0 * n2 + 5.0 / 16.0 * n3 + 41.0 / 180.0 * n4 - 127.0 / 288.0 * n5 + 7891.0 / 37800.0 * n6,
        13.0 / 48.0 * n2 - 3.0 / 5.0 * n3 + 557.0 / 1440.0 * n4 + 281.0 / 630.0 * n5 - 1983433.0 / 1935360.0 * n6,
        61.0 / 240.0 * n3 - 103.0 / 140.0 * n4 + 15061.0 / 26880.0 * n5 + 167603.0 / 181440.0 * n6,
        4956.01 / 161280.0 * n4 - 179.0 / 168.0 * n5 + 6601661.0 / 7257600.0 * n6,
        4729.0 / 80640.0 * n5 - 3418889.0 / 1995840.0 * n6,
        212378941.0 / 319334400.0 * n6];

    let mut xi = xip;
    for j in 1..7 {
        xi += alpha[j] * (2.0 * (j as f64) * xip).sin() * (2.0 * (j as f64) * etap).cosh();
    }
    let mut eta = etap;
    for j in 1..7 {
        eta += alpha[j] * (2.0 * (j as f64) * xip).cos() * (2.0 * (j as f64) * etap).sinh();
    }
    let mut easting = utm_grid::SCALE_FACTOR_CENTERAL_MERIDIAN * a_maj * eta;
    let mut northing = utm_grid::SCALE_FACTOR_CENTERAL_MERIDIAN * a_maj * xi;

    //Determine convergence
    let mut pp = 1.0;
    for j in 1..7 {
        pp += (2.0 * (j as f64) * alpha[j]) * (2.0 * (j as f64) * xip).cos() * (2.0 * (j as f64) * etap).cosh();
    }
    let mut qp = 0.0;
    for j in 1..7 {
        qp += (2.0 * (j as f64) * alpha[j]) * (2.0 * (j as f64) * xip).sin() * (2.0 * (j as f64) * etap).sinh();
    }

    let yp = (taup / (1.0 + taup.powi(2)).sqrt() * t_lambda).atan();
    let ypp = qp.atan2(pp);
    let mut y = yp + ypp;

    //Determine scale
    let sphi = ll_vec.x.sin();
    let kp = (1.0 - e.powi(2) * sphi.powi(2)).sqrt() * (1.0 + ll_vec.x.tan().powi(2)).sqrt() / (taup.powi(2) + c_lambda.powi(2)).sqrt();
    let kpp = a_maj / ellipsoid.get_semi_major_axis() * (pp.powi(2) + qp.powi(2)).sqrt();
    let mut k = utm_grid::SCALE_FACTOR_CENTERAL_MERIDIAN * kp * kpp;

    //Shift northing and easting to false origins
    easting += utm_grid::FALSE_EASTING;
    if northing < 0.0 {
        northing += utm_grid::FALSE_NORTHING;
    }

    //Round to correct precision - easting and northing are nanometer precision
    easting = (easting * 1000000.0).round() / 1000000.0;
    northing = (northing * 1000000.0).round() / 1000000.0;
    y = (y * 1000000000.0).round() / 1000000000.0;
    k = (k * 1000000000000.0).round() / 1000000000000.0;

    //Set values in structure
    ret_utm.set_zone(zone as u32);
    ret_utm.set_easting(easting);
    ret_utm.set_northing(northing);
    ret_utm.set_convergence(y);
    ret_utm.set_scale(k);
    if ll_vec.x >= 0.0 {
        ret_utm.set_hem(utm_grid::hemisphere::NORTH);
    }else{
        ret_utm.set_hem(utm_grid::hemisphere::SOUTH);
    }

    ret_utm
}

/// Converts UTM Grid coordinates to 2-d LL (using Karney method) - accruacy to within a few nanometers within 3900km of the central meridian
/// 
/// # Arguments
/// 
/// * `utm` - utm_grid reference to the UTM grid
/// * `ellipsoid` - geo_ellipsoid reference to the ellipsoid
/// 
/// # Return Value
/// 
/// * `nalgebra::Vector2<f64>` - Lat, Lon (radians, radians)
/// 
/// # Notes
/// 
/// * Based on code from here: http://www.movable-type.co.uk/scripts/latlong-utm-mgrs.html 
/// * Based on white paper from here: https://arxiv.org/abs/1002.1417
/// * (c) Chris Veness 2014-2017 MIT Licence
pub fn utm2ll(utm: &utm_grid::utm_grid, ellipsoid: &geo_ellipsoid::geo_ellipsoid) -> Vector2<f64> {
    let mut ret_val: Vector2<f64> = Vector2::new(0.0, 0.0);
    let z = utm.get_zone();
    let h = utm.get_hem();
    let mut x = utm.get_easting();
    let mut y = utm.get_northing();

    let a = ellipsoid.get_semi_major_axis();
    let f = ellipsoid.get_flattening();
    let e = ellipsoid.get_first_ecc();

    //Make x relative to central meridian
    x = x - utm_grid::FALSE_EASTING;

    //Make y relative to equator
    if h == utm_grid::hemisphere::SOUTH {
        y = y - utm_grid::FALSE_NORTHING;
    }

    let n = f / (2.0 - f);
    let n2 = n * n;
    let n3 = n * n2;
    let n4 = n * n3;
    let n5 = n * n4;
    let n6 = n * n5;

    let A = a / (1.0 + n) * (1.0 + 1.0/4.0 * n2 + 1.0 / 64.0 * n4 + 1.0 / 256.0 * n6);
    let eta = x / (utm_grid::SCALE_FACTOR_CENTERAL_MERIDIAN * A);
    let xi = y / (utm_grid::SCALE_FACTOR_CENTERAL_MERIDIAN * A);

    let beta = 
        [0.0,
            1.0 / 2.0 * n - 2.0 / 3.0 * n2 + 37.0 / 96.0 * n3 - 1.0 / 360.0 * n4 - 81.0 / 512.0 * n5 + 96199.0 / 604800.0 * n6,
            1.0 / 48.0 * n2 + 1.0 / 15.0 * n3 - 437.0 / 440.0 * n4 + 46.0 / 105.0 * n5 - 1118711.0 / 3870720.0 * n6,
            17.0 / 480.0 * n3 - 37.0 / 840.0 * n4 - 209.0 / 4480.0 * n5 + 5569.0 / 90720.0 * n6,
            4397.0 / 161280.0 * n4 - 11.0 / 504.0 * n5 - 830251.0 / 7257600.0 * n6,
            4583.0 / 161280.0 * n5 - 108847.0 / 3991680.0 * n6,
            20648693.0 / 638668800.0 * n6];

    let mut xip = xi;
    for j in 1..7 {
        xip -= beta[j] * (2.0 * (j as f64) * xi).sin() * (2.0 * (j as f64) * eta).cosh();
    }

    let mut etap = eta;
    for j in 1..7 {
        etap -= beta[j] * (2.0 * (j as f64) * xi).cos() * (2.0 * (j as f64) * eta).sinh();
    }

    let sinetap = etap.sinh();
    let sinxip = xip.sin();
    let cosxip = xip.cos();

    let taup = sinxip / (sinetap * sinetap + cosxip * cosxip).sqrt();
    let mut taui = taup;
    while {
        let sigmai = (e * (e * taui / (1.0 + taui * taui).sqrt()).atanh()).sinh();
        let tauip = taui * (1.0 + sigmai * sigmai).sqrt() - sigmai * (1.0 + taui * taui).sqrt();
        let delataui = (taup - tauip) / (1.0 + tauip * tauip).sqrt() * (1.0 + (1.0 - e * e) * taui * taui) / ((1.0 - e * e) * (1.0 + taui * taui).sqrt());
        taui += delataui;
        delataui.abs() > 1e-12
    } {}

    let phi = taui.atan();
    let mut lambda = sinetap.atan2(cosxip);

    let lambda0 = (((z - 1) * 6 - 180 + 3) as f64).to_radians();
    lambda += lambda0;

    ret_val.x = phi;
    ret_val.y = lambda;

    ret_val
}

/// Converts 3-d LLA origin coordinates plus 3-d LLA coordinates and an ellipsoid to 3-d local NED cartesian coordinates
/// 
/// # Arguments
/// 
/// * `lla_origin` - Vector3 reference to the LLA origin vector (lat, long, alt) (radians, radians, meters)
/// * `lla_vec` - Vector3 reference to the LLA vector (lat, long, alt) (radians, radians, meters)
/// * `ellipsoid` - geo_ellipsoid reference to the ellipsoid
/// 
/// # Return Value
/// 
/// * `nalgebra::Vector3<f64>` - x, y, z
/// 
pub fn lla2ned(lla_origin: &Vector3<f64>, lla_vec: &Vector3<f64>, ellipsoid: &geo_ellipsoid::geo_ellipsoid) -> Vector3<f64> {
    let orig = lla2ecef(lla_origin, ellipsoid);
    let actual = lla2ecef(lla_vec, ellipsoid);
    let ned_rot = Matrix3::new(-lla_origin.x.sin() * lla_origin.y.cos(), -lla_origin.x.sin()  * lla_origin.y.sin(), lla_origin.x.cos(),
                                -lla_origin.y.sin(), lla_origin.y.cos(), 0.0,
                                -lla_origin.x.cos() * lla_origin.y.cos(), -lla_origin.x.cos() * lla_origin.y.sin(), -lla_origin.x.sin());
    ned_rot * (actual - orig)
}

/// Converts 3-d LLA origin coordinates plus 3-d LLA coordinates and an ellipsoid to 3-d local ENU cartesian coordinates
/// 
/// # Arguments
/// 
/// * `lla_origin` - Vector3 reference to the LLA origin vector (lat, long, alt) (radians, radians, meters)
/// * `lla_vec` - Vector3 reference to the LLA vector (lat, long, alt) (radians, radians, meters)
/// * `ellipsoid` - geo_ellipsoid reference to the ellipsoid
/// 
/// # Return Value
/// 
/// * `nalgebra::Vector3<f64>` - x, y, z
/// 
pub fn lla2enu(lla_origin: &Vector3<f64>, lla_vec: &Vector3<f64>, ellipsoid: &geo_ellipsoid::geo_ellipsoid) -> Vector3<f64> {
    let orig = lla2ecef(lla_origin, ellipsoid);
    let actual = lla2ecef(lla_vec, ellipsoid);
    let enu_rot = Matrix3::new(-lla_origin.y.sin(), lla_origin.y.cos(), 0.0,
                                    -lla_origin.y.cos() * lla_origin.x.sin(), -lla_origin.y.sin() * lla_origin.x.sin(), lla_origin.x.cos(),
                                    lla_origin.y.cos() * lla_origin.x.cos(), lla_origin.y.sin() * lla_origin.x.cos(), lla_origin.x.sin());
    enu_rot * (actual - orig)
}

/// Converts 3-d local cartesian NED coordinates plus 3-d LLA origin coordinates and an ellipsoid to 3-d LLA coordinates
/// 
/// # Arguments
/// 
/// * `lla_origin` - Vector3 reference to the LLA origin vector (lat, long, alt) (radians, radians, meters)
/// * `ned_vec` - Vector3 reference to the local NED cartesian vector (x, y, z)
/// * `ellipsoid` - geo_ellipsoid reference to the ellipsoid
/// 
/// # Return Value
/// 
/// * `nalgebra::Vector3<f64>` - (lat, long, alt) (radians, radians, meters)
/// 
pub fn ned2lla(lla_origin: &Vector3<f64>, ned_vec: &Vector3<f64>, ellipsoid: &geo_ellipsoid::geo_ellipsoid) -> Vector3<f64> {
    let orig = lla2ecef(lla_origin, ellipsoid);
    let ned_rot = Matrix3::new(-lla_origin.x.sin() * lla_origin.y.cos(), -lla_origin.y.sin(), -lla_origin.x.cos() * lla_origin.y.cos(),
                                -lla_origin.x.sin() * lla_origin.y.sin(), lla_origin.y.cos(), -lla_origin.x.cos() * lla_origin.y.sin(),
                                lla_origin.x.cos(), 0.0, -lla_origin.x.sin());
    let actual = ned_rot * ned_vec;
    let total = actual + orig;
    let ret_vec = ecef2lla(&total, ellipsoid);
    ret_vec
}

/// Converts 3-d local cartesian ENU coordinates plus 3-d LLA origin coordinates and an ellipsoid to 3-d LLA coordinates
/// 
/// # Arguments
/// 
/// * `lla_origin` - Vector3 reference to the LLA origin vector (lat, long, alt) (radians, radians, meters)
/// * `enu_vec` - Vector3 reference to the local NED cartesian vector (x, y, z)
/// * `ellipsoid` - geo_ellipsoid reference to the ellipsoid
/// 
/// # Return Value
/// 
/// * `nalgebra::Vector3<f64>` - (lat, long, alt) (radians, radians, meters)
/// 
pub fn enu2lla(lla_origin: &Vector3<f64>, enu_vec: &Vector3<f64>, ellipsoid: &geo_ellipsoid::geo_ellipsoid) -> Vector3<f64> {
    let orig = lla2ecef(lla_origin, ellipsoid);
    let enu_rot = Matrix3::new(-lla_origin.y.sin(), -lla_origin.y.cos() * lla_origin.x.sin(), lla_origin.y.cos() * lla_origin.x.cos(), 
                                lla_origin.y.cos(), -lla_origin.y.sin() * lla_origin.x.sin(), lla_origin.y.sin() * lla_origin.x.cos(), 
                                0.0, lla_origin.x.cos(), lla_origin.x.sin());
    let actual = enu_rot * enu_vec;
    let total = actual + orig;
    let ret_vec = ecef2lla(&total, ellipsoid);
    ret_vec
}

//Unit tests
#[cfg(test)]
mod tests {
	use super::*;
    use float_cmp::ApproxEqRatio;
    #[test]
    fn test_enu2ned() {
        let enu_vec: Vector3<f64> = Vector3::new(3.0, 4.0, 5.0);
        let ned_vec = enu2ned(&enu_vec);

        let test_x = 4.0;
        let test_y = 3.0;
        let test_z = -5.0;
        assert!(ned_vec.x.approx_eq_ratio(&test_x, 0.00025));
        assert!(ned_vec.y.approx_eq_ratio(&test_y, 0.00025));
        assert!(ned_vec.z.approx_eq_ratio(&test_z, 0.00025));
    }
    #[test]
    fn test_ned2enu() {
        let ned_vec: Vector3<f64> = Vector3::new(3.0, 4.0, 5.0);
        let enu_vec = ned2enu(&ned_vec);

        let test_x = 4.0;
        let test_y = 3.0;
        let test_z = -5.0;
        assert!(enu_vec.x.approx_eq_ratio(&test_x, 0.00025));
        assert!(enu_vec.y.approx_eq_ratio(&test_y, 0.00025));
        assert!(enu_vec.z.approx_eq_ratio(&test_z, 0.00025));
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
        assert!(ecef_vec.x.approx_eq_ratio(&test_x, 0.00025));
        assert!(ecef_vec.y.approx_eq_ratio(&test_y, 0.00025));
        assert!(ecef_vec.z.approx_eq_ratio(&test_z, 0.00025));
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
        assert!(lla_vec.x.approx_eq_ratio(&test_x, 0.00025));
        assert!(lla_vec.y.approx_eq_ratio(&test_y, 0.00025));
        assert!(lla_vec.z.approx_eq_ratio(&test_z, 0.00025));
    }
    #[test]
    fn test_ll2utm() {
        let ellipsoid = geo_ellipsoid::geo_ellipsoid::new(geo_ellipsoid::WGS84_SEMI_MAJOR_AXIS_METERS,
                                            geo_ellipsoid::WGS84_FLATTENING);
        let ll_vec: Vector2<f64> = Vector2::new(1.3804121468, 0.20555336013);
        let utm = ll2utm(&ll_vec, &ellipsoid);
        
        let test_zone = 33;
        let test_hem = utm_grid::hemisphere::NORTH;
        let test_easting = 431952.612166;
        let test_northing = 8782098.22289;
        let test_convergence = -0.055232079;
        let test_scale = 0.999656581563;

        assert_eq!(utm.get_zone(), test_zone);
        assert_eq!(utm.get_hem(), test_hem);
        assert!(utm.get_easting().approx_eq_ratio(&test_easting, 0.00025));
        assert!(utm.get_northing().approx_eq_ratio(&test_northing, 0.00025));
        assert!(utm.get_convergence().approx_eq_ratio(&test_convergence, 0.00025));
        assert!(utm.get_scale().approx_eq_ratio(&test_scale, 0.00025));
    }
    #[test]
    fn test_utm2ll() {
        let ellipsoid = geo_ellipsoid::geo_ellipsoid::new(geo_ellipsoid::WGS84_SEMI_MAJOR_AXIS_METERS,
                                            geo_ellipsoid::WGS84_FLATTENING);
        let mut utm = utm_grid::utm_grid::new(0, utm_grid::hemisphere::NORTH, 0.0, 0.0, 0.0, 0.0);
        utm.set_zone(33);
        utm.set_hem(utm_grid::hemisphere::NORTH);
        utm.set_easting(431952.612166);
        utm.set_northing(8782098.22289);
        utm.set_convergence(-0.055232079);
        utm.set_scale(0.999656581563);
        let ll_vec = utm2ll(&utm, &ellipsoid);

        let test_lat = 1.3804121468;
        let test_lon = 0.20555336013;
        assert!(ll_vec.x.approx_eq_ratio(&test_lat, 0.00025));
        assert!(ll_vec.y.approx_eq_ratio(&test_lon, 0.00025));
    }
    #[test]
    fn test_lla2ned() {
        let ellipsoid = geo_ellipsoid::geo_ellipsoid::new(geo_ellipsoid::WGS84_SEMI_MAJOR_AXIS_METERS,
    										geo_ellipsoid::WGS84_FLATTENING);
        let lla_orig_vec: Vector3<f64> = Vector3::new(0.8527087756759584, 0.04105401863784606, 1000.000000000);
        let lla_vec: Vector3<f64> = Vector3::new(0.8527087756759584, 0.042799347889836060477, 1000.000000000);
        let ned_vec = lla2ned(&lla_orig_vec, &lla_vec, &ellipsoid);

        let test_x = 4.8231982231937990945;
        let test_y = 7339.3050417820732036;
        let test_z = 4.2139798876589225073;
        assert!(ned_vec.x.approx_eq_ratio(&test_x, 0.00025));
        assert!(ned_vec.y.approx_eq_ratio(&test_y, 0.00025));
        assert!(ned_vec.z.approx_eq_ratio(&test_z, 0.00025));
    }
    #[test]
    fn test_lla2enu() {
        let ellipsoid = geo_ellipsoid::geo_ellipsoid::new(geo_ellipsoid::WGS84_SEMI_MAJOR_AXIS_METERS,
    										geo_ellipsoid::WGS84_FLATTENING);
        let lla_orig_vec: Vector3<f64> = Vector3::new(0.8527087756759584, 0.04105401863784606, 1000.000000000);
        let lla_vec: Vector3<f64> = Vector3::new(0.8527087756759584, 0.042799347889836060477, 1000.000000000);
        let enu_vec = lla2enu(&lla_orig_vec, &lla_vec, &ellipsoid);

        let test_x = 7339.3050417820732036;
        let test_y = 4.8231982231937990945;
        let test_z = -4.2139798876589225073;
        assert!(enu_vec.x.approx_eq_ratio(&test_x, 0.00025));
        assert!(enu_vec.y.approx_eq_ratio(&test_y, 0.00025));
        assert!(enu_vec.z.approx_eq_ratio(&test_z, 0.00025));
    }
    #[test]
    fn test_ned2lla() {
        let ellipsoid = geo_ellipsoid::geo_ellipsoid::new(geo_ellipsoid::WGS84_SEMI_MAJOR_AXIS_METERS,
    										geo_ellipsoid::WGS84_FLATTENING);
        let lla_orig_vec: Vector3<f64> = Vector3::new(0.8527087756759584, 0.04105401863784606, 1000.000000000);
        let ned_vec: Vector3<f64> = Vector3::new(4.8231982231937990945, 7339.3050417820732036, 4.2139798876589225073);
        let lla_vec = ned2lla(&lla_orig_vec, &ned_vec, &ellipsoid);

        let test_x = 0.8527087756759584;
        let test_y = 0.042799347889836060477;
        let test_z = 1000.000000000;
        assert!(lla_vec.x.approx_eq_ratio(&test_x, 0.00025));
        assert!(lla_vec.y.approx_eq_ratio(&test_y, 0.00025));
        assert!(lla_vec.z.approx_eq_ratio(&test_z, 0.00025));  
    }
    #[test]
    fn test_enu2lla() {
        let ellipsoid = geo_ellipsoid::geo_ellipsoid::new(geo_ellipsoid::WGS84_SEMI_MAJOR_AXIS_METERS,
    										geo_ellipsoid::WGS84_FLATTENING);
        let lla_orig_vec: Vector3<f64> = Vector3::new(0.8527087756759584, 0.04105401863784606, 1000.000000000);
        let enu_vec: Vector3<f64> = Vector3::new(7339.3050417820732036, 4.8231982231937990945, -4.2139798876589225073);
        let lla_vec = enu2lla(&lla_orig_vec, &enu_vec, &ellipsoid);

        let test_x = 0.8527087756759584;
        let test_y = 0.042799347889836060477;
        let test_z = 1000.000000000;
        assert!(lla_vec.x.approx_eq_ratio(&test_x, 0.00025));
        assert!(lla_vec.y.approx_eq_ratio(&test_y, 0.00025));
        assert!(lla_vec.z.approx_eq_ratio(&test_z, 0.00025));                            
    }
}
