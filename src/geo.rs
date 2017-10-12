use na::Vector3;
use na::Vector2;
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


/// Converts 2-d LL coordinates to UTM Grid (using Karney method) - accruacy to within a few nanometers within 3900km of the central meridian
/// 
/// # Arguments
/// 
/// * `ll_vec` - Vector2 reference to the LL vector (latitude, longitude) (radians, radians)
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
    let mut ret_utm = utm_grid::utm_grid::new();
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
        assert!(utm.get_easting().approx_eq_ratio(&test_easting, 0.0000000001));
        assert!(utm.get_northing().approx_eq_ratio(&test_northing, 0.0000000001));
        assert!(utm.get_convergence().approx_eq_ratio(&test_convergence, 0.0000000001));
        assert!(utm.get_scale().approx_eq_ratio(&test_scale, 0.0000000001));
    }
}
