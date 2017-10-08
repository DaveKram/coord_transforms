use na::Vector3;

/// Converts 3-d spherical coordinates to 3-d cartesian coordinates
/// 
/// # Arguments
/// 
/// * `sphere_vec` - Vector3 reference to the spherical vector (rho, theta, phi) (r, el, az) in radians
/// 
/// # Return Value
/// 
/// * `nalgebra::Vector3<f64>` - x, y, z
/// 
/// # Formula
/// 
/// * x = rho * sin(theta) * cos(phi)
/// * y = rho * sin(theta) * sin(phi)
/// * z = rho * cos(theta)
pub fn spherical2cartesian(sphere_vec: &Vector3<f64>) -> Vector3<f64> {
    let mut ret_vec: Vector3<f64> = Vector3::new(0.0, 0.0, 0.0);
    ret_vec.x = sphere_vec.x * sphere_vec.y.sin() * sphere_vec.z.cos();
    ret_vec.y = sphere_vec.x * sphere_vec.y.sin() * sphere_vec.z.sin();
    ret_vec.z = sphere_vec.x * sphere_vec.y.cos();
    ret_vec
}

/// Converts 3-d cylindrical coordinates to 3-d cartesian coordinates
/// 
/// # Arguments
/// 
/// * `cyl_vec` - Vector3 reference to the cylindrical vector (rho, theta, z) in radians
/// 
/// # Return Value
/// 
/// * `nalgebra::Vector3<f64>` - x, y, z
/// 
/// # Formula
/// 
/// * x = rho * cos(theta)
/// * y = rho * sin(theta)
/// * z = z
pub fn cylindrical2cartesian(cyl_vec: &Vector3<f64>) -> Vector3<f64> {
    let mut ret_vec: Vector3<f64> = Vector3::new(0.0, 0.0, 0.0);
    ret_vec.x = cyl_vec.x * cyl_vec.y.cos();
    ret_vec.y = cyl_vec.x * cyl_vec.y.sin();
    ret_vec.z = cyl_vec.z;
    ret_vec
}

/// Converts 3-d cartesian coordinates to 3-d spherical coordinates
/// 
/// # Arguments
/// 
/// * `cart_vec` - Vector3 reference to the cartesian vector (x, y, z)
/// 
/// # Return Value
/// 
/// * `nalgebra::Vector3<f64>` - rho, theta, phi (in radians)
/// 
/// # Formula
/// 
/// * rho = sqrt( x^2 + y^2 + z^2 )
/// * theta = arctan((sqrt( x2 + y^2 )) / (z))
/// * phi = arctan(y / x)
pub fn cartesian2spherical(cart_vec: &Vector3<f64>) -> Vector3<f64> {
	let mut ret_vec: Vector3<f64> = Vector3::new(0.0, 0.0, 0.0);
	ret_vec.x = (cart_vec.x.powi(2) + cart_vec.y.powi(2) + cart_vec.z.powi(2)).sqrt();
	ret_vec.y = ((cart_vec.x.powi(2) + cart_vec.y.powi(2)).sqrt()).atan2(cart_vec.z); 
	ret_vec.z = cart_vec.y.atan2(cart_vec.x);
	ret_vec
}

/// Converts 3-d cartesian coordinates to 3-d cylindrical coordinates
/// 
/// # Arguments
/// 
/// * `cart_vec` - Vector3 reference to the cartesian vector (x, y, z)
/// 
/// # Return Value
/// 
/// * `nalgebra::Vector3<f64>` - rho, theta, z (in radians)
/// 
/// # Formula
/// 
/// * rho = sqrt( x^2 + y^2 )
/// * theta = arctan(y / x)
/// * z = z
pub fn cartesian2cylindrical(cart_vec: &Vector3<f64>) -> Vector3<f64> {
    let mut ret_vec: Vector3<f64> = Vector3::new(0.0, 0.0, 0.0);
    ret_vec.x = (cart_vec.x.powi(2) + cart_vec.y.powi(2)).sqrt();
    ret_vec.y = cart_vec.y.atan2(cart_vec.x);
    ret_vec.z = cart_vec.z;
    ret_vec
}

//Unit tests
#[cfg(test)]
mod tests {
	use super::*;
    use float_cmp::ApproxEqUlps;
    #[test]
    fn test_spherical2cartesian() {
        let sphere_vec: Vector3<f64> = Vector3::new(3.0, 4.0, 5.0);
        let cart_vec = spherical2cartesian(&sphere_vec);

        let test_x = -0.6440287493492097;
        let test_y = 2.177148851629225;
        let test_z = -1.960930862590836;
        assert!(cart_vec.x.approx_eq_ulps(&test_x, 2));
        assert!(cart_vec.y.approx_eq_ulps(&test_y, 2));
        assert!(cart_vec.z.approx_eq_ulps(&test_z, 2));
    }
    #[test]
    fn test_cylindrical2cartesian() {
        let cyl_vec: Vector3<f64> = Vector3::new(3.0, 4.0, 5.0);
        let cart_vec = cylindrical2cartesian(&cyl_vec);

        let test_x = -1.960930862590836;
        let test_y = -2.2704074859237844;
        let test_z = 5.0;
        assert!(cart_vec.x.approx_eq_ulps(&test_x, 2));
        assert!(cart_vec.y.approx_eq_ulps(&test_y, 2));
        assert!(cart_vec.z.approx_eq_ulps(&test_z, 2));
    }
    #[test]
    fn test_cartesian2spherical() {
        let cart_vec: Vector3<f64> = Vector3::new(3.0, 4.0, 5.0);
        let sphere_vec = cartesian2spherical(&cart_vec);

        let test_x = 7.0710678118654755;
        let test_y = 0.7853981633974483;
        let test_z = 0.9272952180016122;
        assert!(sphere_vec.x.approx_eq_ulps(&test_x, 2));
        assert!(sphere_vec.y.approx_eq_ulps(&test_y, 2));
        assert!(sphere_vec.z.approx_eq_ulps(&test_z, 2));
    }
    #[test]
    fn test_cartesian2cylindrical() {
        let cart_vec: Vector3<f64> = Vector3::new(3.0, 4.0, 5.0);
        let cyl_vec = cartesian2cylindrical(&cart_vec);

        let test_x = 5.0;
        let test_y = 0.9272952180016122;
        let test_z = 5.0;
        assert!(cyl_vec.x.approx_eq_ulps(&test_x, 2));
        assert!(cyl_vec.y.approx_eq_ulps(&test_y, 2));
        assert!(cyl_vec.z.approx_eq_ulps(&test_z, 2));
    }
}
