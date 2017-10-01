use na::Vector3;

/// Converts 3-d cartesian coordinates to 3-d spherical coordinates
/// 
/// # Arguments
/// 
/// * `cart_vec` - Vector3 reference to the cartesian vector (x, y, z)
/// 
/// # Return Value
/// 
/// * nalgebra::Vector3<f64> - rho, theta, phi (in radians)
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

/// Converts 3-d spherical coordinates to 3-d cartesian coordinates
/// 
/// # Arguments
/// 
/// * `sphere_vec` - Vector3 reference to the spherical vector (rho, theta, phi) in radians
/// 
/// # Return Value
/// 
/// * nalgebra::Vector3<f64> - x, y, z
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


/// Converts 3-d ENU coordinates to 3-d NED coordinates
/// 
/// # Arguments
/// 
/// * `enu_vec` - Vector3 reference to the ENU vector (x, y, z)
/// 
/// # Return Value
/// 
/// * nalgebra::Vector3<f64> - x, y, z
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
/// * nalgebra::Vector3<f64> - x, y, z
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

//Unit tests
#[cfg(test)]
mod tests {
	use super::*;
    #[test]
    fn test_cartesian2spherical() {
        let cart_vec: Vector3<f64> = Vector3::new(3.0, 4.0, 5.0);
        let sphere_vec = cartesian2spherical(&cart_vec);
        assert_approx_eq!(sphere_vec.x, 7.0710678118655);
        assert_approx_eq!(sphere_vec.y, 0.78539816339745);
        assert_approx_eq!(sphere_vec.z, 0.92729521800161);
    }
    #[test]
    fn test_spherical2cartesian() {
        let sphere_vec: Vector3<f64> = Vector3::new(3.0, 4.0, 5.0);
        let cart_vec = spherical2cartesian(&sphere_vec);
        assert_approx_eq!(cart_vec.x, -0.6440287493492097);
        assert_approx_eq!(cart_vec.y, 2.177148851629225);
        assert_approx_eq!(cart_vec.z, -1.960930862590836);
    }
    #[test]
    fn test_enu2ned() {
        let enu_vec: Vector3<f64> = Vector3::new(3.0, 4.0, 5.0);
        let ned_vec = enu2ned(&enu_vec);
        assert_eq!(ned_vec.x, 4.0);
        assert_eq!(ned_vec.y, 3.0);
        assert_eq!(ned_vec.z, -5.0);
    }
    #[test]
    fn test_ned2enu() {
        let ned_vec: Vector3<f64> = Vector3::new(3.0, 4.0, 5.0);
        let enu_vec = ned2enu(&ned_vec);
        assert_eq!(enu_vec.x, 4.0);
        assert_eq!(enu_vec.y, 3.0);
        assert_eq!(enu_vec.z, -5.0);
    }
}