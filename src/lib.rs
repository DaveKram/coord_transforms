#[macro_use] extern crate assert_approx_eq;
extern crate nalgebra as na;
use na::Vector3;
use na::Vector2;

/// Converts cartesian coordinates to spherical coordinates
/// 
/// # Arguments
/// 
/// * `cart_vec` - Vector3 reference to the cartesian vector (x, y, z)
/// 
/// # Return Value
/// 
/// * nalgebra::Vector3<f64> - r, theta, phi (in radians)
pub fn cartesian2spherical(cart_vec: &Vector3<f64>) -> Vector3<f64> {
	let mut ret_vec: Vector3<f64> = Vector3::new(0.0, 0.0, 0.0);
	ret_vec.x = (cart_vec.x.powi(2) + cart_vec.y.powi(2) + cart_vec.z.powi(2)).sqrt();
	ret_vec.y = (cart_vec.y / cart_vec.x).atan();
	ret_vec.z = ((cart_vec.x.powi(2) + cart_vec.y.powi(2)).sqrt() / cart_vec.z).atan();
	ret_vec
}


/// Converts cartesian coordinates to log polar coordinates
/// 
/// # Arguments
/// 
/// * `cart_vec` - Vector2 reference to the cartesian vector (x, y)
/// 
/// # Return Value
/// 
/// * nalgebra::Vector2<f64> - rho, theta (in radians)
pub fn cartesian2logpolar(cart_vec: &Vector2<f64>) -> Vector2<f64> {
	let mut ret_vec: Vector2<f64> = Vector2::new(0.0, 0.0);
	ret_vec.x = ((cart_vec.x.powi(2) + cart_vec.y.powi(2)).sqrt()).ln();
	ret_vec.y = (cart_vec.y / cart_vec.x).atan();
	ret_vec
}

/// Converts spherical coordinates to cartesian coordinates
/// 
/// # Arguments
/// 
/// * `sphere_vec` - Vector3 reference to the spherical vector (r, theta, phi) in radians
/// 
/// # Return Value
/// 
/// * nalgebra::Vector3<f64> - x, y, z
pub fn spherical2cartesian(sphere_vec: &Vector3<f64>) -> Vector3<f64> {
	let mut ret_vec: Vector3<f64> = Vector3::new(0.0, 0.0, 0.0);
	ret_vec.x = sphere_vec.x * sphere_vec.z.sin() * sphere_vec.y.cos();
	ret_vec.y = sphere_vec.x * sphere_vec.z.sin() * sphere_vec.y.sin();
	ret_vec.z = sphere_vec.x * sphere_vec.z.cos();
	ret_vec
}

/// Converts log polar coordinates to cartesian coordinates
/// 
/// # Arguments
/// 
/// * `logpol_vec` - Vector2 reference to the log polar vector (rho, theta) in radians
/// 
/// # Return Value
/// 
/// * nalgebra::Vector2<f64> - x, y
pub fn logpolar2cartesian(logpol_vec: &Vector2<f64>) -> Vector2<f64> {
	let mut ret_vec: Vector2<f64> = Vector2::new(0.0, 0.0);
	ret_vec.x = std::f64::consts::E.powf(logpol_vec.x) * logpol_vec.y.cos();
	ret_vec.y = std::f64::consts::E.powf(logpol_vec.x) * logpol_vec.y.sin();
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
        assert_approx_eq!(sphere_vec.y, 0.92729521800161);
        assert_approx_eq!(sphere_vec.z, 0.78539816339745);
    }
    #[test]
    fn test_cartesian2logpolar() {
        let cart_vec: Vector2<f64> = Vector2::new(3.0, 4.0);
        let logpolar_vec = cartesian2logpolar(&cart_vec);
        assert_approx_eq!(logpolar_vec.x, 1.6094379124341003);
        assert_approx_eq!(logpolar_vec.y, 0.9272952180016122);
    }
    #[test]
    fn test_spherical2cartesian() {
        let sphere_vec: Vector3<f64> = Vector3::new(5.0, 0.5, 0.2);
        let cart_vec = spherical2cartesian(&sphere_vec);
        assert_approx_eq!(cart_vec.x, 0.8717437014);
        assert_approx_eq!(cart_vec.y, 0.4762357546);
        assert_approx_eq!(cart_vec.z, 4.900332889);
    }
    #[test]
    fn test_logpolar2cartesian() {
        let logpol_vec: Vector2<f64> = Vector2::new(1.6094379124341003, 0.9272952180016122);
        let cart_vec = logpolar2cartesian(&logpol_vec);
        assert_approx_eq!(cart_vec.x, 3.0);
        assert_approx_eq!(cart_vec.y, 4.0);
    }
}