use na::Vector3;

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
    fn test_spherical2cartesian() {
        let sphere_vec: Vector3<f64> = Vector3::new(5.0, 0.5, 0.2);
        let cart_vec = spherical2cartesian(&sphere_vec);
        assert_approx_eq!(cart_vec.x, 0.8717437014);
        assert_approx_eq!(cart_vec.y, 0.4762357546);
        assert_approx_eq!(cart_vec.z, 4.900332889);
    }
}