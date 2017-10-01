use ::std;
use na::Vector2;

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
    fn test_cartesian2logpolar() {
        let cart_vec: Vector2<f64> = Vector2::new(3.0, 4.0);
        let logpolar_vec = cartesian2logpolar(&cart_vec);
        assert_approx_eq!(logpolar_vec.x, 1.6094379124341003);
        assert_approx_eq!(logpolar_vec.y, 0.9272952180016122);
    }
    #[test]
    fn test_logpolar2cartesian() {
        let logpol_vec: Vector2<f64> = Vector2::new(1.6094379124341003, 0.9272952180016122);
        let cart_vec = logpolar2cartesian(&logpol_vec);
        assert_approx_eq!(cart_vec.x, 3.0);
        assert_approx_eq!(cart_vec.y, 4.0);
    }
}