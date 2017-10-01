use ::std;
use na::Vector2;

/// Converts 2-d polar coordinates to 2-d cartesian coordinates
/// 
/// # Arguments
/// 
/// * `pol_vec` - Vector2 reference to the polar vector (r, theta) in radians
/// 
/// # Return Value
/// 
/// * nalgebra::Vector2<f64> - x, y
pub fn polar2cartesian(logpol_vec: &Vector2<f64>) -> Vector2<f64> {
	let mut ret_vec: Vector2<f64> = Vector2::new(0.0, 0.0);
	ret_vec.x = logpol_vec.x * logpol_vec.y.cos();
	ret_vec.y = logpol_vec.x * logpol_vec.y.sin();
	ret_vec
}

/// Converts 2-d log polar coordinates to 2-d cartesian coordinates
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

/// Converts 2-d bipolar coordinates to 2-d cartesian coordinates
/// 
/// # Arguments
/// 
/// * `bipol_vec` - Vector2 reference to the bipolar vector (sigma, tau) in radians
/// * `a` - f64 value for foci points (-a, 0) and (a, 0)
/// 
/// # Return Value
/// 
/// * nalgebra::Vector2<f64> - x, y
pub fn bipolar2cartesian(logpol_vec: &Vector2<f64>, a: f64) -> Vector2<f64> {
	let mut ret_vec: Vector2<f64> = Vector2::new(0.0, 0.0);
	ret_vec.x = a * ((logpol_vec.y.sinh()) / (logpol_vec.y.cosh() - logpol_vec.x.cos()));
	ret_vec.y = a * ((logpol_vec.x.sin()) / (logpol_vec.y.cosh() - logpol_vec.x.cos()));
	ret_vec
}

/// Converts 2-d cartesian coordinates to 2-d log polar coordinates
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

//Unit tests
#[cfg(test)]
mod tests {
	use super::*;
	#[test]
	fn test_polar2cartesian() {
        let pol_vec: Vector2<f64> = Vector2::new(3.0, 4.0);
        let cart_vec = polar2cartesian(&pol_vec);
        assert_approx_eq!(cart_vec.x, -1.960930862590836);
        assert_approx_eq!(cart_vec.y, -2.2704074859237844);
    }
    #[test]
    fn test_logpolar2cartesian() {
        let logpol_vec: Vector2<f64> = Vector2::new(3.0, 4.0);
        let cart_vec = logpolar2cartesian(&logpol_vec);
        assert_approx_eq!(cart_vec.x, -13.128783081462156);
        assert_approx_eq!(cart_vec.y, -15.20078446306795);
    }
    #[test]
    fn test_bipolar2cartesian() {
        let bipol_vec: Vector2<f64> = Vector2::new(3.0, 4.0);
        let cart_vec = bipolar2cartesian(&bipol_vec, 1.0);
        assert_approx_eq!(cart_vec.x, 0.9643685028429331);
        assert_approx_eq!(cart_vec.y, 0.0049869);
    }
	#[test]
    fn test_cartesian2logpolar() {
        let cart_vec: Vector2<f64> = Vector2::new(3.0, 4.0);
        let logpolar_vec = cartesian2logpolar(&cart_vec);
        assert_approx_eq!(logpolar_vec.x, 1.6094379124341003);
        assert_approx_eq!(logpolar_vec.y, 0.9272952180016122);
    }
}