use ::std;
use na::Vector2;

/// Converts 2-d polar coordinates to 2-d cartesian coordinates
/// 
/// # Arguments
/// 
/// * `pol_vec` - Vector2 reference to the polar vector (rho, theta) - theta in radians
/// 
/// # Return Value
/// 
/// * `nalgebra::Vector2<f64>` - x, y
/// 
/// # Formula
/// 
/// * x = rho * cos(theta)
/// * y = rho * sin(theta)
pub fn polar2cartesian(pol_vec: &Vector2<f64>) -> Vector2<f64> {
	let mut ret_vec: Vector2<f64> = Vector2::new(0.0, 0.0);
	ret_vec.x = pol_vec.x * pol_vec.y.cos();
	ret_vec.y = pol_vec.x * pol_vec.y.sin();
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
/// * `nalgebra::Vector2<f64>` - x, y
/// 
/// # Formula
/// 
/// * x = e^rho * cos(theta)
/// * y = e^rho * sin(theta)
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
/// * `nalgebra::Vector2<f64>` - x, y
/// 
/// # Formula
/// 
/// * x = a * ((sinh(tau)) / (cosh(tau) - cos(sigma)))
/// * y = a * ((sin(sigma)) / (cosh(tau) - cos(sigma)))
pub fn bipolar2cartesian(bipol_vec: &Vector2<f64>, a: f64) -> Vector2<f64> {
	let mut ret_vec: Vector2<f64> = Vector2::new(0.0, 0.0);
	ret_vec.x = a * ((bipol_vec.y.sinh()) / (bipol_vec.y.cosh() - bipol_vec.x.cos()));
	ret_vec.y = a * ((bipol_vec.x.sin()) / (bipol_vec.y.cosh() - bipol_vec.x.cos()));
	ret_vec
}


/// Converts 2-d cartesian coordinates to 2-d polar coordinates
/// 
/// # Arguments
/// 
/// * `cart_vec` - Vector2 reference to the cartesian vector (x, y)
/// 
/// # Return Value
/// 
/// * `nalgebra::Vector2<f64>` - rho, theta (in radians)
/// 
/// # Formula
/// 
/// * r = sqrt( x^2 + y^2 )
/// * theta = arctan(y / x)
pub fn cartesian2polar(cart_vec: &Vector2<f64>) -> Vector2<f64> {
	let mut ret_vec: Vector2<f64> = Vector2::new(0.0, 0.0);
	ret_vec.x = (cart_vec.x.powi(2) + cart_vec.y.powi(2)).sqrt();
	ret_vec.y = cart_vec.y.atan2(cart_vec.x);
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
/// * `nalgebra::Vector2<f64>` - rho, theta (in radians)
/// 
/// # Formula
/// 
/// * r = log(sqrt( x^2 + y^2 ))
/// * theta = arctan(y / x)
pub fn cartesian2logpolar(cart_vec: &Vector2<f64>) -> Vector2<f64> {
	let mut ret_vec: Vector2<f64> = Vector2::new(0.0, 0.0);
	ret_vec.x = ((cart_vec.x.powi(2) + cart_vec.y.powi(2)).sqrt()).ln();
	ret_vec.y = cart_vec.y.atan2(cart_vec.x);
	ret_vec
}

//Unit tests
#[cfg(test)]
mod tests {
	use super::*;
    use float_cmp::ApproxEqUlps;
	#[test]
	fn test_polar2cartesian() {
        let pol_vec: Vector2<f64> = Vector2::new(3.0, 4.0);
        let cart_vec = polar2cartesian(&pol_vec);

        let test_x = -1.960930862590836;
        let test_y = -2.2704074859237844;
        assert!(cart_vec.x.approx_eq_ulps(&test_x, 2));
        assert!(cart_vec.y.approx_eq_ulps(&test_y, 2));
    }
    #[test]
    fn test_logpolar2cartesian() {
        let logpol_vec: Vector2<f64> = Vector2::new(3.0, 4.0);
        let cart_vec = logpolar2cartesian(&logpol_vec);

        let test_x = -13.128783081462156;
        let test_y = -15.20078446306795;
        assert!(cart_vec.x.approx_eq_ulps(&test_x, 2));
        assert!(cart_vec.y.approx_eq_ulps(&test_y, 2));
    }
    #[test]
    fn test_bipolar2cartesian() {
        let bipol_vec: Vector2<f64> = Vector2::new(3.0, 4.0);
        let cart_vec = bipolar2cartesian(&bipol_vec, 1.0);

        let test_x = 0.9643685028429331;
        let test_y = 0.004986885446035738;
        assert!(cart_vec.x.approx_eq_ulps(&test_x, 2));
        assert!(cart_vec.y.approx_eq_ulps(&test_y, 2));
    }
    #[test]
    fn test_cartesian2polar() {
        let cart_vec: Vector2<f64> = Vector2::new(3.0, 4.0);
        let polar_vec = cartesian2polar(&cart_vec);

        let test_x = 5.0;
        let test_y = 0.9272952180016122;
        assert!(polar_vec.x.approx_eq_ulps(&test_x, 2));
        assert!(polar_vec.y.approx_eq_ulps(&test_y, 2));
    }
	#[test]
    fn test_cartesian2logpolar() {
        let cart_vec: Vector2<f64> = Vector2::new(3.0, 4.0);
        let logpolar_vec = cartesian2logpolar(&cart_vec);

        let test_x = 1.6094379124341003;
        let test_y = 0.9272952180016122;
        assert!(logpolar_vec.x.approx_eq_ulps(&test_x, 2));
        assert!(logpolar_vec.y.approx_eq_ulps(&test_y, 2));
    }
}
