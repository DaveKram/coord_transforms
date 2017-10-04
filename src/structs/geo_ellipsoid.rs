pub const WGS84_SEMI_MAJOR_AXIS_METERS: f64 = 6378137.0;
pub const WGS84_FLATTENING: f64 = 298.257223563;

pub struct geo_ellipsoid {
	semi_major_axis: f64,
	flattening: f64,
	semi_minor_axis: f64,
	first_ecc: f64,
	second_ecc: f64
}

impl geo_ellipsoid {
	pub fn new(sma: f64, f: f64) -> geo_ellipsoid {
		let smia = sma * (1.0 - (1.0 / f));
		geo_ellipsoid {
			semi_major_axis: sma,
			flattening: 1.0 / f,
			semi_minor_axis: smia,
			first_ecc: ((sma.powi(2) - smia.powi(2)) / (sma.powi(2))).sqrt(),
			second_ecc: ((sma.powi(2) - smia.powi(2)) / (smia.powi(2))).sqrt(),
		}
	}

	pub fn get_semi_major_axis(&self) -> f64 {
		self.semi_major_axis
	}

	pub fn get_flattening(&self) -> f64 {
		self.flattening
	}	

	pub fn get_semi_minor_axis(&self) -> f64 {
		self.semi_minor_axis
	}

	pub fn get_first_ecc(&self) -> f64 {
		self.first_ecc
	}

	pub fn get_second_ecc(&self) -> f64 {
		self.second_ecc
	}
}
