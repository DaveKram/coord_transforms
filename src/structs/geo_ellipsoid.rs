pub const WGS84_SEMI_MAJOR_AXIS_METERS: f64 = 6378137.0;
pub const WGS84_FLATTENING: f64 = 298.257223563;
pub const WGS72_SEMI_MAJOR_AXIS_METERS: f64 = 6378135.0;
pub const WGS72_FLATTENING: f64 = 298.26;
pub const WGS66_SEMI_MAJOR_AXIS_METERS: f64 = 6378145.0;
pub const WGS66_FLATTENING: f64 = 298.25;
pub const WGS60_SEMI_MAJOR_AXIS_METERS: f64 = 6378165.0;
pub const WGS60_FLATTENING: f64 = 298.3;
pub const SOUTHAMERICAN_1969_SEMI_MAJOR_AXIS_METERS: f64 = 6378160.0;
pub const SOUTHAMERICAN_1969_FLATTENING: f64 = 298.25;
pub const FISCHER_1969_MODIFIED_SEMI_MAJOR_AXIS_METERS: f64 = 6378155.0;
pub const FISCHER_1969_MODIFIED_FLATTENING: f64 = 298.3;
pub const EVEREST_MODIFIED_SEMI_MAJOR_AXIS_METERS: f64 = 6377304.063;
pub const EVEREST_MODIFIED_FLATTENING: f64 = 300.8017;
pub const AIRY_MODIFIED_SEMI_MAJOR_AXIS_METERS: f64 = 6377340.189;
pub const AIRY_MODIFIED_FLATTENING: f64 = 299.3249646;
pub const KRASSOVSKY_SEMI_MAJOR_AXIS_METERS: f64 = 6378245.0;
pub const KRASSOVSKY_FLATTENING: f64 = 298.3;
pub const INTERNATIONAL_SEMI_MAJOR_AXIS_METERS: f64 = 6378388.0;
pub const INTERNATIONAL_FLATTENING: f64 = 297.0;
pub const HOUGH_SEMI_MAJOR_AXIS_METERS: f64 = 6378270.0;
pub const HOUGH_FLATTENING: f64 = 297.0;
pub const HELMERT_1906_SEMI_MAJOR_AXIS_METERS: f64 = 6378200.0;
pub const HELMERT_1906_FLATTENING: f64 = 298.3;
pub const GRS_1980_SEMI_MAJOR_AXIS_METERS: f64 = 6378137.0;
pub const GRS_1980_FLATTENING: f64 = 298.257222101;
pub const GRS_1967_SEMI_MAJOR_AXIS_METERS: f64 = 6378160.0;
pub const GRS_1967_FLATTENING: f64 = 298.247167427;
pub const FISCHER_1968_SEMI_MAJOR_AXIS_METERS: f64 = 6378150.0;
pub const FISCHER_1968_FLATTENING: f64 = 298.3;
pub const FISCHER_1960_MERCURY_SEMI_MAJOR_AXIS_METERS: f64 = 6378166.0;
pub const FISCHER_1960_MERCURY_FLATTENING: f64 = 298.3;
pub const EVEREST_SEMI_MAJOR_AXIS_METERS: f64 = 6377276.345;
pub const EVEREST_FLATTENING: f64 = 300.8017;
pub const CLARKE_1880_SEMI_MAJOR_AXIS_METERS: f64 = 6378249.145;
pub const CLARKE_1880_FLATTENING: f64 = 293.465;
pub const CLARKE_1866_SEMI_MAJOR_AXIS_METERS: f64 = 6378206.4;
pub const CLARKE_1866_FLATTENING: f64 = 294.9786982;
pub const BESSEL_1841_NAMBIA_SEMI_MAJOR_AXIS_METERS: f64 = 6377483.865;
pub const BESSEL_1841_NAMBIA_FLATTENING: f64 = 299.1528128;
pub const BESSEL_1841_SEMI_MAJOR_AXIS_METERS: f64 = 6377397.155;
pub const BESSEL_1841_FLATTENING: f64 = 299.1528128;
pub const AUSTRALIAN_NATIONAL_SEMI_MAJOR_AXIS_METERS: f64 = 6378160.0;
pub const AUSTRALIAN_NATIONAL_FLATTENING: f64 = 298.25;
pub const AIRY_SEMI_MAJOR_AXIS_METERS: f64 = 6377563.396;
pub const AIRY_FLATTENING: f64 = 299.3249646;

#[allow(non_camel_case_types)]
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
