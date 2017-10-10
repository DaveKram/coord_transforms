pub const SCALE_FACTOR_CENTERAL_MERIDIAN: f64 = 0.9996;
pub const FALSE_EASTING: f64 = 500000.0;
pub const FALSE_NORTHING: f64 = 10000000.0;

#[derive(Clone, Debug)]
#[allow(non_camel_case_types)]
pub enum hemisphere {
	NORTH,
	SOUTH
}

#[allow(non_camel_case_types)]
pub struct utm_grid {
	zone: u32,
	hem: hemisphere,
	easting: f64,
	northing: f64,
	convergence: f64,
	scale: f64
}

impl utm_grid {
	pub fn new() -> utm_grid {
		utm_grid {
			zone: 0,
			hem: hemisphere::NORTH,
			easting: 0.0,
			northing: 0.0,
			convergence: 0.0,
			scale: 0.0
		}
	}

	pub fn set_zone(&mut self, z: u32) {
		self.zone = z;
	}

	pub fn set_hem(&mut self, h: hemisphere) {
		self.hem = h;
	}

	pub fn set_easting(&mut self, e: f64) {
		self.easting = e;
	}

	pub fn set_northing(&mut self, n: f64) {
		self.northing = n;
	}

	pub fn set_convergence(&mut self, n: f64) {
		self.convergence = n;
	}

	pub fn set_scale(&mut self, n: f64) {
		self.scale = n;
	}

	pub fn get_zone(&self) -> u32 {
		self.zone
	}

	pub fn get_hem(&self) -> hemisphere {
		self.hem.to_owned()
	}

	pub fn get_easting(&self) -> f64 {
		self.easting
	}

	pub fn get_northing(&self) -> f64 {
		self.northing
	}

	pub fn get_convergence(&self) -> f64 {
		self.convergence
	}

	pub fn get_scale(&self) -> f64 {
		self.scale
	}
}