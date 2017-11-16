#![allow(dead_code)]
#[allow(unused_imports)]
#[macro_use] 
extern crate float_cmp;
extern crate nalgebra as na;

/// Module for 3-dimensional coordinate transformations
pub mod d3;
/// Module for 2-dimensional coordinate transformations
pub mod d2;
/// Module for geographical coordinate transformations
pub mod geo;
/// Module for constants used in the coordinate transformations
pub mod structs;

/// Prelude module for ease of use
pub mod prelude {
	pub use na::Vector3;
	pub use na::Vector2;
	pub use d3;
	pub use d2;
	pub use geo;
	pub use structs;
}