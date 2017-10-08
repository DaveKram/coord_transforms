#![allow(dead_code)]
#[macro_use] 
extern crate float_cmp;
extern crate nalgebra as na;

/// Module for 3-dimensional coordinate transformations
pub mod d3;
/// Module for 2-dimensional coordinate transformations
pub mod d2;
///Module for geographical coordinate transformations
pub mod geo;
/// Module for constants used in the coordinate transformations
pub mod structs;