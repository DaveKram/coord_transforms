extern crate coord_transforms;
extern crate nalgebra as na;
extern crate minifb;
use coord_transforms::prelude::*;
use minifb::{Key, WindowOptions, Window};

const WIDTH: usize = 800;
const HEIGHT: usize = 600;

fn main() {
	let mut buffer: Vec<u32> = vec![0; WIDTH * HEIGHT];

    let mut window = Window::new("Test - Visual - coord_transforms",
                                 WIDTH,
                                 HEIGHT,
                                 WindowOptions::default()).unwrap_or_else(|e| {
        panic!("{}", e);
    });

    let mut line_angle_rads = 0.0;
    while window.is_open() && !window.is_key_down(Key::Escape) {
    	//Circle code
    	let radius: i64 = 200;
    	let x0: i64 = WIDTH as i64 / 2;
    	let y0: i64 = HEIGHT as i64 / 2;
    	let mut x: i64 = radius - 1;
    	let mut y: i64 = 0;
    	let mut dx: i64 = 1;
    	let mut dy: i64 = 1;
    	let mut err: i64 = dx - (radius << 1);
    	let colorPix = 0xFFFFFF;

    	//Line code
    	let mut p1x: i64 = WIDTH as i64 / 2;
    	let mut p1y: i64 = HEIGHT as i64 / 2;
    	let pol_vec: Vector2<f64> = Vector2::new(radius as f64, line_angle_rads);
    	let cart_vec = d2::polar2cartesian(&pol_vec);
    	let mut p2x: i64 = p1x + cart_vec.x as i64;
    	let mut p2y: i64 = p1y + cart_vec.y as i64;
    	let mut psteep = false;
    	if (p2y - p1y).abs() > (p2x - p1x).abs() {
    		psteep = true;
    	}
    	if psteep {
    		let tempx = p1x;
    		p1x = p1y;
    		p1y = tempx;
    		let tempx2 = p2x;
    		p2x = p2y;
    		p2y = tempx2;
    	}
    	if p1x > p2x {
    		let tempx = p1x;
    		p1x = p2x;
    		p2x = tempx;
    		let tempy = p1y;
    		p1y = p2y;
    		p2y = tempy;
    	}
    	let pdx = p2x - p1x;
    	let pdy = (p2y - p1y).abs();
    	let mut perror = pdx / 2;
    	let mut pystep = 0;
    	if p1y < p2y {
    		pystep = 1;
    	}else{
    		pystep = -1;
    	}
    	let mut py0 = p1y;
    	let pmaxx = p2x;

    	//Draw circle
    	while x >= y {
        	buffer[xy2lin(x0 + x, y0 + y) as usize] = colorPix;
        	buffer[xy2lin(x0 + y, y0 + x) as usize] = colorPix;
        	buffer[xy2lin(x0 - y, y0 + x) as usize] = colorPix;
        	buffer[xy2lin(x0 - x, y0 + y) as usize] = colorPix;
        	buffer[xy2lin(x0 - x, y0 - y) as usize] = colorPix;
        	buffer[xy2lin(x0 - y, y0 - x) as usize] = colorPix;
        	buffer[xy2lin(x0 + y, y0 - x) as usize] = colorPix;
        	buffer[xy2lin(x0 + x, y0 - y) as usize] = colorPix;
            if err <= 0 {
            	y += 1;
            	err += dy;
            	dy += 2;
            }
            if err > 0 {
            	x -= 1;
            	dx += 2;
            	err += (-radius << 1) + dx;
            }
		}

		//Draw line
		for x in p1x..pmaxx {
			if psteep {
				buffer[xy2lin(py0, x) as usize] = colorPix;
			}else{
				buffer[xy2lin(x, py0) as usize] = colorPix;
			}
			perror -= pdy;
			if perror < 0 {
				py0 += pystep;
				perror += pdx;
			}
		}

        window.update_with_buffer(&buffer, WIDTH, HEIGHT).unwrap();

    	//Wrap radians
    	if line_angle_rads >= 2.0 * 3.14 {
    		line_angle_rads = 0.0;
    	}else{
    		line_angle_rads += 0.001;
    	}

    	//Clear
    	for x in 0..WIDTH * HEIGHT {
    		buffer[x] = 0x0;
    	}
	}
}

fn xy2lin(x: i64, y: i64) -> i64 {
	WIDTH as i64 * y + x
}
