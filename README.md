## Synopsis

A rust crate use for performing coordinate transformations. Relies on the nalgebra crate for its structures. Currently uses all f64 values. Still very early work in progress.

## Code Example

```
let cart_vec: Vector3<f64> = Vector3::new(3.0, 4.0, 5.0);
let sphere_vec = cartesian2spherical(&cart_vec);
assert_approx_eq!(sphere_vec.x, 7.0710678118655);
assert_approx_eq!(sphere_vec.y, 0.92729521800161);
assert_approx_eq!(sphere_vec.z, 0.78539816339745);
```

## Motivation

Providing simple helper functions for coordinate transforms seemed like a simple crate to create for the Rust community. With unit tested code, the crate will provide a stable baseline for common coordinate transformations.

## Tests

The crate has unit tests built in. If you want to verify this yourself, run:

```
cargo test
```

## License

Copyright (c) 2017 David Kramer

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
