use crate::util::vector::Vector;
use crate::util::xform::Xform;
use crate::engine::bullet_math::{
	fmaxss, fminss, cross, sse_normalize,
};

#[derive(Clone, Copy, PartialEq)]
pub enum GeomShape {
	Ball,
	Cylinder,
	Block,
	Wedge,
	CornerWedge
}

pub struct Shape {
	pub geom: GeomShape,
	pub half: Vector,
	pub margin: f32
}

pub fn make_shape(a1: GeomShape, a2: Vector) -> Shape {
	let x = fmaxss(a2.x, 0.05) * 0.5;
	let mut y = fmaxss(a2.y, 0.05) * 0.5;
	let mut z = fmaxss(a2.z, 0.05) * 0.5;

	if a1 == GeomShape::Ball {
		let bz = fminss(z, fminss(y, x));
		let bx = fminss(x, z);
		let by = fminss(y, z);
		return Shape {
			geom: GeomShape::Ball,
			half: Vector::new(bx, by, bz),
			margin: z
		};
	}

	if a1 == GeomShape::Cylinder {
		y = fminss(z, y);
		z = y;
	}
	return Shape {
		geom: a1,
		half: Vector::new(
			fmaxss(x - 0.05, 0.0),
			fmaxss(y - 0.05, 0.0),
			fmaxss(z - 0.05, 0.0)
		),
		margin: 0.05
	};
}

fn copysign(a1: f32, a2: f32) -> f32 {
	if a1 <= 0.0 && a1.recip() < 0.0 {
		-a2
	} else {
		a2
	}
}

pub fn support(shape: Shape, a2: Vector) -> Vector {
	let half = shape.half;
	let geom = shape.geom;

	match geom {
		GeomShape::Ball => {
			return Vector::new(0.0, 0.0, 0.0);
		},

		GeomShape::Block => {
			return Vector::new(
				copysign(a2.x, half.x),
				copysign(a2.y, half.y),
				copysign(a2.z, half.z)
			);
		},

		GeomShape::Wedge => {
			let yy = a2.y * half.y;
			let zz = a2.z * half.z;
			if yy < 0.0 && zz > 0.0 {
				return Vector::new(
					copysign(a2.x, half.x),
					copysign(a2.y, half.y),
					copysign(a2.z, half.z)
				);
			}
			if yy > -zz {
				return Vector::new(
					copysign(a2.x, half.x), half.y, half.z
				);
			} else {
				return Vector::new(
					copysign(a2.x, half.x), -half.y, -half.z
				);
			}
		},

		GeomShape::Cylinder => {
			let x = copysign(a2.x, half.x);
			let yzlen = (a2.y * a2.y + a2.z * a2.z).sqrt();
			if yzlen > 0.0 {
				let y = a2.y / yzlen * half.y;
				let z = a2.z / yzlen * half.y;
				return Vector::new(x, y, z);
			} else {
				return Vector::new(x, 0.0, 0.0);
			}
		},

		GeomShape::CornerWedge => {
			let xx = a2.x * half.x;
			let yy = a2.y * half.y;
			let zz = a2.z * half.z;

			let mut v11 = false;
			if yy > 0.0 {
				let v12 = zz < -xx;
				if !v12 {
					if yy > zz && !v12 {
						let mut v13 = true;
						if xx <= 0.0 {
							v13 = true;
							if zz >= 0.0 {
								v13 = (zz - xx) < yy;
							}
						}
						v11 = v13;
					}
				} else if -xx < yy || yy > zz && !v12 {
					let mut v13 = true;
					if xx <= 0.0 {
						v13 = true;
						if zz >= 0.0 {
							v13 = (zz - xx) < yy;
						}
					}
					v11 = v13;
				}
			}

			if v11 {
				return Vector::new(half.x, half.y, -half.z);
			} else {
				return Vector::new(
					copysign(a2.x, half.x),
					-half.y,
					copysign(a2.z, half.z)
				);
			}
		}
	}
}

/// a1 type is a guess, cloning a2.p or a2 is probably not a good idea
pub fn collision_transform(geom: GeomShape, a2: &Xform) -> Xform {
	match geom {
		GeomShape::Ball => {
			return Xform {
				r: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
				p: a2.p.clone()
			};
		},

		GeomShape::Cylinder => {
			let v4 = Vector::new(a2.r[0], a2.r[3], a2.r[6]);
			let p = {
				if v4.y.abs() < 0.9 {
					Vector::new(1.0, 0.0, 0.0)
				} else {
					Vector::new(0.0, 1.0, 0.0)
				}
			};
			let v1 = sse_normalize(&cross(&v4, &p));
			let v5 = sse_normalize(&cross(&v1, &v4));
			return Xform {
				r: [
					v4.x, v5.x, v1.x,
					v4.y, v5.y, v1.y,
					v4.z, v5.z, v1.z
				],
				p: a2.p.clone()
			};
		},

		_ => {
			return a2.clone();
		}
	}
}

/// a1 type is a guess
fn extents_size(geom: GeomShape, a2: &Vector) -> Vector {
	let x = fmaxss(a2.x, 0.01);
	let y = fmaxss(a2.y, 0.01);
	let z = fmaxss(a2.z, 0.01);

	match geom {
		GeomShape::Ball => {
			let d = fminss(z, fminss(y, x));
			return Vector::new(d, d, d);
		},
		GeomShape::Cylinder => {
			let d = fminss(z, y);
			return Vector::new(x, d, d);
		},
		_ => {
			return Vector::new(x, y, z);
		}
	}
}

pub fn simulation_extents(geom: GeomShape, a2: &Vector, a3: &Xform)
	-> (Vector, Vector)
{
	let extents = extents_size(geom, a2);
	let hx = extents.x * 0.5;
	let hy = extents.y * 0.5;
	let hz = extents.z * 0.5;

	let r = a3.r;
	let p = a3.p;

	let mut min = [0.0; 3];
	let mut max = [0.0; 3];
	let pos = [p.x, p.y, p.z];
	let half = [hx, hy, hz];

	for i in 0..3 {
		let mut radius = r[i * 3].abs() * half[0];
		radius = radius + (r[i * 3 + 1].abs() * half[1]);
		radius = radius + (r[i * 3 + 2].abs() * half[2]);

		min[i] = pos[i] - radius - 0.01;
		max[i] = pos[i] + radius + 0.01;
	}

	return (
		Vector::new(min[0], min[1], min[2]),
		Vector::new(max[0], max[1], max[2])
	);
}

pub fn fuzzy_overlap(a1: &Vector, a2: &Vector, a3: &Vector, a4: &Vector)
	-> bool
{
	if a4.x >= a1.x && a2.x >= a3.x {
		if a4.y >= a1.y && a2.y >= a3.y {
			if a4.z >= a1.z && a2.z >= a3.z {
				return true;
			}
			return false;
		}
		return false;
	}
	return false;
}