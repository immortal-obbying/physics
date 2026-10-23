use crate::util::transform::Transform;
use crate::util::vector::Vector3;
use crate::engine::bullet_math::{
	fmaxss, dot, length2, length, cross, vneg, mat_apply, mat_apply_t,
	transpose_times, xform, rsqrtss
};
use crate::engine::bullet_shape::{
	Shape,
	support
};
use crate::engine::math::{
	Matrix, Quaternion,
	transpose
};


pub fn det(a1: &Vector3, a2: &Vector3, a3: &Vector3) -> f32 {
	return (a1.y * a2.z * a3.x)
		+ (a1.z * a2.x * a3.y)
		- (a1.x * a2.z * a3.y)
		- (a1.y * a2.x * a3.z)
		+ (a1.x * a2.y * a3.z)
		- (a1.z * a2.y * a3.x);
}

pub fn det_parts(a1: &Vector3, a2: &Vector3, a3: &Vector3) -> (f32, f32) {
	return ((a1.y * a2.z * a3.x)
		+ (a1.z * a2.x * a3.y)
		- (a1.x * a2.z * a3.y)
		- (a1.y * a2.x * a3.z)
		+ (a1.x * a2.y * a3.z),
		(a1.z * a2.y * a3.x));
}

pub struct MinkowskiDiff {
	pub shape0: Shape,
	pub shape1: Shape,
	pub toshape1: Matrix,
	pub toshape0: Transform
}
pub fn new_minkowski_diff(a1: Shape, a2: &Transform, a3: Shape,
	a4: &Transform) -> MinkowskiDiff
{
	let toshape1 = transpose_times(&a4.r, &a2.r);
	let toshape0 = Transform {
		r: transpose(&toshape1),
		p: mat_apply_t(&a2.r, &(a4.p - a2.p))
	};
	return MinkowskiDiff {
		shape0: a1,
		shape1: a3,
		toshape1: toshape1,
		toshape0: toshape0
	};
}

/// could possibly make this public
struct SV {
	pa: Vector3,
	pb: Vector3,
	w: Vector3
}

fn get_support(a1: &MinkowskiDiff, a2: &Vector3) -> SV {
	let v1: Vector3 = a2 * rsqrtss(length2(a2));
	let pa = support(&a1.shape0, v1);
	let pb = mat_apply(&a1.toshape0.r, &support(
		&a1.shape1,
		mat_apply(&a1.toshape1, &vneg(&v1))
	));
	return SV {
		pa: pa,
		pb: pb,
		w: pa - pb
	};
}

pub fn project_origin_2(a1: &Vector3, a2: &Vector3, a3: &mut [f32;2])
	-> (f32, usize)
{
	let v1 = a2 - a1;
	let len = length2(&v1);
	if !(len > 0.0) {
		return (-1.0, 0);
	}

	let v5 = -(dot(a1, &v1)) / len;
	if v5 >= 1.0 {
		a3[0] = 0.0;
		a3[1] = 1.0;
		return (length2(a2), 2);
	} else if v5 <= 0.0 {
		a3[0] = 1.0;
		a3[1] = 0.0;
		return (length2(a1), 1);
	} else {
		a3[1] = v5;
		a3[0] = 1.0 - v5;
		return (length2(&(v1 * v5 + a1)), 3);
	}
}

static U88: [usize; 3] = [2, 3, 1];
pub fn project_origin_3(a1: &Vector3, a2: &Vector3, a3: &Vector3,
	a4: &mut [f32;3]) -> (f32, usize)
{
	let v4 = [a1, a2, a3];
	let v5 = [a1 - a2, a2 - a3, a3 - a1];
	let v6 = cross(&v5[0], &v5[1]);
	let v7 = length2(&v6);

	if !(v7 > 0.0) {
		return (-1.0, 0);
	}

	let mut v8: f32 = -1.0;
	let mut v9: usize = 0;
	let mut v10 = [0.0f32; 2];

	for i in 0..3 {
		let cross_v5_v6 = cross(&v5[i], &v6);
		let v1_val = dot(v4[i], &cross_v5_v6);

		if v1_val > 0.0 {
			let u88_idx = U88[i];
			let target_idx = u88_idx - 1;

			let (v2, v3) = project_origin_2(v4[i], v4[target_idx], &mut v10);

			if v8 < 0.0 || v2 < v8 {
				v8 = v2;
				let bit1 = if v3 % 2 != 1 { 0 } else { 1 << i };
				let bit2 = if v3 < 2 { 0 } else { 1 << target_idx };
				v9 = bit1 + bit2;
				a4[i] = v10[0];
				a4[target_idx] = v10[1];
				a4[U88[target_idx] - 1] = 0.0;
			}
		}
	}

	if v8 < 0.0 {
		let v10_dot = dot(a1, &v6);
		
		let v11 = if v7.is_nan() {
			f32::NAN
		} else if v7 >= 0.0 {
			v7.sqrt()
		} else {
			f32::NAN
		};

		let v12 = v10_dot / v7;
		let v14 = v6 * v12;
		v8 = length2(&v14);
		v9 = 7;

		let v15_1 = length(&cross(&v5[1], &(a2 - &v14))) / v11;
		a4[0] = v15_1;
		let v15_2 = length(&cross(&v5[2], &(a3 - &v14))) / v11;
		a4[1] = v15_2;
		let v16 = a4[1] + a4[0];
		a4[2] = 1.0 - v16;
	}

	return (v8, v9);
}