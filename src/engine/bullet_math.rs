use crate::util::transform::Transform;
use crate::util::vector::Vector3;
use crate::engine::math::{ Matrix, Quaternion };

#[inline]
pub fn fsqrt(a1: f32) -> f32 {
	a1.sqrt()
}

#[inline]
pub fn fmaxss(a1: f32, a2: f32) -> f32 {
	if a2 < a1 { a1 } else { a2 }
}

#[inline]
pub fn fminss(a1: f32, a2: f32) -> f32 {
	if a1 < a2 { a1 } else { a2 }
}

#[inline]
pub fn fneg(a1: f32) -> f32 {
	-a1
}

#[inline]
pub fn dot(a1: &Vector3, a2: &Vector3) -> f32 {
	let v = a1 * a2;
	return v.x + v.y + v.z;
}

#[inline]
pub fn length2(a1: &Vector3) -> f32 {
	let v = a1 * a1;
	return v.x + v.y + v.z;
}

#[inline]
pub fn length(a1: &Vector3) -> f32 {
	let v = a1 * a1;
	return (v.x + v.y + v.z).sqrt();
}

#[inline]
pub fn cross(a1: &Vector3, a2: &Vector3) -> Vector3 {
	// to be concise these are all f32
	let x = (a1.y * a2.z) - (a1.z * a2.y);
	let y = (a1.z * a2.x) - (a1.x * a2.z);
	let z = (a1.x * a2.y) - (a1.y * a2.x);
	return Vector3 { x, y, z };
}

#[inline]
pub fn vneg(a1: &Vector3) -> Vector3 {
	Vector3 {
		x: -a1.x,
		y: -a1.y,
		z: -a1.z
	}
}

#[inline]
pub fn rsqrtss(a1: f32) -> f32 {
	// a1.sqrt().recip()
	let v1 = {
		if a1.is_nan() {
			f32::NAN
		} else if a1 >= 0.0 {
			a1.sqrt()
		} else {
			f32::NAN
		}
	};
	return 1.0 / v1;
}

#[inline]
pub fn sse_normalize(a1: &Vector3) -> Vector3 {
	let sq = a1 * a1;
	let len_sq = sq.x + sq.y + sq.z;
	let len_rsq = rsqrtss(len_sq);
	let sse_mul = len_rsq * (1.5 - (len_sq * 0.5 * len_rsq));
	return a1 * sse_mul;
}

#[inline]
pub fn sse_normalize_len(a1: &Vector3) -> Vector3 {
	let sq = a1 * a1;
	let len_sq = sq.x + sq.y + sq.z;
	let len_rsq = rsqrtss(len_sq);
	let sse_mul = (len_rsq * (len_sq * -0.5 * len_rsq)) + 1.5;
	return a1 * sse_mul;
}

#[inline]
pub fn mat_apply(a1: &Matrix, a2: &Vector3) -> Vector3 {
	let v1 = Vector3 { x: a1[0], y: a1[3], z: a1[6] };
	let v2 = Vector3 { x: a1[1], y: a1[4], z: a1[7] };
	let v3 = Vector3 { x: a1[2], y: a1[5], z: a1[8] };
	return a2.x * v1 + a2.y * v2 + a2.z * v3;
}

#[inline]
pub fn mat_apply_t(a1: &Matrix, a2: &Vector3) -> Vector3 {
	let v1 = Vector3 { x: a1[0], y: a1[1], z: a1[2] };
	let v2 = Vector3 { x: a1[3], y: a1[4], z: a1[5] };
	let v3 = Vector3 { x: a1[6], y: a1[7], z: a1[8] };
	return a2.x * v1 + a2.y * v2 + a2.z * v3;
}

#[inline]
pub fn xform(a1: &Transform, a2: &Vector3) -> Vector3 {
	let v1 = Vector3 { x: a1.r[0], y: a1.r[1], z: a1.r[2] };
	let v2 = Vector3 { x: a1.r[3], y: a1.r[4], z: a1.r[5] };
	let v3 = Vector3 { x: a1.r[6], y: a1.r[7], z: a1.r[8] };
	return a2.x * v1 + a2.y * v2 + a2.z * v3 + a1.p;
}

#[inline]
pub fn inv_xform(a1: &Transform, a2: &Vector3) -> Vector3 {
	let inv = a2 - a1.p;
	let v1 = Vector3 { x: a1.r[0], y: a1.r[1], z: a1.r[2] };
	let v2 = Vector3 { x: a1.r[3], y: a1.r[4], z: a1.r[5] };
	let v3 = Vector3 { x: a1.r[6], y: a1.r[7], z: a1.r[8] };
	return inv.x * v1 + inv.y * v2 + inv.z * v3 + a1.p;
}

#[inline]
pub fn transpose_times(a1: &Matrix, a2: &Matrix) -> Matrix {
	let v1 = Vector3 { x: a2[0], y: a2[1], z: a2[2] };
	let v2 = Vector3 { x: a2[3], y: a2[4], z: a2[5] };
	let v3 = Vector3 { x: a2[6], y: a2[7], z: a2[8] };
	let m1 = a1[0] * v1 + a1[3] * v2 + a1[6] * v3;
	let m2 = a1[1] * v1 + a1[4] * v2 + a1[7] * v3;
	let m3 = a1[2] * v1 + a1[5] * v2 + a1[8] * v3;
	return [
		m1.x, m1.y, m1.z,
		m2.x, m2.y, m2.z,
		m3.x, m3.y, m3.z
	];
}

pub fn mat_get_rotation(m: &Matrix) -> Quaternion {
	// chatgpt helped me here i had no clue what was going on and
	// the decompiler probably had a stroke
	let trace = (m[0] + m[4]) + m[8];
	let mut q: Quaternion = [0.0; 4];

	if !(trace > 0.0) {
		let diag = [m[0], m[4], m[8]];

		let i = if !(m[0] < m[4]) {
			if !(m[0] < m[8]) { 0 } else { 2 }
		} else if !(m[4] < m[8]) {
			1
		} else {
			2
		};

		let j = (i + 1) % 3;
		let k = (i + 2) % 3;
		let at = |row: usize, col: usize| -> f32 { m[row * 3 + col] };
		let s = ((diag[i] - diag[j]) - diag[k]) + 1.0;

		q[i] = s;
		q[3] = at(k, j) - at(j, k);
		q[j] = at(j, i) + at(i, j);
		q[k] = at(k, i) + at(i, k);

		let inv = 0.5 / s.sqrt();
		q[0] *= inv;
		q[1] *= inv;
		q[2] *= inv;
		q[3] *= inv;
	} else {
		let s = trace + 1.0;
		q[3] = s;
		q[0] = m[7] - m[5];
		q[1] = m[2] - m[6];
		q[2] = m[3] - m[1];

		let inv = 0.5 / s.sqrt();
		q[0] *= inv;
		q[1] *= inv;
		q[2] *= inv;
		q[3] *= inv;
	}

	return q;
}

#[inline]
pub fn quat_mul(a1: &Quaternion, a2: &Quaternion) -> Quaternion {
	let v13 = a2[3] * a1[0] + a2[2] * a1[1];
	let v11 = a2[3] * a1[1] + a2[0] * a1[2];
	let v14 = a2[3] * a1[2] + a2[1] * a1[0];

	let v16 = -(a2[0] * a1[0] + a2[1] * a1[1]);
	let v20 = a1[3] * a2[0] - a2[1] * a1[2];
	let v10 = a1[3] * a2[1] - a2[2] * a1[0];
	let v6 = a1[3] * a2[2] - a2[0] * a1[1];
	let v5 = a1[3] * a2[3] - a2[2] * a1[2];

	return [ v20 + v13, v10 + v11, v6 + v14, v5 + v16 ];
}

#[inline]
pub fn dot4(a1: &Quaternion, a2: &Quaternion) -> f32 {
	return a1[0] * a2[0]
		+ a1[1] * a2[1]
		+ a1[2] * a2[2]
		+ a1[3] * a2[3];
}

#[inline]
pub fn acosf(a1: f32) -> f32 {
	a1.acos()
}

// pub static ZERO3: Vector3 = Vector3 { x: 0.0, y: 0.0, z: 0.0 };