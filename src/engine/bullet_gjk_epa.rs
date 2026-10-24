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


const U47: f32 = 0.0001;
const U56: f32 = 0.0001;
const U65: f32 = 0.0001;
const U88: [usize; 3] = [1, 2, 0];
const U113: [Vector3; 3] = [
	Vector3 { x: 1.0, y: 0.0, z: 0.0 },
	Vector3 { x: 0.0, y: 1.0, z: 0.0 },
	Vector3 { x: 0.0, y: 0.0, z: 1.0 },
];

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
#[derive(Clone)]
struct SV {
	pa: Vector3,
	pb: Vector3,
	w: Vector3
}

fn get_support(a1: &MinkowskiDiff, a2: &Vector3, a3: &mut SV) {
	let v1: Vector3 = a2 * rsqrtss(length2(a2));
	let pa = support(&a1.shape0, v1);
	let pb = mat_apply(&a1.toshape0.r, &support(
		&a1.shape1,
		mat_apply(&a1.toshape1, &vneg(&v1))
	));
	a3.pa = pa;
	a3.pb = pb;
	a3.w = pa - pb;
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
	let mut v10: [f32;2] = [0.0; 2];

	for i in 0..3 {
		let v1_val = dot(v4[i], &cross(&v5[i], &v6));
		if v1_val > 0.0 {
			let target_idx = U88[i];
			let (v2, v3) = project_origin_2(v4[i], v4[target_idx], &mut v10);
			if v8 < 0.0 || v2 < v8 {
				v8 = v2;
				let bit1 = if v3 % 2 != 1 { 0 } else { 1 << i };
				let bit2 = if v3 < 2 { 0 } else { 1 << target_idx };
				v9 = bit1 + bit2;
				a4[i] = v10[0];
				a4[target_idx] = v10[1];
				a4[U88[target_idx]] = 0.0;
			}
		}
	}

	if v8 < 0.0 {
		let v10_dot = dot(a1, &v6);
		let v11 = v7.sqrt();
		let v14 = v6 * (v10_dot / v7);
		v8 = length2(&v14);
		v9 = 7;
		a4[0] = length(&cross(&v5[1], &(a2 - &v14))) / v11;
		a4[1] = length(&cross(&v5[2], &(a3 - &v14))) / v11;
		a4[2] = 1.0 - (a4[1] + a4[0]);
	}
	return (v8, v9);
}

pub fn project_origin_4(a1: &Vector3, a2: &Vector3, a3: &Vector3,
	a4: &Vector3, a5: &mut [f32;4]) -> (f32, usize)
{
	let v6 = [a1, a2, a3, a4];
	let v7 = [a1 - a4, a2 - a4, a3 - a4];
	let v8 = det(&v7[0], &v7[1], &v7[2]);
	if v8 == 0.0 {
		return (-1.0, 0);
	}
	let v9 = dot(&a1, &cross(&(a2 - a3), &(a1 - a2))) * v8;
	if 0.0 < v9 {
		return (-1.0, 0);
	}

	let mut v11: f32 = -1.0;
	let mut v12: usize = 0;
	let mut v13: [f32;3] = [0.0; 3];
	for i in 0..3 {
		let target_idx = U88[i];
		let v3 = dot(a4, &cross(&v7[i], &v7[target_idx])) * v8;
		if 0.0 < v3 {
			let (v1, v2) = project_origin_3(&v6[i], &v6[target_idx], a4, &mut v13);
			if v11 < 0.0 || v1 < v11 {
				let bit1 = if v2 % 2 != 1 { 0 } else { 1 << i };
				let bit2 = if (v2 >> 1) % 2 != 1 { 0 } else { 1 << target_idx };
				let v3 = bit1 + bit2;
				v12 = v3 + (if v2 < 4 { 0 } else { 8 });
				a5[i] = v13[0];
				a5[target_idx] = v13[1];
				a5[U88[target_idx]] = 0.0;
				a5[3] = v13[2];
			}
		}
	}

	if v11 < 0.0 {
		v12 = 15;
		a5[0] = det(a3, a2, a4) / v8;
		a5[1] = det(a1, a3, a4) / v8;
		a5[2] = det(a2, a1, a4) / v8;
		a5[3] = 1.0 - (a5[0] + a5[1] + a5[2]);
		v11 = 0.0;
	}
	return (v11, v12);
}

#[derive(Clone)]
struct Simplex {
	rank: usize, // guess
	c: Vec<SV>,
	p: [f32; 4]
}

struct GJK {
	ray: Vector3,
	distance: f32,
	current: usize,
	status: u8, // TODO: replace with enum once enough info is gathered
	shape: MinkowskiDiff,
	simplices: [Simplex;2],
	store: [SV;4],
	free: Vec<SV>,
	nfree: usize,
	simplex: Simplex // TODO: express in a better way (index?)
}

fn append_vertice(a1: &mut GJK, a2: &mut Simplex, a3: &Vector3) {
	let rank = a2.rank + 1;
	a2.p[rank] = 0.0;
	a1.nfree -= 1;
	a2.c[rank] = a1.free[a1.nfree].clone();
	a2.rank = rank;
	get_support(&a1.shape, a3, &mut a1.free[a1.nfree]);
}

fn remove_vertice(a1: &mut GJK, a2: &mut Simplex) {
	let rank = a2.rank;
	a2.rank = rank - 1;
	a1.nfree = a1.nfree + 1;
	a1.free[a1.nfree] = a2.c[rank].clone();
}

pub fn gjk_evaluate(a1: &mut GJK, a2: &Vector3) -> u8 {
	a1.nfree = 4;
	a1.free = a1.store.to_vec();
	a1.status = 0;
	a1.distance = 0.0;
	a1.current = 0;

	let v17 = &mut a1.simplices[a1.current];
	v17.rank = 0;

	let v18_len = length2(a2);
	a1.ray = *a2;

	let mut v19: f32 = 0.0;

	// Push initial support point
	let initial_sv = a1.free.pop().unwrap();
	v17.p[v17.rank] = 0.0;
	v17.c.push(initial_sv);
	v17.rank += 1;

	let initial_support_dir = if v18_len > 0.0 { vneg(a2) } else { Vector3::new(1.0, 0.0, 0.0) };
	get_support(&a1.shape, &initial_support_dir, &mut v17.c[0]);

	v17.p[0] = 1.0;
	a1.ray = v17.c[0].w;

	let mut v18 = [a1.ray, a1.ray, a1.ray, a1.ray];
	let mut v17_idx: usize = 0;
	let mut v13_count: usize = 0;

	loop {
		let current = a1.current;
		let v15 = 1 - current;

		// Split simplices borrow safely
		let (v16, v14) = if current == 0 {
			let (s0, s1) = a1.simplices.split_at_mut(1);
			(&mut s0[0], &mut s1[0])
		} else {
			let (s0, s1) = a1.simplices.split_at_mut(1);
			(&mut s1[0], &mut s0[0])
		};

		let v2_len = length2(&a1.ray);
		let v1 = if v2_len.is_nan() {
			f32::NAN
		} else if v2_len >= 0.0 {
			v2_len.sqrt()
		} else {
			f32::NAN
		};

		if v1 < U56 {
			a1.status = 1;
			a1.simplex = v16.clone();
			a1.distance = if a1.status == 0 { length(&a1.ray) } else { 0.0 };
			return a1.status;
		}

		let v2_neg = vneg(&a1.ray);
		let w_sv = a1.free.pop().unwrap_or_else(|| SV {
			pa: Vector3::new(0.0, 0.0, 0.0),
			pb: Vector3::new(0.0, 0.0, 0.0),
			w: Vector3::new(0.0, 0.0, 0.0),
		});

		let v3_rank_idx = v16.rank;
		v16.p[v3_rank_idx] = 0.0;
		v16.c.push(w_sv);
		v16.rank = v3_rank_idx + 1;

		get_support(&a1.shape, &v2_neg, &mut v16.c[v16.rank - 1]);
		let w = v16.c[v16.rank - 1].w;

		let mut stagnation = false;
		for i in 0..4 {
			if length2(&(w - v18[i])) < U65 {
				stagnation = true;
				break;
			}
		}

		if stagnation {
			v16.rank -= 1;
			let popped = v16.c.pop().unwrap();
			a1.free.push(popped);
			a1.simplex = v16.clone();
			a1.distance = if a1.status == 0 { length(&a1.ray) } else { 0.0 };
			return a1.status;
		}

		v17_idx = (v17_idx + 1) % 4;
		v18[v17_idx] = w;

		let v9 = dot(&w, &a1.ray) / v1;
		v19 = f32::max(v9, v19);

		let v6 = -v1 * U47;
		let v12 = v1 - v19;
		if v6 + v12 <= 0.0 {
			v16.rank -= 1;
			let popped = v16.c.pop().unwrap();
			a1.free.push(popped);
			a1.simplex = v16.clone();
			a1.distance = if a1.status == 0 { length(&a1.ray) } else { 0.0 };
			return a1.status;
		}

		let mut v2_weights = [0.0f32; 4];
		let mut v20_dist: f32 = 0.0;
		let mut v1_mask: usize = 0;

		match v16.rank {
			2 => {
				let (d, mask) = project_origin_2(&v16.c[0].w, &v16.c[1].w, &mut [v2_weights[0], v2_weights[1]]);
				v20_dist = d;
				v1_mask = mask;
			}
			3 => {
				let (d, mask) = project_origin_3(&v16.c[0].w, &v16.c[1].w, &v16.c[2].w, &mut [v2_weights[0], v2_weights[1], v2_weights[2]]);
				v20_dist = d;
				v1_mask = mask;
			}
			4 => {
				let (d, mask) = project_origin_4(&v16.c[0].w, &v16.c[1].w, &v16.c[2].w, &v16.c[3].w, &mut v2_weights);
				v20_dist = d;
				v1_mask = mask;
			}
			_ => {}
		}

		if v20_dist < 0.0 {
			v16.rank -= 1;
			let popped = v16.c.pop().unwrap();
			a1.free.push(popped);
			a1.simplex = v16.clone();
			a1.distance = if a1.status == 0 { length(&a1.ray) } else { 0.0 };
			return a1.status;
		}

		v14.rank = 0;
		v14.c.clear();
		a1.ray = Vector3::new(0.0, 0.0, 0.0);
		a1.current = v15;

		let rank_4 = v16.rank;
		for j in 0..rank_4 {
			if (v1_mask & (1 << j)) == 0 {
				a1.free.push(v16.c[j].clone());
			} else {
				let new_r = v14.rank;
				v14.c.push(v16.c[j].clone());
				v14.p[new_r] = v2_weights[j];
				v14.rank = new_r + 1;
				a1.ray = a1.ray + v16.c[j].w * v2_weights[j];
			}
		}

		if v1_mask == 15 {
			a1.status = 1;
		}

		if v13_count == 127 {
			a1.status = 2;
			a1.simplex = v14.clone();
			return 2;
		}

		v13_count += 1;

		if a1.status != 0 {
			a1.simplex = a1.simplices[a1.current].clone();
			a1.distance = if a1.status == 0 { length(&a1.ray) } else { 0.0 };
			return a1.status;
		}
	}
}

pub fn enclose_origin(a1: &mut GJK) -> bool {
	let rank = a1.simplex.rank;

	if rank == 1 {
		for j in 0..3 {
			let v3 = U113[j];
			
			// Branch 1: positive direction (v3)
			let v4 = a1.simplex.rank;
			a1.simplex.p[v4] = 0.0;
			let v5 = a1.free.pop().unwrap();
			a1.simplex.c.push(v5);
			a1.simplex.rank = v4 + 1;
			
			get_support(&a1.shape, &v3, &mut a1.simplex.c[a1.simplex.rank - 1]);
			
			if enclose_origin(a1) {
				return true;
			}
			
			let rank_2 = a1.simplex.rank;
			a1.simplex.rank = rank_2 - 1;
			let popped1 = a1.simplex.c.pop().unwrap();
			a1.free.push(popped1);

			// Branch 2: negative direction (-v3)
			let v4_neg = vneg(&v3);
			let v5_idx = a1.simplex.rank;
			a1.simplex.p[v5_idx] = 0.0;
			let v6 = a1.free.pop().unwrap();
			a1.simplex.c.push(v6);
			a1.simplex.rank = v5_idx + 1;
			
			get_support(&a1.shape, &v4_neg, &mut a1.simplex.c[a1.simplex.rank - 1]);
			
			if enclose_origin(a1) {
				return true;
			}
			
			let rank_3 = a1.simplex.rank;
			a1.simplex.rank = rank_3 - 1;
			let popped2 = a1.simplex.c.pop().unwrap();
			a1.free.push(popped2);
		}
		return false;
	}

	if rank == 2 {
		let v1 = a1.simplex.c[1].w - a1.simplex.c[0].w;
		for i in 0..3 {
			let v4 = cross(&v1, &U113[i]);
			let v5_len = length2(&v4);
			
			if v5_len > 0.0 {
				// Branch 1: positive direction (v4)
				let v5_idx = a1.simplex.rank;
				a1.simplex.p[v5_idx] = 0.0;
				let v6 = a1.free.pop().unwrap();
				a1.simplex.c.push(v6);
				a1.simplex.rank = v5_idx + 1;
				
				get_support(&a1.shape, &v4, &mut a1.simplex.c[a1.simplex.rank - 1]);
				
				if enclose_origin(a1) {
					return true;
				}
				
				let rank_4 = a1.simplex.rank;
				a1.simplex.rank = rank_4 - 1;
				let popped1 = a1.simplex.c.pop().unwrap();
				a1.free.push(popped1);

				// Branch 2: negative direction (-v4)
				let v5_neg = vneg(&v4);
				let v6_idx = a1.simplex.rank;
				a1.simplex.p[v6_idx] = 0.0;
				let v7 = a1.free.pop().unwrap();
				a1.simplex.c.push(v7);
				a1.simplex.rank = v6_idx + 1;
				
				get_support(&a1.shape, &v5_neg, &mut a1.simplex.c[a1.simplex.rank - 1]);
				
				if enclose_origin(a1) {
					return true;
				}
				
				let rank_5 = a1.simplex.rank;
				a1.simplex.rank = rank_5 - 1;
				let popped2 = a1.simplex.c.pop().unwrap();
				a1.free.push(popped2);
			}
		}
		return false;
	}

	if rank == 4 {
		let (v1, v2) = det_parts(
			&(a1.simplex.c[0].w - a1.simplex.c[3].w),
			&(a1.simplex.c[1].w - a1.simplex.c[3].w),
			&(a1.simplex.c[2].w - a1.simplex.c[3].w),
		);
		return v1 != v2;
	}

	if rank != 3 {
		return false;
	}

	let v1 = cross(
		&(a1.simplex.c[1].w - a1.simplex.c[0].w),
		&(a1.simplex.c[2].w - a1.simplex.c[0].w),
	);
	let v2_len = length2(&v1);

	if v2_len > 0.0 {
		// Branch 1: positive direction (v1)
		let v2_idx = a1.simplex.rank;
		a1.simplex.p[v2_idx] = 0.0;
		let v9 = a1.free.pop().unwrap();
		a1.simplex.c.push(v9);
		a1.simplex.rank = v2_idx + 1;
		
		get_support(&a1.shape, &v1, &mut a1.simplex.c[a1.simplex.rank - 1]);
		
		if enclose_origin(a1) {
			return true;
		}
		
		let rank_6 = a1.simplex.rank;
		a1.simplex.rank = rank_6 - 1;
		let popped1 = a1.simplex.c.pop().unwrap();
		a1.free.push(popped1);

		// Branch 2: negative direction (-v1)
		let v2_neg = vneg(&v1);
		let v9_idx = a1.simplex.rank;
		a1.simplex.p[v9_idx] = 0.0;
		let v3 = a1.free.pop().unwrap();
		a1.simplex.c.push(v3);
		a1.simplex.rank = v9_idx + 1;
		
		get_support(&a1.shape, &v2_neg, &mut a1.simplex.c[a1.simplex.rank - 1]);
		
		if enclose_origin(a1) {
			return true;
		}
		
		let rank_7 = a1.simplex.rank;
		a1.simplex.rank = rank_7 - 1;
		let popped2 = a1.simplex.c.pop().unwrap();
		a1.free.push(popped2);
	}

	return false;
}