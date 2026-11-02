use std::cell::RefCell;
use std::rc::Rc;
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
const U74: f32 = 0.0001;
const U83: f32 = 0.00001;
const U88: [usize; 3] = [1, 2, 0];
const U92: [usize; 3] = [1, 2, 0];
const U96: [usize; 3] = [2, 0, 1];
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

fn new_simplex() -> Simplex {
	return Simplex { rank: 0, c: Vec::new(), p: [0.0, 0.0, 0.0, 0.0] };
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
	simplex: usize
}

fn new_gjk(a1: MinkowskiDiff) -> GJK {
	return GJK {
		ray: Vector3::new(0.0, 0.0, 0.0),
		distance: 0.0,
		nfree: 0,
		current: 0,
		status: 0,
		shape: a1,
		simplices: [new_simplex(), new_simplex()],
		store: [
			SV {
				pa: Vector3::new(0.0, 0.0, 0.0),
				pb: Vector3::new(0.0, 0.0, 0.0),
				w: Vector3::new(0.0, 0.0, 0.0)
			},
			SV {
				pa: Vector3::new(0.0, 0.0, 0.0),
				pb: Vector3::new(0.0, 0.0, 0.0),
				w: Vector3::new(0.0, 0.0, 0.0)
			},
			SV {
				pa: Vector3::new(0.0, 0.0, 0.0),
				pb: Vector3::new(0.0, 0.0, 0.0),
				w: Vector3::new(0.0, 0.0, 0.0)
			},
			SV {
				pa: Vector3::new(0.0, 0.0, 0.0),
				pb: Vector3::new(0.0, 0.0, 0.0),
				w: Vector3::new(0.0, 0.0, 0.0)
			}
		],
		free: Vec::new(),
		simplex: 0
	};
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
			a1.simplex = current;
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
			a1.simplex = current;
			a1.distance = if a1.status == 0 { length(&a1.ray) } else { 0.0 };
			return a1.status;
		}

		v17_idx = (v17_idx + 1) % 4;
		v18[v17_idx] = w;

		let v9 = dot(&w, &a1.ray) / v1;
		v19 = fmaxss(v9, v19);

		let v6 = -v1 * U47;
		let v12 = v1 - v19;
		if v6 + v12 <= 0.0 {
			v16.rank -= 1;
			let popped = v16.c.pop().unwrap();
			a1.free.push(popped);
			a1.simplex = current;
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
			a1.simplex = current;
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
			a1.simplex = v15;
			return 2;
		}

		v13_count += 1;

		if a1.status != 0 {
			a1.simplex = a1.current;
			a1.distance = if a1.status == 0 { length(&a1.ray) } else { 0.0 };
			return a1.status;
		}
	}
}

pub fn enclose_origin(a1: &mut GJK) -> bool {
	let rank = a1.simplices[a1.simplex].rank;

	if rank == 1 {
		for j in 0..3 {
			let v3 = U113[j];
			
			// Branch 1: positive direction (v3)
			let v4 = a1.simplices[a1.simplex].rank;
			a1.simplices[a1.simplex].p[v4] = 0.0;
			let v5 = a1.free.pop().unwrap();
			a1.simplices[a1.simplex].c.push(v5);
			a1.simplices[a1.simplex].rank = v4 + 1;
			
			get_support(&a1.shape, &v3, &mut a1.simplices[a1.simplex].c[a1.simplices[a1.simplex].rank - 1]);
			
			if enclose_origin(a1) {
				return true;
			}
			
			let rank_2 = a1.simplices[a1.simplex].rank;
			a1.simplices[a1.simplex].rank = rank_2 - 1;
			let popped1 = a1.simplices[a1.simplex].c.pop().unwrap();
			a1.free.push(popped1);

			// Branch 2: negative direction (-v3)
			let v4_neg = vneg(&v3);
			let v5_idx = a1.simplices[a1.simplex].rank;
			a1.simplices[a1.simplex].p[v5_idx] = 0.0;
			let v6 = a1.free.pop().unwrap();
			a1.simplices[a1.simplex].c.push(v6);
			a1.simplices[a1.simplex].rank = v5_idx + 1;
			
			get_support(&a1.shape, &v4_neg, &mut a1.simplices[a1.simplex].c[a1.simplices[a1.simplex].rank - 1]);
			
			if enclose_origin(a1) {
				return true;
			}
			
			let rank_3 = a1.simplices[a1.simplex].rank;
			a1.simplices[a1.simplex].rank = rank_3 - 1;
			let popped2 = a1.simplices[a1.simplex].c.pop().unwrap();
			a1.free.push(popped2);
		}
		return false;
	}

	if rank == 2 {
		let v1 = a1.simplices[a1.simplex].c[1].w - a1.simplices[a1.simplex].c[0].w;
		for i in 0..3 {
			let v4 = cross(&v1, &U113[i]);
			let v5_len = length2(&v4);
			
			if v5_len > 0.0 {
				// Branch 1: positive direction (v4)
				let v5_idx = a1.simplices[a1.simplex].rank;
				a1.simplices[a1.simplex].p[v5_idx] = 0.0;
				let v6 = a1.free.pop().unwrap();
				a1.simplices[a1.simplex].c.push(v6);
				a1.simplices[a1.simplex].rank = v5_idx + 1;
				
				get_support(&a1.shape, &v4, &mut a1.simplices[a1.simplex].c[a1.simplices[a1.simplex].rank - 1]);
				
				if enclose_origin(a1) {
					return true;
				}
				
				let rank_4 = a1.simplices[a1.simplex].rank;
				a1.simplices[a1.simplex].rank = rank_4 - 1;
				let popped1 = a1.simplices[a1.simplex].c.pop().unwrap();
				a1.free.push(popped1);

				// Branch 2: negative direction (-v4)
				let v5_neg = vneg(&v4);
				let v6_idx = a1.simplices[a1.simplex].rank;
				a1.simplices[a1.simplex].p[v6_idx] = 0.0;
				let v7 = a1.free.pop().unwrap();
				a1.simplices[a1.simplex].c.push(v7);
				a1.simplices[a1.simplex].rank = v6_idx + 1;
				
				get_support(&a1.shape, &v5_neg, &mut a1.simplices[a1.simplex].c[a1.simplices[a1.simplex].rank - 1]);
				
				if enclose_origin(a1) {
					return true;
				}
				
				let rank_5 = a1.simplices[a1.simplex].rank;
				a1.simplices[a1.simplex].rank = rank_5 - 1;
				let popped2 = a1.simplices[a1.simplex].c.pop().unwrap();
				a1.free.push(popped2);
			}
		}
		return false;
	}

	if rank == 4 {
		let (v1, v2) = det_parts(
			&(a1.simplices[a1.simplex].c[0].w - a1.simplices[a1.simplex].c[3].w),
			&(a1.simplices[a1.simplex].c[1].w - a1.simplices[a1.simplex].c[3].w),
			&(a1.simplices[a1.simplex].c[2].w - a1.simplices[a1.simplex].c[3].w),
		);
		return v1 != v2;
	}

	if rank != 3 {
		return false;
	}

	let v1 = cross(
		&(a1.simplices[a1.simplex].c[1].w - a1.simplices[a1.simplex].c[0].w),
		&(a1.simplices[a1.simplex].c[2].w - a1.simplices[a1.simplex].c[0].w),
	);
	let v2_len = length2(&v1);

	if v2_len > 0.0 {
		// Branch 1: positive direction (v1)
		let v2_idx = a1.simplices[a1.simplex].rank;
		a1.simplices[a1.simplex].p[v2_idx] = 0.0;
		let v9 = a1.free.pop().unwrap();
		a1.simplices[a1.simplex].c.push(v9);
		a1.simplices[a1.simplex].rank = v2_idx + 1;
		
		get_support(&a1.shape, &v1, &mut a1.simplices[a1.simplex].c[a1.simplices[a1.simplex].rank - 1]);
		
		if enclose_origin(a1) {
			return true;
		}
		
		let rank_6 = a1.simplices[a1.simplex].rank;
		a1.simplices[a1.simplex].rank = rank_6 - 1;
		let popped1 = a1.simplices[a1.simplex].c.pop().unwrap();
		a1.free.push(popped1);

		// Branch 2: negative direction (-v1)
		let v2_neg = vneg(&v1);
		let v9_idx = a1.simplices[a1.simplex].rank;
		a1.simplices[a1.simplex].p[v9_idx] = 0.0;
		let v3 = a1.free.pop().unwrap();
		a1.simplices[a1.simplex].c.push(v3);
		a1.simplices[a1.simplex].rank = v9_idx + 1;
		
		get_support(&a1.shape, &v2_neg, &mut a1.simplices[a1.simplex].c[a1.simplices[a1.simplex].rank - 1]);
		
		if enclose_origin(a1) {
			return true;
		}
		
		let rank_7 = a1.simplices[a1.simplex].rank;
		a1.simplices[a1.simplex].rank = rank_7 - 1;
		let popped2 = a1.simplices[a1.simplex].c.pop().unwrap();
		a1.free.push(popped2);
	}

	return false;
}

#[derive(Clone)]
pub struct Face {
	pub n: Vector3,
	pub d: f32,
	pub pass: u8,
	pub c: [SV; 3],
	pub f: [Option<Rc<RefCell<Face>>>; 3],
	pub e: [usize; 3],
	pub l0: Option<Rc<RefCell<Face>>>,
	pub l1: Option<Rc<RefCell<Face>>>,
}

pub struct Hull {
	pub root: Option<Rc<RefCell<Face>>>,
	pub count: usize,
}

pub struct EPA {
	pub status: u8,
	pub normal: Vector3,
	pub depth: f32,
	pub nextsv: usize,
	pub result: Simplex,
	pub sv_store: Vec<SV>,
	pub hull: Hull,
	pub stock: Hull,
}

pub struct ExpandArg {
	pub cf: Option<Rc<RefCell<Face>>>,
	pub ff: Option<Rc<RefCell<Face>>>,
	pub nf: usize,
}

pub fn list_append(hull: &mut Hull, face: &Rc<RefCell<Face>>) {
	{
		let mut f = face.borrow_mut();
		f.l0 = None;
		f.l1 = hull.root.clone();
	}
	if let Some(root) = &hull.root {
		root.borrow_mut().l0 = Some(face.clone());
	}
	hull.root = Some(face.clone());
	hull.count += 1;
}

pub fn list_remove(hull: &mut Hull, face: &Rc<RefCell<Face>>) {
	let (l0, l1) = {
		let f = face.borrow();
		(f.l0.clone(), f.l1.clone())
	};

	if let Some(ref next_face) = l1 {
		next_face.borrow_mut().l0 = l0.clone();
	}
	if let Some(ref prev_face) = l0 {
		prev_face.borrow_mut().l1 = l1.clone();
	}

	if let Some(root) = &hull.root {
		if Rc::ptr_eq(root, face) {
			hull.root = l1.clone();
		}
	}
	hull.count -= 1;
}

pub fn bind(a1: &Rc<RefCell<Face>>, a2: usize, a3: &Rc<RefCell<Face>>, a4: usize) {
	a1.borrow_mut().e[a2] = a4;
	a1.borrow_mut().f[a2] = Some(a3.clone());
	a3.borrow_mut().e[a4] = a2;
	a3.borrow_mut().f[a4] = Some(a1.clone());
}

pub fn new_epa() -> EPA {
	let mut sv_store = Vec::with_capacity(64);
	for _ in 0..64 {
		sv_store.push(SV {
			pa: Vector3::new(0.0, 0.0, 0.0),
			pb: Vector3::new(0.0, 0.0, 0.0),
			w: Vector3::new(0.0, 0.0, 0.0),
		});
	}

	let mut stock = Hull { root: None, count: 0 };
	for _ in 0..128 {
		let face = Rc::new(RefCell::new(Face {
			n: Vector3::new(0.0, 0.0, 0.0),
			d: 0.0,
			pass: 0,
			c: [
				SV { pa: Vector3::new(0.0, 0.0, 0.0), pb: Vector3::new(0.0, 0.0, 0.0), w: Vector3::new(0.0, 0.0, 0.0) },
				SV { pa: Vector3::new(0.0, 0.0, 0.0), pb: Vector3::new(0.0, 0.0, 0.0), w: Vector3::new(0.0, 0.0, 0.0) },
				SV { pa: Vector3::new(0.0, 0.0, 0.0), pb: Vector3::new(0.0, 0.0, 0.0), w: Vector3::new(0.0, 0.0, 0.0) },
			],
			f: [None, None, None],
			e: [0, 0, 0],
			l0: None,
			l1: None,
		}));
		list_append(&mut stock, &face);
	}

	return EPA {
		status: 9,
		normal: Vector3::new(0.0, 0.0, 0.0),
		depth: 0.0,
		nextsv: 0,
		result: Simplex {
			rank: 0,
			c: Vec::new(),
			p: [0.0; 4],
		},
		sv_store,
		hull: Hull { root: None, count: 0 },
		stock,
	};
}

pub fn find_best(a1: &EPA) -> Rc<RefCell<Face>> {
	let root = a1.hull.root.as_ref().expect("Hull root is None").clone();
	let mut best = root.clone();
	let (mut best_d2, mut l1_opt) = {
		let r = best.borrow();
		(r.d * r.d, r.l1.clone())
	};

	while let Some(l1) = l1_opt {
		let (d2, next_l1) = {
			let r = l1.borrow();
			(r.d * r.d, r.l1.clone())
		};
		if d2 < best_d2 {
			best_d2 = d2;
			best = l1.clone();
		}
		l1_opt = next_l1;
	}
	return best;
}

pub fn get_edge_dist(a1: &Face, a2: &SV, a3: &SV) -> (bool, f32) {
	let v1 = a3.w - a2.w;
	if dot(&a2.w, &cross(&v1, &a1.n)) < 0.0 {
		return (false, 0.0);
	}
	if dot(&a2.w, &v1) > 0.0 {
		return (true, length(&a2.w));
	}
	if dot(&a3.w, &v1) < 0.0 {
		return (true, length(&a3.w));
	}
	let v2 = length2(&v1);
	let v3 = dot(&a2.w, &a3.w);
	let v4 = length2(&a2.w) * length2(&a3.w);
	let v5 = v4 - v3 * v3;

	if v2 <= 0.0 {
		return (true, f32::NAN);
	}
	let v9 = v5 / v2;
	if v9.is_nan() || v9 < 0.0 {
		return (true, f32::NAN);
	}
	return (true, v9.sqrt());
}

pub fn epa_new_face(a1: &mut EPA, a2: &SV, a3: &SV, a4: &SV, a5: bool) -> Option<Rc<RefCell<Face>>> {
	let root = match a1.stock.root.take() {
		Some(r) => r,
		None => {
			a1.status = 5;
			return None;
		}
	};

	let stock_l1 = root.borrow().l1.clone();
	a1.stock.root = stock_l1.clone();
	if let Some(ref l1) = stock_l1 {
		l1.borrow_mut().l0 = None;
	}
	a1.stock.count -= 1;

	let hull_root = a1.hull.root.clone();
	{
		let mut r = root.borrow_mut();
		r.l0 = None;
		r.l1 = hull_root.clone();
	}
	if let Some(ref hr) = hull_root {
		hr.borrow_mut().l0 = Some(root.clone());
	}
	a1.hull.root = Some(root.clone());
	a1.hull.count += 1;

	{
		let mut r = root.borrow_mut();
		r.pass = 0;
		r.c[0] = a2.clone();
		r.c[1] = a3.clone();
		r.c[2] = a4.clone();
		r.n = cross(&(a3.w - a2.w), &(a4.w - a2.w));
	}

	let v3 = length2(&root.borrow().n);
	if v3.is_nan() || v3 < 0.0 {
		return None;
	}
	let v2 = v3.sqrt();

	let mut status_val = 2;
	if U74 < v2 {
		let (mut ok, mut dist) = get_edge_dist(&root.borrow(), a2, a3);
		if !ok {
			let res = get_edge_dist(&root.borrow(), a3, a4);
			ok = res.0;
			dist = res.1;
		}
		if !ok {
			let res = get_edge_dist(&root.borrow(), a4, a2);
			ok = res.0;
			dist = res.1;
		}
		{
			let mut r = root.borrow_mut();
			if !ok {
				r.d = dot(&a2.w, &r.n) / v2;
			} else {
				r.d = dist;
			}
			r.n = r.n * (1.0 / v2);
		}

		if a5 || !(root.borrow().d < -U83) {
			return Some(root);
		}
		status_val = 3;
	}

	a1.status = status_val;
	list_remove(&mut a1.hull, &root);
	list_append(&mut a1.stock, &root);
	return None;
}

pub fn expand(a1: &mut EPA, a2: u8, a3: &SV, a4: &Rc<RefCell<Face>>,
	a5: usize, a6: &mut ExpandArg) -> bool
{
	let (pass, n, d) = {
		let f = a4.borrow();
		(f.pass, f.n, f.d)
	};
	if pass == a2 {
		return false;
	}

	let v1 = U92[a5];
	let v2 = dot(&n, &a3.w) - d;

	if v2 < -U83 {
		let c_v1 = a4.borrow().c[v1].clone();
		let c_a5 = a4.borrow().c[a5].clone();
		let v4 = match epa_new_face(a1, &c_v1, &c_a5, a3, false) {
			Some(face) => face,
			None => return false,
		};

		{
			let mut f4 = v4.borrow_mut();
			f4.e[0] = a5;
			f4.f[0] = Some(a4.clone());
		}
		{
			let mut f_a4 = a4.borrow_mut();
			f_a4.e[a5] = 0;
			f_a4.f[a5] = Some(v4.clone());
		}

		if let Some(ref cf) = a6.cf {
			{
				let mut cf_mut = cf.borrow_mut();
				cf_mut.e[1] = 2;
				cf_mut.f[1] = Some(v4.clone());
			}
			{
				let mut v4_mut = v4.borrow_mut();
				v4_mut.e[2] = 1;
				v4_mut.f[2] = Some(cf.clone());
			}
		} else {
			a6.ff = Some(v4.clone());
		}
		a6.cf = Some(v4.clone());
		a6.nf += 1;
		return true;
	}

	let v5 = U96[a5];
	a4.borrow_mut().pass = a2;

	let f_v1_opt = a4.borrow().f[v1].clone();
	let e_v1 = a4.borrow().e[v1];
	let f_v5_opt = a4.borrow().f[v5].clone();
	let e_v5 = a4.borrow().e[v5];

	let res1 = if let Some(ref f) = f_v1_opt {
		expand(a1, a2, a3, f, e_v1, a6)
	} else {
		false
	};

	let res2 = if let Some(ref f) = f_v5_opt {
		expand(a1, a2, a3, f, e_v5, a6)
	} else {
		false
	};

	if res1 && res2 {
		list_remove(&mut a1.hull, a4);
		list_append(&mut a1.stock, a4);
		return true;
	}
	return false;
}

pub fn copy_face(a1: &Face) -> Face {
	return Face {
		n: a1.n,
		d: a1.d,
		pass: a1.pass,
		c: a1.c.clone(),
		f: [None, None, None],
		e: [0, 0, 0],
		l0: None,
		l1: None,
	};
}

pub fn epa_evaluate(a1: &mut EPA, a2: &mut GJK, a3: &Vector3) -> u8 {
	if a2.simplices[a2.simplex].rank > 1 && enclose_origin(a2) {
		while a1.hull.root.is_some() {
			let root = a1.hull.root.as_ref().unwrap().clone();
			list_remove(&mut a1.hull, &root);
			list_append(&mut a1.stock, &root);
		}
		a1.status = 0;
		a1.nextsv = 0;

		let (v4_det, v5_det) = det_parts(
			&(a2.simplices[a2.simplex].c[0].w - a2.simplices[a2.simplex].c[3].w),
			&(a2.simplices[a2.simplex].c[1].w - a2.simplices[a2.simplex].c[3].w),
			&(a2.simplices[a2.simplex].c[2].w - a2.simplices[a2.simplex].c[3].w),
		);

		let mut sim_c = a2.simplices[a2.simplex].c.clone();
		let mut sim_p = a2.simplices[a2.simplex].p;

		if v4_det < v5_det {
			sim_c.swap(0, 1);
			sim_p.swap(0, 1);
		}

		let f1 = epa_new_face(a1, &sim_c[0], &sim_c[1], &sim_c[2], true);
		let f2 = epa_new_face(a1, &sim_c[1], &sim_c[0], &sim_c[3], true);
		let f3 = epa_new_face(a1, &sim_c[2], &sim_c[1], &sim_c[3], true);
		let f4 = epa_new_face(a1, &sim_c[0], &sim_c[2], &sim_c[3], true);

		if a1.hull.count == 4 {
			let faces = [f1.unwrap(), f2.unwrap(), f3.unwrap(), f4.unwrap()];

			bind(&faces[0], 0, &faces[1], 0);
			bind(&faces[0], 1, &faces[2], 0);
			bind(&faces[0], 2, &faces[3], 0);
			bind(&faces[1], 1, &faces[3], 2);
			bind(&faces[1], 2, &faces[2], 2);
			bind(&faces[2], 2, &faces[3], 1);

			a1.status = 0;
			let mut v2_pass: u8 = 0;

			while a1.nextsv < 64 {
				a1.nextsv += 1;
				let nextsv_idx = a1.nextsv - 1;

				let mut v7 = find_best(a1);
				v2_pass = v2_pass.wrapping_add(1);
				v7.borrow_mut().pass = v2_pass;

				let mut support_v = a1.sv_store[nextsv_idx].clone();
				get_support(&a2.shape, &v7.borrow().n, &mut support_v);

				let v13 = dot(&v7.borrow().n, &support_v.w) - v7.borrow().d;
				if !(U74 < v13) {
					a1.status = 7;
				} else {
					let mut expand_arg = ExpandArg {
						cf: None,
						ff: None,
						nf: 0,
					};
					let mut ok_expand = true;
					for i in 0..3 {
						let f_i = v7.borrow().f[i].clone().unwrap();
						let e_i = v7.borrow().e[i];
						if !(ok_expand && expand(a1, v2_pass, &support_v, &f_i, e_i, &mut expand_arg)) {
							ok_expand = false;
							break;
						}
					}

					if !ok_expand || expand_arg.nf < 3 {
						a1.status = 4;
					} else {
						let cf = expand_arg.cf.clone().unwrap();
						let ff = expand_arg.ff.clone().unwrap();
						bind(&cf, 1, &ff, 2);

						list_remove(&mut a1.hull, &v7);
						list_append(&mut a1.stock, &v7);

						v7 = find_best(a1);
						if v2_pass != 255 {
							continue;
						}
					}
				}

				let (n, d) = {
					let r = v7.borrow();
					(r.n, r.d)
				};
				a1.normal = n;
				a1.depth = d;

				let result = &mut a1.result;
				result.rank = 3;
				result.c = vec![v7.borrow().c[0].clone(), v7.borrow().c[1].clone(), v7.borrow().c[2].clone()];

				let v_origin = n * d;
				let len3 = length(&cross(&(result.c[1].w - v_origin), &(result.c[2].w - v_origin)));
				let len18 = length(&cross(&(result.c[2].w - v_origin), &(result.c[0].w - v_origin)));
				let len11 = length(&cross(&(result.c[0].w - v_origin), &(result.c[1].w - v_origin)));
				let total_len = len3 + len18 + len11;

				result.p[0] = len3 / total_len;
				result.p[1] = len18 / total_len;
				result.p[2] = len11 / total_len;

				return a1.status;
			}

			a1.status = 6;
			let v7 = find_best(a1);
			let (n, d) = {
				let r = v7.borrow();
				(r.n, r.d)
			};
			a1.normal = n;
			a1.depth = d;

			let result = &mut a1.result;
			result.rank = 3;
			result.c = vec![v7.borrow().c[0].clone(), v7.borrow().c[1].clone(), v7.borrow().c[2].clone()];

			let v_origin = n * d;
			let len3 = length(&cross(&(result.c[1].w - v_origin), &(result.c[2].w - v_origin)));
			let len18 = length(&cross(&(result.c[2].w - v_origin), &(result.c[0].w - v_origin)));
			let len11 = length(&cross(&(result.c[0].w - v_origin), &(result.c[1].w - v_origin)));
			let total_len = len3 + len18 + len11;

			result.p[0] = len3 / total_len;
			result.p[1] = len18 / total_len;
			result.p[2] = len11 / total_len;

			return a1.status;
		}
	}

	a1.status = 8;
	let v5_len = length2(a3);
	let v4_val = if v5_len.is_nan() || v5_len < 0.0 {
		f32::NAN
	} else {
		v5_len.sqrt()
	};

	let normal = if !(v4_val > 0.0) {
		Vector3::new(1.0, 0.0, 0.0)
	} else {
		vneg(a3) * (1.0 / v4_val)
	};

	a1.normal = normal;
	a1.depth = 0.0;
	a1.result.rank = 1;
	if !a1.result.c.is_empty() {
		a1.result.c[0] = a2.simplices[a2.simplex].c[0].clone();
	} else {
		a1.result.c.push(a2.simplices[a2.simplex].c[0].clone());
	}
	a1.result.p[0] = 1.0;

	return a1.status;
}

pub struct Results {
	pub status: u8,
	pub normal: Vector3,
	pub distance: f32,
	pub witnesses: [Vector3; 2]
}

fn separated(a1: &GJK, a2: &Transform, a3: &mut Results) -> bool {
	let mut v2 = Vector3::new(0.0, 0.0, 0.0);
	let mut v3 = Vector3::new(0.0, 0.0, 0.0);
	let simplex = &a1.simplices[a1.simplex];
	for i in 0..simplex.rank {
		v2 = v2 + simplex.c[i].pa * simplex.p[i];
		v3 = v3 + simplex.c[i].pb * simplex.p[i];
	}
	a3.witnesses[0] = xform(a2, &v2);
	a3.witnesses[1] = xform(a2, &v3);

	let v4 = v2 - v3;
	let dist = length(&v4);
	a3.distance = dist;
	let factor = if !(U56 < dist) { 1.0 } else { 1.0 / dist };
	a3.normal = v4 * factor;
	a3.status = 0;
	return false;
}

pub fn new_results() -> Results {
	return Results {
		status: 0,
		normal: Vector3::new(0.0, 0.0, 0.0),
		distance: 0.0,
		witnesses: [
			Vector3::new(0.0, 0.0, 0.0),
			Vector3::new(0.0, 0.0, 0.0)
		]
	};
}

// unknown external struct
pub struct PenArgBc {
	pub bary: [f32; 4]
}
// unknown external struct
pub struct PenArg {
	pub n: u8,
	pub p: [Vector3; 4],
	pub q: [Vector3; 4],
	pub bc: PenArgBc
}

pub fn penetration_or_distance(a1: Shape, a2: &Transform, a3: Shape,
	a4: &Transform, a5: &Vector3, a6: &mut Results,
	a7_opt: Option<&PenArg>) -> bool
{
	a6.status = 0;
	a6.witnesses[0] = Vector3::new(0.0, 0.0, 0.0);
	a6.witnesses[1] = Vector3::new(0.0, 0.0, 0.0);

	let mut v2 = new_gjk(new_minkowski_diff(a1, a2, a3, a4));
	let mut v3 = false;

	if let Some(a7) = a7_opt {
		if a7.n == 4 {
			let v1_idx = v2.simplex;
			let v1 = &mut v2.simplices[v1_idx];
			v1.rank = 4;
			v1.c.clear();
			for i in 0..4 {
				let mut v4 = v2.store[i].clone();
				v4.pa = mat_apply_t(&a2.r, &(a7.p[i] - a2.p));
				v4.pb = mat_apply_t(&a2.r, &(a7.q[i] - a2.p));
				v4.w = v4.pa - v4.pb;
				v1.c.push(v4);
				v1.p[i] = a7.bc.bary[i];
			}
			if !enclose_origin(&mut v2) {
				v2.ray = Vector3::new(0.0, 0.0, 0.0);
				v2.nfree = 0;
				v2.current = 1;
				v2.status = 2;
				v2.distance = 0.0;
			} else {
				v3 = true;
			}
		}
	}

	if !v3 {
		let v1_status = gjk_evaluate(&mut v2, &vneg(a5));
		if v1_status == 2 {
			a6.status = 2;
			return false;
		}
		if v1_status == 0 {
			return separated(&mut v2, a2, a6);
		}
	}

	let mut v1_epa = new_epa();
	if epa_evaluate(&mut v1_epa, &mut v2, &vneg(a5)) == 9 {
		a6.status = 3;
		if v3 && gjk_evaluate(&mut v2, &vneg(a5)) != 0 {
			return false;
		}
		return separated(&mut v2, a2, a6);
	}

	let mut v5 = Vector3::new(0.0, 0.0, 0.0);
	let rank = v1_epa.result.rank;
	for j in 0..rank {
		v5 = v5 + v1_epa.result.c[j].pa * v1_epa.result.p[j];
	}

	a6.status = 1;
	a6.witnesses[0] = xform(a2, &v5);
	a6.witnesses[1] = xform(a2, &(v5 - v1_epa.normal * v1_epa.depth));
	a6.normal = vneg(&v1_epa.normal);
	a6.distance = -v1_epa.depth;
	return true;
}

pub fn calc_pen_depth(a1_opt: Option<&PenArg>, a2: Shape, a3: Shape,
	a4: &Transform, a5: &Transform
) -> (bool, Option<Vector3>, Option<Vector3>, Option<Vector3>) {
	let v1 = a5.p - a4.p;
	let mut v2 = Results {
		status: 0,
		normal: Vector3::new(0.0, 0.0, 0.0),
		distance: 0.0,
		witnesses: [Vector3::new(0.0, 0.0, 0.0), Vector3::new(0.0, 0.0, 0.0)],
	};
	
	let v3 = penetration_or_distance(a2, a4, a3, a5, &v1, &mut v2, a1_opt);
	
	if v2.status <= 1 {
		(v3, Some(v2.normal), Some(v2.witnesses[0]), Some(v2.witnesses[1]))
	} else {
		(v3, None, None, None)
	}
}