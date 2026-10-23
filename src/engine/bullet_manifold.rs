use crate::util::transform::Transform;
use crate::util::vector::Vector3;
use crate::engine::bullet_math::{
	fminss, dot, length2, xform, inv_xform
};

#[derive(Clone)]
pub struct Point {
	pub local_a: Vector3,
	pub local_b: Vector3,
	pub world_b: Vector3,
	pub world_a: Vector3,
	pub normal_b: Vector3,
	user: Option<u8>,
	pub distance: f32
}
pub struct Manifold {
	pub points: Vec<Point>,
	pub threshold: f32
}

fn area(a1: &Vector3, a2: &Vector3) -> f32 {
	let v6 = a2.y * a1.x - a2.x * a1.y;
	let v3 = a2.z * a1.y - a2.y * a1.z;
	let v7 = a2.x * a1.z - a2.z * a1.x;
	return v7 * v7 + v3 * v3 + v6 * v6;
}

fn max_axis_4(a1: &[f32;4]) -> usize {
	let mut axis: usize = 0;
	let mut axis_val = -1e+18;
	for i in 0..4 {
		let val = a1[i].abs();
		if axis_val < val {
			axis_val = val;
			axis = i;
		}
	}
	return axis;
}

pub fn get_cache_entry(a1: &Manifold, a2: &Point) -> Option<usize> {
	let threshold = a1.threshold * a1.threshold;
	let mut point_i: Option<usize> = None;

	for i in 0..a1.points.len() {
		let p = &a1.points[i];
		let distance = fminss(
			length2(&(p.local_a - a2.local_a)),
			length2(&(p.local_b - a2.local_b))
		);
		if distance < threshold {
			point_i = Some(i);
		}
	}
	return point_i;
}

pub fn sort_cached_points(a1: &Manifold, a2: &Point) -> usize {
	assert_eq!(a1.points.len(), 4);
	let distance = a2.distance;
	let mut v2: isize = -1;
	
	for (i, pt) in a1.points.iter().enumerate() {
		if pt.distance < distance {
			v2 = i as isize;
		}
	}

	let p0 = a1.points[0].local_a;
	let p1 = a1.points[1].local_a;
	let p2 = a1.points[2].local_a;
	let p3 = a1.points[3].local_a;
	let p4 = a2.local_a;

	return max_axis_4(&[
		if v2 == 0 { 0.0 } else { area(&(p4 - p1), &(p3 - p2)) },
		if v2 == 1 { 0.0 } else { area(&(p4 - p0), &(p3 - p2)) },
		if v2 == 2 { 0.0 } else { area(&(p4 - p0), &(p3 - p1)) },
		if v2 == 3 { 0.0 } else { area(&(p4 - p0), &(p2 - p1)) }
	]);
}

pub fn add_manifold_point(a1: &mut Manifold, point: Point) -> usize {
	if a1.points.len() != 4 {
		let i = a1.points.len();
		a1.points.push(point);
		return i;
	}

	let i = sort_cached_points(a1, &point);
	// if i < 0 { i = 0; }
	a1.points[i] = point;
	return i;
}

fn replace_contact_point(a1: &mut Manifold, a2: usize, mut point: Point) {
	point.user = a1.points[a2].user;
	a1.points[a2] = point;
}

pub fn remove_contact_point(a1: &mut Manifold, a2: usize) {
	a1.points[a2].user = None;
	let last = a1.points.len() - 1;
	if last != a2 {
		a1.points[a2] = a1.points[last].clone();
	}
	a1.points.pop();
}

pub fn refresh_contact_points(a1: &mut Manifold, a2: &Transform,
	a3: &Transform, a4: f32, a5: f32)
{
	let len = a1.points.len();
	if len == 0 {
		return;
	}
	
	let mut v9: Vec<f32> = vec![0.0; len];
	for idx in 0..len {
		let i = len - idx - 1;
		v9[i] = a1.points[i].distance;
		a1.points[i].world_a = xform(&a2, &a1.points[i].local_a);
		a1.points[i].world_b = xform(&a3, &a1.points[i].local_b);
		a1.points[i].distance = dot(&(a1.points[i].world_a - a1.points[i].world_b), &a1.points[i].normal_b);
	}

	let threshold = a1.threshold;
	for idx in 0..len {
		let i = len - idx - 1;
		let distance = a1.points[i].distance;
		if !(threshold < distance) {
			let v4 = if !(0.0 < v9[i]) || !(distance < 0.0) { a4 } else { a5 };
			let v3 = length2(&(
				a1.points[i].world_b - (
					a1.points[i].world_a - a1.points[i].normal_b * distance
				)
			));
			if v4 * v4 * threshold * threshold < v3 {
				remove_contact_point(a1, i);
			}
		} else {
			remove_contact_point(a1, i);
		}
	}
}

pub fn add_contact_point(a1: &mut Manifold, a2: &Transform,
	a3: &Transform, a4: &Vector3, a5: &Vector3, a6: f32) -> usize
{
	if a1.threshold < a6 {
		return 0;
	}

	let v1 = a4 * a6 + a5;
	if v1.x.abs() < f32::INFINITY
		&& v1.y.abs() < f32::INFINITY
		&& v1.z.abs() < f32::INFINITY
	{
		let mut p = Point {
			local_a: inv_xform(a2, &v1),
			local_b: inv_xform(a3, a5),
			world_a: v1,
			world_b: a5.clone(),
			normal_b: a4.clone(),
			user: None,
			distance: a6.clone()
		};
		
		let opt = get_cache_entry(a1, &p);
		match opt {
			Some(i) => {
				p.user = a1.points[i].user;
				a1.points[i] = p;
				return i;
			},
			None => {
				return add_manifold_point(a1, p);
			}
		}
	}
	return 0;
}