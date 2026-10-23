use crate::util::vector::Vector;

#[derive(Clone)]
pub struct Xform {
	pub r: [f32;9],
	pub p: Vector
}