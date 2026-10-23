use crate::util::vector::Vector;

#[derive(Copy, Clone)]
pub struct Transform {
	pub r: [f32;9],
	pub p: Vector
}