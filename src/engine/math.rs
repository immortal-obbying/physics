pub type Matrix = [f32;9];
pub type Quaternion = [f32;4];

pub fn transpose(a1: &Matrix) -> Matrix {
	return [a1[0], a1[3], a1[6], a1[1], a1[4], a1[7], a1[2], a1[5], a1[8]];
}