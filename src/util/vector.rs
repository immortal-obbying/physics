use std::ops::{ Add, Sub, Mul, Div };

#[derive(Copy, Clone)]
pub struct Vector3 {
	pub x: f32,
	pub y: f32,
	pub z: f32
}
impl Vector3 {
	pub fn new(x: f32, y: f32, z: f32) -> Vector3 {
		Vector3 { x, y, z }
	}
}

impl Add for Vector3 {
	type Output = Vector3;
	fn add(self, rhs: Vector3) -> Vector3 {
		Vector3 { x: self.x + rhs.x, y: self.y + rhs.y, z: self.z + rhs.z }
	}
}
impl<'a> Add<Vector3> for &'a Vector3 { type Output = Vector3; fn add(self, rhs: Vector3) -> Vector3 { *self + rhs } }
impl<'a> Add<&'a Vector3> for Vector3 { type Output = Vector3; fn add(self, rhs: &'a Vector3) -> Vector3 { self + *rhs } }
impl<'a, 'b> Add<&'b Vector3> for &'a Vector3 { type Output = Vector3; fn add(self, rhs: &'b Vector3) -> Vector3 { *self + *rhs } }

impl Sub for Vector3 {
	type Output = Vector3;
	fn sub(self, rhs: Vector3) -> Vector3 {
		Vector3 { x: self.x - rhs.x, y: self.y - rhs.y, z: self.z - rhs.z }
	}
}
impl<'a> Sub<Vector3> for &'a Vector3 { type Output = Vector3; fn sub(self, rhs: Vector3) -> Vector3 { *self - rhs } }
impl<'a> Sub<&'a Vector3> for Vector3 { type Output = Vector3; fn sub(self, rhs: &'a Vector3) -> Vector3 { self - *rhs } }
impl<'a, 'b> Sub<&'b Vector3> for &'a Vector3 { type Output = Vector3; fn sub(self, rhs: &'b Vector3) -> Vector3 { *self - *rhs } }

impl Mul for Vector3 {
	type Output = Vector3;
	fn mul(self, rhs: Vector3) -> Vector3 {
		Vector3 { x: self.x * rhs.x, y: self.y * rhs.y, z: self.z * rhs.z }
	}
}
impl<'a> Mul<Vector3> for &'a Vector3 { type Output = Vector3; fn mul(self, rhs: Vector3) -> Vector3 { *self * rhs } }
impl<'a> Mul<&'a Vector3> for Vector3 { type Output = Vector3; fn mul(self, rhs: &'a Vector3) -> Vector3 { self * *rhs } }
impl<'a, 'b> Mul<&'b Vector3> for &'a Vector3 { type Output = Vector3; fn mul(self, rhs: &'b Vector3) -> Vector3 { *self * *rhs } }

impl Mul<f32> for Vector3 {
	type Output = Vector3;
	fn mul(self, rhs: f32) -> Vector3 {
		Vector3 { x: self.x * rhs, y: self.y * rhs, z: self.z * rhs }
	}
}
impl<'a> Mul<f32> for &'a Vector3 { type Output = Vector3; fn mul(self, rhs: f32) -> Vector3 { *self * rhs } }
impl Mul<Vector3> for f32 {
	type Output = Vector3;
	fn mul(self, rhs: Vector3) -> Vector3 {
		Vector3 { x: self * rhs.x, y: self * rhs.y, z: self * rhs.z }
	}
}
impl<'a> Mul<&'a Vector3> for f32 { type Output = Vector3; fn mul(self, rhs: &'a Vector3) -> Vector3 { self * *rhs } }

impl Div for Vector3 {
	type Output = Vector3;
	fn div(self, rhs: Vector3) -> Vector3 {
		Vector3 { x: self.x / rhs.x, y: self.y / rhs.y, z: self.z / rhs.z }
	}
}
impl<'a> Div<Vector3> for &'a Vector3 { type Output = Vector3; fn div(self, rhs: Vector3) -> Vector3 { *self / rhs } }
impl<'a> Div<&'a Vector3> for Vector3 { type Output = Vector3; fn div(self, rhs: &'a Vector3) -> Vector3 { self / *rhs } }
impl<'a, 'b> Div<&'b Vector3> for &'a Vector3 { type Output = Vector3; fn div(self, rhs: &'b Vector3) -> Vector3 { *self / *rhs } }

impl Div<f32> for Vector3 {
	type Output = Vector3;
	fn div(self, rhs: f32) -> Vector3 {
		Vector3 { x: self.x / rhs, y: self.y / rhs, z: self.z / rhs }
	}
}
impl<'a> Div<f32> for &'a Vector3 { type Output = Vector3; fn div(self, rhs: f32) -> Vector3 { *self / rhs } }
impl Div<Vector3> for f32 {
	type Output = Vector3;
	fn div(self, rhs: Vector3) -> Vector3 {
		Vector3 { x: self / rhs.x, y: self / rhs.y, z: self / rhs.z }
	}
}
impl<'a> Div<&'a Vector3> for f32 { type Output = Vector3; fn div(self, rhs: &'a Vector3) -> Vector3 { self / *rhs } }
