use std::ops::{ Add, Sub, Mul, Div };

#[derive(Copy, Clone)]
pub struct Vector {
	pub x: f32,
	pub y: f32,
	pub z: f32
}
impl Vector {
	pub fn new(x: f32, y: f32, z: f32) -> Vector {
		Vector { x, y, z }
	}
}

impl Add for Vector {
	type Output = Vector;
	fn add(self, rhs: Vector) -> Vector {
		Vector { x: self.x + rhs.x, y: self.y + rhs.y, z: self.z + rhs.z }
	}
}
impl<'a> Add<Vector> for &'a Vector { type Output = Vector; fn add(self, rhs: Vector) -> Vector { *self + rhs } }
impl<'a> Add<&'a Vector> for Vector { type Output = Vector; fn add(self, rhs: &'a Vector) -> Vector { self + *rhs } }
impl<'a, 'b> Add<&'b Vector> for &'a Vector { type Output = Vector; fn add(self, rhs: &'b Vector) -> Vector { *self + *rhs } }

impl Sub for Vector {
	type Output = Vector;
	fn sub(self, rhs: Vector) -> Vector {
		Vector { x: self.x - rhs.x, y: self.y - rhs.y, z: self.z - rhs.z }
	}
}
impl<'a> Sub<Vector> for &'a Vector { type Output = Vector; fn sub(self, rhs: Vector) -> Vector { *self - rhs } }
impl<'a> Sub<&'a Vector> for Vector { type Output = Vector; fn sub(self, rhs: &'a Vector) -> Vector { self - *rhs } }
impl<'a, 'b> Sub<&'b Vector> for &'a Vector { type Output = Vector; fn sub(self, rhs: &'b Vector) -> Vector { *self - *rhs } }

impl Mul for Vector {
	type Output = Vector;
	fn mul(self, rhs: Vector) -> Vector {
		Vector { x: self.x * rhs.x, y: self.y * rhs.y, z: self.z * rhs.z }
	}
}
impl<'a> Mul<Vector> for &'a Vector { type Output = Vector; fn mul(self, rhs: Vector) -> Vector { *self * rhs } }
impl<'a> Mul<&'a Vector> for Vector { type Output = Vector; fn mul(self, rhs: &'a Vector) -> Vector { self * *rhs } }
impl<'a, 'b> Mul<&'b Vector> for &'a Vector { type Output = Vector; fn mul(self, rhs: &'b Vector) -> Vector { *self * *rhs } }

impl Mul<f32> for Vector {
	type Output = Vector;
	fn mul(self, rhs: f32) -> Vector {
		Vector { x: self.x * rhs, y: self.y * rhs, z: self.z * rhs }
	}
}
impl<'a> Mul<f32> for &'a Vector { type Output = Vector; fn mul(self, rhs: f32) -> Vector { *self * rhs } }
impl Mul<Vector> for f32 {
	type Output = Vector;
	fn mul(self, rhs: Vector) -> Vector {
		Vector { x: self * rhs.x, y: self * rhs.y, z: self * rhs.z }
	}
}
impl<'a> Mul<&'a Vector> for f32 { type Output = Vector; fn mul(self, rhs: &'a Vector) -> Vector { self * *rhs } }

impl Div for Vector {
	type Output = Vector;
	fn div(self, rhs: Vector) -> Vector {
		Vector { x: self.x / rhs.x, y: self.y / rhs.y, z: self.z / rhs.z }
	}
}
impl<'a> Div<Vector> for &'a Vector { type Output = Vector; fn div(self, rhs: Vector) -> Vector { *self / rhs } }
impl<'a> Div<&'a Vector> for Vector { type Output = Vector; fn div(self, rhs: &'a Vector) -> Vector { self / *rhs } }
impl<'a, 'b> Div<&'b Vector> for &'a Vector { type Output = Vector; fn div(self, rhs: &'b Vector) -> Vector { *self / *rhs } }

impl Div<f32> for Vector {
	type Output = Vector;
	fn div(self, rhs: f32) -> Vector {
		Vector { x: self.x / rhs, y: self.y / rhs, z: self.z / rhs }
	}
}
impl<'a> Div<f32> for &'a Vector { type Output = Vector; fn div(self, rhs: f32) -> Vector { *self / rhs } }
impl Div<Vector> for f32 {
	type Output = Vector;
	fn div(self, rhs: Vector) -> Vector {
		Vector { x: self / rhs.x, y: self / rhs.y, z: self / rhs.z }
	}
}
impl<'a> Div<&'a Vector> for f32 { type Output = Vector; fn div(self, rhs: &'a Vector) -> Vector { self / *rhs } }
