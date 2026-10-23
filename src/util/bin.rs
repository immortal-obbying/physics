use std::num::Wrapping;

static TABLE: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

pub fn encode(buffer: Box<[u8]>) -> String {
	let length = buffer.len();
	let capacity = ((length + 2) / 3) * 4;
	let mut result = String::with_capacity(capacity);

	let mut v12 = 0;

	while v12 + 3 <= length {
		let v5 = buffer[v12];
		let v6 = buffer[v12 + 1];
		let v7 = buffer[v12 + 2];

		result.push(TABLE[(v5 >> 2) as usize] as char);
		result.push(TABLE[(((v5 & 3) << 4) | (v6 >> 4)) as usize] as char);
		result.push(TABLE[(((v6 & 15) << 2) | (v7 >> 6)) as usize] as char);
		result.push(TABLE[(v7 & 63) as usize] as char);

		v12 += 3;
	}

	let remaining = length - v12;

	if remaining == 1 {
		let v5 = buffer[v12];

		result.push(TABLE[(v5 >> 2) as usize] as char);
		result.push(TABLE[((v5 & 3) << 4) as usize] as char);
		result.push('=');
		result.push('=');
	} else if remaining == 2 {
		let v6 = buffer[v12];
		let v7 = buffer[v12 + 1];

		result.push(TABLE[(v6 >> 2) as usize] as char);
		result.push(TABLE[(((v6 & 3) << 4) | (v7 >> 4)) as usize] as char);
		result.push(TABLE[((v7 & 15) << 2) as usize] as char);
		result.push('=');
	}

	return result;
}

pub fn decode(input: &str) -> Box<[u8]> {
	let bytes = input.as_bytes();
	let mut v7 = input.len();
	while v7 > 0 {
		if bytes[v7 - 1] != b'=' {
			break;
		}
		v7 -= 1;
	}

	let mut buffer: Vec<u8> = Vec::with_capacity((v7 * 3) >> 2);
	let mut v5: u32 = 0;
	let mut v10 = 0;
	for i in 0..v7 {
		let v1 = TABLE[bytes[i] as usize];
		if v1 != 255 {
			v5 = (v5 << 6) | (v1 as u32);
			v10 += 6;
			if v10 >= 8 {
				v10 -= 8;
				let v4 = ((v5 >> v10) & 255) as u8;
				buffer.push(v4);
			}
		}
	}

	let decoded: Box<[u8]> = buffer.into_boxed_slice();
	return decoded;
}

pub fn hash(data: Box<[u8]>) -> String {
	let mut v1 = Wrapping(2166136261u32);
	for i in 0..data.len() {
		v1 = v1 ^ Wrapping(data[i] as u32);
		v1 = (v1 << 24) + (v1 << 8) + v1 * Wrapping(147u32);
	}
	return format!("{:08x}", v1.0);
}