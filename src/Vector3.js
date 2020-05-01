/**
 * A vector with three components.
 */

export class Vector3 {

	/**
	 * Constructs a new vector.
	 *
	 * @param {Number} [x=0] - The X component.
	 * @param {Number} [y=0] - The Y component.
	 * @param {Number} [z=0] - The Z component.
	 */

	constructor(x = 0, y = 0, z = 0) {

		/**
		 * The X component.
		 *
		 * @type {Number}
		 */

		this.x = x;

		/**
		 * The Y component.
		 *
		 * @type {Number}
		 */

		this.y = y;

		/**
		 * The Z component.
		 *
		 * @type {Number}
		 */

		this.z = z;

	}

	/**
	 * Sets the values of this vector
	 *
	 * @param {Number} x - The X component.
	 * @param {Number} y - The Y component.
	 * @param {Number} z - The Z component.
	 * @return {Vector3} This vector.
	 */

	set(x, y, z) {

		this.x = x;
		this.y = y;
		this.z = z;

		return this;

	}

	/**
	 * Randomizes the values of this vector
	 *
	 * @return {Vector3} This vector.
	 */

	random() {

		this.x = Math.random();
		this.y = Math.random();
		this.z = Math.random();

		return this;

	}

	/**
	 * Copies the values of another vector.
	 *
	 * @param {Vector3} v - A vector.
	 * @return {Vector3} This vector.
	 */

	copy(v) {

		this.x = v.x;
		this.y = v.y;
		this.z = v.z;

		return this;

	}

	/**
	 * Clones this vector.
	 *
	 * @return {Vector3} A clone of this vector.
	 */

	clone() {

		return new this.constructor(this.x, this.y, this.z);

	}

	/**
	 * Copies values from an array.
	 *
	 * @param {Number[]} array - An array.
	 * @param {Number} offset - An offset.
	 * @return {Vector3} This vector.
	 */

	fromArray(array, offset = 0) {

		this.x = array[offset];
		this.y = array[offset + 1];
		this.z = array[offset + 2];

		return this;

	}

	/**
	 * Stores this vector in an array.
	 *
	 * @param {Array} [array] - A target array.
	 * @param {Number} offset - An offset.
	 * @return {Number[]} The array.
	 */

	toArray(array = [], offset = 0) {

		array[offset] = this.x;
		array[offset + 1] = this.y;
		array[offset + 2] = this.z;

		return array;

	}

	/**
	 * Sets the values of this vector based on a spherical description.
	 *
	 * @param {Spherical} s - A spherical description.
	 * @return {Vector3} This vector.
	 */

	setFromSpherical(s) {

		return this.setFromSphericalCoords(s.radius, s.phi, s.theta);

	}

	/**
	 * Sets the values of this vector based on spherical coordinates.
	 *
	 * @param {Number} radius - The radius.
	 * @param {Number} phi - The polar angle.
	 * @param {Number} theta - The angle around the equator of the sphere.
	 * @return {Vector3} This vector.
	 */

	setFromSphericalCoords(radius, phi, theta) {

		const sinPhiRadius = Math.sin(phi) * radius;

		this.x = sinPhiRadius * Math.sin(theta);
		this.y = Math.cos(phi) * radius;
		this.z = sinPhiRadius * Math.cos(theta);

		return this;

	}

	/**
	 * Sets the values of this vector based on a cylindrical description.
	 *
	 * @param {Cylindrical} c - A cylindrical description.
	 * @return {Vector3} This vector.
	 */

	setFromCylindrical(c) {

		return this.setFromCylindricalCoords(c.radius, c.theta, c.y);

	}

	/**
	 * Sets the values of this vector based on cylindrical coordinates.
	 *
	 * @param {Number} radius - The radius.
	 * @param {Number} theta - Theta.
	 * @param {Number} y - The height.
	 * @return {Vector3} This vector.
	 */

	setFromCylindricalCoords(radius, theta, y) {

		this.x = radius * Math.sin(theta);
		this.y = y;
		this.z = radius * Math.cos(theta);

		return this;

	}

	/**
	 * Copies the values of a matrix column.
	 *
	 * @param {Matrix4} m - A 3x3 matrix.
	 * @param {Number} index - A column index of the range [0, 2].
	 * @return {Vector3} This vector.
	 */

	setFromMatrix3Column(m, index) {

		return this.fromArray(m.elements, index * 3);

	}

	/**
	 * Copies the values of a matrix column.
	 *
	 * @param {Matrix4} m - A 4x4 matrix.
	 * @param {Number} index - A column index of the range [0, 3].
	 * @return {Vector3} This vector.
	 */

	setFromMatrixColumn(m, index) {

		return this.fromArray(m.elements, index * 4);

	}

	/**
	 * Extracts the position from a matrix.
	 *
	 * @param {Matrix4} m - A 4x4 matrix.
	 * @return {Vector3} This vector.
	 */

	setFromMatrixPosition(m) {

		const me = m.elements;

		this.x = me[12];
		this.y = me[13];
		this.z = me[14];

		return this;

	}

	/**
	 * Extracts the scale from a matrix.
	 *
	 * @param {Matrix4} m - A 4x4 matrix.
	 * @return {Vector3} This vector.
	 */

	setFromMatrixScale(m) {

		const sx = this.setFromMatrixColumn(m, 0).length();
		const sy = this.setFromMatrixColumn(m, 1).length();
		const sz = this.setFromMatrixColumn(m, 2).length();

		this.x = sx;
		this.y = sy;
		this.z = sz;

		return this;

	}

	/**
	 * Adds a vector to this one.
	 *
	 * @param {Vector3} v - The vector to add.
	 * @return {Vector3} This vector.
	 */

	add(v) {

		this.x += v.x;
		this.y += v.y;
		this.z += v.z;

		return this;

	}

	/**
	 * Adds a scalar to this vector.
	 *
	 * @param {Number} s - The scalar to add.
	 * @return {Vector3} This vector.
	 */

	addScalar(s) {

		this.x += s;
		this.y += s;
		this.z += s;

		return this;

	}

	/**
	 * Sets this vector to the sum of two given vectors.
	 *
	 * @param {Vector3} a - A vector.
	 * @param {Vector3} b - Another vector.
	 * @return {Vector3} This vector.
	 */

	addVectors(a, b) {

		this.x = a.x + b.x;
		this.y = a.y + b.y;
		this.z = a.z + b.z;

		return this;

	}

	/**
	 * Adds a scaled vector to this one.
	 *
	 * @param {Vector3} v - The vector to scale and add.
	 * @param {Number} s - A scalar.
	 * @return {Vector3} This vector.
	 */

	addScaledVector(v, s) {

		this.x += v.x * s;
		this.y += v.y * s;
		this.z += v.z * s;

		return this;

	}

	/**
	 * Subtracts a vector from this vector.
	 *
	 * @param {Vector3} v - The vector to subtract.
	 * @return {Vector3} This vector.
	 */

	sub(v) {

		this.x -= v.x;
		this.y -= v.y;
		this.z -= v.z;

		return this;

	}

	/**
	 * Subtracts a scalar from this vector.
	 *
	 * @param {Number} s - The scalar to subtract.
	 * @return {Vector3} This vector.
	 */

	subScalar(s) {

		this.x -= s;
		this.y -= s;
		this.z -= s;

		return this;

	}

	/**
	 * Sets this vector to the difference between two given vectors.
	 *
	 * @param {Vector3} a - A vector.
	 * @param {Vector3} b - A second vector.
	 * @return {Vector3} This vector.
	 */

	subVectors(a, b) {

		this.x = a.x - b.x;
		this.y = a.y - b.y;
		this.z = a.z - b.z;

		return this;

	}

	/**
	 * Multiplies this vector with another vector.
	 *
	 * @param {Vector3} v - A vector.
	 * @return {Vector3} This vector.
	 */

	multiply(v) {

		this.x *= v.x;
		this.y *= v.y;
		this.z *= v.z;

		return this;

	}

	/**
	 * Multiplies this vector with a given scalar.
	 *
	 * @param {Number} s - A scalar.
	 * @return {Vector3} This vector.
	 */

	multiplyScalar(s) {

		this.x *= s;
		this.y *= s;
		this.z *= s;

		return this;

	}

	/**
	 * Sets this vector to the product of two given vectors.
	 *
	 * @param {Vector3} a - A vector.
	 * @param {Vector3} b - Another vector.
	 * @return {Vector3} This vector.
	 */

	multiplyVectors(a, b) {

		this.x = a.x * b.x;
		this.y = a.y * b.y;
		this.z = a.z * b.z;

		return this;

	}

	/**
	 * Divides this vector by another vector.
	 *
	 * @param {Vector3} v - A vector.
	 * @return {Vector3} This vector.
	 */

	divide(v) {

		this.x /= v.x;
		this.y /= v.y;
		this.z /= v.z;

		return this;

	}

	/**
	 * Divides this vector by a given scalar.
	 *
	 * @param {Number} s - A scalar.
	 * @return {Vector3} This vector.
	 */

	divideScalar(s) {

		this.x /= s;
		this.y /= s;
		this.z /= s;

		return this;

	}

	/**
	 * Sets this vector to the cross product of the given vectors.
	 *
	 * @param {Vector3} a - A vector.
	 * @param {Vector3} b - Another vector.
	 * @return {Vector3} This vector.
	 */

	crossVectors(a, b) {

		const ax = a.x, ay = a.y, az = a.z;
		const bx = b.x, by = b.y, bz = b.z;

		this.x = ay * bz - az * by;
		this.y = az * bx - ax * bz;
		this.z = ax * by - ay * bx;

		return this;

	}

	/**
	 * Calculates the cross product of this vector and the given one.
	 *
	 * @param {Vector3} v - A vector.
	 * @return {Vector3} This vector.
	 */

	cross(v) {

		return this.crossVectors(this, v);

	}

	/**
	 * Applies a matrix to this direction vector.
	 *
	 * @param {Matrix4} m - A matrix.
	 * @return {Vector3} This vector.
	 */

	transformDirection(m) {

		const x = this.x, y = this.y, z = this.z;
		const e = m.elements;

		this.x = e[0] * x + e[4] * y + e[8] * z;
		this.y = e[1] * x + e[5] * y + e[9] * z;
		this.z = e[2] * x + e[6] * y + e[10] * z;

		return this.normalize();

	}

	/**
	 * Applies a matrix to this vector.
	 *
	 * @param {Matrix3} m - A matrix.
	 * @return {Vector3} This vector.
	 */

	applyMatrix3(m) {

		const x = this.x, y = this.y, z = this.z;
		const e = m.elements;

		this.x = e[0] * x + e[3] * y + e[6] * z;
		this.y = e[1] * x + e[4] * y + e[7] * z;
		this.z = e[2] * x + e[5] * y + e[8] * z;

		return this;

	}

	/**
	 * Applies a normal matrix to this vector and normalizes it.
	 *
	 * @param {Matrix3} m - A normal matrix.
	 * @return {Vector3} This vector.
	 */

	applyNormalMatrix(m) {

		return this.applyMatrix3(m).normalize();

	}

	/**
	 * Applies a matrix to this vector.
	 *
	 * @param {Matrix4} m - A matrix.
	 * @return {Vector3} This vector.
	 */

	applyMatrix4(m) {

		const x = this.x, y = this.y, z = this.z;
		const e = m.elements;

		this.x = e[0] * x + e[4] * y + e[8] * z + e[12];
		this.y = e[1] * x + e[5] * y + e[9] * z + e[13];
		this.z = e[2] * x + e[6] * y + e[10] * z + e[14];

		return this;

	}

	/**
	 * Applies a quaternion to this vector.
	 *
	 * @param {Quaternion} q - A quaternion.
	 * @return {Vector3} This vector.
	 */

	applyQuaternion(q) {

		const x = this.x, y = this.y, z = this.z;
		const qx = q.x, qy = q.y, qz = q.z, qw = q.w;

		// Calculate: quaternion * vector.
		const ix = qw * x + qy * z - qz * y;
		const iy = qw * y + qz * x - qx * z;
		const iz = qw * z + qx * y - qy * x;
		const iw = -qx * x - qy * y - qz * z;

		// Calculate: result * inverse quaternion.
		this.x = ix * qw + iw * -qx + iy * -qz - iz * -qy;
		this.y = iy * qw + iw * -qy + iz * -qx - ix * -qz;
		this.z = iz * qw + iw * -qz + ix * -qy - iy * -qx;

		return this;

	}

	/**
	 * Negates this vector.
	 *
	 * @return {Vector3} This vector.
	 */

	negate() {

		this.x = -this.x;
		this.y = -this.y;
		this.z = -this.z;

		return this;

	}

	/**
	 * Calculates the dot product with another vector.
	 *
	 * @param {Vector3} v - A vector.
	 * @return {Number} The dot product.
	 */

	dot(v) {

		return this.x * v.x + this.y * v.y + this.z * v.z;

	}

	/**
	 * Reflects this vector. The given plane normal is assumed to be normalized.
	 *
	 * @param {Vector3} n - A normal.
	 * @return {Vector3} This vector.
	 */

	reflect(n) {

		const nx = n.x;
		const ny = n.y;
		const nz = n.z;

		this.sub(n.multiplyScalar(2 * this.dot(n)));

		// Restore the normal.
		n.set(nx, ny, nz);

		return this;

	}

	/**
	 * Computes the angle to the given vector.
	 *
	 * @param {Vector3} v - A vector.
	 * @return {Number} The angle in radians.
	 */

	angleTo(v) {

		const theta = this.dot(v) / (Math.sqrt(this.lengthSquared() * v.lengthSquared()));

		// Clamp to avoid numerical problems.
		return Math.acos(Math.min(Math.max(theta, -1), 1));

	}

	/**
	 * Calculates the Manhattan length of this vector.
	 *
	 * @return {Number} The length.
	 */

	manhattanLength() {

		return Math.abs(this.x) + Math.abs(this.y) + Math.abs(this.z);

	}

	/**
	 * Calculates the squared length of this vector.
	 *
	 * @return {Number} The squared length.
	 */

	lengthSquared() {

		return this.x * this.x + this.y * this.y + this.z * this.z;

	}

	/**
	 * Calculates the length of this vector.
	 *
	 * @return {Number} The length.
	 */

	length() {

		return Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z);

	}

	/**
	 * Calculates the Manhattan distance to a given vector.
	 *
	 * @param {Vector3} v - A vector.
	 * @return {Number} The distance.
	 */

	manhattanDistanceTo(v) {

		return Math.abs(this.x - v.x) + Math.abs(this.y - v.y) + Math.abs(this.z - v.z);

	}

	/**
	 * Calculates the squared distance to a given vector.
	 *
	 * @param {Vector3} v - A vector.
	 * @return {Number} The squared distance.
	 */

	distanceToSquared(v) {

		const dx = this.x - v.x;
		const dy = this.y - v.y;
		const dz = this.z - v.z;

		return dx * dx + dy * dy + dz * dz;

	}

	/**
	 * Calculates the distance to a given vector.
	 *
	 * @param {Vector3} v - A vector.
	 * @return {Number} The distance.
	 */

	distanceTo(v) {

		return Math.sqrt(this.distanceToSquared(v));

	}

	/**
	 * Normalizes this vector.
	 *
	 * @return {Vector3} This vector.
	 */

	normalize() {

		return this.divideScalar(this.length());

	}

	/**
	 * Sets the length of this vector.
	 *
	 * @param {Number} length - The new length.
	 * @return {Vector3} This vector.
	 */

	setLength(length) {

		return this.normalize().multiplyScalar(length);

	}

	/**
	 * Adopts the min value for each component of this vector and the given one.
	 *
	 * @param {Vector3} v - A vector.
	 * @return {Vector3} This vector.
	 */

	min(v) {

		this.x = Math.min(this.x, v.x);
		this.y = Math.min(this.y, v.y);
		this.z = Math.min(this.z, v.z);

		return this;

	}

	/**
	 * Adopts the max value for each component of this vector and the given one.
	 *
	 * @param {Vector3} v - A vector.
	 * @return {Vector3} This vector.
	 */

	max(v) {

		this.x = Math.max(this.x, v.x);
		this.y = Math.max(this.y, v.y);
		this.z = Math.max(this.z, v.z);

		return this;

	}

	/**
	 * Clamps this vector.
	 *
	 * @param {Vector3} min - The lower bounds. Assumed to be smaller than max.
	 * @param {Vector3} max - The upper bounds. Assumed to be greater than min.
	 * @return {Vector3} This vector.
	 */

	clamp(min, max) {

		this.x = Math.max(min.x, Math.min(max.x, this.x));
		this.y = Math.max(min.y, Math.min(max.y, this.y));
		this.z = Math.max(min.z, Math.min(max.z, this.z));

		return this;

	}

	/**
	 * Floors this vector.
	 *
	 * @return {Vector3} This vector.
	 */

	floor() {

		this.x = Math.floor(this.x);
		this.y = Math.floor(this.y);
		this.z = Math.floor(this.z);

		return this;

	}

	/**
	 * Ceils this vector.
	 *
	 * @return {Vector3} This vector.
	 */

	ceil() {

		this.x = Math.ceil(this.x);
		this.y = Math.ceil(this.y);
		this.z = Math.ceil(this.z);

		return this;

	}

	/**
	 * Rounds this vector.
	 *
	 * @return {Vector3} This vector.
	 */

	round() {

		this.x = Math.round(this.x);
		this.y = Math.round(this.y);
		this.z = Math.round(this.z);

		return this;

	}

	/**
	 * Lerps towards the given vector.
	 *
	 * @param {Vector3} v - The target vector.
	 * @param {Number} alpha - The lerp factor.
	 * @return {Vector3} This vector.
	 */

	lerp(v, alpha) {

		this.x += (v.x - this.x) * alpha;
		this.y += (v.y - this.y) * alpha;
		this.z += (v.z - this.z) * alpha;

		return this;

	}

	/**
	 * Sets this vector to the lerp result of the given vectors.
	 *
	 * @param {Vector3} v1 - A base vector.
	 * @param {Vector3} v2 - The target vector.
	 * @param {Number} alpha - The lerp factor.
	 * @return {Vector3} This vector.
	 */

	lerpVectors(v1, v2, alpha) {

		return this.subVectors(v2, v1).multiplyScalar(alpha).add(v1);

	}

	/**
	 * Checks if this vector equals the given one.
	 *
	 * @param {Vector3} v - A vector.
	 * @return {Boolean} Whether this vector equals the given one.
	 */

	equals(v) {

		return (v.x === this.x && v.y === this.y && v.z === this.z);

	}

}
