/**
 * math-ds v1.1.4 build Thu Jan 23 2020
 * https://github.com/vanruesc/math-ds
 * Copyright 2020 Raoul van Rüschen
 * @license Zlib
 */
/**
 * A vector with three components.
 */

class Vector3 {

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

/**
 * A vector.
 *
 * @type {Vector3}
 * @private
 */

const v = new Vector3();

/**
 * A list of points.
 *
 * @type {Vector3[]}
 * @private
 */

const points = [
	new Vector3(),
	new Vector3(),
	new Vector3(),
	new Vector3(),
	new Vector3(),
	new Vector3(),
	new Vector3(),
	new Vector3()
];

/**
 * A 3D box.
 */

class Box3 {

	/**
	 * Constructs a new box.
	 *
	 * @param {Vector3} [min] - The lower bounds.
	 * @param {Vector3} [max] - The upper bounds.
	 */

	constructor(
		min = new Vector3(Infinity, Infinity, Infinity),
		max = new Vector3(-Infinity, -Infinity, -Infinity)
	) {

		/**
		 * The lower bounds.
		 *
		 * @type {Vector3}
		 */

		this.min = min;

		/**
		 * The upper bounds.
		 *
		 * @type {Vector3}
		 */

		this.max = max;

	}

	/**
	 * Sets the values of this box.
	 *
	 * @param {Vector3} min - The lower bounds.
	 * @param {Vector3} max - The upper bounds.
	 * @return {Box3} This box.
	 */

	set(min, max) {

		this.min.copy(min);
		this.max.copy(max);

		return this;

	}

	/**
	 * Copies the values of a given box.
	 *
	 * @param {Box3} b - A box.
	 * @return {Box3} This box.
	 */

	copy(b) {

		this.min.copy(b.min);
		this.max.copy(b.max);

		return this;

	}

	/**
	 * Clones this box.
	 *
	 * @return {Box3} A clone of this box.
	 */

	clone() {

		return new this.constructor().copy(this);

	}

	/**
	 * Makes this box empty.
	 *
	 * The lower bounds are set to infinity and the upper bounds to negative
	 * infinity to create an infinitely small box.
	 *
	 * @return {Box3} This box.
	 */

	makeEmpty() {

		this.min.x = this.min.y = this.min.z = Infinity;
		this.max.x = this.max.y = this.max.z = -Infinity;

		return this;

	}

	/**
	 * Indicates whether this box is truly empty.
	 *
	 * This is a more robust check for emptiness since the volume can get positive
	 * with two negative axes.
	 *
	 * @return {Box3} This box.
	 */

	isEmpty() {

		return (
			this.max.x < this.min.x ||
			this.max.y < this.min.y ||
			this.max.z < this.min.z
		);

	}

	/**
	 * Computes the center of this box.
	 *
	 * @param {Vector3} [target] - A target vector. If none is provided, a new one will be created.
	 * @return {Vector3} A vector that describes the center of this box.
	 */

	getCenter(target = new Vector3()) {

		return !this.isEmpty() ?
			target.addVectors(this.min, this.max).multiplyScalar(0.5) :
			target.set(0, 0, 0);

	}

	/**
	 * Computes the size of this box.
	 *
	 * @param {Vector3} [target] - A target vector. If none is provided, a new one will be created.
	 * @return {Vector3} A vector that describes the size of this box.
	 */

	getSize(target = new Vector3()) {

		return !this.isEmpty() ?
			target.subVectors(this.max, this.min) :
			target.set(0, 0, 0);

	}

	/**
	 * Computes the bounding box of the given sphere.
	 *
	 * @param {Sphere} sphere - A sphere.
	 * @return {Box3} This box.
	 */

	setFromSphere(sphere) {

		this.set(sphere.center, sphere.center);
		this.expandByScalar(sphere.radius);

		return this;

	}

	/**
	 * Expands this box by the given point.
	 *
	 * @param {Vector3} p - A point.
	 * @return {Box3} This box.
	 */

	expandByPoint(p) {

		this.min.min(p);
		this.max.max(p);

		return this;

	}

	/**
	 * Expands this box by the given vector.
	 *
	 * @param {Vector3} v - A vector.
	 * @return {Box3} This box.
	 */

	expandByVector(v) {

		this.min.sub(v);
		this.max.add(v);

		return this;

	}

	/**
	 * Expands this box by the given scalar.
	 *
	 * @param {Number} s - A scalar.
	 * @return {Box3} This box.
	 */

	expandByScalar(s) {

		this.min.addScalar(-s);
		this.max.addScalar(s);

		return this;

	}

	/**
	 * Defines this box by the given points.
	 *
	 * @param {Vector3[]} points - The points.
	 * @return {Box3} This box.
	 */

	setFromPoints(points) {

		let i, l;

		this.min.set(0, 0, 0);
		this.max.set(0, 0, 0);

		for(i = 0, l = points.length; i < l; ++i) {

			this.expandByPoint(points[i]);

		}

		return this;

	}

	/**
	 * Defines this box by the given center and size.
	 *
	 * @param {Vector3} center - The center.
	 * @param {Number} size - The size.
	 * @return {Box3} This box.
	 */

	setFromCenterAndSize(center, size) {

		const halfSize = v.copy(size).multiplyScalar(0.5);

		this.min.copy(center).sub(halfSize);
		this.max.copy(center).add(halfSize);

		return this;

	}

	/**
	 * Clamps the given point to the boundaries of this box.
	 *
	 * @param {Vector3} point - A point.
	 * @param {Vector3} [target] - A target vector. If none is provided, a new one will be created.
	 * @return {Vector3} The clamped point.
	 */

	clampPoint(point, target = new Vector3()) {

		return target.copy(point).clamp(this.min, this.max);

	}

	/**
	 * Calculates the distance from this box to the given point.
	 *
	 * @param {Vector3} p - A point.
	 * @return {Number} The distance.
	 */

	distanceToPoint(p) {

		const clampedPoint = v.copy(p).clamp(this.min, this.max);

		return clampedPoint.sub(p).length();

	}

	/**
	 * Applies the given matrix to this box.
	 *
	 * @param {Matrix4} m - The matrix.
	 * @return {Box3} This box.
	 */

	applyMatrix4(m) {

		const min = this.min;
		const max = this.max;

		if(!this.isEmpty()) {

			points[0].set(min.x, min.y, min.z).applyMatrix4(m);
			points[1].set(min.x, min.y, max.z).applyMatrix4(m);
			points[2].set(min.x, max.y, min.z).applyMatrix4(m);
			points[3].set(min.x, max.y, max.z).applyMatrix4(m);
			points[4].set(max.x, min.y, min.z).applyMatrix4(m);
			points[5].set(max.x, min.y, max.z).applyMatrix4(m);
			points[6].set(max.x, max.y, min.z).applyMatrix4(m);
			points[7].set(max.x, max.y, max.z).applyMatrix4(m);

			this.setFromPoints(points);

		}

		return this;

	}

	/**
	 * Translates this box.
	 *
	 * @param {Vector3} offset - The offset.
	 * @return {Box3} This box.
	 */

	translate(offset) {

		this.min.add(offset);
		this.max.add(offset);

		return this;

	}

	/**
	 * Intersects this box with the given one.
	 *
	 * @param {Box3} b - A box.
	 * @return {Box3} This box.
	 */

	intersect(b) {

		this.min.max(b.min);
		this.max.min(b.max);

		/* Ensure that if there is no overlap, the result is fully empty to prevent
		subsequent intersections to erroneously return valid values. */
		if(this.isEmpty()) {

			this.makeEmpty();

		}

		return this;

	}

	/**
	 * Expands this box by combining it with the given one.
	 *
	 * @param {Box3} b - A box.
	 * @return {Box3} This box.
	 */

	union(b) {

		this.min.min(b.min);
		this.max.max(b.max);

		return this;

	}

	/**
	 * Checks if the given point lies inside this box.
	 *
	 * @param {Vector3} p - A point.
	 * @return {Boolean} Whether this box contains the point.
	 */

	containsPoint(p) {

		const min = this.min;
		const max = this.max;

		return (
			p.x >= min.x &&
			p.y >= min.y &&
			p.z >= min.z &&
			p.x <= max.x &&
			p.y <= max.y &&
			p.z <= max.z
		);

	}

	/**
	 * Checks if the given box lies inside this box.
	 *
	 * @param {Box3} b - A box.
	 * @return {Boolean} Whether this box contains the given one.
	 */

	containsBox(b) {

		const tMin = this.min;
		const tMax = this.max;
		const bMin = b.min;
		const bMax = b.max;

		return (
			tMin.x <= bMin.x && bMax.x <= tMax.x &&
			tMin.y <= bMin.y && bMax.y <= tMax.y &&
			tMin.z <= bMin.z && bMax.z <= tMax.z
		);

	}

	/**
	 * Checks if this box intersects the given one.
	 *
	 * @param {Box3} b - A box.
	 * @return {Boolean} Whether the boxes intersect.
	 */

	intersectsBox(b) {

		const tMin = this.min;
		const tMax = this.max;
		const bMin = b.min;
		const bMax = b.max;

		return (
			bMax.x >= tMin.x &&
			bMax.y >= tMin.y &&
			bMax.z >= tMin.z &&
			bMin.x <= tMax.x &&
			bMin.y <= tMax.y &&
			bMin.z <= tMax.z
		);

	}

	/**
	 * Checks if this box intersects the given sphere.
	 *
	 * @param {Sphere} s - A sphere.
	 * @return {Boolean} Whether the box intersects the sphere.
	 */

	intersectsSphere(s) {

		// Find the point in this box that is closest to the sphere's center.
		const closestPoint = this.clampPoint(s.center, v);

		// If that point is inside the sphere, it intersects this box.
		return (closestPoint.distanceToSquared(s.center) <= (s.radius * s.radius));

	}

	/**
	 * Checks if this box intersects the given plane.
	 *
	 * Computes the minimum and maximum dot product values. If those values are on
	 * the same side (back or front) of the plane, then there is no intersection.
	 *
	 * @param {Plane} p - A plane.
	 * @return {Boolean} Whether the box intersects the plane.
	 */

	intersectsPlane(p) {

		let min, max;

		if(p.normal.x > 0) {

			min = p.normal.x * this.min.x;
			max = p.normal.x * this.max.x;

		} else {

			min = p.normal.x * this.max.x;
			max = p.normal.x * this.min.x;

		}

		if(p.normal.y > 0) {

			min += p.normal.y * this.min.y;
			max += p.normal.y * this.max.y;

		} else {

			min += p.normal.y * this.max.y;
			max += p.normal.y * this.min.y;

		}

		if(p.normal.z > 0) {

			min += p.normal.z * this.min.z;
			max += p.normal.z * this.max.z;

		} else {

			min += p.normal.z * this.max.z;
			max += p.normal.z * this.min.z;

		}

		return (min <= -p.constant && max >= -p.constant);

	}

	/**
	 * Checks if this box equals the given one.
	 *
	 * @param {Box3} b - A box.
	 * @return {Boolean} Whether this box equals the given one.
	 */

	equals(b) {

		return (b.min.equals(this.min) && b.max.equals(this.max));

	}

}

/**
 * A box.
 *
 * @type {Box3}
 * @private
 */

const box = new Box3();

/**
 * A vector.
 *
 * @type {Vector3}
 * @private
 */

const v$1 = new Vector3();

/**
 * A sphere.
 */

class Sphere {

	/**
	 * Constructs a new sphere.
	 *
	 * @param {Vector3} [center] - The center.
	 * @param {Number} [radius] - The radius.
	 */

	constructor(center = new Vector3(), radius = 0) {

		/**
		 * The center.
		 *
		 * @type {Vector3}
		 */

		this.center = center;

		/**
		 * The radius.
		 *
		 * @type {Number}
		 */

		this.radius = radius;

	}

	/**
	 * Sets the center and the radius.
	 *
	 * @param {Vector3} center - The center.
	 * @param {Number} radius - The radius.
	 * @return {Sphere} This sphere.
	 */

	set(center, radius) {

		this.center.copy(center);
		this.radius = radius;

		return this;

	}

	/**
	 * Copies the given sphere.
	 *
	 * @param {Sphere} s - A sphere.
	 * @return {Sphere} This sphere.
	 */

	copy(s) {

		this.center.copy(s.center);
		this.radius = s.radius;

		return this;

	}

	/**
	 * Clones this sphere.
	 *
	 * @return {Sphere} The cloned sphere.
	 */

	clone() {

		return new this.constructor().copy(this);

	}

	/**
	 * Sets this sphere from points.
	 *
	 * @param {Vector3[]} points - The points.
	 * @param {Vector3} [center] - An optional center.
	 * @return {Sphere} This sphere.
	 */

	setFromPoints(points, center = box.setFromPoints(points).getCenter(this.center)) {

		let maxRadiusSq = 0;
		let i, l;

		for(i = 0, l = points.length; i < l; ++i) {

			maxRadiusSq = Math.max(maxRadiusSq, center.distanceToSquared(points[i]));

		}

		this.radius = Math.sqrt(maxRadiusSq);

		return this;

	}

	/**
	 * Computes the bounding sphere of the given box.
	 *
	 * @param {Box3} box - A box.
	 * @return {Sphere} This pshere.
	 */

	setFromBox(box) {

		box.getCenter(this.center);
		this.radius = box.getSize(v$1).length() * 0.5;

		return this;

	}

	/**
	 * Checks if this sphere is empty.
	 *
	 * @return {Boolean} Whether this sphere is empty.
	 */

	isEmpty() {

		return (this.radius <= 0);

	}

	/**
	 * Translates this sphere.
	 *
	 * @param {Number} offset - An offset.
	 * @return {Sphere} This sphere.
	 */

	translate(offset) {

		this.center.add(offset);

		return this;

	}

	/**
	 * Clamps the given point to this sphere.
	 *
	 * @param {Vector3} p - A point.
	 * @param {Vector3} [target] - A target vector. If none is provided, a new one will be created.
	 * @return {Vector3} The clamped point.
	 */

	clampPoint(p, target = new Vector3()) {

		const deltaLengthSq = this.center.distanceToSquared(p);

		target.copy(p);

		if(deltaLengthSq > (this.radius * this.radius)) {

			target.sub(this.center).normalize();
			target.multiplyScalar(this.radius).add(this.center);

		}

		return target;

	}

	/**
	 * Calculates the distance from this sphere to the given point.
	 *
	 * @param {Vector3} p - A point.
	 * @return {Number} The distance.
	 */

	distanceToPoint(p) {

		return (p.distanceTo(this.center) - this.radius);

	}

	/**
	 * Checks if the given point lies inside this sphere.
	 *
	 * @param {Vector3} p - A point.
	 * @return {Boolean} Whether this sphere contains the point.
	 */

	containsPoint(p) {

		return (p.distanceToSquared(this.center) <= (this.radius * this.radius));

	}

	/**
	 * Checks if the this sphere intersects with the given one.
	 *
	 * @param {Sphere} s - A sphere.
	 * @return {Boolean} Whether this sphere intersects with the given one.
	 */

	intersectsSphere(s) {

		const radiusSum = this.radius + s.radius;

		return s.center.distanceToSquared(this.center) <= (radiusSum * radiusSum);

	}

	/**
	 * Checks if the this sphere intersects with the given box.
	 *
	 * @param {Box3} b - A box.
	 * @return {Boolean} Whether this sphere intersects with the given box.
	 */

	intersectsBox(b) {

		return b.intersectsSphere(this);

	}

	/**
	 * Checks if the this sphere intersects with the given plane.
	 *
	 * @param {Plane} p - A plane.
	 * @return {Boolean} Whether this sphere intersects with the given plane.
	 */

	intersectsPlane(p) {

		return (Math.abs(p.distanceToPoint(this.center)) <= this.radius);

	}

	/**
	 * Checks if this sphere equals the given one.
	 *
	 * @param {Sphere} s - A sphere.
	 * @return {Boolean} Whether the spheres are equal.
	 */

	equals(s) {

		return (s.center.equals(this.center) && (s.radius === this.radius));

	}

}

/**
 * A vector with two components.
 */

class Vector2 {

	/**
	 * Constructs a new vector.
	 *
	 * @param {Number} [x=0] - The X component.
	 * @param {Number} [y=0] - The Y component.
	 */

	constructor(x = 0, y = 0) {

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

	}

	/**
	 * The width. This is an alias for X.
	 *
	 * @type {Number}
	 */

	get width() {

		return this.x;

	}

	/**
	 * Sets the width.
	 *
	 * @type {Number}
	 */

	set width(value) {

		return this.x = value;

	}

	/**
	 * The height. This is an alias for Y.
	 *
	 * @type {Number}
	 */

	get height() {

		return this.y;

	}

	/**
	 * Sets the height.
	 *
	 * @type {Number}
	 */

	set height(value) {

		return this.y = value;

	}

	/**
	 * Sets the values of this vector
	 *
	 * @param {Number} x - The X component.
	 * @param {Number} y - The Y component.
	 * @return {Vector2} This vector.
	 */

	set(x, y) {

		this.x = x;
		this.y = y;

		return this;

	}

	/**
	 * Copies the values of another vector.
	 *
	 * @param {Vector2} v - A vector.
	 * @return {Vector2} This vector.
	 */

	copy(v) {

		this.x = v.x;
		this.y = v.y;

		return this;

	}

	/**
	 * Clones this vector.
	 *
	 * @return {Vector2} A clone of this vector.
	 */

	clone() {

		return new this.constructor(this.x, this.y);

	}

	/**
	 * Copies values from an array.
	 *
	 * @param {Number[]} array - An array.
	 * @param {Number} offset - An offset.
	 * @return {Vector2} This vector.
	 */

	fromArray(array, offset = 0) {

		this.x = array[offset];
		this.y = array[offset + 1];

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

		return array;

	}

	/**
	 * Adds a vector to this one.
	 *
	 * @param {Vector2} v - The vector to add.
	 * @return {Vector2} This vector.
	 */

	add(v) {

		this.x += v.x;
		this.y += v.y;

		return this;

	}

	/**
	 * Adds a scalar to this vector.
	 *
	 * @param {Number} s - The scalar to add.
	 * @return {Vector2} This vector.
	 */

	addScalar(s) {

		this.x += s;
		this.y += s;

		return this;

	}

	/**
	 * Sets this vector to the sum of two given vectors.
	 *
	 * @param {Vector2} a - A vector.
	 * @param {Vector2} b - Another vector.
	 * @return {Vector2} This vector.
	 */

	addVectors(a, b) {

		this.x = a.x + b.x;
		this.y = a.y + b.y;

		return this;

	}

	/**
	 * Adds a scaled vector to this one.
	 *
	 * @param {Vector2} v - The vector to scale and add.
	 * @param {Number} s - A scalar.
	 * @return {Vector2} This vector.
	 */

	addScaledVector(v, s) {

		this.x += v.x * s;
		this.y += v.y * s;

		return this;

	}

	/**
	 * Subtracts a vector from this vector.
	 *
	 * @param {Vector2} v - The vector to subtract.
	 * @return {Vector2} This vector.
	 */

	sub(v) {

		this.x -= v.x;
		this.y -= v.y;

		return this;

	}

	/**
	 * Subtracts a scalar from this vector.
	 *
	 * @param {Number} s - The scalar to subtract.
	 * @return {Vector2} This vector.
	 */

	subScalar(s) {

		this.x -= s;
		this.y -= s;

		return this;

	}

	/**
	 * Sets this vector to the difference between two given vectors.
	 *
	 * @param {Vector2} a - A vector.
	 * @param {Vector2} b - A second vector.
	 * @return {Vector2} This vector.
	 */

	subVectors(a, b) {

		this.x = a.x - b.x;
		this.y = a.y - b.y;

		return this;

	}

	/**
	 * Multiplies this vector with another vector.
	 *
	 * @param {Vector2} v - A vector.
	 * @return {Vector2} This vector.
	 */

	multiply(v) {

		this.x *= v.x;
		this.y *= v.y;

		return this;

	}

	/**
	 * Multiplies this vector with a given scalar.
	 *
	 * @param {Number} s - A scalar.
	 * @return {Vector2} This vector.
	 */

	multiplyScalar(s) {

		this.x *= s;
		this.y *= s;

		return this;

	}

	/**
	 * Divides this vector by another vector.
	 *
	 * @param {Vector2} v - A vector.
	 * @return {Vector2} This vector.
	 */

	divide(v) {

		this.x /= v.x;
		this.y /= v.y;

		return this;

	}

	/**
	 * Divides this vector by a given scalar.
	 *
	 * @param {Number} s - A scalar.
	 * @return {Vector2} This vector.
	 */

	divideScalar(s) {

		this.x /= s;
		this.y /= s;

		return this;

	}

	/**
	 * Applies the given matrix to this vector.
	 *
	 * @param {Matrix3} m - A matrix.
	 * @return {Vector2} This vector.
	 */

	applyMatrix3(m) {

		const x = this.x, y = this.y;
		const e = m.elements;

		this.x = e[0] * x + e[3] * y + e[6];
		this.y = e[1] * x + e[4] * y + e[7];

		return this;

	}

	/**
	 * Calculates the dot product with another vector.
	 *
	 * @param {Vector2} v - A vector.
	 * @return {Number} The dot product.
	 */

	dot(v) {

		return this.x * v.x + this.y * v.y;

	}

	/**
	 * Calculates the cross product with another vector.
	 *
	 * This method calculates a scalar that would result from a regular 3D cross
	 * product of the input vectors, while taking their Z values implicitly as 0.
	 *
	 * @param {Vector2} v - A vector.
	 * @return {Number} The cross product.
	 */

	cross(v) {

		return this.x * v.y - this.y * v.x;

	}

	/**
	 * Calculates the Manhattan length of this vector.
	 *
	 * @return {Number} The length.
	 */

	manhattanLength() {

		return Math.abs(this.x) + Math.abs(this.y);

	}

	/**
	 * Calculates the squared length of this vector.
	 *
	 * @return {Number} The squared length.
	 */

	lengthSquared() {

		return this.x * this.x + this.y * this.y;

	}

	/**
	 * Calculates the length of this vector.
	 *
	 * @return {Number} The length.
	 */

	length() {

		return Math.sqrt(this.x * this.x + this.y * this.y);

	}

	/**
	 * Calculates the Manhattan distance to a given vector.
	 *
	 * @param {Vector2} v - A vector.
	 * @return {Number} The squared distance.
	 */

	manhattanDistanceTo(v) {

		return Math.abs(this.x - v.x) + Math.abs(this.y - v.y);

	}

	/**
	 * Calculates the squared distance to a given vector.
	 *
	 * @param {Vector2} v - A vector.
	 * @return {Number} The squared distance.
	 */

	distanceToSquared(v) {

		const dx = this.x - v.x;
		const dy = this.y - v.y;

		return dx * dx + dy * dy;

	}

	/**
	 * Calculates the distance to a given vector.
	 *
	 * @param {Vector2} v - A vector.
	 * @return {Number} The distance.
	 */

	distanceTo(v) {

		return Math.sqrt(this.distanceToSquared(v));

	}

	/**
	 * Normalizes this vector.
	 *
	 * @return {Vector2} This vector.
	 */

	normalize() {

		return this.divideScalar(this.length());

	}

	/**
	 * Sets the length of this vector.
	 *
	 * @param {Number} length - The new length.
	 * @return {Vector2} This vector.
	 */

	setLength(length) {

		return this.normalize().multiplyScalar(length);

	}

	/**
	 * Adopts the min value for each component of this vector and the given one.
	 *
	 * @param {Vector2} v - A vector.
	 * @return {Vector2} This vector.
	 */

	min(v) {

		this.x = Math.min(this.x, v.x);
		this.y = Math.min(this.y, v.y);

		return this;

	}

	/**
	 * adopts the max value for each component of this vector and the given one.
	 *
	 * @param {Vector2} v - A vector.
	 * @return {Vector2} This vector.
	 */

	max(v) {

		this.x = Math.max(this.x, v.x);
		this.y = Math.max(this.y, v.y);

		return this;

	}

	/**
	 * Clamps this vector.
	 *
	 * @param {Vector2} min - A vector, assumed to be smaller than max.
	 * @param {Vector2} max - A vector, assumed to be greater than min.
	 * @return {Vector2} This vector.
	 */

	clamp(min, max) {

		this.x = Math.max(min.x, Math.min(max.x, this.x));
		this.y = Math.max(min.y, Math.min(max.y, this.y));

		return this;

	}

	/**
	 * Floors this vector.
	 *
	 * @return {Vector2} This vector.
	 */

	floor() {

		this.x = Math.floor(this.x);
		this.y = Math.floor(this.y);

		return this;

	}

	/**
	 * Ceils this vector.
	 *
	 * @return {Vector2} This vector.
	 */

	ceil() {

		this.x = Math.ceil(this.x);
		this.y = Math.ceil(this.y);

		return this;

	}

	/**
	 * Rounds this vector.
	 *
	 * @return {Vector2} This vector.
	 */

	round() {

		this.x = Math.round(this.x);
		this.y = Math.round(this.y);

		return this;

	}

	/**
	 * Negates this vector.
	 *
	 * @return {Vector2} This vector.
	 */

	negate() {

		this.x = -this.x;
		this.y = -this.y;

		return this;

	}

	/**
	 * Computes the angle in radians with respect to the positive X-axis.
	 *
	 * @return {Number} The angle.
	 */

	angle() {

		let angle = Math.atan2(this.y, this.x);

		if(angle < 0) {

			angle += 2 * Math.PI;

		}

		return angle;

	}

	/**
	 * Lerps towards the given vector.
	 *
	 * @param {Vector2} v - The target vector.
	 * @param {Number} alpha - The lerp factor.
	 * @return {Vector2} This vector.
	 */

	lerp(v, alpha) {

		this.x += (v.x - this.x) * alpha;
		this.y += (v.y - this.y) * alpha;

		return this;

	}

	/**
	 * Sets this vector to the lerp result of the given vectors.
	 *
	 * @param {Vector2} v1 - A base vector.
	 * @param {Vector2} v2 - The target vector.
	 * @param {Number} alpha - The lerp factor.
	 * @return {Vector2} This vector.
	 */

	lerpVectors(v1, v2, alpha) {

		return this.subVectors(v2, v1).multiplyScalar(alpha).add(v1);

	}

	/**
	 * Rotates this vector around a given center.
	 *
	 * @param {Vector2} center - The center.
	 * @param {Number} angle - The rotation in radians.
	 * @return {Vector2} This vector.
	 */

	rotateAround(center, angle) {

		const c = Math.cos(angle), s = Math.sin(angle);

		const x = this.x - center.x;
		const y = this.y - center.y;

		this.x = x * c - y * s + center.x;
		this.y = x * s + y * c + center.y;

		return this;

	}

	/**
	 * Checks if this vector equals the given one.
	 *
	 * @param {Vector2} v - A vector.
	 * @return {Boolean} Whether this vector equals the given one.
	 */

	equals(v) {

		return (v.x === this.x && v.y === this.y);

	}

}

/**
 * A vector.
 *
 * @type {Vector2}
 * @private
 */

const v$2 = new Vector2();

/**
 * A 2D box.
 */

class Box2 {

	/**
	 * Constructs a new box.
	 *
	 * @param {Vector2} [min] - The lower bounds.
	 * @param {Vector2} [max] - The upper bounds.
	 */

	constructor(
		min = new Vector2(Infinity, Infinity),
		max = new Vector2(-Infinity, -Infinity)
	) {

		/**
		 * The lower bounds.
		 *
		 * @type {Vector2}
		 */

		this.min = min;

		/**
		 * The upper bounds.
		 *
		 * @type {Vector2}
		 */

		this.max = max;

	}

	/**
	 * Sets the values of this box.
	 *
	 * @param {Vector2} min - The lower bounds.
	 * @param {Vector2} max - The upper bounds.
	 * @return {Box2} This box.
	 */

	set(min, max) {

		this.min.copy(min);
		this.max.copy(max);

		return this;

	}

	/**
	 * Copies the values of a given box.
	 *
	 * @param {Box2} b - A box.
	 * @return {Box2} This box.
	 */

	copy(b) {

		this.min.copy(b.min);
		this.max.copy(b.max);

		return this;

	}

	/**
	 * Clones this box.
	 *
	 * @return {Box2} A clone of this box.
	 */

	clone() {

		return new this.constructor().copy(this);

	}

	/**
	 * Makes this box empty.
	 *
	 * The lower bounds are set to infinity and the upper bounds to negative
	 * infinity to create an infinitely small box.
	 *
	 * @return {Box2} This box.
	 */

	makeEmpty() {

		this.min.x = this.min.y = Infinity;
		this.max.x = this.max.y = -Infinity;

		return this;

	}

	/**
	 * Indicates whether this box is truly empty.
	 *
	 * This is a more robust check for emptiness since the volume can get positive
	 * with two negative axes.
	 *
	 * @return {Box2} This box.
	 */

	isEmpty() {

		return (
			this.max.x < this.min.x ||
			this.max.y < this.min.y
		);

	}

	/**
	 * Computes the center of this box.
	 *
	 * @param {Vector2} [target] - A target vector. If none is provided, a new one will be created.
	 * @return {Vector2} A vector that describes the center of this box.
	 */

	getCenter(target = new Vector2()) {

		return !this.isEmpty() ?
			target.addVectors(this.min, this.max).multiplyScalar(0.5) :
			target.set(0, 0);

	}

	/**
	 * Computes the size of this box.
	 *
	 * @param {Vector2} [target] - A target vector. If none is provided, a new one will be created.
	 * @return {Vector2} A vector that describes the size of this box.
	 */

	getSize(target = new Vector2()) {

		return !this.isEmpty() ?
			target.subVectors(this.max, this.min) :
			target.set(0, 0);

	}

	/**
	 * Computes the bounding sphere of this box.
	 *
	 * @param {Sphere} [target] - A target sphere. If none is provided, a new one will be created.
	 * @return {Sphere} The bounding sphere of this box.
	 */

	getBoundingSphere(target = new Sphere()) {

		this.getCenter(target.center);

		target.radius = this.getSize(v$2).length() * 0.5;

		return target;

	}

	/**
	 * Expands this box by the given point.
	 *
	 * @param {Vector2} p - A point.
	 * @return {Box2} This box.
	 */

	expandByPoint(p) {

		this.min.min(p);
		this.max.max(p);

		return this;

	}

	/**
	 * Expands this box by the given vector.
	 *
	 * @param {Vector2} v - A vector.
	 * @return {Box2} This box.
	 */

	expandByVector(v) {

		this.min.sub(v);
		this.max.add(v);

		return this;

	}

	/**
	 * Expands this box by the given scalar.
	 *
	 * @param {Number} s - A scalar.
	 * @return {Box2} This box.
	 */

	expandByScalar(s) {

		this.min.addScalar(-s);
		this.max.addScalar(s);

		return this;

	}

	/**
	 * Defines this box by the given points.
	 *
	 * @param {Vector2[]} points - The points.
	 * @return {Box2} This box.
	 */

	setFromPoints(points) {

		let i, l;

		this.min.set(0, 0);
		this.max.set(0, 0);

		for(i = 0, l = points.length; i < l; ++i) {

			this.expandByPoint(points[i]);

		}

		return this;

	}

	/**
	 * Defines this box by the given center and size.
	 *
	 * @param {Vector2} center - The center.
	 * @param {Number} size - The size.
	 * @return {Box2} This box.
	 */

	setFromCenterAndSize(center, size) {

		const halfSize = v$2.copy(size).multiplyScalar(0.5);

		this.min.copy(center).sub(halfSize);
		this.max.copy(center).add(halfSize);

		return this;

	}

	/**
	 * Clamps the given point to the boundaries of this box.
	 *
	 * @param {Vector2} point - A point.
	 * @param {Vector2} [target] - A target vector. If none is provided, a new one will be created.
	 * @return {Vector2} The clamped point.
	 */

	clampPoint(point, target = new Vector2()) {

		return target.copy(point).clamp(this.min, this.max);

	}

	/**
	 * Calculates the distance from this box to the given point.
	 *
	 * @param {Vector2} p - A point.
	 * @return {Number} The distance.
	 */

	distanceToPoint(p) {

		const clampedPoint = v$2.copy(p).clamp(this.min, this.max);

		return clampedPoint.sub(p).length();

	}

	/**
	 * Translates this box.
	 *
	 * @param {Vector2} offset - The offset.
	 * @return {Box2} This box.
	 */

	translate(offset) {

		this.min.add(offset);
		this.max.add(offset);

		return this;

	}

	/**
	 * Intersects this box with the given one.
	 *
	 * @param {Box2} b - A box.
	 * @return {Box2} This box.
	 */

	intersect(b) {

		this.min.max(b.min);
		this.max.min(b.max);

		/* Ensure that if there is no overlap, the result is fully empty to prevent
		subsequent intersections to erroneously return valid values. */
		if(this.isEmpty()) {

			this.makeEmpty();

		}

		return this;

	}

	/**
	 * Expands this box by combining it with the given one.
	 *
	 * @param {Box2} b - A box.
	 * @return {Box2} This box.
	 */

	union(b) {

		this.min.min(b.min);
		this.max.max(b.max);

		return this;

	}

	/**
	 * Checks if the given point lies inside this box.
	 *
	 * @param {Vector2} p - A point.
	 * @return {Boolean} Whether this box contains the point.
	 */

	containsPoint(p) {

		const min = this.min;
		const max = this.max;

		return (
			p.x >= min.x &&
			p.y >= min.y &&
			p.x <= max.x &&
			p.y <= max.y
		);

	}

	/**
	 * Checks if the given box lies inside this box.
	 *
	 * @param {Box2} b - A box.
	 * @return {Boolean} Whether this box contains the given one.
	 */

	containsBox(b) {

		const tMin = this.min;
		const tMax = this.max;
		const bMin = b.min;
		const bMax = b.max;

		return (
			tMin.x <= bMin.x && bMax.x <= tMax.x &&
			tMin.y <= bMin.y && bMax.y <= tMax.y
		);

	}

	/**
	 * Checks if this box intersects the given one.
	 *
	 * @param {Box2} b - A box.
	 * @return {Boolean} Whether the boxes intersect.
	 */

	intersectsBox(b) {

		const tMin = this.min;
		const tMax = this.max;
		const bMin = b.min;
		const bMax = b.max;

		return (
			bMax.x >= tMin.x &&
			bMax.y >= tMin.y &&
			bMin.x <= tMax.x &&
			bMin.y <= tMax.y
		);

	}

	/**
	 * Checks if this box equals the given one.
	 *
	 * @param {Box2} b - A box.
	 * @return {Boolean} Whether this box equals the given one.
	 */

	equals(b) {

		return (b.min.equals(this.min) && b.max.equals(this.max));

	}

}

/**
 * A cylindrical coordinate system.
 *
 * For details see: https://en.wikipedia.org/wiki/Cylindrical_coordinate_system
 */

class Cylindrical {

	/**
	 * Constructs a new cylindrical system.
	 *
	 * @param {Number} [radius=1] - The radius of the cylinder.
	 * @param {Number} [theta=0] - Theta.
	 * @param {Number} [y=0] - The height.
	 */

	constructor(radius = 1, theta = 0, y = 0) {

		/**
		 * The distance from the origin to a point in the XZ-plane.
		 *
		 * @type {Number}
		 */

		this.radius = radius;

		/**
		 * The counterclockwise angle in the XZ-plane measured in radians from the
		 * positive Z-axis.
		 *
		 * @type {Number}
		 */

		this.theta = theta;

		/**
		 * The height above the XZ-plane.
		 *
		 * @type {Number}
		 */

		this.y = y;

	}

	/**
	 * Sets the values of this cylindrical system.
	 *
	 * @param {Number} radius - The radius.
	 * @param {Number} theta - Theta.
	 * @param {Number} y - The height.
	 * @return {Cylindrical} This cylindrical system.
	 */

	set(radius, theta, y) {

		this.radius = radius;
		this.theta = theta;
		this.y = y;

		return this;

	}

	/**
	 * Copies the values of the given cylindrical system.
	 *
	 * @param {Cylindrical} c - A cylindrical system.
	 * @return {Cylindrical} This cylindrical system.
	 */

	copy(c) {

		this.radius = c.radius;
		this.theta = c.theta;
		this.y = c.y;

		return this;

	}

	/**
	 * Clones this cylindrical system.
	 *
	 * @return {Cylindrical} The cloned cylindrical system.
	 */

	clone() {

		return new this.constructor().copy(this);

	}

	/**
	 * Sets the values of this cylindrical system.
	 *
	 * @param {Vector3} v - A vector.
	 * @return {Cylindrical} This cylindrical system.
	 */

	setFromVector3(v) {

		return this.setFromCartesianCoords(v.x, v.y, v.z);

	}

	/**
	 * Sets the values of this cylindrical system based on cartesian coordinates.
	 *
	 * @param {Number} x - The X coordinate.
	 * @param {Number} y - The Y coordinate.
	 * @param {Number} z - The Z coordinate.
	 * @return {Cylindrical} This cylindrical system.
	 */

	setFromCartesianCoords(x, y, z) {

		this.radius = Math.sqrt(x * x + z * z);
		this.theta = Math.atan2(x, z);
		this.y = y;

		return this;

	}

}

/**
 * A 3x3 matrix.
 */

class Matrix3 {

	/**
	 * Constructs a new matrix.
	 */

	constructor() {

		/**
		 * The matrix elements.
		 *
		 * @type {Float32Array}
		 */

		this.elements = new Float32Array([

			1, 0, 0,
			0, 1, 0,
			0, 0, 1

		]);

	}

	/**
	 * Sets the values of this matrix.
	 *
	 * @param {Number} m00 - The value of the first row, first column.
	 * @param {Number} m01 - The value of the first row, second column.
	 * @param {Number} m02 - The value of the first row, third column.
	 * @param {Number} m10 - The value of the second row, first column.
	 * @param {Number} m11 - The value of the second row, second column.
	 * @param {Number} m12 - The value of the second row, third column.
	 * @param {Number} m20 - The value of the third row, first column.
	 * @param {Number} m21 - The value of the third row, second column.
	 * @param {Number} m22 - The value of the third row, third column.
	 * @return {Matrix3} This matrix.
	 */

	set(m00, m01, m02, m10, m11, m12, m20, m21, m22) {

		const te = this.elements;

		te[0] = m00; te[3] = m01; te[6] = m02;
		te[1] = m10; te[4] = m11; te[7] = m12;
		te[2] = m20; te[5] = m21; te[8] = m22;

		return this;

	}

	/**
	 * Sets this matrix to the identity matrix.
	 *
	 * @return {Matrix3} This matrix.
	 */

	identity() {

		this.set(

			1, 0, 0,
			0, 1, 0,
			0, 0, 1

		);

		return this;

	}

	/**
	 * Copies the values of a given matrix.
	 *
	 * @param {Matrix3} matrix - A matrix.
	 * @return {Matrix3} This matrix.
	 */

	copy(matrix) {

		const me = matrix.elements;
		const te = this.elements;

		te[0] = me[0]; te[1] = me[1]; te[2] = me[2];
		te[3] = me[3]; te[4] = me[4]; te[5] = me[5];
		te[6] = me[6]; te[7] = me[7]; te[8] = me[8];

		return this;

	}

	/**
	 * Clones this matrix.
	 *
	 * @return {Matrix3} A clone of this matrix.
	 */

	clone() {

		return new this.constructor().fromArray(this.elements);

	}

	/**
	 * Copies the values of a given array.
	 *
	 * @param {Number[]} array - An array.
	 * @param {Number} [offset=0] - An offset into the array.
	 * @return {Matrix3} This matrix.
	 */

	fromArray(array, offset = 0) {

		const te = this.elements;

		let i;

		for(i = 0; i < 9; ++i) {

			te[i] = array[i + offset];

		}

		return this;

	}

	/**
	 * Stores this matrix in an array.
	 *
	 * @param {Number[]} [array] - A target array.
	 * @param {Number} [offset=0] - An offset into the array.
	 * @return {Number[]} The array.
	 */

	toArray(array = [], offset = 0) {

		const te = this.elements;

		let i;

		for(i = 0; i < 9; ++i) {

			array[i + offset] = te[i];

		}

		return array;

	}

	/**
	 * Sets this matrix to the product of the given matrices.
	 *
	 * @param {Matrix3} a - A matrix.
	 * @param {Matrix3} b - A matrix.
	 * @return {Matrix3} This matrix.
	 */

	multiplyMatrices(a, b) {

		const ae = a.elements;
		const be = b.elements;
		const te = this.elements;

		const a11 = ae[0], a12 = ae[3], a13 = ae[6];
		const a21 = ae[1], a22 = ae[4], a23 = ae[7];
		const a31 = ae[2], a32 = ae[5], a33 = ae[8];

		const b11 = be[0], b12 = be[3], b13 = be[6];
		const b21 = be[1], b22 = be[4], b23 = be[7];
		const b31 = be[2], b32 = be[5], b33 = be[8];

		te[0] = a11 * b11 + a12 * b21 + a13 * b31;
		te[3] = a11 * b12 + a12 * b22 + a13 * b32;
		te[6] = a11 * b13 + a12 * b23 + a13 * b33;

		te[1] = a21 * b11 + a22 * b21 + a23 * b31;
		te[4] = a21 * b12 + a22 * b22 + a23 * b32;
		te[7] = a21 * b13 + a22 * b23 + a23 * b33;

		te[2] = a31 * b11 + a32 * b21 + a33 * b31;
		te[5] = a31 * b12 + a32 * b22 + a33 * b32;
		te[8] = a31 * b13 + a32 * b23 + a33 * b33;

		return this;

	}

	/**
	 * Multiplies this matrix with a given one.
	 *
	 * @param {Matrix3} m - A matrix.
	 * @return {Matrix3} This matrix.
	 */

	multiply(m) {

		return this.multiplyMatrices(this, m);

	}

	/**
	 * Multiplies a given matrix with this one.
	 *
	 * @param {Matrix3} m - A matrix.
	 * @return {Matrix3} This matrix.
	 */

	premultiply(m) {

		return this.multiplyMatrices(m, this);

	}

	/**
	 * Multiplies this matrix with a given scalar.
	 *
	 * @param {Number} s - A scalar.
	 * @return {Matrix3} This matrix.
	 */

	multiplyScalar(s) {

		const te = this.elements;

		te[0] *= s; te[3] *= s; te[6] *= s;
		te[1] *= s; te[4] *= s; te[7] *= s;
		te[2] *= s; te[5] *= s; te[8] *= s;

		return this;

	}

	/**
	 * Calculates the determinant of this matrix.
	 *
	 * @return {Number} The determinant.
	 */

	determinant() {

		const te = this.elements;

		const a = te[0], b = te[1], c = te[2];
		const d = te[3], e = te[4], f = te[5];
		const g = te[6], h = te[7], i = te[8];

		return (

			a * e * i -
			a * f * h -
			b * d * i +
			b * f * g +
			c * d * h -
			c * e * g

		);

	}

	/**
	 * Inverts the given matrix and stores the result in this matrix.
	 *
	 * @param {Matrix3} matrix - The matrix that should be inverted.
	 * @return {Matrix3} This matrix.
	 */

	getInverse(matrix) {

		const me = matrix.elements;
		const te = this.elements;

		const n11 = me[0], n21 = me[1], n31 = me[2];
		const n12 = me[3], n22 = me[4], n32 = me[5];
		const n13 = me[6], n23 = me[7], n33 = me[8];

		const t11 = n33 * n22 - n32 * n23;
		const t12 = n32 * n13 - n33 * n12;
		const t13 = n23 * n12 - n22 * n13;

		const det = n11 * t11 + n21 * t12 + n31 * t13;

		let invDet;

		if(det !== 0) {

			invDet = 1.0 / det;

			te[0] = t11 * invDet;
			te[1] = (n31 * n23 - n33 * n21) * invDet;
			te[2] = (n32 * n21 - n31 * n22) * invDet;

			te[3] = t12 * invDet;
			te[4] = (n33 * n11 - n31 * n13) * invDet;
			te[5] = (n31 * n12 - n32 * n11) * invDet;

			te[6] = t13 * invDet;
			te[7] = (n21 * n13 - n23 * n11) * invDet;
			te[8] = (n22 * n11 - n21 * n12) * invDet;

		} else {

			console.error("Can't invert matrix, determinant is zero", matrix);

			this.identity();

		}

		return this;

	}

	/**
	 * Transposes this matrix.
	 *
	 * @return {Matrix3} This matrix.
	 */

	transpose() {

		const me = this.elements;

		let t;

		t = me[1]; me[1] = me[3]; me[3] = t;
		t = me[2]; me[2] = me[6]; me[6] = t;
		t = me[5]; me[5] = me[7]; me[7] = t;

		return this;

	}

	/**
	 * Scales this matrix.
	 *
	 * @param {Number} sx - The X scale.
	 * @param {Number} sy - The Y scale.
	 * @return {Matrix3} This matrix.
	 */

	scale(sx, sy) {

		const te = this.elements;

		te[0] *= sx; te[3] *= sx; te[6] *= sx;
		te[1] *= sy; te[4] *= sy; te[7] *= sy;

		return this;

	}

	/**
	 * Rotates this matrix.
	 *
	 * @param {Number} theta - The rotation.
	 * @return {Matrix3} This matrix.
	 */

	rotate(theta) {

		const c = Math.cos(theta);
		const s = Math.sin(theta);

		const te = this.elements;

		const a11 = te[0], a12 = te[3], a13 = te[6];
		const a21 = te[1], a22 = te[4], a23 = te[7];

		te[0] = c * a11 + s * a21;
		te[3] = c * a12 + s * a22;
		te[6] = c * a13 + s * a23;

		te[1] = -s * a11 + c * a21;
		te[4] = -s * a12 + c * a22;
		te[7] = -s * a13 + c * a23;

		return this;

	}

	/**
	 * Translates this matrix.
	 *
	 * @param {Number} tx - The X offset.
	 * @param {Number} ty - The Y offset.
	 * @return {Matrix3} This matrix.
	 */

	translate(tx, ty) {

		const te = this.elements;

		te[0] += tx * te[2]; te[3] += tx * te[5]; te[6] += tx * te[8];
		te[1] += ty * te[2]; te[4] += ty * te[5]; te[7] += ty * te[8];

		return this;

	}

	/**
	 * Checks if this matrix equals the given one.
	 *
	 * @param {Matrix3} m - A matrix.
	 * @return {Boolean} Whether the matrix are equal.
	 */

	equals(m) {

		const te = this.elements;
		const me = m.elements;

		let result = true;
		let i;

		for(i = 0; result && i < 9; ++i) {

			if(te[i] !== me[i]) {

				result = false;

			}

		}

		return result;

	}

}

/**
 * An enumeration of Euler rotation orders.
 *
 * @type {Object}
 * @property {String} XYZ - X -> Y -> Z.
 * @property {String} YZX - Y -> Z -> X.
 * @property {String} ZXY - Z -> X -> Y.
 * @property {String} XZY - X -> Z -> Y.
 * @property {String} YXZ - Y -> X -> Z.
 * @property {String} ZYX - Z -> Y -> X.
 */

const RotationOrder = {

	XYZ: "XYZ",
	YZX: "YZX",
	ZXY: "ZXY",
	XZY: "XZY",
	YXZ: "YXZ",
	ZYX: "ZYX"

};

/**
 * A quaternion.
 */

class Quaternion {

	/**
	 * Constructs a new quaternion.
	 *
	 * @param {Number} [x=0] - The X component.
	 * @param {Number} [y=0] - The Y component.
	 * @param {Number} [z=0] - The Z component.
	 * @param {Number} [w=0] - The W component.
	 */

	constructor(x = 0, y = 0, z = 0, w = 0) {

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

		/**
		 * The W component.
		 *
		 * @type {Number}
		 */

		this.w = w;

	}

	/**
	 * Sets the components of this quaternion.
	 *
	 * @param {Number} x - The X component.
	 * @param {Number} y - The Y component.
	 * @param {Number} z - The Z component.
	 * @param {Number} w - The W component.
	 * @return {Quaternion} This quaternion.
	 */

	set(x, y, z, w) {

		this.x = x;
		this.y = y;
		this.z = z;
		this.w = w;

		return this;

	}

	/**
	 * Copies the components of the given quaternion.
	 *
	 * @param {Quaternion} q - The quaternion.
	 * @return {Quaternion} This quaternion.
	 */

	copy(q) {

		this.x = q.x;
		this.y = q.y;
		this.z = q.z;
		this.w = q.w;

		return this;

	}

	/**
	 * Clones this quaternion.
	 *
	 * @return {Quaternion} The cloned quaternion.
	 */

	clone() {

		return new this.constructor(this.x, this.y, this.z, this.w);

	}

	/**
	 * Copies values from an array.
	 *
	 * @param {Number[]} array - An array.
	 * @param {Number} offset - An offset.
	 * @return {Quaternion} This quaternion.
	 */

	fromArray(array, offset = 0) {

		this.x = array[offset];
		this.y = array[offset + 1];
		this.z = array[offset + 2];
		this.w = array[offset + 3];

		return this;

	}

	/**
	 * Stores this quaternion in an array.
	 *
	 * @param {Array} [array] - A target array.
	 * @param {Number} offset - An offset.
	 * @return {Number[]} The array.
	 */

	toArray(array = [], offset = 0) {

		array[offset] = this.x;
		array[offset + 1] = this.y;
		array[offset + 2] = this.z;
		array[offset + 3] = this.w;

		return array;

	}

	/**
	 * Sets the components of this quaternion based on the given Euler angles.
	 *
	 * For more details see: https://goo.gl/XRD1kr
	 *
	 * @param {Euler} euler - The euler angles.
	 * @return {Quaternion} This quaternion.
	 */

	setFromEuler(euler) {

		const x = euler.x;
		const y = euler.y;
		const z = euler.z;

		const cos = Math.cos;
		const sin = Math.sin;

		const c1 = cos(x / 2);
		const c2 = cos(y / 2);
		const c3 = cos(z / 2);

		const s1 = sin(x / 2);
		const s2 = sin(y / 2);
		const s3 = sin(z / 2);

		switch(euler.order) {

			case RotationOrder.XYZ:
				this.x = s1 * c2 * c3 + c1 * s2 * s3;
				this.y = c1 * s2 * c3 - s1 * c2 * s3;
				this.z = c1 * c2 * s3 + s1 * s2 * c3;
				this.w = c1 * c2 * c3 - s1 * s2 * s3;
				break;

			case RotationOrder.YXZ:
				this.x = s1 * c2 * c3 + c1 * s2 * s3;
				this.y = c1 * s2 * c3 - s1 * c2 * s3;
				this.z = c1 * c2 * s3 - s1 * s2 * c3;
				this.w = c1 * c2 * c3 + s1 * s2 * s3;
				break;

			case RotationOrder.ZXY:
				this.x = s1 * c2 * c3 - c1 * s2 * s3;
				this.y = c1 * s2 * c3 + s1 * c2 * s3;
				this.z = c1 * c2 * s3 + s1 * s2 * c3;
				this.w = c1 * c2 * c3 - s1 * s2 * s3;
				break;

			case RotationOrder.ZYX:
				this.x = s1 * c2 * c3 - c1 * s2 * s3;
				this.y = c1 * s2 * c3 + s1 * c2 * s3;
				this.z = c1 * c2 * s3 - s1 * s2 * c3;
				this.w = c1 * c2 * c3 + s1 * s2 * s3;
				break;

			case RotationOrder.YZX:
				this.x = s1 * c2 * c3 + c1 * s2 * s3;
				this.y = c1 * s2 * c3 + s1 * c2 * s3;
				this.z = c1 * c2 * s3 - s1 * s2 * c3;
				this.w = c1 * c2 * c3 - s1 * s2 * s3;
				break;

			case RotationOrder.XZY:
				this.x = s1 * c2 * c3 - c1 * s2 * s3;
				this.y = c1 * s2 * c3 - s1 * c2 * s3;
				this.z = c1 * c2 * s3 + s1 * s2 * c3;
				this.w = c1 * c2 * c3 + s1 * s2 * s3;
				break;

		}

		return this;

	}

	/**
	 * Sets the components of this quaternion based on a given axis angle.
	 *
	 * For more information see:
	 *  http://www.euclideanspace.com/maths/geometry/rotations/conversions/angleToQuaternion/index.htm
	 *
	 * @param {Vector3} axis - The axis. Assumed to be normalized.
	 * @param {Number} angle - The angle in radians.
	 * @return {Quaternion} This quaternion.
	 */

	setFromAxisAngle(axis, angle) {

		const halfAngle = angle / 2.0;
		const s = Math.sin(halfAngle);

		this.x = axis.x * s;
		this.y = axis.y * s;
		this.z = axis.z * s;
		this.w = Math.cos(halfAngle);

		return this;

	}

	/**
	 * Sets the components of this quaternion based on a given rotation matrix.
	 *
	 * For more information see:
	 *  http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/index.htm
	 *
	 * @param {Matrix4} m - The rotation matrix. The upper 3x3 is assumed to be a pure rotation matrix (i.e. unscaled).
	 * @return {Quaternion} This quaternion.
	 */

	setFromRotationMatrix(m) {

		const te = m.elements;

		const m00 = te[0], m01 = te[4], m02 = te[8];
		const m10 = te[1], m11 = te[5], m12 = te[9];
		const m20 = te[2], m21 = te[6], m22 = te[10];

		const trace = m00 + m11 + m22;

		let s;

		if(trace > 0) {

			s = 0.5 / Math.sqrt(trace + 1.0);

			this.w = 0.25 / s;
			this.x = (m21 - m12) * s;
			this.y = (m02 - m20) * s;
			this.z = (m10 - m01) * s;

		} else if(m00 > m11 && m00 > m22) {

			s = 2.0 * Math.sqrt(1.0 + m00 - m11 - m22);

			this.w = (m21 - m12) / s;
			this.x = 0.25 * s;
			this.y = (m01 + m10) / s;
			this.z = (m02 + m20) / s;

		} else if(m11 > m22) {

			s = 2.0 * Math.sqrt(1.0 + m11 - m00 - m22);

			this.w = (m02 - m20) / s;
			this.x = (m01 + m10) / s;
			this.y = 0.25 * s;
			this.z = (m12 + m21) / s;

		} else {

			s = 2.0 * Math.sqrt(1.0 + m22 - m00 - m11);

			this.w = (m10 - m01) / s;
			this.x = (m02 + m20) / s;
			this.y = (m12 + m21) / s;
			this.z = 0.25 * s;

		}

		return this;

	}

	/**
	 * Sets the components of this quaternion based on unit vectors.
	 *
	 * @param {Vector3} vFrom - A unit vector. Assumed to be normalized.
	 * @param {Vector3} vTo - A unit vector. Assumed to be normalized.
	 * @return {Quaternion} This quaternion.
	 */

	setFromUnitVectors(vFrom, vTo) {

		let r = vFrom.dot(vTo) + 1;

		if(r < 1e-6) {

			r = 0;

			if(Math.abs(vFrom.x) > Math.abs(vFrom.z)) {

				this.x = -vFrom.y;
				this.y = vFrom.x;
				this.z = 0;
				this.w = r;

			} else {

				this.x = 0;
				this.y = -vFrom.z;
				this.z = vFrom.y;
				this.w = r;

			}

		} else {

			// crossVectors(vFrom, vTo)
			this.x = vFrom.y * vTo.z - vFrom.z * vTo.y;
			this.y = vFrom.z * vTo.x - vFrom.x * vTo.z;
			this.z = vFrom.x * vTo.y - vFrom.y * vTo.x;
			this.w = r;

		}

		return this.normalize();

	}

	/**
	 * Calculates the angle to another quaternion.
	 *
	 * @param {Quaternion} q - A quaternion.
	 * @return {Number} The angle in radians.
	 */

	angleTo(q) {

		return 2.0 * Math.acos(Math.abs(Math.min(Math.max(this.dot(q), -1.0), 1.0)));

	}

	/**
	 * Rotates this quaternion towards the given one by a given step size.
	 *
	 * @param {Quaternion} q - The target quaternion.
	 * @param {Number} step - The step size.
	 * @return {Quaternion} This quaternion.
	 */

	rotateTowards(q, step) {

		const angle = this.angleTo(q);

		if(angle !== 0.0) {

			this.slerp(q, Math.min(1.0, step / angle));

		}

		return this;

	}

	/**
	 * Inverts this quaternion. The quaternion is assumed to have unit length.
	 *
	 * @return {Quaternion} This quaternion.
	 */

	invert() {

		return this.conjugate();

	}

	/**
	 * Conjugates this quaternion.
	 *
	 * @return {Quaternion} This quaternion.
	 */

	conjugate() {

		this.x *= -1;
		this.y *= -1;
		this.z *= -1;

		return this;

	}

	/**
	 * Calculates the squared length of this quaternion.
	 *
	 * @return {Number} The squared length.
	 */

	lengthSquared() {

		return this.x * this.x + this.y * this.y + this.z * this.z + this.w * this.w;

	}

	/**
	 * Calculates the length of this quaternion.
	 *
	 * @return {Number} The length.
	 */

	length() {

		return Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z + this.w * this.w);

	}

	/**
	 * Normalizes this quaternion.
	 *
	 * @return {Quaternion} This quaternion.
	 */

	normalize() {

		const l = this.length();

		let invLength;

		if(l === 0) {

			this.x = 0;
			this.y = 0;
			this.z = 0;
			this.w = 1;

		} else {

			invLength = 1.0 / l;

			this.x = this.x * invLength;
			this.y = this.y * invLength;
			this.z = this.z * invLength;
			this.w = this.w * invLength;

		}

		return this;

	}

	/**
	 * Calculates the dot product with a given vector.
	 *
	 * @param {Vector4} v - A vector.
	 * @return {Number} The dot product.
	 */

	dot(v) {

		return this.x * v.x + this.y * v.y + this.z * v.z + this.w * v.w;

	}

	/**
	 * Multiplies the given quaternions and stores the result in this quaternion.
	 *
	 * For more details see:
	 *  http://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions/code/index.htm
	 *
	 * @param {Quaternion} a - A quaternion.
	 * @param {Quaternion} b - Another quaternion.
	 * @return {Quaternion} This quaternion.
	 */

	multiplyQuaternions(a, b) {

		const qax = a.x, qay = a.y, qaz = a.z, qaw = a.w;
		const qbx = b.x, qby = b.y, qbz = b.z, qbw = b.w;

		this.x = qax * qbw + qaw * qbx + qay * qbz - qaz * qby;
		this.y = qay * qbw + qaw * qby + qaz * qbx - qax * qbz;
		this.z = qaz * qbw + qaw * qbz + qax * qby - qay * qbx;
		this.w = qaw * qbw - qax * qbx - qay * qby - qaz * qbz;

		return this;

	}

	/**
	 * Multiplies this quaternion with the given one and stores the result in
	 * this quaternion.
	 *
	 * @param {Quaternion} q - A quaternion.
	 * @return {Quaternion} This quaternion.
	 */

	multiply(q) {

		return this.multiplyQuaternions(this, q);

	}

	/**
	 * Multiplies the given quaternion with this one and stores the result in
	 * this quaternion.
	 *
	 * @param {Quaternion} q - A quaternion.
	 * @return {Quaternion} This quaternion.
	 */

	premultiply(q) {

		return this.multiplyQuaternions(q, this);

	}

	/**
	 * Performs a spherical linear interpolation towards the given quaternion.
	 *
	 * For more details see:
	 *  http://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions/slerp/
	 *
	 * @param {Quaternion} q - A quaternion.
	 * @param {Number} t - The slerp factor.
	 * @return {Quaternion} This quaternion.
	 */

	slerp(q, t) {

		const x = this.x, y = this.y, z = this.z, w = this.w;

		let cosHalfTheta, sinHalfThetaSquared, sinHalfTheta, halfTheta;
		let s, ratioA, ratioB;

		if(t === 1) {

			this.copy(q);

		} else if(t > 0) {

			cosHalfTheta = w * q.w + x * q.x + y * q.y + z * q.z;

			if(cosHalfTheta < 0.0) {

				this.w = -q.w;
				this.x = -q.x;
				this.y = -q.y;
				this.z = -q.z;

				cosHalfTheta = -cosHalfTheta;

			} else {

				this.copy(q);

			}

			if(cosHalfTheta >= 1.0) {

				this.w = w;
				this.x = x;
				this.y = y;
				this.z = z;

			} else {

				sinHalfThetaSquared = 1.0 - cosHalfTheta * cosHalfTheta;
				s = 1.0 - t;

				if(sinHalfThetaSquared <= Number.EPSILON) {

					this.w = s * w + t * this.w;
					this.x = s * x + t * this.x;
					this.y = s * y + t * this.y;
					this.z = s * z + t * this.z;

					this.normalize();

				} else {

					sinHalfTheta = Math.sqrt(sinHalfThetaSquared);
					halfTheta = Math.atan2(sinHalfTheta, cosHalfTheta);
					ratioA = Math.sin(s * halfTheta) / sinHalfTheta;
					ratioB = Math.sin(t * halfTheta) / sinHalfTheta;

					this.w = (w * ratioA + this.w * ratioB);
					this.x = (x * ratioA + this.x * ratioB);
					this.y = (y * ratioA + this.y * ratioB);
					this.z = (z * ratioA + this.z * ratioB);

				}

			}

		}

		return this;

	}

	/**
	 * Checks if this quaternions equals the given one.
	 *
	 * @param {Quaternion} q - A quaternion.
	 * @return {Boolean} Whether the quaternions are equal.
	 */

	equals(q) {

		return (q.x === this.x) && (q.y === this.y) && (q.z === this.z) && (q.w === this.w);

	}

	/**
	 * Performs a spherical linear interpolation.
	 *
	 * @param {Quaternion} qa - The base quaternion.
	 * @param {Quaternion} qb - The target quaternion.
	 * @param {Quaternion} qr - A quaternion to store the result in.
	 * @param {Number} t - The slerp factor.
	 * @return {Quaternion} The resulting quaternion.
	 */

	static slerp(qa, qb, qr, t) {

		return qr.copy(qa).slerp(qb, t);

	}

	/**
	 * Performs an array-based spherical linear interpolation.
	 *
	 * @param {Number[]} dst - An array to store the result in.
	 * @param {Number} dstOffset - An offset into the destination array.
	 * @param {Number[]} src0 - An array that contains the base quaternion values.
	 * @param {Number} srcOffset0 - An offset into the base array.
	 * @param {Number[]} src1 - An array that contains the target quaternion values.
	 * @param {Number} srcOffset1 - An offset into the target array.
	 * @param {Number} t - The slerp factor.
	 */

	static slerpFlat(dst, dstOffset, src0, srcOffset0, src1, srcOffset1, t) {

		const x1 = src1[srcOffset1];
		const y1 = src1[srcOffset1 + 1];
		const z1 = src1[srcOffset1 + 2];
		const w1 = src1[srcOffset1 + 3];

		let x0 = src0[srcOffset0];
		let y0 = src0[srcOffset0 + 1];
		let z0 = src0[srcOffset0 + 2];
		let w0 = src0[srcOffset0 + 3];

		let s, f;
		let sin, cos, sqrSin;
		let dir, len, tDir;

		if(w0 !== w1 || x0 !== x1 || y0 !== y1 || z0 !== z1) {

			s = 1.0 - t;
			cos = x0 * x1 + y0 * y1 + z0 * z1 + w0 * w1;

			dir = (cos >= 0) ? 1 : -1;
			sqrSin = 1.0 - cos * cos;

			// Skip the Slerp for tiny steps to avoid numeric problems.
			if(sqrSin > Number.EPSILON) {

				sin = Math.sqrt(sqrSin);
				len = Math.atan2(sin, cos * dir);

				s = Math.sin(s * len) / sin;
				t = Math.sin(t * len) / sin;

			}

			tDir = t * dir;

			x0 = x0 * s + x1 * tDir;
			y0 = y0 * s + y1 * tDir;
			z0 = z0 * s + z1 * tDir;
			w0 = w0 * s + w1 * tDir;

			// Normalize in case a lerp has just been performed.
			if(s === 1.0 - t) {

				f = 1.0 / Math.sqrt(x0 * x0 + y0 * y0 + z0 * z0 + w0 * w0);

				x0 *= f;
				y0 *= f;
				z0 *= f;
				w0 *= f;

			}

		}

		dst[dstOffset] = x0;
		dst[dstOffset + 1] = y0;
		dst[dstOffset + 2] = z0;
		dst[dstOffset + 3] = w0;

	}

}

/**
 * Clamps the given value.
 *
 * @private
 * @param {Number} value - The value.
 * @param {Number} min - The lower limit.
 * @param {Number} max - The upper limit.
 * @return {Number} The clamped value.
 */

function clamp(value, min, max) {

	return Math.max(Math.min(value, max), min);

}

/**
 * A matrix.
 *
 * @type {Matrix3}
 * @private
 */

const m = new Matrix3();

/**
 * A quaternion.
 *
 * @type {Quaternion}
 * @private
 */

const q = new Quaternion();

/**
 * Euler angles.
 */

class Euler {

	/**
	 * Constructs a new set of Euler angles.
	 *
	 * @param {Number} [x=0] - The rotation around the X-axis.
	 * @param {Number} [y=0] - The rotation around the Y-axis.
	 * @param {Number} [z=0] - The rotation around the Z-axis.
	 */

	constructor(x = 0, y = 0, z = 0) {

		/**
		 * The rotation around the X-axis.
		 *
		 * @type {Number}
		 */

		this.x = x;

		/**
		 * The rotation around the Y-axis.
		 *
		 * @type {Number}
		 */

		this.y = y;

		/**
		 * The rotation around the Z-axis.
		 *
		 * @type {Number}
		 */

		this.z = z;

		/**
		 * The rotation order.
		 *
		 * @type {RotationOrder}
		 * @default Euler.defaultOrder
		 */

		this.order = Euler.defaultOrder;

	}

	/**
	 * Sets the Euler angles and rotation order.
	 *
	 * @param {Number} x - The rotation around the X-axis.
	 * @param {Number} y - The rotation around the Y-axis.
	 * @param {Number} z - The rotation around the Z-axis.
	 * @param {Number} order - The rotation order.
	 * @return {Euler} This set of Euler angles.
	 */

	set(x, y, z, order) {

		this.x = x;
		this.y = y;
		this.z = z;
		this.order = order;

		return this;

	}

	/**
	 * Copies the values of another set of Euler angles.
	 *
	 * @param {Euler} e - A set of Euler angles.
	 * @return {Euler} This set of Euler angles.
	 */

	copy(e) {

		this.x = e.x;
		this.y = e.y;
		this.z = e.z;
		this.order = e.order;

		return this;

	}

	/**
	 * Clones this set of Euler angles.
	 *
	 * @return {Euler} A clone of this set of Euler angles.
	 */

	clone() {

		return new this.constructor(this.x, this.y, this.z, this.order);

	}

	/**
	 * Copies angles and the rotation order from an array.
	 *
	 * @param {Number[]} array - An array.
	 * @param {Number} offset - An offset.
	 * @return {Euler} This set of Euler angles.
	 */

	fromArray(array, offset = 0) {

		this.x = array[offset];
		this.y = array[offset + 1];
		this.z = array[offset + 2];
		this.order = array[offset + 3];

		return this;

	}

	/**
	 * Stores this set of Euler angles and the rotation order in an array.
	 *
	 * @param {Array} [array] - A target array.
	 * @param {Number} offset - An offset.
	 * @return {Number[]} The array.
	 */

	toArray(array = [], offset = 0) {

		array[offset] = this.x;
		array[offset + 1] = this.y;
		array[offset + 2] = this.z;
		array[offset + 3] = this.order;

		return array;

	}

	/**
	 * Stores this set of Euler angles in a vector.
	 *
	 * @param {Vector3} [vector] - A target vector. If none is provided, a new one will be created.
	 * @return {Vector3} The vector.
	 */

	toVector3(vector = new Vector3()) {

		return vector.set(this.x, this.y, this.z);

	}

	/**
	 * Copies the rotation from a given matrix.
	 *
	 * @param {Matrix4} m - A rotation matrix. The upper 3x3 is assumed to be a pure rotation matrix (i.e. unscaled).
	 * @param {RotationOrder} [order] - An override rotation order.
	 * @return {Euler} This set of Euler angles.
	 */

	setFromRotationMatrix(m, order = this.order) {

		const te = m.elements;
		const m00 = te[0], m01 = te[4], m02 = te[8];
		const m10 = te[1], m11 = te[5], m12 = te[9];
		const m20 = te[2], m21 = te[6], m22 = te[10];

		const THRESHOLD = 1.0 - 1e-7;

		switch(order) {

			case RotationOrder.XYZ: {

				this.y = Math.asin(clamp(m02, -1, 1));

				if(Math.abs(m02) < THRESHOLD) {

					this.x = Math.atan2(-m12, m22);
					this.z = Math.atan2(-m01, m00);

				} else {

					this.x = Math.atan2(m21, m11);
					this.z = 0;

				}

				break;

			}

			case RotationOrder.YXZ: {

				this.x = Math.asin(-clamp(m12, -1, 1));

				if(Math.abs(m12) < THRESHOLD) {

					this.y = Math.atan2(m02, m22);
					this.z = Math.atan2(m10, m11);

				} else {

					this.y = Math.atan2(-m20, m00);
					this.z = 0;

				}

				break;

			}

			case RotationOrder.ZXY: {

				this.x = Math.asin(clamp(m21, -1, 1));

				if(Math.abs(m21) < THRESHOLD) {

					this.y = Math.atan2(-m20, m22);
					this.z = Math.atan2(-m01, m11);

				} else {

					this.y = 0;
					this.z = Math.atan2(m10, m00);

				}

				break;

			}

			case RotationOrder.ZYX: {

				this.y = Math.asin(-clamp(m20, -1, 1));

				if(Math.abs(m20) < THRESHOLD) {

					this.x = Math.atan2(m21, m22);
					this.z = Math.atan2(m10, m00);

				} else {

					this.x = 0;
					this.z = Math.atan2(-m01, m11);

				}

				break;

			}

			case RotationOrder.YZX: {

				this.z = Math.asin(clamp(m10, -1, 1));

				if(Math.abs(m10) < THRESHOLD) {

					this.x = Math.atan2(-m12, m11);
					this.y = Math.atan2(-m20, m00);

				} else {

					this.x = 0;
					this.y = Math.atan2(m02, m22);

				}

				break;

			}

			case RotationOrder.XZY: {

				this.z = Math.asin(-clamp(m01, -1, 1));

				if(Math.abs(m01) < THRESHOLD) {

					this.x = Math.atan2(m21, m11);
					this.y = Math.atan2(m02, m00);

				} else {

					this.x = Math.atan2(-m12, m22);
					this.y = 0;

				}

				break;

			}

		}

		this.order = order;

		return this;

	}

	/**
	 * Copies the rotation from a given quaternion.
	 *
	 * @param {Matrix4} q - A quaternion.
	 * @param {RotationOrder} [order] - An override rotation order.
	 * @return {Euler} This set of Euler angles.
	 */

	setFromQuaternion(q, order) {

		m.makeRotationFromQuaternion(q);

		return this.setFromRotationMatrix(m, order);

	}

	/**
	 * Copies the rotation from a given vector.
	 *
	 * @param {Matrix4} v - A vector.
	 * @param {RotationOrder} [order] - A rotation order.
	 * @return {Euler} This set of Euler angles.
	 */

	setFromVector3(v, order = this.order) {

		return this.set(v.x, v.y, v.z, order);

	}

	/**
	 * Reorder the rotation angles.
	 *
	 * WARNING: this operation discards revolution information!
	 *
	 * @param {RotationOrder} newOrder - The new rotation order.
	 * @return {Euler} This set of Euler angles.
	 */

	reorder(newOrder) {

		q.setFromEuler(this);

		return this.setFromQuaternion(q, newOrder);

	}

	/**
	 * Checks if this set of Euler angles equals the given one.
	 *
	 * @param {Euler} e - Euler angles.
	 * @return {Boolean} Whether this set of Euler angles equals the given one.
	 */

	equals(e) {

		return (e.x === this.x && e.y === this.y && e.z === this.z && e.order === this.order);

	}

	/**
	 * The default rotation order.
	 *
	 * @type {RotationOrder}
	 * @final
	 */

	static get defaultOrder() {

		return RotationOrder.XYZ;

	}

}

/**
 * A vector.
 *
 * @type {Vector3}
 * @private
 */

const a = new Vector3();

/**
 * A vector.
 *
 * @type {Vector3}
 * @private
 */

const b = new Vector3();

/**
 * A plane.
 */

class Plane {

	/**
	 * Constructs a new plane.
	 *
	 * @param {Vector3} [normal] - The normal.
	 * @param {Number} [constant] - The constant.
	 */

	constructor(normal = new Vector3(1, 0, 0), constant = 0) {

		/**
		 * The normal.
		 *
		 * @type {Vector3}
		 */

		this.normal = normal;

		/**
		 * The constant.
		 *
		 * @type {Number}
		 */

		this.constant = constant;

	}

	/**
	 * Sets the normal and the constant.
	 *
	 * @param {Vector3} normal - The normal.
	 * @param {Number} constant - The constant.
	 * @return {Plane} This plane.
	 */

	set(normal, constant) {

		this.normal.copy(normal);
		this.constant = constant;

		return this;

	}

	/**
	 * Sets the components of this plane.
	 *
	 * @param {Number} x - The X component of the normal.
	 * @param {Number} y - The Y component of the normal.
	 * @param {Number} z - The Z component of the normal.
	 * @param {Number} w - The constant.
	 * @return {Plane} This plane.
	 */

	setComponents(x, y, z, w) {

		this.normal.set(x, y, z);
		this.constant = w;

		return this;

	}

	/**
	 * Copies the given plane.
	 *
	 * @param {Plane} p - A plane.
	 * @return {Plane} This plane.
	 */

	copy(p) {

		this.normal.copy(p.normal);
		this.constant = p.constant;

		return this;

	}

	/**
	 * Clones this plane.
	 *
	 * @return {Plane} The cloned plane.
	 */

	clone() {

		return new this.constructor().copy(this);

	}

	/**
	 * Sets this plane from a normal and a coplanar point.
	 *
	 * @param {Vector3} n - The normal.
	 * @param {Vector3} p - The coplanar point.
	 * @return {Sphere} This sphere.
	 */

	setFromNormalAndCoplanarPoint(n, p) {

		this.normal.copy(n);
		this.constant = -p.dot(this.normal);

		return this;

	}

	/**
	 * Sets this plane from three distinct coplanar points.
	 *
	 * @param {Vector3} p0 - A coplanar point.
	 * @param {Vector3} p1 - A coplanar point.
	 * @param {Vector3} p2 - A coplanar point.
	 * @return {Plane} This plane.
	 */

	setFromCoplanarPoints(p0, p1, p2) {

		const normal = a.subVectors(p2, p1).cross(b.subVectors(p0, p1)).normalize();

		this.setFromNormalAndCoplanarPoint(normal, a);

		return this;

	}

	/**
	 * Normalizes this plane.
	 *
	 * @return {Plane} This plane.
	 */

	normalize() {

		const inverseNormalLength = 1.0 / this.normal.length();

		this.normal.multiplyScalar(inverseNormalLength);
		this.constant *= inverseNormalLength;

		return this;

	}

	/**
	 * Negates this plane.
	 *
	 * @return {Plane} This plane.
	 */

	negate() {

		this.normal.negate();
		this.constant = -this.constant;

		return this;

	}

	/**
	 * Calculates the distance from this plane to the given point.
	 *
	 * @param {Vector3} p - A point.
	 * @return {Number} The length.
	 */

	distanceToPoint(p) {

		return this.normal.dot(p) + this.constant;

	}

	/**
	 * Calculates the distance from this plane to the given sphere.
	 *
	 * @param {Sphere} s - A sphere.
	 * @return {Number} The length.
	 */

	distanceToSphere(s) {

		return this.distanceToPoint(s.center) - s.radius;

	}

	/**
	 * Projects the given point on this plane.
	 *
	 * @param {Vector3} p - A point.
	 * @param {Vector3} [target] - A target vector. If none is provided, a new one will be created.
	 * @return {Vector3} The projected point.
	 */

	projectPoint(p, target) {

		return target.copy(this.normal).multiplyScalar(-this.distanceToPoint(p)).add(p);

	}

	/**
	 * Calculates a coplanar point and returns it.
	 *
	 * @param {Vector3} [target] - A target vector. If none is provided, a new one will be created.
	 * @return {Vector3} A coplanar plane.
	 */

	coplanarPoint(target) {

		return target.copy(this.normal).multiplyScalar(-this.constant);

	}

	/**
	 * Translates this plane.
	 *
	 * @param {Vector3} offset - An offset.
	 * @return {Plane} This plane.
	 */

	translate(offset) {

		this.constant -= offset.dot(this.normal);

		return this;

	}

	/**
	 * Finds the point of intersection between this plane and a given line.
	 *
	 * @param {Line3} l - A line.
	 * @param {Vector3} [target] - A target vector. If none is provided, a new one will be created.
	 * @return {Vector3} The intersection point.
	 */

	intersectLine(l, target) {

		const direction = l.delta(a);
		const denominator = this.normal.dot(direction);

		if(denominator === 0) {

			// The line is coplanar, return origin.
			if(this.distanceToPoint(l.start) === 0) {

				target.copy(l.start);

			}

		} else {

			const t = -(l.start.dot(this.normal) + this.constant) / denominator;

			if(t >= 0 && t <= 1) {

				target.copy(direction).multiplyScalar(t).add(l.start);

			}

		}

		return target;

	}

	/**
	 * Checks if this plane intersects with the given line.
	 *
	 * @param {Line3} l - A line.
	 * @return {Boolean} Whether this plane intersects with the given line.
	 */

	intersectsLine(l) {

		const startSign = this.distanceToPoint(l.start);
		const endSign = this.distanceToPoint(l.end);

		return ((startSign < 0 && endSign > 0) || (endSign < 0 && startSign > 0));

	}

	/**
	 * Checks if this plane intersects with the given box.
	 *
	 * @param {Box3} b - A box.
	 * @return {Boolean} Whether this plane intersects with the given box.
	 */

	intersectsBox(b) {

		return b.intersectsPlane(this);

	}

	/**
	 * Checks if this plane intersects with the given sphere.
	 *
	 * @param {Sphere} s - A sphere.
	 * @return {Boolean} Whether this plane intersects with the given sphere.
	 */

	intersectsSphere(s) {

		return s.intersectsPlane(this);

	}

	/**
	 * Checks if this plane equals the given one.
	 *
	 * @param {Plane} p - A plane.
	 * @return {Boolean} Whether this plane equals the given one.
	 */

	equals(p) {

		return (p.normal.equals(this.normal) && (p.constant === this.constant));

	}

}

/**
 * A vector.
 *
 * @type {Vector3}
 * @private
 */

const v$3 = new Vector3();

/**
 * A frustum.
 */

class Frustum {

	/**
	 * Constructs a new frustum.
	 *
	 * @param {Plane} [p0] - A plane.
	 * @param {Plane} [p1] - A plane.
	 * @param {Plane} [p2] - A plane.
	 * @param {Plane} [p3] - A plane.
	 * @param {Plane} [p4] - A plane.
	 * @param {Plane} [p5] - A plane.
	 */

	constructor(
		p0 = new Plane(),
		p1 = new Plane(),
		p2 = new Plane(),
		p3 = new Plane(),
		p4 = new Plane(),
		p5 = new Plane()
	) {

		/**
		 * The six planes that form the frustum.
		 *
		 * @type {Plane[]}
		 */

		this.planes = [p0, p1, p2, p3, p4, p5];

	}

	/**
	 * Sets the planes of this frustum.
	 *
	 * @param {Plane} [p0] - A plane.
	 * @param {Plane} [p1] - A plane.
	 * @param {Plane} [p2] - A plane.
	 * @param {Plane} [p3] - A plane.
	 * @param {Plane} [p4] - A plane.
	 * @param {Plane} [p5] - A plane.
	 * @return {Frustum} This frustum.
	 */

	set(p0, p1, p2, p3, p4, p5) {

		const planes = this.planes;

		planes[0].copy(p0);
		planes[1].copy(p1);
		planes[2].copy(p2);
		planes[3].copy(p3);
		planes[4].copy(p4);
		planes[5].copy(p5);

		return this;

	}

	/**
	 * Clones this frustum.
	 *
	 * @return {Frustum} The cloned frustum.
	 */

	clone() {

		return new this.constructor().copy(this);

	}

	/**
	 * Copies a given frustum.
	 *
	 * @param {Frustum} frustum - A frustum.
	 * @return {Frustum} This frustum.
	 */

	copy(frustum) {

		const planes = this.planes;

		let i;

		for(i = 0; i < 6; ++i) {

			planes[i].copy(frustum.planes[i]);

		}

		return this;

	}

	/**
	 * Sets this frustum based on a given projection matrix.
	 *
	 * @param {Matrix4} m - A matrix.
	 * @return {Frustum} This frustum.
	 * @deprecated Use setFromPerspectiveMatrix instead.
	 */

	setFromMatrix(m) {

		return this.setFromProjectionMatrix(m);

	}

	/**
	 * Sets this frustum based on a given projection matrix.
	 *
	 * @param {Matrix4} m - A matrix.
	 * @return {Frustum} This frustum.
	 */

	setFromProjectionMatrix(m) {

		const planes = this.planes;

		const me = m.elements;
		const me0 = me[0], me1 = me[1], me2 = me[2], me3 = me[3];
		const me4 = me[4], me5 = me[5], me6 = me[6], me7 = me[7];
		const me8 = me[8], me9 = me[9], me10 = me[10], me11 = me[11];
		const me12 = me[12], me13 = me[13], me14 = me[14], me15 = me[15];

		planes[0].setComponents(me3 - me0, me7 - me4, me11 - me8, me15 - me12).normalize();
		planes[1].setComponents(me3 + me0, me7 + me4, me11 + me8, me15 + me12).normalize();
		planes[2].setComponents(me3 + me1, me7 + me5, me11 + me9, me15 + me13).normalize();
		planes[3].setComponents(me3 - me1, me7 - me5, me11 - me9, me15 - me13).normalize();
		planes[4].setComponents(me3 - me2, me7 - me6, me11 - me10, me15 - me14).normalize();
		planes[5].setComponents(me3 + me2, me7 + me6, me11 + me10, me15 + me14).normalize();

		return this;

	}

	/**
	 * Checks if this frustum intersects with the given sphere.
	 *
	 * @param {Sphere} sphere - A sphere.
	 * @return {Boolean} Whether this frustum intersects with the sphere.
	 */

	intersectsSphere(sphere) {

		const planes = this.planes;
		const center = sphere.center;
		const negativeRadius = -sphere.radius;

		let result = true;
		let i, d;

		for(i = 0; i < 6; ++i) {

			d = planes[i].distanceToPoint(center);

			if(d < negativeRadius) {

				result = false;
				break;

			}

		}

		return result;

	}

	/**
	 * Checks if this frustum intersects with the given sphere.
	 *
	 * @param {Box3} box - A box.
	 * @return {Boolean} Whether this frustum intersects with the box.
	 */

	intersectsBox(box) {

		const planes = this.planes;
		const min = box.min, max = box.max;

		let i, plane;

		for(i = 0; i < 6; ++i) {

			plane = planes[i];

			// Corner at max distance.
			v$3.x = (plane.normal.x > 0.0) ? max.x : min.x;
			v$3.y = (plane.normal.y > 0.0) ? max.y : min.y;
			v$3.z = (plane.normal.z > 0.0) ? max.z : min.z;

			if(plane.distanceToPoint(v$3) < 0.0) {

				return false;

			}

		}

		return true;

	}

	/**
	 * Checks if this frustum contains the given point.
	 *
	 * @param {Vector3} point - A point.
	 * @return {Boolean} Whether this frustum contains the point.
	 */

	containsPoint(point) {

		const planes = this.planes;

		let result = true;
		let i;

		for(i = 0; i < 6; ++i) {

			if(planes[i].distanceToPoint(point) < 0) {

				result = false;
				break;

			}

		}

		return result;

	}

}

/**
 * A vector.
 *
 * @type {Vector3}
 * @private
 */

const a$1 = new Vector3();

/**
 * A vector.
 *
 * @type {Vector3}
 * @private
 */

const b$1 = new Vector3();

/**
 * A line.
 */

class Line3 {

	/**
	 * Constructs a new line.
	 *
	 * @param {Vector3} [start] - The starting point. If none is provided, a new vector will be created.
	 * @param {Vector3} [end] - The ending point. If none is provided, a new vector will be created.
	 */

	constructor(start = new Vector3(), end = new Vector3()) {

		/**
		 * The starting point.
		 *
		 * @type {Vector3}
		 */

		this.start = start;

		/**
		 * The ending point.
		 *
		 * @type {Vector3}
		 */

		this.end = end;

	}

	/**
	 * Sets the starting and ending point of this line.
	 *
	 * @param {Vector3} start - The starting point.
	 * @param {Vector3} end - The ending point.
	 * @return {Line3} This line.
	 */

	set(start, end) {

		this.start.copy(start);
		this.end.copy(end);

		return this;

	}

	/**
	 * Copies the values of the given line.
	 *
	 * @param {Line3} l - A line.
	 * @return {Line3} This line.
	 */

	copy(l) {

		this.start.copy(l.start);
		this.end.copy(l.end);

		return this;

	}

	/**
	 * Clones this line.
	 *
	 * @return {Line3} The cloned line.
	 */

	clone() {

		return new this.constructor().copy(this);

	}

	/**
	 * Calculates the center of this line.
	 *
	 * @param {Vector3} [target] - A target vector. If none is provided, a new one will be created.
	 * @return {Vector3} The center of this line.
	 */

	getCenter(target = new Vector3()) {

		return target.addVectors(this.start, this.end).multiplyScalar(0.5);

	}

	/**
	 * Calculates the delta vector of this line.
	 *
	 * @param {Vector3} [target] - A target vector. If none is provided, a new one will be created.
	 * @return {Vector3} The delta vector of this line.
	 */

	delta(target = new Vector3()) {

		return target.subVectors(this.end, this.start);

	}

	/**
	 * Calculates the squared length of this line.
	 *
	 * @return {Vector3} The squared length.
	 */

	lengthSquared() {

		return this.start.distanceToSquared(this.end);

	}

	/**
	 * Calculates the length of this line.
	 *
	 * @return {Vector3} The length.
	 */

	length() {

		return this.start.distanceTo(this.end);

	}

	/**
	 * Adjusts this line to point in the given direction.
	 *
	 * @param {Vector3} d - The direction.
	 * @param {Vector3} [target] - A target vector. If none is provided, a new one will be created.
	 * @return {Vector3} The length.
	 */

	at(d, target) {

		return this.delta(target).multiplyScalar(d).add(this.start);

	}

	/**
	 * Returns a point parameter based on the closest point as projected on the line segement.
	 *
	 * @private
	 * @param {Vector3} p - A point.
	 * @param {Boolean} clampToLine - Whether the point should be clamped to the line.
	 * @return {Vector3} The parameter.
	 */

	closestPointToPointParameter(p, clampToLine) {

		a$1.subVectors(p, this.start);
		b$1.subVectors(this.end, this.start);

		const bb = b$1.dot(b$1);
		const ba = b$1.dot(a$1);

		const t = clampToLine ? Math.min(Math.max(ba / bb, 0), 1) : ba / bb;

		return t;

	}

	/**
	 * Returns the closest point on the line.
	 *
	 * @param {Vector3} p - A point.
	 * @param {Boolean} [clampToLine=false] - Whether the point should be clamped to the line.
	 * @param {Vector3} [target] - A target vector. If none is provided, a new one will be created.
	 * @return {Vector3} The parameter.
	 */

	closestPointToPoint(p, clampToLine = false, target = new Vector3()) {

		const t = this.closestPointToPointParameter(p, clampToLine);

		return this.delta(target).multiplyScalar(t).add(this.start);

	}

	/**
	 * Checks if this line equals the given one.
	 *
	 * @param {Line3} l - A line.
	 * @return {Boolean} Whether the lines are equal.
	 */

	equals(l) {

		return l.start.equals(this.start) && l.end.equals(this.end);

	}

}

/**
 * A vector.
 *
 * @type {Vector3}
 * @private
 */

const a$2 = new Vector3();

/**
 * A vector.
 *
 * @type {Vector3}
 * @private
 */

const b$2 = new Vector3();

/**
 * A vector.
 *
 * @type {Vector3}
 * @private
 */

const c = new Vector3();

/**
 * A 4x4 matrix.
 */

class Matrix4 {

	/**
	 * Constructs a new matrix.
	 */

	constructor() {

		/**
		 * The matrix elements.
		 *
		 * @type {Float32Array}
		 */

		this.elements = new Float32Array([

			1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, 1, 0,
			0, 0, 0, 1

		]);

	}

	/**
	 * Sets the values of this matrix.
	 *
	 * @param {Number} n00 - The value of the first row, first column.
	 * @param {Number} n01 - The value of the first row, second column.
	 * @param {Number} n02 - The value of the first row, third column.
	 * @param {Number} n03 - The value of the first row, fourth column.
	 * @param {Number} n10 - The value of the second row, first column.
	 * @param {Number} n11 - The value of the second row, second column.
	 * @param {Number} n12 - The value of the second row, third column.
	 * @param {Number} n13 - The value of the second row, fourth column.
	 * @param {Number} n20 - The value of the third row, first column.
	 * @param {Number} n21 - The value of the third row, second column.
	 * @param {Number} n22 - The value of the third row, third column.
	 * @param {Number} n23 - The value of the third row, fourth column.
	 * @param {Number} n30 - The value of the fourth row, first column.
	 * @param {Number} n31 - The value of the fourth row, second column.
	 * @param {Number} n32 - The value of the fourth row, third column.
	 * @param {Number} n33 - The value of the fourth row, fourth column.
	 * @return {Matrix4} This matrix.
	 */

	set(n00, n01, n02, n03, n10, n11, n12, n13, n20, n21, n22, n23, n30, n31, n32, n33) {

		const te = this.elements;

		te[0] = n00; te[4] = n01; te[8] = n02; te[12] = n03;
		te[1] = n10; te[5] = n11; te[9] = n12; te[13] = n13;
		te[2] = n20; te[6] = n21; te[10] = n22; te[14] = n23;
		te[3] = n30; te[7] = n31; te[11] = n32; te[15] = n33;

		return this;

	}

	/**
	 * Sets this matrix to the identity matrix.
	 *
	 * @return {Matrix4} This matrix.
	 */

	identity() {

		this.set(

			1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, 1, 0,
			0, 0, 0, 1

		);

		return this;

	}

	/**
	 * Copies the values of a given matrix.
	 *
	 * @param {Matrix4} matrix - A matrix.
	 * @return {Matrix4} This matrix.
	 */

	copy(matrix) {

		const me = matrix.elements;
		const te = this.elements;

		te[0] = me[0]; te[1] = me[1]; te[2] = me[2]; te[3] = me[3];
		te[4] = me[4]; te[5] = me[5]; te[6] = me[6]; te[7] = me[7];
		te[8] = me[8]; te[9] = me[9]; te[10] = me[10]; te[11] = me[11];
		te[12] = me[12]; te[13] = me[13]; te[14] = me[14]; te[15] = me[15];

		return this;

	}

	/**
	 * Clones this matrix.
	 *
	 * @return {Matrix4} A clone of this matrix.
	 */

	clone() {

		return new this.constructor().fromArray(this.elements);

	}

	/**
	 * Copies the values of a given array.
	 *
	 * @param {Number[]} array - An array.
	 * @param {Number} [offset=0] - An offset into the array.
	 * @return {Matrix4} This matrix.
	 */

	fromArray(array, offset = 0) {

		const te = this.elements;

		let i;

		for(i = 0; i < 16; ++i) {

			te[i] = array[i + offset];

		}

		return this;

	}

	/**
	 * Stores this matrix in an array.
	 *
	 * @param {Number[]} [array] - A target array.
	 * @param {Number} [offset=0] - An offset into the array.
	 * @return {Number[]} The array.
	 */

	toArray(array = [], offset = 0) {

		const te = this.elements;

		let i;

		for(i = 0; i < 16; ++i) {

			array[i + offset] = te[i];

		}

		return array;

	}

	/**
	 * Returns the largest scale.
	 *
	 * @return {Number} The largest scale of the three axes.
	 */

	getMaxScaleOnAxis() {

		const te = this.elements;

		const scaleXSq = te[0] * te[0] + te[1] * te[1] + te[2] * te[2];
		const scaleYSq = te[4] * te[4] + te[5] * te[5] + te[6] * te[6];
		const scaleZSq = te[8] * te[8] + te[9] * te[9] + te[10] * te[10];

		return Math.sqrt(Math.max(scaleXSq, scaleYSq, scaleZSq));

	}

	/**
	 * Copies the position values of a given matrix.
	 *
	 * @param {Matrix4} matrix - A matrix.
	 * @return {Matrix4} This matrix.
	 */

	copyPosition(matrix) {

		const te = this.elements;
		const me = matrix.elements;

		te[12] = me[12];
		te[13] = me[13];
		te[14] = me[14];

		return this;

	}

	/**
	 * Sets the position values of this matrix.
	 *
	 * @param {Vector3} p - A position.
	 * @return {Matrix4} This matrix.
	 */

	setPosition(p) {

		const te = this.elements;

		te[12] = p.x;
		te[13] = p.y;
		te[14] = p.z;

		return this;

	}

	/**
	 * Extracts the basis from this matrix.
	 *
	 * @param {Vector3} xAxis - A vector to store the X-axis column in.
	 * @param {Vector3} yAxis - A vector to store the Y-axis column in.
	 * @param {Vector3} zAxis - A vector to store the Z-axis column in.
	 * @return {Matrix4} This matrix.
	 */

	extractBasis(xAxis, yAxis, zAxis) {

		xAxis.setFromMatrixColumn(this, 0);
		yAxis.setFromMatrixColumn(this, 1);
		zAxis.setFromMatrixColumn(this, 2);

		return this;

	}

	/**
	 * Sets the basis of this matrix.
	 *
	 * @param {Vector3} xAxis - The X-axis.
	 * @param {Vector3} yAxis - The Y-axis.
	 * @param {Vector3} zAxis - The Z-axis.
	 * @return {Matrix4} This matrix.
	 */

	makeBasis(xAxis, yAxis, zAxis) {

		this.set(

			xAxis.x, yAxis.x, zAxis.x, 0,
			xAxis.y, yAxis.y, zAxis.y, 0,
			xAxis.z, yAxis.z, zAxis.z, 0,
			0, 0, 0, 1

		);

		return this;

	}

	/**
	 * Extracts the rotation from a given matrix.
	 *
	 * This method does not support reflection matrices.
	 *
	 * @param {Matrix4} m - A matrix.
	 * @return {Matrix4} This matrix.
	 */

	extractRotation(m) {

		const te = this.elements;
		const me = m.elements;

		const scaleX = 1.0 / a$2.setFromMatrixColumn(m, 0).length();
		const scaleY = 1.0 / a$2.setFromMatrixColumn(m, 1).length();
		const scaleZ = 1.0 / a$2.setFromMatrixColumn(m, 2).length();

		te[0] = me[0] * scaleX;
		te[1] = me[1] * scaleX;
		te[2] = me[2] * scaleX;
		te[3] = 0;

		te[4] = me[4] * scaleY;
		te[5] = me[5] * scaleY;
		te[6] = me[6] * scaleY;
		te[7] = 0;

		te[8] = me[8] * scaleZ;
		te[9] = me[9] * scaleZ;
		te[10] = me[10] * scaleZ;
		te[11] = 0;

		te[12] = 0;
		te[13] = 0;
		te[14] = 0;
		te[15] = 1;

		return this;

	}

	/**
	 * Sets the matrix rotation based on the given Euler angles.
	 *
	 * @param {Euler} euler - The euler angles.
	 * @return {Matrix4} This matrix.
	 */

	makeRotationFromEuler(euler) {

		const te = this.elements;

		const x = euler.x;
		const y = euler.y;
		const z = euler.z;

		const a = Math.cos(x), b = Math.sin(x);
		const c = Math.cos(y), d = Math.sin(y);
		const e = Math.cos(z), f = Math.sin(z);

		let ae, af, be, bf;
		let ce, cf, de, df;
		let ac, ad, bc, bd;

		switch(euler.order) {

			case RotationOrder.XYZ: {

				ae = a * e, af = a * f, be = b * e, bf = b * f;

				te[0] = c * e;
				te[4] = -c * f;
				te[8] = d;

				te[1] = af + be * d;
				te[5] = ae - bf * d;
				te[9] = -b * c;

				te[2] = bf - ae * d;
				te[6] = be + af * d;
				te[10] = a * c;

				break;

			}

			case RotationOrder.YXZ: {

				ce = c * e, cf = c * f, de = d * e, df = d * f;

				te[0] = ce + df * b;
				te[4] = de * b - cf;
				te[8] = a * d;

				te[1] = a * f;
				te[5] = a * e;
				te[9] = -b;

				te[2] = cf * b - de;
				te[6] = df + ce * b;
				te[10] = a * c;

				break;

			}

			case RotationOrder.ZXY: {

				ce = c * e, cf = c * f, de = d * e, df = d * f;

				te[0] = ce - df * b;
				te[4] = -a * f;
				te[8] = de + cf * b;

				te[1] = cf + de * b;
				te[5] = a * e;
				te[9] = df - ce * b;

				te[2] = -a * d;
				te[6] = b;
				te[10] = a * c;

				break;

			}

			case RotationOrder.ZYX: {

				ae = a * e, af = a * f, be = b * e, bf = b * f;

				te[0] = c * e;
				te[4] = be * d - af;
				te[8] = ae * d + bf;

				te[1] = c * f;
				te[5] = bf * d + ae;
				te[9] = af * d - be;

				te[2] = -d;
				te[6] = b * c;
				te[10] = a * c;

				break;

			}

			case RotationOrder.YZX: {

				ac = a * c, ad = a * d, bc = b * c, bd = b * d;

				te[0] = c * e;
				te[4] = bd - ac * f;
				te[8] = bc * f + ad;

				te[1] = f;
				te[5] = a * e;
				te[9] = -b * e;

				te[2] = -d * e;
				te[6] = ad * f + bc;
				te[10] = ac - bd * f;

				break;

			}

			case RotationOrder.XZY: {

				ac = a * c, ad = a * d, bc = b * c, bd = b * d;

				te[0] = c * e;
				te[4] = -f;
				te[8] = d * e;

				te[1] = ac * f + bd;
				te[5] = a * e;
				te[9] = ad * f - bc;

				te[2] = bc * f - ad;
				te[6] = b * e;
				te[10] = bd * f + ac;

				break;

			}

		}

		// Bottom row.
		te[3] = 0;
		te[7] = 0;
		te[11] = 0;

		// Last column.
		te[12] = 0;
		te[13] = 0;
		te[14] = 0;
		te[15] = 1;

		return this;

	}

	/**
	 * Sets the matrix rotation based on the given quaternion.
	 *
	 * @param {Quaternion} q - The quaternion.
	 * @return {Matrix4} This matrix.
	 */

	makeRotationFromQuaternion(q) {

		return this.compose(a$2.set(0, 0, 0), q, b$2.set(1, 1, 1));

	}

	/**
	 * Creates a rotation that looks at the given target.
	 *
	 * @param {Vector3} eye - The position of the eye.
	 * @param {Vector3} target - The target to look at.
	 * @param {Vector3} up - The up vector.
	 * @return {Matrix4} This matrix.
	 */

	lookAt(eye, target, up) {

		const te = this.elements;
		const x = a$2, y = b$2, z = c;

		z.subVectors(eye, target);

		if(z.lengthSquared() === 0) {

			// Eye and target are at the same position.
			z.z = 1;

		}

		z.normalize();
		x.crossVectors(up, z);

		if(x.lengthSquared() === 0) {

			// Up and z are parallel.
			if(Math.abs(up.z) === 1) {

				z.x += 1e-4;

			} else {

				z.z += 1e-4;

			}

			z.normalize();
			x.crossVectors(up, z);

		}

		x.normalize();
		y.crossVectors(z, x);

		te[0] = x.x; te[4] = y.x; te[8] = z.x;
		te[1] = x.y; te[5] = y.y; te[9] = z.y;
		te[2] = x.z; te[6] = y.z; te[10] = z.z;

		return this;

	}

	/**
	 * Sets this matrix to the product of the given matrices.
	 *
	 * @param {Matrix4} a - A matrix.
	 * @param {Matrix4} b - A matrix.
	 * @return {Matrix4} This matrix.
	 */

	multiplyMatrices(a, b) {

		const te = this.elements;
		const ae = a.elements;
		const be = b.elements;

		const a00 = ae[0], a01 = ae[4], a02 = ae[8], a03 = ae[12];
		const a10 = ae[1], a11 = ae[5], a12 = ae[9], a13 = ae[13];
		const a20 = ae[2], a21 = ae[6], a22 = ae[10], a23 = ae[14];
		const a30 = ae[3], a31 = ae[7], a32 = ae[11], a33 = ae[15];

		const b00 = be[0], b01 = be[4], b02 = be[8], b03 = be[12];
		const b10 = be[1], b11 = be[5], b12 = be[9], b13 = be[13];
		const b20 = be[2], b21 = be[6], b22 = be[10], b23 = be[14];
		const b30 = be[3], b31 = be[7], b32 = be[11], b33 = be[15];

		te[0] = a00 * b00 + a01 * b10 + a02 * b20 + a03 * b30;
		te[4] = a00 * b01 + a01 * b11 + a02 * b21 + a03 * b31;
		te[8] = a00 * b02 + a01 * b12 + a02 * b22 + a03 * b32;
		te[12] = a00 * b03 + a01 * b13 + a02 * b23 + a03 * b33;

		te[1] = a10 * b00 + a11 * b10 + a12 * b20 + a13 * b30;
		te[5] = a10 * b01 + a11 * b11 + a12 * b21 + a13 * b31;
		te[9] = a10 * b02 + a11 * b12 + a12 * b22 + a13 * b32;
		te[13] = a10 * b03 + a11 * b13 + a12 * b23 + a13 * b33;

		te[2] = a20 * b00 + a21 * b10 + a22 * b20 + a23 * b30;
		te[6] = a20 * b01 + a21 * b11 + a22 * b21 + a23 * b31;
		te[10] = a20 * b02 + a21 * b12 + a22 * b22 + a23 * b32;
		te[14] = a20 * b03 + a21 * b13 + a22 * b23 + a23 * b33;

		te[3] = a30 * b00 + a31 * b10 + a32 * b20 + a33 * b30;
		te[7] = a30 * b01 + a31 * b11 + a32 * b21 + a33 * b31;
		te[11] = a30 * b02 + a31 * b12 + a32 * b22 + a33 * b32;
		te[15] = a30 * b03 + a31 * b13 + a32 * b23 + a33 * b33;

		return this;

	}

	/**
	 * Multiplies this matrix with the given one.
	 *
	 * @param {Matrix4} m - A matrix.
	 * @return {Matrix4} This matrix.
	 */

	multiply(m) {

		return this.multiplyMatrices(this, m);

	}

	/**
	 * Multiplies a given matrix with this one.
	 *
	 * @param {Matrix4} m - A matrix.
	 * @return {Matrix4} This matrix.
	 */

	premultiply(m) {

		return this.multiplyMatrices(m, this);

	}

	/**
	 * Multiplies this matrix with a given scalar.
	 *
	 * @param {Number} s - A scalar.
	 * @return {Matrix4} This matrix.
	 */

	multiplyScalar(s) {

		const te = this.elements;

		te[0] *= s; te[4] *= s; te[8] *= s; te[12] *= s;
		te[1] *= s; te[5] *= s; te[9] *= s; te[13] *= s;
		te[2] *= s; te[6] *= s; te[10] *= s; te[14] *= s;
		te[3] *= s; te[7] *= s; te[11] *= s; te[15] *= s;

		return this;

	}

	/**
	 * Calculates the determinant of this matrix.
	 *
	 * For more details see:
	 *  http://www.euclideanspace.com/maths/algebra/matrix/functions/inverse/fourD/index.htm
	 *
	 * @return {Number} The determinant.
	 */

	determinant() {

		const te = this.elements;

		const n00 = te[0], n01 = te[4], n02 = te[8], n03 = te[12];
		const n10 = te[1], n11 = te[5], n12 = te[9], n13 = te[13];
		const n20 = te[2], n21 = te[6], n22 = te[10], n23 = te[14];
		const n30 = te[3], n31 = te[7], n32 = te[11], n33 = te[15];

		const n00n11 = n00 * n11, n00n12 = n00 * n12, n00n13 = n00 * n13;
		const n01n10 = n01 * n10, n01n12 = n01 * n12, n01n13 = n01 * n13;
		const n02n10 = n02 * n10, n02n11 = n02 * n11, n02n13 = n02 * n13;
		const n03n10 = n03 * n10, n03n11 = n03 * n11, n03n12 = n03 * n12;

		return (

			n30 * (
				n03n12 * n21 -
				n02n13 * n21 -
				n03n11 * n22 +
				n01n13 * n22 +
				n02n11 * n23 -
				n01n12 * n23
			) +

			n31 * (
				n00n12 * n23 -
				n00n13 * n22 +
				n03n10 * n22 -
				n02n10 * n23 +
				n02n13 * n20 -
				n03n12 * n20
			) +

			n32 * (
				n00n13 * n21 -
				n00n11 * n23 -
				n03n10 * n21 +
				n01n10 * n23 +
				n03n11 * n20 -
				n01n13 * n20
			) +

			n33 * (
				-n02n11 * n20 -
				n00n12 * n21 +
				n00n11 * n22 +
				n02n10 * n21 -
				n01n10 * n22 +
				n01n12 * n20
			)

		);

	}

	/**
	 * Inverts the given matrix and stores the result in this matrix.
	 *
	 * For details see:
	 *  http://www.euclideanspace.com/maths/algebra/matrix/functions/inverse/fourD/index.htm
	 *
	 * @param {Matrix4} matrix - The matrix that should be inverted.
	 * @return {Matrix4} This matrix.
	 */

	getInverse(matrix) {

		const te = this.elements;
		const me = matrix.elements;

		const n00 = me[0], n10 = me[1], n20 = me[2], n30 = me[3];
		const n01 = me[4], n11 = me[5], n21 = me[6], n31 = me[7];
		const n02 = me[8], n12 = me[9], n22 = me[10], n32 = me[11];
		const n03 = me[12], n13 = me[13], n23 = me[14], n33 = me[15];

		const t00 = n12 * n23 * n31 - n13 * n22 * n31 + n13 * n21 * n32 - n11 * n23 * n32 - n12 * n21 * n33 + n11 * n22 * n33;
		const t01 = n03 * n22 * n31 - n02 * n23 * n31 - n03 * n21 * n32 + n01 * n23 * n32 + n02 * n21 * n33 - n01 * n22 * n33;
		const t02 = n02 * n13 * n31 - n03 * n12 * n31 + n03 * n11 * n32 - n01 * n13 * n32 - n02 * n11 * n33 + n01 * n12 * n33;
		const t03 = n03 * n12 * n21 - n02 * n13 * n21 - n03 * n11 * n22 + n01 * n13 * n22 + n02 * n11 * n23 - n01 * n12 * n23;

		const det = n00 * t00 + n10 * t01 + n20 * t02 + n30 * t03;

		let invDet;

		if(det !== 0) {

			invDet = 1.0 / det;

			te[0] = t00 * invDet;
			te[1] = (n13 * n22 * n30 - n12 * n23 * n30 - n13 * n20 * n32 + n10 * n23 * n32 + n12 * n20 * n33 - n10 * n22 * n33) * invDet;
			te[2] = (n11 * n23 * n30 - n13 * n21 * n30 + n13 * n20 * n31 - n10 * n23 * n31 - n11 * n20 * n33 + n10 * n21 * n33) * invDet;
			te[3] = (n12 * n21 * n30 - n11 * n22 * n30 - n12 * n20 * n31 + n10 * n22 * n31 + n11 * n20 * n32 - n10 * n21 * n32) * invDet;

			te[4] = t01 * invDet;
			te[5] = (n02 * n23 * n30 - n03 * n22 * n30 + n03 * n20 * n32 - n00 * n23 * n32 - n02 * n20 * n33 + n00 * n22 * n33) * invDet;
			te[6] = (n03 * n21 * n30 - n01 * n23 * n30 - n03 * n20 * n31 + n00 * n23 * n31 + n01 * n20 * n33 - n00 * n21 * n33) * invDet;
			te[7] = (n01 * n22 * n30 - n02 * n21 * n30 + n02 * n20 * n31 - n00 * n22 * n31 - n01 * n20 * n32 + n00 * n21 * n32) * invDet;

			te[8] = t02 * invDet;
			te[9] = (n03 * n12 * n30 - n02 * n13 * n30 - n03 * n10 * n32 + n00 * n13 * n32 + n02 * n10 * n33 - n00 * n12 * n33) * invDet;
			te[10] = (n01 * n13 * n30 - n03 * n11 * n30 + n03 * n10 * n31 - n00 * n13 * n31 - n01 * n10 * n33 + n00 * n11 * n33) * invDet;
			te[11] = (n02 * n11 * n30 - n01 * n12 * n30 - n02 * n10 * n31 + n00 * n12 * n31 + n01 * n10 * n32 - n00 * n11 * n32) * invDet;

			te[12] = t03 * invDet;
			te[13] = (n02 * n13 * n20 - n03 * n12 * n20 + n03 * n10 * n22 - n00 * n13 * n22 - n02 * n10 * n23 + n00 * n12 * n23) * invDet;
			te[14] = (n03 * n11 * n20 - n01 * n13 * n20 - n03 * n10 * n21 + n00 * n13 * n21 + n01 * n10 * n23 - n00 * n11 * n23) * invDet;
			te[15] = (n01 * n12 * n20 - n02 * n11 * n20 + n02 * n10 * n21 - n00 * n12 * n21 - n01 * n10 * n22 + n00 * n11 * n22) * invDet;

		} else {

			console.error("Can't invert matrix, determinant is zero", matrix);

			this.identity();

		}

		return this;

	}

	/**
	 * Transposes this matrix.
	 *
	 * @return {Matrix4} This matrix.
	 */

	transpose() {

		const te = this.elements;

		let t;

		t = te[1]; te[1] = te[4]; te[4] = t;
		t = te[2]; te[2] = te[8]; te[8] = t;
		t = te[6]; te[6] = te[9]; te[9] = t;

		t = te[3]; te[3] = te[12]; te[12] = t;
		t = te[7]; te[7] = te[13]; te[13] = t;
		t = te[11]; te[11] = te[14]; te[14] = t;

		return this;

	}

	/**
	 * Scales this matrix.
	 *
	 * @param {Number} sx - The X scale.
	 * @param {Number} sy - The Y scale.
	 * @param {Number} sz - The Z scale.
	 * @return {Matrix4} This matrix.
	 */

	scale(sx, sy, sz) {

		const te = this.elements;

		te[0] *= sx; te[4] *= sy; te[8] *= sz;
		te[1] *= sx; te[5] *= sy; te[9] *= sz;
		te[2] *= sx; te[6] *= sy; te[10] *= sz;
		te[3] *= sx; te[7] *= sy; te[11] *= sz;

		return this;

	}

	/**
	 * Makes this matrix a scale matrix.
	 *
	 * @param {Number} x - The X scale.
	 * @param {Number} y - The Y scale.
	 * @param {Number} z - The Z scale.
	 * @return {Matrix4} This matrix.
	 */

	makeScale(x, y, z) {

		this.set(

			x, 0, 0, 0,
			0, y, 0, 0,
			0, 0, z, 0,
			0, 0, 0, 1

		);

		return this;

	}

	/**
	 * Makes this matrix a translation matrix.
	 *
	 * @param {Number} x - The X offset.
	 * @param {Number} y - The Y offset.
	 * @param {Number} z - The Z offset.
	 * @return {Matrix4} This matrix.
	 */

	makeTranslation(x, y, z) {

		this.set(

			1, 0, 0, x,
			0, 1, 0, y,
			0, 0, 1, z,
			0, 0, 0, 1

		);

		return this;

	}

	/**
	 * Makes this matrix a rotation matrix.
	 *
	 * @param {Number} theta - The angle in radians.
	 * @return {Matrix4} This matrix.
	 */

	makeRotationX(theta) {

		const c = Math.cos(theta), s = Math.sin(theta);

		this.set(

			1, 0, 0, 0,
			0, c, -s, 0,
			0, s, c, 0,
			0, 0, 0, 1

		);

		return this;

	}

	/**
	 * Makes this matrix a rotation matrix with respect to the Y-axis.
	 *
	 * @param {Number} theta - The angle in radians.
	 * @return {Matrix4} This matrix.
	 */

	makeRotationY(theta) {

		const c = Math.cos(theta), s = Math.sin(theta);

		this.set(

			c, 0, s, 0,
			0, 1, 0, 0,
			-s, 0, c, 0,
			0, 0, 0, 1

		);

		return this;

	}

	/**
	 * Makes this matrix a rotation matrix with respect to the Z-axis.
	 *
	 * @param {Number} theta - The angle in radians.
	 * @return {Matrix4} This matrix.
	 */

	makeRotationZ(theta) {

		const c = Math.cos(theta), s = Math.sin(theta);

		this.set(

			c, -s, 0, 0,
			s, c, 0, 0,
			0, 0, 1, 0,
			0, 0, 0, 1

		);

		return this;

	}

	/**
	 * Makes this matrix a translation matrix with respect to a specific axis.
	 *
	 * For mor einformation see:
	 *  http://www.gamedev.net/reference/articles/article1199.asp
	 *
	 * @param {Vector3} axis - The axis. Assumed to be normalized.
	 * @param {Number} angle - The angle in radians.
	 * @return {Matrix4} This matrix.
	 */

	makeRotationAxis(axis, angle) {

		const c = Math.cos(angle);
		const s = Math.sin(angle);

		const t = 1.0 - c;

		const x = axis.x, y = axis.y, z = axis.z;
		const tx = t * x, ty = t * y;

		this.set(

			tx * x + c, tx * y - s * z, tx * z + s * y, 0,
			tx * y + s * z, ty * y + c, ty * z - s * x, 0,
			tx * z - s * y, ty * z + s * x, t * z * z + c, 0,
			0, 0, 0, 1

		);

		return this;

	}

	/**
	 * Makes this matrix a shear matrix.
	 *
	 * @param {Number} x - The X shear value.
	 * @param {Number} y - The Y shear value.
	 * @param {Number} z - The Z shear value.
	 * @return {Matrix4} This matrix.
	 */

	makeShear(x, y, z) {

		this.set(

			1, y, z, 0,
			x, 1, z, 0,
			x, y, 1, 0,
			0, 0, 0, 1

		);

		return this;

	}

	/**
	 * Sets this matrix based on the given position, rotation and scale.
	 *
	 * @param {Vector3} position - The position.
	 * @param {Quaternion} quaternion - The rotation.
	 * @param {Vector3} scale - The scale.
	 * @return {Matrix4} This matrix.
	 */

	compose(position, quaternion, scale) {

		const te = this.elements;

		const x = quaternion.x, y = quaternion.y, z = quaternion.z, w = quaternion.w;
		const x2 = x + x,	y2 = y + y, z2 = z + z;
		const xx = x * x2, xy = x * y2, xz = x * z2;
		const yy = y * y2, yz = y * z2, zz = z * z2;
		const wx = w * x2, wy = w * y2, wz = w * z2;

		const sx = scale.x, sy = scale.y, sz = scale.z;

		te[0] = (1 - (yy + zz)) * sx;
		te[1] = (xy + wz) * sx;
		te[2] = (xz - wy) * sx;
		te[3] = 0;

		te[4] = (xy - wz) * sy;
		te[5] = (1 - (xx + zz)) * sy;
		te[6] = (yz + wx) * sy;
		te[7] = 0;

		te[8] = (xz + wy) * sz;
		te[9] = (yz - wx) * sz;
		te[10] = (1 - (xx + yy)) * sz;
		te[11] = 0;

		te[12] = position.x;
		te[13] = position.y;
		te[14] = position.z;
		te[15] = 1;

		return this;

	}

	/**
	 * Decomposes this matrix into a position, rotation and scale vector.
	 *
	 * @param {Vector3} position - The target position.
	 * @param {Quaternion} quaternion - The target rotation.
	 * @param {Vector3} scale - The target scale.
	 * @return {Matrix4} This matrix.
	 */

	decompose(position, quaternion, scale) {

		const te = this.elements;

		const n00 = te[0], n10 = te[1], n20 = te[2];
		const n01 = te[4], n11 = te[5], n21 = te[6];
		const n02 = te[8], n12 = te[9], n22 = te[10];

		const det = this.determinant();

		// If the determinant is negative, one scale must be inverted.
		const sx = a$2.set(n00, n10, n20).length() * ((det < 0) ? -1 : 1);
		const sy = a$2.set(n01, n11, n21).length();
		const sz = a$2.set(n02, n12, n22).length();

		const invSX = 1.0 / sx;
		const invSY = 1.0 / sy;
		const invSZ = 1.0 / sz;

		// Export the position.
		position.x = te[12];
		position.y = te[13];
		position.z = te[14];

		// Scale the rotation part.
		te[0] *= invSX; te[1] *= invSX; te[2] *= invSX;
		te[4] *= invSY; te[5] *= invSY; te[6] *= invSY;
		te[8] *= invSZ; te[9] *= invSZ; te[10] *= invSZ;

		// Export the rotation.
		quaternion.setFromRotationMatrix(this);

		// Restore the original values.
		te[0] = n00; te[1] = n10; te[2] = n20;
		te[4] = n01; te[5] = n11; te[6] = n21;
		te[8] = n02; te[9] = n12; te[10] = n22;

		// Export the scale.
		scale.x = sx;
		scale.y = sy;
		scale.z = sz;

		return this;

	}

	/**
	 * Creates a perspective matrix.
	 *
	 * @param {Number} left - The distance to the left plane.
	 * @param {Number} right - The distance to the right plane.
	 * @param {Number} top - The distance to the top plane.
	 * @param {Number} bottom - The distance to the bottom plane.
	 * @param {Number} near - The distance to the near plane.
	 * @param {Number} far - The distance to the far plane.
	 * @return {Matrix4} This matrix.
	 */

	makePerspective(left, right, top, bottom, near, far) {

		const te = this.elements;
		const x = 2 * near / (right - left);
		const y = 2 * near / (top - bottom);

		const a = (right + left) / (right - left);
		const b = (top + bottom) / (top - bottom);
		const c = -(far + near) / (far - near);
		const d = -2 * far * near / (far - near);

		te[0] = x; te[4] = 0; te[8] = a; te[12] = 0;
		te[1] = 0; te[5] = y; te[9] = b; te[13] = 0;
		te[2] = 0; te[6] = 0; te[10] = c; te[14] = d;
		te[3] = 0; te[7] = 0; te[11] = -1; te[15] = 0;

		return this;

	}

	/**
	 * Creates an orthographic matrix.
	 *
	 * @param {Number} left - The distance to the left plane.
	 * @param {Number} right - The distance to the right plane.
	 * @param {Number} top - The distance to the top plane.
	 * @param {Number} bottom - The distance to the bottom plane.
	 * @param {Number} near - The distance to the near plane.
	 * @param {Number} far - The distance to the far plane.
	 * @return {Matrix4} This matrix.
	 */

	makeOrthographic(left, right, top, bottom, near, far) {

		const te = this.elements;
		const w = 1.0 / (right - left);
		const h = 1.0 / (top - bottom);
		const p = 1.0 / (far - near);

		const x = (right + left) * w;
		const y = (top + bottom) * h;
		const z = (far + near) * p;

		te[0] = 2 * w; te[4] = 0; te[8] = 0; te[12] = -x;
		te[1] = 0; te[5] = 2 * h; te[9] = 0; te[13] = -y;
		te[2] = 0; te[6] = 0; te[10] = -2 * p; te[14] = -z;
		te[3] = 0; te[7] = 0; te[11] = 0; te[15] = 1;

		return this;

	}

	/**
	 * Checks if this matrix equals the given one.
	 *
	 * @param {Matrix4} m - A matrix.
	 * @return {Boolean} Whether the matrix are equal.
	 */

	equals(m) {

		const te = this.elements;
		const me = m.elements;

		let result = true;
		let i;

		for(i = 0; result && i < 16; ++i) {

			if(te[i] !== me[i]) {

				result = false;

			}

		}

		return result;

	}

}

/**
 * A list of vectors.
 *
 * @type {Vector3[]}
 * @private
 */

const v$4 = [
	new Vector3(),
	new Vector3(),
	new Vector3(),
	new Vector3()
];

/**
 * A ray.
 */

class Ray {

	/**
	 * Constructs a new ray.
	 *
	 * @param {Vector3} [origin] - The origin.
	 * @param {Vector3} [direction] - The direction.
	 */

	constructor(origin = new Vector3(), direction = new Vector3(0, 0, -1)) {

		/**
		 * The origin.
		 *
		 * @type {Vector3}
		 */

		this.origin = origin;

		/**
		 * The direction.
		 *
		 * @type {Vector3}
		 */

		this.direction = direction;

	}

	/**
	 * Sets the origin and the direction.
	 *
	 * @param {Vector3} origin - The origin.
	 * @param {Vector3} direction - The direction. Should be normalized.
	 * @return {Ray} This ray.
	 */

	set(origin, direction) {

		this.origin.copy(origin);
		this.direction.copy(direction);

		return this;

	}

	/**
	 * Copies the given ray.
	 *
	 * @param {Ray} r - A ray.
	 * @return {Ray} This ray.
	 */

	copy(r) {

		this.origin.copy(r.origin);
		this.direction.copy(r.direction);

		return this;

	}

	/**
	 * Clones this ray.
	 *
	 * @return {Ray} The cloned ray.
	 */

	clone() {

		return new this.constructor().copy(this);

	}

	/**
	 * Computes a point along the ray based on a given scalar t.
	 *
	 * @param {Number} t - The scalar.
	 * @param {Vector3} [target] - A target vector. If none is provided, a new one will be created.
	 * @return {Vector3} The point.
	 */

	at(t, target = new Vector3()) {

		return target.copy(this.direction).multiplyScalar(t).add(this.origin);

	}

	/**
	 * Rotates this ray to look at the given target.
	 *
	 * @param {Vector3} target - A point to look at.
	 * @return {Ray} This ray.
	 */

	lookAt(target) {

		this.direction.copy(target).sub(this.origin).normalize();

		return this;

	}

	/**
	 * Moves the origin along the ray by a given scalar t.
	 *
	 * @param {Number} t - The scalar.
	 * @return {Ray} This ray.
	 */

	recast(t) {

		this.origin.copy(this.at(t, v$4[0]));

		return this;

	}

	/**
	 * Finds the closest point along this ray to a given point.
	 *
	 * @param {Vector3} p - A point.
	 * @param {Vector3} [target] - A target vector. If none is provided, a new one will be created.
	 * @return {Vector3} The point.
	 */

	closestPointToPoint(p, target = new Vector3()) {

		const directionDistance = target.subVectors(p, this.origin).dot(this.direction);

		return (directionDistance >= 0.0) ?
			target.copy(this.direction).multiplyScalar(directionDistance).add(this.origin) :
			target.copy(this.origin);

	}

	/**
	 * Calculates the squared distance from this ray to the given point.
	 *
	 * @param {Vector3} p - The point.
	 * @return {Number} The squared distance.
	 */

	distanceSquaredToPoint(p) {

		const directionDistance = v$4[0].subVectors(p, this.origin).dot(this.direction);

		// Check if the point is behind the ray.
		return (directionDistance < 0.0) ?
			this.origin.distanceToSquared(p) :
			v$4[0].copy(this.direction).multiplyScalar(directionDistance).add(this.origin).distanceToSquared(p);

	}

	/**
	 * Calculates the distance from this ray to the given point.
	 *
	 * @param {Vector3} p - The point.
	 * @return {Number} The distance.
	 */

	distanceToPoint(p) {

		return Math.sqrt(this.distanceSquaredToPoint(p));

	}

	/**
	 * Calculates the distance from this ray to the given plane.
	 *
	 * @param {Plane} p - The plane.
	 * @return {Number} The distance, or null if the denominator is zero.
	 */

	distanceToPlane(p) {

		const denominator = p.normal.dot(this.direction);

		const t = (denominator !== 0.0) ?
			-(this.origin.dot(p.normal) + p.constant) / denominator :
			((p.distanceToPoint(this.origin) === 0.0) ? 0.0 : -1.0);

		return (t >= 0.0) ? t : null;

	}

	/**
	 * Calculates the distance from this ray to a given line segment.
	 *
	 * Based on:
	 *  http://www.geometrictools.com/GTEngine/Include/Mathematics/GteDistRaySegment.h
	 *
	 * @param {Vector3} v0 - The start of the segment.
	 * @param {Vector3} v1 - The end of the segment.
	 * @param {Vector3} [pointOnRay] - If provided, the point on this Ray that is closest to the segment will be stored in this vector.
	 * @param {Vector3} [pointOnSegment] - If provided, the point on the line segment that is closest to this ray will be stored in this vector.
	 * @return {Number} The smallest distance between the ray and the segment defined by v0 and v1.
	 */

	distanceSquaredToSegment(v0, v1, pointOnRay, pointOnSegment) {

		const segCenter = v$4[0].copy(v0).add(v1).multiplyScalar(0.5);
		const segDir = v$4[1].copy(v1).sub(v0).normalize();
		const diff = v$4[2].copy(this.origin).sub(segCenter);

		const segExtent = v0.distanceTo(v1) * 0.5;
		const a01 = -this.direction.dot(segDir);
		const b0 = diff.dot(this.direction);
		const b1 = -diff.dot(segDir);
		const c = diff.lengthSq();
		const det = Math.abs(1.0 - a01 * a01);

		let s0, s1, extDet, invDet, sqrDist;

		if(det > 0.0) {

			// The ray and segment are not parallel.
			s0 = a01 * b1 - b0;
			s1 = a01 * b0 - b1;
			extDet = segExtent * det;

			if(s0 >= 0.0) {

				if(s1 >= -extDet) {

					if(s1 <= extDet) {

						// Region 0.
						// Minimum at interior points of ray and segment.
						invDet = 1.0 / det;
						s0 *= invDet;
						s1 *= invDet;
						sqrDist = s0 * (s0 + a01 * s1 + 2.0 * b0) + s1 * (a01 * s0 + s1 + 2.0 * b1) + c;

					} else {

						// Region 1.
						s1 = segExtent;
						s0 = Math.max(0.0, -(a01 * s1 + b0));
						sqrDist = -s0 * s0 + s1 * (s1 + 2.0 * b1) + c;

					}

				} else {

					// Region 5.
					s1 = -segExtent;
					s0 = Math.max(0.0, -(a01 * s1 + b0));
					sqrDist = -s0 * s0 + s1 * (s1 + 2.0 * b1) + c;

				}

			} else {

				if(s1 <= -extDet) {

					// Region 4.
					s0 = Math.max(0.0, -(-a01 * segExtent + b0));
					s1 = (s0 > 0.0) ? -segExtent : Math.min(Math.max(-segExtent, -b1), segExtent);
					sqrDist = -s0 * s0 + s1 * (s1 + 2.0 * b1) + c;

				} else if(s1 <= extDet) {

					// Region 3.
					s0 = 0.0;
					s1 = Math.min(Math.max(-segExtent, -b1), segExtent);
					sqrDist = s1 * (s1 + 2.0 * b1) + c;

				} else {

					// Region 2.
					s0 = Math.max(0.0, -(a01 * segExtent + b0));
					s1 = (s0 > 0.0) ? segExtent : Math.min(Math.max(-segExtent, -b1), segExtent);
					sqrDist = -s0 * s0 + s1 * (s1 + 2.0 * b1) + c;

				}

			}

		} else {

			// Ray and segment are parallel.
			s1 = (a01 > 0.0) ? -segExtent : segExtent;
			s0 = Math.max(0.0, -(a01 * s1 + b0));
			sqrDist = -s0 * s0 + s1 * (s1 + 2.0 * b1) + c;

		}

		if(pointOnRay !== undefined) {

			pointOnRay.copy(this.direction).multiplyScalar(s0).add(this.origin);

		}

		if(pointOnSegment !== undefined) {

			pointOnSegment.copy(segDir).multiplyScalar(s1).add(segCenter);

		}

		return sqrDist;

	}

	/**
	 * Finds the point where this ray intersects the given sphere.
	 *
	 * @param {Sphere} s - A sphere.
	 * @param {Vector3} [target] - A target vector. If none is provided, a new one will be created.
	 * @return {Vector3} The point of intersection, or null if there is none.
	 */

	intersectSphere(s, target = new Vector3()) {

		const ab = v$4[0].subVectors(s.center, this.origin);
		const tca = ab.dot(this.direction);
		const d2 = ab.dot(ab) - tca * tca;
		const radius2 = s.radius * s.radius;

		let result = null;
		let thc, t0, t1;

		if(d2 <= radius2) {

			thc = Math.sqrt(radius2 - d2);

			// t0 = first intersection point - entrance on front of sphere.
			t0 = tca - thc;

			// t1 = second intersection point - exit point on back of sphere.
			t1 = tca + thc;

			// Check if both t0 and t1 are behind the ray - if so, return null.
			if(t0 >= 0.0 || t1 >= 0.0) {

				/* Check if t0 is behind the ray. If it is, the ray is inside the
				sphere, so return the second exit point scaled by t1 in order to always
				return an intersection point that is in front of the ray. If t0 is in
				front of the ray, return the first collision point scaled by t0. */
				result = (t0 < 0.0) ? this.at(t1, target) : this.at(t0, target);

			}

		}

		return result;

	}

	/**
	 * Determines whether this ray intersects the given sphere.
	 *
	 * @param {Sphere} s - A sphere.
	 * @return {Boolean} Whether this ray intersects the given sphere.
	 */

	intersectsSphere(s) {

		return (this.distanceSqToPoint(s.center) <= (s.radius * s.radius));

	}

	/**
	 * Finds the point where this ray intersects the given plane.
	 *
	 * @param {Plane} p - A plane.
	 * @param {Vector3} [target] - A target vector. If none is provided, a new one will be created.
	 * @return {Vector3} The point of intersection, or null if there is none.
	 */

	intersectPlane(p, target = new Vector3()) {

		const t = this.distanceToPlane(p);

		return (t === null) ? null : this.at(t, target);

	}

	/**
	 * Determines whether this ray intersects the given plane.
	 *
	 * @param {Plane} p - A plane.
	 * @return {Boolean} Whether this ray intersects the given plane.
	 */

	intersectsPlane(p) {

		const distanceToPoint = p.distanceToPoint(this.origin);

		return (distanceToPoint === 0.0 || p.normal.dot(this.direction) * distanceToPoint < 0.0);

	}

	/**
	 * Finds the point where this ray intersects the given box.
	 *
	 * @param {Plane} b - A box.
	 * @param {Vector3} [target] - A target vector. If none is provided, a new one will be created.
	 * @return {Vector3} The point of intersection, or null if there is none.
	 */

	intersectBox(b, target = new Vector3()) {

		const origin = this.origin;
		const direction = this.direction;
		const min = b.min;
		const max = b.max;

		const invDirX = 1.0 / direction.x;
		const invDirY = 1.0 / direction.y;
		const invDirZ = 1.0 / direction.z;

		let result = null;
		let tmin, tmax, tymin, tymax, tzmin, tzmax;

		if(invDirX >= 0.0) {

			tmin = (min.x - origin.x) * invDirX;
			tmax = (max.x - origin.x) * invDirX;

		} else {

			tmin = (max.x - origin.x) * invDirX;
			tmax = (min.x - origin.x) * invDirX;

		}

		if(invDirY >= 0.0) {

			tymin = (min.y - origin.y) * invDirY;
			tymax = (max.y - origin.y) * invDirY;

		} else {

			tymin = (max.y - origin.y) * invDirY;
			tymax = (min.y - origin.y) * invDirY;

		}

		if(tmin <= tymax && tymin <= tmax) {

			/* Handle the case where tmin or tmax is NaN (result of 0 * Infinity).
			Note: x !== x returns true if x is NaN. */
			if(tymin > tmin || tmin !== tmin) {

				tmin = tymin;

			}

			if(tymax < tmax || tmax !== tmax) {

				tmax = tymax;

			}

			if(invDirZ >= 0.0) {

				tzmin = (min.z - origin.z) * invDirZ;
				tzmax = (max.z - origin.z) * invDirZ;

			} else {

				tzmin = (max.z - origin.z) * invDirZ;
				tzmax = (min.z - origin.z) * invDirZ;

			}

			if(tmin <= tzmax && tzmin <= tmax) {

				if(tzmin > tmin || tmin !== tmin) {

					tmin = tzmin;

				}

				if(tzmax < tmax || tmax !== tmax) {

					tmax = tzmax;

				}

				// Return the closest point (positive side).
				if(tmax >= 0.0) {

					result = this.at((tmin >= 0.0) ? tmin : tmax, target);

				}

			}

		}

		return result;

	}

	/**
	 * Determines whether this ray intersects the given box.
	 *
	 * @param {Box3} b - A box.
	 * @return {Boolean} Whether this ray intersects the given box.
	 */

	intersectsBox(b) {

		return (this.intersectBox(b, v$4[0]) !== null);

	}

	/**
	 * Finds the point where this ray intersects the given triangle.
	 *
	 * Based on:
	 *  http://www.geometrictools.com/GTEngine/Include/Mathematics/GteIntrRay3Triangle3.h
	 *
	 * @param {Vector3} a - A triangle vertex.
	 * @param {Vector3} b - A triangle vertex.
	 * @param {Vector3} c - A triangle vertex.
	 * @param {Boolean} [backfaceCulling=false] - Whether backface culling should be considered.
	 * @param {Vector3} [target] - A target vector. If none is provided, a new one will be created.
	 * @return {Vector3} The point of intersection, or null if there is none.
	 */

	intersectTriangle(a, b, c, backfaceCulling, target) {

		const direction = this.direction;

		// Compute the offset origin, edges, and normal.
		const diff = v$4[0];
		const edge1 = v$4[1];
		const edge2 = v$4[2];
		const normal = v$4[3];

		let result = null;
		let DdN, sign, DdQxE2, DdE1xQ, QdN;

		edge1.subVectors(b, a);
		edge2.subVectors(c, a);
		normal.crossVectors(edge1, edge2);

		/* Solve Q + t * D = b1 * E1 + b2 * E2
		 * (Q = kDiff, D = ray direction, E1 = kEdge1, E2 = kEdge2,
		 * N = Cross(E1, E2)):
		 *
		 *   | Dot(D, N) | * b1 = sign(Dot(D, N)) * Dot(D, Cross(Q, E2))
		 *   | Dot(D, N) | * b2 = sign(Dot(D, N)) * Dot(D, Cross(E1, Q))
		 *   | Dot(D, N) | * t = -sign(Dot(D, N)) * Dot(Q, N)
		 */

		DdN = direction.dot(normal);

		// Discard coplanar constellations and cull backfaces.
		if(DdN !== 0.0 && !(backfaceCulling && DdN > 0.0)) {

			if(DdN > 0.0) {

				sign = 1.0;

			} else {

				sign = -1.0;
				DdN = -DdN;

			}

			diff.subVectors(this.origin, a);
			DdQxE2 = sign * direction.dot(edge2.crossVectors(diff, edge2));

			// b1 < 0, no intersection.
			if(DdQxE2 >= 0.0) {

				DdE1xQ = sign * direction.dot(edge1.cross(diff));

				// b2 < 0, or b1 + b2 > 1, no intersection.
				if(DdE1xQ >= 0.0 && DdQxE2 + DdE1xQ <= DdN) {

					// The line intersects the triangle, check if the ray does.
					QdN = -sign * diff.dot(normal);

					// t < 0, no intersection.
					if(QdN >= 0.0) {

						// Ray intersects triangle.
						result = this.at(QdN / DdN, target);

					}

				}

			}

		}

		return result;

	}

	/**
	 * Applies the given matrix to this ray.
	 *
	 * @param {Matrix4} m - A matrix.
	 * @return {Ray} This ray.
	 */

	applyMatrix4(m) {

		this.origin.applyMatrix4(m);
		this.direction.transformDirection(m);

		return this;

	}

	/**
	 * Checks if this ray equals the given one.
	 *
	 * @param {Ray} r - A ray.
	 * @return {Boolean} Whether the rays are equal.
	 */

	equals(r) {

		return (r.origin.equals(this.origin) && r.direction.equals(this.direction));

	}

}

/**
 * A spherical coordinate system.
 *
 * For details see: https://en.wikipedia.org/wiki/Spherical_coordinate_system
 *
 * The poles (phi) are at the positive and negative Y-axis. The equator starts
 * at positive Z.
 */

class Spherical {

	/**
	 * Constructs a new spherical system.
	 *
	 * @param {Number} [radius=1] - The radius of the sphere.
	 * @param {Number} [phi=0] - The polar angle phi.
	 * @param {Number} [theta=0] - The equator angle theta.
	 */

	constructor(radius = 1, phi = 0, theta = 0) {

		/**
		 * The radius of the sphere.
		 *
		 * @type {Number}
		 */

		this.radius = radius;

		/**
		 * The polar angle, up and down towards the top and bottom pole.
		 *
		 * @type {Number}
		 */

		this.phi = phi;

		/**
		 * The angle around the equator of the sphere.
		 *
		 * @type {Number}
		 */

		this.theta = theta;

	}

	/**
	 * Sets the values of this spherical system.
	 *
	 * @param {Number} radius - The radius.
	 * @param {Number} phi - Phi.
	 * @param {Number} theta - Theta.
	 * @return {Spherical} This spherical system.
	 */

	set(radius, phi, theta) {

		this.radius = radius;
		this.phi = phi;
		this.theta = theta;

		return this;

	}

	/**
	 * Copies the values of the given spherical system.
	 *
	 * @param {Spherical} s - A spherical system.
	 * @return {Spherical} This spherical system.
	 */

	copy(s) {

		this.radius = s.radius;
		this.phi = s.phi;
		this.theta = s.theta;

		return this;

	}

	/**
	 * Clones this spherical system.
	 *
	 * @return {Spherical} The cloned spherical system.
	 */

	clone() {

		return new this.constructor().copy(this);

	}

	/**
	 * Restricts phi to `[1e-6, PI - 1e-6]`.
	 *
	 * @return {Spherical} This spherical system.
	 */

	makeSafe() {

		this.phi = Math.max(1e-6, Math.min(Math.PI - 1e-6, this.phi));

		return this;

	}

	/**
	 * Sets the values of this spherical system based on a vector.
	 *
	 * The radius is set to the vector's length while phi and theta are set from
	 * its direction.
	 *
	 * @param {Vector3} v - The vector.
	 * @return {Spherical} This spherical system.
	 */

	setFromVector3(v) {

		return this.setFromCartesianCoords(v.x, v.y, v.z);

	}

	/**
	 * Sets the values of this spherical system based on cartesian coordinates.
	 *
	 * @param {Number} x - The X coordinate.
	 * @param {Number} y - The Y coordinate.
	 * @param {Number} z - The Z coordinate.
	 * @return {Spherical} This spherical system.
	 */

	setFromCartesianCoords(x, y, z) {

		this.radius = Math.sqrt(x * x + y * y + z * z);

		if(this.radius === 0) {

			this.theta = 0;
			this.phi = 0;

		} else {

			// Calculate the equator angle around the positive Y-axis.
			this.theta = Math.atan2(x, z);

			// Calculate the polar angle.
			this.phi = Math.acos(Math.min(Math.max(y / this.radius, -1), 1));

		}

		return this;

	}

}

/**
 * A symmetric 3x3 matrix.
 */

class SymmetricMatrix3 {

	/**
	 * Constructs a new symmetric matrix.
	 */

	constructor() {

		/**
		 * The matrix elements.
		 *
		 * @type {Float32Array}
		 */

		this.elements = new Float32Array([

			1, 0, 0,
			1, 0,
			1

		]);

	}

	/**
	 * Sets the values of this matrix.
	 *
	 * @param {Number} m00 - The value of the first row, first column.
	 * @param {Number} m01 - The value of the first row, second column and the second row, first column.
	 * @param {Number} m02 - The value of the first row, third column and the third row, first column.
	 * @param {Number} m11 - The value of the second row, second column.
	 * @param {Number} m12 - The value of the second row, third column and third row, second column.
	 * @param {Number} m22 - The value of the third row, third column.
	 * @return {SymmetricMatrix3} This matrix.
	 */

	set(m00, m01, m02, m11, m12, m22) {

		const e = this.elements;

		e[0] = m00;
		e[1] = m01; e[3] = m11;
		e[2] = m02; e[4] = m12; e[5] = m22;

		return this;

	}

	/**
	 * Sets this matrix to the identity matrix.
	 *
	 * @return {SymmetricMatrix3} This matrix.
	 */

	identity() {

		this.set(

			1, 0, 0,
			1, 0,
			1

		);

		return this;

	}

	/**
	 * Copies the values of a given symmetric matrix.
	 *
	 * @param {SymmetricMatrix3} m - A matrix.
	 * @return {SymmetricMatrix3} This matrix.
	 */

	copy(m) {

		const me = m.elements;

		this.set(

			me[0], me[1], me[2],
			me[3], me[4],
			me[5]

		);

		return this;

	}

	/**
	 * Clones this matrix.
	 *
	 * @return {SymmetricMatrix3} A clone of this matrix.
	 */

	clone() {

		return new this.constructor().copy(this);

	}

	/**
	 * Copies this symmetric matrix into a given 3x3 matrix.
	 *
	 * @param {Matrix3} m - The target matrix.
	 */

	toMatrix3(m) {

		const me = m.elements;

		m.set(

			me[0], me[1], me[2],
			me[1], me[3], me[4],
			me[2], me[4], me[5]

		);

	}

	/**
	 * Adds the values of a given symmetric matrix to this one.
	 *
	 * @param {SymmetricMatrix3} m - A matrix.
	 * @return {SymmetricMatrix3} This matrix.
	 */

	add(m) {

		const te = this.elements;
		const me = m.elements;

		te[0] += me[0];
		te[1] += me[1]; te[3] += me[3];
		te[2] += me[2]; te[4] += me[4]; te[5] += me[5];

		return this;

	}

	/**
	 * Calculates the Frobenius norm of this matrix.
	 *
	 * @return {Number} The norm of this matrix.
	 */

	norm() {

		const e = this.elements;

		const m01m01 = e[1] * e[1];
		const m02m02 = e[2] * e[2];
		const m12m12 = e[4] * e[4];

		return Math.sqrt(

			e[0] * e[0] + m01m01 + m02m02 +
			m01m01 + e[3] * e[3] + m12m12 +
			m02m02 + m12m12 + e[5] * e[5]

		);

	}

	/**
	 * Calculates the absolute sum of all matrix components except for the main
	 * diagonal.
	 *
	 * @return {Number} The offset of this matrix.
	 */

	off() {

		const e = this.elements;

		return Math.sqrt(2 * (

			// Diagonal = [0, 3, 5].
			e[1] * e[1] + e[2] * e[2] + e[4] * e[4]

		));

	}

	/**
	 * Applies this symmetric matrix to a vector.
	 *
	 * @param {Vector3} v - The vector to modify.
	 * @return {Vector3} The modified vector.
	 */

	applyToVector3(v) {

		const x = v.x, y = v.y, z = v.z;
		const e = this.elements;

		v.x = e[0] * x + e[1] * y + e[2] * z;
		v.y = e[1] * x + e[3] * y + e[4] * z;
		v.z = e[2] * x + e[4] * y + e[5] * z;

		return v;

	}

	/**
	 * Checks if this matrix equals the given one.
	 *
	 * @param {SymmetricMatrix3} m - A matrix.
	 * @return {Boolean} Whether the matrices are equal.
	 */

	equals(m) {

		const te = this.elements;
		const me = m.elements;

		let result = true;
		let i;

		for(i = 0; result && i < 6; ++i) {

			if(te[i] !== me[i]) {

				result = false;

			}

		}

		return result;

	}

	/**
	 * Calculates the linear index of an element from this matrix.
	 *
	 * Let N be the dimension of the symmetric matrix:
	 *
	 *     index = N * (N - 1) / 2 - (N - i) * (N - i - 1) / 2 + j
	 *
	 * @param {Number} i - The row.
	 * @param {Number} j - The column.
	 * @return {Number} The index into the elements of this matrix.
	 */

	static calculateIndex(i, j) {

		return (3 - (3 - i) * (2 - i) / 2 + j);

	}

}

/**
 * A vector with four components.
 */

class Vector4 {

	/**
	 * Constructs a new vector.
	 *
	 * @param {Number} [x=0] - The X component.
	 * @param {Number} [y=0] - The Y component.
	 * @param {Number} [z=0] - The Z component.
	 * @param {Number} [w=0] - The W component.
	 */

	constructor(x = 0, y = 0, z = 0, w = 0) {

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

		/**
		 * The W component.
		 *
		 * @type {Number}
		 */

		this.w = w;

	}

	/**
	 * Sets the values of this vector
	 *
	 * @param {Number} x - The X component.
	 * @param {Number} y - The Y component.
	 * @param {Number} z - The Z component.
	 * @param {Number} w - The W component.
	 * @return {Vector4} This vector.
	 */

	set(x, y, z, w) {

		this.x = x;
		this.y = y;
		this.z = z;
		this.w = w;

		return this;

	}

	/**
	 * Copies the values of another vector.
	 *
	 * @param {Vector4} v - A vector.
	 * @return {Vector4} This vector.
	 */

	copy(v) {

		this.x = v.x;
		this.y = v.y;
		this.z = v.z;
		this.w = v.w;

		return this;

	}

	/**
	 * Clones this vector.
	 *
	 * @return {Vector4} A clone of this vector.
	 */

	clone() {

		return new this.constructor(this.x, this.y, this.z, this.w);

	}

	/**
	 * Copies values from an array.
	 *
	 * @param {Number[]} array - An array.
	 * @param {Number} offset - An offset.
	 * @return {Vector4} This vector.
	 */

	fromArray(array, offset = 0) {

		this.x = array[offset];
		this.y = array[offset + 1];
		this.z = array[offset + 2];
		this.w = array[offset + 3];

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
		array[offset + 3] = this.w;

		return array;

	}

	/**
	 * Stores the axis angle from the given quaternion in this vector.
	 *
	 * For more details see:
	 *  http://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToAngle/index.htm
	 *
	 * @param {Quaternion} q - A quaternion. Assumed to be normalized
	 * @return {Vector4} This vector.
	 */

	setAxisAngleFromQuaternion(q) {

		this.w = 2 * Math.acos(q.w);

		const s = Math.sqrt(1 - q.w * q.w);

		if(s < 1e-4) {

			this.x = 1;
			this.y = 0;
			this.z = 0;

		} else {

			this.x = q.x / s;
			this.y = q.y / s;
			this.z = q.z / s;

		}

		return this;

	}

	/**
	 * Stores the axis angle from the given rotation matrix in this vector.
	 *
	 * For more details see:
	 *  http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToAngle/index.htm
	 *
	 * @param {Matrix4} m - A matrix. The upper 3x3 must be a pure rotation matrix (i.e. unscaled).
	 * @return {Vector4} This vector.
	 */

	setAxisAngleFromRotationMatrix(m) {

		// Margin to allow for rounding errors.
		const E = 0.01;
		// Margin to distinguish between 0 and 180 degrees.
		const H = 0.1;

		const me = m.elements;
		const m00 = me[0], m01 = me[4], m02 = me[8];
		const m10 = me[1], m11 = me[5], m12 = me[9];
		const m20 = me[2], m21 = me[6], m22 = me[10];

		let angle;
		let x, y, z;
		let xx, yy, zz;
		let xy, xz, yz;
		let s;

		if((Math.abs(m01 - m10) < E) && (Math.abs(m02 - m20) < E) && (Math.abs(m12 - m21) < E)) {

			/* Singularity found. First, check for identity matrix which must have +1
			for all terms in the leading diagonal and zero in other terms. */
			if((Math.abs(m01 + m10) < H) && (Math.abs(m02 + m20) < H) && (Math.abs(m12 + m21) < H) && (Math.abs(m00 + m11 + m22 - 3) < H)) {

				// This singularity is the identity matrix. The angle is zero.
				this.set(1, 0, 0, 0);

			} else {

				// The angle is 180.
				angle = Math.PI;

				xx = (m00 + 1) / 2;
				yy = (m11 + 1) / 2;
				zz = (m22 + 1) / 2;
				xy = (m01 + m10) / 4;
				xz = (m02 + m20) / 4;
				yz = (m12 + m21) / 4;

				if((xx > yy) && (xx > zz)) {

					// m00 is the largest diagonal term.
					if(xx < E) {

						x = 0;
						y = 0.707106781;
						z = 0.707106781;

					} else {

						x = Math.sqrt(xx);
						y = xy / x;
						z = xz / x;

					}

				} else if(yy > zz) {

					// m11 is the largest diagonal term.
					if(yy < E) {

						x = 0.707106781;
						y = 0;
						z = 0.707106781;

					} else {

						y = Math.sqrt(yy);
						x = xy / y;
						z = yz / y;

					}

				} else {

					// m22 is the largest diagonal term.
					if(zz < E) {

						x = 0.707106781;
						y = 0.707106781;
						z = 0;

					} else {

						z = Math.sqrt(zz);
						x = xz / z;
						y = yz / z;

					}

				}

				this.set(x, y, z, angle);

			}

		} else {

			// There are no singularities.
			s = Math.sqrt(
				(m21 - m12) * (m21 - m12) +
				(m02 - m20) * (m02 - m20) +
				(m10 - m01) * (m10 - m01)
			);

			// Prevent division by zero.
			if(Math.abs(s) < 0.001) {

				s = 1;

			}

			this.x = (m21 - m12) / s;
			this.y = (m02 - m20) / s;
			this.z = (m10 - m01) / s;
			this.w = Math.acos((m00 + m11 + m22 - 1) / 2);

		}

		return this;

	}

	/**
	 * Adds a vector to this one.
	 *
	 * @param {Vector4} v - The vector to add.
	 * @return {Vector4} This vector.
	 */

	add(v) {

		this.x += v.x;
		this.y += v.y;
		this.z += v.z;
		this.w += v.w;

		return this;

	}

	/**
	 * Adds a scalar to this vector.
	 *
	 * @param {Number} s - The scalar to add.
	 * @return {Vector4} This vector.
	 */

	addScalar(s) {

		this.x += s;
		this.y += s;
		this.z += s;
		this.w += s;

		return this;

	}

	/**
	 * Sets this vector to the sum of two given vectors.
	 *
	 * @param {Vector4} a - A vector.
	 * @param {Vector4} b - Another vector.
	 * @return {Vector4} This vector.
	 */

	addVectors(a, b) {

		this.x = a.x + b.x;
		this.y = a.y + b.y;
		this.z = a.z + b.z;
		this.w = a.w + b.w;

		return this;

	}

	/**
	 * Adds a scaled vector to this one.
	 *
	 * @param {Vector4} v - The vector to scale and add.
	 * @param {Number} s - A scalar.
	 * @return {Vector4} This vector.
	 */

	addScaledVector(v, s) {

		this.x += v.x * s;
		this.y += v.y * s;
		this.z += v.z * s;
		this.w += v.w * s;

		return this;

	}

	/**
	 * Subtracts a vector from this vector.
	 *
	 * @param {Vector4} v - The vector to subtract.
	 * @return {Vector4} This vector.
	 */

	sub(v) {

		this.x -= v.x;
		this.y -= v.y;
		this.z -= v.z;
		this.w -= v.w;

		return this;

	}

	/**
	 * Subtracts a scalar from this vector.
	 *
	 * @param {Number} s - The scalar to subtract.
	 * @return {Vector4} This vector.
	 */

	subScalar(s) {

		this.x -= s;
		this.y -= s;
		this.z -= s;
		this.w -= s;

		return this;

	}

	/**
	 * Sets this vector to the difference between two given vectors.
	 *
	 * @param {Vector4} a - A vector.
	 * @param {Vector4} b - A second vector.
	 * @return {Vector4} This vector.
	 */

	subVectors(a, b) {

		this.x = a.x - b.x;
		this.y = a.y - b.y;
		this.z = a.z - b.z;
		this.w = a.w - b.w;

		return this;

	}

	/**
	 * Multiplies this vector with another vector.
	 *
	 * @param {Vector4} v - A vector.
	 * @return {Vector4} This vector.
	 */

	multiply(v) {

		this.x *= v.x;
		this.y *= v.y;
		this.z *= v.z;
		this.w *= v.w;

		return this;

	}

	/**
	 * Multiplies this vector with a given scalar.
	 *
	 * @param {Number} s - A scalar.
	 * @return {Vector4} This vector.
	 */

	multiplyScalar(s) {

		this.x *= s;
		this.y *= s;
		this.z *= s;
		this.w *= s;

		return this;

	}

	/**
	 * Sets this vector to the product of two given vectors.
	 *
	 * @param {Vector4} a - A vector.
	 * @param {Vector4} b - Another vector.
	 * @return {Vector4} This vector.
	 */

	multiplyVectors(a, b) {

		this.x = a.x * b.x;
		this.y = a.y * b.y;
		this.z = a.z * b.z;
		this.w = a.w * b.w;

		return this;

	}

	/**
	 * Divides this vector by another vector.
	 *
	 * @param {Vector4} v - A vector.
	 * @return {Vector4} This vector.
	 */

	divide(v) {

		this.x /= v.x;
		this.y /= v.y;
		this.z /= v.z;
		this.w /= v.w;

		return this;

	}

	/**
	 * Divides this vector by a given scalar.
	 *
	 * @param {Number} s - A scalar.
	 * @return {Vector4} This vector.
	 */

	divideScalar(s) {

		this.x /= s;
		this.y /= s;
		this.z /= s;
		this.w /= s;

		return this;

	}

	/**
	 * Applies a matrix to this vector.
	 *
	 * @param {Matrix4} m - A matrix.
	 * @return {Vector4} This vector.
	 */

	applyMatrix4(m) {

		const x = this.x, y = this.y, z = this.z, w = this.w;
		const e = m.elements;

		this.x = e[0] * x + e[4] * y + e[8] * z + e[12] * w;
		this.y = e[1] * x + e[5] * y + e[9] * z + e[13] * w;
		this.z = e[2] * x + e[6] * y + e[10] * z + e[14] * w;
		this.w = e[3] * x + e[7] * y + e[11] * z + e[15] * w;

		return this;

	}

	/**
	 * Negates this vector.
	 *
	 * @return {Vector4} This vector.
	 */

	negate() {

		this.x = -this.x;
		this.y = -this.y;
		this.z = -this.z;
		this.w = -this.w;

		return this;

	}

	/**
	 * Calculates the dot product with another vector.
	 *
	 * @param {Vector4} v - A vector.
	 * @return {Number} The dot product.
	 */

	dot(v) {

		return this.x * v.x + this.y * v.y + this.z * v.z + this.w * v.w;

	}

	/**
	 * Calculates the Manhattan length of this vector.
	 *
	 * @return {Number} The length.
	 */

	manhattanLength() {

		return Math.abs(this.x) + Math.abs(this.y) + Math.abs(this.z) + Math.abs(this.w);

	}

	/**
	 * Calculates the squared length of this vector.
	 *
	 * @return {Number} The squared length.
	 */

	lengthSquared() {

		return (
			this.x * this.x +
			this.y * this.y +
			this.z * this.z +
			this.w * this.w
		);

	}

	/**
	 * Calculates the length of this vector.
	 *
	 * @return {Number} The length.
	 */

	length() {

		return Math.sqrt(
			this.x * this.x +
			this.y * this.y +
			this.z * this.z +
			this.w * this.w
		);

	}

	/**
	 * Calculates the Manhattan distance to a given vector.
	 *
	 * @param {Vector4} v - A vector.
	 * @return {Number} The distance.
	 */

	manhattanDistanceTo(v) {

		return (
			Math.abs(this.x - v.x) +
			Math.abs(this.y - v.y) +
			Math.abs(this.z - v.z) +
			Math.abs(this.w - v.w)
		);

	}

	/**
	 * Calculates the squared distance to a given vector.
	 *
	 * @param {Vector4} v - A vector.
	 * @return {Number} The squared distance.
	 */

	distanceToSquared(v) {

		const dx = this.x - v.x;
		const dy = this.y - v.y;
		const dz = this.z - v.z;
		const dw = this.w - v.w;

		return dx * dx + dy * dy + dz * dz + dw * dw;

	}

	/**
	 * Calculates the distance to a given vector.
	 *
	 * @param {Vector4} v - A vector.
	 * @return {Number} The distance.
	 */

	distanceTo(v) {

		return Math.sqrt(this.distanceToSquared(v));

	}

	/**
	 * Normalizes this vector.
	 *
	 * @return {Vector4} This vector.
	 */

	normalize() {

		return this.divideScalar(this.length());

	}

	/**
	 * Sets the length of this vector.
	 *
	 * @param {Number} length - The new length.
	 * @return {Vector4} This vector.
	 */

	setLength(length) {

		return this.normalize().multiplyScalar(length);

	}

	/**
	 * Adopts the min value for each component of this vector and the given one.
	 *
	 * @param {Vector4} v - A vector.
	 * @return {Vector4} This vector.
	 */

	min(v) {

		this.x = Math.min(this.x, v.x);
		this.y = Math.min(this.y, v.y);
		this.z = Math.min(this.z, v.z);
		this.w = Math.min(this.w, v.w);

		return this;

	}

	/**
	 * Adopts the max value for each component of this vector and the given one.
	 *
	 * @param {Vector4} v - A vector.
	 * @return {Vector4} This vector.
	 */

	max(v) {

		this.x = Math.max(this.x, v.x);
		this.y = Math.max(this.y, v.y);
		this.z = Math.max(this.z, v.z);
		this.w = Math.max(this.w, v.w);

		return this;

	}

	/**
	 * Clamps this vector.
	 *
	 * @param {Vector4} min - The lower bounds. Assumed to be smaller than max.
	 * @param {Vector4} max - The upper bounds. Assumed to be greater than min.
	 * @return {Vector4} This vector.
	 */

	clamp(min, max) {

		this.x = Math.max(min.x, Math.min(max.x, this.x));
		this.y = Math.max(min.y, Math.min(max.y, this.y));
		this.z = Math.max(min.z, Math.min(max.z, this.z));
		this.w = Math.max(min.w, Math.min(max.w, this.w));

		return this;

	}

	/**
	 * Floors this vector.
	 *
	 * @return {Vector4} This vector.
	 */

	floor() {

		this.x = Math.floor(this.x);
		this.y = Math.floor(this.y);
		this.z = Math.floor(this.z);
		this.w = Math.floor(this.w);

		return this;

	}

	/**
	 * Ceils this vector.
	 *
	 * @return {Vector4} This vector.
	 */

	ceil() {

		this.x = Math.ceil(this.x);
		this.y = Math.ceil(this.y);
		this.z = Math.ceil(this.z);
		this.w = Math.ceil(this.w);

		return this;

	}

	/**
	 * Rounds this vector.
	 *
	 * @return {Vector4} This vector.
	 */

	round() {

		this.x = Math.round(this.x);
		this.y = Math.round(this.y);
		this.z = Math.round(this.z);
		this.w = Math.round(this.w);

		return this;

	}

	/**
	 * Lerps towards the given vector.
	 *
	 * @param {Vector4} v - The target vector.
	 * @param {Number} alpha - The lerp factor.
	 * @return {Vector4} This vector.
	 */

	lerp(v, alpha) {

		this.x += (v.x - this.x) * alpha;
		this.y += (v.y - this.y) * alpha;
		this.z += (v.z - this.z) * alpha;
		this.w += (v.w - this.w) * alpha;

		return this;

	}

	/**
	 * Sets this vector to the lerp result of the given vectors.
	 *
	 * @param {Vector4} v1 - A base vector.
	 * @param {Vector4} v2 - The target vector.
	 * @param {Number} alpha - The lerp factor.
	 * @return {Vector4} This vector.
	 */

	lerpVectors(v1, v2, alpha) {

		return this.subVectors(v2, v1).multiplyScalar(alpha).add(v1);

	}

	/**
	 * Checks if this vector equals the given one.
	 *
	 * @param {Vector4} v - A vector.
	 * @return {Boolean} Whether this vector equals the given one.
	 */

	equals(v) {

		return (v.x === this.x && v.y === this.y && v.z === this.z && v.w === this.w);

	}

}

export { Box2, Box3, Cylindrical, Euler, Frustum, Line3, Matrix3, Matrix4, Plane, Quaternion, Ray, RotationOrder, Sphere, Spherical, SymmetricMatrix3, Vector2, Vector3, Vector4 };
