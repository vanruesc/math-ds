/**
 * A vector with two components.
 */

export class Vector2 {

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
	 * Randomizes the values of this vector
	 *
	 * @return {Vector2} This vector.
	 */

	random() {

		this.x = Math.random();
		this.y = Math.random();

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

		return Math.atan2(-this.y, -this.x) + Math.PI;

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

		this.x = v1.x + (v2.x - v1.x) * alpha;
		this.y = v1.y + (v2.y - v1.y) * alpha;

		return this;

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
