/**
 * A vector with four components.
 */

export class Vector4 {

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
			if(Math.abs(s) < 0.001) { s = 1; }

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
