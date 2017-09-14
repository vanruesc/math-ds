import { RotationOrder } from "./RotationOrder.js";
import { Vector3 } from "./Vector3.js";

/**
 * A vector.
 *
 * @type {Vector3}
 * @private
 */

const v = new Vector3();

/**
 * A quaternion.
 */

export class Quaternion {

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

				v.set(-vFrom.y, vFrom.x, 0);

			} else {

				v.set(0, -vFrom.z, vFrom.y);

			}

		} else {

			v.crossVectors(vFrom, vTo);

		}

		this.x = v.x;
		this.y = v.y;
		this.z = v.z;
		this.w = r;

		return this.normalize();

	}

	/**
	 * Inverts this quaternion.
	 *
	 * @return {Quaternion} This quaternion.
	 */

	invert() {

		return this.conjugate().normalize();

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
	 * @param {Quaternion} q - A quaternion.
	 * @param {Quaternion} q - Another quaternion.
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

		let cosHalfTheta, sinHalfTheta;
		let halfTheta, ratioA, ratioB;

		if(t === 1) {

			this.copy(q);

		} else if(t > 0) {

			cosHalfTheta = w * q.w + x * q.x + y * q.y + z * q.z;

			if(cosHalfTheta < 0) {

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

				return this;

			}

			sinHalfTheta = Math.sqrt(1.0 - cosHalfTheta * cosHalfTheta);

			if(Math.abs(sinHalfTheta) < 1e-3) {

				this.w = 0.5 * (w + this.w);
				this.x = 0.5 * (x + this.x);
				this.y = 0.5 * (y + this.y);
				this.z = 0.5 * (z + this.z);

				return this;

			}

			halfTheta = Math.atan2(sinHalfTheta, cosHalfTheta);
			ratioA = Math.sin((1.0 - t) * halfTheta) / sinHalfTheta;
			ratioB = Math.sin(t * halfTheta) / sinHalfTheta;

			this.w = (w * ratioA + this.w * ratioB);
			this.x = (x * ratioA + this.x * ratioB);
			this.y = (y * ratioA + this.y * ratioB);
			this.z = (z * ratioA + this.z * ratioB);

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
	 * @param {Quaternion} qr - A quaternion to store the result in. If none is provided, a new one will be created.
	 * @param {Number} t - The slerp factor.
	 * @return {Quaternion} The resulting quaternion.
	 */

	static slerp(qa, qb, qr = new Quaternion(), t) {

		return qr.copy(qa).slerp(qb, t);

	}

	/**
	 * Performs an array-based spherical linear interpolation.
	 *
	 * @param {Number[]} dst - An array to store the result in.
	 * @param {Number} dstOffset - An offset into the destination array.
	 * @param {Number[]} src0 - An array that contains the base quaternion values.
	 * @param {Number} src0Offset - An offset into the base array.
	 * @param {Number[]} src1 - An array that contains the target quaternion values.
	 * @param {Number} src1Offset - An offset into the target array.
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
