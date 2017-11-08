import { Matrix3 } from "./Matrix3.js";
import { Quaternion } from "./Quaternion.js";
import { RotationOrder } from "./RotationOrder.js";
import { Vector3 } from "./Vector3.js";

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

export class Euler {

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
		this.order = z;

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

		const THRESHOLD = 1.0 - 1e-5;

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
