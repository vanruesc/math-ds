/**
 * A 3x3 matrix.
 */

export class Matrix3 {

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

			this.set(0, 0, 0, 0, 0, 0, 0, 0, 0);

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
