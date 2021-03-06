import { RotationOrder } from "./RotationOrder.js";
import { Vector3 } from "./Vector3.js";

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
 * A vector.
 *
 * @type {Vector3}
 * @private
 */

const c = new Vector3();

/**
 * A 4x4 matrix.
 */

export class Matrix4 {

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

		const scaleX = 1.0 / a.setFromMatrixColumn(m, 0).length();
		const scaleY = 1.0 / a.setFromMatrixColumn(m, 1).length();
		const scaleZ = 1.0 / a.setFromMatrixColumn(m, 2).length();

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

		return this.compose(a.set(0, 0, 0), q, b.set(1, 1, 1));

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
		const x = a, y = b, z = c;

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

			this.set(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

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
		const sx = a.set(n00, n10, n20).length() * ((det < 0) ? -1 : 1);
		const sy = a.set(n01, n11, n21).length();
		const sz = a.set(n02, n12, n22).length();

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
