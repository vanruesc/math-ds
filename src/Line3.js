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
 * A line.
 */

export class Line3 {

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

		a.subVectors(p, this.start);
		b.subVectors(this.end, this.start);

		const bb = b.dot(b);
		const ba = b.dot(a);

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
