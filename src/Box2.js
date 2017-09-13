import { Sphere } from "./Sphere.js";
import { Vector2 } from "./Vector2.js";

/**
 * A vector.
 *
 * @type {Vector2}
 * @private
 */

const v = new Vector2();

/**
 * A 2D box.
 */

export class Box2 {

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

		target.radius = this.getSize(v).length() * 0.5;

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

		const halfSize = v.copy(size).multiplyScalar(0.5);

		this.min.copy(center).sub(halfSize);
		this.max.copy(center).add(halfSize);

		return this;

	}

	/**
	 * Clamps the given point to the boundaries of this box.
	 *
	 * @param {Vector2} p - A point.
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

		const clampedPoint = v.copy(p).clamp(this.min, this.max);

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
	 * Expands this box by combining it with the given one.
	 *
	 * @param {Box2} b - A box.
	 * @return {Box2} This box.
	 */

	intersect(b) {

		this.min.max(b.min);
		this.max.min(b.max);

		/* Ensure that if there is no overlap, the result is fully empty to prevent
		subsequent intersections to erroneously return valid values. */
		if(this.isEmpty()) { this.makeEmpty(); }

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

		return !(
			p.x < this.min.x || p.x > this.max.x ||
			p.y < this.min.y || p.y > this.max.y
		);

	}

	/**
	 * Checks if the given box lies inside this box.
	 *
	 * @param {Vector2} b - A box.
	 * @return {Boolean} Whether this box contains the given one.
	 */

	containsBox(b) {

		return (
			this.min.x <= b.min.x && b.max.x <= this.max.x &&
			this.min.y <= b.min.y && b.max.y <= this.max.y
		);

	}

	/**
	 * Checks if this box intersects with the given one.
	 *
	 * @param {Box2} b - A box.
	 * @return {Boolean} Whether the boxes intersect.
	 */

	intersectsBox(b) {

		return !(
			b.max.x < this.min.x || b.min.x > this.max.x ||
			b.max.y < this.min.y || b.min.y > this.max.y
		);

	}

	/**
	 * Checks if this box equals the given one.
	 *
	 * @param {Box2} v - A box.
	 * @return {Boolean} Whether this box equals the given one.
	 */

	equals(b) {

		return (b.min.equals(this.min) && b.max.equals(this.max));

	}

}
