import { Vector3 } from "./Vector3.js";

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

export class Box3 {

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
