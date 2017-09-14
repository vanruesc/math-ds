import { Box3 } from "./Box3.js";
import { Vector3 } from "./Vector3.js";

/**
 * A box.
 *
 * @type {Box3}
 * @private
 */

const box = new Box3();

/**
 * A sphere.
 */

export class Sphere {

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
	 * @param {Sphere} sphere - A sphere.
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
	 * @param {Sphere} [target] - A target vector. If none is provided, a new one will be created.
	 * @return {Sphere} This sphere.
	 */

	setFromPoints(points, target = box.setFromPoints(points).getCenter(this.center)) {

		const center = this.center;

		let maxRadiusSq = 0;
		let i, l;

		for(i = 0, l = points.length; i < l; ++i) {

			maxRadiusSq = Math.max(maxRadiusSq, center.distanceToSquared(points[i]));

		}

		this.radius = Math.sqrt(maxRadiusSq);

		return this;

	}

	/**
	 * Calculates the bounding box of this sphere.
	 *
	 * @param {Box3} [target] - A target sphere. If none is provided, a new one will be created.
	 * @return {Box3} The bounding box.
	 */

	getBoundingBox(target = new Box3()) {

		target.set(this.center, this.center);
		target.expandByScalar(this.radius);

		return target;

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
