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
 * A plane.
 */

export class Plane {

	/**
	 * Constructs a new plane.
	 *
	 * @param {Vector3} [normal] - The normal.
	 * @param {Number} [constant] - The constant.
	 */

	constructor(normal = new Vector3(1, 0, 0), constant = 0) {

		/**
		 * The normal.
		 *
		 * @type {Vector3}
		 */

		this.normal = normal;

		/**
		 * The constant.
		 *
		 * @type {Number}
		 */

		this.constant = constant;

	}

	/**
	 * Sets the normal and the constant.
	 *
	 * @param {Vector3} normal - The normal.
	 * @param {Number} constant - The constant.
	 * @return {Plane} This plane.
	 */

	set(normal, constant) {

		this.normal.copy(normal);
		this.constant = constant;

		return this;

	}

	/**
	 * Sets the components of this plane.
	 *
	 * @param {Number} x - The X component of the normal.
	 * @param {Number} y - The Y component of the normal.
	 * @param {Number} z - The Z component of the normal.
	 * @param {Number} w - The constant.
	 * @return {Plane} This plane.
	 */

	setComponents(x, y, z, w) {

		this.normal.set(x, y, z);
		this.constant = w;

		return this;

	}

	/**
	 * Copies the given plane.
	 *
	 * @param {Plane} p - A plane.
	 * @return {Plane} This plane.
	 */

	copy(p) {

		this.normal.copy(p.normal);
		this.constant = p.constant;

		return this;

	}

	/**
	 * Clones this plane.
	 *
	 * @return {Plane} The cloned plane.
	 */

	clone() {

		return new this.constructor().copy(this);

	}

	/**
	 * Sets this plane from a normal and a coplanar point.
	 *
	 * @param {Vector3} n - The normal.
	 * @param {Vector3} p - The coplanar point.
	 * @return {Sphere} This sphere.
	 */

	setFromNormalAndCoplanarPoint(n, p) {

		this.normal.copy(n);
		this.constant = -p.dot(this.normal);

		return this;

	}

	/**
	 * Sets this plane from three distinct coplanar points.
	 *
	 * @param {Vector3} p0 - A coplanar point.
	 * @param {Vector3} p1 - A coplanar point.
	 * @param {Vector3} p2 - A coplanar point.
	 * @return {Plane} This plane.
	 */

	setFromCoplanarPoints(p0, p1, p2) {

		const normal = a.subVectors(p2, p1).cross(b.subVectors(p0, p1)).normalize();

		this.setFromNormalAndCoplanarPoint(normal, a);

		return this;

	}

	/**
	 * Normalizes this plane.
	 *
	 * @return {Plane} This plane.
	 */

	normalize() {

		const inverseNormalLength = 1.0 / this.normal.length();

		this.normal.multiplyScalar(inverseNormalLength);
		this.constant *= inverseNormalLength;

		return this;

	}

	/**
	 * Negates this plane.
	 *
	 * @return {Plane} This plane.
	 */

	negate() {

		this.normal.negate();
		this.constant = -this.constant;

		return this;

	}

	/**
	 * Calculates the distance from this plane to the given point.
	 *
	 * @param {Vector3} p - A point.
	 * @return {Number} The length.
	 */

	distanceToPoint(p) {

		return this.normal.dot(p) + this.constant;

	}

	/**
	 * Calculates the distance from this plane to the given sphere.
	 *
	 * @param {Sphere} s - A sphere.
	 * @return {Number} The length.
	 */

	distanceToSphere(s) {

		return this.distanceToPoint(s.center) - s.radius;

	}

	/**
	 * Projects the given point on this plane.
	 *
	 * @param {Vector3} p - A point.
	 * @param {Vector3} [target] - A target vector. If none is provided, a new one will be created.
	 * @return {Vector3} The projected point.
	 */

	projectPoint(p, target) {

		return target.copy(this.normal).multiplyScalar(-this.distanceToPoint(p)).add(p);

	}

	/**
	 * Calculates a coplanar point and returns it.
	 *
	 * @param {Vector3} [target] - A target vector. If none is provided, a new one will be created.
	 * @return {Vector3} A coplanar plane.
	 */

	coplanarPoint(target) {

		return target.copy(this.normal).multiplyScalar(-this.constant);

	}

	/**
	 * Translates this plane.
	 *
	 * @param {Vector3} offset - An offset.
	 * @return {Plane} This plane.
	 */

	translate(offset) {

		this.constant -= offset.dot(this.normal);

		return this;

	}

	/**
	 * Finds the point of intersection between this plane and a given line.
	 *
	 * @param {Line3} l - A line.
	 * @param {Vector3} [target] - A target vector. If none is provided, a new one will be created.
	 * @return {Vector3} The intersection point.
	 */

	intersectLine(l, target) {

		const direction = l.delta(a);
		const denominator = this.normal.dot(direction);

		if(denominator === 0) {

			// The line is coplanar, return origin.
			if(this.distanceToPoint(l.start) === 0) {

				target.copy(l.start);

			}

		} else {

			const t = -(l.start.dot(this.normal) + this.constant) / denominator;

			if(t >= 0 && t <= 1) {

				target.copy(direction).multiplyScalar(t).add(l.start);

			}

		}

		return target;

	}

	/**
	 * Checks if this plane intersects with the given line.
	 *
	 * @param {Line3} l - A line.
	 * @return {Boolean} Whether this plane intersects with the given line.
	 */

	intersectsLine(l) {

		const startSign = this.distanceToPoint(l.start);
		const endSign = this.distanceToPoint(l.end);

		return ((startSign < 0 && endSign > 0) || (endSign < 0 && startSign > 0));

	}

	/**
	 * Checks if this plane intersects with the given box.
	 *
	 * @param {Box3} b - A box.
	 * @return {Boolean} Whether this plane intersects with the given box.
	 */

	intersectsBox(b) {

		return b.intersectsPlane(this);

	}

	/**
	 * Checks if this plane intersects with the given sphere.
	 *
	 * @param {Sphere} s - A sphere.
	 * @return {Boolean} Whether this plane intersects with the given sphere.
	 */

	intersectsSphere(s) {

		return s.intersectsPlane(this);

	}

	/**
	 * Checks if this plane equals the given one.
	 *
	 * @param {Plane} p - A plane.
	 * @return {Boolean} Whether this plane equals the given one.
	 */

	equals(p) {

		return (p.normal.equals(this.normal) && (p.constant === this.constant));

	}

}
