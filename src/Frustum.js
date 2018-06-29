import { Plane } from "./Plane.js";
import { Vector3 } from "./Vector3.js";

/**
 * A vector.
 *
 * @type {Vector3}
 * @private
 */

const v = new Vector3();

/**
 * A frustum.
 */

export class Frustum {

	/**
	 * Constructs a new frustum.
	 *
	 * @param {Plane} [p0] - A plane.
	 * @param {Plane} [p1] - A plane.
	 * @param {Plane} [p2] - A plane.
	 * @param {Plane} [p3] - A plane.
	 * @param {Plane} [p4] - A plane.
	 * @param {Plane} [p5] - A plane.
	 */

	constructor(
		p0 = new Plane(),
		p1 = new Plane(),
		p2 = new Plane(),
		p3 = new Plane(),
		p4 = new Plane(),
		p5 = new Plane()
	) {

		/**
		 * The six planes that form the frustum.
		 *
		 * @type {Plane[]}
		 */

		this.planes = [p0, p1, p2, p3, p4, p5];

	}

	/**
	 * Sets the planes of this frustum.
	 *
	 * @param {Plane} [p0] - A plane.
	 * @param {Plane} [p1] - A plane.
	 * @param {Plane} [p2] - A plane.
	 * @param {Plane} [p3] - A plane.
	 * @param {Plane} [p4] - A plane.
	 * @param {Plane} [p5] - A plane.
	 * @return {Frustum} This frustum.
	 */

	set(p0, p1, p2, p3, p4, p5) {

		const planes = this.planes;

		planes[0].copy(p0);
		planes[1].copy(p1);
		planes[2].copy(p2);
		planes[3].copy(p3);
		planes[4].copy(p4);
		planes[5].copy(p5);

		return this;

	}

	/**
	 * Clones this frustum.
	 *
	 * @return {Frustum} The cloned frustum.
	 */

	clone() {

		return new this.constructor().copy(this);

	}

	/**
	 * Copies a given frustum.
	 *
	 * @param {Frustum} frustum - A frustum.
	 * @return {Frustum} This frustum.
	 */

	copy(frustum) {

		const planes = this.planes;

		let i;

		for(i = 0; i < 6; ++i) {

			planes[i].copy(frustum.planes[i]);

		}

		return this;

	}

	/**
	 * Sets this frustm based on a given 4x4 matrix.
	 *
	 * @param {Matrix4} m - A matrix.
	 * @return {Frustum} This frustum.
	 */

	setFromMatrix(m) {

		const planes = this.planes;

		const me = m.elements;
		const me0 = me[0], me1 = me[1], me2 = me[2], me3 = me[3];
		const me4 = me[4], me5 = me[5], me6 = me[6], me7 = me[7];
		const me8 = me[8], me9 = me[9], me10 = me[10], me11 = me[11];
		const me12 = me[12], me13 = me[13], me14 = me[14], me15 = me[15];

		planes[0].setComponents(me3 - me0, me7 - me4, me11 - me8, me15 - me12).normalize();
		planes[1].setComponents(me3 + me0, me7 + me4, me11 + me8, me15 + me12).normalize();
		planes[2].setComponents(me3 + me1, me7 + me5, me11 + me9, me15 + me13).normalize();
		planes[3].setComponents(me3 - me1, me7 - me5, me11 - me9, me15 - me13).normalize();
		planes[4].setComponents(me3 - me2, me7 - me6, me11 - me10, me15 - me14).normalize();
		planes[5].setComponents(me3 + me2, me7 + me6, me11 + me10, me15 + me14).normalize();

		return this;

	}

	/**
	 * Checks if this frustum intersects with the given sphere.
	 *
	 * @param {Sphere} sphere - A sphere.
	 * @return {Boolean} Whether this frustum intersects with the sphere.
	 */

	intersectsSphere(sphere) {

		const planes = this.planes;
		const center = sphere.center;
		const negativeRadius = -sphere.radius;

		let result = true;
		let i, d;

		for(i = 0; i < 6; ++i) {

			d = planes[i].distanceToPoint(center);

			if(d < negativeRadius) {

				result = false;
				break;

			}

		}

		return result;

	}

	/**
	 * Checks if this frustum intersects with the given sphere.
	 *
	 * @param {Box3} box - A box.
	 * @return {Boolean} Whether this frustum intersects with the box.
	 */

	intersectsBox(box) {

		const planes = this.planes;
		const min = box.min, max = box.max;

		let i, plane;

		for(i = 0; i < 6; ++i) {

			plane = planes[i];

			// Corner at max distance.
			v.x = (plane.normal.x > 0.0) ? max.x : min.x;
			v.y = (plane.normal.y > 0.0) ? max.y : min.y;
			v.z = (plane.normal.z > 0.0) ? max.z : min.z;

			if(plane.distanceToPoint(v) < 0.0) {

				return false;

			}

		}

		return true;

	}

	/**
	 * Checks if this frustum contains the given point.
	 *
	 * @param {Vector3} point - A point.
	 * @return {Boolean} Whether this frustum contains the point.
	 */

	containsPoint(point) {

		const planes = this.planes;

		let result = true;
		let i;

		for(i = 0; i < 6; ++i) {

			if(planes[i].distanceToPoint(point) < 0) {

				result = false;
				break;

			}

		}

		return result;

	}

}
