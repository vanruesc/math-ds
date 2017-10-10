import { Vector3 } from "./Vector3.js";

/**
 * A list of vectors.
 *
 * @type {Vector3[]}
 * @private
 */

const v = [
	new Vector3(),
	new Vector3(),
	new Vector3(),
	new Vector3()
];

/**
 * A ray.
 */

export class Ray {

	/**
	 * Constructs a new ray.
	 *
	 * @param {Vector3} [origin] - The origin.
	 * @param {Vector3} [direction] - The direction.
	 */

	constructor(origin = new Vector3(), direction = new Vector3()) {

		/**
		 * The origin.
		 *
		 * @type {Vector3}
		 */

		this.origin = origin;

		/**
		 * The direction.
		 *
		 * @type {Vector3}
		 */

		this.direction = direction;

	}

	/**
	 * Sets the origin and the direction.
	 *
	 * @param {Vector3} origin - The origin.
	 * @param {Vector3} direction - The direction. Should be normalized.
	 * @return {Ray} This ray.
	 */

	set(origin, direction) {

		this.origin.copy(origin);
		this.direction.copy(direction);

		return this;

	}

	/**
	 * Copies the given ray.
	 *
	 * @param {Ray} r - A ray.
	 * @return {Ray} This ray.
	 */

	copy(r) {

		this.origin.copy(r.origin);
		this.direction.copy(r.direction);

		return this;

	}

	/**
	 * Clones this ray.
	 *
	 * @return {Ray} The cloned ray.
	 */

	clone() {

		return new this.constructor().copy(this);

	}

	/**
	 * Computes a point along the ray based on a given scalar t.
	 *
	 * @param {Number} t - The scalar.
	 * @param {Vector3} [target] - A target vector. If none is provided, a new one will be created.
	 * @return {Vector3} The point.
	 */

	at(t, target = new Vector3()) {

		return target.copy(this.direction).multiplyScalar(t).add(this.origin);

	}

	/**
	 * Rotates this ray to look at the given target.
	 *
	 * @param {Vector3} target - A point to look at.
	 * @return {Ray} This ray.
	 */

	lookAt(target) {

		this.direction.copy(target).sub(this.origin).normalize();

		return this;

	}

	/**
	 * Moves the origin along the ray by a given scalar t.
	 *
	 * @param {Number} t - The scalar.
	 * @return {Ray} This ray.
	 */

	recast(t) {

		this.origin.copy(this.at(t, v[0]));

		return this;

	}

	/**
	 * Finds the closest point along this ray to a given point.
	 *
	 * @param {Vector3} p - A point.
	 * @param {Vector3} [target] - A target vector. If none is provided, a new one will be created.
	 * @return {Vector3} The point.
	 */

	closestPointToPoint(p, target = new Vector3()) {

		const directionDistance = target.subVectors(p, this.origin).dot(this.direction);

		if(directionDistance < 0.0) {

			return target.copy(this.origin);

		}

		return target.copy(this.direction).multiplyScalar(directionDistance).add(this.origin);

	}

	/**
	 * Calculates the squared distance from this ray to the given point.
	 *
	 * @param {Vector3} p - The point.
	 * @return {Number} The squared distance.
	 */

	distanceSquaredToPoint(p) {

		const directionDistance = v[0].subVectors(p, this.origin).dot(this.direction);

		// Check if the point is behind the ray.
		return (directionDistance < 0.0) ?
			this.origin.distanceToSquared(p) :
			v[0].copy(this.direction).multiplyScalar(directionDistance).add(this.origin).distanceToSquared(p);

	}

	/**
	 * Calculates the distance from this ray to the given point.
	 *
	 * @param {Vector3} p - The point.
	 * @return {Number} The distance.
	 */

	distanceToPoint(p) {

		return Math.sqrt(this.distanceSquaredToPoint(p));

	}

	/**
	 * Calculates the distance from this ray to the given plane.
	 *
	 * @param {Vector3} p - The point.
	 * @return {Number} The distance, or null if the denominator is zero.
	 */

	distanceToPlane(p) {

		const denominator = p.normal.dot(this.direction);

		const t = (denominator !== 0.0) ?
			-(this.origin.dot(p.normal) + p.constant) / denominator :
			((p.distanceToPoint(this.origin) === 0.0) ? 0.0 : -1.0);

		return (t >= 0.0) ? t : null;

	}

	/**
	 * Calculates the distance from this ray to a given line segment.
	 *
	 * Based on:
	 *  http://www.geometrictools.com/GTEngine/Include/Mathematics/GteDistRaySegment.h
	 *
	 * @param {Vector3} v0 - The start of the segment.
	 * @param {Vector3} v1 - The end of the segment.
	 * @param {Vector3} [pointOnRay] - If provided, the point on this Ray that is closest to the segment will be stored in this vector.
	 * @param {Vector3} [pointOnSegment] - If provided, the point on the line segment that is closest to this ray will be stored in this vector.
	 * @return {Number} The smallest distance between the ray and the segment defined by v0 and v1.
	 */

	distanceSquaredToSegment(v0, v1, pointOnRay, pointOnSegment) {

		const segCenter = v[0].copy(v0).add(v1).multiplyScalar(0.5);
		const segDir = v[1].copy(v1).sub(v0).normalize();
		const diff = v[2].copy(this.origin).sub(segCenter);

		const segExtent = v0.distanceTo(v1) * 0.5;
		const a01 = -this.direction.dot(segDir);
		const b0 = diff.dot(this.direction);
		const b1 = -diff.dot(segDir);
		const c = diff.lengthSq();
		const det = Math.abs(1.0 - a01 * a01);

		let s0, s1, extDet, invDet, sqrDist;

		if(det > 0.0) {

			// The ray and segment are not parallel.
			s0 = a01 * b1 - b0;
			s1 = a01 * b0 - b1;
			extDet = segExtent * det;

			if(s0 >= 0.0) {

				if(s1 >= -extDet) {

					if(s1 <= extDet) {

						// Region 0.
						// Minimum at interior points of ray and segment.
						invDet = 1.0 / det;
						s0 *= invDet;
						s1 *= invDet;
						sqrDist = s0 * (s0 + a01 * s1 + 2.0 * b0) + s1 * (a01 * s0 + s1 + 2.0 * b1) + c;

					} else {

						// Region 1.
						s1 = segExtent;
						s0 = Math.max(0.0, -(a01 * s1 + b0));
						sqrDist = -s0 * s0 + s1 * (s1 + 2.0 * b1) + c;

					}

				} else {

					// Region 5.
					s1 = -segExtent;
					s0 = Math.max(0.0, -(a01 * s1 + b0));
					sqrDist = -s0 * s0 + s1 * (s1 + 2.0 * b1) + c;

				}

			} else {

				if(s1 <= -extDet) {

					// Region 4.
					s0 = Math.max(0.0, -(-a01 * segExtent + b0));
					s1 = (s0 > 0.0) ? -segExtent : Math.min(Math.max(-segExtent, -b1), segExtent);
					sqrDist = -s0 * s0 + s1 * (s1 + 2.0 * b1) + c;

				} else if(s1 <= extDet) {

					// Region 3.
					s0 = 0.0;
					s1 = Math.min(Math.max(-segExtent, -b1), segExtent);
					sqrDist = s1 * (s1 + 2.0 * b1) + c;

				} else {

					// Region 2.
					s0 = Math.max(0.0, -(a01 * segExtent + b0));
					s1 = (s0 > 0.0) ? segExtent : Math.min(Math.max(-segExtent, -b1), segExtent);
					sqrDist = -s0 * s0 + s1 * (s1 + 2.0 * b1) + c;

				}

			}

		} else {

			// Ray and segment are parallel.
			s1 = (a01 > 0.0) ? -segExtent : segExtent;
			s0 = Math.max(0.0, -(a01 * s1 + b0));
			sqrDist = -s0 * s0 + s1 * (s1 + 2.0 * b1) + c;

		}

		if(pointOnRay !== undefined) {

			pointOnRay.copy(this.direction).multiplyScalar(s0).add(this.origin);

		}

		if(pointOnSegment !== undefined) {

			pointOnSegment.copy(segDir).multiplyScalar(s1).add(segCenter);

		}

		return sqrDist;

	}

	/**
	 * Finds the point where this ray intersects the given sphere.
	 *
	 * @param {Sphere} s - A sphere.
	 * @param {Vector3} [target] - A target vector. If none is provided, a new one will be created.
	 * @return {Vector3} The point of intersection, or null if there is none.
	 */

	intersectSphere(s, target = new Vector3()) {

		const ab = v[0].subVectors(s.center, this.origin);
		const tca = ab.dot(this.direction);
		const d2 = ab.dot(ab) - tca * tca;
		const radius2 = s.radius * s.radius;

		let result = null;
		let thc, t0, t1;

		if(d2 <= radius2) {

			thc = Math.sqrt(radius2 - d2);

			// t0 = first intersection point - entrance on front of sphere.
			t0 = tca - thc;

			// t1 = second intersection point - exit point on back of sphere.
			t1 = tca + thc;

			// Check if both t0 and t1 are behind the ray - if so, return null.
			if(t0 >= 0.0 || t1 >= 0.0) {

				/* Check if t0 is behind the ray. If it is, the ray is inside the
				sphere, so return the second exit point scaled by t1 in order to always
				return an intersection point that is in front of the ray. If t0 is in
				front of the ray, return the first collision point scaled by t0. */
				result = (t0 < 0.0) ? this.at(t1, target) : this.at(t0, target);

			}

		}

		return result;

	}

	/**
	 * Determines whether this ray intersects the given sphere.
	 *
	 * @param {Sphere} s - A sphere.
	 * @return {Boolean} Whether this ray intersects the given sphere.
	 */

	intersectsSphere(s) {

		return (this.distanceToPoint(s.center) <= s.radius);

	}

	/**
	 * Finds the point where this ray intersects the given plane.
	 *
	 * @param {Plane} p - A plane.
	 * @param {Vector3} [target] - A target vector. If none is provided, a new one will be created.
	 * @return {Vector3} The point of intersection, or null if there is none.
	 */

	intersectPlane(p, target = new Vector3()) {

		const t = this.distanceToPlane(p);

		return (t === null) ? null : this.at(t, target);

	}

	/**
	 * Determines whether this ray intersects the given plane.
	 *
	 * @param {Plane} p - A plane.
	 * @return {Boolean} Whether this ray intersects the given plane.
	 */

	intersectsPlane(p) {

		const distanceToPoint = p.distanceToPoint(this.origin);

		return (distanceToPoint === 0.0 || p.normal.dot(this.direction) * distanceToPoint < 0.0);

	}

	/**
	 * Finds the point where this ray intersects the given box.
	 *
	 * @param {Plane} b - A box.
	 * @param {Vector3} [target] - A target vector. If none is provided, a new one will be created.
	 * @return {Vector3} The point of intersection, or null if there is none.
	 */

	intersectBox(b, target = new Vector3()) {

		const origin = this.origin;
		const direction = this.direction;
		const min = b.min;
		const max = b.max;

		const invDirX = 1.0 / direction.x;
		const invDirY = 1.0 / direction.y;
		const invDirZ = 1.0 / direction.z;

		let result = null;
		let tmin, tmax, tymin, tymax, tzmin, tzmax;

		if(invDirX >= 0.0) {

			tmin = (min.x - origin.x) * invDirX;
			tmax = (max.x - origin.x) * invDirX;

		} else {

			tmin = (max.x - origin.x) * invDirX;
			tmax = (min.x - origin.x) * invDirX;

		}

		if(invDirY >= 0.0) {

			tymin = (min.y - origin.y) * invDirY;
			tymax = (max.y - origin.y) * invDirY;

		} else {

			tymin = (max.y - origin.y) * invDirY;
			tymax = (min.y - origin.y) * invDirY;

		}

		if(tmin <= tymax && tymin <= tmax) {

			/* Handle the case where tmin or tmax is NaN (result of 0 * Infinity).
			Note: x !== x returns true if x is NaN. */
			if(tymin > tmin || tmin !== tmin) { tmin = tymin; }
			if(tymax < tmax || tmax !== tmax) { tmax = tymax; }

			if(invDirZ >= 0.0) {

				tzmin = (min.z - origin.z) * invDirZ;
				tzmax = (max.z - origin.z) * invDirZ;

			} else {

				tzmin = (max.z - origin.z) * invDirZ;
				tzmax = (min.z - origin.z) * invDirZ;

			}

			if(tmin <= tzmax && tzmin <= tmax) {

				if(tzmin > tmin || tmin !== tmin) { tmin = tzmin; }
				if(tzmax < tmax || tmax !== tmax) { tmax = tzmax; }

				// Return the closest point (positive side).
				if(tmax >= 0.0) {

					result = this.at((tmin >= 0.0) ? tmin : tmax, target);

				}

			}

		}

		return result;

	}

	/**
	 * Determines whether this ray intersects the given box.
	 *
	 * @param {Box3} b - A box.
	 * @return {Boolean} Whether this ray intersects the given box.
	 */

	intersectsBox(b) {

		return (this.intersectBox(b, v[0]) !== null);

	}

	/**
	 * Finds the point where this ray intersects the given triangle.
	 *
	 * Based on:
	 *  http://www.geometrictools.com/GTEngine/Include/Mathematics/GteIntrRay3Triangle3.h
	 *
	 * @param {Vector3} a - A triangle vertex.
	 * @param {Vector3} b - A triangle vertex.
	 * @param {Vector3} c - A triangle vertex.
	 * @param {Boolean} [backfaceCulling=false] - Whether backface culling should be considered.
	 * @param {Vector3} [target] - A target vector. If none is provided, a new one will be created.
	 * @return {Vector3} The point of intersection, or null if there is none.
	 */

	intersectTriangle(a, b, c, backfaceCulling, target) {

		const direction = this.direction;

		// Compute the offset origin, edges, and normal.
		const diff = v[0];
		const edge1 = v[1];
		const edge2 = v[2];
		const normal = v[3];

		let result = null;
		let DdN, sign, DdQxE2, DdE1xQ, QdN;

		edge1.subVectors(b, a);
		edge2.subVectors(c, a);
		normal.crossVectors(edge1, edge2);

		/* Solve Q + t * D = b1 * E1 + b2 * E2
		 * (Q = kDiff, D = ray direction, E1 = kEdge1, E2 = kEdge2,
		 * N = Cross(E1, E2)):
		 *
		 *   | Dot(D, N) | * b1 = sign(Dot(D, N)) * Dot(D, Cross(Q, E2))
		 *   | Dot(D, N) | * b2 = sign(Dot(D, N)) * Dot(D, Cross(E1, Q))
		 *   | Dot(D, N) | * t = -sign(Dot(D, N)) * Dot(Q, N)
		 */

		DdN = direction.dot(normal);

		// Discard coplanar constellations and cull backfaces.
		if(DdN !== 0.0 && !(backfaceCulling && DdN > 0.0)) {

			if(DdN > 0.0) {

				sign = 1.0;

			} else {

				sign = -1.0;
				DdN = -DdN;

			}

			diff.subVectors(this.origin, a);
			DdQxE2 = sign * direction.dot(edge2.crossVectors(diff, edge2));

			// b1 < 0, no intersection.
			if(DdQxE2 >= 0.0) {

				DdE1xQ = sign * direction.dot(edge1.cross(diff));

				// b2 < 0, or b1 + b2 > 1, no intersection.
				if(DdE1xQ >= 0.0 && DdQxE2 + DdE1xQ <= DdN) {

					// The line intersects the triangle, check if the ray does.
					QdN = -sign * diff.dot(normal);

					// t < 0, no intersection.
					if(QdN >= 0.0) {

						// Ray intersects triangle.
						result = this.at(QdN / DdN, target);

					}

				}

			}

		}

		return result;

	}

	/**
	 * Applies the given matrix to this ray.
	 *
	 * @param {Matrix4} m - A matrix.
	 * @return {Ray} This ray.
	 */

	applyMatrix4(m) {

		this.origin.applyMatrix4(m);
		this.direction.transformDirection(m);

		return this;

	}

	/**
	 * Checks if this ray equals the given one.
	 *
	 * @param {Ray} r - A ray.
	 * @return {Boolean} Whether the rays are equal.
	 */

	equals(r) {

		return (r.origin.equals(this.origin) && r.direction.equals(this.direction));

	}

}
