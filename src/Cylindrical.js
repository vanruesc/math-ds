/**
 * A cylindrical coordinate system.
 *
 * For details see: https://en.wikipedia.org/wiki/Cylindrical_coordinate_system
 */

export class Cylindrical {

	/**
	 * Constructs a new cylindrical system.
	 *
	 * @param {Number} [radius=1] - The radius of the cylinder.
	 * @param {Number} [theta=0] - Theta.
	 * @param {Number} [y=0] - The height.
	 */

	constructor(radius = 1, theta = 0, y = 0) {

		/**
		 * The distance from the origin to a point in the XZ-plane.
		 *
		 * @type {Number}
		 */

		this.radius = radius;

		/**
		 * The counterclockwise angle in the XZ-plane measured in radians from the
		 * positive Z-axis.
		 *
		 * @type {Number}
		 */

		this.theta = theta;

		/**
		 * The height above the XZ-plane.
		 *
		 * @type {Number}
		 */

		this.y = y;

	}

	/**
	 * Sets the values of this cylindrical system.
	 *
	 * @param {Number} radius - The radius.
	 * @param {Number} theta - Theta.
	 * @param {Number} y - The height.
	 * @return {Cylindrical} This cylindrical system.
	 */

	set(radius, theta, y) {

		this.radius = radius;
		this.theta = theta;
		this.y = y;

		return this;

	}

	/**
	 * Copies the values of the given cylindrical system.
	 *
	 * @param {Cylindrical} c - A cylindrical system.
	 * @return {Cylindrical} This cylindrical system.
	 */

	copy(c) {

		this.radius = c.radius;
		this.theta = c.theta;
		this.y = c.y;

		return this;

	}

	/**
	 * Clones this cylindrical system.
	 *
	 * @return {Cylindrical} The cloned cylindrical system.
	 */

	clone() {

		return new this.constructor().copy(this);

	}

	/**
	 * Sets the values of this cylindrical system.
	 *
	 * @param {Vector3} v - The vector
	 * @return {Cylindrical} This cylindrical system.
	 */

	setFromVector3(v) {

		this.radius = Math.sqrt(v.x * v.x + v.z * v.z);
		this.theta = Math.atan2(v.x, v.z);
		this.y = v.y;

		return this;

	}

}
