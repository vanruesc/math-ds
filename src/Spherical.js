/**
 * A spherical coordinate system.
 *
 * For details see: https://en.wikipedia.org/wiki/Spherical_coordinate_system
 *
 * The poles (phi) are at the positive and negative Y-axis. The equator starts
 * at positive Z.
 */

export class Spherical {

	/**
	 * Constructs a new spherical system.
	 *
	 * @param {Number} [radius=1] - The radius of the sphere.
	 * @param {Number} [phi=0] - The polar angle phi.
	 * @param {Number} [theta=0] - The equator angle theta.
	 */

	constructor(radius = 1, phi = 0, theta = 0) {

		/**
		 * The radius of the sphere.
		 *
		 * @type {Number}
		 */

		this.radius = radius;

		/**
		 * The polar angle, up and down towards the top and bottom pole.
		 *
		 * @type {Number}
		 */

		this.phi = phi;

		/**
		 * The angle around the equator of the sphere.
		 *
		 * @type {Number}
		 */

		this.theta = theta;

	}

	/**
	 * Sets the values of this spherical system.
	 *
	 * @param {Number} radius - The radius.
	 * @param {Number} phi - Phi.
	 * @param {Number} theta - Theta.
	 * @return {Spherical} This spherical system.
	 */

	set(radius, phi, theta) {

		this.radius = radius;
		this.phi = phi;
		this.theta = theta;

		return this;

	}

	/**
	 * Copies the values of the given spherical system.
	 *
	 * @param {Spherical} s - A spherical system.
	 * @return {Spherical} This spherical system.
	 */

	copy(s) {

		this.radius = s.radius;
		this.phi = s.phi;
		this.theta = s.theta;

		return this;

	}

	/**
	 * Clones this spherical system.
	 *
	 * @return {Spherical} The cloned spherical system.
	 */

	clone() {

		return new this.constructor().copy(this);

	}

	/**
	 * Sets the values of this spherical system.
	 *
	 * @param {Vector3} v - The vector
	 * @return {Spherical} This spherical system.
	 */

	setFromVector3(v) {

		this.radius = v.length();

		if(this.radius === 0) {

			this.theta = 0;
			this.phi = 0;

		} else {

			// Calculate the equator angle around the positive Y-axis.
			this.theta = Math.atan2(v.x, v.z);

			// Calculate the polar angle.
			this.phi = Math.acos(Math.min(Math.max(v.y / this.radius, -1), 1));

		}

		return this.makeSafe();

	}

	/**
	 * Restricts phi to [1e-6, PI - 1e-6].
	 *
	 * @return {Spherical} This spherical system.
	 */

	makeSafe() {

		this.phi = Math.max(1e-6, Math.min(Math.PI - 1e-6, this.phi));

		return this;

	}

}
