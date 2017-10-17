"use strict";

const Sphere = require("../build/math-ds.js").Sphere;

module.exports = {

	"Sphere": {

		"can be instantiated": function(test) {

			const s = new Sphere();

			test.ok(s);
			test.done();

		}

	}

};
