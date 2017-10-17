"use strict";

const Plane = require("../build/math-ds.js").Plane;

module.exports = {

	"Plane": {

		"can be instantiated": function(test) {

			const p = new Plane();

			test.ok(p);
			test.done();

		}

	}

};
