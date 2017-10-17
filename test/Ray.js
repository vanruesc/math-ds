"use strict";

const Ray = require("../build/math-ds.js").Ray;

module.exports = {

	"Ray": {

		"can be instantiated": function(test) {

			const r = new Ray();

			test.ok(r);
			test.done();

		}

	}

};
