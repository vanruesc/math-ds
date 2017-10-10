"use strict";

const Spherical = require("../build/math-ds.js").Spherical;

module.exports = {

	"Spherical": {

		"can be instantiated": function(test) {

			const s = new Spherical();

			test.ok(s);
			test.done();

		}

	}

};
