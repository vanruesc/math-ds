"use strict";

const Euler = require("../build/math-ds.js").Euler;

module.exports = {

	"Euler": {

		"can be instantiated": function(test) {

			const e = new Euler();

			test.ok(e);
			test.done();

		}

	}

};
