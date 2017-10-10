"use strict";

const Line3 = require("../build/math-ds.js").Line3;

module.exports = {

	"Line3": {

		"can be instantiated": function(test) {

			const l = new Line3();

			test.ok(l);
			test.done();

		}

	}

};
