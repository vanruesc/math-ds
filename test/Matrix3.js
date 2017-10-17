"use strict";

const Matrix3 = require("../build/math-ds.js").Matrix3;

module.exports = {

	"Matrix3": {

		"can be instantiated": function(test) {

			const m = new Matrix3();

			test.ok(m);
			test.done();

		}

	}

};
