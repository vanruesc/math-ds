"use strict";

const Matrix4 = require("../build/math-ds.js").Matrix4;

module.exports = {

	"Matrix4": {

		"can be instantiated": function(test) {

			const m = new Matrix4();

			test.ok(m);
			test.done();

		}

	}

};
