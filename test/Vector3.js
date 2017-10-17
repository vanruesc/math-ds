"use strict";

const Vector3 = require("../build/math-ds.js").Vector3;

module.exports = {

	"Vector3": {

		"can be instantiated": function(test) {

			const v = new Vector3();

			test.ok(v);
			test.done();

		}

	}

};
