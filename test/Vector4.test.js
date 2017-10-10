"use strict";

const Vector4 = require("../build/math-ds.js").Vector4;

module.exports = {

	"Vector4": {

		"can be instantiated": function(test) {

			const v = new Vector4();

			test.ok(v);
			test.done();

		}

	}

};
