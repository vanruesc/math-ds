"use strict";

const Vector2 = require("../build/math-ds.js").Vector2;

module.exports = {

	"Vector2": {

		"can be instantiated": function(test) {

			const v = new Vector2();

			test.ok(v);
			test.done();

		}

	}

};
