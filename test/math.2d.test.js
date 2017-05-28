"use strict";

const lib = require("../build/math-ds.js");

module.exports = {

	"Vector2": {

		"can be instantiated": function(test) {

			const v = new lib.Vector2();

			test.ok(v);
			test.done();

		}

	}

};
