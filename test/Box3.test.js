"use strict";

const Box3 = require("../build/math-ds.js").Box3;

module.exports = {

	"Box3": {

		"can be instantiated": function(test) {

			const b = new Box3();

			test.ok(b);
			test.done();

		}

	}

};
