"use strict";

const Box2 = require("../build/math-ds.js").Box2;

module.exports = {

	"Box2": {

		"can be instantiated": function(test) {

			const b = new Box2();

			test.ok(b);
			test.done();

		}

	}

};
