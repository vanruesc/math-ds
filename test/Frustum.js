"use strict";

const Frustum = require("../build/math-ds.js").Frustum;

module.exports = {

	"Frustum": {

		"can be instantiated": function(test) {

			const l = new Frustum();

			test.ok(l);
			test.done();

		}

	}

};
