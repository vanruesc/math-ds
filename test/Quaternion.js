"use strict";

const Quaternion = require("../build/math-ds.js").Quaternion;

module.exports = {

	"Quaternion": {

		"can be instantiated": function(test) {

			const q = new Quaternion();

			test.ok(q);
			test.done();

		}

	}

};
