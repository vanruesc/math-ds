"use strict";

const Cylindrical = require("../build/math-ds.js").Cylindrical;

module.exports = {

	"Cylindrical": {

		"can be instantiated": function(test) {

			const c = new Cylindrical();

			test.ok(c);
			test.done();

		}

	}

};
