"use strict";

const lib = require("../build/math-ds.js");

module.exports = {

	"Box3": {

		"can be instantiated": function(test) {

			const b = new lib.Box3();

			test.ok(b);
			test.done();

		}

	},

	"Matrix3": {

		"can be instantiated": function(test) {

			const m = new lib.Matrix3();

			test.ok(m);
			test.done();

		}

	},

	"SymmetricMatrix3": {

		"can be instantiated": function(test) {

			const sm = new lib.SymmetricMatrix3();

			test.ok(sm);
			test.done();

		}

	},

	"Vector3": {

		"can be instantiated": function(test) {

			const v = new lib.Vector3();

			test.ok(v);
			test.done();

		}

	}

};
