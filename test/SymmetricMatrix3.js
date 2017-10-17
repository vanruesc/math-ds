"use strict";

const lib = require("../build/math-ds.js");

const SymmetricMatrix3 = lib.SymmetricMatrix3;
const Matrix3 = lib.Matrix3;
const Vector3 = lib.Vector3;

module.exports = {

	"SymmetricMatrix3": {

		"can be instantiated": function(test) {

			const sm = new SymmetricMatrix3();

			test.ok(sm);
			test.done();

		},

		"correctly transforms vectors": function(test) {

			const m0 = new SymmetricMatrix3();
			const m1 = new Matrix3();

			const v0 = new Vector3(0, 1, 2);
			const v1 = v0.clone();
			const v2 = v0.clone();
			const v3 = v0.clone();

			m0.identity();
			m1.identity();

			m0.applyToVector3(v0);
			v1.applyMatrix3(m1);

			m0.set(

				1, 2, 3,
				2, 3,
				3

			);

			m1.set(

				1, 2, 3,
				2, 2, 3,
				3, 3, 3

			);

			m0.applyToVector3(v2);
			v3.applyMatrix3(m1);

			test.ok(v0.equals(v1), "should compute the same result as its complete matrix equivalent");
			test.ok(v2.equals(v3), "should compute the same result as its complete matrix equivalent");
			test.done();

		}

	}

};
