import test from "ava";
import { Matrix3, SymmetricMatrix3, Vector3 } from "../build/math-ds.js";

test("can be created", t => {

	const object = new SymmetricMatrix3();

	t.truthy(object);

});

test("correctly transforms vectors", t => {

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

	t.true(v0.equals(v1), "should compute the same result as its complete matrix equivalent");
	t.true(v2.equals(v3), "should compute the same result as its complete matrix equivalent");

});
