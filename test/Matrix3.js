import test from "ava";
import { Matrix3 } from "../build/math-ds.js";

test("can be created", t => {

	const object = new Matrix3();

	t.truthy(object);

});
