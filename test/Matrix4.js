import test from "ava";
import { Matrix4 } from "../build/math-ds.js";

test("can be created", t => {

	const object = new Matrix4();

	t.truthy(object);

});
