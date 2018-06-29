import test from "ava";
import { Quaternion } from "../build/math-ds.js";

test("can be created", t => {

	const object = new Quaternion();

	t.truthy(object);

});
