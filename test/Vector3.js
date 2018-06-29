import test from "ava";
import { Vector3 } from "../build/math-ds.js";

test("can be created", t => {

	const object = new Vector3();

	t.truthy(object);

});
