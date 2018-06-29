import test from "ava";
import { Vector4 } from "../build/math-ds.js";

test("can be created", t => {

	const object = new Vector4();

	t.truthy(object);

});
