import test from "ava";
import { Spherical } from "../build/math-ds.js";

test("can be created", t => {

	const object = new Spherical();

	t.truthy(object);

});
