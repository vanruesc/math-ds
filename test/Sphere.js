import test from "ava";
import { Sphere } from "../build/math-ds.js";

test("can be created", t => {

	const object = new Sphere();

	t.truthy(object);

});
