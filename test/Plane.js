import test from "ava";
import { Plane } from "../build/math-ds.js";

test("can be created", t => {

	const object = new Plane();

	t.truthy(object);

});
