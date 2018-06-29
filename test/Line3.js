import test from "ava";
import { Line3 } from "../build/math-ds.js";

test("can be created", t => {

	const object = new Line3();

	t.truthy(object);

});
