import test from "ava";
import { Box3 } from "../build/math-ds.js";

test("can be created", t => {

	const object = new Box3();

	t.truthy(object);

});
