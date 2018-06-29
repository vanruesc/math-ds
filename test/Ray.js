import test from "ava";
import { Ray } from "../build/math-ds.js";

test("can be created", t => {

	const object = new Ray();

	t.truthy(object);

});
