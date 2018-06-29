import test from "ava";
import { Frustum } from "../build/math-ds.js";

test("can be created", t => {

	const object = new Frustum();

	t.truthy(object);

});
