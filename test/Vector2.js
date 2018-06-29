import test from "ava";
import { Vector2 } from "../build/math-ds.js";

test("can be created", t => {

	const object = new Vector2();

	t.truthy(object);

});
