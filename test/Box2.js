import test from "ava";
import { Box2 } from "../build/math-ds.js";

test("can be created", t => {

	const object = new Box2();

	t.truthy(object);

});
