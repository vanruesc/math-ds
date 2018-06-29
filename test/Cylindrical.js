import test from "ava";
import { Cylindrical } from "../build/math-ds.js";

test("can be created", t => {

	const object = new Cylindrical();

	t.truthy(object);

});
