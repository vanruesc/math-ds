import test from "ava";
import { Euler } from "../build/math-ds.js";

test("can be created", t => {

	const object = new Euler();

	t.truthy(object);

});
