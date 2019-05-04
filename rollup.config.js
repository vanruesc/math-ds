import babel from "rollup-plugin-babel";
import minify from "rollup-plugin-babel-minify";
import resolve from "rollup-plugin-node-resolve";

const pkg = require("./package.json");
const date = (new Date()).toDateString();

const banner = `/**
 * ${pkg.name} v${pkg.version} build ${date}
 * ${pkg.homepage}
 * Copyright ${date.slice(-4)} ${pkg.author.name}, ${pkg.license}
 */`;

const production = (process.env.NODE_ENV === "production");

const lib = {

	esm: {

		input: "src/index.js",
		plugins: [resolve()],
		output: [{
			file: pkg.main,
			format: "esm"
		}].concat(!production ? [] : [
			{
				file: pkg.module,
				format: "esm",
				banner: banner
			}, {
				file: pkg.main.replace(".js", ".min.js"),
				format: "esm"
			}
		])

	},

	umd: {

		input: pkg.main,
		plugins: production ? [babel()] : [],
		output: {
			file: pkg.main,
			format: "umd",
			name: pkg.name.replace(/-/g, "").toUpperCase(),
			banner: banner
		}

	}

};

export default [lib.esm, lib.umd].concat(production ? [{

	input: lib.esm.output[2].file,
	plugins: [babel(), minify({
		bannerNewLine: true,
		comments: false
	})],
	output: {
		file: lib.esm.output[2].file,
		format: "umd",
		name: pkg.name.replace(/-/g, "").toUpperCase(),
		banner: banner
	}

}] : []);

