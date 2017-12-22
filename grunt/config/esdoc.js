module.exports = {

	compile: {
		options: {
			source: "src",
			destination: "docs",
			plugins: [{
				name: "esdoc-standard-plugin"
			}]
		}
	}

};
