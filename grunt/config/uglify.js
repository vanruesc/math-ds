module.exports = {

	lib: {
		options: {
			banner: "<%= banner %>"
		},
		files: {
			"build/<%= package.name %>.min.js": ["build/<%= package.name %>.js"]
		}
	}

};
