/**
 * math-ds v0.0.1 build Aug 02 2017
 * https://github.com/vanruesc/math-ds
 * Copyright 2017 Raoul van Rüschen, Zlib
 */

(function (global, factory) {
  typeof exports === 'object' && typeof module !== 'undefined' ? factory(exports) :
  typeof define === 'function' && define.amd ? define(['exports'], factory) :
  (factory((global.MATHDS = {})));
}(this, (function (exports) { 'use strict';

  var classCallCheck = function (instance, Constructor) {
    if (!(instance instanceof Constructor)) {
      throw new TypeError("Cannot call a class as a function");
    }
  };

  var createClass = function () {
    function defineProperties(target, props) {
      for (var i = 0; i < props.length; i++) {
        var descriptor = props[i];
        descriptor.enumerable = descriptor.enumerable || false;
        descriptor.configurable = true;
        if ("value" in descriptor) descriptor.writable = true;
        Object.defineProperty(target, descriptor.key, descriptor);
      }
    }

    return function (Constructor, protoProps, staticProps) {
      if (protoProps) defineProperties(Constructor.prototype, protoProps);
      if (staticProps) defineProperties(Constructor, staticProps);
      return Constructor;
    };
  }();

  var Vector2 = function () {
  	function Vector2() {
  		var x = arguments.length > 0 && arguments[0] !== undefined ? arguments[0] : 0;
  		var y = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : 0;
  		classCallCheck(this, Vector2);


  		this.x = x;

  		this.y = y;
  	}

  	createClass(Vector2, [{
  		key: "set",
  		value: function set$$1(x, y) {

  			this.x = x;
  			this.y = y;

  			return this;
  		}
  	}, {
  		key: "copy",
  		value: function copy(v) {

  			this.x = v.x;
  			this.y = v.y;

  			return this;
  		}
  	}, {
  		key: "fromArray",
  		value: function fromArray(array) {
  			var offset = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : 0;


  			this.x = array[offset];
  			this.y = array[offset + 1];

  			return this;
  		}
  	}, {
  		key: "toArray",
  		value: function toArray$$1() {
  			var array = arguments.length > 0 && arguments[0] !== undefined ? arguments[0] : [];
  			var offset = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : 0;


  			array[offset] = this.x;
  			array[offset + 1] = this.y;

  			return array;
  		}
  	}, {
  		key: "equals",
  		value: function equals(v) {

  			return v.x === this.x && v.y === this.y;
  		}
  	}, {
  		key: "clone",
  		value: function clone() {

  			return new this.constructor(this.x, this.y);
  		}
  	}, {
  		key: "add",
  		value: function add(v) {

  			this.x += v.x;
  			this.y += v.y;

  			return this;
  		}
  	}, {
  		key: "addScaledVector",
  		value: function addScaledVector(v, s) {

  			this.x += v.x * s;
  			this.y += v.y * s;

  			return this;
  		}
  	}, {
  		key: "addScalar",
  		value: function addScalar(s) {

  			this.x += s;
  			this.y += s;

  			return this;
  		}
  	}, {
  		key: "addVectors",
  		value: function addVectors(a, b) {

  			this.x = a.x + b.x;
  			this.y = a.y + b.y;

  			return this;
  		}
  	}, {
  		key: "sub",
  		value: function sub(v) {

  			this.x -= v.x;
  			this.y -= v.y;

  			return this;
  		}
  	}, {
  		key: "subScalar",
  		value: function subScalar(s) {

  			this.x -= s;
  			this.y -= s;

  			return this;
  		}
  	}, {
  		key: "subVectors",
  		value: function subVectors(a, b) {

  			this.x = a.x - b.x;
  			this.y = a.y - b.y;

  			return this;
  		}
  	}, {
  		key: "multiply",
  		value: function multiply(v) {

  			this.x *= v.x;
  			this.y *= v.y;

  			return this;
  		}
  	}, {
  		key: "multiplyScalar",
  		value: function multiplyScalar(s) {

  			if (isFinite(s)) {

  				this.x *= s;
  				this.y *= s;
  			} else {

  				this.x = 0;
  				this.y = 0;
  			}

  			return this;
  		}
  	}, {
  		key: "multiplyVectors",
  		value: function multiplyVectors(a, b) {

  			this.x = a.x * b.x;
  			this.y = a.y * b.y;

  			return this;
  		}
  	}, {
  		key: "divide",
  		value: function divide(v) {

  			this.x /= v.x;
  			this.y /= v.y;

  			return this;
  		}
  	}, {
  		key: "divideScalar",
  		value: function divideScalar(s) {

  			return this.multiplyScalar(1 / s);
  		}
  	}, {
  		key: "divideVectors",
  		value: function divideVectors(a, b) {

  			this.x = a.x / b.x;
  			this.y = a.y / b.y;

  			return this;
  		}
  	}, {
  		key: "negate",
  		value: function negate() {

  			this.x = -this.x;
  			this.y = -this.y;

  			return this;
  		}
  	}, {
  		key: "dot",
  		value: function dot(v) {

  			return this.x * v.x + this.y * v.y;
  		}
  	}, {
  		key: "lengthSq",
  		value: function lengthSq() {

  			return this.x * this.x + this.y * this.y;
  		}
  	}, {
  		key: "length",
  		value: function length() {

  			return Math.sqrt(this.x * this.x + this.y * this.y);
  		}
  	}, {
  		key: "distanceTo",
  		value: function distanceTo(v) {

  			return Math.sqrt(this.distanceToSquared(v));
  		}
  	}, {
  		key: "distanceToSquared",
  		value: function distanceToSquared(v) {

  			var dx = this.x - v.x;
  			var dy = this.y - v.y;

  			return dx * dx + dy * dy;
  		}
  	}, {
  		key: "normalize",
  		value: function normalize() {

  			return this.divideScalar(this.length());
  		}
  	}, {
  		key: "min",
  		value: function min(v) {

  			this.x = Math.min(this.x, v.x);
  			this.y = Math.min(this.y, v.y);

  			return this;
  		}
  	}, {
  		key: "max",
  		value: function max(v) {

  			this.x = Math.max(this.x, v.x);
  			this.y = Math.max(this.y, v.y);

  			return this;
  		}
  	}, {
  		key: "clamp",
  		value: function clamp(min, max) {

  			this.x = Math.max(min.x, Math.min(max.x, this.x));
  			this.y = Math.max(min.y, Math.min(max.y, this.y));

  			return this;
  		}
  	}]);
  	return Vector2;
  }();

  var Vector3 = function () {
  	function Vector3() {
  		var x = arguments.length > 0 && arguments[0] !== undefined ? arguments[0] : 0;
  		var y = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : 0;
  		var z = arguments.length > 2 && arguments[2] !== undefined ? arguments[2] : 0;
  		classCallCheck(this, Vector3);


  		this.x = x;

  		this.y = y;

  		this.z = z;
  	}

  	createClass(Vector3, [{
  		key: "set",
  		value: function set$$1(x, y, z) {

  			this.x = x;
  			this.y = y;
  			this.z = z;

  			return this;
  		}
  	}, {
  		key: "copy",
  		value: function copy(v) {

  			this.x = v.x;
  			this.y = v.y;
  			this.z = v.z;

  			return this;
  		}
  	}, {
  		key: "fromArray",
  		value: function fromArray(array) {
  			var offset = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : 0;


  			this.x = array[offset];
  			this.y = array[offset + 1];
  			this.z = array[offset + 2];

  			return this;
  		}
  	}, {
  		key: "toArray",
  		value: function toArray$$1() {
  			var array = arguments.length > 0 && arguments[0] !== undefined ? arguments[0] : [];
  			var offset = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : 0;


  			array[offset] = this.x;
  			array[offset + 1] = this.y;
  			array[offset + 2] = this.z;

  			return array;
  		}
  	}, {
  		key: "equals",
  		value: function equals(v) {

  			return v.x === this.x && v.y === this.y && v.z === this.z;
  		}
  	}, {
  		key: "clone",
  		value: function clone() {

  			return new this.constructor(this.x, this.y, this.z);
  		}
  	}, {
  		key: "add",
  		value: function add(v) {

  			this.x += v.x;
  			this.y += v.y;
  			this.z += v.z;

  			return this;
  		}
  	}, {
  		key: "addScaledVector",
  		value: function addScaledVector(v, s) {

  			this.x += v.x * s;
  			this.y += v.y * s;
  			this.z += v.z * s;

  			return this;
  		}
  	}, {
  		key: "addScalar",
  		value: function addScalar(s) {

  			this.x += s;
  			this.y += s;
  			this.z += s;

  			return this;
  		}
  	}, {
  		key: "addVectors",
  		value: function addVectors(a, b) {

  			this.x = a.x + b.x;
  			this.y = a.y + b.y;
  			this.z = a.z + b.z;

  			return this;
  		}
  	}, {
  		key: "sub",
  		value: function sub(v) {

  			this.x -= v.x;
  			this.y -= v.y;
  			this.z -= v.z;

  			return this;
  		}
  	}, {
  		key: "subScalar",
  		value: function subScalar(s) {

  			this.x -= s;
  			this.y -= s;
  			this.z -= s;

  			return this;
  		}
  	}, {
  		key: "subVectors",
  		value: function subVectors(a, b) {

  			this.x = a.x - b.x;
  			this.y = a.y - b.y;
  			this.z = a.z - b.z;

  			return this;
  		}
  	}, {
  		key: "multiply",
  		value: function multiply(v) {

  			this.x *= v.x;
  			this.y *= v.y;
  			this.z *= v.z;

  			return this;
  		}
  	}, {
  		key: "multiplyScalar",
  		value: function multiplyScalar(s) {

  			if (isFinite(s)) {

  				this.x *= s;
  				this.y *= s;
  				this.z *= s;
  			} else {

  				this.x = 0;
  				this.y = 0;
  				this.z = 0;
  			}

  			return this;
  		}
  	}, {
  		key: "multiplyVectors",
  		value: function multiplyVectors(a, b) {

  			this.x = a.x * b.x;
  			this.y = a.y * b.y;
  			this.z = a.z * b.z;

  			return this;
  		}
  	}, {
  		key: "divide",
  		value: function divide(v) {

  			this.x /= v.x;
  			this.y /= v.y;
  			this.z /= v.z;

  			return this;
  		}
  	}, {
  		key: "divideScalar",
  		value: function divideScalar(s) {

  			return this.multiplyScalar(1 / s);
  		}
  	}, {
  		key: "divideVectors",
  		value: function divideVectors(a, b) {

  			this.x = a.x / b.x;
  			this.y = a.y / b.y;
  			this.z = a.z / b.z;

  			return this;
  		}
  	}, {
  		key: "negate",
  		value: function negate() {

  			this.x = -this.x;
  			this.y = -this.y;
  			this.z = -this.z;

  			return this;
  		}
  	}, {
  		key: "dot",
  		value: function dot(v) {

  			return this.x * v.x + this.y * v.y + this.z * v.z;
  		}
  	}, {
  		key: "lengthSq",
  		value: function lengthSq() {

  			return this.x * this.x + this.y * this.y + this.z * this.z;
  		}
  	}, {
  		key: "length",
  		value: function length() {

  			return Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z);
  		}
  	}, {
  		key: "distanceTo",
  		value: function distanceTo(v) {

  			return Math.sqrt(this.distanceToSquared(v));
  		}
  	}, {
  		key: "distanceToSquared",
  		value: function distanceToSquared(v) {

  			var dx = this.x - v.x;
  			var dy = this.y - v.y;
  			var dz = this.z - v.z;

  			return dx * dx + dy * dy + dz * dz;
  		}
  	}, {
  		key: "normalize",
  		value: function normalize() {

  			return this.divideScalar(this.length());
  		}
  	}, {
  		key: "min",
  		value: function min(v) {

  			this.x = Math.min(this.x, v.x);
  			this.y = Math.min(this.y, v.y);
  			this.z = Math.min(this.z, v.z);

  			return this;
  		}
  	}, {
  		key: "max",
  		value: function max(v) {

  			this.x = Math.max(this.x, v.x);
  			this.y = Math.max(this.y, v.y);
  			this.z = Math.max(this.z, v.z);

  			return this;
  		}
  	}, {
  		key: "clamp",
  		value: function clamp(min, max) {

  			this.x = Math.max(min.x, Math.min(max.x, this.x));
  			this.y = Math.max(min.y, Math.min(max.y, this.y));
  			this.z = Math.max(min.z, Math.min(max.z, this.z));

  			return this;
  		}
  	}, {
  		key: "applyMatrix3",
  		value: function applyMatrix3(m) {

  			var x = this.x,
  			    y = this.y,
  			    z = this.z;
  			var e = m.elements;

  			this.x = e[0] * x + e[3] * y + e[6] * z;
  			this.y = e[1] * x + e[4] * y + e[7] * z;
  			this.z = e[2] * x + e[5] * y + e[8] * z;

  			return this;
  		}
  	}, {
  		key: "applyMatrix4",
  		value: function applyMatrix4(m) {

  			var x = this.x,
  			    y = this.y,
  			    z = this.z;
  			var e = m.elements;

  			this.x = e[0] * x + e[4] * y + e[8] * z + e[12];
  			this.y = e[1] * x + e[5] * y + e[9] * z + e[13];
  			this.z = e[2] * x + e[6] * y + e[10] * z + e[14];

  			return this;
  		}
  	}]);
  	return Vector3;
  }();

  var Box3 = function () {
  	function Box3() {
  		var min = arguments.length > 0 && arguments[0] !== undefined ? arguments[0] : new Vector3(Infinity, Infinity, Infinity);
  		var max = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : new Vector3(-Infinity, -Infinity, -Infinity);
  		classCallCheck(this, Box3);


  		this.min = min;

  		this.max = max;
  	}

  	createClass(Box3, [{
  		key: "set",
  		value: function set$$1(min, max) {

  			this.min.copy(min);
  			this.max.copy(max);

  			return this;
  		}
  	}, {
  		key: "copy",
  		value: function copy(b) {

  			this.min.copy(b.min);
  			this.max.copy(b.max);

  			return this;
  		}
  	}, {
  		key: "clone",
  		value: function clone() {

  			return new this.constructor().copy(this);
  		}
  	}, {
  		key: "expandByPoint",
  		value: function expandByPoint(p) {

  			this.min.min(p);
  			this.max.max(p);

  			return this;
  		}
  	}, {
  		key: "union",
  		value: function union(b) {

  			this.min.min(b.min);
  			this.max.max(b.max);

  			return this;
  		}
  	}, {
  		key: "setFromPoints",
  		value: function setFromPoints(points) {

  			var i = void 0,
  			    l = void 0;

  			for (i = 0, l = points.length; i < l; ++i) {

  				this.expandByPoint(points[i]);
  			}

  			return this;
  		}
  	}, {
  		key: "setFromCenterAndSize",
  		value: function setFromCenterAndSize(center, size) {

  			var halfSize = size.clone().multiplyScalar(0.5);

  			this.min.copy(center).sub(halfSize);
  			this.max.copy(center).add(halfSize);

  			return this;
  		}
  	}, {
  		key: "intersectsBox",
  		value: function intersectsBox(box) {

  			return !(box.max.x < this.min.x || box.min.x > this.max.x || box.max.y < this.min.y || box.min.y > this.max.y || box.max.z < this.min.z || box.min.z > this.max.z);
  		}
  	}]);
  	return Box3;
  }();

  var Matrix3 = function () {
  	function Matrix3() {
  		classCallCheck(this, Matrix3);


  		this.elements = new Float32Array([1, 0, 0, 0, 1, 0, 0, 0, 1]);
  	}

  	createClass(Matrix3, [{
  		key: "set",
  		value: function set$$1(m00, m01, m02, m10, m11, m12, m20, m21, m22) {

  			var te = this.elements;

  			te[0] = m00;te[3] = m01;te[6] = m02;
  			te[1] = m10;te[4] = m11;te[7] = m12;
  			te[2] = m20;te[5] = m21;te[8] = m22;

  			return this;
  		}
  	}, {
  		key: "identity",
  		value: function identity() {

  			this.set(1, 0, 0, 0, 1, 0, 0, 0, 1);

  			return this;
  		}
  	}, {
  		key: "copy",
  		value: function copy(m) {

  			var me = m.elements;

  			return this.set(me[0], me[3], me[6], me[1], me[4], me[7], me[2], me[5], me[8]);
  		}
  	}, {
  		key: "clone",
  		value: function clone() {

  			return new this.constructor().copy(this);
  		}
  	}]);
  	return Matrix3;
  }();

  var SymmetricMatrix3 = function () {
  	function SymmetricMatrix3() {
  		classCallCheck(this, SymmetricMatrix3);


  		this.elements = new Float32Array([1, 0, 0, 1, 0, 1]);
  	}

  	createClass(SymmetricMatrix3, [{
  		key: "set",
  		value: function set$$1(m00, m01, m02, m11, m12, m22) {

  			var e = this.elements;

  			e[0] = m00;
  			e[1] = m01;e[3] = m11;
  			e[2] = m02;e[4] = m12;e[5] = m22;

  			return this;
  		}
  	}, {
  		key: "identity",
  		value: function identity() {

  			this.set(1, 0, 0, 1, 0, 1);

  			return this;
  		}
  	}, {
  		key: "copy",
  		value: function copy(m) {

  			var me = m.elements;

  			this.set(me[0], me[1], me[2], me[3], me[4], me[5]);

  			return this;
  		}
  	}, {
  		key: "clone",
  		value: function clone() {

  			return new this.constructor().copy(this);
  		}
  	}, {
  		key: "toMatrix3",
  		value: function toMatrix3(m) {

  			var me = m.elements;

  			m.set(me[0], me[1], me[2], me[1], me[3], me[4], me[2], me[4], me[5]);
  		}
  	}, {
  		key: "add",
  		value: function add(m) {

  			var te = this.elements;
  			var me = m.elements;

  			te[0] += me[0];
  			te[1] += me[1];te[3] += me[3];
  			te[2] += me[2];te[4] += me[4];te[5] += me[5];

  			return this;
  		}
  	}, {
  		key: "norm",
  		value: function norm() {

  			var e = this.elements;

  			var m01m01 = e[1] * e[1];
  			var m02m02 = e[2] * e[2];
  			var m12m12 = e[4] * e[4];

  			return Math.sqrt(e[0] * e[0] + m01m01 + m02m02 + m01m01 + e[3] * e[3] + m12m12 + m02m02 + m12m12 + e[5] * e[5]);
  		}
  	}, {
  		key: "off",
  		value: function off() {

  			var e = this.elements;

  			return Math.sqrt(2 * (e[1] * e[1] + e[2] * e[2] + e[4] * e[4]));
  		}
  	}, {
  		key: "applyToVector3",
  		value: function applyToVector3(v) {

  			var x = v.x,
  			    y = v.y,
  			    z = v.z;
  			var e = this.elements;

  			v.x = e[0] * x + e[1] * y + e[2] * z;
  			v.y = e[1] * x + e[3] * y + e[4] * z;
  			v.z = e[2] * x + e[4] * y + e[5] * z;

  			return v;
  		}
  	}], [{
  		key: "calculateIndex",
  		value: function calculateIndex(i, j) {

  			return 3 - (3 - i) * (2 - i) / 2 + j;
  		}
  	}]);
  	return SymmetricMatrix3;
  }();

  exports.Vector2 = Vector2;
  exports.Box3 = Box3;
  exports.Matrix3 = Matrix3;
  exports.SymmetricMatrix3 = SymmetricMatrix3;
  exports.Vector3 = Vector3;

  Object.defineProperty(exports, '__esModule', { value: true });

})));
