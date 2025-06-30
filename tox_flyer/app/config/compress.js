"use strict";

const MAX_CHAR_CODE = 2 ** 16 - 1;
const MAX_ENC_CHAR_CODE = 2 ** 15 - 1;

const SEPARATORS = new Map();
[
  "sign",
  "float"
].forEach((sep, i) => SEPARATORS.set(sep, String.fromCharCode(MAX_ENC_CHAR_CODE + i + 1)))

const CHAR_CODE_AFTER_SEP = MAX_ENC_CHAR_CODE + 1 + SEPARATORS.size;


function getKeyEncoder(defaults) {
  const offsets = new Map();
  let offset = 0;
  for (const [category, settings] of Object.entries(defaults)) {
    offsets.set(category, offset);
    offset += Object.keys(settings).length;
  }

  function encode(category, key) {
    const keys = Object.keys(defaults[category]).sort();
    return String.fromCharCode(CHAR_CODE_AFTER_SEP + offsets.get(category) + keys.indexOf(key));
  }

  function decode(str, idx) {
    const charCodePlusOffset = str.charCodeAt(idx) - CHAR_CODE_AFTER_SEP;
    for (const [category, offset] of Array.from(offsets).reverse()) {
      if (charCodePlusOffset >= offset) {
        const keys = Object.keys(defaults[category]).sort();
        const keyIdx = charCodePlusOffset - offset;
        const decoded = {category, key: keys[keyIdx]};
        idx++;
        return { decoded, idx };
      }
    }
  }

  return {encode, decode};
}

function encodeInt(int) {
  let encoded = "";
  if (int < 0) {
    encoded = SEPARATORS.get("sign");
    int = -int;
  }

  const shift = 2 ** MAX_ENC_CHAR_CODE.toString(2).length;
  do {
    encoded += String.fromCharCode(int & MAX_ENC_CHAR_CODE);
    int = Number.parseInt(int / shift);
  } while (int > 0)
  return encoded;
}

function decodeInt(str, idx) {
  let sign = 1;
  if (str.charAt(idx) === SEPARATORS.get("sign")) {
    sign = -1;
    idx++;
  }
  if (str.charCodeAt(idx) > MAX_ENC_CHAR_CODE) throw Error("Corrupted int encoding");
  let decoded = 0;
  const shift = MAX_ENC_CHAR_CODE.toString(2).length;

  let i = 0;
  while (str.charCodeAt(idx) <= MAX_ENC_CHAR_CODE) {
    decoded += str.charCodeAt(idx) * (2 ** (shift * i++));
    idx++;
  }
  return { decoded: sign * decoded, idx };
}

function encodeFloat(float) {
  const buffer = new ArrayBuffer(8); // 64 bits
  const floatView = new Float64Array(buffer);
  const intView = new Int32Array(buffer);

  // Write float value
  floatView[0] = float;

  // Read raw bits as two 32-bit integers
  return encodeInt(intView[0]) + SEPARATORS.get("float") + encodeInt(intView[1]);
}

function decodeFloat(str, idx) {
  const buffer = new ArrayBuffer(8); // 64 bits
  const floatView = new Float64Array(buffer);
  const intView = new Int32Array(buffer);

  const firstPart = decodeInt(str, idx);
  idx = firstPart.idx;
  intView[0] = firstPart.decoded;

  if (str.charAt(idx) !== SEPARATORS.get("float")) throw Error("Corrupted float encoding");
  idx++;

  const secondPart = decodeInt(str, idx);
  idx = secondPart.idx;
  intView[1] = secondPart.decoded;

  const decoded = floatView[0];
  return { decoded, idx };
}

function encodeColor(hexCode) {
  const code = hexCode.slice(1).padEnd(8, "F");
  const codeInt = Number.parseInt(code, 16);
  return encodeInt(codeInt);
}

function decodeColor(str, idx) {
  let decoded = "#";
  const decodedInt = decodeInt(str, idx);
  decoded += decodedInt.decoded.toString(16);
  return { decoded, idx: decodedInt.idx };
}

function test() {
  function testGetMapping() {
    const defaults = {allModes: {a: 1}, lightMode: {b: 1}, darkMode: {b: 1}};
    const encoder = getKeyEncoder(defaults);
    for (const [category, settings] of Object.entries(defaults)) {
      for (const key of Object.keys(settings)) {
        const encoded = encoder.encode(category, key);
        const decoded = encoder.decode(encoded, 0);
        if (decoded.idx !== 1) {
          throw Error(`Mapping returned wrong index, expected 1, got ${decoded.idx}`);
        }
        if (decoded.decoded.category !== category) {
          throw Error(`Mapping returned wrong category, expected ${category}, got ${decoded.decoded.category}`);
        }
        if (decoded.decoded.key !== key) {
          throw Error(`Mapping returned wrong key, expected ${key}, got ${decoded.decoded.key}`);
        }
      }
    }
  }

  function testInt() {
    const ints = [MAX_CHAR_CODE * 12, 12, 0, MAX_CHAR_CODE, MAX_CHAR_CODE - 1, 100000000, -100000000, 2882400255, Number.MAX_SAFE_INTEGER, Number.MIN_SAFE_INTEGER];
    for (const int of ints) {
      const encodedInt = encodeInt(int);
      const decodedInt = decodeInt(encodedInt, 0);
      if (decodedInt.decoded !== int) {
        throw Error(`Integer decoding failed, expected ${int}, got ${decodedInt.decoded}`);
      }
      if (decodedInt.idx !== encodedInt.length) {
        throw Error(`Integer decoding for ${int} returned wrong index, expected ${encodedInt.length}, got ${decodedInt.idx}`);
      }
    }
  }

  function testFloat() {
    const floats = [MAX_CHAR_CODE * 12, 12, 0, MAX_CHAR_CODE, MAX_CHAR_CODE - 1, 100000000, -100000000, 2882400255].map(f => f * Math.PI);

    for (const float of [...floats, Number.MAX_VALUE, Number.MIN_VALUE]) {
      const encodedFloat = encodeFloat(float);
      const decodedFloat = decodeFloat(encodedFloat, 0);
      if (decodedFloat.decoded !== float) {
        throw Error(`Float decoding failed, expected ${float}, got ${decodedFloat.decoded}`);
      }
      if (decodedFloat.idx !== encodedFloat.length) {
        throw Error(`Float decoding for ${float} returned wrong index, expected ${encodedFloat.length}, got ${decodedFloat.idx}`);
      }
    }
  }

  function testColor() {
    const colors = ["#ABCDEF", "#ABCDEF10"];
    for (const color of colors) {
      const encodedColor = encodeColor(color);
      const decodedColor = decodeColor(encodedColor, 0);

      const expected = color.padEnd(9, "F").toLowerCase();
      if (decodedColor.decoded !== expected) {
        throw Error(`Color encoding failed, expected '${expected}', got '${decodedColor.decoded}'`);
      }
      if (decodedColor.idx !== encodedColor.length) {
        throw Error(`Color decoding returned wrong index, expected ${encodedColor.length}, got ${decodedColor.idx}`);
      }
    }
  }

  try {
    testGetMapping();
    testInt();
    testFloat();
    testColor();
    console.log("All tests passed");
  } catch (err) {
    console.error(err);
  }
}

test()