"use strict";

const MAX_CHAR_CODE = 2 ** 16 - 1;
const MAX_ENC_CHAR_CODE = 2 ** 15 - 1;

const SEPARATORS = new Map();
[
  "sign",
  "float",
  "familyKey_allModes",
  "familyKey_lightMode",
  "familyKey_darkMode",
  "familyKeyStop",
  "picked",
  "null",
  "array",
].forEach((sep, i) => SEPARATORS.set(sep, String.fromCharCode(MAX_ENC_CHAR_CODE + i + 1)));

const SEPARATORS_REV = new Map();
SEPARATORS.entries().forEach(([k, v]) => SEPARATORS_REV.set(v, k));

function getKeyEncoder(defaults, familyKeyTypes) {
  // add builtin keys as separators
  for (const [category, settings] of Object.entries(defaults)) {
    for (const key of Object.keys(settings)) {
      const sepType = `${key}:${category}`;
      const sep = String.fromCharCode(MAX_ENC_CHAR_CODE + SEPARATORS.size + 1);
      SEPARATORS.set(sepType, sep);
      SEPARATORS_REV.set(sep, sepType);
    }
  }
  const CHAR_CODE_AFTER_SEP = MAX_ENC_CHAR_CODE + 1 + SEPARATORS.size;

  function encode(category, key) {
    if (key.includes("_")) {
      let encoded = SEPARATORS.get("familyKey_" + category);
      const [family, keyType, gene] = key.match(/^(\d+)_(\D+)(:\d+)?$/).slice(1);

      encoded += encodeInt(Number.parseInt(family));

      const keyTypeIdx = familyKeyTypes.indexOf(keyType);
      if (keyTypeIdx === -1) throw new Error(`Missing family key type: '${keyType}'`);
      encoded += String.fromCharCode(CHAR_CODE_AFTER_SEP + keyTypeIdx);

      if (gene !== undefined) {
        encoded += encodeInt(Number.parseInt(gene.slice(1)));
      }
      
      encoded += SEPARATORS.get("familyKeyStop");
      return encoded;
    } else {
      return SEPARATORS.get(`${key}:${category}`);
    }
  }

  function decode(str, idx) {
    const sepType = SEPARATORS_REV.get(str.charAt(idx));
    if (sepType === undefined) throw new Error("Corrupted key encoding: Missing start");

    if (!sepType.startsWith("familyKey")) {
      const [key, category] = sepType.split(":");
      return { decoded: { category, key }, idx: idx + 1 };
    }

    // if this is reached, a family key was encoded
    const category = sepType.split("_")[1];
    if (category === undefined) throw new Error("Corrupted family key encoding: Unknown category");

    idx++;
    const decodedFamily = decodeInt(str, idx);
    idx = decodedFamily.idx;

    const keyTypeIdx = str.charCodeAt(idx) - CHAR_CODE_AFTER_SEP;
    const keyType = familyKeyTypes[keyTypeIdx];
    if (keyType === undefined) throw new Error("Corrupted family key encoding: Unknown key type");
    const key = `${decodedFamily.decoded}_${familyKeyTypes[keyTypeIdx]}`;

    idx++;
    if (str.charAt(idx) === SEPARATORS.get("familyKeyStop")) {
      return { decoded: { category, key }, idx: idx + 1 };
    } else {
      const decodedGene = decodeInt(str, idx);
      if (str.charAt(decodedGene.idx) !== SEPARATORS.get("familyKeyStop")) throw new Error("Corrupted family key encoding: Missing stop");

      return { decoded: { category, key: `${key}:${decodedGene.decoded}` }, idx: decodedGene.idx + 1 };
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
  if (str.charCodeAt(idx) > MAX_ENC_CHAR_CODE) throw new Error("Corrupted int encoding: Unexpected separator");
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

  if (str.charAt(idx) !== SEPARATORS.get("float")) throw new Error("Corrupted float encoding: Missing separator");
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
    const defaults = {allModes: {a: 1, "123_ShiftVector:12": true}, lightMode: {b: 1, "123_OutlierColor": "#FFFFFF"}, darkMode: {b: 1}};
    const encoder = getKeyEncoder(defaults, ["ShiftVector", "OutlierColor"]);
    for (const [category, settings] of Object.entries(defaults)) {
      for (const key of Object.keys(settings)) {
        const encoded = encoder.encode(category, key);
        const decoded = encoder.decode(encoded, 0);
        if (decoded.decoded.key !== key) {
          throw new Error(`Mapping returned wrong key, expected ${key}, got ${decoded.decoded.key}`);
        }
        if (decoded.decoded.category !== category) {
          throw new Error(`Mapping returned wrong category, expected ${category}, got ${decoded.decoded.category}`);
        }
        if (decoded.idx !== encoded.length) {
          throw new Error(`Mapping returned wrong index, expected ${encoded.length}, got ${decoded.idx}`);
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
        throw new Error(`Integer decoding failed, expected ${int}, got ${decodedInt.decoded}`);
      }
      if (decodedInt.idx !== encodedInt.length) {
        throw new Error(`Integer decoding for ${int} returned wrong index, expected ${encodedInt.length}, got ${decodedInt.idx}`);
      }
    }
  }

  function testFloat() {
    const floats = [MAX_CHAR_CODE * 12, 12, 0, MAX_CHAR_CODE, MAX_CHAR_CODE - 1, 100000000, -100000000, 2882400255].map(f => f * Math.PI);

    for (const float of [...floats, Number.MAX_VALUE, Number.MIN_VALUE]) {
      const encodedFloat = encodeFloat(float);
      const decodedFloat = decodeFloat(encodedFloat, 0);
      if (decodedFloat.decoded !== float) {
        throw new Error(`Float decoding failed, expected ${float}, got ${decodedFloat.decoded}`);
      }
      if (decodedFloat.idx !== encodedFloat.length) {
        throw new Error(`Float decoding for ${float} returned wrong index, expected ${encodedFloat.length}, got ${decodedFloat.idx}`);
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
        throw new Error(`Color encoding failed, expected '${expected}', got '${decodedColor.decoded}'`);
      }
      if (decodedColor.idx !== encodedColor.length) {
        throw new Error(`Color decoding returned wrong index, expected ${encodedColor.length}, got ${decodedColor.idx}`);
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