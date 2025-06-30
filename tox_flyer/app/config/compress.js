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

function getEncoder(defaults, familyKeyTypes) {
  // add builtin keys as separators
  for (const [category, settings] of Object.entries(defaults)) {
    for (const key of Object.keys(settings)) {
      const sepType = `${key}:${category}`;
      const sep = String.fromCharCode(MAX_ENC_CHAR_CODE + SEPARATORS.size + 1);
      SEPARATORS.set(sepType, sep);
      SEPARATORS_REV.set(sep, sepType);
    }
  }
  for (const keyType in familyKeyTypes) {
    const sep = String.fromCharCode(MAX_ENC_CHAR_CODE + SEPARATORS.size + 1);
    SEPARATORS.set(keyType, sep);
    SEPARATORS_REV.set(sep, keyType);
  }


  function encode(category, key, value) {
    const encoders = {
      boolean: (encodedKey, value, defaultValue) => {
        if (value !== defaultValue) {
          return encodedKey;
        } else {
          return "";
        }
      },
      number: (encodedKey, value, defaultValue) => {
        if (value !== defaultValue) {
          if (Number.isInteger(value)) {
            return encodedKey + encodeInt(value);
          } else {
            return encodedKey + encodeFloat(value);
          }
        } else {
          return "";
        }
      },
      color: (encodedKey, value, defaultValue) => {
        if (value !== defaultValue) {
          return encodedKey + encodeColor(value);
        } else {
          return "";
        }
      },
      string: (encodedKey, value, defaultValue) => {
        console.log("Encoding strings not yet supported");
        return "";
      },
      object: (encodedKey, value, defaultValue) => {
        console.log("Encoding objects not yet supported");
        return "";
      }
    }

    let type;
    let defaultValue;
    let encodedKey;
    if (key.includes("_")) {
      encodedKey = SEPARATORS.get("familyKey_" + category);
      const [family, keyType, gene] = key.match(/^(\d+)_(\D+)(:\d+)?$/).slice(1);

      encodedKey += encodeInt(Number.parseInt(family));

      if (!SEPARATORS.has(keyType)) throw new Error(`Missing family key type: '${keyType}'`);
      encodedKey += SEPARATORS.get(keyType);

      if (gene !== undefined) {
        encodedKey += encodeInt(Number.parseInt(gene.slice(1)));
      }
      
      encodedKey += SEPARATORS.get("familyKeyStop");
      
      defaultValue = familyKeyTypes[keyType].default;
      type = familyKeyTypes[keyType].type
    } else {
      encodedKey = SEPARATORS.get(`${key}:${category}`);
      defaultValue = defaults[category][key];
    }
    let encoder;
    if (category !== "allModes") {
      encoder = encoders.color;
    } else {
      encoder = encoders[typeof value];
    }
    return encoder(encodedKey, value, defaultValue);
  }

  function decode(str, idx) {
    const decoders = {
      number: (str, idx) => {
        const decoded = decodeNumber(str, idx);
        return { value: decoded.decoded, idx: decoded.idx };
      },
      color: (str, idx) => {
        const decoded = decodeColor(str, idx);
        return { value: decoded.decoded, idx: decoded.idx };
      },
      string: (str, idx) => {
        console.log("Decoding strings not yet supported");
        return "";
      },
      object: (str, idx) => {
        console.log("Decoding objects not yet supported");
        return "";
      }
    }

    const sepType = SEPARATORS_REV.get(str.charAt(idx));
    if (sepType === undefined) throw new Error("Corrupted key encoding: Missing start");
    idx++;

    let category;
    let key;
    let type;
    let defaultValue;
    if (!sepType.startsWith("familyKey")) {
      [key, category] = sepType.split(":");
      type = typeof defaults[category][key];
      if (type === "boolean") {
        defaultValue = defaults[category][key];
      }
    } else {
      category = sepType.split("_")[1];
      if (category === undefined) throw new Error("Corrupted family key encoding: Unknown category");

      const decodedFamily = decodeInt(str, idx);
      idx = decodedFamily.idx;

      const keyType = SEPARATORS_REV.get(str.charAt(idx));
      if (keyType === undefined) throw new Error("Corrupted family key encoding: Unknown key type");
      key = `${decodedFamily.decoded}_${keyType}`;

      idx++;
      type = familyKeyTypes[keyType].type;
      if (type === "boolean") {
        defaultValue = familyKeyTypes[keyType].default;
      }
      if (str.charAt(idx) === SEPARATORS.get("familyKeyStop")) {
        idx++;
      } else {
        const decodedGene = decodeInt(str, idx);
        if (str.charAt(decodedGene.idx) !== SEPARATORS.get("familyKeyStop")) throw new Error("Corrupted family key encoding: Missing stop");

        idx = decodedGene.idx + 1;
        key += ":" + decodedGene.decoded;
      }
    }

    if (type === "boolean") {
      return { category, key, value: !defaultValue, idx };
    } else {
      let decoder;
      if (category !== "allModes") {
        decoder = decoders.color;
      } else {
        decoder = decoders[type];
      }

      return { category, key, ...decoder(str, idx) };
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

function decodeNumber(str, idx) {
  const firstPart = decodeInt(str, idx);
  idx = firstPart.idx;

  let decoded = firstPart.decoded;
  if (str.charAt(idx) === SEPARATORS.get("float")) {
    const buffer = new ArrayBuffer(8); // 64 bits
    const floatView = new Float64Array(buffer);
    const intView = new Int32Array(buffer);
    
    intView[0] = firstPart.decoded;
    idx++;

    const secondPart = decodeInt(str, idx);
    idx = secondPart.idx;
    intView[1] = secondPart.decoded;
  
    decoded = floatView[0];
  }

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

function encode(encoder, values) {
  let encoded = "";
  for (const [category, settings] of Object.entries(values)) {
    for (const [key, value] of Object.entries(settings)) {
      encoded += encoder(category, key, value);
    }
  }
  return encoded;
}

function encodeBase64(encoder, values) {
  const encoded = encode(encoder, values);
  return btoa(unescape(encodeURIComponent(encoded)));
}

function decode(decoder, encodedStr) {
  const values = {
    allModes: {},
    lightMode: {},
    darkMode: {}
  };

  let idx = 0;
  while (idx < encodedStr.length) {
    const decoded = decoder(encodedStr, idx);
    values[decoded.category][decoded.key] = decoded.value;
    idx = decoded.idx;
  }

  return values;
}

function decodeBase64(decoder, base64) {
  const encoded = decodeURIComponent(escape(atob(base64)));
  return decode(decoder, encoded);
}

function test() {
  const tests = {
    testEncodeDecode() {
      const defaults = {
        allModes: {
          changedInt: 1,
          unchangedFloat: 1.00000001,
          unchangedBool: true
        },
        lightMode: {
          changedColor: "#ABCDEF05"
        }, 
        darkMode: {
          unchangedColor: "#123456"
        }
      };
      const custom = {
        allModes: {
          changedInt: 2,
          unchangedFloat: 1.00000001,
          unchangedBool: true,
          "123_ChangedFamilySettingWithGene:12": false
        },
        lightMode: {
          changedColor: "#66666666",
          "123_UnchangedFamilySettingWithoutGene": "#FFFFFF"
        },
        darkMode: {
          unchangedColor: "#123456"
        }
      };
      const expected = {
        allModes: {
          changedInt: 2,
          "123_ChangedFamilySettingWithGene:12": false
        }, lightMode: {
          changedColor: "#66666666"
        },
        darkMode: {}
      };

      const familyKeyTypes = {
        ChangedFamilySettingWithGene: { type: "boolean", default: true },
        UnchangedFamilySettingWithoutGene: { type: "string", default: "#FFFFFF" },
      }

      const encoder = getEncoder(defaults, familyKeyTypes);
      const encoded = encodeBase64(encoder.encode, custom);
      const decoded = decodeBase64(encoder.decode, encoded);

      if (JSON.stringify(decoded) !== JSON.stringify(expected)) {
        throw new Error("Mismatched decoding");
      }
    },
    testInt() {
      const ints = [MAX_CHAR_CODE * 12, 12, 0, MAX_CHAR_CODE, MAX_CHAR_CODE - 1, 100000000, -100000000, 2882400255, Number.MAX_SAFE_INTEGER, Number.MIN_SAFE_INTEGER];
      for (const int of ints) {
        const encodedInt = encodeInt(int);
        const decodedInt = decodeNumber(encodedInt, 0);
        if (decodedInt.decoded !== int) {
          throw new Error(`Integer decoding failed, expected ${int}, got ${decodedInt.decoded}`);
        }
        if (decodedInt.idx !== encodedInt.length) {
          throw new Error(`Integer decoding for ${int} returned wrong index, expected ${encodedInt.length}, got ${decodedInt.idx}`);
        }
      }
    },
    testFloat() {
      const floats = [MAX_CHAR_CODE * 12, 12, 0, MAX_CHAR_CODE, MAX_CHAR_CODE - 1, 100000000, -100000000, 2882400255].map(f => f * Math.PI);

      for (const float of [...floats, Number.MAX_VALUE, Number.MIN_VALUE]) {
        const encodedFloat = encodeFloat(float);
        const decodedFloat = decodeNumber(encodedFloat, 0);
        if (decodedFloat.decoded !== float) {
          throw new Error(`Float decoding failed, expected ${float}, got ${decodedFloat.decoded}`);
        }
        if (decodedFloat.idx !== encodedFloat.length) {
          throw new Error(`Float decoding for ${float} returned wrong index, expected ${encodedFloat.length}, got ${decodedFloat.idx}`);
        }
      }
    },
    testColor() {
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
  }

  for (const [test, func] of Object.entries(tests)) {
    try {
      func();
      console.log(`${test}: SUCCESS`)
    } catch (err) {
      console.error(`${test}: FAILED: ${err.message}`);
    }
  }
}

test()