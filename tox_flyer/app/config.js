"use strict";

const DEFAULTS = {
  allModes: {
    orbitMode: false,
    darkMode: true,
    x: 0,
    y: 0,
    z: 0,
    rotationX: 0,
    rotationY: 0,
    orbitModeTargetDistance: 10,
    mouseSensibility: 2000,  // the higher, the slower, greater zero
    movementSpeed: 0.5,
    scale: 100,
    chunkDiameter: 50,
    chunkLoadRange: 2,
    shownFamilies: null,
    tissueX: "Liver",
    tissueY: "Heart",
    tissueZ: "Lung",
    defaultDiameter: 0.25
  },
  lightMode: {
    selectedDataPointColor: "#FFFF00FF",

    backgroundColor: "#FFFFFFFF",

    xAxisColor: "#FF0000FF",
    yAxisColor: "#00FF00FF",
    zAxisColor: "#0000FFFF",
  },
  darkMode: {
    selectedDataPointColor: "#FFFF00FF",

    backgroundColor: "#1B1A1FFF",

    xAxisColor: "#DE0000FF",
    yAxisColor: "#19CF00FF",
    zAxisColor: "#0092FFFF",
  }
};

export async function setupConfig() {

  const values = {
    allModes: {},
    lightMode: {},
    darkMode: {}
  };

  let familyKeyTypes;
  
  const callbacks = {};
  const validate = getValidator();

  const triggersChunkReload = [
    "tissueX",
    "tissueY",
    "tissueZ",
    "chunkDiameter",
    "chunkLoadRange",
    "scale",
    "darkMode",
    "shownFamilies",
    "defaultDiameter"
  ];

  const config = {
    get(key) {
      let value = values.allModes[key] ?? DEFAULTS.allModes[key];
      if (value === undefined) {
        if (values.allModes.darkMode ?? DEFAULTS.allModes.darkMode) {
          value = values.darkMode[key] ?? DEFAULTS.darkMode[key];
        } else {
          value = values.lightMode[key] ?? DEFAULTS.lightMode[key];
        }
      }

      return value ?? familyDefault(key, familyKeyTypes);
    },
    set(key, value, runCallback=true) {
      validate(key, value);
      if (DEFAULTS.allModes[key] !== undefined || !key.endsWith("Color")) values.allModes[key] = value;
      else if (this.get("darkMode")) values.darkMode[key] = value;
      else values.lightMode[key] = value;

      if (runCallback) {
        callbacks[key]?.(value);

        // when changing family related stuff (like <familyname>_Color) or other things that need to trigger a chunk reload
        if (key.includes("_") || triggersChunkReload.includes(key)) {
          document.dispatchEvent(new CustomEvent("chunkReload", {
            detail: { setting: key }
          }));
        }
      }
    },
    setSetterCallback(key, callback) {
      if (callbacks[key] === undefined) {
        callbacks[key] = (value) => callback(value); // wrapping the callback to avoid this-context on the private callbacks object
        callback(this.get(key));
      } else {
        throw new Error(`Another callback function has been already registered for '${key}' in the past.`);
      }
    },
    runCallbacks() {
      config.set("darkMode", config.get("darkMode"));
    },
    async asURL() {
      const currentURL = new URL(document.URL);
      const encode = await getCompressor(familyKeyTypes, values);
      const base64 = await encode(true);
      return `${currentURL.origin}${currentURL.pathname}?config=${base64}`;
    },
    async asFile(filename="tox_flyer.conf", compressed=true) {
      let content;

      if (compressed) {
        const encode = await getCompressor(familyKeyTypes, values);
        content = await encode();
      } else {
        content = JSON.stringify(values);
      }

      // Create a blob from the string
      const blob = new Blob([content], { type: "text/plain" });

      // Create a temporary URL for the blob
      const url = URL.createObjectURL(blob);

      // Create a link and trigger the download
      const a = document.createElement("a");
      a.href = url;
      a.download = filename;
      document.body.appendChild(a);
      a.click();

      // Clean up
      document.body.removeChild(a);
      URL.revokeObjectURL(url);
    }
  }

  Object.freeze(config);

  familyKeyTypes = {
    ShiftVector: { type: "boolean", default: () => false },
    Centroid: { type: "boolean", default: () => false },
    Hull: { type: "boolean", default: () => false },
    Color: { type: "string", default: (family) => dataHandler.getColor(family) },
    OutlierColor: { type: "string", default: (family) => config.get(`${family}_Color`) },
    Diameter: { type: "number", default: () => config.get("defaultDiameter") },
    OutlierDiameter: { type: "number", default: () => config.get("defaultDiameter") },
    PickedGene: { type: "boolean", default: () => false },
    PickedShiftVector: { type: "boolean", default: () => false },
    PickedCentroid: { type: "boolean", default: () => false },
  };

  config.set("darkMode", window.matchMedia('(prefers-color-scheme: dark)').matches);

  // on dark mode switch, all callbacks need to be triggered
  config.setSetterCallback("darkMode", (enable) => {
    for (const [key, callback] of Object.entries(callbacks)) {
      if (key !== "darkMode") {
        callback(config.get(key));
      }
    }
    const event = new CustomEvent("chunkReload", {
      detail: { setting: "darkMode" }
    });
    document.dispatchEvent(event);
  })

  try {
    const currentURL = new URL(document.URL);
    const configArg = currentURL.searchParams.get("config");
    if (configArg) {
      const decode = await getCompressor(familyKeyTypes);
      const importingConfig = await decode(configArg, true);

      values.allModes = importingConfig.allModes;
      values.lightMode = importingConfig.lightMode;
      values.darkMode = importingConfig.darkMode;
    }
  } catch (err) {
    console.error("Could not import config from URL");
  }

  return config;
}

function familyDefault(key, familyKeyTypes) {
  const [family, keyType, gene] = key.match(/^(\d+)_(\D+)(:\d+)?$/)?.slice(1) ?? [];
  if (keyType !== undefined) {
    return familyKeyTypes[keyType].default(family);
  } else if (key === "shownFamilies") {
    return dataHandler.families;
  }
}

async function getCompressor(familyKeyTypes, values=null) {
  const {
    getEncoder,
    encode,
    encodeBase64,
    decode,
    decodeBase64
  } = await import("./config/compress.js");

  const defaultsCopy = JSON.parse(JSON.stringify(DEFAULTS));
  defaultsCopy.allModes.tissueX = dataHandler.tissues.indexOf(DEFAULTS.allModes.tissueX);
  defaultsCopy.allModes.tissueY = dataHandler.tissues.indexOf(DEFAULTS.allModes.tissueY);
  defaultsCopy.allModes.tissueZ = dataHandler.tissues.indexOf(DEFAULTS.allModes.tissueZ);
  const encoder = getEncoder(defaultsCopy, familyKeyTypes);

  if (values) {
    const valuesCopy = JSON.parse(JSON.stringify(values));
    for (const tissue of ["tissueX", "tissueY", "tissueZ"]) {
      const val = valuesCopy.allModes[tissue];
      if (val !== undefined) {
        valuesCopy.allModes[tissue] = dataHandler.tissues.indexOf(val);
      }
    }

    function compress(asBase64=false) {
      if (asBase64) {
        return encodeBase64(encoder.encode, valuesCopy);
      } else {
        return encode(encoder.encode, valuesCopy);
      }
    }
    return compress;
  } else {
    async function decompress(encoded, isBase64=false) {
      let decodedValues;
      if (isBase64) {
        decodedValues = await decodeBase64(encoder.decode, encoded);
      } else {
        decodedValues = await decode(encoder.decode, encoded);
      }

      for (const tissue of ["tissueX", "tissueY", "tissueZ"]) {
        const val = decodedValues.allModes[tissue];
        if (val !== undefined) {
          decodedValues.allModes[tissue] = dataHandler.tissues[val];
        }
      }

      return decodedValues;
    }

    return decompress;
  }
}

function getValidator() {
  const validators = {};
  {
    const asArray = [
      [
        ["orbitMode", "darkMode", "ShiftVector", "Centroid", "Hull", "PickedGene", "PickedShiftVector", "PickedCentroid"],
        v => {
          if (typeof v !== "boolean") throw new Error(`Expecting boolean value, got: ${typeof v}`);
        }
      ],
      [
        ["x", "y", "z", "rotationX", "rotationY"],
        v => {
          if (typeof v !== "number" || !Number.isFinite(v)) throw new Error(`Expecting number, got: ${typeof v}`);
        }
      ],
      [
        ["orbitModeTargetDistance", "mouseSensibility", "movementSpeed", "scale", "defaultDiameter", "Diameter", "OutlierDiameter"],
        v => {
          if (typeof v !== "number" || v <= 0 || !Number.isFinite(v)) throw new Error(`Expecting true positive number, got: ${v} (${typeof v})`);
        }
      ],
      [
        ["chunkDiameter"],
        v => {
          if (!Number.isInteger(v) || v <= 0 || v % 2 === 1) throw new Error(`Expecting true positive even integer, got: ${v} (${typeof v})`);
        }
      ],
      [
        ["chunkLoadRange"],
        v => {
          if (!Number.isInteger(v) || v <= 0) throw new Error(`Expecting true positive integer, got: ${v} (${typeof v})`);
        }
      ],
      [
        ["shownFamilies"],
        v => {
          if (v !== null && !(v instanceof Array)) throw new Error(`Expecting either null or Array of family names, got: ${typeof v}`);
        }
      ],
      [["tissueX", "tissueY", "tissueZ"], () => {}],
      [
        ["selectedDataPointColor", "backgroundColor", "xAxisColor", "yAxisColor", "zAxisColor", "Color"],
        v => {
          if (!/^#[A-Fa-f0-9]{6}(?:[A-Fa-f0-9]{2})?$/.test(v)) throw new Error(`Expecting RGB(A) hex color code, got: ${v}`);
        }
      ],
    ]
    for (const [keys, validator] of asArray) {
      for (const key of keys) {
        validators[key] = validator;
      }
    }
  }
  function validate(key, value) {
    const [family, keyType, gene] = key.match(/^(\d+)_(\D+)(:\d+)?$/)?.slice(1) ?? [];
    let validator;
    if (keyType !== undefined) {
      validator = validators[keyType];
    } else if (key[0].toUpperCase() !== key[0]) {
      validator = validators[key];
    }
    if (validator !== undefined) {
      try {
        validator(value);
        return true;
      } catch (err) {
        throw new Error(`${key}: ${err.message}`);
        return;
      }
    } else {
      throw new Error(`Unknown key: ${key}`);
    }
  }

  return validate;
}
