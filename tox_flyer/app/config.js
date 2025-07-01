"use strict";

export function setupConfig() {
  const defaults = {
    allModes: {
      orbitMode: false,
      darkMode: window.matchMedia('(prefers-color-scheme: dark)').matches,
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

  const values = {
    allModes: {},
    lightMode: {},
    darkMode: {}
  };
  
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
      let value = values.allModes[key] ?? defaults.allModes[key];
      if (value === undefined) {
        if (values.allModes.darkMode ?? defaults.allModes.darkMode) {
          value = values.darkMode[key] ?? defaults.darkMode[key];
        } else {
          value = values.lightMode[key] ?? defaults.lightMode[key];
        }
      }

      return value;
    },
    set(key, value, runCallback=true) {
      if (validate(key, value) !== undefined) {
        if (!key.includes("_") && this.get(key) === undefined) {
          console.error(`${key}: Unknown key`);
          return;
        }
        if (defaults.allModes[key] !== undefined || !key.endsWith("Color")) values.allModes[key] = value;
        else if (this.get("darkMode")) values.darkMode[key] = value;
        else values.lightMode[key] = value;

        if (runCallback) {
          callbacks[key]?.(value);

          // when changing family related stuff (like <familyname>_Color) or other things that need to trigger a chunk reload
          if (key.includes("_") || triggersChunkReload.includes(key) || /_ShiftVector:\d+/.test(key)) {
            document.dispatchEvent(new CustomEvent("chunkReload", {
              detail: { setting: key }
            }));
          }
        }
      }
    },
    setSetterCallback(key, callback) {
      if (callbacks[key] === undefined) {
        callbacks[key] = (value) => callback(value); // wrapping the callback to avoid this-context on the private callbacks object
        callback(this.get(key));
      } else {
        throw Error(`Another callback function has been already registered for '${key}' in the past.`);
      }
    },
    async asURL() {
      const currentURL = new URL(document.URL);
      const picked = await new Promise(resolve => {
        document.dispatchEvent(new CustomEvent(
          "feedConfig",
          { detail: {meshSelectedPoints: resolve } }
        ));
      });
      const encode = await getCompressor(defaults, values, picked);
      const base64 = encode(true);
      return `${currentURL.origin}${currentURL.pathname}?config=${base64}`;
    },
    async asFile() {
      const picked = await new Promise(resolve => {
        document.dispatchEvent(new CustomEvent(
          "feedConfig",
          { detail: {meshSelectedPoints: resolve } }
        ));
      });
      const encode = await getCompressor(defaults, values, picked);
      const content = encode();
      const filename = "tox_flyer.conf";

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
      getCompressor(defaults).then(decode => {
        const importingConfig = decode(configArg, true);
        values.allModes = importingConfig.values.allModes;
        values.lightMode = importingConfig.values.lightMode;
        values.darkMode = importingConfig.values.darkMode;
        document.addEventListener("initializePicked", evt => {
          evt.detail(importingConfig.picked);
        }, { once: true })
      })
    }
  } catch (err) {
    console.error("Could not import config from URL");
  }

  return config;
}

async function getCompressor(defaults, values=null, picked=null) {
  const {
    getEncoder,
    encode,
    encodeBase64,
    decode,
    decodeBase64
  } = await import("./config/compress.js");

  const familyKeyTypes = {
    ShiftVector: { type: "boolean", default: () => false },
    Centroid: { type: "boolean", default: () => false },
    Hull: { type: "boolean", default: () => false }, 
    Color: { type: "string", default: (family) => dataHandler.getColor(family) },
    OutlierColor: { type: "string", default: (family) => dataHandler.getFamily(family) },
    Diameter: { type: "number", default: () => defaults.allModes.defaultDiameter },
    OutlierDiameter: { type: "number", default: () => defaults.allModes.defaultDiameter },
    Picked: { type: "object", default: () => [] }
  }

  const defaultsCopy = JSON.parse(JSON.stringify(defaults));
  defaultsCopy.allModes.tissueX = dataHandler.tissues.indexOf(defaults.allModes.tissueX);
  defaultsCopy.allModes.tissueY = dataHandler.tissues.indexOf(defaults.allModes.tissueY);
  defaultsCopy.allModes.tissueZ = dataHandler.tissues.indexOf(defaults.allModes.tissueZ);
  const encoder = getEncoder(defaultsCopy, familyKeyTypes);

  if (values) {
    const valuesCopy = JSON.parse(JSON.stringify(values));
    for (const tissue of ["tissueX", "tissueY", "tissueZ"]) {
      const val = valuesCopy.allModes[tissue];
      if (val !== undefined) {
        valuesCopy.allModes[tissue] = dataHandler.tissues.indexOf(val);
      }
    }

    for (const [family, genes] of Object.entries(picked)) {
      valuesCopy.allModes[`${family}_Picked`] = genes;
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
    function decompress(encoded, isBase64=false) {
      let decodedValues;
      if (isBase64) {
        decodedValues = decodeBase64(encoder.decode, encoded);
      } else {
        decodedValues = decode(encoder.decode, encoded);
      }

      const decodedPicked = {};
      for (const [key, values] of Object.entries(decodedValues.allModes)) {
        const family = key.match(/^(.*)_Picked/)?.[1];
        if (family) {
          decodedPicked[family] = decodedValues.allModes[`${family}_Picked`];
          delete decodedValues.allModes[`${family}_Picked`];
        }
      }

      for (const tissue of ["tissueX", "tissueY", "tissueZ"]) {
        const val = decodedValues.allModes[tissue];
        if (val !== undefined) {
          decodedValues.allModes[tissue] = dataHandler.tissues[val];
        }
      }

      return { values: decodedValues, picked: decodedPicked };
    }

    return decompress;
  }
}

function getValidator() {
  const validators = {};
  {
    const asArray = [
      [
        ["orbitMode", "darkMode"],
        v => {
          if (typeof v !== "boolean") throw Error(`Expecting boolean value, got: ${typeof v}`);
        }
      ],
      [
        ["x", "y", "z", "rotationX", "rotationY"],
        v => {
          if (typeof v !== "number" || !Number.isFinite(v)) throw Error(`Expecting number, got: ${typeof v}`);
        }
      ],
      [
        ["orbitModeTargetDistance", "mouseSensibility", "movementSpeed", "scale"],
        v => {
          if (typeof v !== "number" || v <= 0 || !Number.isFinite(v)) throw Error(`Expecting true positive number, got: ${v} (${typeof v})`);
        }
      ],
      [
        ["chunkDiameter"],
        v => {
          if (!Number.isInteger(v) || v <= 0 || v % 2 === 1) throw Error(`Expecting true positive even integer, got: ${v} (${typeof v})`);
        }
      ],
      [
        ["chunkLoadRange"],
        v => {
          if (!Number.isInteger(v) || v <= 0) throw Error(`Expecting true positive integer, got: ${v} (${typeof v})`);
        }
      ],
      [
        ["shownFamilies"],
        v => {
          if (v !== null && !(v instanceof Array)) throw Error(`Expecting either null or Array of family names, got: ${typeof v}`)
        }
      ],
      [["tissueX", "tissueY", "tissueZ"], () => {}]
    ]
    for (const [keys, validator] of asArray) {
      for (const key of keys) {
        validators[key] = validator;
      }
    }
  }
  function validate(key, value) {
    const validator = validators[key];
    if (validator !== undefined) {
      try {
        validator(value);
        return true;
      } catch (err) {
        console.error(`${key}: ${err.message}`);
        return;
      }
    } else {
      if (key.endsWith("Diameter")) {
        if (typeof value !== "number" || value <= 0 || !Number.isFinite(value)) {
          console.error(`${key}: Expecting true positive number, got: ${value} (${typeof value})`);
          return;
        }
        return true;
      } else if (key.endsWith("Color")) {
        if (!/^#[A-Fa-f0-9]{6}(?:[A-Fa-f0-9]{2})?$/.test(value)) {
          console.error(`${key}: Expecting RGB(A) hex color code, got: ${value}`);
          return;
        }
        return true;
      } else if (key.endsWith("_Centroid") || key.endsWith("_Hull") || /_ShiftVector:\d+/.test(key)) {
        if (typeof value !== "boolean") {
          console.error(`${key}: Expecting boolean value, got: ${typeof value}`);
          return;
        }
        return true;
      }
      console.error(`Unknown key: ${key}`);
    }
  }

  return validate;
}
