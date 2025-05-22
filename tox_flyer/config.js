"use strict";

function setupConfig() {
  const defaults = {
    allModes: {
      orbitMode: false,
      darkMode: window.matchMedia('(prefers-color-scheme: dark)').matches,
      outlierDataPointDiameter: 0.25,
      x: 0,
      y: 0,
      z: 0,
      rotationX: 0,
      rotationY: 0,
      orbitModeTargetDistance: 10,
      mouseSensibility: 2000,  // the higher, the slower
      movementSpeed: 0.5,
      tissueX: "Liver",
      tissueY: "Heart",
      tissueZ: "Lung",
      chunkDiameter: 50,
      chunkLoadRange: 2,
      scale: 100,
      shownFamilies: null
    },
    lightMode: {
      selectedDataPointColor: "#FFFF00FF",
      outlierDataPointColor: "#0000FFFF",

      backgroundColor: "#FFFFFFFF",

      xAxisColor: "#FF0000FF",
      yAxisColor: "#00FF00FF",
      zAxisColor: "#0000FFFF",
    },
    darkMode: {
      selectedDataPointColor: "#FFFF00FF",
      outlierDataPointColor: "#F15829FF",

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

  const triggersChunkReload = [
    "tissueX",
    "tissueY",
    "tissueZ",
    "chunkDiameter",
    "chunkLoadRange",
    "scale",
    "darkMode",
    "shownFamilies"
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
      if (defaults.allModes[key] !== undefined) values.allModes[key] = value;
      else if (config.get("darkMode")) values.darkMode[key] = value;
      else values.lightMode[key] = value;

      if (runCallback) {
        callbacks[key]?.(value);

        // when changing family related stuff (like <familyname>_Color) or other things that need to trigger a chunk reload
        if (key.includes("_") || triggersChunkReload.includes(key)) {
          const event = new CustomEvent("chunkReload", {
            detail: { setting: key }
          });
          document.dispatchEvent(event);
        }
      }
    },
    setSetterCallback(key, callback) {
      if (callbacks[key] === undefined) {
        callbacks[key] = (value) => callback(value); // wrapping the callback to avoid this-context on the private callbacks object
        callback(config.get(key));
      } else {
        throw Error(`Another callback function has been already registered for '${key}' in the past.`);
      }
    },
    asURL() {
      const currentURL = new URL(document.URL);
      return `${currentURL.origin}${currentURL.pathname}?config=${btoa(JSON.stringify(values))}`;
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
      const importingConfig = JSON.parse(atob(configArg));
      if (importingConfig.allModes) values.allModes = importingConfig.allModes;
      if (importingConfig.darkMode) values.darkMode = importingConfig.darkMode;
      if (importingConfig.lightMode) values.lightMode = importingConfig.lightMode;
    }
  } catch (err) {
    console.log("Could not import config from URL");
  }

  return config;
}

export const config = setupConfig();
