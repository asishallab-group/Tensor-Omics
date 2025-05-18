"use strict";

function setupConfig() {
  const values = {
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
      tissueZ: "Lung"
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
    },
    callbacks: {}
  }

  const config = {
    get(key) {
      return values.allModes[key] ?? (values.allModes.darkMode ? values.darkMode[key] : values.lightMode[key]);
    },
    set(key, value, runCallback=true) {
      if (values.allModes[key] !== undefined) values.allModes[key] = value;
      else if (values.allModes.darkMode) values.darkMode[key] = value;
      else values.lightMode[key] = value;

      if (runCallback) values.callbacks[key]?.(value);
    },
    setSetterCallback(key, callback) {
      if (values.callbacks[key] === undefined) {
        values.callbacks[key] = (value) => callback(value); // wrapping the callback to avoid this-context on the private callbacks object
        callback(config.get(key));
      } else {
        throw Error(`Another callback function has been already registered for '${key}' in the past.`);
      }
    },
    asURL() {
      const currentURL = new URL(document.URL);
      return `${currentURL.origin}${currentURL.pathname}?config=${btoa(JSON.stringify({
        allModes: values.allModes,
        lightMode: values.lightMode,
        darkMode: values.darkMode
      }))}`;
    }
  }

  Object.freeze(config);

  // on dark mode switch, all callbacks need to be triggered
  config.setSetterCallback("darkMode", (enable) => {
    const entries = enable ? Object.entries({...values.darkMode}) : Object.entries({...values.lightMode});
    for (const [key, value] of entries) {
      config.set(key, value);
    }
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

Object.defineProperty(window, "config", {
    value: setupConfig(),
    writable: false, // Prevents modification
    configurable: false // Prevents deletion
});
