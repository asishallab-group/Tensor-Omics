"use strict";

import { dataHandler } from "./dataHandler.js";

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
      movementSpeed: 0.5
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

/***************************************************************
 * Function: configureCanvas
 * Purpose: Retrieve and configure the canvas element.
 * - Sets focus to the canvas on page load to enable keyboard controls immediately.
 * - Throws an error if the canvas element is not found.
 ***************************************************************/
function configureCanvas() {
  // Retrieve the canvas element by its id "view".
  const canvas = document.getElementById("view");
  if (!canvas) {
    throw new Error("Canvas element with id 'view' not found.");
  }

  // Set focus to the canvas element to enable keyboard controls on load.
  canvas.tabIndex = 1; // Ensure the canvas is focusable.
  canvas.focus();

  return canvas;
}


/***************************************************************
 * Function: initializeEngine
 * Purpose: Set up Babylon.js with the experimental WebGPU engine.
 * - Checks for WebGPU support.
 * - Creates and initializes the WebGPU engine.
 * - Wraps any initialization errors in a try/catch.
 ***************************************************************/
async function initializeEngine(canvas) {
  // Check if the browser supports WebGPU. navigator.gpu is defined only if WebGPU is available.
  if (!navigator.gpu || typeof BABYLON === "undefined") {
    throw new Error("WebGPU is not supported on this browser.");
  }
  // Create a new WebGPU engine.
  // Babylon.js automatically detects that we want to use WebGPU based on this engine.
  const engine = new BABYLON.WebGPUEngine(canvas);
  try {
    // Asynchronously initialize the engine. This prepares the WebGPU adapter.
    await engine.initAsync();
  } catch (err) {
    console.error("Failed to initialize WebGPU engine: ", err);
    throw err;
  }
  // Disable offline support for a faster startup (optional setting)
  engine.enableOfflineSupport = false;

  return engine;
}

/**
 * Function: setupCamera
 * Purpose: Configure a UniversalCamera and ArcRotationCamera for movement and rotation controls.
 * - Sets up WASD and arrow keys for navigation in 3D space.
 * - Q/E for rotating left/right, Space/Shift for upward/downward movement.
 * - Attaches camera controls to the canvas for mouse-based view rotation.
 * - Binds mouse wheel to apply zoom effect to the view.
 */
function setupCamera(scene, canvas) {  
  // Create a UniversalCamera placed initially above the ground and away from the origin.
  // UniversalCamera is suited for first-person style movement and rotation in 3D space.
  const camera = new BABYLON.UniversalCamera("camera", new BABYLON.Vector3(0, 0, 0), scene);
  scene.switchActiveCamera(camera);

  // Customize key bindings for movement (WASD, Arrow keys, etc.).
  // Movement forwa rd/backward is controlled by W/S (87/83) and ArrowUp/ArrowDown keys.
  camera.keysUp = [87, 38]; // W (87) and ArrowUp (38)
  camera.keysDown = [83, 40]; // S (83) and ArrowDown (40)

  // Movement left/right is controlled by A/D (65/68) and ArrowLeft/ArrowRight keys.
  camera.keysLeft = [65, 37]; // A (65) and ArrowLeft (37)
  camera.keysRight = [68, 39]; // D (68) and ArrowRight (39)

  // Add bindings for upward and downward movement.
  // Space (32) for upward movement, Shift (16) for downward movement.
  camera.keysUpward = [32]; // Space
  camera.keysDownward = [16]; // Shift

  // Add bindings for rotation controls.
  // Q (81) for rotating left, E (69) for rotating right.
  camera.keysRotateLeft = [81]; // Q
  camera.keysRotateRight = [69]; // E

  // Set the movement speed and mouse sensitivity for a smooth experience.
  config.setSetterCallback("movementSpeed", (speed) => {
    camera.speed = speed; // Controls the speed of movement for WASD and arrow keys.
  });
  config.setSetterCallback("mouseSensibility", (sensibility) => {
    camera.angularSensibility = sensibility; // Controls mouse drag sensitivity for view rotation.
  })

  // // Listen to the mouse wheel event on the canvas to simulate zooming.
  // canvas.addEventListener("wheel", event => {
  //   // event.deltaY is positive when scrolling down (zoom out) and negative when scrolling up (zoom in).
  //   const delta = event.deltaY * 0.0005;
  //   // Adjust the camera's field of view (fov) to simulate zoom changes.
  //   camera.fov += delta;
  //   // Clamp the FOV value to keep the zoom within sensible limits.
  //   camera.fov = Math.min(Math.max(camera.fov, 0.1), 1.5);
  // });

  // create an ArcRotateCamera for orbit view
  const orbitCam = new BABYLON.ArcRotateCamera("orbitCamera", null, null, 10, new BABYLON.Vector3.Zero(), scene);
  const meshSelectedPoints = createSphereMesh(scene, "meshSelectedPoints", "selectedDataPointColor");
  setupOrbitView(scene, meshSelectedPoints);

  config.setSetterCallback("rotationX", (radians) => {
    if (config.get("orbitMode")) {
      orbitCam.beta = -radians + Math.PI / 2;
    } else {
      camera.rotation.x = radians;
    }
  });
  config.setSetterCallback("rotationY", (radians) => {
    if (config.get("orbitMode")) {
      orbitCam.alpha = radians + 1.5 * Math.PI;
    } else {
      camera.rotation.y = -radians;
    }
  });

  for (const axis of "xyz") {
    config.setSetterCallback(axis, (position) => {
      if (config.get("orbitMode")) {
        const newPosition = new BABYLON.Vector3(config.get("x"), config.get("y"), config.get("z"));
        const newTarget = getOrbitTargetFromPosition(scene, newPosition, orbitCam.radius);
        orbitCam.setTarget(newTarget);
        orbitCam.position = newPosition;
      } else {
        camera.position[axis] = position;
      }
    })
  }

  config.setSetterCallback("orbitModeTargetDistance", (radius) => {
    if (config.get("orbitMode")) {
      const newTarget = getOrbitTargetFromPosition(scene, orbitCam.position, radius);
      orbitCam.setTarget(newTarget);
    }
  })

  scene.registerBeforeRender(() => {
    // will work for both cameras
    // disabling callback function to run, as it would just set the camera to its current position
    config.set("x", scene.activeCamera.position.x, false);
    config.set("y", scene.activeCamera.position.y, false);
    config.set("z", scene.activeCamera.position.z, false);

    if (config.get("orbitMode")) {
      config.set("rotationX", (-scene.activeCamera.beta + Math.PI / 2) % (2 * Math.PI), false);
      config.set("rotationY", (scene.activeCamera.alpha - 1.5 * Math.PI) % (2 * Math.PI), false);
    } else {
      config.set("rotationX", scene.activeCamera.rotation.x, false);
      config.set("rotationY", -scene.activeCamera.rotation.y, false);
    }

    config.set("orbitModeTargetDistance", orbitCam.radius, false);
  });
}

function getOrbitTargetFromPosition(scene, position, radius) {
  const forward = scene.activeCamera.getDirection(BABYLON.Vector3.Forward());
  forward.normalize();
  const newTarget = position.add(forward.scale(radius));
  return newTarget;
}

function setupOrbitView(scene, meshSelectedPoints) {
  config.setSetterCallback("orbitMode", (enable) => {
    if (!enable && scene.activeCamera.name !== "camera") {
      const camera = scene.getCameraByName("camera");
      const position = scene.activeCamera.position;
      scene.switchActiveCamera(camera);
      camera.position = position;
    }
    else if (enable && scene.activeCamera.name !== "orbitCamera") {
      const orbitCamera = scene.getCameraByName("orbitCamera");

      let target;
      let radius;

      if (meshSelectedPoints.instances.length === 0) {
        radius = 10;
        target = getOrbitTargetFromPosition(scene, scene.activeCamera.position, radius);
      } else {
        // calculate mid point of all selected points and set as target,
        // set distance to this point as radius

        // Initialize variables to calculate the sum of positions
        let sumX = 0;
        let sumY = 0;
        let sumZ = 0;

        // Loop through all instances and sum up their positions
        meshSelectedPoints.instances.forEach(instance => {
            const position = instance.position;
            sumX += position.x;
            sumY += position.y;
            sumZ += position.z;
        });

        // Calculate the average position
        const numInstances = meshSelectedPoints.instances.length;
        target = new BABYLON.Vector3(
            sumX / numInstances,
            sumY / numInstances,
            sumZ / numInstances
        );
        radius = BABYLON.Vector3.Distance(target, scene.activeCamera.position);
      }
      scene.switchActiveCamera(orbitCamera);
      orbitCamera.setTarget(target);
      orbitCamera.radius = radius;
    }
  })
  scene.getEngine().getRenderingCanvas().addEventListener("keydown", (evt) => {
    const key = evt.key.toLowerCase();
    if (key === "f") {
      config.set("orbitMode", !config.get("orbitMode"));
    }
  })
}

/***************************************************************
 * Function: setupScene
 * Purpose: Create the Babylon scene, camera, lighting and input controls.
 * - Uses a UniversalCamera for game-like WASD movement and mouse-click drag.
 * - Sets up basic keyboard events and mouse wheel for zooming.
 * - Updates the global view state (position and rotation) every frame.
 ***************************************************************/
function setupScene(engine, canvas) {
  // Create a new Babylon scene.
  const scene = new BABYLON.Scene(engine);

  // set background color
  config.setSetterCallback("backgroundColor", (hexColorCode) => {
    setBackgroundColor(scene, hexColorCode);
  });

  create3DGrid(scene);
  setupCamera(scene, canvas);

  // Create a basic hemispheric light to illuminate the scene.
  const light = new BABYLON.HemisphericLight("light", new BABYLON.Vector3(0, 1, 0), scene);
  light.intensity = 0.8;
  
  return scene;
}

function setBackgroundColor(scene, hexColorCode) {
  scene.clearColor = BABYLON.Color4.FromHexString(hexColorCode);
}

/**
 * Function: create3DGrid
 * Purpose: Creates a 3D grid with customizable size and density.
 * - Useful for visualizing spatial boundaries in a 3D plot.
 * @param {BABYLON.Scene} scene - The Babylon.js scene where the grid will be added.
 * @param {number} size - The total size of the grid (length of each axis).
 * @param {number} step - The spacing between grid lines.
 */
function create3DGrid(scene, size = 100, step = 10) {
  // Create a parent node to group all grid lines for easy management.
  const gridParent = new BABYLON.TransformNode("gridParent", scene);

  // Helper function to create a single line.
  function createLine(start, end, color) {
    const line = BABYLON.MeshBuilder.CreateLines("line", { points: [start, end] }, scene);
    const material = new BABYLON.StandardMaterial("lineMat", scene);
    material.emissiveColor = color; // Use emissive color to make lines bright.
    material.disableLighting = true; // Lines are unaffected by scene lighting.
    material.alpha = 0.2;
    line.material = material;
    line.parent = gridParent; // Attach the line to the parent node.
  }

  // Create grid lines parallel to each axis.
  const halfSize = size / 2;
  const color = BABYLON.Color3.Gray(); // Default color for grid lines.

  // Lines along the X-axis.
  for (let a = -halfSize; a <= halfSize; a += step) {
    for (let b = -halfSize; b <= halfSize; b += step) {
      // along x axis
      createLine(
        new BABYLON.Vector3(-halfSize, a, b), // Start point
        new BABYLON.Vector3(halfSize, a, b),  // End point
        color
      );
      // along y axis
      createLine(
        new BABYLON.Vector3(a, -halfSize, b), // Start point
        new BABYLON.Vector3(a, halfSize, b),  // End point
        color
      );
      // along z axis
      createLine(
        new BABYLON.Vector3(a, b, -halfSize), // Start point
        new BABYLON.Vector3(a, b, halfSize),  // End point
        color
      );
    }
  }

  return gridParent; // Return the parent node containing all grid lines.
}

/***************************************************************
 * Function: plotData
 * Purpose: Visualize the data points and add axis ordinates.
 * - Creates a base blue sphere mesh; for performance reasons instances are used for each data point.
 * - Registers pointer events for tooltips (on hover) and for clicking to light up a point.
 * - Also creates 3D arrows along the X, Y, and Z axes (in black) to serve as ordinates.
 ***************************************************************/
function plotData(scene, data) {
  // Create a base sphere that serves as a template.
  // Using instances is more performant than creating compconstely separate meshes.
  // Mesh for basic spheres
  const meshOutliers = createSphereMesh(scene, "meshOutliers", "outlierDataPointColor", "outlierDataPointDiameter");

  // Loop through the mock data and create an instance for each point.
  for (const family of dataHandler.families) {
    dataHandler.iterGenes(family, "Liver", "Heart", "Lung").forEach(
      ({coordinates, ...metaData}, i) => {
        meshOutliers.dataPoint("dataPoint_" + i, new BABYLON.Vector3(...coordinates).scale(1000));
      }
    )
  }
}

function createSphereMesh(scene, name, configColorAttribute, configDiameterAttribute) {
    const mesh = BABYLON.MeshBuilder.CreateSphere(name, { diameter: 1, segments: 16 }, scene);

    // Create and assign a blue material for the spheres.
    const sphereMaterial = new BABYLON.StandardMaterial(name + "Mat", scene);
    mesh.material = sphereMaterial;

    // set size
    if (configDiameterAttribute) {
      config.setSetterCallback(configDiameterAttribute, (diameter) => {
        setSphereSize(mesh, diameter);
        for (const instance of mesh.instances) {
          setSphereSize(instance, diameter);
        }
      })
    }

    // set color
    config.setSetterCallback(configColorAttribute, (hexColorCode) => {
      setSphereColor(mesh, hexColorCode);
    });

    // Hide the original sphere since we will use instances.
    mesh.isVisible = false;

    mesh.dataPoint = function (name, position) {
      const instance = this.createInstance(name);
      instance.position = position;
      instance.actionManager = this.actionManager;
      return instance;
    }

    // Enable pointer interactions by attaching an ActionManager to each instance.
    mesh.actionManager = new BABYLON.ActionManager(scene);


    // Register a hover action to display a tooltip with the data values.
    const datapointDiv = document.getElementById("datapoint");
    mesh.actionManager.registerAction(
      new BABYLON.ExecuteCodeAction(BABYLON.ActionManager.OnPointerOverTrigger, function (evt) {
        const dataPoint = evt.source;
        if (dataPoint) {
          datapointDiv.style.display = "block";
          document.body.style.cursor = "pointer";
          // Format the tooltip content with two decimal places.
          datapointDiv.innerHTML = "x: " + dataPoint.position.x.toFixed(2) + 
                                 "<br>y: " + dataPoint.position.y.toFixed(2) + 
                                 "<br>z: " + dataPoint.position.z.toFixed(2);
        }
      })
    );

    // Hide the tooltip when the pointer leaves the sphere.
    mesh.actionManager.registerAction(
      new BABYLON.ExecuteCodeAction(BABYLON.ActionManager.OnPointerOutTrigger, function () {
        datapointDiv.style.display = "none";
        document.body.style.cursor = "unset";
      })
    );

    // Register a click action to select a sphere.
    mesh.actionManager.registerAction(
      new BABYLON.ExecuteCodeAction(BABYLON.ActionManager.OnPickTrigger, function (evt) {
        const instance = evt.source;
        const meshSelectedPoints = scene.getMeshByName("meshSelectedPoints");
        if (instance.TOX_unselectedInstance === undefined) {
          // Create an instance of meshSelectedPoints
          const selectedInstance = meshSelectedPoints.dataPoint(instance.name + "_selected", instance.position);
          selectedInstance.scaling = instance.scaling;

          // Hide instance
          instance.setEnabled(false);
          selectedInstance.TOX_unselectedInstance = instance;
        } else {
          instance.TOX_unselectedInstance.setEnabled(true);

          // sync scaling, as it may have changed during disabled state
          instance.TOX_unselectedInstance.scaling = instance.scaling;

          instance.dispose();
        }
      })
    );

    return mesh;
}

function setSphereSize(sphere, diameter) {
  sphere.scaling = new BABYLON.Vector3(diameter, diameter, diameter);
}

function setSphereColor(sphere, color) {
  sphere.material.diffuseColor = BABYLON.Color4.FromHexString(color);
  sphere.material.alpha = sphere.material.diffuseColor.a;
}

/***************************************************************
 * Function: setupTooltipFollow
 * Purpose: Make the tooltip div follow the mouse pointer.
 * - Listens for mousemove events on the canvas.
 * - Offsets the tooltip by a few pixels from the pointer for better visibility.
 ***************************************************************/
function setupTooltipFollow(canvas) {
  const datapointDiv = document.getElementById("datapoint");
  canvas.addEventListener("mousemove", function (evt) {
    datapointDiv.style.left = (evt.clientX + 10) + "px";
    datapointDiv.style.top = (evt.clientY + 10) + "px";
  });
}

/**
 * Function: showPositionOverlay
 * Purpose: Creates a GUI overlay that displays the camera's current position.
 * - Updates dynamically as the camera moves.
 * - Position is displayed in the bottom-right corner.
 * @param {BABYLON.Scene} scene - The Babylon.js scene for GUI integration.
 */
function showPositionOverlay(scene, xAxis, yAxis, zAxis) {
  // Create a fullscreen GUI overlay using Babylon.js's GUI library.
  const advancedTexture = BABYLON.GUI.AdvancedDynamicTexture.CreateFullscreenUI("UI", true, scene);

  // Create a text block to display the position.
  const xPosition = new BABYLON.GUI.TextBlock();
  xPosition.fontSize = "3%"; // Font size
  xPosition.fontStyle = "bold"; // Font style
  xPosition.textHorizontalAlignment = BABYLON.GUI.Control.HORIZONTAL_ALIGNMENT_RIGHT; // Align text to the right
  xPosition.textVerticalAlignment = BABYLON.GUI.Control.VERTICAL_ALIGNMENT_BOTTOM; // Align text to the bottom
  xPosition.left = "-1%"; // Add some padding from the right
  xPosition.top = "-31%"; // Add some padding from the top
  const yPosition = xPosition.clone();
  yPosition.top = "-28%"; // Add some padding from the top
  const zPosition = yPosition.clone();
  zPosition.top = "-25%"; // Add some padding from the top

  function setColorCallback(attribute, textfield, axis) {
    config.setSetterCallback(attribute, (hexColorCode) => {
      setTextfieldColor(textfield, hexColorCode);
      axis.material.diffuseColor = BABYLON.Color4.FromHexString(hexColorCode);
      axis.material.alpha = axis.material.diffuseColor.a;
    })
  }
  setColorCallback("xAxisColor", xPosition, xAxis);
  setColorCallback("yAxisColor", yPosition, yAxis);
  setColorCallback("zAxisColor", zPosition, zAxis);

  // Add the text blocks to the GUI overlay.
  advancedTexture.addControl(xPosition);
  advancedTexture.addControl(yPosition);
  advancedTexture.addControl(zPosition);

  // Update the position text dynamically as the camera moves.
  scene.registerBeforeRender(() => {
    xPosition.text = `X: ${scene.activeCamera.position.x.toFixed(2)}`;
    yPosition.text = `Y: ${scene.activeCamera.position.y.toFixed(2)}`;
    zPosition.text = `Z: ${scene.activeCamera.position.z.toFixed(2)}`;
  });
}

function setTextfieldColor(textfield, color) {
  textfield.color = color;
}

/**
 * Function: add3DCompass
 * Purpose: Adds a 3D compass fixed to the bottom-right corner of the canvas.
 * - A mini coordinate system (X, Y, Z axes) is always visible.
 * - Rotates synchronously (inverted) with the camera to show the correct orientation of the axes.
 * @param {BABYLON.Scene} mainScene - The primary scene of your application.
 * @param {BABYLON.Engine} engine - The Babylon.js engine used for rendering.
 * @returns {BABYLON.Scene} - The mini scene containing the interactive compass.
 */
function add3DCompass(mainScene, engine) {
  // Create a new scene dedicated to the compass visualization.
  const compassScene = new BABYLON.Scene(engine);

  // Prevent the compass scene from clearing the canvas to maintain visibility of the main scene.
  compassScene.autoClear = false;

  // Create an orthographic ArcRotateCamera for the compass.
  const compassCamera = new BABYLON.ArcRotateCamera(
    "compassCamera",
    Math.PI / 2, Math.PI / 2, 5, BABYLON.Vector3.Zero(), compassScene
  );
  compassCamera.mode = BABYLON.Camera.ORTHOGRAPHIC_CAMERA; // Fixed scaling.
  compassCamera.orthoLeft = -1;
  compassCamera.orthoRight = 1;
  compassCamera.orthoBottom = -1;
  compassCamera.orthoTop = 1;
  compassCamera.viewport = new BABYLON.Viewport(0.85, 0, 0.15, 0.25); // Bottom-right corner.

  // Add light to the compass scene to illuminate the axes.
  new BABYLON.HemisphericLight("compassLight", new BABYLON.Vector3(0, 1, 0), compassScene);

  // Axis and arrowhead size settings.
  const axisSize = 0.8;    // Length of the axis lines.
  const axisRadius = 0.0375;
  const cross = new BABYLON.TransformNode("cross", compassScene);

  const createAxis = (direction) => {
    // Create the axis line.
    const axis = BABYLON.MeshBuilder.CreateTube(
      `${direction}Axis`,
      { path: [BABYLON.Vector3.Zero(), direction.scale(axisSize)], radius: axisRadius, cap: BABYLON.Mesh.CAP_END },
      compassScene
    );
    const axisMaterial = new BABYLON.StandardMaterial(`${direction}AxisMat`, compassScene);
    axis.material = axisMaterial;
    axis.parent = cross;
    return axis;
  };

  // Create the X, Y, and Z axes with their respective colors and show current position above
  showPositionOverlay(
    mainScene,
    createAxis(new BABYLON.Vector3(-1, 0, 0)),
    createAxis(new BABYLON.Vector3(0, 1, 0)),
    createAxis(new BABYLON.Vector3(0, 0, -1))
  )
  // add origin
  BABYLON.MeshBuilder.CreateSphere("origin", { diameter: 5 * axisRadius }, compassScene);

  // Synchronize the compass with the main camera's rotation.
  mainScene.onBeforeRenderObservable.add(() => {
    const activeCamera = mainScene.activeCamera;
    let alpha;
    let beta;

    if (activeCamera instanceof BABYLON.ArcRotateCamera) {
      // Compute quaternion from alpha and beta for ArcRotateCamera.
      alpha = activeCamera.alpha + Math.PI / 2;
      beta = -activeCamera.beta + Math.PI / 2;
    } else {
      // Fallback for UniversalCamera or other camera types
      alpha = -activeCamera.rotation.y;
      beta = activeCamera.rotation.x;
    }
    cross.rotation.x = beta;
    cross.rotation.y = alpha;
  });


  // Return the compass scene for further customization or control.
  return compassScene;
}


/***************************************************************
 * Function: main
 * Purpose: Entry point of the application.
 * - Calls each setup function in order.
 * - Wraps the initialization steps within a try/catch block for error handling.
 * - Starts the render loop and handles browser window resizes.
 ***************************************************************/
async function main() {
  try {
    // Step 1: Configure the canvas
    const canvas = configureCanvas();

    // Step 2: Initialize Babylon's experimental WebGPU engine.
    const engine = await initializeEngine(canvas);

    // Step 3: Setup the scene, including the camera, lighting, and basic controls.
    const scene = setupScene(engine, canvas);
    console.log(scene)
    const compassScene = add3DCompass(scene, engine);

    // Step 4: Setup the tooltip follow behavior so tooltips stay near the pointer.
    setupTooltipFollow(canvas);

    // Step 5: Plot the data points and draw the ordinate arrows.
    plotData(scene);

    // Step 6: Run the render loop to continuously update the scene.
    engine.runRenderLoop(() => {
      scene.render();
      compassScene.render();
    });

    // Handle browser window resize events to adjust the canvas dimension accordingly.
    window.addEventListener("resize", () => {
      engine.resize();
    });

    // Enable snapshot rendering for the engine
    // engine.snapshotRendering = true;

  } catch (err) {
    // Log any errors during initialization to the console.
    console.log(err);
    document.body.innerHTML = `Error during initialization: ${err}`;
  }
}

// Start the application by calling the main function.
main();
