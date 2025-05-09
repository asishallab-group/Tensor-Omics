"use strict";

import { plotData, createSphereMesh } from "./plotData.js";

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

/***************************************************************
 * Function: configureCanvas
 * Purpose: Retrieve and configure the canvas element.
 * - Sets focus to the canvas on page load to enable keyboard controls immediately.
 * - Throws an error if the canvas element is not found.
 ***************************************************************/
function configureCanvas(id) {
  // Retrieve the canvas element by its id "view".
  const canvas = document.getElementById(id);
  if (!canvas) {
    throw new Error(`Canvas element with id ${id} not found.`);
  }

  // Set focus to the canvas element to enable keyboard controls on load.
  canvas.tabIndex = 1; // Ensure the canvas is focusable.
  canvas.focus();

  return canvas;
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
    const canvas = configureCanvas("view");
    const engine = await initializeEngine(canvas);
    const scene = setupScene(engine, canvas);
    console.log(scene)
    const compassScene = add3DCompass(scene, engine);

    plotData(scene);

    // Run the render loop to continuously update the scene.
    engine.runRenderLoop(() => {
      scene.render();
      compassScene.render();
    });

    // Handle browser window resize events to adjust the canvas dimension accordingly.
    window.addEventListener("resize", () => {
      engine.resize();
    });

  } catch (err) {
    // Log any errors during initialization to the console.
    console.log(err);
    document.body.innerHTML = `Error during initialization: ${err}`;
  }
}

// Start the application by calling the main function.
main();
