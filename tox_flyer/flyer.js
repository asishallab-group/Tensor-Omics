"use strict";
class Config {
  #values = {
    selectedDataPointColor: "#FFFF00FF",
    outlierDataPointColor: "#0000FFFF",
    outlierDataPointDiameter: 0.25,

    backgroundColor: "#FFFFFFFF",

    xAxisColor: "#FF0000FF",
    yAxisColor: "#00FF00FF",
    zAxisColor: "#0000FFFF",

    x: 0,
    y: 0,
    z: 0
  }
  #callbacks = {}

  constructor() {
    for (const [key, value] of Object.entries(this.#values)) {
      Object.defineProperty(this, key, {
        get() {
          return this.#values[key];
        },
        set(value) {
          this.#values[key] = value;
          this.#callbacks[key]?.(value);
        }
      })
      this[key] = value;
    }
  }

  setSetterCallback(key, callback) {
    if (this.#callbacks[key] === undefined) {
      this.#callbacks[key] = (value) => callback(value); // wrapping the callback to avoid this-context on the private callbacks object
      callback(this[key]);
    } else {
      throw Error(`Another callback function has been already registered for '${key}' in the past.`);
    }
  }
}
const config = new Config();

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
  // create an ArcRotateCamera for orbit view
  const orbitCam = new BABYLON.ArcRotateCamera("orbitCamera", -Math.PI / 2, Math.PI / 2, 10, new BABYLON.Vector3(0, 0, 0), scene);
  const meshSelectedPoints = createSphereMesh(scene, "meshSelectedPoints", "selectedDataPointColor");
  setupOrbitView(scene, meshSelectedPoints);

  // Create a UniversalCamera placed initially above the ground and away from the origin.
  // UniversalCamera is suited for first-person style movement and rotation in 3D space.
  // camera = new BABYLON.UniversalCamera("camera", new BABYLON.Vector3(0, 5, -20), scene);
  const camera = new BABYLON.UniversalCamera("camera", new BABYLON.Vector3(0, 0, 0), scene);

  // set position
  for (const axis of "xyz") {
    config.setSetterCallback(axis, (value) => {
      orbitCam.position[axis] = value;
      camera.position[axis] = value;
    })
  }

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
  camera.speed = 0.5; // Controls the speed of movement for WASD and arrow keys.
  camera.angularSensibility = 2000; // Controls mouse drag sensitivity for view rotation.

  // Listen to the mouse wheel event on the canvas to simulate zooming.
  canvas.addEventListener("wheel", event => {
    // event.deltaY is positive when scrolling down (zoom out) and negative when scrolling up (zoom in).
    const delta = event.deltaY * 0.0005;
    // Adjust the camera's field of view (fov) to simulate zoom changes.
    camera.fov += delta;
    // Clamp the FOV value to keep the zoom within sensible limits.
    camera.fov = Math.min(Math.max(camera.fov, 0.1), 1.5);
  });

  // Add function to easily switch cameras
  scene.TOX_switchCamera = (target, radius) => {
    if (scene.activeCamera.name === "camera") {
      scene.setActiveCameraByName("orbitCamera");
      if (target !== undefined) {
        scene.activeCamera.setTarget(target);
      }
      if (radius !== undefined) {
        scene.activeCamera.radius = radius;
      }
    } else {
      scene.setActiveCameraByName("camera");
      scene.activeCamera.positionQ = target;
    }
    // Attach the camera controls to the canvas to enable mouse and keyboard usage.
    scene.activeCamera.attachControl(canvas);
  }

  scene.TOX_switchCamera();
}

function setupOrbitView(scene, meshSelectedPoints) {
  scene.getEngine().getRenderingCanvas().addEventListener("keydown", (evt) => {
    const key = evt.key.toLowerCase();
    if (key === "f") {
      const camera = scene.activeCamera;
      if (camera.name === "orbitCamera") {
        scene.TOX_switchCamera();
      }
      else if (meshSelectedPoints.instances.length === 0) {
        // Calculate the forward direction vector of the camera
        const forward = camera.getDirection(BABYLON.Vector3.Forward());

        // Scale the forward vector by the orbitCamera's radius
        const orbitCam = scene.getCameraByName("orbitCamera");
        const radius = 10;
        const orbitTarget = camera.position.add(forward.scale(radius));

        scene.TOX_switchCamera(orbitTarget, radius);
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
        const middlePoint = new BABYLON.Vector3(
            sumX / numInstances,
            sumY / numInstances,
            sumZ / numInstances
        );
        const radius = BABYLON.Vector3.Distance(middlePoint, scene.activeCamera.position)
        scene.TOX_switchCamera(middlePoint, radius);
      }
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
  config.setSetterCallback("backgroundColor", (value) => {
    setBackgroundColor(scene, value);
  });

  create3DGrid(scene);
  setupCamera(scene, canvas);

  // Create a basic hemispheric light to illuminate the scene.
  const light = new BABYLON.HemisphericLight("light", new BABYLON.Vector3(0, 1, 0), scene);
  light.intensity = 0.8;
  
  return scene;
}

function setBackgroundColor(scene, color) {
  scene.clearColor = BABYLON.Color4.FromHexString(color);
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
 * Function: createMockData
 * Purpose: Generate three arrays (tissueX, tissueY, tissueZ) with random numbers.
 * - These arrays simulate the data that will eventually be provided via WASM.
 * - All arrays are of the same length to represent 3D points.
 * @param {number} numPoints - The number of points to generate.
 ***************************************************************/
function createMockData(numPoints = 100) {
  const tissueX = new Float32Array(numPoints);
  const tissueY = new Float32Array(numPoints);
  const tissueZ = new Float32Array(numPoints);
  
  // Fill each array with random values between -25 and 25.
  for (let i = 0; i < numPoints; i++) {
    tissueX[i] = Math.random() * 50 - 25;
    tissueY[i] = Math.random() * 50 - 25;
    tissueZ[i] = Math.random() * 50 - 25;
  }
  return { tissueX, tissueY, tissueZ };
}

/***************************************************************
 * Function: plotData
 * Purpose: Visualize the data points and add axis ordinates.
 * - Creates a base blue sphere mesh; for performance reasons instances are used for each data point.
 * - Registers pointer events for tooltips (on hover) and for clicking to light up a point.
 * - Also creates 3D arrows along the X, Y, and Z axes (in black) to serve as ordinates.
 ***************************************************************/
function plotData(scene, data) {
  const { tissueX, tissueY, tissueZ } = data;
  const numPoints = tissueX.length;
  
  // Create a base sphere that serves as a template.
  // Using instances is more performant than creating compconstely separate meshes.
  // Mesh for basic spheres
  const meshOutliers = createSphereMesh(scene, "meshOutliers", "outlierDataPointColor", "outlierDataPointDiameter");

  // Loop through the mock data and create an instance for each point.
  for (let i = 0; i < numPoints; i++) {
    meshOutliers.dataPoint("dataPoint_" + i, new BABYLON.Vector3(tissueX[i], tissueY[i], tissueZ[i]));
  }
}

function createSphereMesh(scene, name, configColorAttribute, configDiameterAttribute) {
    const mesh = BABYLON.MeshBuilder.CreateSphere(name, { diameter: 1, segments: 16 }, scene);

    // Create and assign a blue material for the spheres.
    const sphereMaterial = new BABYLON.StandardMaterial(name + "Mat", scene);
    mesh.material = sphereMaterial;

    // set size
    if (configDiameterAttribute) {
      config.setSetterCallback(configDiameterAttribute, (value) => {
        setSphereSize(mesh, value);
        for (const instance of mesh.instances) {
          setSphereSize(instance, value);
        }
      })
    }

    // set color
    config.setSetterCallback(configColorAttribute, (value) => {
      setSphereColor(mesh, value);
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
          instance.dispose();
        }
      })
    );

    return mesh;
}

function setSphereSize(sphere, diameter) {
  sphere.scaling.x = diameter;
  sphere.scaling.y = diameter;
  sphere.scaling.z = diameter;
}

function setSphereColor(sphere, color) {
  sphere.material.diffuseColor = BABYLON.Color4.FromHexString(color);
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
    config.setSetterCallback(attribute, (value) => {
      setTextfieldColor(textfield, value);
      axis.material.diffuseColor = BABYLON.Color4.FromHexString(value);
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
    let camQuat;

    if (activeCamera instanceof BABYLON.ArcRotateCamera) {
      // Compute quaternion from alpha and beta for ArcRotateCamera.
      const alpha = -activeCamera.alpha - Math.PI / 2;
      const beta = activeCamera.beta - Math.PI / 2;

      camQuat = BABYLON.Quaternion.RotationYawPitchRoll(alpha, beta, 0);
    } else {
      // Fallback for UniversalCamera or other camera types
      camQuat = BABYLON.Quaternion.RotationYawPitchRoll(
          activeCamera.rotation.y, // Yaw
          -activeCamera.rotation.x, // Pitch
          activeCamera.rotation.z  // Roll
      );
    }
    cross.rotationQuaternion = camQuat.invert();
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
    const compassScene = add3DCompass(scene, engine);

    // Step 4: Setup the tooltip follow behavior so tooltips stay near the pointer.
    setupTooltipFollow(canvas);

    // Step 5: Generate mock data for plotting (simulate WASM-provided data).
    const mockData = createMockData(200); // Increase the number of points for a denser plot.

    // Step 6: Plot the data points and draw the ordinate arrows.
    plotData(scene, mockData);

    // Step 7: Run the render loop to continuously update the scene.
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
