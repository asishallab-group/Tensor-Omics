/***************************************************************
 * Global Variables and State
 ***************************************************************/
// Global object to store the current view state (position and rotation)
// This makes the view reproducible as you can always use 'currentViewState'
// to retrieve the camera's current position and orientation.
let currentViewState = {
  position: { x: 0, y: 0, z: 0 },
  rotation: { x: 0, y: 0, z: 0 }
};

// Babylon engine, scene and camera will be stored in these global variables.
let engine;
let scene;
let camera;
// Reference to the tooltip element for displaying hovered point data.
let datapointDiv = document.getElementById("datapoint");

/***************************************************************
 * Function: configureCanvas
 * Purpose: Retrieve and configure the canvas element.
 * - Sets focus to the canvas on page load to enable keyboard controls immediately.
 * - Throws an error if the canvas element is not found.
 ***************************************************************/
function configureCanvas() {
  // Retrieve the canvas element by its id "view".
  let canvas = document.getElementById("view");
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
  engine = new BABYLON.WebGPUEngine(canvas);
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
 * Purpose: Configure a UniversalCamera for movement and rotation controls.
 * - Sets up WASD and arrow keys for navigation in 3D space.
 * - Q/E for rotating left/right, Space/Shift for upward/downward movement.
 * - Attaches camera controls to the canvas for mouse-based view rotation.
 * - Binds mouse wheel to apply zoom effect to the view.
 */
function setupCamera(scene, canvas) {
  // Create a UniversalCamera placed initially above the ground and away from the origin.
  // UniversalCamera is suited for first-person style movement and rotation in 3D space.
  camera = new BABYLON.UniversalCamera("camera", new BABYLON.Vector3(0, 5, -20), scene);

  // Customize key bindings for movement (WASD, Arrow keys, etc.).
  // Movement forward/backward is controlled by W/S (87/83) and ArrowUp/ArrowDown keys.
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

  // Attach the camera controls to the canvas to enable mouse and keyboard usage.
  camera.attachControl(canvas, true);

  // Set the movement speed and mouse sensitivity for a smooth experience.
  camera.speed = 0.5; // Controls the speed of movement for WASD and arrow keys.
  camera.angularSensibility = 2000; // Controls mouse drag sensitivity for view rotation.

  // Listen to the mouse wheel event on the canvas to simulate zooming.
  canvas.addEventListener("wheel", event => {
    // event.deltaY is positive when scrolling down (zoom out) and negative when scrolling up (zoom in).
    let delta = event.deltaY * 0.0005;
    // Adjust the camera's field of view (fov) to simulate zoom changes.
    camera.fov += delta;
    // Clamp the FOV value to keep the zoom within sensible limits.
    camera.fov = Math.min(Math.max(camera.fov, 0.1), 1.5);
  });
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
  let scene = new BABYLON.Scene(engine);

  // set white background
  scene.clearColor = BABYLON.Color3.White();

  create3DGrid(scene);
  setupCamera(scene, canvas);

  // Create a basic hemispheric light to illuminate the scene.
  let light = new BABYLON.HemisphericLight("light", new BABYLON.Vector3(0, 1, 0), scene);
  light.intensity = 0.8;
  
  showPositionOverlay(scene, camera);

  // Update the global view state every frame. This ensures reproducibility of the view.
  scene.registerBeforeRender(() => {
    currentViewState.position = {
      x: camera.position.x,
      y: camera.position.y,
      z: camera.position.z
    };
    // The camera rotation is stored in Euler angles.
    currentViewState.rotation = {
      x: camera.rotation.x,
      y: camera.rotation.y,
      z: camera.rotation.z
    };
  });

  return scene;
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
  let gridParent = new BABYLON.TransformNode("gridParent", scene);

  // Helper function to create a single line.
  function createLine(start, end, color) {
    let line = BABYLON.MeshBuilder.CreateLines("line", { points: [start, end] }, scene);
    let material = new BABYLON.StandardMaterial("lineMat", scene);
    material.emissiveColor = color; // Use emissive color to make lines bright.
    material.disableLighting = true; // Lines are unaffected by scene lighting.
    material.alpha = 0.2;
    line.material = material;
    line.parent = gridParent; // Attach the line to the parent node.
  }

  // Create grid lines parallel to each axis.
  let halfSize = size / 2;
  let color = BABYLON.Color3.Gray(); // Default color for grid lines.

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
  let tissueX = new Float32Array(numPoints);
  let tissueY = new Float32Array(numPoints);
  let tissueZ = new Float32Array(numPoints);
  
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
  let { tissueX, tissueY, tissueZ } = data;
  let numPoints = tissueX.length;
  
  // Create a base sphere that serves as a template.
  // Using instances is more performant than creating completely separate meshes.
  let baseSphere = BABYLON.MeshBuilder.CreateSphere("baseSphere", { diameter: 0.25, segments: 16 }, scene);
  
  // Create and assign a blue material for the spheres.
  let sphereMaterial = new BABYLON.StandardMaterial("sphereMat", scene);
  sphereMaterial.diffuseColor = new BABYLON.Color3.Blue(); // Blue color
  baseSphere.material = sphereMaterial;
  
  // Hide the original sphere since we will use instances.
  baseSphere.isVisible = false;
  
  // Loop through the mock data and create an instance for each point.
  for (let i = 0; i < numPoints; i++) {
    let sphereInstance = baseSphere.createInstance("sphere_" + i);
    sphereInstance.position = new BABYLON.Vector3(tissueX[i], tissueY[i], tissueZ[i]);

    // Store the coordinate data in the instance's metadata.
    sphereInstance.metadata = {
      x: tissueX[i],
      y: tissueY[i],
      z: tissueZ[i]
    };

    // Enable pointer interactions by attaching an ActionManager to each instance.
    sphereInstance.actionManager = new BABYLON.ActionManager(scene);

    // Register a hover action to display a tooltip with the data values.
    sphereInstance.actionManager.registerAction(
      new BABYLON.ExecuteCodeAction(BABYLON.ActionManager.OnPointerOverTrigger, function (evt) {
        const pickResult = evt.source;
        if (pickResult && pickResult.metadata) {
          datapointDiv.style.display = "block";
          // Format the tooltip content with two decimal places.
          datapointDiv.innerHTML = "x: " + pickResult.metadata.x.toFixed(2) + 
                                 "<br>y: " + pickResult.metadata.y.toFixed(2) + 
                                 "<br>z: " + pickResult.metadata.z.toFixed(2);
        }
      })
    );

    // Hide the tooltip when the pointer leaves the sphere.
    sphereInstance.actionManager.registerAction(
      new BABYLON.ExecuteCodeAction(BABYLON.ActionManager.OnPointerOutTrigger, function () {
        datapointDiv.style.display = "none";
      })
    );

    // Register a click action to temporarily highlight the sphere.
    sphereInstance.actionManager.registerAction(
      new BABYLON.ExecuteCodeAction(BABYLON.ActionManager.OnPickTrigger, function (evt) {
        let mesh = evt.source;
        // Toggle a highlight effect: here, setting emissiveColor to Yellow.
        if (mesh.material?.emissiveColor?.equals(BABYLON.Color3.Yellow())) {
          mesh.material.emissiveColor = BABYLON.Color3.Black();
        } else {
          mesh.material.emissiveColor = BABYLON.Color3.Yellow();
        }
      })
    );
  }
}

/***************************************************************
 * Function: setupTooltipFollow
 * Purpose: Make the tooltip div follow the mouse pointer.
 * - Listens for mousemove events on the canvas.
 * - Offsets the tooltip by a few pixels from the pointer for better visibility.
 ***************************************************************/
function setupTooltipFollow(canvas) {
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
 * @param {BABYLON.Camera} camera - The camera whose position is displayed.
 */
function showPositionOverlay(scene, camera) {
  // Create a fullscreen GUI overlay using Babylon.js's GUI library.
  let advancedTexture = BABYLON.GUI.AdvancedDynamicTexture.CreateFullscreenUI("UI");

  // Create a text block to display the position.
  let xPosition = new BABYLON.GUI.TextBlock();
  xPosition.fontSize = 24; // Font size
  xPosition.fontStyle = "bold"; // Font style
  xPosition.textHorizontalAlignment = BABYLON.GUI.Control.HORIZONTAL_ALIGNMENT_RIGHT; // Align text to the right
  xPosition.textVerticalAlignment = BABYLON.GUI.Control.VERTICAL_ALIGNMENT_BOTTOM; // Align text to the bottom
  xPosition.paddingRight = 10; // Add some padding from the right
  xPosition.paddingBottom = 300; // Add some padding from the top
  xPosition.color = "red"; // Text color

  yPosition = xPosition.clone();
  yPosition.color = "green"; // Text color
  yPosition.paddingBottom = 270;
  zPosition = yPosition.clone();
  zPosition.color = "blue"; // Text color
  zPosition.paddingBottom = 240;

  // Add the text blocks to the GUI overlay.
  advancedTexture.addControl(xPosition);
  advancedTexture.addControl(yPosition);
  advancedTexture.addControl(zPosition);

  // Update the position text dynamically as the camera moves.
  scene.registerBeforeRender(() => {
    xPosition.text = `X: ${camera.position.x.toFixed(2)}`;
    yPosition.text = `Y: ${camera.position.y.toFixed(2)}`;
    zPosition.text = `Z: ${camera.position.z.toFixed(2)}`;
  });
}

/**
 * Function: add3DCompass
 * Purpose: Adds a 3D compass fixed to the bottom-right corner of the canvas.
 * - A mini coordinate system (X, Y, Z axes) is always visible.
 * - Rotates synchronously (inverted) with the camera to show the correct orientation of the axes.
 * @param {BABYLON.Scene} mainScene - The primary scene of your application.
 * @param {BABYLON.Camera} mainCamera - The main camera whose rotation is used to synchronize the compass.
 * @param {BABYLON.Engine} engine - The Babylon.js engine used for rendering.
 * @returns {BABYLON.Scene} - The mini scene containing the interactive compass.
 */
function add3DCompass(mainScene, mainCamera, engine) {
  // Create a new scene dedicated to the compass visualization.
  const compassScene = new BABYLON.Scene(engine);

  // Prevent the compass scene from clearing the canvas to maintain visibility of the main scene.
  compassScene.autoClear = false;
  compassScene.clearColor = new BABYLON.Color4(0, 0, 0, 0); // Transparent background.

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

  const createAxis = (direction, color) => {
    // Create the axis line.
    const axis = BABYLON.MeshBuilder.CreateTube(
      `${direction}Axis`,
      { path: [BABYLON.Vector3.Zero(), direction.scale(axisSize)], radius: axisRadius, cap: BABYLON.Mesh.CAP_END },
      compassScene
    );
    const axisMaterial = new BABYLON.StandardMaterial(`${direction}AxisMat`, compassScene);
    axisMaterial.diffuseColor = color; // Match the arrowhead's color with the axis.
    axis.material = axisMaterial;
    axis.parent = cross;
  };

  // Create the X, Y, and Z axes with their respective colors.
  createAxis(new BABYLON.Vector3(-1, 0, 0), new BABYLON.Color3.Red())
  createAxis(new BABYLON.Vector3(0, 1, 0), new BABYLON.Color3.Green())
  createAxis(new BABYLON.Vector3(0, 0, -1), new BABYLON.Color3.Blue())
  BABYLON.MeshBuilder.CreateSphere("origin", { diameter: 5 * axisRadius }, compassScene);

  // Synchronize the compass with the main camera's rotation.
  mainScene.onBeforeRenderObservable.add(() => {
    let camQuat = mainCamera.rotationQuaternion; // Get the camera's rotation as a quaternion.
    if (!camQuat) {
      // If the camera does not use quaternions, convert its Euler angles to a quaternion.
      camQuat = BABYLON.Quaternion.RotationYawPitchRoll(
        mainCamera.rotation.y, // Yaw (horizontal rotation).
        -mainCamera.rotation.x, // Pitch (vertical tilt).
        mainCamera.rotation.z  // Roll (bank rotation).
      );
    }
    const inverseQuat = camQuat.clone().invert(); // Compute the inverse rotation.

    // Apply the inverse rotation to the axes (arrowheads automatically follow).
    cross.rotationQuaternion = inverseQuat; // Arrowheads inherit this rotation from their parent.

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
    await initializeEngine(canvas);

    // Step 3: Setup the scene, including the camera, lighting, and basic controls.
    scene = setupScene(engine, canvas);
    compassScene = add3DCompass(scene, camera, engine);

    // Step 4: Setup the tooltip follow behavior so tooltips stay near the pointer.
    setupTooltipFollow(canvas);

    // Step 5: Generate mock data for plotting (simulate WASM-provided data).
    let mockData = createMockData(200); // Increase the number of points for a denser plot.

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
    document.body.innerHTML = `Error during initialization: ${err}`;
  }
}

// Start the application by calling the main function.
main();
