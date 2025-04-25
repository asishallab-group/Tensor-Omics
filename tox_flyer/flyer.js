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
 * - Throws an error if the canvas element is not found.
 * - The styling is handled via CSS.
 ***************************************************************/
function configureCanvas() {
  // Retrieve the canvas element by its id "view"
  let canvas = document.getElementById("view");
  if (!canvas) {
    throw new Error("Canvas element with id 'view' not found.");
  }
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
  if (!navigator.gpu) {
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
  setupCamera(scene, canvas);

  // Create a basic hemispheric light to illuminate the scene.
  let light = new BABYLON.HemisphericLight("light", new BABYLON.Vector3(0, 1, 0), scene);
  light.intensity = 0.8;
  
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
  sphereMaterial.diffuseColor = new BABYLON.Color3(0, 0, 1); // Blue color
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

  // Helper function to create an axis arrow.
  // Each arrow is created by combining a cylinder (shaft) and a cone (head).
  function createArrow(axis, color, length) {
    // Create a parent transform node to group the arrow's parts.
    let arrow = new BABYLON.TransformNode("arrow_" + axis, scene);

    // Create the shaft of the arrow as a cylinder.
    let shaft = BABYLON.MeshBuilder.CreateCylinder("shaft_" + axis, { diameter: 0.2, height: length }, scene);
    let shaftMat = new BABYLON.StandardMaterial("shaftMat_" + axis, scene);
    shaftMat.diffuseColor = color;
    shaft.material = shaftMat;
    shaft.parent = arrow;
    // Position the shaft so that its base is at the origin.
    shaft.position.y = length / 2;
    
    // Create the head of the arrow as a cone (a cylinder with diameterTop = 0).
    let head = BABYLON.MeshBuilder.CreateCylinder("head_" + axis, { diameterTop: 0, diameterBottom: 0.5, height: 1.5, tessellation: 20 }, scene);
    let headMat = new BABYLON.StandardMaterial("headMat_" + axis, scene);
    headMat.diffuseColor = color;
    head.material = headMat;
    head.parent = arrow;
    // Position the head on top of the shaft.
    head.position.y = length + .75;
    
    // Rotate the arrow so that it points in the correct axis direction.
    if (axis === "x") {
      arrow.rotation.z = -Math.PI / 2; // Point along the X-axis.
    } else if (axis === "z") {
      arrow.rotation.x = Math.PI / 2;  // Point along the Z-axis.
    }
    return arrow;
  }
  
  // Create the three coordinate axes arrows in black.
  createArrow("x", new BABYLON.Color3(0, 0, 0), Math.max(...tissueX) * 1.1);
  createArrow("y", new BABYLON.Color3(0, 0, 0), Math.max(...tissueY) * 1.1);
  createArrow("z", new BABYLON.Color3(0, 0, 0), Math.max(...tissueZ) * 1.1);
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

    // Step 4: Setup the tooltip follow behavior so tooltips stay near the pointer.
    setupTooltipFollow(canvas);

    // Step 5: Generate mock data for plotting (simulate WASM-provided data).
    let mockData = createMockData(200); // Increase the number of points for a denser plot.

    // Step 6: Plot the data points and draw the ordinate arrows.
    plotData(scene, mockData);

    // Step 7: Run the render loop to continuously update the scene.
    engine.runRenderLoop(() => {
      scene.render();
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


