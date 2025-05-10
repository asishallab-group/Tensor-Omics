"use strict";

import { handler as dataHandler } from "./dataHandler.js";

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

/**
 * Plots data points in the scene by instancing a base sphere mesh onto different chunks.
 * 
 * Data points are grouped into "chunks" based on spatial positions determined by a chunk diameter.
 * Chunks that fall outside a defined sight range relative to the config coordinates are disabled,
 * and an update function dynamically loads/unloads chunks as the config position (e.g. camera)
 * changes. A 3D grid is also created for visual reference.
 *
 * @param {BABYLON.Scene} scene - The BabylonJS scene in which to plot the data.
 */
function plotData(scene) {
  // Create a base sphere mesh that serves as a template for outlier data points.
  // Using mesh instances is more performant than creating completely separate meshes.
  const meshOutliers = createSphereMesh(
    scene,
    "meshOutliers",
    "outlierDataPointColor",
    "outlierDataPointDiameter"
  );

  // Object to store TransformNodes representing spatial chunks.
  // Each key in this object corresponds to a chunk's centroid in string form.
  const chunks = {};

  // Diameter of each chunk in world units.
  const chunkDiameter = 50;

  // Range (in number of chunks) around the current chunk that should be loaded.
  const chunkLoadRange = 2;

  // Define the "sight" range: the distance threshold (in world units) used to decide whether a chunk is visible.
  // The value includes an extra 0.5 chunk diameter margin.
  const sight = (chunkLoadRange + 0.5) * chunkDiameter;

  // Retrieve the initial position from the configuration.
  const posX = config.get("x");
  const posY = config.get("y");
  const posZ = config.get("z");

  // Loop through each family available in the data handler.
  // For each family, iterate through the genes for specific organs ("Liver", "Heart", "Lung")
  // and create an instance of the outlier mesh for each data point.
  for (const family of dataHandler.families) {
    dataHandler.iterGenes(family, "Liver", "Heart", "Lung").forEach(({ coordinates, ...metaData }, i) => {
      // Convert the coordinate array into a BabylonJS Vector3 and scale it.
      const position = new BABYLON.Vector3(...coordinates).scale(1000);

      // Determine the centroid of the chunk this position falls into.
      // getChunkCentroid is assumed to return an array-like coordinate (e.g. [x, y, z])
      // which is also used as a key in the `chunks` object.
      const chunk = getChunkCentroid(position, chunkDiameter);

      // If this chunk has not been created yet, create a TransformNode for it.
      if (chunks[chunk] === undefined) {
        chunks[chunk] = new BABYLON.TransformNode(`chunk_${chunk}`, scene);

        // Disable the chunk if it lies outside the sight range relative to the config position.
        // This is done by checking if the absolute difference in any axis is greater than or equal to sight.
        if (
          Math.abs(chunk[0] - posX) >= sight ||
          Math.abs(chunk[1] - posY) >= sight ||
          Math.abs(chunk[2] - posZ) >= sight
        ) {
          chunks[chunk].setEnabled(false);
        }
      }

      // Create a new instance of the sphere mesh for the data point.
      // Name it uniquely using the index and position.
      const instance = meshOutliers.dataPoint("dataPoint_" + i, position);

      // Parent the instance to the corresponding chunk's transform node.
      instance.parent = chunks[chunk];
    });
  }

  // Determine the initial chunk centroid based on the config position.
  // This represents the "active" chunk coordinates in which data is loaded.
  let chunkCentroid = getChunkCentroid({ x: posX, y: posY, z: posZ }, chunkDiameter);

  // The maximum distance (in world units) to load/unload chunks.
  const lastChunkDist = chunkLoadRange * chunkDiameter;

  // Temporary array used to compute the centroid for chunks that need updating.
  // This array is reused within the render loop for efficiency.
  const triggeredChunkCentroid = [0, 0, 0];

  /**
   * Register a callback that fires before every render.
   * This callback compares the current chunk centroid (from config)
   * with the previous one and determines which neighboring chunks need to be loaded/unloaded.
   */
  scene.registerBeforeRender(() => {
    // Get the current chunk centroid from the config position.
    const currentChunkCentroid = getChunkCentroid(
      { x: config.get("x"), y: config.get("y"), z: config.get("z") },
      chunkDiameter
    );

    // Loop through each axis (x, y, z) and detect any change in the chunk centroid.
    // If a change is detected along an axis, adjust chunks along that axis.
    chunkCentroid.forEach((axis, i) => {
      const currentAxis = currentChunkCentroid[i];
      // Only proceed if the coordinate along this axis has changed.
      if (currentAxis !== axis) {
        // Determine the direction of movement on the changed axis.
        // If currentAxis > axis then we are moving positively (compare = 1) else negatively (compare = -1).
        const compare = currentAxis > axis ? 1 : -1;

        // Loop over all possible offsets on the remaining two axes (j and k)
        // covering the range from -lastChunkDist to lastChunkDist (in steps of chunkDiameter).
        for (let aAxis = -lastChunkDist; aAxis <= lastChunkDist; aAxis += chunkDiameter) {
          // Calculate the index for the first non-changing axis.
          const j = (i + 1) % 3;
          triggeredChunkCentroid[j] = chunkCentroid[j] + aAxis;
          for (let bAxis = -lastChunkDist; bAxis <= lastChunkDist; bAxis += chunkDiameter) {
            // Calculate the index for the second non-changing axis.
            const k = (i + 2) % 3;
            triggeredChunkCentroid[k] = chunkCentroid[k] + bAxis;

            // For the changing axis, first compute the position for the chunk that should be loaded,
            // i.e. the new chunk in range along this axis.
            triggeredChunkCentroid[i] = currentAxis + compare * lastChunkDist;
            // Call loadChunk with a flag "true" to enable the chunk.
            loadChunk(chunks[triggeredChunkCentroid], true);

            // Next, compute the position for the chunk that should be disabled,
            // i.e. the chunk that is no longer within range, so on the complementary side.
            triggeredChunkCentroid[i] = axis - compare * lastChunkDist;
            // Call loadChunk with the flag "false" to disable the chunk.
            loadChunk(chunks[triggeredChunkCentroid], false);
          }
        }
      }
    });
    // Update the chunkCentroid to the current value for use in the next frame.
    chunkCentroid = currentChunkCentroid;
  });

  // Create a 3D grid in the scene for improving interpretability.
  create3DGrid(scene);
}

function getChunkCentroid({ x, y, z }, diameter) {
  function trim(a) {
    return Math.floor((a + diameter / 2) / diameter) * diameter;
  }
  return [trim(x), trim(y), trim(z)];
}

function loadChunk(chunk, state=true) {
  chunk?.setEnabled(state);
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

  setupTooltip(scene, mesh);

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
function setupTooltip(scene, mesh) {
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

  scene.getEngine().getRenderingCanvas().addEventListener("mousemove", function (evt) {
    datapointDiv.style.left = (evt.clientX + 10) + "px";
    datapointDiv.style.top = (evt.clientY + 10) + "px";
  });
}

export { plotData, createSphereMesh };