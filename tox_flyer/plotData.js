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

/***************************************************************
 * Function: plotData
 * Purpose: Visualize the data points and add axis ordinates.
 * - Creates a base blue sphere mesh; for performance reasons instances are used for each data point.
 * - Registers pointer events for tooltips (on hover) and for clicking to light up a point.
 * - Also creates 3D arrows along the X, Y, and Z axes (in black) to serve as ordinates.
 ***************************************************************/
function plotData(scene) {
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

  create3DGrid(scene);
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