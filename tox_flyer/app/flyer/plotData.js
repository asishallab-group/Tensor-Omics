"use strict";

import {
  Mesh,
  Color,
  Vector,
  TransformNode,
  calcVectorDistance,
  Material
} from "./babylon.js";

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
  // get the initial position
  const posX = scene.activeCamera.position.x;
  const posY = scene.activeCamera.position.y;
  const posZ = scene.activeCamera.position.z;

  // Diameter of each chunk in world units.
  let chunkDiameter = config.get("chunkDiameter");

  // Range (in number of chunks) around the current chunk that should be loaded.
  let chunkLoadRange = config.get("chunkLoadRange");

  // The maximum distance (in world units) to load/unload chunks.
  let lastChunkDist = chunkLoadRange * chunkDiameter;

  // Each key in this object corresponds to a chunk's centroid in string form.
  let [ chunks, activeChunks ] = calculateChunks(posX, posY, posZ, chunkDiameter, chunkLoadRange);

  for (const chunk of activeChunks) {
    loadChunk(scene, chunks, chunk);
  }

  setupSelectionMesh(scene);

  document.addEventListener("chunkReload", (evt) => {
    [chunks, activeChunks, chunkDiameter, chunkLoadRange] = reloadChunks(scene, chunks, activeChunks);
    lastChunkDist = chunkLoadRange * chunkDiameter;
    scene.getMeshByName("meshSelectedPoints").TOX_update();
  })

  // Determine the initial chunk centroid based on the config position.
  // This represents the "active" chunk coordinates in which data is loaded.
  let chunkCentroid = getChunkCentroid([posX, posY, posZ], chunkDiameter);

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
      [ scene.activeCamera.position.x, scene.activeCamera.position.y, scene.activeCamera.position.z ],
      chunkDiameter
    );

    // Loop through each axis (x, y, z) and detect any change in the chunk centroid.
    // If a change is detected along an axis, adjust chunks along that axis.
    chunkCentroid.forEach((axis, i) => {
      const currentAxis = currentChunkCentroid[i];
      // Only proceed if the coordinate along this axis has changed.
      if (currentAxis !== axis) {

        // Determine the direction of movement on the changed axis.
        // If currentAxis > axis then we are moving positively (direction = 1) else negatively (direction = -1).
        const direction = currentAxis > axis ? 1 : -1;

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
            triggeredChunkCentroid[i] = currentAxis + direction * lastChunkDist;
            // Call loadChunk with a flag "true" to enable the chunk.
            loadChunk(scene, chunks, triggeredChunkCentroid.toString(), true);
            activeChunks.push(triggeredChunkCentroid.toString());

            // Next, compute the position for the chunk that should be disabled,
            // i.e. the chunk that is no longer within range, so on the complementary side.
            triggeredChunkCentroid[i] = axis - direction * lastChunkDist;
            // Call loadChunk with the flag "false" to disable the chunk.
            loadChunk(scene, chunks, triggeredChunkCentroid.toString(), false);
            const chunkIndex = activeChunks.indexOf(triggeredChunkCentroid.toString());
            if (chunkIndex !== -1) {
              activeChunks.splice(chunkIndex, 1);
            }
          }
        }
      }
      if (activeChunks.length !== (chunkLoadRange*2 + 1)**3) {
        console.error("Less active chunks than expected: ", activeChunks.length)
      }
    });
    // Update the chunkCentroid to the current value for use in the next frame.
    chunkCentroid = currentChunkCentroid;
  });

  let picked = null;
  scene.onPointerObservable.add((evt) => {
    switch (evt.type) {
      case BABYLON.PointerEventTypes.POINTERTAP:
        picked = pickFromMeshes(scene, chunks, activeChunks);
        if (picked !== null) {
          const selected = scene.getMeshByName("meshSelectedPoints");
          if (!selected.TOX_unpick(picked.family, picked.geneIndex)) {
            selected.TOX_pick(picked.family, picked.geneIndex);
          }
        }
        break;
      case BABYLON.PointerEventTypes.POINTERDOUBLETAP:
        if (picked !== null) {
          const selected = scene.getMeshByName("meshSelectedPoints");
          const wasSelected = selected.TOX_unpick(picked.family, picked.geneIndex);
          if (!wasSelected) {
            selected.TOX_unpick(picked.family);
          } else {
            selected.TOX_pick(picked.family);
          }
        }
        break;
    }
  });
}

function create3DGridFromChunk(scene, chunk) {
  const [chunkX, chunkY, chunkZ] = chunk.split(",").map(Number);
  const chunkRadius = config.get("chunkDiameter") / 2;
  const gridStep = 10;
  const lines = [];

  const left = chunkX - chunkRadius;
  const right = chunkX + chunkRadius;
  const bottom = chunkY - chunkRadius;
  const top = chunkY + chunkRadius;
  const back = chunkZ - chunkRadius;
  const front = chunkZ + chunkRadius;
  for (let z = back + Math.abs(back % gridStep); z <= front; z+=gridStep) {
    for (let y = bottom + Math.abs(bottom % gridStep); y <= top; y+=gridStep) {
      lines.push([Vector(left, y, z), Vector(right, y, z)]);
    }
  }
  for (let z = back + Math.abs(back % gridStep); z <= front; z+=gridStep) {
    for (let x = left + Math.abs(left % gridStep); x <= right; x+=gridStep) {
      lines.push([Vector(x, bottom, z), Vector(x, top, z)]);
    }
  }
  for (let y = bottom + Math.abs(bottom % gridStep); y <= top; y+=gridStep) {
    for (let x = left + Math.abs(left % gridStep); x <= right; x+=gridStep) {
      lines.push([Vector(x, y, back), Vector(x, y, front)]);
    }
  }

  const grid = BABYLON.MeshBuilder.CreateLineSystem(null, { lines }, scene);
  grid.color = Color(.5, .5, .5);
  return grid;
}

function setupSelectionMesh(scene) {
  const highlightLayer = new BABYLON.HighlightLayer("highlight", scene);
  const meshSelectedPoints = Mesh.Sphere(scene, "meshSelectedPoints");

  config.setSetterCallback("selectedDataPointColor", hexColorCode => {
    highlightLayer.removeMesh(meshSelectedPoints);
    highlightLayer.addMesh(meshSelectedPoints, Color(hexColorCode));
  })
  highlightLayer.setEffectIntensity(meshSelectedPoints, 0.7);

  // disable as long as spheres are picked
  highlightLayer.isEnabled = false;

  meshSelectedPoints.material = Material(scene, null, Color(0, 0, 0, 0));
  Mesh.setSize(meshSelectedPoints, 0); // hide initial instance
  meshSelectedPoints.TOX_pick = function (family, geneIndex) {
    const scale = config.get("scale");
    const inlierDiameter = config.get(`${family}_Diameter`) ?? config.get("defaultDiameter");
    const outlierDiameter = config.get(`${family}_OutlierDiameter`) ?? config.get("defaultDiameter");
    const tissues = [config.get("tissueX"), config.get("tissueY"), config.get("tissueZ")];

    function pickOne(geneIndex) {
      const { is_outlier, coordinates } = dataHandler.getGeneData(family, geneIndex, tissues, ["is_outlier"]);
      const instance = this.createInstance();
      instance.position = Vector(...coordinates.map(v => v * scale));
      instance.TOX_family = family;
      instance.TOX_geneIndex = geneIndex;
      const diameter = is_outlier ? outlierDiameter : inlierDiameter;
      Mesh.setSize(instance, diameter + 0.001); // slightly larger so the highlightLayer can truly distinguish it from the actual sphere
    }
    this.TOX_unpick(family, geneIndex);

    if (geneIndex !== undefined) {
      pickOne.call(this, geneIndex);
    } else {
      const geneCount = dataHandler.getGeneCount(family);
      for (let geneIndex = 0; geneIndex < geneCount; geneIndex++) {
        pickOne.call(this, geneIndex);
      }
    }
    highlightLayer.isEnabled = true;
  }
  meshSelectedPoints.TOX_unpick = function (family, geneIndex) {
    const instances = this.instances.filter(i => i.TOX_family === family && ((geneIndex === undefined) || (i.TOX_geneIndex === geneIndex)))
    for (const instance of instances) {
      if (this.instances.length === 1) {
        highlightLayer.isEnabled = false;
      }
      instance.dispose();
    }
    return instances.length;
  }
  meshSelectedPoints.TOX_update = function () {
    for (const instance of [...this.instances]) {
      this.TOX_pick(instance.TOX_family, instance.TOX_geneIndex);
    }
  }

  // initially fetch picked instances from config and set them up
  new Promise(resolve => {
    document.dispatchEvent(new CustomEvent("initializePicked", { detail: resolve }));
  }).then(picked => {
    try {
      if (typeof picked === "object") {
        for (const [family, genes] of Object.entries(picked)) {
          for (const geneIndex of genes) {
            meshSelectedPoints.TOX_pick(family, geneIndex);
          }
        }
      }
    } catch {
      console.error("Could not restore picked elements");;
    }
  })

  // send picked instances to config
  document.addEventListener("feedConfig", (evt) => {
    const picked = {};
    for (const instance of meshSelectedPoints.instances) {
      picked[instance.TOX_family] ??= [];
      picked[instance.TOX_family].push(instance.TOX_geneIndex);
    }
    evt.detail.meshSelectedPoints(picked);
  })
}

function pickFromMeshes(scene, chunks, activeChunks) {
  const pickRay = scene.createPickingRay(
    scene.pointerX, scene.pointerY,
    BABYLON.Matrix.Identity(),    // you can pass other transforms if you want
    scene.activeCamera
  );


  let closestDist = Infinity;
  let picked = null;
  for (const chunk of activeChunks) {
    const meshes = chunks[chunk]?.slice(-2) ?? [];
    for (const is_outlier in meshes) {
      const mesh = meshes[is_outlier];
      if (mesh) {
        const sphereMatrices = mesh.thinInstanceGetWorldMatrices(); 

        for (const i in sphereMatrices) {
          const matrix = sphereMatrices[i].m;
          // extract translation
          const tx = matrix[12];
          const ty = matrix[13];
          const tz = matrix[14];
          const spherePosition = Vector(tx, ty, tz);

          const intersects = pickRay.intersectsSphere(
            { center: spherePosition, radius: matrix[0] / 2 }
          );
          if (intersects) {
            const distance = calcVectorDistance(scene.activeCamera.position, spherePosition);
            if (distance < closestDist) {
              closestDist = distance;
              picked = {
                position: spherePosition,
                diameter: matrix[0],
                index: i,
                is_outlier,
                chunk
              }
            }
          }
        }
      }
    }
  }

  if (picked) {
    const genes = chunks[picked.chunk][0];
    let pickedIndex = picked.index;
    for (const [family, members] of genes.entries()) {
      pickedIndex -= members[picked.is_outlier].length;
      if (pickedIndex < 0) {
        picked.family = family;
        picked.geneIndex = members[picked.is_outlier][pickedIndex + members[picked.is_outlier].length];
        break;
      }
    }
  }
  return picked;
}

function reloadChunks(scene, chunks, activeChunks) {
  clearChunks(scene, chunks, activeChunks);
  const chunkDiameter = config.get("chunkDiameter");
  const chunkLoadRange = config.get("chunkLoadRange");
  const [newChunks, newActiveChunks] = calculateChunks(
    scene.activeCamera.position.x,
    scene.activeCamera.position.y,
    scene.activeCamera.position.z,
    chunkDiameter,
    chunkLoadRange
  );
  for (const chunk of newActiveChunks) {
    loadChunk(scene, newChunks, chunk);
  }
  return [newChunks, newActiveChunks, chunkDiameter, chunkLoadRange];
}

function clearChunks(scene, chunks, activeChunks) {
  for (const chunk of activeChunks) {
    loadChunk(scene, chunks, chunk, false);
  }
}

function calculateChunks(posX, posY, posZ, chunkDiameter, chunkLoadRange) {
  // Each key in this object corresponds to a chunk's centroid in string form.
  const chunks = {};

  // Define the "sight" range: the distance threshold (in world units) used to decide whether a chunk is visible.
  // The value includes an extra 0.5 chunk diameter margin.
  const sight = (chunkLoadRange + 0.5) * chunkDiameter;

  chunks.scale = config.get("scale");
  chunks.tissues = [config.get("tissueX"), config.get("tissueY"), config.get("tissueZ")]

  const familiesToShow = config.get("shownFamilies") ?? dataHandler.families;

  // Loop through each family available in the data handler.
  // For each family, iterate through the genes for specific tissues
  // and create an instance of the outlier mesh for each data point.
  for (const family of familiesToShow) {
    const geneCount = dataHandler.getGeneCount(family);
    for (let geneIndex = 0; geneIndex < geneCount; geneIndex++) {
      const { coordinates, is_outlier } = dataHandler.getGeneData(family, geneIndex, chunks.tissues, ["is_outlier"]);
      const scaled = coordinates.map((v) => v*chunks.scale);

      // Determine the centroid of the chunk this position falls into.
      // getChunkCentroid is assumed to return an array-like coordinate (e.g. [x, y, z])
      // which is also used as a key in the `chunks` object.
      const chunk = getChunkCentroid(scaled, chunkDiameter);

      if (chunks[chunk] === undefined) {
        chunks[chunk] = [new Map(), 0, 0, null, null, null];
      }

      const genes = chunks[chunk][0];
      if (!genes.has(family)) {
        genes.set(family, [[], []]);
      }
      genes.get(family)[is_outlier ? 1 : 0].push(geneIndex);
      chunks[chunk][1] += !is_outlier;
      chunks[chunk][2] += is_outlier;
    }
  }

  const activeChunks = [];
  for (let x = posX - sight; x <= posX + sight; x+=chunkDiameter) {
    for (let y = posY - sight; y <= posY + sight; y+=chunkDiameter) {
      for (let z = posZ - sight; z <= posZ + sight; z+=chunkDiameter) {
        const chunk = getChunkCentroid([x, y, z], chunkDiameter);
        if (
          Math.abs(chunk[0] - posX) < sight &&
          Math.abs(chunk[1] - posY) < sight &&
          Math.abs(chunk[2] - posZ) < sight
        ) {
          activeChunks.push(chunk.toString());
        }
      }
    }
  }

  return [ chunks, activeChunks ]
}

function getChunkCentroid([ x, y, z ], diameter) {
  function trim(a) {
    return Math.floor((a + diameter / 2) / diameter) * diameter;
  }
  return [trim(x), trim(y), trim(z)];
}

function loadChunk(scene, chunks, chunk, state=true) {
  const chunkData = chunks[chunk];
  if (chunkData) {
    if (state) {
      const [genes, inlierCount, outlierCount, ...meshes] = chunkData;

      for (const mesh of meshes) {
        mesh?.dispose();
      }

      chunkData[3] = create3DGridFromChunk(scene, chunk);

      // data points -- spheres
      const sphereDimensionsBuffer = new Float32Array(16 * inlierCount); // the translation buffer for one position takes 16 entries (it is a 4x4 rotation matrix)
      const sphereColorBuffer = new Float32Array(4 * inlierCount); // rgba
      if (sphereColorBuffer.length > 0) { // false if all members are outliers
        chunkData[4] = Mesh.Sphere(scene);
      }

      // outliers -- octahedrons
      const octDimensionsBuffer = new Float32Array(16 * outlierCount); // the translation buffer for one position takes 16 entries (it is a 4x4 rotation matrix)
      const octColorBuffer = new Float32Array(4 * outlierCount); // rgba
      if (outlierCount > 0) {
        chunkData[5] = Mesh.Octahedron(scene)
        chunkData[5].enableEdgesRendering();
        chunkData[5].edgesWidth = config.get("defaultDiameter") * 12;
        chunkData[5].edgesColor = Color(0, 0, 0, 1); // Black edges
        chunkData[5].edgesShareWithThinInstances = true;
      }

      let inlierIndex = 0;
      let outlierIndex = 0;
      for (const [family, [inliers, outliers]] of genes.entries()) {
        const familyColor = Color.FromHexString(config.get(`${family}_Color`) ?? dataHandler.getColor(family));
        for (const geneIndex of inliers) {
          const diameter = config.get(`${family}_Diameter`) ?? config.get("defaultDiameter");
          const pointData = dataHandler.getGeneData(family, geneIndex, chunks.tissues, []);
          fillThinInstanceBuffers(
            sphereDimensionsBuffer, inlierIndex * 16,
            sphereColorBuffer, inlierIndex * 4,
            diameter,
            pointData.coordinates.map(v => v * chunks.scale),
            familyColor
          );
          inlierIndex++;
        }
        const outlierColorHex = config.get(`${family}_OutlierColor`);
        const outlierColor = outlierColorHex === undefined ? familyColor : Color.FromHexString(outlierColorHex);
        for (const geneIndex of outliers) {
          const pointData = dataHandler.getGeneData(family, geneIndex, chunks.tissues, []);
          const diameter = config.get(`${family}_OutlierDiameter`) ?? config.get("defaultDiameter");
          fillThinInstanceBuffers(
            octDimensionsBuffer, outlierIndex * 16,
            octColorBuffer, outlierIndex * 4,
            diameter,
            pointData.coordinates.map(v => v * chunks.scale),
            outlierColor
          );
          outlierIndex++;
        }
      }
      if (inlierIndex !== inlierCount || outlierIndex !== outlierCount) {
        console.error("Misfilled buffers", inlierIndex, inlierCount, outlierIndex, outlierCount);
      }

      chunkData[4]?.thinInstanceSetBuffer("matrix", sphereDimensionsBuffer, 16);
      chunkData[4]?.thinInstanceSetBuffer("color", sphereColorBuffer, 4);
      chunkData[5]?.thinInstanceSetBuffer("matrix", octDimensionsBuffer, 16);
      chunkData[5]?.thinInstanceSetBuffer("color", octColorBuffer, 4);
    } else {
      if (chunkData[0].size > 0) {
        for (let i = chunkData.length - 3; i < chunkData.length; i++) {
          chunkData[i]?.dispose();
          chunkData[i] = null;
        }
      } else {
        delete chunks[chunk];
      }
    }
  } else {
    chunks[chunk] = [new Map(), 0, 0, create3DGridFromChunk(scene, chunk), null, null];
  }
}

function fillThinInstanceBuffers(dimensionsBuffer, dIndex, colorBuffer, cIndex, diameter, [x, y, z], color) {
  dimensionsBuffer[dIndex] = diameter; // set x scale
  dimensionsBuffer[dIndex + 5] = diameter; // set y scale
  dimensionsBuffer[dIndex + 10] = diameter; // set z scale

  dimensionsBuffer[dIndex + 12] = x;
  dimensionsBuffer[dIndex + 13] = y;
  dimensionsBuffer[dIndex + 14] = z;

  dimensionsBuffer[dIndex + 15] = 1;
  // the unchanged indices affect the rotation of the sphere -> zero 

  // setting color
  colorBuffer[cIndex++] = color.r;
  colorBuffer[cIndex++] = color.g;
  colorBuffer[cIndex++] = color.b;
  colorBuffer[cIndex] = color.a;
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

export { plotData };