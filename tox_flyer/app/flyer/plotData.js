"use strict";

import {
  getChunks,
  createDynamicThinInstance,
  removeDynamicThinInstance,
  setupDynamicThinInstanceMesh,
  dynamicThinInstanceBufferUpdated
} from "./chunks.js";
import { createTooltip, removeTooltip } from "./gui.js";
import {
  Mesh,
  Color,
  Vector,
  TransformNode,
  Material,
  fillThinInstanceBuffers,
  decomposeMatrix,
  getInstanceMatrix,
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
  setupSelectionMeshes(scene);
  setupFamilyHullMesh(scene);
  setupShiftVectorMesh(scene);

  document.dispatchEvent(new CustomEvent("initialTrigger", { detail: [
    "Hull",
    "ShiftVector",
    "PickedGene",
    "PickedCentroid",
    "PickedShiftVector"
  ]}));

  const chunks = getChunks(scene);

  let picked = null;
  scene.onPointerObservable.add((evt) => {
    switch (evt.type) {
      case BABYLON.PointerEventTypes.POINTERDOWN:
        removeTooltip();
        break;
      case BABYLON.PointerEventTypes.POINTERTAP: {
        function handlePick(picked, onPick) {
          const keyType = `Picked${picked.type}`;
          const configAttr = `${picked.family}_${keyType}` + (picked.geneIndex === undefined ? "" : `:${picked.geneIndex}`);
          const isPicked = config.get(configAttr);
          if (!isPicked) {
            onPick?.();
          }
          
          config.set(configAttr, !isPicked);
        }
        picked = pickFromMeshes(chunks);
        if (picked !== null) {
          switch (picked.type) {
            case "ShiftVector": {
              handlePick(picked);
              break;
            }
            case "Gene": {
              handlePick(picked, () => {
                const geneData = dataHandler.getGeneData(picked.family, picked.geneIndex, chunks.tissues, ["genes", "species", "is_outlier"]);
                createTooltip(evt.event.clientX, evt.event.clientY,
                  `<center>${geneData.is_outlier ? "Outlier" : "Inlier"}</center><table><tbody>` +
                  `<tr><td>Gene:</td><td>${geneData.genes}</td></tr>` +
                  `<tr><td>Species:</td><td>${geneData.species}</td></tr>` +
                  `<tr><td>${chunks.tissues[0]}:</td><td>${geneData.coordinates[0].toFixed(2)}</td></tr>` +
                  `<tr><td>${chunks.tissues[1]}:</td><td>${geneData.coordinates[1].toFixed(2)}</td></tr>` +
                  `<tr><td>${chunks.tissues[2]}:</td><td>${geneData.coordinates[2].toFixed(2)}</td></tr>` +
                  "</tbody></table>"
                );
              });
              break;
            }
            case "Centroid": {
              handlePick(picked, () => {
                const { centroid } = dataHandler.getFamilyData(picked.family, ...chunks.tissues);
                createTooltip(evt.event.clientX, evt.event.clientY,
                  "Centroid<table><tbody>" +
                  `<tr><td>Family:</td><td>${dataHandler.getFamilyIDs(picked.family)[0]}</td></tr>` +
                  `<tr><td>${chunks.tissues[0]}:</td><td>${centroid[0].toFixed(2)}</td></tr>` +
                  `<tr><td>${chunks.tissues[1]}:</td><td>${centroid[1].toFixed(2)}</td></tr>` +
                  `<tr><td>${chunks.tissues[2]}:</td><td>${centroid[2].toFixed(2)}</td></tr>` +
                  "</tbody></table>"
                );
              });
              break;
            }
          }
        } else {
          removeTooltip();
        }
        break;
      }
      case BABYLON.PointerEventTypes.POINTERDOUBLETAP:
        // // should pick/unpick all instances of a family, TODO, not yet working
        // if (picked !== null) {
        //   // if unpick was successful, so it was picked already, the original pick was initiated by the simultaneously triggered POINTERTAP
        //   const wasUnselected = unpickInstance(picked);
        //   if (!wasUnselected) {
        //     unpickInstance(picked);
        //   } else {
        //     pickInstance(picked);
        //     createTooltip(evt.event.clientX, evt.event.clientY,
        //       "Family<table><tbody>" +
        //       `<tr><td>Family:</td><td>${dataHandler.getFamilyIDs(picked.family)[0]}</td></tr>` +
        //       `<tr><td>Members:</td><td>${dataHandler.getGeneCount(picked.family)}</td></tr>` +
        //       "</tbody></table>"
        //     );
        //   }
        // }
        break;
    }
  });
}

function setupFamilyHullMesh(scene) {
  const hull = BABYLON.MeshBuilder.CreateCapsule("hull", {height: 1, radius: 1/3, subdivisions: 2, capSubdivisions: 3}, scene);
  setupDynamicThinInstanceMesh(hull);

  hull.material = Material(scene, null, {wireframe: true});

  function createHull(evt) {
    const { family, gene, value } = evt.detail;
    if (value) {
      const scale = config.get("scale");
      const familyData = dataHandler.getFamilyData(family, config.get("tissueX"), config.get("tissueY"), config.get("tissueZ"));
      const centroid = Vector(...familyData.centroid.map(v => v*scale));
      const stdDevs = familyData.stdDevs.map(v => v*scale);
      const color = Color(config.get(family + "_Color")).scale(2);
      color.a /= 2;

      const instanceMatrix = getInstanceMatrix(
        centroid,
        Vector(stdDevs[0] * 3, stdDevs[1] * 2, stdDevs[2] * 3)
      );

      createDynamicThinInstance(hull, family, undefined, instanceMatrix, color);
    } else {
      removeDynamicThinInstance(hull, family);
    }
  }

  document.addEventListener("Hull", createHull);
  document.addEventListener("HullUpdated", () => dynamicThinInstanceBufferUpdated(hull));

  function recreate() {
    for (const { family } of hull.TOX_metadata) {
      config.set(`${family}_Hull`, true, false);
    }
    document.dispatchEvent(new CustomEvent("HullUpdated"));
  }
  for (const setting of ["tissueX", "tissueY", "tissueZ", "scale", "Color"]) {
    document.addEventListener(setting, recreate);
  }
}

function createVectorPartsInstanceMatrices(family, geneIndex, grow=0) {
  const scale = config.get("scale");
  const tissues = [config.get("tissueX"), config.get("tissueY"), config.get("tissueZ")];

  const centroid = Vector(...dataHandler.getFamilyData(family, ...tissues).centroid.map(v => v*scale));
  const sphereDiameter = config.get(`${family}_Diameter`);

  const { coordinates } = dataHandler.getGeneData(family, geneIndex, tissues, []);
  const genePos = Vector(...coordinates.map(v => v*scale));
  const direction = genePos.subtract(centroid);
  const vectorLength = direction.length() - sphereDiameter / 2;
  const shaftLengthScale = 1 - 2 * sphereDiameter / vectorLength;
  const shaftPosition = centroid.add(direction.scale(shaftLengthScale / 2));
  const headPosition = centroid.add(direction.scale(shaftLengthScale + sphereDiameter / vectorLength / 2));

  return {
    shaft: getInstanceMatrix(
      shaftPosition,
      Vector(sphereDiameter / 2 + grow, vectorLength * shaftLengthScale + grow, sphereDiameter / 2 + grow),
      genePos
    ),
    head: getInstanceMatrix(
      headPosition,
      Vector(sphereDiameter + grow, sphereDiameter * 2 + grow, sphereDiameter + grow),
      genePos
    )
  };
}

function setupShiftVectorMesh(scene) {
  const shiftVectorShaft = Mesh.Cylinder(scene, "shiftVectorShaft");
  setupDynamicThinInstanceMesh(shiftVectorShaft);
  const shiftVectorHead = Mesh.Cone(scene, "shiftVectorHead");
  setupDynamicThinInstanceMesh(shiftVectorHead);

  document.addEventListener("ShiftVector", evt => {
    const { family, gene, value } = evt.detail;
    if (value) {
      const matrices = createVectorPartsInstanceMatrices(family, gene);

      let color = Color(config.get(family + "_Color"));
      const colorScale = 1 / Math.max(color.r, color.g, color.b);
      color = color.scale(colorScale);
      color.a /= colorScale;
      createDynamicThinInstance(shiftVectorShaft, family, gene, matrices.shaft, color);
      createDynamicThinInstance(shiftVectorHead, family, gene, matrices.head, color);
    } else {
      removeDynamicThinInstance(shiftVectorHead, family, gene);
      removeDynamicThinInstance(shiftVectorShaft, family, gene);
    }
  });

  document.addEventListener("ShiftVectorUpdated", () => {
    dynamicThinInstanceBufferUpdated(shiftVectorHead);
    dynamicThinInstanceBufferUpdated(shiftVectorShaft);
  });

  function recreate() {
    for (const { family, geneIndex } of shiftVectorHead.TOX_metadata) {
      config.set(`${family}_ShiftVector:${geneIndex}`, true, false);
    }
    document.dispatchEvent(new CustomEvent("ShiftVectorUpdated"));
  }
  for (const setting of ["tissueX", "tissueY", "tissueZ", "scale", "Diameter", "defaultDiameter", "Color"]) {
    document.addEventListener(setting, recreate);
  }
}

function setupSelectionMeshes(scene) {
  const meshes = [
    Mesh.Sphere(scene, "pickedSphere"),
    Mesh.Octahedron(scene, "pickedOctahedron"),
    Mesh.Cylinder(scene, "pickedVectorShaft"),
    Mesh.Cone(scene, "pickedVectorHead")
  ];

  const highlightLayer = new BABYLON.HighlightLayer("highlight", scene);

  {
    function setHighlightColor(evt) {
      const color = Color(evt.detail);
      for (const mesh of meshes) {
        highlightLayer.removeMesh(mesh);
        highlightLayer.addMesh(mesh, color);
      }
    }
    document.addEventListener("selectedDataPointColor", setHighlightColor);
    setHighlightColor({ detail: config.get("selectedDataPointColor") });
  }

  const material =  Material(scene, null, {color: Color(0, 0, 0, 0)});

  for (const mesh of meshes) {
    setupDynamicThinInstanceMesh(mesh, false);
    mesh.material = material;
    highlightLayer.setEffectIntensity(mesh, 0.7);
  }

  setupGenePicking(scene);
  setupCentroidPicking(scene);
  setupVectorPicking(scene);
  {
    function repickGenesAndCentroids() {
      for (const meshName of ["pickedSphere", "pickedOctahedron"]) {
        const mesh = scene.getMeshByName(meshName);
        for (const { family, geneIndex } of mesh.TOX_metadata) {
          if (geneIndex !== undefined) {
            config.set(`${family}_PickedGene:${geneIndex}`, true, false);
          } else {
            config.set(`${family}_PickedCentroid`, true, false);
          }
        }
      }
      document.dispatchEvent(new CustomEvent("PickedCentroidUpdated"));
      document.dispatchEvent(new CustomEvent("PickedGeneUpdated"));
    }
    for (const setting of ["tissueX", "tissueY", "tissueZ", "scale", "OutlierDiameter", "Diameter", "defaultDiameter"]) {
      document.addEventListener(setting, repickGenesAndCentroids);
    }
  }

}

function selectionMeshPick(selectionMesh, family, gene, type, instanceMatrix, dispatchEvent=true) {
  if (instanceMatrix !== undefined) {
    const instanceCount = selectionMesh.TOX_instanceCount;
    createDynamicThinInstance( selectionMesh, family, gene, instanceMatrix );
    if (instanceCount !== selectionMesh.TOX_instanceCount && dispatchEvent) {
      document.dispatchEvent(new CustomEvent("pick", { detail: { family, gene, type } }));
    }
  } else {
    if (removeDynamicThinInstance(selectionMesh, family, gene) && dispatchEvent) {
      document.dispatchEvent(new CustomEvent("unpick", { detail: { family, gene, type } }));
    }
  }
}

function setupGenePicking(scene) {
  function pick(evt) {
    const { family, gene, value } = evt.detail;
    const { coordinates, is_outlier } = dataHandler.getGeneData(family, gene, [config.get("tissueX"), config.get("tissueY"), config.get("tissueZ")], ["is_outlier"]);
    const mesh = scene.getMeshByName(`picked${is_outlier ? "Octahedron" : "Sphere"}`);
    if (value) {
      const diameter = (config.get(is_outlier ? `${family}_OutlierDiameter` : `${family}_Diameter`)) + .001;
      const scale = config.get("scale");
      selectionMeshPick(
        mesh,
        family,
        gene,
        "Gene",
        getInstanceMatrix(
          Vector(...coordinates.map(v => v*scale)),
          Vector(diameter, diameter, diameter)
        )      
      );
    } else {
      selectionMeshPick(mesh, family, gene, "Gene");
    }
  }
  document.addEventListener("PickedGene", pick);
  for (const meshName of ["pickedSphere", "pickedOctahedron"]) {
    document.addEventListener("PickedGeneUpdated", () => dynamicThinInstanceBufferUpdated(scene.getMeshByName(meshName)));
  }
}

function setupCentroidPicking(scene) {
  function pick(evt) {
    const { family, value } = evt.detail;
    const mesh = scene.getMeshByName("pickedSphere");
    if (value) {
      const { centroid } = dataHandler.getFamilyData(family, config.get("tissueX"), config.get("tissueY"), config.get("tissueZ"));
      const diameter = config.get(`${family}_Diameter`) * 4 + .001;
      const scale = config.get("scale");
      selectionMeshPick(
        mesh,
        family,
        undefined,
        "Centroid",
        getInstanceMatrix(
          Vector(...centroid.map(v => v*scale)),
          Vector(diameter, diameter, diameter)
        )      
      );
    } else {
      selectionMeshPick(mesh, family, undefined, "Centroid");
    }
  }
  document.addEventListener("PickedCentroid", pick);
  document.addEventListener("PickedCentroidUpdated", () => dynamicThinInstanceBufferUpdated(scene.getMeshByName("pickedSphere")));
}

function setupVectorPicking(scene) {
  function pick(evt) {
    const { family, gene, value } = evt.detail;
    if (value) {
      const matrices = createVectorPartsInstanceMatrices(family, gene, .001);
      selectionMeshPick(scene.getMeshByName("pickedVectorShaft"), family, gene, "ShiftVector", matrices.shaft);
      selectionMeshPick(scene.getMeshByName("pickedVectorHead"), family, gene, "ShiftVector", matrices.head, false);
    } else {
      selectionMeshPick(scene.getMeshByName("pickedVectorShaft"), family, gene, "ShiftVector");
      selectionMeshPick(scene.getMeshByName("pickedVectorHead"), family, gene, "ShiftVector", undefined, false);
    }
  }
  document.addEventListener("PickedShiftVector", pick);
  document.addEventListener("PickedShiftVectorUpdated", () => {
    dynamicThinInstanceBufferUpdated(scene.getMeshByName("pickedVectorShaft"));
    dynamicThinInstanceBufferUpdated(scene.getMeshByName("pickedVectorHead"));
  });

  function repick() {
    for (const meshName of ["pickedVectorShaft", "pickedVectorHead"]) {
      const mesh = scene.getMeshByName(meshName);
      for (const { family, geneIndex } of mesh.TOX_metadata) {
        config.set(`${family}_PickedShiftVector:${geneIndex}`, true, false);
      }
    }
    document.dispatchEvent(new CustomEvent("PickedShiftVectorUpdated"));
  }
  for (const setting of ["tissueX", "tissueY", "tissueZ", "scale", "Diameter", "defaultDiameter"]) {
    document.addEventListener(setting, repick);
  }
}

function intersectsNonSpherical(
  worldMatrix,    // Matrix of the instance
  worldRay,       // Ray in world space
  scaling
) {
  // 1. Build local‐space AABB
  const max = Vector(.5, .5, .5);
  const min = max.negate();

  // 2. Transform ray into local space
  const inv    = BABYLON.Matrix.Invert(worldMatrix);
  const localR = new BABYLON.Ray(
    BABYLON.Vector3.TransformCoordinates(worldRay.origin, inv),
    BABYLON.Vector3.TransformNormal(worldRay.direction, inv).normalize(),
    worldRay.length
  );

  // 3. Test box
  return localR.intersectsBoxMinMax(min, max);
}

function meshHit(ray, mesh, maxDistance=Infinity) {
  const sphereMatrices = mesh.thinInstanceGetWorldMatrices(); 

  let picked = null;

  for (const i in sphereMatrices) {
    const { position, scaling, rotation } = decomposeMatrix(sphereMatrices[i]);

    const distance = Vector.Distance(ray.origin, position);
    if (distance < (picked?.distance ?? maxDistance)) {
      let intersects;
      if (["sphere", "octahedron"].includes(mesh.TOX_shape)) {
        intersects = ray.intersectsSphere(
          { center: position, radius: scaling.x / 2 }
        );
      } else {
        intersects = intersectsNonSpherical(sphereMatrices[i], ray, scaling);
      }
      if (intersects) {
        picked = {
          index: i,
          distance
        }
      }
    }
  }

  return picked;
}

function pickFromMeshes(chunks) {
  const pickRay = chunks.scene.createPickingRay(
    chunks.scene.pointerX, chunks.scene.pointerY,
    BABYLON.Matrix.Identity(),    // you can pass other transforms if you want
    chunks.scene.activeCamera
  );

  let picked = null;
  for (const chunkCentroid of chunks.active) {
    const meshes = chunks.chunks.get(chunkCentroid)[2];

    for (const [meshType, mesh] of meshes) {
      const hit = meshHit(pickRay, mesh, picked?.distance);
      if (hit !== null) {
        picked = hit;
        picked.meshType = meshType;
        picked.chunkCentroid = chunkCentroid;
      }
    }
  }

  const unchunkedMeshes = {
    arrow: ["shiftVectorShaft", "shiftVectorHead"]
  };
  for (const [meshType, meshNames] of Object.entries(unchunkedMeshes)) {
    let hit = null;
    for (const meshName of meshNames) {
      const mesh = chunks.scene.getMeshByName(meshName);
      if (mesh) {
        hit = meshHit(pickRay, mesh, picked?.distance);
        if (hit !== null) {
          const { family, geneIndex } = mesh.TOX_metadata[hit.index];
          picked = { ...hit, meshType, family, geneIndex };
          break;
        }
      }
    }
  }

  if (picked) {
    let pickedIndex = picked.index;
    switch (picked.meshType) {
      case "arrow": {
        picked.type = "ShiftVector";
        break;
      }
      case "centroids": {
        const [genes] = chunks.chunks.get(picked.chunkCentroid);
        picked.type = "Centroid";
        for (const [family, members] of genes) {
          pickedIndex -= Boolean(members.centroids);
          if (pickedIndex === -1) {
            picked.family = family;
            break;
          }
        }
        break;
      }
      case "inliers":
      case "outliers": {
        const [genes] = chunks.chunks.get(picked.chunkCentroid);
        picked.type = "Gene";
        for (const [family, members] of genes) {
          const indices = members[picked.meshType];
          if (indices) {
            pickedIndex -= indices.length;
            if (pickedIndex < 0) {
              picked.family = family;
              picked.geneIndex = indices[pickedIndex + indices.length];
              break;
            }
          }
        }
        break;
      }
    }
  }
  return picked;
}

export { plotData };