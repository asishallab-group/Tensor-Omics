"use strict";

import {
  getChunks,
  pickInstance,
  unpickInstance,
  createDynamicThinInstance,
  removeDynamicThinInstance,
  setupDynamicThinInstanceMesh
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

  const chunks = getChunks(scene);

  document.addEventListener("chunkReload", (evt) => {
    chunks.recalculate();
    chunks.load();
  })

  let picked = null;
  scene.onPointerObservable.add((evt) => {
    switch (evt.type) {
      case BABYLON.PointerEventTypes.POINTERDOWN:
        removeTooltip();
        break;
      case BABYLON.PointerEventTypes.POINTERTAP: {
        function handlePick(picked, onPick, dispatchEvent=true) {
          const configAttr = `${picked.family}_Picked${picked.type}` + (picked.geneIndex === undefined ? "" : `:${picked.geneIndex}`);
          if (!unpickInstance(picked, dispatchEvent)) {
            pickInstance(picked, dispatchEvent);
            onPick?.();
            config.set(configAttr, true, false);
          } else {
            config.set(configAttr, false, false);
          }
        }
        picked = pickFromMeshes(chunks);
        if (picked !== null) {
          switch (picked.type) {
            case "ShiftVector": {
              let dispatchEvent = true;
              for (const pickedPart of picked.parts) {
                handlePick(pickedPart, null, dispatchEvent);
                dispatchEvent = false;
              }
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
  hull.TOX_create = function (family, tissues, scale) {
    const familyData = dataHandler.getFamilyData(family, ...tissues);
    const centroid = Vector(...familyData.centroid.map(v => v*scale));
    const stdDevs = familyData.stdDevs.map(v => v*scale);
    const color = Color(config.get(family + "_Color")).scale(2);
    color.a /= 2;

    const instanceMatrix = getInstanceMatrix(
      centroid,
      Vector(stdDevs[0] * 3, stdDevs[1] * 2, stdDevs[2] * 3)
    );

    createDynamicThinInstance(this, family, undefined, instanceMatrix, color);
  }

  hull.TOX_remove = function (family) {
    removeDynamicThinInstance(this, family);
  }
}

function setupShiftVectorMesh(scene) {
  const shiftVectorShaft = Mesh.Cylinder(scene, "shiftVectorShaft");
  setupDynamicThinInstanceMesh(shiftVectorShaft);
  const shiftVectorHead = Mesh.Cone(scene, "shiftVector");
  setupDynamicThinInstanceMesh(shiftVectorHead);

  shiftVectorHead.TOX_create = function create(family, geneIndex, tissues, scale) {
    let color = Color(config.get(family + "_Color"));
    const colorScale = 1 / Math.max(color.r, color.g, color.b);
    color = color.scale(colorScale);
    color.a /= colorScale;

    const centroid = Vector(...dataHandler.getFamilyData(family, ...tissues).centroid.map(v => v*scale));
    const sphereDiameter = config.get(`${family}_Diameter`);

    const { coordinates } = dataHandler.getGeneData(family, geneIndex, tissues, []);
    const genePos = Vector(...coordinates.map(v => v*scale));
    const direction = genePos.subtract(centroid);
    const vectorLength = direction.length() - sphereDiameter / 2;
    const shaftLengthScale = 1 - 2 * sphereDiameter / vectorLength;
    const shaftPosition = centroid.add(direction.scale(shaftLengthScale / 2));
    const headPosition = centroid.add(direction.scale(shaftLengthScale + sphereDiameter / vectorLength / 2));

    // create shaft
    const shaftInstanceMatrix = getInstanceMatrix(
      shaftPosition,
      Vector(sphereDiameter / 2, vectorLength * shaftLengthScale, sphereDiameter / 2),
      genePos
    );
    // create head
    const headInstanceMatrix = getInstanceMatrix(
      headPosition,
      Vector(sphereDiameter, sphereDiameter * 2, sphereDiameter),
      genePos
    );

    for (const mesh of [this, shiftVectorShaft]) {
      const isThis = mesh === this;
      const matrix = isThis ? headInstanceMatrix : shaftInstanceMatrix;
      createDynamicThinInstance(mesh, family, geneIndex, matrix, color);
      if (config.get(`${family}_PickedShiftVector:${geneIndex}`)) {
        pickInstance({
          mesh,
          family,
          geneIndex,
          ...decomposeMatrix(matrix),
          type: "ShiftVector"
        }, isThis);  // fire pick event only once for shaft
      }
    }
  }

  shiftVectorHead.TOX_remove = function (family, gene) {
    removeDynamicThinInstance(this, family, gene);
    removeDynamicThinInstance(shiftVectorShaft, family, gene);
  }
}

function setupSelectionMeshes(scene) {
  const highlightLayer = new BABYLON.HighlightLayer("highlight", scene);

  const material =  Material(scene, null, {color: Color(0, 0, 0, 0)});

  function setupMesh(mesh) {
    setupDynamicThinInstanceMesh(mesh, false);
    mesh.material = material;
    highlightLayer.setEffectIntensity(mesh, 0.7);
  }

  const meshes = [
    Mesh.Sphere(scene, "picked_sphere"),
    Mesh.Octahedron(scene, "picked_octahedron"),
    Mesh.Cylinder(scene, "picked_cylinder"),
    Mesh.Cone(scene, "picked_cone")
  ];

  for (const mesh of meshes) {
    setupMesh(mesh);
  }

  config.setSetterCallback("selectedDataPointColor", hexColorCode => {
    const color = Color(hexColorCode);
    for (const mesh of meshes) {
      highlightLayer.removeMesh(mesh);
      highlightLayer.addMesh(mesh, color);
    }
  });
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
          position,
          rotation,
          scaling,
          index: i,
          mesh,
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
    arrow: ["shiftVectorShaft", "shiftVector"]
  };
  for (const [meshType, meshNames] of Object.entries(unchunkedMeshes)) {
    let hit = null;
    for (const meshName of meshNames) {
      const mesh = chunks.scene.getMeshByName(meshName);
      if (mesh) {
        hit = meshHit(pickRay, mesh, picked?.distance);
        if (hit !== null) {
          picked = { meshType, parts: [] };
          for (const meshName of meshNames) {
            const mesh = chunks.scene.getMeshByName(meshName);
            const worldMatrices = mesh.thinInstanceGetWorldMatrices();
            const decomposed = decomposeMatrix(worldMatrices[hit.index]);
            picked.parts.push({ ...hit, ...decomposed, mesh });
          }
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
        for (const pickedPart of picked.parts) {
          pickedPart.type = "ShiftVector";
          pickedPart.family = 76;
          pickedPart.geneIndex = 21;
        }
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