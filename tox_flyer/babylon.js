"use strict"

export function SphereMesh(scene, name) {
  return BABYLON.MeshBuilder.CreateSphere(name, { diameter: 1, segments: 16 }, scene);
}

export function Octahedron(scene, name) {
  return BABYLON.MeshBuilder.CreatePolyhedron(name, { type: 2, size: 0.5, flat: false }, scene);
}

export function Vector(x, y, z) {
  return new BABYLON.Vector3(x, y, z);
}

export function OrbitCam(scene, name) {
  const cam = new BABYLON.ArcRotateCamera(name, null, null, 10, Vector(0, 0, 0), scene);
  cam.lowerRadiusLimit = 1;
  return cam;
}

export function UniversalCam(scene, name) {
  return new BABYLON.UniversalCamera(name, Vector(0, 0, 0), scene);
}

export function calcVectorDistance(v1, v2) {
  return BABYLON.Vector3.Distance(v1, v2);
}

export function WebGPUEngine(canvas) {
  return new BABYLON.WebGPUEngine(canvas);
}

export function WebGLEngine(canvas) {
  return new BABYLON.WebGPUEngine(canvas);
}

export function Scene(engine) {
  return new BABYLON.Scene(engine);
}
Scene.FOGMODE_LINEAR = BABYLON.Scene.FOGMODE_LINEAR;

export function Color(r, g, b, a) {
  if (a === undefined) {
    return new BABYLON.Color3(r, g, b);
  } else {
    return new BABYLON.Color4(r, g, b, a);
  }
}
Color.FromHexString = (hexString) => {
    return BABYLON.Color4.FromHexString(hexString);
}

export function Viewport(left, top, width, height) {
  return new BABYLON.Viewport(left, top, width, height);
}

export function Light(scene, name, direction) {
  return new BABYLON.HemisphericLight(name, direction, scene);
}

export function TransformNode(scene, name) {
  return new BABYLON.TransformNode(name, scene);
}