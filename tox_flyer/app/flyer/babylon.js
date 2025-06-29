"use strict";

export function Vector(x=0, y=0, z=0) {
  return new BABYLON.Vector3(x, y, z);
}

export const Mesh = {
  Sphere(scene, name, options={}) {
    const mesh = BABYLON.MeshBuilder.CreateSphere(name, { diameter: 1, segments: 16, ...options }, scene);
    mesh.isPickable = false;
    mesh.freezeWorldMatrix();
    return mesh;
  },
  Octahedron(scene, name) {
    const mesh = BABYLON.MeshBuilder.CreatePolyhedron(name, { type: 2, size: 0.5, flat: false }, scene);
    mesh.isPickable = false;
    mesh.freezeWorldMatrix();
    return mesh;
  },
  setColor(sphere, color) {
    sphere.material.unfreeze();
    sphere.material.diffuseColor = typeof color === "string" ? Color.FromHexString(color) : color;
    sphere.material.alpha = sphere.material.diffuseColor.a;
    sphere.material.freeze();
  },
  setSize(sphere, diameter) {
    sphere.unfreezeWorldMatrix();
    sphere.scaling = Vector(diameter, diameter, diameter);
    sphere.freezeWorldMatrix();
  }
};
Object.freeze(Mesh);

export function OrbitCam(scene, name) {
  const cam = new BABYLON.ArcRotateCamera(name, null, null, 10, Vector(), scene);
  cam.lowerRadiusLimit = 1;
  return cam;
}

export function UniversalCam(scene, name) {
  return new BABYLON.UniversalCamera(name, Vector(), scene);
}

export function calcVectorDistance(v1, v2) {
  return BABYLON.Vector3.Distance(v1, v2);
}

export function WebGPUEngine(canvas, options) {
  return new BABYLON.WebGPUEngine(canvas, options);
}

export function WebGLEngine(canvas, options) {
  return new BABYLON.Engine(canvas, options.antialias ?? true, options);
}

export function Scene(engine) {
  return new BABYLON.Scene(engine);
}
Scene.FOGMODE_LINEAR = BABYLON.Scene.FOGMODE_LINEAR;

export function Color(r, g, b, a) {
  if (typeof r === "string") {
    return Color.FromHexString(r);
  } else if (r instanceof BABYLON.Color3 || r instanceof BABYLON.Color4) {
    return r;
  } else if (a === undefined) {
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

export function Material(scene, name, options={}) {
  const material = new BABYLON.StandardMaterial(name, scene);
  if (options.color) {
    Mesh.setColor({material}, options.color);
  }
  if (options.wireframe) {
    material.wireframe = true;
  }
  material.freeze();
  return material;
}
