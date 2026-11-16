// Nube de células 3D: puntos redondos con color por estado, empaquetados en base64.
import { THREE } from "nicegui-scene";

function decode(b64) {
  const bin = atob(b64);
  const out = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) out[i] = bin.charCodeAt(i);
  return out;
}

let discTexture = null;
function disc() {
  if (discTexture) return discTexture;
  const c = document.createElement("canvas");
  c.width = c.height = 64;
  const g = c.getContext("2d");
  g.beginPath();
  g.arc(32, 32, 29, 0, 2 * Math.PI);
  g.fillStyle = "#ffffff";
  g.fill();
  discTexture = new THREE.CanvasTexture(c);
  return discTexture;
}

export default class CellCloud {
  mesh;

  create_mesh(payload, palette, size) {
    this.palette = palette.map((hex) => new THREE.Color(hex || "#000000"));
    const geometry = new THREE.BufferGeometry();
    const material = new THREE.PointsMaterial({
      size: size,
      vertexColors: true,
      map: disc(),
      alphaTest: 0.5,
      sizeAttenuation: true,
    });
    this.mesh = new THREE.Points(geometry, material);
    this.set_cells(payload);
    return this.mesh;
  }

  set_cells(payload) {
    const pos = decode(payload.positions);
    const states = decode(payload.states);
    const n = states.length;
    const [cz, cy, cx] = payload.center;
    const position = new Float32Array(n * 3);
    const color = new Float32Array(n * 3);
    for (let i = 0; i < n; i++) {
      // modelo (z, y, x) -> escena (x, z, y): el eje del conducto queda horizontal
      position[i * 3] = pos[i * 3] - cz;
      position[i * 3 + 1] = pos[i * 3 + 2] - cx;
      position[i * 3 + 2] = pos[i * 3 + 1] - cy;
      const c = this.palette[states[i]] || this.palette[0];
      color[i * 3] = c.r;
      color[i * 3 + 1] = c.g;
      color[i * 3 + 2] = c.b;
    }
    const geometry = this.mesh.geometry;
    geometry.setAttribute("position", new THREE.BufferAttribute(position, 3));
    geometry.setAttribute("color", new THREE.BufferAttribute(color, 3));
    geometry.computeBoundingSphere();
  }

  set_size(size) {
    this.mesh.material.size = size;
  }

  apply_material({ opacity }) {
    // el color viene de cada punto; solo se respeta la opacidad
    this.mesh.material.transparent = opacity < 1;
    this.mesh.material.opacity = opacity;
  }
}
