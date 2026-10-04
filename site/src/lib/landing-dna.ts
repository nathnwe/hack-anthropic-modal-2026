import * as THREE from "three";

// Decorative, idealised DNA: not sequence-derived atomic coordinates or a
// molecular-dynamics simulation. Workbench data and geometry are independent.
const PAIRS = 480;
const HALF_STEPS = PAIRS * 2;
const UP = new THREE.Vector3(0, 1, 0);
const BLUE = new THREE.Color("#6fb3e8");
const COPPER = new THREE.Color("#d98a6a");
const BASE_BLUE = new THREE.Color("#bfd3e1");
const BASE_COPPER = new THREE.Color("#dbc3b3");
const OPEN = [
  [-65, -8, -2],
  [-44, -5, 1],
  [-31, -9, 2],
  [-23, -8, 1],
  [-25, 0, 2],
  [-22, 8, -3],
  [-12, 13, -4],
  [0, 15, -2],
  [14, 11, 2],
  [24, 3, 4],
  [26, -5, 3],
  [23, -9, 1],
  [34, -10, -1],
  [47, -6, -2],
  [65, -8, 0],
];
const CLOSED = [
  [-65, -7, -2],
  [-43, -6, 0],
  [-22, -13, -1],
  [-1.8, -10, 0],
  [-10, -2, -1],
  [-19, 8, -4],
  [-13, 18, -3],
  [0, 22, 0],
  [16, 17, 4],
  [21, 6, 6],
  [12, -3, 4],
  [1.8, -10, 1],
  [22, -13, 3],
  [43, -7, 0],
  [65, -6, -1],
];
const ease = (t: number) => t * t * (3 - 2 * t);

export class LandingDNA {
  private renderer: THREE.WebGLRenderer;
  private scene = new THREE.Scene();
  private camera = new THREE.OrthographicCamera(-55, 55, 26, -20, 0.1, 250);
  private model = new THREE.Group();
  private material = new THREE.MeshPhongMaterial({
    color: 0xffffff,
    shininess: 65,
    specular: 0x6b777e,
  });
  private spheres = new THREE.InstancedMesh(
    new THREE.SphereGeometry(1, 12, 8),
    this.material,
    PAIRS * 10 + 2,
  );
  private bonds = new THREE.InstancedMesh(
    new THREE.CylinderGeometry(1, 1, 1, 8),
    this.material,
    PAIRS * 11,
  );
  private curve = new THREE.CatmullRomCurve3(
    OPEN.map((p) => new THREE.Vector3(...p)),
    false,
    "centripetal",
  );
  private glow: THREE.Sprite;
  private glowTexture: THREE.CanvasTexture;
  private light = new THREE.PointLight(0xd8b75a, 0, 28, 2);
  private matrix = new THREE.Matrix4();
  private quaternion = new THREE.Quaternion();
  private scale = new THREE.Vector3();
  private delta = new THREE.Vector3();
  private middle = new THREE.Vector3();
  private observer: ResizeObserver;
  private abort = new AbortController();
  private playing = false;
  private reduced: boolean;
  private lost = false;
  private disposed = false;
  private time = 0;
  private last = 0;
  private request = 0;
  private width = 0;
  private height = 0;

  constructor(
    private canvas: HTMLCanvasElement,
    private panel: HTMLElement,
    reduced: boolean,
  ) {
    this.reduced = reduced;
    this.renderer = new THREE.WebGLRenderer({
      canvas,
      alpha: true,
      antialias: true,
      powerPreference: "low-power",
    });
    this.renderer.setPixelRatio(Math.min(devicePixelRatio || 1, 1.75));
    this.renderer.setClearColor(0x000000, 0);
    this.renderer.outputColorSpace = THREE.SRGBColorSpace;
    this.camera.position.set(0, 3, 100);
    this.camera.lookAt(0, 3, 0);
    this.scene.add(new THREE.HemisphereLight(0xd5e8fa, 0x292032, 1.6));
    const key = new THREE.DirectionalLight(0xf3f7ff, 2.6);
    key.position.set(-20, 30, 60);
    const rim = new THREE.DirectionalLight(0x8fb9e8, 1.6);
    rim.position.set(25, 12, -20);
    const warm = new THREE.DirectionalLight(0xe9a482, 0.75);
    warm.position.set(-35, -10, 25);
    this.scene.add(key, rim, warm, this.light);
    this.spheres.instanceMatrix.setUsage(THREE.DynamicDrawUsage);
    this.bonds.instanceMatrix.setUsage(THREE.DynamicDrawUsage);
    this.spheres.frustumCulled = this.bonds.frustumCulled = false;
    this.model.add(this.spheres, this.bonds);
    this.scene.add(this.model);
    this.curve.arcLengthDivisions = 1800;

    const glowCanvas = document.createElement("canvas");
    glowCanvas.width = glowCanvas.height = 128;
    const context = glowCanvas.getContext("2d")!;
    const gradient = context.createRadialGradient(64, 64, 0, 64, 64, 64);
    gradient.addColorStop(0, "rgba(216,183,90,0.48)");
    gradient.addColorStop(0.28, "rgba(216,183,90,0.18)");
    gradient.addColorStop(1, "rgba(216,183,90,0)");
    context.fillStyle = gradient;
    context.fillRect(0, 0, 128, 128);
    this.glowTexture = new THREE.CanvasTexture(glowCanvas);
    this.glow = new THREE.Sprite(
      new THREE.SpriteMaterial({
        map: this.glowTexture,
        transparent: true,
        depthWrite: false,
        blending: THREE.AdditiveBlending,
      }),
    );
    this.glow.scale.set(28, 23, 1);
    this.glow.position.set(0, -10, -7);
    this.scene.add(this.glow);
    canvas.addEventListener(
      "webglcontextlost",
      (event) => {
        event.preventDefault();
        this.lost = true;
        cancelAnimationFrame(this.request);
        this.request = 0;
        panel.classList.remove("has-molecular-scene");
      },
      { signal: this.abort.signal },
    );
    canvas.addEventListener(
      "webglcontextrestored",
      () => {
        this.lost = false;
        this.render();
        this.setMotion(this.playing);
      },
      { signal: this.abort.signal },
    );
    this.observer = new ResizeObserver(() => this.resize());
    this.observer.observe(canvas);
    this.resize();
  }

  private resize() {
    this.width = this.canvas.clientWidth;
    this.height = this.canvas.clientHeight;
    if (!this.width || !this.height) return;
    this.renderer.setSize(this.width, this.height, false);
    const horizontal = this.width < 600 ? 54 : 110;
    const vertical = Math.max(45, (horizontal * this.height) / this.width);
    const worldWidth = (vertical * this.width) / this.height;
    this.camera.left = -worldWidth / 2;
    this.camera.right = worldWidth / 2;
    this.camera.top = vertical / 2;
    this.camera.bottom = -vertical / 2;
    this.camera.updateProjectionMatrix();
    this.render();
  }

  setReducedMotion(value: boolean) {
    this.reduced = value;
    if (value) this.setMotion(false);
    this.render();
  }

  setMotion(active: boolean) {
    this.playing = active && !this.reduced;
    this.panel.dataset.dnaMotion = this.playing ? "playing" : "paused";
    if (!this.playing) {
      cancelAnimationFrame(this.request);
      this.request = 0;
      return;
    }
    if (this.request || this.lost || this.disposed) return;
    this.last = performance.now();
    this.request = requestAnimationFrame(this.frame);
  }

  private frame = (now: number) => {
    this.request = 0;
    if (!this.playing || this.lost || this.disposed) return;
    const elapsed = now - this.last;
    if (elapsed >= 1000 / 30) {
      this.time += Math.min(elapsed / 1000, 0.08);
      this.last = now;
      this.render();
    }
    this.request = requestAnimationFrame(this.frame);
  };

  private render() {
    if (!this.width || !this.height || this.lost || this.disposed) return;
    const t = this.reduced ? 0 : this.time;
    // Slow approach, held contact, then relaxation, with no abrupt loop restart.
    const cycle = (1 - Math.cos((t * Math.PI * 2) / 22)) / 2;
    const contact = this.reduced ? 1 : ease(Math.min(1, cycle * 1.22));
    const narrow = this.width < 600;
    for (let i = 0; i < this.curve.points.length; i++) {
      const a = OPEN[i],
        b = CLOSED[i];
      const anchor = i === 3 || i === 11;
      const flexibility = anchor ? 0.16 * (1 - contact) : 1;
      const x = a[0] + (b[0] - a[0]) * contact;
      const y = a[1] + (b[1] - a[1]) * contact;
      const z = a[2] + (b[2] - a[2]) * contact;
      // Correlated travelling bends, not rigid pieces translating together.
      const wave = this.reduced ? 0 : flexibility;
      this.curve.points[i].set(
        x * (narrow ? 0.58 : 1) +
          (narrow ? Math.sign(x) * Math.exp(-Math.abs(x) / 9) * contact : 0) +
          wave * 0.6 * Math.sin(i * 0.8 - t * 0.55),
        y * (narrow ? 1.35 : 1) +
          wave *
            (0.85 * Math.sin(i * 0.72 - t * 0.7) +
              0.35 * Math.sin(i * 1.7 + t * 0.4)),
        z + wave * 0.85 * Math.sin(i * 0.93 + t * 0.46),
      );
    }
    this.curve.updateArcLengths();
    const points = this.curve.getSpacedPoints(HALF_STEPS);
    // Parallel transport avoids sudden helix-frame flips through tight bends.
    const normal = new THREE.Vector3(0, 0, 1);
    let previousA: THREE.Vector3 | undefined;
    let previousB: THREE.Vector3 | undefined;
    let sphereCount = 0,
      bondCount = 0;
    const sphere = (
      position: THREE.Vector3,
      radius: number,
      color: THREE.Color,
    ) => {
      this.matrix.makeScale(radius, radius, radius).setPosition(position);
      this.spheres.setMatrixAt(sphereCount, this.matrix);
      this.spheres.setColorAt(sphereCount++, color);
    };
    const bond = (
      a: THREE.Vector3,
      b: THREE.Vector3,
      radius: number,
      color: THREE.Color,
    ) => {
      this.delta.subVectors(b, a);
      const length = this.delta.length();
      this.quaternion.setFromUnitVectors(
        UP,
        this.delta.multiplyScalar(1 / length),
      );
      this.middle.addVectors(a, b).multiplyScalar(0.5);
      this.scale.set(radius, length, radius);
      this.matrix.compose(this.middle, this.quaternion, this.scale);
      this.bonds.setMatrixAt(bondCount, this.matrix);
      this.bonds.setColorAt(bondCount++, color);
    };
    for (let i = 0; i <= HALF_STEPS; i++) {
      const center = points[i];
      const tangent = new THREE.Vector3()
        .subVectors(
          points[Math.min(i + 1, HALF_STEPS)],
          points[Math.max(0, i - 1)],
        )
        .normalize();
      normal.addScaledVector(tangent, -normal.dot(tangent)).normalize();
      const binormal = new THREE.Vector3()
        .crossVectors(tangent, normal)
        .normalize();
      const phase = (i / 2 / 10.5) * Math.PI * 2;
      const offset = (angle: number) =>
        center
          .clone()
          .addScaledVector(normal, Math.cos(angle))
          .addScaledVector(binormal, Math.sin(angle));
      const a = offset(phase),
        b = offset(phase + Math.PI * 0.86);
      const rung = i % 2 === 0;
      sphere(a, rung ? 0.31 : 0.25, BLUE);
      sphere(b, rung ? 0.31 : 0.25, COPPER);
      if (previousA && previousB) {
        bond(previousA, a, 0.18, BLUE);
        bond(previousB, b, 0.18, COPPER);
      }
      if (rung && i < HALF_STEPS) {
        let previous = a;
        for (let j = 1; j <= 6; j++) {
          const point = a.clone().lerp(b, j / 7);
          const color = j <= 3 ? BASE_BLUE : BASE_COPPER;
          sphere(point, 0.21, color);
          bond(previous, point, 0.105, color);
          previous = point;
        }
        bond(previous, b, 0.105, BASE_COPPER);
      }
      previousA = a;
      previousB = b;
    }
    this.spheres.count = sphereCount;
    this.bonds.count = bondCount;
    this.spheres.instanceMatrix.needsUpdate =
      this.bonds.instanceMatrix.needsUpdate = true;
    this.spheres.instanceColor!.needsUpdate =
      this.bonds.instanceColor!.needsUpdate = true;
    this.model.rotation.y = this.reduced ? 0.06 : 0.06 * Math.sin(t * 0.19);
    (this.glow.material as THREE.SpriteMaterial).opacity =
      ease(Math.max(0, (contact - 0.6) / 0.4)) * 0.7;
    this.light.position.set(0, -10, 5);
    this.light.intensity = contact * 7;
    this.renderer.render(this.scene, this.camera);
    if (!this.panel.classList.contains("has-molecular-scene"))
      this.panel.classList.add("has-molecular-scene");
  }

  dispose() {
    this.disposed = true;
    cancelAnimationFrame(this.request);
    this.observer.disconnect();
    this.abort.abort();
    this.spheres.geometry.dispose();
    this.bonds.geometry.dispose();
    this.material.dispose();
    this.glow.material.dispose();
    this.glowTexture.dispose();
    this.renderer.dispose();
  }
}
