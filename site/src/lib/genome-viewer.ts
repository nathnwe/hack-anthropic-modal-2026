import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";
import {
  displayFrames,
  locusWindow,
  displayCenter,
  helixDisplay,
  helixSamples,
  type HelixDisplay,
  type LocusWindow,
} from "./dna-geometry";
import {
  formatNumber,
  regionLabel,
  type CRE,
  type GeneRecord,
} from "./records";

const STEPS = 600;
const GENE = new THREE.Color("#2877a5");
const ELEMENT = new THREE.Color("#ad5641");
const STRAND_A = new THREE.Color("#9baea0");
const STRAND_B = new THREE.Color("#c8b4ab");
const BASE_A = new THREE.Color("#e0e3da");
const BASE_B = new THREE.Color("#c1ccd0");
const FADE = new THREE.Color("#f9f8f2");
const Y_AXIS = new THREE.Vector3(0, 1, 0);
const byId = <T extends HTMLElement = HTMLElement>(id: string) =>
  document.getElementById(id) as T;

export class GenomeViewer {
  private root = byId("genome-viewer");
  private stage = byId("dna-stage");
  private scene = new THREE.Scene();
  private camera = new THREE.OrthographicCamera(-30, 30, 22, -22, 0.1, 500);
  private renderer: THREE.WebGLRenderer;
  private controls: OrbitControls;
  private model = new THREE.Group();
  private spheres: THREE.InstancedMesh;
  private bonds: THREE.InstancedMesh;
  private contact: THREE.Line;
  private omission: THREE.Line;
  private sphereRegions: Array<"gene" | "element" | null> = [];
  private frame = 0;
  private pending = false;
  private visible = true;
  private playing = false;
  private progress = 0;
  private lastTime = 0;
  private locus: LocusWindow | null = null;
  private display: HelixDisplay | null = null;
  private record: GeneRecord | null = null;
  private cre: CRE | null = null;
  private selectedKey = "";
  private observer: ResizeObserver;
  private visibility: IntersectionObserver;
  private abort = new AbortController();
  private reducedMotion = window.matchMedia("(prefers-reduced-motion: reduce)");
  private expanded = false;
  private oldOverflow = "";
  private disposed = false;
  private lostContext = false;
  private matrix = new THREE.Matrix4();
  private quaternion = new THREE.Quaternion();
  private scale = new THREE.Vector3();
  private raycaster = new THREE.Raycaster();
  private down = { x: 0, y: 0 };

  constructor() {
    this.renderer = new THREE.WebGLRenderer({
      antialias: true,
      alpha: true,
      powerPreference: "low-power",
    });
    this.renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
    this.renderer.setClearColor(0xf9f8f2, 0);
    this.renderer.outputColorSpace = THREE.SRGBColorSpace;
    this.renderer.setSize(
      this.stage.clientWidth || 600,
      this.stage.clientHeight || 420,
      false,
    );
    const canvas = this.renderer.domElement;
    canvas.tabIndex = 0;
    canvas.setAttribute("role", "img");
    canvas.setAttribute(
      "aria-label",
      "Interactive schematic DNA. Drag to rotate, scroll to zoom. Arrow keys pan; reset view restores the camera.",
    );
    this.stage.prepend(canvas);
    this.camera.position.set(8, 14, 85);
    this.controls = new OrbitControls(this.camera, canvas);
    this.controls.target.set(0, 2, 0);
    this.controls.enableDamping = true;
    this.controls.dampingFactor = 0.1;
    this.controls.minZoom = 0.65;
    this.controls.maxZoom = 5;
    this.controls.listenToKeyEvents(canvas);
    this.controls.update();
    this.controls.saveState();
    this.controls.addEventListener("change", () => this.requestFrame());
    this.scene.add(new THREE.HemisphereLight(0xffffff, 0x96988a, 2.1));
    const key = new THREE.DirectionalLight(0xffffff, 2.3);
    key.position.set(-20, 35, 60);
    this.scene.add(key);
    const fill = new THREE.DirectionalLight(0xe9f0ff, 0.8);
    fill.position.set(30, -8, -30);
    this.scene.add(fill);
    const material = new THREE.MeshPhongMaterial({
      color: 0xffffff,
      shininess: 40,
      specular: 0x555c58,
    });
    this.spheres = new THREE.InstancedMesh(
      new THREE.SphereGeometry(1, 14, 10),
      material,
      1800,
    );
    this.bonds = new THREE.InstancedMesh(
      new THREE.CylinderGeometry(1, 1, 1, 8),
      material,
      2000,
    );
    this.spheres.instanceMatrix.setUsage(THREE.DynamicDrawUsage);
    this.bonds.instanceMatrix.setUsage(THREE.DynamicDrawUsage);
    // The deformed model's bounds change. Culling stale instance bounds can hide it.
    this.spheres.frustumCulled = this.bonds.frustumCulled = false;
    this.model.add(this.spheres, this.bonds);
    this.contact = new THREE.Line(
      new THREE.BufferGeometry(),
      new THREE.LineDashedMaterial({
        color: 0x52683f,
        dashSize: 0.4,
        gapSize: 0.3,
        transparent: true,
        opacity: 0.7,
      }),
    );
    this.contact.visible = false;
    this.model.add(this.contact);
    this.omission = new THREE.Line(
      new THREE.BufferGeometry(),
      new THREE.LineDashedMaterial({
        color: 0x83887d,
        dashSize: 0.18,
        gapSize: 0.45,
        transparent: true,
        opacity: 0.45,
      }),
    );
    this.omission.visible = false;
    this.model.add(this.omission);
    this.model.rotation.z = 0.28;
    this.scene.add(this.model);
    const signal = this.abort.signal;
    const on = (id: string, action: () => void) =>
      byId(id).addEventListener("click", action, { signal });
    on("dna-reset", () => {
      this.controls.reset();
      this.requestFrame();
    });
    on("dna-zoom-in", () => this.zoom(1.25));
    on("dna-zoom-out", () => this.zoom(0.8));
    on("dna-expand", () => this.setExpanded(!this.expanded));
    on("dna-play", () => this.togglePlayback());
    on("view-3d", () => this.switchView(true));
    on("view-contacts", () => this.switchView(false));
    on("dna-fallback-map", () => this.switchView(false));
    byId<HTMLInputElement>("dna-progress").addEventListener(
      "input",
      (event) => {
        this.playing = false;
        this.setProgress(
          Number((event.target as HTMLInputElement).value) / 100,
        );
      },
      { signal },
    );
    canvas.addEventListener(
      "pointerdown",
      (event) => {
        this.down = { x: event.clientX, y: event.clientY };
      },
      { signal },
    );
    canvas.addEventListener(
      "pointerup",
      (event) => {
        if (
          Math.hypot(event.clientX - this.down.x, event.clientY - this.down.y) <
          5
        )
          this.pick(event);
      },
      { signal },
    );
    canvas.addEventListener(
      "webglcontextlost",
      (event) => {
        event.preventDefault();
        this.lostContext = true;
        this.playing = false;
        this.showFallback(
          "The graphics context was interrupted. The contact map remains available; reload to retry 3D.",
        );
      },
      { signal },
    );
    this.root.addEventListener("keydown", (event) => this.handleKeys(event), {
      signal,
    });
    document.addEventListener(
      "visibilitychange",
      () => {
        this.lastTime = 0;
        if (!document.hidden) this.requestFrame();
      },
      { signal },
    );
    this.reducedMotion.addEventListener(
      "change",
      () => {
        this.playing = false;
        this.controls.enableDamping = !this.reducedMotion.matches;
        this.updatePlaybackUI();
      },
      { signal },
    );
    this.controls.enableDamping = !this.reducedMotion.matches;
    this.observer = new ResizeObserver(() => this.resize());
    this.observer.observe(this.stage);
    this.visibility = new IntersectionObserver(([entry]) => {
      this.visible = entry.isIntersecting;
      this.lastTime = 0;
      if (this.visible) this.requestFrame();
    });
    this.visibility.observe(this.stage);
    byId("dna-loading").hidden = true;
    this.stage.dataset.state = "ready";
  }

  update(record: GeneRecord, cre: CRE | null) {
    const key = JSON.stringify([
      record.gene.symbol,
      record.gene.chrom,
      record.gene.tss,
      record.locus.tad,
      cre,
      record.illustrative,
    ]);
    if (key === this.selectedKey) return;
    this.selectedKey = key;
    this.record = record;
    this.cre = cre;
    this.locus = cre ? locusWindow(record, cre) : null;
    this.display = this.locus ? helixDisplay(this.locus) : null;
    this.playing = false;
    this.progress = 0;
    byId("dna-region-info").hidden = true;
    if (!this.locus || this.lostContext) {
      this.model.visible = false;
      this.showFallback(
        this.lostContext
          ? "Reload to retry 3D, or use the contact map."
          : "This record does not provide two distinct in-domain positions for a contact preview.",
      );
      return;
    }
    this.model.visible = true;
    this.stage.dataset.state = "ready";
    byId("dna-fallback").hidden = true;
    byId("dna-controls").hidden = false;
    byId("dna-label-gene").textContent = `${record.gene.symbol} · TSS`;
    byId("dna-label-element").textContent = cre!.id;
    const display = this.display!;
    byId("dna-distance").hidden = false;
    byId("dna-separation").textContent =
      `${formatNumber(display.separationBp)} bp apart`;
    byId("dna-coordinate-kind").textContent = record.illustrative
      ? "Placeholder coordinates"
      : "Record coordinates";
    byId("dna-gap-label").hidden = !display.gap;
    byId("dna-gap-label").textContent =
      `⋯ ${formatNumber(display.omittedBp)} bp omitted`;
    byId("dna-resolution").textContent =
      `${formatNumber(display.shownBp)} bp shown${display.gap ? ` · ${formatNumber(display.omittedBp)} bp omitted from this view` : " · no internal omission"}`;
    this.stage.dataset.separationBp = String(display.separationBp);
    this.stage.dataset.omittedBp = String(display.omittedBp);
    this.stage.dataset.element = cre!.id;
    this.renderer.domElement.setAttribute(
      "aria-label",
      `Interactive 3D DNA cutaway: ${record.gene.symbol} transcription start and ${cre!.id} midpoint, ${formatNumber(display.separationBp)} base pairs apart along DNA. ${formatNumber(display.omittedBp)} base pairs omitted. Three-dimensional spacing and sequence are illustrative.`,
    );
    this.setProgress(0);
    this.resize();
  }

  private region(position: number): "gene" | "element" | null {
    if (!this.locus) return null;
    // Position markers only, not a claimed promoter or full element footprint.
    if (Math.abs(position - this.locus.tss) <= 4) return "gene";
    if (Math.abs(position - this.locus.midpoint) <= 4) return "element";
    return null;
  }

  private rebuild() {
    if (!this.display) return;
    const { shape } = this.display;
    const progress = this.progress;
    const { normals, binormals } = displayFrames(this.display, progress, STEPS);
    let sphereCount = 0,
      bondCount = 0;
    const sphere = (
      point: THREE.Vector3,
      radius: number,
      color: THREE.Color,
      region: "gene" | "element" | null,
    ) => {
      this.matrix.makeScale(radius, radius, radius).setPosition(point);
      this.spheres.setMatrixAt(sphereCount, this.matrix);
      this.spheres.setColorAt(sphereCount, color);
      this.sphereRegions[sphereCount++] = region;
    };
    const bond = (
      a: THREE.Vector3,
      b: THREE.Vector3,
      radius: number,
      color: THREE.Color,
    ) => {
      const delta = b.clone().sub(a);
      this.quaternion.setFromUnitVectors(Y_AXIS, delta.clone().normalize());
      this.scale.set(radius, delta.length(), radius);
      this.matrix.compose(
        a.clone().add(b).multiplyScalar(0.5),
        this.quaternion,
        this.scale,
      );
      this.bonds.setMatrixAt(bondCount, this.matrix);
      this.bonds.setColorAt(bondCount++, color);
    };
    for (const segment of helixSamples(this.display)) {
      let previousA: THREE.Vector3 | null = null;
      let previousB: THREE.Vector3 | null = null;
      for (const sample of segment) {
        const { t, position, fade, rung } = sample;
        const i = Math.round(t * STEPS);
        const center = displayCenter(t, progress, this.display);
        const phase = ((position % 10.5) / 10.5) * Math.PI * 2;
        const offset = (angle: number) =>
          center
            .clone()
            .addScaledVector(normals[i], Math.cos(angle) * 1.05)
            .addScaledVector(binormals[i], Math.sin(angle) * 1.05);
        const a = offset(phase),
          b = offset(phase + Math.PI * 0.86);
        const region = this.region(position);
        const faded = (color: THREE.Color) =>
          color.clone().lerp(FADE, 1 - fade);
        const aColor = faded(
          region === "gene" ? GENE : region === "element" ? ELEMENT : STRAND_A,
        );
        const bColor = faded(
          region === "gene" ? GENE : region === "element" ? ELEMENT : STRAND_B,
        );
        const taper = 0.35 + 0.65 * fade;
        sphere(a, (rung ? 0.3 : 0.24) * taper, aColor, region);
        sphere(b, (rung ? 0.3 : 0.24) * taper, bColor, region);
        if (previousA && previousB) {
          bond(previousA, a, 0.19 * taper, aColor);
          bond(previousB, b, 0.19 * taper, bColor);
        }
        previousA = a;
        previousB = b;
        if (rung) {
          let previous = a;
          for (let j = 1; j <= 6; j++) {
            const p = a.clone().lerp(b, j / 7);
            const color =
              region === "gene"
                ? GENE
                : region === "element"
                  ? ELEMENT
                  : j <= 3
                    ? BASE_A
                    : BASE_B;
            sphere(p, 0.2 * taper, faded(color), region);
            bond(previous, p, 0.1 * taper, faded(color));
            previous = p;
          }
          bond(previous, b, 0.1 * taper, bColor);
        }
      }
    }
    this.spheres.count = sphereCount;
    this.bonds.count = bondCount;
    this.spheres.instanceMatrix.needsUpdate =
      this.bonds.instanceMatrix.needsUpdate = true;
    this.spheres.instanceColor!.needsUpdate =
      this.bonds.instanceColor!.needsUpdate = true;
    this.spheres.computeBoundingSphere();
    const a = displayCenter(shape.promoter, progress, this.display);
    const b = displayCenter(shape.element, progress, this.display);
    this.contact.geometry.dispose();
    this.contact.geometry = new THREE.BufferGeometry().setFromPoints([a, b]);
    this.contact.computeLineDistances();
    this.contact.visible = this.progress > 0.8;
    const gap = this.display.gap;
    this.omission.visible = Boolean(gap);
    this.omission.geometry.dispose();
    this.omission.geometry = new THREE.BufferGeometry().setFromPoints(
      gap
        ? Array.from({ length: 33 }, (_, i) =>
            displayCenter(
              gap.start + ((gap.end - gap.start) * i) / 32,
              progress,
              this.display!,
            ),
          )
        : [],
    );
    if (gap) this.omission.computeLineDistances();
  }

  private setProgress(value: number) {
    this.progress = Math.max(0, Math.min(1, value));
    this.rebuild();
    this.updatePlaybackUI();
    this.requestFrame();
  }

  private updatePlaybackUI() {
    const label =
      this.progress < 0.02
        ? "Open model"
        : this.progress > 0.98
          ? "In proximity"
          : "Approaching";
    const slider = byId<HTMLInputElement>("dna-progress");
    slider.value = String(Math.round(this.progress * 100));
    slider.setAttribute("aria-valuetext", label);
    byId("dna-progress-label").textContent = label;
    const button = byId("dna-play");
    const action = this.playing
      ? "Pause"
      : this.progress >= 0.999
        ? "Replay contact"
        : "Preview contact";
    button.setAttribute("aria-label", action);
    button.innerHTML = `${this.playing ? "Ⅱ" : "▷"} <span>${action}</span>`;
  }

  private togglePlayback() {
    if (!this.locus || this.lostContext) return;
    if (this.reducedMotion.matches) {
      this.playing = false;
      this.setProgress(this.progress > 0.5 ? 0 : 1);
      return;
    }
    if (this.playing) this.playing = false;
    else {
      if (this.progress >= 0.999) this.setProgress(0);
      this.playing = true;
      this.lastTime = 0;
      this.requestFrame();
    }
    this.updatePlaybackUI();
  }

  private zoom(factor: number) {
    this.camera.zoom = THREE.MathUtils.clamp(
      this.camera.zoom * factor,
      this.controls.minZoom,
      this.controls.maxZoom,
    );
    this.camera.updateProjectionMatrix();
    this.requestFrame();
  }

  private resize() {
    const width = this.stage.clientWidth,
      height = this.stage.clientHeight;
    if (!width || !height) return;
    const halfHeight = this.display?.gap
      ? Math.max(18, (18 * height) / width)
      : Math.max(22, (28 * height) / width);
    this.camera.left = (-halfHeight * width) / height;
    this.camera.right = (halfHeight * width) / height;
    this.camera.top = halfHeight;
    this.camera.bottom = -halfHeight;
    this.camera.updateProjectionMatrix();
    this.renderer.setSize(width, height, false);
    this.requestFrame();
  }

  private requestFrame() {
    if (
      this.pending ||
      this.disposed ||
      this.lostContext ||
      document.hidden ||
      !this.visible ||
      byId("dna-panel").hidden
    )
      return;
    this.pending = true;
    this.frame = requestAnimationFrame((time) => this.render(time));
  }

  private render(time: number) {
    this.pending = false;
    if (this.disposed || this.lostContext) return;
    if (this.playing) {
      const dt = this.lastTime
        ? Math.min(0.05, (time - this.lastTime) / 1000)
        : 0;
      this.lastTime = time;
      this.progress = Math.min(1, this.progress + dt / 4);
      if (this.progress === 1) this.playing = false;
      this.rebuild();
      this.updatePlaybackUI();
    }
    const changed = this.controls.update();
    this.renderer.render(this.scene, this.camera);
    this.positionLabels();
    if (this.playing || changed) this.requestFrame();
  }

  private positionLabels() {
    if (!this.display || !this.model.visible) return;
    const { shape, gap } = this.display;
    const progress = this.progress;
    const width = this.stage.clientWidth,
      height = this.stage.clientHeight;
    for (const [id, t, side] of [
      ["dna-label-gene", shape.promoter, -1],
      ["dna-label-element", shape.element, 1],
    ] as const) {
      const point = displayCenter(t, progress, this.display);
      this.model.localToWorld(point);
      point.project(this.camera);
      const label = byId(id);
      label.hidden = point.z < -1 || point.z > 1;
      const x = (point.x * 0.5 + 0.5) * width;
      const y = (-point.y * 0.5 + 0.5) * height;
      const half = (label.offsetWidth || 110) / 2;
      label.style.left = `${Math.max(half + 8, Math.min(width - half - 8, x + side * 30))}px`;
      label.style.top = `${Math.max(18, Math.min(height - 26, y + (side === -1 ? -40 : 38)))}px`;
    }
    const label = byId("dna-gap-label");
    label.hidden = !gap;
    if (gap) {
      const point = displayCenter(
        (gap.start + gap.end) / 2,
        progress,
        this.display,
      );
      this.model.localToWorld(point);
      point.project(this.camera);
      const half = label.offsetWidth / 2;
      label.style.left = `${Math.max(half + 8, Math.min(width - half - 8, (point.x * 0.5 + 0.5) * width))}px`;
      label.style.top = `${Math.max(20, Math.min(height - 25, (-point.y * 0.5 + 0.5) * height - 28))}px`;
    }
  }

  private pick(event: PointerEvent) {
    if (!this.locus || !this.record || !this.cre) return;
    const rect = this.renderer.domElement.getBoundingClientRect();
    this.raycaster.setFromCamera(
      new THREE.Vector2(
        ((event.clientX - rect.left) / rect.width) * 2 - 1,
        -((event.clientY - rect.top) / rect.height) * 2 + 1,
      ),
      this.camera,
    );
    const hit = this.raycaster.intersectObject(this.spheres)[0];
    const region =
      hit?.instanceId !== undefined ? this.sphereRegions[hit.instanceId] : null;
    const panel = byId("dna-region-info");
    panel.hidden = !region;
    if (region === "gene")
      panel.textContent = `${this.record.gene.symbol} transcription start · ${regionLabel({ chrom: this.record.gene.chrom, start: this.record.gene.tss, end: this.record.gene.tss })}`;
    if (region === "element")
      panel.textContent = `${this.cre.id} · ${regionLabel({ chrom: this.record.gene.chrom, start: this.cre.start, end: this.cre.end })}`;
  }

  private switchView(threeD: boolean) {
    byId("dna-panel").hidden = !threeD;
    byId("contact-panel").hidden = threeD;
    for (const id of ["dna-reset", "dna-zoom-in", "dna-zoom-out"])
      byId(id).hidden = !threeD;
    byId("view-3d").setAttribute("aria-pressed", String(threeD));
    byId("view-contacts").setAttribute("aria-pressed", String(!threeD));
    if (!threeD) {
      this.playing = false;
      this.updatePlaybackUI();
    } else {
      this.resize();
      this.requestFrame();
    }
  }

  private showFallback(message: string) {
    byId("dna-loading").hidden = true;
    byId("dna-fallback").hidden = false;
    byId("dna-fallback-message").textContent = message;
    byId("dna-controls").hidden = true;
    byId("dna-distance").hidden = true;
    byId("dna-gap-label").hidden = true;
    byId("dna-resolution").textContent = "";
    byId("dna-label-gene").hidden = byId("dna-label-element").hidden = true;
    this.stage.dataset.state = "unavailable";
  }

  private setExpanded(expanded: boolean) {
    this.expanded = expanded;
    this.root.classList.toggle("is-expanded", expanded);
    const button = byId("dna-expand");
    button.setAttribute("aria-pressed", String(expanded));
    button.setAttribute(
      "aria-label",
      expanded ? "Close expanded viewer" : "Expand 3D viewer",
    );
    button.textContent = expanded ? "×" : "⤢";
    if (expanded) {
      this.oldOverflow = document.body.style.overflow;
      document.body.style.overflow = "hidden";
      this.root.setAttribute("role", "dialog");
      this.root.setAttribute("aria-modal", "true");
      this.root.setAttribute("aria-label", "Interactive genome viewer");
      button.focus();
    } else {
      document.body.style.overflow = this.oldOverflow;
      this.root.removeAttribute("role");
      this.root.removeAttribute("aria-modal");
      this.root.removeAttribute("aria-label");
      button.focus();
    }
    this.resize();
  }

  private handleKeys(event: KeyboardEvent) {
    if (!this.expanded) return;
    if (event.key === "Escape") {
      event.preventDefault();
      this.setExpanded(false);
    }
    if (event.key !== "Tab") return;
    const focusable = Array.from(
      this.root.querySelectorAll<HTMLElement>(
        "button, input, canvas[tabindex]",
      ),
    ).filter(
      (e) => e.getClientRects().length && !(e as HTMLButtonElement).disabled,
    );
    const first = focusable[0],
      last = focusable.at(-1);
    if (event.shiftKey && document.activeElement === first) {
      event.preventDefault();
      last?.focus();
    } else if (!event.shiftKey && document.activeElement === last) {
      event.preventDefault();
      first?.focus();
    }
  }

  dispose() {
    this.disposed = true;
    if (this.expanded) document.body.style.overflow = this.oldOverflow;
    cancelAnimationFrame(this.frame);
    this.abort.abort();
    this.observer.disconnect();
    this.visibility.disconnect();
    this.controls.dispose();
    this.spheres.geometry.dispose();
    this.bonds.geometry.dispose();
    (this.spheres.material as THREE.Material).dispose();
    this.contact.geometry.dispose();
    (this.contact.material as THREE.Material).dispose();
    this.omission.geometry.dispose();
    (this.omission.material as THREE.Material).dispose();
    this.renderer.dispose();
    this.renderer.domElement.remove();
  }
}
