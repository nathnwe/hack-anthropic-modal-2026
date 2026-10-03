import * as THREE from "three";
import {
  foldedDNA,
  foldedFrame,
  foldedAnchor,
  hazeCenter,
  type FoldedDNA,
} from "./folded-dna";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";
import {
  dnaCenter,
  dnaFrame,
  locusWindow,
  detailIntervals,
  RISE_NM,
  BP_PER_TURN,
  type LocusWindow,
} from "./dna-geometry";
import {
  formatNumber,
  regionLabel,
  type CRE,
  type GeneRecord,
} from "./records";

const STEPS = 4096;
const GENE = new THREE.Color("#2877a5");
const ELEMENT = new THREE.Color("#ad5641");
const STRAND_A = new THREE.Color("#9baea0");
const STRAND_B = new THREE.Color("#c8b4ab");
const BASE_A = new THREE.Color("#e0e3da");
const BASE_B = new THREE.Color("#c1ccd0");

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
  private haze = new THREE.Group();
  private hazeTexture: THREE.CanvasTexture;
  private cutaway: FoldedDNA | null = null;
  private viewMode: "folded" | "full" = "folded";
  private sphereFade = new THREE.InstancedBufferAttribute(
    new Float32Array(24000).fill(1),
    1,
  );
  private bondFade = new THREE.InstancedBufferAttribute(
    new Float32Array(28000).fill(1),
    1,
  );
  private spheres: THREE.InstancedMesh;
  private bonds: THREE.InstancedMesh;
  private contact: THREE.Line;
  private trace: THREE.Line;
  private geneMarker: THREE.Mesh;
  private elementMarker: THREE.Mesh;
  private focus: "gene" | "element" | "junction" | null = null;
  private geometryDirty = true;
  private overviewTracking = true;
  private viewWidth = 30;
  private viewHeight = 30;
  private sphereRegions: Array<"gene" | "element" | null> = [];
  private frame = 0;
  private pending = false;
  private visible = true;
  private playing = false;
  private progress = 0;
  private lastTime = 0;
  private locus: LocusWindow | null = null;
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
    this.controls.maxZoom = 1_000_000;
    this.controls.listenToKeyEvents(canvas);
    this.controls.update();
    this.controls.saveState();
    this.controls.addEventListener("change", () => {
      this.geometryDirty = true;
      this.requestFrame();
    });
    this.controls.addEventListener("start", () => {
      this.focus = null;
      this.overviewTracking = false;
    });
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
    material.alphaHash = true;
    material.onBeforeCompile = (shader) => {
      shader.vertexShader =
        "attribute float instanceFade; varying float vFade;\n" +
        shader.vertexShader;
      shader.vertexShader = shader.vertexShader.replace(
        "#include <begin_vertex>",
        "#include <begin_vertex>\nvFade = instanceFade;",
      );
      shader.fragmentShader = "varying float vFade;\n" + shader.fragmentShader;
      shader.fragmentShader = shader.fragmentShader.replace(
        "#include <alphahash_fragment>",
        "diffuseColor.a *= vFade;\n#include <alphahash_fragment>",
      );
    };
    this.spheres = new THREE.InstancedMesh(
      new THREE.SphereGeometry(1, 14, 10),
      material,
      24000,
    );
    this.bonds = new THREE.InstancedMesh(
      new THREE.CylinderGeometry(1, 1, 1, 8),
      material,
      28000,
    );
    this.spheres.geometry.setAttribute("instanceFade", this.sphereFade);
    this.bonds.geometry.setAttribute("instanceFade", this.bondFade);
    this.sphereFade.setUsage(THREE.DynamicDrawUsage);
    this.bondFade.setUsage(THREE.DynamicDrawUsage);
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
    this.trace = new THREE.Line(
      new THREE.BufferGeometry(),
      new THREE.LineBasicMaterial({
        vertexColors: true,
        transparent: true,
        opacity: 0.85,
      }),
    );
    this.model.add(this.trace);
    this.geneMarker = new THREE.Mesh(
      new THREE.SphereGeometry(1, 16, 12),
      new THREE.MeshBasicMaterial({ color: GENE }),
    );
    this.elementMarker = new THREE.Mesh(
      new THREE.SphereGeometry(1, 16, 12),
      new THREE.MeshBasicMaterial({ color: ELEMENT }),
    );
    this.model.add(this.geneMarker, this.elementMarker);
    this.scene.add(this.model);
    const smoke = document.createElement("canvas");
    smoke.width = smoke.height = 192;
    const context = smoke.getContext("2d")!;
    // Soft overlapping lobes give a quiet, smoky boundary rather than a hard mask.
    for (let i = 0; i < 13; i++) {
      const angle = i * 2.39996,
        r = 12 + (18 * (i % 3)) / 2;
      const x = 96 + Math.cos(angle) * r,
        y = 96 + Math.sin(angle) * r * 0.7;
      const gradient = context.createRadialGradient(
        x,
        y,
        0,
        x,
        y,
        55 + (i % 4) * 5,
      );
      gradient.addColorStop(0, "rgba(128, 141, 119, 0.10)");
      gradient.addColorStop(0.45, "rgba(145, 154, 136, 0.045)");
      gradient.addColorStop(1, "rgba(160, 167, 150, 0)");
      context.fillStyle = gradient;
      context.fillRect(0, 0, 192, 192);
    }
    this.hazeTexture = new THREE.CanvasTexture(smoke);
    for (let i = 0; i < 3; i++) {
      const sprite = new THREE.Sprite(
        new THREE.SpriteMaterial({
          map: this.hazeTexture,
          transparent: true,
          depthWrite: false,
          depthTest: false,
          opacity: 0.55,
        }),
      );
      sprite.renderOrder = 3 + i;
      sprite.userData.lobe = i;
      this.haze.add(sprite);
    }
    this.scene.add(this.haze);
    const signal = this.abort.signal;
    const on = (id: string, action: () => void) =>
      byId(id).addEventListener("click", action, { signal });
    on("dna-reset", () => {
      this.focusView("folded");
      this.requestFrame();
    });
    for (const view of [
      "folded",
      "full",
      "gene",
      "element",
      "junction",
    ] as const)
      on(`dna-focus-${view}`, () => this.focusView(view));
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
    this.cutaway = this.locus ? foldedDNA(this.locus) : null;
    this.viewMode = this.cutaway ? "folded" : "full";
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
    byId("dna-focus-gene").textContent = `${record.gene.symbol} / TSS`;
    byId("dna-label-element").textContent = cre!.id;
    const locus = this.locus!;
    byId("dna-distance").hidden = false;
    const separation = Math.abs(locus.midpoint - locus.tss);
    byId("dna-separation").textContent = `${formatNumber(separation)} bp apart`;
    byId("dna-coordinate-kind").textContent = record.illustrative
      ? "Placeholder coordinates"
      : record.caseStudy
        ? `${record.locus.cell_type} · ${record.caseStudy.assembly}`
        : "Record coordinates";
    this.stage.dataset.separationBp = String(separation);
    this.stage.dataset.element = cre!.id;
    byId<HTMLButtonElement>("dna-play").disabled = locus.overlap;
    byId<HTMLInputElement>("dna-progress").disabled = locus.overlap;
    byId<HTMLButtonElement>("dna-focus-junction").disabled = locus.overlap;
    this.setProgress(0);
    this.focusView("folded");
    this.resize();
  }

  private region(position: number): "gene" | "element" | null {
    if (!this.locus || !this.record) return null;
    if (this.viewMode === "folded") {
      // Local position patches, not a claim that a full gene fits in this window.
      if (Math.abs(position - this.locus.tss) <= 10) return "gene";
      if (Math.abs(position - this.locus.midpoint) <= 10) return "element";
      return null;
    }
    if (Math.abs(position - this.locus.tss) < 4) return "gene";
    if (position >= this.locus.elementStart && position < this.locus.elementEnd)
      return "element";
    const ref = this.record.caseStudy?.reference;
    if (ref && position >= ref.start && position < ref.end) return "gene";
    return null;
  }

  private center(position: number) {
    return this.viewMode === "folded" && this.cutaway
      ? foldedAnchor(position, this.progress, this.cutaway)
      : dnaCenter(position, this.progress, this.locus!);
  }

  private updateScaleDescription() {
    if (!this.locus || !this.record || !this.cre) return;
    const folded = this.viewMode === "folded" && this.cutaway;
    const omitted = folded ? folded.omittedBp : 0;
    const shown = folded ? folded.shownBp : this.locus.end - this.locus.start;
    this.root.dataset.view = this.viewMode;
    this.stage.dataset.view = this.viewMode;
    this.stage.dataset.omittedBp = String(omitted);
    this.stage.dataset.totalBp = String(
      folded ? folded.end - folded.start : shown,
    );
    byId("dna-resolution").textContent = folded
      ? `${formatNumber(shown)} bp in two local windows · gap compressed, not to scale`
      : `${formatNumber(shown)} bp in the full span · 0 bp omitted`;
    byId("dna-fold-note").textContent = folded
      ? "Illustrative folding, not a measured TAD or atomic structure. Coloured patches locate the reference TSS and element midpoint."
      : "Full DNA contour at one scale; folding remains illustrative, not a measured TAD or atomic structure.";
    byId("dna-omitted-count").textContent = `${formatNumber(omitted)} bp`;
    byId("dna-omission").hidden = !folded;
    this.renderer.domElement.setAttribute(
      "aria-label",
      `${this.record.gene.symbol} and ${this.cre.id}: ${folded ? `folded molecular close-up, ${formatNumber(omitted)} base pairs omitted in the haze` : "full continuous genomic span"}. Three-dimensional folding is illustrative. Drag to rotate, scroll to zoom.`,
    );
    for (const view of ["folded", "full"])
      byId(`dna-focus-${view}`).setAttribute(
        "aria-pressed",
        String(view === this.viewMode),
      );
  }

  private rebuildFolded() {
    if (!this.cutaway || !this.locus) return;
    this.trace.visible =
      this.geneMarker.visible =
      this.elementMarker.visible =
        false;
    this.haze.visible = true;
    const center = hazeCenter(this.cutaway, this.progress);
    this.haze.position.copy(center);
    this.haze.children.forEach((child, i) => {
      const width = this.cutaway!.hazeWidth;
      child.position.set(
        (i - 1) * width * 0.18,
        Math.sin(i * 2) * width * 0.1,
        0,
      );
      child.scale.set(width * (1.65 - i * 0.12), width * (1.2 + i * 0.06), 1);
      (child as THREE.Sprite).material.opacity = 0.6;
    });
    const a = this.center(this.locus.tss),
      b = this.center(this.locus.midpoint);
    this.contact.geometry.dispose();
    this.contact.geometry = new THREE.BufferGeometry().setFromPoints([a, b]);
    this.contact.computeLineDistances();
    this.contact.visible = this.progress > 0.96;
    this.stage.dataset.anchorSeparationNm = String(a.distanceTo(b));
    // This mode intentionally has no whole-locus physical contour measurement.
    delete this.stage.dataset.contourNm;
    if (this.focus) this.moveFocus(this.focus);
    this.geometryDirty = true;
  }

  private rebuild() {
    if (!this.locus) return;
    if (this.viewMode === "folded" && this.cutaway) {
      this.rebuildFolded();
      return;
    }
    this.haze.visible = false;
    const locus = this.locus;
    const positions = Array.from(
      { length: STEPS + 1 },
      (_, i) => locus.start + ((locus.end - locus.start) * i) / STEPS,
    );
    // Include annotation boundaries exactly, even when narrower than an overview segment.
    positions.push(
      locus.tss,
      locus.midpoint,
      locus.elementStart,
      locus.elementEnd,
    );
    const ref = this.record?.caseStudy?.reference;
    if (ref)
      positions.push(
        Math.max(locus.start, ref.start),
        Math.min(locus.end, ref.end),
      );
    positions.sort((a, b) => a - b);
    const colors = positions.flatMap((p) =>
      (this.region(p) === "gene"
        ? GENE
        : this.region(p) === "element"
          ? ELEMENT
          : STRAND_A
      ).toArray(),
    );
    this.trace.geometry.dispose();
    const points = positions.map((p) => dnaCenter(p, this.progress, locus));
    this.trace.geometry = new THREE.BufferGeometry().setFromPoints(points);
    const box = new THREE.Box3().setFromPoints(points);
    if (this.overviewTracking) {
      const target = box.getCenter(new THREE.Vector3());
      this.camera.position.add(target.clone().sub(this.controls.target));
      this.controls.target.copy(target);
    }
    const size = box.getSize(new THREE.Vector3());
    const scale = Math.min(126 / Math.max(size.x, 1), 66 / Math.max(size.y, 1));
    const xy = (v: THREE.Vector3) => [
      80 + (v.x - (box.min.x + box.max.x) / 2) * scale,
      42 - (v.y - (box.min.y + box.max.y) / 2) * scale,
    ];
    const path = points
      .filter((_, i) => i % 32 === 0 || i === points.length - 1)
      .map((v, i) => `${i ? "L" : "M"}${xy(v).join(" ")}`)
      .join(" ");
    const gene = xy(dnaCenter(locus.tss, this.progress, locus)),
      element = xy(dnaCenter(locus.midpoint, this.progress, locus));
    byId("dna-minimap").innerHTML =
      `<svg viewBox="0 0 160 88" role="img" aria-label="Full continuous DNA span"><path d="${path}" fill="none" stroke="#748C63" stroke-width="1.4"/><circle cx="${gene[0]}" cy="${gene[1]}" r="3" fill="#2877a5"/><circle cx="${element[0]}" cy="${element[1]}" r="3" fill="#ad5641"/></svg><span>Full DNA span</span>`;
    this.trace.geometry.setAttribute(
      "color",
      new THREE.Float32BufferAttribute(colors, 3),
    );
    const a = dnaCenter(locus.tss, this.progress, locus),
      b = dnaCenter(locus.midpoint, this.progress, locus);
    this.geneMarker.position.copy(a);
    this.elementMarker.position.copy(b);
    this.contact.geometry.dispose();
    this.contact.geometry = new THREE.BufferGeometry().setFromPoints([a, b]);
    this.contact.computeLineDistances();
    this.contact.visible = this.progress > 0.999;
    this.stage.dataset.contourNm = String((locus.end - locus.start) * RISE_NM);
    this.stage.dataset.anchorSeparationNm = String(a.distanceTo(b));
    this.geometryDirty = true;
    if (this.focus) this.moveFocus(this.focus);
  }

  private detailGeometry() {
    if (!this.locus || !this.stage.clientWidth || !this.stage.clientHeight)
      return;
    const height = (this.camera.top - this.camera.bottom) / this.camera.zoom;
    const width =
      (height * this.stage.clientWidth) / Math.max(1, this.stage.clientHeight);
    const folded = this.viewMode === "folded" && this.cutaway;
    const ranges = folded
      ? folded.arms
      : height < 130
        ? detailIntervals(
            this.locus,
            this.progress,
            this.controls.target,
            Math.hypot(width, height) * 0.6,
          )
        : [];
    const bpCount = ranges.reduce((n, r) => n + r.end - r.start, 0);
    const detailed = bpCount > 0 && bpCount <= 2100;
    this.spheres.visible = this.bonds.visible = detailed;
    // A continuous centreline remains at all scales, including beyond the viewport.
    this.trace.visible = !folded;
    (this.trace.material as THREE.LineBasicMaterial).opacity = detailed
      ? 0.12
      : 0.9;
    this.geneMarker.visible = this.elementMarker.visible = !folded && !detailed;
    const markerSize = (height * 5) / Math.max(1, this.stage.clientHeight);
    this.geneMarker.scale.setScalar(markerSize);
    this.elementMarker.scale.setScalar(markerSize);
    byId("dna-minimap").hidden = !!folded || !detailed;
    this.stage.dataset.resolution = detailed ? "base-pair" : "centreline";
    this.stage.dataset.detailBp = detailed ? String(bpCount) : "0";
    byId("dna-scale-mode").textContent = folded
      ? "Molecular close-up · one rung per local bp"
      : detailed
        ? "Molecular detail · one rung per bp"
        : "Continuous DNA · zoom in to resolve base pairs";
    byId("dna-scale-line").hidden = byId("dna-scale-label").hidden = !!folded;
    const raw = (width * 100) / Math.max(1, this.stage.clientWidth);
    const power = 10 ** Math.floor(Math.log10(raw));
    const nm =
      [1, 2, 5, 10]
        .map((n) => n * power)
        .filter((n) => n <= raw)
        .at(-1) || power;
    byId("dna-scale-line").style.width =
      `${(nm / width) * this.stage.clientWidth}px`;
    byId("dna-scale-label").textContent =
      nm >= 1000 ? `${formatNumber(nm / 1000)} µm` : `${formatNumber(nm)} nm`;
    if (!detailed) return;
    let sphereCount = 0,
      bondCount = 0;
    const sphere = (
      point: THREE.Vector3,
      radius: number,
      color: THREE.Color,
      region: "gene" | "element" | null,
      fade = 1,
    ) => {
      this.matrix.makeScale(radius, radius, radius).setPosition(point);
      this.spheres.setMatrixAt(sphereCount, this.matrix);
      this.spheres.setColorAt(sphereCount, color);
      this.sphereFade.setX(sphereCount, fade);
      this.sphereRegions[sphereCount++] = region;
    };
    const bond = (
      a: THREE.Vector3,
      b: THREE.Vector3,
      radius: number,
      color: THREE.Color,
      fade = 1,
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
      this.bondFade.setX(bondCount, fade);
      this.bonds.setColorAt(bondCount++, color);
    };
    for (const range of ranges) {
      let previousA: THREE.Vector3 | null = null,
        previousB: THREE.Vector3 | null = null;
      for (
        let position = range.start + 0.5;
        position < range.end;
        position += 0.5
      ) {
        const { center, normal, binormal, fade } = folded
          ? foldedFrame(position, this.progress, folded)
          : { ...dnaFrame(position, this.progress, this.locus), fade: 1 };
        const phase = ((position % BP_PER_TURN) / BP_PER_TURN) * Math.PI * 2;
        const offset = (angle: number) =>
          center
            .clone()
            .addScaledVector(normal, Math.cos(angle))
            .addScaledVector(binormal, Math.sin(angle));
        const a = offset(phase),
          b = offset(phase + Math.PI * 0.86);
        const region = this.region(position),
          rung = position % 1 === 0.5;
        const aColor =
          region === "gene" ? GENE : region === "element" ? ELEMENT : STRAND_A;
        const bColor =
          region === "gene" ? GENE : region === "element" ? ELEMENT : STRAND_B;
        sphere(a, rung ? 0.3 : 0.24, aColor, region, fade);
        sphere(b, rung ? 0.3 : 0.24, bColor, region, fade);
        if (previousA && previousB) {
          bond(previousA, a, 0.17, aColor, fade);
          bond(previousB, b, 0.17, bColor, fade);
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
            sphere(p, 0.2, color, region, fade);
            bond(previous, p, 0.1, color, fade);
            previous = p;
          }
          bond(previous, b, 0.1, bColor, fade);
        }
      }
    }
    this.sphereFade.needsUpdate = this.bondFade.needsUpdate = true;
    this.spheres.count = sphereCount;
    this.bonds.count = bondCount;
    this.spheres.instanceMatrix.needsUpdate =
      this.bonds.instanceMatrix.needsUpdate = true;
    this.spheres.instanceColor!.needsUpdate =
      this.bonds.instanceColor!.needsUpdate = true;
    this.spheres.computeBoundingSphere();
  }

  private moveFocus(view: "gene" | "element" | "junction") {
    if (!this.locus) return;
    const gene = this.center(this.locus.tss),
      element = this.center(this.locus.midpoint);
    const target =
      view === "gene"
        ? gene
        : view === "element"
          ? element
          : gene.add(element).multiplyScalar(0.5);
    this.camera.position.add(target.clone().sub(this.controls.target));
    this.controls.target.copy(target);
  }

  private focusView(view: "folded" | "full" | "gene" | "element" | "junction") {
    if (!this.locus) return;
    this.viewMode = view !== "full" && this.cutaway ? "folded" : "full";
    this.updateScaleDescription();
    if (this.viewMode === "folded" && this.cutaway) {
      this.focusFolded(view);
      return;
    }
    this.playing = false;
    this.focus = view === "full" || view === "folded" ? null : view;
    this.overviewTracking = view === "full";
    const length = (this.locus.end - this.locus.start) * RISE_NM;
    this.camera.near = 0.01;
    this.camera.far = Math.max(1000, length * 8);
    const bounds = new THREE.Box3();
    for (const p of [0, 0.5, 1])
      for (let i = 0; i <= 128; i++)
        bounds.expandByPoint(
          dnaCenter(
            this.locus.start + ((this.locus.end - this.locus.start) * i) / 128,
            p,
            this.locus,
          ),
        );
    const center = bounds.getCenter(new THREE.Vector3()),
      size = bounds.getSize(new THREE.Vector3());
    this.viewWidth = Math.max(20, size.x * 0.63);
    this.viewHeight = Math.max(20, size.y * 0.63);
    this.controls.target.copy(center);
    this.camera.position
      .copy(center)
      .add(new THREE.Vector3(0, 0, Math.max(300, length * 3)));
    this.camera.zoom = 1;
    this.resize();
    if (view === "junction") this.setProgress(1);
    if (view !== "full" && view !== "folded") {
      this.moveFocus(view);
      this.camera.zoom = (this.camera.top - this.camera.bottom) / 18;
    }
    this.rebuild();
    this.controls.update();
    this.updatePlaybackUI();
    this.resize();
    this.geometryDirty = true;
    this.requestFrame();
  }

  private focusFolded(
    view: "folded" | "full" | "gene" | "element" | "junction",
  ) {
    if (!this.cutaway) return;
    this.playing = false;
    this.overviewTracking = false;
    this.focus =
      view === "gene" || view === "element" || view === "junction"
        ? view
        : null;
    const bounds = new THREE.Box3();
    for (const progress of [0, 1])
      for (const arm of this.cutaway.arms)
        for (let i = 0; i <= 64; i++)
          bounds.expandByPoint(
            foldedAnchor(
              arm.start + ((arm.end - arm.start) * i) / 64,
              progress,
              this.cutaway,
            ),
          );
    const center = bounds.getCenter(new THREE.Vector3()),
      size = bounds.getSize(new THREE.Vector3());
    // Keep the same molecular magnification across the three source examples.
    // Their compressed gaps can then look different without shrinking the DNA.
    this.viewWidth = Math.max(26, size.x * 0.56 + 1.2);
    this.viewHeight = Math.max(12, size.y * 0.65 + 2);
    this.camera.near = 0.01;
    this.camera.far = 500;
    this.camera.zoom = 1;
    this.controls.target.copy(center);
    this.camera.position.copy(center).add(new THREE.Vector3(0, 2, 110));
    this.resize();
    if (view === "junction") this.setProgress(1);
    if (this.focus) {
      this.moveFocus(this.focus);
      this.camera.zoom = (this.camera.top - this.camera.bottom) / 18;
    }
    this.rebuild();
    this.controls.update();
    this.updatePlaybackUI();
    this.resize();
    this.requestFrame();
  }

  private setProgress(value: number) {
    this.progress = Math.max(0, Math.min(1, value));
    this.rebuild();
    this.updatePlaybackUI();
    this.requestFrame();
  }

  private updatePlaybackUI() {
    const label = this.locus?.overlap
      ? "Overlaps TSS"
      : this.progress < 0.02
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
    const halfHeight = Math.max(
      this.viewHeight,
      (this.viewWidth * height) / width,
    );
    this.camera.left = (-halfHeight * width) / height;
    this.camera.right = (halfHeight * width) / height;
    this.camera.top = halfHeight;
    this.camera.bottom = -halfHeight;
    this.camera.updateProjectionMatrix();
    this.renderer.setSize(width, height, false);
    this.geometryDirty = true;
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
    if (this.geometryDirty) {
      this.detailGeometry();
      this.geometryDirty = false;
    }
    this.renderer.render(this.scene, this.camera);
    this.positionLabels();
    if (this.playing || changed) this.requestFrame();
  }

  private positionLabels() {
    if (!this.locus || !this.model.visible) return;
    const width = this.stage.clientWidth,
      height = this.stage.clientHeight;
    for (const [id, position, side] of [
      ["dna-label-gene", this.locus.tss, -1],
      ["dna-label-element", this.locus.midpoint, 1],
    ] as const) {
      const point = this.center(position).project(this.camera);
      const label = byId(id);
      label.hidden =
        Math.abs(point.x) > 1.1 ||
        Math.abs(point.y) > 1.1 ||
        Math.abs(point.z) > 1;
      const x = (point.x * 0.5 + 0.5) * width,
        y = (-point.y * 0.5 + 0.5) * height;
      const half = (label.offsetWidth || 100) / 2;
      label.style.left = `${Math.max(half + 8, Math.min(width - half - 8, x + side * (half + 10)))}px`;
      label.style.top = `${Math.max(22, Math.min(height - 65, y + side * 34))}px`;
    }
    const label = byId("dna-omission");
    if (this.viewMode === "folded" && this.cutaway) {
      const point = hazeCenter(this.cutaway, this.progress).project(
        this.camera,
      );
      label.hidden =
        Math.abs(point.x) > 1.1 ||
        Math.abs(point.y) > 1.1 ||
        Math.abs(point.z) > 1;
      const half = (label.offsetWidth || 150) / 2;
      label.style.left = `${Math.max(half + 6, Math.min(width - half - 6, (point.x * 0.5 + 0.5) * width))}px`;
      label.style.top = `${Math.max(35, Math.min(height - 70, (-point.y * 0.5 + 0.5) * height))}px`;
    } else label.hidden = true;
  }

  private pick(event: PointerEvent) {
    if (!this.locus || !this.record || !this.cre || !this.spheres.visible)
      return;
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
    byId("dna-resolution").textContent = "";
    byId("dna-omission").hidden = true;
    this.haze.visible = false;
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
    this.trace.geometry.dispose();
    (this.trace.material as THREE.Material).dispose();
    for (const marker of [this.geneMarker, this.elementMarker]) {
      marker.geometry.dispose();
      (marker.material as THREE.Material).dispose();
    }
    this.hazeTexture.dispose();
    this.haze.children.forEach((child) =>
      (child as THREE.Sprite).material.dispose(),
    );
    this.renderer.dispose();
    this.renderer.domElement.remove();
  }
}
