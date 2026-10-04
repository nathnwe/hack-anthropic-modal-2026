// Drawing geometry and scroll reveal adapted from feat/site's LoopPanel.
// This is a schematic curve: neither contour length nor 3D conformation is measured.
function initLoopPanel() {
  const wrap = document.querySelector<HTMLElement>(".loop-panel-wrap");
  const film = document.querySelector<HTMLElement>(".loop-panel");
  const canvas = document.querySelector<HTMLCanvasElement>(".loop-canvas");
  const replay = document.querySelector<HTMLButtonElement>(".loop-replay");
  if (!wrap || !film || !canvas || !replay) return;
  const context = canvas.getContext("2d");
  if (!context) return;
  const media = matchMedia("(prefers-reduced-motion: reduce)");
  let width = 0;
  let height = 0;
  let queued = false;
  let replayStart: number | null = null;

  function resize() {
    const dpr = Math.min(devicePixelRatio || 1, 2);
    const rect = canvas!.getBoundingClientRect();
    width = rect.width;
    height = rect.height;
    canvas!.width = Math.round(width * dpr);
    canvas!.height = Math.round(height * dpr);
    context!.setTransform(dpr, 0, 0, dpr, 0, 0);
  }

  function draw(progress: number) {
    if (!width || !height) return;
    const ctx = context!;
    ctx.clearRect(0, 0, width, height);
    const midY = height * 0.78;
    const initialReach = width * 0.3;
    const reach = initialReach * (1 - progress);
    const ax = width * 0.5 - reach;
    const bx = width * 0.5 + reach;
    const lift = height * 0.38 * progress;
    const spread = initialReach * 0.95 * progress;

    function strand(offset: number, colour: string, thickness: number) {
      ctx.beginPath();
      ctx.moveTo(-40, midY + offset);
      ctx.lineTo(ax, midY + offset);
      ctx.bezierCurveTo(
        ax - spread,
        midY - lift + offset,
        bx + spread,
        midY - lift + offset,
        bx,
        midY + offset,
      );
      ctx.lineTo(width + 40, midY + offset);
      ctx.strokeStyle = colour;
      ctx.lineWidth = thickness;
      ctx.lineCap = "round";
      ctx.stroke();
    }
    strand(0, "rgba(111,179,232,0.55)", 1.6);
    strand(9, "rgba(217,138,106,0.4)", 1.2);

    function dot(x: number, y: number, colour: string) {
      ctx.beginPath();
      ctx.arc(x, y, 4.5, 0, Math.PI * 2);
      ctx.fillStyle = colour;
      ctx.fill();
    }
    dot(ax, midY, "rgba(165,76,56,0.95)");
    dot(bx, midY, "rgba(111,179,232,0.95)");
    if (progress > 0.6) {
      const glow = (progress - 0.6) / 0.4;
      const gradient = ctx.createRadialGradient(
        width * 0.5,
        midY,
        0,
        width * 0.5,
        midY,
        150 * glow,
      );
      gradient.addColorStop(0, `rgba(216,183,90,${0.5 * glow})`);
      gradient.addColorStop(1, "rgba(216,183,90,0)");
      ctx.fillStyle = gradient;
      ctx.fillRect(0, midY - 200, width, 400);
    }
    // Text makes the schematic readable without relying on colour alone.
    ctx.font = "12px 'DM Sans', sans-serif";
    ctx.textAlign = "center";
    ctx.fillStyle = "#d9dbdd";
    if (bx - ax < 170) {
      ctx.fillText("Enhancer · promoter", width * 0.5, midY + 34);
    } else {
      ctx.fillText("Enhancer", ax, midY + 34);
      ctx.fillText("Promoter", bx, midY + 34);
    }
  }

  function frame(now: number) {
    queued = false;
    const rect = film!.getBoundingClientRect();
    const scrollProgress = media.matches
      ? 1
      : Math.max(
          0,
          Math.min(1, (innerHeight - rect.top) / (innerHeight * 0.95)),
        );
    const inset =
      Math.max(
        0,
        (innerWidth -
          Math.min(1400, innerWidth - (innerWidth >= 640 ? 80 : 48))) /
          2,
      ) *
      (1 - scrollProgress);
    film!.style.clipPath = `inset(0 ${inset.toFixed(1)}px round ${(16 * (1 - scrollProgress)).toFixed(1)}px)`;
    let progress = scrollProgress;
    if (replayStart !== null && !media.matches) {
      const elapsed = Math.min(1, (now - replayStart) / 2800);
      progress = elapsed * elapsed * (3 - 2 * elapsed);
      if (elapsed < 1) queueFrame();
      else {
        replayStart = null;
        replay!.disabled = false;
      }
    }
    draw(progress);
  }

  function queueFrame() {
    if (queued) return;
    queued = true;
    requestAnimationFrame(frame);
  }
  function stopReplay() {
    replayStart = null;
    replay!.disabled = false;
  }
  replay.hidden = media.matches;
  replay.addEventListener("click", () => {
    if (media.matches) return;
    replayStart = performance.now();
    replay.disabled = true;
    queueFrame();
  });
  addEventListener(
    "scroll",
    () => {
      stopReplay();
      queueFrame();
    },
    { passive: true },
  );
  addEventListener("resize", () => {
    resize();
    queueFrame();
  });
  media.addEventListener("change", () => {
    stopReplay();
    replay.hidden = media.matches;
    queueFrame();
  });
  // Fonts and browser zoom can change element sizes without a window resize.
  new ResizeObserver(() => {
    resize();
    queueFrame();
  }).observe(film);
  resize();
  queueFrame();
}

initLoopPanel();
