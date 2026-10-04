// Preserve the scroll expansion; load the molecular scene only near the panel.
function initLoopPanel() {
  const film = document.querySelector<HTMLElement>(".loop-panel");
  const canvas = document.querySelector<HTMLCanvasElement>(".loop-canvas");
  const toggle = document.querySelector<HTMLButtonElement>(".loop-motion");
  if (!film || !canvas || !toggle) return;
  const media = matchMedia("(prefers-reduced-motion: reduce)");
  const abort = new AbortController();
  let scene: import("../lib/landing-dna").LandingDNA | undefined;
  let loading = false;
  let inView = false;
  let manuallyPaused = false;
  let queued = 0;
  let disposed = false;

  function expand() {
    queued = 0;
    const progress = media.matches
      ? 1
      : Math.max(
          0,
          Math.min(
            1,
            (innerHeight - film!.getBoundingClientRect().top) /
              (innerHeight * 0.95),
          ),
        );
    const inset =
      Math.max(
        0,
        (innerWidth -
          Math.min(1400, innerWidth - (innerWidth >= 640 ? 80 : 48))) /
          2,
      ) *
      (1 - progress);
    film!.style.clipPath = `inset(0 ${inset.toFixed(1)}px round ${(16 * (1 - progress)).toFixed(1)}px)`;
  }
  function queueExpansion() {
    if (!queued) queued = requestAnimationFrame(expand);
  }
  function syncMotion() {
    const paused = manuallyPaused || media.matches;
    toggle!.setAttribute(
      "aria-label",
      paused ? "Play DNA animation" : "Pause DNA animation",
    );
    toggle!.setAttribute(
      "title",
      paused ? "Play DNA animation" : "Pause DNA animation",
    );
    toggle!.dataset.paused = String(paused);
    toggle!.hidden = media.matches || !scene;
    scene?.setMotion(
      inView &&
        !document.hidden &&
        !film!.classList.contains("is-entering") &&
        !paused,
    );
  }
  const observer = new IntersectionObserver(
    async ([entry]) => {
      inView = entry.isIntersecting;
      if (inView && !loading) {
        loading = true;
        try {
          const { LandingDNA } = await import("../lib/landing-dna");
          if (disposed) return;
          scene = new LandingDNA(canvas!, film!, media.matches);
        } catch (error) {
          // Keep the static molecular illustration visible without WebGL.
          console.warn("Landing DNA uses its static illustration:", error);
        }
      }
      syncMotion();
    },
    { threshold: 0 },
  );
  observer.observe(film);
  const filmObserver = new MutationObserver(syncMotion);
  filmObserver.observe(film, { attributes: true, attributeFilter: ["class"] });
  toggle.addEventListener(
    "click",
    () => {
      manuallyPaused = !manuallyPaused;
      syncMotion();
    },
    { signal: abort.signal },
  );
  media.addEventListener(
    "change",
    () => {
      scene?.setReducedMotion(media.matches);
      syncMotion();
      queueExpansion();
    },
    { signal: abort.signal },
  );
  addEventListener("scroll", queueExpansion, {
    passive: true,
    signal: abort.signal,
  });
  addEventListener("resize", queueExpansion, { signal: abort.signal });
  document.addEventListener("visibilitychange", syncMotion, {
    signal: abort.signal,
  });
  addEventListener(
    "pagehide",
    (event) => {
      scene?.setMotion(false);
      if (event.persisted) return;
      disposed = true;
      cancelAnimationFrame(queued);
      observer.disconnect();
      filmObserver.disconnect();
      abort.abort();
      scene?.dispose();
    },
    { signal: abort.signal },
  );
  addEventListener("pageshow", syncMotion, { signal: abort.signal });
  queueExpansion();
}

initLoopPanel();
