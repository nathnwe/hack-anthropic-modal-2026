function initLandingVideo() {
  const launch = document.querySelector<HTMLAnchorElement>(".loop-launch");
  const panel = document.querySelector<HTMLElement>(".loop-panel");
  const dialog = document.querySelector<HTMLDialogElement>(".cell-film-dialog");
  const video = document.querySelector<HTMLVideoElement>(".cell-film-video");
  const closeButton =
    document.querySelector<HTMLButtonElement>(".cell-film-close");
  const playButton =
    document.querySelector<HTMLButtonElement>(".cell-film-play");
  const loading = document.querySelector<HTMLElement>(".cell-film-loading");
  const error = document.querySelector<HTMLElement>(".cell-film-error");
  if (
    !launch ||
    !panel ||
    !dialog ||
    !video ||
    !closeButton ||
    !playButton ||
    !loading ||
    !error
  )
    return;
  // The link remains a direct-video fallback if modal dialogs are unsupported.
  if (typeof dialog.showModal !== "function") return;

  const motion = matchMedia("(prefers-reduced-motion: reduce)");
  let closing = false;
  let closeTimer: number | undefined;
  let previousOverflow = "";
  let generation = 0;

  const active = () => dialog.open && !closing;
  const revealFrame = () => {
    if (!active()) return;
    dialog.classList.add("is-ready");
    loading.hidden = true;
  };

  const tryPlay = () => {
    const attempt = generation;
    // Keep play() within the click gesture, before waiting for the fade, so
    // Safari and mobile browsers can honour the user's playback request.
    video.play().catch(() => {
      if (!active() || attempt !== generation || video.error) return;
      loading.hidden = true;
      playButton.hidden = false;
      if (video.readyState >= HTMLMediaElement.HAVE_CURRENT_DATA) revealFrame();
    });
  };

  const finishClose = () => {
    window.clearTimeout(closeTimer);
    closeTimer = undefined;
    video.pause();
    if (dialog.open) dialog.close();
  };

  const close = () => {
    if (!active()) return;
    closing = true;
    generation++;
    video.pause();
    dialog.classList.remove("is-open");
    panel.classList.remove("is-entering");
    if (motion.matches) finishClose();
    else closeTimer = window.setTimeout(finishClose, 550);
  };

  launch.addEventListener("click", (event) => {
    // Preserve open-in-new-tab and other normal link behaviour.
    if (event.metaKey || event.ctrlKey || event.shiftKey || event.altKey)
      return;
    event.preventDefault();
    if (dialog.open) return;
    generation++;
    closing = false;
    previousOverflow = document.documentElement.style.overflow;
    document.documentElement.style.overflow = "hidden";
    dialog.classList.remove("is-open", "is-ready");
    loading.hidden = false;
    error.hidden = true;
    playButton.hidden = true;
    dialog.showModal();
    closeButton.focus({ preventScroll: true });
    panel.classList.add("is-entering");
    // Establish the transparent dialog before starting the crossfade.
    dialog.getBoundingClientRect();
    requestAnimationFrame(() => {
      if (active()) dialog.classList.add("is-open");
    });

    // No video bytes are requested by the landing page until this point.
    if (!video.getAttribute("src") || video.error) {
      video.src = video.dataset.src!;
      video.load();
    } else {
      video.currentTime = 0;
      if (video.readyState >= HTMLMediaElement.HAVE_CURRENT_DATA) revealFrame();
    }
    tryPlay();
  });

  closeButton.addEventListener("click", close);
  dialog.addEventListener("cancel", (event) => {
    event.preventDefault();
    close();
  });
  dialog.addEventListener("close", () => {
    generation++;
    video.pause();
    window.clearTimeout(closeTimer);
    closeTimer = undefined;
    closing = false;
    dialog.classList.remove("is-open", "is-ready");
    panel.classList.remove("is-entering");
    document.documentElement.style.overflow = previousOverflow;
    launch.focus({ preventScroll: true });
  });
  playButton.addEventListener("click", () => {
    playButton.hidden = true;
    tryPlay();
  });
  video.addEventListener("loadeddata", revealFrame);
  video.addEventListener("playing", () => {
    revealFrame();
    playButton.hidden = true;
  });
  video.addEventListener("error", () => {
    if (!active()) return;
    loading.hidden = true;
    playButton.hidden = true;
    error.hidden = false;
  });
  // If the preference changes during the closing fade, release the modal now.
  motion.addEventListener("change", () => {
    if (motion.matches && closing) finishClose();
  });
  addEventListener("pagehide", () => {
    video.pause();
    if (dialog.open) finishClose();
  });
}

initLandingVideo();
