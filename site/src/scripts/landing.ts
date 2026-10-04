const story = document.querySelector<HTMLElement>(".story");
const thread = document.querySelector<SVGPathElement>(".thread-progress");
const media = matchMedia("(prefers-reduced-motion: reduce)");
const pending = new Set(
  document.querySelectorAll<HTMLElement>(".story .landing-reveal"),
);
let scheduled = false;

function draw() {
  scheduled = false;
  // Reveal past anchors too: a direct #team visit must not leave earlier copy
  // hidden. Content is visible by default if JavaScript is unavailable.
  for (const element of pending) {
    if (
      media.matches ||
      element.getBoundingClientRect().top < innerHeight * 0.92
    ) {
      element.classList.remove("reveal-pending");
      pending.delete(element);
    }
  }
  if (!story || !thread) return;
  const rect = story.getBoundingClientRect();
  const progress = media.matches
    ? 1
    : Math.max(0, Math.min(1, (innerHeight * 0.65 - rect.top) / rect.height));
  thread.style.strokeDashoffset = String(1 - progress);
}

function queueDraw() {
  if (scheduled) return;
  scheduled = true;
  requestAnimationFrame(draw);
}

if (!media.matches) {
  for (const element of pending) element.classList.add("reveal-pending");
  // Establish the starting style before removing it on the next frame.
  // This mirrors the 26px / 1.1s fade-and-rise used by the original landing.
  story?.getBoundingClientRect();
}

addEventListener("scroll", queueDraw, { passive: true });
addEventListener("resize", queueDraw);
media.addEventListener("change", queueDraw);
// A keyboard user should never focus an invisible link between scroll frames.
document.addEventListener("focusin", draw);
queueDraw();
