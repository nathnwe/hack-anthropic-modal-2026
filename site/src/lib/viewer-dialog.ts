// Shared expansion works even when the 3D renderer is unavailable or not loaded.
export function viewerDialog() {
  const root = document.getElementById("genome-viewer")!;
  const button = document.getElementById("dna-expand")!;
  const controller = new AbortController();
  let expanded = false,
    oldOverflow = "";
  function setExpanded(value: boolean) {
    expanded = value;
    root.classList.toggle("is-expanded", value);
    button.setAttribute("aria-pressed", String(value));
    button.setAttribute(
      "aria-label",
      value ? "Close expanded viewer" : "Expand viewer",
    );
    button.textContent = value ? "×" : "⤢";
    if (value) {
      oldOverflow = document.body.style.overflow;
      document.body.style.overflow = "hidden";
      root.setAttribute("role", "dialog");
      root.setAttribute("aria-modal", "true");
      root.setAttribute("aria-label", "Interactive genome viewer");
    } else {
      document.body.style.overflow = oldOverflow;
      root.removeAttribute("role");
      root.removeAttribute("aria-modal");
      root.removeAttribute("aria-label");
    }
    button.focus();
  }
  button.addEventListener("click", () => setExpanded(!expanded), {
    signal: controller.signal,
  });
  root.addEventListener(
    "keydown",
    (event) => {
      if (!expanded) return;
      if (event.key === "Escape") {
        event.preventDefault();
        setExpanded(false);
      }
      if (event.key !== "Tab") return;
      const focusable = Array.from(
        root.querySelectorAll<HTMLElement>(
          'button, input, select, a[href], [tabindex="0"]',
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
    },
    { signal: controller.signal },
  );
  return () => {
    if (expanded) setExpanded(false);
    controller.abort();
  };
}
