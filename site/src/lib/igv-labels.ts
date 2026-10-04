/** IGV 3.8 uses an open shadow root and draws its default names over the data.
 * Keep native labels hidden and mirror their row geometry in an external gutter.
 * Measuring the rendered viewports keeps names aligned after zoom/reflow/loading.
 */
export function trackNameGutter(
  host: HTMLElement,
  gutter: HTMLElement,
  metadata: (name: string) => { color?: string; detail?: string; url?: string },
) {
  const root = host.shadowRoot!;
  const style = document.createElement("style");
  style.textContent = `
    .igv-container { font-family: 'DM Sans', sans-serif; border: 0; background: #fffefb; }
    .igv-track-label { display: none !important; }
    .igv-viewport { border-bottom: 1px solid #efeee6; }
    .igv-trackgear-container, .igv-track-drag-handle { opacity: .55; }
  `;
  root.append(style);
  let frame = 0;
  const observed = new Set<Element>();
  const labels = new Map<Element, HTMLElement>();
  const resize = new ResizeObserver(schedule);
  function render() {
    frame = 0;
    if (!host.isConnected || !host.offsetWidth) return;
    const origin = gutter.getBoundingClientRect().top;
    const rows = [
      ...root.querySelectorAll<HTMLElement>(".igv-column .igv-viewport"),
    ];
    for (const [row, label] of labels)
      if (!rows.includes(row as HTMLElement)) {
        label.remove();
        labels.delete(row);
        resize.unobserve(row);
        observed.delete(row);
      }
    for (const row of rows) {
      const name = row.querySelector(".igv-track-label")?.textContent?.trim();
      if (!name) continue;
      if (!observed.has(row)) {
        resize.observe(row);
        observed.add(row);
      }
      let label = labels.get(row);
      const meta = metadata(name);
      if (!label) {
        label = document.createElement("div");
        label.className = "igv-track-name";
        gutter.append(label);
        labels.set(row, label);
      }
      if (label.dataset.name !== name) {
        label.replaceChildren();
        label.dataset.name = name;
        const title = document.createElement(meta.url ? "a" : "span");
        title.textContent = name;
        if (title instanceof HTMLAnchorElement && meta.url) {
          title.href = meta.url;
          title.target = "_blank";
          title.rel = "noopener noreferrer";
        }
        label.append(title);
        if (meta.detail) {
          const detail = document.createElement("small");
          detail.textContent = meta.detail;
          label.append(detail);
        }
      }
      label.style.setProperty("--track-color", meta.color || "#626657");
      label.classList.toggle("is-anchor", name.startsWith("Anchor "));
      const box = row.getBoundingClientRect();
      label.style.top = `${box.top - origin}px`;
      label.style.height = `${box.height}px`;
    }
  }
  function schedule() {
    if (!frame) frame = requestAnimationFrame(render);
  }
  const mutation = new MutationObserver(schedule);
  mutation.observe(root, {
    childList: true,
    subtree: true,
    attributes: true,
    attributeFilter: ["style"],
  });
  resize.observe(host);
  schedule();
  return () => {
    mutation.disconnect();
    resize.disconnect();
    cancelAnimationFrame(frame);
    style.remove();
    gutter.replaceChildren();
  };
}
