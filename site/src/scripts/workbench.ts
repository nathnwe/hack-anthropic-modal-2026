import { suggestGenes, type GeneSuggestion } from "../lib/live-genes";
import { escapeHTML as esc } from "../lib/records";
const form = document.querySelector<HTMLFormElement>("#query-form")!;
const input = document.querySelector<HTMLInputElement>("#gene")!;
const source = document.querySelector<HTMLInputElement>("#gene-source")!;
const list = document.querySelector<HTMLElement>("#gene-suggestions")!;
const status = document.querySelector<HTMLElement>("#gene-help")!;
let items: GeneSuggestion[] = [],
  active = -1,
  timer: ReturnType<typeof setTimeout>,
  controller: AbortController | null = null;
function close() {
  list.hidden = true;
  input.setAttribute("aria-expanded", "false");
  input.removeAttribute("aria-activedescendant");
  active = -1;
}
function cancel() {
  clearTimeout(timer);
  controller?.abort();
  close();
  if (status.textContent === "Searching…") status.textContent = "";
}
function choose(index: number) {
  input.value = items[index].symbol;
  source.value = "reference";
  close();
  status.textContent = items[index].name;
  input.focus();
}
function highlight() {
  list
    .querySelectorAll("[role=option]")
    .forEach((e, i) => e.setAttribute("aria-selected", String(i === active)));
  if (active >= 0)
    input.setAttribute("aria-activedescendant", `gene-option-${active}`);
  else input.removeAttribute("aria-activedescendant");
}
input.addEventListener("input", () => {
  source.value = "reference";
  clearTimeout(timer);
  controller?.abort();
  close();
  const query = input.value.trim();
  if (query.length < 2) {
    status.textContent = "";
    return;
  }
  timer = setTimeout(async () => {
    const request = new AbortController();
    controller = request;
    status.textContent = "Searching…";
    try {
      const result = await suggestGenes(query, request.signal);
      if (request.signal.aborted || input.value.trim() !== query) return;
      items = result;
      active = -1;
      list.innerHTML = items
        .map(
          (g, i) =>
            `<li id="gene-option-${i}" role="option" aria-selected="false" data-option="${i}"><strong>${esc(g.symbol)}</strong><span>${esc(g.name)}</span></li>`,
        )
        .join("");
      list.hidden = !items.length;
      input.setAttribute("aria-expanded", String(items.length > 0));
      status.textContent = items.length
        ? `${items.length} matches · ↑ ↓ to select, Enter to confirm.`
        : "No suggestions found. Try another symbol or name.";
    } catch {
      if (!request.signal.aborted)
        status.textContent = "Gene search unavailable. Try an example dataset.";
    }
  }, 250);
});
list.addEventListener("pointerdown", (e) => e.preventDefault());
list.addEventListener("click", (e) => {
  const option = (e.target as HTMLElement).closest<HTMLElement>(
    "[data-option]",
  );
  if (option) choose(Number(option.dataset.option));
});
input.addEventListener("keydown", (e) => {
  if (e.key === "Escape") cancel();
  if (list.hidden) return;
  if (e.key === "ArrowDown" || e.key === "ArrowUp") {
    e.preventDefault();
    active =
      active < 0
        ? e.key === "ArrowDown"
          ? 0
          : items.length - 1
        : (active + (e.key === "ArrowDown" ? 1 : -1) + items.length) %
          items.length;
    highlight();
  }
  if (e.key === "Enter" && active >= 0) {
    e.preventDefault();
    choose(active);
  }
});
input.addEventListener("blur", cancel);
form.addEventListener("submit", () => {
  controller?.abort();
  input.value = input.value.trim().toUpperCase();
});
document.querySelectorAll<HTMLButtonElement>("[data-gene]").forEach((button) =>
  button.addEventListener("click", () => {
    controller?.abort();
    clearTimeout(timer);
    input.value = button.dataset.gene!;
    source.value = "sample";
    close();
    form.querySelector<HTMLInputElement>(
      `input[value="${button.dataset.mode}"]`,
    )!.checked = true;
    status.textContent =
      input.value === "MYC"
        ? "MYC · rE2G predictions"
        : `${input.value} · illustrative data`;
  }),
);
