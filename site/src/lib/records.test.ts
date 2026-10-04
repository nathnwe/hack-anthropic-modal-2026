import test from "node:test";
import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import {
  type GeneRecord,
  contactsFor,
  contactStrength,
  candidatesFor,
  sortElements,
  selectedElementFor,
  evidenceLink,
  escapeHTML,
} from "./records.ts";

const fixture = (): GeneRecord =>
  JSON.parse(
    readFileSync(
      new URL("../../public/data/myc.json", import.meta.url),
      "utf8",
    ),
  );

test("a missing contact is distinct from a supplied zero; zero sorts before missing", () => {
  const record = fixture();
  record.locus.contacts = [
    { a: record.gene.tss, b: record.locus.cres[1].start, strength: 0 },
  ];
  assert.equal(contactStrength(record, record.locus.cres[0]), null);
  assert.equal(contactStrength(record, record.locus.cres[1]), 0);
  assert.equal(sortElements(record, "contact")[0].id, record.locus.cres[1].id);
});
test("contacts must match the exact TSS and element interval, in either direction", () => {
  const record = fixture(),
    cre = record.locus.cres[0];
  record.locus.contacts = [
    { a: record.gene.tss + 1, b: cre.start, strength: 9 },
    { a: record.gene.tss, b: cre.start - 1, strength: 8 },
    { a: cre.end, b: record.gene.tss, strength: 0.4 },
    { a: cre.start, b: record.gene.tss, strength: 0.7 },
  ];
  assert.equal(contactsFor(record, cre).length, 2);
  assert.equal(contactStrength(record, cre), 0.7);
});
test("direction filters do not relabel incompatible staples or reuse their tiers", () => {
  const record = fixture();
  assert.equal(candidatesFor(record, "down").length, 1);
  assert.equal(candidatesFor(record, "up").length, 0);
  assert.equal(candidatesFor(record, "off").length, 0);
});
test("unsupported rankings preserve record order; supported sort does not mutate input", () => {
  const record = fixture();
  record.locus.cres.reverse();
  const order = record.locus.cres.map((c) => c.id);
  assert.deepEqual(
    sortElements(record, "composite").map((c) => c.id),
    order,
  );
  assert.notEqual(sortElements(record, "re2g")[0].id, order[0]);
  assert.deepEqual(
    record.locus.cres.map((c) => c.id),
    order,
  );
});
test("source links retain supplied URLs and PMIDs; placeholders never become citations", () => {
  assert.equal(
    evidenceLink("Study PMID: 25126789"),
    "https://pubmed.ncbi.nlm.nih.gov/25126789/",
  );
  assert.equal(
    evidenceLink("Source: https://example.org/paper."),
    "https://example.org/paper",
  );
  assert.equal(evidenceLink("PLACEHOLDER PMID: 123456"), null);
  assert.equal(evidenceLink("javascript:alert(1)"), null);
});
test("record content cannot inject markup into result tables", () => {
  assert.equal(
    escapeHTML('<img src=x onerror="alert(1)">'),
    "&lt;img src=x onerror=&quot;alert(1)&quot;&gt;",
  );
});

test("viewer defaults follow each ranking, including opposite orders and an empty locus", () => {
  const record = fixture();
  const [first, second] = record.locus.cres;
  // Test-only values exercise diverging rankings without changing a fixture.
  first.score_rE2G = 0.1;
  second.score_rE2G = 0.9;
  assert.equal(selectedElementFor(record, "contact")?.id, first.id);
  assert.equal(selectedElementFor(record, "re2g")?.id, second.id);
  record.locus.cres = [];
  assert.equal(selectedElementFor(record, "contact"), null);
});
