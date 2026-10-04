"""Real K562 inputs for one 2 Mb window, range-read from public servers.

Nothing here downloads a whole file. The Hi-C is 33.8 GB; we ask the server
for the bytes of one 2 Mb x 2 Mb block. Everything is cached under data/raw/
(git-ignored) and a manifest records exactly what was read and when.

Sources
-------
Hi-C        ENCODE ENCFF621AIY, K562 intact Hi-C (ENCSR479XDG), GRCh38,
            mapq-thresholded. Balanced with Juicer SCALE: this file carries no
            KR vectors at 5 kb, so this is the nearest published equivalent and
            is a deviation from ABC's KR normalisation.
DNase       ENCODE ENCSR000EKS (K562). Peaks ENCFF070TML; signal ENCFF972GVB.
Genes       UCSC ncbiRefSeqCurated, hg38. TSS = the one shared by the most
            isoforms of the gene (the rE2G convention).
CTCF        UCSC jaspar2022 track, TFName == CTCF, restricted to accessible
            sites. (Track is JASPAR 2022; the method doc cites MA0139.1.)
Dosage      ClinGen gene dosage curation list, GRCh38.

NOT covered, and not faked: COSMIC Cancer Gene Census (needs a login) and
DepMap common-essential (download is behind a verification page). The
oncogene and essential-gene axes of the safety gate are therefore empty in
this run. Only ClinGen haploinsufficiency / triplosensitivity is real.
"""

from __future__ import annotations

import collections
import csv
import datetime as dt
import http.client
import io
import json
import time
import urllib.request
from pathlib import Path

import numpy as np

from pipeline.locus import BIN_KB, Locus

REPO = Path(__file__).resolve().parents[1]
RAW = REPO / "data" / "raw" / "real"

HIC_URL = "https://www.encodeproject.org/files/ENCFF621AIY/@@download/ENCFF621AIY.hic"
PEAK_URL = "https://www.encodeproject.org/files/ENCFF070TML/@@download/ENCFF070TML.bigBed"
SIGNAL_URL = "https://www.encodeproject.org/files/ENCFF972GVB/@@download/ENCFF972GVB.bigWig"
UCSC = "https://api.genome.ucsc.edu/getData/track?genome=hg38;track={track};chrom={chrom};start={start};end={end}"
RNA_URL = "https://www.encodeproject.org/files/ENCFF222UVT/@@download/ENCFF222UVT.tsv"
CLINGEN_URL = "https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv"

BIN = BIN_KB * 1000
PROMOTER_HALF = 2500          # an element within this of a TSS is a promoter, not an enhancer


def _get(url: str, timeout: int = 120, tries: int = 4) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "chromatin-staple-poc/0.1"})
    for attempt in range(tries):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except (OSError, http.client.HTTPException):    # incl. IncompleteRead on big transfers
            if attempt == tries - 1:
                raise
            time.sleep(2 ** attempt)


def _cached_json(name: str, url: str) -> dict:
    path = RAW / name
    if not path.exists():
        RAW.mkdir(parents=True, exist_ok=True)
        path.write_bytes(_get(url))
    return json.loads(path.read_text())


def fetch_hic(chrom: str, start: int, end: int) -> np.ndarray:
    """SCALE-balanced observed contacts at 5 kb for [start, end)^2."""
    cache = RAW / f"hic_{chrom}_{start}_{end}.npy"
    if cache.exists():
        return np.load(cache)
    import hicstraw

    hic = hicstraw.HiCFile(HIC_URL)
    zoom = hic.getMatrixZoomData(chrom, chrom, "observed", "SCALE", "BP", BIN)
    m = np.asarray(zoom.getRecordsAsMatrix(start, end - 1, start, end - 1), dtype=float)
    m = np.nan_to_num(m, nan=0.0)
    m = np.triu(m) + np.triu(m, 1).T          # straw returns the upper triangle
    RAW.mkdir(parents=True, exist_ok=True)
    np.save(cache, m)
    return m


def fetch_genes(chrom: str, start: int, end: int) -> dict[str, dict]:
    """symbol -> {tss, strand, start, end}. TSS shared by the most isoforms."""
    d = _cached_json(f"genes_{chrom}_{start}_{end}.json",
                     UCSC.format(track="ncbiRefSeqCurated", chrom=chrom, start=start, end=end))
    by_gene: dict[str, list[dict]] = collections.defaultdict(list)
    for t in d["ncbiRefSeqCurated"]:
        if t["name"].startswith(("NM_", "NR_")):
            by_gene[t["name2"]].append(t)
    out = {}
    for sym, ts in by_gene.items():
        tss_of = lambda t: t["txStart"] if t["strand"] == "+" else t["txEnd"]
        tss = collections.Counter(tss_of(t) for t in ts).most_common(1)[0][0]
        out[sym] = {"tss": tss, "strand": ts[0]["strand"],
                    "start": min(t["txStart"] for t in ts),
                    "end": max(t["txEnd"] for t in ts),
                    "coding": any(t["name"].startswith("NM_") for t in ts)}
    return out


def fetch_peaks(chrom: str, start: int, end: int) -> list[tuple[int, int]]:
    cache = RAW / f"peaks_{chrom}_{start}_{end}.json"
    if cache.exists():
        return [tuple(x) for x in json.loads(cache.read_text())]
    import pybigtools

    bb = pybigtools.open(PEAK_URL)
    peaks = [(int(r[0]), int(r[1])) for r in bb.records(chrom, start, end)]
    RAW.mkdir(parents=True, exist_ok=True)
    cache.write_text(json.dumps(peaks))
    return peaks


def fetch_signal(chrom: str, spans: list[tuple[int, int]]) -> np.ndarray:
    """Mean DNase signal over each span."""
    lo, hi = min(s for s, _ in spans), max(e for _, e in spans)
    cache = RAW / f"signal_{chrom}_{lo}_{hi}.npy"
    if cache.exists():
        track = np.load(cache)
    else:
        import pybigtools

        bw = pybigtools.open(SIGNAL_URL)       # only touch the network on a cache miss
        track = np.nan_to_num(np.asarray(bw.values(chrom, lo, hi), dtype=float), nan=0.0)
        RAW.mkdir(parents=True, exist_ok=True)
        np.save(cache, track)
    return np.array([track[s - lo:e - lo].mean() for s, e in spans])


def fetch_ctcf(chrom: str, start: int, end: int) -> list[int]:
    """CTCF motif centres in the window. The UCSC endpoint returns the whole
    chromosome (130+ MB), so only the filtered centres are cached."""
    cache = RAW / f"ctcf_sites_{chrom}_{start}_{end}.json"
    if cache.exists():
        return json.loads(cache.read_text())
    d = json.loads(_get(UCSC.format(track="jaspar2022", chrom=chrom, start=start, end=end), timeout=300))
    sites = [(r["chromStart"] + r["chromEnd"]) // 2 for r in d.get("jaspar2022", [])
             if r.get("TFName") == "CTCF" and start <= r["chromStart"] < end]
    RAW.mkdir(parents=True, exist_ok=True)
    cache.write_text(json.dumps(sites))
    return sites


def fetch_expression(chrom: str, start: int, end: int, genes: dict[str, dict]) -> dict[str, float]:
    """symbol -> K562 TPM (ENCODE polyA+ RNA-seq ENCSR000AEM rep 1).

    ENCODE quantifies by Ensembl ID, our genes are RefSeq symbols, so they are
    joined by TSS: a RefSeq gene takes the TPM of any GENCODE gene whose
    transcript starts within 500 bp on the same strand (highest wins).
    """
    path = RAW / "k562_rna_ENCFF222UVT.tsv"
    if not path.exists():
        RAW.mkdir(parents=True, exist_ok=True)
        path.write_bytes(_get(RNA_URL, timeout=300))
    # transcript (ENST, unversioned) -> gene TPM. UCSC names annotated genes by
    # symbol, not Ensembl ID, so the join goes through transcripts instead.
    tx_tpm = {}
    for row in csv.DictReader(path.open(), delimiter="\t"):
        if row["gene_id"].startswith("ENSG"):
            for tx in row["transcript_id(s)"].split(","):
                tx_tpm[tx.split(".")[0]] = float(row["TPM"])
    d = _cached_json(f"gencode_{chrom}_{start}_{end}.json",
                     UCSC.format(track="knownGene", chrom=chrom, start=start, end=end))
    starts = [(t["chromStart"] if t["strand"] == "+" else t["chromEnd"], t["strand"],
               tx_tpm.get(t["name"].split(".")[0], 0.0)) for t in d["knownGene"]]
    out = {}
    for sym, g in genes.items():
        hits = [v for tss, st, v in starts if st == g["strand"] and abs(tss - g["tss"]) <= 500]
        out[sym] = max(hits) if hits else 0.0
    return out


def fetch_dosage() -> dict[str, str]:
    """symbol -> 'haploinsufficient' | 'triplosensitive' (ClinGen score 3 = sufficient evidence)."""
    path = RAW / "clingen_dosage.tsv"
    if not path.exists():
        RAW.mkdir(parents=True, exist_ok=True)
        path.write_bytes(_get(CLINGEN_URL))
    flags = {}
    raw = path.read_text().splitlines()
    # The header row itself begins with '#', so find it rather than dropping every '#' line.
    head = next(i for i, l in enumerate(raw) if l.startswith("#Gene Symbol"))
    for row in csv.DictReader(io.StringIO("\n".join(raw[head:])), delimiter="\t"):
        sym = row["#Gene Symbol"]
        if row["Haploinsufficiency Score"].strip() == "3":
            flags[sym] = "haploinsufficient"
        elif row["Triplosensitivity Score"].strip() == "3":
            flags[sym] = "triplosensitive"
    return flags


def load_k562_locus(chrom: str, start: int, end: int) -> tuple[Locus, dict]:
    """Assemble a Locus from the real inputs. Returns (locus, provenance)."""
    assert start % BIN == 0 and end % BIN == 0
    n = (end - start) // BIN
    to_bin = lambda pos: int(np.clip((pos - start) // BIN, 0, n - 1))

    contact = fetch_hic(chrom, start, end)
    assert contact.shape == (n, n), contact.shape
    np.fill_diagonal(contact, 0.0)

    genes_all = fetch_genes(chrom, start, end)
    peaks = fetch_peaks(chrom, start, end)
    signal = fetch_signal(chrom, peaks)

    gene_bins = {g: to_bin(v["tss"]) for g, v in genes_all.items() if v["coding"]}
    tss_bins = np.array(sorted(set(gene_bins.values())))

    accessible = np.zeros(n, dtype=bool)
    el_bins, el_act = [], []
    for (s, e), sig in zip(peaks, signal):
        b = to_bin((s + e) // 2)
        accessible[b] = True
        near_tss = any(abs((s + e) // 2 - genes_all[g]["tss"]) < PROMOTER_HALF
                       for g in gene_bins)
        if not near_tss and sig > 0:
            el_bins.append(b)
            el_act.append(sig)
    # one element per bin: keep the strongest
    best: dict[int, float] = {}
    for b, a in zip(el_bins, el_act):
        best[b] = max(best.get(b, 0.0), a)
    elements = np.array(sorted(best))
    activities = np.array([best[b] for b in elements])

    ctcf = sorted({to_bin(p) for p in fetch_ctcf(chrom, start, end) if accessible[to_bin(p)]})

    expression = fetch_expression(chrom, start, end, genes_all)
    dosage = fetch_dosage()
    flags = {g: dosage[g] for g in gene_bins if g in dosage}

    bodies = {g: (to_bin(v["start"]), to_bin(v["end"])) for g, v in genes_all.items() if g in gene_bins}
    loc = Locus(contact=contact, elements=elements, activities=activities,
                genes=gene_bins, gene_bodies=bodies, ctcf=ctcf, tad_edges=[0, n],
                accessible=accessible, flags=flags, illustrative=False,
                expression={g: expression[g] for g in gene_bins})
    prov = {
        "retrieved": dt.date.today().isoformat(),
        "window": f"{chrom}:{start}-{end}", "bins": n, "bin_bp": BIN,
        "hic": {"file": "ENCFF621AIY", "experiment": "ENCSR479XDG", "cell": "K562",
                "assembly": "GRCh38", "balancing": "Juicer SCALE (no KR at 5 kb)"},
        "dnase": {"experiment": "ENCSR000EKS", "peaks": "ENCFF070TML", "signal": "ENCFF972GVB"},
        "genes": "UCSC ncbiRefSeqCurated hg38", "ctcf": "UCSC jaspar2022, CTCF, accessible only",
        "rna": "ENCODE polyA+ RNA-seq ENCSR000AEM rep 1 (ENCFF222UVT), joined to RefSeq by TSS +/-500 bp",
        "dosage": "ClinGen gene dosage curation list GRCh38 (score 3 only)",
        "not_covered": ["COSMIC Cancer Gene Census (login)", "DepMap common essential (verification page)"],
        "counts": {"dnase_peaks": len(peaks), "candidate_elements": int(len(elements)),
                   "coding_genes": len(gene_bins), "ctcf_sites": len(ctcf),
                   "flagged_genes": len(flags)},
    }
    return loc, prov
