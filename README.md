# 3D Genome Engineering: target shortlisting for chromatin "stapling"

Can we reversibly tune disease-relevant gene expression by reshaping 3D genome contacts, using existing regulatory machinery?
Given a dosage-sensitive gene and cell type, this tool proposes which two genomic regions to staple (or unstaple), predicts the direction of effect, flags secondary risks, and ranks the options for a drug-discovery researcher.

- Project context: [docs/background.md](docs/background.md)
- Delivery plan and team split: [docs/plan.md](docs/plan.md)
- Data contract: [contracts/](contracts/)

## Quick start
```bash
uv sync
uv run python contracts/validate.py
```

## Credits
- UniversalEPI: Grover et al., 2026 (see `docs/`)
- ENCODE-rE2G: Gschwind et al., 2026 (see `docs/`)
- ABC model: Fulco et al., 2019, Nat Genet 51:1664-1669
- Data (real-data run): ENCODE K562 Hi-C (ENCFF621AIY), DNase (ENCSR000EKS), RNA-seq (ENCSR000AEM); UCSC Genome Browser API (RefSeq, GENCODE, JASPAR 2022); ClinGen gene dosage curation
- Libraries: hic-straw, pybigtools, NumPy
Add every further library, model and dataset here as it is used.
