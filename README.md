[![Rewire Bio — visit the website](docs/assets/rewire-bio-banner.png)](https://leandre-tappenden.github.io/rewire-bio-site/)

# 3D Genome Engineering: target shortlisting for chromatin "stapling"

[**Open the Rewire Bio website →**](https://leandre-tappenden.github.io/rewire-bio-site/)

The website is published separately from this source repository. See
[publishing instructions](site/README.md#public-deployment) to deploy updates from `main`.

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
- Animation renderer: NumPy, SciPy, Pillow, FFmpeg (via imageio-ffmpeg); Segoe UI font (Windows system font); DM Sans font (SIL OFL 1.1, `animation/assets/fonts/`)
Add every further library, model and dataset here as it is used.
