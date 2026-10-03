# animation/  (Operator C, branch `feat/animation`)
30 s explainer following `storyboard.md`: cell -> nucleus -> one TAD; enhancer and repressor change
transcription; an engineered dCas9-dCas9 stapler holds enhancer and promoter together.

`procedural.py` is a 2.5D renderer (NumPy, SciPy, Pillow; FFmpeg via imageio-ffmpeg). No Blender.
Everything is illustrative and not to scale; no structure comes from PDB.

## Render
```
pip install -r animation/requirements.txt
python animation/procedural.py --width 1920 --out animation/renders/procedural_1080p
python animation/procedural.py --width 1280 --out animation/renders/preview --frames 0,150,600   # stills
```
Frames are cached in `<out>/frames/` (delete to re-render). Renders go to `renders/` (git-ignored); share finals via release/link.
