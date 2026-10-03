# animation/  (Operator C, branch `feat/animation`)
Blender explainer: chromatin loop forms, enhancer contacts promoter, gene expression increases (gene dosage up).
Keep it a small parametric `bpy` script (`scene.py`). Renders go to `renders/` (git-ignored); share finals via release/link.

## Rendering on Modal
One-time setup (Python 3.11+, venv is git-ignored):
```
py -3.11 -m venv animation/.venv
animation/.venv/Scripts/pip install modal
animation/.venv/Scripts/modal setup      # browser login
```
Then `animation/.venv/Scripts/modal run animation/modal_render.py --still 1` (one frame) or without `--still` for the full film.
Presets (`draft` 720p / `review` 1080p / `final` 4K) are in `scene.py`. `scene.py` is currently a pipeline smoke-test scene, not the storyboard.
