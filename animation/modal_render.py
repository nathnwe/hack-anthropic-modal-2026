"""Render the Blender animation on Modal GPUs, fanning frame chunks out across containers.

Usage (from repo root, with animation/.venv active or via its Scripts/modal.exe):
    modal run animation/modal_render.py --still 1                 # one draft frame, to eyeball the look
    modal run animation/modal_render.py                           # full 720-frame draft -> animation/renders/
    modal run animation/modal_render.py --preset review
    modal run animation/modal_render.py --preset final --chunk 4  # 4K, slow; only for finals

Presets live in scene.py. Draft is tuned so the whole film renders in ~15 min wall clock when
chunks run in parallel; the real limit is how many GPU containers your Modal plan lets you run.
"""
import os
import pathlib
import subprocess
import time

import modal

BLENDER_VERSION = "4.2.3"
BLENDER_URL = f"https://download.blender.org/release/Blender4.2/blender-{BLENDER_VERSION}-linux-x64.tar.xz"
GPU = os.environ.get("RENDER_GPU", "L4")
HERE = pathlib.Path(__file__).parent

app = modal.App("tad-animation")
frames_vol = modal.Volume.from_name("tad-animation-frames", create_if_missing=True)
assets_vol = modal.Volume.from_name("tad-animation-assets", create_if_missing=True)  # PDB files etc.

image = (
    modal.Image.debian_slim(python_version="3.11")
    .apt_install(
        "wget", "xz-utils", "ffmpeg", "libxi6", "libxrender1", "libxkbcommon0", "libsm6",
        "libxxf86vm1", "libxfixes3", "libgl1", "libegl1", "libxext6", "libx11-6",
    )
    .run_commands(
        f"wget -q {BLENDER_URL} -O /tmp/blender.tar.xz",
        "mkdir -p /opt/blender && tar -xf /tmp/blender.tar.xz -C /opt/blender --strip-components=1",
        "rm /tmp/blender.tar.xz",
    )
    .add_local_file(HERE / "scene.py", "/app/scene.py")
)

VOLUMES = {"/frames": frames_vol, "/assets": assets_vol}


@app.function(image=image, gpu=GPU, volumes=VOLUMES, timeout=60 * 30, retries=1)
def render_chunk(run: str, preset: str, start: int, end: int) -> str:
    out = f"/frames/{run}"
    os.makedirs(out, exist_ok=True)
    t0 = time.time()
    subprocess.run(
        ["/opt/blender/blender", "-b", "-noaudio", "--python", "/app/scene.py", "--",
         "--preset", preset, "--start", str(start), "--end", str(end), "--out", out],
        check=True,
    )
    frames_vol.commit()
    return f"frames {start}-{end}: {time.time() - t0:.0f}s"


@app.function(image=image, volumes=VOLUMES, timeout=60 * 20)
def encode(run: str, fps: int = 24) -> bytes:
    frames_vol.reload()
    out = f"/tmp/{run}.mp4"
    subprocess.run(
        ["ffmpeg", "-y", "-framerate", str(fps), "-i", f"/frames/{run}/frame_%04d.png",
         "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "16", out],
        check=True,
    )
    return pathlib.Path(out).read_bytes()


@app.function(image=image, volumes=VOLUMES)
def read_frame(run: str, frame: int) -> bytes:
    frames_vol.reload()
    return pathlib.Path(f"/frames/{run}/frame_{frame:04d}.png").read_bytes()


@app.local_entrypoint()
def main(preset: str = "draft", start: int = 1, end: int = 720, chunk: int = 12, still: int = 0):
    out_dir = HERE / "renders"
    out_dir.mkdir(exist_ok=True)
    run = f"{preset}_{int(time.time())}"

    if still:
        print(render_chunk.remote(run, preset, still, still))
        path = out_dir / f"{run}_f{still:04d}.png"
        path.write_bytes(read_frame.remote(run, still))
        print(f"Wrote {path}")
        return

    chunks = [(run, preset, s, min(s + chunk - 1, end)) for s in range(start, end + 1, chunk)]
    print(f"{run}: {len(chunks)} chunks of <= {chunk} frames")
    t0 = time.time()
    for line in render_chunk.starmap(chunks):
        print(line)
    print(f"Rendered in {time.time() - t0:.0f}s; encoding")
    path = out_dir / f"{run}.mp4"
    path.write_bytes(encode.remote(run))
    print(f"Wrote {path} (total {time.time() - t0:.0f}s)")
