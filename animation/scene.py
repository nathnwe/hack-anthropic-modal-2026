"""Parametric Blender scene for the explainer animation (30 s, 24 fps = 720 frames).

Run by Modal (see modal_render.py), or locally:
    blender -b -noaudio --python animation/scene.py -- --preset draft --start 1 --end 1 --out animation/renders/test/

STATUS: pipeline smoke-test scene only (pastel glossy membrane + placeholder blobs).
It is NOT the storyboard and uses NO PDB structures yet. Replace build_scene() once the
storyboard is agreed.
"""
import argparse
import math
import sys

import bpy

FPS = 24
DURATION_S = 30
TOTAL_FRAMES = FPS * DURATION_S  # 720

# Draft must stay cheap: whole 720-frame run in <= ~15 min wall clock when fanned out.
PRESETS = {
    "draft": dict(res=(1280, 720), samples=24, denoise=True),
    "review": dict(res=(1920, 1080), samples=64, denoise=True),
    "final": dict(res=(3840, 2160), samples=256, denoise=True),
}

# Soft pastel palette (reference style: pink / lilac / aqua / peach on a blue-grey backdrop)
PASTEL = {
    "pink": (0.96, 0.62, 0.70, 1),
    "lilac": (0.72, 0.60, 0.90, 1),
    "aqua": (0.50, 0.85, 0.82, 1),
    "peach": (0.98, 0.72, 0.55, 1),
    "membrane": (0.90, 0.95, 0.98, 1),
    "backdrop": (0.45, 0.62, 0.72, 1),
}


def parse_args():
    argv = sys.argv[sys.argv.index("--") + 1:] if "--" in sys.argv else []
    p = argparse.ArgumentParser()
    p.add_argument("--preset", default="draft", choices=PRESETS)
    p.add_argument("--start", type=int, default=1)
    p.add_argument("--end", type=int, default=TOTAL_FRAMES)
    p.add_argument("--out", default="//renders/")
    return p.parse_args(argv)


def glossy(name, color, roughness=0.25, coat=0.6):
    mat = bpy.data.materials.new(name)
    mat.use_nodes = True
    bsdf = mat.node_tree.nodes["Principled BSDF"]
    bsdf.inputs["Base Color"].default_value = color
    bsdf.inputs["Roughness"].default_value = roughness
    for key, val in (("Coat Weight", coat), ("Subsurface Weight", 0.15)):
        if key in bsdf.inputs:
            bsdf.inputs[key].default_value = val
    return mat


def build_scene():
    bpy.ops.wm.read_factory_settings(use_empty=True)
    scene = bpy.context.scene

    # Membrane slab
    bpy.ops.mesh.primitive_plane_add(size=14, location=(0, 0, 0))
    mem = bpy.context.object
    mem.data.materials.append(glossy("membrane", PASTEL["membrane"], roughness=0.12, coat=1.0))

    # Placeholder blobs that bob gently (stand-ins for proteins)
    for i, (col, x, y) in enumerate([("pink", -2, 0), ("lilac", 0, 1), ("aqua", 2, -0.5), ("peach", 0.8, -1.8)]):
        bpy.ops.mesh.primitive_uv_sphere_add(radius=0.7, location=(x, y, 1.0), segments=48, ring_count=24)
        s = bpy.context.object
        bpy.ops.object.shade_smooth()
        s.data.materials.append(glossy(f"blob_{col}", PASTEL[col]))
        for f, z in ((1, 1.0), (TOTAL_FRAMES // 2, 1.5 + 0.1 * i), (TOTAL_FRAMES, 1.0)):
            s.location.z = z
            s.keyframe_insert("location", index=2, frame=f)

    # Camera: slow push-in, framing the whole scene
    cam_data = bpy.data.cameras.new("cam")
    cam_data.lens = 50
    cam = bpy.data.objects.new("cam", cam_data)
    scene.collection.objects.link(cam)
    scene.camera = cam
    for f, r in ((1, 11), (TOTAL_FRAMES, 9)):
        cam.location = (0, -r * 0.9, r * 0.55)
        cam.rotation_euler = (math.radians(65), 0, 0)
        cam.keyframe_insert("location", frame=f)

    # Lights + world
    bpy.ops.object.light_add(type="AREA", location=(-4, -4, 8))
    key = bpy.context.object
    key.data.energy = 900
    key.data.size = 6
    world = bpy.data.worlds.new("world")
    scene.world = world
    world.use_nodes = True
    bg = world.node_tree.nodes["Background"]
    bg.inputs["Color"].default_value = PASTEL["backdrop"]
    bg.inputs["Strength"].default_value = 0.8


def configure_render(args):
    preset = PRESETS[args.preset]
    scene = bpy.context.scene
    scene.render.engine = "CYCLES"
    cy = scene.cycles
    cy.device = "GPU"
    cy.samples = preset["samples"]
    cy.use_denoising = preset["denoise"]
    cy.denoiser = "OPTIX"
    cy.use_adaptive_sampling = True
    scene.render.resolution_x, scene.render.resolution_y = preset["res"]
    scene.render.resolution_percentage = 100
    scene.render.fps = FPS
    scene.frame_start, scene.frame_end = args.start, args.end
    scene.render.image_settings.file_format = "PNG"
    scene.render.filepath = args.out.rstrip("/") + "/frame_"

    prefs = bpy.context.preferences.addons["cycles"].preferences
    for dev_type in ("OPTIX", "CUDA"):
        try:
            prefs.compute_device_type = dev_type
            prefs.get_devices()
            gpus = [d for d in prefs.devices if d.type == dev_type]
            if gpus:
                for d in prefs.devices:
                    d.use = d.type == dev_type
                print(f"[scene] Cycles device type: {dev_type} ({len(gpus)} GPU)")
                return
        except TypeError:
            continue
    print("[scene] WARNING: no GPU found, rendering on CPU")
    cy.device = "CPU"
    cy.use_denoising = True
    cy.denoiser = "OPENIMAGEDENOISE"


if __name__ == "__main__":
    args = parse_args()
    build_scene()
    configure_render(args)
    bpy.ops.render.render(animation=True)
