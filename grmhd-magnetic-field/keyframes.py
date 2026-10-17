"""Plan reviewed snapshots and interpolate SETTINGS over real snapshot times."""
import argparse
import copy
import glob
import os
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator
from common import ROOT, header, read_json, settings, validate, write_json

# Seed topology/rotation policy, trace count/selection, contours, sampling resolution,
# output dimensions, and smoothing policy stay fixed to avoid temporal popping.
CONTINUOUS = {
    "geometry.density_gain", "field.line_radius", "field.line_opacity",
    *["render."+key for key in (
        "inclination","azimuth","camera_roll","scale","warm_radius","warm_glow",
        "warm_crisp","cool_glow","gas_glow","color_gain","ambient","diffuse",
        "key","rim","fill","specular","specular_power","opacity","fade_start",
        "fade_end","horizon_contrast","sharpness","outer_silver","outer_blue")]
}
ANGLES = {"render.azimuth","render.camera_roll"}


def value(cfg,path):
    group,key = path.split(".")
    return cfg[group][key]


def overrides(base,changes):
    cfg = copy.deepcopy(base)
    for path,v in changes.items():
        if path not in CONTINUOUS: raise ValueError("Non-interpolatable/unknown setting: "+path)
        if isinstance(v,bool) or not isinstance(v,(int,float)): raise ValueError("Expected numeric override: "+path)
        group,key = path.split(".")
        cfg[group][key] = v
    validate(cfg)
    return cfg


def interpolate(base,frames,times):
    if len(frames)<2: raise ValueError("At least two keyframes are required for interpolation")
    anchors = np.array([f["time"] for f in frames],dtype=float)
    if not np.isfinite(anchors).all() or not np.all(np.diff(anchors)>0):
        raise ValueError("Keyframe times must be finite and strictly increasing")
    times = np.asarray(times,dtype=float)
    if not np.isfinite(times).all() or np.any(times<anchors[0]) or np.any(times>anchors[-1]):
        raise ValueError("Schedule times must lie inside the keyframe interval")
    configs = [overrides(base,f.get("overrides",{})) for f in frames]
    paths = sorted(set().union(*(f.get("overrides",{}).keys() for f in frames)))
    curves = {}
    for path in paths:
        points = np.array([value(c,path) for c in configs],dtype=float)
        if path in ANGLES: points = np.rad2deg(np.unwrap(np.deg2rad(points)))
        curves[path] = PchipInterpolator(anchors,points)(times)
    result = []
    for i in range(len(times)):
        cfg = overrides(base,{path:float(values[i]) for path,values in curves.items()})
        result.append(cfg)
    return result


def catalog(pattern):
    paths = sorted(glob.glob(pattern))
    if not paths: raise ValueError("No snapshots match: "+pattern)
    entries = []
    for path in paths:
        h = header(path)
        if not np.isfinite(h["time"]): raise ValueError("Nonfinite simulation time")
        if not {"dens","bcc1","bcc2","bcc3"}.issubset(h["variables"]):
            raise ValueError("Snapshot lacks required density/B fields: "+path)
        entries.append(dict(source=h["source"],time=h["time"],spin=h["spin"]))
    entries.sort(key=lambda x:x["time"])
    if len({e["spin"] for e in entries})!=1: raise ValueError("Mixed black-hole spins: use snapshots from one simulation")
    if not np.all(np.diff([e["time"] for e in entries])>0):
        raise ValueError("Duplicate/non-increasing times: remove duplicate outputs or mixed simulations")
    return entries


def choose(entries,count):
    if count<2 or count>len(entries): raise ValueError("Need 2 <= keyframe count <= available snapshots")
    times = np.array([e["time"] for e in entries])
    result,previous = [],-1
    for i,target in enumerate(np.linspace(times[0],times[-1],count)):
        nearest = int(np.argmin(np.abs(times-target)))
        index = int(np.clip(nearest,previous+1,len(entries)-(count-i)))
        result.append(entries[index])
        previous = index
    return result


def load_plan(path):
    path = Path(path).resolve()
    plan = read_json(path)
    if plan.get("schema_version")!=1: raise ValueError("Unsupported keyframe-plan schema")
    base = settings(path.parent/plan["base_settings"])
    frames = copy.deepcopy(plan["keyframes"])
    for frame in frames:
        frame["source"] = str((path.parent/frame["source"]).resolve())
    return plan,base,frames


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command",required=True)
    p = commands.add_parser("plan")
    p.add_argument("--glob",required=True)
    p.add_argument("--count",type=int,default=100)
    p.add_argument("--settings",default=str(ROOT/"settings.json"))
    p.add_argument("--out",required=True)
    p = commands.add_parser("schedule")
    p.add_argument("--plan",required=True)
    p.add_argument("--glob")
    p.add_argument("--keyframes-only",action="store_true")
    p.add_argument("--require-reviewed",action="store_true")
    p.add_argument("--out",required=True)
    args = parser.parse_args()
    if args.command=="plan":
        settings(args.settings)
        entries = catalog(args.glob)
        frames = choose(entries,args.count)
        destination = Path(args.out).resolve()
        frames = [dict(e,source=os.path.relpath(e["source"],destination.parent),
                       overrides={},reviewed=False,notes="") for e in frames]
        write_json(destination,dict(schema_version=1,
            base_settings=os.path.relpath(Path(args.settings).resolve(),destination.parent),
            interpolation="PCHIP over simulation time; shortest-path unwrapped angles; no extrapolation",
            keyframes=frames))
        print("Selected",len(frames),"snapshots, approximately uniform in simulation time")
    else:
        plan,base,frames = load_plan(args.plan)
        if args.require_reviewed and not all(f.get("reviewed",False) for f in frames):
            raise ValueError("Some keyframes are not reviewed; see QUALITY_CRITERIA.txt")
        key_headers = [header(f["source"]) for f in frames]
        for frame,h in zip(frames,key_headers):
            if frame["time"]!=h["time"]: raise ValueError("Plan time does not match snapshot header")
        if len({h["spin"] for h in key_headers})!=1: raise ValueError("Mixed keyframe spins")
        if args.keyframes_only:
            entries = [dict(source=f["source"],time=f["time"],spin=h["spin"]) for f,h in zip(frames,key_headers)]
        else:
            if not args.glob: parser.error("schedule requires --glob or --keyframes-only")
            entries = [e for e in catalog(args.glob) if frames[0]["time"]<=e["time"]<=frames[-1]["time"]]
        if not entries: raise ValueError("No snapshots lie in the keyframe interval")
        if any(e["spin"]!=key_headers[0]["spin"] for e in entries): raise ValueError("Mixed schedule/keyframe spins")
        configs = interpolate(base,frames,[e["time"] for e in entries])
        write_json(args.out,dict(schema_version=1,keyframe_plan=str(Path(args.plan).resolve()),
            keyframes_reviewed=all(f.get("reviewed",False) for f in frames),
            frames=[dict(e,index=i,settings=c) for i,(e,c) in enumerate(zip(entries,configs))]))
        print("Scheduled",len(entries),"real snapshots; rendering settings interpolated in simulation time")


if __name__=="__main__":
    main()
