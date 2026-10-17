"""Portable AthenaK GRMHD binary -> accepted square scientific rendering."""
import argparse
import gc
import hashlib
import os
from pathlib import Path
import sys
from common import ROOT, digest, header, read_json, settings, validate, versions, write_json


def render_work(work,cfg,output,passes=False):
    from prepare import prepared_signature
    from renderer import read_poly,render
    work = Path(work)
    if Path(output).suffix.lower()!=".png": raise ValueError("Output must be a .png path")
    meta = read_json(work/"prepared.json")
    samples = read_json(work/"samples.json")
    if meta["signature"]!=prepared_signature(samples,cfg):
        raise ValueError("Prepared geometry/settings mismatch: rerun prepare or run")
    info = render(read_poly(work/"disk.vtk"),read_poly(work/"fields.vtk"),samples["spin"],cfg,output,passes)
    info.update(settings=cfg,input=samples["input"],versions=versions(),
                seed_motion=meta["seed_motion"],seed_positions=meta["seed_positions"],
                prepared_signature=meta["signature"],
                png_sha256=hashlib.sha256(Path(output).read_bytes()).hexdigest())
    if info["missing_ids"]: print("WARNING: missing selected traces:",info["missing_ids"],flush=True)
    write_json(Path(output).with_suffix(".json"),info)
    print("Saved",output,flush=True)
    return info


def run(source,work,cfg,output,passes=False):
    from prepare import sample,prepare
    sample(source,work,cfg)
    prepare(work,cfg)
    return render_work(work,cfg,output,passes)


def doctor(render_test=False):
    import vtk
    from renderer import capture,horizon_actor,scene
    info = dict(python=sys.version,versions=versions(),DISPLAY=os.environ.get("DISPLAY"),
                window_backend=vtk.vtkRenderWindow().GetClassName())
    if render_test:
        r = settings()["render"]
        r.update(width=128,height=128,scale=4.)
        renderer,window = scene(r)
        renderer.SetBackground(0,0,0)
        shell = horizon_actor(.9375,white=True)
        shell.SetScale(2.5,2.5,2.5)
        shell.GetProperty().SetOpacity(.25)
        renderer.AddActor(shell)
        renderer.AddActor(horizon_actor(.9375))
        pixels = capture(window,r)
        info["depth_peeling_used"] = bool(renderer.GetLastRenderingUsedDepthPeeling())
        info["image_range"] = [int(pixels.min()),int(pixels.max())]
        window.Finalize()
        if not info["depth_peeling_used"] or pixels.max()==pixels.min():
            raise RuntimeError("Rendering/depth-peeling environment test failed")
    print(__import__("json").dumps(info,indent=2))
    return info


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command",required=True)
    p = commands.add_parser("doctor")
    p.add_argument("--render",action="store_true")
    p = commands.add_parser("seed-info")
    p.add_argument("source")
    p.add_argument("--settings")
    p.add_argument("--points",action="store_true")
    for command in ("sample","prepare","render","run"):
        p = commands.add_parser(command)
        if command in ("sample","run"): p.add_argument("source")
        p.add_argument("--work",required=True)
        p.add_argument("--settings")
        if command in ("run","render"):
            p.add_argument("--out",required=True)
            p.add_argument("--width",type=int)
            p.add_argument("--height",type=int)
            p.add_argument("--passes",action="store_true")
    p = commands.add_parser("batch")
    p.add_argument("--schedule",required=True)
    p.add_argument("--work-root",required=True)
    p.add_argument("--out-dir",required=True)
    p.add_argument("--start",type=int,default=0)
    p.add_argument("--stop",type=int)
    p.add_argument("--width",type=int)
    p.add_argument("--height",type=int)
    p.add_argument("--passes",action="store_true")
    args = parser.parse_args()
    if args.command=="doctor":
        doctor(args.render)
        return
    if args.command=="seed-info":
        from prepare import seed_motion,seed_positions
        cfg = settings(args.settings)
        h = header(args.source)
        info = seed_motion(cfg,h["time"],h["spin"])
        if args.points: info["seed_positions"] = seed_positions(cfg,h["time"],h["spin"]).tolist()
        print(__import__("json").dumps(info,indent=2))
        return
    if args.command=="batch":
        schedule = read_json(args.schedule)
        if schedule.get("schema_version")!=1: raise ValueError("Unsupported schedule schema")
        if args.start<0 or (args.stop is not None and args.stop<args.start): raise ValueError("Invalid batch slice")
        for frame in schedule["frames"][args.start:args.stop]:
            cfg = frame["settings"]
            if args.width: cfg["render"]["width"] = args.width
            if args.height: cfg["render"]["height"] = args.height
            validate(cfg)
            if header(frame["source"])["time"]!=frame["time"]: raise ValueError("Scheduled time does not match bin header")
            index = f'{frame["index"]:06d}'
            output = Path(args.out_dir)/(index+".png")
            # Reuse only a matching completed render, never just an existing PNG.
            if output.exists() and output.with_suffix(".json").exists():
                saved = read_json(output.with_suffix(".json"))
                if (saved["settings"]==cfg and saved["input"]==header(frame["source"])
                    and saved["versions"]==versions()
                    and saved["png_sha256"]==hashlib.sha256(output.read_bytes()).hexdigest()):
                    print("Keeping matching completed frame",index,flush=True)
                    continue
            run(frame["source"],Path(args.work_root)/index,cfg,output,args.passes)
            gc.collect()
        return
    cfg = settings(args.settings)
    for key in ("width","height"):
        v = getattr(args,key,None)
        if v is not None: cfg["render"][key] = v
    validate(cfg)
    if args.command=="sample":
        from prepare import sample
        sample(args.source,args.work,cfg)
    elif args.command=="prepare":
        from prepare import prepare
        prepare(args.work,cfg)
    elif args.command=="render":
        render_work(args.work,cfg,args.out,args.passes)
    else:
        run(args.source,args.work,cfg,args.out,args.passes)


if __name__=="__main__":
    main()
