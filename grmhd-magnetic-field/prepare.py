"""AthenaK bin -> stretched Cartesian samples -> display surfaces/native traces."""
import gc
from pathlib import Path
import numpy as np
from scipy.ndimage import gaussian_filter
import vtk
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy
from common import digest, header, read_json, versions, write_json
from renderer import stable_curve_ids, taper

VARIABLES = ["dens","bcc1","bcc2","bcc3"]
TRACE_KEYS = ("grid_stride","seed_z","seed_radii","seed_phi_count","initial_step",
              "minimum_step","maximum_step","maximum_propagation","maximum_steps","terminal_speed",
              "seed_motion","seed_reference_time","seed_phase","seed_omega_fraction")


def axes(cfg):
    c = cfg["sampling"]
    n,stretch = c["n"],c["stretch"]
    result = [(stretch*np.sinh(np.linspace(-np.arcsinh(e/stretch),np.arcsinh(e/stretch),n))).astype(np.float32)
              for e in c["extents"]]
    if not all(np.all(np.diff(a)>0) for a in result): raise ValueError("Nonmonotonic axes")
    return result


def sample(source, work, cfg):
    from wongutils.grmhd.athenak import AthenaKSnapshot
    work = Path(work)
    work.mkdir(parents=True,exist_ok=True)
    meta = header(source)
    if not set(VARIABLES).issubset(meta["variables"]): raise ValueError("Snapshot needs dens,bcc1,bcc2,bcc3")
    for i,extent in enumerate(cfg["sampling"]["extents"],1):
        mesh = meta["parameters"]["mesh"]
        if int(mesh[f"nx{i}"])<=1: raise ValueError("A full 3D snapshot is required")
        if -extent<float(mesh[f"x{i}min"]) or extent>float(mesh[f"x{i}max"]):
            raise ValueError("Sampling box exceeds native domain on axis " + str(i))
    environment = versions()
    signature = digest(dict(schema=1,input=meta,sampling=cfg["sampling"],
                            reader=environment["reader_sha256"],meshblocks=environment["meshblocks_sha256"]))
    if (work/"samples.json").exists() and (work/"samples.npy").exists():
        existing = read_json(work/"samples.json")
        if existing["signature"]==signature:
            print("Keeping matching sampled cache",flush=True)
            return existing
    print("Loading",source,flush=True)
    snapshot = AthenaKSnapshot(str(source),variable_mapping=VARIABLES)
    if not np.isclose(float(snapshot.header["coord"]["a"]),meta["spin"]):
        raise ValueError("Header spin mismatch")
    x,y,z = axes(cfg)
    n = cfg["sampling"]["n"]
    temp = work/"samples.tmp.npy"
    cache = np.lib.format.open_memmap(temp,mode="w+",dtype=np.float32,shape=(n,n,n,4))
    for start in range(0,n,cfg["sampling"]["slab"]):
        stop = min(start+cfg["sampling"]["slab"],n)
        X,Y,Z = np.broadcast_arrays(x[None,None,:],y[None,:,None],z[start:stop,None,None])
        values = snapshot.get_primitives_at(X,Y,Z,interpolation="linear")
        if values.shape != X.shape+(4,) or not np.isfinite(values).all() or values[...,0].min()<0:
            raise ValueError("Invalid interpolated primitives in slab " + str(start))
        cache[start:stop] = values
        if start%40==0 or stop==n: print(f"Sampled {stop}/{n} z planes",flush=True)
    cache.flush()
    del cache,snapshot,values
    gc.collect()
    temp.replace(work/"samples.npy")
    result = dict(signature=signature,input=meta,spin=meta["spin"],time=meta["time"],
                  sampling=cfg["sampling"],shape=[n,n,n,4],axis_order="z,y,x,component",
                  components=VARIABLES,coordinates="Cartesian Kerr-Schild, code length units rg",
                  interpolation="wongutils trilinear, native ghost zones populated",versions=environment)
    write_json(work/"samples.json",result)
    return result


def grid(cfg, stride=1):
    result = vtk.vtkRectilinearGrid()
    coords = [a[::stride] for a in axes(cfg)]
    result.SetDimensions(*[len(a) for a in coords])
    for setter,axis in zip((result.SetXCoordinates,result.SetYCoordinates,result.SetZCoordinates),coords):
        setter(numpy_to_vtk(axis.copy(),deep=True))
    return result


def write_poly(mesh,path):
    path = Path(path)
    temp = path.with_name(path.stem+".tmp.vtk")
    writer = vtk.vtkPolyDataWriter()
    writer.SetFileName(str(temp))
    writer.SetFileVersion(42)
    writer.SetFileTypeToBinary()
    writer.SetInputData(mesh)
    if writer.Write()!=1: raise RuntimeError("VTK write failed: " + str(path))
    temp.replace(path)


def disk(cache,cfg):
    g = cfg["geometry"]
    x,y,_ = axes(cfg)
    radius = np.hypot(y[:,None],x[None,:])
    start,end = g["blend_middle"]
    t = np.clip((radius-start)/(end-start),0,1)
    middle = (t*t*(3-2*t))[None]
    outer = (1-taper(radius,*g["blend_outer"]))[None]
    print("Building radially smoothed display density",flush=True)
    density = gaussian_filter(np.asarray(cache[...,0]),g["sigma_inner"],mode="nearest")
    density *= 1-middle
    medium = gaussian_filter(np.asarray(cache[...,0]),g["sigma_middle"],mode="nearest")
    density += medium*middle
    del medium
    soft = gaussian_filter(np.asarray(cache[...,0]),g["sigma_outer"],mode="nearest")
    density *= 1-outer
    density += soft*outer
    density *= g["density_gain"]
    del soft
    if not np.isfinite(density).all() or density.min()<0: raise ValueError("Invalid display density")
    data = grid(cfg)
    values = numpy_to_vtk(density.reshape(-1),deep=False)
    values.SetName("DisplayDensity")
    data.GetPointData().SetScalars(values)
    append = vtk.vtkAppendPolyData()
    for index,level in enumerate(g["levels"]):
        contour = vtk.vtkContourFilter()
        contour.SetInputData(data)
        contour.SetValue(0,level)
        contour.ComputeNormalsOff()
        smooth = vtk.vtkWindowedSincPolyDataFilter()
        smooth.SetInputConnection(contour.GetOutputPort())
        configure_sinc(smooth,g["inner_iterations"],g["inner_passband"])
        smooth.Update()
        mesh = vtk.vtkPolyData()
        mesh.DeepCopy(smooth.GetOutput())
        if not mesh.GetNumberOfPoints():
            print("WARNING: empty contour",level,flush=True)
            continue
        stronger = vtk.vtkWindowedSincPolyDataFilter()
        stronger.SetInputData(mesh)
        configure_sinc(stronger,g["outer_iterations"],g["outer_passband"])
        stronger.Update()
        xyz = vtk_to_numpy(mesh.GetPoints().GetData())
        xyz_soft = vtk_to_numpy(stronger.GetOutput().GetPoints().GetData())
        weight = (1-taper(np.hypot(xyz[:,0],xyz[:,1]),*g["blend_outer"]))[:,None]
        xyz[:] = xyz*(1-weight)+xyz_soft*weight
        mesh.GetPoints().Modified()
        decimate = vtk.vtkQuadricDecimation()
        decimate.SetInputData(mesh)
        decimate.SetTargetReduction(g["decimation"])
        decimate.VolumePreservationOn()
        normals = vtk.vtkPolyDataNormals()
        normals.SetInputConnection(decimate.GetOutputPort())
        normals.SplittingOff()
        normals.ConsistencyOn()
        normals.AutoOrientNormalsOn()
        normals.Update()
        result = vtk.vtkPolyData()
        result.DeepCopy(normals.GetOutput())
        shell = numpy_to_vtk(np.full(result.GetNumberOfPoints(),index,dtype=np.float32),deep=True)
        shell.SetName("ShellLevel")
        result.GetPointData().AddArray(shell)
        append.AddInputData(result)
        print("Shell",level,result.GetNumberOfPoints(),"points",flush=True)
    append.Update()
    if not append.GetOutput().GetNumberOfPoints(): raise ValueError("No density surfaces at requested levels")
    mesh = vtk.vtkPolyData()
    mesh.DeepCopy(append.GetOutput())
    return mesh


def configure_sinc(smooth,iterations,passband):
    smooth.SetNumberOfIterations(iterations)
    smooth.SetPassBand(passband)
    smooth.BoundarySmoothingOff()
    smooth.FeatureEdgeSmoothingOff()
    smooth.NonManifoldSmoothingOn()
    smooth.NormalizeCoordinatesOn()


def seed_motion(cfg,time,spin):
    """Original Omega_F prescription, anchored to retain the accepted frame."""
    f = cfg["field"]
    if not np.isfinite([time,spin]).all() or not -1<=spin<=1:
        raise ValueError("Invalid snapshot time/spin for footpoint rotation")
    horizon = 1+np.sqrt(1-spin**2)
    omega_h = spin/(2*horizon)
    omega_f = f["seed_omega_fraction"]*omega_h if f["seed_motion"]=="rotate" else 0.
    phase = f["seed_phase"]+omega_f*(time-f["seed_reference_time"])
    return dict(mode=f["seed_motion"],time=time,spin=spin,omega_h=omega_h,omega_f=omega_f,
                reference_time=f["seed_reference_time"],phase_at_reference=f["seed_phase"],
                phase_unwrapped_radians=phase,phase_radians=float(np.remainder(phase,2*np.pi)),
                prescription="phi=phi0+phase_at_reference+fraction*Omega_H*(header_time-reference_time)")


def seed_positions(cfg,time,spin):
    f = cfg["field"]
    phase = seed_motion(cfg,time,spin)["phase_radians"]
    points = []
    for sign in (-1,1):
        for radius in f["seed_radii"]:
            for phi in np.arange(f["seed_phi_count"])*2*np.pi/f["seed_phi_count"]+phase:
                points.append((radius*np.cos(phi),radius*np.sin(phi),sign*f["seed_z"]))
    return np.asarray(points,dtype=np.float32)


def fields(cache,cfg,time,spin):
    f = cfg["field"]
    stride = f["grid_stride"]
    data = grid(cfg,stride)
    vectors = np.ascontiguousarray(cache[::stride,::stride,::stride,1:4]).reshape(-1,3)
    magnetic = numpy_to_vtk(vectors,deep=False)
    magnetic.SetName("B-field")
    data.GetPointData().SetVectors(magnetic)
    points = vtk.vtkPoints()
    for point in seed_positions(cfg,time,spin):
        points.InsertNextPoint(*point)
    seeds = vtk.vtkPolyData()
    seeds.SetPoints(points)
    tracer = vtk.vtkStreamTracer()
    tracer.SetInputData(data)
    tracer.SetSourceData(seeds)
    tracer.SetIntegratorTypeToRungeKutta45()
    tracer.SetIntegrationDirectionToBoth()
    tracer.SetIntegrationStepUnit(vtk.vtkStreamTracer.LENGTH_UNIT)
    tracer.SetInitialIntegrationStep(f["initial_step"])
    tracer.SetMinimumIntegrationStep(f["minimum_step"])
    tracer.SetMaximumIntegrationStep(f["maximum_step"])
    tracer.SetMaximumPropagation(f["maximum_propagation"])
    tracer.SetMaximumNumberOfSteps(f["maximum_steps"])
    tracer.SetTerminalSpeed(f["terminal_speed"])
    tracer.Update()
    result = vtk.vtkPolyData()
    result.DeepCopy(tracer.GetOutput())
    if result.GetNumberOfLines()==0: raise ValueError("No field traces integrated")
    ids = stable_curve_ids(result,points.GetNumberOfPoints())
    identity = numpy_to_vtk(np.asarray(ids,dtype=np.int32),deep=True)
    identity.SetName("StableCurveId")
    result.GetCellData().AddArray(identity)
    print("Traced",result.GetNumberOfLines(),"seed/direction branches",flush=True)
    return result


def prepared_signature(sample_meta,cfg):
    return digest(dict(schema=2,sample=sample_meta["signature"],sampling=cfg["sampling"],
                       geometry=cfg["geometry"],trace={k:cfg["field"][k] for k in TRACE_KEYS},
                       seed_motion=seed_motion(cfg,sample_meta["time"],sample_meta["spin"]),
                       vtk=vtk.vtkVersion.GetVTKVersion()))


def prepare(work,cfg):
    work = Path(work)
    meta = read_json(work/"samples.json")
    if meta["sampling"]!=cfg["sampling"]: raise ValueError("Sampling settings changed: resample the bin file first")
    signature = prepared_signature(meta,cfg)
    if (work/"prepared.json").exists() and all((work/f).exists() for f in ("disk.vtk","fields.vtk")):
        existing = read_json(work/"prepared.json")
        if existing["signature"]==signature:
            print("Keeping matching prepared geometry",flush=True)
            return existing
    vtk.vtkMultiThreader.SetGlobalDefaultNumberOfThreads(cfg["render"]["threads"])
    cache = np.load(work/"samples.npy",mmap_mode="r")
    n = cfg["sampling"]["n"]
    if cache.shape!=(n,n,n,4) or cache.dtype!=np.float32: raise ValueError("Unexpected sample cache schema")
    mesh = disk(cache,cfg)
    write_poly(mesh,work/"disk.vtk")
    stats = dict(points=mesh.GetNumberOfPoints(),cells=mesh.GetNumberOfCells())
    del mesh
    gc.collect()
    lines = fields(cache,cfg,meta["time"],meta["spin"])
    write_poly(lines,work/"fields.vtk")
    result = dict(signature=signature,samples=meta,geometry=cfg["geometry"],
                  trace_settings={k:cfg["field"][k] for k in TRACE_KEYS},disk=stats,
                  seed_motion=seed_motion(cfg,meta["time"],meta["spin"]),
                  seed_positions=seed_positions(cfg,meta["time"],meta["spin"]).tolist(),
                  trace_count=lines.GetNumberOfLines(),versions=versions())
    write_json(work/"prepared.json",result)
    return result
