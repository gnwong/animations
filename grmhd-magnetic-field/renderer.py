"""Accepted scientific surface/field renderer and conventional bloom.

No VisIt, local sessions, external image assets, or image synthesis required.
"""
import os
from pathlib import Path
import numpy as np
from PIL import Image, ImageFilter
from scipy.ndimage import gaussian_filter
import vtk
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy


def taper(radius, start, end):
    t = np.clip((radius-start)/(end-start), 0, 1)
    return np.clip(1-(6*t**5-15*t**4+10*t**3), 0, 1)


def read_poly(path):
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(str(path))
    reader.ReadAllScalarsOn()
    reader.ReadAllVectorsOn()
    reader.Update()
    mesh = reader.GetOutput()
    if mesh.GetNumberOfPoints() == 0: raise ValueError("Empty geometry: " + str(path))
    return mesh


def disk_actor(source, cfg):
    mesh = vtk.vtkPolyData()
    mesh.ShallowCopy(source)
    xyz = vtk_to_numpy(mesh.GetPoints().GetData())
    shell = np.rint(vtk_to_numpy(mesh.GetPointData().GetArray("ShellLevel"))).astype(int)
    radius = np.hypot(xyz[:, 0], xyz[:, 1])
    cool = np.array(((36,49,94),(59,77,124),(114,135,171),(159,149,152),(205,168,121)), dtype=float)
    warm = np.array(((169,66,22),(211,87,22),(244,120,29),(255,168,46),(255,225,135)), dtype=float)
    silver = cool.copy()
    silver[:3] = ((50,60,84),(88,104,131),(150,162,184))
    point_cool = cool[shell]
    blend = (1-taper(radius,12,40))[:,None]*cfg["outer_silver"]
    point_cool = point_cool*(1-blend)+silver[shell]*blend
    if cfg["outer_blue"]:
        target = point_cool*np.array([.62,.84,1.30])
        weights = np.array([.2126,.7152,.0722])
        target *= ((point_cool@weights)/np.maximum(target@weights,1e-6))[:,None]
        blend = (1-taper(radius,18,60))[:,None]*cfg["outer_blue"]
        point_cool = point_cool*(1-blend)+target*blend
    heat = np.exp(-((radius/cfg["warm_radius"])**2))
    rgb = (point_cool*(1-heat[:,None])+warm[shell]*heat[:,None])*cfg["color_gain"]
    alpha = np.array((.055,.085,.14,.40,.88))[shell]*cfg["opacity"]*(.85+.35*heat)
    alpha *= taper(radius,cfg["fade_start"],cfg["fade_end"])
    outside = 1-taper(radius,cfg["feather_start"],cfg["feather_end"])
    normals = mesh.GetPointData().GetNormals()
    if normals is None: raise ValueError("Prepared disk lacks normals")
    theta, phi = np.deg2rad([cfg["inclination"],cfg["azimuth"]])
    view = np.array([np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),-np.cos(theta)])
    facing = np.clip(np.abs(vtk_to_numpy(normals)@view),0,1)
    edge = .2+.8*facing**.7
    strength = np.where(shell<=2,outside,0)
    alpha *= 1-strength+strength*edge*np.asarray(cfg["feather_shell_weights"])[shell]
    rgba = np.column_stack((np.clip(rgb,0,255),np.clip(alpha*255,0,255))).astype(np.uint8)
    colors = numpy_to_vtk(rgba,deep=True,array_type=vtk.VTK_UNSIGNED_CHAR)
    colors.SetName("DisplayRGBA")
    mesh.GetPointData().SetScalars(colors)
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(mesh)
    mapper.SetColorModeToDirectScalars()
    mapper.SetScalarModeToUsePointData()
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    prop = actor.GetProperty()
    prop.SetInterpolationToPhong()
    prop.SetAmbient(cfg["ambient"])
    prop.SetDiffuse(cfg["diffuse"])
    prop.SetSpecular(cfg["specular"])
    prop.SetSpecularPower(cfg["specular_power"])
    return actor


def stable_curve_ids(poly, seed_count):
    """IDs remain seed+direction based even if other seeds fail to integrate."""
    seeds = poly.GetCellData().GetArray("SeedIds")
    times = poly.GetPointData().GetArray("IntegrationTime")
    if seeds is None or times is None: raise ValueError("Field traces lack seed/time identities")
    result = []
    for i in range(poly.GetNumberOfLines()):
        seed = int(seeds.GetTuple1(i))
        cell = poly.GetCell(i)
        final_time = times.GetTuple1(cell.GetPointId(cell.GetNumberOfPoints()-1))
        if not 0 <= seed < seed_count or final_time == 0:
            raise ValueError("Invalid seed/direction identity")
        result.append(seed+(seed_count if final_time<0 else 0))
    if len(result) != len(set(result)): raise ValueError("Duplicate seed/direction traces")
    return result


def field_actor(source, cfg):
    seed_count = 2*len(cfg["seed_radii"])*cfg["seed_phi_count"]
    identities = stable_curve_ids(source,seed_count)
    selected_ids = np.random.default_rng(cfg["selection_seed"]).permutation(2*seed_count)[:cfg["line_count"]].tolist()
    wanted = set(selected_ids)
    lines = vtk.vtkCellArray()
    source_lines = source.GetLines()
    source_lines.InitTraversal()
    ids = vtk.vtkIdList()
    index = 0
    while source_lines.GetNextCell(ids):
        if identities[index] in wanted: lines.InsertNextCell(ids)
        index += 1
    selected = vtk.vtkPolyData()
    selected.SetPoints(source.GetPoints())
    selected.SetLines(lines)
    if lines.GetNumberOfCells()==0:
        raise ValueError("None of the fixed selected traces integrated; review seed placement")
    spline = vtk.vtkSplineFilter()
    spline.SetInputData(selected)
    spline.SetSubdivideToLength()
    spline.SetLength(cfg["spline_length"])
    tube = vtk.vtkTubeFilter()
    tube.SetInputConnection(spline.GetOutputPort())
    tube.SetRadius(cfg["line_radius"])
    tube.SetNumberOfSides(cfg["tube_sides"])
    tube.CappingOff()
    tube.Update()
    mesh = vtk.vtkPolyData()
    mesh.DeepCopy(tube.GetOutput())
    xyz = vtk_to_numpy(mesh.GetPoints().GetData())
    alpha = cfg["line_opacity"]*taper(np.abs(xyz[:,2]),*cfg["fade_z"])
    alpha *= taper(np.hypot(xyz[:,0],xyz[:,1]),*cfg["fade_radius"])
    rgba = np.empty((len(xyz),4),dtype=np.uint8)
    rgba[:,:3] = (215,229,255)
    rgba[:,3] = np.clip(alpha*255,0,255)
    colors = numpy_to_vtk(rgba,deep=True,array_type=vtk.VTK_UNSIGNED_CHAR)
    colors.SetName("FieldRGBA")
    mesh.GetPointData().SetScalars(colors)
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(mesh)
    mapper.SetColorModeToDirectScalars()
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetAmbient(1)
    actor.GetProperty().SetDiffuse(0)
    actor.GetProperty().SetSpecular(0)
    missing = sorted(wanted-set(identities))
    return actor, dict(requested_ids=selected_ids,missing_ids=missing,visible_trace_count=len(wanted)-len(missing))


def horizon_actor(spin, white=False):
    rh = 1+np.sqrt(1-spin**2)
    sphere = vtk.vtkSphereSource()
    sphere.SetRadius(1)
    sphere.SetThetaResolution(96)
    sphere.SetPhiResolution(64)
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(sphere.GetOutputPort())
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.SetScale(np.sqrt(rh**2+spin**2),np.sqrt(rh**2+spin**2),rh)
    actor.GetProperty().SetColor((1,1,1) if white else (.001,.002,.005))
    actor.GetProperty().SetAmbient(1)
    actor.GetProperty().SetDiffuse(0)
    actor.GetProperty().SetSpecular(0)
    return actor


def scene(cfg):
    vtk.vtkMultiThreader.SetGlobalDefaultNumberOfThreads(cfg["threads"])
    window = vtk.vtkRenderWindow()
    if window.IsA("vtkXOpenGLRenderWindow") and not os.environ.get("DISPLAY"):
        raise RuntimeError("This VTK needs X11. Run under xvfb-run or an accessible GPU X server; see README.txt.")
    window.SetOffScreenRendering(1)
    window.SetSize(cfg["width"]*cfg["supersampling"],cfg["height"]*cfg["supersampling"])
    window.SetAlphaBitPlanes(1)
    window.SetMultiSamples(0)
    renderer = vtk.vtkRenderer()
    renderer.SetUseDepthPeeling(True)
    renderer.SetMaximumNumberOfPeels(cfg["peels"])
    renderer.SetOcclusionRatio(cfg["occlusion_ratio"])
    renderer.AutomaticLightCreationOff()
    window.AddRenderer(renderer)
    theta,phi = np.deg2rad([cfg["inclination"],cfg["azimuth"]])
    normal = np.array([np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),-np.cos(theta)])
    right = np.cross([0.,0.,1.],normal)
    right /= np.linalg.norm(right)
    up = np.cross(normal,right)
    camera = renderer.GetActiveCamera()
    camera.SetPosition(*(600*normal))
    camera.SetFocalPoint(0,0,0)
    camera.SetViewUp(*up)
    camera.ParallelProjectionOn()
    camera.SetParallelScale(cfg["scale"])
    camera.SetClippingRange(200,1100)
    camera.Roll(cfg["camera_roll"])
    for direction,color,intensity in (
        (.85*normal-.65*right+.7*up,(1.,.94,.85),cfg["key"]),
        (.5*normal+.85*right+.1*up,(.55,.68,1.),cfg["fill"]),
        (-.45*normal+.4*right+.8*up,(1.,.77,.57),cfg["rim"]),
    ):
        light = vtk.vtkLight()
        light.SetLightTypeToSceneLight()
        light.SetPosition(*(600*direction))
        light.SetFocalPoint(0,0,0)
        light.SetColor(*color)
        light.SetIntensity(intensity)
        renderer.AddLight(light)
    return renderer,window


def capture(window, cfg):
    window.Render()
    buffer = vtk.vtkWindowToImageFilter()
    buffer.SetInput(window)
    buffer.SetInputBufferTypeToRGB()
    buffer.ReadFrontBufferOff()
    buffer.Update()
    data = buffer.GetOutput()
    w,h,_ = data.GetDimensions()
    rgb = np.flipud(vtk_to_numpy(data.GetPointData().GetScalars()).reshape(h,w,3)).copy()
    return np.asarray(Image.fromarray(rgb).resize((cfg["width"],cfg["height"]),Image.Resampling.BOX)).copy()


def linear(rgb):
    rgb = rgb.astype(np.float32)/255
    return np.where(rgb<=.04045,rgb/12.92,((rgb+.055)/1.055)**2.4)


def display(light):
    light = np.where(light<=.72,light,.72+.28*(1-np.exp(-(light-.72)/.28)))
    light = np.clip(light,0,1)
    rgb = np.where(light<=.0031308,12.92*light,1.055*light**(1/2.4)-.055)
    return np.uint8(np.clip(rgb*255+.5,0,255))


def composite(disk, beauty, mask, cfg):
    """Preserve the 3D base; horizon protection affects added glow only.

    Within the projected opaque horizon, the 3D renderer has already hidden
    all material behind it. Visible emission there is foreground emission.
    Retain its bloom while optionally blocking spill from surrounding pixels.
    This is an illustrative bloom control, not a physical PSF calculation.
    """
    base = linear(beauty)
    warm = linear(disk)
    weight = np.clip((warm[...,0]-1.15*warm[...,2])/np.maximum(warm[...,0],.001),0,1)
    warm *= weight[...,None]
    field = np.maximum(linear(beauty)-linear(disk),0)
    unit = (cfg["height"]/cfg["effect_reference_height"])*(cfg["effect_reference_scale"]/cfg["scale"])*(cfg["effect_reference_width"]/720)
    blur = lambda rgb,r: gaussian_filter(rgb,sigma=(r*unit,r*unit,0),mode="constant")
    gas = np.maximum(linear(disk)-linear(np.array([[cfg["background"]]],dtype=np.uint8)),0)*(1-weight[...,None])
    def glow(warm,field,gas):
        warm_glow = cfg["warm_glow"]*(.30*warm+.8*blur(warm,3)+1.4*blur(warm,11)+1.7*blur(warm,28))
        crisp = cfg["warm_crisp"]*np.maximum(warm-blur(warm,1.2),0)
        cool_glow = cfg["cool_glow"]*(.5*blur(field,2)+1.1*blur(field,9)+1.8*blur(field,24))
        cool_glow *= np.array([.26,.45,1.],dtype=np.float32)
        gas_glow = cfg["gas_glow"]*(.6*blur(gas,5)+1.1*blur(gas,16))
        return warm_glow+crisp+cool_glow+gas_glow
    added = glow(warm,field,gas)
    if cfg["horizon_glow_protection"]:
        protect = mask.max(axis=2,keepdims=True).astype(float)/255
        added = added*(1-protect)
        yy,xx = np.nonzero(protect[...,0])
        if len(yy):
            region = np.s_[yy.min():yy.max()+1,xx.min():xx.max()+1]
            coverage = protect[region]
            # Only this small region has nonzero foreground sources. Constant
            # zero boundaries give the same Gaussian values inside it as a
            # full-frame convolution of the masked sources, without that cost.
            front = glow(warm[region]*coverage,field[region]*coverage,gas[region]*coverage)
            added[region] += front*coverage
    light = base+added
    result = Image.fromarray(display(light))
    return result.filter(ImageFilter.UnsharpMask(radius=unit,percent=round(cfg["sharpness"]),threshold=3))


def render(disk, fields, spin, cfg, output, passes=False):
    r = cfg["render"]
    renderer,window = scene(r)
    renderer.SetBackground(*(np.asarray(r["background"])/255))
    renderer.AddActor(disk_actor(disk,r))
    renderer.AddActor(horizon_actor(spin))
    disk_image = capture(window,r)
    field,info = field_actor(fields,cfg["field"])
    renderer.AddActor(field)
    beauty = capture(window,r)
    info["depth_peeling_used"] = bool(renderer.GetLastRenderingUsedDepthPeeling())
    if not info["depth_peeling_used"]:
        window.Finalize()
        raise RuntimeError("Depth peeling unavailable: refusing to substitute incorrect occlusion")
    info["window_backend"] = window.GetClassName()
    info["horizon_mask_policy"] = "Base 3D image never masked/darkened; optional glow-spill protection retains visible foreground glow"
    renderer.RemoveAllViewProps()
    renderer.SetBackground(0,0,0)
    renderer.AddActor(horizon_actor(spin,white=True))
    mask = capture(window,r)
    window.Finalize()
    result = composite(disk_image,beauty,mask,r)
    rotate = lambda im: im.transpose(Image.Transpose.ROTATE_180) if r["screen_rotation"]==180 else im
    output = Path(output)
    output.parent.mkdir(parents=True,exist_ok=True)
    temp = output.with_name(output.stem+".tmp.png")
    rotate(result).save(temp,dpi=(300,300))
    temp.replace(output)
    if passes:
        for suffix,image in (("disk",disk_image),("raw",beauty),("horizon_mask",mask)):
            rotate(Image.fromarray(image)).save(output.with_name(output.stem+"_"+suffix+".png"))
    pixels = np.asarray(result)
    info["pixel_metrics"] = dict(max=int(pixels.max()),std=float(pixels.std()),
        saturated_pixel_fraction=float(np.mean(np.any(pixels>=250,axis=2))))
    return info
