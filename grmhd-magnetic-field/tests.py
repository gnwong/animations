"""Small portability tests: real binary layout, settings interpolation, trace IDs."""
import tempfile
import unittest
import subprocess
import sys
from pathlib import Path
import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk
from common import ROOT, header, read_json, settings, write_json
from keyframes import choose, interpolate, overrides
from renderer import composite,stable_curve_ids
from prepare import fields,prepared_signature,seed_motion,seed_positions


def fixture(path,time=12.5):
    text = """<mesh>
nx1=8
nx2=8
nx3=8
nghost=2
x1min=-4
x1max=4
x2min=-4
x2max=4
x3min=-4
x3max=4
<meshblock>
nx1=8
nx2=8
nx3=8
<coord>
a=0.9375
""".encode()
    axis = np.arange(8)-3.5
    z,y,x = np.meshgrid(axis,axis,axis,indexing="ij")
    data = np.array([10+.2*x+.3*y+.1*z,np.ones_like(x),.5*np.ones_like(x),-np.ones_like(x)],dtype=np.float32)
    with path.open("wb") as stream:
        stream.write(b"Athena binary output version=1.1\n")
        stream.write(f"pheader_count=5\ntime={time}\ncycle=0\nsize of location=8\nsize of variable=4\n".encode())
        stream.write(b"nvars=4\nvariables: dens bcc1 bcc2 bcc3\n")
        stream.write(f"header offset={len(text)}\n".encode())
        stream.write(text)
        stream.write(np.array([2,9,2,9,2,9],dtype=np.int32).tobytes())
        stream.write(np.array([0,0,0,0],dtype=np.int32).tobytes())
        stream.write(np.array([-4,4,-4,4,-4,4],dtype=np.float64).tobytes())
        stream.write(data.tobytes())


class WorkflowTests(unittest.TestCase):
    def test_horizon_mask_never_darkens_foreground_base(self):
        cfg = settings()["render"]
        cfg.update(width=32,height=32,sharpness=0,warm_glow=0,warm_crisp=0,cool_glow=0,gas_glow=0)
        raw = np.full((32,32,3),[40,60,90],dtype=np.uint8)
        raw[14:18,14:18] = [180,210,240]
        masked = composite(raw,raw,np.full_like(raw,255),cfg)
        unmasked = composite(raw,raw,np.zeros_like(raw),cfg)
        np.testing.assert_array_equal(masked,unmasked)

    def test_foreground_bloom_survives_horizon_protection(self):
        cfg = settings()["render"]
        cfg.update(width=64,height=64,sharpness=0)
        disk = np.full((64,64,3),cfg["background"],dtype=np.uint8)
        beauty = disk.copy()
        disk[30:32,30:32] = [180,100,20]
        beauty[:] = disk
        beauty[32:34,33:35] = [210,225,250]
        mask = np.zeros_like(disk)
        mask[24:40,24:40] = 255
        protected = composite(disk,beauty,mask,cfg)
        cfg["horizon_glow_protection"] = False
        unprotected = composite(disk,beauty,mask,cfg)
        np.testing.assert_array_equal(protected,unprotected)

    def test_rotating_footpoints_and_reference_anchor(self):
        cfg = settings()
        spin = .9375
        t0 = cfg["field"]["seed_reference_time"]
        motion = seed_motion(cfg,t0,spin)
        omega = spin/(4*(1+np.sqrt(1-spin**2)))
        self.assertEqual(motion["phase_radians"],0.)
        self.assertAlmostEqual(motion["omega_f"],omega)
        points = seed_positions(cfg,t0,spin)
        quarter = seed_positions(cfg,t0+np.pi/(2*omega),spin)
        expected = np.column_stack((-points[:,1],points[:,0],points[:,2]))
        np.testing.assert_allclose(quarter,expected,atol=1e-6)
        np.testing.assert_allclose(seed_positions(cfg,t0+2*np.pi/omega,spin),points,atol=1e-6)
        self.assertLess(seed_motion(cfg,t0+1,-spin)["omega_f"],0)
        np.testing.assert_array_equal(seed_positions(cfg,t0+100,0),points)
        fixed = settings()
        fixed["field"]["seed_motion"] = "fixed"
        np.testing.assert_array_equal(seed_positions(fixed,t0+100,spin),points)
        meta = dict(signature="fixture",time=t0,spin=spin)
        self.assertNotEqual(prepared_signature(meta,cfg),prepared_signature(meta,fixed))

    def test_rotated_seeds_reach_the_actual_stream_tracer(self):
        cfg = settings()
        cfg["sampling"]["n"] = 9
        cfg["field"]["maximum_propagation"] = 4.
        cache = np.zeros((9,9,9,4),dtype=np.float32)
        cache[...,0] = 1
        cache[...,3] = 1  # constant native Cartesian Bz: traces stay vertical
        t0 = cfg["field"]["seed_reference_time"]
        omega = seed_motion(cfg,t0,.9375)["omega_f"]
        t = t0+np.pi/(2*omega)
        traces = fields(cache,cfg,t,.9375)
        expected = seed_positions(cfg,t,.9375)
        seed_ids = traces.GetCellData().GetArray("SeedIds")
        for i in range(traces.GetNumberOfLines()):
            first = traces.GetCell(i).GetPointId(0)
            np.testing.assert_allclose(traces.GetPoint(first),expected[int(seed_ids.GetTuple1(i))],atol=1e-6)

    def test_plan_schedule_cli(self):
        with tempfile.TemporaryDirectory() as folder:
            folder = Path(folder)
            for i,time in enumerate([0.,5.,10.]): fixture(folder/f"snap{i}.bin",time)
            plan = folder/"keyframes.json"
            schedule = folder/"schedule.json"
            cli = [sys.executable,str(ROOT/"keyframes.py")]
            subprocess.run(cli+["plan","--glob",str(folder/"*.bin"),"--count","2","--out",str(plan)],check=True,capture_output=True)
            data = read_json(plan)
            for frame,angle,scale in zip(data["keyframes"],[359.,1.],[82.,90.]):
                frame.update(reviewed=True,overrides={"render.azimuth":angle,"render.scale":scale})
            write_json(plan,data)
            subprocess.run(cli+["schedule","--plan",str(plan),"--glob",str(folder/"*.bin"),"--require-reviewed","--out",str(schedule)],check=True,capture_output=True)
            frames = read_json(schedule)["frames"]
            self.assertEqual([f["time"] for f in frames],[0.,5.,10.])
            self.assertAlmostEqual(frames[1]["settings"]["render"]["azimuth"],360.)
            self.assertAlmostEqual(frames[1]["settings"]["render"]["scale"],86.)

    def test_native_reader_component_and_axis_order(self):
        from wongutils.grmhd.athenak import AthenaKSnapshot
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/"affine.bin"
            fixture(path)
            self.assertEqual(header(path)["time"],12.5)
            snap = AthenaKSnapshot(str(path),variable_mapping=["dens","bcc1","bcc2","bcc3"])
            x = np.array([-.75,.25,1.1])
            y = np.array([.2,-.5,.9])
            z = np.array([.1,.8,-.2])
            got = snap.get_primitives_at(x,y,z)
            np.testing.assert_allclose(got[:,0],10+.2*x+.3*y+.1*z,rtol=1e-6)
            np.testing.assert_allclose(got[:,1:],np.tile([1.,.5,-1.],(3,1)),atol=1e-7)

    def test_shortest_angle_and_physical_time(self):
        frames = [dict(time=0.,overrides={"render.azimuth":359.,"render.scale":82.}),
                  dict(time=10.,overrides={"render.azimuth":1.,"render.scale":90.})]
        result = interpolate(settings(),frames,[0.,2.5,5.,10.])
        self.assertAlmostEqual(result[2]["render"]["azimuth"],360.)
        self.assertAlmostEqual(result[1]["render"]["scale"],84.)
        self.assertEqual(result[-1]["render"]["scale"],90.)
        with self.assertRaises(ValueError): interpolate(settings(),frames,[11.])
        with self.assertRaises(ValueError): overrides(settings(),{"field.line_count":80})

    def test_time_selection_is_unique_and_includes_endpoints(self):
        entries = [dict(time=t,source=str(t)) for t in [0,1,2,3,30,60,90,100]]
        selected = choose(entries,5)
        self.assertEqual(selected[0],entries[0])
        self.assertEqual(selected[-1],entries[-1])
        self.assertEqual(len({e["time"] for e in selected}),5)

    def test_trace_identity_survives_missing_branches(self):
        points = vtk.vtkPoints()
        for point in [(0,0,0),(0,0,1),(1,0,0),(1,0,-1)]: points.InsertNextPoint(*point)
        lines = vtk.vtkCellArray()
        for pair in [(0,1),(2,3)]:
            lines.InsertNextCell(2)
            for p in pair: lines.InsertCellPoint(p)
        poly = vtk.vtkPolyData()
        poly.SetPoints(points)
        poly.SetLines(lines)
        seeds = numpy_to_vtk(np.array([0,1],dtype=np.int32),deep=True)
        seeds.SetName("SeedIds")
        poly.GetCellData().AddArray(seeds)
        times = numpy_to_vtk(np.array([0.,1.,0.,-1.]),deep=True)
        times.SetName("IntegrationTime")
        poly.GetPointData().AddArray(times)
        self.assertEqual(stable_curve_ids(poly,4),[0,5])


if __name__=="__main__":
    unittest.main()
