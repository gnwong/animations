"""Configuration, cache provenance, and fast AthenaK 1.1 header inspection."""
import copy
import hashlib
import importlib.metadata
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def read_json(path):
    return json.loads(Path(path).read_text())


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".tmp")
    temp.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temp.replace(path)


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False).encode()).hexdigest()


def settings(path=None):
    defaults = read_json(ROOT / "settings.json")
    if path is not None:
        def merge(base, overrides, prefix=""):
            for key, value in overrides.items():
                if key not in base:
                    raise ValueError("Unknown setting: " + prefix + key)
                if isinstance(base[key], dict):
                    if not isinstance(value, dict):
                        raise ValueError("Expected an object for " + prefix + key)
                    merge(base[key], value, prefix + key + ".")
                else:
                    base[key] = copy.deepcopy(value)
        merge(defaults, read_json(path))
    validate(defaults)
    return defaults


def validate(cfg):
    import math
    def finite(value):
        if isinstance(value, dict):
            for v in value.values(): finite(v)
        elif isinstance(value, list):
            for v in value: finite(v)
        elif isinstance(value, (int, float)) and not math.isfinite(value):
            raise ValueError("Settings contain a nonfinite number")
    finite(cfg)
    if cfg["schema_version"] != 1:
        raise ValueError("Unsupported settings schema")
    sample, geo, field, render = [cfg[k] for k in ("sampling", "geometry", "field", "render")]
    if sample["coordinates"] != "cartesian-kerr-schild":
        raise ValueError("This pipeline supports Cartesian Kerr-Schild coordinates only")
    for value in (sample["n"], sample["slab"], field["grid_stride"], field["seed_phi_count"],
                  field["maximum_steps"], field["line_count"], field["tube_sides"],
                  geo["inner_iterations"], geo["outer_iterations"], render["width"],
                  render["height"], render["supersampling"], render["threads"], render["peels"]):
        if not isinstance(value, int) or isinstance(value, bool) or value <= 0:
            raise ValueError("Grid sizes, counts, and iterations must be positive integers")
    if sample["n"] < 5 or sample["n"] % 2 != 1 or (sample["n"]-1) % field["grid_stride"]:
        raise ValueError("Use an odd sampling n>=5 with (n-1) divisible by field.grid_stride")
    if len(sample["extents"]) != 3 or min(sample["extents"]) <= 0 or sample["stretch"] <= 0:
        raise ValueError("Invalid sampling extents/stretch")
    if len(geo["levels"]) != 5 or sorted(set(geo["levels"])) != geo["levels"] or min(geo["levels"]) <= 0:
        raise ValueError("The accepted palette requires five increasing positive contour levels")
    for pair in (geo["blend_middle"], geo["blend_outer"], field["fade_z"], field["fade_radius"]):
        if len(pair) != 2 or not 0 <= pair[0] < pair[1]: raise ValueError("Invalid radius/fade pair")
    for key in ("sigma_inner", "sigma_middle", "sigma_outer", "density_gain", "inner_passband", "outer_passband"):
        if geo[key] <= 0: raise ValueError("Invalid geometry." + key)
    if not 0 <= geo["decimation"] < 1: raise ValueError("Invalid decimation")
    if field["line_count"] > 4*len(field["seed_radii"])*field["seed_phi_count"]:
        raise ValueError("Requested more traces than the fixed seed/direction inventory")
    if min(field["seed_radii"]) <= 0 or field["seed_z"] <= 0:
        raise ValueError("Invalid magnetic seed layout")
    if not {"seed_motion","seed_reference_time","seed_phase","seed_omega_fraction"}.issubset(field):
        raise ValueError("This schedule predates moving footpoints; regenerate it with the current package")
    if field["seed_motion"] not in ("fixed","rotate"):
        raise ValueError("field.seed_motion must be fixed or rotate")
    for key in ("seed_reference_time","seed_phase","seed_omega_fraction"):
        if isinstance(field[key],bool) or not isinstance(field[key],(int,float)):
            raise ValueError("Invalid field."+key)
    if not 0 < field["minimum_step"] <= field["initial_step"] <= field["maximum_step"]:
        raise ValueError("Invalid integration step limits")
    for key in ("maximum_propagation", "terminal_speed", "line_radius", "spline_length"):
        if field[key] <= 0: raise ValueError("Invalid field."+key)
    if not 0 <= field["line_opacity"] <= 1: raise ValueError("Invalid line opacity")
    if not 0 < render["inclination"] < 180 or render["scale"] <= 0 or render["warm_radius"] <= 0:
        raise ValueError("Invalid camera inclination/scale or warm radius")
    if render["screen_rotation"] not in (0, 180): raise ValueError("Screen rotation must be 0 or 180")
    if not 0 <= render["fade_start"] < render["fade_end"] or not 0 <= render["feather_start"] < render["feather_end"]:
        raise ValueError("Invalid disk fade/feather radii")
    if len(render["feather_shell_weights"]) != 5 or any(not 0 <= x <= 1 for x in render["feather_shell_weights"]):
        raise ValueError("Invalid feather weights")
    if any(not 0 <= render[k] <= 1 for k in ("outer_silver", "outer_blue", "horizon_contrast", "occlusion_ratio")):
        raise ValueError("Invalid display blend/occlusion ratio")
    for key in ("warm_glow", "warm_crisp", "cool_glow", "gas_glow", "color_gain", "ambient", "diffuse", "key", "rim", "fill", "specular", "opacity", "sharpness"):
        if render[key] < 0: raise ValueError("Invalid render."+key)
    if len(render["background"]) != 3 or any(not 0 <= x <= 255 for x in render["background"]):
        raise ValueError("Invalid background RGB")
    if min(render["effect_reference_width"], render["effect_reference_height"], render["effect_reference_scale"]) <= 0:
        raise ValueError("Invalid effect reference scale")


def header(path):
    """Read metadata without loading meshblock payloads (for snapshot catalogs)."""
    path = Path(path).resolve()
    with path.open("rb") as stream:
        first = stream.readline()
        if not first.startswith(b"Athena ") or first.split()[-1].split(b"=")[-1] != b"1.1":
            raise ValueError("Expected an AthenaK binary format 1.1 snapshot: " + str(path))
        count = int(stream.readline().split(b"=")[-1])
        physical = {}
        for _ in range(count-1):
            key, value = stream.readline().decode().split("=", 1)
            physical[key.strip()] = value.strip()
        nvars = int(stream.readline().split(b"=")[-1])
        variables = stream.readline().decode().split()[1:]
        size = int(stream.readline().split(b"=")[-1])
        raw = stream.read(size)
        end = stream.tell()
        stream.seek(0)
        header_hash = hashlib.sha256(stream.read(end)).hexdigest()
    if len(raw) != size or len(variables) != nvars:
        raise ValueError("Truncated/inconsistent binary header")
    groups, group = {}, None
    for line in raw.decode().splitlines():
        line = line.split("#")[0].strip()
        if line.startswith("<") and line.endswith(">"):
            group = groups.setdefault(line[1:-1], {})
        elif "=" in line and group is not None:
            key, value = line.split("=", 1)
            group[key.strip()] = value.strip()
    spin = float(groups["coord"]["a"])
    if not -1 <= spin <= 1: raise ValueError("Invalid Kerr spin")
    stat = path.stat()
    return dict(source=str(path), size_bytes=stat.st_size, mtime_ns=stat.st_mtime_ns,
                header_sha256=header_hash, format="Athena binary 1.1", time=float(physical["time"]),
                spin=spin, variables=variables, parameters=groups)


def versions():
    import wongutils.grmhd.athenak as reader
    import wongutils.grmhd.meshblocks as meshblocks
    packages = {name: importlib.metadata.version(name) for name in ("numpy", "scipy", "Pillow", "vtk", "tqdm")}
    packages["reader_sha256"] = hashlib.sha256(Path(reader.__file__).read_bytes()).hexdigest()
    packages["meshblocks_sha256"] = hashlib.sha256(Path(meshblocks.__file__).read_bytes()).hexdigest()
    return packages
