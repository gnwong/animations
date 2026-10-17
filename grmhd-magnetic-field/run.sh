#!/usr/bin/env bash
set -eu
package_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
runtime_root=${GRMHD_RENDER_RUNTIME:-"$package_root/runtime"}
mkdir -p "$runtime_root"/{tmp,cache,config,xdg}
chmod 700 "$runtime_root/xdg"
export TMPDIR="$runtime_root/tmp"
export XDG_CACHE_HOME="$runtime_root/cache"
export XDG_CONFIG_HOME="$runtime_root/config"
export XDG_RUNTIME_DIR="$runtime_root/xdg"
export MESA_SHADER_CACHE_DIR="$runtime_root/cache/mesa"
export __GL_SHADER_DISK_CACHE=0
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}
exec "${PYTHON:-python}" "$package_root/pipeline.py" "$@"
