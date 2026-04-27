#!/usr/bin/env python3
"""
CAC Model Editor  V2
======================
Integrated GUI for editing, running, and visualizing CAC simulations.

Changes in v1.2:
  1. write_tecplot parser  — reads "write_tecplot <file>" for TecPlot auto-load
  2. Run command now adds  -in  before the input script: ./cac -in script.cac
  3. Variable panel        — scans script for all `variable` definitions & values
  4. Inline command hints  — shows syntax + example from manual when you type a command

Dependencies:
    - tkinter       (sudo apt install python3-tk  or  sudo dnf install python3-tkinter)
    - pytecplot     (pip install pytecplot)   <- optional, only for TecPlot auto-load
"""

import os, sys, re, subprocess, threading, queue, time, glob
from pathlib import Path
from typing import Optional

import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext

try:
    import tecplot as tp
    TECPLOT_AVAILABLE = True
except ImportError:
    TECPLOT_AVAILABLE = False

# ─────────────────────────────────────────────────────────────────────────────
#  PALETTE & FONTS
# ─────────────────────────────────────────────────────────────────────────────
P = {
    "bg":          "#E8E8E8",
    "bg2":         "#D0D0D0",
    "bg3":         "#F5F5F5",
    "border":      "#A8A8A8",
    "accent":      "#1A6FA5",
    "accent2":     "#2E8B57",
    "accent3":     "#C07800",
    "danger":      "#C0392B",
    "text":        "#1A1A1A",
    "text_dim":    "#5A5A5A",
    "text_bright": "#000000",
    "syn_kw":      "#7B2D8B",
    "syn_comment": "#808080",
    "syn_string":  "#2E7D32",
    "syn_number":  "#B85000",
    "syn_section": "#1A6FA5",
    "syn_var":     "#8B4513",
    "hint_bg":     "#DDEEFF",
    "hint_title":  "#1A6FA5",
    "hint_syntax": "#7B2D8B",
    "hint_ex":     "#2E7D32",
    "term_bg":     "#1A1A2E",
    "term_fg":     "#FFD700",
}

FM  = ("Courier", 11)
FMS = ("Courier", 10)
FU  = ("Arial", 10)
FUB = ("Arial", 10, "bold")
FT  = ("Arial", 12, "bold")

# ─────────────────────────────────────────────────────────────────────────────
#  CAC COMMAND HINTS  (built from the CAC manual)
# ─────────────────────────────────────────────────────────────────────────────
CAC_HINTS = {
    "add_atoms": {
        "syntax":  "add_atoms group-ID dim nlayers",
        "args":    "group-ID: group of elements\ndim: +x -x +y -y +z -z (face to add atoms)\nnlayers: number of layers",
        "example": "add_atoms all +x 2",
    },
    "add_etype": {
        "syntax":  "add_etype type_id element_shape apc ulx uly ulz intx inty intz",
        "args":    "type_id: index\nelement_shape: Quad/Tri/Hex/Tet/Oct/Pyr/Wedge\napc: atoms per unit cell\nulx,uly,ulz: size in unit cells\nintx,inty,intz: integration points (use 4 for Hex)",
        "example": "add_etype 1 Hex 2 4 4 4 4 4 4",
    },
    "atom_style": {
        "syntax":  "atom_style atomic",
        "args":    "Only 'atomic' style currently supported",
        "example": "atom_style atomic",
    },
    "atomic_strain": {
        "syntax":  "atomic_strain input_file output_file [keywords values]",
        "args":    "keywords: out/format (plt/szplt/ascii/atom), cutoff, out/config,\n  def/gradient, rotation, strain/tensor, stretch/tensor,\n  wrap, average, no/atom, no/element, region, step",
        "example": "atomic_strain ref.dat strain_out.plt out/format plt cutoff 6",
    },
    "balance": {
        "syntax":  "balance keyword args",
        "args":    "See LAMMPS. Atoms weight=1; elements weight=num integration points",
        "example": "balance 1.1 rcb",
    },
    "boundary": {
        "syntax":  "boundary x y z",
        "args":    "p=periodic  f=fixed  s=shrink-wrapped  m=shrink-wrapped+min",
        "example": "boundary p p p",
    },
    "box": {
        "syntax":  "box keyword value",
        "args":    "tilt: small / large / none",
        "example": "box tilt large",
    },
    "change_box": {
        "syntax":  "change_box group-ID parameter args",
        "args":    "Similar to LAMMPS change_box",
        "example": "change_box all z final 0.0 50.0",
    },
    "charge": {
        "syntax":  "charge atom_type charge_value",
        "args":    "atom_type: atom type index\ncharge_value: charge in electron units",
        "example": "charge 1 1.84  # Sr\ncharge 2 2.36  # Ti\ncharge 3 -1.40 # O",
    },
    "clear": {
        "syntax":  "clear",
        "args":    "Clears all atoms, elements, and settings",
        "example": "clear",
    },
    "coarse_graining": {
        "syntax":  "coarse_graining group_ID cutoff structure_type element_etype atomic_type [keywords]",
        "args":    "structure_type: fcc/bcc/cubic/diamond\nkeywords: compress, migrate, wrap, tol, a1, a2, a3, centroid, seeds",
        "example": "coarse_graining all 2.6 fcc 1 1 tol 0.01 wrap yes &\n    a1 0 0 1 a2 0.8660254 0 0.5 a3 0.2886751 0.81649658 0.5",
    },
    "compute": {
        "syntax":  "compute ID group-ID style args",
        "args":    "styles: centroid/stress/atom, cna/atom, displace/atom, ids/atom,\n  ke, msd, pe, pressure, stress/atom, stress/mech, temp, vacf",
        "example": "compute 1 all pe\ncompute cna all cna/atom adaptive\ncompute disp all displace/atom",
    },
    "create_atoms": {
        "syntax":  "create_atoms type style args [keywords]",
        "args":    "style: box | region region-ID | single x y z\nSee LAMMPS for keywords",
        "example": "create_atoms 1 box",
    },
    "create_box": {
        "syntax":  "create_box N region-ID",
        "args":    "N: number of atom types\nregion-ID: simulation box region",
        "example": "create_box 3 simbox",
    },
    "create_elements": {
        "syntax":  "create_elements etype ctype style args [keywords]",
        "args":    "etype: element type (from add_etype)\nctype: atom/chemical type\nstyle: box | region region-ID | single x y z",
        "example": "create_elements 1 1 box\ncreate_elements 1 1 region myregion",
    },
    "delete_atoms": {
        "syntax":  "delete_atoms style args",
        "args":    "style: group group-ID  |  region region-ID",
        "example": "delete_atoms region delregion",
    },
    "delete_elements": {
        "syntax":  "delete_elements style args",
        "args":    "style: group group-ID  |  region region-ID",
        "example": "delete_elements region delregion",
    },
    "dimension": {
        "syntax":  "dimension N",
        "args":    "N: 2 or 3",
        "example": "dimension 3",
    },
    "disc_elements": {
        "syntax":  "disc_elements region region_ID [migrate yes/no] [compress yes/no] [style atom/element args]",
        "args":    "Discretize elements into atoms or smaller elements\nmigrate: yes/no (default yes)\ncompress: yes/no (default yes)\nstyle atom: fully discrete atoms\nstyle element: split shape (Hex->Hex/HexQuad/Wedge/PyrTet)",
        "example": "disc_elements region disloc_region migrate yes compress yes",
    },
    "displace": {
        "syntax":  "displace group-ID style args [keywords]",
        "args":    "style: move | ramp | rotate | dislocation | multidisl | disp/file\ndislocation args: N ddim dim1 dim2 sdim b angle v theta\nmultidisl args: N ddim dim1 dim2 sdim scale b angle v start_pos",
        "example": "displace all dislocation 1 x 0 0 +y 2.55 90 0.33 0\ndisplace all multidisl 5 x 0 0 +y 100 1.1 2.55 90 0.33 0",
    },
    "dump": {
        "syntax":  "dump ID group-ID style N file [keywords values]",
        "args":    "style: tecplot | atom | structure/atom | stress/mech | cac/data\ntecplot keywords: stress, force, velocity, pe, displace, average, format\nformat: ascii(.dat) / plt / szplt",
        "example": "dump 1 all tecplot 100 output.plt stress yes velocity yes\ndump 2 all atom 100 traj.dat",
    },
    "echo": {
        "syntax":  "echo style",
        "args":    "style: none | screen | log | both",
        "example": "echo both",
    },
    "element_modify": {
        "syntax":  "element_modify keyword values ...",
        "args":    "keywords: id, map, max/ucell(4000), max/gcell(64), max/sub/ucell(500),\n  max/sub/elem(125), nodal/force/style(lumped/consistent), rigid",
        "example": "element_modify max/ucell 8000 max/gcell 64",
    },
    "element_style": {
        "syntax":  "element_style cac max_npe max_apc [wscale nodew edgew surfacew innerw]",
        "args":    "max_npe: max nodes per element\nmax_apc: max atoms per unit cell",
        "example": "element_style cac 8 2",
    },
    "fix": {
        "syntax":  "fix ID group-ID style args",
        "args":    "styles: adaptive, addforce, aveforce, balance, deform, enforce2d,\n  langevin, momentum, move, nve, print, press/berendsen,\n  setforce, temp/press, temp/rescale, viscous",
        "example": "fix 1 all nve\nfix 2 all temp/rescale 1 300 300 10 1.0\nfix 3 top setforce 0 0 NULL",
    },
    "fix adaptive": {
        "syntax":  "fix ID group_ID adaptive N d slipplane region_ID compute_ID perfect_crystal_ID",
        "args":    "N: Nevery (>0)\nd: width between slip planes\nslipplane: xy/yz/xz\ncompute_ID: cna/atom or ids/atom",
        "example": "fix 1 all adaptive 100 2.0 xy slipregion cna 1",
    },
    "fix deform": {
        "syntax":  "fix ID group-ID deform N keyword args ...",
        "args":    "N: every N steps\nkeywords: x/y/z/xy/xz/yz + final/delta/scale/vel/erate/trate/wiggle",
        "example": "fix 1 all deform 1 z erate 1e-4 units box remap v",
    },
    "fix nve": {
        "syntax":  "fix ID group-ID nve",
        "args":    "Constant NVE (microcanonical) integration",
        "example": "fix 1 all nve",
    },
    "fix setforce": {
        "syntax":  "fix ID group-ID setforce fx fy fz [group element/node]",
        "args":    "fx,fy,fz: force components (NULL to skip)",
        "example": "fix 1 top setforce 0 0 NULL",
    },
    "group": {
        "syntax":  "group group-ID style args",
        "args":    "style: region | subtract | union | intersect | atom | element\nregion options: center (default) | allnode | onenode | oneatom",
        "example": "group top region topreg center\ngroup mobile subtract all frozen",
    },
    "jump": {
        "syntax":  "jump file [label]",
        "args":    "Jump to label in file for loop control",
        "example": "jump SELF loop",
    },
    "label": {
        "syntax":  "label ID",
        "args":    "Define a label for jump/loop",
        "example": "label loop",
    },
    "lattice": {
        "syntax":  "lattice style scale [keywords]",
        "args":    "CAC built-in styles: pri/fcc, pri/bcc, pri/hcp, pri/cubic/diamond\nscale: lattice parameter",
        "example": "lattice fcc 3.615\nlattice pri/bcc 2.87",
    },
    "log": {
        "syntax":  "log file [append]",
        "args":    "Redirect thermo output to file",
        "example": "log run.log",
    },
    "mass": {
        "syntax":  "mass type value",
        "args":    "type: atom type index\nvalue: mass in current units",
        "example": "mass 1 87.62  # Sr\nmass 2 47.87  # Ti\nmass 3 16.00  # O",
    },
    "neighbor": {
        "syntax":  "neighbor skin style",
        "args":    "skin: extra distance\nstyle: bin (only option currently)",
        "example": "neighbor 2.0 bin",
    },
    "neigh_modify": {
        "syntax":  "neigh_modify keyword values",
        "args":    "delay, every, check, once, include, exclude, page, one",
        "example": "neigh_modify delay 0 every 1 check yes",
    },
    "newton": {
        "syntax":  "newton flag",
        "args":    "flag: on or off",
        "example": "newton on",
    },
    "next": {
        "syntax":  "next variable",
        "args":    "Increment loop variable",
        "example": "next i",
    },
    "pair_coeff": {
        "syntax":  "pair_coeff I J args",
        "args":    "I, J: atom type indices\nargs depend on pair_style",
        "example": "pair_coeff * * SrTiO3.eam Sr Ti O",
    },
    "pair_style": {
        "syntax":  "pair_style style args",
        "args":    "Supported: buck, coul_wolf, deepmd, eam, eam/alloy, eam/fs,\n  gauss_cut, hybrid, hybrid/overlay, lj/cut, lj/mdf, meam, morse, sw, table",
        "example": "pair_style eam/alloy\npair_style hybrid eam/alloy coul_wolf 0.2 12.0",
    },
    "processors": {
        "syntax":  "processors Nx Ny Nz [style file outfile]",
        "args":    "Use * for auto in any direction",
        "example": "processors * * *",
    },
    "read_data": {
        "syntax":  "read_data file [keywords]",
        "args":    "keywords: add(append/merge), offset toff aoff, shift Sx Sy Sz,\n  extra/atom/types N, extra/element/types N, group groupID,\n  reference(cur/frame | ref/frame | off)",
        "example": "read_data bilayer.data\nread_data grain2.data add append group g2 offset 1 0",
    },
    "region": {
        "syntax":  "region region-ID style args [units box/lattice] [side in/out]",
        "args":    "style:\n  block: xlo xhi ylo yhi zlo zhi\n  sphere: x y z radius\n  cylinder: dim c1 c2 radius lo hi\n  plane: px py pz nx ny nz\n  prism: xlo xhi ylo yhi zlo zhi xy xz yz\nUnits default to box",
        "example": "region simbox block 0 100 0 100 0 50\nregion top block INF INF INF INF 45 50 units box\nregion sph sphere 50 50 25 10.0",
    },
    "reset_ids": {
        "syntax":  "reset_ids",
        "args":    "Re-number atom/element IDs sequentially",
        "example": "reset_ids",
    },
    "reset_timestep": {
        "syntax":  "reset_timestep N",
        "args":    "Set timestep counter to N",
        "example": "reset_timestep 0",
    },
    "run": {
        "syntax":  "run N",
        "args":    "N: number of timesteps (no keyword options in CAC)",
        "example": "run 10000",
    },
    "set_atoms": {
        "syntax":  "set_atoms style ID keyword values",
        "args":    "Operates on atoms only. See LAMMPS set command.",
        "example": "set_atoms group mobile type 2",
    },
    "set_elements": {
        "syntax":  "set_elements style ID keyword values",
        "args":    "Operates on elements only\nstyle: group/ctype — select by group + basis ctype\nkeywords: etype, ctype",
        "example": "set_elements group/ctype noforce 1 ctype 10",
    },
    "thermo": {
        "syntax":  "thermo N",
        "args":    "Print thermodynamic output every N steps",
        "example": "thermo 100",
    },
    "thermo_style": {
        "syntax":  "thermo_style style [args]",
        "args":    "style: one | custom\ncustom keywords: step, elapsed, dt, cpu, elements, atoms, press, pe,\n  vol, lx,ly,lz, pxx,pyy,pzz,pxy,pxz,pyz, c_ID, f_ID, v_name",
        "example": "thermo_style custom step pe press lx ly lz\nthermo_style one",
    },
    "timestep": {
        "syntax":  "timestep dt",
        "args":    "dt: timestep size in time units",
        "example": "timestep 0.001",
    },
    "units": {
        "syntax":  "units style",
        "args":    "style: metal (Bar) | metalgpa (GPa) | si | cgs | real | lj",
        "example": "units metal\nunits metalgpa",
    },
    "variable": {
        "syntax":  "variable name style args",
        "args":    "style: index/loop/world/string/equal/atom/file/python...\nequal: arithmetic using v_, c_, f_ references",
        "example": "variable a equal 3.905\nvariable i loop 10\nvariable str string hello",
    },
    "velocity": {
        "syntax":  "velocity group-ID style args [keywords]",
        "args":    "style: create | set | scale | ramp | zero\ncreate args: temp seed\nkeywords: dist, sum, mom, rot, loop\ngroup: element or node (units=box by default)",
        "example": "velocity all create 300 12345 dist gaussian mom yes rot no",
    },
    "write_data": {
        "syntax":  "write_data file [keywords]",
        "args":    "keywords: noinit yes/no, novelocity yes/no, noelement yes/no, noatom yes/no",
        "example": "write_data final.data\nwrite_data output.data noelement yes",
    },
    "write_tecplot": {
        "syntax":  "write_tecplot file [format ascii/plt/szplt] [noinit yes/no]",
        "args":    "file: output filename (.dat / .plt / .szplt)\nformat: ascii/plt/szplt — default=plt, overridden by file suffix\nnoinit: yes/no (default no)\nRequires TECPLOT package compiled into CAC",
        "example": "write_tecplot bilayer.plt\nwrite_tecplot output.dat format ascii",
    },
    "quit": {
        "syntax":  "quit [N]",
        "args":    "Exit CAC. Optional N = exit code.",
        "example": "quit",
    },
}

# ─────────────────────────────────────────────────────────────────────────────
#  SYNTAX PATTERNS
# ─────────────────────────────────────────────────────────────────────────────
_KW_PAT = r'\b(' + '|'.join(re.escape(k) for k in CAC_HINTS) + r')\b'
SYNTAX_PATTERNS = [
    ("comment",  r'#[^\n]*'),
    ("section",  r'^\s*\[.*?\]'),
    ("string",   r'"[^"]*"'),
    ("number",   r'\b[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?\b'),
    ("variable", r'\$\{?\w+\}?'),
    ("keyword",  _KW_PAT),
]

# ─────────────────────────────────────────────────────────────────────────────
#  PARSE HELPERS
# ─────────────────────────────────────────────────────────────────────────────
def parse_write_tecplot(text):
    # type: (str) -> Optional[str]
    m = re.search(r'^\s*write_tecplot\s+"?([^\s"#]+)"?',
                  text, re.IGNORECASE | re.MULTILINE)
    return m.group(1) if m else None

def parse_write_data(text):
    # type: (str) -> Optional[str]
    m = re.search(r'^\s*write_data\s+"?([^\s"#]+)"?',
                  text, re.IGNORECASE | re.MULTILINE)
    return m.group(1) if m else None

# LAMMPS/CAC built-in runtime variables — skip evaluation
_LAMMPS_BUILTINS = {
    "lx","ly","lz","xlo","xhi","ylo","yhi","zlo","zhi",
    "xy","xz","yz","xlat","ylat","zlat",
    "step","elapsed","elaplong","dt","time","cpu","hcpu",
    "tpcpu","spcpu","cpuremain","hcpuremain","part","timeremain",
    "imbalance","elements","atoms","press","pe","vol",
    "pxx","pyy","pzz","pxy","pxz","pyz",
    "nbuild","ndanger","temp","ke","etotal",
}

import math as _math

# ── safe math context — complete LAMMPS equal-style function set ──────────────
# Source: https://docs.lammps.org/variable.html  (30Mar2026 release)
#
# LAMMPS math operators: + - * / ^ % == != < <= > >= && || |^ !
#   ^ = exponentiation (converted to ** before eval)
#   % = modulo (native Python)
#
# LAMMPS constants: PI, version, on/off/true/false/yes/no
#
# LAMMPS math functions (pure, evaluable):
#   sqrt exp ln log abs sign
#   sin cos tan asin acos atan atan2
#   ceil floor round
#
# LAMMPS runtime functions (need live simulation state → show as "(runtime)"):
#   random normal ternary
#   ramp stagger logfreq logfreq2 logfreq3 stride stride2
#   vdisplace swiggle cwiggle
#   bound                              ← group boundary coordinate
#   count mass charge                  ← group properties
#   xcm vcm fcm                        ← group center of mass / force
#   gyration ke angmom torque inertia omega
#   sum min max ave trap slope sort rsort   ← vector reductions
#   gmask rmask grmask                 ← atom masks
#   next                               ← file variable
#   is_file is_os is_timeout           ← system queries
#   is_available is_active is_defined  ← feature queries
#   label2type is_typelabel            ← type label queries
#   extract_setting                    ← LAMMPS setting query
#   py_varname(...)                    ← python function wrapper

def _lammps_ln(x):
    """LAMMPS ln(x) = natural logarithm."""
    return _math.log(x)

def _lammps_log(x):
    """LAMMPS log(x) = base-10 logarithm."""
    return _math.log10(x)

def _lammps_round(x):
    """LAMMPS round(x) = round half-up, matching C behaviour."""
    return float(_math.floor(x + 0.5))

def _lammps_sign(x):
    """LAMMPS sign(x) = 1 if x>0, -1 if x<0, 0 if x==0."""
    if x > 0: return 1.0
    if x < 0: return -1.0
    return 0.0

# Sentinel: evaluable but depends on simulation runtime state
_RUNTIME_SENTINEL = object()
def _runtime_fn(*a, **kw):
    return _RUNTIME_SENTINEL

_SAFE_CTX = {
    "__builtins__": {},

    # ── constants ─────────────────────────────────────────────────────────────
    "PI":      _math.pi,
    "pi":      _math.pi,
    "e":       _math.e,
    # LAMMPS also supports: version (integer), on/off/true/false/yes/no (1/0)
    "version": 0,          # placeholder — actual value is runtime
    "on":  1.0, "off": 0.0,
    "true": 1.0, "false": 0.0,
    "yes": 1.0, "no":  0.0,

    # ── pure math functions ───────────────────────────────────────────────────
    "sqrt":  _math.sqrt,
    "exp":   _math.exp,
    "ln":    _lammps_ln,       # natural log
    "log":   _lammps_log,      # base-10 log
    "abs":   abs,
    "sign":  _lammps_sign,
    "sin":   _math.sin,
    "cos":   _math.cos,
    "tan":   _math.tan,
    "asin":  _math.asin,
    "acos":  _math.acos,
    "atan":  _math.atan,
    "atan2": _math.atan2,
    "ceil":  _math.ceil,
    "floor": _math.floor,
    "round": _lammps_round,
    "int":   int,              # truncation toward zero

    # ── runtime-only functions → return sentinel ──────────────────────────────
    # Time-dependent / simulation-state-dependent math functions
    "random":    _runtime_fn,  # random(lo,hi,seed)
    "normal":    _runtime_fn,  # normal(mean,sigma,seed)
    "ternary":   _runtime_fn,  # ternary(cond,yes,no) — conditional, may crash
    "ramp":      _runtime_fn,  # ramp(x,y) — timestep-interpolated
    "stagger":   _runtime_fn,  # stagger(x,y)
    "logfreq":   _runtime_fn,  # logfreq(x,y,z)
    "logfreq2":  _runtime_fn,
    "logfreq3":  _runtime_fn,
    "stride":    _runtime_fn,  # stride(x,y,z)
    "stride2":   _runtime_fn,
    "vdisplace": _runtime_fn,  # vdisplace(x,y) — velocity-based displacement
    "swiggle":   _runtime_fn,  # sinusoidal wiggle
    "cwiggle":   _runtime_fn,  # cosinusoidal wiggle

    # Group functions — require atom data
    "bound":     _runtime_fn,  # bound(group,dir[,region]) → box boundary coord
    "count":     _runtime_fn,  # count(group[,region])
    "mass":      _runtime_fn,
    "charge":    _runtime_fn,
    "xcm":       _runtime_fn,  # xcm(group,dim[,region])
    "vcm":       _runtime_fn,
    "fcm":       _runtime_fn,
    "gyration":  _runtime_fn,
    "ke":        _runtime_fn,
    "angmom":    _runtime_fn,
    "torque":    _runtime_fn,
    "inertia":   _runtime_fn,
    "omega":     _runtime_fn,

    # Special / vector reduction functions
    "sum":       _runtime_fn,
    "min":       _runtime_fn,
    "max":       _runtime_fn,
    "ave":       _runtime_fn,
    "trap":      _runtime_fn,
    "slope":     _runtime_fn,
    "sort":      _runtime_fn,
    "rsort":     _runtime_fn,
    "gmask":     _runtime_fn,
    "rmask":     _runtime_fn,
    "grmask":    _runtime_fn,
    "next":      _runtime_fn,

    # System / feature query functions
    "is_file":        _runtime_fn,
    "is_os":          _runtime_fn,
    "is_timeout":     _runtime_fn,
    "is_available":   _runtime_fn,
    "is_active":      _runtime_fn,
    "is_defined":     _runtime_fn,
    "label2type":     _runtime_fn,
    "is_typelabel":   _runtime_fn,
    "extract_setting":_runtime_fn,
}

def _all_refs(expr):
    """Return set of all lower-case variable names referenced in an expression.
    Catches v_name, ${name}, $name, and bare LAMMPS builtin names."""
    refs = set()
    for m in re.finditer(r'\bv_(\w+)\b', expr, re.IGNORECASE):
        refs.add(m.group(1).lower())
    for m in re.finditer(r'\$\{(\w+)\}', expr):
        refs.add(m.group(1).lower())
    for m in re.finditer(r'\$(\w+)', expr):
        refs.add(m.group(1).lower())
    # Also catch bare LAMMPS builtin names used directly (e.g. xhi, xlo, lx)
    for m in re.finditer(r'\b(\w+)\b', expr):
        if m.group(1).lower() in _LAMMPS_BUILTINS:
            refs.add(m.group(1).lower())
    return refs

def _substitute(expr, env):
    """Replace all v_name, ${name}, $name with their numeric values from env.
    Returns (substituted_expr, set_of_still_unresolved_names)."""
    unresolved = set()
    def _rep(vn_lower):
        if vn_lower in env:
            return "(%s)" % env[vn_lower]
        unresolved.add(vn_lower)
        return "__UNRESOLVED__"
    out = re.sub(r'\bv_(\w+)\b',  lambda m: _rep(m.group(1).lower()), expr)
    out = re.sub(r'\$\{(\w+)\}',  lambda m: _rep(m.group(1).lower()), out)
    out = re.sub(r'\$(\w+)',        lambda m: _rep(m.group(1).lower()), out)
    return out, unresolved

# Set of runtime function names for fast lookup
_RUNTIME_FN_NAMES = {
    k for k, v in _SAFE_CTX.items()
    if v is _runtime_fn and k != "__builtins__"
}

def _contains_runtime_fn(expr):
    """Return True if expr contains a call to a runtime-only LAMMPS function."""
    for m in re.finditer(r'\b(\w+)\s*\(', expr):
        if m.group(1).lower() in _RUNTIME_FN_NAMES:
            return True
    return False

def _safe_eval(expr_raw, env):
    """Evaluate a LAMMPS equal-style expression.
    Returns (value_string_or_'(runtime)', ok_bool).

    Handles:
    - Surrounding quotes stripped  (e.g. variable x equal "PI/180.0")
    - LAMMPS ^ exponent converted to Python **
    - v_name / ${name} / $name variable substitution
    - Full LAMMPS math function set (sqrt floor ceil ln log sin cos tan
      asin acos atan atan2 abs sign round exp ...)
    - Runtime-only group/special functions (bound count xcm ramp ...) -> '(runtime)'
    - % modulo operator (native Python)
    """
    # Strip surrounding quotes
    expr = expr_raw.strip().strip('"\'')

    # Early-exit: if expression calls any runtime-only function, it is runtime
    if _contains_runtime_fn(expr):
        return "(runtime)", True

    # Convert LAMMPS exponent ^ to Python **
    expr = re.sub(r'(?<![\*])\^(?![\*])', '**', expr)

    # Substitute variable references (v_name, ${name}, $name only)
    expr2, unresolved = _substitute(expr, env)
    if unresolved:
        return None, False  # still has unresolved variable references

    try:
        val = eval(expr2, _SAFE_CTX)
        # Paranoia: sentinel shouldn't appear here (caught above), but just in case
        if val is _RUNTIME_SENTINEL:
            return "(runtime)", True
        if isinstance(val, (int, float)):
            if isinstance(val, float) and val == int(val) and abs(val) < 1e15:
                return str(int(val)), True
            return ("%.6g" % val), True
        return str(val), True
    except Exception:
        return None, False

def parse_variables(text):
    # type: (str) -> list
    """Return list of (name, style, raw_expr, evaluated_value).
    evaluated_value:
      - numeric string  : successfully computed from script context
      - "(runtime)"     : depends on LAMMPS box/thermo builtins (known at run time only)
      - None            : genuinely unresolvable

    Uses multi-pass evaluation so forward-referenced variables (defined later
    in the script) are still resolved — e.g. costheta uses v_theta which is
    declared several lines below it.
    """
    raw = []
    for m in re.finditer(
            r'^\s*variable\s+(\S+)\s+(\S+)\s*(.*?)(?:\s*#.*)?$',
            text, re.MULTILINE):
        raw.append((m.group(1).strip(), m.group(2).strip().lower(), m.group(3).strip()))

    # ── Pass 1: collect all non-equal styles and simple equal values ─────────
    env          = {}           # name -> float, for resolved variables
    runtime_vars = set(_LAMMPS_BUILTINS)
    evaluated_map = {}          # name -> result string (or None / "(runtime)")

    def _process_one(name, style, expr):
        """Try to evaluate one variable. Returns evaluated string or None."""
        nm = name.lower()
        if style == "equal":
            refs = _all_refs(expr)
            if refs & runtime_vars:
                runtime_vars.add(nm)
                return "(runtime)"
            val_str, ok = _safe_eval(expr, env)
            if ok:
                if val_str == "(runtime)":
                    runtime_vars.add(nm)
                    return "(runtime)"
                try:
                    env[nm] = float(val_str)
                except (ValueError, TypeError):
                    pass
                return val_str
            return None   # still has unresolved refs — retry later

        elif style in ("index", "loop"):
            parts = expr.split()
            if parts:
                try: env[nm] = float(parts[0])
                except (ValueError, TypeError): pass
                return parts[0]
            return None

        elif style == "string":
            return expr.strip('"\'')

        elif style in ("file","atomfile","python","universe",
                       "world","format","uloop","getenv"):
            runtime_vars.add(nm)
            return "(runtime)"

        return None

    # First pass
    for name, style, expr in raw:
        evaluated_map[name.lower()] = _process_one(name, style, expr)

    # ── Passes 2-N: retry failed equal-style vars until no more progress ────
    # This resolves forward references (e.g. costheta uses v_theta defined later)
    MAX_PASSES = 6
    for _pass in range(MAX_PASSES):
        progress = False
        for name, style, expr in raw:
            nm = name.lower()
            if style == "equal" and evaluated_map.get(nm) is None:
                result = _process_one(name, style, expr)
                if result is not None:
                    evaluated_map[nm] = result
                    progress = True
        if not progress:
            break   # nothing new resolved — stop

    # ── Build final result list in declaration order ─────────────────────────
    results = []
    for name, style, expr in raw:
        results.append((name, style, expr, evaluated_map.get(name.lower())))
    return results

def current_line_keyword(text_widget):
    # type: (tk.Text) -> Optional[str]
    idx  = text_widget.index("insert linestart")
    line = text_widget.get(idx, idx + " lineend").strip()
    if not line or line.startswith('#'):
        return None
    tokens = line.split()
    if not tokens:
        return None
    if len(tokens) >= 2:
        two = tokens[0].lower() + " " + tokens[1].lower()
        if two in CAC_HINTS:
            return two
    return tokens[0].lower() if tokens[0].lower() in CAC_HINTS else None

# ─────────────────────────────────────────────────────────────────────────────
#  SYNTAX HIGHLIGHTER
# ─────────────────────────────────────────────────────────────────────────────
class SyntaxHighlighter(object):
    TAG_MAP = {
        "comment":  P["syn_comment"],
        "section":  P["syn_section"],
        "string":   P["syn_string"],
        "number":   P["syn_number"],
        "variable": P["syn_var"],
        "keyword":  P["syn_kw"],
    }
    def __init__(self, widget):
        self.w = widget
        for tag, color in self.TAG_MAP.items():
            self.w.tag_configure(tag, foreground=color)
        self.w.tag_configure("keyword",
            foreground=P["syn_kw"], font=(FM[0], FM[1], "bold"))
        self.w.tag_configure("section",
            foreground=P["syn_section"], font=(FM[0], FM[1], "bold"))
        self._aid = None

    def schedule(self, _=None):
        if self._aid:
            self.w.after_cancel(self._aid)
        self._aid = self.w.after(150, self.run)

    def run(self):
        content = self.w.get("1.0", "end-1c")
        for tag in self.TAG_MAP:
            self.w.tag_remove(tag, "1.0", "end")
        for tag, pat in SYNTAX_PATTERNS:
            for m in re.finditer(pat, content, re.MULTILINE):
                s = self._ofs(content, m.start())
                e = self._ofs(content, m.end())
                self.w.tag_add(tag, s, e)

    def _ofs(self, text, offset):
        before = text[:offset]
        ln  = before.count('\n') + 1
        col = offset - before.rfind('\n') - 1
        return "%d.%d" % (ln, col)

# ─────────────────────────────────────────────────────────────────────────────
#  LINE NUMBERS
# ─────────────────────────────────────────────────────────────────────────────
class LineNumbers(tk.Canvas):
    def __init__(self, parent, tw, **kw):
        super(LineNumbers, self).__init__(parent, **kw)
        self.tw = tw
        self.configure(bg=P["bg2"], highlightthickness=0, width=46)
        for ev in ("<Configure>","<<Modified>>","<KeyRelease>",
                   "<Button-4>","<Button-5>","<MouseWheel>"):
            self.tw.bind(ev, self._draw, add="+")

    def _draw(self, _=None):
        self.delete("all")
        i = self.tw.index("@0,0")
        while True:
            dl = self.tw.dlineinfo(i)
            if dl is None:
                break
            self.create_text(38, dl[1], anchor="ne",
                text=str(i).split(".")[0],
                fill=P["text_dim"], font=FMS)
            ni = self.tw.index("%s+1line" % i)
            if ni == self.tw.index("%s+0lines" % ni):
                break
            i = ni

# ─────────────────────────────────────────────────────────────────────────────
#  HINT PANEL
# ─────────────────────────────────────────────────────────────────────────────
def _parse_manual_hints(path):
    """Parse Manual.md and return dict of {command: {syntax,args,example}}.
    Falls back to built-in CAC_HINTS for any command not found in the file."""
    if not path or not os.path.isfile(path):
        return {}
    try:
        with open(path, encoding="utf-8", errors="replace") as f:
            text = f.read()
    except Exception:
        return {}

    hints = {}
    # Split on ### headings
    sections = re.split(r'(?m)^###\s+', text)
    for sec in sections[1:]:
        lines = sec.strip().splitlines()
        if not lines:
            continue
        cmd = lines[0].strip().lower()
        body = "\n".join(lines[1:])
        # Extract syntax block (indented lines after "syntax:")
        syntax_m = re.search(r'(?:command syntax|syntax)[:\s]*\n((?:[ \t]+[^\n]+\n?)+)',
                              body, re.IGNORECASE)
        syntax_str = ""
        if syntax_m:
            syntax_str = "\n".join(l.strip() for l in
                                    syntax_m.group(1).strip().splitlines()
                                    if l.strip())
        # Extract example block
        ex_m = re.search(r'(?:example[s]?)[:\s]*\n((?:[ \t]+[^\n]+\n?)+)',
                         body, re.IGNORECASE)
        ex_str = ""
        if ex_m:
            ex_str = "\n".join(l.strip() for l in
                                ex_m.group(1).strip().splitlines()
                                if l.strip())
        # Body text (strip markdown links/tags, first 400 chars)
        clean = re.sub(r'<[^>]+>', '', body)
        clean = re.sub(r'\[.*?\]\(.*?\)', '', clean).strip()
        # trim to first ~400 chars sensibly
        if len(clean) > 400:
            clean = clean[:400].rsplit("\n", 1)[0] + "..."

        if cmd:
            hints[cmd] = {
                "syntax":  syntax_str or "(see manual)",
                "args":    clean,
                "example": ex_str or "",
                "_from_file": True,
            }
    return hints


class HintPanel(tk.Frame):
    def __init__(self, parent, manual_path_var, **kw):
        super(HintPanel, self).__init__(parent, bg=P["hint_bg"], **kw)
        self._cur = None
        self._manual_path_var = manual_path_var   # tk.StringVar from app
        self._manual_cache = {}    # path -> parsed hints dict
        self._last_mtime  = 0

        hdr = tk.Frame(self, bg=P["hint_bg"])
        hdr.pack(fill="x", padx=8, pady=(4,0))
        self.title_lbl = tk.Label(hdr, text="",
            fg=P["hint_title"], bg=P["hint_bg"],
            font=(FM[0], 11, "bold"), anchor="w")
        self.title_lbl.pack(side="left")
        self.src_lbl = tk.Label(hdr, text="COMMAND REFERENCE",
            fg=P["text_dim"], bg=P["hint_bg"], font=FU)
        self.src_lbl.pack(side="right")

        self.txt = tk.Text(self,
            bg=P["hint_bg"], fg=P["text"], font=FMS,
            relief="flat", bd=0, height=5,
            wrap="word", padx=10, pady=4, state="disabled")
        self.txt.tag_configure("label",  foreground=P["text_dim"],
                               font=(FMS[0], FMS[1], "bold"))
        self.txt.tag_configure("syntax", foreground=P["hint_syntax"])
        self.txt.tag_configure("ex",     foreground=P["hint_ex"])
        self.txt.tag_configure("body",   foreground=P["text"])
        self.txt.tag_configure("file",   foreground=P["accent3"])
        self.txt.pack(fill="both", expand=True)

    def _get_hints(self):
        """Return merged hint dict: file hints override built-ins if manual is set."""
        path = self._manual_path_var.get().strip()
        if path and os.path.isfile(path):
            try:
                mtime = os.path.getmtime(path)
            except OSError:
                mtime = 0
            if path not in self._manual_cache or mtime != self._last_mtime:
                self._manual_cache = {path: _parse_manual_hints(path)}
                self._last_mtime   = mtime
                self.src_lbl.config(text="MANUAL: %s" % os.path.basename(path),
                                    fg=P["accent3"])
            file_hints = self._manual_cache.get(path, {})
            # merge: file takes priority, fall back to built-in
            merged = dict(CAC_HINTS)
            merged.update(file_hints)
            return merged
        else:
            self.src_lbl.config(text="COMMAND REFERENCE (built-in)",
                                fg=P["text_dim"])
            return CAC_HINTS

    def show(self, key):
        # type: (Optional[str]) -> None
        hints = self._get_hints()
        # force refresh if same key but manual may have changed
        self._cur = None
        if key == self._cur:
            return
        self._cur = key
        self.txt.config(state="normal")
        self.txt.delete("1.0", "end")
        if key and key in hints:
            h = hints[key]
            from_file = h.get("_from_file", False)
            self.title_lbl.config(text="  " + key)
            self.txt.insert("end", "Syntax:  ", "label")
            self.txt.insert("end", h["syntax"] + "\n", "syntax")
            self.txt.insert("end", "Args:    ", "label")
            self.txt.insert("end", h.get("args","") + "\n", "body")
            if h.get("example"):
                self.txt.insert("end", "Example: ", "label")
                self.txt.insert("end", h["example"] + "\n",
                                "file" if from_file else "ex")
        else:
            self.title_lbl.config(text="  —")
            self.txt.insert("end",
                "Start typing a CAC command on any line to see its syntax here.",
                "body")
        self.txt.config(state="disabled")

# ─────────────────────────────────────────────────────────────────────────────
#  VARIABLE PANEL
# ─────────────────────────────────────────────────────────────────────────────
class VariablePanel(tk.Frame):
    def __init__(self, parent, **kw):
        super(VariablePanel, self).__init__(parent, bg=P["bg2"], **kw)

        hdr = tk.Frame(self, bg=P["bg2"], pady=4)
        hdr.pack(fill="x", padx=8)
        tk.Label(hdr, text="VARIABLES", fg=P["text_dim"],
                 bg=P["bg2"], font=FUB).pack(side="left")
        self.count_lbl = tk.Label(hdr, text="",
                                  fg=P["text_dim"], bg=P["bg2"], font=FU)
        self.count_lbl.pack(side="right")

        cols = ("Name", "Style", "Expression", "Value")
        self.tree = ttk.Treeview(self, columns=cols, show="headings", height=10)
        for c, w in zip(cols, (110, 70, 200, 120)):
            self.tree.heading(c, text=c)
            self.tree.column(c, width=w, anchor="w")
        # tag colours
        self.tree.tag_configure("resolved",  foreground=P["accent2"])
        self.tree.tag_configure("runtime",   foreground=P["accent3"])
        self.tree.tag_configure("unresolved",foreground=P["text_dim"])

        vsb = ttk.Scrollbar(self, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=vsb.set)
        vsb.pack(side="right", fill="y")
        self.tree.pack(fill="both", expand=True)

    def refresh(self, script_text):
        # type: (str) -> None
        self.tree.delete(*self.tree.get_children())
        vars_ = parse_variables(script_text)
        for name, style, expr, evaluated in vars_:
            if evaluated is None:
                disp  = "—"
                tag   = "unresolved"
            elif evaluated == "(runtime)":
                disp  = "(runtime)"
                tag   = "runtime"
            else:
                disp  = evaluated
                tag   = "resolved"
            self.tree.insert("", "end",
                values=(name, style, expr, disp), tags=(tag,))
        n = len(vars_)
        self.count_lbl.config(
            text="%d variable%s" % (n, "s" if n != 1 else ""))

# ─────────────────────────────────────────────────────────────────────────────
#  SCRIPT EDITOR TAB
# ─────────────────────────────────────────────────────────────────────────────
class ScriptEditorTab(ttk.Frame):
    def __init__(self, parent, app, **kw):
        super(ScriptEditorTab, self).__init__(parent, **kw)
        self.app = app
        self.current_file = None
        self._unsaved = False
        self._build()

    def _build(self):
        # toolbar
        tb = tk.Frame(self, bg=P["bg2"], pady=6, padx=8)
        tb.pack(fill="x")

        def btn(text, cmd, color=P["accent"]):
            tk.Button(tb, text=text, command=cmd,
                bg=P["bg3"], fg=color, relief="flat",
                activebackground=P["border"], font=FUB,
                padx=10, pady=4, cursor="hand2", bd=0
            ).pack(side="left", padx=3)

        btn("New",    self.new_file)
        btn("Open",   self.open_file)
        btn("Save",   self.save_file)
        btn("Save As",self.save_as)
        self.fname_lbl = tk.Label(tb, text="(no file)",
            fg=P["text_dim"], bg=P["bg2"], font=FU)
        self.fname_lbl.pack(side="left", padx=12)
        self.out_lbl = tk.Label(tb, text="",
            fg=P["accent2"], bg=P["bg2"], font=FU)
        self.out_lbl.pack(side="right", padx=8)
        def _set_manual():
            p = filedialog.askopenfilename(
                title="Select CAC Manual",
                filetypes=[("Markdown","*.md *.MD"),("All","*.*")])
            if p:
                self.app.manual_path_var.set(p)
                manual_btn.config(
                    text="Manual: %s" % os.path.basename(p),
                    fg=P["accent2"])
        manual_btn = tk.Button(tb, text="Manual...",
            command=_set_manual,
            bg=P["bg3"], fg=P["text_dim"], relief="flat",
            activebackground=P["border"], font=("Arial", 9),
            padx=6, pady=4, cursor="hand2", bd=0)
        manual_btn.pack(side="right", padx=(0,4))
        tk.Checkbutton(tb,
            text="Default",
            variable=self.app.save_manual_var,
            command=self.app._save_default_manual,
            bg=P["bg2"], fg=P["text_dim"],
            selectcolor=P["bg3"], activebackground=P["bg2"],
            font=("Arial", 8), cursor="hand2"
        ).pack(side="right", padx=(0,2))

        # horizontal split: editor | variable panel
        pane = tk.PanedWindow(self, orient="horizontal",
                              bg=P["bg"], sashwidth=4, sashrelief="flat")
        pane.pack(fill="both", expand=True)

        # left: editor + hint below
        left = tk.Frame(pane, bg=P["bg"])
        pane.add(left, minsize=400)

        edit_row = tk.Frame(left, bg=P["bg"])
        edit_row.pack(fill="both", expand=True)

        self.text = tk.Text(edit_row,
            bg=P["bg3"], fg=P["text"],
            insertbackground=P["accent"],
            selectbackground=P["border"],
            font=FM, relief="flat", bd=0,
            tabs="    ", undo=True, wrap="none",
            padx=10, pady=8)

        self.lnums = LineNumbers(edit_row, self.text)
        self.lnums.pack(side="left", fill="y")

        vs = ttk.Scrollbar(edit_row, orient="vertical",   command=self.text.yview)
        hs = ttk.Scrollbar(left,     orient="horizontal", command=self.text.xview)
        self.text.configure(yscrollcommand=vs.set, xscrollcommand=hs.set)
        vs.pack(side="right", fill="y")
        hs.pack(side="bottom", fill="x")
        self.text.pack(side="left", fill="both", expand=True)

        # hint panel at the bottom of editor
        self.hint = HintPanel(left, self.app.manual_path_var)
        self.hint.pack(fill="x", side="bottom")
        self.hint.show(None)

        # right: variable panel
        right = tk.Frame(pane, bg=P["bg2"])
        pane.add(right, minsize=220)
        self.varpanel = VariablePanel(right)
        self.varpanel.pack(fill="both", expand=True)

        # status bar
        self._sb = tk.Label(self, text="Ln 1, Col 1",
            bg=P["bg2"], fg=P["text_dim"], font=FU, anchor="w", padx=8)
        self._sb.pack(fill="x", side="bottom")

        # bindings
        self.hl = SyntaxHighlighter(self.text)
        self.text.bind("<KeyRelease>", self._on_key)
        self.text.bind("<<Modified>>", self._on_mod)
        self.text.bind("<ButtonRelease>", self._update_status)

    # events
    def _on_key(self, _=None):
        self.hl.schedule()
        self._update_status()
        self._update_out_label()
        self._update_hint()
        self._refresh_vars()

    def _on_mod(self, _=None):
        self._unsaved = True
        self.text.edit_modified(False)

    def _update_status(self, _=None):
        idx = self.text.index("insert")
        ln, col = idx.split(".")
        self._sb.config(text="Ln %s, Col %d" % (ln, int(col)+1))

    def _update_out_label(self):
        txt = self.text.get("1.0", "end-1c")
        tp_f = parse_write_tecplot(txt)
        wd_f = parse_write_data(txt)
        if tp_f:
            self.out_lbl.config(text="-> TecPlot out: %s" % tp_f)
            self.app.detected_output_file = tp_f
        elif wd_f:
            self.out_lbl.config(text="-> data out: %s" % wd_f)
            self.app.detected_output_file = wd_f
        else:
            self.out_lbl.config(text="")
            self.app.detected_output_file = None

    def _update_hint(self):
        self.hint.show(current_line_keyword(self.text))

    def _refresh_vars(self):
        self.varpanel.refresh(self.text.get("1.0", "end-1c"))

    # file ops
    def new_file(self):
        if self._unsaved:
            if not messagebox.askyesno("Unsaved changes","Discard?"):
                return
        self.text.delete("1.0","end")
        self.current_file = None
        self.fname_lbl.config(text="(no file)")
        self._unsaved = False

    def open_file(self, path=None):
        if path is None:
            path = filedialog.askopenfilename(
                title="Open CAC Script",
                filetypes=[("All files","*.*")])
        if not path: return
        try:
            with open(path) as f:
                content = f.read()
            self.text.delete("1.0","end")
            self.text.insert("1.0", content)
            self.current_file = path
            self.fname_lbl.config(text=os.path.basename(path))
            self.app.workdir_var.set(os.path.dirname(os.path.abspath(path)))
            self.hl.run()
            self._update_out_label()
            self._refresh_vars()
            self._unsaved = False
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def save_file(self):
        if self.current_file: self._write(self.current_file)
        else: self.save_as()

    def save_as(self):
        path = filedialog.asksaveasfilename(
            defaultextension=".cac",
            filetypes=[("All files","*.*")])
        if path:
            self._write(path)
            self.current_file = path
            self.fname_lbl.config(text=os.path.basename(path))

    def _write(self, path):
        try:
            with open(path,"w") as f:
                f.write(self.text.get("1.0","end-1c"))
            self._unsaved = False
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def get_script_text(self):  return self.text.get("1.0","end-1c")
    def get_current_file(self): return self.current_file

# ─────────────────────────────────────────────────────────────────────────────
#  RUN TAB
# ─────────────────────────────────────────────────────────────────────────────
class RunTab(ttk.Frame):
    def __init__(self, parent, app, **kw):
        super(RunTab, self).__init__(parent, **kw)
        self.app = app
        self._process = None
        self._log_q   = queue.Queue()
        self._running = False
        self._build()
        self._poll()

    def _build(self):
        # ── Shared: Working Dir + Input Script (always enabled) ───────────────
        top = tk.LabelFrame(self, text="  INPUT  ",
            bg=P["bg2"], fg=P["text_bright"], font=FUB,
            padx=10, pady=8, relief="groove", bd=2)
        top.pack(fill="x", padx=8, pady=(8,4))
        top.columnconfigure(1, weight=1)

        def trow(label, var, r, browse=None):
            tk.Label(top, text=label, fg=P["text_dim"], bg=P["bg2"],
                font=FU, anchor="w", width=14
            ).grid(row=r, column=0, sticky="w", pady=3)
            tk.Entry(top, textvariable=var, bg=P["bg3"], fg=P["text"],
                insertbackground=P["accent"], relief="flat", font=FMS, bd=3
            ).grid(row=r, column=1, sticky="ew", padx=4, pady=3)
            if browse:
                tk.Button(top, text="...", command=browse,
                    bg=P["bg3"], fg=P["accent"], relief="flat",
                    font=FU, bd=0, padx=6, cursor="hand2"
                ).grid(row=r, column=2, padx=2)

        trow("Working Dir",  self.app.workdir_var, 0, self._br_dir)
        trow("Input Script", self.app.script_var,  1, self._br_script)

        # ── Run method selector ───────────────────────────────────────────────
        sel = tk.Frame(self, bg=P["bg2"], padx=12, pady=6)
        sel.pack(fill="x", padx=8)
        tk.Label(sel, text="Run method:", fg=P["text_dim"],
            bg=P["bg2"], font=FUB).pack(side="left", padx=(0,12))

        self._run_mode = tk.StringVar(value="regular")
        for val, label in [("regular", "Regular"), ("sif", "Singularity")]:
            tk.Radiobutton(sel, text=label,
                variable=self._run_mode, value=val,
                command=self._toggle_run_mode,
                bg=P["bg2"], fg=P["text"],
                selectcolor=P["bg3"], activebackground=P["bg2"],
                font=FUB, cursor="hand2"
            ).pack(side="left", padx=6)

        # ── Split panels ──────────────────────────────────────────────────────
        panels = tk.Frame(self, bg=P["bg"])
        panels.pack(fill="x", padx=8, pady=(4,4))
        panels.columnconfigure(0, weight=1)
        panels.columnconfigure(1, weight=1)

        # helper: create entry row inside a given frame, return (label,entry,[btn])
        def prow(parent, bg, label, var, r, browse=None):
            lbl = tk.Label(parent, text=label, fg=P["text_dim"], bg=bg,
                font=FU, anchor="w", width=14)
            lbl.grid(row=r, column=0, sticky="w", pady=2)
            ent = tk.Entry(parent, textvariable=var, bg=P["bg3"], fg=P["text"],
                insertbackground=P["accent"], relief="flat", font=FMS, bd=3)
            ent.grid(row=r, column=1, sticky="ew", padx=4, pady=2)
            widgets = [lbl, ent]
            if browse:
                btn = tk.Button(parent, text="...", command=browse,
                    bg=P["bg3"], fg=P["accent"], relief="flat",
                    font=FU, bd=0, padx=5, cursor="hand2")
                btn.grid(row=r, column=2, padx=2)
                widgets.append(btn)
            return widgets

        # ── LEFT: Regular ─────────────────────────────────────────────────────
        self._reg_frame = tk.LabelFrame(panels, text="  REGULAR RUN  ",
            bg=P["bg2"], fg=P["accent"], font=FUB,
            padx=10, pady=8, relief="groove", bd=2)
        self._reg_frame.grid(row=0, column=0, sticky="nsew", padx=(0,4))
        self._reg_frame.columnconfigure(1, weight=1)

        self._reg_widgets = []
        self._reg_widgets += prow(self._reg_frame, P["bg2"],
            "CAC Executable", self.app.cac_exe_var, 0, self._br_exe)
        cb_exe = tk.Checkbutton(self._reg_frame, text="Set as default",
            variable=self.app.save_exe_var,
            command=self.app._save_default_exe,
            bg=P["bg2"], fg=P["text_dim"], selectcolor=P["bg3"],
            activebackground=P["bg2"], font=("Arial",9), cursor="hand2")
        cb_exe.grid(row=1, column=1, sticky="w", pady=(0,2))
        self._reg_widgets.append(cb_exe)
        self._reg_widgets += prow(self._reg_frame, P["bg2"],
            "MPI ranks", self.app.mpi_var, 2)
        self._reg_widgets += prow(self._reg_frame, P["bg2"],
            "module load", self.app.module_var, 3)
        tk.Label(self._reg_frame,
            text="(e.g. CAC openmpi — blank=skip)",
            fg=P["text_dim"], bg=P["bg2"],
            font=("Arial",8), anchor="w"
        ).grid(row=4, column=1, sticky="w")

        self.run_btn = tk.Button(self._reg_frame, text="RUN",
            command=self._run_local,
            bg=P["accent2"], fg=P["bg"],
            font=(FUB[0],11,"bold"),
            relief="flat", bd=0, padx=18, pady=6, cursor="hand2")
        self.run_btn.grid(row=5, column=0, columnspan=2,
                          sticky="w", pady=(8,2))
        self._reg_widgets.append(self.run_btn)

        # ── RIGHT: Singularity ────────────────────────────────────────────────
        self._sif_frame = tk.LabelFrame(panels, text="  SINGULARITY RUN  ",
            bg=P["bg2"], fg=P["accent3"], font=FUB,
            padx=10, pady=8, relief="groove", bd=2)
        self._sif_frame.grid(row=0, column=1, sticky="nsew", padx=(4,0))
        self._sif_frame.columnconfigure(1, weight=1)

        self._sif_widgets = []
        self._sif_widgets += prow(self._sif_frame, P["bg2"],
            ".sif path", self.app.sif_var, 0, self._br_sif)
        cb_sif = tk.Checkbutton(self._sif_frame, text="Set as default",
            variable=self.app.save_sif_var,
            command=self.app._save_default_sif,
            bg=P["bg2"], fg=P["text_dim"], selectcolor=P["bg3"],
            activebackground=P["bg2"], font=("Arial",9), cursor="hand2")
        cb_sif.grid(row=1, column=1, sticky="w", pady=(0,2))
        self._sif_widgets.append(cb_sif)
        self._sif_widgets += prow(self._sif_frame, P["bg2"],
            "Bind dir", self.app.sif_bind_var, 2, self._br_sif_bind)
        tk.Label(self._sif_frame,
            text="(outermost bind dir, e.g. /data1)",
            fg=P["text_dim"], bg=P["bg2"],
            font=("Arial",8), anchor="w"
        ).grid(row=3, column=1, sticky="w")
        self._sif_widgets += prow(self._sif_frame, P["bg2"],
            "MPI ranks", self.app.sif_mpi_var, 4)
        self._sif_widgets += prow(self._sif_frame, P["bg2"],
            "module load", self.app.sif_module_var, 5)
        tk.Label(self._sif_frame,
            text="(leave blank to skip)",
            fg=P["text_dim"], bg=P["bg2"],
            font=("Arial",8), anchor="w"
        ).grid(row=6, column=1, sticky="w")

        self.sif_run_btn = tk.Button(self._sif_frame, text="RUN SIF",
            command=self._run_singularity,
            bg=P["accent3"], fg=P["bg"],
            font=(FUB[0],11,"bold"),
            relief="flat", bd=0, padx=18, pady=6, cursor="hand2")
        self.sif_run_btn.grid(row=7, column=0, columnspan=2,
                               sticky="w", pady=(8,2))
        self._sif_widgets.append(self.sif_run_btn)

        # Initial state: regular active, sif greyed
        self._toggle_run_mode()

        # ── Kill / status / progress ──────────────────────────────────────────
        bf = tk.Frame(self, bg=P["bg"], pady=4)
        bf.pack(fill="x", padx=8)
        self.kill_btn = tk.Button(bf, text="KILL",
            command=self._kill,
            bg=P["danger"], fg=P["text_bright"],
            font=FUB, relief="flat", bd=0,
            padx=12, pady=6, cursor="hand2", state="disabled")
        self.kill_btn.pack(side="left")
        self.status_lbl = tk.Label(bf, text="Idle",
            fg=P["text_dim"], bg=P["bg"], font=FU)
        self.status_lbl.pack(side="left", padx=12)
        self.progress = ttk.Progressbar(bf, mode="indeterminate", length=180)
        self.progress.pack(side="left")

        # ── Terminal log ──────────────────────────────────────────────────────
        self._log_font_sz = 10
        lh = tk.Frame(self, bg=P["bg2"], padx=10, pady=4)
        lh.pack(fill="x")
        tk.Label(lh, text="  TERMINAL OUTPUT",
            fg=P["text_dim"], bg=P["bg2"], font=FUB).pack(side="left")
        tk.Button(lh, text="+", command=lambda: self._resize_log(+1),
            bg=P["bg3"], fg=P["accent"], relief="groove",
            font=FUB, bd=1, padx=6, pady=1, cursor="hand2"
        ).pack(side="right", padx=(2,0))
        tk.Button(lh, text="-", command=lambda: self._resize_log(-1),
            bg=P["bg3"], fg=P["accent"], relief="groove",
            font=FUB, bd=1, padx=6, pady=1, cursor="hand2"
        ).pack(side="right", padx=(0,2))
        tk.Button(lh, text="Clear", command=self._clear,
            bg=P["bg3"], fg=P["text_dim"], relief="flat",
            font=FU, bd=0, padx=8, cursor="hand2"
        ).pack(side="right", padx=(0,8))

        lvs  = ttk.Scrollbar(self, orient="vertical")
        lhs_ = ttk.Scrollbar(self, orient="horizontal")
        self.log = tk.Text(self,
            bg=P["term_bg"], fg=P["term_fg"],
            font=("Courier", self._log_font_sz),
            relief="flat", bd=0, state="disabled",
            wrap="none", padx=10, pady=8,
            yscrollcommand=lvs.set, xscrollcommand=lhs_.set)
        lvs.config(command=self.log.yview)
        lhs_.config(command=self.log.xview)
        self.log.tag_configure("err",  foreground=P["danger"])
        self.log.tag_configure("warn", foreground=P["accent3"])
        self.log.tag_configure("ok",   foreground=P["accent2"])
        self.log.tag_configure("info", foreground=P["accent"])
        lvs.pack(side="right", fill="y")
        lhs_.pack(side="bottom", fill="x")
        self.log.pack(fill="both", expand=True)

    def _toggle_run_mode(self):
        """Grey out the non-selected panel, enable the selected one."""
        mode = self._run_mode.get()
        # regular panel
        reg_state = "normal" if mode == "regular" else "disabled"
        sif_state = "normal" if mode == "sif"     else "disabled"
        for w in self._reg_widgets:
            try: w.config(state=reg_state)
            except Exception: pass
        for w in self._sif_widgets:
            try: w.config(state=sif_state)
            except Exception: pass
        # dim panel labels via fg colour
        dim = P["text_dim"]
        act = P["text_bright"]
        try:
            self._reg_frame.config(
                fg=P["accent"] if mode=="regular" else dim)
            self._sif_frame.config(
                fg=P["accent3"] if mode=="sif" else dim)
        except Exception:
            pass

    def _resize_log(self, delta):
        self._log_font_sz = max(7, min(24, self._log_font_sz + delta))
        self.log.config(font=("Courier", self._log_font_sz))

    # ── toggle helpers ────────────────────────────────────────────────────────
    def _set_group_state(self, tag, state):
        """Enable or disable all widgets in a named group."""
        for w in self._row_widgets.get(tag, []):
            try:
                w.config(state=state)
            except Exception:
                pass

    def _toggle_sif(self):
        if self.app.use_sif_var.get():
            self._set_group_state("regular", "disabled")
            self._set_group_state("sif",     "normal")
        else:
            self._set_group_state("regular", "normal")
            self._set_group_state("sif",     "disabled")

    # ── browse helpers ────────────────────────────────────────────────────────
    def _br_exe(self):
        p = filedialog.askopenfilename(title="Select CAC Executable")
        if p: self.app.cac_exe_var.set(p)

    def _br_dir(self):
        p = filedialog.askdirectory(title="Select Working Directory")
        if p: self.app.workdir_var.set(p)

    def _br_script(self):
        p = filedialog.askopenfilename(title="Select Input Script",
            filetypes=[("All files", "*.*")])
        if p: self.app.script_var.set(p)

    def _br_sif(self):
        p = filedialog.askopenfilename(title="Select .sif file",
            filetypes=[("Singularity image", "*.sif"), ("All files", "*.*")])
        if p:
            self.app.sif_var.set(p)
            # Auto-detect bind dir from sif path
            parts = p.strip("/").split("/")
            if parts:
                self.app.sif_bind_var.set("/" + parts[0])

    def _br_sif_bind(self):
        p = filedialog.askdirectory(title="Select Singularity Bind Directory")
        if p: self.app.sif_bind_var.set(p)

    # ── run ───────────────────────────────────────────────────────────────────
    def _run(self):
        self._run_local()

    def _run_local(self):
        exe     = self.app.cac_exe_var.get().strip()
        wdir    = self.app.workdir_var.get().strip()
        script  = self.app.script_var.get().strip()
        mpi     = self.app.mpi_var.get().strip()
        modules = self.app.module_var.get().strip()

        if not exe:
            messagebox.showwarning("Missing", "Set the CAC executable path.")
            return
        if not script:
            ef = self.app.editor_tab.get_current_file()
            if ef:
                script = ef
                self.app.script_var.set(script)
        if not script:
            messagebox.showwarning("Missing", "Specify an input script.")
            return
        if self.app.editor_tab._unsaved and script == self.app.editor_tab.get_current_file():
            self.app.editor_tab.save_file()
        if not wdir:
            wdir = os.path.dirname(os.path.abspath(script))
            self.app.workdir_var.set(wdir)

        mod_init = ("source /etc/profile.d/modules.sh 2>/dev/null || "
                    "source /usr/share/lmod/lmod/init/bash 2>/dev/null || true\n")
        mods = "".join("module load %s\n" % m for m in modules.split() if m)
        if mpi and mpi != "0":
            cac_cmd = "mpirun -np %s %s -in %s" % (mpi, exe, script)
        else:
            cac_cmd = "%s -in %s" % (exe, script)
        shell = mod_init + mods + cac_cmd + "\n"

        self._log_append("[LOCAL] dir: %s\n" % wdir, "info")
        if modules:
            self._log_append("[LOCAL] modules: %s\n" % modules, "info")
        self._log_append("[LOCAL] %s\n" % cac_cmd, "info")
        self._log_append("-" * 64 + "\n", "info")
        self._start_job(["bash", "-c", shell], wdir)

    def _run_singularity(self):
        sif     = self.app.sif_var.get().strip()
        bind    = self.app.sif_bind_var.get().strip()
        wdir    = self.app.workdir_var.get().strip()
        script  = self.app.script_var.get().strip()
        mpi     = self.app.sif_mpi_var.get().strip()
        modules = self.app.sif_module_var.get().strip()

        if not sif:
            messagebox.showwarning("Missing", "Set the .sif file path.")
            return
        if not os.path.isfile(sif):
            messagebox.showerror("Not found", "Singularity image not found:\n%s" % sif)
            return
        if not wdir:
            messagebox.showwarning("Missing", "Set the working directory.")
            return
        if not script:
            ef = self.app.editor_tab.get_current_file()
            if ef:
                script = ef
                self.app.script_var.set(script)
        if not script:
            messagebox.showwarning("Missing", "Specify an input script.")
            return
        if self.app.editor_tab._unsaved and script == self.app.editor_tab.get_current_file():
            self.app.editor_tab.save_file()

        # bind dir defaults to working dir if not set
        if not bind:
            bind = wdir

        mod_init = ("source /etc/profile.d/modules.sh 2>/dev/null || "
                    "source /usr/share/lmod/lmod/init/bash 2>/dev/null || true\n")
        mods = "".join("module load %s\n" % m for m in modules.split() if m)
        n = mpi if (mpi and mpi != "0") else "8"
        cac_cmd = ("mpirun -n %s singularity run --bind %s %s -in %s"
                   % (n, bind, sif, script))
        shell = mod_init + mods + cac_cmd + "\n"

        self._log_append("[SIF] dir:  %s\n" % wdir, "info")
        self._log_append("[SIF] bind: %s\n" % bind, "info")
        self._log_append("[SIF] sif:  %s\n" % sif, "info")
        if modules:
            self._log_append("[SIF] modules: %s\n" % modules, "info")
        self._log_append("[SIF] %s\n" % cac_cmd, "info")
        self._log_append("-" * 64 + "\n", "info")
        self._start_job(["bash", "-c", shell], wdir)

    # ── shared job infrastructure ─────────────────────────────────────────────
    def _start_job(self, cmd, cwd):
        self._running = True
        self.run_btn.config(state="disabled")
        try: self.sif_run_btn.config(state="disabled")
        except Exception: pass
        self.kill_btn.config(state="normal")
        self.status_lbl.config(text="Running...", fg=P["accent3"])
        self.progress.start(12)
        threading.Thread(target=self._worker, args=(cmd, cwd), daemon=True).start()

    def _worker(self, cmd, cwd):
        try:
            kw = dict(stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                      text=True, bufsize=1)
            if cwd:
                kw["cwd"] = cwd
            self._process = subprocess.Popen(cmd, **kw)
            for line in self._process.stdout:
                self._log_q.put(("line", line))
            self._process.wait()
            self._log_q.put(("done", self._process.returncode))
        except Exception as e:
            self._log_q.put(("err", str(e)))
            self._log_q.put(("done", -1))

    def _kill(self):
        if self._process:
            self._process.terminate()
            self._log_append("\n[KILL] Terminated by user.\n", "warn")

    def _poll(self):
        try:
            while True:
                kind, data = self._log_q.get_nowait()
                if kind == "line":
                    lo = data.lower()
                    tag = ("err"  if "error" in lo or "fatal" in lo else
                           "warn" if "warning" in lo else None)
                    self._log_append(data, tag)
                elif kind == "err":
                    self._log_append(data + "\n", "err")
                elif kind == "done":
                    self._on_done(data)
        except queue.Empty:
            pass
        self.after(80, self._poll)

    def _on_done(self, rc):
        self._running = False
        self.run_btn.config(state="normal")
        try: self.sif_run_btn.config(state="normal")
        except Exception: pass
        self.kill_btn.config(state="disabled")
        self.progress.stop()
        if rc == 0:
            self.status_lbl.config(text="Finished", fg=P["accent2"])
            self._log_append("\n[DONE] Completed successfully.\n", "ok")
            self.app.tecplot_tab.on_run_complete()
        else:
            self.status_lbl.config(text="Exit (%d)" % rc, fg=P["danger"])
            self._log_append("\n[DONE] Exit code %d.\n" % rc, "err")

    def _log_append(self, text, tag=None):
        self.log.config(state="normal")
        if tag:
            self.log.insert("end", text, tag)
        else:
            self.log.insert("end", text)
        self.log.see("end")
        self.log.config(state="disabled")

    def _clear(self):
        self.log.config(state="normal")
        self.log.delete("1.0", "end")
        self.log.config(state="disabled")

# ─────────────────────────────────────────────────────────────────────────────
#  TECPLOT TAB
# ─────────────────────────────────────────────────────────────────────────────
class TecplotTab(ttk.Frame):
    def __init__(self, parent, app, **kw):
        super(TecplotTab, self).__init__(parent, **kw)
        self.app = app
        self._build()

    def _build(self):
        hdr = tk.Frame(self, bg=P["bg2"], padx=14, pady=10)
        hdr.pack(fill="x")
        tk.Label(hdr, text="TecPlot Integration",
            fg=P["text_bright"], bg=P["bg2"], font=FT).pack(side="left")
        tk.Button(hdr, text="Open TecPlot 360",
            command=self._open_tecplot,
            bg=P["accent"], fg=P["bg"],
            font=FUB, relief="flat", bd=0,
            padx=12, pady=5, cursor="hand2"
        ).pack(side="right", padx=(8,0))
        avail = ("pytecplot available" if TECPLOT_AVAILABLE
                 else "pytecplot not installed  (pip install pytecplot)")
        tk.Label(hdr, text=avail,
            fg=P["accent2"] if TECPLOT_AVAILABLE else P["danger"],
            bg=P["bg2"], font=FU).pack(side="right")

        cfg = tk.Frame(self, bg=P["bg"], padx=14, pady=10)
        cfg.pack(fill="x")
        cfg.columnconfigure(1, weight=1)

        self.datafile_var = tk.StringVar()
        tk.Label(cfg, text="Output Data File", fg=P["text_dim"], bg=P["bg"],
                 font=FU, width=20, anchor="w").grid(row=0, column=0, sticky="w", pady=4)
        tk.Entry(cfg, textvariable=self.datafile_var,
                 bg=P["bg3"], fg=P["text"], insertbackground=P["accent"],
                 relief="flat", font=FMS, bd=4
                 ).grid(row=0, column=1, sticky="ew", padx=4)
        tk.Button(cfg, text="Browse", command=self._br_data,
                  bg=P["bg3"], fg=P["accent"], relief="flat",
                  font=FU, bd=0, padx=8, cursor="hand2"
                  ).grid(row=0, column=2, padx=4)

        self.autoload_var = tk.BooleanVar(value=True)
        tk.Checkbutton(cfg,
            text="Auto-load after run completes",
            variable=self.autoload_var,
            bg=P["bg"], fg=P["text"], selectcolor=P["bg3"],
            activebackground=P["bg"], font=FU
            ).grid(row=1, column=1, sticky="w", pady=2)

        self.port_var = tk.StringVar(value="7600")
        tk.Label(cfg, text="TecPlot Port", fg=P["text_dim"], bg=P["bg"],
                 font=FU, width=20, anchor="w").grid(row=2, column=0, sticky="w", pady=4)
        tk.Entry(cfg, textvariable=self.port_var, width=10,
                 bg=P["bg3"], fg=P["text"], insertbackground=P["accent"],
                 relief="flat", font=FMS, bd=4
                 ).grid(row=2, column=1, sticky="w", padx=4)

        # ── TecPlot connection how-to box ────────────────────────────────────
        hint_frame = tk.Frame(self, bg=P["hint_bg"], padx=14, pady=8)
        hint_frame.pack(fill="x")

        tk.Label(hint_frame, text="HOW TO CONNECT",
            fg=P["hint_title"], bg=P["hint_bg"], font=FUB).pack(anchor="w")

        if TECPLOT_AVAILABLE:
            step4 = "4.  Run this command in a terminal to verify the connection:"
        else:
            step4 = "4.  pytecplot NOT installed. Run the command below to install it:"
        steps = (
            "1.  Open TecPlot 360 and load your data file.",
            "2.  In TecPlot menu:  Scripting  ->  PyTecplot Connections...",
            "3.  Check  \"Accept connections\"  and click OK.",
            step4,
        )
        for s in steps:
            tk.Label(hint_frame, text=s,
                fg=P["text"], bg=P["hint_bg"], font=FMS, anchor="w"
            ).pack(fill="x")

        # copyable terminal command
        cmd_frame = tk.Frame(hint_frame, bg=P["bg3"], padx=8, pady=4)
        cmd_frame.pack(fill="x", pady=(2,6))
        if TECPLOT_AVAILABLE:
            cmd_txt = "python3 -c \"import tecplot as tp; tp.session.connect(); print(\'Connected\')\""
        else:
            cmd_txt = "pip install pytecplot"
        self._tp_cmd_var = tk.StringVar(value=cmd_txt)
        cmd_entry = tk.Entry(cmd_frame, textvariable=self._tp_cmd_var,
            bg=P["bg3"], fg=P["hint_ex"], font=FMS,
            relief="flat", bd=0, state="readonly",
            readonlybackground=P["bg3"])
        cmd_entry.pack(side="left", fill="x", expand=True)

        def _copy_cmd():
            self.clipboard_clear()
            self.clipboard_append(cmd_txt)
            copy_btn.config(text="Copied!")
            self.after(1500, lambda: copy_btn.config(text="Copy"))

        copy_btn = tk.Button(cmd_frame, text="Copy", command=_copy_cmd,
            bg=P["accent"], fg=P["bg"], font=FUB,
            relief="flat", bd=0, padx=8, cursor="hand2")
        copy_btn.pack(side="right", padx=(4,0))

        tk.Label(hint_frame,
            text="5.  Then click  Connect  below (port 7600 is default).",
            fg=P["text"], bg=P["hint_bg"], font=FMS, anchor="w"
        ).pack(fill="x")

        # ── connection / load buttons ────────────────────────────────────────
        bf = tk.Frame(self, bg=P["bg"], padx=14, pady=6)
        bf.pack(fill="x")

        def ab(text, cmd, color=P["accent"]):
            return tk.Button(bf, text=text, command=cmd,
                bg=color, fg=P["bg"] if color != P["bg3"] else P["text"],
                font=FUB, relief="flat", bd=0,
                padx=14, pady=7, cursor="hand2")

        ab("Connect",     self._connect,     P["accent"]).pack(side="left", padx=(0,8))
        ab("Load Now",    self._load_now,    P["accent2"]).pack(side="left", padx=(0,8))
        ab("Close Model", self._close_model, P["danger"]).pack(side="left")

        # ── Zone Style Presets ───────────────────────────────────────────────
        zf_outer = tk.Frame(self, bg=P["bg2"], padx=14, pady=10)
        zf_outer.pack(fill="x")

        hrow = tk.Frame(zf_outer, bg=P["bg2"])
        hrow.pack(fill="x", pady=(0,8))
        tk.Label(hrow, text="ZONE STYLE PRESETS",
            fg=P["text_dim"], bg=P["bg2"], font=FUB).pack(side="left")
        tk.Label(hrow, text="Applies zone display style to TecPlot based on model type",
            fg=P["text_dim"], bg=P["bg2"],
            font=("Arial", 8)).pack(side="left", padx=10)

        zf = tk.Frame(zf_outer, bg=P["bg2"])
        zf.pack(fill="x")

        def zpreset(parent, label, cmd, color, desc):
            f = tk.Frame(parent, bg=P["bg2"])
            f.pack(side="left", padx=(0,16))
            tk.Button(f, text=label, command=cmd,
                bg=color, fg=P["bg"],
                font=FUB, relief="flat", bd=0,
                padx=18, pady=8, cursor="hand2"
            ).pack()
            tk.Label(f, text=desc,
                fg=P["text_dim"], bg=P["bg2"],
                font=("Arial", 8), justify="left"
            ).pack(anchor="w", pady=(3,0))

        zpreset(zf, "Apply CG", self._zone_cg, P["accent"],
                "Elements only\nAll zones: boundary surface")
        zpreset(zf, "Apply CAC", self._zone_cac, P["accent2"],
                "Elements + atoms\nElements: surface  |  Atoms: scatter")

        # ── Axis ratio checkbox ───────────────────────────────────────────────
        ax_row = tk.Frame(zf_outer, bg=P["bg2"])
        ax_row.pack(fill="x", pady=(8,0))
        self.axis_equal_var = tk.BooleanVar(value=False)
        tk.Checkbutton(ax_row,
            text="Lock XYZ aspect ratio 1:1:1  (XYZDependent axes)",
            variable=self.axis_equal_var,
            command=self._apply_axis_ratio,
            bg=P["bg2"], fg=P["text"],
            selectcolor=P["bg3"],
            activebackground=P["bg2"],
            font=FU, cursor="hand2"
        ).pack(side="left")
        tk.Label(ax_row,
            text="applies immediately when checked",
            fg=P["text_dim"], bg=P["bg2"],
            font=("Arial", 8)
        ).pack(side="left", padx=8)


        self._tp_log_font_sz = 10

        lhdr = tk.Frame(self, bg=P["bg2"], padx=10, pady=4)
        lhdr.pack(fill="x", pady=(8,0))
        tk.Label(lhdr, text="  LOG",
            fg=P["text_dim"], bg=P["bg2"], font=FUB).pack(side="left")
        tk.Button(lhdr, text="+",
            command=lambda: self._resize_tp_log(+1),
            bg=P["bg3"], fg=P["accent"], relief="groove",
            font=FUB, bd=1, padx=6, pady=1, cursor="hand2"
        ).pack(side="right", padx=(2,0))
        tk.Button(lhdr, text="-",
            command=lambda: self._resize_tp_log(-1),
            bg=P["bg3"], fg=P["accent"], relief="groove",
            font=FUB, bd=1, padx=6, pady=1, cursor="hand2"
        ).pack(side="right", padx=(0,2))

        self.tp_log = scrolledtext.ScrolledText(self,
            bg=P["term_bg"], fg=P["term_fg"],
            font=("Courier", self._tp_log_font_sz), relief="flat", bd=0,
            state="disabled", padx=10, pady=8)
        self.tp_log.pack(fill="both", expand=True, pady=(0,0))
        self.tp_log.tag_configure("err",  foreground=P["danger"])
        self.tp_log.tag_configure("ok",   foreground=P["accent2"])
        self.tp_log.tag_configure("warn", foreground=P["accent3"])

        self._log("TecPlot Integration ready.\n")
        if not TECPLOT_AVAILABLE:
            self._log("pip install pytecplot\n"
                      "Then enable: Scripting -> Allow PyTecplot Connections\n", "warn")

    def _resize_tp_log(self, delta):
        self._tp_log_font_sz = max(7, min(24, self._tp_log_font_sz + delta))
        self.tp_log.config(font=("Courier", self._tp_log_font_sz))

    def _br_data(self):
        p = filedialog.askopenfilename(
            title="Select TecPlot Data File",
            filetypes=[("TecPlot","*.dat *.plt *.szplt"),("All","*.*")])
        if p: self.datafile_var.set(p)

    def _log(self, msg, tag=None):
        self.tp_log.config(state="normal")
        line = "[%s] %s" % (time.strftime("%H:%M:%S"), msg)
        if tag: self.tp_log.insert("end", line, tag)
        else:   self.tp_log.insert("end", line)
        self.tp_log.see("end")
        self.tp_log.config(state="disabled")

    def _resolve(self):
        # type: () -> Optional[str]
        ex = self.datafile_var.get().strip()
        if ex and os.path.isfile(ex): return ex
        wdir   = self.app.workdir_var.get().strip()
        detect = self.app.detected_output_file
        if detect:
            for base in [wdir, os.path.dirname(self.app.script_var.get())]:
                if base:
                    c = os.path.join(base, detect)
                    if os.path.isfile(c): return c
        if wdir:
            for ext in ("*.plt","*.dat","*.szplt"):
                hits = sorted(glob.glob(os.path.join(wdir, ext)),
                              key=os.path.getmtime, reverse=True)
                if hits:
                    self._log("Auto-detected: %s\n" % hits[0], "warn")
                    return hits[0]
        return None

    def _connect(self):
        if not TECPLOT_AVAILABLE:
            self._log("pytecplot not installed.\n","err"); return
        try:
            tp.session.connect(port=int(self.port_var.get()))
            self._log("Connected on port %s.\n" % self.port_var.get(),"ok")
        except Exception as e:
            self._log("Connection failed: %s\n" % e,"err")

    def _load_now(self):
        if not TECPLOT_AVAILABLE:
            self._log("pytecplot not installed.\n","err"); return
        f = self._resolve()
        if not f: self._log("No data file found.\n","err"); return
        self._do_load(f)

    def _do_load(self, path):
        # type: (str) -> None
        self._log("Loading: %s\n" % path)
        try:
            ds = tp.data.load_tecplot(path)
            self._log("Loaded: %d zone(s), %d var(s).\n" % (
                ds.num_zones, ds.num_variables), "ok")
            try:
                if ds.num_variables >= 3:
                    tp.active_frame().plot(
                        tp.constant.PlotType.Cartesian3D).activate()
                    self._log("Plot type -> Cartesian3D.\n","ok")
            except Exception:
                pass
        except Exception as e:
            self._log("Load failed: %s\n" % e,"err")

    def _open_tecplot(self):
        for c in ["tec360","tecplot","tec360evo",
                  "/opt/tecplot/360ex_2024r1/bin/tec360",
                  "/usr/local/tecplot/360ex/bin/tec360"]:
            try:
                subprocess.Popen([c], close_fds=True)
                self._log("Launched: %s\n" % c,"ok"); return
            except FileNotFoundError:
                continue
        p = filedialog.askopenfilename(title="Locate TecPlot Executable")
        if p:
            subprocess.Popen([p], close_fds=True)
            self._log("Launched: %s\n" % p,"ok")

    # ── Zone Style Presets ──────────────────────────────────────────────────
    def _require_tp(self):
        if not TECPLOT_AVAILABLE:
            self._log("pytecplot not installed.\n","err")
            return False
        return True

    # ── correct pytecplot API zone style helpers ──────────────────────────────
    # API refs:
    #   plot.fieldmap(zone)          — access fieldmap by zone object (0-based index)
    #   fmap.show = True/False       — show/hide zone in active fieldmaps
    #   fmap.surfaces.surfaces_to_plot = SurfacesToPlot.BoundaryFaces
    #   fmap.mesh.show               — wire mesh overlay
    #   fmap.scatter.show            — scatter layer on this fieldmap
    #   fmap.scatter.size            — scatter symbol size (% of frame height)
    #   fmap.scatter.symbol().shape  — GeomShape.Circle / Square / Diamond etc.
    #   plot.show_scatter = True     — must enable scatter layer at PLOT level too
    #   plot.show_mesh    = True     — enable mesh layer at plot level

    def _get_plot_3d(self):
        """Return active Cartesian3D plot, switching if needed."""
        frame = tp.active_frame()
        plot  = frame.plot(tp.constant.PlotType.Cartesian3D)
        plot.activate()
        return frame, plot

    def _surface_on(self, fmap):
        """Configure fieldmap as boundary surface (solid faces, no scatter)."""
        fmap.show = True
        fmap.surfaces.surfaces_to_plot = tp.constant.SurfacesToPlot.BoundaryFaces
        fmap.mesh.show    = False
        fmap.scatter.show = False
        fmap.contour.show = False
        try: fmap.shade.show = False
        except Exception: pass

    def _coarse_on(self, fmap):
        """Coarse Element (zone 2): boundary surface WITH shade."""
        self._surface_on(fmap)
        try:
            fmap.shade.show                = True
            fmap.shade.use_lighting_effect  = True
        except Exception:
            pass

    def _scatter_on(self, fmap):
        """Discrete Atoms: scatter points, no shade, no surface."""
        fmap.show = True
        fmap.scatter.show = True
        fmap.scatter.size = 0.4
        try:
            fmap.scatter.symbol().shape = tp.constant.GeomShape.Point
        except Exception:
            pass
        fmap.surfaces.surfaces_to_plot = tp.constant.SurfacesToPlot.None_
        fmap.mesh.show    = False
        fmap.contour.show = False
        try: fmap.shade.show = False
        except Exception: pass

    def _hide(self, fmap):
        fmap.show = False

    def _simcell_on(self, fmap):
        """Simulation Cell (zone 3): boundary surface + overlap mesh, black, 0.5%."""
        fmap.show = True
        fmap.surfaces.surfaces_to_plot = tp.constant.SurfacesToPlot.BoundaryFaces
        fmap.mesh.show = True
        try:
            fmap.mesh.line_pattern  = tp.constant.LinePattern.Solid
            fmap.mesh.color         = tp.constant.Color.Black
            fmap.mesh.line_thickness = 0.5
            fmap.mesh.overlay_mode  = tp.constant.MeshOverlayMode.Overlay
        except Exception:
            pass
        fmap.scatter.show = False
        fmap.contour.show = False
        try:
            fmap.shade.show = False
            pass
        except Exception:
            pass

    def _apply_preset(self, surface_kws, atom_kws, fallback="surface"):
        """Core preset logic. Returns (surface_count, scatter_count)."""
        frame, plot = self._get_plot_3d()
        ds = frame.dataset

        # enable relevant plot-level layers
        plot.show_scatter = True
        plot.show_mesh    = False

        s_count = a_count = other_count = 0
        zones_list = list(ds.zones())

        for z in zones_list:
            nm  = z.name.lower()
            fm  = plot.fieldmap(z)          # pass zone object directly
            if surface_kws and any(k in nm for k in surface_kws):
                self._surface_on(fm);  s_count += 1
            elif atom_kws and any(k in nm for k in atom_kws):
                self._scatter_on(fm);  a_count += 1
            else:
                if fallback == "surface":
                    self._surface_on(fm);  other_count += 1
                elif fallback == "scatter":
                    self._scatter_on(fm);  other_count += 1
                else:
                    self._hide(fm)

        # if nothing matched at all, show everything as fallback
        if s_count + a_count == 0:
            for z in zones_list:
                fm = plot.fieldmap(z)
                if fallback == "scatter":
                    self._scatter_on(fm)
                else:
                    self._surface_on(fm)
            self._log("No keyword-matched zones — applied %s to all %d zone(s).\n" % (
                fallback, len(zones_list)), "warn")
        return s_count, a_count

    def _zone_cg(self):
        """CG preset: Zone1=Discrete Atoms(hidden), Zone2=Coarse Element(surface),
        Zone3=Simulation Cell(surface+mesh). Extra zones -> surface."""
        if not self._require_tp(): return
        self._log("Zone preset: CG\n")
        try:
            frame, plot = self._get_plot_3d()
            plot.show_scatter = False
            plot.show_mesh    = True
            plot.show_shade   = True
            zones = list(frame.dataset.zones())
            for i, z in enumerate(zones):
                fm = plot.fieldmap(z)
                if i == 0:    # Zone 1 — Discrete Atoms: hide in CG
                    self._hide(fm)
                elif i == 2:  # Zone 3 — Simulation Cell
                    self._simcell_on(fm)
                else:         # Zone 2 — Coarse Element + any extras
                    self._coarse_on(fm)
            self._log("  -> CG: %d zone(s). Zone 3 = sim cell mesh.\n" % len(zones), "ok")
        except Exception as e:
            self._log("Zone preset error: %s\n" % e, "err")

    def _zone_cac(self):
        """CAC preset: Zone1=Discrete Atoms(scatter), Zone2=Coarse Element(surface),
        Zone3=Simulation Cell(surface+mesh). Extra zones -> surface."""
        if not self._require_tp(): return
        self._log("Zone preset: CAC\n")
        try:
            frame, plot = self._get_plot_3d()
            plot.show_scatter = True
            plot.show_mesh    = True
            plot.show_shade   = True
            zones = list(frame.dataset.zones())
            for i, z in enumerate(zones):
                fm = plot.fieldmap(z)
                if i == 0:    # Zone 1 — Discrete Atoms: scatter
                    self._scatter_on(fm)
                elif i == 2:  # Zone 3 — Simulation Cell: surface + mesh
                    self._simcell_on(fm)
                else:         # Zone 2 — Coarse Element + extras: surface
                    self._coarse_on(fm)
            self._log("  -> CAC: %d zone(s). Zone1=scatter Zone2=surface Zone3=mesh.\n" % len(zones), "ok")
        except Exception as e:
            self._log("Zone preset error: %s\n" % e, "err")


    # ── Axis ratio ───────────────────────────────────────────────────────────
    def _apply_axis_ratio(self):
        """Lock X:Y=1:1 and X:Z=1:1 using XYZDependent + explicit ratios."""
        if not self._require_tp(): return
        try:
            frame, plot = self._get_plot_3d()
            if self.axis_equal_var.get():
                # Step 1: set dependency mode
                plot.axes.axis_mode = tp.constant.AxisMode.XYZDependent
                # Step 2: force ratios to exactly 100% (= 1:1:1)
                # xy_ratio and xz_ratio are in percent:
                # 100 means X:Y = 1:1, X:Z = 1:1
                plot.axes.xy_ratio = 1
                plot.axes.xz_ratio = 1
                try:
                    plot.view.fit()
                except Exception:
                    pass
                self._log("Axis -> XYZDependent: X:Y=1:1, X:Z=1:1.\n", "ok")
            else:
                plot.axes.axis_mode = tp.constant.AxisMode.Independent
                self._log("Axis -> Independent.\n", "ok")
        except Exception as e:
            self._log("Axis ratio error: %s\n" % e, "err")

    def _close_model(self):
        """Clear current TecPlot data — equivalent to File -> New Layout."""
        if not self._require_tp(): return
        try:
            if not messagebox.askyesno(
                    "Close Model",
                    "This will clear all data in TecPlot.\nContinue?"):
                return
            tp.new_layout()
            self._log("TecPlot layout cleared (new layout).\n", "ok")
        except Exception as e:
            self._log("Close model error: %s\n" % e, "err")

    def on_run_complete(self):
        if self.autoload_var.get():
            self._log("Run complete - auto-loading...\n","ok")
            self.after(1500, self._auto_load)

    def _auto_load(self):
        f = self._resolve()
        if f:
            self.datafile_var.set(f)
            if TECPLOT_AVAILABLE: self._do_load(f)
            else: self._log("File ready: %s\n" % f,"warn")
        else:
            self._log("Could not locate output file.\n","warn")

# ─────────────────────────────────────────────────────────────────────────────
#  MAIN WINDOW
# ─────────────────────────────────────────────────────────────────────────────
class CACModelEditor(tk.Tk):
    def __init__(self):
        super(CACModelEditor, self).__init__()
        self.title("CAC Model Editor  V2")
        self.geometry("1280x860")
        self.minsize(960, 640)
        self.configure(bg=P["bg"])
        self._theme()

        self.save_exe_var    = tk.BooleanVar(value=False)
        self.save_sif_var    = tk.BooleanVar(value=False)
        self.save_manual_var = tk.BooleanVar(value=False)
        _cfg = self._load_config()

        self.cac_exe_var = tk.StringVar(
            value=_cfg.get("default_exe",
                "/mnt/yan_2/CAC-main/V3.1.1/src/cac.3.1.1"))
        self.script_var           = tk.StringVar()
        self.mpi_var              = tk.StringVar(value="4")
        self.module_var           = tk.StringVar(value="CAC")
        self.workdir_var          = tk.StringVar(value=os.getcwd())
        self.detected_output_file = None
        self.manual_path_var      = tk.StringVar(
            value=_cfg.get("default_manual", ""))

        self.use_sif_var     = tk.BooleanVar(value=False)
        self.sif_var         = tk.StringVar(
            value=_cfg.get("default_sif", ""))
        self.sif_bind_var    = tk.StringVar(value="")
        self.sif_mpi_var     = tk.StringVar(value="4")
        self.sif_module_var  = tk.StringVar(value="")

        if _cfg.get("save_exe"):    self.save_exe_var.set(True)
        if _cfg.get("save_sif"):    self.save_sif_var.set(True)
        if _cfg.get("save_manual"): self.save_manual_var.set(True)

        self._build()
        self.protocol("WM_DELETE_WINDOW", self._close)

    def _theme(self):
        s = ttk.Style(self)
        s.theme_use("clam")
        s.configure(".",
            background=P["bg"], foreground=P["text"],
            fieldbackground=P["bg3"], troughcolor=P["border"],
            bordercolor=P["border"],
            darkcolor=P["bg2"], lightcolor=P["bg3"],
            relief="groove")
        s.configure("TNotebook",
            background=P["bg2"], borderwidth=1,
            tabmargins=[2, 2, 0, 0])
        s.configure("TNotebook.Tab",
            background=P["bg2"], foreground=P["text_dim"],
            padding=[14, 5], font=FUB, borderwidth=1)
        s.map("TNotebook.Tab",
            background=[("selected", P["bg3"]), ("active", "#C8C8C8")],
            foreground=[("selected", P["text_bright"]), ("active", P["text"])],
            relief=[("selected", "flat"), ("!selected", "raised")])
        s.configure("TScrollbar",
            background=P["bg2"], troughcolor=P["border"],
            arrowcolor=P["text_dim"], borderwidth=1, relief="raised")
        s.map("TScrollbar",
            background=[("active", P["bg3"])])
        s.configure("TProgressbar",
            troughcolor=P["border"], background=P["accent2"], borderwidth=1)
        s.configure("TFrame", background=P["bg"])
        s.configure("Treeview",
            background=P["bg3"], foreground=P["text"],
            fieldbackground=P["bg3"], borderwidth=1,
            relief="sunken", rowheight=20)
        s.configure("Treeview.Heading",
            background=P["bg2"], foreground=P["text_dim"],
            font=FUB, relief="raised", borderwidth=1)
        s.map("Treeview",
            background=[("selected", P["accent"])],
            foreground=[("selected", "#FFFFFF")])
        s.configure("TCombobox",
            background=P["bg3"], foreground=P["text"],
            fieldbackground=P["bg3"], arrowcolor=P["text_dim"],
            borderwidth=1)
        s.map("TCombobox",
            fieldbackground=[("readonly", P["bg3"])],
            background=[("active", P["bg2"])])

    def _build(self):
        ban = tk.Frame(self, bg=P["bg2"], height=46)
        ban.pack(fill="x")
        ban.pack_propagate(False)
        tk.Label(ban, text="  CAC MODEL EDITOR",
            fg=P["text_bright"], bg=P["bg2"],
            font=(FT[0], 14, "bold")).pack(side="left", padx=12)
        tk.Label(ban, text="V2  |  Linux",
            fg=P["text_dim"], bg=P["bg2"], font=FU).pack(side="right", padx=12)

        nb = ttk.Notebook(self)
        nb.pack(fill="both", expand=True)

        self.editor_tab  = ScriptEditorTab(nb, self)
        self.run_tab     = RunTab(nb, self)
        self.tecplot_tab = TecplotTab(nb, self)

        nb.add(self.editor_tab,  text="  Script Editor  ")
        nb.add(self.run_tab,     text="  Run & Monitor  ")
        nb.add(self.tecplot_tab, text="  TecPlot  ")

        sb = tk.Frame(self, bg=P["bg2"], height=24)
        sb.pack(fill="x", side="bottom")
        sb.pack_propagate(False)
        tk.Label(sb, text="Ready", fg=P["text_dim"], bg=P["bg2"],
                 font=FU, anchor="w", padx=10).pack(side="left")

        self.bind("<Control-o>", lambda _: self.editor_tab.open_file())
        self.bind("<Control-s>", lambda _: self.editor_tab.save_file())
        self.bind("<Control-r>", lambda _: self.run_tab._run())
        self.bind("<F5>",        lambda _: self.run_tab._run())

    def _config_path(self):
        return os.path.join(os.path.expanduser("~"), ".cac_editor.conf")

    def _load_config(self):
        cfg = {}
        try:
            with open(self._config_path()) as f:
                for line in f:
                    line = line.strip()
                    if "=" in line and not line.startswith("#"):
                        k, _, v = line.partition("=")
                        cfg[k.strip()] = v.strip()
            for key in ("save_exe","save_sif","save_manual"):
                cfg[key] = (cfg.get(key) == "true")
        except (IOError, OSError):
            pass
        return cfg

    def _write_config(self):
        try:
            with open(self._config_path(), "w") as f:
                s_exe = self.save_exe_var.get()
                s_sif = self.save_sif_var.get()
                f.write("save_exe=%s\n" % ("true" if s_exe else "false"))
                if s_exe:
                    f.write("default_exe=%s\n" % self.cac_exe_var.get().strip())
                f.write("save_sif=%s\n" % ("true" if s_sif else "false"))
                if s_sif:
                    f.write("default_sif=%s\n" % self.sif_var.get().strip())
                s_man = self.save_manual_var.get()
                f.write("save_manual=%s\n" % ("true" if s_man else "false"))
                if s_man:
                    f.write("default_manual=%s\n" % self.manual_path_var.get().strip())
        except (IOError, OSError):
            pass

    def _save_default_exe(self): self._write_config()
    def _save_default_sif(self):    self._write_config()
    def _save_default_manual(self): self._write_config()

    def _close(self):
        if any([self.save_exe_var.get(), self.save_sif_var.get(),
                self.save_manual_var.get()]):
            self._write_config()
        if hasattr(self, "run_tab") and self.run_tab._running:
            if not messagebox.askyesno("Job running",
                    "A job is still running. Kill and exit?"):
                return
            self.run_tab._kill()
        self.destroy()
if __name__ == "__main__":
    app = CACModelEditor()
    if len(sys.argv) > 1 and os.path.isfile(sys.argv[1]):
        app.editor_tab.open_file(sys.argv[1])
    app.mainloop()
