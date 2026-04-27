# CAC Model Editor  V2

Desktop GUI for editing, running, and visualizing CAC simulations.


# V2 version update
1. Add singularity run pannal (supported CAC version V3.1.2)
2. Add tecplot simulation box (zone style) preset for CAC/CG
3. fix bug: **Default** button, to restore default on next launch
4. fix bug:  Tecplot x-y-z 1:1:1 ratio checkbox not working
---

## Requirements

| | Install |
|---|---|
| **Python 3.6+** | Already on most Linux systems |
| **tkinter** | `sudo apt install python3-tk` · `sudo dnf install python3-tkinter` · `conda install -c conda-forge tk` |
| **pytecplot** *(optional)* | `pip install pytecplot` — only needed for TecPlot integration |

---

## Launch

```bash
python3 CACModelEditorV2.py
python3 CACModelEditorV2.py my_model.cac   # open script directly
```

---

## Tabs

### Script Editor
- Open / save CAC scripts (any extension)
- Syntax highlighting + line numbers
- **Variable panel** — evaluates all `variable equal` expressions live, resolves chains and forward references; marks LAMMPS runtime values as `(runtime)`
- **Command hint panel** — shows syntax, args, and example for the CAC command under your cursor
- **Manual...** button — point to `Manual.md` for live hints; check **Default** to restore on next launch

### Run & Monitor

**Shared (top):** Working Dir and Input Script — apply to both run modes.

**Run method** radio switch:
- **Regular** — `./cac -in script` or `mpirun -np N ./cac -in script`
  - CAC executable (with **Set as default** checkbox)
  - MPI ranks (default 4), module load
- **Singularity** — `mpirun -n N singularity run --bind /dir /path/cac.sif -in script`
  - .sif path (with **Set as default** checkbox), Bind dir, MPI ranks, module load

Non-selected panel is greyed out. `module load` is optional in both modes — leave blank to skip.

Shared: **KILL** button, status, progress bar, terminal log (yellow on black, font size `+`/`−`).

### TecPlot
- Step-by-step connection instructions with copyable terminal command
- **Open TecPlot 360** button (top right)
- Auto-loads output file (`write_tecplot` or `write_data`) when run completes
- **Apply CG** — Zone 2 (Coarse Element): boundary surface + shade; Zone 3 (Sim Cell): mesh only
- **Apply CAC** — Zone 1 (Atoms): scatter points; Zone 2 (Coarse Element): surface + shade; Zone 3 (Sim Cell): mesh only
- **Lock XYZ 1:1:1** checkbox — sets `XYZDependent` with `xy_ratio=1`, `xz_ratio=1`
- **Close Model** — clears current TecPlot data (File → New Layout)

---

## Defaults (saved to `~/.cac_editor.conf`)

| Field | Default | Save checkbox |
|---|---|---|
| CAC executable | `/mnt/yan_2/CAC-main/V3.1.1/src/cac.3.1.1` | in Regular panel |
| .sif path | — | in Singularity panel |
| Manual path | — | next to Manual... button |

---

## Connecting TecPlot

1. Open TecPlot 360
2. **Scripting → PyTecplot Connections…** → check **Accept connections** → OK
3. Click **Connect** in the TecPlot tab (port 7600)

---

## Keyboard Shortcuts

| Key | Action |
|---|---|
| `Ctrl+O` | Open script |
| `Ctrl+S` | Save script |
| `Ctrl+R` / `F5` | Run (active mode) |
