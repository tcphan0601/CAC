# CAC Model Editor  v1.0 by Chunxu Yan (ychunxu@ncsu.edu)

A desktop GUI for editing, running, and visualizing CAC simulations — replacing the write-script → terminal → TecPlot cycle with a single integrated tool.

---

## Requirements

| Dependency | How to install |
|---|---|
| **Python 3.6+** | Already on most Linux systems |
| **tkinter** | `sudo apt install python3-tk` (Ubuntu/Debian) |
| | `sudo dnf install python3-tkinter` (RHEL/Rocky) |
| | `conda install -c conda-forge tk` (no sudo needed) |
| **pytecplot** *(optional)* | `pip install pytecplot` — only needed for TecPlot auto-load |

Verify tkinter works: `python3 -c "import tkinter; print('OK')"`

---

## Launch

```bash
python3 CACModelEditor.py

# or open a script directly:
python3 CACModelEditor.py my_model.cac
```

No installation. Single file.

---

## Tabs

### Script Editor
- Open / Save CAC input scripts (any file extension)
- Syntax highlighting for all CAC/LAMMPS commands
- **Variable panel** (right side) — evaluates all `variable equal` expressions live, resolves chains and forward references, marks LAMMPS runtime values (lx, bound, etc.) as `(runtime)`
- **Command hint panel** (bottom) — shows syntax, arguments, and example for whichever CAC command your cursor is on
- **Manual...** button (toolbar) — point to `Manual.md` for live hints from the actual manual

### Run & Monitor
- Set CAC executable, working directory, input script, MPI ranks, module names
- **Set as default executable** checkbox — saves exe path to `~/.cac_editor.conf`, restored on next launch
- **RUN CAC** — runs `./cac -in script.cac` (or with `mpirun`) via bash, sourcing the module system automatically
- Live terminal output — yellow text on black background
- `+` / `−` buttons adjust terminal font size
- **KILL** button to terminate a running job

### TecPlot
- Step-by-step connection instructions with copyable terminal command
- **Open TecPlot 360** button (top right)
- Auto-loads output file (from `write_tecplot` or `write_data` in script) when run completes
- **Apply CG** — all zones as boundary surface (elements-only model)
- **Apply CAC** — element zones as surface, atom zones as scatter points
- **Lock XYZ 1:1:1** checkbox — sets `AxisMode.XYZDependent` immediately

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
| `Ctrl+R` / `F5` | Run CAC |

---

## Default Settings

Edit these in the **Run & Monitor** tab:

| Setting | Default |
|---|---|
| CAC executable | `/mnt/yan_2/CAC-main/V3.1.1/src/cac.3.1.1` |
| MPI ranks | `4` |
| module load | `CAC` |
| TecPlot port | `7600` |

The exe path is saved to `~/.cac_editor.conf` when the **Set as default** checkbox is on.
