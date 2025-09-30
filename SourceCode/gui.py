"""
GUI skeleton for running Molecular Dynamics simulations in this project.

How this file connects to the rest of the project:
- settings.json: SettingsModel loads and saves the JSON file used by PreProcessing.readSettings.
- POSCAR: The GUI lets the user select a POSCAR path. If not provided, PreProcessing currently
  generates a default FCC Cu structure; the GUI mirrors that behavior and passes the path if present.
- PreProcessing.py: The GUI constructs a PreProcessing instance with (settings_path, poscar_path, flags)
  where flags emulates command-line pairs (e.g., ["-T", "300", "-E", "NVT"]).
- MDBase.py: The PreProcessing.createMD returns an MDBase instance; the GUI calls runMD and visualizeTraj
  on that instance to execute and visualize the simulation.
- MolecularDynamics.py: That file demonstrates the CLI entry-point path; this GUI uses the same flow
  programmatically rather than via sys.argv.
- data.traj / data.log: MDBase writes these (as configured by output_file setting). The GUI provides a
  button to visualize the trajectory file produced by the simulation.

This is a skeleton focusing on structure and wiring. Validation and error handling are basic and
can be extended as needed.
"""
from __future__ import annotations

import json
import os
import threading
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from typing import Any, Dict, List, Optional

# Project-local imports
from PreProcessing import PreProcessing  # Uses settings.json + POSCAR and maps flags to settings
# MDBase is used indirectly through PreProcessing.createMD()


class SettingsModel:
    """Load/save settings.json and provide helpers for GUI binding.

    Connections:
    - Reads and writes to settings.json in the project root so PreProcessing.readSettings gets the same values.
    - Translates GUI fields to the flag pairs expected by PreProcessing.expected_keys ("-T" -> Temperature, "-E" -> Ensemble).
    """

    DEFAULT_SETTINGS_PATH = os.path.join(os.path.dirname(__file__), "settings.json")

    def __init__(self, settings_path: Optional[str] = None):
        self.settings_path = settings_path or self.DEFAULT_SETTINGS_PATH
        self.data: Dict[str, Any] = {}
        self.load()

    def load(self):
        try:
            with open(self.settings_path, "r") as f:
                self.data = json.load(f)
        except FileNotFoundError:
            # Minimal defaults consistent with project usage
            self.data = {
                "Temperature": 293,
                "Ensemble": "NVT",
                # Additional MDBase init fields: these are optional and can be extended
                "Timestep_fs": 10,
                "Steps": 200,
                "Interval": 10,
                "OutputFile": "data",
                "Friction": 0.01,
                "Potential": "LJ",
                "Integrator": "NVT",
                "Pressure": 1.0e6,
                "Compressibility": 1.0e-11,
            }

    def save(self):
        with open(self.settings_path, "w") as f:
            json.dump(self.data, f, indent=2)

    def to_flags(self) -> List[str]:
        """Return flag pairs consumed by PreProcessing.readTerminalInput.

        PreProcessing.expected_keys = {"-T": "Temperature", "-E": "Ensemble"}
        So we only produce flags for these keys; all other keys are already persisted in settings.json
        and can be picked up by a future extension of PreProcessing/MDBase.
        """
        flags: List[str] = []
        if "Temperature" in self.data:
            flags += ["-T", str(self.data["Temperature"])]
        if "Ensemble" in self.data:
            flags += ["-E", str(self.data["Ensemble"])]
        return flags


class MDGUI(tk.Tk):
    """Tkinter-based GUI skeleton.

    - Edits settings.json fields relevant for current PreProcessing.
    - Selects POSCAR file path (optional; if none, PreProcessing uses FCC Cu).
    - Runs simulation in background thread to avoid blocking the UI.
    - Offers a button to visualize the produced trajectory (MDBase.visualizeTraj via MD instance).
    """

    def __init__(self):
        super().__init__()
        self.title("Molecular Dynamics GUI - Skeleton")
        self.geometry("720x520")

        self.settings = SettingsModel()
        self.poscar_path: Optional[str] = None
        self.md_instance = None  # type: Optional[Any]  # MDBase returned from PreProcessing.createMD

        self._build_widgets()
        self._load_settings_into_form()

    # UI BUILDING
    def _build_widgets(self):
        # Settings frame
        settings_frame = ttk.LabelFrame(self, text="Settings (settings.json)")
        settings_frame.pack(fill="x", padx=10, pady=10)

        # Temperature
        ttk.Label(settings_frame, text="Temperature (K)").grid(row=0, column=0, sticky="w", padx=5, pady=5)
        self.var_temp = tk.StringVar()
        ttk.Entry(settings_frame, textvariable=self.var_temp, width=10).grid(row=0, column=1, padx=5, pady=5)

        # Ensemble
        ttk.Label(settings_frame, text="Ensemble").grid(row=0, column=2, sticky="w", padx=5, pady=5)
        self.var_ensemble = tk.StringVar()
        ttk.Combobox(settings_frame, textvariable=self.var_ensemble, values=["NVE", "NVT", "NPT"], width=8).grid(row=0, column=3, padx=5, pady=5)

        # Optional auxiliary parameters that are currently not wired by PreProcessing flags but live in settings.json
        # (These can be used by future extensions where PreProcessing/MDBase are expanded to read them.)
        ttk.Label(settings_frame, text="Timestep (fs)").grid(row=1, column=0, sticky="w", padx=5, pady=5)
        self.var_dt = tk.StringVar()
        ttk.Entry(settings_frame, textvariable=self.var_dt, width=10).grid(row=1, column=1, padx=5, pady=5)

        ttk.Label(settings_frame, text="Steps").grid(row=1, column=2, sticky="w", padx=5, pady=5)
        self.var_steps = tk.StringVar()
        ttk.Entry(settings_frame, textvariable=self.var_steps, width=10).grid(row=1, column=3, padx=5, pady=5)

        ttk.Label(settings_frame, text="Output file base").grid(row=2, column=0, sticky="w", padx=5, pady=5)
        self.var_output = tk.StringVar()
        ttk.Entry(settings_frame, textvariable=self.var_output, width=20).grid(row=2, column=1, padx=5, pady=5)

        ttk.Label(settings_frame, text="Potential").grid(row=2, column=2, sticky="w", padx=5, pady=5)
        self.var_potential = tk.StringVar()
        ttk.Combobox(settings_frame, textvariable=self.var_potential, values=["LJ", "EMT"], width=8).grid(row=2, column=3, padx=5, pady=5)

        # POSCAR selection
        poscar_frame = ttk.LabelFrame(self, text="Structure")
        poscar_frame.pack(fill="x", padx=10, pady=10)

        self.lbl_poscar = ttk.Label(poscar_frame, text="POSCAR: <not selected, will use default FCC Cu>")
        self.lbl_poscar.pack(side="left", padx=5, pady=5)

        ttk.Button(poscar_frame, text="Choose POSCAR...", command=self._choose_poscar).pack(side="left", padx=5)
        ttk.Button(poscar_frame, text="Clear", command=self._clear_poscar).pack(side="left", padx=5)

        # Actions
        action_frame = ttk.Frame(self)
        action_frame.pack(fill="x", padx=10, pady=10)

        ttk.Button(action_frame, text="Save Settings", command=self._save_settings_from_form).pack(side="left", padx=5)
        ttk.Button(action_frame, text="Run Simulation", command=self._run_simulation_async).pack(side="left", padx=5)
        ttk.Button(action_frame, text="View Trajectory", command=self._view_traj).pack(side="left", padx=5)

        # Log/Status
        log_frame = ttk.LabelFrame(self, text="Status / Log (see also data.log)")
        log_frame.pack(fill="both", expand=True, padx=10, pady=10)
        self.txt_log = tk.Text(log_frame, height=12)
        self.txt_log.pack(fill="both", expand=True)

    def _load_settings_into_form(self):
        self.var_temp.set(str(self.settings.data.get("Temperature", 293)))
        self.var_ensemble.set(str(self.settings.data.get("Ensemble", "NVT")))
        self.var_dt.set(str(self.settings.data.get("Timestep_fs", 10)))
        self.var_steps.set(str(self.settings.data.get("Steps", 200)))
        self.var_output.set(str(self.settings.data.get("OutputFile", "data")))
        self.var_potential.set(str(self.settings.data.get("Potential", "LJ")))

    # UI HANDLERS
    def _choose_poscar(self):
        path = filedialog.askopenfilename(title="Select POSCAR file", filetypes=[("POSCAR/CONTCAR", "POSCAR CONTCAR *"), ("All files", "*")])
        if path:
            self.poscar_path = path
            self.lbl_poscar.config(text=f"POSCAR: {path}")

    def _clear_poscar(self):
        self.poscar_path = None
        self.lbl_poscar.config(text="POSCAR: <not selected, will use default FCC Cu>")

    def _save_settings_from_form(self):
        # Basic validation and save
        try:
            self.settings.data["Temperature"] = float(self.var_temp.get())
        except ValueError:
            messagebox.showerror("Invalid temperature", "Temperature must be a number (K).")
            return
        ens = self.var_ensemble.get().strip().upper()
        if ens not in ("NVE", "NVT", "NPT"):
            messagebox.showerror("Invalid ensemble", "Ensemble must be one of NVE, NVT, NPT.")
            return
        self.settings.data["Ensemble"] = ens

        # Save auxiliary fields back to JSON (even if not used yet by PreProcessing)
        for key, var, caster in [
            ("Timestep_fs", self.var_dt, float),
            ("Steps", self.var_steps, int),
            ("OutputFile", self.var_output, str),
            ("Potential", self.var_potential, str),
        ]:
            try:
                self.settings.data[key] = caster(var.get())
            except Exception:
                # Keep previous value on failure
                pass

        self.settings.save()
        self._log("Settings saved to settings.json")

    def _run_simulation_async(self):
        # Save current settings to make sure PreProcessing picks them up
        self._save_settings_from_form()
        self._log("Starting simulation... (running in background)")
        t = threading.Thread(target=self._run_simulation, daemon=True)
        t.start()

    def _run_simulation(self):
        try:
            settings_path = self.settings.settings_path
            poscar_path = self.poscar_path or "POSCAR"  # If file exists, PreProcessing.readAtomicStructure can use it

            # Emulate CLI flags for PreProcessing.readTerminalInput
            flags = self.settings.to_flags()

            pp = PreProcessing(settings_path, poscar_path, flags)
            md = pp.createMD()  # MDBase instance configured by ensemble/temperature
            self.md_instance = md

            # Run simulation using atoms prepared in PreProcessing
            md.runMD(pp.atoms)
            self._log("Simulation completed. Output written to data.traj and data.log (or custom output base).")
        except Exception as e:
            self._log(f"Error during simulation: {e}")

    def _view_traj(self):
        if self.md_instance is None:
            # Attempt to construct a minimal instance for visualization if data.traj exists
            if not os.path.exists("data.traj"):
                messagebox.showinfo("No trajectory", "No MD instance available and data.traj not found.")
                return
            try:
                # Build a minimal PreProcessing/MD to reuse visualizeTraj method
                pp = PreProcessing(self.settings.settings_path, self.poscar_path or "POSCAR", self.settings.to_flags())
                md = pp.createMD()
                md.visualizeTraj()
            except Exception as e:
                self._log(f"Error opening trajectory: {e}")
        else:
            try:
                self.md_instance.visualizeTraj()
            except Exception as e:
                self._log(f"Error opening trajectory: {e}")

    def _log(self, msg: str):
        self.txt_log.insert("end", msg + "\n")
        self.txt_log.see("end")


def main():
    app = MDGUI()
    app.mainloop()


if __name__ == "__main__":
    main()