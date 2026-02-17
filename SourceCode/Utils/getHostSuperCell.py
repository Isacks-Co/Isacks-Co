# MIT License
#
# Copyright (c) 2025 Isacks-Co contributors
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.


import httk
import httk.db
import os
from classes import DefectCell,HostSuperCell
from httk.atomistic import Structure
from httk.external import ase_glue
from scipy.constants import precision
from mace.calculators import MACECalculator
import torch
import json
from ase.io import read



def get_defect_cif(material, defect_name, db_path="../../../defect/defects.sqlite", out_dir="cif_out"):

    os.makedirs(out_dir, exist_ok=True)
    print(db_path)
    backend = httk.db.backend.Sqlite(db_path)
    store = httk.db.store.SqlStore(backend)
    search = store.searcher()

    search_supercell = search.variable(HostSuperCell)
    #search.add(search_supercell.material == material)
    search.output(search_supercell, "host_materials")

    saved_paths = []
    i = 0

    dict = {}
    for match, _ in search:
        material_key = match[0].material
        httk.save( match[0].host_supercell, f"{material_key}.cif" )
        atoms = read(f"{material_key}.cif")
        total_energy = atomObject_to_energy(atoms, "../../../mace-mpa-0-medium.model")
        dict[material_key] = total_energy

    with open("HostEnergyCell.json", "w") as o:
        json.dump(dict, o)



    backend.close()
    return saved_paths

def atomObject_to_energy(atom_struct, model_path ):
    device = "cuda" if torch.cuda.is_available() else "cpu"
    atom_struct.calc = MACECalculator(device = device, model_paths=model_path, default_dtype="float64",
        head="default")
    energy = atom_struct.get_total_energy()

    return energy



if __name__ == "__main__":
    get_defect_cif(
        material="PtSe2",
        defect_name="C_int",
        out_dir="cif_out"
    )












