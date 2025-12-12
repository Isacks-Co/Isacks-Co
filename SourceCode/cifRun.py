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
import httk.task
import os
from classes import DefectCell
from httk.atomistic import Structure


def get_defect_cif(
        host_name,
        defect_name,
        db_path="../../defect/defects.sqlite",
        out_root="Runs",
        template_name="template",  # httk-task-template
):
    """
    Create HTTK run folders and write CIF files for matching defect structures.

    Parameters
    ----------
    host_name : str
        Base host identifier without the ``"_PBE"`` suffix (e.g. ``"As2"``).
    defect_name : str
        Defect name prefix to match (e.g. ``"C_"``).
    db_path : str, optional
        Path to the SQLite defect database.
    out_root : str, optional
        Root directory where run folders will be created.
    template_name : str, optional
        HTTK task template name passed to ``httk.task.create_batch_task``.

    Returns
    -------
    None

    Side Effects
    ------------
    - Creates directories under `out_root`.
    - Writes one CIF file per matching defect inside its run directory.
    - Creates HTTK batch task directories via ``create_batch_task``.
    """
    os.makedirs(out_root, exist_ok=True)

    backend = httk.db.backend.Sqlite(db_path)
    store = httk.db.store.SqlStore(backend)
    search = store.searcher()

    search_defect = search.variable(DefectCell)

    host_name_pbe = host_name + "_PBE"
    search.add(search_defect.host_name == host_name_pbe)

    search.output(search_defect, "defect")

    for match, header in search:
        defect_cell = match[0]

        if not defect_cell.defect_name.startswith(defect_name):
            continue

        struct = defect_cell.defect_structure

        run_name = f"{host_name_pbe}_{defect_cell.defect_name}"

        dir_path = httk.task.create_batch_task(
            out_root,
            template_name,
            {"structure": struct},
            name=run_name,
            overwrite_head_dir=True,
        )

        # CIF-file gets put in folder that create_batch_task made
        cif_filename = f"{run_name}.cif"
        full_path = os.path.abspath(os.path.join(dir_path, cif_filename))

        httk.save(struct, full_path)

    backend.close()


if __name__ == "__main__":
    get_defect_cif(
        host_name="As2",
        defect_name="C_",
        db_path="../../defect/defects.sqlite",
        out_root="Runs",
        template_name="template",
    )
