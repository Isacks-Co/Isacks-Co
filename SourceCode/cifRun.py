#!/usr/bin/env python
from __future__ import print_function

import os
import httk
import httk.db
import httk.task
from httk.atomistic import Structure
from classes import DefectCell


def get_defect_cif(
    host_name,
    defect_name,
    db_path="../../defect/defects.sqlite",
    out_root="Runs",
    template_name="template",  #  httk-task-template
    ):

    os.makedirs(out_root, exist_ok=True)

    backend = httk.db.backend.Sqlite(db_path)
    store = httk.db.store.SqlStore(backend)
    search = store.searcher()

    search_defect = search.variable(DefectCell)

    host_name_pbe = host_name + "_PBE"
    search.add(search_defect.host_name == host_name_pbe)

    search.output(search_defect, "defect")
    saved_paths = []

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
        print(f"CIF saved: {full_path}")

        saved_paths.append(full_path)

    backend.close()
    return saved_paths


if __name__ == "__main__":
    get_defect_cif(
        host_name="As2",
        defect_name="C_int",
        db_path="../../defect/defects.sqlite",
        out_root="Runs",      
        template_name="template",
    )
