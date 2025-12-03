#!/usr/bin/env python
import os
import sys

import httk
import httk.db
import httk.task
from classes import DefectCell
from abad_classes import AbadParameters
from httk.atomistic import Structure


def get_defect_cif(
        host_name,
        defect_name,
        db_path="../../defect/defects.sqlite",
        out_root="Runs",
        template_name="template",  # httk-task-template
):
    os.makedirs(out_root, exist_ok=True)

    backend = httk.db.backend.Sqlite(db_path)
    store = httk.db.store.SqlStore(backend)
    search = store.searcher()

    search_defect = search.variable(DefectCell)
    search_abad = search.variable(AbadParameters)
    host_name_pbe = host_name + "_PBE"
    search.add(search_defect.host_name == host_name_pbe)
    search.add(search_defect.key == search_abad.key)

    search.output(search_defect, "defect")
    search.output(search_abad, "abad")

    for match, header in search:
        defect_cell = match[0]
        abad_parameters = match[1]

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

        key_file_path = dir_path + "/key"
        defect_info_file_path = dir_path + "/defect_info"

        with open(defect_info_file_path, "w", encoding="utf-8") as f:
            f.write(str(abad_parameters.defect_index))

        with open(key_file_path, "w", encoding="utf-8") as f:
            f.write(str(defect_cell.key))

        
        httk.save(struct, full_path)

    backend.close()


if __name__ == "__main__":
    get_defect_cif(
        host_name="MoTe2",
        defect_name="C_",
        db_path="../../defect/defects.sqlite",
        out_root="../../MD_runs/Runs",
        template_name="template",
    )
