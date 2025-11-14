import os
import httk
import httk.db
from httk.atomistic import Structure
from classes import DefectCell

def get_defect_cif(host_name, defect_name, db_path="../../defect/defects.sqlite", out_dir="cif_out"):
    os.makedirs(out_dir, exist_ok=True)
    backend = httk.db.backend.Sqlite(db_path)
    store = httk.db.store.SqlStore(backend)
    search = store.searcher()

    search_defect = search.variable(DefectCell)
    host_name = host_name+"_PBE"
    search.add(search_defect.host_name == host_name)
    #search.add(search_defect.defect_name.startswith(defect_name) )

    search.output(search_defect, "defect")

    saved_paths = []
    i = 0

    for match, _ in search:
        defect_cell = match[0]

        if not defect_cell.defect_name.startswith(defect_name):
            continue

        struct = defect_cell.defect_structure

        # here we use the actual defect name (C_int0, C_int1,osv)
        filename = f"{host_name}_{defect_cell.defect_name}.cif"
        full_path = os.path.abspath(os.path.join(out_dir, filename))

        httk.save(struct, full_path)
        print(f"CIF saved: {full_path}")
        saved_paths.append(full_path) #dont need now, maybe later
        i += 1

    backend.close()
    return saved_paths


if __name__ == "__main__":
    get_defect_cif(
        host_name="As2",
        defect_name="C_int",
        db_path="../../defect/defects.sqlite",
        out_dir="cif_out"
    )