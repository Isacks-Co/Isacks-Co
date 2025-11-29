import httk
import httk.db
from classes import DefectCell
from newDefect import newScreenResult


def getAllDefectKeys(store):

    search = store.searcher()
    search_defectcell = search.variable(DefectCell)

    search.output(search_defectcell.host_name, "host")
    search.output(search_defectcell.defect_name, "defect")

    all_keys = set()

    for match, header in search:
        host_name = match[0]
        defect_name = match[1]

        all_keys.add((host_name, defect_name))

    return all_keys


def getMDKeys(store):
    search = store.searcher()
    search_defectcell = search.variable(newScreenResult)

    search.output(search_defectcell.host_name, "host")
    search.output(search_defectcell.defect_name, "defect")

    all_md_keys = set()

    for match, header in search:
        host_name = match[0]
        defect_name = match[1]

        all_md_keys.add((host_name, defect_name))

    return all_md_keys

def allDoneMDRuns(store):
    defect_keys = getAllDefectKeys(store)
    md_keys = getMDKeys(store)

    missing = defect_keys - md_keys #ger tillbaka lista/set med de nycklar som saknas
    if missing:
        print("MISSING MD FOR: ")
        for host, defect in sorted(missing):
            print(host, defect)
        return False
    
    return True






