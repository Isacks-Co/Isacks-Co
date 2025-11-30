import httk
import httk.db
from classes import DefectCell
from newDefect import newScreenResult


def getAllDefectKeys(store):
    search = store.searcher()
    search_defectcell = search.variable(DefectCell)

    #smarter to use keys
    #search.output(search_defectcell.host_name, "host")
    #search.output(search_defectcell.defect_name, "defect")
    search.output(search_defectcell.key, "defect_key")

    all_keys = set()
    for match, header in search:
        #host_name = match[0]
        #defect_name = match[1]
        key = match[0]

        #all_keys.add((host_name, defect_name))
        all_keys.add(key)

    return all_keys


def getMDKeys(store):
    search = store.searcher()
    search_MDcell = search.variable(newScreenResult)

    #search.output(search_defectcell.host_name, "host")
    #search.output(search_defectcell.defect_name, "defect")
    search.output(search_MDcell.defect_key, "defect_key")

    all_md_keys = set()
    for match, header in search:
        #host_name = match[0]
        #defect_name = match[1]
        key = match[0]

        all_md_keys.add(key)

    return all_md_keys

def allDoneMDRuns(store):
    defect_keys = getAllDefectKeys(store)
    md_keys = getMDKeys(store)

    missing_keys = defect_keys - md_keys #ger tillbaka lista/set med de nycklar som saknas
    if missing_keys:
        print(f"These keys are missing: {sorted(missing_keys)} ")
        return False

    return True






