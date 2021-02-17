from hhsuite import HHResult# import hhResult class
import os
from pathlib import Path


def analyze(file_path):
    good_matches=[]
    hhresult_object = HHResult(file_path)
    hhresult_object.parse_result()
    for match in hhresult_object.matches:
        if float(match.probability)> 0.9:
            good_matches.append(match.query_id)
    return good_matches




def main():
    og=Path.cwd()
    p = og/"HHResult_package"/"HHResult_package"/"D29_HHsearch_results"
    for files in p.iterdir():
        analyze(files)

if __name__ == '__main__':
    main()
