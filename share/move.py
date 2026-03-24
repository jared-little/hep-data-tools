#!/usr/bin/env python3

from pathlib import Path
from glob import glob
import shutil


def main(ntuple = False, histogram = True) -> None:
    """Move ntuples and histograms from the "Exports" and "Hists" folders
       to the "ntuples" or "histograms" folders."""

    campaigns = ["mc23a", "mc23d", "mc23e"]  # mc23a, mc23d, mc23e

    for campaign in campaigns:
        print(f"Processing campaign: {campaign}")
        if ntuple == True:
            in_path_glob = f"/Users/jlittle/work/HHbbVV/bbvv-presel-loose/{campaign}_sysv03_nominal/Exports/*"
            out_path = Path("/Users/jlittle/work/HHbbVV/bbvv-presel-loose/ntuples/")
        elif histogram == True:
            in_path_glob = f"/Users/jlittle/work/HHbbVV/bbvv-presel-loose/{campaign}_sysv03_nominal/Hists/*"
            out_path = Path("/Users/jlittle/work/HHbbVV/bbvv-presel-loose/histograms/")
        else:
            print("Please specify whether you want to move ntuples or histograms by setting the appropriate flag to True.")

        for src_path in glob(in_path_glob):
            src = Path(src_path)
            name = src.name

            dest = out_path / f"{campaign}_{name}"

            if str(src) != str(dest):
                shutil.copy2(src, dest)
                # print(f"source: {src}")
                # print(f"destination: {dest}")
            else:
                print(f"source and destination are the same: {src}")


if __name__ == "__main__":
    main(ntuple = False, histogram = True)
