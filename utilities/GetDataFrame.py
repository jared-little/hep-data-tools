import ROOT
import os
from dotenv import load_dotenv

load_dotenv("histograms.env")
input_trees = os.getenv("INPUT_TREES")


def get_signal_df(Process="XHS_X4000_S2000", campaigns=None):
    """Get the dataframe for a given signal sample."""

    files = []
    for campaign in campaigns:
        sample_name = f"{campaign}_bbVV_{Process}_bbWW_allhad_PreselectionLoose.root"
        full_path = os.path.join(input_trees, sample_name)
        if not os.path.exists(full_path):
            raise FileNotFoundError(f"File not found: {full_path}")
        files.append(full_path)

    dataframe = ROOT.RDataFrame(f"bbVV_data", files)

    return dataframe


def get_background_df(Process="ttbar", campaigns=None):
    """Get the dataframe for a given background sample."""

    files = []
    for campaign in campaigns:
        sample_name = f"{campaign}_bbVV_{Process}_PreselectionLoose.root"
        full_path = os.path.join(input_trees, sample_name)
        if not os.path.exists(full_path):
            raise FileNotFoundError(f"File not found: {full_path}")
        files.append(full_path)

    dataframe = ROOT.RDataFrame(f"bbVV_data", files)

    return dataframe


if __name__ == "__main__":

    campaigns = ["mc23a", "mc23d", "mc23e"]  # mc23a, mc23d, mc23e
    df = get_signal_df(campaigns=campaigns)
    print(df.GetColumnNames())
