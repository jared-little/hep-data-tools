import ROOT
from utilities.ComputeSignificance import get_Zn_histogram, get_SB_histogram
import os
from dotenv import load_dotenv
ROOT.gROOT.SetStyle("ATLAS")

load_dotenv("histograms.env")
input_trees = os.getenv("NOTRIG_TREES")
# SCRIPT FOR CHECKING SIGNAL EFFICIENCY
# INPUT IS EASYJET NTUPLES

def get_signal_df(Process="X4000_S2000", campaigns: list[str] | None = None):
    """Get the dataframe for a given signal sample."""

    if campaigns is None:
        raise ValueError("campaigns must be provided")

    files = []
    for campaign in campaigns:
        sample_name = f"ejOutput_PHYS_bbVV_0lep_splitboosted_{Process}_{campaign}.root"
        full_path = os.path.join(input_trees, sample_name)
        if not os.path.exists(full_path):
            raise FileNotFoundError(f"File not found: {full_path}")
        files.append(full_path)

    dataframe = ROOT.RDataFrame(f"AnalysisMiniTree", files)

    return dataframe


def make_canvas(name ="canvas"):
    """Make a canvas for plotting."""
    left_margin = 0.12
    right_margin = 0.03

    canvas = ROOT.TCanvas(name, name, 700, 700)

    return canvas


def set_bins():
    """Set the binning for the histograms, based on the variable being plotted."""
    bins = {
        "largeRjetpt_1": (100, 0, 2000), # Trigger turn-on is around 500 GeV, so start there
    }
    return bins


def new_columns(df):
    """Define new columns in the RDataFrame for the variables we want to plot."""
    df = df.Define("largeRjetpt_1", "recojet_antikt10UFO_pt___NOSYS[0] / 1000")

    return df


def make_histogram(df, Var, selections=None):
    """Make a histogram of a given variable from a given RDataFrame.
    Apply selections here as well for optimization."""

    bins = set_bins()[Var]
    df = new_columns(df)

    if selections is not None:
        for sel in selections:
            df = df.Filter(sel)

    hist = df.Histo1D((f"{Var}", f"{Var}", bins[0], bins[1], bins[2]), Var)

    return hist


def calculate_yields(df, Var, selections=None):
    """Calculate the yields for a given variable and selections."""
    df = new_columns(df)

    if selections is not None:
        for sel in selections:
            df = df.Filter(sel)

    yield_ = df.Count().GetValue()

    return yield_


def combine_histograms():
    """Combine histograms for different signals to check efficiency."""
    sig_names = ["X1000_S500", "X2000_S1000", "X4000_S2000"]
    campaigns = ["mc23a", "mc23e"]
    for sig in sig_names:
        for campaign in campaigns:
            canvas = make_canvas(f"canvas_{sig}_{campaign}")
            df = get_signal_df(sig, [campaign])
            total_events = calculate_yields(df, "largeRjetpt_1")
            passing_events = calculate_yields(df, "largeRjetpt_1", ["largeRjetpt_1 > 520"])
            passing_rate = (passing_events / total_events) if total_events else 0.0
            print(
                f"{sig} {campaign}: {passing_events}/{total_events} events "
                f"with largeRjetpt_1 > 520 GeV ({passing_rate:.2%})"
            )
            hist = make_histogram(df, "largeRjetpt_1")
            # hist.Scale(1/hist.Integral())
            hist.SetLineColor(ROOT.kRed if sig == sig_names[0] else ROOT.kBlue)
            hist.SetLineWidth(2)
            hist.Draw("hist same")
            canvas.SaveAs(f"largeRjetpt_1_{sig}_{campaign}.pdf")


if __name__ == "__main__":
    ROOT.gROOT.SetBatch(True)

    combine_histograms()
