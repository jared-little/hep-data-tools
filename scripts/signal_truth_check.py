import ROOT
from utilities.ComputeSignificance import get_Zn_histogram, get_SB_histogram
import os
from dotenv import load_dotenv
ROOT.gROOT.SetStyle("ATLAS")

load_dotenv("histograms.env")
input_trees = os.getenv("INPUT_TREES")
# SCRIPT FOR CHECKING SIGNAL TRUTHMX

def get_signal_df(Process="X4000_S2000"):
    """Get the dataframe for a given signal sample."""

    files = []

    sample_name = f"bbVV_XHS_{Process}_bbWW_allhad_All.root"
    full_path = os.path.join(input_trees, sample_name)
    if not os.path.exists(full_path):
        raise FileNotFoundError(f"File not found: {full_path}")
    files.append(full_path)

    dataframe = ROOT.RDataFrame(f"bbVV_data", files)

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
        "largeRjetpt_1": (100, 500, 2500),
        "largeRjetpt_2": (100, 0, 2000),
        "largeRjetpt_3": (100, 0, 2000),
        "largeRjeteta_1": (10, -2.0, 2.0),
        "largeRjeteta_2": (10, -2.0, 2.0),
        "largeRjeteta_3": (10, -2.0, 2.0),
        "largeRjetphi_1": (10, -3.14, 3.14),
        "largeRjetphi_2": (10, -3.14, 3.14),
        "largeRjetphi_3": (10, -3.14, 3.14),
        "largeRjetm_1": (20, 0, 300),
        "largeRjetm_2": (20, 0, 300),
        "largeRjetm_3": (20, 0, 300),
        "truth_mX_gev": (50, 0, 5000),
    }
    return bins


def new_columns(df):
    """Define new columns in the RDataFrame for the variables we want to plot."""
    df = df.Define("largeRjetpt_1", "largeRjetpt[0] / 1000")
    df = df.Define("largeRjetpt_2", "largeRjetpt[1] / 1000")
    df = df.Define("largeRjetpt_3", "largeRjetpt[2] / 1000")
    df = df.Define("largeRjeteta_1", "largeRjeteta[0]")
    df = df.Define("largeRjeteta_2", "largeRjeteta[1]")
    df = df.Define("largeRjeteta_3", "largeRjeteta[2]")
    df = df.Define("largeRjetphi_1", "largeRjetphi[0]")
    df = df.Define("largeRjetphi_2", "largeRjetphi[1]")
    df = df.Define("largeRjetphi_3", "largeRjetphi[2]")
    df = df.Define("largeRjetm_1", "largeRjetm[0] / 1000")
    df = df.Define("largeRjetm_2", "largeRjetm[1] / 1000")
    df = df.Define("largeRjetm_3", "largeRjetm[2] / 1000")
    df = df.Define("truth_mX_gev", "truth_mX / 1000")

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


def make_histograms(variable = "largeRjetpt_1", signal = "X4000_S2000"):
    """Make histograms for different signals."""

    canvas = make_canvas(f"canvas_{signal}")
    df = get_signal_df(signal)
    total_events = calculate_yields(df, variable)

    hist = make_histogram(df, variable)

    hist.SetLineColor(ROOT.kRed)
    hist.SetLineWidth(2)
    hist.Draw("hist same")
    canvas.SaveAs(f"plots/signal_plots/{variable}_{signal}.pdf")


def make_histograms_truth_comparison(variable = "largeRjetpt_1", signal = "X4000_S2000"):
    """Make histograms for different signals."""

    canvas = make_canvas(f"canvas_{signal}")
    df = get_signal_df(signal)
    legend = ROOT.TLegend(0.6, 0.7, 0.9, 0.9)
    legend.SetBorderSize(0)

    hist_notruth = make_histogram(df, variable, selections=["truth_mX == -666"])
    hist_truth = make_histogram(df, variable, selections=["truth_mX != -666"])

    hist_notruth.SetLineColor(ROOT.kRed)
    hist_notruth.SetLineWidth(2)
    if hist_notruth.Integral() != 0: hist_notruth.Scale(1 / hist_notruth.Integral()) # Normalize to unit area
    hist_truth.SetLineColor(ROOT.kBlue)
    hist_truth.SetLineWidth(2)
    if hist_truth.Integral() != 0: hist_truth.Scale(1 / hist_truth.Integral()) # Normalize to unit area
    legend.AddEntry(hist_notruth.GetPtr(), "Truth mX does not exist", "l")
    legend.AddEntry(hist_truth.GetPtr(), "Truth mX exists", "l")
    maximum = 0
    for i in range(1, hist_truth.GetNbinsX() + 1):
        if max(hist_truth.GetBinContent(i), hist_notruth.GetBinContent(i)) > maximum:
            maximum = max(hist_truth.GetBinContent(i), hist_notruth.GetBinContent(i))
    hist_truth.SetMaximum(2 * maximum)
    hist_truth.GetXaxis().SetTitle(variable)
    hist_truth.Draw("hist same")
    hist_notruth.Draw("hist same")
    legend.Draw()
    canvas.SaveAs(f"plots/signal_plots/truth_comparison_{variable}_{signal}.pdf")


if __name__ == "__main__":
    ROOT.gROOT.SetBatch(True)

    signal_names = ["X4000_S2000"]
    variables = ["largeRjetpt_1", "largeRjetpt_2", "largeRjetpt_3", "truth_mX_gev",
                 "largeRjeteta_1", "largeRjeteta_2", "largeRjeteta_3",
                 "largeRjetphi_1", "largeRjetphi_2", "largeRjetphi_3",
                 "largeRjetm_1", "largeRjetm_2", "largeRjetm_3"]
    for var in variables:
        for signal in signal_names:
            make_histograms_truth_comparison(var, signal)
