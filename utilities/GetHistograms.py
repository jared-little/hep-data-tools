import ROOT
import os
from dotenv import load_dotenv

load_dotenv("histograms.env")

input_folder = os.getenv("INPUT_FOLDER")

def _get_detached_histogram(file_path, hist_path, clone_suffix):
    """Load a histogram from a ROOT file and detach it from file ownership."""

    root_file = ROOT.TFile.Open(file_path)
    if not root_file or root_file.IsZombie():
        raise FileNotFoundError(f"Could not open ROOT file: {file_path}")

    hist = root_file.Get(hist_path)
    if not hist:
        root_file.Close()
        raise KeyError(f"Histogram not found: {hist_path} in {file_path}")

    detached = hist.Clone(f"{hist.GetName()}_{clone_suffix}")
    detached.SetDirectory(0)
    root_file.Close()

    return detached


def get_signal_histogram(Signal="XHS_X4000_S2000", Var="NN_score", Region="Preselection", Rebin=1, campaigns=None):
    """Get the signal histogram for a given signal, variable, region, and rebinning factor."""

    if campaigns is None:
        campaigns = ["mc23a"]
    if not input_folder:
        raise RuntimeError("INPUT_FOLDER is not set. Check histograms.env")

    hists = []
    for campaign in campaigns:
        file_path = os.path.join(input_folder, f"{campaign}_{Signal}_bbWW_allhad.root")
        hist_path = f"{Region}/bbVVSplitHadAnalysis_13p6TeV_{Signal}_bbWW_allhad/{Var}"
        hists.append(_get_detached_histogram(file_path, hist_path, campaign))
    
    # Sum histograms from different campaigns
    if not hists:
        raise RuntimeError(f"No histograms loaded for signal {Signal}")

    sig_histogram = hists[0].Clone(f"{Signal}_{Region}_{Var}_sum")
    sig_histogram.SetDirectory(0)
    for hist in hists[1:]:
        sig_histogram.Add(hist)

    if "2000" in Signal: sig_histogram.SetLineColor(ROOT.kOrange)
    if "4000" in Signal: sig_histogram.SetLineColor(ROOT.kCyan)
    if "6000" in Signal: sig_histogram.SetLineColor(ROOT.kViolet)

    sig_histogram.SetLineWidth(4)
    sig_histogram.SetLineStyle(2)
    sig_histogram.SetDirectory(0)
    sig_histogram.Rebin(Rebin)

    return sig_histogram


def get_bkg_histogram(Bkg="dijet", Var="NN_score", Region="Preselection", Rebin=1, campaigns=None):
    """Get the background histogram for a given background, variable, region, and rebinning factor."""

    if campaigns is None:
        campaigns = ["mc23a"]
    if not input_folder:
        raise RuntimeError("INPUT_FOLDER is not set. Check histograms.env")

    hists = []
    for campaign in campaigns:
        file_path = os.path.join(input_folder, f"{campaign}_{Bkg}.root")
        hist_path = f"{Region}/bbVVSplitHadAnalysis_13p6TeV_{Bkg}/{Var}"
        hists.append(_get_detached_histogram(file_path, hist_path, campaign))

    if not hists:
        raise RuntimeError(f"No histograms loaded for background {Bkg}")

    bkg_histogram = hists[0].Clone(f"{Bkg}_{Region}_{Var}_sum")
    bkg_histogram.SetDirectory(0)
    for hist in hists[1:]:
        bkg_histogram.Add(hist)

    bkg_histogram.Rebin(Rebin)
    bkg_histogram.GetYaxis().SetTitle("Events")
    bkg_histogram.GetYaxis().SetLabelSize(0.05)
    bkg_histogram.GetYaxis().SetTitleSize(0.1)
    bkg_histogram.GetYaxis().SetTitleOffset(0.6)
    bkg_histogram.SetDirectory(0)

    return bkg_histogram


def get_data_histogram(Var="NN_score", region="Preselection", rebin=1, campaigns=None):
    """Get the data histogram for a given variable, region, rebinning factor, and campaigns."""
    
    if campaigns is None:
        campaigns = ["22"]
    if not input_folder:
        raise RuntimeError("INPUT_FOLDER is not set. Check histograms.env")

    hists = []
    for campaign in campaigns:
        file_path = os.path.join(input_folder, f"data{campaign}.root")
        hist_path = f"{region}/bbVVSplitHadAnalysis_13p6TeV_data/{Var}"
        hists.append(_get_detached_histogram(file_path, hist_path, campaign))

    if not hists:
        raise RuntimeError("No histograms loaded for data")
    
    data_histogram = hists[0].Clone(f"data_{region}_{Var}_sum")
    data_histogram.SetDirectory(0)
    for hist in hists[1:]:
        data_histogram.Add(hist)

    data_histogram.Rebin(rebin)
    data_histogram.SetMarkerStyle(20)
    # data_histogram.SetMarkerSize(1.2)
    data_histogram.SetDirectory(0)

    return data_histogram


def get_var_name(Var):
    """
    Get the x-axis title for a given variable.
    """
    var_name = {
        "NN_score": "NN Score",
        "largeRjetpt_1": "Leading Large-R UFO Jet p_{T}",
        "largeRjetpt_2": "Subleading Large-R UFO Jet p_{T}",
        "largeRjetpt_3": "Third Leading Large-R UFO Jet p_{T}",
        "largeRjetm_1": "Leading Large-R UFO Jet Mass",
        "largeRjetm_2": "Subleading Large-R UFO Jet Mass",
        "largeRjetm_3": "Third Leading Large-R UFO Jet Mass",
    }

    return var_name[Var]
