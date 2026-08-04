"""Data-driven ABCD estimate of the dijet background.

In every region the dijet background is taken as (data - non-dijet MC).
The CRs give a per-bin transfer factor, which predicts the SR and is validated in the VRs:

    TF  = CR1 / CR2
    SR  = CR0 x TF          (prediction)
    VR1 = VR2 x TF          (closure test)

Outputs a ROOT file with one directory per region (TRExFitter/HHARD layout) plus
transfer-factor and closure plots in plots/ABCD/.
"""

import os
import ctypes

import ROOT

from utilities.GetHistograms import get_var_name
from utilities.DijetEstimate import (
    ABCD_REGIONS as REGIONS,
    ANALYSIS_BINNING,
    PREDICTED_REGIONS,
    apply_transfer_factor,
    check_positive_bins,
    get_dijet_histogram,
    get_transfer_factor,
)
from utilities.DijetEstimate import rebin_histogram as _rebin_histogram

ROOT.gROOT.SetStyle("ATLAS")


def rebin_histogram(hist, binning, name):
    """Rebin a histogram to the analysis binning, refusing to go on if the subtraction
    left a negative bin behind."""

    hist_rebin = _rebin_histogram(hist, binning, name)
    check_positive_bins(hist_rebin, name)

    return hist_rebin


def make_canvas_pads(can_name="ABCD"):
    """Make a canvas with two pads, one for the main plot and one for the ratio plot."""

    c = ROOT.TCanvas(can_name, can_name, 700, 600)
    c.cd()
    pad1 = ROOT.TPad(can_name + "_pad1", can_name + "_pad1", 0., 0.305, .99, 1)
    pad1.SetLeftMargin(0.12)
    pad1.SetRightMargin(0.03)
    pad1.SetBottomMargin(0.0)
    pad1.SetFillColor(ROOT.kWhite)
    pad1.SetTickx()
    pad1.SetTicky()
    pad1.Draw()

    pad2 = ROOT.TPad(can_name + "_pad2", can_name + "_pad2", 0., 0.01, .99, 0.295)
    pad2.SetLeftMargin(0.12)
    pad2.SetRightMargin(0.03)
    pad2.SetBottomMargin(0.38)
    pad2.SetFillColor(ROOT.kWhite)
    pad2.SetTickx()
    pad2.SetTicky()
    pad2.Draw()

    return c, pad1, pad2


def make_ratio_plot(hist1, hist2, Var):
    """Make a ratio plot of hist1 / hist2, with appropriate axis labels and styling."""

    h_ratio = hist1.Clone(f"h_ratio_{Var}")
    h_ratio.SetDirectory(0)
    h_ratio.Divide(hist2)

    h_ratio.GetXaxis().SetTitle(get_var_name(Var))
    h_ratio.GetXaxis().SetLabelSize(0.13)
    h_ratio.GetXaxis().SetLabelOffset(0.02)
    h_ratio.GetXaxis().SetTitleSize(0.15)
    h_ratio.GetYaxis().SetRangeUser(0.0, 2.0)
    h_ratio.GetYaxis().SetNdivisions(505)
    h_ratio.GetYaxis().SetLabelSize(0.13)
    h_ratio.GetYaxis().SetTitleSize(0.17)
    h_ratio.GetYaxis().SetTitleOffset(0.36)

    return h_ratio


def plot_transfer_factor(hist_tf, Var, outdir):
    """Plot the per-bin transfer factor."""

    c = ROOT.TCanvas(f"TF_{Var}", f"TF_{Var}", 700, 600)
    hist_tf.GetXaxis().SetTitle(get_var_name(Var))
    hist_tf.GetYaxis().SetTitle("Transfer factor (CR1 / CR2)")
    hist_tf.SetMarkerStyle(20)
    hist_tf.SetLineColor(ROOT.kBlack)
    hist_tf.Draw("E")

    c.SaveAs(f"{outdir}/TF_{Var}.pdf")


def plot_closure(hist_vr1, hist_prediction, Var, outdir):
    """Plot VR1 against its prediction VR2 x TF, with a ratio panel underneath."""

    hist_vr1.SetLineColor(ROOT.kBlack)
    hist_vr1.SetMarkerStyle(20)
    hist_prediction.SetLineColor(ROOT.kRed)
    hist_prediction.SetMarkerStyle(24)

    c, pad1, pad2 = make_canvas_pads(f"Closure_{Var}")

    pad1.cd()
    hist_vr1.GetYaxis().SetTitle("Events")
    hist_vr1.GetYaxis().SetLabelSize(0.05)
    hist_vr1.GetYaxis().SetTitleSize(0.06)
    hist_vr1.GetYaxis().SetTitleOffset(0.9)
    hist_vr1.SetMinimum(0.0)
    hist_vr1.SetMaximum(1.6 * max(hist_vr1.GetMaximum(), hist_prediction.GetMaximum()))
    hist_vr1.Draw("E")
    hist_prediction.Draw("E same")

    leg = ROOT.TLegend(0.6, 0.7, 0.9, 0.9, "")
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.AddEntry(hist_vr1, "#font[42]{VR1}", "lep")
    leg.AddEntry(hist_prediction, "#font[42]{VR2 #times TF}", "lep")
    leg.Draw()

    pad2.cd()
    h_ratio = make_ratio_plot(hist_vr1, hist_prediction, Var)
    h_ratio.GetYaxis().SetTitle("VR1 / pred.")
    h_ratio.Draw("E")
    line = ROOT.TLine(h_ratio.GetXaxis().GetXmin(), 1, h_ratio.GetXaxis().GetXmax(), 1)
    line.SetLineStyle(7)
    line.Draw("same")
    ROOT.gPad.RedrawAxis()

    c.SaveAs(f"{outdir}/Closure_{Var}.pdf")


def print_closure_yields(hist_vr1, hist_prediction, Var):
    """Print the observed and predicted VR yields, including overflow."""

    error_vr1 = ctypes.c_double(0.0)
    error_prediction = ctypes.c_double(0.0)
    yield_vr1 = hist_vr1.IntegralAndError(0, hist_vr1.GetNbinsX() + 1, error_vr1)
    yield_prediction = hist_prediction.IntegralAndError(0, hist_prediction.GetNbinsX() + 1, error_prediction)

    print(f"  {Var} VR1        : {yield_vr1:10.2f} +/- {error_vr1.value:.2f}")
    print(f"  {Var} VR2 x TF   : {yield_prediction:10.2f} +/- {error_prediction.value:.2f}")


def write_histogram(outfile, hist, Region, name, Syst=None):
    """Write one histogram into <Region>/ (or systs/<Syst>/<Region>/) of the output file,
    following the HHARD/TRExFitter layout."""

    parts = ["systs", Syst, Region] if Syst else [Region]
    directory = outfile
    for part in parts:
        # outfile.mkdir(path, "", True) is supposed to return an existing directory instead of
        # failing, but that only works for single-level names -- for a repeated nested path
        # (e.g. writing several histograms into the same systs/<Syst>/<Region>/) it errors out
        # and returns a null pointer instead. Build the path one level at a time to avoid that.
        next_directory = directory.GetDirectory(part)
        if not next_directory:
            next_directory = directory.mkdir(part)
        directory = next_directory

    directory.cd()
    hist.Write(name)
    outfile.cd()


def run_abcd_variable(Var, binning, campaigns, outfile, outdir, Syst=None):
    """Run the ABCD estimate for one variable, either nominal (Syst=None) or one systematic variation."""

    print(f"\n=== {Var}" + (f" [{Syst}]" if Syst else "") + " ===")

    # 1. dijet = data - MC in every region, at input binning and at analysis binning
    hists = {Region: get_dijet_histogram(Var, Region, campaigns=campaigns, Syst=Syst) for Region in REGIONS}
    hists_rebin = {Region: rebin_histogram(hist, binning, f"{Var}_{Region}_rebin")
                   for Region, hist in hists.items()}

    for Region in REGIONS:
        write_histogram(outfile, hists[Region], Region, Var, Syst)
        write_histogram(outfile, hists_rebin[Region], Region, f"{Var}_rebin", Syst)

    # 2. transfer factor from the control regions
    hist_tf = get_transfer_factor(hists_rebin["CR1"], hists_rebin["CR2"], f"TF_{Var}")

    # 3. regions the TF can predict as well as measure, written alongside the measurement so a
    # plot can put one against the other. See PREDICTED_REGIONS for what each one is worth.
    predictions = {Region: apply_transfer_factor(hists[source], hist_tf, f"{Var}_{Region}_prediction")
                   for Region, source in PREDICTED_REGIONS.items()}
    predictions_rebin = {Region: apply_transfer_factor(hists_rebin[source], hist_tf, f"{Var}_{Region}_prediction_rebin")
                         for Region, source in PREDICTED_REGIONS.items()}

    for Region in PREDICTED_REGIONS:
        write_histogram(outfile, predictions[Region], Region, f"{Var}_prediction", Syst)
        write_histogram(outfile, predictions_rebin[Region], Region, f"{Var}_prediction_rebin", Syst)

    # 4. closure test in the validation regions
    hist_prediction = predictions_rebin["VR1"]
    print_closure_yields(hists_rebin["VR1"], hist_prediction, Var)

    # Plots are only useful for the nominal estimate
    if Syst is None:
        plot_transfer_factor(hist_tf, Var, outdir)
        plot_closure(hists_rebin["VR1"], hist_prediction, Var, outdir)

    # 5. signal region prediction, written at both binnings so that <Var> means the input
    # binning and <Var>_rebin the analysis binning in every region, the SR included. The
    # coarse TF is applied bin by bin to the finer CR0, since measuring it at the input
    # binning would just be noise.
    hist_sr = apply_transfer_factor(hists["CR0"], hist_tf, f"{Var}_SR")
    hist_sr_rebin = apply_transfer_factor(hists_rebin["CR0"], hist_tf, f"{Var}_SR_rebin")
    write_histogram(outfile, hist_sr, "SR", Var, Syst)
    write_histogram(outfile, hist_sr_rebin, "SR", f"{Var}_rebin", Syst)


def run_abcd(variables, campaigns, outfile, outdir, systematics=None):
    """Run the full ABCD estimate for each variable, nominal plus any requested systematic variations."""

    for Var, binning in variables.items():
        run_abcd_variable(Var, binning, campaigns, outfile, outdir)

        for Syst in (systematics or []):
            run_abcd_variable(Var, binning, campaigns, outfile, outdir, Syst=Syst)


if __name__ == "__main__":

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(False)

    campaigns = ["mc23a", "mc23d", "mc23e"]

    # Systematic variations to also run, e.g. ["JET_JER_up", "JET_JER_down"].
    # Each name must have a matching systs/<name>/<Region>/... directory in every input file.
    systematics = []

    # Variable -> analysis binning, shared with anything that reads these histograms back.
    variables = ANALYSIS_BINNING

    outdir = "plots/ABCD"
    os.makedirs(outdir, exist_ok=True)

    outfile_name = f"{outdir}/dijet_ABCD.root"
    outfile = ROOT.TFile(outfile_name, "RECREATE")

    run_abcd(variables, campaigns, outfile, outdir, systematics=systematics)

    outfile.Close()
    print(f"\nWrote {outfile_name}")
