import ROOT
import sys,os
import math
import array
from optparse import OptionParser
from utilities.GetHistograms import get_bkg_histogram, get_data_histogram, get_var_name
from utilities.DijetEstimate import (
    ABCD_FILE,
    ANALYSIS_BINNING,
    campaigns_to_data_years,
    get_dijet_from_abcd_file,
    has_dijet_estimate,
    rebin_histogram,
)


def make_canvas_pads(can_name="DataMC"):
    """Make a canvas with two pads, one for the main plot and one for the ratio plot."""

    c = ROOT.TCanvas(can_name,can_name, 700, 600)
    c.cd()
    pad1 = ROOT.TPad(can_name+"_pad1", can_name+"_pad1", 0., 0.305, .99, 1)
    pad1.SetLogy(1)
    pad1.SetLeftMargin(0.12)
    pad1.SetRightMargin(0.03)
    pad1.SetBottomMargin(0.0)
    pad1.SetFillColor(ROOT.kWhite)
    pad1.SetTickx()
    pad1.SetTicky()
    pad1.Draw()

    pad2 = ROOT.TPad(can_name+"_pad2", can_name+"_pad2", 0., 0.01, .99, 0.295)
    pad2.SetLeftMargin(0.12)
    pad2.SetRightMargin(0.03)
    pad2.SetBottomMargin(0.38)
    pad2.SetFillColor(ROOT.kWhite)
    pad2.SetTickx()
    pad2.SetTicky()
    pad2.Draw()

    return c, pad1, pad2


def make_ratio_plot(hist1, hist2, Var, Region, y_title="Data / MC", y_range=(0.5, 1.5)):
    """Make a ratio plot of hist1 / hist2, with appropriate axis labels and styling."""

    h_ratio = hist1.Clone(f"h_ratio_{Var}_{Region}")
    h_ratio.SetDirectory(0)
    h_ratio.Divide(hist2)

    h_ratio.GetXaxis().SetLabelSize(0.13)
    h_ratio.GetXaxis().SetLabelOffset(0.02)
    h_ratio.GetXaxis().SetTitleSize(0.15)
    h_ratio.GetYaxis().SetRangeUser(*y_range)
    h_ratio.GetYaxis().SetNdivisions(505)
    h_ratio.GetYaxis().SetTitle(y_title)
    h_ratio.GetYaxis().SetLabelSize(0.13)
    h_ratio.GetYaxis().SetTitleSize(0.17)
    h_ratio.GetYaxis().SetTitleOffset(0.36)


    return h_ratio


def _yield_and_error(hist):
    """Return (integral, stat error) over the full histogram range."""
    err = array.array("d", [0.])
    integral = hist.IntegralAndError(1, hist.GetNbinsX(), err)
    return integral, err[0]


def print_yield_table(hist_data, hist_bkgs, bkg_histo, Region, has_data):
    """Print a simple per-process yield table for a region, with each background's % of total background."""

    total_bkg, total_bkg_err = _yield_and_error(bkg_histo)

    print(f"\nYields for {Region}")
    print("-" * 58)
    print(f"{'Process':<12}{'Yield':>12}{'+/- Err':>12}{'% of Bkg':>14}")
    print("-" * 58)
    for name, hist in hist_bkgs.items():
        y, err = _yield_and_error(hist)
        pct = 100 * y / total_bkg if total_bkg else 0
        print(f"{name:<12}{y:>12.1f}{err:>12.1f}{pct:>13.1f}%")
    print("-" * 58)
    print(f"{'Total Bkg':<12}{total_bkg:>12.1f}{total_bkg_err:>12.1f}{100.0:>13.1f}%")

    if has_data:
        data_y, data_err = _yield_and_error(hist_data)
        ratio = data_y / total_bkg if total_bkg else float("nan")
        print(f"{'Data':<12}{data_y:>12.1f}{data_err:>12.1f}{'':>14}")
        print(f"{'Data/Bkg':<12}{ratio:>12.2f}")
    print("-" * 58)


def plot_data_mc(Var, Region, rebin=1, campaigns=["mc23a"], print_yields=False, dijet_method="mc",
                 blind_data=True):
    """Make a data/MC comparison plot for a given variable, region, rebinning factor, and campaigns.

    dijet_method picks where the dijet background comes from. Run ABCD.py first for the
    data-driven options.

        mc              the dijet MC sample
        abcd            the data-driven estimate ABCD.py wrote, so the plot shows the same
                        histogram the fit is given
        abcd-closure    the transfer-factor prediction wherever ABCD.py has one, falling back
                        to the measurement elsewhere

    Whatever the ABCD does not cover keeps the dijet MC: Preselection is not an ABCD region,
    and NN_score and Hbb_bjR_mass are the variables the regions are defined in.

    Under "abcd" the dijet is (data - non-dijet MC), so the stack reproduces data by construction
    and the ratio panel is exactly 1 in every region. "abcd-closure" is what makes the ratio
    mean something: VR1 becomes VR2 x TF tested against VR1 data, and CR1 becomes CR2 x TF,
    which is closed-loop and should land on 1 -- a check that the TF machinery is sound."""

    bkg_names = ["dijet", "ttbar","Vjets", "VV", "top"]
    years = campaigns_to_data_years(campaigns)

    # A prediction is only comparable to a measurement in the binning the transfer factor was
    # measured in. At a finer binning one TF value gets spread flat over several bins, and the
    # CR1 closed loop stops being exact, so put data and MC in the analysis binning as well.
    analysis_binning = ANALYSIS_BINNING.get(Var) if dijet_method == "abcd-closure" else None
    fetch_rebin = 1 if analysis_binning else rebin

    hist_data = get_data_histogram(Var, Region, fetch_rebin, campaigns=years, blind_data=blind_data)
    hist_bkgs = {name: get_bkg_histogram(name, Var, Region, fetch_rebin, campaigns) for name in bkg_names}
    if analysis_binning:
        hist_data = rebin_histogram(hist_data, analysis_binning, f"data_{Region}_{Var}_analysis")
        hist_bkgs = {name: rebin_histogram(hist, analysis_binning, f"{name}_{Region}_{Var}_analysis")
                     for name, hist in hist_bkgs.items()}

    # Keyed off the method rather than off what was found, so that a plot which fell back to the
    # dijet MC still lands in its own file instead of overwriting the --dijet mc output.
    labels = {name: name for name in bkg_names}
    suffix = {"mc": "", "abcd": "_ABCD", "abcd-closure": "_ABCD_closure"}[dijet_method]

    if dijet_method in ("abcd", "abcd-closure"):
        prediction = (dijet_method == "abcd-closure"
                      and has_dijet_estimate(Var, Region, prediction=True))
        if prediction or has_dijet_estimate(Var, Region):
            hist_bkgs["dijet"] = get_dijet_from_abcd_file(Var, Region, rebin, prediction=prediction,
                                                          rebinned=bool(analysis_binning))
            labels["dijet"] = "dijet (CR x TF)" if prediction else "dijet (ABCD)"
        else:
            print(f"  [dijet] no ABCD estimate for {Var} in {Region} in {ABCD_FILE}, keeping the dijet MC")

    colors = [ROOT.kBlue, ROOT.kRed, ROOT.kGreen+2, ROOT.kMagenta, ROOT.kCyan+1]

    leg = ROOT.TLegend(0.55, 0.65, 0.85, 0.9, "")

    stack = ROOT.THStack()
    for id, (k, v) in enumerate(reversed(list(hist_bkgs.items()))):
        # For calculating zn or s/b
        if id == 0: bkg_histo = v.Clone()
        else: bkg_histo.Add(v)
        v.SetFillColor(colors[id])
        v.SetLineWidth(1)
        v.SetName("h"+k)
        leg.AddEntry(v,"#font[42]{"+labels[k]+"}","f")
        # print(f"{k} integral: {v.Integral()}")
        stack.Add(v)

    # bkg_histo.GetYaxis().SetTitle("Events")
    c, pad1, pad2 = make_canvas_pads(f"DataMC_{Region}_{Var}{suffix}")
    # c.cd()
    pad1.cd()

    leg.SetFillStyle(0)
    leg.SetFillColor(0)
    leg.SetBorderSize(0)
    leg.AddEntry(hist_data, "#font[42]{Data}", "p")

    stack.SetMinimum(0.5)
    stack.SetMaximum(10 ** 5)
    # bkg_histo.Draw("E2 same")
    stack.Draw("HIST")
    if "SR" not in Region: hist_data.Draw("eX0 same")
    leg.Draw()

    pad2.cd()
    y_title = "Data / Pred." if suffix else "Data / MC"
    # A closure test scatters well outside the usual band, and ROOT just drops the points that
    # fall off the axis, which would read as missing bins rather than as bad closure.
    y_range = (0.0, 3.0) if dijet_method == "abcd-closure" else (0.5, 1.5)
    if "SR" not in Region: h_ratio = make_ratio_plot(hist_data, bkg_histo, Var, Region, y_title, y_range)
    else: h_ratio = make_ratio_plot(bkg_histo, bkg_histo, Var, Region, y_title, y_range)
    h_ratio.GetXaxis().SetTitle(get_var_name(Var))
    h_ratio.Draw()
    line = ROOT.TLine(hist_data.GetXaxis().GetXmin(), 1, hist_data.GetXaxis().GetXmax(), 1)
    line.Draw("same")
    ROOT.gPad.RedrawAxis()

    c.SaveAs(f"plots/DataMC/DataMC_{Region}_{Var}{suffix}.pdf")

    if print_yields:
        print_yield_table(hist_data, hist_bkgs, bkg_histo, Region, has_data="SR" not in Region)


if __name__ == "__main__":

    parser = OptionParser()
    parser.add_option("--dijet", dest="dijet_method", type="choice", default="mc",
                      choices=["mc", "abcd", "abcd-closure"],
                      help="dijet background: 'mc' for the dijet MC sample, 'abcd' for the "
                           "data-driven estimate, 'abcd-closure' for the transfer-factor "
                           "prediction where there is one [default: %default]")
    (options, args) = parser.parse_args()

    ROOT.gROOT.SetStyle("ATLAS")
    ROOT.gROOT.SetBatch(True)
    os.makedirs("plots/DataMC", exist_ok=True)

    Variable = ["NN_score", "largeRjetpt", "largeRjetm", "NLargeRjets", "mS", "mX"]
    # Variable = ["NN_score", "Hbb_bjR_mass"]

    # Regions = ["Preselection"]
    Regions = ["Preselection", "CR0", "CR1", "CR2", "VR1", "VR2"]
    campaigns = ["mc23a", "mc23d", "mc23e"]

    # Per-variable rebinning; anything not listed uses the default.
    rebin = 4
    rebins = {"NLargeRjets": 1, "NN_score": 5, "mS": 10, "mX": 20}

    for Var in Variable:
        for Region in Regions:
            plot_data_mc(Var, Region, rebin=rebins.get(Var, rebin), campaigns=campaigns,
                         print_yields=(Var == Variable[0]), dijet_method=options.dijet_method)
