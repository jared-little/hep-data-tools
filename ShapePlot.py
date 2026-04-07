import ROOT
from utilities.ComputeSignificance import get_Zn_histogram, get_SB_histogram
from utilities.GetDataFrame import get_signal_df, get_background_df
from utilities.GetHistograms import get_var_name

ROOT.gROOT.SetStyle("ATLAS")

def make_canvas(name ="canvas"):
    """Make a canvas for plotting."""
    left_margin = 0.15
    right_margin = 0.03

    canvas = ROOT.TCanvas(name, name, 800, 600)
    pad1 = ROOT.TPad("pad1", "pad1", 0., 0.01, .99, 1)
    pad1.SetLeftMargin(left_margin)
    pad1.SetRightMargin(right_margin)
    pad1.SetBottomMargin(0.0)
    pad1.SetTopMargin(0.05)
    pad1.SetFillColor(ROOT.kWhite)
    pad1.SetTickx()
    pad1.SetTicky()
    pad1.SetLogy(1)
    pad1.SetBottomMargin(0.15)

    return canvas, pad1


def set_bins():
    """Set the binning for the histograms, based on the variable being plotted."""
    bins = {
        "NN_score": (10, 0.8, 1),
        "Hbb_bjR_mass": (15, 110, 140),
        # "NN_score": (50, 0, 1),
        # "Hbb_bjR_mass": (50, 0, 300),
        "largeRjetpt": (50, 0, 1000),
        "largeRjetpt_1": (50, 500, 1000), # Trigger turn-on is around 500 GeV, so start there
        "largeRjetpt_2": (50, 0, 1000),
        "largeRjetpt_3": (50, 0, 1000),
        "largeRjetm": (70, 30, 100),
        "largeRjetm_1": (70, 30, 100),
        "largeRjetm_2": (120, 30, 150),
        "largeRjetm_3": (70, 30, 100)
    }
    return bins


def new_columns(df):
    """Define new columns in the RDataFrame for the variables we want to plot."""
    df = df.Define("largeRjetpt_1", "largeRjetpt[0] / 1000")
    df = df.Define("largeRjetpt_2", "largeRjetpt[1] / 1000")
    df = df.Define("largeRjetpt_3", "largeRjetpt[2] / 1000")
    df = df.Define("largeRjetm_1", "largeRjetm[0] / 1000")
    df = df.Define("largeRjetm_2", "largeRjetm[1] / 1000")
    df = df.Define("largeRjetm_3", "largeRjetm[2] / 1000")

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


def plot_shape(Var, optimize, selections=None, plot_all_signals=False):
    """Make plots of the Zn or S/B as a function of the cut value on a given variable, for a given region and rebinning factor."""

    bkg_names = ["dijet", "ttbar", "VV", "Vjets", "top"]
    sig_names = ["XHS_X2000_S1000", "XHS_X3000_S1500", "XHS_X4000_S2000"]
    if plot_all_signals:
        sig_names.extend(["XHS_X1000_S500", "XHS_X2000_S1500", "XHS_X2500_S1500",
                         "XHS_X3000_S2000", "XHS_X3000_S2500", "XHS_X4000_S2000"])

    campaigns = ["mc23a", "mc23d", "mc23e"] # "mc23a, mc23d, mc23e"
    
    colors = [ROOT.kBlue, ROOT.kRed, ROOT.kGreen+2, ROOT.kMagenta, ROOT.kCyan+1]

    hist_bkgs = {}
    for bkg in bkg_names:
        df = get_background_df(bkg, campaigns)
        bkg_hist = make_histogram(df, Var, selections)
        bkg_hist.SetFillColor(colors[bkg_names.index(bkg)])
        bkg_hist.SetLineWidth(1)
        hist_bkgs[bkg] = bkg_hist

    hist_sigs = {}
    for sig in sig_names:
        df = get_signal_df(sig, campaigns)
        sig_hist = make_histogram(df, Var, selections)
        hist_sigs[sig] = sig_hist

    for id, (k, v) in enumerate(reversed(list(hist_bkgs.items()))):
        if id == 0: bkg_histo_total = v.Clone()
        else: bkg_histo_total.Add(v.GetPtr())
    
    index = 0
    for sig_name, sig_hist in hist_sigs.items():
        # automate colors and line styles
        sig_hist.SetLineWidth(4)
        if "X2000_S1000" in sig_name:
            sig_hist.SetLineColor(ROOT.kOrange)
            sig_hist.SetLineStyle(2)
        elif "X3000_S1500" in sig_name:
            sig_hist.SetLineColor(ROOT.kCyan)
            sig_hist.SetLineStyle(2)
        elif "X4000_S2000" in sig_name:
            sig_hist.SetLineColor(ROOT.kViolet)
            sig_hist.SetLineStyle(2)
        else:
            sig_hist.SetLineColor(colors[index % len(colors)])
            sig_hist.SetLineStyle(1)
            sig_hist.SetLineWidth(2)
            index += 1

    canvas, pad1 = make_canvas()

    canvas.cd()
    pad1.Draw()
    pad1.cd()

    leg = ROOT.TLegend(0.55, 0.65, 0.85, 0.9, "")
    leg.SetFillStyle(0)
    leg.SetFillColor(0)
    leg.SetBorderSize(0)

    bkg_histo_total.SetFillColor(0)
    bkg_histo_total.SetLineColor(ROOT.kBlack)
    bkg_histo_total.SetLineWidth(2)
    bkg_histo_total.Draw("HIST")

    bkg_histo_total.GetYaxis().SetTitle("Entries")
    bkg_histo_total.SetMinimum(0.5)
    bkg_histo_total.SetMaximum(10**9)
    bkg_histo_total.GetXaxis().SetTitle(get_var_name(Var))

    for sig_name, sig_hist in hist_sigs.items():
        leg.AddEntry(sig_hist.GetPtr(), "#font[42]{"+sig_name+"}", "l")
        sig_hist.Draw("HIST same")

    leg.AddEntry(bkg_histo_total, "#font[42]{Total Bkg}", "l")
    leg.Draw()

    canvas.cd()

    ROOT.gPad.RedrawAxis()

    canvas.SaveAs(f"plots/Shape/SigBkg-{Var}.pdf")


if __name__ == "__main__":

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(False)
    Optimize = "Zn" # "Zn" or "SB"

    Vars = ["NN_score", "Hbb_bjR_mass"]
    # Vars = ["NN_score",
    #         "largeRjetpt_1", "largeRjetpt_2", "largeRjetpt_3",
    #         "largeRjetm_1", "largeRjetm_2", "largeRjetm_3"]

    # Preselections
    selections = ["largeRjetm_1 > 60", "largeRjetm_2 > 70", "largeRjetm_3 > 70"]
    selections.extend(["largeRjetpt_1 > 500", "largeRjetpt_2 > 350", "largeRjetpt_3 > 200"])
    selections.extend(["Hbb_bjR_mass < 140", "Hbb_bjR_mass > 110"]) # Signal window to study NN_score optimization

    for var in Vars:
        plot_shape(var, Optimize, selections, plot_all_signals=True)
