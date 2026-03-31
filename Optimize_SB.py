import ROOT
from utilities.ComputeSignificance import get_Zn_histogram, get_SB_histogram
from utilities.GetDataFrame import get_signal_df, get_background_df

ROOT.gROOT.SetStyle("ATLAS")


def make_canvas(name ="canvas", left_margin=0.12, right_margin=0.03):
    """Make a canvas for plotting."""

    canvas = ROOT.TCanvas(name, name, 800, 600)
    pad1 = ROOT.TPad("pad1", "pad1", 0., 0.305, .99, 1)
    pad1.SetLeftMargin(left_margin)
    pad1.SetRightMargin(right_margin)
    pad1.SetBottomMargin(0.025)
    pad1.SetTopMargin(0.05)
    pad1.SetFillColor(ROOT.kWhite)
    pad1.SetTickx()
    pad1.SetTicky()
    pad1.SetLogy(1)

    pad2 = ROOT.TPad("pad2", "pad2", 0., 0.01, .99, 0.295)
    pad2.SetTopMargin(0.05)
    pad2.SetLeftMargin(left_margin)
    pad2.SetRightMargin(right_margin)
    pad2.SetBottomMargin(0.38)
    pad2.SetFillColor(ROOT.kWhite)
    pad2.SetTickx()
    pad2.SetTicky()

    return canvas, pad1, pad2


def set_bins():
    """Set the binning for the histograms, based on the variable being plotted."""
    bins = {
        "NN_score": (50, 0, 1),
        "Hbb_bjR_mass": (50, 0, 300),
        "largeRjetpt_1": (50, 500, 1000), # Trigger turn-on is around 500 GeV, so start there
        "largeRjetpt_2": (50, 0, 1000),
        "largeRjetpt_3": (50, 0, 1000),
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


def make_Zn_plots(Var, optimize, selections=None):
    """Make plots of the Zn or S/B as a function of the cut value on a given variable, for a given region and rebinning factor."""

    bkg_names = ["dijet", "ttbar", "VV", "Vjets", "top"]
    sig_names = ["XHS_X2000_S1000", "XHS_X3000_S1500", "XHS_X4000_S2000"]
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

    stack = ROOT.THStack()
    for id, (k, v) in enumerate(reversed(list(hist_bkgs.items()))):
        if id == 0: bkgHisto = v.Clone()
        else: bkgHisto.Add(v.GetPtr())
        v.SetFillColor(colors[id])
        v.SetLineWidth(1)
        v.SetName("h"+k)
        stack.Add(v.GetPtr())
    
    for sig_name, sig_hist in hist_sigs.items():
        sig_hist.SetLineWidth(4)
        sig_hist.SetLineStyle(2)
        if "X2000" in sig_name: sig_hist.SetLineColor(ROOT.kOrange)
        if "X3000" in sig_name: sig_hist.SetLineColor(ROOT.kCyan)
        if "X4000" in sig_name: sig_hist.SetLineColor(ROOT.kViolet)

    canvas, pad1, pad2 = make_canvas()


    canvas.cd()
    pad1.Draw()
    pad1.cd()

    leg = ROOT.TLegend(0.55, 0.65, 0.85, 0.9, "")
    leg.SetFillStyle(0)
    leg.SetFillColor(0)
    leg.SetBorderSize(0)

    bkgHisto.Draw("E2 same")
    stack.Draw("HIST")
    stack.GetYaxis().SetTitle("Entries")
    stack.SetMinimum(0.01)
    stack.SetMaximum(10**6)

    for sig_name, sig_hist in hist_sigs.items():
        leg.AddEntry(sig_hist.GetPtr(), "#font[42]{"+sig_name+"}", "l")
        sig_hist.Draw("HIST same")

    for k, v in hist_bkgs.items():
       leg.AddEntry(v.GetPtr(),"#font[42]{"+k+"}","f")

    leg.Draw()

    canvas.cd()
    pad2.Draw()
    pad2.cd()

    if optimize == "Zn": hZnUpper, ymax = get_Zn_histogram(hist_sigs, bkgHisto,"upper")
    else: hZnUpper, ymax = get_SB_histogram(hist_sigs, bkgHisto, "upper")

    for h in hZnUpper:
        h.GetXaxis().SetTitle(Var)
        h.GetYaxis().SetRangeUser(0.,ymax+0.5)

    hZnUpper[0].Draw()
    for k in range(1,len(hZnUpper)): hZnUpper[k].Draw("same")
    ROOT.gPad.RedrawAxis()

    if optimize == "Zn": canvas.SaveAs(f"plots/Optimize/ZnOptimizer-{Var}.pdf")
    else: canvas.SaveAs(f"plots/Optimize/SBOptimizer-{Var}.pdf")


if __name__ == "__main__":

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(False)
    Optimize = "SB" # "Zn" or "SB"

    Vars = ["NN_score", "Hbb_bjR_mass"]
    # Vars = ["NN_score",
    #         "largeRjetpt_1", "largeRjetpt_2", "largeRjetpt_3",
    #         "largeRjetm_1", "largeRjetm_2", "largeRjetm_3"]

    # Preselections
    selections = ["largeRjetm_1 > 60", "largeRjetm_2 > 70", "largeRjetm_3 > 70"]
    selections.extend(["largeRjetpt_1 > 500", "largeRjetpt_2 > 350", "largeRjetpt_3 > 200"])

    for var in Vars:
        make_Zn_plots(var, Optimize, selections)
