import ROOT
from utilities.GetHistograms import get_signal_histogram, get_bkg_histogram, get_var_name

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


def print_yields(Region, hist_bkgs, hist_sigs):
    """Print the event yield for each background and signal process in a given region."""

    print(f"\n=== Yields in region: {Region} ===")
    print("-- Background --")
    total_bkg = 0.0
    for name, hist in hist_bkgs.items():
        yield_ = hist.Integral()
        total_bkg += yield_
        print(f"  {name:12s}: {yield_:10.2f}")
    print(f"  {'Total bkg':12s}: {total_bkg:10.2f}")

    print("-- Signal --")
    for name, hist in hist_sigs.items():
        print(f"  {name:22s}: {hist.Integral():10.2f}")


def plot_shape(Var, Region, campaigns, rebin=1, plot_all_signals=False):
    """Make a shape plot of signal vs. stacked background for a given variable and region."""

    bkg_names = ["dijet", "ttbar", "VV", "Vjets", "top"]
    sig_names = ["XHS_X2000_S1000", "XHS_X3000_S1500", "XHS_X4000_S2000"]
    if plot_all_signals:
        sig_names.extend(["XHS_X1000_S500", "XHS_X2000_S1500", "XHS_X2500_S1500",
                         "XHS_X3000_S2000", "XHS_X3000_S2500", "XHS_X4000_S2000"])

    colors = [ROOT.kBlue, ROOT.kRed, ROOT.kGreen+2, ROOT.kMagenta, ROOT.kCyan+1]

    hist_bkgs = {name: get_bkg_histogram(name, Var, Region, rebin, campaigns) for name in bkg_names}
    hist_sigs = {name: get_signal_histogram(name, Var, Region, rebin, campaigns) for name in sig_names}

    print_yields(Region, hist_bkgs, hist_sigs)

    leg = ROOT.TLegend(0.55, 0.65, 0.85, 0.9, "")
    leg.SetFillStyle(0)
    leg.SetFillColor(0)
    leg.SetBorderSize(0)

    stack = ROOT.THStack()
    for id, (k, v) in enumerate(reversed(list(hist_bkgs.items()))):
        v.SetFillColor(colors[id])
        v.SetLineColor(ROOT.kBlack)
        v.SetLineWidth(1)
        leg.AddEntry(v, "#font[42]{"+k+"}", "f")
        stack.Add(v)

    index = 0
    for sig_name, sig_hist in hist_sigs.items():
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
            index += 1

    canvas, pad1 = make_canvas()

    canvas.cd()
    
    pad1.Draw()
    pad1.cd()
    stack.SetMinimum(0.1)

    stack.Draw("HIST")
    stack.GetXaxis().SetTitle(get_var_name(Var))
    stack.GetYaxis().SetTitle("Events")
    if Region == "SR" and Var == "Hbb_bjR_mass":
        stack.GetXaxis().SetRangeUser(110, 140)
        stack.SetMaximum(1e3)
        # pad1.SetLogy(0)

    for sig_name, sig_hist in hist_sigs.items():
        leg.AddEntry(sig_hist, "#font[42]{"+sig_name+"}", "l")
        sig_hist.Draw("HIST same")

    leg.Draw()

    canvas.cd()

    ROOT.gPad.RedrawAxis()

    canvas.SaveAs(f"plots/Shape/SigBkg-{Region}-{Var}_GN3X.pdf")


if __name__ == "__main__":

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(False)

    Vars = ["Hbb_bjR_mass"]
    # Vars = ["NN_score"]
    # Vars = ["NN_score",
    #         "largeRjetpt", "largeRjetpt", "largeRjetpt",
    #         "largeRjetm", "largeRjetm", "largeRjetm"]

    Regions = ["SR"]
    # Regions = ["Preselection", "CR0", "CR1", "CR2", "VR1", "VR2","SR"]

    campaigns = ["mc23a", "mc23d", "mc23e"]

    for region in Regions:
        for var in Vars:
            plot_shape(var, region, campaigns, rebin=1, plot_all_signals=False)
