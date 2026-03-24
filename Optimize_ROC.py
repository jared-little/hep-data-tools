import ROOT
import time
from utilities.ComputeSignificance import get_efficiency_selection
# from utilities.GetHistograms import get_signal_histogram, get_bkg_histogram
from utilities.GetDataFrame import get_signal_df, get_background_df
from array import array

ROOT.gROOT.SetStyle("ATLAS")


def set_bins():
    """Set the binning for the histograms, based on the variable being plotted."""
    bins = {
        "NN_score": (50, 0, 1),
        "largeRjetpt_1": (50, 500, 1000), # Trigger turn-on is around 500 GeV, so start there
        "largeRjetpt_2": (50, 0, 1000),
        "largeRjetpt_3": (50, 0, 1000),
        "largeRjetm_1": (50, 0, 500),
        "largeRjetm_2": (50, 0, 500),
        "largeRjetm_3": (50, 0, 500)
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


def MakeROC(Var="largeRjetpt", VarName="", direction="upper", bg="total"):
    
    bkg_names = ["dijet", "ttbar", "VV", "Vjets", "top"]
    sig_names = ["XHS_X2000_S1000", "XHS_X3000_S1500", "XHS_X4000_S2000"]
    campaigns = ["mc23a", "mc23d", "mc23e"] # "mc23a, mc23d, mc23e"
    selections = ["NN_score > 0.0"] # No selection for now, but can add later for optimization

    hist_bkgs = {}
    hist_bkg_total = None
    for bkg in bkg_names:
        df = get_background_df(bkg, campaigns)
        bkg_hist = make_histogram(df, Var, selections)
        if bkg == bkg_names[0]: hist_bkg_total = bkg_hist.Clone("hist_bkg_total")
        else: hist_bkg_total.Add(bkg_hist.GetPtr())
        hist_bkgs[bkg] = bkg_hist
    hist_bkg_total.Scale(1/hist_bkg_total.Integral())
    
    hist_sigs = {}
    for sig in sig_names:
        df = get_signal_df(sig, campaigns)
        sig_hist = make_histogram(df, Var, selections)
        sig_hist.Scale(1/sig_hist.Integral())
        hist_sigs[sig] = sig_hist

    # hist_bkg_total.GetYaxis().SetTitleOffset(1.3)
    hist_bkg_total.GetYaxis().SetTitleSize(0.06)

    canName = "Canvas_"+Var
    c = ROOT.TCanvas(canName,canName, 700, 700)
    c.cd()

    hSelection_bkg = get_efficiency_selection(hist_bkg_total, direction)
    hist_XHS_2000_1000 = hist_sigs["XHS_X2000_S1000"]
    hist_XHS_3000_1500 = hist_sigs["XHS_X3000_S1500"]
    hist_XHS_4000_2000 = hist_sigs["XHS_X4000_S2000"]

    hSelection_XHS_2000_1000 = get_efficiency_selection(hist_XHS_2000_1000, direction)
    hSelection_XHS_3000_1500 = get_efficiency_selection(hist_XHS_3000_1500, direction)
    hSelection_XHS_4000_2000 = get_efficiency_selection(hist_XHS_4000_2000, direction)

    n = 0
    nSteps = 30
    if "jetm" in Var: nSteps = 60
    sig2000,sig3000,sig4000 = array('d'),array('d'),array('d')
    bkgRejection, ttbarRejection, dijetRejection = array('d'),array('d'),array('d')

    for sel in range(1,nSteps):
      n+=1
      print(
        "Variable:", Var, "\tbin:", hSelection_XHS_2000_1000.GetBinLowEdge(sel), 
        "\tbkg Rejection:", 1-hSelection_bkg.GetBinContent(sel), 
        "\tEff. 2000:", hSelection_XHS_2000_1000.GetBinContent(sel),
        "\tEff. 3000:", hSelection_XHS_3000_1500.GetBinContent(sel), 
        "\tEff. 4000:", hSelection_XHS_4000_2000.GetBinContent(sel))
      value = hSelection_XHS_2000_1000.GetBinContent(sel)
      if hSelection_XHS_2000_1000.GetBinContent(sel)>=1: sig2000.append(1.0)
      else: sig2000.append(value)
      value = hSelection_XHS_3000_1500.GetBinContent(sel)
      if hSelection_XHS_3000_1500.GetBinContent(sel)>=1: sig3000.append(1.0)
      else: sig3000.append(value)
      value = hSelection_XHS_4000_2000.GetBinContent(sel)
      if hSelection_XHS_4000_2000.GetBinContent(sel)>=1: sig4000.append(1.0)
      else: sig4000.append(value)
      value = hist_bkg_total.Integral()-hSelection_bkg.GetBinContent(sel)
      if hSelection_bkg.GetBinContent(sel)<=0: bkgRejection.append(0.)
      else: bkgRejection.append(value)

    gROC_total = ROOT.TGraph(n,sig2000,bkgRejection)
    gROC_total.GetXaxis().SetTitle('Signal Efficiency')
    if bg == "total": gROC_total.GetYaxis().SetTitle('Background Rejection')
    elif bg == "ttbar": gROC_total.GetYaxis().SetTitle('t#bar{t} Rejection')
    else: gROC_total.GetYaxis().SetTitle('jj Rejection')
    gROC_total.Draw("AL")
    ROOT.gPad.RedrawAxis()
    c.SaveAs("plots/ROC/ROC_total_"+Var+"_"+direction+"_"+bg+".pdf")

    return gROC_total


def combine_ROCs(bkg="total"):

    canName = "Canvas"
    c = ROOT.TCanvas(canName,canName, 700, 700)
    c.cd()

    ROC_jet1pt = MakeROC(Var="largeRjetpt_1", VarName="p_{T}^{Jet 1}", direction="upper",bg=bkg)
    ROC_jet2pt = MakeROC(Var="largeRjetpt_2", VarName="p_{T}^{Jet 2}", direction="upper",bg=bkg)
    ROC_jet3pt = MakeROC(Var="largeRjetpt_3", VarName="p_{T}^{Jet 3}", direction="upper",bg=bkg)
    ROC_jet1m = MakeROC(Var="largeRjetm_1", VarName="m_{Jet 1}", direction="upper",bg=bkg)
    ROC_jet2m = MakeROC(Var="largeRjetm_2", VarName="m_{Jet 2}", direction="upper",bg=bkg)
    ROC_jet3m = MakeROC(Var="largeRjetm_3", VarName="m_{Jet 3}", direction="upper",bg=bkg)

    multiGraph = ROOT.TMultiGraph()
    multiGraph.Add(ROC_jet1pt)
    multiGraph.Add(ROC_jet2pt)
    multiGraph.Add(ROC_jet3pt)
    multiGraph.Add(ROC_jet1m)
    multiGraph.Add(ROC_jet2m)
    multiGraph.Add(ROC_jet3m)

    ROC_jet1pt.SetLineWidth(2)
    ROC_jet2pt.SetLineColor(11)
    ROC_jet2pt.SetLineWidth(2)
    ROC_jet3pt.SetLineColor(3)
    ROC_jet3pt.SetLineWidth(2)
    ROC_jet1m.SetLineColor(4)
    ROC_jet1m.SetLineWidth(2)
    ROC_jet2m.SetLineColor(5)
    ROC_jet2m.SetLineWidth(2)
    ROC_jet3m.SetLineColor(6)
    ROC_jet3m.SetLineWidth(2)
    multiGraph.GetXaxis().SetTitle('Signal Efficiency')
    if bkg=="total": multiGraph.GetYaxis().SetTitle('Background Rejection')
    elif bkg=="ttbar": multiGraph.GetYaxis().SetTitle('t#bar{t} Rejection')
    else: multiGraph.GetYaxis().SetTitle('jj Rejection')
    leg = ROOT.TLegend(0.75, 0.65, 0.95, 0.9, "")
    leg.SetFillStyle(0)
    leg.SetFillColor(0)
    leg.SetBorderSize(0)
    leg.SetTextSize(0.025)
    leg.AddEntry(ROC_jet1pt, "p_{T}^{Jet 1}", "L")
    leg.AddEntry(ROC_jet2pt, "p_{T}^{Jet 2}", "L")
    leg.AddEntry(ROC_jet3pt, "p_{T}^{Jet 3}", "L")
    leg.AddEntry(ROC_jet1m, "m^{Jet 1}", "L")
    leg.AddEntry(ROC_jet2m, "m^{Jet 2}", "L")
    leg.AddEntry(ROC_jet3m, "m^{Jet 3}", "L")

    # multiGraph.Draw("apc")
    multiGraph.Draw("AL")
    leg.Draw()

    ROOT.gPad.RedrawAxis()
    c.SaveAs("plots/ROC/ROC_Combined_ " + bkg + ".pdf")

if __name__ == "__main__":


    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(False)

    start_time = time.time()
    combine_ROCs("total")
    # combine_ROCs("ttbar")
    # combine_ROCs("dijet")

    # MakeROC(Var="largeRjetpt_1", VarName="p_{T}^{Jet 1}", direction="upper",bg="total")
    # MakeROC(Var="largeRjetpt_2", VarName="p_{T}^{Jet 2}", direction="upper",bg="total")
    # MakeROC(Var="largeRjetpt_3", VarName="p_{T}^{Jet 1}", direction="upper",bg="total")
    # MakeROC(Var="largeRjetm_1", VarName="p_{T}^{Jet 1}", direction="upper",bg="total")
    # MakeROC(Var="largeRjetm_2", VarName="p_{T}^{Jet 2}", direction="upper",bg="total")
    # MakeROC(Var="largeRjetm_3", VarName="p_{T}^{Jet 1}", direction="upper",bg="total")

    print("--- %s seconds ---" % (time.time() - start_time))
