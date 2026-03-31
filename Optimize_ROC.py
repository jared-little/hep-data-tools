import ROOT
import time
from utilities.ComputeSignificance import get_efficiency_selection
from utilities.GetDataFrame import get_signal_df, get_background_df
from array import array

ROOT.gROOT.SetStyle("ATLAS")


def set_bins():
    """Set the binning for the histograms, based on the variable being plotted."""
    bins = {
        "NN_score": (50, 0, 1),
        "largeRjetpt_1": (50, 500, 1000), # Trigger turn-on is around 500 GeV, so start there
        "largeRjetpt_2": (80, 200, 1000),
        "largeRjetpt_3": (80, 200, 1000),
        "largeRjetm_1": (60, 40, 100),
        "largeRjetm_2": (110, 40, 150),
        "largeRjetm_3": (60, 40, 100)
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


def MakeROC(Var="largeRjetpt", VarName="", direction="upper", bg="total", selections=None, print_efficiencies=False):
    
    bkg_names = ["dijet", "ttbar", "VV", "Vjets", "top"]
    sig_names = ["XHS_X2000_S1000", "XHS_X3000_S1500", "XHS_X4000_S2000"]
    campaigns = ["mc23a", "mc23d", "mc23e"] # "mc23a, mc23d, mc23e"

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

    h_eff_bkg = get_efficiency_selection(hist_bkg_total, direction)
    h_2000_1000 = hist_sigs["XHS_X2000_S1000"]
    h_3000_1500 = hist_sigs["XHS_X3000_S1500"]
    h_4000_2000 = hist_sigs["XHS_X4000_S2000"]

    h_eff_2000_1000 = get_efficiency_selection(h_2000_1000, direction)
    h_eff_3000_1500 = get_efficiency_selection(h_3000_1500, direction)
    h_eff_4000_2000 = get_efficiency_selection(h_4000_2000, direction)

    n = 0
    print(f"bkg bins: {h_eff_bkg.GetNbinsX()}, sig bins: {h_eff_2000_1000.GetNbinsX()}")
    # nSteps = 30
    nSteps = h_eff_bkg.GetNbinsX()
    # if "jetm" in Var: nSteps = 800
    sig2000,sig3000,sig4000 = array('d'),array('d'),array('d')
    bkgRejection, ttbarRejection, dijetRejection = array('d'),array('d'),array('d')

    for sel in range(1,nSteps):
      n+=1
      if print_efficiencies:
        print(
            "Variable:", Var, "\tbin:", h_eff_2000_1000.GetBinLowEdge(sel), 
            "\tbkg Rejection:", 1-h_eff_bkg.GetBinContent(sel), 
            "\tEff. 2000:", h_eff_2000_1000.GetBinContent(sel),
            "\tEff. 3000:", h_eff_3000_1500.GetBinContent(sel), 
            "\tEff. 4000:", h_eff_4000_2000.GetBinContent(sel))
      value = h_eff_2000_1000.GetBinContent(sel)
      if h_eff_2000_1000.GetBinContent(sel)>=1: sig2000.append(1.0)
      else: sig2000.append(value)
      value = h_eff_3000_1500.GetBinContent(sel)
      if h_eff_3000_1500.GetBinContent(sel)>=1: sig3000.append(1.0)
      else: sig3000.append(value)
      value = h_eff_4000_2000.GetBinContent(sel)
      if h_eff_4000_2000.GetBinContent(sel)>=1: sig4000.append(1.0)
      else: sig4000.append(value)
      value = hist_bkg_total.Integral()-h_eff_bkg.GetBinContent(sel)
      if h_eff_bkg.GetBinContent(sel)<=0: bkgRejection.append(0.)
      else: bkgRejection.append(value)

    gROC_total = ROOT.TGraph(n,sig2000,bkgRejection)
    gROC_total.GetXaxis().SetTitle('Signal Efficiency')
    if bg == "total": gROC_total.GetYaxis().SetTitle('Background Rejection')
    elif bg == "ttbar": gROC_total.GetYaxis().SetTitle('t#bar{t} Rejection')
    else: gROC_total.GetYaxis().SetTitle('jj Rejection')
    gROC_total.Draw("AL")
    ROOT.gPad.RedrawAxis()
    c.SaveAs("plots/ROC/ROC_"+Var+"_"+direction+"_"+bg+".pdf")

    return gROC_total


def combine_ROCs(bkg="total", variables = ["largeRjetpt_1"], selections=None):

    style = {
        "largeRjetpt_1": [ROOT.kBlue, "Jet 1 - P_{T}"],
        "largeRjetpt_2": [ROOT.kRed, "Jet 2 - P_{T}"],
        "largeRjetpt_3": [ROOT.kGreen+2, "Jet 3 - P_{T}"],
        "largeRjetm_1": [ROOT.kMagenta, "Jet 1 - M"],
        "largeRjetm_2": [ROOT.kCyan+1, "Jet 2 - M"],
        "largeRjetm_3": [ROOT.kOrange+1, "Jet 3 - M"]
    }

    canName = "Canvas"
    c = ROOT.TCanvas(canName,canName, 700, 700)
    c.cd()
    roc_dict = {}
    multiGraph = ROOT.TMultiGraph()
    leg = ROOT.TLegend(0.75, 0.65, 0.95, 0.9, "")
    leg.SetFillStyle(0)
    leg.SetFillColor(0)
    leg.SetBorderSize(0)

    for var in variables:
        roc = MakeROC(Var=var, VarName=var, direction="upper", bg=bkg, selections = selections)
        roc_dict[var] = roc
        multiGraph.Add(roc)
        leg.AddEntry(roc, style[var][1], "L")
        ROOT.gPad.RedrawAxis()

    for var, roc in roc_dict.items():
        print(f"Index: {variables.index(var)}, Variable: {var}")
        roc.SetLineWidth(2)
        roc.SetLineColor(style[var][0])
        roc.Draw("L")

    multiGraph.GetXaxis().SetTitle('Signal Efficiency')
    multiGraph.GetYaxis().SetTitle('Background Rejection')
    # multiGraph.GetXaxis().SetLimits(0.8,1.0)
    # multiGraph.GetYaxis().SetRangeUser(0.0, 0.5)
    multiGraph.Draw("AL")
    leg.Draw()

    ROOT.gPad.RedrawAxis()
    c.SaveAs("plots/ROC/ROC_Combined_" + bkg + ".pdf")


if __name__ == "__main__":

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(False)

    start_time = time.time()
    selections = ["NN_score > 0.5", "largeRjetm_1 > 75", "largeRjetm_2 > 75", "largeRjetm_3 > 70"]

    variables = ["largeRjetpt_1", "largeRjetpt_2", "largeRjetpt_3", "largeRjetm_1", "largeRjetm_2", "largeRjetm_3"]
    combine_ROCs("total", variables, selections)
    # combine_ROCs("ttbar")
    # combine_ROCs("dijet")

    # MakeROC(Var="largeRjetpt_1", VarName="p_{T}^{Jet 1}", direction="upper",bg="total", print_efficiencies=True)
    # MakeROC(Var="largeRjetpt_2", VarName="p_{T}^{Jet 2}", direction="upper",bg="total")
    # MakeROC(Var="largeRjetpt_3", VarName="p_{T}^{Jet 1}", direction="upper",bg="total")
    # MakeROC(Var="largeRjetm_1", VarName="p_{T}^{Jet 1}", direction="upper",bg="total")
    # MakeROC(Var="largeRjetm_2", VarName="p_{T}^{Jet 2}", direction="upper",bg="total")
    # MakeROC(Var="largeRjetm_3", VarName="p_{T}^{Jet 1}", direction="upper",bg="total")

    print("--- %s seconds ---" % (time.time() - start_time))
