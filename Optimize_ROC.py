import ROOT
import time
from utilities.ComputeSignificance import get_efficiency_selection
from utilities.GetHistograms import get_signal_histogram, get_bkg_histogram
from array import array

ROOT.gROOT.SetStyle("ATLAS")


def MakeROC(Var="largeRjetpt", VarName="", direction="upper", bg="total"):

    if bg == "total":
      hist_bkg = get_bkg_histogram("ttbar",Var, Region="Preselection", Rebin=1)
      hist_bkg.Add(get_bkg_histogram("dijet",Var, Region="Preselection", Rebin=1))
    else:
      hist_bkg = get_bkg_histogram(bg,Var, Region="Preselection", Rebin=1)

    hist_XHS_2000_1000 = get_signal_histogram("XHS_X2000_S1000",Var, Region="Preselection", Rebin=1)
    hist_XHS_3000_1500 = get_signal_histogram("XHS_X3000_S1500",Var, Region="Preselection", Rebin=1)
    hist_XHS_4000_2000 = get_signal_histogram("XHS_X4000_S2000",Var, Region="Preselection", Rebin=1)

    hist_bkg.Scale(1/hist_bkg.Integral())
    hist_XHS_2000_1000.Scale(1/hist_XHS_2000_1000.Integral())
    hist_XHS_3000_1500.Scale(1/hist_XHS_3000_1500.Integral())
    hist_XHS_4000_2000.Scale(1/hist_XHS_4000_2000.Integral())

    # hist_bkg.GetYaxis().SetTitleOffset(1.3)
    hist_bkg.GetYaxis().SetTitleSize(0.06)

    canName = "Canvas_"+Var
    c = ROOT.TCanvas(canName,canName, 700, 700)
    c.cd()

    hSelection_bkg = get_efficiency_selection(hist_bkg, direction)

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
      value = hist_bkg.Integral()-hSelection_bkg.GetBinContent(sel)
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
