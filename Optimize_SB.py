import ROOT
# from utilities.ComputeSignificance import computeSignificance
from utilities.ComputeSignificance import GetZnHisto, GetSBHisto
from utilities.GetHistograms import getSignalHistogram, getBkgHistogram

ROOT.gROOT.SetStyle("ATLAS")


def MakeZnPlots(Var, Region, optimize, Rebin=1):
  """Make plots of the Zn or S/B as a function of the cut value on a given variable, for a given region and rebinning factor."""

  # fOutput = ROOT.TFile("ZnOptimizer-XHS.root","UPDATE")

  bkg_names = ["dijet", "ttbar", "VV", "Vjets", "top"]
  campaigns = ["mc23a", "mc23e"] # "mc23a, mc23d, mc23e"

  hist_bkgs = {name: getBkgHistogram(name, Var, Region, Rebin) for name in bkg_names}
  colors = [ROOT.kBlue, ROOT.kRed, ROOT.kGreen+2, ROOT.kMagenta, ROOT.kCyan+1]

  hist_XHS_X2000_S1000_mc23a = getSignalHistogram("XHS_X2000_S1000", Var, Region,Rebin, campaigns)
  hist_XHS_X3000_S1500_mc23a = getSignalHistogram("XHS_X3000_S1500", Var, Region,Rebin, campaigns)
  hist_XHS_X4000_S2000_mc23a = getSignalHistogram("XHS_X4000_S2000", Var, Region,Rebin, campaigns)

  sigHistoDict = {
    "X2000_S1000": hist_XHS_X2000_S1000_mc23a,
    "X3000_S1500": hist_XHS_X3000_S1500_mc23a,
    "X4000_S2000": hist_XHS_X4000_S2000_mc23a
  }

  Stack = ROOT.THStack()
  for id, (k, v) in enumerate(reversed(list(hist_bkgs.items()))):
    # For calculating zn or s/b
    if id == 0: bkgHisto = v.Clone()
    else: bkgHisto.Add(v)
    v.SetFillColor(colors[id])
    v.SetLineWidth(1)
    v.SetName("h"+k)
    Stack.Add(v)

  canName = Var+"_Upper"
  c = ROOT.TCanvas(canName,canName, 700, 600)
  c.cd()
  pad1 = ROOT.TPad(canName+"_pad1", canName+"_pad1", 0., 0.305, .99, 1)
  pad1.SetLogy(1)
  pad1.SetLeftMargin(0.12)
  pad1.SetRightMargin(0.03)
  pad1.SetBottomMargin(0.025)
  pad1.SetTopMargin(0.05)
  pad1.SetFillColor(ROOT.kWhite)
  pad1.SetTickx()
  pad1.SetTicky()
  pad1.Draw()
  pad1.cd()

  leg = ROOT.TLegend(0.55, 0.65, 0.85, 0.9, "")
  leg.SetFillStyle(0)
  leg.SetFillColor(0)
  leg.SetBorderSize(0)

  bkgHisto.Draw("E2 same")
  Stack.Draw("HIST")
  Stack.GetXaxis().SetLabelOffset(0.2)
  Stack.GetYaxis().SetTitle("Entries")
  Stack.GetYaxis().SetTitleOffset(1)
  Stack.GetYaxis().SetLabelSize(0.055)
  Stack.GetYaxis().SetTitleSize(0.06)
  Stack.SetMinimum(0.01)
  Stack.SetMaximum(10**6)

  for sigName, sigHist in sigHistoDict.items():
    sigHist.Draw("HIST same")
    leg.AddEntry(sigHist, "#font[42]{"+sigName+"}", "l")

  for k, v in hist_bkgs.items():
    leg.AddEntry(v,"#font[42]{"+k+"}","f")

  leg.Draw()

  c.cd()
  pad2 = ROOT.TPad(canName+"_pad2", canName+"_pad2", 0., 0.01, .99, 0.295)
  pad2.SetTopMargin(0.05)
  pad2.SetLeftMargin(0.12)
  pad2.SetRightMargin(0.03)
  pad2.SetBottomMargin(0.38)
  pad2.SetFillColor(ROOT.kWhite)
  pad2.SetTickx()
  pad2.SetTicky()
  pad2.Draw()
  pad2.cd()

  if optimize == "Zn": hZnUpper, ymax = GetZnHisto(sigHistoDict,bkgHisto,"upper")
  else: hZnUpper, ymax = GetSBHisto(sigHistoDict,bkgHisto,"upper")

  for h in hZnUpper:
    h.GetXaxis().SetTitle(Var)
    h.GetXaxis().SetLabelSize(0.13)
    h.GetXaxis().SetLabelOffset(0.02)
    h.GetXaxis().SetTitleSize(0.15)

    h.GetYaxis().SetRangeUser(0.,ymax+0.5)
    h.GetYaxis().SetNdivisions(505)
    if optimize == "Zn": h.GetYaxis().SetTitle("Zn")
    else: h.GetYaxis().SetTitle("S / B")
    h.GetYaxis().SetLabelSize(0.13)
    h.GetYaxis().SetTitleSize(0.17)
    h.GetYaxis().SetTitleOffset(0.36)

  hZnUpper[0].Draw()
  for k in range(1,len(hZnUpper)): hZnUpper[k].Draw("same")
  ROOT.gPad.RedrawAxis()

  if optimize == "Zn": c.SaveAs(f"plots/ZnOptimizer-{Region}-{Var}.pdf")
  else: c.SaveAs(f"plots/SBOptimizer-{Region}-{Var}.pdf")


def CalculateAndSaveSignalEfficiency(sig, var):
    '''
    Calculate the signal efficiency for a given signal and variable,
    and save the results to a CSV file.
    '''
    hsig = getSignalHistogram(sig, var)
    
    total_signal_events = hsig.Integral(0, hsig.GetNbinsX() + 1)
    if total_signal_events == 0:
      return []

    efficiencies = []
    num_bins = hsig.GetNbinsX()
    
    for i in range(1, num_bins + 1):
      passing_events = hsig.Integral(i, num_bins + 1)
      efficiency = passing_events / total_signal_events
      cut_value = hsig.GetBinLowEdge(i)
      efficiencies.append((cut_value, efficiency))

    output = f"SignalEfficiency_{sig}.csv"
    with open(output, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['NN_score cut value', 'Sig Efficiency'])
        for cut_value, efficiency in efficiencies:
            writer.writerow([f"{cut_value:.2f}", f"{efficiency:.4f}"])
            
    print(f"Signal efficiencies successfully written to {output}")


if __name__ == "__main__":

  ROOT.gROOT.SetBatch(True)
  ROOT.gStyle.SetOptStat(False)
  Rebin=1

  MakeZnPlots("largeRjetpt_1", "All", optimize="SB", Rebin=Rebin)
  MakeZnPlots("largeRjetpt_2", "All", optimize="SB", Rebin=Rebin)
  MakeZnPlots("largeRjetpt_3", "All", optimize="SB", Rebin=Rebin)

  MakeZnPlots("largeRjetm_1", "All", optimize="SB", Rebin=Rebin)
  MakeZnPlots("largeRjetm_2", "All", optimize="SB", Rebin=Rebin)
  MakeZnPlots("largeRjetm_3", "All", optimize="SB", Rebin=Rebin)

  #CalculateAndSaveSignalEfficiency("XHS_X2000_S1000", "NN_score")
  #CalculateAndSaveSignalEfficiency("XHS_X3000_S1500", "NN_score")
  #CalculateAndSaveSignalEfficiency("XHS_X4000_S2000", "NN_score")
