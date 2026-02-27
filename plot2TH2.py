import ROOT
import sys,os
import math
from optparse import OptionParser
import numpy as np

# ROOT.gROOT.LoadMacro("atlasstyle-00-04-02/AtlasStyle.C")
# ROOT.gROOT.SetStyle("ATLAS")


# def computeSignificance(s,b,sigma):
#   n = s + b
#   sigma = b * sigma
#   expression = (2 * (n * math.log (n * (b + sigma * sigma)/(b * b + n * sigma * sigma)) - b * b/(sigma * sigma) * math.log((b * b + n * sigma * sigma)/(b * (b + sigma * sigma))) ))

#   if expression < 0:
#     z = 0
#   elif s > 0 and b > 0:
#     z = math.sqrt(expression)
#   else:
#     z = 0
#   return z

def computeSignificance(n, b, sigma):
    if sigma <= 0:
        # No uncertainty case
        if n == 0 or b == 0:
            return 0.0
        t0 = 2 * (n * np.log(n / b) - (n - b))
    else:
        term1 = n * np.log((n * (b + sigma**2)) / (b**2 + n * sigma**2)) if n > 0 else 0
        term2 = (b**2 / sigma**2) * np.log(1 + (sigma**2 * (n - b)) / (b * (b + sigma**2)))
        t0 = 2 * (term1 - term2)

    return np.sqrt(t0) if n >= b else -np.sqrt(t0)

def getSignalHistogram(Signal="XHS_X4000_S2000", Var="NNscoreVSmH", Region="Preselection", Rebin=1):

  file = ROOT.TFile.Open(inputFolder + Signal + "_bbWW_allhad.root")
  sigHistogram = file.Get(f"{Region}/bbVVSplitHadAnalysis_{Signal}_bbWW_allhad_{Region}_{Var}")
  sigHistogram.RebinX(Rebin)
  sigHistogram.RebinY(Rebin)
  sigHistogram.SetDirectory(0)

  return sigHistogram


def getBkgHistogram(Bkg="dijet", Var="NNscoreVSmH", Region="Preselection", Rebin=1):

  file = ROOT.TFile.Open(inputFolder + Bkg + ".root")
  bkgHistogram = file.Get(f"{Region}/bbVVSplitHadAnalysis_{Bkg}_{Region}_{Var}")
  bkgHistogram.RebinX(Rebin)
  bkgHistogram.RebinY(Rebin)
  bkgHistogram.SetDirectory(0)

  return bkgHistogram


def Make2Dplot(signal, Var, Region, optimize, Rebin=1):

  hist_ttbar = getBkgHistogram("ttbar",Var,Region,Rebin)
  hist_dijet = getBkgHistogram("dijet",Var,Region,Rebin)
  sigHist    = getSignalHistogram(signal,Var,Region,Rebin)

  bkgHisto = hist_dijet.Clone("bkgHisto")
  bkgHisto.Add(hist_ttbar)

  canName = "c" + Var
  c = ROOT.TCanvas(canName,canName, 700, 700)
  c.cd()

  hist_SB = sigHist.Clone("hist_SB")
 
  hist_SB.Divide(bkgHisto)
  hist_SB.GetXaxis().SetTitle("m_{H_{bb}} [GeV]")
  hist_SB.GetYaxis().SetTitle("NN Score")
  hist_SB.Draw("COLZ")

  c.cd()
  ROOT.gPad.RedrawAxis()

  c.SaveAs(f"plots/NNscoreVSmH.pdf")

if __name__ == "__main__":

  ROOT.gROOT.SetBatch(True)
  ROOT.gStyle.SetOptStat(False)
  Rebin=2
  # inputFolder = "/data/jlittle/HHARDout/Out_SplitHad/Hists/"
  inputFolder = "/home/jlittle/runHHARD/run/Out_SplitHad/Hists/"

  Make2Dplot("XHS_X2000_S1000", "NNscoreVSmH", "Preselection", optimize="SB", Rebin=Rebin)

