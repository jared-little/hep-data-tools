import ROOT
import sys,os
import math
from optparse import OptionParser

ROOT.gROOT.SetStyle("ATLAS")

def getSignalHistogram(Signal="XHS_X4000_S2000", Var="nn_score", region="Preselection", rebin=1):

  file = ROOT.TFile.Open(inputFolder + Signal + "_bbWW_allhad.root")
  sigHistogram = file.Get(region+"/bbVVSplitHadAnalysis_13p6TeV_"+Signal+"_bbWW_allhad/"+Var)
  sigHistogram.Rebin(rebin)

  if "2000" in Signal: sigHistogram.SetLineColor(ROOT.kOrange)
  if "4000" in Signal: sigHistogram.SetLineColor(ROOT.kCyan)
  if "6000" in Signal: sigHistogram.SetLineColor(ROOT.kViolet)

  sigHistogram.SetLineWidth(4)
  sigHistogram.SetLineStyle(2)
  sigHistogram.SetDirectory(0)

  return sigHistogram


def getBkgHistogram(Bkg="dijet", Var="nn_score", region="Preselection", rebin=1):

  file = ROOT.TFile.Open(inputFolder + Bkg + ".root")
  bkgHistogram = file.Get(region+"/bbVVSplitHadAnalysis_13p6TeV_"+Bkg+"/"+Var)
  bkgHistogram.Rebin(rebin)
  bkgHistogram.SetDirectory(0)

  return bkgHistogram


def getDataHistogram(Var="nn_score", region="Preselection", rebin=1):

  file = ROOT.TFile.Open(inputFolder + "data.root")
  dataHistogram = file.Get(region+"/bbVVSplitHadAnalysis_13p6TeV_data/"+Var)
  dataHistogram.Rebin(rebin)
  dataHistogram.SetDirectory(0)

  return dataHistogram


def getVarName(Var):

  if Var == "leadinglargeRjetpt": xTitle = "Large-R UFO Jet p_{T}"
  elif Var == "leadinglargeRjeteta": xTitle = "Large-R UFO Jet #eta"
  elif Var == "leadinglargeRjetphi": xTitle = "Large-R UFO Jet #phi"
  elif Var == "nn_score": xTitle = "NN Score"
  else: xTitle = ""

  return xTitle


def PlotDataMC(Var, Region, rebin=1):

  hist_ttbar = getBkgHistogram("ttbar",Var, Region, rebin)
  hist_dijet = getBkgHistogram("dijet",Var, Region, rebin)
  hist_top = getBkgHistogram("top",Var, Region, rebin)
  hist_Vjets = getBkgHistogram("Vjets",Var, Region, rebin)
  hist_VV = getBkgHistogram("VV",Var, Region, rebin)
  hist_data = getDataHistogram(Var, Region, rebin)
  hist_XHS_X4000_S2000 = getSignalHistogram("XHS_X4000_S2000",Var, Region, rebin)

  bkgHisto = hist_dijet.Clone("bkgHisto")
  bkgHisto.Add(hist_ttbar)
  bkgHisto.Add(hist_top)
  bkgHisto.Add(hist_Vjets)
  bkgHisto.Add(hist_VV)

  histoList = [hist_dijet, hist_ttbar, hist_top, hist_Vjets, hist_VV]

  Stack = ROOT.THStack()
  for histo in reversed(histoList):
    Stack.Add(histo)
  hist_dijet.SetFillColor(ROOT.kBlue)
  hist_dijet.SetLineWidth(1)
  hist_ttbar.SetFillColor(ROOT.kRed)
  hist_ttbar.SetLineWidth(1)
  hist_top.SetFillColor(ROOT.kGreen)
  hist_top.SetLineWidth(1)
  hist_Vjets.SetFillColor(ROOT.kMagenta)
  hist_Vjets.SetLineWidth(1)
  hist_VV.SetFillColor(ROOT.kYellow)
  hist_VV.SetLineWidth(1)
  bkgHisto.GetXaxis().SetLabelOffset(0.2)
  bkgHisto.GetYaxis().SetTitle("Entries")
  bkgHisto.GetYaxis().SetTitleOffset(1)
  bkgHisto.GetYaxis().SetLabelSize(0.055)
  bkgHisto.GetYaxis().SetTitleSize(0.06)
  hist_data.SetMarkerStyle(20)
  hist_data.SetMarkerColor(1)

  canName = Var+"_Upper"
  c = ROOT.TCanvas(canName,canName, 700, 600)
  c.cd()
  pad1 = ROOT.TPad(canName+"_pad1", canName+"_pad1", 0., 0.305, .99, 1)
  # pad1.SetLogy(1)
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
  Stack.Draw("HIST same")
  hist_data.Draw("eX0 same")

  # bkgHisto.SetMinimum(0.01)
  # bkgHisto.SetMaximum(10**6)
  # for sig in sigHistoList:
  #   sig.Scale(bkgHisto.Integral()/sig.Integral())
  #   sig.Draw("HIST same")
  #   legendEntry = sig.GetName().split("_")[2]+", "+sig.GetName().split("_")[3]
  #   leg.AddEntry(sig, "#font[42]{"+legendEntry+"}", "l")
  leg.AddEntry(hist_data, "#font[42]{Data}", "p")
  for histo in histoList:
    if "ttbar" in histo.GetName(): leg.AddEntry(histo,"#font[42]{t#bar{t}}","f")
    elif "dijet" in histo.GetName(): leg.AddEntry(histo,"#font[42]{Dijet}","f")
    elif "top" in histo.GetName(): leg.AddEntry(histo,"#font[42]{Top}","f")
    elif "Vjets" in histo.GetName(): leg.AddEntry(histo,"#font[42]{V+Jets}","f")
    elif "VV" in histo.GetName(): leg.AddEntry(histo,"#font[42]{VV}","f")
    else: leg.AddEntry(histo,"#font[42]{"+histo.GetName()[1:]+"}","f")

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

  hist_data.GetXaxis().SetTitle(Var)
  hist_data.GetXaxis().SetLabelSize(0.13)
  hist_data.GetXaxis().SetLabelOffset(0.02)
  hist_data.GetXaxis().SetTitleSize(0.15)

  hist_data.GetYaxis().SetRangeUser(0.5,1.5)
  hist_data.GetYaxis().SetNdivisions(505)
  hist_data.GetYaxis().SetTitle("Data / MC")
  hist_data.GetYaxis().SetLabelSize(0.13)
  hist_data.GetYaxis().SetTitleSize(0.17)
  hist_data.GetYaxis().SetTitleOffset(0.36)

  hRatio = hist_data.Clone("hRatio")
  hRatio.Divide(bkgHisto)
  hRatio.Draw()
  ROOT.gPad.RedrawAxis()

  c.SaveAs("plots/DataMC_"+Region+".pdf")

if __name__ == "__main__":

  ROOT.gROOT.SetBatch(True)
  ROOT.gStyle.SetOptStat(False)
  inputFolder = "/data/jlittle/HHARDout/Out_SplitHad/Hists/"

  PlotDataMC("Hbb_bjR_mass_gev", "Preselection", rebin=1)
  PlotDataMC("Hbb_bjR_mass_gev", "CR0", rebin=1)
  PlotDataMC("Hbb_bjR_mass_gev", "CR1", rebin=1)
  PlotDataMC("Hbb_bjR_mass_gev", "CR2", rebin=1)
