import ROOT
import sys,os
import math
from optparse import OptionParser
from utilities.GetHistograms import get_bkg_histogram, get_data_histogram, get_var_name


def make_canvas_pads():
    """Make a canvas with two pads, one for the main plot and one for the ratio plot."""

    can_name = "DataMC"
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


def make_ratio_plot(hist1, hist2):
    """Make a ratio plot of hist1 / hist2, with appropriate axis labels and styling."""

    h_ratio = hist1.Clone("h_ratio")
    h_ratio.Divide(hist2)

    h_ratio.GetXaxis().SetTitle(get_var_name(Var))
    h_ratio.GetXaxis().SetLabelSize(0.13)
    h_ratio.GetXaxis().SetLabelOffset(0.02)
    h_ratio.GetXaxis().SetTitleSize(0.15)
    h_ratio.GetYaxis().SetRangeUser(0.5,1.5)
    h_ratio.GetYaxis().SetNdivisions(505)
    h_ratio.GetYaxis().SetTitle("Data / MC")
    h_ratio.GetYaxis().SetLabelSize(0.13)
    h_ratio.GetYaxis().SetTitleSize(0.17)
    h_ratio.GetYaxis().SetTitleOffset(0.36)


    return h_ratio



def _campaigns_to_data_years(campaigns):
    """Convert a list of campaign names to corresponding run-3 data years."""
    years = []
    for campaign in campaigns:
        if "a" in campaign:
            years.append("22")
        elif "d" in campaign:
            years.append("23")
        elif "e" in campaign:
            years.append("24")
    return years

def plot_data_mc(Var, Region, rebin=1, campaigns=["mc23a"]):
    """Make a data/MC comparison plot for a given variable, region, rebinning factor, and campaigns."""

    bkg_names = ["dijet", "ttbar","Vjets", "VV", "top"]
    years = _campaigns_to_data_years(campaigns)

    hist_data = get_data_histogram(Var, Region, rebin, campaigns = years)
    hist_bkgs = {name: get_bkg_histogram(name, Var, Region, rebin, campaigns) for name in bkg_names}
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
        leg.AddEntry(v,"#font[42]{"+k+"}","f")
        # print(f"{k} integral: {v.Integral()}")
        stack.Add(v)

    # bkg_histo.GetYaxis().SetTitle("Events")
    c, pad1, pad2 = make_canvas_pads()
    # c.cd()
    pad1.cd()

    leg.SetFillStyle(0)
    leg.SetFillColor(0)
    leg.SetBorderSize(0)
    leg.AddEntry(hist_data, "#font[42]{Data}", "p")

    bkg_histo.Draw("E2 same")
    stack.Draw("HIST same")
    hist_data.Draw("eX0 same")
    leg.Draw()

    pad2.cd()
    h_ratio = make_ratio_plot(hist_data, bkg_histo)
    h_ratio.Draw()
    line = ROOT.TLine(hist_data.GetXaxis().GetXmin(), 1, hist_data.GetXaxis().GetXmax(), 1)
    line.Draw("same")
    ROOT.gPad.RedrawAxis()

    c.SaveAs(f"plots/DataMC/DataMC_{Region}_{Var}.pdf")


if __name__ == "__main__":

    ROOT.gROOT.SetStyle("ATLAS")
    ROOT.gROOT.SetBatch(True)

    inputFolder = "/data/jlittle/HHARDout/Out_SplitHad/Hists/"
    # Variable = ["NN_score", "largeRjetpt", "largeRjetm", "NLargeRjets"]
    Variable = ["NN_score"]

    Regions = ["Preselection", "Preselection_CR0", "Preselection_VR2", "Preselection_CR2"]
    campaigns = ["mc23a", "mc23d", "mc23e"]
    rebin = 2

    for Var in Variable:
        for Region in Regions:
            plot_data_mc(Var, Region, rebin=1, campaigns=campaigns) if Var == "NLargeRjets" else plot_data_mc(Var, Region, rebin=rebin, campaigns=campaigns)
