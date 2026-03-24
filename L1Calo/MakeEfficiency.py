import ROOT

ROOT.gROOT.SetStyle("ATLAS")
# ROOT.gStyle.SetOptStat(0)
# ROOT.gStyle.SetOptTitle(1)


def get_df(fileName="l1calo_hist_noTimingCuts"):

    """Loads the ROOT file and creates a RDataFrame from the gFEX ntuple, defining necessary variables."""

    df = ROOT.RDataFrame("gFEX_ntuple", f"gFEX_ntuples/{fileName}.root")
    df = df.Define("gFEX_largeRjet_et_1", "gFEX_largeRjet_et[0]")
    df = df.Define("reco_leadingLRjet_pt", "AntiKt10UFOCSSKJets_pt[0]")

    return df


def Make_Hist(filename="l1calo_hist_noTimingCuts", rebin_factor=1):

    """Creates histograms for the passed and reconstructed events from the given ROOT file."""

    print(f"Making histograms for {filename}...")
    df = get_df(filename)

    hist_passed = df.Filter("gFEX_largeRjet_et_1 > 50").Histo1D(("hist_passed", "hist_passed", 60, 0, 600), "reco_leadingLRjet_pt").GetPtr()
    hist_reco = df.Histo1D(("hist_reco", "hist_reco", 60, 0, 600), "reco_leadingLRjet_pt").GetPtr()
    hist_reco.GetXaxis().SetTitle("Leading Offline Large-R Jet p_{T} [GeV]")
    hist_reco.SetTitle("gFEX Large-R Jet Trigger, Emulated Threshold: 50 GeV")
    hist_passed.Rebin(rebin_factor)
    hist_reco.Rebin(rebin_factor)

    return hist_passed, hist_reco


def get_legend_name(filename):

    """Generates a legend name based on the filename, replacing specific substrings for better readability."""
    if "data" in filename:
        return "Data"
    elif "noTimingCuts" in filename:
        return "No Timing Cuts"
    elif "withTimingCuts" in filename:
        return "With Timing Cuts"
    elif "TimingWithoutHECwithoutOL" in filename:
        return "Timing Cuts excluding All HEC"
    elif "TimingWithoutHEC" in filename:
        return "Timing Cuts excluding HEC"
    elif "TimingWithoutOL" in filename:
        return "Timing Cuts excluding OL"
    else:
        return filename.split('_')[-1]  # Default to last part of filename if no match


def Make_Efficiency_Plot(filename="l1calo_hist_noTimingCuts", write=False, index=0):

    """Creates an efficiency plot for a given ROOT file containing the gfex ntuple."""

    print(f"Making efficiency plot for {filename}...")
    data = True if "data" in filename else False
    color_array = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen, ROOT.kMagenta, ROOT.kCyan, ROOT.kOrange]

    hist_passed, hist_reco = Make_Hist(filename, rebin_factor=1)

    efficiency = ROOT.TEfficiency(hist_passed, hist_reco)
    if not data:
        efficiency.SetLineColor(color_array[index % len(color_array)])
        efficiency.SetMarkerColor(color_array[index % len(color_array)])
        efficiency.SetMarkerStyle(20)
    else:
        efficiency.SetLineColor(ROOT.kBlack)
        efficiency.SetMarkerColor(ROOT.kBlack)
        efficiency.SetMarkerStyle(21)

    if write:
        file = ROOT.TFile(f"efficiency_histograms.root", "UPDATE")
        efficiency.Write(f"efficiency_{filename.split('_')[-1]}")
        hist_passed.Write(f"hist_passed_{filename.split('_')[-1]}")
        hist_reco.Write(f"hist_reco_{filename.split('_')[-1]}")
        file.Close()
    
    hist_passed.Delete()
    hist_reco.Delete()

    return efficiency


def make_ratio_plot(efficiency1, efficiency2):
    """Creates a ratio plot comparing two efficiency plots."""

    nbin = efficiency1.GetTotalHistogram().GetNbinsX()
    bin_low = efficiency1.GetTotalHistogram().GetBinLowEdge(1)
    bin_high = efficiency1.GetTotalHistogram().GetBinLowEdge(efficiency1.GetTotalHistogram().GetNbinsX() + 1)

    ratio = ROOT.TH1D("ratio", "ratio", nbin, bin_low, bin_high)
    for bin in range(1, efficiency1.GetTotalHistogram().GetNbinsX() + 1):
        eff1 = efficiency1.GetEfficiency(bin)
        eff2 = efficiency2.GetEfficiency(bin)
        ratio.SetBinContent(bin, eff1 / eff2 if eff2 > 0 else 0)
    
    ratio.GetXaxis().SetTitle("Leading Offline Large-R Jet p_{T} [GeV]")
    ratio.SetLineColor(efficiency2.GetLineColor())
    ratio.SetMarkerColor(efficiency2.GetMarkerColor())
    ratio.SetMarkerStyle(efficiency2.GetMarkerStyle())
    ratio.SetLineStyle(efficiency2.GetLineStyle())

    return ratio

def Compare_Timing(plots = []):

    """Compares the efficiency plots for different timing cut configurations.
       Args:     plots (list): List of plot names to compare. 
                 Each name should correspond to a ROOT file containing the gfex ntuple.
    """

    canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)
    canvas.cd()

    legend = ROOT.TLegend(0.40, 0.20, 0.65, 0.28)
    legend.SetBorderSize(0)

    plot_objects = {}  # Keep references to prevent garbage collection
    for name in plots:
        plot = Make_Efficiency_Plot(name, write=False, index=plots.index(name))
        plot_objects[name] = plot  # Store reference
        plot.Draw("AP SAME") if plots.index(name) == 0 else plot.Draw("P SAME")
        legend.AddEntry(plot, get_legend_name(name), "l")

    legend.Draw()
    title = ROOT.TLatex()
    title.SetNDC()
    title.DrawLatex(0.2,  0.96, "Threshold: 50 GeV")
    legend.Draw()

    canvas.Modified()
    canvas.Update()
    canvas.SaveAs("efficiency_comparison.pdf")


def main():
    # plots = ["l1calo_hist_data", "l1calo_hist_noTimingCuts", "l1calo_hist_withTimingCuts",
    #          "l1calo_hist_TimingWithoutOL", "l1calo_hist_TimingWithoutHEC",
    #          "l1calo_hist_TimingWithoutHECwithoutOL"]

    plots = ["l1calo_hist_data", "l1calo_hist_withTimingCuts", "l1calo_hist_noTimingCuts",
             "l1calo_hist_TimingWithoutHECwithoutOL"]
    
    Compare_Timing(plots)


if __name__=="__main__":
    main()
