import ROOT

ROOT.gROOT.SetStyle("ATLAS")
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(1)


def get_df(fileName="l1calo_hist_noTimingCuts"):

    """Loads the ROOT file and creates a RDataFrame from the gFEX ntuple, defining necessary variables."""

    df = ROOT.RDataFrame("gFEX_ntuple", f"gFEX_ntuples/{fileName}.root")
    df = df.Define("gFEX_largeRjet_et_1", "gFEX_largeRjet_et[0]")
    df = df.Define("reco_leadingLRjet_pt", "AntiKt10UFOCSSKJets_pt[0]")

    return df


def Make_Hist(filename="l1calo_hist_noTimingCuts"):

    """Creates histograms for the passed and reconstructed events from the given ROOT file."""

    print(f"Making histograms for {filename}...")
    df = get_df(filename)

    hist_passed = df.Filter("gFEX_largeRjet_et_1 > 50").Histo1D(("hist_passed", "hist_passed", 60, 0, 600), "reco_leadingLRjet_pt").GetPtr()
    hist_reco = df.Histo1D(("hist_reco", "hist_reco", 60, 0, 600), "reco_leadingLRjet_pt").GetPtr()
    hist_reco.GetXaxis().SetTitle("Leading Offline Large-R Jet p_{T} [GeV]")
    hist_reco.SetTitle("gFEX Large-R Jet Trigger, Emulated Threshold: 50 GeV")

    return hist_passed, hist_reco


def Make_Efficiency_Plot(filename="l1calo_hist_noTimingCuts", write=False, index=0):

    """Creates an efficiency plot for a given ROOT file containing the gfex ntuple."""

    print(f"Making efficiency plot for {filename}...")
    color_array = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen, ROOT.kMagenta, ROOT.kCyan, ROOT.kOrange]

    hist_passed, hist_reco = Make_Hist(filename)

    efficiency = ROOT.TEfficiency(hist_passed, hist_reco)
    efficiency.SetLineColor(color_array[index % len(color_array)])
    efficiency.SetMarkerColor(color_array[index % len(color_array)])
    efficiency.SetMarkerStyle(20)

    if write:
        file = ROOT.TFile(f"efficiency_histograms.root", "UPDATE")
        efficiency.Write(f"efficiency_{filename.split('_')[-1]}")
        hist_passed.Write(f"hist_passed_{filename.split('_')[-1]}")
        hist_reco.Write(f"hist_reco_{filename.split('_')[-1]}")
        file.Close()
    
    hist_passed.Delete()
    hist_reco.Delete()

    return efficiency


def Compare_Timing(plots = []):

    """Compares the efficiency plots for different timing cut configurations.
       Args:     plots (list): List of plot names to compare. 
                 Each name should correspond to a ROOT file containing the gfex ntuple.
    """

    canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)
    legend = ROOT.TLegend(0.40, 0.20, 0.65, 0.28)
    legend.SetBorderSize(0)
    canvas.cd()
    title = ROOT.TLatex()
    title.SetNDC()

    plot_objects = []  # Keep references to prevent garbage collection
    for name in plots:
        plot = Make_Efficiency_Plot(name, write=False, index=plots.index(name))
        plot_objects.append(plot)  # Store reference
        plot.Draw("AP SAME") if plots.index(name) == 0 else plot.Draw("P SAME")
        legend.AddEntry(plot, name.split('_')[-1].replace("noTimingCuts", "No Timing Cuts").replace("withTimingCuts", "With Timing Cuts"), "l")

    legend.Draw()
    title.DrawLatex(0.2,  0.96, "Threshold: 50 GeV");
    canvas.Modified()
    canvas.Update()
    canvas.SaveAs("efficiency_comparison.pdf")


def main():
    plots = ["l1calo_hist_noTimingCuts", "l1calo_hist_withTimingCuts"]
    Compare_Timing(plots)


if __name__=="__main__":
    main()
