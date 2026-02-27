import ROOT


# inputFolder = "/data/hslien/bbVV_0lep/hist/ABCDRegion/"
inputFolder = "/data/jlittle/HHARDout/all_mc23_nominal/"


def getSignalHistogram(Signal="XHS_X4000_S2000", Var="NN_score", Region="All", Rebin=1, campaigns=["mc23a"]):
    """Get the signal histogram for a given signal, variable, region, and rebinning factor."""

    hists = []
    for campaign in campaigns:
        file = ROOT.TFile.Open(inputFolder + campaign + "_" + Signal + "_bbWW_allhad.root")
        hist = file.Get(f"{Region}/bbVVSplitHadAnalysis_13p6TeV_{Signal}_bbWW_allhad/{Var}")
        # hist = file.Get(f"{Region}/bbVVSplitHadAnalysis_{Signal}_bbWW_allhad_{Region}_{Var}")
        hists.append(hist)
    
    # Sum histograms from different campaigns
    sigHistogram = hists[0].Clone()
    for hist in hists[1:]:
        sigHistogram.Add(hist)

    if "2000" in Signal: sigHistogram.SetLineColor(ROOT.kOrange)
    if "4000" in Signal: sigHistogram.SetLineColor(ROOT.kCyan)
    if "6000" in Signal: sigHistogram.SetLineColor(ROOT.kViolet)

    sigHistogram.SetLineWidth(4)
    sigHistogram.SetLineStyle(2)
    sigHistogram.SetDirectory(0)
    sigHistogram.Rebin(Rebin)

    return sigHistogram


def getBkgHistogram(Bkg="dijet", Var="NN_score", Region="All", Rebin=1, campaigns=["mc23a"]):
    """Get the background histogram for a given background, variable, region, and rebinning factor."""

    hists = []
    for campaign in campaigns:
        file = ROOT.TFile.Open(inputFolder + campaign + "_" + Bkg + ".root")

        hist = file.Get(f"{Region}/bbVVSplitHadAnalysis_13p6TeV_{Bkg}/{Var}")
        # hist = file.Get(f"{Region}/bbVVSplitHadAnalysis_{Bkg}_{Region}_{Var}")
        hists.append(hist)

    bkgHistogram = hists[0].Clone()
    for hist in hists[1:]:
        bkgHistogram.Add(hist)

    bkgHistogram.Rebin(Rebin)
    bkgHistogram.SetDirectory(0)

    return bkgHistogram


def getVarName(Var):
    """
    Get the x-axis title for a given variable.
    CURRENTLY UNUSED, SHOULD BE MOVED TO A DICTIONARY
    """

    if Var == "leadinglargeRjetpt": xTitle = "Large-R UFO Jet p_{T}"
    elif Var == "leadinglargeRjeteta": xTitle = "Large-R UFO Jet #eta"
    elif Var == "leadinglargeRjetphi": xTitle = "Large-R UFO Jet #phi"
    else: xTitle = ""

    return xTitle
