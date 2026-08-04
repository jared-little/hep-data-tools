"""Data-driven (ABCD) estimate of the dijet background.

In every ABCD region the dijet background is whatever is left of the data once the
non-dijet MC has been subtracted:

    dijet(Region) = data(Region) - sum(MC_BACKGROUNDS)(Region)

That subtraction cannot be used in a region one wants to keep unbiased by data (the SR),
so those regions are predicted instead, with the per-bin transfer factor measured in the
control regions:

    TF  = CR1 / CR2
    SR  = CR0 x TF          (prediction)
    VR1 = VR2 x TF          (closure test)

This module holds the pieces ABCD.py uses to build the estimate, plus the reader other
scripts use to pick it back up out of ABCD.py's output file, so that a plot shows the same
histogram the fit is given rather than a recomputation that can drift.
"""

import math
import os
from array import array

import ROOT

from utilities.GetHistograms import get_bkg_histogram, get_data_histogram, get_detached_histogram

# Regions where the dijet background is measured directly as data - MC.
ABCD_REGIONS = ["CR0", "CR1", "CR2", "VR1", "VR2"]

# The SR is not measured, it is predicted as CR0 x TF. It is still written to the output file.
SR_REGION = "SR"

# Regions that are measured, but that the transfer factor can also predict from another
# region: <region to predict>: <source region>. ABCD.py writes both so a plot can put one
# against the other.
#   VR1 = VR2 x TF   the real closure test, in the inner mass sideband
#   CR1 = CR2 x TF   closed-loop by construction, since TF is itself CR1/CR2. This predicts
#                    nothing physical; it is a check that the TF, the bin mapping and the
#                    file round-trip are all sound, and it should come out at exactly 1.
PREDICTED_REGIONS = {"VR1": "VR2", "CR1": "CR2"}

# Everything that is subtracted from data; what is left over is the dijet background.
MC_BACKGROUNDS = ["ttbar", "top", "Vjets", "VV"]

# Where ABCD.py writes its histograms, and what the reader below defaults to.
ABCD_FILE = "plots/ABCD/dijet_ABCD.root"

# Variable -> analysis binning: either a list of bin edges or an integer group size. This is the
# binning ABCD.py measures the transfer factor in, so anything comparing a prediction against a
# measurement has to use it too. Adding a variable needs an entry in GetHistograms.get_var_name,
# and a binning coarse enough that data - MC stays positive in every region (check_positive_bins
# will say so otherwise).
#
# The two variables that define the ABCD plane are deliberately absent. The regions are cuts on
# Hbb_bjR_mass (window vs. sidebands) and on NN_score (CR1/VR1 above 0.95, the rest below), so
# CR1 and CR2 do not overlap in either one and TF = CR1/CR2 is degenerate -- predicting them
# this way is circular.
ANALYSIS_BINNING = {
    "mX": [1000, 1900, 2100, 2300, 2500, 2700, 2900, 3100, 3300, 3500, 3800, 5000],
    "mS": [0, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2600, 3000],
    "Hbb_bjR_pt": [40, 285, 530, 775, 1020, 1265, 2000],
    "largeRjetpt": [0, 200, 400, 600, 800, 1000, 1200, 2000],
    "largeRjetm": [0, 50, 100, 150, 200, 250, 300, 600],
    "NLargeRjets": [0, 1, 2, 3, 4, 10],
}


def campaigns_to_data_years(campaigns):
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


def is_abcd_region(Region):
    """Whether the data-driven dijet estimate is defined for this region."""

    return Region in ABCD_REGIONS or Region == SR_REGION


def get_dijet_histogram(Var, Region, Rebin=1, campaigns=None, Syst=None, blind_data=False):
    """Get the dijet estimate (data - non-dijet MC) for one variable in one region.
    Data has no systematic variations, so Syst only affects the MC backgrounds that get subtracted.

    Note that blind_data has to match whatever the data histogram in the plot uses: blinding data
    without blinding this estimate (or the other way round) breaks the data = dijet + MC identity
    this method is built on."""

    if campaigns is None:
        campaigns = ["mc23a"]

    years = campaigns_to_data_years(campaigns)
    hist = get_data_histogram(Var, Region, rebin=Rebin, campaigns=years, blind_data=blind_data)
    hist.SetName(f"{Var}_{Region}")
    if not hist.GetSumw2N():
        hist.Sumw2()

    for bkg in MC_BACKGROUNDS:
        hist_bkg = get_bkg_histogram(bkg, Var, Region, Rebin=Rebin, campaigns=campaigns, Syst=Syst)
        if hist_bkg.GetNbinsX() != hist.GetNbinsX():
            raise ValueError(f"Binning mismatch between {bkg} and data for {Var} in {Region}")
        hist.Add(hist_bkg, -1.0)

    return hist


def rebin_histogram(hist, binning, name):
    """Rebin a histogram, with either a group size (int) or a list of bin edges."""

    if isinstance(binning, int):
        return hist.Rebin(binning, name)

    edges = array("d", binning)

    return hist.Rebin(len(edges) - 1, name, edges)


def check_positive_bins(hist, name=None):
    """Raise if any bin is negative. data - MC can undershoot, and a negative bin would
    silently poison the transfer factor and the SR prediction."""

    name = name or hist.GetName()
    negative_bins = [b for b in range(1, hist.GetNbinsX() + 1) if hist.GetBinContent(b) < 0]
    if negative_bins:
        raise ValueError(f"Negative bin content in {name}, bins {negative_bins}. Widen the binning.")


def zero_negative_bins(hist):
    """Clamp negative bins to zero and return how many there were. The bin errors are left
    alone, so a bin that got clamped still carries the uncertainty of the subtraction."""

    n_clamped = 0
    for b in range(0, hist.GetNbinsX() + 2):
        if hist.GetBinContent(b) < 0:
            hist.SetBinContent(b, 0.0)
            n_clamped += 1

    return n_clamped


def get_transfer_factor(hist_cr1, hist_cr2, name="TF"):
    """Get the per-bin transfer factor TF = CR1 / CR2."""

    hist_tf = hist_cr1.Clone(name)
    hist_tf.SetDirectory(0)
    hist_tf.Divide(hist_cr2)

    return hist_tf


def apply_transfer_factor(hist, hist_tf, name):
    """Extrapolate a region with the transfer factor, e.g. SR = CR0 x TF.

    hist_tf may be binned more coarsely than hist, in which case every bin of hist picks up
    the transfer factor of the TF bin containing its centre."""

    hist_extrapolated = hist.Clone(name)
    hist_extrapolated.SetDirectory(0)

    if hist_extrapolated.GetNbinsX() == hist_tf.GetNbinsX():
        hist_extrapolated.Multiply(hist_tf)
        return hist_extrapolated

    for b in range(1, hist_extrapolated.GetNbinsX() + 1):
        tf_bin = hist_tf.GetXaxis().FindBin(hist_extrapolated.GetXaxis().GetBinCenter(b))
        tf, tf_error = hist_tf.GetBinContent(tf_bin), hist_tf.GetBinError(tf_bin)
        content, error = hist_extrapolated.GetBinContent(b), hist_extrapolated.GetBinError(b)
        hist_extrapolated.SetBinContent(b, content * tf)
        # Same error propagation as TH1::Multiply, treating the two as independent.
        hist_extrapolated.SetBinError(b, math.sqrt((error * tf) ** 2 + (content * tf_error) ** 2))

    return hist_extrapolated


def abcd_histogram_path(Var, Region, Syst=None, rebinned=False, prediction=False):
    """Path of one histogram inside ABCD.py's output file, following the HHARD/TRExFitter layout.

    <Var> is always the native input binning and <Var>_rebin the analysis binning, in every
    region including the SR. Regions in PREDICTED_REGIONS also carry <Var>_prediction, the
    transfer factor applied to their source region instead of the direct measurement."""

    name = Var
    if prediction:
        name += "_prediction"
    if rebinned:
        name += "_rebin"

    return f"systs/{Syst}/{Region}/{name}" if Syst else f"{Region}/{name}"


def has_dijet_estimate(Var, Region, Syst=None, rebinned=False, prediction=False, path=ABCD_FILE):
    """Whether ABCD.py has actually written this histogram. Cheaper than catching the read,
    and it lets a caller fall back quietly for a variable the ABCD does not cover."""

    if not is_abcd_region(Region) or not os.path.exists(path):
        return False
    if prediction and Region not in PREDICTED_REGIONS:
        return False

    root_file = ROOT.TFile.Open(path)
    found = bool(root_file and not root_file.IsZombie()
                 and root_file.Get(abcd_histogram_path(Var, Region, Syst, rebinned, prediction)))
    if root_file:
        root_file.Close()

    return found


def get_dijet_from_abcd_file(Var, Region, Rebin=1, Syst=None, rebinned=False, prediction=False,
                             path=ABCD_FILE, clip_negative=True, verbose=True):
    """Read the data-driven dijet estimate that ABCD.py wrote, so that a plot shows exactly
    what the fit is given.

    Rebin groups the native binning by that factor, and is ignored when rebinned=True asks
    for the analysis binning instead. clip_negative clamps bins the subtraction pushed below
    zero, which a stacked plot needs.

    prediction=True takes the transfer-factor extrapolation for regions that have one, rather
    than the direct data - MC measurement. Note that it only reproduces the measurement bin for
    bin at the analysis binning the TF was measured in; asking for a finer binning applies one
    coarse TF value across each group of fine bins, which is exact only where CR1 and CR2 have
    the same shape within that group.

    Run ABCD.py first; it does not cover NN_score or Hbb_bjR_mass, which define the regions."""

    if not is_abcd_region(Region):
        raise ValueError(f"No data-driven dijet estimate for region {Region}, "
                         f"expected one of {ABCD_REGIONS + [SR_REGION]}")
    if prediction and Region not in PREDICTED_REGIONS:
        raise ValueError(f"No transfer-factor prediction for region {Region}, "
                         f"expected one of {list(PREDICTED_REGIONS)}")
    if not os.path.exists(path):
        raise FileNotFoundError(f"{path} not found. Run ABCD.py to produce it.")

    hist = get_detached_histogram(path, abcd_histogram_path(Var, Region, Syst, rebinned, prediction), "abcd")
    hist.SetDirectory(0)

    if not rebinned and Rebin != 1:
        hist = hist.Rebin(Rebin, f"dijet_{Region}_{Var}")
    hist.SetName(f"dijet_{Region}_{Var}")

    if clip_negative:
        n_clamped = zero_negative_bins(hist)
        if n_clamped and verbose:
            print(f"  [dijet] {Var} in {Region}: clamped {n_clamped} negative bin(s) to zero")

    return hist
