import ROOT

from utilities.GetHistograms import get_bkg_histogram, get_signal_histogram

# Same background set and region names as ABCD.py, so the two stay consistent.
BKG_NAMES = ["dijet", "ttbar", "top", "Vjets", "VV"]

CR_REGIONS = ["CR0", "CR1", "CR2"]
VR_REGIONS = ["VR1", "VR2"]


def _yield(hist):
    """Full-range yield, including under/overflow so the result doesn't depend on Var's binning range."""
    return hist.Integral(0, hist.GetNbinsX() + 1)


def get_region_yields(regions, signal, var="NN_score", campaigns=None):
    """Return {process: {region: yield}} for signal + each background, in the given regions."""

    if campaigns is None:
        campaigns = ["mc23a", "mc23d", "mc23e"]

    yields = {signal: {r: _yield(get_signal_histogram(signal, var, r, campaigns=campaigns)) for r in regions}}
    for bkg in BKG_NAMES:
        yields[bkg] = {r: _yield(get_bkg_histogram(bkg, var, r, campaigns=campaigns)) for r in regions}

    return yields


def print_region_table(regions, signal, var="NN_score", campaigns=None, title=""):
    """Print signal and background yields and their fractional contribution to each region's total,
    in the style of Tables 7.2/7.3 of Hsuan Chu's thesis."""

    yields = get_region_yields(regions, signal, var, campaigns)
    totals = {r: sum(yields[p][r] for p in yields) for r in regions}

    label_width = 16
    col_width = 22
    rule = "-" * (label_width + col_width * len(regions))

    print(f"\n{title}")
    print(rule)
    print(f"{'Process':<{label_width}}" + "".join(f"{r:>{col_width}}" for r in regions))
    print(rule)

    signal_label = signal.replace("XHS_", "")
    for process in yields:
        label = signal_label if process == signal else process
        row = f"{label:<{label_width}}"
        for r in regions:
            y = yields[process][r]
            pct = 100 * y / totals[r] if totals[r] else 0
            row += f"{f'{y:.2f} ({pct:.2f}%)':>{col_width}}"
        print(row)

    print(rule)
    print(f"{'Total':<{label_width}}" + "".join(f"{totals[r]:>{col_width}.2f}" for r in regions))


if __name__ == "__main__":

    ROOT.gROOT.SetBatch(True)

    campaigns = ["mc23a", "mc23d", "mc23e"]
    signal = "XHS_X4000_S2000"

    print_region_table(CR_REGIONS, signal, campaigns=campaigns,
                        title="Signal and background yields and fractional contributions in each CR")
    print_region_table(VR_REGIONS, signal, campaigns=campaigns,
                        title="Signal and background yields and fractional contributions in each VR")
