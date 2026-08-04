"""Compare the NN classifier ROC against the cut-based tagger working points,
using full-statistics HHARD output (no train/test/validation split).

GN3X has per-event ntuples (NN_score, combinedWeight, and per-jet GN2X_count/
Wtagger_count arrays), so its DNN curve and cut-based WP markers are both
computed directly from them. GN2X only has histograms (no ntuple), and its
GN2X_count histogram is filled per-jet -- the per-event jet combinatorics
needed for the WP markers are already lost there -- so GN2X only gets a DNN
curve, built from its binned NN_score histogram.

Both taggers' "Preselection" region already requires passing a loosest-WP
H-tag cut (has_GN2X == true, using each config's own tagger), so most of a
tagger's real improvement is baked into Preselection itself, before the NN
score is even evaluated -- and it's invisible if efficiencies are computed
relative to each tagger's own Preselection sample, since that cancels out
exactly like a uniform weight rescale. Efficiencies are instead computed
relative to the shared "All"-region weight (identical for both taggers,
since it's pre-tagger-cut), so the Preselection-level suppression shows up
in the plotted curves.
"""

import glob
import os
from types import SimpleNamespace

import uproot
import awkward as ak
import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)  # no GUI event loop, we only ever write to file

# -------------------------
# CONFIG
# -------------------------
region = "Preselection"
all_region = "All"  # pre-tagger-cut reference stage, shared denominator for both taggers

gn3_ntuple_dir = "/data/jlittle/HHARDout/gn3/ntuples"
gn2_hist_dir = "/data/jlittle/HHARDout/gn2/histograms"
gn3_hist_dir = "/data/jlittle/HHARDout/gn3/histograms"  # only used to look up the All-region reference weight

bkg_processes = ["dijet", "ttbar", "Vjets", "VV", "top"]
# Signal mass points included in NN training/testing (config.py drops X1000_S500 as uncalibrated)
signal_masses = ["X2000_S1000", "X2000_S1500", "X2500_S1500", "X3000_S1500",
                  "X3000_S2000", "X3000_S2500", "X4000_S2000"]

nn_branch = "NN_score"
gn_branch = "GN2X_count"  # branch name is unchanged between the GN2X and GN3X flattenings
wtag_branch = "Wtagger_count"
weight_branch = "combinedWeight"

thresholds = np.linspace(0.0001, 0.9999, 9998)

# The *_count branches hold, per jet, how many of the tagger's nested working
# points that jet passes. The WPs are nested, so "count >= k" selects the k-th
# loosest WP. The list below therefore runs loosest -> tightest, and entry k-1
# names the WP that "count >= k" corresponds to.
gn3x_wps = ["0p6", "0p4", "0p25", "0p15", "0p08"]  # GN3XPV01_hbb
wtag_wps = ["50", "60", "70", "80", "90"]  # ANN W-tagger

# Per-WP color, and a filled/open marker pair (index picks the shape, tagger picks fill).
colors = [ROOT.kRed, ROOT.kOrange + 1, ROOT.kSpring - 6, ROOT.kGreen + 2,
          ROOT.kAzure + 1, ROOT.kViolet - 3, ROOT.kPink + 1, ROOT.kGray + 2]
marker_pairs = [(20, 24), (21, 25), (22, 26), (23, 32), (33, 27), (34, 28), (29, 30)]

# GN2X has no ntuple (only histograms), so it gets a DNN curve only -- no cut-based
# WP markers, since those need per-jet combinatorics a 1D histogram can't provide.
taggers = [
    dict(label="GN2X", kind="histogram", line_style=ROOT.kDashed, gn_wps=[]),
    dict(label="GN3X", kind="ntuple", line_style=ROOT.kSolid, gn_wps=gn3x_wps, filled=True),
]


# -------------------------
# GN3X: per-event ntuples
# -------------------------

def _gn3_files(process=None, mass=None):
    """All campaign files (mc23a/d/e) for one background process or signal mass point."""
    if mass is not None:
        pattern = f"*_bbVV_XHS_{mass}_bbWW_allhad_{region}.root"
    else:
        pattern = f"*_bbVV_{process}_{region}.root"
    files = sorted(glob.glob(os.path.join(gn3_ntuple_dir, pattern)))
    if not files:
        raise RuntimeError(f"No GN3X ntuples matching {pattern}")
    return files


def load_gn3_sample(names, is_signal):
    """Concatenate one class (signal or background) across all processes/mass points and campaigns."""
    branches = [nn_branch, gn_branch, wtag_branch, weight_branch]
    paths = []
    for name in names:
        paths.extend(_gn3_files(mass=name) if is_signal else _gn3_files(process=name))

    data = uproot.concatenate([f"{p}:bbVV_data" for p in paths], filter_name=branches, library="ak")
    # Sum weights in float64; float32 loses precision over millions of events.
    weight = ak.to_numpy(data[weight_branch]).astype(np.float64)
    return SimpleNamespace(
        score=ak.to_numpy(data[nn_branch]).astype(np.float64),  # one per event
        weight=weight,
        total_weight=weight.sum(),
        gn=data[gn_branch],  # jagged: one entry per large-R jet
        wtag=data[wtag_branch],
    )


def load_gn3():
    print(f"Loading GN3X ntuples from {gn3_ntuple_dir}...")
    sig = load_gn3_sample(signal_masses, is_signal=True)
    bkg = load_gn3_sample(bkg_processes, is_signal=False)
    sig_ref = region_total_weight(gn3_hist_dir, signal_masses, True, all_region)
    bkg_ref = region_total_weight(gn3_hist_dir, bkg_processes, False, all_region)
    print(f"Preselection weighted Signal: {sig.total_weight:.2f} (of {sig_ref:.2f} in {all_region}), "
          f"Bkg: {bkg.total_weight:.2f} (of {bkg_ref:.2f} in {all_region})")
    return sig, bkg, sig_ref, bkg_ref


def weight_fraction_above(sample, thresholds, reference_weight):
    """Fraction of reference_weight with score > t, evaluated for every t at once."""
    order = np.argsort(sample.score)
    score, weight = sample.score[order], sample.weight[order]

    weight_upto = np.concatenate([[0.0], np.cumsum(weight)])
    n_at_or_below = np.searchsorted(score, thresholds, side="right")

    return (sample.total_weight - weight_upto[n_at_or_below]) / reference_weight


def cut_efficiency(sample, gn_wp, wtag_wp, reference_weight):
    """Fraction of reference_weight passing: one H-tagged jet, plus >=2 W-tagged jets among the rest."""
    gn_pass = sample.gn >= gn_wp
    has_h = ak.any(gn_pass, axis=1)

    # H candidate is the highest-scoring jet, which passes whenever any jet does.
    h_idx = ak.where(has_h, ak.argmax(sample.gn, axis=1), -1)
    is_h = (ak.local_index(sample.gn) == h_idx) & has_h

    w_count = ak.sum((sample.wtag >= wtag_wp) & ~is_h, axis=1)
    event_pass = ak.to_numpy(has_h & (w_count >= 2))

    return sample.weight[event_pass].sum() / reference_weight


def nn_roc(sig, bkg, sig_ref, bkg_ref):
    """(Signal efficiency, background rejection) along the NN score scan, relative to the All-region weight."""
    eff_s = weight_fraction_above(sig, thresholds, sig_ref)
    eff_b = weight_fraction_above(bkg, thresholds, bkg_ref)
    keep = eff_b > 0
    return eff_s[keep], 1 / eff_b[keep]


def wp_points(sig, bkg, gn_wp, sig_ref, bkg_ref):
    """One (sig eff, bkg rejection) point per W-tagger WP, at a fixed GN working point."""
    sig_eff, bkg_rej = [], []
    for wtag_wp in range(1, len(wtag_wps) + 1):
        bkg_eff = cut_efficiency(bkg, gn_wp, wtag_wp, bkg_ref)
        if bkg_eff > 0:
            sig_eff.append(cut_efficiency(sig, gn_wp, wtag_wp, sig_ref))
            bkg_rej.append(1 / bkg_eff)
    return sig_eff, bkg_rej


# -------------------------
# GN2X: histograms only
# -------------------------

def _sum_histogram(names, is_signal, hist_name, hist_dir, region_name):
    """Sum one 1D histogram over all campaigns (mc23a/d/e) and processes/mass points."""
    total_vals, edges = None, None
    for name in names:
        if is_signal:
            file_pattern = f"*_XHS_{name}_bbWW_allhad.root"
            sample_name = f"XHS_{name}_bbWW_allhad"
        else:
            file_pattern = f"*_{name}.root"
            sample_name = name

        files = sorted(glob.glob(os.path.join(hist_dir, file_pattern)))
        if not files:
            raise RuntimeError(f"No histograms matching {file_pattern} in {hist_dir}")

        for path in files:
            with uproot.open(path) as f:
                hist_path = f"{region_name}/bbVVSplitHadAnalysis_13p6TeV_{sample_name}/{hist_name}"
                vals, e = f[hist_path].to_numpy()
            if total_vals is None:
                total_vals, edges = vals.copy(), e
            else:
                total_vals += vals
    return total_vals, edges


def region_total_weight(hist_dir, names, is_signal, region_name):
    """Total weight in a given region, summed across processes/mass points and campaigns."""
    vals, _ = _sum_histogram(names, is_signal, "NLargeRjets", hist_dir, region_name)
    return float(vals.sum())


def load_gn2():
    print(f"Loading GN2X histograms from {gn2_hist_dir}...")
    sig_vals, _ = _sum_histogram(signal_masses, True, nn_branch, gn2_hist_dir, region)
    bkg_vals, _ = _sum_histogram(bkg_processes, False, nn_branch, gn2_hist_dir, region)
    sig_ref = region_total_weight(gn2_hist_dir, signal_masses, True, all_region)
    bkg_ref = region_total_weight(gn2_hist_dir, bkg_processes, False, all_region)
    print(f"Preselection weighted Signal: {sig_vals.sum():.2f} (of {sig_ref:.2f} in {all_region}), "
          f"Bkg: {bkg_vals.sum():.2f} (of {bkg_ref:.2f} in {all_region})")
    return sig_vals, bkg_vals, sig_ref, bkg_ref


def hist_roc(sig_vals, bkg_vals, sig_ref, bkg_ref):
    """(Signal efficiency, background rejection) from binned NN-score histograms, relative to the All-region weight."""
    def tail_weight(vals):
        cum = np.concatenate([[0.0], np.cumsum(vals)])
        return vals.sum() - cum[:-1]  # weight with score >= each bin's left edge

    eff_s = tail_weight(sig_vals) / sig_ref
    eff_b = tail_weight(bkg_vals) / bkg_ref
    keep = eff_b > 0
    return eff_s[keep], 1 / eff_b[keep]


# -------------------------
# Plotting
# -------------------------

def plot_roc(taggers, output_name):
    canvas = ROOT.TCanvas("c", "ROC", 900, 700)
    canvas.SetLogy()
    canvas.SetGrid()

    # One compact legend per tagger, placed side by side, so each tagger's
    # entries sit in their own column instead of interleaving row-major.
    col_width, col_gap, row_height = 0.16, 0.02, 0.035
    x_start, y_top = 0.58, 0.88

    legends = []
    for tagger_i, tagger in enumerate(taggers):
        n_rows = 2 + len(tagger["gn_wps"])  # header + DNN entry + WP entries
        x1 = x_start + tagger_i * (col_width + col_gap)
        leg = ROOT.TLegend(x1, y_top - n_rows * row_height, x1 + col_width, y_top)
        leg.SetHeader(tagger["label"], "C")
        leg.SetTextSize(0.028)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        legends.append(leg)

    # Compute every curve/marker set up front so the y-axis can be sized to fit
    # all of them, rather than auto-scaling off whichever graph draws first.
    curves = []
    all_rej = []
    for tagger in taggers:
        if tagger["kind"] == "ntuple":
            sig, bkg, sig_ref, bkg_ref = load_gn3()
            eff_s, rej_b = nn_roc(sig, bkg, sig_ref, bkg_ref)
        else:
            sig_vals, bkg_vals, sig_ref, bkg_ref = load_gn2()
            eff_s, rej_b = hist_roc(sig_vals, bkg_vals, sig_ref, bkg_ref)
        all_rej.append(rej_b)

        wps = []
        for i, wp_name in enumerate(tagger["gn_wps"]):
            sig_eff, bkg_rej = wp_points(sig, bkg, i + 1, sig_ref, bkg_ref)
            if sig_eff:
                wps.append((wp_name, i, sig_eff, bkg_rej))
                all_rej.append(np.asarray(bkg_rej))

        curves.append((tagger, eff_s, rej_b, wps))

    rej_min = min(r.min() for r in all_rej)
    rej_max = max(r.max() for r in all_rej)
    axis_min = 10 ** np.floor(np.log10(rej_min))
    axis_max = 10 ** np.ceil(np.log10(rej_max))

    graphs = []  # keep references alive, ROOT will not draw garbage-collected graphs

    for tagger_i, (tagger, eff_s, rej_b, wps) in enumerate(curves):
        legend = legends[tagger_i]

        nn_graph = ROOT.TGraph(len(eff_s), np.asarray(eff_s, "f8"), np.asarray(rej_b, "f8"))
        nn_graph.SetLineColor(ROOT.kBlue)
        nn_graph.SetLineWidth(3)
        nn_graph.SetLineStyle(tagger["line_style"])
        if tagger_i == 0:
            nn_graph.GetXaxis().SetLimits(0.1, 1.0)
            nn_graph.SetTitle(";#varepsilon_{s};1 / #varepsilon_{b}")
            nn_graph.SetMinimum(axis_min)
            nn_graph.SetMaximum(axis_max)
            nn_graph.Draw("AL")
        else:
            nn_graph.Draw("L SAME")
        legend.AddEntry(nn_graph, "DNN Classifier", "l")
        graphs.append(nn_graph)

        for wp_name, i, sig_eff, bkg_rej in wps:
            filled_style, open_style = marker_pairs[i % len(marker_pairs)]
            graph = ROOT.TGraph(len(sig_eff), np.array(sig_eff, "f8"), np.array(bkg_rej, "f8"))
            graph.SetMarkerStyle(filled_style if tagger["filled"] else open_style)
            graph.SetMarkerColor(colors[i % len(colors)])
            graph.Draw("P SAME")
            legend.AddEntry(graph, wp_name, "p")
            graphs.append(graph)

    for legend in legends:
        legend.Draw()
    canvas.SaveAs(output_name)


plot_roc(taggers, "weighted_roc.pdf")
