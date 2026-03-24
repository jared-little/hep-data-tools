import ROOT
import numpy as np

def compute_significance(n, b, sigma):
    """Compute the significance using the Asimov formula, accounting for background uncertainty."""

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


def get_Zn_histogram(sigHistoDict,bkgHisto,versus):
    """Calculate the Zn histogram for each signal histogram compared to the background histogram, for either upper or lower cuts."""

    hZn = []
    max = 0 # For setting y-axis range
    for sigName, sig in sigHistoDict.items():
      h = ROOT.TH1D("hZn_"+sigName+"_"+versus,"",sig.GetNbinsX(),sig.GetXaxis().GetXmin(),sig.GetXaxis().GetXmax())
      if versus=="upper":
        for i in range(1,sig.GetNbinsX()+1):
          if bkgHisto.Integral(i,sig.GetNbinsX())>0:
            error = 0.0
            # print(sigName+":\t" "Cut:", h.GetBinLowEdge(i), "\tSignal Events:", sig.Integral(i,sig.GetNbinsX()), "\tBackground Events:", bkgHisto.Integral(i,sig.GetNbinsX()))
            zn = compute_significance(sig.Integral(i,sig.GetNbinsX()),bkgHisto.Integral(i,sig.GetNbinsX()),error)
            if zn > max: max = zn
            h.SetBinContent(i, zn)
          else:
            h.SetBinContent(i, 100)
      else:
        for i in range(0,sig.GetNbinsX()):
          if bkgHisto.Integral(1,sig.GetNbinsX()-i)>0:
            error = 0.0
            zn = compute_significance(sig.Integral(1,sig.GetNbinsX()-i),bkgHisto.Integral(1,sig.GetNbinsX()-i),error)
            if zn > max: max = zn
            h.SetBinContent(sig.GetNbinsX()-i, zn)
          else:
            h.SetBinContent(sig.GetNbinsX()-i, 100)

      if "X2000" in sigName: h.SetLineColor(ROOT.kOrange)
      if "X3000" in sigName: h.SetLineColor(ROOT.kCyan)
      if "X4000" in sigName: h.SetLineColor(ROOT.kViolet)
      h.SetLineWidth(4)
      h.SetLineStyle(2)
      hZn.append(h)

    return hZn, max


def get_SB_histogram(sigHistoDict,bkgHisto,versus):
    """Calculate the S/B histogram for each signal histogram compared to the background histogram, for either upper or lower cuts."""

    hSB = []
    max = 0 # For setting y-axis range
    for sigName, sig in sigHistoDict.items():
      h = ROOT.TH1D("hSB_"+sigName+"_"+versus,"",sig.GetNbinsX(),sig.GetXaxis().GetXmin(),sig.GetXaxis().GetXmax())
      if versus=="upper":
        for i in range(1,sig.GetNbinsX()+1):
          if bkgHisto.Integral(i,sig.GetNbinsX())>0:
            # print(sigName+":\t" "Cut:", h.GetBinLowEdge(i), "\tSignal Events:", sig.Integral(i,sig.GetNbinsX()), "\tBackground Events:", bkgHisto.Integral(i,sig.GetNbinsX()))
            sb = sig.Integral(i,sig.GetNbinsX()) / bkgHisto.Integral(i,sig.GetNbinsX())
            if sb > max: max = sb
            h.SetBinContent(i, sb)
          else:
            h.SetBinContent(i, 100)
      else:
        for i in range(0,sig.GetNbinsX()):
          if bkgHisto.Integral(1,sig.GetNbinsX()-i)>0:
            sb = sig.Integral(1,sig.GetNbinsX()-i) / bkgHisto.Integral(1,sig.GetNbinsX()-i)
            if sb > max: max = sb
            h.SetBinContent(sig.GetNbinsX()-i, sb)
          else:
            h.SetBinContent(sig.GetNbinsX()-i, 100)

      if "X2000" in sigName: h.SetLineColor(ROOT.kOrange)
      if "X3000" in sigName: h.SetLineColor(ROOT.kCyan)
      if "X4000" in sigName: h.SetLineColor(ROOT.kViolet)
      h.SetLineWidth(4)
      h.SetLineStyle(2)
      hSB.append(h)

    return hSB, max


def get_efficiency_selection(Histo,versus):

    h = ROOT.TH1D("hEff_"+Histo.GetName()+"_"+versus,"",Histo.GetNbinsX(),Histo.GetXaxis().GetXmin(),Histo.GetXaxis().GetXmax())
    if versus=="upper":
      for i in range(1,Histo.GetNbinsX()+1):
        # print(Histo.GetName()+":\t" "Cut:", h.GetBinLowEdge(i), "\tSelected Events:", Histo.Integral(i,Histo.GetNbinsX()))
        h.SetBinContent(i,Histo.Integral(i,Histo.GetNbinsX()))
    else:
      for i in range(0,Histo.GetNbinsX()):
        # print(Histo.GetName()+":\t" "Cut:", h.GetBinLowEdge(i), "\tSelected Events:", Histo.Integral(i,Histo.GetNbinsX()))
        h.SetBinContent(Histo.GetNbinsX()-i, Histo.Integral(1,Histo.GetNbinsX()-i))

    return h


def compute_significance_deprecated(s,b,sigma):
    """Compute the significance using the Asimov formula, accounting for background uncertainty."""
    import math
    n = s + b
    sigma = b * sigma
    if s > 0 and b > 0:
        z = math.sqrt(2 * (n * math.log (n * (b + sigma * sigma)/(b * b + n * sigma * sigma)) - b * b/(sigma * sigma) * math.log((b * b + n * sigma * sigma)/(b * (b + sigma * sigma))) ) )
    else:
        z = 0
    return z
