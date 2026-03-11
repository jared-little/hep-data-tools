import ROOT
import csv


def CalculateAndSaveSignalEfficiency(sig, var):
    '''
    Calculate the signal efficiency for a given signal and variable,
    and save the results to a CSV file.
    '''
    hsig = getSignalHistogram(sig, var)
    
    total_signal_events = hsig.Integral(0, hsig.GetNbinsX() + 1)
    if total_signal_events == 0:
      return []

    efficiencies = []
    num_bins = hsig.GetNbinsX()
    
    for i in range(1, num_bins + 1):
      passing_events = hsig.Integral(i, num_bins + 1)
      efficiency = passing_events / total_signal_events
      cut_value = hsig.GetBinLowEdge(i)
      efficiencies.append((cut_value, efficiency))

    output = f"SignalEfficiency_{sig}.csv"
    with open(output, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['NN_score cut value', 'Sig Efficiency'])
        for cut_value, efficiency in efficiencies:
            writer.writerow([f"{cut_value:.2f}", f"{efficiency:.4f}"])
            
    print(f"Signal efficiencies successfully written to {output}")