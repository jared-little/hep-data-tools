from utilities.GetHistograms import get_data_histogram


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


def get_yield(var, region, rebin=1, campaigns=None, selection=None):
    if campaigns is None:
        campaigns = ["mc23a"]

    years = _campaigns_to_data_years(campaigns)
    hist_data = get_data_histogram(var, region, rebin, campaigns=years, blind_data=False)

    if selection is None:
        print(f"Yield: {hist_data.Integral()}")
        return

    if not isinstance(selection, (tuple, list)) or len(selection) != 2:
        raise ValueError("selection must be a tuple/list like (x_min, x_max)")

    x_min, x_max = selection
    if x_min > x_max:
        raise ValueError("selection must satisfy x_min <= x_max")

    x_axis = hist_data.GetXaxis()
    bin_min = x_axis.FindBin(x_min)
    bin_max = x_axis.FindBin(x_max)
    selected_yield = hist_data.Integral(bin_min, bin_max)

    # print(f"Yield [{x_min}, {x_max}]: {selected_yield}")

    return selected_yield, x_min, x_max


def make_table(regions, campaigns, selection=(0.0, 1.0)):
    """Make a table of yields for different regions and selections."""

    for region in regions:
        events, x_min, x_max = get_yield("NN_score", region, rebin=2, campaigns=campaigns, selection=selection)
        print(f"{region}:\t{events} events in [{x_min}, {x_max}]")


if __name__ == "__main__":
    campaigns = ["mc23a", "mc23d", "mc23e"]
    regions = ["Preselection", "Preselection_CR0", "Preselection_CR2", "Preselection_VR2"]
    # selection = (0.5, 0.95) # For example, to get yields in the NN_score < 0.95 region
    for x in [0.0, 0.2, 0.3, 0.4, 0.5]:
        make_table(regions, campaigns, selection=(x, .95))
    # make_table(regions, campaigns, selection)
