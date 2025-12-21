import matplotlib.pyplot as plt

# Plotting the changes in P0 and T0 along the stations
def plot_T0P0_vs_Stations(thermo):
    # Consolidating T0 and P0 data
    station_T0 = [a.T0 for a in vars(thermo).values()]
    station_P0 = [a.P0 for a in vars(thermo).values()]

    # Setting up the plot
    fig, axLeft = plt.subplots()
    plt.title("Station Total Temperatures and Pressures")
    plt.subplots_adjust(left=0.15, bottom=0.1, right=0.85, top=0.9, wspace=0, hspace=0)

    # Plotting the left axis
    axLeft.plot([0,1.5,2,2.5,3,4,4.5,5,6,7,8], station_T0, '-ob', label="Temperatures")
    axLeft.set_ylabel("Total Temperatures")
    axLeft.set_xlabel("Station Number")
    plt.legend(loc="upper left")

    # Plotting the right axis
    axRight = axLeft.twinx()
    axRight.plot([0,1.5,2,2.5,3,4,4.5,5,6,7,8], station_P0, '-or', label="Pressures")
    axRight.set_ylabel("Total Pressure")
    plt.legend(loc="upper right")

    plt.show()
