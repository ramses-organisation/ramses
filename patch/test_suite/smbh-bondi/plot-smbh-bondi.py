"""Plot results of sink accretion test"""
from matplotlib import pyplot as plt
from matplotlib import ticker


def main():
    """Main function"""

    # read in the data
    times = []
    rates = []
    with open('data4plot.dat', 'r') as datafile:
        for line in datafile:
            times.append(float(line.split(' ')[0]))
            rates.append(float(line.split(' ')[1]))

    # generating a simple plot
    fig = plt.figure(figsize=(8, 7), dpi=100)
    axis = fig.add_subplot(1, 1, 1)

    # plot as line and scatter
    axis.plot(times, rates, color='C0')
    axis.scatter(times, rates, color='C0', marker='+')

    # labels and eye-candy
    axis.set_xlabel('Time [kyr]', fontsize=20)
    axis.set_ylabel(r'dot(M$_\mathrm{Bondi}$) [Msol/yr]', fontsize=20)
    axis.set_xlim([0, 60])
    axis.set_ylim([1e-3, 3e-3])
    axis.yaxis.set_minor_formatter(ticker.ScalarFormatter())
    axis.set_yscale('log')
    axis.tick_params(axis='both', color='k', labelcolor='k',
                     which='both', labelsize=16,
                     bottom=True, top=True, left=True, right=True)
    axis.xaxis.label.set_color('k')
    axis.yaxis.label.set_color('k')

    # adjusting plot
    plt.subplots_adjust(left=0.22, right=0.95,
                        bottom=0.15, top=0.95,
                        wspace=0.25, hspace=0.1)

    plt.savefig('smbh-bondi.pdf')


if __name__ == '__main__':
    main()
