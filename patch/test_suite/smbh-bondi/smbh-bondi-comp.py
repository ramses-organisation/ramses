"""Extract values for plotting"""


def extract_comparison_data():
    """Extract data necesary for comparison
    i.e. the final accretion rate from output_00002"""

    with open('output_00002/sink_00002.csv') as infile:
        for line in infile:
            if line.split(' ')[1] != '#':
                with open('data.dat', 'w') as outfile:
                    outfile.write('{}\n'.format(line.split(',')[12][4:]))

    return


def extract_acc_rates():
    """Extract accretion rates for plotting scripts"""

    # containers
    times = []
    rates = []

    with open('smbh-bondi.nml') as nml:
        for line in nml:
            # loading unit of time from the namelist
            if line[0:10] == 'units_time':
                unit_t = float(line.split(' ')[-4].replace('d', 'e'))

    with open('log') as logfile:
        for line in logfile:
            # reading in the timestep
            # and scaling to kyr
            if line[0:5] == ' Fine':
                times.append(float(line[22:33])*unit_t/365./86400./1.e3)
            if line[0:3] == '  1':
                # reading in the accretion rate
                rates.append(float(line[20:33]))

    with open('data4plot.dat', 'w') as outfile:
        for i in xrange(len(rates)):
            outfile.write('{:e} {:e}\n'.format(times[i+1], rates[i]))

    return


def main():
    """Main function"""
    extract_comparison_data()
    extract_acc_rates()


if __name__ == '__main__':
    main()
