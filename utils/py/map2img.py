"""Script to create image out of amr2map/part2map outputs"""
import os
import sys
from argparse import ArgumentParser
import numpy
import fortranfile
from matplotlib import pyplot as plt


def main():

    # Parse command line
    parser = ArgumentParser(description='Script to create image out of amr2map/part2map outputs')
    parser.add_argument('-l', '--logscale', dest='logscale', action='store_true', \
                        help='use log color scaling',
                        default=False)
    parser.add_argument('-i', '--input', dest='infile', metavar='FILE', \
                        help='input binary image file',
                        type=str, default=None)
    parser.add_argument('-o', '--output', dest='outfile', metavar='FILE', \
                        help='output image file',
                        type=str, default=None)
    parser.add_argument('-m', '--min', dest='min', metavar='VALUE', \
                        help='min value',
                        type=float, default=None)
    parser.add_argument('-M', '--max', dest='max', metavar='VALUE', \
                        help='max value',
                        type=float, default=None)
    parser.add_argument('-a', '--autorange', dest='autorange', action='store_true', \
                        help='use automatic dynamic range (overrides min & max)',
                        default=False)
    parser.add_argument('-c', '--colormap', dest='cmap_str', metavar='VALUE', \
                        help='matplotlib color map to use',
                        type=str, default='gray')

    args = parser.parse_args()

    if args.infile is None or not os.path.isfile(args.infile):
        print 'Incorrect input. Exiting'
        sys.exit()

    # Read image data
    fileobj = fortranfile.FortranFile(args.infile)
    [time, delta_x, delta_y, delta_z] = fileobj.readReals('d')
    del time, delta_x, delta_y, delta_z
    [nx, ny] = fileobj.readInts()
    dat = fileobj.readReals()
    fileobj.close()

    rawmin = numpy.amin(dat)
    rawmax = numpy.amax(dat)

    # Bounds
    if args.min is None:
        plotmin = rawmin
    else:
        plotmin = float(args.min)

    if args.max is None:
        plotmax = rawmax
    else:
        plotmax = float(args.max)

    # Log scale?
    if args.logscale:
        dat = numpy.log10(dat)
        rawmin = numpy.log10(rawmin)
        rawmax = numpy.log10(rawmax)
        plotmin = numpy.log10(plotmin)
        plotmax = numpy.log10(plotmax)

    # Auto-adjust dynamic range?
    if args.autorange:
        # Overrides any provided bounds
        nbins = 200
        # Compute histogram
        (hist, bins) = numpy.histogram(dat, nbins, (rawmin, rawmax), normed=True)
        chist = numpy.cumsum(hist)
        chist = chist / numpy.amax(chist)
        # Compute black and white point
        clip_k = chist.searchsorted(0.05)
        plotmin = bins[clip_k]
        plotmax = rawmax


    # Reshape data to 2d
    dat = dat.reshape(ny, nx)

    # Plotting
    fig = plt.figure(frameon=False)
    fig.set_size_inches(nx/100, ny/100)
    axis = plt.Axes(fig, [0., 0., 1., 1.])
    axis.set_axis_off()
    fig.add_axes(axis)
    axis.imshow(dat, interpolation='nearest', cmap=args.cmap_str,
                vmin=plotmin, vmax=plotmax, aspect='auto')

    plt.axis('off') # removes axis
    plt.xlim(0, nx) # trims image to borders
    plt.ylim(0, ny)

    # corrects window extent
    axis.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    plt.savefig(args.outfile, dpi=100)
    plt.close(fig)


if __name__ == '__main__':
    main()
