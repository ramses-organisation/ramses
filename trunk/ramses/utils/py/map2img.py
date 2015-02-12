#!/usr/bin/env python
 
import sys
import numpy
import os
import warnings
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
import fortranfile
from optparse import OptionParser


def label(xy, text):
	y = xy[1] + 15 # shift y-value for label so that it's below the artist
	plt.text(xy[0], y, text, ha="center",  size=14, color='white')

def main():

	# Parse command line arguments
	parser = OptionParser()
	parser.usage = "%prog [options] map_file"
	parser.add_option('-l','--logscale',dest='logscale', action='store_true', \
	       help='use log color scaling', default=False)
	parser.add_option("-o","--output",  dest="outfile", metavar="FILE", \
			help='output image file [default: <map_file>.png]', default=None)
	parser.add_option("-m","--min",  dest="min", metavar="VALUE", \
			help='min value', default=None)
	parser.add_option("-M","--max",  dest="max", metavar="VALUE", \
			help='max value', default=None)
	parser.add_option("-f","--fmin",  dest="fmin", metavar="VALUE", \
			help='frame min value', default=0, type=int)
	parser.add_option("-F","--fmax",  dest="fmax", metavar="VALUE", \
			help='frame max value', default=0, type=int)
	parser.add_option("-d","--dir", dest="dir", \
			help='map directory', default=None)
	parser.add_option("-i","--iter", dest="iter", \
			help="iterator index", default=2)
	parser.add_option("-p","--proj", dest="proj", \
			help="proj_nd", default=2)
	parser.add_option("-s","--step", dest="step", \
			help="framing step", default=5, type=int)
	parser.add_option('-k','--kind', dest="kind", \
			help="kind of plot [temp, dens, metal]", default='dens')	
	parser.add_option('-a','--autorange',dest='autorange', action='store_true', \
	       help='use automatic dynamic range (overrides min & max)', default=False)
	parser.add_option('--big-endian',dest='big_endian', action='store_true', \
	       help='input binary data is stored as big endian', default=False)
	parser.add_option('-c','--colormap',dest='cmap_str', metavar='CMAP', \
	       help='matplotlib color map to use', default="jet")
	parser.add_option('-b','--barlen',dest='barlen', metavar='VALUE', \
	       help='length of the bar (specify unit!)', default=50)
	parser.add_option('-B','--barlen_unit',dest='barlen_unit', metavar='VALUE', \
	       help='unit of the bar length [pc/kpc/Mpc]', default='pc')
	parser.add_option('-t','--time_unit',dest='time_unit', metavar='VALUE', \
	       help='unit of time [Myr/Gyr]', default='Myr')



	(opts,args)=parser.parse_args()

	if opts.barlen_unit == 'pc':
		scale_l = 1e0
	elif opts.barlen_unit == 'kpc':
		scale_l = 1e3
	elif opts.barlen_unit == 'Mpc':
		scale_l = 1e3
	else:
		print "Wrong length unit!"
		sys.exit()

	if opts.time_unit == 'Myr':
		scale_t = 1e6
	elif opts.time_unit == 'Gyr':
		scale_t = 1e9
	else:
		print "Wrong time unit!"
		sys.exit()



	
	proj_ind = int(opts.proj)-1
	
	# Searching for namelist
	for file in os.listdir(opts.dir):
		if file.endswith(".nml"):
			namelist = opts.dir + file
	try:
		nmlf = open(namelist)
	except IOError:
		print "No namelist found! Aborting!"
		sys.exit()
	

	sink_flag = False
	# Loading parameters from the namelist
        for i, line in enumerate(nmlf):
		if line.split('=')[0] == 'xcentre_frame':
			xcentre_frame = float(line.split('=')[1].split(',')[0+4*proj_ind])
		if line.split('=')[0] == 'ycentre_frame':
			ycentre_frame = float(line.split('=')[1].split(',')[0+4*proj_ind])
		if line.split('=')[0] == 'zcentre_frame':
			zcentre_frame = float(line.split('=')[1].split(',')[0+4*proj_ind])
		if line.split('=')[0] == 'deltax_frame':
			deltax_frame = float(line.split('=')[1].split(',')[0+2*proj_ind])
		if line.split('=')[0] == 'deltay_frame':
			deltay_frame = float(line.split('=')[1].split(',')[0+2*proj_ind])
		if line.split('=')[0] == 'deltaz_frame':
			deltaz_frame = float(line.split('=')[1].split(',')[0+2*proj_ind])
		if line.split('=')[0] == 'proj_axis':
			proj_axis = line.split('=')[1]
		if line.split('=')[0] == 'imovout':
			max_iter = int(line.split('=')[1])
		if (line.split('=')[0] == 'sink') and (line.split('=')[1][:-1] == '.true.'):
			sink_flag = True
	
	print 'Projection axis: %s' % (proj_axis[proj_ind+1])
	
	# Progressbar imports/inits
	try:
		from widgets import Percentage, Bar, ETA
		from progressbar import ProgressBar
		progressbar_avail = True
	except ImportError:
		progressbar_avail = False
  
	if (int(opts.fmax) > 0):
		max_iter=int(opts.fmax)
	
	if progressbar_avail:
	  widgets = ['Working...', Percentage(), Bar(marker='#'),ETA()]
	  pbar = ProgressBar(widgets=widgets, maxval = max_iter+1).start()
	else:
		print 'Working!'
        
	# Looping over movie snapshots
	for i in xrange(int(opts.fmin)+int(opts.step),max_iter+1,int(opts.step)):

		infile = "%smovie%d/%s_%05d.map" % (opts.dir, int(opts.proj), opts.kind, i)		
	
		infof = open("%smovie%d/info_%05d.txt" % (opts.dir, int(opts.proj), i))
		for j, line in enumerate(infof):
			if j == 7:
				boxlen = float(line.split()[2])
			if j == 8:
				time = float(line.split()[2])
			if j == 15:
				unit_l = float(line.split()[2])
			if j == 16:
				unit_d = float(line.split()[2])
			if j == 17:
				unit_t = float(line.split()[2])
			if j> 18:
				break
    
		unit_m = unit_d*unit_l**3/2e33 # in MSun
		
		if sink_flag:
			sink_file = "%smovie1/sink_%05d.txt" % (opts.dir, i)
			try:
				with warnings.catch_warnings():
					warnings.simplefilter("ignore")
					sink_m,sink_x,sink_y,sink_z = numpy.loadtxt(sink_file, delimiter=',',usecols=(1,2,3,4),unpack=True)/float(boxlen)
					sink_m *= unit_m*float(boxlen)
					sink_m = numpy.log10(sink_m)
					plot_sinks = True
			except ValueError:
				plot_sinks = False
			except IOError:
				print "No sink file" 

		if(opts.outfile==None):
			if not os.path.exists("%smovie%d/pngs/" % (opts.dir, int(opts.proj))):
    				os.makedirs("%smovie%d/pngs/" % (opts.dir, int(opts.proj)))
			outfile="%smovie%d/pngs/%s_%05d.png" % (opts.dir, int(opts.proj), opts.kind, i/int(opts.step)-int(opts.fmin))
		else:
			outfile=opts.outfile
		
		# Read image data
		f = fortranfile.FortranFile(infile)
		[time, dx, dy, dz] = f.readReals('d')
		[nx,ny] = f.readInts()
		dat = f.readReals()
		f.close()
	
		if(opts.logscale):
			dat = numpy.array(dat)

		rawmin = numpy.amin(dat)
		rawmax = numpy.amax(dat)

		# Bounds
		if opts.min==None:
			plotmin = rawmin
		else:
			plotmin = float(opts.min)

		if opts.max==None:
			plotmax = rawmax
		else:
			plotmax = float(opts.max)

		# Log scale?
		if(opts.logscale):
			dat = numpy.log10(dat)
			rawmin = numpy.log10(rawmin)
			rawmax = numpy.log10(rawmax)
			plotmin = numpy.log10(plotmin)
			plotmax = numpy.log10(plotmax)

		# Auto-adjust dynamic range?
		if(opts.autorange):
			# Overrides any provided bounds
			NBINS = 200
			# Compute histogram
			(hist,bins) = numpy.histogram(dat, NBINS, (rawmin,rawmax), normed=True)
			chist = numpy.cumsum(hist); chist = chist / numpy.amax(chist)
			# Compute black and white point
			clip_k = chist.searchsorted(0.05)
			plotmin = bins[clip_k]
			plotmax = rawmax

		#if(plotmax-plotmin>0):
		#	dat = numpy.clip((dat-plotmin)/(plotmax-plotmin), 0.0, 1.0)
		#else:
		#	dat = 0.5*dat/plotmax
    
		axis = proj_axis[proj_ind+1]
	
		# Reshape data to 2d
		dat = dat.reshape(ny,nx)
	
		# Plotting
		fig = plt.figure(figsize=(8,8*ny/nx),frameon=False)
		
		ax = fig.add_subplot(1,1,1)
		ax.set_axis_off()
		fig.add_axes(ax)	
		ax.imshow(dat, interpolation = 'nearest', cmap = opts.cmap_str,\
				vmin = plotmin, vmax = plotmax)

		# Plotting sink
		if (sink_flag and plot_sinks):
			if axis == 'x':
				ax.scatter((sink_y-ycentre_frame/boxlen)/(deltay_frame/boxlen/2)*nx/2+nx/2,\
					(sink_z-zcentre_frame/boxlen)/(deltaz_frame/boxlen/2)*ny/2+ny/2,\
					marker='+',c='r',s=4*sink_m**2)
			elif axis == 'y':
				ax.scatter((sink_x-xcentre_frame/boxlen)/(deltax_frame/boxlen/2)*nx/2+nx/2,\
					(sink_z-zcentre_frame/boxlen)/(deltaz_frame/boxlen/2)*ny/2+ny/2,\
					marker='+',c='r',s=4*sink_m**2)
			else:
				ax.scatter((sink_x-xcentre_frame/boxlen)/(deltax_frame/boxlen/2)*nx/2+nx/2,\
					(sink_y-ycentre_frame/boxlen)/(deltay_frame/boxlen/2)*ny/2+ny/2,\
					marker='+',c='r',s=4*sink_m**2)

		patches = []
		rect = mpatches.Rectangle((nx/100,ny/100), opts.barlen*nx/(dx/boxlen*unit_l*3.24e-19/scale_l),10)
		label([nx/100+opts.barlen*nx/(dx/boxlen*unit_l*3.24e-19/scale_l)/2,ny/100],"%d %s" % (opts.barlen, opts.barlen_unit))
		patches.append(rect)

		plt.axis('off') # removes axis
		plt.xlim(0,nx) # trims image to borders
		plt.ylim(0,ny)
		ax.text(0.95, 0.95, '%.1f %s' % (time*unit_t/86400/365.25/scale_t, opts.time_unit),
				        verticalalignment='bottom', horizontalalignment='right',
				        transform=ax.transAxes,
				        color='white', fontsize=14)
		collection = PatchCollection(patches, facecolor='white')
		ax.add_collection(collection)


		# corrects window extent
		extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
		plt.savefig(outfile,bbox_inches=extent,dpi=100)
		plt.close(fig)
		
		if progressbar_avail:
			pbar.update(i) # updates progressbar

	if progressbar_avail:
		pbar.finish()
	else:
		print 'Finished!'
if __name__ == '__main__':
	main()
