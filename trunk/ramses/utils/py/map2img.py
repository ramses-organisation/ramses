#!/usr/bin/env python
 
import sys
import numpy
import pylab
import os
from matplotlib import pyplot as plt
import fortranfile
from optparse import OptionParser

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
	parser.add_option("-d","--dir", dest="dir", \
			help='map directory', default=None)
	parser.add_option("-i","--iter", dest="iter", \
			help="iterator index", default=2)
	parser.add_option("-p","--proj", dest="proj", \
			help="proj_nd", default=2)
	parser.add_option("-s","--step", dest="step", \
			help="framing step", default=5)
	parser.add_option("-z","--zero", dest="zero", \
			help="Starting frame", default=0, type=int)
	parser.add_option('-k','--kind', dest="kind", \
			help="kind of plot [temp, dens, metal]", default='dens')	
	parser.add_option('-a','--autorange',dest='autorange', action='store_true', \
	       help='use automatic dynamic range (overrides min & max)', default=False)
	parser.add_option('--big-endian',dest='big_endian', action='store_true', \
	       help='input binary data is stored as big endian', default=False)
	parser.add_option('-c','--colormap',dest='cmap_str', metavar='CMAP', \
	       help='matplotlib color map to use', default="jet")
	(opts,args)=parser.parse_args()

	
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
		if (line.split('=')[0] == 'sink' and line.split('=')[1] == '.true.'):
			sink_flag = True
			
	
	print 'Projection axis: %s' % (proj_axis[proj_ind+1])
	
	# Progressbar imports/inits	
	from widgets import Percentage, Bar, ETA
	from progressbar import ProgressBar
	
	widgets = ['Working...', Percentage(), Bar(marker='#'),ETA()]
	pbar = ProgressBar(widgets=widgets, maxval = max_iter+1).start()
	
	# Looping over movie snapshots
	for i in xrange(opts.zero,int(opts.step),max_iter+1,int(opts.step)):

		infile = "%smovie%d/%s_%05d.map" % (opts.dir, int(opts.proj), opts.kind, i)		
	
		infof = open("%smovie%d/info_%05d.txt" % (opts.dir, int(opts.proj), i))
		for j, line in enumerate(infof):
			if j == 7:
				boxlen = float(line.split()[2])
			if j == 8:
				time = float(line.split()[2])
			if j> 18:
				break	
		if sink_flag:
			try:
				sinkf = open("%smovie%d/sink_%05d.txt" % (opts.dir, int(opts.proj), i))
				sink = sinkf.readline().split(',')
				sink_x = float(sink[2])/float(boxlen)
				sink_y = float(sink[3])/float(boxlen)
				sink_z = float(sink[4])/float(boxlen)
				#print sink_x, sink_y, sink_z
			except IOError:
				print "No sink file" 

	
		if(opts.outfile==None):
			if not os.path.exists("%smovie%d/pngs/" % (opts.dir, int(opts.proj))):
    				os.makedirs("%smovie%d/pngs/" % (opts.dir, int(opts.proj)))
			outfile="%smovie%d/pngs/%s_%05d.png" % (opts.dir, int(opts.proj), opts.kind, i/int(opts.step)-opts.zero)
		else:
			outfile=opts.outfile
		
		# Read image data
		f = fortranfile.FortranFile(infile)
		[time, dx, dy, dz] = f.readReals('d')
		[nx,ny] = f.readInts()
		dat = f.readReals()
		f.close()
	
		if(opts.logscale):
			dat = numpy.array(dat)+1e-12

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
		fig = plt.figure(figsize=(8,8),frameon=False)
		
		ax = fig.add_subplot(1,1,1)
		ax.set_axis_off()
		fig.add_axes(ax)	
		ax.imshow(dat, interpolation = 'nearest', cmap = opts.cmap_str,\
				vmin = plotmin, vmax = plotmax)

		# Plotting sink
		if sink_flag:
			if axis == 'x':
				ax.plot((sink_y-ycentre_frame/boxlen)/(deltay_frame/boxlen/2)*nx/2+nx/2,\
					(sink_z-zcentre_frame/boxlen)/(deltaz_frame/boxlen/2)*ny/2+ny/2,\
					marker='+',c='r',mew=1,ms=6)
			elif axis == 'y':
				ax.plot((sink_x-xcentre_frame/boxlen)/(deltax_frame/boxlen/2)*nx/2+nx/2,\
					(sink_z-zcentre_frame/boxlen)/(deltaz_frame/boxlen/2)*ny/2+ny/2,\
					marker='+',c='r',mew=1,ms=6)
			else:
				ax.plot((sink_x-xcentre_frame/boxlen)/(deltax_frame/boxlen/2)*nx/2+nx/2,\
					(sink_y-ycentre_frame/boxlen)/(deltay_frame/boxlen/2)*ny/2+ny/2,\
					marker='+',c='r',mew=1,ms=6)

		plt.axis('off') # removes axis
		plt.xlim(0,ny) # trims image to borders
		plt.ylim(0,nx)
		# corrects window extent
		extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
		plt.savefig(outfile,bbox_inches=extent,dpi=100)
		plt.close(fig)
		
		pbar.update(i) # updates progressbar

	pbar.finish()
		

if __name__ == '__main__':
	main()
