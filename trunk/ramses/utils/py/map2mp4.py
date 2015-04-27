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
import subprocess


def label(xy, text):
	y = xy[1] + 15 # shift y-value for label so that it's below the artist
	plt.text(xy[0], y, text, ha="center",  size=14, color='white')

def load_map(opts,k,i):
	kind = [item for item in opts.kind.split(' ')]
	# define map path
	map_file = "%s/movie%d/%s_%05d.map" % (opts.dir, int(opts.proj), kind[k], i)
		
	# read image data
	f = fortranfile.FortranFile(map_file)
	[time, dx, dy, dz] = f.readReals('d')
	[nx,ny] = f.readInts()
	dat = f.readReals()
	f.close()

	return dat

def load_sink(opts,i):
	# setting dummy values
	sink_id, sink_m, sink_x, sink_y, sink_z = [-1 for a in xrange(5)]
	# defnining sink path
	sink_file = "%s/movie1/sink_%05d.txt" % (opts.dir, i)
	try:
		with warnings.catch_warnings(): # load sink id, mass and position
			warnings.simplefilter("ignore")
			sink_id, sink_m,sink_x,sink_y,sink_z = numpy.loadtxt(sink_file, delimiter=',',usecols=(0,1,2,3,4),unpack=True)
			plot_sinks = True
	except ValueError: # catch if no sinks exist (for sink creation)
		plot_sinks = False
	except IOError: # sink file missing
		print "No sink file"
		plot_sinks = False
	
	return plot_sinks, sink_id, sink_m, [sink_x, sink_y, sink_z]

def load_namelist_info(opts):
	proj_list = [int(item) for item in opts.proj.split(' ')]
	proj_ind = int(proj_list[0])-1

	namelist = opts.dir + '/output_00001/namelist.txt'
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
		if line.split('=')[0] == 'boxlen':
			boxlen = float(line.split('=')[1])
		if line.split('=')[0] == 'proj_axis':
			proj_axis = line.split('=')[1]
		if line.split('=')[0] == 'nw_frame':
			nx = int(line.split('=')[1])
		if line.split('=')[0] == 'nh_frame':
			ny = int(line.split('=')[1])
		if line.split('=')[0] == 'imovout':
			max_iter = int(line.split('=')[1])
		if (line.split('=')[0] == 'sink') and (line.split('=')[1][:-1] == '.true.'):
			sink_flag = True


	return [xcentre_frame, ycentre_frame, zcentre_frame],[deltax_frame, deltay_frame, deltaz_frame], \
			boxlen, proj_axis, nx, ny, max_iter, sink_flag

def load_units(opts):
	proj_list = [int(item) for item in opts.proj.split(' ')]
	proj_ind = int(proj_list[0])-1

	infof = open("%s/movie%d/info_00001.txt" % (opts.dir, proj_list[0]))
	for j, line in enumerate(infof):
		if j == 15:
			unit_l = float(line.split()[2])
		if j == 16:
			unit_d = float(line.split()[2])
		if j == 17:
			unit_t = float(line.split()[2])
		if j> 18:
			break
	unit_m = unit_d*unit_l**3/2e33 # in MSun

	return unit_l, unit_d, unit_t, unit_m


def main():

	# Parse command line arguments
	parser = OptionParser()
	parser.usage = "%prog [options] map_file"
	parser.add_option('-l','--logscale',dest='logscale', action='store', \
	    help='use log color scaling', default=False)
	parser.add_option("-m","--min",  dest="min", metavar="VALUE", \
			help='min value', default=None)
	parser.add_option("-M","--max",  dest="max", metavar="VALUE", \
			help='max value', default=None)
	parser.add_option("-f","--fmin",  dest="fmin", metavar="VALUE", \
			help='frame min value', default=0, type=int)
	parser.add_option("-F","--fmax",  dest="fmax", metavar="VALUE", \
			help='frame max value', default=-1, type=int)
	parser.add_option("-d","--dir", dest="dir", \
			help='map directory', default=os.environ['PWD'], metavar="VALUE")
	parser.add_option("-p","--proj", dest="proj", default=1, \
			help="projection index")
	parser.add_option("-s","--step", dest="step", \
			help="framing step", default=1, type=int)
	parser.add_option('-k','--kind', dest="kind", \
			help="kind of plot [temp, dens, metal]", default='dens')	
	parser.add_option('-a','--autorange',dest='autorange', action='store_true', \
	    help='use automatic dynamic range (overrides min & max)', default=False)
	parser.add_option('--clean_plot',dest='clean_plot', action='store_true', \
	    help='do not annotate plot with bar and timestamp', default=False)
	parser.add_option('--big-endian',dest='big_endian', action='store_true', \
	    help='input binary data is stored as big endian', default=False)
	parser.add_option('-c','--colormap',dest='cmap_str', metavar='CMAP', \
	    help='matplotlib color map to use', default="bone")
	parser.add_option('-b','--barlen',dest='barlen', metavar='VALUE', \
	    help='length of the bar (specify unit!)', default=5)
	parser.add_option('-B','--barlen_unit',dest='barlen_unit', metavar='VALUE', \
	    help='unit of the bar length [pc/kpc/Mpc]', default='kpc')
	parser.add_option('-t','--time_unit',dest='time_unit', metavar='VALUE', \
	    help='unit of time [Myr/Gyr]', default='Myr')
	parser.add_option('-g','--geometry',dest='geometry', metavar='VALUE', \
	    help=' [montage geometry (x,y)]', default=(1,1))
	parser.add_option("-o","--output",  dest="outfile", metavar="FILE", \
			help='output image file [default: <map_file>.png]', default=None)

	(opts,args) = parser.parse_args()

	proj_list = [int(item) for item in opts.proj.split(' ')]
	geo = [int(item) for item in opts.geometry.split(' ')]
	kind = [item for item in opts.kind.split(' ')]
	
	if opts.barlen_unit == 'pc':
		scale_l = 1e0
	elif opts.barlen_unit == 'kpc':
		scale_l = 1e3
	elif opts.barlen_unit == 'Mpc':
		scale_l = 1e6
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

	proj_ind = int(proj_list[0])-1
	sink_flag = False
	
	frame_centre, frame_delta, boxlen, proj_axis, nx, ny, max_iter, sink_flag = load_namelist_info(opts)
	unit_l, unit_d, unit_t, unit_m = load_units(opts)

	
	if (int(opts.fmax) > 0):
		max_iter=int(opts.fmax)


	# Progressbar imports/inits
	try:
		from widgets import Percentage, Bar, ETA
		from progressbar import ProgressBar
		progressbar_avail = True
	except ImportError:
		progressbar_avail = False

	if progressbar_avail:
	  widgets = ['Working...', Percentage(), Bar(marker='#'),ETA()]
	  pbar = ProgressBar(widgets=widgets, maxval = max_iter+1).start()
	else:
		print 'Working!'
        
	# Looping over movie snapshots
	for i in xrange(int(opts.fmin)+int(opts.step),max_iter+1,int(opts.step)):
		fig = plt.figure(frameon=False)
		fig.set_size_inches(nx/100*geo[1],ny/100*geo[0])
		for p in xrange(len(proj_list)):
			opts.proj = proj_list[p]
			axis = proj_axis[opts.proj]
			dat = load_map(opts,p,i)
			if sink_flag:
				plot_sinks, sink_id, sink_m, sink_pos = load_sink(opts,i)
				if plot_sinks:
					sink_m = numpy.log10(sink_m*unit_m)
					sink_pos = [x/boxlen for x in sink_pos]
		
			infof = open("%s/movie%d/info_%05d.txt" % (opts.dir, int(opts.proj), i))
			for j, line in enumerate(infof):
				if j == 8:
					time = float(line.split()[2])
				if j> 8:
					break
			
			dx_min=boxlen/2.**14
			
			if(opts.outfile==None):
				if not os.path.exists("%s/pngs/" % (opts.dir)):
							os.makedirs("%s/pngs/" % (opts.dir))
				outfile="%s/pngs/%s_%05d.png" % (opts.dir, kind[p], i/int(opts.step)-int(opts.fmin))
				if sum(geo)>2:
					outfile="%s/pngs/multi_%05d.png" % (opts.dir, i/int(opts.step)-int(opts.fmin))
			else:
				outfile=opts.outfile
			
			if(opts.logscale):
				dat = numpy.array(dat)
				if(kind[p] == 'stars'):
					dat += 1e-16
			# Reshape data to 2d
			dat = dat.reshape(ny,nx)

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
			if(opts.logscale and (kind[p][0] != 'v')):
			#if(opts.logscale and (kind[p][0] != 'v' and kind[p] != 'stars' and kind[p]  != 'dm')):
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
			
		
			# Plotting
			
			ax = fig.add_subplot(geo[0],geo[1],p+1)
			ax.set_axis_off()
			ax.axis([0,nx,0,ny])
			fig.add_axes(ax)

			if kind[p] == 'temp':
				cmap = 'jet'
			elif kind[p][0] == 'v':
				cmap = 'RdBu_r'
			else:
				cmap=opts.cmap_str
			ax.imshow(dat, interpolation = 'nearest', cmap = cmap,\
					vmin = plotmin, vmax = plotmax, aspect='auto')
			
			# Plotting sink
			#if (sink_flag and plot_sinks):
			if axis == 'x':
				w=1 # y
				h=2	# z
			elif axis == 'y':
				w=0 # x
				h=2 # z
			else:
				w=0	# x
				h=1 # y

			if (sink_flag and plot_sinks):
				ax.scatter((sink_pos[w]-frame_centre[w]/boxlen)/(frame_delta[w]/boxlen/2)*nx/2+nx/2,\
						(sink_pos[h]-frame_centre[h]/boxlen)/(frame_delta[h]/boxlen/2)*ny/2+ny/2,\
						marker='+',c='r',s=6*sink_m**2)

			if not opts.clean_plot:
				patches = []
				barlen_px = opts.barlen*scale_l*nx/(float(boxlen)*unit_l*3.24e-19*frame_delta[w]/float(boxlen))
				rect = mpatches.Rectangle((0.5*nx,0.05*ny), barlen_px, 10)
				ax.text(0.5+barlen_px/nx/2., 0.05+15./ny,"%d %s" % (opts.barlen, opts.barlen_unit),
									verticalalignment='bottom', horizontalalignment='right',
									transform=ax.transAxes,
									color='white', fontsize=14)
				patches.append(rect)

				ax.text(0.95, 0.95, '%.1f %s' % (time*unit_t/86400/365.25/scale_t, opts.time_unit),
									verticalalignment='bottom', horizontalalignment='right',
									transform=ax.transAxes,
									color='white', fontsize=14)
				collection = PatchCollection(patches, facecolor='white')
				ax.add_collection(collection)
		
		plt.axis('off') # removes axis

		# corrects window extent
		plt.subplots_adjust(left=0., bottom=0., right=1., top=1., wspace=0., hspace=0.)
		plt.savefig(outfile,dpi=100)
		plt.close(fig)
		
		if progressbar_avail:
			pbar.update(i) # updates progressbar

	if progressbar_avail:
		pbar.finish()
	else:
		print 'Finished!'
	
	# movie name for montage
	if sum(geo) > 2:
		frame = "{dir}/pngs/multi_%05d.png".format(dir=opts.dir)
		mov = "{dir}/multi.mp4".format(dir=opts.dir)
	else:
		frame = "{dir}/pngs/{kind}_%05d.png".format(dir=opts.dir,kind=opts.kind)
		mov = "{dir}/{kind}{proj}.mp4".format(dir=opts.dir, kind=opts.kind, proj=opts.proj)
	
	print 'Calling ffmpeg!'
	subprocess.call("ffmpeg -loglevel quiet -i {input} -y -vcodec h264 -pix_fmt yuv420p  -r 25 -qp 15 {output}".format(input=frame, output=mov), shell=True)
	print 'Movie created! Cleaning up!'
	if sum(geo) > 2:
		subprocess.call("rm {dir}/pngs -r".format(dir=opts.dir), shell=True)
	else:
		subprocess.call("rm {dir}/movie{proj}/pngs -r".format(dir=opts.dir, proj=opts.proj), shell=True)
	subprocess.call("chmod a+r {mov}".format(mov=mov), shell=True)

if __name__ == '__main__':
	main()
	

