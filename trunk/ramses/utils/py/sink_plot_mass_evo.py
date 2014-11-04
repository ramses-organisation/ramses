from matplotlib import pyplot as plt
from numpy import loadtxt

from argparse import ArgumentParser 

def main():
	# parsing arguments
	parser = ArgumentParser(description = "Analyze mass evolution of sink")
	parser.add_argument('--dir', dest='dir', type=str, \
			default=None, help='Directory')

	args = parser.parse_args()
	
	if args.dir == None:
		from sys import exit
		exit()
	print args.dir
	path = "/zbox/user/biernack/"  
	dir = path + args.dir + "movie1/"

	# loading sink data and times of snapshots
	data = loadtxt(dir+"all_sink.txt",delimiter=",")
	times = loadtxt(dir+"all_time.txt")

	data = data.T
	times =  times.T

	# getting units from one of the info files
	infof = open("%sinfo_00002.txt" % (dir))
	for i, line in enumerate(infof):
	  if i == 17:
	    unit_t = float(line.split()[2]) # in seconds
	  if i == 16:
	    unit_d = float(line.split()[2]) 
	  if i == 15:
	    unit_l = float(line.split()[2])
	  if i > 18:
	    break

	unit_m = unit_d*unit_l**3
	unit_t /= (86400.*365.25*1.e9) # in Gyr

	# plotting
	plt.plot(times*unit_t,data[1]*unit_m/1.988435e33)
	plt.xlabel("Time [Gyr]")
	plt.ylabel(r"Mass [$M_{\odot}$]")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	#plt.show()
	plt.savefig(path + args.dir + "M_sink_evo.png")

if __name__ == '__main__':
	main()
