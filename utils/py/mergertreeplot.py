#!/usr/bin/python2

#===============================================================================
# -----------------------------
#       TREEPLOT.PY
# -----------------------------
#
# This python script creates a plot for the merger tree of a given clump/halo
# as the root. By default, this script looks for the  clump/halo in the last 
# available output_XXXXX directory. You need to call this script from the
# directory where the output_XXXXX directories are stored.
#
# -------------
#   Usage:
# -------------
#       treeplot.py <halo-ID>
#
# -------------
#   Options:
# -------------
#       -lc, --label-clumps     if creating particle plots, label the clumps you
#                               plotted in the scatterplot
#       -pp, --plotparticles    also create a plot of particles at each output,
#                               marking the clumps in the tree with their colour
#                               ! Needs the unbinding_XXXXX.txtYYYYY  files for 
#                               this function ! 
#
#       -v, --verbose           verbose: print more details of what you're doing
#
#       -s <output_nr>          don't start at last available output_XXXXX 
#       --start_at <output_nr>  directory, but at output_<output_nr>. 
#                               <output_nr> does not need to be formatted, e.g. 5 
#                               is good enough, no need to type 00005. Naturally, 
#                               the <halo-ID> must be a halo/clump in this output 
#                               directory.
#
#       -t, --use-time          use time instead of redshift for y label.
#                               For non-cosmo runs, you NEED to use this, unless
#                               you want z=0 everywhere. (In case you forgot, 
#                               you will be warned.)
#                               if you used the -pp flag, not using the -t flag
#                               will be interpreted as plotting for a cosmo run
#                               (=> periodic boundary conditions are assumed)
#
#
# ---------------  
#  More details
# ---------------  
#
# Each drawn branch will have a color pair: an inner and an outer color. This
# is to help identifying them uniquely. If you nevertheless need/want to add
# more colors, add them in class global_params.colorlist.
# If a drawn clump doesn't have a direct progenitor, this will be signified by
# a progenitor with clump ID 0. It's supposed to signify that the clump has 
# formed in the time between the last and this output. All other nodes will be 
# labelled by the drawn clump ID at that output step.
# If a clump re-emerges more than one timestep later, there will be some clump
# ID's missing in the branch, as it skips a few output steps. Provided the
# descenant that it appears to have been merged into and later re-emerges from
# is in the tree, this will be symbolised by a dashed line. But it is possible
# that such a descendant is not in this tree, so in some cases, there may be
# gaps without the corresponding dashed lines.
#
# By default, this script creates an image called 
#   "merger_tree_halo_<last_output_dir_nr>_<halo_nr>.png"
#
# If you chose the --plot-particles flag, instead it will first make a directory 
# called
#   "merger_tree_halo_<last_output_dir_nr>_<halo_nr>/
# and save 
#   "particleplot_XXXXX.png" 
# images inside the directory.
#
# If you furthermore also chose the --label-clumps flag, instead it will create
#   "particleplot-with_labels_XXXXX.png"
# images inside the directory, so the others won't be overwritten.
# 
#
# ---------------
#  Requirements
# ---------------
# This script was tested with python 2.7.12, matplotlib 1.5.1, numpy 1.14.2, and 
# mpl_scatter_density 0.3 and fortranfile.py as it is in the ramses repository
# as of 03/2018. mpl_scatter_density and fortranfile.py are used to plot the 
# particles, but the tree plot itself doesn't require them. mpl_scatter_density
# It can be found on https://github.com/astrofrog/mpl-scatter-density and/or 
# installed via pip: $ pip install --user mpl-scatter-density
#
#===============================================================================






#===================
class _branch_x:
#=================== 
    """
    This class makes objects to store in lists to 
    dynamically allocate and track the x coordinate
    of branches, because in lists, objects are pointed to,
    not copied.
    """

    def __init__(self):
        self.x = 0
        return

    def set_x(self, x):
        self.x = x
        return





#======================
class global_params:
#======================
    """
    An object to store all global parameters in, so you can only
    pass 1 object to functions that need it.

    """


    #======================
    def __init__(self):
    #======================
        """
        Initialises object. 
        """
        self.halo = 0               # the root halo for which to plot for
        self.lastdir = ''           # last output directory
        self.noutput = 0            # number of outputs to work with
        self.ncpu = 0               # how many processors were used for simulation run

        self.workdir = ''           # current working directory
        self.prefix=''              # filename prefix
        self.outputfilename = ''    # filename for finished tree image

        self.verbose = False        # whether to print more details of what this program is doing
        self.plotparticles = False  # whether to also create plots with particles
        self.use_t = False          # whether to use time instead of redshift for y-axis labels
        self.label_clumps = False   # whether to label clumps in particle plots
        self.start = 0              # which output directory to start with

        # dictionnary of accepted keyword command line arguments
        self.accepted_args = {
            '--label-clumps' : self.set_clumplabels,
            '-lc': self.set_clumplabels,
            '--plotparticles' : self.set_plotparticles,
            '-pp': self.set_plotparticles,
            '-s' : self.set_start,
            '--start_at': self.set_start,
            '-t' : self.set_yaxis_labels,
            '--use-time' : self.set_yaxis_labels,
            '-v' : self.set_verbose,
            '--verbose' : self.set_verbose
            } 


        # Define a long list of colors for the branches
        self.colorlist=[ 
                'black',
                'red', 
                'green', 
                'gold', 
                'cornflowerblue',
                'lime', 
                'magenta', 
                'orange', 
                'mediumpurple', 
                'deeppink',
                'lightgreen',
                'saddlebrown', 
                'orchid',
                'mediumseagreen']

        return






    #=============================
    def read_cmdlineargs(self):
    #=============================
        """
        Reads in the command line arguments and stores them in the
        global_params object.
        """
        from sys import argv 

        nargs = len(argv)
        i = 1 # first cmdlinearg is filename of this file, so skip it

        while i < nargs:
            arg  = argv[i]
            arg = arg.strip()
            if arg in self.accepted_args.keys():
                if arg in ['-s', '-start_at']:
                    startnr = argv[i+1]
                    try:
                        startnr = int(startnr)
                        self.accepted_args[arg](startnr)
                    except ValueError:
                        print '"'+argv[i+1]+'"', 
                        "is not a valid start number for an output directory."
                        quit()
                    i += 1
                else:
                    self.accepted_args[arg]()
            else:
                try:
                    self.halo = int(arg)
                except ValueError:
                    print "I didn't recognize the argument '", arg, "'"
                    quit()

            i+= 1


        # defensive programming
        if self.halo <= 0:
            print "No or wrong halo given. Halo ID must be > 0"
            quit()


        return





    #==========================
    def get_output_info(self):
    #==========================
        """
        Read in the output info based on the files in the current
        working directory.
        Reads in last directory, ncpu, noutputs. 
        """

        from os import getcwd
        from os import listdir

        self.workdir = getcwd()
        filelist = listdir(self.workdir)
        
        outputlist = []
        for filename in filelist:
            if "output_" in filename:
                outputlist.append(filename)

        outputlist.sort()


        self.lastdir = outputlist[-1]
        self.lastdirnr = int(self.lastdir[-5:])
        self.noutput = len(outputlist)

        if (self.start > 0):
            startnrstr = str(self.start).zfill(5)
            dir_found = False
            for i,out in enumerate(outputlist):
                if startnrstr in out:
                    dir_found = True
                    self.lastdir = out 
                    self.lastdirnr = self.start
                    self.noutput = i+1
                    break

            if not dir_found:
                print "Didn't find specified starting directory output_"+startnrstr
                quit()
                
                    

        infofile = self.lastdir+'/'+'info_'+self.lastdir[-5:]+'.txt'
        f = open(infofile, 'r')
        ncpuline = f.readline()
        line = ncpuline.split()
        
        self.ncpu = int(line[-1])

        return


    

    #========================
    def print_params(self):
    #========================
        """
        Prints the current parameters to screen.
        Parameters:
            none
        Returns:
            nothing
        """

        print "==============================================="
        print "Working parameters are:"
        print "halo:", self.halo
        print "starting from output directory:", self.lastdir
        print "number of outputs:", self.noutput
        print "ncpus used for sim:", self.ncpu
        
        if self.use_t:
            labelname = 'time'
        else:
            labelname = 'redshift'

        print "y-axis label will be:", labelname
        print "Particles will be plotted?", self.plotparticles
        if self.plotparticles:
            print "Clumps in particle plots will be labelled?", self.label_clumps
        print "==============================================="

        return



    #========================
    # Setter methods
    #========================

    def set_clumplabels(self):
        self.label_clumps = True
        return

    def set_outputfilename(self):
        if self.plotparticles:
            self.set_prefix()

        # output filename for image
        self.outputfilename = self.prefix+"merger_tree_"+self.lastdir[-5:]+"_halo_"+str(self.halo)
        return
        

    def set_plotparticles(self):
        self.plotparticles = True
        return

    def set_prefix(self):
        import  os
        # get subdirectory name
        self.prefix = 'merger_tree_'+str(self.lastdirnr).zfill(5)+'_halo_'+str(self.halo)+'/'
        # create subdirectory if it doesn't exist
        if not os.path.exists(self.prefix):
            os.makedirs(self.prefix)
        return

    def set_start(self, startnr):
        self.start = startnr
        return

    def set_verbose(self):
        self.verbose = True
        return

    def set_yaxis_labels(self):
        self.use_t = True
        return





#==================
class Node:
#==================
    """
    A class for each node to be drawn. 
    Parameters:
        id:         (sub)halo ID
        y:          y axis value for plot (time/redshift)
        desc_id:    the ID of the "main descendant". Needed to distinguish
                    jumps across multiple snapshots.
        
    """

    #-------------------------------------
    def __init__(self, ID, y, desc_id):
    #-------------------------------------
        """
        Arguments:
            id:         (sub)halo ID [any number]
            y:          y axis value for plot (time/redshift) [any number]
            desc_id:    the ID of the "main descendant". Needed to distinguish
                        jumps across multiple snapshots. [any number]
            
        """
        
        self.id = ID                # clump ID at timestep y
        self.main_desc_id = desc_id # id of "main descendant" to check for jumpers
        self.y = y                  # value for y axis: outputnumber where it was active clump
        self.x = None               # will be set later
        self.is_jumper = False      # whether this clump will re-emerge later

        self.progs = []             # list of progenitor nodes
        self.is_main_prog = []      # wether the i-th progenitor of the progs list is the main progenitor

        self.branches = []          # list of all branches at this node
        self.branches_level = 0     # how many times this branch will split
        self.branches_tot = 0       # total number of branches for this branch
                                    # Needed for nice plotting
        self.start_of_branch = None # Node at which the branch this clump is in starts
        self.walked = False         # whether this node has been walked over already
        

        return



    #-----------------------------------------------
    def add_progenitor(self, prog_node, is_main):
    #-----------------------------------------------
        """
        Adds a progenitor to this node's list.
        Arguments:
            prog_node:  node of the progenitor to add to list [class Node object]
            is_main:    whether this progenitor is the main progenitor of this node [boolean] 
        """
        self.progs.append(prog_node)
        self.is_main_prog.append(is_main)
        return





#=========================================
def _draw_tree(node, ax, colorlist):
#=========================================
    """
    Plots all connections of Node node to its progenitors (by calling
    _draw_tree_on_plot), then calls itself recursively.

    Arguments:
        node:       class Node object whose progenitors are to be plotted
        ax:         axis object of the plot
        colorlist:  list of colors for different branches

    returns:
        nothing
    """


    node.walked = True
    
    for i, prog in enumerate(node.progs):
        # call actual drawing function
        _draw_tree_on_plot(node, prog, ax, colorlist)
        # call yourself recursively for the main progenitor
        if (not prog.walked) or prog.is_jumper:
            _draw_tree(prog, ax, colorlist)


    # if you reached the leaves, just add the progenitor numbers
    if len(node.progs) == 0:
        _draw_tree_on_plot(node, node, ax, colorlist)



    return

    



#====================================================
def _draw_tree_on_plot(node, prog, ax, colorlist):
#====================================================
    """
    The actual drawing function of the tree. Draws a line between Node node 
    and its progenitor prog.
    If the progenitor re-imerges at a later timestep, draw a dotted line 
    instead of a solid line to signify where it looks like it  merged into.

    Arguments:
        node:       class Node object of descendant
        prog:       class Node object of progenitor
        ax:         axis object of the plot
        colorlist:  list of colors for different branches

    returns:
        nothing
    """

    # get x and y values for plot
    x = [node.x.x, prog.x.x]
    y = [node.y, prog.y]


    # get line colors
    # color index is same as x!
    inner, outer = _get_plotcolors(prog.x.x, colorlist)
    myfacecolor = colorlist[inner]
    myedgecolor = colorlist[outer]

    # Determine line-style
    # Dashed for mergers that will re-emerge later
    linestyle='-'
    linewidth = 4
    outerlinewidth = 8
    alp = 1
    zorder = 1

    if node.id != prog.id: 
        # if not just adding label
        if  node.id != prog.main_desc_id: 
            linestyle = '--'
            linewidth = 3
            outerlinewidth = 4
            alp = 0.2
            zorder = 0



    #---------------
    # Plot the line
    #---------------

    # plot outer line
    ax.plot(x, y,
            color = myedgecolor,
            lw=outerlinewidth,
            ls='-',
            alpha=alp,
            zorder=zorder)

    # plot inner line
    ax.plot(x, y, 
            color = myfacecolor, 
            lw=linewidth,
            ls=linestyle,
            zorder=zorder)



    #---------------
    # Annotation
    #---------------
 
    # take here node data instead of prog!
    inner, outer = _get_plotcolors(node.x.x, colorlist)
    myfacecolor = colorlist[inner]
    myedgecolor = colorlist[outer]

    # Annotate the dots with the clump ID

    bbox_props = dict(boxstyle="round,pad=0.15", 
            fc=myfacecolor, 
            ec=myedgecolor,
            lw=1.5, 
            alpha=1)

    ax.text(node.x.x, node.y, str(node.id),
            size=6,
            bbox=bbox_props,
            horizontalalignment = 'center',
            verticalalignment = 'center',
            rotation=45,
            zorder=2)




    return





#==============================================
def _get_plotcolors(color_index, colorlist):
#==============================================
    """
    Gets inner and outer plot line / scatterpoint colors.
    NEVER returns zero, zeroth index is black, reserved for special cases.
    (= for particles not in the tree)

    Parameters:
        color_index:    index for this color-pair
        colorlist:      list of colors to choose from

    Returns:
        inner, outer:   indices for inner / outer colors from colorlist
    """

    shortlen = len(colorlist) - 1
    multiples = color_index//shortlen
    inner = color_index - (multiples) * shortlen + 1
    while multiples >= shortlen:
        multiples = multiples//shortlen
    outer = -(multiples+1)
    if outer == inner:
        outer += 1

    if inner==0:
        print "got inner 0"

    return inner, outer





#=============================================================
def _get_x(node, x_values, straight):
#=============================================================
    """
    Assign x values for the plot to a node and its progenitors.

    Function summary:
    if straight = True:
        recursively descend down the line of main
        progenitors, then return to where you left off.
        Main progenitors inherit their descendants x object.
    else:
        if this node's progenitor hasn't been walked over already,
        assign a new x-value to it. Then call this function 
        recursively for all appropriate progenitors of this node 
        with straight=True
    For all appropriate progenitors:
        call this function recursively for every appropriate
        progenitor of this node with straight=False 

    Arguments:
        node:       class Node object to check for
        x_values:   list of _branch_x objects 
        straight:   whether to go straight first

    returns:
        x_values:   list of _branch_x objects 

    """
 

    # Loop over progenitors
    for i, prog in enumerate(node.progs):
        
        if straight:
            # Only go down main progenitors
            if not prog.walked and node.is_main_prog[i]:
                # Prog inherits descendant's x
                prog.x = node.x
                prog.walked = True
                if prog.id > 0:
                    # call yourself recursively
                    x_values = _get_x(prog, x_values, straight=True)
                # when you're done, return! Need to go down branches 
                # from top to bottom!
                return x_values


        else:
            if (not node.is_main_prog[i]) and (not prog.is_jumper):

                # if the progenitor doesn't have an x already, give it one
                if not prog.walked:
                    ids = (branch.id for branch in prog.start_of_branch.branches)
                    ids = list(ids)

                    bind = ids.index(prog.id)
                    new_x = _branch_x()
                    
                    # figure out which way to go
                    x = x_values.index(node.x)
                    if bind%2==0: # go right
                        x_values = x_values[:x+1]+[new_x]+x_values[x+1:]
                        prog.x = x_values[x+1]
                    else: # go left
                        x_values = x_values[:x]+[new_x]+x_values[x:]
                        prog.x = x_values[x]

                    prog.walked = True

                # descend down branch along main progenitors only 
                if prog.id > 0:
                    x_values = _get_x(prog, x_values, straight=True)


    # When you're done going down straight lines, now check for progenitor's progenitors
    # IMPORTANT: DO THIS ONLY AFTER PREVIOUS FOR LOOP IS DONE
    for i,prog in enumerate(node.progs):
        descend = (not prog.is_jumper) or node.is_main_prog[i]
        if prog.id > 0 and descend:
            x_values = _get_x(prog, x_values, straight=False)



    return x_values





#=======================================================================
def make_tree(progenitors, descendants, progenitor_outputnrs, params):
#=======================================================================
    """
    makes a tree out of read in lists. 
    The tree is a list containing lists of class Node objects
    for each output step.

    parameters:
        progenitors:            list of lists of progenitors for each timestep
        descendants:            list of lists of descendants for each timestep
        progenitor_outputnrs:   list of lists of the output numbers when 
                                any given progenitor was an active clump
        params:                 class global_params object

    returns:
        tree:                   list of lists of nodes constituting the tree
    """

    
    #------------------
    # Setup
    #------------------

    halo = params.halo
    lastdirnr = params.lastdirnr

    nout_present = len(progenitors) 
    output_start = lastdirnr - nout_present

    if params.verbose:
        print "Creating tree."
        print 



    #---------------------
    # initialise tree
    #---------------------

    tree = []

    # create empty list for each output
    for i in range(nout_present+1):
        tree.append([])

    # enter root
    rootnode = Node(halo, lastdirnr, 0)
    tree[0]=[rootnode]

    

    #---------------------
    # Make tree
    #---------------------

    print 
    print "----------------------------------------------------------------------"

    # Loop over all snapshots
    for out in range(nout_present):
        found_one = False

        # for each branch of the tree at that snapshots:
        for branch in tree[out]:

            # find for which progenitor the descendant matches
            for i in range(len(progenitors[out])) :

                snapnr = progenitor_outputnrs[out][i]       # snapshot nr for progenitor
                ind = nout_present + output_start - snapnr  # index in tree / tree level
                progid = abs(progenitors[out][i])           # progenitor ID

                if abs(descendants[out][i]) == branch.id:
                    found_one = True
                    is_main = descendants[out][i] == branch.id
                    # is a main progenitor if desc ID == branch ID
                    # is a merger if desc ID == -branch ID

                    # Check first if progenitor is already in list
                    prog_not_in_list = True

                    if (progid > 0):    # always add a new 0!
                        if (len(tree[ind]) > 0):
                            for j, candidate in enumerate(tree[ind]):
                                if (candidate.id == progid):
                                    # Then this progenitor is already in the tree.
                                    branch.add_progenitor(tree[ind][j], is_main)
                                    prog_not_in_list = False
                                    break

                    if prog_not_in_list:
                        # create new node for progenitor
                        newnode = Node(progid, snapnr, branch.id)
                        if branch.y-snapnr > 1 :
                            newnode.is_jumper = True

                        # add new node to tree
                        tree[ind].append(newnode)

                        # link to its progenitor:
                        # you know in which outputnr it will be
                        # since it's the last added element, it's index will
                        # be len(tree at that outputnr) - 1
                        branch.add_progenitor(tree[ind][-1], is_main)


                    # Print informations
                    output_string = 'Adding progenitor '
                    output_string += '{0:7d}'.format(progid)+' '
                    output_string += 'for descendant ' 
                    output_string += '{0:7d}'.format(branch.id)+' '
                    output_string += "|| snapshot "  
                    output_string += '{0:3d}'.format(snapnr)+' '
                    output_string += "->"  
                    output_string += '{0:3d}'.format(branch.y)

                    if (is_main):
                        if branch.y-snapnr > 1 :
                            output_string += "   found jumper"
                        print output_string 

                    else:
                        if branch.y-snapnr > 1 :
                            output_string += "   found jumper and merger"
                        else:
                            output_string += "   found merger"
                        print output_string

        # make new line after every snapshot read-in
        if found_one:
            print "----------------------------------------------------------------------"

    print "\n\n"


    # remove empty lists at the end of the tree:
    # clumps might not have been around since the first output

    list_is_empty = (len(tree[-1]) == 0)

    while list_is_empty:
        del tree[-1]
        list_is_empty = len(tree[-1]) == 0

            


    return tree





#============================================================================
def _plot_particles(x, y, z, clumpid, time, clumps_in_tree, colors, params):
#============================================================================
    """
    This function plots the particles for this output. Is meant
    to be called for every output by the plot_treeparticles function.

    Parameters:
        x, y, z:        numpy arrays of particle positions
        clumpid:        numpy arrays of particle clump IDs [list of integers]
        time:           time or redshift of output
        clumps_in_tree: list of clumps in tree at this output
        colors:         the colors of clumps in trees [list of integers]

    Returns:
        nothing
    """

    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    import numpy as np
    import matplotlib.colors as mc
    import mpl_scatter_density


    if params.verbose:
        print "Creating particle figure"



    #---------------------------
    # Set up figure
    #---------------------------

    fig = plt.figure(facecolor='white', figsize=(15,6), dpi=150)

    # add subplots on fixed positions
    gs = gridspec.GridSpec(1, 3,
                       width_ratios=[1, 1, 1], height_ratios=[1],
                       left=0.04, bottom=0.13, right=0.98, top=0.9, wspace=0.12,)
    ax1 = fig.add_subplot(gs[0], aspect='equal', projection='scatter_density')
    ax2 = fig.add_subplot(gs[1], aspect='equal', projection='scatter_density')
    ax3 = fig.add_subplot(gs[2], aspect='equal', projection='scatter_density')



    #------------------------------
    # Setup plot region
    #------------------------------

    to_plot_all = np.array(x.shape, dtype='bool')
    to_plot_all = False
    colorarray = np.zeros(x.shape, dtype='int')

    for i,clump in enumerate(clumps_in_tree):
        is_this_clump = (clumpid == clump)
        to_plot_all = to_plot_all | is_this_clump
        fc, ec = _get_plotcolors(colors[i], params.colorlist)
        colorarray[is_this_clump] = fc
   
    borders = _setup_plot_region(x, y, z, to_plot_all, params)
    xmin, xmax, ymin, ymax, zmin, zmax = borders



    # reset to_plot here!
    # plot only particles within the limits
    to_plot_all = ((x <= xmax) & (x >= xmin))
    to_plot_all = to_plot_all & ((y <= ymax) & (y >= ymin))
    to_plot_all = to_plot_all & ((z <= zmax) & (z >= zmin))





    #---------------------------
    # Actual plotting
    #---------------------------

    # set up colormap
    bounds=np.linspace(0, len(params.colorlist), len(params.colorlist)+1)
    cmap=mc.ListedColormap(params.colorlist, name='My colormap')
    norm=mc.BoundaryNorm(bounds, len(params.colorlist))



    # first plot particles not in tree
    to_plot = to_plot_all & (clumpid==0)

    try:
        ax1.scatter_density(x[to_plot], y[to_plot], c=colorarray[to_plot], 
            cmap=cmap, norm=norm, dpi=150, alpha=1, zorder=0)
        ax2.scatter_density(y[to_plot], z[to_plot], c=colorarray[to_plot], 
            cmap=cmap, norm=norm, dpi=150, alpha=1, zorder=0)
        ax3.scatter_density(x[to_plot], z[to_plot], c=colorarray[to_plot], 
            cmap=cmap, norm=norm, dpi=150, alpha=1, zorder=0)
    except ValueError:
        # in case there were no particles to plot:
        print "No particles outside of tree were found to be plotted."
    


    # now plot particles in tree
    for i, c in enumerate(clumps_in_tree):
        try:
            to_plot = to_plot_all & (clumpid==c)
            ax1.scatter_density(x[to_plot], y[to_plot], c=colorarray[to_plot], 
                cmap=cmap, norm=norm, dpi=72, alpha=1, zorder=1)
            ax2.scatter_density(y[to_plot], z[to_plot], c=colorarray[to_plot], 
                cmap=cmap, norm=norm, dpi=72, alpha=1, zorder=1)
            ax3.scatter_density(x[to_plot], z[to_plot], c=colorarray[to_plot], 
                cmap=cmap, norm=norm, dpi=72, alpha=1, zorder=1)


            if params.label_clumps: 
                #------------------------
                # create clump labels
                #------------------------
                inner, outer = _get_plotcolors(colors[i], params.colorlist)
                fc = params.colorlist[inner]
                ec = params.colorlist[outer]

                bbox_props = dict(boxstyle="round,pad=0.15", 
                        fc=fc, 
                        ec=ec,
                        lw=1.5, 
                        alpha=1)

                xc = np.mean(x[to_plot])
                yc = np.mean(y[to_plot])
                zc = np.mean(z[to_plot])

                ax1.text(xc, yc, str(c),
                    size=6,
                    bbox=bbox_props,
                    horizontalalignment='center',
                    verticalalignment='center',
                    )
                ax2.text(yc, zc, str(c),
                    size=6,
                    bbox=bbox_props,
                    horizontalalignment='center',
                    verticalalignment='center',
                    )
                ax3.text(xc, zc, str(c),
                    size=6,
                    bbox=bbox_props,
                    horizontalalignment='center',
                    verticalalignment='center',
                    )

        except ValueError:
            #in case there were no particles to plot:
            print "WARNING: ValueError while plotting clump", c, "raised."
            print "WARNING: Maybe clump has no particles for some reason?"
            print "WARNING: Or are you playing around with the code?"
            acceptable_answer = False
            while not acceptable_answer:
                answer = raw_input("Do you still want to continue? (y/n) ")
                if (answer == 'y'):
                    acceptable_answer = True
                    break
                elif (answer == 'n'):
                    print "Exiting."
                    quit()
                else:
                    print "Please answer with 'y' or 'n'"







    #-------------------------------
    # Tweak and save figure
    #-------------------------------
    
    _tweak_particleplot(fig, time, borders, params)


    _save_fig(fig, params)

    plt.close()
    
    return






#=======================================================
def plot_tree(tree, yaxis_int, yaxis_phys, params):
#=======================================================
    """
    The main function for plotting the tree.
    After some preparation, it distributes x-values to branches
    and then plots them.

    Arguments:
        tree:       list of lists of class Node objects, constituting the tree
        yaxis_int:  array for y axis ticks containing integers
                    (= output numbers)
        yaxis_phys: array for y axis ticks containing physical data
                    (= redshift or time)
        params:     class global_params object

    returns:
        nothing
    """

    import matplotlib.pyplot as plt


    #-----------------
    # Preparation
    #-----------------

    if params.verbose:
        print "Preparing to plot tree."

    # First find the number of branches for each branch.
    tree[0][0].start_of_branch = tree[0][0]
    _walk_tree(tree[0][0], tree[0][0])

    # Now recursively sum up the number of branches
    # to know how much space to leave between them for plotting
    if params.verbose:
        print "Sorting branches."
    _sum_branches(tree[0][0])

    # lastly, sort the root's branches
    _sort_branch(tree[0][0]) 

    # reset whether nodes have been walked for _get_x()
    for level in tree:
        for branch in level:
            branch.walked = False




    #--------------------------------
    # start distributing x values      
    #--------------------------------

    if params.verbose:
        print "Assigning positions to nodes and branches."

    # give initial values for root and borders for axes
    x_values = list([_branch_x()])
    tree[0][0].x = x_values[0]

    x_values = _get_x(tree[0][0], x_values, straight=True)
    x_values = _get_x(tree[0][0], x_values, straight=False)



    # once you're done, assign actual integer values for nodes
    for i,x in enumerate(x_values):
        x.set_x(i) 

    # get x-axis borders for plot
    borders = [-len(x_values)/20, len(x_values)*(1+1/20)]





    #---------------------
    # Plot the tree
    #---------------------

    if params.verbose:
        print "Plotting tree with", tree[0][0].branches_tot+1, \
                "branches in total over", tree[0][0].y - tree[-1][0].y + 1, \
                "output steps."
        print "Creating figure."

    # create figure
    fig = plt.figure(facecolor='white')
    ax = fig.add_subplot(111)


    # reset whether nodes have been walked for _draw_tree()
    for level in tree:
        for branch in level:
            branch.walked = False

    # draw the tree
    _draw_tree(tree[0][0], ax, params.colorlist)

    # Tweak the plot
    _tweak_treeplot(fig, yaxis_int[:len(tree)+1], yaxis_phys[:len(tree)+1], borders, params)
    
    # Save the figure
    _save_fig(fig, params)




    return 





#=================================================
def plot_treeparticles(tree, yaxis_phys, params):
#=================================================
    """
    The main function to plot the particles currently in trees.
    Calls _plot_particles for every output in the tree.
    NOTE: This function assumes that the colors for the tree have
    already been assigned and takes the same color of the node
    for the particles / clump labels..

    Parameters:
        tree:       list of list of class Node objects containing the tree
        yaxis_phys: array for y axis ticks containing physical data
                    (= redshift or time)
        params:     class global_params object, containing global parameters

    Returns:
        nothing
    """

    if params.verbose:
        print "Started plotting tree particles."


    # for every output in the tree:
    for i,out in enumerate(tree):

        # find which clumps are in the tree and their colors
        clumps_in_tree = []
        colors = []
        for branch in out:
            if branch.id != 0:
                clumps_in_tree.append(branch.id)
                colors.append(branch.x.x)

        # if there is work to be done
        if len(clumps_in_tree) > 0:

            outnr = params.lastdirnr - i
            srcdir = 'output_'+str(outnr).zfill(5)

            # read in particles of this output
            x, y, z, clumpid = _read_particle_data(srcdir, params)

            # get time/redshift
            time = yaxis_phys[i]

            # reset outputfilename
            if params.label_clumps:
                params.outputfilename = params.prefix+'particleplot-with-labels_'+str(outnr).zfill(5)
            else:
                params.outputfilename = params.prefix+'particleplot_'+str(outnr).zfill(5)

            # plot the particles
            _plot_particles(x, y, z, clumpid, time, clumps_in_tree, colors, params)

            del x, y, z, clumpid, clumps_in_tree, colors
   

    return





#===================================
def read_mergertree_data(params):
#===================================
    """
    reads in mergertree data as written by the mergertree patch.
    NOTE: requires mergertree.txtYYYYY files.

    parameters:
        params:     class global_params object, containing the global parameters

    returns:
        progenitors :           lists of lists of progenitors and descendants,
        descendants :           starting with the last output step.
        progenitor_outputnrs:   the output number at which the progenitor is
        outputnrs   :           the output number at which descendants were taken from
        time :                  list of times correspondig to each output step

    """ 

    import numpy as np
    import warnings
    import gc

    noutput = params.noutput
    ncpu = params.ncpu
    use_t = params.use_t

    if params.verbose:
        print "Reading in mergertree data."

    


    # create lists where to store stuff
    fname = 'mergertree.txt'

    progenitors = []
    descendants = []
    time = []
    progenitor_outputnrs = []
    outputnrs = []

    startnr=params.lastdirnr
    dir_template = params.lastdir[:-5] # = 'output_'


    #---------------------------
    # Loop over directories
    #---------------------------

    for output in range(noutput):
        # loop through every output: Progenitor data only starts at output_00002,
        # but you'll need time/redshift data from output_00001!

        # Start with last directory (e.g. output_00060),
        # work your way to first directory (e.g. output_00001)
        dirnr =  str(startnr - output).zfill(5)
        srcdir = dir_template + dirnr

        if output < noutput-1: # don't try to read progenitor stuff from output_00001
            #------------------------------
            # Read in progenitor data 
            #------------------------------
            progs_snapshot = []
            descs_snapshot = []
            prog_outputnr_snapshot = []

            # Stop early if you reach a directory that has no mergertree.txt* files
            # (Can happen if there are no halos in the simulation yet)
            try:
                # loop over files
                for cpu in range(ncpu):
                    filenr = str(cpu + 1).zfill(5)          # get 00001, 00002, ...
                    fileloc = srcdir + '/' + fname + filenr


                    # ignore "empty files" warnings
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        raw_data = np.loadtxt(fileloc, dtype='int', skiprows=1, usecols=([0, 1, 2]))

                    #intercept index error if there is only one halo
                    try:
                        dshape = raw_data.shape[1] #try if the shape is 2dimensional

                        desc = raw_data[:,0].tolist()
                        prog = raw_data[:,1].tolist()
                        prog_outputnr = raw_data[:,2].tolist()

                    except IndexError: 

                        # in this case, the read file has only 1 halo, or none.
                        # Check if there is 1:

                        if (len(raw_data)>0):
                            desc = [raw_data[0]]
                            prog = [raw_data[1]]
                            prog_outputnr = [raw_data[2]]
                        else:
                            continue

                    # append read list for each CPU to list for entire snapshot
                    for i in range(len(prog)):
                        progs_snapshot.append(prog[i])
                        descs_snapshot.append(desc[i])
                        prog_outputnr_snapshot.append(prog_outputnr[i])
                

                # append entire list here!
                descendants.append(descs_snapshot)
                progenitors.append(progs_snapshot)
                progenitor_outputnrs.append(prog_outputnr_snapshot)
                outputnrs.append(startnr - output)

            except IOError: # If file doesn't exist
                print "Didn't find any progenitor data in", srcdir



        
        try:
            #-------------------------------------
            # get time, even for output_00001
            #-------------------------------------
            fileloc = srcdir+'/info_'+dirnr+'.txt'
            infofile = open(fileloc)
            for i in range(8):
                infofile.readline() # skip first 8 lines
            
            if not use_t:
                infofile.readline() # skip another line for redshift

            timeline = infofile.readline()
            timestring, equal, timeval = timeline.partition("=")
            timefloat = float(timeval)

            if not use_t:
                timefloat = 1.0/timefloat - 1
                
            time.append(timefloat)
    
        except IOError: # If file doesn't exist
            print "Didn't find any info data in ", srcdir
            break



    
    #----------------------------------------------------------
    # print warning if -t cmdline arg might've been forgotten
    #----------------------------------------------------------

    if (len(time)>1 and (not use_t)):
        if (time[0] == time[1]):
            print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            print "WARNING: The first two elements in the physical y-axis data are the same."
            print "WARNING: If you are trying to plot a non-cosmo simulation, you should use"
            print "WARNING: the -t flag. Otherwise, you'll always have z = 0 on the y-axis."
            print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            acceptable_answer = False
            while not acceptable_answer:
                answer = raw_input("Do you still want to continue? (y/n) ")
                if (answer == 'y'):
                    acceptable_answer = True
                    break
                elif (answer == 'n'):
                    print "Exiting."
                    quit()
                else:
                    print "Please answer with 'y' or 'n'"

    # collect garbage
    gc.collect()

    return descendants, progenitors, progenitor_outputnrs, outputnrs, time
   




#========================================
def _read_particle_data(srcdir, params):
#========================================
    """
    Reads in the particle data from directory srcdir.
    NOTE: requires part_XXXXX.outYYYYY and unbinding_XXXXX.outYYYYY files

    Parameters:
        srcdir:     String of directory where to read data from
        params:     class global_params object

    returns:
        x,y,z:      numpy arrays of particle positions
        clumpid:    numpy arrays of particle clump IDs
    """

    import numpy as np
    from os import listdir
    import fortranfile as ff

    if params.verbose:
        print "Reading in particles of output", int(srcdir[-5:])

    srcdirlist = listdir(srcdir)

    if 'unbinding_'+srcdir[-5:]+'.out00001' not in srcdirlist:
        print "Couldn't find unbinding_"+srcdir[-5:]+".out00001 in", srcdir
        print "To plot particles, I require the unbinding output."
        quit()






    #-----------------------
    # First read headers
    #-----------------------
    nparts = np.zeros(params.ncpu, dtype='int')
    partfiles = [0]*params.ncpu 

    for cpu in range(params.ncpu):
        srcfile = srcdir+'/part_'+srcdir[-5:]+'.out'+str(cpu+1).zfill(5)
        partfiles[cpu] = ff.FortranFile(srcfile)

        ncpu = partfiles[cpu].readInts()
        ndim = partfiles[cpu].readInts()
        nparts[cpu] = partfiles[cpu].readInts()
        localseed = partfiles[cpu].readInts()
        nstar_tot = partfiles[cpu].readInts()
        mstar_tot = partfiles[cpu].readReals('d')
        mstar_lost = partfiles[cpu].readReals('d')
        nsink = partfiles[cpu].readInts()

        del ncpu, ndim, localseed, nstar_tot, mstar_tot, mstar_lost, nsink


    
    #-------------------
    # Allocate arrays
    #-------------------
    nparttot = nparts.sum()
    x = np.zeros(nparttot, dtype='float')
    y = np.zeros(nparttot, dtype='float')
    z = np.zeros(nparttot, dtype='float')
    clumpid = np.zeros(nparttot, dtype='int')


    #----------------------
    # Read particle data
    #----------------------

    start_ind = np.zeros(params.ncpu, dtype='int')
    for cpu in range(params.ncpu-1):
        start_ind[cpu+1] = nparts[cpu] + start_ind[cpu]

    for cpu in range(params.ncpu):
        x[start_ind[cpu]:start_ind[cpu]+nparts[cpu]] = partfiles[cpu].readReals('d')
        y[start_ind[cpu]:start_ind[cpu]+nparts[cpu]] = partfiles[cpu].readReals('d')
        z[start_ind[cpu]:start_ind[cpu]+nparts[cpu]] = partfiles[cpu].readReals('d')
        
        srcfile = srcdir+'/part_'+srcdir[-5:]+'.out'+str(cpu+1).zfill(5)
        unbfile = srcdir+'/unbinding_'+srcdir[-5:]+'.out'+str(cpu+1).zfill(5)
        unbffile = ff.FortranFile(unbfile)

        clumpid[start_ind[cpu]:start_ind[cpu]+nparts[cpu]] = unbffile.readInts()


    clumpid = np.absolute(clumpid) 


    return x, y, z, clumpid





#============================================
def _save_fig(fig, params):
#============================================
    """
    Save figure as png.
    this_name:  name to save figure with
    fig:        figure object
    params:     class global_params object

    returns:
        nothing

    """
    
    import matplotlib.pyplot as plt

    fig_path = params.workdir+'/'+params.outputfilename+'.png'
    if params.verbose:
        print "saving figure as "+fig_path

    plt.savefig(fig_path, format='png', facecolor=fig.get_facecolor(), transparent=False, dpi=150)
    print "\nSaved", fig_path, "\n"


    plt.close()
   
    return




#======================================================
def _setup_plot_region(x, y, z, to_plot, params):
#======================================================
    """
    This function sets up the plot region: 
        - find borders for the plot: min/max in x,y,z direction
        - account for periodic scenarios
    Parameters:
        x, y, z :   numpy arrays of x, y, z positions of particles
        to_plot :   numpy boolean array of particles which are to be plotted
        params :    class global_params object

    Returns:
        borders :   tuple of borders:  (xmin, xmax, ymin, ymax, zmin, zmax)
    """

      
    xmin = x[to_plot].min()
    xmax = x[to_plot].max()
    ymin = y[to_plot].min()
    ymax = y[to_plot].max()
    zmin = z[to_plot].min()
    zmax = z[to_plot].max()

    xc = (xmax + xmin)/2
    dx = xmax - xmin
    yc = (ymax + ymin)/2
    dy = ymax - ymin
    zc = (zmax + zmin)/2
    dz = zmax - zmin


    if not params.use_t: 
        # if not use_t, then assume cosmo run
        # find out whether you need to shift for periodicity
        #----------------------------
        # Check for periodicity
        #----------------------------
        while dx > 0.5:
            if (xc >= 0.5) :
                x[x<0.5] += 1
            else:
                x[x>=0.5] -= 1
            xmin = x[to_plot].min()
            xmax = x[to_plot].max()
            xc = (xmax + xmin)/2
            dx = xmax - xmin

        while dy > 0.5:
            if (yc >= 0.5) :
                y[y<0.5] += 1
            else:
                y[y>=0.5] -= 1
            ymin = y[to_plot].min()
            ymax = y[to_plot].max()
            yc = (ymax + ymin)/2
            dy = ymax - ymin

        while dz > 0.5:
            if (zc >= 0.5) :
                z[z<0.5] += 1
            else:
                z[z>=0.5] -= 1
            zmin = z[to_plot].min()
            zmax = z[to_plot].max()
            zc = (zmax + zmin)/2
            dz = zmax - zmin


    maxd = max(dx, dy, dz)
    xmin = xc - maxd*0.55
    xmax = xc + maxd*0.55
    ymin = yc - maxd*0.55
    ymax = yc + maxd*0.55
    zmin = zc - maxd*0.55
    zmax = zc + maxd*0.55

    
    borders = (xmin, xmax, ymin, ymax, zmin, zmax)



    return borders





#===========================
def _sort_branch(branch):
#===========================
    """
    Sort the list of branches of a given start of a branch (class Node object)
    by increasing time of appearance. The earliest branches (lowest on y-axis)
    come first. If there are multiple mergers at the same time, put the one 
    with more own branches further out.

    Parameters:
        branch:     class Node objecs whose .branches list is to be sorted

    returns:
        nothing
    """
            
    import numpy as np
    import gc

    if len(branch.branches) > 0:
        branchlist = (bbranch for bbranch in branch.branches)
        branchlist = list(branchlist)

        nbranch = len(branchlist)
        times = np.zeros(nbranch, dtype='int')
        all_branches = np.zeros(nbranch, dtype='int')
        branchlist_id = np.zeros(nbranch, dtype='int')

        for b in range(nbranch):
            times[b]=branchlist[b].y
            all_branches[b]=branchlist[b].branches_tot
            branchlist_id[b] = branchlist[b].id

        # find out if you have multiple mergers at same timestep
        needs_branchcount_sorting = False

        for t in times:
            temp = times[times==t]
            if ((temp/temp).sum()>1):
                needs_branchcount_sorting = True
                break

        sort_ind = times.argsort()
        times = times[sort_ind]
        branchlist_id = branchlist_id[sort_ind]


        if needs_branchcount_sorting:

            all_branches = all_branches[sort_ind]
            
            # at least the next element is the same
            i = 0
            while i < (len(times)-1):
                if (times[i+1]==times[i]):
                    # if the following is the same:
                    startind = i
                    endind = i+1

                    # find where the same end
                    while times[startind]==times[endind]:
                        if endind == len(times)-1:
                            break
                        else:
                            endind += 1

                    # sort all of them by ascending branch_tot
                    branchtots = all_branches[startind:endind+1]
                    sort_ind2 = branchtots.argsort()

                    branchlist_id[startind:endind+1] = branchlist_id[startind+sort_ind2]


                    i = endind

                else:
                    i+=1


        # overwrite branches' branch list
        # work with branch.ids to not copy entire objects all the time
        for i in range(len(times)):
            for bbranch in branchlist:
                if bbranch.id == branchlist_id[i]:
                    branch.branches[i] = bbranch
                    break

        # release memory
        del branchlist
        del branchlist_id
        del times
        del all_branches
        del sort_ind

    gc.collect()


    return





#============================
def _sum_branches(node):
#============================
    """
    Recursively sum up the total number of branches of each "root"

    Arguments:
        node:   class Node object to check for

    returns:
        nothing
    """


    # each branch root has itself as a branch.
    # If the tree splits up somewhere else along the way,
    # then there must be > 1 branch in the root node.
    
    if len(node.branches) > 0:
        node.branches_level = 1
        for i, branch in enumerate(node.branches):
            # don't call yourself again
            _sum_branches(branch)
            node.branches_tot += branch.branches_tot
            node.branches_tot += 1 # count yourself too :)
            
            node.branches_level += branch.branches_level
            _sort_branch(branch)

    return





#==========================================================================
def _tweak_particleplot(fig, time, borders, params):
#==========================================================================
    """
    Tweak the particle plot. Set ticks, background color, 
    label axes, get title, and some other cosmetics.

    Parameters:
        fig:            figure object
        time:           time or redshift of output
        borders:        tuple of plot edges: xmin, xmax, ..., zmin, zmax
        params:         class global_params object

    Returns:
        nothing
    """

    ax1, ax2, ax3 = fig.axes


    #------------------------------------------
    # set tick params (especially digit size)
    #------------------------------------------

    ax1.tick_params(axis='both', which='major', labelsize=8,top=5)
    ax2.tick_params(axis='both', which='major', labelsize=8,top=5)
    ax3.tick_params(axis='both', which='major', labelsize=8,top=5)
    



    #--------------
    # label axes
    #--------------

    ax1.set_xlabel(r'x', labelpad=4, family='serif',size=12)
    ax1.set_ylabel(r'y', labelpad=4, family='serif',size=12)

    ax2.set_xlabel(r'y', labelpad=4, family='serif',size=12)
    ax2.set_ylabel(r'z', labelpad=4, family='serif',size=12)

    ax3.set_xlabel(r'x', labelpad=4, family='serif',size=12)
    ax3.set_ylabel(r'z', labelpad=4, family='serif',size=12)




    #--------------
    # Add title
    #--------------

    title = "Particles in tree at output "+str(int(params.outputfilename[-5:]))+"; "

    if params.use_t:
        title += "t = "
    else:
        title += "z = "

   
    # pad the time/redshift with zeros at the end
    num = str(round(time, 3))
    if (round(time, 3)) < 0:
        # one longer for negative z
        end = 7
    else:
        end = 6
    while len(num) < end:
        num += '0'
    title += num

    fig.suptitle(title, family='serif', size=18)




    #-------------------
    # set axes limits
    #-------------------

    xmin, xmax, ymin, ymax, zmin, zmax = borders

    ax1.set_xlim([xmin, xmax])
    ax1.set_ylim([ymin, ymax])

    ax2.set_xlim([ymin, ymax])
    ax2.set_ylim([zmin, zmax])

    ax3.set_xlim([xmin, xmax])
    ax3.set_ylim([zmin, zmax])





    #------------------------
    # Set background color
    #------------------------
    ax1.set_axis_bgcolor('aliceblue')
    ax2.set_axis_bgcolor('aliceblue')
    ax3.set_axis_bgcolor('aliceblue')



    return





#====================================================================
def _tweak_treeplot(fig, yaxis_int, yaxis_phys, borders, params):
#====================================================================
    """
    tweaks the plot. Removes x-axis ticks, labels y-axis
    ticks nicely, adds right y axis.

    parameters:
        fig:        pyplot figure object
        yaxis_int:  array for y axis ticks containing integers
                    (= output numbers)
        yaxis_phys: array for y axis ticks containing physical data
                    (= redshift or time)
        borders:    borders for x axis
        params:     class global_params object

    returns:
        nothing
    """

    import matplotlib.pyplot as plt
    from copy import deepcopy

    #-----------------
    # preparation
    #-----------------

    if params.verbose:
        print "Tweaking the plot."

    halo = params.halo 
    lastdirnr = params.lastdirnr



    # add additional yaxis points for padding
    outputnr = deepcopy(yaxis_int)
    outputnr = [yaxis_int[0]+1] + outputnr + [yaxis_int[-1]-1]
    noutput = len(outputnr)






    #------------------
    #  Prepare y axis
    #------------------

    # determine how many y axis ticks you want.
    # find step to go through loop

    nyticks_step = int(noutput/10+0.5) 
    if nyticks_step < 2:
        nyticks_step = 2

    yticks = []
    yticks_right = []
    yticks_right_labels=[]
    
    ind = 1
    while ind < noutput - 1:
        
        if ind % nyticks_step == 0:
            yticks.append(outputnr[ind])
            yticks_right.append(outputnr[ind])
            yticks_right_labels.append(round(yaxis_phys[ind-1],3))
    
        ind += 1
    



    #---------------------
    # Set ticks and title
    #---------------------

    # left y ticks
    ax = fig.axes[0]
    ax.set_yticks(yticks)
    ax.set_ylim([outputnr[-1]-0.5, outputnr[0]+0.5])
    ax.set_ylabel('output number', size=14)
    ax.tick_params(axis='both', labelsize=10)

    # right y ticks
    ax2 = ax.twinx()
    ax2.set_yticks(yticks_right)
    ax2.set_ylim([outputnr[-1]-0.5, outputnr[0]+0.5])
    ax2.set_yticklabels(yticks_right_labels)
    if params.use_t:
        ax2.set_ylabel("t [code units]", size=14)
    else:
        ax2.set_ylabel("redshift z", size=14)
    ax2.tick_params(axis='both', labelsize=10)


    # x axis and ticks
    ax.set_xlim(borders)
    ax.set_xticks([]) # no x ticks

    # title
    title = "Merger tree for halo "+str(halo)+" at output "+str(lastdirnr)
    ax.set_title(title, size=18, y=1.01)

    



    #-------------------------
    # Other cosmetics
    #-------------------------

    # add grid
    ax.grid(which='both')

    
    # set figure size
    height = 0.3*len(yaxis_int)
    if height < 12:
        height = 12

    dx = borders[1]-borders[0] - 2
    width = 0.4*dx
    if width < 10:
        width = 10
    
    fig.set_size_inches(width, height)


    # cleaner layout
    plt.tight_layout()
    
    return





#=================================
def _walk_tree(node, root):
#=================================
    """ 
    Walk the tree and count the branches. Add branches to the 
    root/start of that branch.

    Arguments:
        node:   class Node object to check whether a new branch starts here
        root:   class Node object which is current source of the branch.

    returns:
        nothing
    """

    # mark node as walked.
    node.walked = True

    # First check out new branches to mark possible jumps
    # over multiple timesteps as "walked"

    for i, prog in enumerate(node.progs):
        # check only if node hasn't been walked over already:
        if not prog.walked:
            if not node.is_main_prog[i]:
                # a new branch starts!
                # this progenitor will be the root for the new branch.
                root.branches.append(prog)
                prog.start_of_branch = root

                _walk_tree(prog, prog)


    # then just resume where you left off
    for i, prog in enumerate(node.progs):
        if not prog.walked:
            if node.is_main_prog[i]:
                prog.start_of_branch = root
                _walk_tree(prog, root)

    return

#!/usr/bin/python2


#==============
def main():
#==============
    """
    Main function. Calls all the rest.
    """

    #-----------------------
    # Set up
    #-----------------------
    params = global_params()
    params.read_cmdlineargs()
    params.get_output_info()
    params.set_outputfilename()

    if params.verbose:
        params.print_params()



    #----------------
    # read in data
    #----------------
    descendants, progenitors, progenitor_outputnrs, outputnrs, t = read_mergertree_data(params)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #--------------------
    # Debug mode:
    #--------------------

    # In case you want to test stuff out, uncomment the following few lines. They contain
    # dummy arrays for testing. The data doesn't need to be read in anymore, but will
    # be overwritten in case you still read it in.

    #  progenitors =           [ [202], [201,  150], [ 186, 263,   9], [182, 166], \
    #                          [   1,    2,   3,   4], [5,  6,  7, 8,  9, 10, 11, 12], \
    #                          [13, 14, 15, 16, 17, 18, 19, 20] ]
    #  descendants =           [ [203], [202, -202], [-201, 201, 150], [186, 263], \
    #                          [-166, -166, 182, 166], [4, -2, -1, 1, -3,  2, -4,  3], \
    #                          [5,  6,  7, 8,  9, 10, 11, 12] ]
    #
    #  params.noutput = len(progenitors)
    #  progenitor_outputnrs =  [ [  7], [  6,    6], [   5,   5,   2], [  4,   4], \
    #                          [   3,    3,   3,   3], [2,  2,  2, 2,  2,  2,  2,  2], \
    #                          [1,  1,  1, 1,  1,  1,  1,  1]]
    #  t = [ params.noutput-x for x in range(params.noutput) ]
    #  outputnrs = t
    #  params.halo= 203
    #  params.lastdir = 'output_00008'
    #  params.lastdirnr = 8
    #  params.plotparticles = False

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   

    #------------------
    # Make tree
    #------------------
    tree = make_tree(progenitors, descendants, progenitor_outputnrs, params)
    del progenitors
    del descendants
    del progenitor_outputnrs

    #------------------
    # Plot the tree
    #------------------
    plot_tree(tree, outputnrs, t, params)

    #------------------------------
    # If needed: Plot particles
    #------------------------------
    if params.plotparticles:
        plot_treeparticles(tree, t, params)


    return






#===============================
if __name__ == "__main__":
#===============================
    
    main()

