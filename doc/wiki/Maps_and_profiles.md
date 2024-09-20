

In the following, we will assume that a RAMSES output named `output_00001` exists in the current workind directory

# 1. Maps
## 1.1. Fortran tools

## 1.2. Python tools

### 1.2.1. pymses

### 1.2.2. pynbody

### 1.2.3. osiris

### 1.2.4. yt

The following script will extract a 2D map from the output:

    import yt
    # Load a dataset
    ds = yt.load('output_00001/info_00001.txt')
    # This is to make a slice plot and a projection plot
    p = yt.SlicePlot(ds, 'z', 'density')
    p.save()  # save it as a png

    # Do a projection plot (actually integrate over the l.o.s.)
    p = yt.ProjectionPlot(ds, 'z', 'density', weight_field='ones')
    p.save()  # save it as a png

For more information, see [the YT wiki](http://yt-project.org/doc/visualizing/plots.html).

# 2. Profiles

## 2.1. Fortran tools

## 2.2. Python tools

### 2.2.1. pymses

### 2.2.2. pynbody

### 2.2.3. osiris

### 1.2.4. yt