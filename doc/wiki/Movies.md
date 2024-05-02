

# Movies

The block named `&MOVIE_PARAMS` contains the parameters related to the movies. One can produce up to 5 movies, which center at different parts of the simulation box and have various camera behavior (see example blocks in the bottom). 

The movie routine is very useful for creating on the fly visualizations of a simulation with a small time step between frames, without having to write and process many huge snapshots for that. For turning the binary images created by this routine into a movie that can be watched on any media player, check out [RAM](https://bitbucket.org/biernacki/ram).

| Variable name, syntax, default value | Fortran type  | Description               |
|:---------------------------- |:------------- |:------------------------- |
| `movie`           | `boolean`   | Turns movies on and off |
| `imov`            | `integer `  | Number of starting frame |
| `imovout`         | `integer`   | Total number of frames |
| `tstartmov`       | `float`     | Start time of the movie |
| `astartmov`       | `float`     | Start aexp of the movie |
| `tendmov`         | `float`     | End time of the movie |
| `aendmov`         | `float`     | End aexp of the movie |
| `levelmax_frame`  | `integer`   | Maximum level of the frame |
| `nw_frame`        | `integer`   | Width of frame in px |
| `nh_frame`        | `integer`   | Height of frame in px |
| `xcentre_frame`   | 5*`float, float, float, float` | Four (**for each projection**) coordinates for z position of frame centre [code units] - spline coefficients |
| `ycentre_frame`   | 5*`float, float, float, float` | Four (**for each projection**) coordinates for y position of frame centre [code units] - spline coefficients |
| `zcentre_frame`   | 5*`float, float, float, float` | Four (**for each projection**) coordinates for z position of frame centre [code units] - spline coefficients |
| `deltax_frame`    | `float, float` | Two coordinates for x delta of frame [code units] |
| `deltay_frame`    | `float, float` | Two coordinates for y delta of frame [code units] |
| `deltaz_frame`    | `float, float` | Two coordinates for z delta of frame [code units] |
| `zoom_only_frame` | 5*`boolean`    | Consider only particles in the zoom region of cosmological runs |
| `proj_axis`       | 5*`character`  | Letter code of projection axis; x, y, z - maximum 5 line of sights |
| `movie_vars_txt`  | `strings` | Determines which variables to save - `dens`, `temp`, `pres`, `vx`, `vy`, `vz`, `varX`, `FpX` (RT-only), `stars` (star particles), `dm` (dark matter particles), `lum` (stellar luminosity) |
| `theta_camera`    | 5*`float` | Azimuthal angle of the camera with respect to the line of sight [degrees] |
| `phi_camera`      | 5*`float `    | Polar angle of the camera with respect to the line of sight [degrees] |
| `dtheta_camera`   | 5*`float `    | Azimuthal angle rotation completed between `tstartmov` and `tendmov` [degrees] |
| `dphi_camera`     | 5*`float `    | Polar angle rotation completed between `tstartmov` and `tendmov` [degrees] |
| `focal_camera`    | 5*`float `    | Distance of the focal plane of the camera [code units]. Camera distance is set to `boxlen` and `focal_camera` is set to `boxlen` by default. |
| `dist_camera`    | 5*`float `    | Distance of the camera [code units]. Camera distance is set to `boxlen` and `focal_camera` is set to `boxlen` by default. |
| `ddist_camera`    | 5*`float `    | Motion of the camera between `tstartmov` and `tendmov` [code units]|
| `perspective_camera`| 5*`boolean` | Perspective corrections for the projected cells |
| `shader_frame`| 5*`strings` | Shader applied for the leaf cells projection [`square`, `sphere`, `cube`] |
| `ìvar_frame`| 5*`integer` | Index of hydro variable to use for selecting cells to project (default=-1) |
| `varmin_frame`| 5*`integer` | Only project cells within `varmin_frame<uold(*,ivar_frame)<varmax_frame`. Density is expressed in cm^-3, velocities in km/s, temperature in K/mu, all other hydro variables in code units. |
| `varmax_frame`| 5*`integer` | Only project cells within `varmin_frame<uold(*,ivar_frame)<varmax_frame`. Density is expressed in cm^-3, velocities in km/s, temperature in K/mu, all other hydro variables in code units. |
| `method_frame` | 5*`strings` | Available projection methods: `mean_mass` (default), `mean_dens`, `mean_vol`, `sum`, `min`, `max` |

Examples:

* cosmological simulation
```
#!fortran
&MOVIE_PARAMS
movie=.true.
imov=0
aendmov=1.
imovout=1000
nw_frame=1920
nh_frame=1080
levelmax_frame=18
xcentre_frame=0.49944825,-0.00391524,0.02661462,-0.01714339
ycentre_frame=0.49512072,0.00498905,-0.01826436,0.01097619
zcentre_frame=0.483284613,0.0246328776,-0.0149075526,-1.06750114e-04
deltax_frame=0.0,0.001
deltay_frame=0.0,0.0005625
deltaz_frame=0.0,0.001
proj_axis='z'
movie_vars_txt='dm','stars','temp','dens','var6'
/
```

This will result in 1000 frames 1920x1080 (Full HD) between starting redshift and z=0. Values for `xcentre_frame`, `ycentre_frame` and `zcentre_frame` come from fitting the spline coefficients with `utils/py/get_spline_coeffs.py` (it requires previous run of the same setup with clump finder; e.g. DM-only). Size of the frame will be reflecting physical units (e.g. delta_x=deltax_frame[0]+deltax_frame[1]/a; where a is the explansion factor). Frames are going to be projected along the z axis and projections of dark matter density, stellar density, temperature, gas density and var6 (extra variable, here metallicity) will be saved.

* non-cosmological
```
#!fortran
&MOVIE_PARAMS
movie=.true.
imov=0
tendmov=10.
imovout=1000
nw_frame=1920
nh_frame=1080
levelmax_frame=14
xcentre_frame=5.0,0.,0.,0.,5.0,0.,0.,0.
ycentre_frame=5.0,0.,0.,0.,5.0,0.,0.,0.
zcentre_frame=5.0,0.,0.,0.,5.0,0.,0.,0.
deltax_frame=1.0,0.,1.0,0.
deltay_frame=1.0,0.,1.0,0.
deltaz_frame=1.0,0.,1.0,0.
proj_axis='zy'
movie_vars_txt='dm','stars','temp','dens','var6'
/
```
Same variables as for the cosmological run are going to be saved with the same resolution. Here, the camera is always centred at the centre of the box (assuming boxlen=10) and encompasses volume of 10% x 10% x 10% of the box. Tendmovie is set in the code units of time.