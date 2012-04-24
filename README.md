USAGE
------

`main.m` details typical usage: loading data, preparing calculations, 
running calculations and making some plots.


INPUTS
------

place the following two files containing your data in a folder 
such as `peach/`:

- `stas.mat` should contain a struct with fields:
      - `stas(i).spikes` the list of spike times of cell i
      - `stas(i).spatial` the spatial component of the STA of cell i

- `cone_params.mat` contains a struct with fields:
      - `stimulus_variance` : variance of each stim pixel color chanel (usually 1)
      - `supersample`       : the integer number of cone positions per pixel 
                              width/height  (usually 4)
      - `colors`            : 3x3 matrix of cone color sensitivities
      - `support_radius`    : radius of cone receptive field filter  (usually 3.0)
      - `repulsion_radii`

see `fullstas2stas.m` for an example of how to make `stas(i).spatial` from raw data


OUTPUTS
-------

the current `main.m` saves results from GREEDY.m, MCMC.m and CAST.m
to the current folder.

a new folder is created by `make_plots.m` reflecting the name of the 
dataset;  all appropriate plots are placed in it.

some plots are in .svg format; in order for these plots to be automatically 
converted to pdf, there must be a folder called `batik/` containing a working
version of [Apache Batik](http://xmlgraphics.apache.org/batik/download.cgi) 
in the same directory as `main.m`.

one of the plots for MCMC and CAST output is `dancing_cones_movie.svg`; 
this movie can be viewed by opening it in a modern browser.


DOCUMENTATION
-------------

most important files are commented; [M2HTML](http://www.artefact.tk/software/matlab/m2html/) is your friend