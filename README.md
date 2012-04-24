USAGE
------

`main.m` details typical usage: loading data, preparing calculations, 
running calculations and making some plots.


INPUTS
------

place the following two files containing your data in a folder 
such as `peach/`:

- `stas.mat` should contain a struct with fields:

        * `stas(i).spikes`    : the list of spike times of cell i

        * `stas(i).spatial`   : the spatial component of the STA of cell i

- `cone_params.mat` contains a struct with fields:

        * `stimulus_variance` : variance of each stim pixel color chanel (usually 1)

        * `supersample`       : the integer number of cone positions per pixel width/height  (usually 4)

        * `colors`            : 3x3 matrix of cone color sensitivities

        * `support_radius`    : radius of cone receptive field filter  (usually 3.0)

        * `repulsion_radii`

see `fullstas2stas.m` for an example of how to make `stas(i).spatial` from raw data