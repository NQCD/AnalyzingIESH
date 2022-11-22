# AnalyzingIESH

[![DOI](https://zenodo.org/badge/567320255.svg)](https://zenodo.org/badge/latestdoi/567320255)

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> AnalyzingIESH

It is authored by James Gardner.

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> ]
   pkg> registry add https://github.com/NQCD/NQCRegistry
   pkg> registry add https://github.com/jamesgardner1421/JamesRegistry
   pkg> add DrWatson
   pkg> activate path/to/this/project
   pkg> instantiate
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.
