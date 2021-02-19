# LIMahl.jl: Line-intensity modelling approximating halo luminosities, in Julia

LIMahl.jl is a proof-of-concept equivalent of the [lim](https://github.com/pcbreysse/lim) package for Python. The goal of the package is to produce real- and redshift-space calculations of line-intensity power spectra.

Most dependencies are available through the official Julia registry, except for these remote packages:
* Eiichiro Komatsu's [HaloMF.jl](https://github.com/komatsu5147/HaloMF.jl)
* machakann's [DoubleExponentialFormulas.jl](https://github.com/machakann/DoubleExponentialFormulas.jl), which we use for some integrals instead of QuadGK (following Eiichiro Komatsu's MatterPower.jl code)

The current emphasis in development is on getting the package to a working state, rather than performance. Preliminary tests suggest it performs around 2x faster than lim (time-to-P(k) is 0.33 seconds for non-first-run timing in the Julia REPL versus 0.65 seconds with lim in the Python REPL), but there is surely still significant room for improvement here down the line.
