# LuaLB
## A lattice Boltzmann simulation with Lua scripting

**Automatically exported from [code.google.com/p/lualb](https://code.google.com/p/lualb/)** 


**LuaLB** is an open-source parallelizable [lattice Boltzmann simulation](http://en.wikipedia.org/wiki/Lattice_Boltzmann_methods) for computational fluid dynamics with [Lua scripting](http://lua.org)


**News**
* **August 19, 2008** - *Experimentals*: with Lua 5.1.x (compiled with make linux) and HDF5 1.8.x (with configure flag --with-default-api-version=v16) please check out here: svn checkout http://lualb.googlecode.com/svn/trunk/lualb-lua-5.1.x.
* **July 2, 2008** - No binary released now. You need to build from source (see instructions below).

**System requirements**
* Linux or Windows ([Cygwin](http://www.cygwin.com/))
* GCC compatible compiler
* [HDF5 1.6.x](http://hdf.ncsa.uiuc.edu/HDF5/)
* [Lua 5.0.x](http://www.lua.org/)
* [libpng](http://www.libpng.org/pub/png/libpng.html) (optional for graphics)
* [SWIG >= 1.3.31](http://www.swig.org/) ( optional for generating wrapper files)
* [MPICH](http://www-unix.mcs.anl.gov/mpi/mpich1/) (optional for parallel version)


To build the software, please see [the installation page](https://github.com/xderhan/LuaLB/wiki/Installation).


**Credit**
- Firstly developed by Waipot Ngamsaad (ngamsaad.waipot[at]gmail.com)
- This is derived work from the original authors V. Babin and R. Holyst (http://pepe.ichf.edu.pl/lb/index.html)

**References**

*Here is the list of related articles:*

- Shiyi Chen and ­ Gary D.Doolen, " Lattice Boltzmann Method for Fluid Flows", Annual Review of Fluid Mechanics Vol. 30: 329-364 (1998). [link](http://arjournals.annualreviews.org/doi/abs/10.1146/annurev.fluid.30.1.329)
- Swift M. R., Orlandini E., Osborn W. R., Yeomans J. M., "Lattice Boltzmann simulations of liquid-gas and binary fluid systems" Phys. Rev. E 54(5), pp.5041-5052 (1996). [link](http://prola.aps.org/abstract/PRE/v54/i5/p5041_1)
- R.R. Nourgaliev, T.N. Dinh, T.G.Theofanous and D. Joseph, "The Lattice Boltzmann Equation Method: Theoretical Interpretation, Numerics and Implications", Intern. J. Multiphase Flows, 29, pp.117-169 (2003). [link](http://www.sciencedirect.com/science?_ob=ArticleURL&_udi=B6V45-47DM4F8-2&_user=206209&_coverDate=01%2F31%2F2003&_rdoc=7&_fmt=full&_orig=browse&_srch=doc-info(%23toc%235749%232003%23999709998%231%23FLA%23display%23Volume)&_cdi=5749&_sort=d&_docanchor=&view=c&_ct=8&_acct=C000014278&_version=1&_urlVersion=0&_userid=206209&md5=99173a6846c184b139e48643ba534127)
