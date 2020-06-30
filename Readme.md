# Introduction

HiPIMS standards for **Hi**gh-**P**erformance **I**ntegrated hydrodynamic
**M**odelling **S**ystem. It uses state-of-art numerical schemes
(Godunov-type finite volume) to solve the 2D shallow water equations for flood simulations. To support high-resolution flood simulations, HiPIMS is implemented on multiple
GPUs (Graphics Processing Unit) using CUDA/C++ languages to achieve high-performance computing. Since HiPIMS has a modular and flexible structure, it has a great potential to be further developed for other applications in hydrological science as long as the problem can be solved on a uniform rectangular grid. To find out how to use the model, please see the **[wiki](https://github.com/HEMLab/hipims/wiki)**.

# Using the Python APIs

The easiest way to use the model is to use the **[Python APIs](https://pypi.org/project/pypims/)** for hipims, which provides a user-friendly toolkit for data pre-processing, running hipims and result visualisation.

# Contact

The official version of HiPIMS is maintained by the [Hydro-Environmental Modelling Labarotory](http://www.hemlab.org) at Loughborough University. Please contact Qiuhua Liang (q.liang@lboro.ac.uk) for more information.

# License

The code is licensed under GPLv3. Please see LICENCE for more information.
