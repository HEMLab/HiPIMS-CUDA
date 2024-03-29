# Introduction ![GitHub all releases](https://img.shields.io/github/downloads/HEMLab/HiPIMS-CUDA/total?logo=GitHub)

HiPIMS standards for **Hi**gh-**P**erformance **I**ntegrated hydrodynamic **M**odelling **S**ystem. It uses state-of-art numerical schemes (Godunov-type finite volume) to solve the 2D shallow water equations for flood simulations. To support high-resolution flood simulations, HiPIMS is implemented on multiple GPUs (Graphics Processing Units) using CUDA/C++ languages to achieve high-performance computing. Since HiPIMS has a modular and flexible structure, it has great potential to be further developed for other applications in hydrological science as long as the problem can be solved on a uniform rectangular grid.

# Using HiPIMS

HiPIMS has been embedded in **[Pypims](https://pypims.readthedocs.io/en/latest/)**, the Python APIs for HiPIMS, which provides a user-friendly integrated toolchain for preparing inputs, running HiPIMS and visualising outputs.

# Contributing and Acknowledgement

HiPIMS is developed and maintained by the **[Hydro-Environmental Modelling Labarotory](http://www.hemlab.org)**, a research hub for technological innovation and interdisciplinary collaboration. The development of HiPIMS is led by Professor [Qiuhua Liang](https://www.lboro.ac.uk/departments/abce/staff/qiuhua-liang/) and represents the effort of nearly two decades and many researchers.

# References

The implementation of this code is documented in the following papers

Xia X, Liang Q, Ming X (2019) A full-scale fluvial flood modelling framework based on a high-performance integrated hydrodynamic modelling system (HiPIMS). Advances in Water Resources, 132: 103392.

Xia X, Liang Q (2018) A new efficient implicit scheme for discretising the stiff friction terms in the shallow water equations. Advances in Water Resources, 117: 87-97.

Xia X, Liang Q, Ming X, Hou J (2017) An efficient and stable hydrodynamic model with novel source term discretization schemes for overland flow and flood simulations. Water Resources Research, 53: 3730-3759.

Xia X, Liang Q (2018) A new depth-averaged model for flow-like landslides over complex terrains with curvatures and steep slopes. Engineering Geology, 234: 174-191.

A full list of related publications can be found [here](https://github.com/HEMLab/hipims/wiki/References). Please cite appropriate papers if you use hipims for your project.

# License

The code is licensed under GPLv3. Please see LICENCE for more information.

# Contact us

For more information, please contact [Xiaoli Su](mailto:x.su@lboro.ac.uk), [Jinghua Jiang](mailto:j.jiang3@lboro.ac.uk), [Huili Chen](mailto:H.Chen2@lboro.ac.uk), or [Xue Tong](mailto:x.tong2@lboro.ac.uk).
