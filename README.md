
ellShiftedPkg
=============

ETNA package attached to the paper:

**Range Restricted Iterative Methods for the 
Solution of Linear Discrete Ill-posed Problems**
_Alessandro Buccini, Lucas Onisk, and Lothar Reichel_
ETNA Volume XX, Pages yy-zz, year
[link](https://etna.math.kent.edu/url_available_after_acceptance)

---

DESCRIPTION
-----------

This package provides six functions capable of solving large linear and block
linear discrete ill-posed problems when the matrix involved is nonsymmetric.
The functions GMRES.m and ellShiftGMRES.m are capable of solving the former and
the functions BGMRES.m, ellShiftBGMRES.m, glGMRES.m, and ellShiftGlGMRES.m are
capable of solving the latter. Two demos are provided in the package for the
linear case and the block linear case.

Release: 1.0, May 2022

Programming language: Matlab 9.9 (R2020b)

License: (see license.md)

---

INSTALLATION
------------

Download and extract the ellShiftedPkg package. Additionally, visit
https://www.netlib.org/numeralgo/ to download the na33 package. 
Upon extracting this second package, place the functions phillips_alt.m and
shaw_alt.m in the same directory as the extracted ellShiftedPkg package. These
are all the files that are needed to run the demos in ellShiftedPkg package.
The linearSys_demo.m and blockLinearSys_demo.m demos are scripts that do not
need any additional information to run. As long as all necessary files are in
the working directory, the demos may be executed by selecting "Run". See the
documentation in the demo files for expectations of a successful run.

---

PACKAGE USE
------------

As long as all necessary files are in the working directory (see above), the
linearSys_demo.m and blockLinearSys_demo.m demos may be executed by selecting
"Run" without entry of  additional information (i.e. the demos are already
set-up). The linearSys_demo.m demo partially replicates partially replicates
the results and figures of the Shaw example in the paper. The
blockLinearSys_demo.m demo provides a basic example for the block linear case
using two right-hand sides. For more information, see the codePrimer.pdf file
in this package.

---

PACKAGE STRUCTURE
-----------------

The following is a list with a brief description of the contents of the
ellShiftedPkg package. For a more detailed description see the codePrimer.pdf
file in this package.

* codePrimer.pdf        : A detailed overview of the ellShiftedPkg package
                          contents specifically directed towards the functions
                          provided.
* license.md            : Markdown file containing license for package use. 
* README.md             : This file.
* GMRES.m               : The GMRES algorithm for linear discrete ill-posed
                          problems with a square nonsymmetric matrix.
* ellShiftGMRES.m       : The ellShiftGMRES algorithm for linear discrete
                          ill-posed problems with a square nonsymmetric matrix.
                          This version allows the user to select the level of
                          range restriction.
* linearSys_demo.m      : A demo to showcase the use of the GMRES and
                          ellShiftGMRES algorithms on a linear discrete
                          ill-posed problem.
* BGMRES.m              : The BGMRES algorithm for block linear discrete
                          ill-posed problems with a square nonsymmetric matrix. 
* ellShiftBGMRES.m      : The ellShiftBGMRES algorithm for linear discrete
                          ill-posed problems with a square nonsymmetric matrix.
                          This version allows the user to select the level of
                          range restriction.
* glGMRES.m             : The glGMRES algorithm for block linear discrete
                          ill-posed problems with a square nonsymmetric matrix.
* ellShiftGlGMRES.m     : The ellShiftGlGMRES algorithm for block linear
                          discrete ill-posed problems with a square
                          nonsymmetric matrix. This version allows the user to
                          select the level of range restriction.
* blockLinearSys_demo.m : A demo to showcase the use of the BGMRES,
                          ellShiftBGMRES, glGMRES, and ellShiftGlGMRES
                          algorithms on a block linear discrete ill-posed
                          problem.

