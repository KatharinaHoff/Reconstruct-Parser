# Reconstruct-Parser User Guide

Author and Contact Information
------------------------------

Katharina J. Hoff

University of Greifswald,
Institute for Mathematics and Computer Science,
Walther-Rathenau-Str. 47,
17489 Greifswald

University of Greifswald,
Center for Functional Genomics of Microbes,
Felix-Hausdorff-Str. 8,
17489 Greifswald

katharina.hoff@uni-greifswald.de

Contents
========

-   [What is the Reconstruct-Parser?](#what-is-reconstruct-parser)
-   [Installation](#installation)
    -   [Dependencies](#dependencies)
    -   [Reconstruct-Parser](#reconstruct-parser)
-   [Data preparation](#data-preparation)
-   [Running Reconstruct-Parser](#running-reconstruct-parser)
-   [Example data](#example-data)
-   [Output of Reconstruct-Parser](#output-of-reconstruct-parser)
-   [Bug reporting](#bug-reporting)
-   [Citing Reconstruct-Parser](#citing-reconstruct-parser)
-   [License](#license)

What is the Reconstruct-Parser?
===============================

Reconstruct-Parser is a script that allows extraction of 3D coordinates from Reconstruct [1, 2] XML-files that contain marked contour information of microscopy images. The output is a text file in flat file format than can easily be imported into R (or other software) for statistical analysis.

For all contour types except for one, Reconstruct-Parser will return the 3D coordinates of the contour itself. However, Reconstruct-Parser was originally developed for statistical analysis of the spatial distribution of gold particles in electron microscopy images. Any contour that translates to ```gold``` in the configuration file is therefore not reported as a list of coordinates, but as a single averaged coordinate (the arithmetic mean of x- and y-coordinates is reported).

Installation
============

Reconstruct-Parser is a plain Perl script (parse_coordinates.pl) for Linux systems that usually does not require specific installation if Perl itself is available.

Dependencies
------------

Reconstruct-Parser requires Perl (it was tested with version 5.26.1) and the Perl modules strict, warnings, Getopt::Long, and File::Basename. On Ubuntu, missing modules can for example be installed with CPANminus (```sudo apt-get install cpanminus```): ```sudo cpanm Module::Name```, e.g. ```sudo cpanm File::Basename```.

Reconstruct-Parser
------------------

We recommend obtaining the parser from github with ```git clone https://github.com/KatharinaHoff/Reconstruct-Parser.git```. You can subsequently execute the script with full path to the script or add the path to the script to your $PATH-variable (e.g. in ~/.bashrc file). For the latter, add the following to the foot of your ~/.bashrc file (adapt path to the actual path on your system!):

```
PATH=/home/user/Reconstruct-Parser:$PATH
```
Load the new bash configuration in your current session with ```source ~/.bashrc```. New bash sessions will automatically load the new configuration.

If you transferred the script via e.g.~a USB stick, the executability flag might get lost. To be safe, make the script exectuable with:

```
chmod u+x parse_coordinates.pl
```

(Alternatively, you also have the option to preceed each call with the executable perl binary.)

Data Preparation
================

Process the contours of interest in Reconstruct to your liking. Pay attention to the naming/colours of contours. If for example you want to mark a nucleus-like structure, follow the same naming-pattern in all slice files (e.g. red followed by a number, which translates to the regular expression pattern red\d+).

Reconstruct parser requires a two-column tabulator separated configuration file that contains contour name patterns and the actual biological meaning assignment, e.g. for above example, you'd state

```
red\d+\tnucleus-like
```
Export the data to XML files. The files for a single cell are typically arranged in one image per slice, and the image ends on a the index number of the slice, e.g. cell1.1, cell1.2, cell1.3, ...

Have a look at the default configuration file ```translate.cfg``` that resides in the same folder as ```parse_coordinates.pl```, but remember that you'll need to adapt the file to the actual contour names that you use in your own project.

Running Reconstruct-Parser
==========================

The Reconstruct-Parser ```parse_coordinates.pl``` has two obligatory arguments:

 * ```--in_stem=STRING``` where STRING is the path to the reconstruct slice files including the slice name stem, e.g. /home/user/files/cell1. (don't forget the dot that preceeds the image index number).
 * ```--out=STRING``` where STRING is the name of the output file.

The argument ```--cfg=STRING``` is optional. STRING should be the name of the configuration file (by default, ```translate.cfg``` is used).

Run the parser as follows:

```
parse_reconstruct.pl --in_stem=example/cell1. --out=test.txt
```

Beware:  Reconstruct-parser was developed for batch-processing, i.e. for processing more than one pile of image slices at a time. Therefore, the output is always appended to the output file if it already exists!

For batch parsing, a script ```batch_parsing.sh``` is provided. You need to edit the script in order to point it to the correct location of your files and possibly the correct configuration file prior execution. (It is executable with the example data by default, though, e.g. try ```batch_parsing.sh```, results will be appended to ```test.txt```.

Example Data
============

Three image slice files are provided as example data in directory example. The slice files start with the in_stem cell.:

```
ls example/
cell1.1  cell1.2  cell1.3

```
Executing the parser with this input set will result in an output file of 471 lines.

Output of Reconstruct-Parser
============================

The output file is a tabulator separated 5-columns text file. Example:

```x       y       z       type    cell
3.29023 2.14294 0.05    cell_wall       example_cell1.
3.30097 2.06652 0.05    cell_wall       example_cell1.
3.30575 2.04742 0.05    cell_wall       example_cell1.

```

The first column contains the x-coordinate, second column contains the y-coordinate (both are extracted from contours in XML-files), the third column contains the z-coordinate (computed from slice thickness in XML file). The fourth column contains the biological label of a coordinate (e.g. cell_wall, membrane, ...) and the fifth column contains the name stem of a cell (the slash between directory name and actual file name_stem is replaced by an underscore, e.g. original example/cell1. is converted to example_cell1.).

Bug reporting
=============

Before reporting bugs, please check that you are using the most recent versions of RECONSTRUCT-Parser. Also, check the open and closed issues on github at <https://github.com/KatharinaHoff/Reconstruct-Parser/issues> for possible solutions to your problem.

Reporting bugs on github
------------------------

If you found a bug, please open an issue at <https://github.com/KatharinaHoff/Reconstruct-Parser/issues>  (or contact katharina.hoff@uni-greifswald.de).


Citing Reconstruct-Parser
=========================

Petersen et al. (2019) Manuscript in preparation.

License
=======

All source code is under GNU public license 3.0 (see
<https://www.gnu.org/licenses/gpl-3.0.de.html>).

References
==========

[1] Fiala, J.C. (2005) Reconstruct: a free ditor for serial section microscopy. Journal of Microscopy 218(1):52-61.

[2] <https://synapseweb.clm.utexas.edu/software-0>, accessed on June 13th 2019.