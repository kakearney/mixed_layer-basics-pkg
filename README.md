# The mixed_layer model

## Introduction

The mixed_layer model is a Matlab-based suite of functions that simulates advection and diffusion in a one-dimensional water column with prescribed surface forcing.  It was originally written by Charlie Stock and then rewritten and extended by me (Kelly Kearney), and is intended to be a simple physical environment in which to develop and test novel biogeochemical models.  This code (or earlier versions of it) was used to perform the experiments underlying several of my end-to-end model publications, including the following:

 - Kearney KA (2012) An analysis of marine ecosystem dynamics through development of a coupled physical-biogeochemical-fisheries food web model. Princeton University
 - Kearney KA, Stock C, Aydin K, Sarmiento JL (2012) Coupling planktonic ecosystem and fisheries food web models for a pelagic ecosystem: Description and validation for the subarctic Pacific. Ecol Modell 237-238:43–62
   [DOI: 10.1016/j.ecolmodel.2012.04.006](http://dx.doi.org/10.1016/j.ecolmodel.2012.04.006)
 - Kearney KA, Stock C, Sarmiento JL (2013) Amplification and attenuation of increased primary production in a marine food web. Mar Ecol Prog Ser 491:1–14
   [DOI: 10.3354/meps10484](http://dx.doi.org/10.3354/meps10484)
 - Kearney KA, Tommasi D, Stock C (2015) Simulated ecosystem response to volcanic iron fertilization in the subarctic Pacific ocean. Fish Oceanogr 24:395–413
   [DOI: 10.1111/fog.12118](http://dx.doi.org/10.1111/fog.12118)

Variations of the code have also been used for the following studies:

 - Song H, Ji R, Stock C, Kearney K, Wang Z (2011) Interannual variability in phytoplankton blooms and plankton productivity over the Nova Scotian Shelf and in the Gulf of Maine. Mar Ecol Prog Ser 426:105–118
   [DOI: 10.3354/meps09002](http://dx.doi.org/10.3354/meps09002)
 - Bianchi D, Stock C, Galbraith ED, Sarmiento JL (2013) Diel vertical migration: Ecological controls and impacts on the biological pump in a one-dimensional ocean model. Global Biogeochem Cycles 27:478–491
   [DOI: 10.1002/gbc.20031](http://dx.doi.org/10.1002/gbc.20031)
 
## Getting Started

### Prerequisites

This software requires [Matlab](http://www.mathworks.com/products/matlab/).  It was primarily designed using R2012a, but I try to keep it up to date with later versions (let me know if you encounter imcompatibilities).

The full package requires additional Mathworks toolboxes:

- [Curve Fitting Toolbox](http://www.mathworks.com/products/curvefitting/)
- [Image Processing Toolbox](http://www.mathworks.com/products/image/)
- [Parallel Computing Toolbox](http://www.mathworks.com/products/parallel-computing/) (for parallel runs only... can be run without)
- [Statistics Toolbox](http://www.mathworks.com/products/statistics/)
 
(The dependencies on the Curve Fitting and Image Processing toolboxes could probably be removed by altering one or two lines of code; if those dependencies cause you major headaches, let me know and I'll try to modify the code accordingly).
 
### Downloading

Git users can clone (or fork, then clone) directly from this repository.

Alternatively, you may download a zipped version of this code by clicking on "Clone or Download > Download ZIP" in the upper right of the GitHub respository page.
 
### Installation

Once downloaded (and unzipped, if necessary), the following subfolders need to be added to your Matlab path:

```matlab
% Replace with location of the downloaded folder
pth = './mixed_layer-basics-pkg/';

addpath('mixed_layer-basics-pkg/mixed_layer');
addpath('mixed_layer-basics-pkg/mixed_layer/biomodules');
addpath('mixed_layer-basics-pkg/mixed_layer/biomodules_subs');
addpath('mixed_layer-basics-pkg/seawater_ver3_2');
addpath('mixed_layer-basics-pkg/mergestruct');
```

Additional subfolders may need to be added to your path as well, depending on the biological modules you choose to run.  See the [the examples document](https://rawgit.com/kakearney/mixed_layer-basics-pkg/master/mixedlayer_examples.html) for further information.

 
### Contents

The primary functions in this package (`mixed_layer.m` and `runmixedlayer.m`) are found in the `mixed_layer` subfolder.  

The remaining folders hold generic utility functions that the `mixed_layer` function calls during its calculations.  Some of these are my own functions, and some are third-party ones.  License details, where available, are in the appropriate subfolder.  

Below is a brief description of the third-party functions/toolboxes used:

- [wstress](http://woodshole.er.usgs.gov/operations/sea-mat/RPSstuff-html/index.html "RPSstuff"): Rich Signell's function to calculate wind stress
- [mergestruct](http://www.mathworks.com/matlabcentral/fileexchange/7842-catstruct) Originally called catstruct (but renamed due to a filename clash with some other code), concatenates structures.
- [seawater\_ver3\_2](http://www.marine.csiro.au/datacentre/processing.htm): Phil Morgan's seawater toolbox
- [cprintf](http://www.mathworks.com/matlabcentral/fileexchange/24093-cprintf-display-formatted-colored-text-in-the-command-window): utility to print colored text to the command window, which is used when displaying Ecopath models (related to the wce biological module)
- [ConsoleProgressBar](http://www.mathworks.com/matlabcentral/fileexchange/30297-consoleprogressbar): a little utility to show a text progress bar in the command window
- [gridxy](http://www.mathworks.com/matlabcentral/fileexchange/9973-gridxy--v2-2-feb-2008-): plots reference lines (used by Ecopath growth curve plots related to the wce module)
- [plots](http://www.mathworks.com/matlabcentral/fileexchange/10242-plots-m--plotses-m): Plots lines using multiple x/y axes (used by Ecopath growth curve plots related to the wce module)

See [the examples document](https://rawgit.com/kakearney/mixed_layer-basics-pkg/master/mixedlayer_examples.html) for a description of exactly which folders need to be added to your path based on the biology you intend to use.

### The biological modules

The code here comes with several different options for including biogeochemistry, ranging from simple toy models that can be run more or less out of the box to an end-to-end ecosystem model that requires some pretty substantial setup.  From simplest to most complex:

 - **tracer**: The simplest possible model, this simply adds a tracer variable with no sources or sinks. 
 - **np**: A classical nutrient-phytoplankton model, based on the box models described in Chapter 4 of Sarmiento & Gruber (2006).
 - **npz**: A classical nutrient-phytoplankton-zooplankton model, based on the box models described in Chapter 4 of Sarmiento & Guber (2006).
 - **nemurokak**: A variant on the NEMURO model (Kishi et al., 2007), a more realistic biogeochemical model for the subarctic north Pacific Ocean. This variant includes modifications for explicit iron cycling and limitation and some modification to the grazing functional responses tailored to my own work.
 - **cobalt_fweb**: The Carbon, Ocean Biogeochemistry and Lower Trophics (COBALT) model (Stock et al., 2014), a biogeochemical model designed to capture primary and secondary production dynamics in global climate models.
 - **wce**: The water column ecosystem model, which couples traditional biogeochemistry and plankton dynamics to an Ecopath-based upper trophic level fisheries model.  Not for the faint of heart... I'm including this one in the repository primarily for the purpose of scientific transparency.  Setting up the input variables for even a single simulation with this module is complex;  I have a lot of additional code dedicated to parameterizing the multiple-ensemble-member simulations that underly the publications I listed above.  If you're interested in working with this module, please contact me directly.
 
#### References
 
  - Kishi M, Kashiwai M, Ware D, Megrey B, Eslinger D, Werner F, Noguchiaita M, Azumaya T, Fujii M, Hashimoto S (2007) NEMURO- a lower trophic level model for the North Pacific marine ecosystem. Ecol Modell 202:12–25
  - Sarmiento J, Gruber N (2006) Ocean Biogeochemical Dynamics. Princeton University Press
  - Stock CA, Dunne JP, John J (2014) Global-scale carbon and energy flows through the planktonic food web: an analysis with a coupled physical-biological model. Prog Oceanogr 120

## Usage

### Syntax

To run a single instance of the model:

```matlab
mixed_layer(outputfile, param1, val1, param2, val2, ...)
```

Please see the function help (standard Matlab function headers) for a description of the many possible input parameters.  Additional parameters can be found in the function help for the various biological modules.

To run an ensemble run (i.e. multiple simulations using the same spatial and temporal setup and biological module, but varying the other input parameters):

```matlab
runmixedlayer(In, p1, v1, ...)
```

Again, see the function help for a full description of use.

  
### Examples

[The examples document](https://rawgit.com/kakearney/mixed_layer-basics-pkg/master/mixedlayer_examples.html) (and the .m file of the same name, which created it) go through some simple examples of use in detail.




 
 
 

