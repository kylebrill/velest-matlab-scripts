# velest-matlab-scripts

These are wrapper functions to facilitate using Velest 1-D velocity modeling with MATLAB. Some of these scripts call on BRTT's Antelope Software, AWK, and GMT. It has been a while since I've used these scripts, so some trial and error might be needed. Some also might be obsolete based on more recent changes. I originally wrote these for use on a Linux machine, but I am pretty sure that I've also used them all in Windows, possibly calling some bash files using Windows Subsystem for Linux.

Reference:
Kissling, E., W.L. Ellsworth, D. Eberhart-Phillips, and U. Kradolfer: Initial reference models in local earthquake tomography, J. Geophys. Res., 99, 19635-19646, 1994. doi: 10.1029/93JB03138

Chet J. Hopp, Gregory P. Waite, Characterization of seismicity at Volcán Barú, Panama: May 2013 through April 2014, Journal of Volcanology and Geothermal Research, Volume 328, 2016, Pages 187-197, ISSN 0377-0273, doi: 10.1016/j.jvolgeores.2016.11.002.

Brill KA, Waite GP and Chigna G (2018) Foundations for Forecasting: Defining Baseline Seismicity at Fuego Volcano, Guatemala. Front. Earth Sci. 6:87. doi: 10.3389/feart.2018.00087


## Getting Started

Read VelestInstructions.docx, it should get you through everything you need. More details on Velest itself can be found in the above reference or here: 
https://ethz.ch/content/dam/ethz/special-interest/erdw/geophysics/seismology-and-geodynamics-dam/document/velest_guide.pdf

### Prerequisites

Depending on what you want to do, other than Velest (https://seg.ethz.ch/software/velest.html), you might want to use AWK (https://www.gnu.org/software/gawk/) and GMT (https://www.soest.hawaii.edu/gmt/). 

The Plots require that a topography file be specified as well. 


## Authors

* **Kyle Brill** - *Work on the MATLAB wrapper* 
* **Chet Hopp** - *Initial random model setup*
* **Greg Waite** - *Advisor* [MichiganTechnologicalUniversity](https://github.com/MichiganTech/volcanoseismo)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
