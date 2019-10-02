# Gannet 3.1

Open-source, MATLAB-based software for automated data processing and quantification of edited magnetic resonance spectroscopy (MRS) data.

Visit our [website](http://www.gabamrs.com/) for the latest news on Gannet and our developments in edited MRS methodology.

## Getting Started

### Prerequisites

Gannet runs in MATLAB (we recommend using the latest release if possible). Additionally, Gannet requires that the following MATLAB toolboxes are installed:

* Image Processing
* Optimization
* Signal Processing
* Statistics and Machine Learning

To run the voxel co-registration and structural image segmentation modules, [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) must be installed.

### Installing

The simplest way to install Gannet 3.1 is to download the zipped [master folder](https://github.com/richardedden/Gannet3.1/archive/master.zip), unzip it and then move the Gannet3.1-master folder to your MATLAB directory. To add the folder to your search path, run the following line of code in the command window of MATLAB:

macOS
```
addpath(genpath('~/Documents/MATLAB/Gannet3.1-master'));
```

Alternatively, open the Set Path dialog box from the MATLAB menu, click 'Add with Subfolders', find the Gannet3.1-master folder and then select it. When done, press 'Save' to permanently save the Gannet folder to MATLAB's search path.

## Built With

* [MATLAB](https://www.mathworks.com/products/matlab.html)

## Bugs, Contributions, Queries

For bug reporting, contribution requests and/or queries, please contact us at: gabamrs@gmail.com

## Versioning

Semantic versioning is used when updates are made to Gannet using the style 3.1.x. Versioning is also conducted on a module-specific basis using the style YYMMDD. That is, each Gannet module has its own release version.

## Developers

* **Richard Edden** (The Johns Hopkins University)
* **Mark Mikkelsen** (The Johns Hopkins University)
* **Georg Oeltzschner** (The Johns Hopkins University)
* **Sofie Tapper** (The Johns Hopkins University)
* **Muhammad Saleh** (The Johns Hopkins University)
* **C. John Evans** (Cardiff University)
* **Ashley Harris** (University of Calgary)
* **Nicolaas Puts** (The Johns Hopkins University)

## License

This software is open-source and does not have a specific license, but should you publish material that made use of Gannet please cite the following publications (as appropriate):

* Edden RAE, Puts NAJ, Harris AD, Barker PB, Evans CJ. Gannet: A batch-processing tool for the quantitative analysis of gamma-aminobutyric acid-edited MR spectroscopy spectra. J. Magn. Reson. Imaging 2014;40:1445–1452. doi: [10.1002/jmri.24478](http://doi.wiley.com/10.1002/jmri.24478)
* Ashburner J, Friston KJ. Unified segmentation. Neuroimage 2005;26:839–851 doi: [10.1016/j.neuroimage.2005.02.018](https://doi.org/10.1016/j.neuroimage.2005.02.018) (if you perform tissue segmentation)
* Harris AD, Puts NAJ, Edden RAE. Tissue correction for GABA-edited MRS: Considerations of voxel composition, tissue segmentation, and tissue relaxations. J. Magn. Reson. Imaging 2015;42:1431–1440. doi: [10.1002/jmri.24903](http://doi.wiley.com/10.1002/jmri.24903) (if you report tissue-corrected, water-referenced measurements based on the Harris et al. method)
* Gasparovic C, Song T, Devier D, et al. Use of tissue water as a concentration reference for proton spectroscopic imaging. Magn. Reson. Med. 2006;55:1219–1226 doi: [10.1002/mrm.20901](http://doi.wiley.com/10.1002/mrm.20901) (if you report tissue-corrected, water-referenced measurements based on the Gasparovic et al. method)

## Acknowledgments

We wish to thank the following individuals for their contributions to the development of Gannet and shared data processing code:

* Ralph Noeske (GE Berlin)
* Jamie Near (McGill University)
* Peter Barker (The Johns Hopkins University)
* Robin de Graaf (Yale School of Medicine)
* Philipp Ehses (Max Planck Institute for Biological Cybernetics, Tübingen)
* Wouter Potters (UMC Amsterdam)
