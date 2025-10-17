# xshelfGeostrophicTransportNWES

This repository contains code for a manuscript titled **"Persistent ocean-shelf transport across the North West European Shelf edge"**, by A. Palóczy, J. Hopkins, A. Wise, J. Huthnance, submitted to the Journal of Geophysical Research: Oceans. This [Jupyter notebook](https://github.com/apaloczy/xshelfGeostrophicTransportNWES/blob/main/index.ipynb) provides an overview of the contents.

The directory `plot_figs/` contains Jupyter notebooks used to produce the figures in the manuscript (Figures 1-7), and the vector graphic file for Figure 8. These notebooks depend on the data files in the `data/` directory. Scripts to generate derived data files are in the `proc/` directory. *TODO: Add code for table transport entries*.

## Python environment

To run the code on your machine, first clone the repository. Then create a mamba/conda environment from the `environment.yml` file:

```
mamba env create -f environment.yml
```

## Abstract

Transport mechanisms between the deep ocean and adjacent continental shelf seas play an important role in the spatial distribution of nutrient delivery to the coastal ocean and in the temporal variability of shelf biogeochemical processes. Along the North West European Shelf (NWES) edge, nutrient-rich waters of oceanic origin are found below the mixed layer, representing a potential nutrient source for fueling new production on the shelf. We find persistent cross-isobath geostrophic transport across the NWES edge in hydrographic climatologies and altimetric sea surface height gradients. This transport is O(1 cm/s), has little vertical structure, and is onshore along the entire extent of the 200 m isobath, except along the southern rim of the Norwegian Trench. Despite strong temporal variability in the shelf-edge hydrography on seasonal to decadal timescales, changes in the ocean-shelf geostrophic transport are subtle. This is due to a persistent large-scale steric sea surface slope along the shelf edge. The geostrophic flow induces local depth-integrated cross-isobath nitrate fluxes of O(1-10 mmol/m/s). This is similar in magnitude to the winter wind-driven nitrate transport, but is much less variable at seasonal and inter-annual time scales. Variability in the geostrophic advection of nitrate is thus determined by the ocean-shelf nitrate gradient's variability, rather than by the cross-isobath flow's variability. Geostrophic transport may therefore be an important baseline component of the nutrient and carbon budgets on the NWES and other continental shelves, and should be considered in their long-term response to climate-scale forcing.

## Authors
* [André Palóczy](https://noc.ac.uk/n/Andre%20Paloczy) (<apaloczy@noc.ac.uk>)
* [Joanne Hopkins](https://noc.ac.uk/n/Joanne%20Hopkins) (<jeh200@noc.ac.uk>)
* [Anthony Wise](https://noc.ac.uk/n/Anthony%20Wise) (<anwise@noc.ac.uk>)
* [John Huthnance](https://noc.ac.uk/n/John%20Huthnance) (<jmh@noc.ac.uk>)

## Acknowledgments

Palóczy and Hopkins were funded by the NERC National Capability Programme Atlantic Climate and Environment Strategic Science (AtlantiS), NE/Y005589/1. Thanks to ICDC, CEN, University of Hamburg for data support (NSBC climatology). This study has been conducted using E.U. Copernicus Marine Service Information; [https://doi.org/10.48670/mds-00337](https://doi.org/10.48670/mds-00337), [https://doi.org/10.48670/moi-00141](https://doi.org/10.48670/moi-00141). EN.4.2.2 data were obtained from [https://www.metoffice.gov.uk/hadobs/en4/](https://www.metoffice.gov.uk/hadobs/en4/) and are ©British Crown Copyright, Met Office, 2025, provided under a Non-Commercial Government Licence [http://www.nationalarchives.gov.uk/doc/non-commercial-government-licence/version/2/](http://www.nationalarchives.gov.uk/doc/non-commercial-government-licence/version/2/).

## Data Availability Statement

Code necessary to reproduce all results is available from [https://github.com/apaloczy/xshelfGeostrophicTransportNWES](https://github.com/apaloczy/xshelfGeostrophicTransportNWES), archived under DOI (to be added). The World Ocean Atlas 2023 temperature, salinity, and nitrate datasets are available from the National Centers for Environmental Information/NOAA website, [https://www.ncei.noaa.gov/access/world-ocean-atlas-2023/](https://www.ncei.noaa.gov/access/world-ocean-atlas-2023/), (Locarnini _et al._, 2024; Reagan _et al._, 2024; Garcia _et al._, 2024). The North Sea Biogeochemical Climatology is available from the University of Hamburg website, [https://www.cen.uni-hamburg.de/en/icdc/data/ocean/nsbc.html](https://www.cen.uni-hamburg.de/en/icdc/data/ocean/nsbc.html) (Hinrichs _et al._, 2017). The EN4 temperature and salinity profile dataset includes the bias correction from Gouretski and Reseghetti (2010) and is available from UKMO (2025). The regional Absolute Dynamic Topography and Mean Dynamic Topography products for European seas are distributed by the Copernicus Marine Environment Monitoring Service (CMEMS, 2024a, 2024b). The SRTM15+ dataset is available from the Institute of Geophysics and Planetary Physics (IGPP, 2025).
