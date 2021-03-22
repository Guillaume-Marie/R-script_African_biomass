# R-script_African_biomass
R scripts used to carry out all the analysis and figures of the African biomass paper under review. [R2OpenBUGS](https://cran.r-project.org/web/packages/R2OpenBUGS/R2OpenBUGS.pdf) is an R package that allowed to use BUGS framework into [R](https://cran.r-project.org/). It also means that R2OpenBUGS is mandatory to be able to run our R script.

# Files description

- BA_optim_CWT.r : the R script use to run OpenBUGS 
- model_OPENBUGS.txt : the OpenBUGS model framework use to estimate cover fractions and reference biomasses
- CWT_25-975_refine.csv : the estimated Cross-walking table. Cover fraction are available for the 2.5% and 97.5% CI from the estimated distribution obtained in BA_optim_CWT.r. 
- ybio.tar.xz: a data table with all pure land cover type pixels of 1 km² associated with sum annual precipitation (mm/year) and biomass (t/ha). This dataset is used to feed the OpenBUGS model in order to estimate cover fractions and reference biomasses. 
- Hybrid_biomass_map.tif : it is a biomass map spatially discretized in pixels of 1km^2. It map is obtain according to the method describe in 2.2.3 and use to obtain the data set ybio.tar.xz 

# Dependencies
The scripts use formatted data from original NetCDF,GIF and SHP files. Here are the link of the original files :
- Biomass map of Africa created by CESBIO can be downloaded on demand. It consists of a GIF file in which Africa is spatially discretized in pixels of 1km^2. The unit is a tonne of dry mass per hectare (t/ha). [Contact person](thuy.letoan@cesbio.cnes.fr) 
- [Land cover map](http://www.esa-landcover-cci.org/) consist of a GIF file in which Africa is spatially discretized in pixels of 300m^2.
- [Ecoregion](https://databasin.org/datasets/68635d7c77f1475f9b6c1d1dbe0a4c4c/) map use follows the work of Olson et al. 2001. It is used in the figure 3-6.

