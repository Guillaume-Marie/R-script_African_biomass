# R-script_African_biomass
[R](https://cran.r-project.org/) and [Qgis](https://www.qgis.org/fr/site/) are used to carry out all the analysis and to create figures of the African biomass paper under review. [R2OpenBUGS](https://cran.r-project.org/web/packages/R2OpenBUGS/R2OpenBUGS.pdf) is an R package that allowed to use [OPENBUGS framework](http://www.openbugs.net/w/FrontPage) into [R](https://cran.r-project.org/). It also means that [R2OpenBUGS](https://cran.r-project.org/web/packages/R2OpenBUGS/R2OpenBUGS.pdf) is mandatory to be able to run [our R script](https://github.com/volarex84/R-script_African_biomass/blob/main/BA_optim_CWT.r).

# Files description

- [BA_optim_CWT.r](https://github.com/volarex84/R-script_African_biomass/blob/main/BA_optim_CWT.r) : the R script use to run OpenBUGS 
- [model_OPENBUGS.txt](https://github.com/volarex84/R-script_African_biomass/blob/main/model_OPENBUGS.txt) : the OpenBUGS model framework use to estimate cover fractions and reference biomasses
- [CWT_25-975_refine.csv](https://github.com/volarex84/R-script_African_biomass/blob/main/CWT_25-975_refine.csv) : the estimated Cross-walking table. Cover fraction are available for the 2.5% and 97.5% CI from the estimated distribution obtained in BA_optim_CWT.r. 
- [ybio.tar.xz](https://github.com/volarex84/R-script_African_biomass/blob/main/ybio.tar.xz): a data table with all pure land cover type pixels of 1x1km associated with sum annual precipitation (mm/year) and biomass (t/ha). This dataset is used to feed the [OpenBUGS](http://www.openbugs.net/w/FrontPage) model in order to estimate cover fractions and reference biomasses. 
- [Hybrid_biomass_map.tif](https://github.com/volarex84/R-script_African_biomass/blob/main/Hybrid_biomass_map.tif): it is a biomass map spatially discretized in pixels of 1x1km. It map is obtain according to the method describe in 2.2.3 and use to obtain the data set ybio.tar.xz 

# Dependencies
The scripts use formatted data from original NetCDF,GIF and SHP files. Here are the link of the original files :
- Biomass map of Africa created by [CESBIO](https://www.cesbio.cnrs.fr/) can be downloaded on demand. It consists of a GIF file in which Africa is spatially discretized in pixels of 1x1km. The unit is a tonne of dry mass per hectare (t/ha). [Contact person](thuy.letoan@cesbio.cnes.fr) 
- [Land cover map](http://www.esa-landcover-cci.org/) consist of a GIF file in which Africa is spatially discretized in pixels of 300x300m.
- [Ecoregion](https://databasin.org/datasets/68635d7c77f1475f9b6c1d1dbe0a4c4c/) map use follows the work of Olson et al. 2001. It is used in the figure 3-6.

