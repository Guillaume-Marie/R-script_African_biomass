# R-script_African_biomass
R scripts used to carry out all the analysis and figures of the African biomass paper under review. R2OpenBUGS is an R package that allowed us to BUGS framework into R script. It also mean that R2OpenBUGS is mandatory to be able to run our R script. You can found more detzail on R2OpenBUGS on the link : 
chrome-extension://oemmndcbldboiebfnladdacbdfmadadm/https://cran.r-project.org/web/packages/R2OpenBUGS/R2OpenBUGS.pdf

# Files description

- BA_optim_CWT.r : the R script use to run OpenBUGS 
- model_OPENBUGS.txt : the OpenBUGS model framework use to estimate cover fractions and reference biomasses
- CWT_25-975_refine.csv : the estimated Cross-walking table
- ybio.tar.xz: a data table with all pure land cover type pixels of 1 km² associated wit sum annual precipitation (mm/year) and biomass (t/ha). This dataset is used to feed the OpenBUGS model in order tho esrtimate cover fractions and reference biomasses. 

# Dependencies
The scripts use formated data from original NetCDF,GIF and SHP files. Here are the link of the original files :
- Biomass map of Africa created by CESBIO can be download on demand. It consist of a GIF file in which Africa is spatially discretize in pixel of 1km^2. The unit is tonne of dry mass per hectare (t/ha). Contact person : thuy.letoan@cesbio.cnes.fr 
- Land cover map is freely available here : http://www.esa-landcover-cci.org/. It consist of a GIF file in which Africa is spatially discretize in pixel of 300m^2.
- Ecoregion map use follow the work of Olson et al. 2001. It is used in the figure 3-6. This map is freely available here : https://databasin.org/datasets/68635d7c77f1475f9b6c1d1dbe0a4c4c/
