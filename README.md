# Growth_variability
Scripts used to calculate the datasets underlying the study **“Drivers of increasing tree radial growth variability in Central European forests.”**

## Project Structure

### `Calculated_datasets/`
Recalculated and processed datasets used in statistical analyses.

- `Bootstrapped_crn_variability/` – Site chronologies and bootstrapped variability metrics  
- `Calculated_models/` – Fitted LME models - the full statistical models used in the manuscript are **not included** in this repository. They cannot be reliably estimated on the small demo subset of the data provided here, and therefore are not part of this GitHub version.  
- `Depozitions/` – Annual sulphur and nitrogen deposition values for each site  
- `Finalized_datasets/` – Final datasets used for graphing and model calculations  
- `Recalculated_climate/` – Reconstructed climate variables (temperature, precipitation) for each site  
- `rwi_data/` – Detrended/standardized tree-ring indices for all trees at each site  
- `Site_climate/` – Site-level temperature and precipitation  
- `Soil_data/` – Soil characteristics for each site  

---

### `Data_preparation/`
Scripts used for data cleaning, harmonization, and preprocessing.

### `Graphs/`
Scripts used to generate figures included in the study and analyses related to figure content.

### `Input_data/`
Contains only metadata and small example datasets.  
**Most raw data are not included** due to licensing and file size limitations (details below).

- `Clima_grids/` – netCDF files with temperature and precipitation data.  
  *Not included due to size.*  
- `Database_file/` – CSV with tree metadata.  
  *Only example data provided.*  
- `Depozitions/` – Rasters with annual S and N deposition.  
  *Not included due to size.*  
  **Source:** Oulehle, F. et al. (2016). *Predicting sulphur and nitrogen deposition using a simple statistical method.* Atmospheric Environment 140: 456–468.  
- `Hillshade/` – Elevation model of Central Europe.  
  *Not included.*  
  **Source:** EU-DEM (https://www.eea.europa.eu/data-and-maps/data/eu-dem)  
- `Site_meta/` – CSV with site-level metadata.  
  *Only example data provided.*  
- `Soil_data/` – Detailed soil information.  
  *Not included.*  
  **Source:** Treml, V. et al. (2025). *Rapid trend towards larger and more moisture-limited trees in Central-European temperate forests.* Environmental Research Letters 20: 084033.  
- `Species_distribution/` – Shapefiles of species distributions.  
  *Not included.*  
  **Source:** Caudullo, G., Welk, E., San-Miguel-Ayanz, J. (2017). *Chorological maps for the main European woody species.* Data in Brief 12: 662–666. https://doi.org/10.1016/j.dib.2017.05.007

---

## Tree and Climate Data
Most tree-ring, climate, and site-level data used in this study are stored in the **treedataclim database** (https://treedataclim.cz).  
These datasets have been used in the following publications:

- **Kašpar J et al. (2024)**. *Major tree species of Central European forests differ in their proportion of positive, negative, and nonstationary growth trends.* Global Change Biology 30: e17146.  
- **Mašek J et al. (2024)**. *Shifting climatic responses of tree rings and NDVI along environmental gradients.* Science of The Total Environment 908: 168275.  
- **Treml V et al. (2025)**. *Rapid trend towards larger and more moisture-limited trees in Central-European temperate forests.* Environmental Research Letters 20: 084033.  
- **Tumajer J et al. (2025)**. *Longer growing seasons will not offset growth loss in drought-prone temperate forests of Central–Southeast Europe.* Nature Communications 16: 9535.


