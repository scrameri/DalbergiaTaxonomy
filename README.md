This README_CRAMERI_etal_2022.txt file was generated on 2022-03-15 by S. Crameri

**GENERAL INFORMATION**
1. Title of Dataset: 
```
Taxonomic studies on Malagasy Dalbergia (Fabaceae). III. Two New Species from
Southeastern Madagascar and an Emended Description of the Rosewood Species
Dalbergia maritima
```

2. Author Information:

```
A. First Author Contact Information (contact regarding analyses)
	Name:         Simon Crameri
	Institution:  Institute of Integrative Biology, ETH Zurich, Zurich, Switzerland
	Address:      Universitätstrasse 16, 8092 Zürich
	Email:        simon.crameri@usys.ethz.ch ; sfcrameri@gmail.com

B. Last Author Contact Information
	Name:         Alex Widmer
	Institution:  Institute of Integrative Biology, ETH Zurich, Zurich, Switzerland
	Address:      Universitätstrasse 16, 8092 Zürich
	Email:        alex.widmer@usys.ethz.ch
```

3. Date of data collection:
```
2018 - 2021 (plant measurements)
```

4. Geographic location of data collection: 
```
Madagascar
```

5. Information about funding sources that supported the collection of the data:
```
This work was supported by ETH Zurich and a grant from the Rübel Foundation to AW.
The participation of PBP, NW, and PPL was supported by a grant from the Fondation
Franklinia to PPL. Field work and NR’s participation were funded by the Délégation de
l'Union Européenne à Madagascar (DEUM).
```


**SHARING/ACCESS INFORMATION**
1. Licenses/restrictions placed on the data:
```
GNU GENERAL PUBLIC LICENCE version 3
```

2. Links to publications that cite or use the data: 
	[https://doi.org/10.5061/dryad.3n5tb2rhg](https://doi.org/10.5061/dryad.3n5tb2rhg)
	
3. Was data derived from another source?
```
The ecological characteristics were extracted from available spatial raster/vector files:
```					
- bioclimate raster data:	CHELSA version 1.2 Bioclim database (Karger et al. 2017).
- vegetation classes: 		raster of the Atlas of the Vegetation of Madagascar by Moat and Smith (2007)
- surface lithology classes: 	SERVIR database available at http://geoportal.rcmrd.org/data/africa_surface_lithology.zip
- surface boundaries: 		GADM database available through the R raster package (Hijmans and van Etten 2012) version 3.4.5


**References**
- Hijmans, R. J. and J. van Etten. 2012. raster: geographic analysis and modeling with
  raster data. R package version 3.4-5. http://CRAN.R-project.org/package=raster.
- Karger, D. N., O. Conrad, J. Böhner, T. Kawohl, H. Kreft, R. W. Soria-Auza, N. E. 
  Zimmermann, H. P. Linder, and M. Kessler. 2017. Data from: Climatologies at high 
  resolution for the earth’s land surface areas. Dryad Digital Repository.
  https://doi.org/10.5061/dryad.kd1d4.
- Moat, J. and P. Smith (eds.). 2007. Atlas of the Vegetation of Madagascar. Kew: Royal Botanic Gardens.

The spatial raster data from which 17 ecological characteristics were extracted:
```
r.eco_utm_30.tif
```
The extracted data on 17 ecological characteristics:
```
Supplementary_Material.xlsx
```

4. Recommended citation for this dataset: 
```
Crameri, S., P. B. Phillipson, N. Rakotonirina, N. Wilding, R. L. Andriamiarisoa,
P. P. Lowry II, and A. Widmer. 2022. Data from: Taxonomic studies on Malagasy 
Dalbergia (Fabaceae). III. Two new species from southeastern Madagascar and an emended
description of the rosewood species Dalbergia maritima. Dryad Digital Repository.
https://doi.org/10.5061/dryad.3n5tb2rhg.
```

5. Instructions for reproduction of Manuscript Figure 3:
```
1) download the Dryad dataset (https://doi.org/10.5061/dryad.3n5tb2rhg) or the
   Github repository (https://github.com/scrameri/DalbergiaTaxonomy) and place the
   entire unzipped folder in your project working directory.
   
2) open the R script [SupplementaryMaterial_2.R](https://github.com/scrameri/DalbergiaTaxonomy/blob/main/SupplementaryMaterial_2.R) in RStudio, and set the working dir
   to Source File Location (Session/Set Working Directory/To Source File Location).
   
3) run the script line by line using Ctrl/Cmd + Enter	

NOTE:
- the R script should work using R version 4.0.2 and higher, but it depends
- lines that are 'commented out' in the R script are for reference, but cannot be 
  fully reproduced without precise geo-coordinates
  - the extraction of ecological characteristics from spatial raster data (script
    lines 789-951) is shown for reference, but it cannot be reproduced due to non-public
    co-ordinates of endangered rosewood species.	
  - the distribution map (script lines 974 ff) is also shown for reference, but it 
    cannot be reproduced due to non-public co-ordinates of endangered rosewood species.	
- in case of specific questions, please contact Simon Crameri (sfcrameri@gmail.com)
```

**DATA & FILE OVERVIEW**
1. File List:

	[Supplementary_Material.xlsx](https://github.com/scrameri/DalbergiaTaxonomy) (with 4 tabs for Supplementary Materials 1, 3-5):	
	
		Supplementary Material 1	Malagasy Dalbergia collections of six taxa for which 
									morphological measurements were made (n = 424), and of
									the undescribed SAVA material (n = 9).  The last four
									columns denote individuals that were included in
									Figure 3. See column 'TropicosLink' for further
									collection details.												
		
		Supplementary Material 3	13 leaf/leaflet characters of 57 Malagasy Dalbergia
									collections. Continuous values represent medians of
									repeated measurements per individual. See Table 1 for
									details on character codes. NA = not available.														
		
		Supplementary Material 4	11 inflorescence/flower characters of 17 Malagasy
									Dalbergia collections. Continuous values represent
									medians of repeated measurements per individual. See
									Table 1 for details on character codes.	
																				
		Supplementary Material 5	17 ecological characterisitcs of 330 Malagasy
									Dalbergia collections. Latitude, Longitude and
									Elevation are included for reproducibility but
									centered and scaled to unit variance to prevent
									localization of endangered rosewoods. See Table 2 for
									details on character codes. NA = not available.																		
									
	r.eco_utm_30.tif			Spatial raster data of ecological variables for Madagascar
	
	[Tables1-2.docx](https://github.com/scrameri/DalbergiaTaxonomy)					Copy of manuscript Tables 1-2, with details on 
									morphological variables and ecological characteristics
									in Supplementary_Material.xlsx.
	
	[SupplementaryMaterial_2.R](https://github.com/scrameri/DalbergiaTaxonomy/blob/main/SupplementaryMaterial_2.R)		R Script that reads	Supplementary_Material.xlsx
									and generates Manuscript Figure 3. Tested on R 4.0.2
									and 4.1.1.


2. Relationship between files, if important: 
```
SupplementaryMaterial_2.R reads data from the 'data' subfolder. The Rscript and 
subfolder must be located in the same working directory when running the R script.
```

**METHODOLOGICAL INFORMATION**
1. Methods for collection / processing the data: 

Morphological Measurements

```
We assessed 13 leaf and leaflet characters (Supplementary Material 3), along with
11 inflorescence and flower characters (Supplementary Material 4). Continuous and
discrete characters (see Manuscript Table 1), which were measured several times on
a given collection, were recorded as sample medians.
```

Ecological Characteristics
```
We assessed 17 potentially relevant ecological characteristics (Supplementary 
Material 5) from available spatial raster or vector data for Madagascar (see
Manuscript Table 2).
		
We used R (R Core Team 2020) version 4.0.2 and the elevatr package (Hollister et
al. 2020) version 0.3.1, the terra package (Hijmans 2021) version 1.0.10, and the
fasterRaster package (Smith 2020) version 0.6.0 to download high-resolution
elevation data (3 arc seconds, ca. 90 m resolution) and to perform raster
calculations.
		
All rasters were projected to Universal Transverse Mercator (UTM) zone 38S and
re-sampled to the resolution of the highest-resolved raster (30 m) when needed.
We then extracted the ecological characteristics of the selected collections,
resulting in an ecological dataset for multivariate analysis (see 
SupplementaryMaterial_2.R).
		
The projected and re-sampled spatial rasters are in r.eco_utm_30.tif.

```

Multivariate analyses
```
Commented R code documenting the extraction of ecological characteristics from
occurrence data as well as all panels of Figure 3 is available in
SupplementaryMaterial_2.R.
```

2. Instrument- or software-specific information needed to interpret the data: 

```r
> sessionInfo()
R version 4.0.2 (2020-06-22)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Catalina 10.15.7

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] XML_3.99-0.8     rgdal_1.5-27     cleangeo_0.2-4   maptools_1.1-2   rgeos_0.5-8      ggrepel_0.9.1   
[7] ggforce_0.3.3    ggplot2_3.3.5    rmapshaper_0.4.5 concaveman_1.1.0 tmap_3.3-2       alphahull_2.2   
[13] ConR_1.3.0       raster_3.5-2     sp_1.4-6         terra_1.4-22     ape_5.5          cluster_2.1.2   
[19] adegenet_2.1.5   ade4_1.7-18      tabplot_1.4.1    ffbase_0.13.3    ff_4.0.5         bit_4.0.4   

```

3. Describe any quality-assurance procedures performed on the data: 
```
	Some missing values in the raw data were imputed using class means before subjecting
	it to multivariate analyses (PCA, PCoA), see R Script. Missing data code is "NA".
```

4. People involved with sample collection, processing, analysis and/or submission: 
```
Simon Crameri, Nivohenintsoa Rakotonirina
```


