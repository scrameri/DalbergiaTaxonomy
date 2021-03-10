##################################################################
### Multivariate Analyses on Morphological and Ecological data ###
##################################################################

### Date and Author
# author: simon.crameri@env.ethz.ch
# date: March 2021


# library(devtools)
# install_github("mtennekes/tabplot")
library(tabplot) # tableplot()
library(adegenet) # transp() funky()
library(cluster) # daisy()
library(ape) # pcoa()
library(terra) # rast() extract()

### Create spatial rasters
# library(fasterRaster)
# library(enmSdm)
# library(elevatr)
# proj.MDG <- "+proj=utm +zone=38 +south +datum=WGS84 +units=m +no_defs"
# ## Administrative border of Madagascar
# for2017 <- raster("../../Ecology/MadagascarForestModel/RasterData/forest/for2017.tif") # from Extension to Vieilledent et al. (2018)
# grassDir <- c("/Applications/GRASS-7.8.app/Contents/Resources") # on SC's Mac
# mdg <- raster::getData("GADM",country="Madagascar",level=0,download=T) # SpatialPolygonsDataFrame from GADM
# mdg_utm <- spTransform(mdg, CRSobj = proj.MDG)
# MDG_border_utm_30 <- fasterRasterize(vect = mdg_utm, rast = for2017, grassDir = grassDir, outGrassName = "MDG_border_utm_30")
# MDG_border_utm_30 <- writeRaster(MDG_border_utm_30, filename = "../../Ecology/MadagascarForestModel/RasterData/terrain/MDG_border_utm_30.tif", format = "GTiff")
# ## Long/Lat rasters
# ll <- enmSdm::longLatRasters(MDG_border_utm_30)
# ll <- ll * MDG_border_utm_30
# names(ll) <- c("longitude", "latitude")
# writeRaster(ll, filename = "RasterData/terrain/MDG_LongLat_utm_30", datatype = "INT4S")
# ## Elevation
# alt <- get_elev_raster(locations = mdg, prj = "+proj=longlat +ellps=WGS84 +no_defs", z = 15, src = "gl3") # raster of elevation (3 arcseconds, c. 80-90m))
# MDG_alt_utm_30_orig <- projectRaster(from = alt, to = MDG_border_utm_30, method = "bilinear",filename = "RasterData/terrain/MDG_alt_utm_30_orig", format = "GTiff") # original altitude (to compute terrain variables)
# MDG_alt_utm_30 <- MDG_alt_utm_30_orig * MDG_border_utm_30
# MDG_alt_utm_30 <- round(100 * MDG_alt_utm_30)
# writeRaster(MDG_alt_utm_30, file.path(folder, paste0("MDG_alt_utm_30.tif")), datatype = "INT4S") # rounded altitude*100 to save disk space
# ## Aspect, Slope (and other terrain variables derived from the elevation raster)
# for (opt in c("aspect","flowdir","roughness","slope","TPI","TRI")) {
#   cat("processing", opt, "...\n")
#   r <- terrain(x = MDG_alt_utm_30_orig, opt = opt, unit = "degrees", neighbors = 8)
#   r <- r * MDG_border_utm_30
#   if (opt != "aspect") {r <- round(100 * r)}
#   writeRaster(r, file.path("../../Ecology/MadagascarForestModel/RasterData/terrain", paste0("MDG_", opt, "_utm_30.tif")), datatype = "INT4S")
# }
# ## Climate
# library(terra)
# tdir <- "../../Ecology/DataRastersVectors/Climate/CHELSA_Madagascar/terra"
# re-write input rasters for terra (projection is not working on the original rasters cropped and written with the *raster* package)
# MDG_border_utm_30 <- rast("../../Ecology/MadagascarForestModel/RasterData/terrain/MDG_border_utm_30.tif")
# terra::writeRaster(MDG_border_utm_30, filename = file.path(tdir, paste0("DMG_border_utm_30.tif")),wopt = list(filetype = "GTiff", gdal = c("COMPRESS=LZW", "TFW=NO"), datatype = "INT1U"), overwrite = TRUE)
# p.bioclim <- list.files("../../Ecology/DataRastersVectors/Climate/CHELSA_Madagascar/cropped", pattern = ".tif$", full.names = TRUE) # downloaded from CHELSA, unzipped and cropped to extent(c(42,52,-26,-11))
# bioclim <- do.call(what = "sds", args = list(p.bioclim))
# for (i in 1:length(bioclim)) {
#   rast <- bioclim[[i]]
#   oname <- gsub("CHELSA_bio10_", "CHELSA_Bio", names(rast))
#   r30 <- project(x = rast, y = MDG_border_utm_30, method = "bilinear")
#   names(r30) <- paste0(oname, "_utm_30")
#   r30 <- round(r30)
#   writeRaster(r30, filename = file.path(tdir, paste0(oname, "_utm_30.tif")),wopt = list(filetype = "GTiff", gdal = c("COMPRESS=LZW", "TFW=NO"), datatype = "INT2U"), overwrite = TRUE)
# }
# ## Distance to Coast
# out <- fasterRastDistance(MDG_border_utm_30, metric = "euclidean", meters = TRUE,fillNAs = FALSE, grassDir = grassDir)
# out <- out * MDG_border_utm_30
# out <- round(out)
# names(out) <- "MDG_dist2Coast_utm_30"
# writeRaster(out, "RasterData/terrain/MDG_dist2Coast_utm_30", datatype = "INT4S")
# ## Distance to Lakes
# vect_ll <- shapefile("../../Ecology/MadagascarForestModel/RasterData/DIVA-GIS/MDG_wat/MDG_water_areas_dcw.shp")
# vect <- sp::spTransform(vect_ll, CRSobj = proj.MDG)
# vect <- crop(vect, MDG_border_utm_30) # nothing is cropped
# out <- fasterVectToRastDistance(vect = vect, rast = MDG_border_utm_30, metric = "euclidean", meters = TRUE, grassDir = grassDir)
# out <- out * MDG_border_utm_30
# out <- round(out)
# names(out) <- "MDG_dist2Lakes_utm_30"
# r.dist2Lakes <- writeRaster(out, "../../Ecology/MadagascarForestModel/RasterData/terrain/MDG_dist2Lakes_utm_30.tif", datatype = "INT4S")
# ## Distance to Rivers
# vect_ll <- shapefile("../../Ecology/MadagascarForestModel/RasterData/DIVA-GIS/MDG_wat/MDG_water_lines_dcw.shp")
# vect <- sp::spTransform(vect_ll, proj.MDG)
# vect <- crop(vect, MDG_border_utm_30) # nothing is cropped
# out <- fasterVectToRastDistance(vect = vect, rast = MDG_border_utm_30, metric = "euclidean", meters = TRUE, grassDir = grassDir)
# out <- out * MDG_border_utm_30
# out <- round(out)
# names(out) <- "MDG_dist2Rivers_utm_30"
# r.dist2Rivers <- writeRaster(out, "RasterData/terrain/MDG_dist2Rivers_utm_30.tif", datatype = "INT4S")
# ## Distance to Inland Waters
# water <- raster::stack(r.dist2Lakes, r.dist2Rivers)
# beginCluster(7)
# out <- calc(water, min)
# endCluster()
# names(out) <- "MDG_dist2InlandWaters_utm_30"
# writeRaster(out, "RasterData/terrain/MDG_dist2InlandWaters_utm_30", datatype = "INT4S")
## Geology
# r.geol.ll <- crop(raster("../../Ecology/DataRastersVectors/Geology_Lithology/Africa_Surface_Lithology/Africa_Surface_Lithology.tif"), extent(c(42,52,-26,-11))) # longlat
# r.geol.origres.utm <- projectRaster(r.geol.ll, crs = CRS(proj.MDG), method = "nearest", grassDir = grassDir)
# r.geol <- writeRaster(r.geol.origres.utm, filename = "../../Ecology/DataRastersVectors/Geology_Lithology/Africa_Surface_Lithology/Madagascar_Surface_Lithology_utm.tif", datatype = "INT1U", format = "GTiff") # utm, original resolution
# r.geol_utm <- fasterProjectRaster(r.geol, MDG_border_utm_30, method = "nearest", grassDir = grassDir)
# writeRaster(r.geol_utm, filename = "../../Ecology/MadagascarForestModel/RasterData/terrain/MDG_geol_utm_30.tif", datatype = "INT1U", format = "GTiff") # utm, 30x30 resolution
## Vegetation
# r.vege <- rast("../../Ecology/DataRastersVectors/Vegetation/AtlasOfTheVegetationOfMadagascar/AtlasOfTheVegetationOfMadagascar2007.tif") # merc, c. 30x30 m
# r.vege_utm <- project(r.vege, MDG_border_utm_30, method = "near")
# terra::writeRaster(r.vege_utm, filename = "../../Ecology/DataRastersVectors/Vegetation/AtlasOfTheVegetationOfMadagascar/AtlasOfTheVegetationOfMadagascar2007_utm_30.tif",wopt = list(filetype = "GTiff", gdal = c("COMPRESS=LZW", "TFW=NO"), datatype = "INT1U"), overwrite = TRUE)
            
### Parameters
# paths
# p.leaves <- "../../Pictures/ZT_digitalization/Chapter2/results_shape1-721/d.morph.708.rda" # leaf shape data
# p.shape <- "../../THESIS/CHAPTER_3/Morphology/shapes.xlsx" # leaf shape data
p.morpho <- "../../THESIS/CHAPTER_3/Morphology/charactertable.rda" # flower, fruit, leaf data
p.data <- "specimens.txt" # collections for ecological data extractions

# spatial rasters (UTM, 30m resolution)
r.longlat <- rast("../../Ecology/MadagascarForestModel/RasterData/terrain/MDG_LongLat_utm_30.tif") # two layers
r.prec.annual <- rast("../../Ecology/DataRastersVectors/Climate/CHELSA_Madagascar/terra/CHELSA_Bio12_Annual_Precipitation_utm_30.tif")
r.prec.seasonality <- rast("../../Ecology/DataRastersVectors/Climate/CHELSA_Madagascar/terra/CHELSA_Bio15_Precipitation_Seasonality_utm_30.tif")
r.temp.annual <- rast("../../Ecology/DataRastersVectors/Climate/CHELSA_Madagascar/terra/CHELSA_Bio01_Annual_Mean_Temperature_utm_30.tif")
r.temp.seasonality <- rast("../../Ecology/DataRastersVectors/Climate/CHELSA_Madagascar/terra/CHELSA_Bio04_Temperature_Seasonality_utm_30.tif")
r.isothermality <- rast("../../Ecology/DataRastersVectors/Climate/CHELSA_Madagascar/terra/CHELSA_Bio03_Isothermality_utm_30.tif")
r.elevation <- rast("../../Ecology/MadagascarForestModel/RasterData/terrain/MDG_alt_utm_30.tif")
r.aspect <- rast("../../Ecology/MadagascarForestModel/RasterData/terrain/MDG_aspect_utm_30.tif")
r.slope <- rast("../../Ecology/MadagascarForestModel/RasterData/terrain/MDG_slope_utm_30.tif")
r.dist2Coast <- rast("../../Ecology/MadagascarForestModel/RasterData/terrain/MDG_dist2Coast_utm_30.tif")
r.dist2InlandWaters <- rast("../../Ecology/MadagascarForestModel/RasterData/terrain/MDG_dist2InlandWaters_utm_30.tif")
r.geol <- rast("../../Ecology/MadagascarForestModel/RasterData/terrain/MDG_geol_utm_30.tif")
r.vege <- rast("../../Ecology/DataRastersVectors/Vegetation/AtlasOfTheVegetationOfMadagascar/AtlasOfTheVegetationOfMadagascar2007_utm_30.tif")

# categorical variables
categorical <- c("geol","vege")
min.frac.dummy <- 0.05 # at least this fraction of observations need to have a categorical level for it to be kept
na.strings <- c("NA","",NA)

# variable names
geovars <- c("LongitudeDecimal","LatitudeDecimal")
shpvars <- c("shp1_orbicular","shp3_ovate","shp4_elliptic")
ecovars <- c("prec.annual","prec.seasonality","temp.annual","temp.seasonality", "isothermality","elevation","aspect","slope","dist2Coast","dist2InlandWaters","geol","vege")

idvar <- "Collection"
classvar <- "Species"
n.meta <- c(idvar,classvar)

# collection filter
identifyers <- "Crameri|Phillipson|Wilding"
accurate.geo <- c("GPS","google earth","map","post facto, high", "corrected post facto")

# classes
spec <- c("chapelieri","louvelii","maritima_pubescens","maritima_maritima","normandii","occulta","pseudomaritima","razakamalalae")
cols.all <- adegenet::funky(length(spec))

# plotting
plot.size <- 10

# pcoa types (binary values would be handled automatically by cluster::daisy)
# ordratio = "^PUB_|^MAR_leaflet" # grep ordinal factors
# asymm = "^CORI_|^VEN_high_low|^DENS_infl|^REL_infl_len|TYPE_smooth_pericarp|^TYPE_ven_pericarp|^FISS_pericarp|^PERS_calyx|^vege_|^geol_|^eco_" # grep binary factors
ordratio = "^PUB_|^MAR_leaflet" # grep ordinal factors
asymm = "^CORI_|^TYPE_infl|^vege_|^geol_" # grep binary nominal factors (dummy variables)


#####################################################################################################

### Helperfunctions
# numerize
numerize <- function(x) {
  t <- suppressWarnings(as.numeric(as.character(x)))
  if (all(is.na(t))) {
    t <- rep(as.numeric(NA), length(x))
  }
  return(t)
}

# factor2dummy
factor2dummy <- function(df, factorname, keep.levels = levels(df[,factorname]), na.method = "zero") {
  
  dfac <- data.frame(factor(df[,factorname]))
  names(dfac) <- factorname
  keep.levels <- levels(dfac[,factorname])[levels(dfac[,factorname]) %in% keep.levels]

  m <- model.matrix(~dfac[,factorname])
  m <- m[match(rownames(dfac), rownames(m)), , drop = FALSE]
  m[,1] <- ifelse(dfac[,factorname] == levels(dfac[,factorname])[1], 1, 0) # first column is a dummy intercept (first level not returned in model.matrix)
  colnames(m) <- paste0(factorname, "_", levels(dfac[,factorname]))
  rownames(m) <- rownames(dfac)
  if (na.method == "zero") m[is.na(m)] <- 0
  
  stopifnot(all.equal(sum(m), sum(table(dfac[,factorname]))))
  df2 <- df[,!colnames(df) %in% factorname]
  pos <- which(colnames(df) == factorname)
  if (pos == 1) {
    df2 <- cbind(m[,paste0(factorname, "_", keep.levels),drop=FALSE], df2)
  } else if (pos == ncol(df)) {
    df2 <- cbind(df2, m[,paste0(factorname, "_", keep.levels),drop=FALSE])
  } else {
    df2 <- cbind(df2[,1:(pos-1),drop=FALSE], m[,paste0(factorname, "_", keep.levels),drop=FALSE], df2[,pos:ncol(df2),drop=FALSE])
  }
  names(df2)[pos:(pos + length(keep.levels) - 1)] <- paste0(factorname, "_", keep.levels)
  rownames(df2) <- rownames(df)
  return(df2)
}

# imputation
impute <- function(df, fun = "mean", classes = NULL, summary = TRUE) {
  
  ## Author: simon.crameri@env.ethz.ch
  
  ## Check input ##
  stopifnot(inherits(df, "data.frame"),
            fun %in% c("mean", "zero", "class.mean"))
  if (!is.null(classes)) {
    stopifnot(length(classes) == nrow(df),
              is.factor(classes))
  }
  
  ## Get missingness summary ##
  missingness <- list()
  
  if (summary) {
    xmis <- df
  
    # total missingness (displayed in summary(genind))
    dfmis <- apply(xmis, 2, function(x) {gsub("NA", NA, x)})
    df.perc <- sum(is.na(dfmis)) / (nrow(dfmis)*ncol(dfmis))
    missingness[["total"]] <- df.perc
    
    # missingness per variable
    df.var.perc <- apply(dfmis, 2, function(x) {sum(is.na(x)) / length(x)})
    missingness[["per.variable"]] <- df.var.perc
    
    # missingness per individual
    df.ind.perc <- apply(dfmis, 1, function(x) {sum(is.na(x)) / length(x)})
    missingness[["per.individual"]] <- df.ind.perc
    
  }
  
  ## Define impute function ##
  impute <- function(x, fun, ...) {
    
    # check input
    stopifnot(is.numeric(x))
    if (any(is.na(x))) {
      
      # impute missing data
      switch(fun,
             mean = {x[is.na(x)] <- mean(x, na.rm = TRUE)},
             zero = {x[is.na(x)] <- 0}
      )
    }
    
    # return imputed vector
    return(x) 
  }
  
  ## Impute missing data ##
  switch(fun, 
         mean = {
           missingness[["imputed"]] <- data.frame(apply(df, MARGIN = 2, FUN = impute, fun = "mean"))
         }, 
         zero = {
           missingness[["imputed"]] <- data.frame(apply(df, MARGIN = 2, FUN = impute, fun = "zero"))
         },
         class.mean = {
           tab <- array(NA, dim = c(0, ncol(df)), dimnames = list(NULL, colnames(df)))
           for (i in levels(classes)) {
             dfsub <- df[classes == i,]
             tmp <- apply(dfsub, MARGIN = 2, FUN = impute, fun = "mean", simplify = F)
             if (!is.array(tmp)) tmp <- data.frame(t(tmp))
             rownames(tmp) <- rownames(dfsub)
             colnames(tmp) <- colnames(tab)
             tab <- rbind(tab, tmp)
           }
           stopifnot(all(rownames(tab) %in% rownames(df)))
           tab <- tab[rownames(df),]
           missingness[["imputed"]] <- data.frame(tab)
         }
  )
  
  
  ## Get imputation info ##
  switch(fun, 
         mean = {
           mean.imp <- apply(df, 2, mean, na.rm = TRUE)
           missingness[["imputation"]] <- mean.imp
         },
         zero = {
           zero.imp <- rep(0, ncol(df))
           names(zero.imp) <- colnames(df)
           missingness[["imputation"]] <- zero.imp
         },
         class.mean = {
           class.mean.imp <- apply(df, 2, tapply, classes, mean, na.rm = TRUE)
           missingness[["imputation"]] <- class.mean.imp
         }
  )
  
  ## return results
  return(missingness)
}

# remove constant variables
strip.constant <- function(df, verbose = FALSE) {
  kept <- apply(df, 2, function(x) {length(unique(x[!x%in%na.strings])) > 1})
  if ((length(kept) < ncol(df)) & verbose) {
    message("removed the following ", ncol(df)-length(kept), " constant variables:\n", paste(colnames(df)[!colnames(df) %in% kept], collapse = ","))
  }
  return(df[,kept])
}

# multivariate analyses
ordplot <- function(df, impute.fun = "class.mean", x = 1, y = 2, classes = factor(rep("", nrow(df))),
                    cols = adegenet::funky(nlevels(droplevels(classes))),
                    scale = TRUE, center = TRUE, pca = TRUE, pcoa = TRUE,
                    hulls = TRUE, voronoi = FALSE, biplot = TRUE, biplot.quantile = 0, 
                    biplot.scaling = 0.8, voronoi.scaling = 0.1,
                    pcoa.metric = "gower", ordratio = NULL, asymm = NULL, title = substitute(df), verbose = TRUE) {
  
  # author: simon.crameri@env.ethz.ch, Mar 2021
  library(ggforce)
  library(ggrepel)
  
  # helperfunction
  strip.constant <- function(df, verbose = FALSE) {
    kept <- apply(df, 2, function(x) {length(unique(x[!x%in%na.strings])) > 1})
    if ((length(kept) < ncol(df)) & verbose) {
      message("removed the following ", ncol(df)-length(kept), " constant variables:\n", paste(colnames(df)[!colnames(df) %in% kept], collapse = ","))
    }
    return(df[,kept])
  }
  
  # subset data rows
  sel <- seq(nrow(df))
  
  # subset data columns
  df.sub <- df
  df.sub <- strip.constant(df.sub, verbose = verbose)
  
  # remove non-numeric variables
  d.num <- df.sub[,sapply(df.sub, is.numeric)]
  
  # fill NAs with means
  classes <- droplevels(classes[sel])
  d.imp <- impute(df = d.num, fun = impute.fun, classes = classes)
  if (impute.fun == "class.mean") {
    # if there is no class mean, use global mean
    d.imp <- impute(df = d.imp$imputed, fun = "mean")
  }
  
  # prepare datasets
  d.pca <- d.imp$imputed
  colnames(d.pca) <- colnames(d.num) # imputed numerical variables
  # d.pcoa <- cbind(d.pca, df.sub[,!colnames(df.sub) %in% colnames(d.pca)]) ## all variables, possibly including (ordinal / nominal) factors with NAs
  d.pcoa <- d.pca # only (imputed) numeric variables

  #########
  ## PCA ##
  #########
  if (pca) {
    d.pca.res <- prcomp(x = d.pca, scale = scale, center = center)
    res.pca <- data.frame(d.pca.res$x)
    res.pca$Class <- classes
    
    # eigenvectors / rotation matrix
    rot <- d.pca.res$rotation
    
    # eigenvalues / % variance
    dvar <- round(100*(d.pca.res$sdev^2 / sum(d.pca.res$sdev^2)), 2)
    
    # biplot scaling
    scaler <- min(max(abs(res.pca[, x]))/max(abs(rot[,x])),
                  max(abs(res.pca[, y]))/max(abs(rot[,y]))) * biplot.scaling
    
    # biplot filtering
    qt <- apply(abs(rot), 2, quantile, biplot.quantile, na.rm = TRUE)
    torm <- abs(rot[,x]) < qt[x] & abs(rot[,y]) < qt[y]
    if ((sum(torm) > 0) & verbose) {
      message("removed the following ", sum(torm), " variables from PCA biplot:\n", paste(names(torm[torm]), collapse = ","))
    }
    rot <- data.frame(rot[!torm,])
    
    # biplot
    p1 <- ggplot(data = res.pca, aes_string(x = paste0("PC", x), y = paste0("PC", y), colour = "Class"))
    if (hulls) {
      p1 <- p1 +
        geom_mark_hull(aes(fill = Class))
    }
    if (voronoi) {
      p1 <- p1 +
        geom_voronoi_tile(aes(fill = Class, group = -1L), 
                          max.radius = diff(range(res.pca[,x]))*voronoi.scaling, colour = NA,
                          alpha = 0.5, show.legend = TRUE)
    }
    p1 <- p1 +
      geom_point()
    if (biplot) {
      p1 <- p1 +
        geom_segment(data = rot, aes(x = rep(0, nrow(rot)), y = rep(0, nrow(rot)),
                                     xend = scaler*rot[,x], yend = scaler*rot[,y]),
                     arrow = arrow(angle = 20, length = unit(5, "mm")), colour = "tomato") +
        geom_label_repel(data = rot, aes(x = scaler*rot[,x], y = scaler*rot[,y], label = rownames(rot)),
                         label.size = 0.2, inherit.aes = FALSE)
    }
    p1 <- p1 +
      scale_colour_manual(values = cols) +
      scale_fill_manual(values = cols) +
      coord_fixed() +
      labs(x = paste0("PC ", x, " (", dvar[x], "%)"), y = paste0("PC ", y, " (", dvar[y], "%)"), fill = "", colour = "") +
      theme_bw() +
      ggtitle(paste0("PCA of ", title, " (nind = ", nrow(d.pca), " ; nvar = ", ncol(d.pca), ")"))
    print(p1)
  }
      
  
  ##########
  ## PCoA ##
  ##########
  if (pcoa) {
    type <- list()
    if (!is.null(ordratio)) {
      # check if variable is truly a nominal factor
      d.ord <- d.pcoa[,grep(ordratio, colnames(d.pcoa)), drop = FALSE]
      lu <- apply(d.ord, 2, function(x) length(unique(x)))
      ordratio <- names(lu)[lu > 2]
      asymm.add <- names(lu)[lu == 2]
      type[["ordratio"]] <- ordratio
      if (verbose) message("treating as ordinal variables (ratio scaled):\n", paste0(ordratio, collapse = ","))
    } else {
      asymm.add <- character()
    }
    if (!is.null(asymm)) {
      # check if variable is truly a binary factor
      # d.asym <- d.pcoa[,grep(asymm, colnames(d.pcoa)), drop = FALSE]
      d.asym <- d.pcoa
      as <- apply(d.asym, 2, function(x) length(unique(x)))
      asymm <- c(names(as)[as == 2], asymm.add)
      asymm <- asymm[!duplicated(asymm)]
      
      # rescale binary variables to [0,1]
      d.pcoa[,asymm] <- apply(d.pcoa[,asymm,drop=F], 2, scales::rescale)
      
      if (verbose) message("treating as binary variables (asymmetric binary scaled):\n", paste0(asymm, collapse = ","))
      type[["asymm"]] <- asymm
    }
    weights <- rep(1, ncol(d.pcoa))
    d.pcoa.dist <- cluster::daisy(x = d.pcoa, metric = pcoa.metric, stand = scale, weights = weights, type = type)
    d.pcoa.res <- ape::pcoa(D = d.pcoa.dist, correction = "cailliez")
    
    res.pcoa <- data.frame(d.pcoa.res$vectors)
    res.pcoa$Class <- classes
    
    # eigenvectors / rotation matrix
    plot.axes <- c(x,y)
    pr.coo <- d.pcoa.res$vectors
    pr.ev <- d.pcoa.res$values
    Y <- scale(d.pcoa, center = center, scale = scale)
    points.stand <- scale(pr.coo[,plot.axes])
    # points.stand <- pr.coo[,plot.axes]
    n <- nrow(Y)
    S <- cov(Y, points.stand, use = "na.or.complete")
    rot <- S %*% diag((pr.ev$Eigenvalues[plot.axes]/(n - 1))^(-0.5))
    
    # eigenvalues / % variance
    dvar <- round(100*(d.pcoa.res$values$Eigenvalues / sum(d.pcoa.res$values$Eigenvalues)), 2)
    
    # biplot scaling
    scaler <- min(max(abs(pr.coo[, x]))/max(abs(rot[,1])),
                  max(abs(pr.coo[, y]))/max(abs(rot[,2]))) * biplot.scaling
    
    # biplot filtering
    qt <- apply(abs(rot), 2, quantile, biplot.quantile, na.rm = TRUE)
    torm <- abs(rot[,x]) < qt[x] & abs(rot[,y]) < qt[y]
    if ((sum(torm) > 0) & verbose) {
      message("removed the following ", sum(torm), " variables from PCoA biplot:\n", paste(names(torm[torm]), collapse = ","))
    }
    rot <- data.frame(rot[!torm,])
    
    # biplot
    # biplot(d.pcoa.res, Y = Y, rn = rep("o", nrow(d.pcoa.res$vectors))) # compare with biplot()
    p2 <- ggplot(data = res.pcoa, aes_string(x = paste0("Axis.", x), y = paste0("Axis.", y), colour = "Class"))
    if (hulls) {
      p2 <- p2 +
        geom_mark_hull(aes(fill = Class))
    }
    if (voronoi) {
      p2 <- p2 +
        geom_voronoi_tile(aes(fill = Class, group = -1L), 
                          max.radius = diff(range(res.pcoa[,x]))*voronoi.scaling, colour = NA,
                          alpha = 0.5, show.legend = TRUE)
    }
    p2 <- p2 +
      geom_point()
    if (biplot) {
      p2 <- p2 +
        geom_segment(data = rot, aes(x = rep(0, nrow(rot)), y = rep(0, nrow(rot)),
                                     xend = scaler*rot[,x], yend = scaler*rot[,y]),
                     arrow = arrow(angle = 20, length = unit(5, "mm")), colour = "tomato") +
        geom_label_repel(data = rot, aes(x = scaler*rot[,x], y = scaler*rot[,y], label = rownames(rot)),
                         label.size = 0.2, inherit.aes = FALSE)
    }
    p2 <- p2 +
      scale_colour_manual(values = cols) +
      scale_fill_manual(values = cols) +
      coord_fixed() +
      labs(x = paste0("PCo ", x, " (", dvar[x], "%)"), y = paste0("PCo ", y, " (", dvar[y], "%)"), fill = "", colour = "") +
      theme_bw() +
      ggtitle(paste0("PCoA of ", title, " (nind = ", nrow(d.pcoa), " ; nvar = ", ncol(d.pcoa), ")"))
    print(p2)
  }
    
  
  ## Return results¨
  res <- list(df = df.sub, df.imp = d.imp, classes = classes)
  if (pca) {
    res[["df.pca"]] <- d.pca
    res[["res.pca"]] <- d.pca.res
    res[["plot.pca"]] <- p1
  }
  if (pcoa) {
    res[["df.pcoa"]] <- d.pcoa
    res[["res.pcoa"]] <- d.pcoa.res
    res[["type.pcoa"]] <- type
    res[["plot.pcoa"]] <- p2
  }
  invisible(res)
}

# subset.data
subset.data <- function(df, grep.keep = NULL, grep.rm = NULL, var.full = NULL, 
                        na.strings = c(NA,"NA",""), frac.ind = 0.75, frac.var = 0.75, 
                        tableplot = TRUE, verbose = TRUE) {
  # helperfunctions
  strip.constant <- function(df, verbose = FALSE) {
    kept <- apply(df, 2, function(x) {length(unique(x[!x%in%na.strings])) > 1})
    if ((length(kept) < ncol(df)) & verbose) {
      message("removed the following ", ncol(df)-length(kept), " constant variables:\n", paste(colnames(df)[!colnames(df) %in% kept], collapse = ","))
    }
    return(df[,kept])
  }
  if (tableplot) library(tabplot)
  
  # keep variables in grep.keep
  if (!is.null(grep.keep)) {
    dd <- df[,grep(paste0(grep.keep, collapse = "|"), colnames(df))]
    if ((ncol(dd) < ncol(df)) & verbose) {
      message("retained the following ", ncol(dd), " variables:\n", paste(colnames(dd), collapse = ","))
    }
  } else {
    dd <- df
  }
  
  # remove variables in grep.rm
  if (!is.null(grep.rm)) {
    dd2 <- dd[,grep(paste0(grep.rm, collapse = "|"), colnames(dd), invert = TRUE)]
    if ((ncol(dd2) < ncol(dd)) & verbose) {
      if ((ncol(dd) < ncol(df)) & verbose) {
        message("removed the following ", ncol(dd) - ncol(dd2), " variables:\n", paste(colnames(dd)[!colnames(dd) %in% colnames(dd2)], collapse = ","))
      }
    }
    dd <- dd2 ; rm(dd2)
  }
  
  # remove constant (or all-missing) variables
  dd <- strip.constant(dd, verbose = verbose)
  if (tableplot) t1 <- tableplot(dat = dd, nBins = min(c(100, nrow(dd))))
  
  # remove individuals with missingness in var.full
  if (!is.null(var.full)) {
    dd <- dd[! dd[,var.full] %in% na.strings,]
  }
  
  # remove individuals with high missingness
  comp.ind <- apply(dd, 1, function(x) {sum(!x %in% na.strings)/length(x)})
  dd <- dd[comp.ind >= frac.ind,]
  
  comp.var <- apply(dd, 2, function(x) {sum(!x %in% na.strings)/length(x)})
  dd <- dd[,comp.var >= frac.var]
  
  # finish
  if (tableplot) t2 <- tableplot(dat = dd, nBins = min(c(100, nrow(dd))))
  return(list(df = dd, p1 = t1, p2 = t2))
}

#####################################################################################################

### Load data
# vegetative traits
# load(p.leaves) # d.morph
load(p.morpho) # dd


### Fix data
# fix column names
names(dd) <- gsub(" ","",gsub("__","_",gsub(",","_",gsub("\n","_",gsub("\r","_",names(dd))))))

# fix species
dd$Species <- gsub("^maritima$", "maritima_maritima", dd$Species)

# remove rows with no identifyer
dd <- dd[-which(dd[,idvar] %in% na.strings),]
rownames(dd) <- dd[,idvar]

### Combine data
# # add shape
# dshape <- data.frame(readxl::read_excel(p.shape))
# dshape <- data.frame(apply(dshape[,!names(dshape) %in% c("Collection","Collection_momocs")], 2, numerize), row.names = dshape$Collection)
# dshape$Collection <- rownames(dshape)
# dd <- merge(dd[,!duplicated(colnames(dd))], dshape[,c("Collection",shpvars)], by = "Collection", all.x = TRUE, all.y = FALSE, sort = FALSE)
dd <- dd[,!duplicated(colnames(dd))] # if shape is not used
rownames(dd) <- dd[,idvar]


### Subset data
# select species
dd <- subset(dd, Species %in% spec)

# select variables
dd <- dd[,grep("\\(", colnames(dd), invert = TRUE)] # removes variables with multiple measurements

# create ordered factors
# ordvar <- names(dd)[grep("^PUB_|^MAR_", names(dd))]
# for (i in ordvar) dd[,i] <- ordered(dd[,i], levels = sort(unique(dd[,i])))

# make numeric (Quantitative Continuous and Discrete and Binary, means)
numvar <- names(dd)[grep("^LEN_|^WID_|^RATIO_|^NB_|^CORI_|^DENS_|^VEN_high_low|^MAR_|^PUB_", names(dd))]
for (i in numvar) dd[,i] <- suppressWarnings(as.numeric(as.character(dd[,i])))

# fix NA strings
for (i in seq(ncol(dd))) {
  dd[dd[,i] %in% na.strings,i] <- NA
  if (is.factor(dd[,i])) dd[,i] <- droplevels(dd[,i])
}


##############
### FLOWER ###
##############

# var.flower <- c("infl","flower","calyx","standard","wing","keel","pedicel","gyn","LEN_style") # 26 variables
# var.flower <- c("^TYPE_infl", "PUB_infl_axis", "LEN_pedicel", "^LEN_calyx$","^LEN_calyx_low$", "LEN_standard", "WID_standard", "RATIO_standard", "LEN_wing","WID_wing", "^PUB_gyn", "LEN_style") # 13 variables
var.flower <- c("^TYPE_infl", "PUB_infl_axis", "LEN_pedicel", "^LEN_calyx$","^LEN_calyx_low$", "LEN_standard", "WID_standard", "RATIO_standard", "LEN_wing","WID_wing", "^PUB_gyn") # 12 variables

var.flower.excl <- c("claw","NB_infl_units","LEN_infl_comp")
# var.flower.full <- "LEN_style"
# var.flower.fact <- c("TYPE_infl","COL_gyn") # 26 variables
var.flower.fact <- "TYPE_infl"

ss.flower <- subset.data(df = dd, grep.keep = var.flower, grep.rm = var.flower.excl, frac.ind = 0.75)
dd.flower <- ss.flower$df
for (i in var.flower.fact) {dd.flower <- factor2dummy(df = dd.flower, factorname = i)}
dd.flower.spec <- factor(as.character(dd[rownames(dd.flower),classvar]))
dd.flower.cols <- cols.all[spec %in% dd.flower.spec]

pdf("DATASET_flower.pdf", width = 20, height = 20)
plot(ss.flower$p2)
graphics.off()

pdf("PCA_PCOA_flower.pdf", width = plot.size, height = plot.size)
res.flower <- ordplot(df = dd.flower, x = 1, y = 2, classes = dd.flower.spec, cols = dd.flower.cols, 
                      ordratio = ordratio, asymm = asymm,
                      pca = TRUE, pcoa = TRUE, hulls = FALSE, voronoi = TRUE, voronoi.scaling = 0.1,
                      title = "FLOWER characters")
graphics.off()


##############
### LEAVES ###
##############

# var.leaf <- c("leaf","leaflet","petiole","petiolule","rachis","VEN_high_low","^shp") # 17 variables
var.leaf <- c("PUB_petiolule","PUB_rachis","PUB_leaflet_low","PUB_leaflet_up","MAR_leaflet","NB_leaflets","LEN_leaf$","LEN_leaflet_", "WID_leaflet_", "RATIO_leaflet_dist", "LEN_petiole","LEN_petiolule", "CORI_leaflet") # 14 variables
var.leaf.excl <- NULL
# var.leaf.full <- "LEN_petiolule" #"shp1_orbicular"
var.leaf.fact <- NULL

ss.leaf <- subset.data(df = dd, grep.keep = var.leaf, grep.rm = var.leaf.excl, frac.ind = 0.65)
dd.leaf <- ss.leaf$df
for (i in var.leaf.fact) {dd.leaf <- factor2dummy(df = dd.leaf, factorname = i)}
dd.leaf.spec <- factor(as.character(dd[rownames(dd.leaf),classvar]))
dd.leaf.cols <- cols.all[spec %in% dd.leaf.spec]

pdf("DATASET_leaf.pdf", width = 20, height = 20)
plot(ss.leaf$p2)
graphics.off()

pdf("PCA_PCOA_leaf.pdf", width = plot.size, height = plot.size)
res.leaf <- ordplot(df = dd.leaf, x = 1, y = 2, classes = dd.leaf.spec, cols = dd.leaf.cols,
                    ordratio = ordratio, asymm = asymm,
                    pca = TRUE, pcoa = TRUE, hulls = FALSE, voronoi = TRUE, voronoi.scaling = 0.1,
                    title = "LEAF characters")
graphics.off()


######################
### LEAVES & FRUIT ###
######################

var.lfruit <- c("leaf","leaflet","petiole","petiolule","rachis","VEN_high_low","^shp", "fruit")
var.lfruit.excl <- NULL
# var.leaf.full <- "LEN_petiolule"
var.lfruit.fact <- "SHP_fruit_base"

ss.lfruit <- subset.data(df = dd, grep.keep = var.lfruit, grep.rm = var.lfruit.excl, frac.ind = 0.5, frac.var = 0.15)
dd.lfruit <- ss.lfruit$df
for (i in var.lfruit.fact) {dd.lfruit <- factor2dummy(df = dd.lfruit, factorname = i)}
dd.lfruit.spec <- factor(as.character(dd[rownames(dd.lfruit),classvar]))
dd.lfruit.cols <- cols.all[spec %in% dd.lfruit.spec]

table(!is.na(dd.lfruit$LEN_fruit_1), dd.lfruit.spec) # no fruit measurements for maritima_pubescens and chapelieri

pdf("DATASET_leaf_fruit.pdf", width = 20, height = 20)
plot(ss.lfruit$p2)
graphics.off()

pdf("PCA_PCOA_leaf_fruit.pdf", width = plot.size, height = plot.size)
res.lfruit <- ordplot(df = dd.lfruit, x = 1, y = 2, classes = dd.lfruit.spec, cols = dd.lfruit.cols,
                      ordratio = ordratio, asymm = asymm,
                      pca = TRUE, pcoa = TRUE, hulls = FALSE, voronoi = TRUE, voronoi.scaling = 0.1,
                      title = "LEAF & FRUIT characters")
graphics.off()


###############
### ECOLOGY ###
###############


### Subset rows
# load("../../Sampling/data.rda") # data: specimen data (with coordinates) downloaded from Tropicos.org
# dsub <- subset(data, Species %in% spec)[,c("Species","SpecimenID","DeterminationBy","LatitudeDecimal","LongitudeDecimal","Coordinate.Method")]
# dsub <- subset(dsub, Coordinate.Method %in% accurate.geo)
# dsub$ID_Lab <- rownames(dsub)
# dsub <- dsub[,c("ID_Lab","Species", names(dsub)[!names(dsub)%in%c("ID_Lab","Species")])]
# write.table(dsub, file = "specimens.txt", sep = "\t", row.names = FALSE, quote = FALSE)
dsub <- read.delim(p.data)
dsub <- dsub[grep(identifyers, dsub$DeterminationBy),]
dsub <- dsub[!is.na(dsub[,geovars[1]]),]
dsub <- dsub[!duplicated(dsub[,c(classvar,geovars)]),]

### Extract ecological variables
# spatial vector of sample coordinates
v.eco <- vect(x = dsub[,geovars], geom = geovars, type = "points", crs = "+init=epsg:4326")
v.eco <- terra::project(v.eco, r.longlat) # same projection as rasters (UTM)

# extract
deco <- dsub
for (ecovar in ecovars) {
  deco[,ecovar] <- extract(x = get(paste0("r.", ecovar)), y = v.eco, method = ifelse(ecovar %in% categorical, "simple", "bilinear"))[,2]
}


### Combine geology levels
legend.geo <- c("0 unmapped",
                "1 Carbonate [Marble (Cipolin)]",
                "2 Karst [Mesozoic Limestones incl Tsingy]",
                "3 NonCarbonate [Sandstones]",
                "4 Metasedimentary [Basement Rocks]",
                "5 Alkaline Intrusive Volcanic",
                "6 Silicic [Basement Rocks]",
                "7 Metaigneous [Basement Rocks]",
                "8 Ultramafic [Ultrabasics]",
                "9 Extrusive Volcanic [Lavas incl Basalts & Gabbros]",
                "11 Hydric Organic [Mangroves]",
                "14 Alluvium Fluvial [Alluvial & Lake deposits]",
                "15 Alluvium Beach Strand Coastal Dune",
                "18 Alluvium Other [Unconsolidated Sands]",
                "20 Water")
# table(deco$geol, useNA="a")
deco$geol[deco$geol %in% c(14,15,18)] <- 18


### Combine vegetation levels
legend.vege <- c("0 unmapped / Clouds",
                 "1 Water Bodies",
                 "2 Bare soil/rock",
                 "3 Mangroves",
                 "4 Cultivation",
                 "5 Western dry forest [NW/W]",
                 "6 Plateau grassland-wooded grassland mosaic",
                 "7 Wooded grassland-bushland mosaic",
                 "9 Western humid forest",
                 "10 Western dry forest [SW]",
                 "11 Degraded south western dry spiny forest",
                 "12 South western dry spiny forest-thicket",
                 "13 Wetlands",
                 "14 Humid forest",
                 "15 Littoral forest",
                 "16 Degraded humid Forest",
                 "18 unmapped / Clouds",
                 "19 South western coastal bushland",
                 "22 Western sub-humid forest",
                 "23 Tapia forest",
                 "25 Sea")
# table(deco$vege, useNA="a")
# deco$vege[deco$vege %in% c(5,10)] <- 10
deco$vege[deco$vege %in% c(6,7)] <- 7

### Create dummy variables for categorical variables
for (i in categorical) {
  tab.dummy <- sort(table(deco[,i], useNA = "always"), decreasing = TRUE)
  dummy <- names(tab.dummy[tab.dummy > sum(tab.dummy) * min.frac.dummy])
  deco <- factor2dummy(df = deco, factorname = i, keep.levels = dummy)
}


### Fix column names
colnames(deco) <- gsub("Decimal$", "", colnames(deco))


### PCA Biplot
var.eco <- c("^prec","^temp","^isothermality","^elevation","^aspect","^slope","^dist2","^geol","^vege","^Latitude","Longitude")
var.eco.excl <- NULL
# var.eco.full <- "Latitude"
var.eco.fact <- NULL

ss.eco <- subset.data(df = deco, grep.keep = var.eco, grep.rm = var.eco.excl)
dd.eco <- ss.eco$df
for (i in var.eco.fact) {dd.eco <- factor2dummy(df = dd.eco, factorname = i)}
dd.eco.spec <- factor(as.character(deco[rownames(dd.eco),classvar]))
dd.eco.cols <- cols.all[spec %in% dd.eco.spec]

pdf("DATASET_eco.pdf", width = 20, height = 20)
plot(ss.eco$p2)
graphics.off()

pdf("PCA_PCOA_eco.pdf", width = plot.size, height = plot.size)
res.eco <- ordplot(df = dd.eco, x = 1, y = 2, classes = dd.eco.spec, cols = dd.eco.cols,
                   ordratio = ordratio, asymm = asymm,
                   pca = TRUE, pcoa = TRUE, hulls = FALSE, voronoi = TRUE, voronoi.scaling = 0.1, biplot.quantile = 0.15,
                   title = "ECOLOGICAL variables")
graphics.off()

