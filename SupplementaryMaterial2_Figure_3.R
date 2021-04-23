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
library(ConR) # EOO.computing()
library(alphahull) # EOO.computing()
library(tmap) # tm_shape()
library(concaveman) # geom_mark_hull()
library(rmapshaper) # tm_shape(simplify = 0.01)

### Parameters
# paths
p.data.eco <- "SpecimenData/specimens_343.txt" # includes precise geo-coordinates [not public]
# p.morph <- "SpecimenData/charactertable.xlsm" # includes repeated measurements on all organs and individuals [median values reported]

p.specimens <- "SpecimenData/SupplementaryMaterial1_343.txt"
p.morph.leaf <- "SpecimenData/SupplementaryMaterial_3.txt"
p.morph.flower <- "SpecimenData/SupplementaryMaterial_4.txt"
p.eco <- "SpecimenData/SupplementaryMaterial_5.txt"

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

for2017 <- raster("../../Ecology/MadagascarForestModel/RasterData/forest/for2017.tif")
eco2017 <- readRDS("../../Ecology/DataRastersVectors/Ecoregions/Madagascar_ecoregions2017.rds")
alt <- raster("../../Ecology/DataRastersVectors/DIVA-GIS/MDG_alt/MDG_alt.grd")
border <- readRDS("../../Ecology/DataRastersVectors/AdministrativeBoundaries/Madagascar_border.rds")
regions <- readRDS("../../Ecology/DataRastersVectors/AdministrativeBoundaries/MDG_regions.rds")
pas <- readRDS("../../Ecology/DataRastersVectors/ProtectedAreas/MDG_PA_202102_124.rds")

# categorical variables
categorical <- c("geol","vege")
min.frac.dummy <- 0.075 # minimum level frequency: at least this fraction of observations need to be at this categorical level for that level to be kept
na.strings <- c("NA","",NA)

# variable names
geovars <- c("LongitudeDecimal","LatitudeDecimal")
ecovars <- c("prec.annual","prec.seasonality","temp.annual","temp.seasonality", "isothermality","elevation","aspect","slope","dist2Coast","dist2InlandWaters","geol","vege")

idvar <- "Collection"
classvar <- "Species"
n.meta <- c(idvar,classvar)

# collection filter
accurate.geo <- c("GPS","google earth","map","post facto, high", "corrected post facto")
identifyers <- "Crameri|Phillipson|Wilding"

# classes
spec.morpho <- c("chapelieri","aff._chapelieri_coriaceous","aff._chapelieri_coriaceous_large","aff._chapelieri_Manombo","aff._chapelieri_N","aff._chapelieri_S","cf._chapelieri", # chapelieri sensu lato
                 "louvelii","aff._louvelii_Mananara","aff._louvelii_PAL", # louvelii sensu lato
                 "maritima_pubescens", "maritima_maritima",
                 "pseudomaritima", # southeastern paniculate
                 "razakamalalae") # southeastern racemose
s.lat <- c("chapelieri","louvelii") # will be merged to form sensu lato groups

spec.map <- c(spec.morpho,
              "racemosa") # northeastern racemose = SAVA material

# class colour
col.pal <- c("#A6CEE3","#A99099","#4F9F3B","#E93E3F","#FDAC4F","#B15928","#FFFF99")

# plotting
plot.size <- 10

# pcoa types (binary values would be handled automatically by cluster::daisy)
ordratio = "^IND_|^MAR_leaflet" # grep ordinal factors
asymm = "^TEX_|^TYPE_infl|^VEGE_|^LITH_" # grep binary nominal factors (dummy variables)


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

# add N
addN <- function(x) {
  t <- table(x)
  n <- as.character(x)
  for (i in names(t)) {
    n[x == i] <- paste0(n[x == i], " (n=", t[i],")")
  }
  if (is.factor(x)) {
    levs <- n[sapply(levels(x), function(z) {which(as.character(x)==z)[1]})]
    levs[is.na(levs)] <- paste0(levels(x)[is.na(levs)], " (n=0)")
    n <- factor(n, levels = levs)
  }
  return(n)
}

# factor2dummy
factor2dummy <- function(df, factorname, keep.levels = levels(factor(df[,factorname])), na.method = "zero") {
  
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
ordplot <- function(df, impute.fun = "class.mean", classes = factor(rep("", nrow(df))),
                    x = 1, y = 2, pca.flipx = 1, pca.flipy = 1, pcoa.flipx = 1, pcoa.flipy = 1,
                    cols = adegenet::funky(nlevels(droplevels(classes))),
                    scale = TRUE, center = TRUE, pca = TRUE, pcoa = TRUE, labels = FALSE, label.size.ind = 5,
                    hulls = TRUE, concavity = 2, radius = unit(2.5, "mm"), expand = unit(5, "mm"),
                    voronoi = FALSE, voronoi.scaling = 0.1,
                    biplot = TRUE, label.size.var = 5, biplot.quantile = 0, biplot.scaling = 0.8, 
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
    d.imp2 <- impute(df = d.imp$imputed, fun = "mean")
    d.imp2$total <- d.imp$total
    d.imp2$per.variable <- d.imp$per.variable
    d.imp2$per.individual <- d.imp$per.individual
    d.imp2$imputation[!is.na(d.imp$imputation)] <- d.imp$imputation[!is.na(d.imp$imputation)]
    d.imp <- d.imp2 ; rm(d.imp2)
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
    
    # axis flipping
    res.pca[,paste0("PC", x)] <- res.pca[,paste0("PC", x)] * pca.flipx
    res.pca[,paste0("PC", y)] <- res.pca[,paste0("PC", y)] * pca.flipy
    rot[,x] <- rot[,x] * pca.flipx
    rot[,y] <- rot[,y] * pca.flipy
    
    # biplot
    p1 <- ggplot(data = res.pca, aes_string(x = paste0("PC", x), y = paste0("PC", y), colour = "Class"))
    if (hulls) {
      p1 <- p1 +
        # geom_mark_rect(aes(fill = Class), radius = radius, expand = expand)
        # geom_mark_ellipse(aes(fill = Class), radius = radius, expand = expand)
        geom_mark_hull(aes(fill = Class), concavity = concavity, radius = radius, expand = expand)
      
    }
    if (voronoi) {
      p1 <- p1 +
        geom_voronoi_tile(aes(fill = Class, group = -1L), 
                          max.radius = diff(range(res.pca[,x]))*voronoi.scaling, colour = NA,
                          alpha = 0.5, show.legend = TRUE)
    }
    p1 <- p1 +
      geom_point()
    if (labels) {
      p1 <- p1 +
        geom_label_repel(label = rownames(res.pca), size = label.size.ind)
    }
    if (biplot) {
      p1 <- p1 +
        geom_segment(data = rot, aes(x = rep(0, nrow(rot)), y = rep(0, nrow(rot)),
                                     xend = scaler*rot[,x], yend = scaler*rot[,y]),
                     arrow = arrow(angle = 20, length = unit(5, "mm")), colour = "tomato") +
        geom_label_repel(data = rot, aes(x = scaler*rot[,x], y = scaler*rot[,y], label = rownames(rot)),
                         size = label.size.var, inherit.aes = FALSE)
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
    colnames(rot) <- c(x,y)
    
    # eigenvalues / % variance
    dvar <- round(100*(d.pcoa.res$values$Eigenvalues / sum(d.pcoa.res$values$Eigenvalues)), 2)
    
    # biplot scaling
    scaler <- min(max(abs(pr.coo[, x]))/max(abs(rot[,1])),
                  max(abs(pr.coo[, y]))/max(abs(rot[,2]))) * biplot.scaling
    
    # biplot filtering
    qt <- apply(abs(rot), 2, quantile, biplot.quantile, na.rm = TRUE)
    torm <- abs(rot[,1]) < qt[1] & abs(rot[,2]) < qt[2]
    if ((sum(torm) > 0) & verbose) {
      message("removed the following ", sum(torm), " variables from PCoA biplot:\n", paste(names(torm[torm]), collapse = ","))
    }
    rot <- data.frame(rot[!torm,])
    
    # axis flipping
    res.pcoa[,paste0("Axis.", x)] <- res.pcoa[,paste0("Axis.", x)] * pcoa.flipx
    res.pcoa[,paste0("Axis.", y)] <- res.pcoa[,paste0("Axis.", y)] * pcoa.flipy
    rot[,1] <- rot[,1] * pcoa.flipx
    rot[,2] <- rot[,2] * pcoa.flipy
    
    # biplot
    # biplot(d.pcoa.res, Y = Y, rn = rep("o", nrow(d.pcoa.res$vectors))) # compare with biplot()
    p2 <- ggplot(data = res.pcoa, aes_string(x = paste0("Axis.", x), y = paste0("Axis.", y), colour = "Class"))
    if (hulls) {
      p2 <- p2 +
        # geom_mark_rect(aes(fill = Class), radius = radius, expand = expand)
        # geom_mark_ellipse(aes(fill = Class), radius = radius, expand = expand)
        geom_mark_hull(aes(fill = Class), concavity = concavity, radius = radius, expand = expand)
    }
    if (voronoi) {
      p2 <- p2 +
        geom_voronoi_tile(aes(fill = Class, group = -1L), 
                          max.radius = diff(range(res.pcoa[,x]))*voronoi.scaling, colour = NA,
                          alpha = 0.5, show.legend = TRUE)
    }
    p2 <- p2 +
      geom_point()
    if (labels) {
      p2 <- p2 +
        geom_label_repel(label = rownames(res.pcoa), size = label.size.ind)
    }
    if (biplot) {
      p2 <- p2 +
        geom_segment(data = rot, aes(x = rep(0, nrow(rot)), y = rep(0, nrow(rot)),
                                     xend = scaler*rot[,1], yend = scaler*rot[,2]),
                     arrow = arrow(angle = 20, length = unit(5, "mm")), colour = "tomato") +
        geom_label_repel(data = rot, aes(x = scaler*rot[,1], y = scaler*rot[,2], label = rownames(rot)),
                         size = label.size.var, inherit.aes = FALSE)
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
    
  
  ## Return results
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

# ### Load data (including precise [not public] coordinates)
# ## vegetative and fertile traits (measured in P or on scans)
# dd.all <- data.frame(readxl::read_excel(p.morph, na = c("", "NA")))
# 
# # subset rows
# dd <- dd.all[!is.na(dd.all$SpecimenID),]
# 
# # subset columns
# dd <- dd[,1:grep("^Notes$", names(dd))]
# 
# # fix col names
# names(dd) <- gsub("[.][.][.][.]|[.][.][.]|[.][.]", ".", names(dd))
# 
# # fix row names
# rownames(dd) <- dd[,idvar]
# 
# # numerize
# isnum <- apply(dd, 2, function(x) {sum(!is.na(suppressWarnings(as.numeric(x)))) > 1 & length(grep(",", x)) == 0})
# dd[,isnum] <- sapply(dd[,isnum], as.numeric)
# 
# # min, max, mean, median
# numstrings <- names(dd)[grep("^LEN_|^WID_|^NB_", names(dd))-1] # identify cells with comma-separated repeated measurements
# numstrings <- numstrings[grep("^TYPE_", numstrings, invert = TRUE)]
# 
# for (numstring in numstrings) {
# 
#   basename <- names(dd)[which(names(dd) == numstring)+1]
# 
#   ls <- sapply(strsplit(as.character(dd[,numstring]), split = ","), as.numeric)
#   ls.min <- as.numeric(suppressWarnings(sapply(ls, FUN = min)))
#   ls.max <- as.numeric(suppressWarnings(sapply(ls, FUN = max)))
#   ls.mean <- as.numeric(suppressWarnings(sapply(ls, FUN = mean)))
#   ls.median <- as.numeric(suppressWarnings(sapply(ls, FUN = median)))
#   # ls.sd <- as.numeric(suppressWarnings(sapply(ls, FUN = sd)))
# 
#   dd[,paste0(basename, "_min")] <- ls.min
#   dd[,paste0(basename, "_max")] <- ls.max
#   dd[,paste0(basename, "_mean")] <- ls.mean
#   dd[,paste0(basename, "_median")] <- ls.median
#   # dd[,paste0(basename, "_sd")] <- ls.sd
#   dd[,basename] <- NULL
# }
# 
# ### Subset data
# # select species
# dd <- subset(dd, Species %in% spec.morpho)
# 
# # combine sensu lato
# sp.morpho <- spec.morpho
# sp.map <- spec.map
# for (i in s.lat) {
#   dd$Species[grep(i, dd$Species)] <- paste(i, "s.lat.")
#   sp.morpho[grep(i, sp.morpho)] <- paste(i, "s.lat.")
#   sp.map[grep(i, sp.map)] <- paste(i, "s.lat.")
# }
# sp.morpho <- sort(unique(paste0("D. ", gsub("_", " subsp. ", sp.morpho))))
# sp.map[grep("racemosa", sp.map)] <- "sp. (SAVA material)"
# sp.map <- sort(unique(paste0("D. ", gsub("_", " subsp. ", sp.map))))

scan()

###################
### LEAVES PCOA ###
###################

# var.leaf <- c("LEN_petiole","LEN_leaf_","NB_leaflets","IND_rachis",
#               "IND_petiolule","IND_leaflet",
#               "LEN_petiolule",
#               # "LEN_leaflet_prox","WID_leaflet_prox",
#               "LEN_leaflet_dist","WID_leaflet_dist",
#               "RATIO_leaflet_dist",
#               "TEX_leaflet","MAR_leaflet")
# var.leaf.excl <- c("_min$","_mean$", "_max$", "_sd$")
# # var.leaf.full <- "LEN_petiolule" #"shp1_orbicular"
# 
# ss.leaf <- subset.data(df = dd, grep.keep = var.leaf, grep.rm = var.leaf.excl, frac.ind = 0.8) # frac.ind = 0.65
# message("kept the following ", ncol(ss.leaf$df), " variables:\n", paste(colnames(ss.leaf$df), collapse = "\n"))
# 
# pdf("DATASET_leaf.pdf", width = 20, height = 20)
# plot(ss.leaf$p2)
# graphics.off()
# 
# dd.leaf <- ss.leaf$df
# 
# dd.leaf.spec <- factor(addN(paste("D.", gsub("s.lat.", "s.l.", gsub("_", " ssp. ", as.character(dd[rownames(dd.leaf),classvar]))))))
# dd.leaf.cols <- col.pal[sp.map %in% sp.morpho]
# 
# dd.leaf$TAXON <- dd.leaf.spec
# 
# write.table(dd.leaf, file = p.morph.leaf, row.names = TRUE, quote = FALSE, sep = "\t")

dd.leaf <- read.delim(p.morph.leaf, row.names = 1)
dd.leaf.spec <- factor(dd.leaf$TAXON) ; dd.leaf <- dd.leaf[,!names(dd.leaf) %in% "COLOR"]
dd.leaf.cols <- c("#A6CEE3","#A99099","#4F9F3B","#E93E3F","#FDAC4F","#B15928")

var.leaf.fact <- NULL
for (i in var.leaf.fact) {dd.leaf <- factor2dummy(df = dd.leaf, factorname = i)}
names(dd.leaf) <- gsub("_median$", "", names(dd.leaf))


pdf("PCOA_leaf.pdf", width = plot.size, height = plot.size)
res.leaf <- ordplot(df = dd.leaf, x = 1, y = 2, classes = dd.leaf.spec, cols = dd.leaf.cols,
                    ordratio = ordratio, asymm = asymm, labels = F,
                    pca = FALSE, pcoa = TRUE, hulls = TRUE, concavity = 10, radius = unit(5, "mm"), expand = unit(5, "mm"), voronoi = FALSE,
                    title = "LEAF characters")
graphics.off()


scan()

###################
### FLOWER PCOA ###
###################

# var.flower <- c("TYPE_infl", "IND_infl_axis", "LEN_pedicel",
#                 "LEN_flower","LEN_pedicel",
#                 "LEN_calyx",
#                 "LEN_calyx_up",#"WID_calyx_up",
#                 # "LEN_calyx_mid","WID_calyx_mid",
#                 "LEN_calyx_low",#"WID_calyx_low",
#                 # "LEN_standard",
#                 "LEN_standard_lamina",
#                 # "LEN_standard_claw",
#                 "WID_standard",
#                 "RATIO_standard",
#                 "RATIO_flower",
#                 # "LEN_wing","LEN_wing_lamina","LEN_wing_claw","WID_wing",
#                 # "LEN_keel","LEN_keel_lamina","LEN_keel_claw","WID_keel",
#                 # "LEN_stamens","LEN_stamens_free_min","LEN_stamens_free_max",
#                 # "LEN_gyn"
#                 "IND_gyn"
#                 # "LEN_ovary",
#                 # "LEN_style"
#                 # "NB_ovules"
# )
# 
# var.flower.excl <- c("claw","NB_infl_units","LEN_infl_comp","LEN_calyx_mid","LEN_calyx_low", "WID_calyx", "_min$","_mean$", "_max$", "_sd$")
# # var.flower.full <- "LEN_style"
# 
# ss.flower <- subset.data(df = dd[,grep(paste(var.flower.excl, collapse = "|"), names(dd), invert = TRUE)], grep.keep = var.flower, grep.rm = var.flower.excl, frac.var = 0.1, frac.ind = 0.75)
# message("kept the following ", ncol(ss.flower$df), " variables:\n", paste(colnames(ss.flower$df), collapse = "\n"))
# 
# pdf("DATASET_flower.pdf", width = 20, height = 20)
# plot(ss.flower$p2) # Ramamonjiarisoa 4 has missing measurements -> ask photograph of dissection?
# graphics.off()
# 
# dd.flower <- ss.flower$df
# 
# dd.flower.spec <- factor(addN(paste("D.", gsub("s.lat.", "s.l.", gsub("_", " ssp. ", as.character(dd[rownames(dd.flower),classvar]))))))
# dd.flower.cols <- col.pal[sp.map %in% sp.morpho]
# 
# dd.flower$TAXON <- dd.flower.spec
# 
# write.table(dd.flower, file = p.morph.flower, row.names = TRUE, quote = FALSE, sep = "\t")

dd.flower <- read.delim(p.morph.flower, row.names = 1)
dd.flower.spec <- factor(dd.flower$TAXON) ; dd.flower <- dd.flower[,!names(dd.flower) %in% "COLOR"]
dd.flower.cols <- c("#A6CEE3","#A99099","#4F9F3B","#E93E3F","#FDAC4F","#B15928")

var.flower.fact <- "TYPE_infl"
for (i in var.flower.fact) {dd.flower <- factor2dummy(df = dd.flower, factorname = i)}
names(dd.flower) <- gsub("_median$", "", names(dd.flower))
names(dd.flower) <- gsub("_lamina$", "", names(dd.flower))

pdf("PCOA_flower.pdf", width = plot.size, height = plot.size)
res.flower <- ordplot(df = dd.flower, x = 1, y = 2, classes = dd.flower.spec, cols = dd.flower.cols,
                      ordratio = ordratio, asymm = asymm, labels = F,
                      pca = FALSE, pcoa = TRUE, hulls = TRUE, concavity = 10, radius = unit(5, "mm"), expand = unit(5, "mm"), voronoi = FALSE,
                      title = "FLOWER characters")
graphics.off()


scan()


###################
### ECOLOGY PCA ###
###################

# ### Read tropicos data
# load("../../Sampling/data.rda") # data: specimen data (with coordinates) downloaded from Tropicos.org
# dsub <- data[!duplicated(data$Collection),c("Collection","ID_Lab","SpecimenID","Genus","Species","CurrentDetermination","DeterminationBy","LatitudeDecimal","LongitudeDecimal","Minimum.Year", "Coordinate.Method")] # remove duplicated entries
# dsub <- subset(dsub, !SpecimenID %in% c("NA",NA,"")) # remove entries not on tropicos
# dsub <- subset(dsub, !LatitudeDecimal %in% c("NA",NA,"")) # remove entries without any coordinates
# 
# # select taxa
# dsub <- subset(dsub, Genus == "Dalbergia" & Species %in% spec.map)
# 
# # handle sensu lato
# for (i in s.lat) {
#   dsub$Species[grep(i, dsub$Species)] <- paste(i, "s.lat.")
# }
# 
# # handle undescribed SAVA material (D. racemosa ined.)
# dsub$Species[dsub$Species %in% "racemosa"] <- "sp. (SAVA material)"
# 
# dsub <- dsub[!duplicated(dsub[,c(classvar,geovars)]),] # remove conspecific specimens georeferenced to the same coordinates
# dim(dsub) # 343, with precise coordinates
# write.table(dsub, file = paste0("SpecimenData/specimens_", nrow(dsub), ".txt"), sep = "\t", row.names = FALSE, quote = FALSE)


# ### Read tropicos data
# dsub343 <- read.delim(p.data.eco) # 343, with precise coordinates
# 
# ### Subset tropicos data
# # select specimens with accurate georeferencing
# dsub <- subset(dsub343, !is.na(LatitudeDecimal) & Coordinate.Method %in% accurate.geo)
# 
# # select specimens with identification by authors
# # dsub <- dsub[grep(identifyers, dsub$DeterminationBy),] # select collections identified by authors
# 
# # filter out old collections (< 2000)
# dsub <- subset(dsub, Minimum.Year >= 2000)
# 
# # reorder columns
# dsub <- dsub[,c("ID_Lab","Species", names(dsub)[!names(dsub)%in%c("ID_Lab","Species")])]
# rownames(dsub) <- dsub[,idvar]
# dim(dsub) # 263
# 
# ### Extract ecological variables
# # spatial vector of sample coordinates
# v.eco <- vect(x = dsub[,geovars], geom = geovars, type = "points", crs = "+init=epsg:4326")
# v.eco <- terra::project(v.eco, r.longlat) # same projection as rasters (UTM)
# 
# # extract
# deco <- dsub
# for (ecovar in ecovars) {
#   deco[,ecovar] <- terra::extract(x = get(paste0("r.", ecovar)), y = v.eco, method = ifelse(ecovar %in% categorical, "simple", "bilinear"))[,2]
# }
# 
# 
# ### Combine geology levels
# legend.geo <- c("0 unmapped",
#                 "1 Carbonate [Marble (Cipolin)]",
#                 "2 Karst [Mesozoic Limestones incl Tsingy]",
#                 "3 NonCarbonate [Sandstones]",
#                 "4 Metasedimentary [Basement Rocks]",
#                 "5 Alkaline Intrusive Volcanic",
#                 "6 Silicic [Basement Rocks]",
#                 "7 Metaigneous [Basement Rocks]",
#                 "8 Ultramafic [Ultrabasics]",
#                 "9 Extrusive Volcanic [Lavas incl Basalts & Gabbros]",
#                 "11 Hydric Organic [Mangroves]",
#                 "14 Alluvium Fluvial [Alluvial & Lake deposits]",
#                 "15 Alluvium Beach Strand Coastal Dune",
#                 "18 Alluvium Other [Unconsolidated Sands]",
#                 "20 Water")
# # table(deco$geol, useNA="a")
# deco$geol[deco$geol %in% c(14,15,18)] <- 18
# 
# 
# ### Combine vegetation levels
# legend.vege <- c("0 unmapped / Clouds",
#                  "1 Water Bodies",
#                  "2 Bare soil/rock",
#                  "3 Mangroves",
#                  "4 Cultivation",
#                  "5 Western dry forest [NW/W]",
#                  "6 Plateau grassland-wooded grassland mosaic",
#                  "7 Wooded grassland-bushland mosaic",
#                  "9 Western humid forest",
#                  "10 Western dry forest [SW]",
#                  "11 Degraded south western dry spiny forest",
#                  "12 South western dry spiny forest-thicket",
#                  "13 Wetlands",
#                  "14 Humid forest",
#                  "15 Littoral forest",
#                  "16 Degraded humid Forest",
#                  "18 unmapped / Clouds",
#                  "19 South western coastal bushland",
#                  "22 Western sub-humid forest",
#                  "23 Tapia forest",
#                  "25 Sea")
# # table(deco$vege, useNA="a")
# # deco$vege[deco$vege %in% c(5,10)] <- 10
# deco$vege[deco$vege %in% c(6,7)] <- 7
# 
# ### Save Supplementary Material 1
# dd.SM1 <- merge(dsub343, deco, by = intersect(names(dsub343), names(deco)), all.x = TRUE, all.y = TRUE, sort = FALSE)
# names(dd.SM1) <- gsub("temp.", "TEMP_", gsub("prec.", "PRECIP_", gsub("^geol$", "LITHology", gsub("^vege$", "VEGEtation", names(dd.SM1)))))
# names(dd.SM1) <- gsub("annual", "Annual", gsub("seasonality", "Seasonality", gsub("isothermality", "Isothermality", gsub("elevation", "Elevation", gsub("slope", "Slope", gsub("aspect", "Aspect", gsub("dist2","DIST_", names(dd.SM1))))))))
# dd.SM1 <- dd.SM1[,!names(dd.SM1) %in% c("ID_Lab","Aspect","DIST_InlandWaters","PRECIP_Seasonality")]
# dd.SM1 <- dd.SM1[order(dd.SM1$Species, dd.SM1$Collection),]
# for (i in c("LatitudeDecimal","LongitudeDecimal","Elevation")) {dd.SM1[,i] <- scale(dd.SM1[,i])} # standardize geo-coordinates
# write.table(dd.SM1[,c("SpecimenID","Collection","Genus","Species","CurrentDetermination","DeterminationBy","Minimum.Year","Coordinate.Method")], file = p.specimens, sep = "\t", quote = FALSE, row.names = FALSE)
# 
# 
# ### Create dummy variables for categorical variables
# for (i in categorical) {
#   tab.dummy <- sort(table(deco[,i], useNA = "always"), decreasing = TRUE)
#   dummy <- names(tab.dummy[tab.dummy > sum(tab.dummy) * min.frac.dummy])
#   deco <- factor2dummy(df = deco, factorname = i, keep.levels = dummy)
# }
# 
# 
# ### Fix column names
# colnames(deco) <- gsub("geol_4", "LITH_Metasedimentary", gsub("geol_6", "LITH_Silicic", gsub("geol_7", "LITH_Metaigneous", gsub("geol_9", "LITH_Volcanic", gsub("geol_18", "LITH_Alluvium",
#                    gsub("vege_14", "VEGE_Humid", gsub("vege_15", "VEGE_Littoral", gsub("vege_16", "VEGE_DegradedHumid",
#                     gsub("dist2Coast", "DIST_Coast", gsub("dist2InlandWaters", "DIST_InlandWater",
#                      gsub("elevation", "Elevation", gsub("aspect", "Aspect", gsub("slope", "Slope",
#                       gsub("isothermality", "Isothermality", gsub("temp.annual", "TEMP_Annual", gsub("temp.seasonality", "TEMP_Seasonality",
#                        gsub("prec.annual", "PRECIP_Annual", gsub("prec.seasonality", "PRECIP_Seasonality", gsub("Decimal$", "", colnames(deco))))))))))))))))))))
# 
# ### Fix row names
# rownames(deco) <- deco$Collection
# 
# ### PCA Biplot
# var.eco <- c("Latitude","Longitude","Elevation","Slope","^DIST_Coast","^LITH","^VEGE","^PRECIP_Annual","^TEMP","^Iso")
# 
# var.eco.excl <- NULL
# # var.eco.full <- "Latitude"
# dd.eco.spec.excl <- "sp. (SAVA material)"
# 
# ss.eco <- subset.data(df = deco, grep.keep = var.eco, grep.rm = var.eco.excl)
# message("kept the following ", ncol(ss.eco$df), " variables:\n", paste(colnames(ss.eco$df), collapse = "\n"))
# 
# pdf("DATASET_eco.pdf", width = 20, height = 20)
# plot(ss.eco$p2)
# graphics.off()
# 
# dd.eco <- ss.eco$df
# dd.eco <- dd.eco[!deco$Species %in% dd.eco.spec.excl,] # exclude SAVA population
# 
# # scale geo-coordinates
# for (i in c("Latitude","Longitude","Elevation")) {dd.eco[,i] <- scale(dd.eco[,i])} # standardize geo-coordinates
# dd.eco$TAXON <- factor(addN(paste("D.", gsub("s.lat.", "s.l.", gsub("_", " ssp. ", as.character(deco[rownames(dd.eco),classvar]))))))
# write.table(dd.eco, file = p.eco, sep = "\t", quote = FALSE, row.names = TRUE) # 257, scaled coordinates

dd.eco <- read.delim(p.eco, row.names = 1) # 257, scaled coordinates
dd.eco.spec <- factor(dd.eco$TAXON) ; dd.eco <- dd.eco[,!names(dd.eco) %in% "TAXON"]
dd.eco.cols <- col.pal

var.eco.fact <- NULL
for (i in var.eco.fact) {dd.eco <- factor2dummy(df = dd.eco, factorname = i)}

pdf("PCA_PCOA_eco.pdf", width = plot.size, height = plot.size)
res.eco <- ordplot(df = dd.eco, x = 1, y = 2, classes = dd.eco.spec, cols = dd.eco.cols,
                   ordratio = ordratio, asymm = asymm,
                   pca = TRUE, pcoa = FALSE, hulls = TRUE, concavity = 1E06, 
                   voronoi = FALSE, # biplot.quantile = 0.15,
                   pca.flipx = -1, pca.flipy = -1, #biplot.scaling = 1.5,
                   title = "ECOLOGICAL variables")
graphics.off()

scan()


########################
### Distribution map ###
########################

# read data
p.data.map <- p.data.eco
dsub.map <- read.delim(p.data.map) # 343, with precise coordinates (not public)

# filter out doubtful specimens
dsub.map <- subset(dsub.map, !Collection %in% c("Service Forestier Madagascar 38-R-118 SF"))
dim(dsub.map) # 342

# calculate EOO
sp.map <- SpatialPointsDataFrame(coords = dsub.map[,geovars], data = dsub.map, 
                                 proj4string = CRS("+proj=utm +zone=38 +south +datum=WGS84 +units=m +no_defs"))
sp.map.spec <- factor(addN(paste("D.", gsub("s.lat.", "s.l.", gsub("_", " ssp. ", as.character(dsub.map[,classvar]))))))

# relevel
sp.order <- c(1,2,4,3,6,5,7)
sp.map.spec <- factor(sp.map.spec, levels = levels(sp.map.spec)[sp.order]) # relevel such that small hulls are on top
sp.map.spec.lev <- levels(factor(dsub.map[,classvar]))[sp.order]

if (exists("sp.eoo")) rm(sp.eoo)
eoo <- list()
for (i in sp.map.spec.lev) {
  cat("processing", i, "\n")
  
  switch(i,
         "chapelieri s.lat." = {alpha = 5; buff.alpha = 300/3600}, 
         "louvelii s.lat." = {alpha = 1; buff.alpha = 300/3600}, # alpha = 10
         "maritima_maritima" = {alpha = 10; buff.alpha = 300/3600},
         "maritima_pubescens" = {alpha = 10; buff.alpha = 300/3600},
         "pseudomaritima" = {alpha = 10; buff.alpha = 100/3600}, # 100 arc seconds ~= 3 km
         "razakamalalae" = {alpha = 10; buff.alpha = 300/3600}, # 300 arc seconds ~= 9 km
         "sp. (SAVA material)" = {alpha = 10; buff.alpha = 300/3600},
         )
  
  repeat{
    d.sp.map <- subset(sp.map, Species %in% c(i))
    d.eoo <- try(EOO.computing(XY = data.frame(ddlat = d.sp.map$LatitudeDecimal, ddlon = d.sp.map$LongitudeDecimal, tax = i),
                              exclude.area = TRUE, country_map = border, export_shp = TRUE, write_shp = FALSE, 
                              alpha = alpha, buff.alpha = buff.alpha, method.range = "alpha.hull", write_results = FALSE), silent = TRUE)
    if(!inherits(d.eoo, "try-error")) break()
  }
  eoo <- c(eoo, d.eoo)
  sp <- d.eoo[[paste("spatial.polygon_1")]]
  crs(sp) <- "+proj=longlat +ellps=WGS84 +no_defs"
  
  sp@polygons[[1]]@ID <- i
  if (!exists("sp.eoo")) sp.eoo <- sp else sp.eoo <- maptools::spRbind(sp.eoo, sp)
}
sp.eoo <- SpatialPolygonsDataFrame(sp.eoo, data = data.frame(Taxon = sp.map.spec.lev), match.ID=F)

# display tmap map
tmap_mode("view")
tmap_mode("plot")

# place names (coordinates taken from Google Maps)
dp <- data.frame("Place" = c("Antananarivo","Makiro-\nvana","Antalaha","Antongil\nBay","Makira","Masoala",
                             "Mahavelona","Toamasina","Betampona","Tampina",
                             "Manombo","Tolagnaro","Mandena","Sainte Luce","Tsitongambarika","Sahafina","Ambila-Lemaitso","Manakara","Analalava","Vangaindrano"),
                 "Latitude" = c(-18.57378, # Antananarivo
                                -14.12185, # Makirovana
                                -14.87895, # Antalaha
                                -15.80000, # Antongil
                                -15.40817, # Makira
                                -15.45000, # Masoala
                                -17.68017, # Mahavelona
                                -18.11662, # Toamasina
                                -17.88687, # Betampona
                                -18.52551, # Tampina
                                -22.99505, # Manombo
                                -25.00722, # Tolagnaro
                                -24.95071, # Mandena
                                -24.76165, # Sainte Luce
                                -24.524713, # Tsitongambarika
                                -18.81357, # Sahafina
                                -18.858472, # Ambila-Lemaitso
                                -22.143357, # Manakara
                                -17.701283, # Analalava
                                -23.350128), # Vangaindrano
                 "Longitude" = c(47.52856, # Antananarivo
                                 49.95169, # Makirovana
                                 50.27981, # Antalaha
                                 49.84036, # Antongil
                                 49.36581, # Makira
                                 50.20000, # Masoala
                                 49.51421, # Mahavelona
                                 49.39595, # Toamasina
                                 49.22591, # Betampona
                                 49.27662, # Tampina
                                 47.71408, # Manombo
                                 46.92000, # Tolagnaro
                                 47.02000, # Mandena
                                 47.13000, # Sainte Luce
                                 47.116828, # Tsitongambarika
                                 48.95837, # Sahafina
                                 49.145987, # Ambila-Lemaitso
                                 48.006673, # Manakara
                                 49.456564, # Analalava
                                 47.598809)) # Vangaindrano

# add
dp$Type <- c("Place")
dp$Type[dp$Place %in% c("Betampona","Sahafina","Masoala","Makira","Makiro-\nvana","Manombo","Mandena","Sainte Luce","Tsitongambarika","Analalava")] <- "Protected Area"
dp$Type <- factor(dp$Type)

# filter
dp <- subset(dp, !Place %in% c("Antalaha"))

# create spatial object
sp.place <- SpatialPointsDataFrame(coords = dp[,c(3,2)], data = dp, proj4string = crs(sp.eoo))
sp.place.lab <- sp.place

# adjust
addLong <- rep(0, nrow(dp))
addLat <- rep(0, nrow(dp))

addLong[dp$Place %in% c("Antananarivo")] <- -0.5
addLong[dp$Place %in% c("Masoala","Makira","Antongil\nBay")] <- -0.3
addLong[dp$Place %in% c("Betampona")] <- - 0.92
addLong[dp$Place %in% c("Sahafina")] <- - 0.75

addLat[dp$Place %in% c("Antananarivo","Tolagnaro")] <- c(-0.15,-0.1)
addLat[dp$Place %in% c("Makira","Masoala")] <- c(-0.15,-0.15)

sp.place.lab@coords[,1] <- sp.place.lab@coords[,1] + addLong
sp.place.lab@coords[,2] <- sp.place.lab@coords[,2] + addLat

# plotting aestethics
legend.show <- F #TRUE
shape.omit <- "Antongil\nBay" #c("Makira","Masoala","Antongil\nBay")

# plot
pdf("FIG3D.pdf", width = 12, height = 12)

tm_shape(alt, raster.downsample = T) +
  tm_raster(breaks = c(0, 100, 250, 500, 750, 1000, 1250, 1500, 2000, 2800), legend.show = legend.show,
            palette = gray.colors(10, rev = T), alpha = 0.5, title = "Elevation [m a.s.l.]") +
  tm_shape(for2017, raster.downsample = T) +
  tm_raster(palette = "forestgreen", alpha = 0.75, legend.show = legend.show,
            title = "Forest Cover (2017, Vieilledent et al. 2018b)") +
  tm_shape(shp = regions, simplify = 0.01) +
  tm_polygons(col = NA, alpha = 0, border.col = "black", border.alpha = 1, lwd = 1) +
  tm_grid(labels.inside.frame = TRUE, n.x = 7, n.y = 11) +  # SCALE
  tm_shape(shp = pas, simplify = 1) +
  tm_polygons(col = NA, alpha = 0, border.col = "#0096FF", border.alpha = 1, lwd = 2) +
  tm_shape(shp = sp.eoo, simplify = 1) + 
  tm_polygons(col = "Taxon", border.col = col.pal[sp.order], legend.show = legend.show, 
              palette = col.pal, lwd = 3, alpha = 0.6) +
  tm_shape(sp.place[!sp.place$Place %in% shape.omit,]) +
  tm_symbols(col = "red", border.col = "black", shape = "Type", shapes = c(21,22), size = 0.35, alpha = 0.65, legend.shape.show = legend.show) +
  tm_shape(sp.place.lab) +
  tm_text(text = "Place", size = 0.7, xmod = 0.5, just = "left") +
  tm_compass() + # type="radar", position=c("left", "top"), show.labels = 3
  tm_scale_bar(position=c("right", "bottom"), text.size = 0.7)
  
graphics.off()
