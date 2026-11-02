# HortControl / PlantEye trait display names (User Manual §4).
# Shared by biology readout + DR loadings / biplot axes.
#
# !!! THIS CODE SHOULD BE CHECKED BY A HUMAN PLZ !!!
# Bin interval maps and polarity heuristics are easy to get wrong;
# verify against HortControl User Manual §4.3 before trusting biology text.

`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Interval string for a sanitized spectral-bin column stem.
#' e.g. gli_0_0_15 → "[0, 0.15]"; gli_0_2_0 → "[-0.2, 0]".
.sf_spectral_bin_interval <- function(base_name) {
  s <- tolower(as.character(base_name)[[1L]])
  s <- sub("^(gli|ndvi|psri|npci)_", "", s)
  s <- sub("_percent$", "", s)
  known <- c(
    "1_0" = "[-1, 0]",
    "1_0_2" = "[-1.0, -0.2]",
    "0_2_0" = "[-0.2, 0]",
    "0_0_15" = "[0, 0.15]",
    "0_15_0_3" = "[0.15, 0.3]",
    "0_3_0_45" = "[0.3, 0.45]",
    "0_3_0_6" = "[0.3, 0.6]",
    "0_45_0_6" = "[0.45, 0.6]",
    "0_6_1" = "[0.6, 1]"
  )
  if (s %in% names(known)) return(unname(known[[s]]))
  NA_character_
}

#' Morphologic / technical base names → HortControl-style titles.
.sf_hortcontrol_base_label <- function(base) {
  b <- tolower(as.character(base)[[1L]])
  known <- c(
    digital_biomass_mm3 = "Digital Biomass",
    digital_biomass = "Digital Biomass",
    x3d_leaf_area_mm2 = "3D Leaf Area",
    projected_leaf_area_mm2 = "Projected Leaf Area",
    plant_height_averaged_mm = "Plant Height Average",
    plant_height_max_mm = "Plant Height Max",
    canopy_light_penetration_depth_mm = "Canopy Light Penetration Depth",
    leaf_area_index = "Leaf Area Index",
    leaf_inclination = "Leaf Inclination",
    convex_hull_area_mm2 = "Convex Hull Area",
    convex_hull_area_coverage_percent = "Convex Hull Area Coverage",
    convex_hull_circumference_mm = "Convex Hull Circumference",
    convex_hull_maximum_width_mm = "Convex Hull Maximum Width",
    convex_hull_aspect_ratio = "Convex Hull Aspect Ratio",
    voxel_volume_total_mm3 = "Voxel Volume Total",
    surface_angle_average = "Surface Angle Average",
    hue_average = "Hue average",
    lightness_average_percent = "Lightness average",
    saturation_average_percent = "Saturation average",
    ndvi_average = "NDVI average",
    gli_average = "GLI average",
    psri_average = "PSRI average",
    npci_average = "NPCI average"
  )
  if (b %in% names(known)) return(unname(known[[b]]))

  idx <- if (grepl("^gli", b)) {
    "GLI"
  } else if (grepl("^ndvi", b)) {
    "NDVI"
  } else if (grepl("^psri", b)) {
    "PSRI"
  } else if (grepl("^npci", b)) {
    "NPCI"
  } else {
    NA_character_
  }
  if (!is.na(idx) && grepl("_percent$", b)) {
    iv <- .sf_spectral_bin_interval(b)
    if (!is.na(iv)) return(paste0(idx, " bin ", iv))
    return(paste0(idx, " bin"))
  }
  # Fallback: readable snake_case
  gsub("_", " ", b)
}

#' Full HortControl-style trait label, including optional GS metric suffix.
#' Examples: "GLI bin [0, 0.15] · auc", "Digital Biomass · last", "NDVI average".
.sf_hortcontrol_trait_label <- function(trait_name) {
  nm <- as.character(trait_name)[[1L]]
  base <- sub("__(auc|peak|growth_rate|last|mean)$", "", tolower(nm))
  metric <- if (grepl("__(auc|peak|growth_rate|last|mean)$", tolower(nm))) {
    sub("^.*__", "", tolower(nm))
  } else {
    NA_character_
  }
  lab <- .sf_hortcontrol_base_label(base)
  if (!is.na(metric)) paste0(lab, " · ", metric) else lab
}
