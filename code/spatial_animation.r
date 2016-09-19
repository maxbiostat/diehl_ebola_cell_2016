library(maptools)
library(ggplot2)
#######
Map <- readShapePoly("../data/gis/locations.shp")
Country_info <- read.csv("../data/location_country.csv", header = TRUE)
getCountry <- function(loc) as.character(Country_info$country[match(loc, Country_info$location)])
#######
plotGenotypeMap <- function(MutDat, title = "", index = NULL){
  freqV <- MutDat$V/(MutDat$A + MutDat$V + 1E-5)
  names(freqV) <- MutDat$location
  genotype_table <- data.frame(location = names(freqV),
                               genotype = ifelse(as.numeric(freqV) > .5, "V", "A"))
  ## fix to make it include both genotypes in the legend
  genotype_table$genotype <- factor(genotype_table$genotype, levels = c("A", "V"))
  ## end fix
  Map@data$location <- gsub("GrandCapeMo","GrandCapeMount", gsub("WesternRura", "WesternRural", gsub("WesternUrba", "WesternUrban", Map@data$location)))
  Map.f <- fortify(Map, region = "location")
  Map.f$location <- Map.f$id
  merge.shp <- merge(x = Map.f, y = genotype_table, by = "location", all = TRUE)
  merge.shp <- merge(merge.shp, Map@data)
  final.plot <- merge.shp[order(merge.shp$order), ]
  final.plot$country <- getCountry(final.plot$location)
  ##########################
  ##########################
  Map@data$country <- getCountry(Map@data$location)
  RedMap <- subset(Map, country == "GIN" | country == "LBR" | country == "SLE")
  CUnion <- maptools::unionSpatialPolygons(RedMap, RedMap@data$country, avoidUnaryUnion = TRUE)
  newDt <- data.frame(country = unique(RedMap@data$country))
  row.names(newDt) <- as.character(newDt$country)
  CUnion <- SpatialPolygonsDataFrame(CUnion, data = newDt)
  Map.U <- fortify(CUnion, region = "country")
  p0 <- ggplot() + 
    ggtitle(title) +
    geom_polygon(data = final.plot, 
                 aes(x = long, y = lat, group = group, colour = "black", fill = genotype),
                 size = .2) + 
    geom_polygon(data = Map.U, 
                 aes(x = long, y = lat, group = group, colour = "black"),
                 size = .7, alpha = 0.00001) + 
    scale_fill_manual(values = c( "A" = "black", "V" = "red"), drop = FALSE) +
    scale_color_manual(values = c("black", "black")) +
    labs(fill = "Predominant genotype") +
    guides(colour = FALSE) + 
    theme_bw() +
    coord_map()
  png(file = paste("../plots/GP82_animation/", LETTERS[index], "_", title, ".png", sep = ""), height = 4*480,
      width = 4*480, units = "px", res = 300)
  print(p0)
  dev.off()
}
####
library(zoo)
library(lubridate)
#
All.Data <- read.csv("../data/GP82_raw_data.csv", header = TRUE)
All.Data$date <- as.numeric(as.Date(All.Data$date))
NoAmbiguity <- subset(All.Data, location != "?" & GP82 != "X")
NoAmbiguity$location <- as.character(NoAmbiguity$location)
NoAmbiguity$GP82 <- as.character(NoAmbiguity$GP82)
summary(NoAmbiguity)
Windows <- as.Date(seq(ymd('2014-03-01'), ymd('2015-11-30'), by = '1 month'))
Titles <- as.yearmon(Windows)
ws <- as.numeric(Windows)
SplitData <- vector(length(ws)-1, mode = "list")
for(i in 2:length(ws)) SplitData[[i-1]] <- subset(NoAmbiguity, ws[i-1] < date & date < ws[i])
getMutperLoc <- function(dt){
  res <- xtabs(~ location + GP82, dt)
  if(ncol(res) == 1){
    missing <- setdiff(c("A", "V"), colnames(res))
    dt <- data.frame(rownames(res), as.vector(res), rep(0, nrow(res)))
    names(dt) <- c("location", colnames(res), missing)
  }else{
    dt <- data.frame(location = rownames(res), A = res[, "A"], V = res[, "V"])
  }
  return(dt)
} 
MutsPerLoc <- lapply(SplitData, getMutperLoc)
# plotGenotypeMap(MutsPerLoc[[13]], title = Titles[13], index = 13)
lapply(seq_along(MutsPerLoc), function(i) plotGenotypeMap(MutsPerLoc[[i]], title = Titles[i], index = i))