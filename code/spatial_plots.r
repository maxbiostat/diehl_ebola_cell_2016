library(maptools)
library(ggplot2)
#######
All <- TRUE ## should use all data or only complete data?
if(All){
  MutDat <- read.csv("../data/mutations_per_location_GP82.csv", header = TRUE)
  freqV <- MutDat$V/(MutDat$A + MutDat$V)
  names(freqV) <- MutDat$location
  freqV
}else{
  EVD_GP82 <- read.csv("../data/EVD_case_fatality_GP82.csv", header = TRUE)
  EVD_GP82 <- subset(EVD_GP82, location != "?")
  EVD_GP82$location <- paste(EVD_GP82$location)
  EVD_GP82 <- subset(EVD_GP82, country == "GIN")
  (MutDat <- xtabs(~GP82 + location, EVD_GP82) )
  ( freqV <- MutDat[2, ]/colSums(MutDat) )
}
prop_table <- data.frame(location = names(freqV), prop_V = as.numeric(freqV))
Map <- readShapePoly("../data/gis/locations.shp")
Map@data$location <- gsub("GrandCapeMo","GrandCapeMount", gsub("WesternRura", "WesternRural", gsub("WesternUrba", "WesternUrban", Map@data$location)))
Map.f <- fortify(Map, region = "location")
Map.f$location <- Map.f$id
merge.shp <- merge(x = Map.f, y = prop_table, by = "location", all = TRUE)
merge.shp <- merge(merge.shp, Map@data)
final.plot <- merge.shp[order(merge.shp$order), ]
Country_info <- read.csv("../data/location_country.csv", header = TRUE)

getCountry <- function(loc) as.character(Country_info$country[match(loc, Country_info$location)])
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

library(RColorBrewer)
Cols <- brewer.pal(9, "Blues")

p0 <- ggplot() + 
  ggtitle("") +
  geom_polygon(data = final.plot, 
               aes(x = long, y = lat, group = group, colour = "black", fill = prop_V),
               size = .2) + 
  geom_polygon(data = Map.U, 
               aes(x = long, y = lat, group = group, colour = "red"),
               size = .7, alpha = 0.00001) + 
  scale_color_manual(values = c("black", "red")) +
  labs(fill = "GP82V proportion") +
  scale_fill_continuous(low = Cols[5], high = Cols[9], 
                         na.value = "white") + 
  guides(colour = FALSE) + 
  theme_bw() +
  coord_map()

pdf("../plots/GP82V_proportion_map.pdf")
p0
dev.off()