library(dplyr)
library("cdcfluview", lib.loc="~/Library/R/3.4/library")
library(mapdata)
library(maps)
library(ggmap)
library(ggplot2)
usa <- map_data("usa")
county_map_data <- map_data("county")
state_map_data <- map_data("state")
county_names <- as.list(unique(county_map_data$subregion))
county_data <- read.csv("~/Box Sync/Heterogeneity_overview_paper/sources/map_sources/2016_US_County_Level_Presidential_Results.csv")
stroke_data <- read.csv("~/Box Sync/Heterogeneity_overview_paper/sources/map_sources/stroke_all.csv")
drink_data <- read.csv("~/Box Sync/Heterogeneity_overview_paper/sources/map_sources/drinks_processed.csv")
county_data$county_name <- as.character(county_data$county_name)
stroke_data$county.Name <- as.character(stroke_data$County.Name)
for (i in 1:length(county_map_data$subregion)) {
  if (length(as.numeric(unique(county_data$per_dem[grep(county_map_data$subregion[i],county_data$county_name,ignore.case = TRUE)]))) == 0) {
    county_map_data$per_dem[i] = NA
    county_map_data$per_gop[i] = NA
    county_map_data$diff[i] = NA
  } else {
    county_map_data$per_dem[i] <- as.numeric(unique(county_data$per_dem[grep(county_map_data$subregion[i],county_data$county_name,ignore.case = TRUE)]))
    county_map_data$per_gop[i] <- as.numeric(unique(county_data$per_gop[grep(county_map_data$subregion[i],county_data$county_name,ignore.case = TRUE)]))
    county_map_data$diff[i] <- as.numeric(county_map_data$per_dem[i] - county_map_data$per_gop[i])
    county_map_data$fips[i] <- as.integer(unique(county_data$combined_fips[grep(county_map_data$subregion[i],county_data$county_name,ignore.case = TRUE)]))
  }
  if (length(as.numeric(unique(stroke_data$Adults.35.[grep(county_map_data$subregion[i],stroke_data$County.Name,ignore.case = TRUE)]))) == 0) {
    county_map_data$stroke_mortality[i] = NA
  } else {
    county_map_data$stroke_mortality[i] <- as.numeric(unique(stroke_data$Adults.35.[grep(county_map_data$subregion[i],stroke_data$County.Name,ignore.case = TRUE)]))
  }
  if (length(as.numeric(unique(drink_data$pop.coke.or.soda[grep(county_map_data$subregion[i],drink_data$County_Name,ignore.case = TRUE)]))) == 0) {
    county_map_data$drink_name[i] = NA
  } else {
    county_map_data$drink_name[i] <- as.numeric(unique(drink_data$pop.coke.or.soda[grep(county_map_data$subregion[i],drink_data$County_Name,ignore.case = TRUE)]))
  }   
}
no_axes <- theme(
  axis.text = element_blank(),
  axis.line = element_blank(),
  axis.ticks = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.title = element_blank()
)
county_map_object <- ggplot(data = usa,mapping = aes(x = long, y = lat, group = group)) + coord_fixed(1.3) + geom_polygon(color="black",fill="gray")
county_map_object <- county_map_object + theme_nothing() + geom_polygon(data=state_map_data,fill=NA,color="black") + geom_polygon(color="black",fill=NA)
county_map_object
per_dem_map <- county_map_object + geom_polygon(data = county_map_data,aes(fill = per_dem)) + geom_polygon(color= "black",fill=NA) + theme_bw() + no_axes
per_dem_map + scale_fill_gradient(trans = "log10") + geom_polygon(data=state_map_data,fill=NA,color="black") + geom_polygon(color="black",fill=NA)
per_gop_map <- county_map_object + geom_polygon(data = county_map_data,aes(fill = per_gop)) + geom_polygon(color= "black",fill=NA) + theme_bw() + no_axes
per_gop_map + scale_fill_gradient(trans = "log10") + geom_polygon(data=state_map_data,fill=NA,color="black") + geom_polygon(color="black",fill=NA)
per_diff_map <- county_map_object + geom_polygon(data = county_map_data,aes(fill = diff)) + geom_polygon(color= "black",fill=NA) + theme_bw() + no_axes
per_diff_map + scale_fill_gradient2(low = "#ff0000",mid="#D3D3D3",high="#0000ff") + geom_polygon(data=state_map_data,fill=NA,color="black") + geom_polygon(color="black",fill=NA)
stroke_map <- county_map_object + geom_polygon(data = county_map_data,aes(fill = stroke_mortality)) + geom_polygon(color= "black",fill=NA) + theme_bw() + no_axes
stroke_map + scale_fill_gradientn(colours = c("#F0DFEF","#D2A9D5","#9E62A3","#6A2770","#630D66"),values = c(0,0.32,0.36,0.4,1)) + geom_polygon(data=state_map_data,fill=NA,color="black") + geom_polygon(color="black",fill=NA)
dialect_map <- county_map_object + geom_polygon(data = county_map_data,aes(fill = drink_name)) + geom_polygon(color = "black",fill=NA) + theme_bw() + no_axes
dialect_map + scale_fill_gradientn(colours = c("grey50","blue","blue","green","green","red","red"),values=c(0,2/7,3/7,4/7,5/7,6/7,1)) + geom_polygon(data=state_map_data,fill=NA,color="black") + geom_polygon(color="black",fill=NA)
