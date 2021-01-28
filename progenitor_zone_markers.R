rm(list=ls())
library(readr)
library(stringr)
library(ggsignif)
library(ggplot2)
library(EnvStats)

plot_PZ = T # will plot PZ, if set to TRUE.
y_axis_label = "Distance from DTC (cell diameter"

dir_path1 <- "/run/media/areeba/New Volume/Acad/Image_analysis/2021/01.12.2021_PX3525_EdU_pSUN-1_CYE-1/"
#dir_path2 <- "/run/media/areeba/New Volume/Acad/Image_analysis/2021/01.12.2021_PX3525_EdU_pSUN-1_CYE-1/PHX3525/"

color_list <- c("deepskyblue", "black", "lawngreen", "blue", "indianred1", "purple", "royalblue", "olive", "cyan", "darkturquoise", "aquamarine4", "darkgreen", "navajowhite4", "cornflowerblue", "magenta", "violet", "sienna", "red", "darkorange", "slategray", "mediumslateblue", "lightseagreen", "deeppink", "rosybrown", "darkgoldenrod", "dodgerblue", "olivedrab", "darkviolet", "forestgreen", "orange")

# colors will be used for plotting, make changes or create new list
colors <- color_list

colors <- c("dodgerblue",  "black", "orange", "indianred1", "sienna")


##########################################################
## Function  getMarkersData
##########################################################

getMarkersData <- function(markers_dirs){
  raw_data <- data.frame()
  for (markers_dir in markers_dirs) {
    genotype <- getGenotype(markers_dir)
    genotype_markers <- getMarkers(markers_dir, genotype, file_type = "csv")
    raw_data <- rbind(raw_data, genotype_markers)
  }
  colnames(raw_data) <- c("genotype", "gonad", "marker1", "marker2")
  return(raw_data)
}  

##########################################################
getMarkers <- function(markers_dir, genotype, file_type = "csv") {
  csv_files <- list.files(path = markers_dir, full.names = TRUE, recursive = FALSE, pattern = file_type, include.dirs = FALSE)
  genotype_markers <- data.frame()
  for (csv_file in csv_files){
    gonad <- stringr::str_replace(basename(csv_file), "csv", "tiff")
    temp_data <- as.data.frame(suppressMessages(read_csv(file = csv_file, col_names = TRUE, progress = F)))
    
    temp_data <- as.data.frame(table(temp_data$Counter), stringsAsFactors =FALSE)
    temp_data$Var1 <- as.numeric(temp_data$Var1)
    temp_data <- temp_data[order(temp_data$Var1),]
    gonad_markers <- c(genotype, gonad, temp_data$Freq[temp_data$Var1 == 1], ifelse(max(temp_data$Var1) > 2, temp_data$Freq[temp_data$Var1 == 2], NA))
    genotype_markers <- rbind(genotype_markers, gonad_markers)
  }
  colnames(genotype_markers) <- c("genotype", "gonad", "marker1", "marker2")
  return(genotype_markers)
}

########################################################

# function getMarkerDirs

getMarkerDirs <- function(dirs){
  markers_dirs <- NULL
  for (dir_path in dirs){
    markers_dir <- list.files(path = dir_path, full.names = TRUE, recursive = TRUE, pattern = "markers", include.dirs = TRUE)
    markers_dir <- markers_dir[file.info(markers_dir)$isdir]
    if (length(markers_dir)==0){
      if (stringr::str_detect(string = dir_path, "markers")){
        markers_dirs <- c(markers_dirs, stringr::str_extract(dir_path, paste0("^.*markers.+?", .Platform$file.sep))) 
      }else{
        markers_dirs <- c(markers_dirs, dir_path)
        print("No markerss were found")
      }
    }else{
      markers_dirs <- c(markers_dirs, markers_dir)
    }
  }
  print(markers_dirs)
  return(markers_dirs)
}



########################################################################
# for progenitor zone length plotting
###########################################################

geom_mean <- function(size = 8, color = "black", width = 0.2) {
  list(
    stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = 0.3, linetype = "solid", color = color, size = 2),
    stat_summary(fun.data = "data_summary", geom = "errorbar", width = 0.2, color = color)
  )
}

##############################################################################################

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}
#######################################################################


getGenotype <- function(dir_name){
  while (TRUE) {
    genotype <- basename(dir_name)
    if (genotype %in% c("stitched_images", "tiff", "markers", "tif")){
      dir_name <- dirname(dir_name)
    }else{
      return(genotype)
    }
  }
}

##########################################################

prog_zone_plot <- function(genotypeToPlot, plot_data) {
  row_plot_data <- plot_data[plot_data$genotype %in% genotypeToPlot,]
  row_plot_data$genotype[row_plot_data$genotype=="N2"] <- "WT"
  genotypeToPlot[genotypeToPlot=="N2"] <- "WT"
  
  toCompare <- list(genotypeToPlot[c(1,2)]) #, genotypeToPlot[c(1,3)], genotypeToPlot[c(2,3)])
  signifMapping = c("***"=0.0001, "**"=0.001, "*"=0.01)
  #colorList <- darkColorPicker(length(genotypeToPlot))
  colorList <- c("orange", "lawngreen", "royalblue")
  
  #row_plot_data$prog_length <- ifelse(is.na(row_plot_data$marker2), row_plot_data$marker1, row_plot_data$marker1 + row_plot_data$marker2)
  
  maxY <- max(row_plot_data$marker1 + row_plot_data$marker2)
  print(maxY)
  maxY <- 10 * (round(maxY/10,0) +1)
  if (plot_PZ) {
    y = "marker1 + marker2"
  }else{
    y = "marker2"
  }
  row_plot <- ggplot(row_plot_data, aes(x = factor(genotype, levels = genotypeToPlot[length(genotypeToPlot):1])  , y = (eval(parse(text=y))), color = genotype, fill = genotype)) + 
    geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.5, dotsize = 1) +
    geom_mean(color = "black") +
    theme_bw(base_size = 9, base_family = "arial") +
    theme(axis.text.x = element_text(size=9, hjust=.5,vjust=.5,face="plain"), legend.position = "none") + theme(axis.text = element_text(face = "italic")) +
    expand_limits(y = c(0, maxY)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(labels = function(x) str_wrap(str_replace_all(x, "\\_", " "), width = 10),
                     expand = c(0, 0)) +
    coord_flip() + theme(panel.grid = element_blank(), axis.ticks = element_line(size = 1), 
                         axis.title = element_text(size = 9), 
                         axis.text = element_text(size = 9)) + labs(x = NULL, y = "Distance in cell diameter") +
    theme(axis.text = element_text(colour = "black")) +
    scale_fill_manual(values = colorList) +
    scale_color_manual(values = colorList) +
    geom_signif(aes(textsize = 15), comparisons = toCompare, test = "t.test", step_increase = 0.05) +
    stat_n_text(y.pos = 2)
  
  print(row_plot)
  Sys.sleep(1)
  
  #  savePlots(fileName = "dotplot.svg", width = 9, autoSave = F)
  #  savePlots(fileName = "dotplot.pdf", width = 9, autoSave = F)
}




all_dir <- ls()[grep("dir_path[0-9]?", ls())]
dirs <- NULL
for (dir in all_dir) dirs <- c(dirs, get(dir))
print(dirs)

#script_dir <- getSrcDirectory(function(x) {x})  # to get path of the current script
markers_dirs <- getMarkerDirs(dirs)
rm(markers_data)
markers_data <- getMarkersData(markers_dirs)
markers_data[is.na(markers_data)] <- 0
markers_data$marker1 <- as.numeric(markers_data$marker1)
markers_data$marker2 <- as.numeric(markers_data$marker2)

frequency <- as.data.frame(table(markers_data$genotype), stringsAsFactors =FALSE)
frequency$Freq <- as.character(frequency$Freq)

genotypeToPlot <- c("N2",  "PROM-1_HA", "PROM-1_HA_puf-8") #, "gld-1_LAG-1_HA", "gld-2_LAG-1_HA")
genotypeToPlot <- unique(markers_data$genotype)
prog_zone_plot(genotypeToPlot, markers_data)

# prog_zone_plot(genotypeToPlot, plot_data)
# 
# row_plot_data <- plot_data[plot_data$genotype %in% genotypeToPlot,]
# row_plot_data$genotype[row_plot_data$genotype=="N2"] <- "WT"
# 
# genotypeToPlot[genotypeToPlot=="N2"] <- "WT"


#ggsave(filename = "~/Dropbox/Acad/projects/LAG-1/Manuscript/figure_plots/progenitor_zone_4h.svg", width = 3, height = 2)


#write.csv(x = markers_data, file="/run/media/areeba/New Volume/Acad/Image_analysis/2021/01.12.2021_PX3525_EdU_pSUN-1_CYE-1/row_data.csv", row.names = F)
