# Functions ot read in structure_f files
# TO DO:
# include into plot_struc:
# 1. sure palette colors =< k ; and give guidance
# 2. Reverse plot order
# 3. Plot by Pop | clustering
# Evaluate & include other clustering methods

install.packages("reshape")
install.packages("ggplot2")
install.packages("seriation")
install.packages("sata.table")
library(reshape)
library(ggplot2)
library(seriation)  # Different clsutering methods for ordering series
library(data.table)

read_structure_f <- function(structure_file){
  i <- grep("Inferred ancestry",readLines(structure_file, warn=FALSE))
  j <- grep("Estimated Allele",readLines(structure_file, warn=FALSE))
  struc_q <- read.table(structure_file, skip = i+1, nrows = j-i-4, row.names=1) # data frame containg sturcture q data per K
  # struc_q <- struc_q[1:10,] # temprorary subset for dev, remove for final
  struc_q[4] <- NULL
  return(struc_q)
}

#' @importFrom seriation seriate get_order
cluster_structure <- function(q_matrix, method = 1){
  METHODS = c("HC", "MDS")
  if (all((1:2) != method)) {
    stop("The method selection is not valid")
  }
  k <- ncol(q_matrix) # Column of last K in structure file
  d <- dist(as.matrix(q_matrix[,4:k])) # infer sort order by cluster analysis
  if (method == 1) {
    sort_order <- get_order(seriate(d, method="HC"))
  }
  else if (method == 2) {
    sort_order <- get_order(seriate(d, method="MDS"))
  }
  return(sort_order)
}

#' @importFrom reshape melt
format_for_plot <- function(q_matrix){
  q_matrix <- q_matrix[sort_order,]
  a <- melt(q_matrix[1:nrow(q_matrix),], id=c("V2"), measure.vars = c(4:ncol(q_matrix)))
  a$V2 <- factor(a$V2, levels = rev(unique(a$V2)))
  a <- data.table(a)
  #recaculate all proportions to sum to exactly one:
  a[, sum := sum(value), by=list(V2)] 
  a[, proportion := value/sum]
  #delete unneeded rows
  a <- a[, c("value","sum") := NULL]
  plot_data <- a
  return(plot_data)
}

#' @importFrom ggplot2 ggplot
plot_struc <-function(plot_data, method = 1, palette = 1){
  METHODS = c("full", "reduced", "stripped")
  if (all((1:3) != palette)) {
    stop("The palette selection is not valid")
  }
  if (palette == 1) {
    pal = "Set1"
  }
  else if (palette == 2) {
    pal = "Set3"
  }
  else if (palette == 3) {
    pal = "Paired"
  }
  if (all((1:3) != method)) {
    stop("The method selection is not valid")
  }
  if (method == 1) {
    q_plot <- ggplot(data = plot_data, aes(plot_data$V2, proportion, fill = variable, width=1)) + 
      geom_bar(stat = "identity") + 
      ylab("Proportion") +
      xlab( "Individual") +
      scale_fill_brewer(palette=pal, name="K", labels=c(1:length(unique(plot_data$variable)))) +
      theme(
        # axis.text = element_blank(), 
        # axis.ticks = element_blank(), 
        axis.text.x = element_text(angle =90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
  }
  else if (method == 2) {
    q_plot <- ggplot(data = plot_data, aes(plot_data$V2, proportion, fill = variable, width=1)) + 
      geom_bar(stat = "identity") + 
      ylab("Proportion") +
      xlab( "Individual") +
      scale_fill_brewer(palette=pal, name="K", labels=c(1:length(unique(plot_data$variable)))) +
      theme(
        # axis.text = element_blank(), 
        # axis.ticks = element_blank(), 
        axis.text.x = element_text(angle =90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
  }
  else if (method == 3) {
    q_plot <- ggplot(data = plot_data, aes(plot_data$V2, proportion, fill = variable, width=1)) + 
      geom_bar(stat = "identity") + 
      scale_fill_brewer(palette=pal, name="K", labels=c(1:length(unique(plot_data$variable)))) +
      theme(
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
  }
  return(q_plot)
}

structure_file <- "casia_f" # Read STRUCTURE file in xxx_f format and extract data for plotting
structure_file <- "MX_two_run_40_f" # Read STRUCTURE file in xxx_f format and extract data for plotting
structure_file <- "SA_New_run_19_f" # Read STRUCTURE file in xxx_f format and extract data for plotting

q_matrix <- read_structure_f(structure_file)
sort_order <- cluster_structure(q_matrix, method = 1)
plot_data <- format_for_plot(q_matrix)
plot_struc(plot_data, method = 3, palette = 1)

###############################################################
# Nik's sandbox

q_plot <- ggplot(data = plot_data, aes(plot_data$V2, proportion, fill = variable, width=1)) + 
  geom_bar(stat = "identity") + 
  scale_fill_brewer(palette="Set1", name="K", labels=c(1:length(unique(plot_data$variable)))) +
  theme(
    axis.text = element_blank(), 
    axis.ticks = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank())
