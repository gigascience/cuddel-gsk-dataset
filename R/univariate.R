# Title     : univariate.R
# Objective : Analyse GSK dataset with univariate statistics
# Created by: peterli
# Created on: 9/10/2018

library(ggfortify)
library(ggplot2)
library(cowplot)

source("functions.R")

# Location of GSK data set
datadir = "/home/peter/"
# File path for positive files
pos_dir = paste0(datadir, "gsk/raw/esi_pos/netcdf")
# Output directory
output_path <- paste0(pos_dir, "/output")

# Read in data
pos_dataMatrix <- read.table(file=paste(output_path, "dataMatrix.tsv", sep="/"), header=TRUE, row.names=1)
rownames(pos_dataMatrix) <- pos_dataMatrix[,1]
pos_dataMatrix <- pos_dataMatrix[, 2:length(pos_dataMatrix)]
pos_variableMetadata <- read.table(file=paste(output_path, "variableMetadata.tsv", sep="/"), header=TRUE, row.names=1)
pos_sampleMetadata <- read.table(file=paste(output_path, "sampleMetadata.tsv", sep="/"), header=TRUE, row.names=1)

#' Calculate means for multiple groups
#'
#' @param x A vector of values for mean calculation
#' @param tp5RegARegB_regimen_fac A vector of values representing factors associated with x vector
#' @return The matrix of meanvalues for each factor represented in tp5RegARegB_regimen_fac
#' @examples
#' calcMean(10, 1)
calcMean <- function(x, fac) {
    tapply(x, fac, mean)
}

#' Calculate means for multiple groups
#'
#' @param x A vector of values for mean calculation
#' @param tp5RegARegB_regimen_fac A vector of values representing factors associated with x vector
#' @return The matrix of mean values calculated for each factor represented in tp5RegARegB_regimen_fac
#' @examples
#' calcMean(10, 1)
calcStdDev <- function(x, fac) {
    tapply(x, fac, sd)
}

#' Calculate t-test for multiple groups
#'
#' @param x A vector of values for mean calculation
#' @param tp5RegARegB_regimen_fac A vector of values representing factors associated with x vector
#' @return The matrix of mean values calculated for each factor represented in tp5RegARegB_regimen_fac
#' @examples
#' calcMean(10, 1)
performTTest <- function(x, fac) {
    # Need to split data based on the 2 factor groups
    grp_data <- split(x, fac)
    ttest <- t.test(grp_data[[1]], grp_data[[2]])
    return(ttest$p.value)
}

# Get timepoints as a factor
tp_fac <- factor(pos_sampleMetadata[, "timepoint"])
# Get the levels of the timepoint factor
tp_levs <- levels(tp_fac)
# Remove -1 timepoint factor level as this corresponds to QCs
tp_levs <- tp_levs[2:length(tp_levs)]

# Transpose
t_pos_dataMatrix <- t(pos_dataMatrix)

# Create 5 matrices to hold results from analysis of regimens A and B - 23 rows * 8 TP levels
matrix_names <- c("tp_RegA_means", "tp_RegA_stdevs", "tp_RegB_means", "tp_RegB_stdevs", "tp_pvalues")
for(i in 1:length(matrix_names)) {
    m <- matrix(0, nrow=23, ncol=8)
    colnames(m) <- c("TP0", "TP1", "TP4", "TP5", "TP8", "TP12", "TP16", "TP24")
    rownames(m) <- colnames(t_pos_dataMatrix)
    assign(matrix_names[i], m)
}

# Loop over timepoint levels to add data into 5 matrices
for(i in 1:length(tp_levs)) {
    # Extract sample IDs data based on timepoint 5
    tp_sample_meta <- pos_sampleMetadata[pos_sampleMetadata[, "timepoint"] == tp_levs[i], ]
    # Get regimen factor
    # Extract sample IDs for regimen A and B data from tp_sample_meta
    tp_reg_sample_meta <- tp_sample_meta[tp_sample_meta[, "regimen"] == "A" | tp_sample_meta[, "regimen"] == "B", ]
    # Get regimen factor
    regimen_fac <- factor(tp_reg_sample_meta[, "regimen"])

    # Extract peak feature values for Reg A and B
    data <- t_pos_dataMatrix[rownames(tp_reg_sample_meta), ]

    # Call function to do t-tests
    if(nrow(data) != length(regimen_fac)) {
        stop("The number of rows in data needs to equal length of fac!!")
    } else {
        means <- apply(data, 2, calcMean, fac=regimen_fac)
        t_means <- t(means)
        tp_RegA_means[, i] <- t_means[, "A"]
        tp_RegB_means[, i] <- t_means[, "B"]

        stdevs <- apply(data, 2, calcStdDev, fac=regimen_fac)
        t_stdevs <- t(stdevs)
        tp_RegA_stdevs[, i] <- t_stdevs[, "A"]
        tp_RegB_stdevs[, i] <- t_stdevs[, "B"]

        pvalues <- apply(data, 2, performTTest, fac=regimen_fac)
        t_pvalues <- t(pvalues)
        tp_pvalues[, i] <- t_pvalues
    }
}

# Transform data into tidyr format for plotting as line graph using ggplot2
tidyr_data <- data.frame(peak=character(),
    regimen=character(),
    tp=character(),
    mean=numeric(),
    sd=numeric(),
    stringsAsFactors=TRUE)

for(i in 1:length(rownames(tp_RegA_means))) {
    peak_name <- rownames(tp_RegA_means)[i]

    for(x in 1:length(colnames(tp_RegA_means))) {
        newRow <- data.frame(peak=peak_name, regimen="A", timepoint=colnames(tp_RegA_means)[x], mean=tp_RegA_means[i, x], sd=tp_RegA_stdevs[i, x])
        tidyr_data <- rbind(tidyr_data, newRow)
    }
}

for(i in 1:length(rownames(tp_RegB_means))) {
    peak_name <- rownames(tp_RegB_means)[i]

    for(x in 1:length(colnames(tp_RegB_means))) {
        newRow <- data.frame(peak=peak_name, regimen="B", timepoint=colnames(tp_RegB_means)[x], mean=tp_RegB_means[i, x], sd=tp_RegA_stdevs[i, x])
        tidyr_data <- rbind(tidyr_data, newRow)
    }
}

# Error bars overlap so use position_dodge to move them horizontally
pd <- position_dodge(0.1)  # Move them .05 to the left and right

# Use subset function to display only P210 data
p1 <- ggplot(subset(tidyr_data, peak %in% c("P210")), aes(x=timepoint, y=mean, color=regimen, group=regimen)) +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), color="black", width=.1, position=pd) +
    geom_line(position=pd) +
    geom_point(position=pd, size=3, shape=21, fill="white") +
    xlab('Time Point') +
    ylab('P210') +
    scale_colour_hue(name="Regimen",        # Legend label, use darker colors
    breaks=c("A", "B"),
    labels=c("A", "B"),
    l=40) +                                 # Use darker colors, lightness=40
    ggtitle("Compare regimens A and B") +
    expand_limits(y=0.055) +                # Expand y range
    scale_y_continuous(breaks=0:20*4) +     # Set tick every 4
    theme_bw() +
    theme(legend.justification=c(1,0),
    legend.position=c(1,0))                 # Position legend in bottom right

p2 <- ggplot(subset(tidyr_data, peak %in% c("P257")), aes(x=timepoint, y=mean, color=regimen, group=regimen)) +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), color="black", width=.1, position=pd) +
    geom_line(position=pd) +
    geom_point(position=pd, size=3, shape=21, fill="white") +
    xlab('Time Point') +
    ylab('P210') +
    scale_colour_hue(name="Regimen",        # Legend label, use darker colors
    breaks=c("A", "B"),
    labels=c("A", "B"),
    l=40) +                                 # Use darker colors, lightness=40
    ggtitle("Compare regimens A and B") +
    expand_limits(y=0.055) +                # Expand y range
    scale_y_continuous(breaks=0:20*4) +     # Set tick every 4
    theme_bw() +
    theme(legend.justification=c(1,0),
    legend.position=c(1,0))                 # Position legend in bottom right

# Save last ggplot object
ggsave("ggplot.pdf")

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    library(grid)

    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)

    numPlots = length(plots)

    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
        ncol = cols, nrow = ceiling(numPlots/cols))
    }

    if (numPlots==1) {
        print(plots[[1]])

    } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
            layout.pos.col = matchidx$col))
        }
    }
}

pdf("plots.pdf")
multiplot(p1, p2, cols=4)
dev.off()

peaks <- unique(tidyr_data$peak)
plot_names <- paste0(peaks, "_plot")

for(i in 1:length(peaks)) {
    thePeak <- peaks[i]
    # Create object name for this plot
    plot_name <- paste0(thePeak, "_plot")

    aPlot <- ggplot(subset(tidyr_data, peak %in% c(as.character(thePeak))),
        aes(x=timepoint, y=mean, color=regimen, group=regimen)) +
        geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), color="black", width=.1, position=pd) +
        geom_line(position=pd) +
        geom_point(position=pd, size=2, shape=21, fill="white") +
        xlab('Time Point') +
        ylab(thePeak) +
        scale_colour_hue(name="Regimen",        # Legend label, use darker colors
            breaks=c("A", "B"),
            labels=c("A", "B"),
            l=40) +                                 # Use darker colors, lightness=40
        # ggtitle("Compare regimens A and B") +
        expand_limits(y=0.055) +                # Expand y range
        scale_y_continuous(breaks=0:20*4) +     # Set tick every 4
        theme_bw() +
        # theme(legend.justification=c(1,0),
        #       legend.title = element_blank(),
        #       legend.text = element_text(size=8),
        #       legend.box = "horizontal",
        #       legend.position=c(1,0))         # Position legend in bottom right
        theme(legend.position="none",           # Do not display legend
              axis.title.x=element_blank(),     # Do not display x axis title
              axis.text.x = element_text(size=6))
    assign(plot_name, aPlot)
}

# Get the objects assocoiated with the object string names
l <- mget(plot_names)

# Plot graphs
pdf("plots3.pdf")
multiplot(plotlist=l, cols=3)
dev.off()

# Using cowplot to display multiple graphs
graphs <- list(P210_plot, P257_plot, P439_plot, P478_plot, P727_plot, P815_plot, P816_plot, P884_plot, P988_plot, P1051_plot, P2554_plot, P2617_plot, P4141_plot, P4390_plot, P4402_plot, P4514_plot, P4518_plot, P4625_plot, P4727_plot, P4732_plot, P4736_plot, P4993_plot, P6073_plot)
# plot2by2 <- plot_grid(P210_plot, P257_plot, P439_plot, P478_plot, ncol=2)
plot2by2 <- plot_grid(plotlist=graphs, ncol=2, scale = c(.5, .5, .5, .5))
save_plot("plot2by2.png", plot2by2,
    ncol = 2, # we're saving a grid plot of 2 columns
    nrow = 2, # and 2 rows
    # each individual subplot should have an aspect ratio of 1.3
    base_aspect_ratio = 1.3
)