rm( list=ls(all=TRUE ) )
# Load and install necessary packages
load_required_packages <- function(packages){
    new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
    if(length(new.packages)) install.packages(new.packages)
    invisible(lapply(packages, library, character.only = TRUE))
}

required.packages <- c('reticulate', 'stringr', 'svglite', 'quantmod', 'xts', 'zoo')
load_required_packages(required.packages)

use_condaenv("myenv", required = TRUE)
# Import pandas
pd <- reticulate::import("pandas")

# Metadata settings
plotsave <- TRUE
spn <- 'dspn'
sim <- '6A'
model <- 1
path <- paste0('/Documents/Repositories/SPNfinal/', spn, '/model', model, '/physiological/simulations/sim', sim,'/')

# Set working directory
cd <- sub("/Documents.*", "", getwd())
wd <- paste0(cd, path)
setwd(wd)

# Extract folder data
folders <- dir()

# function to read data
read_data_from_folders <- function(folders, wd) {
	Vtree.path <- list()
  	Vsoma <- list()
  	n.folders <- length(folders)
	for (iii in 1:n.folders){
		XX = folders[iii]
		setwd(XX)
		V_soma <- pd$read_pickle(paste0(XX, '_v_data.pickle'))
		V_dend <- pd$read_pickle(paste0(XX, '_vdend_data.pickle'))
		time <- pd$read_pickle(paste0(XX, '_t_data.pickle'))
		m <- length(time)
		if (iii==1){
			t <- time[[1]]
			gt <- rep(0, n.folders)
			n <- length(time[[1]])
			mat <- matrix(rep(0, n.folders*n), ncol=n.folders)
			Vsoma <- Vdend <- rep(list(mat), m)
		}
		gt[iii] <- pd$read_pickle(paste0(XX, '_gt.pickle'))[[1]][[1]]
		for (ii in 1:m){
			Vsoma[[ii]][,iii] = V_soma[[ii]]
			Vdend[[ii]][,iii] = V_dend[[ii]]
		}
		setwd(wd)
	}
  	# Return data
  	return(list(Vtree.path = Vtree.path, Vsoma = Vsoma, Vdend = Vdend, gt = gt, t = t, m = m, n.folders = n.folders))
}

# read data
data <- read_data_from_folders(folders, wd)
list2env(data, envir = .GlobalEnv)

# Clean environment
to_keep <- c('sim', 'n.folders', 't', 'm', 'gt','Vdend','Vsoma','Vtree.path', 'plotsave','wd')
to_remove <- setdiff(ls(), to_keep)
rm(list = to_remove)

# Sort and reorder based on 'gt' values
ind <- sort(gt, index.return=TRUE)$ix
gt <- gt[ind]
Vsoma <- lapply(1:m, function(ii) Vsoma[[ii]][,ind])
Vdend <- lapply(1:m, function(ii) Vdend[[ii]][,ind])

# Color palette function
fun_color_range <- colorRampPalette(c('slateblue', 'indianred', 'aquamarine4', 'darkorchid', 'sienna1'))

# Plotting function for graphical representation
grph1 <- function(mat, t, x.lim, y.lim, main='', tcl=-0.3, lwd=1, xaxis=TRUE, y.int=10){
	n.folders <- dim(mat)[2]
	color_range <- c('slateblue', 'indianred', 'black')
	ind <- sapply(1:length(x.lim), function(ii) which.min(abs(t - x.lim[ii])) )
	if (xaxis) x.lab='time (ms)' else x.lab=''
	for (ii in 1:n.folders){
		if (ii==1){
			plot(t[ind[1]:ind[2]], mat[ind[1]:ind[2],ii], xlab=x.lab, ylim=y.lim,  xlim=x.lim, ylab='mV', main=main, type="l", bty="n", col=color_range[ii], 					axes=FALSE, lwd=lwd, lty=1, frame=FALSE) 
			axis(2, at=seq(y.lim[1], y.lim[2], y.int), labels=seq(y.lim[1], y.lim[2], y.int), las=1, tcl=tcl, lwd=lwd)
		}else{
			lines(t[ind[1]:ind[2]], mat[ind[1]:ind[2],ii], col=color_range[ii], lwd=lwd, lty=1)
		}
	}
	if (xaxis) axis(1, at=seq(x.lim[1], x.lim[2], 50), labels=seq(x.lim[1], x.lim[2], 50) - x.lim[1], las=1, tcl=tcl, lwd=lwd)
}

# Main plotting script
x.lim = c(135,385)
y.lim1 = c(-90,-50)
y.lim2 = c(-100,-20)
lwd=1
##
dev.new(width=6 ,height=10/3,noRStudioGD=TRUE)
par(mar=c(1, 1, 1, 1), mfrow=c(1,2), oma = c(2, 2, 2, 0), ps=10, cex = 0.9, cex.main = 0.9)
grph1(Vdend[[1]], t, x.lim, y.lim2, main='dendritic voltage', lwd=lwd, xaxis=TRUE, y.int=20)
grph1(Vsoma[[1]], t, x.lim, y.lim1, main='somatic voltage', lwd=lwd, xaxis=TRUE, y.int=10)

if (plotsave) {	
    # Determine the save directory
    save_dir <- str_replace(wd, 'simulations', 'images')
    setwd(save_dir)
    
    # Create a timestamp string for file names
    timestamp <- gsub(':', '-', Sys.time())
    
    # Save the combined plots as a PDF for records
    quartz.save(paste0('sim', sim, ' ', timestamp, '.pdf'), type='pdf')
    
    # Save dendritic voltage plot as an SVG
    svglite(paste0('sim', sim, 'a', ' ', timestamp, '.svg'), width=3.5, height=3.75, pointsize=10)
    grph1(Vdend[[1]], t, x.lim, y.lim2, main='dendritic voltage', lwd=lwd, xaxis=TRUE, y.int=20)
	# Add reference lines to the plot
    abline(h=-20, lty=3, lwd=lwd)
	abline(h=-85, lty=3, lwd=lwd)
	dev.off()
    
    # Save somatic voltage plot as an SVG
    svglite(paste0('sim', sim, 'b', ' ',  timestamp, '.svg'), width=3.5, height=3.75, pointsize=10)
    grph1(Vsoma[[1]], t, x.lim, y.lim1, main='somatic voltage', lwd=lwd, xaxis=TRUE, y.int=10)
	# Add reference lines to the plot
    abline(h=-50, lty=3, lwd=lwd)
    abline(h=-85, lty=3, lwd=lwd)
    dev.off()
    
    # Return to original working directory
    setwd(wd)
}

