rm( list=ls(all=TRUE ) )
# Load required packages
required.packages <- c('reticulate', 'stringr', 'svglite', 'quantmod', 'xts', 'zoo')
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(required.packages, library, character.only=TRUE)
# install miniconda if necessary
if (!py_available(initialize = TRUE)) {
	reticulate::install_miniconda()
}
pd <- import("pandas")

# metadata
spn <- 'dspn'
plotsave <- TRUE
sim <- 13
GABA <- 0 
num_gluts <- 18
XX = paste0('upstate vs distance with GABA_', GABA)
path <- paste0('/Documents/GitHub/SPNfinal/dspn/model1/physiological/simulations/sim13/', XX)

# # Set working directory
cd <- sub("/Documents.*", "", getwd())
wd <- paste0(cd, path)
setwd(wd)

# Load Data function
load_data <- function(file_suffix) {
    data <- pd$read_pickle(paste0(XX, file_suffix))
    if(class(data) == "list") sapply(1:length(data), function(ii) data[[ii]]) else unlist(data)
}

# Load data
t <- load_data('_t_data.pickle')
Vsoma <- load_data('_v_data.pickle')
Vdend <- load_data('_vdend_data.pickle')
up <- load_data('_upstates.pickle')
dist <- load_data('_distances.pickle')
mids <- load_data('_midpoints.pickle')
sp <- load_data('_spines.pickle')
den <- load_data('_name_data.pickle')


# Data restructuring
m <- length(t)
upstates <- matrix(up, ncol=1)
distances <- matrix(dist, ncol=1)
midpoints <- matrix(mids, ncol=1)
spines <- matrix(sp, ncol=1)
dend <- matrix(den, ncol=1)
t <- t[,1]

# Clean environment
to_keep <- c('n.folders', 'm', 't', 'Vdend','Vsoma', 'upstates', 'distances', 'midpoints', 'spines', 'dend', 'num_gluts', 'wd', 'plotsave', 'sim')
rm(list=ls()[!(ls() %in% to_keep)])

###### analyse and make graphs ######

# sort by increasing distance from the soma
ind <- sort(distances, index.return=TRUE)$ix
distances <- sapply(1:dim(distances)[2], function(ii) distances[ind,1])
midpoints <- sapply(1:dim(midpoints)[2], function(ii) midpoints[ind,1])
upstates <- sapply(1:dim(distances)[2], function(ii) upstates[ind,1])
spines <- sapply(1:dim(distances)[2], function(ii) spines[ind,1])
dend <- sapply(1:dim(distances)[2], function(ii) dend[ind,1])
Vsoma <- Vsoma[,ind]
Vdend <- Vdend[,ind]

# routine to obtain quarterdrop time
stim.time <- 150
baseline.time = 50
dt <- t[2]-t[1]

idx <-  c((stim.time - baseline.time)/dt, stim.time/dt) # indices for baseline
mu.baseline <- apply(Vsoma[idx[1]:idx[2],], 2, mean)
ind.max <- apply(Vsoma, 2, which.max)
peak.v <- sapply(1:dim(Vsoma)[2], function(ii) Vsoma[ind.max[ii],ii])- mu.baseline
target.v = peak.v - peak.v/4 + mu.baseline
# find index for quarter drop
ind.drop <- sapply(1:dim(Vsoma)[2], function(iii){ 
	which.min( abs( Vsoma[ind.max[iii]:length(t),iii] - target.v[iii] ) ) + ind.max[iii]
})

upstates.check <- sapply(1:dim(Vsoma)[2], function(iii) t[ind.drop[iii]]) - stim.time - num_gluts
cbind(distances, upstates, upstates.check)

ind.max <- apply(Vdend, 2, which.max)			
ind.drop <- sapply(1:dim(Vdend)[2], function(iii){ 
	which.min( abs(Vdend[ind.max[iii]:length(t),iii] - target.v[iii] ) ) + ind.max[iii]
})
upstates.dend <- sapply(1:dim(Vdend)[2], function(iii) t[ind.drop[iii]]) - stim.time - num_gluts


# plot all in Vsoma[[ii]]
fun_color_range <- colorRampPalette(c('slateblue', 'indianred', 'aquamarine4', 'darkorchid', 'sienna1'))   # Apply colorRampPalette
grph1 <- function(mat, t, x.lim, y.lim, iii=1, main='', tcl=-0.3, lwd=1, xaxis=TRUE, y.int=10){
	n <- dim(mat)[2]
	color_range <- fun_color_range(n)  
	ind <- sapply(1:length(x.lim), function(ii) which.min(abs(t - x.lim[ii])) )
	if (xaxis) x.lab='time (ms)' else x.lab=''
	for (ii in 1:n){
		if (ii==1){
			plot(t[ind[1]:ind[2]], mat[ind[1]:ind[2],ii], xlab=x.lab, ylim=y.lim,  xlim=x.lim, ylab='mV', main=main, type="l", bty="n", col=color_range[ii], axes=FALSE, lwd=lwd, lty=1, frame=FALSE) 
			axis(2, at=seq(y.lim[1], y.lim[2], y.int), labels=seq(y.lim[1], y.lim[2], y.int), las=1, tcl=tcl, lwd=lwd)
		}else{
			lines(t[ind[1]:ind[2]], mat[ind[1]:ind[2],ii], col=color_range[ii], lwd=lwd, lty=1)
		}
	}
	if (xaxis) axis(1, at=seq(x.lim[1], x.lim[2], 50), labels=seq(x.lim[1], x.lim[2], 50) - x.lim[1], las=1, tcl=tcl, lwd=lwd)
}

############################################## basic graph ##############################################
lwd=1.5
plt1 <- function(){
	main=''
	tcl=-0.3
	lwd=1.5
	xaxis=TRUE
	x.lim <- c(25, 275)
	x.lab = 'distance from soma (um)'

	y.lim <- c(25*floor(min(upstates)/25),  25*ceiling(max(upstates)/25))
	y.lab = 'time (ms)'
	plot(midpoints[,1], upstates[,1], pch=16, xlim=x.lim, ylim=y.lim, ylab='', xlab='', bty="n", col='slateblue', axes=FALSE, lwd=lwd, frame=FALSE)
	axis(1, at=seq(x.lim[1], x.lim[2], 50), labels=seq(x.lim[1], x.lim[2], 50), las=1, tcl=tcl, lwd=lwd)
	axis(2, at=seq(y.lim[1], y.lim[2], 25), labels=seq(y.lim[1], y.lim[2], 25), las=1, tcl=tcl, lwd=lwd)
	mtext(text = y.lab, side = 2, line = 2.5)
	mtext(text = x.lab, side = 1, line = 2.5)
}

plt2 <- function(){
	x.lab = 'time (ms)'
	y.lab = 'mV'
	x.lim = c(125,375)
	y.lim1 = c(-90,-60)
	y.lim2 = c(-100,-20)
	# plot":
	grph1(Vsoma[,c(11,33)], t, x.lim, y.lim1, iii=iii, main='somatic voltage', lwd=lwd, xaxis=TRUE, y.int=10)
}

plt2_dend <- function(){
	x.lab = 'time (ms)'
	y.lab = 'mV'
	x.lim = c(125,375)
	y.lim1 = c(-90,-60)
	y.lim2 = c(-100,-20)
	# plot
	grph1(Vdend[,c(11,33)], t, x.lim, y.lim2, iii=iii, main='dendritic voltage', lwd=lwd, xaxis=TRUE, y.int=10)
}

if (plotsave) {	
    # Determine the directory to save plots
    save_dir <- str_replace(wd, 'simulations', 'images')
    setwd(save_dir)  # Set working directory to save_dir
    
    # Create a new graphics device for the combined plot
    dev.new(width=6, height=10/3, noRStudioGD=TRUE)
    
    # Set graphical parameters for the plots
    par(mar=c(1, 1, 1, 1), mfrow=c(1,2), oma = c(2, 2, 2, 0), ps=10, cex = 0.9, cex.main = 0.9)
    
    # Generate the plots using plt1 and plt2 functions
    plt1()
    plt2()
    
    # Save the combined plot as a PDF using the current timestamp
    quartz.save(paste0('sim13a ', gsub(':', '-', Sys.time()), '.pdf'), type='pdf')
    
    # Save individual plots as SVG files
    # Create an SVG file for plt2 and add reference lines
    svglite(paste0('sim', sim, ' ', gsub(':', '-', Sys.time()), '.svg'), width=3, height=10/3, pointsize=10)
    plt2()
    abline(h=-60, lty=3, lwd=lwd)  # Reference line at -60
    abline(h=-85, lty=3, lwd=lwd)  # Reference line at -85
    dev.off()  # Close the current graphics device
    
    # Pause execution for a second to ensure file operations don't overlap
    Sys.sleep(1)
    
    # Create another SVG file for plt1
    svglite(paste0('sim', sim, ' ', gsub(':', '-', Sys.time()), '.svg'), width=3, height=10/3, pointsize=10)
    plt1()
    dev.off()  # Close the current graphics device

    # Return to the original working directory
    setwd(wd)
}

if (plotsave) {	
    # Determine the directory to processed data
    save_dir <- str_replace(wd, 'simulations', 'images')
    setwd(save_dir)  # Set working directory to save_dir
	data <- cbind(distance=midpoints[,1], qdrop=upstates[,1])
	# Write to a CSV file
	write.csv(data, file = paste0('sim', sim, ' ', gsub(':', '-', Sys.time()), '.csv'), row.names = FALSE)
	# Return to the original working directory
    setwd(wd)
}

