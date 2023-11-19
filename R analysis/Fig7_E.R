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
sim <- '320f'
model <- 1

# set appropriate path
path <- paste0('/Documents/Repositories/SPNfinal/', spn, '/model', model, '/physiological/simulations/sim', sim,'/')

# Set working directory
cd <- sub("/Documents.*", "", getwd())
wd <- paste0(cd, path)
setwd(wd)

names <- dir()

Vtree.path <- pd$read_pickle('v_tree.pickle')
V_soma <- pd$read_pickle('v_data.pickle')
V_dend <- pd$read_pickle('vdend_data.pickle')
time <- pd$read_pickle('t_data.pickle')
nsims <- length(time)

t <- time[[1]]
gt <- pd$read_pickle('gt.pickle')
gt <- sapply(1:length(gt), function(ii) gt[[ii]][[1]])
Vsoma <- sapply(1:length(V_soma), function(ii) V_soma[[ii]])
Vdend <- sapply(1:length(V_dend), function(ii) V_dend[[ii]])

# isolate 0 GABA and 0 Glut traces
Vsoma_GABA0 <- Vsoma[,1]
Vdend_GABA0 <- Vdend[,1]

gt <- gt[2:length(gt)]
Vsoma <- Vsoma[,2:dim(Vsoma)[2]]
Vdend <- Vdend[,2:dim(Vdend)[2]]

# sort columns by gt
ind <- sort(gt, index.return=TRUE)$ix
gt <- gt[ind]
Vsoma <- Vsoma[,ind]
Vdend <- Vdend[,ind]

# i_mechs <- pd$read_pickle('i_mechs.pickle')
# lst.names <- names(i_mechs)
# g_gaba <- as.numeric(gsub("[^0-9.e-]", "", lst.names))

g_gaba <- c(0e+00, 1e-06, 3e-06, 1e-05, 3e-05, 1e-04, 3e-04, 1e-03, 3e-03, 1e-02, 3e-02)

# Clean environment
to_keep <- c('g_gaba',  't', 'gt', 'Vdend', 'Vsoma', 'Vdend_GABA0', 'Vsoma_GABA0', 'nGlut', 'wd', 'plotsave', 'names', 'sim')
rm(list = ls()[!ls() %in% to_keep])


nGlut=15
dt <- t[2]-t[1]
# remove first 200 ms (burntime)
burntime <- 200
idx <- 200/dt

Vsoma <- Vsoma[idx:dim(Vsoma)[1],]
Vdend <- Vdend[idx:dim(Vdend)[1],]
Vsoma_GABA0 <- Vsoma_GABA0[idx:length(Vsoma_GABA0)]
Vdend_GABA0 <- Vdend_GABA0[idx:length(Vdend_GABA0)]
t <- t[idx:length(t)] - t[idx]

fun_color_range <- colorRampPalette(c('slateblue', 'indianred', 'aquamarine4', 'darkorchid', 'sienna1'))   # Apply colorRampPalette

timing <- gt
x.lim=c(min(timing)-10,min(timing)+190)

# Final plot:
stim.time <- 150
baseline <- 25
offset <- stim.time - baseline
t.max <- 325 # max x lim value

t1 <- t[t < t.max]

ind1 <- offset/dt
ind2 <- length(t1)
t2 <- t[ind1:ind2] 

ind1 <- offset/dt
ind2 <- length(t1)
t2 <- t[ind1:ind2] 


plot.timing <- gt
ind <- gt %in% plot.timing
n.traces <- length(plot.timing)

Vdend.plot <- Vdend[ind1:ind2,ind]
Vsoma.plot <- Vsoma[ind1:ind2,ind]

x.lim = c(0,t.max)
y.lim1 = c(-90,-60)
y.lim2 = c(-96,-30)
lwd=1.5
tcl=-0.3
lwd=1
xaxis=TRUE
y.int=10
x.int=50

funplot1 <- function(V_GABA0, V, y.lim, main='dendritic voltage', ab1=-30, ab2=-85){
	main = 'dendritic voltage'
	plot(t2-offset, V_GABA0[ind1:ind2], xlab='time (ms)', ylim=y.lim,  xlim=x.lim, ylab='mV', main=main, type="l", bty="n", col='black', axes=FALSE, lwd=lwd, lty=1, frame=FALSE) 
	axis(2, at=seq(y.lim[1], y.lim[2], y.int), labels=seq(y.lim[1], y.lim[2], y.int), las=1, tcl=tcl, lwd=lwd)
	if (xaxis) axis(1, at=seq(x.lim[1], x.lim[2], x.int), labels=seq(x.lim[1], x.lim[2], x.int) - x.lim[1], las=1, tcl=tcl, lwd=lwd)
	color_range <- fun_color_range(n.traces) 
	for (ii in 1:n.traces){
		lines(t2, V[,ii], col=color_range[ii], lwd=lwd, lty=1)
		}
	abline(h=ab1, lty=3, lwd=lwd)
	abline(h=ab2, lty=3, lwd=lwd)
}
	
dev.new(width=6 ,height=10/3,noRStudioGD=TRUE)
par(mar=c(1, 1, 1, 1), mfrow=c(1,2), oma = c(2, 2, 2, 0), ps=10, cex = 0.9, cex.main = 0.9)
funplot1(Vdend_GABA0, Vdend.plot, y.lim2, main='dendritic voltage', ab1=-30, ab2=-85)
funplot1(Vsoma_GABA0, Vsoma.plot, y.lim1, main='somatic voltage', ab1=-60, ab2=-85)

if (plotsave) {	
	save_dir <- str_replace(wd, 'simulations', 'images')
	setwd(save_dir)
	quartz.save(paste0('sim320f_a ', gsub(':', '-', Sys.time()), '.pdf'), type='pdf')
	dev.off()
	setwd(wd)
}

if (plotsave) {	
	save_dir <- str_replace(wd, 'simulations', 'images')
	setwd(save_dir)
	svglite(paste0('sim320f_a1 ', gsub(':', '-', Sys.time()), '.svg'), width=3.5 ,height=3.75, pointsize=10)
	funplot1(Vdend_GABA0, Vdend.plot, y.lim2, main='dendritic voltage', ab1=-30, ab2=-85)
	dev.off()
	svglite(paste0('sim320f_a2 ', gsub(':', '-', Sys.time()), '.svg'), width=3.5 ,height=3.75, pointsize=10)
	funplot1(Vsoma_GABA0, Vsoma.plot, y.lim1, main='somatic voltage', ab1=-60, ab2=-85)
	dev.off()
	setwd(wd)
}

### graph for figure ####
x.lim = c(0,600)
y.lim1 = c(-90,-60)
y.lim2 = c(-96,-30)

funplot1a <- function(V_GABA0, V, y.lim, main='dendritic voltage', ab1=-30, ab2=-85, ab3=-60){
	plot(t2-offset, V_GABA0[ind1:ind2], xlab='time (ms)', ylim=y.lim,  xlim=x.lim, ylab='mV', main=main, type="l", bty="n", col='black', axes=FALSE, lwd=lwd, lty=1, frame=FALSE) 
	axis(2, at=seq(y.lim[1], y.lim[2], y.int), labels=seq(y.lim[1], y.lim[2], y.int), las=1, tcl=tcl, lwd=lwd)
	if (xaxis) axis(1, at=seq(x.lim[1], x.lim[2], x.int), labels=seq(x.lim[1], x.lim[2], x.int) - x.lim[1], las=1, tcl=tcl, lwd=lwd)
	color_range <- fun_color_range(n.traces) 
	for (ii in 1:n.traces){
		lines(t2 + (ii-1)*50 -75, V[,ii], col=color_range[ii], lwd=lwd, lty=1)
		}
	abline(h=ab1, lty=3, lwd=lwd)
	abline(h=ab3, lty=3, lwd=lwd, col='indianred')
	abline(h=ab2, lty=3, lwd=lwd)
}
	
dev.new(width=6 ,height=10/3,noRStudioGD=TRUE)
par(mar=c(1, 1, 1, 1), mfrow=c(1,2), oma = c(2, 2, 2, 0), ps=10, cex = 0.9, cex.main = 0.9)
funplot1a(Vdend_GABA0, Vdend.plot, y.lim2, main='dendritic voltage', ab1=-30, ab2=-85, ab3=-60)
funplot1a(Vsoma_GABA0, Vsoma.plot, y.lim1, main='somatic voltage', ab1=-60, ab2=-85, ab3=-60)


if (plotsave) {	
	save_dir <- str_replace(wd, 'simulations', 'images')
	setwd(save_dir)
	quartz.save(paste0('sim320_c ', gsub(':', '-', Sys.time()), '.pdf'), type='pdf')
	dev.off()
	setwd(wd)
}

if (plotsave) {	
	save_dir <- str_replace(wd, 'simulations', 'images')
	setwd(save_dir)
	svglite(paste0('sim320_c1 ', gsub(':', '-', Sys.time()), '.svg'), width=3.5 ,height=3.75, pointsize=10)
	funplot1a(Vdend_GABA0, Vdend.plot, y.lim2, main='dendritic voltage', ab1=-30, ab2=-85, ab3=-60)
	dev.off()
	svglite(paste0('sim320_c2 ', gsub(':', '-', Sys.time()), '.svg'), width=3.5 ,height=3.75, pointsize=10)
	funplot1a(Vsoma_GABA0, Vsoma.plot, y.lim1, main='somatic voltage', ab1=-60, ab2=-85, ab3=-60)
	dev.off()
	setwd(wd)
}

#### Analysis ####
analyse <- function(V, V_GABA0){
	ind <- which.max(V_GABA0)
	Vpeak <- V_GABA0[ind] - mean(V_GABA0[((stim.time - 10)/dt):(stim.time/dt)])

	n <- dim(V)[2]
	dtpeak <- ind*dt 
	idx.peak <- dtpeak/dt 
	idx.window <- c(idx.peak - 0.5/dt, idx.peak + 0.5/dt) 
	idx <- c((stim.time-10)/dt, stim.time/dt) 
	
	baseline0 <- mean(V_GABA0[idx[1]:idx[2]])
	
	# indices for baseline
	if (n == 1){
		baseline <- mean(V[idx[1]:idx[2]])
	}else{
		baseline <- apply(V[idx[1]:idx[2],], 2, mean)
	}

	peaks  <- sapply(1:dim(V)[2], function(ii) max(V[,ii][idx.window[1]:idx.window[2]]) - baseline[ii])
	abs.peaks  <- sapply(1:dim(V)[2], function(ii) max(V[,ii]) - baseline[ii])
	
	rel.peaks = peaks / Vpeak
	rel.abs.peaks = abs.peaks / Vpeak
	list(Vpeak0 = Vpeak, baseline0=baseline0, peaks=peaks, rel.peaks=rel.peaks, abs.peaks=abs.peaks, rel.abs.peaks=rel.abs.peaks, baseline=baseline)
}

out1 <- analyse(V=Vdend, V_GABA0=Vdend_GABA0)
dV <- c(0, out1$baseline - out1$baseline0)


######
# functions for graphs
log10Tck <- function(side, type){
   lim <- switch(side, 
     x = par('usr')[1:2],
     y = par('usr')[3:4],
     stop("side argument must be 'x' or 'y'"))
   at <- floor(lim[1]) : ceiling(lim[2])
   return(switch(type, 
     minor = outer(1:9, 10^(min(at):max(at))),
     major = 10^at,
     stop("type argument must be 'major' or 'minor'")
   ))
}

roundUp <- function(x) 10^ceiling(log10(x))
roundDown <- function(x) 10^floor(log10(x))

# make graph
funplot2 = function(X,Y, y.lim, x.lim, y.int, color_range='black', main=''){
	plot(X, Y, log="x", pch=20, xlim=x.lim, ylim=y.lim, xlab='',ylab='', main=main, bty="n", col=color_range, axes=FALSE, lwd=lwd, frame=FALSE, cex=0.5)
	axis(2, at=seq(y.lim[1], y.lim[2], y.int), labels=seq(y.lim[1], y.lim[2], y.int), las=1, tcl=-0.3, lwd=lwd)
	axis(1, at=log10Tck('x','major'), tcl= -0.3) # bottom
	axis(1, at=log10Tck('x','minor'), tcl= -0.1, labels=NA) # bottom
}


X = g_gaba[2:length(g_gaba)]
Y1 = out1$rel.peaks
x.lim <- c( roundDown(min(X)), roundUp(max(X)) )
y.lim1=c(0,1.2)	
y.int1 = 0.2

Y2 = out1$baseline - out1$baseline0
y.lim2=c(0,25)	
y.int2 = 5

Y3 <- out1$rel.abs.peaks
y.lim3=c(0,2)	
y.int3 = 0.5

dev.new(width=12 ,height=4,noRStudioGD=TRUE)
par(mar=c(1, 1, 1, 1), mfrow=c(1,2), oma = c(2, 2, 2, 0), ps=10, cex = 0.9, cex.main = 0.9)
funplot2(X=X, Y=Y1, y.lim=y.lim1, x.lim=x.lim, y.int=y.int1, main='dendritic P2/P1')
funplot2(X=X, Y=Y2, y.lim=y.lim2, y.int=y.int2, x.lim=x.lim, main='dendritic depolarisation (mV)')

if (plotsave) {	
	save_dir <- str_replace(wd, 'simulations', 'images')
	setwd(save_dir)
	quartz.save(paste0('sim320f_b1 ', gsub(':', '-', Sys.time()), '.pdf'), type='pdf')
	dev.off()
	setwd(wd)
}



dev.new(width=12 ,height=4,noRStudioGD=TRUE)
par(mar=c(1, 1, 1, 1), mfrow=c(1,2), oma = c(2, 2, 2, 0), ps=10, cex = 0.9, cex.main = 0.9)
funplot2(X=X, Y=Y1, y.lim=y.lim3, x.lim=x.lim, y.int=y.int3, main='dendritic P2/P1')
funplot2(X=X, Y=Y3, y.lim=y.lim3, x.lim=x.lim, y.int=y.int3, main='dendritic P3/P1')

if (plotsave) {	
	save_dir <- str_replace(wd, 'simulations', 'images')
	setwd(save_dir)
	quartz.save(paste0('sim320f_b2 ', gsub(':', '-', Sys.time()), '.pdf'), type='pdf')
	dev.off()
	setwd(wd)
}

if (plotsave) {	
	save_dir <- str_replace(wd, 'simulations', 'images')
	setwd(save_dir)
		
	svglite(paste0('sim320f_b3 ', gsub(':', '-', Sys.time()), '.svg'), width=3.625 ,height=2.875, pointsize=10)
	funplot2(X=X, Y=Y1, y.lim=y.lim3, x.lim=x.lim, y.int=y.int3, main='dendritic P2/P1')
	dev.off()

	svglite(paste0('sim320f_b4 ', gsub(':', '-', Sys.time()), '.svg'), width=3.625 ,height=2.875, pointsize=10)
	funplot2(X=X, Y=Y3, y.lim=y.lim3, x.lim=x.lim, y.int=y.int3, main='dendritic P3/P1')
	dev.off()
	
	setwd(wd)
}

# Save the processed data as CSV to image directory if saving figures
if (plotsave) {	
    # Determine the directory to processed data
    save_dir <- str_replace(wd, 'simulations', 'images')
    setwd(save_dir)  # Set working directory to save_dir
	data <- cbind(gGABA=X, 'dendritic P3/P1'=Y3, 'dendritic depolarisation'=Y2)
	# Write to a CSV file
	write.csv(data, file = paste0('sim', sim, ' ', gsub(':', '-', Sys.time()), '.csv'), row.names = FALSE)
	# Return to the original working directory
    setwd(wd)
}

