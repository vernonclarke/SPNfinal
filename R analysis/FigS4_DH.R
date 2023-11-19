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
sim <- 312
model <- 1
nGlut <- 15

path <- paste0('/Documents/Repositories/SPNfinal/', spn, '/model', model, '/physiological/simulations/sim', sim,'/')

# Set working directory
cd <- sub("/Documents.*", "", getwd())
wd <- paste0(cd, path)
setwd(wd)

# Extract folder data
folders1 <- dir()

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
Vsoma_Glut0 <- Vsoma[,1]
Vdend_Glut0 <- Vdend[,1]
Vsoma_GABA0 <- Vsoma[,2]
Vdend_GABA0 <- Vdend[,2]


gt <- gt[3:length(gt)]
Vsoma <- Vsoma[,3:dim(Vsoma)[2]]
Vdend <- Vdend[,3:dim(Vdend)[2]]

# sort columns by gt
ind <- sort(gt, index.return=TRUE)$ix
gt <- gt[ind]
Vsoma <- Vsoma[,ind]
Vdend <- Vdend[,ind]


# Clean environment
to_keep <- c('sim', 'n.folders', 't', 'm', 'gt','Vdend','Vsoma','Vsoma_GABA0','Vdend_Glut0','Vsoma_Glut0','Vdend_GABA0', 'Vtree.path', 'plotsave', 'wd', 'nGlut')
to_remove <- setdiff(ls(), to_keep)
rm(list = to_remove)

dt <- t[2]-t[1]
fun_color_range <- colorRampPalette(c('slateblue', 'indianred', 'aquamarine4', 'darkorchid', 'sienna1'))   # Apply colorRampPalette

timing <- gt
x.lim <- c(min(timing)-10,min(timing)+190)

# Final plot:
stim.time <- 150
baseline <- 25
offset <- stim.time - baseline
t.max <- 325 # max x lim value

t1 <- t[t < t.max]

ind1 <- offset/dt
ind2 <- length(t1)
t2 <- t[ind1:ind2] 

plot.timing <- c(140, 150, 160, 170, 180) 
ind <- gt %in% plot.timing
n.traces <- length(plot.timing)

Vdend.plot <- Vdend[ind1:ind2,ind]
Vsoma.plot <- Vsoma[ind1:ind2,ind]
x.lim <- c(0,t.max)
y.lim1 <- c(-90,-60)
y.lim2 <- c(-96,-30)
lwd <- 1.5
tcl <- -0.3
lwd <- 1
xaxis <- TRUE
y.int <- 10
x.int <- 50

main <- 'dendritic voltage'
ab1 <- -30
ab2 <- -85

funplot1 = function(V, V_Glut0, V_GABA0, main='', y.lim, x.lim, y.int, x.int, ab1, ab2){
	plot(t2, V_Glut0[ind1:ind2], xlab='time (ms)', ylim=y.lim, xlim=x.lim, ylab='mV', main=main, type="l", bty="n", col='gray', axes=FALSE, lwd=lwd, lty=1, frame=FALSE) 
	axis(2, at=seq(y.lim[1], y.lim[2], y.int), labels=seq(y.lim[1], y.lim[2], y.int), las=1, tcl=tcl, lwd=lwd)
	if (xaxis) axis(1, at=seq(x.lim[1], x.lim[2], x.int), labels=seq(x.lim[1], x.lim[2], x.int) - x.lim[1], las=1, 	tcl=tcl, lwd=lwd)
	lines(t2-offset, V_GABA0[ind1:ind2], col='black', lwd=lwd)
	color_range <- fun_color_range(n.traces) 
	for (ii in 1:n.traces){
		lines(t2, V[,ii], col=color_range[ii], lwd=lwd, lty=1)
		}
	abline(h=ab1, lty=3, lwd=lwd)
	abline(h=ab2, lty=3, lwd=lwd)
}

	
dev.new(width=6 ,height=10/3,noRStudioGD=TRUE)
par(mar=c(1, 1, 1, 1), mfrow=c(1,2), oma = c(2, 2, 2, 0), ps=10, cex = 0.9, cex.main = 0.9)
funplot1(Vdend.plot, Vdend_Glut0, Vdend_GABA0, main='dendritic voltage', y.lim=y.lim2, x.lim=x.lim, y.int=10, x.int=50, ab1=-30, ab2=-85)
funplot1(Vsoma.plot, Vsoma_Glut0, Vsoma_GABA0, main='somatic voltage', y.lim=y.lim1, x.lim=x.lim, y.int=5, x.int=50, ab1=-60, ab2=-85)

if (plotsave) {	
	save_dir <- str_replace(wd, 'simulations', 'images')
	setwd(save_dir)
	quartz.save(paste0('sim11 ', gsub(':', '-', Sys.time()), '.pdf'), type='pdf')
	dev.off()
	setwd(wd)
}

if (plotsave) {	
	save_dir <- str_replace(wd, 'simulations', 'images')
	setwd(save_dir)
	svglite(paste0('sim111 ', gsub(':', '-', Sys.time()), '.svg'), width=3.5 ,height=3.75, pointsize=10)
	funplot1(Vdend.plot, Vdend_Glut0, Vdend_GABA0, main='dendritic voltage', y.lim=y.lim2, x.lim=x.lim, y.int=10, x.int=50, ab1=-30, ab2=-85)
	dev.off()
	svglite(paste0('sim112 ', gsub(':', '-', Sys.time()), '.svg'), width=3.5 ,height=3.75, pointsize=10)
	funplot1(Vsoma.plot, Vsoma_Glut0, Vsoma_GABA0, main='somatic voltage', y.lim=y.lim1, x.lim=x.lim, y.int=5, x.int=50, ab1=-60, ab2=-85)
	dev.off()
	setwd(wd)
}


# timing vs P2/P1 & P3/P1

grph2bw <- function(XX, rel.x, x.lim, y.lim, x.lab='relative timing (ms)', y.lab='relative PSP amplitude', main='', tcl=-0.3, lwd=1, xaxis=TRUE, yaxis=TRUE, y.int=1, x.int=10, abl=TRUE){
	n <- length(rel.x)
	color_range <- 'black'
	plot(rel.x, XX, pch=20, xlim=x.lim, ylim=y.lim, , xlab='',ylab='', bty="n", col=color_range, axes=FALSE, lwd=lwd, frame=FALSE, cex=0.5)
	if (xaxis){
		axis(1, at=seq(x.lim[1], x.lim[2], x.int), labels=seq(x.lim[1], x.lim[2], x.int), las=1, tcl=tcl, lwd=lwd)
		mtext(text = x.lab, side = 1, line = 2)
	}
	if (yaxis){
		axis(2, at=seq(y.lim[1], y.lim[2], y.int), labels=seq(y.lim[1], y.lim[2], y.int), las=1, tcl=tcl, lwd=lwd)
		mtext(text = y.lab, side = 2, line = 3)
	}
	if (abl) abline(h=1, col='gray', lty=3, lwd=lwd)
	
}


vary_GABA=FALSE

# Function to normalize the data 
# finds absolute max relative to either baseline (for P3) or relative to GABA envelope (ie P2)
normalise.fun3 <- function(V, V_Glut0, V_GABA0, timing, stim.time, dt){ ## XX is overall matrix eg Vsoma; X is GABA only trace for subtraction;
	m <- dim(V)[1]
	n <- dim(V)[2]
	
	# find time to peak for V_GABA0
	ind <- which.max(V_GABA0)
	Vpeak <- V_GABA0[ind] - mean(V_GABA0[((stim.time - 10)/dt):(stim.time/dt)])
	idx <- if (!vary_GABA) c((min(timing)-10)/dt, min(timing)/dt) else c((stim.time-10)/dt, stim.time/dt) 

	# indices for baseline
	if (n == 1){
		baseline <- mean(V[idx[1]:idx[2]])
	}else{
		baseline <- apply(V[idx[1]:idx[2],], 2, mean)
	}
	time2 <- (1:1:(m-idx[1]+1))*dt - dt
	time.range <- if (!vary_GABA)  timing else rep(stim.time, n)

	idxs <- sapply(1:dim(V)[2], function(ii) which.max(V[,ii]))
	# abs.peaks <- sapply(1:dim(V)[2], function(ii) max(V[,ii])) - baseline	
	abs.peaks <- sapply(1:dim(V)[2], function(ii) V[idxs[ii],ii]) - baseline
	
	abs.peaks.rel.gaba <- sapply(1:dim(V)[2], function(ii) V[idxs[ii],ii] - V_Glut0[idxs[ii]])

	list(abs.peaks=abs.peaks, rel.abs.peaks=abs.peaks.rel.gaba, Vpeak = Vpeak)
}	

# abs.peaks are absolute relative to baseline
# rel.abs.peaks are relative to underlying GABA response

out <- normalise.fun3(V=Vdend, V_Glut0=Vdend_Glut0, V_GABA0=Vdend_GABA0, timing=gt, stim.time, dt)

Vdend.absolute <- out$abs.peaks
Vdend.absolute.GABA <- out$rel.abs.peaks
Vpeak <- out$Vpeak 

Vdend.absolute.rel <- Vdend.absolute/Vpeak
Vdend.absolute.GABA <- Vdend.absolute.GABA/Vpeak


offset <- 150
x.lim <- c(120,210) - offset
timing <- gt - offset
y.lim <- c(0,2)
x.lab <- 'relative Glut timing (ms)'
y.lab <- 'relative potential amplitude'


# Set up for plots
dev.new(width=12 ,height=4,noRStudioGD=TRUE)
par(mar=c(1, 1, 1, 1), mfrow=c(1,2), oma = c(2, 2, 2, 0), ps=10, cex = 0.9, cex.main = 0.9)
# Generate plots using the grph2bw function
grph2bw(Vdend.absolute.GABA, timing, x.lim=x.lim , y.lim=y.lim, x.lab=x.lab, y.lab=y.lab, main='', tcl=-0.3, lwd=1, xaxis=TRUE, yaxis=TRUE, y.int=0.5, x.int=10)
grph2bw(Vdend.absolute.rel, timing, x.lim=x.lim, y.lim=y.lim, x.lab=x.lab, y.lab=y.lab, main='', tcl=-0.3, lwd=1, xaxis=TRUE, yaxis=FALSE, y.int=0.5, x.int=10)

# Save the plot as PDF if the 'plotsave' flag is TRUE
if (plotsave) {	
	save_dir <- str_replace(wd, 'simulations', 'images')
	setwd(save_dir)
	quartz.save(paste0('sim11c ', gsub(':', '-', Sys.time()), '.pdf'), type='pdf')
	setwd(wd)
}

# Save the plot as SVG if the 'plotsave' flag is TRUE
if (plotsave) {	
	save_dir <- str_replace(wd, 'simulations', 'images')
	setwd(save_dir)
	svglite(paste0('sim11c1 ', gsub(':', '-', Sys.time()), '.svg'), width=3.625 ,height=2.875, pointsize=10)
	grph2bw(Vdend.absolute.GABA, timing, x.lim=x.lim , y.lim=y.lim, x.lab=x.lab, y.lab=y.lab, main='', tcl=-0.3, lwd=1, xaxis=TRUE, yaxis=TRUE, y.int=0.5, x.int=10)

	dev.off()
	svglite(paste0('sim11c2 ', gsub(':', '-', Sys.time()), '.svg'), width=3.625 ,height=2.875, pointsize=10)
	grph2bw(Vdend.absolute.rel, timing, x.lim=x.lim, y.lim=y.lim, x.lab=x.lab, y.lab=y.lab, main='', tcl=-0.3, lwd=1, xaxis=TRUE, yaxis=FALSE, y.int=0.5, x.int=10)
	dev.off()
	setwd(wd)
}

# Save the processed data as CSV to image directory if saving figures
if (plotsave) {	
    # Determine the directory to processed data
    save_dir <- str_replace(wd, 'simulations', 'images')
    setwd(save_dir)  # Set working directory to save_dir
	data <- cbind(timing=timing,Vdend.GABA=Vdend.absolute.GABA, Vdend.abs=Vdend.absolute.rel)
	# Write to a CSV file
	write.csv(data, file = paste0('sim', sim, ' ', gsub(':', '-', Sys.time()), '.csv'), row.names = FALSE)
	# Return to the original working directory
    setwd(wd)
}


