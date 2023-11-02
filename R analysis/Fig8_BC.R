rm( list=ls(all=TRUE ) )
# Load and install necessary packages
load_required_packages <- function(packages){
	new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  	if(length(new.packages)) install.packages(new.packages)
  	lapply(packages, library, character.only=TRUE)
}

required.packages <- c('reticulate', 'stringr', 'svglite', 'quantmod', 'xts', 'zoo')
load_required_packages(required.packages)

# Install miniconda if necessary
if (!reticulate::py_available(initialize = TRUE)) {
  	reticulate::install_miniconda()
}

# Import necessary python modules
pd <- reticulate::import("pandas")

# Metadata settings
plotsave <- TRUE
spn <- 'dspn'
sim <- 9
model <- 1
path <- paste0('/Documents/GitHub/SPNfinal/', spn, '/model', model, '/physiological/simulations/sim', sim,'/')

# Set working directory
cd <- sub("/Documents.*", "", getwd())
wd <- paste0(cd, path)
setwd(wd)

offsite <- FALSE

folders <- dir()
Vtree.path <- list()
Vsoma <- list()
n.folders <- length(folders)
imp <- TRUE
return.currents <- TRUE

for (iii in 1:n.folders){
	XX = folders[iii]
	setwd(XX)
	V_soma <- pd$read_pickle(paste0(XX, '_v_data.pickle'))
	V_dend <- pd$read_pickle(paste0(XX, '_vdend_data.pickle'))
	time <- pd$read_pickle(paste0(XX, '_t_data.pickle'))
	if (imp){
		Z_soma <- pd$read_pickle(paste0(XX, '_imp_soma_data.pickle'))
		Z_dend <- pd$read_pickle(paste0(XX, '_imp_dend_data.pickle'))
	}
	m <- length(time)
	if (iii==1){
		t <- time[[1]]
		gt <- rep(0, n.folders)
		n <- length(time[[1]])
		mat <- matrix(rep(0, n.folders*n), ncol=n.folders)
		Vsoma <- Vdend <- rep(list(mat), m)
		if (imp){
			Zsoma <- Zdend  <- rep(list(mat), m)
		}
	}
	gt[iii] <- pd$read_pickle(paste0(XX, '_gt.pickle'))[[1]][1]
	for (ii in 1:m){
		Vsoma[[ii]][,iii] <- V_soma[[ii]]
		Vdend[[ii]][,iii] <- V_dend[[ii]]
		if (imp){
			Zsoma[[ii]][,iii] <- Z_soma[[ii]]
			Zdend[[ii]][,iii] <- Z_dend[[ii]]
		}
	}
	
	setwd(wd)
}


if (return.currents){ 
	for (iii in 1:n.folders){
		XX = folders[iii]
		setwd(XX)
		i_mechs <- pd$read_pickle(paste0(XX, '_i_mechanisms.pickle'))
		gGABA <- pd$read_pickle(paste0(XX, '_g_gaba.pickle'))
		if (iii==1){
			mechanisms <- i_mechanisms <- g_GABA <- list()
		}
		i_mechanisms[[iii]] <- i_mechs
		g_GABA[[iii]] <- gGABA
		setwd(wd)
	}
}


if (sum(grepl('GABA0', folders))==1){ # contains 0 GABA folder
	XX = folders[grepl('GABA0', folders)]
	setwd(XX)
	Vsoma_GABA0 <- pd$read_pickle(paste0(XX, '_v_data.pickle'))
	Vdend_GABA0 <- pd$read_pickle(paste0(XX, '_vdend_data.pickle'))
	if (imp){
		Zsoma_GABA0 <- pd$read_pickle(paste0(XX, '_imp_soma_data.pickle'))
		Zdend_GABA0 <- pd$read_pickle(paste0(XX, '_imp_dend_data.pickle'))
	}
	i_mechs_GABA0 <- pd$read_pickle(paste0(XX, '_i_mechanisms.pickle'))
	setwd(wd)
}


Vsoma_GABA0 <- Vsoma_GABA0[[1]]
Vdend_GABA0 <- Vdend_GABA0[[1]]
i_mechs_GABA0 <- i_mechs_GABA0[[1]]

if (imp){
	Zsoma_GABA0 <- Zsoma_GABA0[[1]]
	Zdend_GABA0 <- Zdend_GABA0[[1]]
	}

# sort by index
nGABA <- sapply(1:length(folders), function(ii) as.numeric(gsub('.*GABA(.+)_GABA.*', '\\1', folders[ii])) )
ind <- sort(nGABA, index.return=TRUE)$ix
nGABA <- nGABA[ind]
gt <- unlist(gt[ind])

Vsoma <- lapply(1:m, function(ii) Vsoma[[ii]][,ind])
Vdend <- lapply(1:m, function(ii) Vdend[[ii]][,ind])
if (imp){
	Zsoma <- lapply(1:m, function(ii) Zsoma[[ii]][,ind])
	Zdend <- lapply(1:m, function(ii) Zdend[[ii]][,ind])
}

if (return.currents) {
	i_mechanisms <- i_mechanisms[ind]
	mechanisms <- mechanisms[ind]
}
# in these examples m = 1 so reduce Vsoma from list to matrix:
if (m==1){
	Vsoma <- Vsoma[[1]]
	Vdend <- Vdend[[1]]
	if (imp){
		Zsoma <- Zsoma[[1]]
		Zdend <- Zdend[[1]]
	}
}

to_keep <- c('imp', 'g_GABA', 'i_GABA', 'n.folders', 't', 'm', 'gt', 'i_mechanisms', 'mechanisms', 'Vdend','Vsoma', 'Zdend','Zsoma', 'mechs_GABA0', 'i_mechs_GABA0', 'Vdend_GABA0', 'Vsoma_GABA0', 'Vdend_Glut0', 'Vsoma_Glut0', 'Zdend_GABA0', 'Zsoma_GABA0', 'Zdend_Glut0', 'Zsoma_Glut0', 'nGABA', 'offsite', 'wd', 'plotsave', 'sim')

rm(list=ls()[!(ls() %in% to_keep)])


dt <- t[2]-t[1]
fun_color_range <- colorRampPalette(c('slateblue', 'indianred', 'aquamarine4', 'darkorchid', 'sienna1'))   # Apply colorRampPalette
stim.time <- 150
timing <- gt

######## FUNCTIONS ########
grph1a <- function(mat, t, x.lim, y.lim, x.offsets=NULL, iii=1, main='', tcl=-0.3, lwd=1, xaxis=TRUE, y.int=10, x.int=50){
	n.folders <- dim(mat)[2]
	color_range <- fun_color_range(n.folders) 
	ind <- sapply(1:length(x.lim), function(ii) which.min(abs(t - x.lim[ii])) )
	if (xaxis) x.lab='time (ms)' else x.lab=''
	if (is.null(x.offsets)) x.offsets <- rep(0, n.folders)
	for (ii in 1:n.folders){
		if (ii==1){
			plot(t[ind[1]:ind[2]]+ x.offsets[ii], mat[ind[1]:ind[2],ii], xlab=x.lab, ylim=y.lim,  xlim=x.lim, ylab='mV', main=main, type='l', bty='n', col=color_range[ii], axes=FALSE, lwd=lwd, lty=1, frame=FALSE) 
			axis(2, at=seq(y.lim[1], y.lim[2], y.int), labels=seq(y.lim[1], y.lim[2], y.int), las=1, tcl=tcl, lwd=lwd)
		}else{
			lines(t[ind[1]:ind[2]]+ x.offsets[ii], mat[ind[1]:ind[2],ii], col=color_range[ii], lwd=lwd, lty=1)
		}
	}
	if (xaxis) axis(1, at=seq(x.lim[1], x.lim[2], x.int), labels=seq(x.lim[1], x.lim[2], x.int) - x.lim[1], las=1, tcl=tcl, lwd=lwd)
}

grph1 <- function(mat, time, x.lim, y.lim, y.lab='PSP (mV)', main='', tcl=-0.3, lwd=1, xaxis=TRUE, x.int=25, y.int=1, offset=NULL){
	n <- dim(mat)[2]
	if (is.null(offset)) offset <- rep(0, n-1)
	color_range <- fun_color_range(n)  
	if (xaxis) x.lab='time (ms)' else x.lab=''

	for (ii in 1:n){
		if (ii==1){
			plot(time + offset[ii], mat[,ii], xlab=x.lab, ylim=y.lim,  xlim=x.lim, ylab=y.lab, main=main, type='l', bty='n', col=color_range[ii], axes=FALSE, lwd=lwd, lty=1, frame=FALSE) 
			axis(2, at=seq(y.lim[1], y.lim[2], y.int), labels=seq(y.lim[1], y.lim[2], y.int), las=1, tcl=tcl, lwd=lwd)
			mtext(text = 'PSP (mV)', side = 2, line = 2)
		}else{
			lines(time + offset[ii], mat[,ii], col=color_range[ii], lwd=lwd, lty=1)
		}
	}
	if (xaxis){
		axis(1, at=seq(x.lim[1], x.lim[2], x.int), labels=seq(x.lim[1], x.lim[2], x.int) - x.lim[1], las=1, tcl=tcl, lwd=lwd)
		mtext(x.lab, side=1, line=2,  outer = FALSE)
	}
}

grph2 <- function(XX, rel.x, x.lim, y.lim, x.lab='relative timing (ms)', y.lab='relative PSP amplitude', main='', tcl=-0.3, lwd=1, xaxis=TRUE, yaxis=TRUE, y.int=1, x.int=10){
	n <- length(rel.x)
	color_range <- fun_color_range(n) 
	plot(rel.x, XX, pch=20, xlim=x.lim, ylim=y.lim, , xlab='',ylab='', bty='n', col=color_range, axes=FALSE, lwd=lwd, frame=FALSE)
	if (xaxis){
		axis(1, at=seq(x.lim[1], x.lim[2], x.int), labels=seq(x.lim[1], x.lim[2], x.int), las=1, tcl=tcl, lwd=lwd)
		mtext(text = x.lab, side = 1, line = 2)
	}
	if (yaxis){
		axis(2, at=seq(y.lim[1], y.lim[2], y.int), labels=seq(y.lim[1], y.lim[2], y.int), las=1, tcl=tcl, lwd=lwd)
		mtext(text = y.lab, side = 2, line = 3)
	}
	abline(h=1, col='gray', lty=3, lwd=lwd)	
}

grph2bw <- function(XX, rel.x, x.lim, y.lim, x.lab='relative timing (ms)', y.lab='relative PSP amplitude', main='', tcl=-0.3, lwd=1, col = 'slateblue', xaxis=TRUE, yaxis=TRUE, y.int=1, x.int=10, abl=TRUE){
	n <- length(rel.x)
	color_range <- col
	plot(rel.x, XX, pch=20, xlim=x.lim, ylim=y.lim, , xlab='',ylab='', bty='n', col=color_range, axes=FALSE, lwd=lwd, frame=FALSE)
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

# Define the exp.fit1 function to model an exponential curve and compute model selection criteria
exp.fit1 <- function(x,y){  
  	# Initial parameter estimation
  	theta.0 <- max(y) * 1.1
  	model.0 <- lm(log(- y + theta.0) ~ x)
  	alpha.0 <- -exp(coef(model.0)[1])
  	beta.0 <- coef(model.0)[2]
  
  	# Provide starting values for the non-linear model
  	start <- list(alpha = alpha.0, beta = beta.0, theta = theta.0)
  
  	# Fit the non-linear model
  	model <- nls(y ~ alpha * exp(beta * x) + theta, start = start)
  
  	# Calculate the residual sum of squares (RSS)
  	RSS <-  sum(residuals(model)^2) 
  
  	# Calculate the number of data points and the number of coefficients
  	n <- length(y)
  	k <- length(coef(model))
  
  	# Compute the log-likelihood
  	logL <- 0.5 * (- n * (log(2 * pi) + 1 - log(n) + log(RSS)))
  
  	# Compute Akaike's Information Criterion (AIC)
  	aic <- 2*(k+1)-2 * logL 
  
  	# Compute Bayesian Information Criterion (BIC)
  	bic <- 2*(k-1)*log(n) - 2*logL
  
  	# Create a sequence of x values for prediction
  	xfit <- seq(min(x), max(x), length.out=1000)
  	
  	# Predict the y values based on the model
  	yfit <- predict(model, list(x = xfit))
  
  	# Return the predictions, coefficients, and model selection criteria
  	list(fits=cbind(xfit=xfit,yfit=yfit), coeffs=coef(model), AIC=aic, BIC=bic, logL=logL)
}

# Define the exp.fit2 function to model an exponential curve and compute model selection criteria
exp.fit2 <- function(x,y){
	# Initial parameter estimation
	theta.0 <- min(y) * 0.5 
	model.0 <- lm(log(y - theta.0) ~ x)  
	alpha.0 <- exp(coef(model.0)[1])
	beta.0 <- coef(model.0)[2]
	
	# Provide starting values for the non-linear model
  	start <- list(alpha = alpha.0, beta = beta.0, theta = theta.0)
	
	# Fit the model
	model <- nls(y ~ alpha * exp(beta * x) + theta, start = start)
	
	# calculate model selection criteria
	RSS <-  sum(residuals(model)^2) 
	
	# Calculate the number of data points and the number of coefficients
  	n <- length(y)
	k <- length(coef(model))
	
	# Compute the log-likelihood
  	logL <- 0.5 * (- n * (log(2 * pi) + 1 - log(n) + log(RSS)))
	
	# Compute Akaike's Information Criterion (AIC)
  	aic <- 2*(k+1)-2 * logL 
	
	# Compute Bayesian Information Criterion (BIC)
  	bic <- 2*(k-1)*log(n) - 2*logL
	
	# Create a sequence of x values for prediction
  	xfit <- seq(min(x), max(x), length.out=1000)
	
	# Predict the y values based on the model
  	yfit <- predict(model, list(x = xfit))
	
	# Return the predictions, coefficients, and model selection criteria
  	list(fits=cbind(xfit=xfit,yfit=yfit), coeffs=coef(model), AIC=aic, BIC=bic, logL=logL)
	}
	
linear.fit <- function(x,y){
	# fit linear model
	model <- lm(y ~ x)
	
	# calculate model selection criteria
	RSS <-  sum(residuals(model)^2) 
	
	# Calculate the number of data points and the number of coefficients
  	n <- length(y)
	k <- length(coef(model))
	
	# Compute the log-likelihood
  	logL <- 0.5 * (- n * (log(2 * pi) + 1 - log(n) + log(RSS)))
	
	# Compute Akaike's Information Criterion (AIC)
  	aic <- 2*(k+1)-2 * logL 
	
	# Compute Bayesian Information Criterion (BIC)
  	bic <- 2*(k-1)*log(n) - 2*logL
	
	# Create a sequence of x values for prediction
  	xfit <- seq(min(x), max(x), length.out=1000)
	
	# Predict the y values based on the model
  	yfit <- predict(model, list(x = xfit))
	
	# Return the predictions, coefficients, and model selection criteria
  	list(fits=cbind(xfit=xfit,yfit=yfit), coeffs=coef(model), AIC=aic, BIC=bic, logL=logL)
}



# wrapper for impedance plots
fun.imp.plot <- function(V, Z, t, x.lim, y.lim1, y.int1, y.lim2, y.int2){
	plot(t, Z, type='l', xlim=x.lim, ylim =y.lim1, col='indianred', bty='n', lwd=1, lty=1, axes=FALSE, frame=FALSE)
	axis(2, yaxp=c(y.lim1[1], y.lim1[2], y.int1), las=1, tck= -0.03)
	par(new=TRUE)
	plot(t, V, type='l', xlim=x.lim, ylim =y.lim2, col='slateblue', bty='n', lwd=1, lty=1, axes=FALSE, frame=FALSE)
	axis(4, yaxp=c(y.lim2[1], y.lim2[2], y.int2), las=1, tck= -0.03)
	mtext('', side = 4, line = 3)
}

# finds min impedance / peak PSP amplitudes etc
normalise.fun <- function(V, stim.time, dt){ ## XX is overall matrix eg Vsoma; X is GABA only trace for subtraction;
	m <- dim(V)[1]
	n <- dim(V)[2]

	# find time to peak for V_GABA0
	idx <-  c((min(stim.time)-10)/dt, min(stim.time)/dt) # indices for baseline
	if (n == 1){
		baseline = mean(V[idx[1]:idx[2]])
	}else{
		baseline <- apply(V[idx[1]:idx[2],], 2, mean)
	}
	time2 = (1:1:(m-idx[1]+1))*dt - dt
	ind <-  c(stim.time/dt, (stim.time+100)/dt)
	peaks <- sapply(1:length(timing), function(ii){
		dV <- V[,ii] - baseline[ii]
		dV <- dV[ind[1]:ind[2]]
		idx.peak <- which.max(dV)
		# expected glutamate peak
		V.peak  <- dV[idx.peak]
		peak.time <-idx.peak*dt
		c(V.peak, peak.time)	
	})
	list(V = V[idx[1]:m,], time=time2, peaks=peaks[1,], 'time2peak'=peaks[2,])
}	

Z.GABA <- function(Z, stim.time, dt){
	m <- dim(Z)[1]
	n <- dim(Z)[2]
	
	# find time to peak for V_GABA0
	idx <-  c((min(stim.time)-10)/dt, min(stim.time)/dt) # indices for baseline
	if (n == 1){
		baseline = mean(Z[idx[1]:idx[2]])
	}else{
		baseline <- apply(Z[idx[1]:idx[2],], 2, mean)
	}
	time2 = (1:1:(m-idx[1]+1))*dt - dt
	ind <-  c(stim.time/dt, (stim.time+100)/dt)
	troughs <- sapply(1:length(timing), function(ii){
		dZ <- Z[,ii] 
		dZ <- dZ[ind[1]:ind[2]]
		idx.trough <- which.min(dZ)
		dZ.trough  <- dZ[idx.trough] 
		rel.change <- (baseline[ii]-dZ.trough) / baseline[ii]
		min.time <-idx.trough*dt
		c(dZ.trough, rel.change, min.time)	
	})
	list(Z = Z[idx[1]:m,], Zbaseline=baseline, Zmin=troughs[1,], rel.change=troughs[2,], 'time2min'=troughs[3,])
}
#########################################################################################################################


# do analysis
out <- normalise.fun(Vsoma, stim.time, dt)
Vsoma2 <- out$V	
time2 <- out$time
Vsoma.peak <- out$peaks	# GABA peak
Vsoma.peak.time <- out$'time2peak' 
Vsoma.peak.time[nGABA==0] <- NA

out <- normalise.fun(Vdend, stim.time, dt)
Vdend2 <- out$V
time2 <- out$time
Vdend.peak <- out$peaks	# GABA peak
Vdend.peak.time <- out$'time2peak' 
Vdend.peak.time[nGABA==0] <- NA


out <- Z.GABA(Zsoma, stim.time,dt)
Zsoma2 <- out$Z
Zsoma.baseline <- out$Zbaseline
Zsoma.min <- out$Zmin
Zsoma.rel <- out$rel.change
Zsoma.min.time <- out$time2min
Zsoma.min.time[nGABA==0] <- NA
dZsoma <- Zsoma.min  - Zsoma.baseline

out <- Z.GABA(Zdend, stim.time,dt)
Zdend2 <- out$Z
Zdend.baseline<- out$Zbaseline
Zdend.min <- out$Zmin
Zdend.rel <- out$rel.change
Zdend.min.time <- out$time2min
Zdend.min.time[nGABA==0] <- NA
dZdend <- Zdend.min  - Zdend.baseline

# generate plots of GABAergic PSPs
# choose subset for main plot
# want timings 170, 190, 200, 210, 220, 230
plot.timing <- c(1, 3, 6, 12, 24)
ind <- nGABA %in% plot.timing

Vdend.plot <- Vdend2[, ind]
Vsoma.plot <- Vsoma2[, ind]
x.lim <- c(0,150)
y.lim1 <- c(-90,-60)
y.lim2 <- c(-90,-60)
lwd <- 1.5
	
offset = seq(0, (dim(Vdend.plot)[2]-1)*10, length.out=dim(Vdend.plot)[2])  # offset for clarity
dev.new(width=6 ,height=10/3,noRStudioGD=TRUE)
par(mar=c(1, 1, 1, 1), mfrow=c(1,2), oma = c(2, 2, 2, 0), ps=10, cex = 0.9, cex.main = 0.9)
grph1a(Vdend.plot, time2,  x.lim, y.lim2, x.offsets=offset, iii=iii, main='dendritic voltage', lwd=lwd, xaxis=TRUE, y.int=10)
abline(h=-60, lty=3, lwd=lwd)
abline(h=-85, lty=3, lwd=lwd)
	
grph1a(Vsoma.plot, time2, x.lim, y.lim1, x.offsets=offset, iii=iii, main='somatic voltage', lwd=lwd, xaxis=TRUE, y.int=10)
abline(h=-60, lty=3, lwd=lwd)
abline(h=-85, lty=3, lwd=lwd)


if (plotsave) {	
	save_dir <- str_replace(wd, 'simulations', 'images')
	setwd(save_dir)
	quartz.save(paste0('sim9a ', gsub(':', '-', Sys.time()), '.pdf'), type='pdf')
	# save as individuals for plots
	svglite(paste0('sim', sim, 'a', ' ', gsub(':', '-', Sys.time()), '.svg'), width=3.5 ,height=2.65, pointsize=10)
	grph1a(Vdend.plot, time2,  x.lim, y.lim2, x.offsets=offset, iii=iii, main='dendritic voltage', lwd=lwd, xaxis=TRUE, y.int=10)
	abline(h=-60, lty=3, lwd=lwd)
	abline(h=-85, lty=3, lwd=lwd)
	dev.off()
	svglite(paste0('sim', sim, 'b', ' ',  gsub(':', '-', Sys.time()), '.svg'), width=3.5 ,height=2.65, pointsize=10)
	grph1a(Vsoma.plot, time2, x.lim, y.lim1, x.offsets=offset, iii=iii, main='somatic voltage', lwd=lwd, xaxis=TRUE, y.int=10)
	abline(h=-60, lty=3, lwd=lwd)
	abline(h=-85, lty=3, lwd=lwd)	
	dev.off()
	setwd(wd)
}

# generate plots of GABAergic impedance
Zdend.plot <- Zdend2[, ind]
Zsoma.plot <- Zsoma2[, ind]

y.lim1 <- c(50, 300)
y.lim2 <- c(50, 80)
y.int1 <- 50
y.int2 <- 5	
lwd <- 1.5
	
dev.new(width=6 ,height=10/3,noRStudioGD=TRUE)
par(mar=c(1, 1, 1, 1), mfrow=c(1,2), oma = c(2, 2, 2, 0), ps=10, cex = 0.9, cex.main = 0.9)
grph1a(Zdend.plot, time2,  x.lim, y.lim1,  x.offsets=offset, iii=iii, main='dendritic impedance', lwd=lwd, xaxis=TRUE, y.int=y.int1)
# abline(h=-60, lty=3, lwd=lwd)
# abline(h=-85, lty=3, lwd=lwd)
	
grph1a(Zsoma.plot, time2, x.lim, y.lim2,  x.offsets=offset, iii=iii, main='somatic impedance', lwd=lwd, xaxis=TRUE, y.int=y.int2)
abline(h=-60, lty=3, lwd=lwd)
abline(h=-85, lty=3, lwd=lwd)

if (plotsave) {	
	save_dir <- str_replace(wd, 'simulations', 'images')
	setwd(save_dir)
	quartz.save(paste0('sim9b ', gsub(':', '-', Sys.time()), '.pdf'), type='pdf')
	# save as individuals for plots
	svglite(paste0('sim', sim, 'c', ' ', gsub(':', '-', Sys.time()), '.svg'), width=3.5 ,height=2.2, pointsize=10)
	grph1a(Zdend.plot, time2,  x.lim, y.lim1, x.offsets=offset, iii=iii, main='dendritic impedance', lwd=lwd, xaxis=TRUE, y.int=50)
	abline(h=50, lty=3, lwd=lwd)
	abline(h=300, lty=3, lwd=lwd)
	dev.off()
	svglite(paste0('sim', sim, 'd', ' ',  gsub(':', '-', Sys.time()), '.svg'), width=3.5 ,height=2.2, pointsize=10)
	grph1a(Zsoma.plot, time2, x.lim, y.lim2, x.offsets=offset, iii=iii, main='somatic impedance', lwd=lwd, xaxis=TRUE, y.int=10)
	abline(h=50, lty=3, lwd=lwd)
	abline(h=80, lty=3, lwd=lwd)	
	dev.off()
	setwd(wd)
	
}


# Prepare a good inital state
# plot vs relative peak time
y.lim1 <- c(0, ceiling(max(Vdend.peak/10))*10)
y.lim2 <- c(0, ceiling(max(Vsoma.peak/10))*10)
x.lim <- c(min(nGABA), ceiling(max(nGABA/5))*5) 
y.int1 <- 5
y.int2 <- 1

dev.new(width=6 ,height=10/3,noRStudioGD=TRUE)
par(mar=c(1, 1, 1, 1), mfrow=c(1,2), oma = c(2, 2, 2, 0), ps=10, cex = 0.9, cex.main = 0.9)
grph2bw(Vdend.peak, nGABA, x.lim=x.lim, y.lim=y.lim1, x.lab='# GABA synapses/dendrite', y.lab='mV', main='', tcl=-0.3, lwd=1, col='slateblue', xaxis=TRUE, yaxis=TRUE, y.int=y.int1, x.int=5, abl=F)
fits.Vdend <- exp.fit1(x=nGABA,y=Vdend.peak)
lines(fits.Vdend$fits[,1],fits.Vdend$fits[,2], col = 'slateblue', lty=3, lwd = lwd)
	
grph2bw(Vsoma.peak, nGABA, x.lim=x.lim, y.lim=y.lim2, x.lab='# GABA synapses/dendrite', y.lab='mV', main='', tcl=-0.3, lwd=1, col='slateblue', xaxis=TRUE, yaxis=FALSE, y.int=y.int2, x.int=5, abl=F)
fits.Vsoma <- exp.fit1(x=nGABA,y=Vsoma.peak)
lines(fits.Vsoma$fits[,1],fits.Vsoma$fits[,2], col = 'slateblue', lty=3, lwd = lwd)

if (plotsave) {	
	save_dir <- str_replace(wd, 'simulations', 'images')
	setwd(save_dir)
	quartz.save(paste0('sim9c ', gsub(':', '-', Sys.time()), '.pdf'), type='pdf')
	# save as individuals for plots
	svglite(paste0('sim', sim, 'e', ' ', gsub(':', '-', Sys.time()), '.svg'), width=2.65 ,height=3.4, pointsize=10)
	grph2bw(Vdend.peak, nGABA, x.lim=x.lim, y.lim=y.lim1, x.lab='# GABA synapses/dendrite', y.lab='mV', main='', tcl=-0.3, lwd=1, col='slateblue', xaxis=TRUE, yaxis=TRUE, y.int=y.int1, x.int=5, abl=F)
	lines(fits.Vdend$fits[,1],fits.Vdend$fits[,2], col = 'slateblue', lty=3, lwd = lwd)
	dev.off()
	svglite(paste0('sim', sim, 'f', ' ',  gsub(':', '-', Sys.time()), '.svg'), width=2.65 ,height=3.4, pointsize=10)
	grph2bw(Vsoma.peak, nGABA, x.lim=x.lim, y.lim=y.lim2, x.lab='# GABA synapses/dendrite', y.lab='mV', main='', tcl=-0.3, lwd=1, col='slateblue', xaxis=TRUE, yaxis=FALSE, y.int=y.int2, x.int=5, abl=F)
	fits <- exp.fit1(x=nGABA,y=Vsoma.peak)
	lines(fits.Vsoma$fits[,1],fits.Vsoma$fits[,2], col = 'slateblue', lty=3, lwd = lwd)
	dev.off()
	setwd(wd)
}


y.lim1 <- c(0, ceiling(max(Zdend.min/50))*50)
y.lim2 <- c(50, ceiling(max(Zsoma.min/5))*5)

y.lim1 <- c(0, 250)
y.lim2 <- c(50, 75)
y.int1 <- 50
y.int2 <- 5


dev.new(width=6 ,height=10/3,noRStudioGD=TRUE)
par(mar=c(1, 1, 1, 1), mfrow=c(1,2), oma = c(2, 2, 2, 0), ps=10, cex = 0.9, cex.main = 0.9)
iii= 1 # on-site distal
grph2bw(Zdend.min, nGABA, x.lim=x.lim, y.lim=y.lim1, x.lab='# GABA synapses/dendrite', y.lab='mV', main='', tcl=-0.3, lwd=1, col='indianred', xaxis=TRUE, yaxis=TRUE, y.int=y.int1, x.int=5, abl=F)	
fits.Zdend <- if (offsite) linear.fit(x=nGABA,y=Zdend.min) else exp.fit2(x=nGABA,y=Zdend.min)
lines(fits.Zdend$fits[,1],fits.Zdend$fits[,2], col = 'indianred', lty=3, lwd = lwd)
	
grph2bw(Zsoma.min, nGABA, x.lim=x.lim, y.lim=y.lim2, x.lab='# GABA synapses/dendrite', y.lab='mV', main='', tcl=-0.3, lwd=1, col='indianred', xaxis=TRUE, yaxis=TRUE, y.int=y.int2, x.int=5, abl=F)
fits.Zsoma <- exp.fit2(x=nGABA,y=Zsoma.min)
lines(fits.Zsoma$fits[,1],fits.Zsoma$fits[,2], col = 'indianred', lty=3, lwd = lwd)

if (plotsave) {	
	save_dir <- str_replace(wd, 'simulations', 'images')
	setwd(save_dir)
	quartz.save(paste0('sim9d ', gsub(':', '-', Sys.time()), '.pdf'), type='pdf')
	
	# save as individuals for plots
	svglite(paste0('sim', sim, 'g', ' ', gsub(':', '-', Sys.time()), '.svg'), width=2.65 ,height=3.4, pointsize=10)
	grph2bw(Zdend.min, nGABA, x.lim=x.lim, y.lim=y.lim1, x.lab='# GABA synapses/dendrite', y.lab='mV', main='', tcl=-0.3, lwd=1, col='indianred', xaxis=TRUE, yaxis=TRUE, y.int=y.int1, x.int=5, abl=F)	
	lines(fits.Zdend$fits[,1],fits.Zdend$fits[,2], col = 'indianred', lty=3, lwd = lwd)
	dev.off()
	
	svglite(paste0('sim', sim, 'h', ' ',  gsub(':', '-', Sys.time()), '.svg'), width=2.65 ,height=3.4, pointsize=10)
	grph2bw(Zsoma.min, nGABA, x.lim=x.lim, y.lim=y.lim2, x.lab='# GABA synapses/dendrite', y.lab='mV', main='', tcl=-0.3, lwd=1, col='indianred', xaxis=TRUE, yaxis=TRUE, y.int=y.int2, x.int=5, abl=F)
	lines(fits.Zsoma$fits[,1],fits.Zsoma$fits[,2], col = 'indianred', lty=3, lwd = lwd)
	dev.off()

	setwd(wd)
}

# Save the processed data as CSV to the image directory if saving figures
if (plotsave) {	
	# Determine the directory to save processed data
	save_dir <- str_replace(wd, 'simulations', 'images')
	setwd(save_dir)  # Set working directory to save_dir
    
	# Prepare data for export
	data1 <- cbind(nGABA = nGABA, PSP.dend = Vdend.peak, PSP.soma = Vsoma.peak)
	data2 <- cbind(nGABA = nGABA, Z.dend = Zdend.min, Z.soma = Zsoma.min)
	data3 <- cbind(nGABA.fits = fits.Vdend$fits[, 1], PSP.dend.fits = fits.Vdend$fits[, 2], PSP.soma.fits = fits.Vsoma$fits[, 2])
	data4 <- cbind(nGABA.fits = fits.Zdend$fits[, 1], PSP.dend.fits = fits.Zdend$fits[, 2], PSP.soma.fits = fits.Zsoma$fits[, 2])
    
	# Generate timestamp for file names
	timestamp <- gsub(':', '-', Sys.time())
    
	# Define a function for writing data to CSV
	write_data_to_csv <- function(data, index) {
		filename <- paste0('sim', sim, '_', index, ' ', timestamp, '.csv')
		write.csv(data, file = filename, row.names = FALSE)
	}
    
	# Write data to CSV files
	write_data_to_csv(data1, '1')
	write_data_to_csv(data2, '2')
	write_data_to_csv(data3, '3')
	write_data_to_csv(data4, '4')
    
	# Return to the original working directory
	setwd(wd)
}


