###################### Multi-session SCR on polecat data 2021 #########################

rm(list=ls())
dir <- "D:/OneDrive - Sveriges lantbruksuniversitet/R/Polecam/"
setwd(dir)

## Load packages
library(oSCR)      # Handling SCR data and models
library(sf)        # Spatial data manipulation
library(ggplot2)   # Improved plotting

#=====================================================================================
## Step 1: Read in capture history & trap detection files ----

# enter the detection history and trap detection:

# detection history - ALL SESSIONS TOGETHER IN ONE FILE.
pc_edf <- read.csv2("polecat.edf.csv")

# check data
head(pc_edf)
str(pc_edf)

# trap array
ba.tdf <- read.csv2("tdf.baldringe.csv",header = TRUE) #baldringe autumn
ba.tdf2 <- read.csv2("tdf.baldringe2.csv",header = TRUE) #baldringe spring
ch.tdf <- read.csv2("tdf.christinehof.csv",header = TRUE) #christinehof autumn
ho.tdf <- read.csv2("tdf.hogestad.csv",header = TRUE) #hogerstad spring
vi.tdf <- read.csv2("tdf.vitemolla.csv",header = TRUE) #vittemolla autumn

# check tdfS
# tdf baldringe spring
head(ba.tdf[,1:10])
tail(ba.tdf[,1:10])

# tdf baldringe autumn
head(ba.tdf2[,1:10])
tail(ba.tdf2[,1:10])

# tdf christinehof
head(ch.tdf[,1:10])
tail(ch.tdf[,1:10])

# tdf hogestad
head(ho.tdf[,1:10])

# tdf vitemolla
head(vi.tdf[,1:10])
tail(vi.tdf[,1:10])

# Note: coordinates on m scale so convert to km scale
ba.tdf$X <- ba.tdf$X/1000
ba.tdf$Y <- ba.tdf$Y/1000
ba.tdf2$X <- ba.tdf2$X/1000
ba.tdf2$Y <- ba.tdf2$Y/1000
ch.tdf$X <- ch.tdf$X/1000
ch.tdf$Y <- ch.tdf$Y/1000
ho.tdf$X <- ho.tdf$X/1000
ho.tdf$Y <- ho.tdf$Y/1000
vi.tdf$X <- vi.tdf$X/1000
vi.tdf$Y <- vi.tdf$Y/1000

## take away readable columns and give right name to trap ID column
ba.tdf <- ba.tdf[,-1]
ba.tdf2 <- ba.tdf2[,-1]
ch.tdf <- ch.tdf[,-1]
ho.tdf <- ho.tdf[,-1]
vi.tdf <- vi.tdf[,-1]

colnames(ba.tdf)[1] <- "Trap_ID"
colnames(ba.tdf2)[1] <- "Trap_ID"
colnames(ch.tdf)[1] <- "Trap_ID"
colnames(ho.tdf)[1] <- "Trap_ID"
colnames(vi.tdf)[1] <- "Trap_ID"

## Add season to tdfs
ba.tdf$covs <- "A"
ba.tdf2$covs <- "A"
ch.tdf$covs <- "A"
ho.tdf$covs <- "A"
vi.tdf$covs <- "A"

ba.tdf$season <- 2
ba.tdf2$season <- 1
ch.tdf$season <- 2
ho.tdf$season <- 1
vi.tdf$season <- 2

#=====================================================================================
## Step 2: Format data, make scrFrame and state space ----

pc.data <- 
  data2oscr(edf = pc_edf,
            tdf = list(ho.tdf,ba.tdf,ba.tdf2,ch.tdf,vi.tdf), # 5 tdf files - one for each session
            sess.col = which(colnames(pc_edf)%in%"SESSION2"), 
            id.col = which(colnames(pc_edf)%in%"ANIMAL_ID2"),
            occ.col = which(colnames(pc_edf)%in%"SO"), 
            trap.col = which(colnames(pc_edf)%in%"TRAP_ID2"),
            K = rep(81,5), # no of sampling occ. per session
            ntraps = c(nrow(ho.tdf),nrow(ba.tdf),nrow(ba.tdf2),nrow(ch.tdf),nrow(vi.tdf)),
            remove.zeros = T,
            trapcov.names = "season")

# make scrFrame
pc.sf <- pc.data$scrFrame

# get summary
pc.sf

par(mfrow=c(2,2),mar=c(1,1,1,1),oma=c(0,0,0,0))
plot(pc.sf, ax = F) #plot a summary

############################# MAKE STATE SPACE ############################

# buffer should be ~3/4x sigma. 
# If you don't know sigma use HMMDM as a guide
pc.sf$mmdm 

# MMDM: 1.0km so HMMDM ~0.5km
# will use a 2km buffer here (2x MMDM)

ss.buffer <- make.ssDF(pc.sf,res=0.1,buff=2) 
## Used res=0.01 in model which is (10m x 10m) = 0.0001km2 pixels, buffer 2km
# now using res=0.1, which is 100m x 100m -> more realistic and faster running models

# plot state space
#plot(ss.buffer)

# plot state space and detections
plot(ss.buffer, pc.sf, spider=TRUE)

#=====================================================================================
## Step 3: Model fitting ----

# null model
m0 <- oSCR.fit(list(D~1,p0~1,sig~1), pc.sf, ss.buffer)
save(m0,file="polecat-m0.RData")

# session specific density
mf <- oSCR.fit(list(D~session,p0~season,sig~1), pc.sf, ss.buffer)
save(mf,file="polecat-mf.RData")

## Load the previously run models
load("polecat-m0.RData")
load("polecat-mf.RData")

## checking model results
m0
mf

#### Getting estimates for the different parameters
# detection
pred.df <- data.frame(session = factor(c(1)),season=c(1,2)) # only one session as not session-specific
pred.det <- get.real(model = mf, type = "det", newdata = pred.df)
pred.det

# sigma
pred_sig <- data.frame(session = factor(c(1))) # only one session and no season
pred.sig <- get.real(model = mf, type = "sig", newdata = pred_sig)
pred.sig

# density
# define the values we want to make prediction for
pred.df.dens <- data.frame(session = factor(c(1,2,3,4,5))) # different from above as this IS session-specific

# make predictions on the real scale
(pred.dens <- get.real(model = mf, type = "dens", newdata = pred.df.dens))

# scale up prediction from resolution of 0.0001 km2 / 100 m2
# to resolution of 1000 ha (10 km2)
(pred.dens <- get.real(model = mf, type = "dens", newdata = pred.df.dens, d.factor = 100000))

#=====================================================================================
## Step 6: Data visualisation ----

library(ggplot2)

# plot session estimates, by sex, from session-specific model (m1)

# make a table with our session-specific and constant density estimates
# estimates and confidence intervals bind into table
est.tab <- rbind(pred.dens[1:5,c(1,3,4)])

# see table
est.tab

# add season column
est.tab$season <- factor(c("Spring","Autumn","Spring","Autumn","Autumn"))

# add location column
est.tab$location <- factor(c("Högestad","Baldringe","Baldringe","Christinehof","Vitemölla"))

# see table with columns added
est.tab

# plot density per session and by sex
# add a line to show where the lion management plan started

gd <- ggplot(est.tab,aes(x=location,y=estimate,fill=season)) +
  geom_errorbar(aes(ymin=lwr,ymax=upr,width=0)) +
  geom_point(aes(color=season,shape=season), size=3.5) +
  scale_shape_manual(values=c(15,19)) + 
  scale_color_manual(values=c("black", "black")) +
  facet_grid(season~.) +
  theme_bw() +
  ylab(expression(bold(Density~~"(per 10 "*km^"2"*")"))) +
  xlab("Location") +
  theme(axis.title.y = element_text(vjust = 2.5)) +
  scale_y_continuous(limits=c(0,16),breaks=c(0,2,4,6,8,10,12,14,16)) +
  theme(axis.text=element_text(colour="black")) +
  theme(legend.position = "none") +    
  theme(axis.title = element_text(face="bold"))
gd

# make a density map
#pred <- predict.oSCR(scrFrame = pc.sf, mf, ssDF = ss.buffer)
save(pred,file="mf-map-prediction.RData")
load("mf-map-prediction.RData")

pdf("figureS2.pdf",width=9,height=6)
par(mfrow=c(2,3),mar=c(2,2,2,2))
image(pred$r[[1]],main="Högestad Spring")
points(ho.tdf[,2:3], pch=20, lwd=0.5)
image(pred$r[[2]],main="Baldringe Autumn")
points(ba.tdf[,2:3], pch=20, lwd=0.5)
image(pred$r[[3]],main="Baldringe Spring")
points(ba.tdf2[,2:3], pch=20, lwd=0.5)
image(pred$r[[4]],main="Christinehof Autumn")
points(ch.tdf[,2:3], pch=20, lwd=0.5)
image(pred$r[[5]],main="Vitemölla Autumn")
points(vi.tdf[,2:3], pch=20, lwd=0.5)
dev.off()

png("figureS2.png",width=900,height=600)
par(mfrow=c(2,3),mar=c(2,2,2,2))
image(pred$r[[1]],main="Högestad Spring")
points(ho.tdf[,2:3], pch=20, lwd=0.5)
image(pred$r[[2]],main="Baldringe Autumn")
points(ba.tdf[,2:3], pch=20, lwd=0.5)
image(pred$r[[3]],main="Baldringe Spring")
points(ba.tdf2[,2:3], pch=20, lwd=0.5)
image(pred$r[[4]],main="Christinehof Autumn")
points(ch.tdf[,2:3], pch=20, lwd=0.5)
image(pred$r[[5]],main="Vitemölla Autumn")
points(vi.tdf[,2:3], pch=20, lwd=0.5)
dev.off()
