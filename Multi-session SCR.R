###################### Multi-session SCR on polecat data 2021 #########################

## Load package
library(oSCR)      # Handling SCR data and models

## Thanks to Rob Davis for providing the script on which this script was based.
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


############################# MAKE STATE SPACE ############################

# buffer should be ~3/4x sigma. 
# If you don't know sigma use HMMDM as a guide
pc.sf$mmdm 

# MMDM: 1.0km so HMMDM ~0.5km
# will use a 2km buffer here (2x MMDM)

ss.buffer <- make.ssDF(pc.sf,res=0.1,buff=2) 
# Using res=0.1 in model which is (100m x 100m) = 0.01km2 pixels, buffer 2km

#=====================================================================================
## Step 3: Model fitting ----

# session specific density with seasonal estimate for p0
mf <- oSCR.fit(list(D~session,p0~season,sig~1), pc.sf, ss.buffer)

## checking model results
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

# make predictions at a resolution of 1000 ha (10 km2)
(pred.dens <- get.real(model = mf, type = "dens", newdata = pred.df.dens, d.factor = 100000))
