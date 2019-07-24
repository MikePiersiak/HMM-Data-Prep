# HMM-Data-Prep
preparation of data for HMM imputation using momentuHMM

require(momentuHMM)

final.crabs_pt1 <- read.csv("FourHrTags_Final_1.csv")


str(final.crabs_pt1)

final.crabs_pt1$release.time <- as.POSIXct(final.crabs_pt1$release.time, format = "%m/%d/%Y %H:%M", tz = "UTC")
final.crabs_pt1$time <- as.POSIXct(final.crabs_pt1$time, format = "%m/%d/%Y %H:%M", tz = "UTC")



##specify columns
final.crabs_pt1 <- final.crabs_pt1[,c(4,2,7,8)]
head(final.crabs_pt1)

##regularize tracks so total nDetections can be determined and <100 tags can be removed

crwOut1 <- crawlWrap(obsData = final.crabs_pt1, timeStep = "30 mins", theta = c(2, -0.001), fixPar = c(NA,NA))
crabs1 <- prepData(crwOut1) 


str(crabs1)


all.tags <- NULL

for (i in levels(crabs1$ID)){
  
  Tag <- subset(crabs1, ID==i)
  if (nrow(Tag) < 100) {
    Tag$TrackStatus <- "Bad"
  } else {
    Tag$TrackStatus <- "Good"
  }
  
  all.tags <- rbind(all.tags, Tag)
  
  
}

head(all.tags)
all.tags$TrackStatus <- as.factor(all.tags$TrackStatus)
str(all.tags)

bad.tags <- subset(all.tags, TrackStatus=="Bad")

##Identifies bad tags for removal
bad.tags$ID <- droplevels(bad.tags$ID)
a <- as.factor(levels(bad.tags$ID))

final.crabs_pt1 <- read.csv("FourHrTags_Final_1.csv")
final.crabs_pt1 <- final.crabs_pt1[,c(4,2,7,8,20,18,13,15,16)]
final.crabs_pt1$ID <- as.factor(final.crabs_pt1$ID)

##Removes bad tags
str(final.crabs_pt1)
final.crabs_pt1 <- subset(final.crabs_pt1, ! ID %in% a) 
final.crabs_pt1$ID <- droplevels(final.crabs_pt1$ID)
str(final.crabs_pt1)

##Interpolate covariate values for crawlMerge
##requires prediction of covariate values at same time step used in crawlWrap
##momentuHMM specifically requires columns named 'ID' 'time' 'x' and 'y' for crawlWrap

colnames(final.crabs_pt1) <-  c("ID","time","x","y","Bottom.Temp","Exposure","Location","Year","Release.Time")
final.crabs_pt1$Release.Time <- as.POSIXct(final.crabs_pt1$Release.Time, format = "%m/%d/%Y %H:%M", tz = "UTC")
final.crabs_pt1$time <- as.POSIXct(final.crabs_pt1$time, format = "%m/%d/%Y %H:%M", tz = "UTC")

##This loop is for interpolating environmental data for the crawl merge function

pt1.sheet <- NULL
for (i in levels(final.crabs_pt1$ID)){
  
  Tag <- subset(final.crabs_pt1, ID==i)
  
#Create a list of times at the desired time interval. Here 30-mins is used.
time.test <- seq(Tag$time[1], Tag$time[length(Tag$time)], 
                 by = 30 * 60)

##approx lat/long locations at time intervals stipulated above.
location <- as.data.frame(cbind(lon = approx(Tag$time, Tag$x, xout = time.test)$y, 
                                lat = approx(Tag$time, Tag$y, xout = time.test)$y))

##approximate temperature values
temp <- as.data.frame(cbind(temp = approx(Tag$time, Tag$Bottom.Temp, xout = time.test)$y))

##bind new interpolated track data
all.dat <- cbind(time = time.test, location, temp)

##Add year and time since release (TSR)
all.dat$year <- substring(all.dat$time,1,4)
all.dat$TSR <- as.numeric(all.dat$time - Tag$Release.Time[1]) 

##rename columns. MomentuHMM requires column names of 'time', 'x', and 'y' for time and lat/long respectively
colnames(all.dat) <- c("time", 
                       "x", 
                       "y", 
                       "Bottom.Temp",
                       "Year",
                       "TSR")

##add exposure treatments to new data frame
a <- max(Tag$time[which(Tag$Exposure == "Pre")])
b <- min(Tag$time[which(Tag$Exposure == "Post")])

all.dat$exp1 <- ifelse((all.dat$time<=a), "Pre","Other")
all.dat$exp2 <- ifelse((all.dat$time>a) & 
                         (all.dat$time<b), "During",all.dat$exp1)
all.dat$Exposure <- ifelse((all.dat$time>=b), "Post",all.dat$exp2)

##add location and ID to data frame
all.dat$location <- unique(Tag$location)
all.dat$ID <- unique(Tag$ID)
head(all.dat)
final.dat <- all.dat[,c(10,1,4,5,6,9)]
pt1.sheet <- rbind(pt1.sheet, final.dat)

}

str(pt1.sheet)
##add column for Hour of the day (HoD) for use as model covariate in momentuHMM
pt1.sheet$HoD <- as.integer(strftime(pt1.sheet$time, format = "%H", tz="UTC"))

pt1.sheet$Location <- Tag$Location[1]
head(pt1.sheet)

final.crabs_pt1 <- read.csv("FourHrTags_Final_1.csv")
final.crabs_pt1 <- final.crabs_pt1[,c(4,2,7,8)]
final.crabs_pt1$ID <- as.factor(final.crabs_pt1$ID)
a <- as.factor(levels(bad.tags$ID))
final.crabs_pt1 <- subset(final.crabs_pt1, ! ID %in% a) 
final.crabs_pt1$ID <- droplevels(final.crabs_pt1$ID)
final.crabs_pt1$time <- as.POSIXct(final.crabs_pt1$time, format = "%m/%d/%Y %H:%M", tz = "UTC")
str(final.crabs_pt1)

crwOut1 <- crawlWrap(obsData = final.crabs_pt1, timeStep = "30 mins", theta = c(2, -0.001), fixPar = c(NA,NA))
crawl.pt1 <- crawlMerge(crwOut1, pt1.sheet, "time")
pt1Data <- prepData(crawl.pt1, covNames = c("Bottom.Temp",
                                            "Exposure",
                                            "Location",
                                            "Year",
                                            "TSR",
                                            "HoD"))
View(pt1Data)


###Now repeat the above process for the second half of the data set and then merge the outputs from prepData
###This will amalgamate the data for imputation in momentuHMM

final.crabs_pt2 <- read.csv("FourHrTags_Final_2.csv")

str(final.crabs_pt2)

final.crabs_pt2$release.time <- as.POSIXct(final.crabs_pt2$release.time, format = "%m/%d/%Y %H:%M", tz = "UTC")
final.crabs_pt2$time <- as.POSIXct(final.crabs_pt2$time, format = "%m/%d/%Y %H:%M", tz = "UTC")

##specify columns
final.crabs_pt2 <- final.crabs_pt2[,c(4,2,7,8)]
head(final.crabs_pt2)

##regularize tracks so total nDetections can be determined and <100 tags can be removed
crwOut2 <- crawlWrap(obsData = final.crabs_pt2, timeStep = "30 mins", theta = c(2, -0.001), fixPar = c(NA,NA))
crabs2 <- prepData(crwOut2) 

str(crabs2)

all.tags <- NULL

for (i in levels(crabs2$ID)){
  
  Tag <- subset(crabs2, ID==i)
  if (nrow(Tag) < 100) {
    Tag$TrackStatus <- "Bad"
  } else {
    Tag$TrackStatus <- "Good"
  }
  
  all.tags <- rbind(all.tags, Tag)
  
  
}

head(all.tags)
all.tags$TrackStatus <- as.factor(all.tags$TrackStatus)
str(all.tags)

bad.tags <- subset(all.tags, TrackStatus=="Bad")

##Identifies bad tags for removal
bad.tags$ID <- droplevels(bad.tags$ID)
a <- as.factor(levels(bad.tags$ID))

final.crabs_pt2 <- read.csv("FourHrTags_Final_2.csv")
final.crabs_pt2 <- final.crabs_pt2[,c(4,2,7,8,20,18,13,15,16)]
final.crabs_pt2$ID <- as.factor(final.crabs_pt2$ID)

##Removes bad tags
str(final.crabs_pt2)
final.crabs_pt2 <- subset(final.crabs_pt2, ! ID %in% a) 
final.crabs_pt2$ID <- droplevels(final.crabs_pt2$ID)
str(final.crabs_pt2)

##Interpolate covariate values for crawlMerge
##requires prediction of covariate values at same time step used in crawlWrap
##momentuHMM specifically requires columns named 'ID' 'time' 'x' and 'y' for crawlWrap

colnames(final.crabs_pt2) <-  c("ID","time","x","y","Bottom.Temp","Exposure","Location","Year","Release.Time")
final.crabs_pt2$Release.Time <- as.POSIXct(final.crabs_pt2$Release.Time, format = "%m/%d/%Y %H:%M", tz = "UTC")
final.crabs_pt2$time <- as.POSIXct(final.crabs_pt2$time, format = "%m/%d/%Y %H:%M", tz = "UTC")


##This loop is for interpolating environmental data for the crawl merge function
pt2.sheet <- NULL
for (i in levels(final.crabs_pt2$ID)){
  
  Tag <- subset(final.crabs_pt2, ID==i)
  
  #Create a list of times at the desired time interval. Here 30-mins is used.
  time.test <- seq(Tag$time[1], Tag$time[length(Tag$time)], 
                   by = 30 * 60)
  
  ##approx lat/long locations at time intervals stipulated above.
  location <- as.data.frame(cbind(lon = approx(Tag$time, Tag$x, xout = time.test)$y, 
                                  lat = approx(Tag$time, Tag$y, xout = time.test)$y))
  
  ##approximate temperature values
  temp <- as.data.frame(cbind(temp = approx(Tag$time, Tag$Bottom.Temp, xout = time.test)$y))
  
  ##bind new interpolated track data
  all.dat <- cbind(time = time.test, location, temp)
  
  ##Add year and time since release (TSR)
  all.dat$year <- substring(all.dat$time,1,4)
  all.dat$TSR <- as.numeric(all.dat$time - Tag$Release.Time[1]) 
  
  ##rename columns. MomentuHMM requires column names of 'time', 'x', and 'y' for time and lat/long respectively
  colnames(all.dat) <- c("time", 
                         "x", 
                         "y", 
                         "Bottom.Temp",
                         "Year",
                         "TSR")
  
  ##add exposure treatments to new data frame
  a <- max(Tag$time[which(Tag$Exposure == "Pre")])
  b <- min(Tag$time[which(Tag$Exposure == "Post")])
  
  all.dat$exp1 <- ifelse((all.dat$time<=a), "Pre","Other")
  all.dat$exp2 <- ifelse((all.dat$time>a) & 
                           (all.dat$time<b), "During",all.dat$exp1)
  all.dat$Exposure <- ifelse((all.dat$time>=b), "Post",all.dat$exp2)
  
  ##add location and ID to data frame
  all.dat$location <- unique(Tag$location)
  all.dat$ID <- unique(Tag$ID)
  head(all.dat)
  final.dat <- all.dat[,c(10,1,4,5,6,9)]
  pt2.sheet <- rbind(pt2.sheet, final.dat)
  
}

str(pt2.sheet)
##add column for Hour of the day (HoD) for use as model covariate in momentuHMM
pt2.sheet$HoD <- as.integer(strftime(pt2.sheet$time, format = "%H", tz="UTC"))

pt2.sheet$Location <- Tag$Location[1]
head(pt2.sheet)

final.crabs_pt2 <- read.csv("FourHrTags_Final_2.csv")
final.crabs_pt2 <- final.crabs_pt2[,c(4,2,7,8)]
final.crabs_pt2$ID <- as.factor(final.crabs_pt2$ID)
a <- as.factor(levels(bad.tags$ID))
final.crabs_pt2 <- subset(final.crabs_pt2, ! ID %in% a) 
final.crabs_pt2$ID <- droplevels(final.crabs_pt2$ID)
final.crabs_pt2$time <- as.POSIXct(final.crabs_pt2$time, format = "%m/%d/%Y %H:%M", tz = "UTC")
str(final.crabs_pt2)

crwOut2 <- crawlWrap(obsData = final.crabs_pt2, timeStep = "30 mins", theta = c(2, -0.001), fixPar = c(NA,NA))
plot(crwOut2)
crawl.pt2 <- crawlMerge(crwOut2, pt2.sheet, "time")
pt2Data <- prepData(crawl.pt2, covNames = c("Bottom.Temp",
                                            "Exposure",
                                            "Location",
                                            "Year",
                                            "TSR",
                                            "HoD"))
View(pt2Data)

all.crabs <- rbind(pt1Data, pt2Data)
str(all.crabs)

###data frame 'all.crabs' is ready for momentuHMM
