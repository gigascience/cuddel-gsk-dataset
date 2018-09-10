# Title     : TODO
# Objective : TODO
# Created by: peterli
# Created on: 10/9/2018

# Create input data
mydata <- read.table("/Users/peterli/PUTMEDID_LCMS_v1.01/Study_posdata.txt", sep='\t', header=TRUE)
indata <- read.table("/Users/peterli/PUTMEDID_LCMS_v1.01/Study_pospeaks.txt", sep='\t', header=TRUE)
calccorr <- "P"
corrcheck <- "0.7"
rtread <- "5"
filetype <- "txt"

rtdiff <- as.double(rtread)
corrlim <- as.double(corrcheck)
tcorr <- 0.0
tempStr <- ""
OKdata <- TRUE

# Sort out indata object
iline <- colnames(indata)
# String iline=indata.get(0).trim();
# String[] iitems=iline.split("\t");
# // Number of peaks in npeaks
# int npeaks=indata.size();
npeaks <- nrow(indata)
# // To hold retention times
# double[] rt = new double[npeaks];
rt <- vector(mode="double", length=npeaks)
# // To hold peak numbers
# int [] PeakNo = new int[npeaks];
PeakNo <- vector(mode="integer", length=npeaks)
# // To hold peak numbers as a double object
# double [] PeakNoDble = new double[npeaks]; //DCW
PeakNoDble <- vector(mode="double", length=npeaks)
# for(int i=1;i<npeaks;i++) {
    # iline=indata.get(i).trim();
    # iitems=iline.split("\t");
    # PeakNo[i]=Integer.parseInt(iitems[0].trim());
    # PeakNoDble[i] = (double)(PeakNo[i]);
    # rt[i] = Double.parseDouble(iitems[2].trim());
# }
PeakNo <- indata[, 1]
PeakNoDble <- as.double(PeakNo)
rt <- indata[, 3]

# Sort out mydata object
dimensionality <- nrow(mydata)
ncols_mydata <- ncol(mydata)
# Check if mydata has correct data structure
if (calccorr.substring(0,1).equals("P") || calccorr.substring(0,1).equals("p") || calccorr.substring(0,1).equals("S") || calccorr.substring(0,1).equals("s")) {
    // data needs to be n rows (peaks) x ncols columns (samples)
    if ((ncols+1 == dimensionality) || (dimensionality+1 != npeaks)) {
        OKdata=false;
        tempStr="Incorrect data to calculate correlation - program aborted";
    }
}
else {
    if ((ncols != dimensionality) || (dimensionality+1 != npeaks)) {
        OKdata=false;
        tempStr="Data should be a square matrix (n x n) with n = no. of peaks - program aborted";
    }
}
