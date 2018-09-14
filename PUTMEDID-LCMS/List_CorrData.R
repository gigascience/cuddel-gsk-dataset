# Title     : TODO
# Objective : TODO
# Created by: peterli
# Created on: 10/9/2018

library('pracma')
library('Hmisc')

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

listofdata <- data.frame(PeakNo=integer(), PeakNo=integer(), tcorr=integer())

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

# Find peak number order
order <- shellSort(PeakNoDble)

# Sort out mydata object
dimensionality <- nrow(mydata)
ncols <- ncol(mydata)
# Check if mydata has correct data structure
if (calccorr == "P" | calccorr == "p" | calccorr == "S" | calccorr == "s") {
    # Data needs to be n rows (peaks) x ncols columns (samples)
    if (ncols+1 == dimensionality | dimensionality+1 != npeaks) {
        OKdata <- FALSE
        tempStr <- "Incorrect data to calculate correlation - program aborted"
    }
} else {
    if (ncols != dimensionality | dimensionality+1 != npeaks) {
        OKdata <- FALSE
        tempStr <- "Data should be a square matrix (n x n) with n = no. of peaks - program aborted"
    }
}

if (OKdata == TRUE) {
    # Create matrix from mydata object with ncols rows and dimensionality columns
    # double[][] values=new double [ncols][dimensionality];
    # for(int i=0;i<dimensionality;i++) {
    #     line=mydata.get(i).trim();
    #     items=line.split(myseparator);
    #     for(int j=0;j<ncols;j++) {
    #         values[j][i] = Double.parseDouble(items[j].trim());
    #     }
    # }
    values <- t(mydata)

    rtgap <- 0
    for(i in 1:npeaks) {
        j <- i + 1
        for(j in 1:npeaks) {
            # uses Peak No. order
            rtgap=rt[order[i]]-rt[order[j]]
            if (rtgap<0) {
                rtgap <- -rtgap
            }
            if (rtgap<rtdiff) {
                if (calccorr == "P" | calccorr == "p") {
                    tcorr <- getPearsonCorrelationMatrix(order[i]-1, order[j]-1, values)
                    cor.test()
                }
                else if (calccorr == "S" | calccorr == "s") {
                    tcorr <- getSpearmanCorrelationMatrix(order[i]-1, order[j]-1, values)
                }
                else {
                    tcorr <- values[order[i]-1][order[j]-1];
                }
            }
            if (tcorr>corrlim) {
                if (PeakNo[order[i]]<PeakNo[order[j]]) {
                    listofdata.append(PeakNo[order[i]]+"\t"+PeakNo[order[j]]+"\t"+tcorr+"\n");
                }
                else {
                    listofdata.append(PeakNo[order[j]]+"\t"+PeakNo[order[i]]+"\t"+tcorr+"\n");
                }
            }
        }
    }
} else {
        listofdata <- rbind(listofdata, tempStr)
}

corrlist <- listofdata
