# Title     : annotate_Massmatch.R
# Objective : Annotating masses to features
# Created by: peterli
# Created on: 24/9/2018

# Load input files
# Select Isotopes FileTav_Isotopes.txt
indata <- read.table("/Users/peterli/PUTMEDID_LCMS_v1.01/Tav_Isotopes.txt", sep='\t', header=TRUE)
# Select Adducts File Tav_Adducts.txt
adddata <- read.table("/Users/peterli/PUTMEDID_LCMS_v1.01/Tav_Adducts.txt", sep='\t', header=TRUE)
# Select Reference Mass File trimMMD_sortAmass.txt
refdata <- read.table("/Users/peterli/PUTMEDID_LCMS_v1.01/trimMMD_sortAmass.txt", sep='\t', header=TRUE)
# Select Mass File to Check Study_pospeaks.txt
mydata <- read.table("/Users/peterli/PUTMEDID_LCMS_v1.01/Study_pospeaks.txt", sep='\t', header=TRUE)
# Select Correlations File Study_pos_CorrListP.txt
corrdata <- read.table("/Users/peterli/PUTMEDID_LCMS_v1.01/Study_pos_CorrListP.txt", sep='\t', header=FALSE)

# Define parameters
# Ion_mode
elabel <- "pos"
# Mass tolerance
lowval <- "3"
# min_massdiff
limtol <- "0.00075"
# rtlim
rtlim <- "2.5"
# rtmin
rtmin <- "30"
# rtmax
rtmax <- "1200"
# Corrlim
mycorrlim <- "0.9"

# Process inputs
# PeakName <- new int[mydata.size()];
PeakName <- nrow(mydata)

# int[] sortPeaks=new int[mydata.size()];
sortPeaks <- nrow(mydata)

# String [] WorkLabel= new String[mydata.size()];
WorkLabel <- vector(mode="char", length=nrow(mydata))

# double[] massdata=new double[mydata.size()];
massdata <- vector(mode="double", length=nrow(mydata))

# double[] matchMF=new double[mydata.size()];
matchMF <- vector(mode="double", length=nrow(mydata))

# double[] RTdata=new double[mydata.size()];
RTdata <- vector(mode="double", length=nrow(mydata))

# double[] medval=new double[mydata.size()];
medval <- vector(mode="double", length=nrow(mydata))

# int[] useLabel=new int[mydata.size()];
useLabel <- vector(mode="integer", length=nrow(mydata))

# Integer headLen;
headLen <- 0

# Integer rheadLen;
rheadLen <- 0

# double rtlim=new Double(rtlim).doubleValue();
rtlim <- as.double(rtlim)

# double rtmin=new Double(rtmin).doubleValue();
rtmin <- as.double(rtmin)

# double rtmax=new Double(rtmax).doubleValue();
rtmax <- as.double(rtmax)

# double tol=new Double(lowval).doubleValue(); // use tol*2 for adduct/mass diffs
tol <- as.double(lowval)

# double masslim=new Double(limtol).doubleValue();
masslim <- as.double(limtol)

# double corrlim=new Double(mycorrlim).doubleValue();
corrlim <- as.double(mycorrlim)

# double[] MFdata=new double[refdata.size()];
MFdata <- vector(mode="double", length=nrow(refdata))

# String [] MFformula= new String[refdata.size()];
MFformula <- vector(mode="char", length=nrow(refdata))

# String ppmStr;
ppmStr <- ""

# int nsym;
nsym <- 0

# double kcount=0;
kcount <- 0

# int tcount=0;
tcount <- 0

# boolean sortedMF=true;
sortedMF <- TRUE

# boolean rsortedRefMF=true;
rsortedRefMF <- TRUE

# String tempStr1="";
tempStr1 <- ""

# String tempStr3="";
tempStr3 <- ""

# UseLabel is 0 for A, 1 for I, F, .... 9 for E (invalid RT)
# double etol=1; // for positive ion mode
etol <- as.double(1)

if (substr(elabel, 0, 1) == "n" | substr(elabel, 0, 1) == "N") {
    etol <- -1
}

# Outputs
# StringBuffer listofdata = new StringBuffer(); // annotated data
listofdata <- character(0)

# StringBuffer listofcorrs = new StringBuffer();
listofcorrs <- character(0)

# StringBuffer listofsortcorrs = new StringBuffer();
listofcorrs <- character(0)

# Integer gno=0;
gno <- 0

# Integer regroup=0;
regroup <- 0

# int[] PeakMatch=new int[mydata.size()];
PeakMatch <- numeric(nrow(mydata))

# int[] PeakIndex=new int[mydata.size()];
PeakIndex <- numeric(nrow(mydata))

# int[] Gcount=new int[mydata.size()];
Gcount <- numeric(nrow(mydata))

# int[] SGcount=new int[mydata.size()];
SGcount <- numeric(nrow(mydata))

# int[] Gorig=new int[mydata.size()];
Gorig <- numeric(nrow(mydata))

# double[] PCharge=new double[mydata.size()];
PCharge <- double(nrow(mydata))

# double[] Inum=new double[mydata.size()];
Inum <- double(nrow(mydata))

# double[] MassMatch=new double[mydata.size()];
MassMatch <- double(nrow(mydata))

# String[] PeakLabel=new String[mydata.size()];
PeakLabel <- char(nrow(mydata))

# String[] SaltLabel=new String[mydata.size()];
SaltLabel <- char(nrow(mydata))

# String[] ILabel=new String[mydata.size()];
ILabel <- char(nrow(mydata))

# int[] Pidx1=new int[corrdata.size()];
Pidx1 <- numeric(nrow(corrdata))

# int[] Pidx2=new int[corrdata.size()];
Pidx2 <- numeric(nrow(corrdata))

# int[] totCorrs=new int[corrdata.size()];
totCorrs <- numeric(nrow(corrdata))

# double[] usecorr=new double[corrdata.size()];
usecorr <- double(nrow(corrdata))

# Integer idx1;
idx1 <- 0

# Integer idx2;
idx2 <- 0

# Initialise - read sorted mass and RT data into an Array for fast search
# double compMF=0;
compMF <- double(0)

for (i in i:mydata) {
    inStr <- mydata[i]
    tempdata <- strsplit(inStr, "\t")

    print(tempdata[0])

    # int pval = Integer.parseInt(tempdata[0]);
    pval <- as.integer(tempdata[0])

    # double massval = new Double(tempdata[1]).doubleValue();
    massval <- as.double(tempdata[1])
    # double tval = new Double(tempdata[2]).doubleValue();
    tval <- as.double(tempdata[2])

    # double respval = new Double(tempdata[3]).doubleValue();
    respval <- as.double(tempdata[3])

    # PeakName[i]=pval;
    PeakName[i] <- pval

    # sortPeaks[i]=PeakName[i];
    sortPeaks[i] <- PeakName[i]

    # massdata[i]=massval;
    massdata[i] <- massval

    # MassMatch[i]=massval;
    MassMatch[i] <- massval

    # matchMF[i]=0;
    matchMF[i] <- 0

    if (massval<compMF) {
        sortedMF <- FALSE
        tempStr1 <- "Check file is not sorted by Accurate Mass"
    }

    compMF <- massval
    RTdata[i] <- tval
    medval[i] <- respval
    PeakMatch[i] <- pval # default - initialise Peak to itself before matching
    PeakIndex[i] <- i
    Gcount[i] <- 0
    Gorig[i] <- 0
    SGcount[i] <- 0
    PCharge[i] <- 1
    Inum[i] <- 0
    PeakLabel[i] <- "A" # default - initialise all Peaks to A - any/all
    WorkLabel[i] <- "A"
    SaltLabel[i] <- ""
    ILabel[i] <- "N"
    useLabel[i] <- 0
    if (RTdata[i]<rtmin | RTdata[i]>rtmax){
        useLabel[i] <- 9
        PeakLabel[i] <- "E"
    }
}