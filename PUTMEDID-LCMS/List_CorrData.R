# Title     : TODO
# Objective : TODO
# Created by: peterli
# Created on: 10/9/2018

library('pracma')

# Read input data from files
mydata <- read.table("/Users/peterli/PUTMEDID_LCMS_v1.01/Study_posdata.txt", sep='\t', header=FALSE)
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
# iline <- colnames(indata)
# String iline=indata.get(0).trim();
# String[] iitems=iline.split("\t");
# // Number of peaks in npeaks
# int npeaks=indata.size();
npeaks <- nrow(indata)
print(paste0("npeaks = ", npeaks))
# // To hold retention times
# double[] rt = new double[npeaks];
rt <- vector(mode="double", length=npeaks)
rt <- indata[, 3]
# // To hold peak numbers
# int [] PeakNo = new int[npeaks];
PeakNo <- vector(mode="integer", length=npeaks)
PeakNo <- indata[, 1]
# // To hold peak numbers as a double object
# double [] PeakNoDble = new double[npeaks]; //DCW
PeakNoDble <- vector(mode="double", length=npeaks)
PeakNoDble <- as.double(PeakNo)
# for(int i=1;i<npeaks;i++) {
    # iline=indata.get(i).trim();
    # iitems=iline.split("\t");
    # PeakNo[i]=Integer.parseInt(iitems[0].trim());
    # PeakNoDble[i] = (double)(PeakNo[i]);
    # rt[i] = Double.parseDouble(iitems[2].trim());
# }

# Find peak number order
order <- shellSort(PeakNoDble)

# Sort out mydata object
dimensionality <- nrow(mydata)
print(paste0("dimensionality = ", dimensionality))
ncols <- ncol(mydata)
print(paste0("ncols = ", ncols))
# Check if mydata has correct data structure
if (calccorr == "P" | calccorr == "p" | calccorr == "S" | calccorr == "s") {
    # Data needs to be n rows (peaks) x ncols columns (samples)
    # In beanshell mydata has header so has extra row, it is ((ncols+1 == dimensionality) || (dimensionality+1 != npeaks))
    if (ncols+1 == dimensionality-1 | dimensionality != npeaks) {
        OKdata <- FALSE
        tempStr <- "Incorrect data to calculate correlation - program aborted"
    }
} else {
    # In beanshell mydata has header so has extra row, it is (ncols != dimensionality) || (dimensionality+1 != npeaks)
    if (ncols != dimensionality-1 | dimensionality != npeaks) {
        OKdata <- FALSE
        tempStr <- "Data should be a square matrix (n x n) with n = no. of peaks - program aborted"
    }
}

# To hold correlation results between peaks
listofdata <- mat.or.vec(0, 3)

if (OKdata == TRUE) {
    values <- t(mydata)

    rtgap <- 0
    for(i in 1:npeaks) {
        # if(i == 3)
        #     break
        j <- i + 1
        #print(paste0("i: ", i))
        for(j in i:npeaks) {
            #print(paste0("j: ", j))
            tcorr <- 0
            # uses Peak No. order
            rtgap=rt[order[i]]-rt[order[j]]
            #print(paste0("order[i]: ", order[i]))

            #print(paste0("order[j]: ", order[j]))
            if (rtgap<0) {
                #print(paste0("rtgap is less than zero"))
                rtgap <- -rtgap
            }
            if (rtgap<rtdiff) {
                #print(paste0("rtgap is less than rtdiff"))
                if (calccorr == "P" | calccorr == "p") {
                    # tcorr <- getPearsonCorrelationMatrix(order[i]-1, order[j]-1, values)
                    #print(paste0("Doing Pearson correlation test"))
                    tcorr <- cor.test(values[, i], values[, j], method = "pearson")
                    tcorr <- tcorr$estimate
                    #print(paste0("tcorr result: ", tcorr))
                }
                else if (calccorr == "S" | calccorr == "s") {
                    # tcorr <- getSpearmanCorrelationMatrix(order[i]-1, order[j]-1, values)
                    tcorr <- cor.test(values[, order[i]-1], values[, order[j]-1], method = "spearman")
                    tcorr <- tcorr$estimate
                }
                else {
                    tcorr <- values[order[i]-1][order[j]-1]
                }
            }

            if (tcorr>corrlim) {
                print(paste0("tcorr is greater than corrlim threshold: ", corrlim))
                # print(paste0("i: ", i))
                # print(paste0(values[, i]))
                # print(paste0("j: ", j))
                # print(paste0(values[, j]))
                print(paste0("tcorr result: ", tcorr))
                if (PeakNo[order[i]]<PeakNo[order[j]]) {
                    #listofdata.append(PeakNo[order[i]]+"\t"+PeakNo[order[j]]+"\t"+tcorr+"\n");
                    row <- c(order[i], order[j], tcorr)
                    listofdata <- rbind(listofdata, row)
                }
                else {
                    #listofdata.append(PeakNo[order[j]]+"\t"+PeakNo[order[i]]+"\t"+tcorr+"\n");
                    row <- c(order[j], order[i], tcorr)
                    rbind(listofdata, row)
                }
            }
        }
    }
} else {
    listofdata <- rbind(listofdata, tempStr)
}

corrlist <- listofdata
#print(corrlist)

write.table(corrlist,
            file = "Study_pos_CorrListP_Routput.txt",
            append = TRUE,
            sep = "\t",
            row.names=FALSE,
            col.names=FALSE,
            quote=FALSE)
