# Title     : annotate_Massmatch.R
# Objective : Annotating masses to features
# Created by: peterli
# Created on: 24/9/2018

library('pracma')  # For shell sort but can use sort(PeakName, method = "shell")
library('Rfast')   # For match function instead of binary search
library('useful')  # For shift.column

annotateMassmatch <- function() {

    # Load input files
    indata <- read.table("/Users/peterli/PUTMEDID_LCMS_v1.01/Tav_Isotopes.txt", sep = '\t', header = TRUE)
    adddata <- read.table("/Users/peterli/PUTMEDID_LCMS_v1.01/Tav_Adducts.txt", sep = '\t', header = TRUE)
    refdata <- read.table("/Users/peterli/PUTMEDID_LCMS_v1.01/trimMMD_sortAmass.txt", sep = '\t', header = TRUE)
    mydata <- read.table("/Users/peterli/PUTMEDID_LCMS_v1.01/Study_pospeaks.txt", sep = '\t', header = TRUE)
    corrdata <- read.table("/Users/peterli/PUTMEDID_LCMS_v1.01/Study_pos_CorrListP.txt", sep = '\t', header = FALSE)

    # Define parameters
    elabel <- "pos"      # Ion_mode
    lowval <- "3"        # Mass tolerance
    limtol <- "0.00075"  # min_massdiff
    rtlim <- "2.5"       # rtlim
    rtmin <- "30"        # rtmin
    rtmax <- "1200"      # rtmax
    mycorrlim <- "0.9"   # Corrlim

    # Process inputs
    PeakName <- nrow(mydata)
    sortPeaks <- nrow(mydata)
    WorkLabel <- character(nrow(mydata))

    massdata <- vector(mode = "double", length = nrow(mydata))
    matchMF <- vector(mode = "double", length = nrow(mydata))
    RTdata <- vector(mode = "double", length = nrow(mydata))
    medval <- vector(mode = "double", length = nrow(mydata))
    useLabel <- vector(mode = "integer", length = nrow(mydata))

    rtlim <- as.double(rtlim)
    rtmin <- as.double(rtmin)
    rtmax <- as.double(rtmax)
    tol <- as.double(lowval)
    masslim <- as.double(limtol)
    corrlim <- as.double(mycorrlim)

    MFdata <- vector(mode = "double", length = nrow(refdata))

    MFformula <- character(nrow(refdata))

    ppmStr <- ""

    headLen <- 0
    rheadLen <- 0
    nsym <- 0
    kcount <- 0
    tcount <- 0
    gno <- 0
    regroup <- 0

    sortedMF <- TRUE
    rsortedRefMF <- TRUE

    tempStr1 <- ""
    tempStr3 <- ""

    # UseLabel is 0 for A, 1 for I, F, .... 9 for E (invalid RT)
    # double etol=1 // for positive ion mode
    etol <- as.double(1)

    # Set etol value based on elabel positive or negative ion mode
    if (substr(elabel, 0, 1) == "n" | substr(elabel, 0, 1) == "N") {
        etol <- -1
    }

    # Outputs
    check <- matrix(ncol=7)
    colnames(check) <- c("tdiff", "rtlim", "useLabel[oidx1]", "useLabel[oidx2]", "ncorr", "oidx1", "oidx2")

    # StringBuffer listofdata = new StringBuffer() // annotated data
    listofdata <- character(0)
    # StringBuffer listofcorrs = new StringBuffer()
    listofcorrs <- character(0)
    # StringBuffer listofsortcorrs = new StringBuffer()
    listofsortcorrs <- character(0)

    PeakMatch <- numeric(nrow(mydata))
    PeakIndex <- numeric(nrow(mydata))
    Gcount <- numeric(nrow(mydata))
    SGcount <- numeric(nrow(mydata))
    Gorig <- numeric(nrow(mydata))
    PCharge <- double(nrow(mydata))
    Inum <- double(nrow(mydata))
    MassMatch <- double(nrow(mydata))
    PeakLabel <- character(nrow(mydata))
    SaltLabel <- character(nrow(mydata))
    ILabel <- character(nrow(mydata))

    Pidx1 <- numeric(nrow(corrdata))
    Pidx2 <- numeric(nrow(corrdata))
    totCorrs <- numeric(nrow(corrdata))
    usecorr <- double(nrow(corrdata))

    idx1 <- 0
    idx2 <- 0

    # Sort mydata peak data by accurate mass
    # mydata <- mydata[order(mydata$mz), ]

    # Initialise - read sorted mass and RT data into an Array for fast search
    compMF <- as.double(0.0)
    # for (int i=1;i<mydata.size();i++) {
    # String inStr=mydata.get(i);
    # String [] tempdata = inStr.split("\t");
    #
    # System.out.println(tempdata[0]);
    #
    # int pval = Integer.parseInt(tempdata[0]);
    # double massval = new Double(tempdata[1]).doubleValue();
    # double tval = new Double(tempdata[2]).doubleValue();
    # double respval = new Double(tempdata[3]).doubleValue();
    # PeakName[i]=pval;
    PeakName <- mydata[, 1]
    # sortPeaks[i]=PeakName[i];
    sortPeaks <- PeakName
    # massdata[i]=massval;
    massdata <- mydata[, 2]
    # MassMatch[i]=massval;
    MassMatch <- massdata
    # matchMF[i] <- 0
    matchMF <- double(nrow(mydata))
    # Check if data in file has been sorted by accurate mass by comparing
    # mass of one peak with the mass of the preceding peak - might not need
    # this because mydata has been sorted by mz value above
    m <- shift.column(data=mydata, columns="mz")
    m <- within(m, "massCheck" <- FALSE)
    m[(m$mz.Shifted < m$mz), "massCheck"] <- TRUE
    if(any(m$massCheck)) {
        tempStr1="Check file is not sorted by Accurate Mass"
        stop("Check file is not sorted by Accurate Mass")
    }
    rm(m)

    # compMF=massval;
    compMF <- mydata[, 2]
    # RTdata[i]=tval;
    RTdata <- mydata[, 3]
    # medval[i]=respval;
    medval <- mydata[, 4]
    # PeakMatch[i]=pval; // default - initialise Peak to itself before matching
    PeakMatch <- mydata[, 1]
    # PeakIndex[i]=i;
    PeakIndex <- seq(1, nrow(mydata))
    # Gcount already initialised with zeros
    # Gcount[i]=0;
    # Gcount <- integer(nrow(mydata))
    # Gorig[i]=0;
    Gorig <- integer(nrow(mydata))
    # SGcount[i]=0;
    SGcount <- integer(nrow(mydata))
    # PCharge[i]=1;
    PCharge <- integer(nrow(mydata))
    # Inum[i]=0;
    Inum <- integer(nrow(mydata))
    # PeakLabel[i]="A"; // default - initialise all Peaks to A - any/all
    PeakLabel <- character(nrow(mydata))
    PeakLabel <- sapply(PeakLabel, function(x) x <- 'A')
    # WorkLabel[i]="A";
    WorkLabel <- character(nrow(mydata))
    WorkLabel <- sapply(WorkLabel, function(x) x <- 'A')
    # SaltLabel[i]="";
    # Already created above
    # ILabel[i]="N";
    ILabel <- character(nrow(mydata))
    ILabel <- sapply(ILabel, function(x) x <- 'N')
    # useLabel[i]=0;
    useLabel <- rep(0, nrow(mydata))

    for(i in 1:length(RTdata)) {
        if (RTdata[i]<rtmin | RTdata[i]>rtmax){
            useLabel[i] <- 9
            PeakLabel[i] <- "E"
        }
    }
    # }

    write.table(PeakLabel, file="PeakLabel_R.csv", sep=",", col.names=FALSE, row.names=FALSE, quote=FALSE)

    # Find number of characters in first line of mydata file
    # headLen=mydata.get(0).length();
    conn <- file("/Users/peterli/PUTMEDID_LCMS_v1.01/Study_pospeaks.txt", "r")
    first_line <- readLines(conn, n=1)
    headLen <- nchar(first_line)
    close(conn)
    # int[] ordpname = ShellSortValues(PeakName, true);
    ordpname <- order(PeakName, method = "shell")
    write.table(ordpname, file="ordpname_R.csv", sep=",", col.names=FALSE, row.names=FALSE, quote=FALSE)
    # Arrays.sort(sortPeaks);
    # Sort sortPeaks in ascending order
    sortPeaks <- sort(sortPeaks)
    write.table(sortPeaks, file="sortPeaks_R.csv", sep=",", col.names=FALSE, row.names=FALSE, quote=FALSE)

    ccount <- 0
    oidx1 <- 0
    oidx2 <- 0
    compMF <- 0
    ncorr <- 0
    # ncorr <- 1
    tdiff <- 0
    counter <- 1
    # Process peak correlation data
    for (i in 1:nrow(corrdata)) {
        # if(i == 3)
        #     break
        pval1 <- as.integer(corrdata[i, 1])
        write(pval1, file="pval1_R.csv", append=TRUE)
        pval2 <- as.integer(corrdata[i, 2])
        tempcorr <- as.double(corrdata[i, 3])
    #     idx1=Arrays.binarySearch(sortPeaks,pval1);
        # Need to check mydata[,1] has been sorted in ascending order
        idx1 <- match(sortPeaks, pval1)
        idx1 <- which(idx1 == 1)
        # write.table(idx1, file="idx1_R.csv", sep=",", col.names=FALSE, row.names=FALSE, quote=FALSE, append=TRUE)
        # idx2=Arrays.binarySearch(sortPeaks,pval2);
        idx2 <- match(sortPeaks, pval2)
        idx2 <- which(idx2 == 1)
        oidx1 <- ordpname[idx1]
        # write.table(oidx1, file="oidx1_R.csv", sep=",", col.names=FALSE, row.names=FALSE, quote=FALSE, append=TRUE)
        oidx2 <- ordpname[idx2]
        tdiff <- RTdata[oidx1]-RTdata[oidx2]
        if (tdiff<0) {
            tdiff <- -tdiff
            # print(paste0("tdiff less than zero: ", tdiff))
        }
        # write.table(tdiff, file="tdiff_R.csv", sep=",", col.names=FALSE, row.names=FALSE, quote=FALSE, append=TRUE)
        if (tdiff<rtlim & useLabel[oidx1]<9 & useLabel[oidx2]<9) {
            ncorr <- ncorr+1
            stuff <- c(tdiff, rtlim, useLabel[oidx1], useLabel[oidx2], ncorr, oidx1, oidx2)
            check <- rbind(check, stuff)

            Pidx1[ncorr] <- oidx1
            Pidx2[ncorr] <- oidx2
            usecorr[ncorr] <- tempcorr
            if (oidx1>oidx2) {
                Pidx1[ncorr] <- oidx2
                Pidx2[ncorr] <- oidx1
            }
            totCorrs[ncorr] <- ((Pidx1[ncorr]*50000)+Pidx2[ncorr])
            # ncorr <- ncorr+1
        }
        counter <- counter + 1
    }

    # })

    # write.table(pval1, file="pval1_R.csv", sep = ",", quote=FALSE, col.names = TRUE, row.names = FALSE)
    write.table(Pidx1, file="Pidx1_R.csv", sep = ",", quote=FALSE, col.names = FALSE, row.names = FALSE)
    write.table(totCorrs, file="totCorrs_R.csv", sep = ",", quote=FALSE, col.names = FALSE, row.names = FALSE)

    # int[] finCorrs=new int[ncorr];
    # for (int i=0;i<ncorr;i++) {
    #     finCorrs[i]=totCorrs[i];
    # }
    finCorrs <- numeric(ncorr)
    for (i in 1:ncorr) {
        finCorrs[i] <- totCorrs[i]
    }
    write.table(finCorrs, file="finCorrs_R.csv", sep = "\t", quote=FALSE, col.names = FALSE, row.names = FALSE)

    # ShellSortValues function sorts values and returns the indices of this
    # sorted order. It does not return the actual values sorted in requested
    # order!
    # int[] ordcorr = ShellSortValues(finCorrs, true);
    ordcorr <- order(finCorrs, method = "shell")
    write.table(ordcorr, file="ordcorr_R.csv", sep = ",", quote=FALSE, col.names = FALSE, row.names = FALSE)

    # Sort MF reference data by accurate mass
    refdata <- refdata[order(refdata$Amass), ]

    # Read sorted MF Reference data into an Array for fast search
    rcompMF <- as.double(0.0)
    # for (i in 1 : nrow(refdata)) {
    #     tcount <- tcount + 1
    #     rheadLen <- length(refdata[i, ])
    #     inStr <- substr(refdata[i, ], 0, rheadLen - 1)
    #     tempdata <- strsplit(inStr, "\t")
    #     tval <- as.double(refdata[i, 1])
    #     MFdata[i] <- tval
    #     MFformula[i] <- refdata[i, 2]
    #     This part is not required because refdata has been sorted above
    #     if (tval < rcompMF) {
    #         tempStr3 <- "Reference file data needs to be sorted by Accurate Mass"
    #         rsortedRefMF <- FALSE
    #     }
    #     rcompMF <- tval
    # }
    MFdata <- as.double(refdata[, 1])
    MFformula <- as.double(refdata[, 2])

    print(paste0("sortedMF: ", sortedMF))
    print(paste0("rsortedRefMF: ", rsortedRefMF))

    print(paste0("Gcount size: ", length(Gcount)))
    if (sortedMF & rsortedRefMF) {
        # Create framework of correlated Peaks
        tnum <- 0
        hno <- 0
        tno <- 0
        print(paste0("ncorr: ", ncorr))
        # loop thru correlations
        for (j in 1:ncorr) {
            # if(j == 10)
            #     break
            # print(paste0("j: ", j))
            # print(paste0("ordcorr[j]] ", ordcorr[j]))
            # print(paste0("usecorr[ordcorr[j]]: ", usecorr[ordcorr[j]]))
            if (usecorr[ordcorr[j]]>0.94) {
                # print(paste0("usecorr[ordcorr[j]]: ", usecorr[ordcorr[j]]))
                i <- Pidx1[ordcorr[j]]
                # print(paste0("i: ", i))
                tnum <- Pidx2[ordcorr[j]]
                if (Gcount[i]==0 & Gcount[tnum]==0) {
                    gno <- gno+1
                    Gcount[i] <- gno
                    Gcount[tnum] <- gno
                }
                if (Gcount[i]>0 & Gcount[tnum]==0) {
                    Gcount[tnum] <- Gcount[i]
                }
                if (Gcount[i]==0 & Gcount[tnum]>0) {
                    Gcount[i] <- Gcount[tnum]
                }
                if (Gcount[i]>0 & Gcount[tnum]>0) {
                    if (Gcount[i] != Gcount[tnum]) {
                        if (Gcount[i]>Gcount[tnum]) {
                            hno <- Gcount[i]
                            tno <- Gcount[tnum]
                        }
                        if (Gcount[i]<Gcount[tnum]) {
                            hno <- Gcount[tnum]
                            tno <- Gcount[i]
                        }
                        for (k in 1:nrow(mydata)) {
                            if (Gcount[k]==hno) {
                                Gcount[k] <- tno
                            }
                        }
                    }
                }
            }
        }
        write.table(Gcount, file="Gcount_R.csv", sep = ",", col.names = FALSE, row.names = FALSE)

        ogno <- as.integer(gno)
        print(paste0("ogno: ", ogno))
        for (i in 1:nrow(mydata)) {
            Gorig[i] <- Gcount[i]
        }

        cbind(listofdata, substr(mydata[1], 0, headLen - 1), "MetGroup", "PeakType", "Salt", "useLabel", "MatchedMF", "Adduct", "ppm", "Ionmode")

        # Initialise for Fourier/Isotope Peaks handling
        holdval <- as.double(0.0)
        iratio <- as.double(0.0)
        # tdiff <- as.double(0.0)
        tmass <- vector(mode = "double", length = nrow(indata))
        tlabel <- vector(mode = "character", length = nrow(indata))
        llabel <- vector(mode = "character", length = nrow(indata))
        ulabel <- vector(mode = "character", length = nrow(indata))
        llim <- vector(mode = "double", length = nrow(indata))
        ulim <- vector(mode = "double", length = nrow(indata))
        tcharge <- vector(mode = "double", length = nrow(indata))
        tperc <- vector(mode = "double", length = nrow(indata))
        tempval <- 0
        # Initialise and sort out isotope info
        # for (i in 1 < nrow(indata)) {
        # tempdata <- strsplit(indata[i], "\t")
        # clabel <- tempdata[1]
        clabel <- indata[,1]
        # cmass <- as.double(tempdata[2])
        cmass <- as.double(indata[,2])
        # cllabel <- as.character(tempdata[3])
        cllabel <- as.character(indata[,3])
        # culabel <- as.character(tempdata[4])
        culabel <- as.character(indata[,4])
        # cllim <- as.double(tempdata[5])
        cllim <- as.double(indata[,5])
        # culim <- as.double(tempdata[6])
        culim <- as.double(indata[,6])
        # ccharge <- as.double(tempdata[7])
        ccharge <- as.double(indata[,7])
        # cperc <- as.double(tempdata[8])
        cperc <- as.double(indata[,8])
        tlabel <- clabel
        tmass <- cmass
        llabel <- cllabel
        ulabel <- culabel
        llim <- cllim
        ulim <- culim
        tcharge <- ccharge
        tperc <- cperc
        # }


        # Test for Fourier peaks - set peaks to F and useLabel=8
        holdPeak <- 0
        mtol <- 0.6
        for (m in 1:ncorr) {
            i <- Pidx1[ordcorr[m]]
            k <- Pidx2[ordcorr[m]]
            if (i > holdPeak) {
                if (useLabel[i] < 2 & useLabel[k] < 2) {
                    if (Gcount[i] == 0 |
                        Gcount[k] == 0 |
                        (Gcount[i] == Gcount[k])) {
                        # Check whether mass differences are within mass diff range
                        if (massdata[k] > massdata[i] + mtol) {
                            holdPeak <- i
                        }
                        else {
                            if (usecorr[ordcorr[m]] > 0.8) {
                                iratio <- 100 * medval[k] / medval[i]  # peak 2 % of peak 1
                                if (iratio > 1000) {
                                    useLabel[i] <- 8
                                    PeakLabel[i] <- "F"
                                    PeakMatch[i] <- PeakName[k]
                                    PeakIndex[i] <- k
                                    MassMatch[i] <- massdata[k]
                                    useLabel[k] <- 1
                                    PeakLabel[k] <- "M"
                                }
                                if (iratio < 1) {
                                    useLabel[i] <- 1
                                    PeakLabel[i] <- "M"
                                    PeakLabel[k] <- "F"
                                    PeakMatch[k] <- PeakName[i]
                                    PeakIndex[k] <- i
                                    MassMatch[k] <- massdata[i]
                                    useLabel[k] <- 8
                                }
                            }
                        }
                    }
                }
            }
        }
        write.table(MassMatch, file="MassMatch_R.csv", sep = ",", col.names = FALSE, row.names = FALSE)

        # check for isotopes
        # set tol to acceptable standard dependent on small diffs between masdiffs
        calcc <- as.double(0.0)
        calcdiff <- as.double(0.0)
        lnum <- 0
        hnum <- 0
        for (j in 1:nrow(indata)) {
            holdPeak <- 0
            if (tmass[j] < 0.26) {
                mtol <- 0.00075
            }
            else {
                mtol <- 0.0015
            }
            for (m in 1:ncorr) {
                i <- Pidx1[ordcorr[m]]
                k <- Pidx2[ordcorr[m]]
                if (i > holdPeak) {
                    if (useLabel[i] < 8 & useLabel[k] < 2) {
                        if (Gcount[i] == 0 |
                            Gcount[k] == 0 |
                            (Gcount[i] == Gcount[k])) {
                            if (massdata[k] > massdata[i] + tmass[j] + mtol) {
                                holdPeak <- i
                            }
                            else {
                                if ((massdata[k] - massdata[i] > tmass[j] - mtol) & usecorr[ordcorr[m]] > corrlim) {
                                    iratio <- 100 * medval[k] / medval[i]  # peak 2 % of peak 1
                                    if ((iratio / tperc[j]) > llim[j] & (iratio / tperc[j]) < ulim[j] * tcharge[j]) {
                                        lnum <- i
                                        hnum <- k
                                        if (useLabel[lnum] == 0) {
                                            PeakLabel[lnum] <- llabel[j]
                                            if (tcharge[j] == 1) {
                                                useLabel[lnum] <- 1
                                            }
                                            else {
                                                useLabel[lnum] <- 7
                                            }
                                        }
                                        useLabel[hnum] <- 7
                                        PeakLabel[hnum] <- ulabel[j]
                                        PeakMatch[hnum] <- PeakName[i]
                                        PeakIndex[hnum] <- lnum
                                        MassMatch[hnum] <- massdata[lnum]
                                        PCharge[lnum] <- tcharge[j]
                                        PCharge[hnum] <- tcharge[j]
                                        ILabel[lnum] <- tlabel[j]
                                        ILabel[hnum] <- tlabel[j]
                                        Inum[lnum] <- iratio / tperc[j]
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        write.table(Inum, file="Inum_R.csv", sep = ",", col.names = FALSE, row.names = FALSE)


        # Simplified Dimer Code - H only
        holdPeak <- 0
        for (m in 1:ncorr) {
            i <- Pidx1[ordcorr[m]]
            k <- Pidx2[ordcorr[m]]
            if (i > holdPeak) {
                if (useLabel[i] < 2 && useLabel[k] < 2) {  # then A or M
                    if (Gcount[i] == 0 |
                        Gcount[k] == 0 |
                        (Gcount[i] == Gcount[k])) {
                        mtol <- (massdata[i] * 2 * tol) / 1000000
                        if (mtol < masslim) {
                            mtol <- masslim
                        }
                        newval <- massdata[i] - 1.00727 * etol
                        newval <- (newval * 2) + 1.00727 * etol - mtol
                        if (massdata[k] > newval + 2 * mtol) {
                            holdPeak <- i
                        }
                        else {
                            if ((massdata[k] - newval > 0) &
                                (massdata[k] - newval < 2 * mtol) &
                                usecorr[ordcorr[m]] > 0.94) {
                                useLabel[i] <- 2
                                PeakLabel[i] <- "M"
                                PeakLabel[k] <- "D"
                                useLabel[k] <- 2
                                PeakMatch[k] <- PeakName[i]
                                PeakIndex[k] <- i
                                MassMatch[k] <- massdata[i]
                            }
                        }
                    }
                }
            }
        }


        # Identify and annotate Adduct Peaks
        # initialise adduct info and look for simple adducts and then combined adducts
        cmode <- 0
        Amode <- vector(mode = "integer", length = nrow(adddata))
        ALabel <- vector(mode = "character", length = nrow(adddata))
        iLabel <- vector(mode = "character", length = nrow(adddata))
        Amass <- vector(mode = "double", length = nrow(adddata))
        # for (i in 1:nrow(adddata)) {
        # tempdata <- strsplit(adddata[i,], "\t")
        # cmode <- as.integer(tempdata[1])
        cmode <- as.integer(adddata[,1])
        # calabel <- as.character(tempdata[2])
        calabel <- as.character(adddata[,2])
        # cilabel <- as.character(tempdata[3])
        cilabel <- as.character(adddata[,3])
        # camass <- as.double(tempdata[4])
        camass <- as.double(adddata[,4])

        for(i in 1:length(cmode)) {
            if (cmode[i] == 3) {
                cmode[i] <- 5
            }
            if (cmode[i] == 1 & etol > 0) {
                if (calabel[i] == "NH3") {
                    cmode[i] <- 3
                }
                if (calabel[i] == "HCOOH") {
                    cmode[i] <- 4
                }
            }
            Amode[i] <- cmode[i]
        }

        ALabel <- calabel
        iLabel <- cilabel
        Amass <- camass

        write.table(Amode, file="Amode_R.csv", sep = ",", col.names = FALSE, row.names = FALSE)

        tempLabel <- ""
        for (j in 1:nrow(adddata)) {
            holdPeak <- 0
            for (m in 1:3) {
                for (c in 1:ncorr) {
                    i <- Pidx1[ordcorr[c]]
                    k <- Pidx2[ordcorr[c]]
                    mtol <- (massdata[i] * 2 * tol) / 1000000
                    if (mtol < masslim) {
                        mtol <- masslim
                    }
                    if (i > holdPeak & Amode[j] == m) {
                        if (useLabel[i] < 4 & useLabel[k] < 3) {  # then A or M or M adduct
                            if (Gcount[i] == 0 |
                                Gcount[k] == 0 |
                                useLabel[i] < 2 |
                                useLabel[k] < 2 |
                                (Gcount[i] == Gcount[k])) {
                                if (massdata[k] > massdata[i] + Amass[j] + mtol) {
                                    holdPeak <- i
                                }
                                else {
                                    if ((massdata[k] - massdata[i] > Amass[j] - mtol) & usecorr[ordcorr[c]] > corrlim) {
                                        lnum <- i
                                        hnum <- k
                                        tempLabel <- ALabel[j]
                                        if (useLabel[lnum] < 2) {
                                            useLabel[lnum] <- 2
                                            PeakLabel[lnum] <- "M"
                                        }
                                        useLabel[hnum] <- 3
                                        PeakLabel[hnum] <- paste0(PeakLabel[lnum], "_", tempLabel)
                                        PeakMatch[hnum] <- PeakName[lnum]
                                        PeakIndex[hnum] <- lnum
                                        MassMatch[hnum] <- massdata[lnum]
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        write.table(useLabel, file="useLabel_R.csv", sep = ",", col.names = FALSE, row.names = FALSE)

        # Relook for special cases
        for (j in 1:nrow(adddata)) {
            if (Amode[j] > 2) {
                holdPeak <- 0
                for (c in 1:ncorr) {
                    i <- Pidx1[ordcorr[c]]
                    k <- Pidx2[ordcorr[c]]
                    mtol <- (massdata[i] * 2 * tol) / 1000000
                    if (mtol < masslim) {
                        mtol <- masslim
                    }
                    if (i > holdPeak) {
                        if (useLabel[i] < 4 & useLabel[k] < 3) {  # then A or M or M adduct
                            if (Gcount[i] == 0 |
                                Gcount[k] == 0 |
                                useLabel[i] < 2 |
                                useLabel[k] < 2 |
                                (Gcount[i] == Gcount[k])) {
                                if (massdata[k] > massdata[i] + Amass[j] + mtol) {
                                    holdPeak <- i
                                }
                                else {
                                    if ((massdata[k] - massdata[i] > Amass[j] - mtol) & usecorr[ordcorr[c]] > corrlim) {
                                        lnum <- i
                                        hnum <- k
                                        tempLabel <- ALabel[j]
                                        if (Amode[j] == 3) {
                                            if (useLabel[i] < 2 &
                                                useLabel[k] < 3 &
                                                (usecorr[ordcorr[c]] > 0.97 | (medval[i] < medval[k] / 2))) {
                                                lnum <- k
                                                hnum <- i
                                                tempLabel <- paste0("[", ALabel[j], "]")
                                            }
                                            useLabel[lnum] <- 2
                                            PeakLabel[lnum] <- "M"
                                            useLabel[hnum] <- 3
                                            PeakLabel[hnum] <- paste0(PeakLabel[lnum], "_", tempLabel)
                                            PeakMatch[hnum] <- PeakName[lnum]
                                            PeakIndex[hnum] <- lnum
                                            MassMatch[hnum] <- massdata[lnum]
                                        }
                                        if (Amode[j] == 4 & useLabel[k] < 3 & useLabel[i] < 2) {
                                            tempLabel <- paste0("[", ALabel[j], "]")
                                            lnum <- k
                                            hnum <- i
                                            useLabel[lnum] <- 2
                                            PeakLabel[lnum] <- "M"
                                            useLabel[hnum] <- 3
                                            PeakLabel[hnum] <- paste0(PeakLabel[lnum], "_", tempLabel)
                                            PeakMatch[hnum] <- PeakName[lnum]
                                            PeakIndex[hnum] <- lnum
                                            MassMatch[hnum] <- massdata[lnum]
                                        }
                                        if (Amode[j] == 5 &
                                            etol < 0 &
                                            useLabel[k] < 3 &
                                            useLabel[i] < 3) {
                                            useLabel[i] <- 2
                                            PeakLabel[i] <- paste0("M_", iLabel[j])
                                            useLabel[k] <- 3
                                            PeakLabel[k] <- paste0("M_", tempLabel)
                                            PeakMatch[k] <- PeakName[i]
                                            PeakIndex[k] <- i
                                            MassMatch[k] <- massdata[i]
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        # RETRY with greater tolerance for isotopes for unmatched features
    #     for (j in 1: nrow(indata)) {
    #         holdPeak <- 0
    #         if (tmass[j] < 0.26) {
    #             mtol <- 0.0015
    #         }
    #         else {
    #             mtol <- 0.002
    #         }
    #         for (m in 1:ncorr) {
    #             i <- Pidx1[ordcorr[m]]
    #             k <- Pidx2[ordcorr[m]]
    #             if (i > holdPeak) {
    #                 if (useLabel[i] < 8 & useLabel[k] < 2) {
    #                     if (Gcount[i] == 0 |
    #                         Gcount[k] == 0 |
    #                         (Gcount[i] == Gcount[k])) {
    #                         if (massdata[k] > massdata[i] + tmass[j] + 1.5 * mtol) {
    #                             holdPeak <- i
    #                         }
    #                         else {
    #                             if ((massdata[k] - massdata[i] > tmass[j] - mtol) & usecorr[ordcorr[m]] > 0.75) {
    #                                 iratio <- 100 * medval[k] / medval[i]
    #                                 if ((iratio / tperc[j]) > llim[j] && (iratio / tperc[j]) < ulim[j] * tcharge[j]) {
    #                                     lnum <- i
    #                                     hnum <- k
    #                                     if ((ILabel[lnum] == tlabel[j]) | ILabel[lnum] == "N") {
    #                                         if (useLabel[lnum] == 0) {
    #                                             PeakLabel[lnum] <- llabel[j]
    #                                             if (tcharge[j] == 1) {
    #                                                 useLabel[lnum] <- 1
    #                                             }
    #                                             else {
    #                                                 useLabel[lnum] <- 6
    #                                             }
    #                                         }
    #                                     }
    #                                     useLabel[hnum] <- 6
    #                                     PeakLabel[hnum] <- ulabel[j]
    #                                     PeakMatch[hnum] <- PeakName[i]
    #                                     PeakIndex[hnum] <- lnum
    #                                     MassMatch[hnum] <- massdata[lnum]
    #                                     PCharge[lnum] <- tcharge[j]
    #                                     PCharge[hnum] <- tcharge[j]
    #                                     ILabel[lnum] <- tlabel[j]
    #                                     ILabel[hnum] <- tlabel[j]
    #                                     Inum[lnum] <- iratio / tperc[j]
    #                                 }
    #                             }
    #                         }
    #                     }
    #                 }
    #             }
    #         }
    #     }
    #     write.table(PCharge, file="PCharge_R.csv", sep = ",", col.names = FALSE, row.names = FALSE)
    #
    #     # Write out annotated group info into existing groups
    #     tnum <- 0
    #     hno <- 0
    #     tno <- 0
    #     for (i in 1:nrow(mydata)) {
    #         WorkLabel[i] <- PeakLabel[i]
    #         if (i != PeakIndex[i]) {
    #             tnum <- PeakIndex[i]
    #             if (Gcount[i] == 0 & Gcount[tnum] == 0) {
    #                 gno <- gno + 1
    #                 Gcount[i] <- gno
    #                 Gcount[tnum] <- gno
    #             }
    #             if (Gcount[i] > 0 & Gcount[tnum] == 0) {
    #                 Gcount[tnum] <- Gcount[i]
    #             }
    #             if (Gcount[i] == 0 & Gcount[tnum] > 0) {
    #                 Gcount[i] <- Gcount[tnum]
    #             }
    #             if (Gcount[i] > 0 & Gcount[tnum] > 0) {
    #                 if (Gcount[i] != Gcount[tnum]) {
    #                     if (Gcount[i] > Gcount[tnum]) {
    #                         hno <- Gcount[i]
    #                         tno <- Gcount[tnum]
    #                     }
    #                     if (Gcount[i] < Gcount[tnum]) {
    #                         hno <- Gcount[tnum]
    #                         tno <- Gcount[i]
    #                     }
    #                     for (k in 1:nrow(mydata)) {
    #                         if (Gcount[k] == hno) {
    #                             Gcount[k] <- tno
    #                         }
    #                     }
    #                 }
    #             }
    #         }
    #     }
    #     write.table(WorkLabel, file="WorkLabel_R.csv", sep = ",", col.names = FALSE, row.names = FALSE)
    #
    #     # SubGroup related Peaks
    #     tnum <- 0
    #     hno <- 0
    #     tno <- 0
    #     Sgno <- 0
    #     for (i in 1:nrow(mydata)) {
    #         if (i != PeakIndex[i]) {
    #             tnum <- PeakIndex[i]
    #             if (SGcount[i] == 0 & SGcount[tnum] == 0) {
    #                 Sgno <- Sgno + 1
    #                 SGcount[i] <- Sgno
    #                 SGcount[tnum] <- Sgno
    #             }
    #             if (SGcount[i] > 0 & SGcount[tnum] == 0) {
    #                 SGcount[tnum] <- SGcount[i]
    #             }
    #             if (SGcount[i] == 0 & SGcount[tnum] > 0) {
    #                 SGcount[i] <- SGcount[tnum]
    #             }
    #             if (SGcount[i] > 0 & SGcount[tnum] > 0) {
    #                 if (SGcount[i] != SGcount[tnum]) {
    #                     if (SGcount[i] > SGcount[tnum]) {
    #                         hno <- SGcount[i]
    #                         tno <- SGcount[tnum]
    #                     }
    #                     if (SGcount[i] < SGcount[tnum]) {
    #                         hno <- SGcount[tnum]
    #                         tno <- SGcount[i]
    #                     }
    #                     for (k in 1:nrow(mydata)) {
    #                         if (SGcount[k] == hno) {
    #                             SGcount[k] <- tno
    #                         }
    #                     }
    #                 }
    #             }
    #         }
    #     }
    #     write.table(SGcount, file="SGcount_R.csv", sep = ",", col.names = FALSE, row.names = FALSE)
    #
    #     # annotate Cl info
    #     for (k in 1:nrow(mydata)) {
    #         if (useLabel[k] < 6) {
    #             if (ILabel[k] == "Cl") {
    #                 if (Inum[k] > 0.56) {
    #                     Inum[k] <- Inum[k] - 0.4
    #                     tfloor <- 1 + as.integer(Inum[k])
    #                     if (PeakLabel[k] == "M" & ILabel[k] == "Cl") {
    #                         if (tfloor == 2 | tfloor == 3) {
    #                             PeakLabel[k] <- paste0("M_MCl", "x", tfloor)
    #                         }
    #                         else {
    #                             PeakLabel[k] <- "M_MCl"
    #                         }
    #                         useLabel[k] <- 2
    #                     }
    #                     PeakLabel[k] <- paste0(PeakLabel[k], " ", ILabel[k], " ", tfloor)
    #                     Inum[k] <- Inum[k] + 0.4
    #                 }
    #             }
    #         }
    #     }
    #
    #
    #     # Write out annotated Peak data
    #     myStr <- ""
    #     UseStr <- ""
    #     tnum <- 0
    #     igno <- as.integer(Sgno + 1)
    #     # print("OK here!")
    #     # glist <- integer[igno]
    #     glist <- vector(mode = "integer", length = igno)
    #     # print("OK here!")
    #     mygcount <- 0
    #
    #     for (j in 1:igno) {
    #         glist[j] <- mygcount
    #         mygcount <- mygcount + 1
    #     }
    #     # print("OK here!")
    #     # Calculate salt mass
    #     # SaltMass <- double[200]
    #     SaltMass <- vector(mode = "double", length = 200)
    #     # SaltType <- character[200]
    #     SaltType <- vector(mode = "character", length = 200)
    #     # print("OK here!")
    #     NaCl <- as.double(57.95862)
    #     Kdiff <- as.double(15.97394)
    #     Fodiff <- as.double(10.0288)
    #     Cldiff <- as.double(1.9972)
    #     # print("OK here!")
    #     scount <- 0
    #     tempmass <- as.double(0.0)
    #     holdmass <- as.double(0.0)
    #     Kmass <- as.double(0.0)
    #     massion <- as.double(0.0)
    #     saltout <- ""
    #     formout <- ""
    #     salt <- "(NaCl)"
    #     iso <- ""
    #     niso <- "Cl-"
    #     piso <- "Na+"
    #     massion <- as.double(0.0)
    #     istart <- 0
    #     iend <- 0
    #     if (etol > 0) {
    #         iso <- paste0(iso, piso)
    #         massion <- massion + 22.98977
    #         istart <- istart + 2
    #     }
    #     else {
    #         iso <- iso + niso
    #         massion <- 34.96885
    #         istart <- istart + 1
    #     }
    #     # print("OK here!")
    #     for (i in 1:6) {
    #         # print("OK here!")
    #         print(i)
    #         tempmass <- (NaCl * i) + massion
    #         if (i < 2) {
    #             saltout <- paste0(salt, iso)
    #             formout <- paste0("(HCOONa)", iso)
    #             print("i is less than 2")
    #         }
    #         else {
    #             saltout <- paste0(salt, i, iso)
    #             formout <- paste0("HCOONa_", salt, (i - 1), iso)
    #         }
    #
    #         # print("OK here!")
    #         SaltMass[scount] <- tempmass
    #         SaltType[scount] <- saltout
    #         scount <- scount + 1
    #         holdmass <- tempmass
    #         Kmass <- tempmass
    #         iend <- i + 1
    #         # print("OK here!")
    #         if (i > 3) {
    #             print("i is more than 3")
    #             iend <- 4
    #         }
    #
    #         for (j in 1:iend) {
    #             tempmass <- tempmass + Cldiff
    #             SaltMass[scount] <- tempmass
    #             # print("OK here!")
    #             SaltType[scount] <- paste0(saltout, " Isotope")
    #             scount <- scount + 1
    #
    #         }
    #         # add in HCOONa + isotopes
    #         holdmass <- holdmass + Fodiff
    #         SaltMass[scount] <- holdmass
    #         SaltType[scount] <- formout
    #         scount <- scount + 1
    #         j <- istart
    #         for (j in 1:iend) {
    #             holdmass <- holdmass + Cldiff
    #             SaltMass[scount] <- holdmass
    #             SaltType[scount] <- paste0(formout, " Isotope")
    #             scount <- scount + 1
    #         }
    #         # add in K adducts
    #         Kmass <- Kmass + Kdiff
    #         if (etol > 0 & i == 1) {
    #             print("In: etol > 0 & i == 1")
    #             SaltMass[scount] <- Kmass
    #             SaltType[scount] <- "(NaCl)K+"
    #             scount <- scount + 1
    #             SaltMass[scount] <- Kmass + Cldiff
    #             SaltType[scount] <- "(NaCl)K+ Isotope"
    #             scount <- scount + 1
    #             SaltMass[scount] <- Kmass + Fodiff
    #             SaltType[scount] <- "HCOONa_K+"
    #             scount <- scount + 1
    #             SaltMass[scount] <- Kmass + Kdiff
    #             SaltType[scount] <- "(KCl)K+"
    #             scount <- scount + 1
    #             SaltMass[scount] <- Kmass + Kdiff + Cldiff
    #             SaltType[scount] <- "(KCl)K+ Isotope"
    #             scount <- scount + 1
    #         }
    #         if (etol > 0 & i == 2) {
    #             print("In: etol > 0 & i == 2")
    #             SaltMass[scount] <- Kmass
    #             SaltType[scount] <- "(NaCl)2K+"
    #             scount <- scount + 1
    #             SaltMass[scount] <- Kmass + Cldiff
    #             SaltType[scount] <- "(NaCl)2K+ Isotope"
    #             scount <- scount + 1
    #             SaltMass[scount] <- Kmass + (Cldiff) * 2
    #             SaltType[scount] <- "(NaCl)2K+ Isotope"
    #             scount <- scount + 1
    #             Kmass <- Kmass + Kdiff
    #             SaltMass[scount] <- Kmass
    #             SaltType[scount] <- "(KCl)2Na+"
    #             scount <- scount + 1
    #             SaltMass[scount] <- Kmass + Cldiff
    #             SaltType[scount] <- "(KCl)2Na+ Isotope"
    #             scount <- scount + 1
    #             SaltMass[scount] <- Kmass + (Cldiff) * 2
    #             SaltType[scount] <- "(KCl)2Na+ Isotope"
    #             scount <- scount + 1
    #             Kmass <- Kmass + Kdiff
    #             SaltMass[scount] <- Kmass
    #             SaltType[scount] <- "(KCl)2K+"
    #             scount <- scount + 1
    #             SaltMass[scount] <- Kmass + Cldiff
    #             SaltType[scount] <- "(KCl)2K+ Isotope"
    #             scount <- scount + 1
    #             SaltMass[scount] <- Kmass + (Cldiff) * 2
    #             SaltType[scount] <- "(KCl)2K+ Isotope"
    #             scount <- scount + 1
    #         }
    #         if (etol < 0 & i == 1) {
    #             print("In: etol < 0 & i == 1")
    #             SaltMass[scount] <- Kmass
    #             SaltType[scount] <- "(KCl)Cl-"
    #             scount <- scount + 1
    #             SaltMass[scount] <- Kmass + Cldiff
    #             SaltType[scount] <- "(KCl)Cl- Isotope"
    #             scount <- scount + 1
    #             SaltMass[scount] <- Kmass + (Cldiff) * 2
    #             SaltType[scount] <- "(KCl)Cl- Isotope"
    #             scount <- scount + 1
    #             SaltMass[scount] <- Kmass + Fodiff
    #             SaltType[scount] <- "HCOOK_Cl-"
    #             scount <- scount + 1
    #             SaltMass[scount] <- Kmass + Fodiff + Cldiff
    #             SaltType[scount] <- "HCOOK_Cl- Isotope"
    #             scount <- scount + 1
    #         }
    #         if (etol < 0 & i == 2) {
    #             print("In: etol < 0 & i == 2")
    #             SaltMass[scount] <- Kmass
    #             SaltType[scount] <- "(NaKCl2)Cl-"
    #             scount <- scount + 1
    #             SaltMass[scount] <- Kmass + Cldiff
    #             SaltType[scount] <- "(NaKCl2)Cl- Isotope"
    #             scount <- scount + 1
    #             SaltMass[scount] <- Kmass + (Cldiff) * 2
    #             SaltType[scount] <- "(NaKCl2)Cl- Isotope"
    #             scount <- scount + 1
    #             SaltMass[scount] <- Kmass + Kdiff
    #             SaltType[scount] <- "(KCl)2Cl-"
    #             scount <- scount + 1
    #             SaltMass[scount] <- Kmass + Kdiff + Cldiff
    #             SaltType[scount] <- "(KCl)2Cl- Isotope"
    #             scount <- scount + 1
    #             SaltMass[scount] <- Kmass + Kdiff + (Cldiff) * 2
    #             SaltType[scount] <- "(KCl)2Cl- Isotope"
    #             scount <- scount + 1
    #             SaltMass[scount] <- Kmass + Kdiff + (Cldiff) * 3
    #             SaltType[scount] <- "(KCl)2Cl- Isotope"
    #             scount <- scount + 1
    #         }
    #     }
    #
    #     mtol <- 0.0015
    #     for (j in 1:scount) {
    #         newval <- SaltMass[j] - mtol
    #         print(paste0("newval: ", newval))
    #         print(paste0("SaltMass[j]: ", SaltMass[j]))
    #         print(paste0("mtol: ", mtol))
    #         idx1 <- 0
    #         idx2 <- 0
    #         idx1 <- match(massdata, newval)
    #         # print(paste0("massdata: ", massdata))
    #         print(paste0("newval: ", newval))
    #         # print(paste0("idx1: ", idx1))
    #         idx1 <- which(idx1 == 1)
    #         print(paste0("idx1: ", idx1))
    #         idx2 <- - idx1
    #         print(paste0("idx2: ", idx2))
    #         k <- idx2 - 1
    #         if (idx2 > 1) {
    #             for (k in 1:nrow(mydata)) {
    #                 if (massdata[k] > newval + 2 * mtol) {
    #                     break
    #                 }
    #                 else {
    #                     if (useLabel[k] < 9) {
    #                         SaltLabel[k] <- SaltType[j]
    #                         print(paste0("SaltType[j]: ", SaltType[j]))
    #                         useLabel[k] <- 5
    #                     }
    #                 }
    #             }
    #         }
    #     }
    #
    #     # Rest of code goes here
    }

    print(paste0("Size of totCorrs: ", length(totCorrs)))
    print(paste0("Size of ncorr: ", ncorr))
    print("Size of ncorr is used to create size of finCorrs")
    print(paste0("Size of corrdata: ", nrow(corrdata)))
    print(paste0("Size of finCorrs: ", length(finCorrs)))
    print(paste0("Size of ordcorr: ", length(ordcorr)))
    print(paste0("totCorrs[1]: ", totCorrs[1]))
    print(paste0("totCorrs[37424]: ", totCorrs[37424]))
    print(paste0("finCorrs[2]: ", finCorrs[2]))
    print(paste0("ordcorr[37424]: ", ordcorr[37424]))
    print(paste0("useLabel[3802]: ", useLabel[3802]))

    # Output for performing testthat tests
    return(rbind(elabel,
    lowval,
    limtol,
    rtlim,
    rtmin,
    rtmax,
    mycorrlim,
    etol,
    PeakLabel[3269],
    PeakLabel[3270],
    ordpname[3802],
    sortPeaks[3802],
    tdiff, MFdata[1],
    ncorr,
    length(totCorrs),
    totCorrs[1],
    length(finCorrs),
    finCorrs[2],
    ogno,
    MassMatch[3802],
    Inum[3801],
    Amode[15],
    useLabel[3802],
    PeakLabel[3802],
    PCharge[3787]))

}




#
#
#     # MASS MATCH AND MF DETERMINATION
#     nsym <- 18
#     iondata <- double[nsym]
#     ionsym <- character[nsym]
#     Apriority <- integer[nsym]
#     clflag <- integer[nsym]
#     # check for ions and adducts
#     ionsym[0] <- "H"
#     ionsym[1] <- "Na"
#     ionsym[2] <- "K"
#     ionsym[3] <- "NH3"
#     ionsym[4] <- "Na_Na"
#     ionsym[5] <- "NaCl"
#     ionsym[6] <- "KCl"
#     ionsym[7] <- "HCOONa"
#     ionsym[8] <- "HCOOK"
#     ionsym[9] <- "Na_HCOONa"
#     ionsym[10] <- "NaClx2"
#     ionsym[11] <- "NaCl_HCOONa"
#     ionsym[12] <- "HCOOH"
#     ionsym[13] <- "MCl"
#     ionsym[14] <- "Cl"
#     ionsym[15] <- "[NH3]"
#     ionsym[16] <- "[HCOOH]"
#     ionsym[17] <- "HCOONa_HCOONa"
#
#     # put in Alternative for Adduct if appropriate
#     Apriority[0] <- 1
#     Apriority[1] <- 4
#     Apriority[2] <- 100
#     Apriority[3] <- 0
#     Apriority[4] <- 100
#     Apriority[5] <- 100
#     Apriority[6] <- 100
#     Apriority[7] <- 9
#     Apriority[8] <- 100
#     Apriority[9] <- 100
#     Apriority[10] <- 100
#     Apriority[11] <- 100
#     Apriority[12] <- 0
#     Apriority[13] <- 14
#     Apriority[14] <- 5
#     Apriority[15] <- 0
#     Apriority[16] <- 0
#     Apriority[17] <- 100
#
#     clflag[0] <- 0
#     clflag[1] <- 0
#     clflag[2] <- 0
#     clflag[3] <- 0
#     clflag[4] <- 0
#     clflag[5] <- 1
#     clflag[6] <- 1
#     clflag[7] <- 0
#     clflag[8] <- 0
#     clflag[9] <- 0
#     clflag[10] <- 1
#     clflag[11] <- 1
#     clflag[12] <- 0
#     clflag[13] <- 1
#     clflag[14] <- 1
#     clflag[15] <- 0
#     clflag[16] <- 0
#     clflag[17] <- 0
#
#     # adjust for ion mode
#     H <- 1.007825
#     e <- 0.000549
#     Hp <- 1.00728
#     mode_adjust <- as.double(0.0)
#     eval <- 1
#     if (etol < 0) {
#         Hp <- - Hp
#     }
#
#     iondata[0] <- Hp
#     iondata[1] <- 21.98194 + Hp
#     iondata[2] <- 37.95588 + Hp
#     iondata[3] <- 17.02655 + Hp
#     iondata[4] <- iondata[1] * 2 - Hp
#     iondata[5] <- 57.95862 + Hp
#     iondata[6] <- 73.93256 + Hp
#     iondata[7] <- 67.98742 + Hp
#     iondata[8] <- 83.96136 + Hp
#     iondata[9] <- iondata[1] + iondata[7] - Hp
#     iondata[10] <- iondata[5] * 2 - Hp
#     iondata[11] <- iondata[5] + iondata[7] - Hp
#     iondata[12] <- 46.00548 + Hp
#     iondata[13] <- 57.95862 + Hp
#     iondata[14] <- 35.97668 + Hp
#     iondata[15] <- - 17.02655 + Hp
#     iondata[16] <- - 46.00548 + Hp
#     iondata[17] <- iondata[7] * 2 - Hp
#
#     multi <- 0
#     # Match for metabolites (not isotope, Fourier peaks etc)
#     # with no adducts
#     for (k in k < nrow(mydata[1,])) {
#         if (useLabel[k] < 2) {
#             myStr <- mydata[k]
#             headLen <- length(mydata[k,])
#             MFStr <- ""
#             AdductStr <- ""
#             idx1 <- 0
#             idx2 <- 0
#             multi <- 1
#             for (j in 1 < 13) {
#                 if (j < 12 | (j == 12 & etol < 1 & RTdata[k] > 80)) {
#                     mval <- massdata[k] - iondata[j]
#                     ppmStr <- ""
#                     idx1 <- binary_search(MFdata, mval - 0.1, index = TRUE)
#                     idx2 <- - idx1
#                     if (idx2 > 0) {
#                         i <- idx2 - 1
#                         for (i in 1 < nrow(refdata[1,])) {
#                             nlim <- as.double(0.0)
#                             nlim <- tol * massdata[k] / 1000000
#                             if (nlim < masslim) {
#                                 nlim <- masslim
#                             }
#                             if (MFdata[i] > (mval + 0.1))
#                             break
#                             if (mval >= (MFdata[i] - nlim) & mval <= (MFdata[i] + nlim)) {
#                                 diffval <- as.double(MFdata[i] - mval)
#                                 ppmcalc <- (MFdata[i] - mval) * 1000000 / massdata[k]
#                                 kcount = kcount + 1
#                                 if (ppmcalc < 0) {
#                                     ppmcalc <- - ppmcalc
#                                 }
#                                 ppmStr <- as.character(ppmcalc)
#                                 cbind(listofdata, substr(mydata[k], 0, headLen - 1), SGcount[k], PeakLabel[k], SaltLabel[k], useLabel[k], MFformula[i], ionsym[j], ppmStr, elabel)
#                                 matchMF[k] <- matchMF[k] + 1
#                             }
#                         }
#                     }
#                 }
#             }
#         }
#     }
#
#     # Match for metabolites with adducts (not isotope or Fourier peaks,
#     OK_Cl <- 1
#     for (k in 1 < nrow(mydata[1,])) {
#         if (useLabel[k] > 1 & useLabel[k] < 5) {
#             multi <- 1
#             OK_Cl <- 1
#             ddata <- strsplit(PeakLabel[k], "_")
#             if (length(ddata) < 4) {
#                 # look to split off Cl and multiplier x
#                 Cldata <- strsplit(PeakLabel[k], "Cl")
#                 mdata <- strsplit(PeakLabel[k], "x")
#                 if (length(Cldata) == 2) {
#                     WorkLabel[k] <- Cldata[1]
#                     multi <- as.integer(Cldata[2])
#                     if (multi > 1) {
#                         if (length(mdata) == 2) {
#                             mdata <- strsplit(Cldata[1], "x")
#                             WorkLabel[k] <- mdata[1]
#                         }
#                     }
#                 }
#                 if (length(Cldata) < 2 &
#                     length(mdata) == 2 &
#                     length(ddata) < 3) {
#                     WorkLabel[k] <- mdata[1]
#                     multi <- as.integer(mdata[2])
#                 }
#                 if (multi < 4) {
#                     for (j in 1 : nsym) {
#                         if (j == 0) {
#                             tempLabel <- "M"
#                         }
#                         else {
#                             tempLabel <- paste0("M_", ionsym[j])
#                         }
#                         if (WorkLabel[k] == tempLabel) {
#                             myStr <- mydata[k]
#                             headLen <- length(mydata[k])
#                             if ((length(Cldata) == 2 & clflag[j] == 0) | (length(mdata) == 2 & clflag[j] == 0)) {
#                                 OK_Cl <- - 1
#                             }
#                             AdductStr <- ""
#                             if (OK_Cl > 0) {
#                                 idx1 <- 0
#                                 idx2 <- 0
#                                 if (multi == 1) {
#                                     mval <- massdata[k] - iondata[j]
#                                 }
#                                 else {
#                                     mval <- massdata[k] - (iondata[j] * multi) + (multi - 1) * Hp
#                                 }
#                                 ppmStr <- ""
#                                 idx1 <- binary_search(MFdata, mval - 0.1, index = TRUE)
#                                 idx2 <- - idx1
#                                 if (idx2 > 0) {
#                                     i <- idx2 - 1
#                                     for (i in 1 : nrow(refdata[1,])) {
#                                         nlim <- as.double(0.0)
#                                         nlim <- tol * massdata[k] / 1000000
#                                         if (nlim < 0.00075) {
#                                             nlim <- 0.00075
#                                         }
#                                         if (MFdata[i] > (mval + 0.1))
#                                         break
#                                         if (mval >= (MFdata[i] - nlim) & mval <= (MFdata[i] + nlim)) {
#                                             diffval <- as.double(MFdata[i] - mval)
#                                             ppmcalc <- (MFdata[i] - mval) * 1000000 / massdata[k]
#                                             AdductStr <- ionsym[j]
#                                             if (AdductStr == "MCl") {
#                                                 AdductStr <- "NaCl"
#                                             }
#                                             if (multi > 1) {
#                                                 AdductStr <- paste0(AdductStr, "*", multi)
#                                             }
#                                             if (ppmcalc < 0) {
#                                                 ppmcalc <- - ppmcalc
#                                             }
#                                             ppmStr <- as.character(ppmcalc)
#                                             cbind(listofdata, substr(mydata[k], 0, headLen - 1), SGcount[k], PeakLabel[k], SaltLabel[k], useLabel[k], MFformula[i], AdductStr, ppmStr, elabel)
#                                             matchMF[k] <- matchMF[k] + 1
#                                         }
#                                     }
#                                 }
#                             }
#                         }
#                     }
#                 }
#             }
#         }
#     }
#
#     for (k in 1 < nrow(mydata)) {
#         if (useLabel[k] > 1 & useLabel[k] < 5) {
#             multi <- 1
#             OK_Cl <- 1
#             ddata <- strsplit(PeakLabel[k], "_")
#             if (length(ddata) < 3) {
#                 # look to split off Cl and multiplier x
#                 Cldata <- strsplit(PeakLabel[k], "Cl")
#                 mdata <- strsplit(PeakLabel[k], "x")
#                 if (length(Cldata) == 2) {
#                     WorkLabel[k] <- Cldata[1]
#                     multi <- as.integer(Cldata[2])
#                     if (multi > 1) {
#                         if (length(mdata) == 2) {
#                             mdata <- strsplit(Cldata[1], "x")
#                             WorkLabel[k] <- mdata[1]
#                         }
#                     }
#                 }
#                 if (length(Cldata) < 2 &
#                     length(mdata) == 2 &
#                     length(ddata) < 3) {
#                     WorkLabel[k] <- mdata[1]
#                     multi <- as.integer(mdata[2])
#                 }
#                 if (multi < 4) {
#                     for (j in 1 < nsym) {
#                         if (j == 0) {
#                             tempLabel <- "M"
#                         }
#                         else {
#                             tempLabel <- paste0("M_", ionsym[j])
#                         }
#                         if (WorkLabel[k] == tempLabel) {
#                             # need to check that multi is valid (Cl in adduct)
#                             myStr <- mydata[k]
#                             headLen <- length(mydata[k])
#                             if ((length(Cldata) == 2 & clflag[j] == 0) | (length(mdata) == 2 & clflag[j] == 0)) {
#                                 OK_Cl <- - 1
#                             }
#                             if (j == 0 & etol < 0) {
#                                 tnum <- 12
#                             }
#                             else {
#                                 tnum <- Apriority[j]
#                             }
#                             if (tnum < 100) {
#                                 AdductStr <- ""
#                                 if (OK_Cl > 0) {
#                                     idx1 <- 0
#                                     idx2 <- 0
#                                     if (multi == 1) {
#                                         mval <- massdata[k] - iondata[tnum]
#                                     }
#                                     else {
#                                         mval <- massdata[k] - (iondata[tnum] * multi) + (multi - 1) * Hp
#                                     }
#                                     ppmStr <- ""
#                                     idx1 <- binary_search(MFdata, mval - 0.1, index = TRUE)
#                                     idx2 <- - idx1
#                                     if (idx2 > 0) {
#                                         i <- idx2 - 1
#                                         for (i in 1 < nrow(refdata)) {
#                                             nlim <- as.double(0.0)
#                                             nlim <- tol * massdata[k] / 1000000
#                                             if (MFdata[i] > (mval + 0.1))
#                                             break
#                                             if (mval >= (MFdata[i] - nlim) & mval <= (MFdata[i] + nlim)) {
#                                                 diffval <- as.double(MFdata[i] - mval)
#                                                 ppmcalc <- (MFdata[i] - mval) * 1000000 / massdata[k]
#                                                 AdductStr <- ionsym[tnum]
#                                                 if (multi > 1) {
#                                                     AdductStr <- cbind(AdductStr, "*", multi)
#                                                 }
#                                                 if (ppmcalc < 0) {
#                                                     ppmcalc <- - ppmcalc
#                                                 }
#                                                 ppmStr <- as.character(ppmcalc)
#                                                 cbind(listofdata, substr(mydata[k](0, headLen - 1), SGcount[k], PeakLabel[k], SaltLabel[k], useLabel[k], MFformula[i], AdductStr, ppmStr, elabel))
#                                                 matchMF[k] <- matchMF[k] + 1
#                                             }
#                                         }
#                                     }
#                                 }
#                             }
#                         }
#                     }
#                 }
#             }
#         }
#     }
#
#     for (k in 1 < nrow(mydata)) {
#         if (matchMF[k] < 1) {
#             myStr <- mydata[k]
#             checkdata <- strsplit(myStr, "\t")
#             temptab <- length(checkdata)
#             headLen <- length(mydata[k])
#             UseStr <- substr(mydata[k], 0, headLen - 1)
#             if (useLabel[k] == 5) {
#                 PeakLabel[k] <- "S"
#             }
#             cbind(listofdata, UseStr, SGcount[k], PeakLabel[k], SaltLabel[k], useLabel[k], elabel)
#         }
#     }
# }
# # else {
#     # cbind(listofdata, "Data not processed - ", tempStr1, tempStr3)
#     print(paste0("Data not processed - ", tempStr1, tempStr3))
# }
#
# outtxt <- listofdata
#
# print(listofdata)
