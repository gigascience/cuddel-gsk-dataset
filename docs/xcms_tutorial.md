# LC-MS data analysis using xcms

XCMS seems to be a very popular tool for this. Testing it from the R
command line using this [documentation](https://bioconductor.org/packages/devel/bioc/vignettes/xcms/inst/doc/xcmsPreprocess.pdf). 

There are five steps involved in processing LC-MS data:

1. Filter and identify peaks
2. Match peaks across samples
3. Correct retention time
4. Fill in missing peaks
5. Arrange peaks into a peak intensity table for statistical analysis

### Preparation of raw data files

Before the first initial step, the MS data has to be loaded into the 
xcms package. These data has to be in the form of NetCDF, mzXML or 
mzData files. Raw data in vendor-specific formats need to be converted
into one of these file formats.

The directory structure that the data files are stored in needs to be 
considered. xcms uses sample class information during pre-processing to 
decide which peaks are significant. If organised into sub-folders, the
samples will be automatically assigned to separate classes based on
their location. Samples may be separated into classes based on tissue 
type, mutation, gender, disease or time. For example, if you are 
analyzing the longitudinal effect of a drug in two patient groups, you 
may wish to put the groups into two directories “GroupA” and “GroupB”. 
Within each of these directories, you could further separate the 
samples by the time they were taken, such as “Day1”, “Day2”, etc. In 
xcms, they will be automatically assigned class names “GroupA/Day1”,
“GroupA/Day2”, etc.

The test data we used is the 12 spectra (NetCDF format) that come 
together with XCMS package. They are from 12 mice spinal cord samples 
collected by LC-MS (Saghatelian et al, PMID: 15533037). Group 1- 
wild-type (WT) or FAAH(+/+); group 2 – knock-out (KO) or FAAH (-/-). 
FAAH is the abbreviation of fatty acid amide hydrolase.

To load the test data:

```
> library(faahKO)
> cdfpath <- system.file("cdf", package = "faahKO")
> list.files(cdfpath, recursive = TRUE)
 [1] "KO/ko15.CDF" "KO/ko16.CDF" "KO/ko18.CDF" "KO/ko19.CDF" "KO/ko21.CDF"
 [6] "KO/ko22.CDF" "WT/wt15.CDF" "WT/wt16.CDF" "WT/wt18.CDF" "WT/wt19.CDF"
[11] "WT/wt21.CDF" "WT/wt22.CDF"
```

## Filtration and peak identification

The class object used for pre-processing analyte data from multiple 
LC/MS files is xcmsSet. It stores peak lists and provides methods for 
grouping and aligning those peaks. To create an xcmsSet object from a 
set of NetCDF files, use the xcmsSet() constructor function. It handles 
batch peak picking and generation of the xcmsSet object.

During peak identification, xcms uses a separate line for each sample 
to report on the status of processing. It outputs pairs of numbers 
separated by a colon. The first number is the m/z it is currently 
processing. The second number is the number of peaks that have been 
identified so far. It is important to note that the number may be 
significantly larger than the final number of peaks as a vicinity 
elimination post-processing step removes duplicate peaks corresponding 
to the same ion.

```
> library(xcms)
> cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)
> xset <- xcmsSet(cdffiles)
> xset
An "xcmsSet" object with 12 samples

Time range: 2506.1-4147.7 seconds (41.8-69.1 minutes)
Mass range: 200.1-599.3338 m/z
Peaks: 4721 (about 393 per sample)
Peak Groups: 0 
Sample classes: KO, WT 

Peak picking was performed on MS1.
Profile settings: method = bin
                  step = 0.1

Memory usage: 0.714 MB
```

The default parameters for xcmsSet should be acceptable but a number of
parameters may need to be optimised for a particular instrument or set
of samples. Information on these parameters can be found in the xcmsSet
function documentation and the findPeaks method. For example, findPeaks
can use different algorithms for peak detection and these are
findPeaks.matchedFilter, findPeaks.centWave, and 
findPeaks.centWaveWithPredictedIsotopeROIs.

## Matching peaks across samples

After peak identification, peaks representing the same analyte across 
samples must be placed into groups. This is accomplished with the group 
method, which returns a new xcmsSet object with the additional group 
information. The grouping process is non-destructive and does not 
affect other data stored in the xcmsSet object. The grouping 
algorithm processes the peak lists in order of increasing mass and will 
regularly output the mass it is currently working on.

```R
> xset <- group(xset)
262 325 387 450 512 575 
```

There are several grouping parameters to consider optimizing for your 
chromatography and mass spectrometer. Please consult the `group` 
documentation for more details. To see what the algorithm is doing while 
running, use the sleep argument to specify a time (in seconds) to pause 
and plot each iteration. That can be quite useful for visualizing 
parameter effects.

## Retention time correction

xcms can use these peak groups to identify and correct correlated drifts 
in retention time from run to run. The aligned peaks can then be used 
for a second pass of peak grouping which will be more accurate than the 
first. The whole process can be repeated in an iterative fashion, 
although we will only demonstrate a single pass of retention time 
alignment here.

Not all peak groups will be helpful for identifying retention time 
drifts. Some groups may be missing peaks from a large fraction of 
samples and thus provide an incomplete picture of the drift at that 
time point. Still others may contain multiple peaks from the same 
sample, which is a sign of improper grouping. xcms ignores these groups 
by only considering “well-behaved” peak groups which are missing at 
most one sample and have at most one extra peak. (Those values can be 
changed with the missing and extra arguments.)

For each of those well-behaved groups, the algorithm calculates a median 
retention time and, for every sample, a deviation from that median. 
Within a sample, the observed deviation generally changes over time in 
a nonlinear fashion. Those changes are approximated using a local 
polynomial regression technique implemented in the loess function. By 
default, the curve fitting is done using least-squares on all data 
points. However, it is possible to enable outlier detection and removal 
by setting the family argument to "symmetric", as shown here.

Retention time correction is performed by the retcor method, which 
returns an xcmsSet object with corrected retention times. Because it 
changes the retention times of all peaks, it is important to store the 
new object under a new variable name. That will allow you to backtrack 
and repeat retention time correction if necessary.

```
> xset2 <- retcor(xset, family = "symmetric", plottype = "mdevden")
Retention Time Correction Groups: 133 
```

The above command uses the plottype argument to produce a plot, shown in 
Figure 2, which is useful for supervising the algorithm. It includes the 
data points used for loess regression and the resulting deviation 
profiles. It additionally shows the distribution of peak groups across 
retention time.

After retention time correction, the initial peak grouping becomes 
invalid and is discarded. Therefore, the resulting object needs to be 
regrouped. Here, we decrease the inclusiveness of the grouping using the 
bw argument (default 30 seconds).

```R
> xset2 <- group(xset2, bw = 10)
262 325 387 450 512 575 
```

## Filling in missing peak data

After the second pass of peak grouping, there will still be peak groups 
which are missing peaks from some of the samples. That can occur because 
peaks were missed during peak identification or because an analyte was 
not present in a sample. In any case, those missing data points can be 
filled in by rereading the raw data files and integrating them in the 
regions of the missing peaks. This is performed using the fillPeaks 
method, which returns a xcmsSet object with the filled in peak data. 
While running, it outputs the name of the sample it is currently 
processing.

```R
> xset3 <- fillPeaks(xset2)
/opt/local/Library/Frameworks/R.framework/Versions/3.3/Resources/library/faahKO/cdf/KO/ko21.CDF
/opt/local/Library/Frameworks/R.framework/Versions/3.3/Resources/library/faahKO/cdf/KO/ko22.CDF
/opt/local/Library/Frameworks/R.framework/Versions/3.3/Resources/library/faahKO/cdf/WT/wt15.CDF
/opt/local/Library/Frameworks/R.framework/Versions/3.3/Resources/library/faahKO/cdf/WT/wt16.CDF
/opt/local/Library/Frameworks/R.framework/Versions/3.3/Resources/library/faahKO/cdf/WT/wt21.CDF
/opt/local/Library/Frameworks/R.framework/Versions/3.3/Resources/library/faahKO/cdf/WT/wt22.CDF
/opt/local/Library/Frameworks/R.framework/Versions/3.3/Resources/library/faahKO/cdf/KO/ko15.CDF
/opt/local/Library/Frameworks/R.framework/Versions/3.3/Resources/library/faahKO/cdf/KO/ko16.CDF
/opt/local/Library/Frameworks/R.framework/Versions/3.3/Resources/library/faahKO/cdf/WT/wt18.CDF
/opt/local/Library/Frameworks/R.framework/Versions/3.3/Resources/library/faahKO/cdf/WT/wt19.CDF
/opt/local/Library/Frameworks/R.framework/Versions/3.3/Resources/library/faahKO/cdf/KO/ko18.CDF
/opt/local/Library/Frameworks/R.framework/Versions/3.3/Resources/library/faahKO/cdf/KO/ko19.CDF
> xset3
An "xcmsSet" object with 12 samples

Time range: 2501.2-4148 seconds (41.7-69.1 minutes)
Mass range: 200.1-599.3338 m/z
Peaks: 6040 (about 503 per sample)
Peak Groups: 400
Sample classes: KO, WT

Peak picking was performed on MS1.
Profile settings: method = bin
                  step = 0.1

Memory usage: 1.74 MB
```

## Analysing and visualising results

A report showing the most statistically significant differences in 
analyte intensities can be generated with the diffreport method. It will 
automatically generate extracted ion chro- matograms for a given number 
of them, in this case 10. Several of those chromatograms are shown in 
Figure 3.

```
reporttab <- diffreport(xset3, "WT", "KO", "example", 10, metlin = 0.15, h=480, w=640)
ko15 Applying retention time correction to /opt/local/Library/Frameworks/R.framework/Versions/3.3/Resources/library/faahKO/cdf/KO/ko15.CDF
ko16 Applying retention time correction to /opt/local/Library/Frameworks/R.framework/Versions/3.3/Resources/library/faahKO/cdf/KO/ko16.CDF
ko18 Applying retention time correction to /opt/local/Library/Frameworks/R.framework/Versions/3.3/Resources/library/faahKO/cdf/KO/ko18.CDF
ko19 Applying retention time correction to /opt/local/Library/Frameworks/R.framework/Versions/3.3/Resources/library/faahKO/cdf/KO/ko19.CDF
ko21 Applying retention time correction to /opt/local/Library/Frameworks/R.framework/Versions/3.3/Resources/library/faahKO/cdf/KO/ko21.CDF
ko22 Applying retention time correction to /opt/local/Library/Frameworks/R.framework/Versions/3.3/Resources/library/faahKO/cdf/KO/ko22.CDF
wt15 Applying retention time correction to /opt/local/Library/Frameworks/R.framework/Versions/3.3/Resources/library/faahKO/cdf/WT/wt15.CDF
wt16 Applying retention time correction to /opt/local/Library/Frameworks/R.framework/Versions/3.3/Resources/library/faahKO/cdf/WT/wt16.CDF
wt18 Applying retention time correction to /opt/local/Library/Frameworks/R.framework/Versions/3.3/Resources/library/faahKO/cdf/WT/wt18.CDF
wt19 Applying retention time correction to /opt/local/Library/Frameworks/R.framework/Versions/3.3/Resources/library/faahKO/cdf/WT/wt19.CDF
wt21 Applying retention time correction to /opt/local/Library/Frameworks/R.framework/Versions/3.3/Resources/library/faahKO/cdf/WT/wt21.CDF
wt22 Applying retention time correction to /opt/local/Library/Frameworks/R.framework/Versions/3.3/Resources/library/faahKO/cdf/WT/wt22.CDF

> reporttab[1:4,]
       name      fold    tstat       pvalue    mzmed    mzmin    mzmax    rtmed
1 M300T3390  5.693594 14.44368 5.026336e-08 300.1898 300.1706 300.2000 3390.324
2 M301T3390  5.876588 15.57570 6.705719e-08 301.1879 301.1659 301.1949 3389.627
3 M298T3187  3.870918 11.93891 3.310025e-07 298.1508 298.1054 298.1592 3186.803
4 M491T3397 24.975703 16.83986 4.463361e-06 491.2000 491.1877 491.2063 3397.160
     rtmin    rtmax npeaks KO WT
1 3386.765 3396.335     12  6  6
2 3386.765 3392.101      7  6  1
3 3184.124 3191.312      4  4  0
4 3367.123 3424.681      6  6  0
                                                                     metlin
1 http://metlin.scripps.edu/metabo_list.php?mass_min=299.04&mass_max=299.34
2 http://metlin.scripps.edu/metabo_list.php?mass_min=300.04&mass_max=300.34
3     http://metlin.scripps.edu/metabo_list.php?mass_min=297&mass_max=297.3
4 http://metlin.scripps.edu/metabo_list.php?mass_min=490.05&mass_max=490.35
              ko15             ko16             ko18             ko19
1 4534353.62273683 4980914.48421051 5290739.13866664 4564262.89684209
2 962353.429578945 1047934.14136842 1109303.04472222 946943.392842103
3 180780.817277777 203926.952354905 191015.910842105  190626.84952381
4 432037.001363632     332159.07255     386966.75145 334951.452952381
              ko21             ko22             wt15             wt16
1 4733236.07999997      3931592.586  349660.88536842 491793.181333331
2 984787.204999993 806171.472899999 86450.4116463556 120096.519533855
3 156869.080488805      220288.6218 16269.0960107969  43677.783977158
4     294816.23565 373577.607619048 7643.13806685429 10519.9421688158
              wt18             wt19             wt21             wt22
1 645526.704947367 634108.848947367 1438254.44559999      1364627.844
2 143007.948300124 137319.686021111  218483.14259183 291392.971409092
3 54739.1288875568 76318.0076842747 54726.1153827847 49679.9424723104
4 26472.2930413739 33598.3227594968 8030.46734256427                0
>
```
If the metlin argument is set to a numeric value, the report will 
include links to the Metlin Metabolite Database 
(http://metlin.scripps.edu/) showing potential metabolite identities. A 
positive value indicates the data was acquired in positive ion mode and 
the neutral mass is calculated assuming all ions are M+H. A negative 
value does the opposite. The value itself indicates the uncertainty in 
mass accuracy. For instance, the first and third metabolites in the 
report produce the following URLs:

* http://metlin.scripps.edu/metabo_list.php?mass_min=299.04&mass_max=299.34   
* http://metlin.scripps.edu/metabo_list.php?mass_min=297&mass_max=297.3

### Selecting and visualising peaks

It is also possible to generate extracted ion chromatograms for 
arbitrary peak groups selected using various criteria. Here we generate 
EICs for two analytes eluting at different times. They are shown using 
both unaligned and aligned retention times. The resulting plots are 
shown in Figure 4.

```
> gt <- groups(xset3)
Error: unexpected '>' in ">"
> gt <- groups(xset3)
> colnames(gt)
[1] "mzmed"  "mzmin"  "mzmax"  "rtmed"  "rtmin"  "rtmax"  "npeaks" "KO"
[9] "WT"

> groupidx1 <- which(gt[,"rtmed"] > 2600 & gt[,"rtmed"] < 2700 & gt[,"npeaks"] == 12)
> groupidx2 <- which(gt[,"rtmed"] > 3600 & gt[,"rtmed"] < 3700 & gt[,"npeaks"] == 12)
> eiccor <- getEIC(xset3, groupidx = c(groupidx1, groupidx2))
> eicraw <- getEIC(xset3, groupidx = c(groupidx1, groupidx2), rt = "raw")
> plot(eicraw, xset3, groupidx = 1)
> plot(eicraw, xset3, groupidx = 2)
> plot(eiccor, xset3, groupidx = 1)
> plot(eiccor, xset3, groupidx = 2)
> cat("These are the warning")
These are the warning
> warnings()
NULL
```


## GSK dataset


Eva's GSK data set files are already converted into NetCDF format.

How to organise the GSK data set?




