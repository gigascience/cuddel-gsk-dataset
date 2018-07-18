# Title     : TODO
# Objective : TODO
# Created by: peterli
# Created on: 17/7/2018

library(MAIT)

# Load data for test analysis
library(faahKO)
cdfFiles <- system.file("cdf", package="faahKO", mustWork=TRUE)
MAIT <- sampleProcessing(dataDir = cdfFiles, project = "MAIT_Demo", snThres=2, rtStep=0.03)

summary(MAIT)
## A MAIT object built of 12 samples
## The object contains 6 samples of class KO
## The object contains 6 samples of class WT
##
## Parameters of the analysis:
## Value
## dataDir             "/usr/lib64/R/library/faahKO/cdf"
## snThres             "2"
## Sigma               "2.12332257516562"
## mzSlices            "0.3"
## retcorrMethod       "loess"
## groupMethod         "density"
## bwGroup             "3"
## mzWidGroup          "0.25"
## filterMethod        "centWave"
## rtStep              "0.03"
## nSlaves             "0"
## project             "MAIT_Demo"
## ppm                 "10"
## minfrac             "0.5"
## fwhm                "30"
## family1             "gaussian"
## family2             "symmetric"
## span                "0.2"
## centWave peakwidth1 "5"
## centWave peakwidth2 "20"

# Do peak annotation
MAIT <- peakAnnotation(MAIT.object = MAIT, corrWithSamp = 0.7, corrBetSamp = 0.75, perfwhm = 0.6)
## WARNING: No input adduct/fragment table was given. Selecting default MAIT table for positive polarity...
## Set adductTable equal to negAdducts to use the default MAIT table for negative polarity
## Start grouping after retention time.
## Created 321 pseudospectra.
## Spectrum build after retention time done
## Generating peak matrix!
## Run isotope peak annotation
## % finished: 10  20  30  40  50  60  70  80  90  100
## Found isotopes: 122
## Isotope annotation done
## Start grouping after correlation.
## Generating EIC's ..
##
## Calculating peak correlations in 321 Groups...
##  % finished: 10  20  30  40  50  60  70  80  90  100
##
## Calculating peak correlations across samples.
##  % finished: 10  20  30  40  50  60  70  80  90  100
##
## Calculating isotope assignments in 321 Groups...
##  % finished: 10  20  30  40  50  60  70  80  90  100
## Calculating graph cross linking in 321 Groups...
##  % finished: 10  20  30  40  50  60  70  80  90  100
## New number of ps-groups:  751
## xsAnnotate has now 751 groups, instead of 321
## Spectrum number increased after correlation done
## Generating peak matrix for peak annotation!
## Found and use user-defined ruleset!
## Calculating possible adducts in 751 Groups...
##  % finished: 10  20  30  40  50  60  70  80  90  100
## Adduct/fragment annotation done

# MAIT object has an xsAnnotated object containing the information about peaks, spectra and their annotation.
rawData(MAIT)
## $xsaFA
## An "xsAnnotate" object!
## With 751 groups (pseudospectra)
## With 12 samples and 1331 peaks
## Polarity mode is set to:  positive
## Using automatic sample selection
## Annotated isotopes: 122
## Annotated adducts & fragments: 81
## Memory usage: 4.45 MB

# Do statistical analysis to identify features different between classes using spectralSigFeatures function
MAIT <- spectralSigFeatures(MAIT.object = MAIT, pvalue = 0.05, p.adj = "none", scale = FALSE)

# Summarise analysis
summary(MAIT)
## A MAIT object built of 12 samples and 1331 peaks. No peak aggregation technique has been applied
## 63 of these peaks are statistically significant
## The object contains 6 samples of class KO
## The object contains 6 samples of class WT

## Parameters of the analysis:
## Value
## dataDir                              "/usr/lib64/R/library/faahKO/cdf"
## snThres                              "2"
## Sigma                                "2.12332257516562"
## mzSlices                             "0.3"
## retcorrMethod                        "loess"
## groupMethod                          "density"
## bwGroup                              "3"
## mzWidGroup                           "0.25"
## filterMethod                         "centWave"
## rtStep                               "0.03"
## nSlaves                              "0"
## project                              "MAIT_Demo"
## ppm                                  "10"
## minfrac                              "0.5"
## fwhm                                 "30"
## family1                              "gaussian"
## family2                              "symmetric"
## span                                 "0.2"
## centWave peakwidth1                  "5"
## centWave peakwidth2                  "20"
## corrWithSamp                         "0.7"
## corrBetSamp                          "0.75"
## perfwhm                              "0.6"
## sigma                                "6"
## peakAnnotation pvalue                "0.05"
## calcIso                              "TRUE"
## calcCiS                              "TRUE"
## calcCaS                              "TRUE"
## graphMethod                          "hcs"
## annotateAdducts                      "TRUE"
## peakAggregation method               "None"
## peakAggregation PCAscale             "FALSE"
## peakAggregation PCAcenter            "FALSE"
## peakAggregation scale                "FALSE"
## peakAggregation RemoveOnePeakSpectra "FALSE"
## Welch pvalue                         "0.05"
## Welch p.adj                          "none"

# In the spectralSigFeatures analysis, a table called signifcantFeatures.csv is created in the Tables subfolder. This
# table shows the characteristics of the statistically significant features. This table can be retrieved using:
signTable <- sigPeaksTable(MAIT.object = MAIT, printCSVfile = FALSE)

head(signTable)
## mz mzmin mzmax    rt rtmin rtmax npeaks KO WT        ko15       ko16       ko18       ko19       ko21       ko22         wt15        wt16
## 249 328.2 328.1 328.2 56.31 56.27 56.45      4  4  0    43851.29   88615.33   41311.31   35218.11   40095.58   47006.34     1907.784    1165.861
## 884 496.2 496.2 496.2 56.27 56.12 56.44      7  3  3 11275649.72 3795994.86 2624223.82 3630452.88 8335183.79 5624245.57 36141998.610 3377994.510
## 891 497.2 497.2 497.2 56.23 56.12 56.41      6  3  3  9239784.03 2603425.00  798816.25 1219126.64 2340526.64 1560252.53  9219730.435 1993433.398
## 896 498.2 498.2 498.2 56.21 56.06 56.23      4  1  3   853625.83   42028.88  150286.13  217955.75  198805.08  356467.11  1837432.274   96019.739
## 899 499.2 499.1 499.2 56.19 55.93 56.34      8  3  4    86962.36   13390.91   53202.17   17416.74   36593.84   54465.13   124631.905   12083.597
## 953 508.2 508.1 508.2 56.19 56.15 56.20      5  4  1   164074.60  139084.68  271571.75  155172.88       0.00  199981.96    99744.641  196969.877
## wt18        wt19        wt21        wt22   isotopes adduct pcgroup    P.adjust           p Fisher.Test Mean Class KO Mean Class WT
## 249   12587.95    7536.961    6619.609    6892.127                         4 0.002514234 0.002514234          NA      49349.66      6118.382
## 884 2885261.59 5243788.829 7706061.535 7808922.868   [62][M]+              4 0.421686461 0.421686461          NA    5880958.44  10527337.990
## 891  978322.88 1472473.239 2221906.241 2241262.797 [62][M+1]+              4 0.973650701 0.973650701          NA    2960321.85   3021188.165
## 896  202991.67  243572.435  628464.047  357991.803 [62][M+2]+              4 0.404976403 0.404976403          NA     303194.80    561078.661
## 899   14849.16   94691.930   41735.420   42589.910 [62][M+3]+              4 0.610170147 0.610170147          NA      43671.86     55096.988
## 953  153144.41  245829.063  496743.520   31504.200                         4 0.534531870 0.534531870          NA     154980.98    203989.286
## Median Class KO Median Class WT
## 249        42581.30        6755.868
## 884      4710120.21     6474925.182
## 891      1950389.59     2107669.820
## 896       208380.41      300782.119
## 899        44898.01       42162.665
## 953       159623.74      175057.144

# Summarise results
MAIT
## A MAIT object built of 12 samples and 1331 peaks. No peak aggregation technique has been applied
## 63 of these peaks are statistically significant
## The object contains 6 samples of class KO
## The object contains 6 samples of class WT

# Do statistical plots. MAIT objects are created wutg PCA and PLS models saved inside them
plotBoxplot(MAIT)
plotHeatmap(MAIT)

# The plot output figures are saved in the subfolders in the project folder
MAIT <- plotPCA(MAIT, plot3d=FALSE)
MAIT <- plotPLS(MAIT, plot3d=FALSE)

PLSmodel <- model(MAIT, type = "PLS")
PCAmodel <- model(MAIT, type = "PCA")

PLSmodel
## Partial Least Squares
##
## 12 samples
## 63 predictors
## 2 classes: 'KO', 'WT'
##
## No pre-processing
## Resampling: Bootstrapped (25 reps)
## Summary of sample sizes: 12, 12, 12, 12, 12, 12, ...
## Resampling results across tuning parameters:
##
## ncomp  Accuracy  Kappa
## 1      1         1
## 2      1         1
## 3      1         1
##
## Accuracy was used to select the optimal model using the largest value.
## The final value used for the model was ncomp = 1.

pcaScores(MAIT)
## PC1         PC2        PC3        PC4         PC5        PC6          PC7         PC8         PC9        PC10        PC11          PC12
## [1,] -8.823400  6.31998966 -0.5526132 -1.2916266 -1.39191835  0.6247246 -0.186452381  0.17513153 -0.03013042  0.07214817 -0.13531712  2.466777e-15
## [2,] -8.202057  0.08372414  2.0231200  3.4017405  1.26354572 -0.4580168  1.896542591 -0.57535747  0.71556592  0.27847053  0.37346246  6.869505e-16
## [3,] -6.647327 -3.00939301  3.6390245 -2.6080710 -0.47162186 -0.9921232 -1.539667265 -0.36479469  0.51448150 -0.11756634 -0.19180616 -6.678685e-16
## [4,] -3.554215 -2.11458548  0.2277483 -0.4454514  0.93503858  1.8743863  1.000764293  2.36685217 -1.32658830 -0.45467997 -0.15334779 -8.083811e-16
## [5,] -4.394232 -3.25707904 -4.0828290 -0.5586029 -0.59580779  1.7252921  0.037541531 -1.35815385  0.59885276  0.56612633  0.16379772 -9.783840e-16
## [6,] -3.165613 -0.45844726 -2.8491469  1.4610838  0.38569573 -2.9636763 -1.383997315  0.16386948 -1.20232530 -0.45472153 -0.16106150 -4.440892e-16
## [7,]  5.305730  1.64665278 -0.9821137 -1.8563330  3.87308588 -0.1983468 -0.007830579 -0.08448479  1.11571569 -0.43147651 -0.09692091 -1.169204e-15
## [8,]  5.487610  0.31429828  0.9351818  2.6945683 -0.73282854  1.7723295 -1.748001788 -0.43156894  0.44464047 -0.96597551 -0.92322133 -3.400058e-16
## [9,]  6.313361  0.39820865  1.2046214 -0.9395212 -0.02530098 -0.1413301  1.580728764 -1.87584484 -1.75244866  0.53327213 -0.66087917  7.109114e-16
## [10,]  5.641545  0.73284993  1.0706125  0.4904065  0.21710858  0.7084426 -1.784234044  0.09591630 -0.68029637  0.77764390  1.64880640 -6.730727e-16
## [11,]  6.093996 -0.28217905 -0.2864657  0.3521655 -1.13341090 -1.0305014  0.315105967  1.71724384  0.99906322  1.59307934 -0.67748775  1.377370e-15
## [12,]  5.944601 -0.37403961 -0.3471400 -0.7003586 -2.32358607 -0.9211803  1.819500227  0.17119127  0.60346949 -1.39632054  0.81397514 -3.538836e-16

# Before identifying metabolites, peak annotation can be improved using the Biotransformations function. Use the default
# MAIT table for biotransformations:
Biotransformations(MAIT.object = MAIT, peakPrecision = 0.005)

# A user-defined biotransformations table can be used - see MAIT manual

# Metabolite identification
MAIT <- identifyMetabolites(MAIT.object = MAIT, peakTolerance = 0.005)
## WARNING: No input database table was given. Selecting default MAIT database...
## Metabolite identification initiated
##
## % Metabolite identification in progress: 10  20  30  40  50  60  70  80  90  100
## Metabolite identification finished

# A table is created containing the possible metabolite identifications
metTable <- metaboliteTable(MAIT)
# View some results in the table
metTable[1:5,1:ncol(metTable)]
## Query Mass Database Mass (neutral mass)    rt    Isotope Adduct    Name spectra Biofluid   ENTRY        p.adj            p Fisher.Test
## 1      328.2                      Unknown 56.31                   Unknown       4  unknown unknown 2.514234e-03 2.514234e-03          NA
## 2      454.1                      Unknown 55.88   [45][M]+        Unknown       7  unknown unknown 6.124024e-01 6.124024e-01          NA
## 3      549.1                      Unknown  53.1 [91][M+1]+        Unknown     215  unknown unknown 6.646044e-01 6.646044e-01          NA
## 4      411.2                      Unknown  65.6                   Unknown     221  unknown unknown 6.761483e-07 6.761483e-07          NA
## 5      324.2                      Unknown 54.59                   Unknown     233  unknown unknown 3.990481e-02 3.990481e-02          NA
## Mean Class KO Mean Class WT Median Class KO Median Class WT KO WT      ko15      ko16       ko18      ko19      ko21      ko22       wt15
## 1      49349.66      6118.382        42581.30        6755.868  4  0  43851.29  88615.33   41311.31  35218.11  40095.58  47006.34   1907.784
## 2     710915.10    882177.176       595445.77      769797.473  6  4 381784.88 878240.29 1511025.72 797661.72 303548.19 393229.82 377503.040
## 3      59572.67     55136.583        57388.44       54346.972  2  5  61471.41  79003.52   85458.18  52338.98  25858.49  53305.46  73945.095
## 4     148772.97     33519.391       145294.01       23729.675  6  4 172806.65 133447.55  167487.24 128308.37 141713.88 148874.14  47465.487
## 5     109873.48      9594.587        84825.35        7569.202  5  0  96901.74 287813.15   79667.89  55325.88  49549.44  89982.80  19287.067
## wt16        wt18        wt19       wt21        wt22
## 1   1165.861   12587.948    7536.961   6619.609    6892.127
## 2 492109.145 1419864.480 1815052.071 141048.523 1047485.800
## 3  51594.920   52613.735   39506.860  56080.210   57078.680
## 4  64596.940   22508.083   24951.268  20562.535   21032.035
## 5      0.000    8493.043    6645.361  20716.752    2425.299

# Validation to check predictive value of significant features
MAIT <- Validation(Iterations = 20, trainSamples= 3, MAIT.object = MAIT)
summary(MAIT)

# External peak data can be analysed by MAIT using hte MAITbuilder function to import peak data and analyse it using
# MAIT statistical functions. Consider that we have the following data:
peaks <- scores(MAIT)
masses <- getPeaklist(MAIT)$mz
rt <- getPeaklist(MAIT)$rt/60

# Use the MAITbuilder to annotate and identify metabolites on these data - check parameter values are suitable!!!
importMAIT <- MAITbuilder(data=peaks, masses=masses, rt=rt, significantFeatures=TRUE, spectraEstimation=TRUE, rtRange=0.2, corThresh=0.7)

# Perform biotransformations - to run negative annotation, sert adductTable = negAdducts
importMAIT <- Biotransformations(MAIT.object = importMAIT, adductAnnotation = TRUE, peakPrecision = 0.005, adductTable = NULL)

# Identify metabolites
importMAIT <- identifyMetabolites(MAIT.object = importMAIT, peakTolerance=0.005, polarity="positive")
## WARNING: No input database table was given. Selecting default MAIT database...
## Metabolite identification initiated
##
## % Metabolite identification in progress: 0  10  20  30  40  50  60  70  80  90  100
## Metabolite identification finished
