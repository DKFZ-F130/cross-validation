Introduction
============

The following R-script was used to perform a performance evaluation of
frequently used online available MHC class-I prediction algorithms for
calculation and cross-validation of individual thresholds.

The performance evaluation is described in the article "Performance
evaluation of MHC class-I binding prediction tools based on an
experimentally validated MHC-peptide binding dataset" by Bonsack, Hoppe
et al. in Cancer Immunology Research (in Review).

From Material and Methods: "Calculation and validation of recommended
decision thresholds Based on our experimental binding affinity data we
calculated new threshold values for each analyzed prediction algorithm,
HLA type and peptide length. These recommended threshold values were
selected based on the following criteria: 1) specificity ≥0.66 (equal to
FPR ≤0.33), 2) TPR ≥2xFPR and 3) the threshold yielding the highest
possible sensitivity within the limits defined in 1) and 2). Lacking a
validation dataset (a second set of similarly obtained binding affinity
data) we used a bootstrapping algorithm to statistically validate our
recommended threshold values and their respective sensitivity,
specificity and accuracy. In brief, per HLA allele the dataset was
randomly split into 2/3 training data and 1/3 test data. Applying the
mentioned criteria, the optimal threshold for each prediction method was
calculated based on the training data. The calculated optimal threshold
was applied to the test data and sensitivity, specificity and accuracy
was calculated. This was repeated 100 times. From these 100 runs of
resampling, the median optimal threshold and the confidence intervals
for sensitivity, specificity and accuracy were calculated (data for
confidence intervals based on this bootstrapping are not shown). To
check the reliability of the validated threshold on arbitrary data sets
a second bootstrapping was performed. Again, our dataset was randomly
split into 2/3 training data and 1/3 test data in 100 runs of
resampling. The mean optimal threshold from the first bootstrapping as
well as the strong, intermediate and low binding affinity threshold
indicated by the methods were applied to the test data. In each run and
for each applied threshold, sensitivity, specificity and accuracy was
calculated. After 100 runs, the confidence intervals of sensitivity,
specificity and accuracy were calculated (shown in Figs. 4, 5, S4 and
S5). Differences of mean sensitivity, specificity and accuracy of
applied thresholds were compared by one-way ANOVA for repeated measures
followed by Dunnett's multiple comparisons test."

Here, we provide the R-scipt used to perform this statistical analysis,
the dataset per HLA-type and a description how to operate the script.

Cross Validation Example
========================

get cross-validation function from wd

    source("crossvalidate-data.R")

get data from csv file in wd

    a1_dat <- read.csv("A1.csv",sep=";",dec=".")

    head(a1_dat)

    ##   ba.result NetMHCpan.4.0 NetMHC.4.0 NetMHCpan.3.0 NetMHcpan.2.8
    ## 1    binder          56.3      32.89         107.9         81.20
    ## 2    binder         678.9    1142.27         413.9        267.28
    ## 3        nb        9369.8    8074.03        7907.1       6440.89
    ## 4    binder        1230.1     330.74        1377.0       2365.85
    ## 5        nb        5056.9    6457.31        4559.1       2286.75
    ## 6    binder        1323.5    3675.65        2245.9       2342.88
    ##   NetMHC.3.4 NetMHCcons.1.1 NetMHCcons.pickpocket.1.1 IEDB.recommended
    ## 1         55          66.92                   1145.59             0.25
    ## 2        250         257.38                   2145.70             0.45
    ## 3       5773        6095.65                  10641.85             0.70
    ## 4       2027        2192.64                   1967.79             0.55
    ## 5       5882        3665.80                   2443.19             0.95
    ## 6       1275        1728.19                   4288.48             0.95
    ##   IEDB.consensus IEDB.smmpmbec IEDB.smm MHC.flurry.1.2 ba.result.1
    ## 1           0.25        292.70   299.61        39.7664 0.898833333
    ## 2           0.35        384.18   314.63       312.5597 1.026666667
    ## 3           0.40        498.36   545.51     10190.3497          nb
    ## 4           0.45        764.79   866.56      1453.3475      15.656
    ## 5           0.65       1206.20  1001.27      2477.1102          nb
    ## 6           0.70       1041.21   941.46       193.5056        1.62
    ##     Sequence
    ## 1  ISEYRHYCY
    ## 2 HGDTPTLHEY
    ## 3 QPETTDLYCY
    ## 4 KISEYRHYCY
    ## 5  VCDKCLKFY
    ## 6 DLQPETTDLY

set number of cross-validations (min=100)

    b <- 100

give variable more a vector that contains individual threshold for each
predictor aka column

    more <- c(8641, 6706, 7097, 5684, 9154, 7693, 5211, 0.95, 0.7, 1042, 315)

maxgrid gives the maximum of the window to look for optimal threshold,
default is IC50=10000nM, consider increasing when optimal threshold is
expected to be higher than 10000nM

    val(b, a1_dat, more, maxgrid = 10000)
