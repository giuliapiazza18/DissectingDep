Dissecting Depression Symptomatology: Brain Phenotypes as Mediators
Between Polygenic Scores and Individual Symptoms in the UK Biobank
================
Giulia Piazza
2025-01-16

# Data cleaning

## Load data on the RAP

``` r
system('dx download "/output/pgs/DEP_05June2024/DEP_pred_auto.txt"')
system('dx download "/output/pgs/ADHD_17June2024/ADHD_pred_auto.txt"')
system('dx download "/output/pgs/SCZ_26July2024/SCZ_pred_auto.txt"')
system('dx download "/output/pgs/BD_26July2024/BD_pred_auto.txt"')
system('dx download "/data/brain_data.csv"')
system('dx download "/data/demo_pheno_data.csv"')
system("dx download /scripts/mediation/functions/make_multiple_model_cov.R")
system("dx download /scripts/mediation/functions/single_model.R")
system("dx download /scripts/mediation/functions/adjust_multiple_df.R")
system("dx download /scripts/mediation/functions/create_single_label.R")
system("dx download /scripts/mediation/functions/split_data.R")
system("dx download /scripts/mediation/functions/create_path_label.R")
system("dx download /scripts/mediation/functions/create_path_label_tot.R")
system("dx download /output/merged/qc/imputed.allchr.hapmap.subset.QC.pruned.pca.eigenvec")

source("make_multiple_model_cov.R") # function to create model
source("adjust_multiple_df.R") # function to adjust data frame for plots and multiple comparisons
source("single_model.R")
source("split_data.R")
source("create_path_label.R")
source("create_single_label.R")
source("create_path_label_tot.R")
```

## Install packages and load

``` r
install.packages(c("lavaan", "psych", "corrplot", "corpcor", "mice", "naniar", "miceadds"), repos = "http://cran.us.r-project.org")

remotes::install_github("TDJorgensen/lavaan.mi")
devtools::install_github("slowkow/ggrepel")
```

``` r
library(dplyr)
library(lavaan)
library(psych)
library(corpcor)
library(ggplot2)
library(corrplot)
library(naniar)
library(mice)
library(lavaan.mi)
library(ggrepel)
library(miceadds)
library(forcats)
```

## Load data in R

``` r
dep_pgs = read.table("DEP_pred_auto.txt", header = TRUE) %>%
  dplyr::select(sample.ID, final_pred_auto) %>% rename(., dep.pgs = final_pred_auto)
adhd_pgs = read.table("ADHD_pred_auto.txt", header = TRUE) %>%
  dplyr::select(sample.ID, final_pred_auto)  %>% rename(., adhd.pgs = final_pred_auto)
scz_pgs = read.table("SCZ_pred_auto.txt", header = TRUE) %>%
  dplyr::select(sample.ID, final_pred_auto)  %>% rename(., scz.pgs = final_pred_auto)
bd_pgs = read.table("BD_pred_auto.txt", header = TRUE) %>%
  dplyr::select(sample.ID, final_pred_auto)  %>% rename(., bd.pgs = final_pred_auto)

pheno_data = read.csv("demo_pheno_data.csv")
brain_data = read.csv("brain_data.csv")
```

## Keep individuals with polygenic scores

``` r
data.pgs = merge(dep_pgs, pheno_data, by.x = "sample.ID", by.y = "Participant.ID", all.x = TRUE) %>%
  merge(., brain_data, by.x = "sample.ID", by.y = "Participant.ID", all.x = T) %>% 
  merge(., adhd_pgs, by = "sample.ID", all.x = T) %>%
  merge(., scz_pgs, by = "sample.ID", all.x = T) %>%
  merge(., bd_pgs, by = "sample.ID", all.x = T)
```

## Select and rename variables for mediation

``` r
data.med = data.pgs %>%
  rename(
    id = sample.ID,
    age = Age.at.recruitment,
    dep.pgs = dep.pgs,
    adhd.pgs = adhd.pgs,
    scz.pgs = scz.pgs,
    bd.pgs = bd.pgs,
    tow = Townsend.deprivation.index.at.recruitment,
    sex = Sex,
    yob = Year.of.birth,
    mob = Month.of.birth,
    emp = Current.employment.status...Instance.2,
    qual = Qualifications...Instance.2,
    eth = Ethnic.background...Instance.2,
    anx = Recent.feelings.or.nervousness.or.anxiety,
    wor = Recent.inability.to.stop.or.control.worrying,
    psy = Recent.changes.in.speed.amount.of.moving.or.speaking,
    dep = Recent.feelings.of.depression,
    ina = Recent.feelings.of.inadequacy,
    tir = Recent.feelings.of.tiredness.or.low.energy,
    int = Recent.lack.of.interest.or.pleasure.in.doing.things,
    app = Recent.poor.appetite.or.overeating,
    sui = Recent.thoughts.of.suicide.or.self.harm,
    con = Recent.trouble.concentrating.on.things,
    sle = Trouble.falling.or.staying.asleep..or.sleeping.too.much,
    irr = Recent.easy.annoyance.or.irritability,
    fore = Recent.feelings.of.foreboding,
    res = Recent.restlessness,
    rel = Recent.trouble.relaxing,
    wor.t = Recent.worrying.too.much.about.different.things,
    mofc.left = Volume.of.medialorbitofrontal..left.hemisphere....Instance.2,
    mofc.right = Volume.of.medialorbitofrontal..right.hemisphere....Instance.2,
    fusi.left = Volume.of.fusiform..left.hemisphere....Instance.2,
    fusi.right = Volume.of.fusiform..right.hemisphere....Instance.2,
    ins.left = Volume.of.insula..left.hemisphere....Instance.2,
    ins.right = Volume.of.insula..right.hemisphere....Instance.2,
    ra.cing.right = Volume.of.rostralanteriorcingulate..right.hemisphere....Instance.2,
    ra.cing.left = Volume.of.rostralanteriorcingulate..left.hemisphere....Instance.2,
    p.cing.right = Volume.of.posteriorcingulate..right.hemisphere....Instance.2,
    p.cing.left = Volume.of.posteriorcingulate..left.hemisphere....Instance.2,
    ca.cing.right = Volume.of.caudalanteriorcingulate..right.hemisphere....Instance.2,
    ca.cing.left = Volume.of.caudalanteriorcingulate..left.hemisphere....Instance.2,
    hip.right = Volume.of.Hippocampus..right.hemisphere....Instance.2,
    hip.left = Volume.of.Hippocampus..left.hemisphere....Instance.2,
    tot.vol=Volume.of.brain..grey.white.matter..normalised.for.head.size....Instance.2) %>%
  dplyr::select(id, age, sex, yob, mob, emp, tow, qual, eth, dep.pgs, adhd.pgs,
         scz.pgs, bd.pgs, anx, wor, psy:wor.t, mofc.left:hip.right)
```

## Keep the cases with neuroimaging data

``` r
data.med %>% dplyr::select(mofc.left:hip.left) %>% gg_miss_upset(nsets = 50) # looks like if they are missing one IDP, they are missing all IDPs

data.med = data.med[complete.cases(data.med$mofc.left),]
```

## Recode data

``` r
# data coding info is here https://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=168

data.reco = data.med %>%
  mutate_at(vars(anx:wor.t), ~recode(., `0` = "0", `1`= "1", `2`= "2", `3`= "3", .default = NA_character_)) %>%
  mutate_at(vars(anx:wor.t), as.numeric) # -3 (prefer not to answer) is now NA
```

## Combine Imaging Derived Phenotypes across hemispheres

``` r
data.reco.brain = data.reco %>%
  rowwise() %>%
  mutate(mofc = mean(c(mofc.left, mofc.right))) %>%
  mutate(fusi = mean(c(fusi.right, fusi.left))) %>%
  mutate(ins = mean(c(ins.right, ins.left))) %>%
  mutate(hip = mean(c(hip.right, hip.left))) %>%
  mutate(cing = mean(c(ca.cing.right, ca.cing.left, ra.cing.left, ra.cing.right)))
```

## Scale brain data and PGS

``` r
# (need to learn a better way to do this)
data.scale = data.reco.brain
data.scale$mofc = scale(data.scale$mofc)[,1]
data.scale$hip = scale(data.scale$hip)[,1]
data.scale$fusi = scale(data.scale$fusi)[,1]
data.scale$ins = scale(data.scale$ins)[,1]
data.scale$cing = scale(data.scale$cing)[,1]
data.scale$dep.pgs = scale(data.scale$dep.pgs)[,1]
data.scale$adhd.pgs = scale(data.scale$adhd.pgs)[,1]
data.scale$scz.pgs = scale(data.scale$scz.pgs)[,1]
data.scale$bd.pgs = scale(data.scale$bd.pgs)[,1]
```

## Add genetic principal components

``` r
princ.comp = read.table("imputed.allchr.hapmap.subset.QC.pruned.pca.eigenvec")
names(princ.comp) = c("IID", "FAM", paste0("PC", c(1:10)))

# Add updated ethnicity and check we're only including Caucasian participants
# data.cov = merge(data.scale, updated.eth, by = "id", all.x = TRUE) # all Caucasian


# Add 10 principal components
data.scale.pc = merge(data.scale, princ.comp, by.x = "id", by.y = "IID")
```

## Describe dataset

``` r
sum.data.scale = data.scale.pc %>%
  rowwise() %>%
  mutate(tot.anx = sum(c(anx, wor, wor.t, rel, res, irr, fore)),
         tot.dep = sum(c(psy, dep, ina, tir, int, app, sui, con, sle))) %>%
  psych::describe()
sum.data.scale.vis = as.data.frame(lapply(sum.data.scale, round, digits = 2))
rownames(sum.data.scale.vis) = rownames(sum.data.scale)
sum.data.scale.vis
```

    ##               vars     n       mean         sd     median    trimmed        mad
    ## id               1 17823 3520596.11 1448625.35 3513615.00 3521024.53 1854645.13
    ## age              2 17823      55.17       7.45      56.00      55.32       8.90
    ## sex              3 17823       0.48       0.50       0.00       0.48       0.00
    ## yob              4 17823    1952.84       7.44    1952.00    1952.69       8.90
    ## mob              5 17823       6.43       3.43       6.00       6.41       4.45
    ## emp              6 17823      15.09       9.15      20.00      14.79       0.00
    ## tow              7 17805      -2.07       2.62      -2.73      -2.38       2.12
    ## qual             8 17823      28.92      22.88      27.00      27.89      34.10
    ## eth              9  4186    1000.29      26.84    1001.00    1001.00       0.00
    ## dep.pgs         10 17823       0.00       1.00       0.00       0.00       1.01
    ## adhd.pgs        11 17823       0.00       1.00       0.01       0.00       1.00
    ## scz.pgs         12 17823       0.00       1.00       0.00       0.00       1.00
    ## bd.pgs          13 17823       0.00       1.00       0.01       0.00       1.01
    ## anx             14 13616       0.32       0.60       0.00       0.20       0.00
    ## wor             15 13611       0.28       0.60       0.00       0.16       0.00
    ## psy             16 13725       0.06       0.32       0.00       0.00       0.00
    ## dep             17 13699       0.26       0.55       0.00       0.14       0.00
    ## ina             18 13675       0.22       0.56       0.00       0.09       0.00
    ## tir             19 13711       0.62       0.81       0.00       0.46       0.00
    ## int             20 13710       0.24       0.56       0.00       0.11       0.00
    ## app             21 13724       0.23       0.60       0.00       0.08       0.00
    ## sui             22 13647       0.04       0.26       0.00       0.00       0.00
    ## con             23 13723       0.23       0.55       0.00       0.10       0.00
    ## sle             24 13707       0.69       0.90       0.00       0.52       0.00
    ## irr             25 13596       0.28       0.55       0.00       0.18       0.00
    ## fore            26 13597       0.21       0.53       0.00       0.08       0.00
    ## res             27 13623       0.12       0.42       0.00       0.00       0.00
    ## rel             28 13616       0.31       0.63       0.00       0.18       0.00
    ## wor.t           29 13606       0.37       0.64       0.00       0.25       0.00
    ## mofc.left       30 17823    5022.88     613.94    4997.00    5008.65     607.87
    ## mofc.right      31 17823    4892.32     578.35    4875.00    4880.38     575.25
    ## fusi.left       32 17823    9109.58    1238.33    9033.00    9067.22    1211.28
    ## fusi.right      33 17823    8915.44    1268.52    8840.00    8871.05    1264.66
    ## ins.left        34 17823    6417.51     692.87    6375.00    6394.92     689.41
    ## ins.right       35 17823    6630.61     715.92    6586.00    6606.19     708.68
    ## ra.cing.right   36 17823    2674.30     530.28    2643.00    2657.39     521.88
    ## ra.cing.left    37 17823    3751.48     642.19    3714.00    3729.66     634.55
    ## p.cing.right    38 17823    3607.48     557.02    3582.00    3595.65     536.70
    ## p.cing.left     39 17823    3598.57     539.73    3569.00    3583.88     524.84
    ## ca.cing.right   40 17823    2230.61     571.49    2220.00    2220.01     538.18
    ## ca.cing.left    41 17823    3105.97     583.94    3069.00    3085.26     561.91
    ## hip.left        42 17823    3942.48     408.83    3930.10    3935.50     401.04
    ## hip.right       43 17823    4111.71     430.92    4101.30    4104.55     416.91
    ## mofc            44 17823       0.00       1.00      -0.04      -0.02       0.99
    ## fusi            45 17823       0.00       1.00      -0.05      -0.03       1.00
    ## ins             46 17823       0.00       1.00      -0.06      -0.03       0.99
    ## hip             47 17823       0.00       1.00      -0.02      -0.02       0.98
    ## cing            48 17823       0.00       1.00      -0.07      -0.04       0.97
    ## FAM             49 17823 3520596.11 1448625.35 3513615.00 3521024.53 1854645.13
    ## PC1             50 17823       0.00       0.01       0.00       0.00       0.01
    ## PC2             51 17823       0.00       0.01       0.00       0.00       0.01
    ## PC3             52 17823       0.00       0.01       0.00       0.00       0.00
    ## PC4             53 17823       0.00       0.01       0.00       0.00       0.01
    ## PC5             54 17823       0.00       0.01       0.00       0.00       0.01
    ## PC6             55 17823       0.00       0.01       0.00       0.00       0.01
    ## PC7             56 17823       0.00       0.01       0.00       0.00       0.00
    ## PC8             57 17823       0.00       0.01       0.00       0.00       0.01
    ## PC9             58 17823       0.00       0.01       0.00       0.00       0.00
    ## PC10            59 17823       0.00       0.01       0.00       0.00       0.01
    ## tot.anx         60 13508       1.88       3.12       0.00       1.19       0.00
    ## tot.dep         61 13520       2.54       3.52       1.00       1.81       1.48
    ##                      min        max      range   skew kurtosis       se
    ## id            1000166.00 6025840.00 5025674.00   0.00    -1.19 10850.90
    ## age                40.00      70.00      30.00  -0.16    -0.95     0.06
    ## sex                 0.00       1.00       1.00   0.07    -2.00     0.00
    ## yob              1937.00    1970.00      33.00   0.16    -0.92     0.06
    ## mob                 1.00      12.00      11.00   0.04    -1.20     0.03
    ## emp                 1.00      50.00      49.00   0.27     0.28     0.07
    ## tow                -6.26       9.23      15.49   1.12     1.01     0.02
    ## qual                1.00      64.00      63.00   0.19    -1.62     0.17
    ## eth                -3.00    1003.00    1006.00 -37.30  1389.66     0.41
    ## dep.pgs            -4.58       4.44       9.03   0.01     0.00     0.01
    ## adhd.pgs           -4.22       4.17       8.39   0.00     0.02     0.01
    ## scz.pgs            -4.00       4.19       8.20  -0.01     0.00     0.01
    ## bd.pgs             -4.55       4.26       8.81  -0.03     0.01     0.01
    ## anx                 0.00       3.00       3.00   2.21     5.52     0.01
    ## wor                 0.00       3.00       3.00   2.50     6.88     0.01
    ## psy                 0.00       3.00       3.00   6.26    44.55     0.00
    ## dep                 0.00       3.00       3.00   2.47     6.92     0.00
    ## ina                 0.00       3.00       3.00   2.98     9.66     0.00
    ## tir                 0.00       3.00       3.00   1.40     1.57     0.01
    ## int                 0.00       3.00       3.00   2.74     8.18     0.00
    ## app                 0.00       3.00       3.00   3.01     9.26     0.01
    ## sui                 0.00       3.00       3.00   7.07    59.44     0.00
    ## con                 0.00       3.00       3.00   2.85     9.05     0.00
    ## sle                 0.00       3.00       3.00   1.27     0.78     0.01
    ## irr                 0.00       3.00       3.00   2.22     5.95     0.00
    ## fore                0.00       3.00       3.00   3.04    10.39     0.00
    ## res                 0.00       3.00       3.00   4.12    19.49     0.00
    ## rel                 0.00       3.00       3.00   2.32     5.74     0.01
    ## wor.t               0.00       3.00       3.00   1.98     4.38     0.01
    ## mofc.left        2219.00    7616.00    5397.00   0.24     0.07     4.60
    ## mofc.right       2800.00    7401.00    4601.00   0.22     0.13     4.33
    ## fusi.left        4063.00   15964.00   11901.00   0.38     0.37     9.28
    ## fusi.right       4758.00   14997.00   10239.00   0.36     0.15     9.50
    ## ins.left         3389.00    9389.00    6000.00   0.33     0.14     5.19
    ## ins.right        3747.00    9855.00    6108.00   0.35     0.19     5.36
    ## ra.cing.right     430.00    5429.00    4999.00   0.33     0.21     3.97
    ## ra.cing.left     1369.00    6633.00    5264.00   0.36     0.25     4.81
    ## p.cing.right     1271.00    6503.00    5232.00   0.22     0.43     4.17
    ## p.cing.left      1168.00    6433.00    5265.00   0.29     0.41     4.04
    ## ca.cing.right     431.00    5754.00    5323.00   0.29     0.70     4.28
    ## ca.cing.left      867.00   10707.00    9840.00   0.51     1.98     4.37
    ## hip.left         1625.10    6498.10    4873.00   0.17     0.43     3.06
    ## hip.right        2196.30    6106.00    3909.70   0.16     0.35     3.23
    ## mofc               -3.61       3.89       7.51   0.22     0.13     0.01
    ## fusi               -3.06       5.61       8.67   0.33     0.15     0.01
    ## ins                -3.57       4.21       7.77   0.35     0.14     0.01
    ## hip                -4.82       4.86       9.68   0.17     0.28     0.01
    ## cing               -3.55       6.50      10.05   0.41     0.37     0.01
    ## FAM           1000166.00 6025840.00 5025674.00   0.00    -1.19 10850.90
    ## PC1                -0.03       0.01       0.04  -0.60    -0.03     0.00
    ## PC2                -0.02       0.01       0.03  -0.14    -0.79     0.00
    ## PC3                -0.02       0.03       0.04   0.91     0.74     0.00
    ## PC4                -0.03       0.02       0.05  -0.35     0.41     0.00
    ## PC5                -0.03       0.02       0.05   0.01     0.28     0.00
    ## PC6                -0.02       0.03       0.04   0.30     0.05     0.00
    ## PC7                -0.03       0.05       0.08   1.10     5.37     0.00
    ## PC8                -0.02       0.03       0.05   0.29     0.21     0.00
    ## PC9                -0.05       0.03       0.07  -0.67     3.32     0.00
    ## PC10               -0.04       0.02       0.06  -1.06     2.99     0.00
    ## tot.anx             0.00      21.00      21.00   2.44     7.34     0.03
    ## tot.dep             0.00      27.00      27.00   2.41     7.67     0.03

``` r
table(data.scale.pc$sex) # 0 = female, 1 = male
```

    ## 
    ##    0    1 
    ## 9217 8606

## Save dataset

``` r
write.csv(data.scale.pc,"data.scale.pc.csv", quote = F, row.names = F)
```

## Explore correlations

``` r
data.cor = data.scale.pc %>% dplyr::select(id, dep.pgs:wor.t, mofc:cing, sex, age, "PC1":"PC10")

cor_all = cor(data.cor[-1], use = "p")
pcor_all = cor2pcor(cor_all)
colnames(pcor_all) = colnames(cor_all)
rownames(pcor_all) = rownames(cor_all)
```

### Partial correlations (everything)

``` r
pcor_plot = corrplot(pcor_all, col = COL2('PuOr'), type = "upper", method = "shade", diag = FALSE)
```

![](mediation_analysis_ukb17122024-114527_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

## Explore missing data

``` r
# percent missing in variables, per case etc 
miss.var = miss_var_summary(data.cor)
cat("The maximum percentage of missing values across variables is",
  max(miss.var$pct_miss)
  )
```

    ## The maximum percentage of missing values across variables is 23.71655

``` r
miss.case = miss_case_table(data.cor)
cat("The maximum percentage of missing values across cases is",
  miss.case %>% filter(n_miss_in_case > 0) %>% dplyr::select(pct_cases) %>% max()
)
```

    ## The maximum percentage of missing values across cases is 22.91421

``` r
data.cor %>% dplyr::select(anx:rel) %>% gg_miss_var()
```

![](mediation_analysis_ukb17122024-114527_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

## Multiple imputation

### Prepare data frame for model (function needs specific names)

``` r
NO = data.cor %>% dplyr::select(anx:wor.t) %>% colnames() %>% length() # number of outcomes
NM = data.cor %>% dplyr::select(mofc:cing) %>% colnames() %>% length() # number of mediators
NE = data.cor %>% dplyr::select(dep.pgs:bd.pgs) %>% colnames() %>% length() # number of exposures
NCV = data.cor %>% dplyr::select(age, sex, "PC1":"PC10") %>% colnames() %>% length() # number of covariates

model_names = c("id",
                paste0("E", 1:NE),
                paste0("O", 1:NO),
                paste0("M", 1:NM),
                paste0("CV", 1:NCV))

names(model_names) = names(data.cor)

data.lav = data.cor
names(data.lav) = model_names 
```

### Impute data

``` r
data.imputed = mice(data.lav,
                    m = 25, 
                    method = "pmm",
                    printFlag = FALSE)
```

# Mediation analysis

## Create multiple PGS model

``` r
mediation.model = make_multiple_model_cov(NO, NE, NM, NCV)
mediation.model
```

    ##   [1] "O1 ~ b1.1*M1 + b2.1*M2 + b3.1*M3 + b4.1*M4 + b5.1*M5 + c1.1*E1 + c2.1*E2 + c3.1*E3 + c4.1*E4 + CV1 + CV2 + CV3 + CV4 + CV5 + CV6 + CV7 + CV8 + CV9 + CV10 + CV11 + CV12"          
    ##   [2] "O2 ~ b1.2*M1 + b2.2*M2 + b3.2*M3 + b4.2*M4 + b5.2*M5 + c1.2*E1 + c2.2*E2 + c3.2*E3 + c4.2*E4 + CV1 + CV2 + CV3 + CV4 + CV5 + CV6 + CV7 + CV8 + CV9 + CV10 + CV11 + CV12"          
    ##   [3] "O3 ~ b1.3*M1 + b2.3*M2 + b3.3*M3 + b4.3*M4 + b5.3*M5 + c1.3*E1 + c2.3*E2 + c3.3*E3 + c4.3*E4 + CV1 + CV2 + CV3 + CV4 + CV5 + CV6 + CV7 + CV8 + CV9 + CV10 + CV11 + CV12"          
    ##   [4] "O4 ~ b1.4*M1 + b2.4*M2 + b3.4*M3 + b4.4*M4 + b5.4*M5 + c1.4*E1 + c2.4*E2 + c3.4*E3 + c4.4*E4 + CV1 + CV2 + CV3 + CV4 + CV5 + CV6 + CV7 + CV8 + CV9 + CV10 + CV11 + CV12"          
    ##   [5] "O5 ~ b1.5*M1 + b2.5*M2 + b3.5*M3 + b4.5*M4 + b5.5*M5 + c1.5*E1 + c2.5*E2 + c3.5*E3 + c4.5*E4 + CV1 + CV2 + CV3 + CV4 + CV5 + CV6 + CV7 + CV8 + CV9 + CV10 + CV11 + CV12"          
    ##   [6] "O6 ~ b1.6*M1 + b2.6*M2 + b3.6*M3 + b4.6*M4 + b5.6*M5 + c1.6*E1 + c2.6*E2 + c3.6*E3 + c4.6*E4 + CV1 + CV2 + CV3 + CV4 + CV5 + CV6 + CV7 + CV8 + CV9 + CV10 + CV11 + CV12"          
    ##   [7] "O7 ~ b1.7*M1 + b2.7*M2 + b3.7*M3 + b4.7*M4 + b5.7*M5 + c1.7*E1 + c2.7*E2 + c3.7*E3 + c4.7*E4 + CV1 + CV2 + CV3 + CV4 + CV5 + CV6 + CV7 + CV8 + CV9 + CV10 + CV11 + CV12"          
    ##   [8] "O8 ~ b1.8*M1 + b2.8*M2 + b3.8*M3 + b4.8*M4 + b5.8*M5 + c1.8*E1 + c2.8*E2 + c3.8*E3 + c4.8*E4 + CV1 + CV2 + CV3 + CV4 + CV5 + CV6 + CV7 + CV8 + CV9 + CV10 + CV11 + CV12"          
    ##   [9] "O9 ~ b1.9*M1 + b2.9*M2 + b3.9*M3 + b4.9*M4 + b5.9*M5 + c1.9*E1 + c2.9*E2 + c3.9*E3 + c4.9*E4 + CV1 + CV2 + CV3 + CV4 + CV5 + CV6 + CV7 + CV8 + CV9 + CV10 + CV11 + CV12"          
    ##  [10] "O10 ~ b1.10*M1 + b2.10*M2 + b3.10*M3 + b4.10*M4 + b5.10*M5 + c1.10*E1 + c2.10*E2 + c3.10*E3 + c4.10*E4 + CV1 + CV2 + CV3 + CV4 + CV5 + CV6 + CV7 + CV8 + CV9 + CV10 + CV11 + CV12"
    ##  [11] "O11 ~ b1.11*M1 + b2.11*M2 + b3.11*M3 + b4.11*M4 + b5.11*M5 + c1.11*E1 + c2.11*E2 + c3.11*E3 + c4.11*E4 + CV1 + CV2 + CV3 + CV4 + CV5 + CV6 + CV7 + CV8 + CV9 + CV10 + CV11 + CV12"
    ##  [12] "O12 ~ b1.12*M1 + b2.12*M2 + b3.12*M3 + b4.12*M4 + b5.12*M5 + c1.12*E1 + c2.12*E2 + c3.12*E3 + c4.12*E4 + CV1 + CV2 + CV3 + CV4 + CV5 + CV6 + CV7 + CV8 + CV9 + CV10 + CV11 + CV12"
    ##  [13] "O13 ~ b1.13*M1 + b2.13*M2 + b3.13*M3 + b4.13*M4 + b5.13*M5 + c1.13*E1 + c2.13*E2 + c3.13*E3 + c4.13*E4 + CV1 + CV2 + CV3 + CV4 + CV5 + CV6 + CV7 + CV8 + CV9 + CV10 + CV11 + CV12"
    ##  [14] "O14 ~ b1.14*M1 + b2.14*M2 + b3.14*M3 + b4.14*M4 + b5.14*M5 + c1.14*E1 + c2.14*E2 + c3.14*E3 + c4.14*E4 + CV1 + CV2 + CV3 + CV4 + CV5 + CV6 + CV7 + CV8 + CV9 + CV10 + CV11 + CV12"
    ##  [15] "O15 ~ b1.15*M1 + b2.15*M2 + b3.15*M3 + b4.15*M4 + b5.15*M5 + c1.15*E1 + c2.15*E2 + c3.15*E3 + c4.15*E4 + CV1 + CV2 + CV3 + CV4 + CV5 + CV6 + CV7 + CV8 + CV9 + CV10 + CV11 + CV12"
    ##  [16] "O16 ~ b1.16*M1 + b2.16*M2 + b3.16*M3 + b4.16*M4 + b5.16*M5 + c1.16*E1 + c2.16*E2 + c3.16*E3 + c4.16*E4 + CV1 + CV2 + CV3 + CV4 + CV5 + CV6 + CV7 + CV8 + CV9 + CV10 + CV11 + CV12"
    ##  [17] "M1 ~ a1.1*E1 + a2.1*E2 + a3.1*E3 + a4.1*E4 + CV1 + CV2 + CV3 + CV4 + CV5 + CV6 + CV7 + CV8 + CV9 + CV10 + CV11 + CV12"                                                            
    ##  [18] "M2 ~ a1.2*E1 + a2.2*E2 + a3.2*E3 + a4.2*E4 + CV1 + CV2 + CV3 + CV4 + CV5 + CV6 + CV7 + CV8 + CV9 + CV10 + CV11 + CV12"                                                            
    ##  [19] "M3 ~ a1.3*E1 + a2.3*E2 + a3.3*E3 + a4.3*E4 + CV1 + CV2 + CV3 + CV4 + CV5 + CV6 + CV7 + CV8 + CV9 + CV10 + CV11 + CV12"                                                            
    ##  [20] "M4 ~ a1.4*E1 + a2.4*E2 + a3.4*E3 + a4.4*E4 + CV1 + CV2 + CV3 + CV4 + CV5 + CV6 + CV7 + CV8 + CV9 + CV10 + CV11 + CV12"                                                            
    ##  [21] "M5 ~ a1.5*E1 + a2.5*E2 + a3.5*E3 + a4.5*E4 + CV1 + CV2 + CV3 + CV4 + CV5 + CV6 + CV7 + CV8 + CV9 + CV10 + CV11 + CV12"                                                            
    ##  [22] "m1.1 := a1.1*b1.1 + a1.2*b2.1 + a1.3*b3.1 + a1.4*b4.1 + a1.5*b5.1"                                                                                                                
    ##  [23] "m1.2 := a1.1*b1.2 + a1.2*b2.2 + a1.3*b3.2 + a1.4*b4.2 + a1.5*b5.2"                                                                                                                
    ##  [24] "m1.3 := a1.1*b1.3 + a1.2*b2.3 + a1.3*b3.3 + a1.4*b4.3 + a1.5*b5.3"                                                                                                                
    ##  [25] "m1.4 := a1.1*b1.4 + a1.2*b2.4 + a1.3*b3.4 + a1.4*b4.4 + a1.5*b5.4"                                                                                                                
    ##  [26] "m1.5 := a1.1*b1.5 + a1.2*b2.5 + a1.3*b3.5 + a1.4*b4.5 + a1.5*b5.5"                                                                                                                
    ##  [27] "m1.6 := a1.1*b1.6 + a1.2*b2.6 + a1.3*b3.6 + a1.4*b4.6 + a1.5*b5.6"                                                                                                                
    ##  [28] "m1.7 := a1.1*b1.7 + a1.2*b2.7 + a1.3*b3.7 + a1.4*b4.7 + a1.5*b5.7"                                                                                                                
    ##  [29] "m1.8 := a1.1*b1.8 + a1.2*b2.8 + a1.3*b3.8 + a1.4*b4.8 + a1.5*b5.8"                                                                                                                
    ##  [30] "m1.9 := a1.1*b1.9 + a1.2*b2.9 + a1.3*b3.9 + a1.4*b4.9 + a1.5*b5.9"                                                                                                                
    ##  [31] "m1.10 := a1.1*b1.10 + a1.2*b2.10 + a1.3*b3.10 + a1.4*b4.10 + a1.5*b5.10"                                                                                                          
    ##  [32] "m1.11 := a1.1*b1.11 + a1.2*b2.11 + a1.3*b3.11 + a1.4*b4.11 + a1.5*b5.11"                                                                                                          
    ##  [33] "m1.12 := a1.1*b1.12 + a1.2*b2.12 + a1.3*b3.12 + a1.4*b4.12 + a1.5*b5.12"                                                                                                          
    ##  [34] "m1.13 := a1.1*b1.13 + a1.2*b2.13 + a1.3*b3.13 + a1.4*b4.13 + a1.5*b5.13"                                                                                                          
    ##  [35] "m1.14 := a1.1*b1.14 + a1.2*b2.14 + a1.3*b3.14 + a1.4*b4.14 + a1.5*b5.14"                                                                                                          
    ##  [36] "m1.15 := a1.1*b1.15 + a1.2*b2.15 + a1.3*b3.15 + a1.4*b4.15 + a1.5*b5.15"                                                                                                          
    ##  [37] "m1.16 := a1.1*b1.16 + a1.2*b2.16 + a1.3*b3.16 + a1.4*b4.16 + a1.5*b5.16"                                                                                                          
    ##  [38] "m2.1 := a2.1*b1.1 + a2.2*b2.1 + a2.3*b3.1 + a2.4*b4.1 + a2.5*b5.1"                                                                                                                
    ##  [39] "m2.2 := a2.1*b1.2 + a2.2*b2.2 + a2.3*b3.2 + a2.4*b4.2 + a2.5*b5.2"                                                                                                                
    ##  [40] "m2.3 := a2.1*b1.3 + a2.2*b2.3 + a2.3*b3.3 + a2.4*b4.3 + a2.5*b5.3"                                                                                                                
    ##  [41] "m2.4 := a2.1*b1.4 + a2.2*b2.4 + a2.3*b3.4 + a2.4*b4.4 + a2.5*b5.4"                                                                                                                
    ##  [42] "m2.5 := a2.1*b1.5 + a2.2*b2.5 + a2.3*b3.5 + a2.4*b4.5 + a2.5*b5.5"                                                                                                                
    ##  [43] "m2.6 := a2.1*b1.6 + a2.2*b2.6 + a2.3*b3.6 + a2.4*b4.6 + a2.5*b5.6"                                                                                                                
    ##  [44] "m2.7 := a2.1*b1.7 + a2.2*b2.7 + a2.3*b3.7 + a2.4*b4.7 + a2.5*b5.7"                                                                                                                
    ##  [45] "m2.8 := a2.1*b1.8 + a2.2*b2.8 + a2.3*b3.8 + a2.4*b4.8 + a2.5*b5.8"                                                                                                                
    ##  [46] "m2.9 := a2.1*b1.9 + a2.2*b2.9 + a2.3*b3.9 + a2.4*b4.9 + a2.5*b5.9"                                                                                                                
    ##  [47] "m2.10 := a2.1*b1.10 + a2.2*b2.10 + a2.3*b3.10 + a2.4*b4.10 + a2.5*b5.10"                                                                                                          
    ##  [48] "m2.11 := a2.1*b1.11 + a2.2*b2.11 + a2.3*b3.11 + a2.4*b4.11 + a2.5*b5.11"                                                                                                          
    ##  [49] "m2.12 := a2.1*b1.12 + a2.2*b2.12 + a2.3*b3.12 + a2.4*b4.12 + a2.5*b5.12"                                                                                                          
    ##  [50] "m2.13 := a2.1*b1.13 + a2.2*b2.13 + a2.3*b3.13 + a2.4*b4.13 + a2.5*b5.13"                                                                                                          
    ##  [51] "m2.14 := a2.1*b1.14 + a2.2*b2.14 + a2.3*b3.14 + a2.4*b4.14 + a2.5*b5.14"                                                                                                          
    ##  [52] "m2.15 := a2.1*b1.15 + a2.2*b2.15 + a2.3*b3.15 + a2.4*b4.15 + a2.5*b5.15"                                                                                                          
    ##  [53] "m2.16 := a2.1*b1.16 + a2.2*b2.16 + a2.3*b3.16 + a2.4*b4.16 + a2.5*b5.16"                                                                                                          
    ##  [54] "m3.1 := a3.1*b1.1 + a3.2*b2.1 + a3.3*b3.1 + a3.4*b4.1 + a3.5*b5.1"                                                                                                                
    ##  [55] "m3.2 := a3.1*b1.2 + a3.2*b2.2 + a3.3*b3.2 + a3.4*b4.2 + a3.5*b5.2"                                                                                                                
    ##  [56] "m3.3 := a3.1*b1.3 + a3.2*b2.3 + a3.3*b3.3 + a3.4*b4.3 + a3.5*b5.3"                                                                                                                
    ##  [57] "m3.4 := a3.1*b1.4 + a3.2*b2.4 + a3.3*b3.4 + a3.4*b4.4 + a3.5*b5.4"                                                                                                                
    ##  [58] "m3.5 := a3.1*b1.5 + a3.2*b2.5 + a3.3*b3.5 + a3.4*b4.5 + a3.5*b5.5"                                                                                                                
    ##  [59] "m3.6 := a3.1*b1.6 + a3.2*b2.6 + a3.3*b3.6 + a3.4*b4.6 + a3.5*b5.6"                                                                                                                
    ##  [60] "m3.7 := a3.1*b1.7 + a3.2*b2.7 + a3.3*b3.7 + a3.4*b4.7 + a3.5*b5.7"                                                                                                                
    ##  [61] "m3.8 := a3.1*b1.8 + a3.2*b2.8 + a3.3*b3.8 + a3.4*b4.8 + a3.5*b5.8"                                                                                                                
    ##  [62] "m3.9 := a3.1*b1.9 + a3.2*b2.9 + a3.3*b3.9 + a3.4*b4.9 + a3.5*b5.9"                                                                                                                
    ##  [63] "m3.10 := a3.1*b1.10 + a3.2*b2.10 + a3.3*b3.10 + a3.4*b4.10 + a3.5*b5.10"                                                                                                          
    ##  [64] "m3.11 := a3.1*b1.11 + a3.2*b2.11 + a3.3*b3.11 + a3.4*b4.11 + a3.5*b5.11"                                                                                                          
    ##  [65] "m3.12 := a3.1*b1.12 + a3.2*b2.12 + a3.3*b3.12 + a3.4*b4.12 + a3.5*b5.12"                                                                                                          
    ##  [66] "m3.13 := a3.1*b1.13 + a3.2*b2.13 + a3.3*b3.13 + a3.4*b4.13 + a3.5*b5.13"                                                                                                          
    ##  [67] "m3.14 := a3.1*b1.14 + a3.2*b2.14 + a3.3*b3.14 + a3.4*b4.14 + a3.5*b5.14"                                                                                                          
    ##  [68] "m3.15 := a3.1*b1.15 + a3.2*b2.15 + a3.3*b3.15 + a3.4*b4.15 + a3.5*b5.15"                                                                                                          
    ##  [69] "m3.16 := a3.1*b1.16 + a3.2*b2.16 + a3.3*b3.16 + a3.4*b4.16 + a3.5*b5.16"                                                                                                          
    ##  [70] "m4.1 := a4.1*b1.1 + a4.2*b2.1 + a4.3*b3.1 + a4.4*b4.1 + a4.5*b5.1"                                                                                                                
    ##  [71] "m4.2 := a4.1*b1.2 + a4.2*b2.2 + a4.3*b3.2 + a4.4*b4.2 + a4.5*b5.2"                                                                                                                
    ##  [72] "m4.3 := a4.1*b1.3 + a4.2*b2.3 + a4.3*b3.3 + a4.4*b4.3 + a4.5*b5.3"                                                                                                                
    ##  [73] "m4.4 := a4.1*b1.4 + a4.2*b2.4 + a4.3*b3.4 + a4.4*b4.4 + a4.5*b5.4"                                                                                                                
    ##  [74] "m4.5 := a4.1*b1.5 + a4.2*b2.5 + a4.3*b3.5 + a4.4*b4.5 + a4.5*b5.5"                                                                                                                
    ##  [75] "m4.6 := a4.1*b1.6 + a4.2*b2.6 + a4.3*b3.6 + a4.4*b4.6 + a4.5*b5.6"                                                                                                                
    ##  [76] "m4.7 := a4.1*b1.7 + a4.2*b2.7 + a4.3*b3.7 + a4.4*b4.7 + a4.5*b5.7"                                                                                                                
    ##  [77] "m4.8 := a4.1*b1.8 + a4.2*b2.8 + a4.3*b3.8 + a4.4*b4.8 + a4.5*b5.8"                                                                                                                
    ##  [78] "m4.9 := a4.1*b1.9 + a4.2*b2.9 + a4.3*b3.9 + a4.4*b4.9 + a4.5*b5.9"                                                                                                                
    ##  [79] "m4.10 := a4.1*b1.10 + a4.2*b2.10 + a4.3*b3.10 + a4.4*b4.10 + a4.5*b5.10"                                                                                                          
    ##  [80] "m4.11 := a4.1*b1.11 + a4.2*b2.11 + a4.3*b3.11 + a4.4*b4.11 + a4.5*b5.11"                                                                                                          
    ##  [81] "m4.12 := a4.1*b1.12 + a4.2*b2.12 + a4.3*b3.12 + a4.4*b4.12 + a4.5*b5.12"                                                                                                          
    ##  [82] "m4.13 := a4.1*b1.13 + a4.2*b2.13 + a4.3*b3.13 + a4.4*b4.13 + a4.5*b5.13"                                                                                                          
    ##  [83] "m4.14 := a4.1*b1.14 + a4.2*b2.14 + a4.3*b3.14 + a4.4*b4.14 + a4.5*b5.14"                                                                                                          
    ##  [84] "m4.15 := a4.1*b1.15 + a4.2*b2.15 + a4.3*b3.15 + a4.4*b4.15 + a4.5*b5.15"                                                                                                          
    ##  [85] "m4.16 := a4.1*b1.16 + a4.2*b2.16 + a4.3*b3.16 + a4.4*b4.16 + a4.5*b5.16"                                                                                                          
    ##  [86] "m1.1.1 := a1.1*b1.1"                                                                                                                                                              
    ##  [87] "m1.1.2 := a1.1*b1.2"                                                                                                                                                              
    ##  [88] "m1.1.3 := a1.1*b1.3"                                                                                                                                                              
    ##  [89] "m1.1.4 := a1.1*b1.4"                                                                                                                                                              
    ##  [90] "m1.1.5 := a1.1*b1.5"                                                                                                                                                              
    ##  [91] "m1.1.6 := a1.1*b1.6"                                                                                                                                                              
    ##  [92] "m1.1.7 := a1.1*b1.7"                                                                                                                                                              
    ##  [93] "m1.1.8 := a1.1*b1.8"                                                                                                                                                              
    ##  [94] "m1.1.9 := a1.1*b1.9"                                                                                                                                                              
    ##  [95] "m1.1.10 := a1.1*b1.10"                                                                                                                                                            
    ##  [96] "m1.1.11 := a1.1*b1.11"                                                                                                                                                            
    ##  [97] "m1.1.12 := a1.1*b1.12"                                                                                                                                                            
    ##  [98] "m1.1.13 := a1.1*b1.13"                                                                                                                                                            
    ##  [99] "m1.1.14 := a1.1*b1.14"                                                                                                                                                            
    ## [100] "m1.1.15 := a1.1*b1.15"                                                                                                                                                            
    ## [101] "m1.1.16 := a1.1*b1.16"                                                                                                                                                            
    ## [102] "m2.1.1 := a2.1*b1.1"                                                                                                                                                              
    ## [103] "m2.1.2 := a2.1*b1.2"                                                                                                                                                              
    ## [104] "m2.1.3 := a2.1*b1.3"                                                                                                                                                              
    ## [105] "m2.1.4 := a2.1*b1.4"                                                                                                                                                              
    ## [106] "m2.1.5 := a2.1*b1.5"                                                                                                                                                              
    ## [107] "m2.1.6 := a2.1*b1.6"                                                                                                                                                              
    ## [108] "m2.1.7 := a2.1*b1.7"                                                                                                                                                              
    ## [109] "m2.1.8 := a2.1*b1.8"                                                                                                                                                              
    ## [110] "m2.1.9 := a2.1*b1.9"                                                                                                                                                              
    ## [111] "m2.1.10 := a2.1*b1.10"                                                                                                                                                            
    ## [112] "m2.1.11 := a2.1*b1.11"                                                                                                                                                            
    ## [113] "m2.1.12 := a2.1*b1.12"                                                                                                                                                            
    ## [114] "m2.1.13 := a2.1*b1.13"                                                                                                                                                            
    ## [115] "m2.1.14 := a2.1*b1.14"                                                                                                                                                            
    ## [116] "m2.1.15 := a2.1*b1.15"                                                                                                                                                            
    ## [117] "m2.1.16 := a2.1*b1.16"                                                                                                                                                            
    ## [118] "m3.1.1 := a3.1*b1.1"                                                                                                                                                              
    ## [119] "m3.1.2 := a3.1*b1.2"                                                                                                                                                              
    ## [120] "m3.1.3 := a3.1*b1.3"                                                                                                                                                              
    ## [121] "m3.1.4 := a3.1*b1.4"                                                                                                                                                              
    ## [122] "m3.1.5 := a3.1*b1.5"                                                                                                                                                              
    ## [123] "m3.1.6 := a3.1*b1.6"                                                                                                                                                              
    ## [124] "m3.1.7 := a3.1*b1.7"                                                                                                                                                              
    ## [125] "m3.1.8 := a3.1*b1.8"                                                                                                                                                              
    ## [126] "m3.1.9 := a3.1*b1.9"                                                                                                                                                              
    ## [127] "m3.1.10 := a3.1*b1.10"                                                                                                                                                            
    ## [128] "m3.1.11 := a3.1*b1.11"                                                                                                                                                            
    ## [129] "m3.1.12 := a3.1*b1.12"                                                                                                                                                            
    ## [130] "m3.1.13 := a3.1*b1.13"                                                                                                                                                            
    ## [131] "m3.1.14 := a3.1*b1.14"                                                                                                                                                            
    ## [132] "m3.1.15 := a3.1*b1.15"                                                                                                                                                            
    ## [133] "m3.1.16 := a3.1*b1.16"                                                                                                                                                            
    ## [134] "m4.1.1 := a4.1*b1.1"                                                                                                                                                              
    ## [135] "m4.1.2 := a4.1*b1.2"                                                                                                                                                              
    ## [136] "m4.1.3 := a4.1*b1.3"                                                                                                                                                              
    ## [137] "m4.1.4 := a4.1*b1.4"                                                                                                                                                              
    ## [138] "m4.1.5 := a4.1*b1.5"                                                                                                                                                              
    ## [139] "m4.1.6 := a4.1*b1.6"                                                                                                                                                              
    ## [140] "m4.1.7 := a4.1*b1.7"                                                                                                                                                              
    ## [141] "m4.1.8 := a4.1*b1.8"                                                                                                                                                              
    ## [142] "m4.1.9 := a4.1*b1.9"                                                                                                                                                              
    ## [143] "m4.1.10 := a4.1*b1.10"                                                                                                                                                            
    ## [144] "m4.1.11 := a4.1*b1.11"                                                                                                                                                            
    ## [145] "m4.1.12 := a4.1*b1.12"                                                                                                                                                            
    ## [146] "m4.1.13 := a4.1*b1.13"                                                                                                                                                            
    ## [147] "m4.1.14 := a4.1*b1.14"                                                                                                                                                            
    ## [148] "m4.1.15 := a4.1*b1.15"                                                                                                                                                            
    ## [149] "m4.1.16 := a4.1*b1.16"                                                                                                                                                            
    ## [150] "m1.2.1 := a1.2*b2.1"                                                                                                                                                              
    ## [151] "m1.2.2 := a1.2*b2.2"                                                                                                                                                              
    ## [152] "m1.2.3 := a1.2*b2.3"                                                                                                                                                              
    ## [153] "m1.2.4 := a1.2*b2.4"                                                                                                                                                              
    ## [154] "m1.2.5 := a1.2*b2.5"                                                                                                                                                              
    ## [155] "m1.2.6 := a1.2*b2.6"                                                                                                                                                              
    ## [156] "m1.2.7 := a1.2*b2.7"                                                                                                                                                              
    ## [157] "m1.2.8 := a1.2*b2.8"                                                                                                                                                              
    ## [158] "m1.2.9 := a1.2*b2.9"                                                                                                                                                              
    ## [159] "m1.2.10 := a1.2*b2.10"                                                                                                                                                            
    ## [160] "m1.2.11 := a1.2*b2.11"                                                                                                                                                            
    ## [161] "m1.2.12 := a1.2*b2.12"                                                                                                                                                            
    ## [162] "m1.2.13 := a1.2*b2.13"                                                                                                                                                            
    ## [163] "m1.2.14 := a1.2*b2.14"                                                                                                                                                            
    ## [164] "m1.2.15 := a1.2*b2.15"                                                                                                                                                            
    ## [165] "m1.2.16 := a1.2*b2.16"                                                                                                                                                            
    ## [166] "m2.2.1 := a2.2*b2.1"                                                                                                                                                              
    ## [167] "m2.2.2 := a2.2*b2.2"                                                                                                                                                              
    ## [168] "m2.2.3 := a2.2*b2.3"                                                                                                                                                              
    ## [169] "m2.2.4 := a2.2*b2.4"                                                                                                                                                              
    ## [170] "m2.2.5 := a2.2*b2.5"                                                                                                                                                              
    ## [171] "m2.2.6 := a2.2*b2.6"                                                                                                                                                              
    ## [172] "m2.2.7 := a2.2*b2.7"                                                                                                                                                              
    ## [173] "m2.2.8 := a2.2*b2.8"                                                                                                                                                              
    ## [174] "m2.2.9 := a2.2*b2.9"                                                                                                                                                              
    ## [175] "m2.2.10 := a2.2*b2.10"                                                                                                                                                            
    ## [176] "m2.2.11 := a2.2*b2.11"                                                                                                                                                            
    ## [177] "m2.2.12 := a2.2*b2.12"                                                                                                                                                            
    ## [178] "m2.2.13 := a2.2*b2.13"                                                                                                                                                            
    ## [179] "m2.2.14 := a2.2*b2.14"                                                                                                                                                            
    ## [180] "m2.2.15 := a2.2*b2.15"                                                                                                                                                            
    ## [181] "m2.2.16 := a2.2*b2.16"                                                                                                                                                            
    ## [182] "m3.2.1 := a3.2*b2.1"                                                                                                                                                              
    ## [183] "m3.2.2 := a3.2*b2.2"                                                                                                                                                              
    ## [184] "m3.2.3 := a3.2*b2.3"                                                                                                                                                              
    ## [185] "m3.2.4 := a3.2*b2.4"                                                                                                                                                              
    ## [186] "m3.2.5 := a3.2*b2.5"                                                                                                                                                              
    ## [187] "m3.2.6 := a3.2*b2.6"                                                                                                                                                              
    ## [188] "m3.2.7 := a3.2*b2.7"                                                                                                                                                              
    ## [189] "m3.2.8 := a3.2*b2.8"                                                                                                                                                              
    ## [190] "m3.2.9 := a3.2*b2.9"                                                                                                                                                              
    ## [191] "m3.2.10 := a3.2*b2.10"                                                                                                                                                            
    ## [192] "m3.2.11 := a3.2*b2.11"                                                                                                                                                            
    ## [193] "m3.2.12 := a3.2*b2.12"                                                                                                                                                            
    ## [194] "m3.2.13 := a3.2*b2.13"                                                                                                                                                            
    ## [195] "m3.2.14 := a3.2*b2.14"                                                                                                                                                            
    ## [196] "m3.2.15 := a3.2*b2.15"                                                                                                                                                            
    ## [197] "m3.2.16 := a3.2*b2.16"                                                                                                                                                            
    ## [198] "m4.2.1 := a4.2*b2.1"                                                                                                                                                              
    ## [199] "m4.2.2 := a4.2*b2.2"                                                                                                                                                              
    ## [200] "m4.2.3 := a4.2*b2.3"                                                                                                                                                              
    ## [201] "m4.2.4 := a4.2*b2.4"                                                                                                                                                              
    ## [202] "m4.2.5 := a4.2*b2.5"                                                                                                                                                              
    ## [203] "m4.2.6 := a4.2*b2.6"                                                                                                                                                              
    ## [204] "m4.2.7 := a4.2*b2.7"                                                                                                                                                              
    ## [205] "m4.2.8 := a4.2*b2.8"                                                                                                                                                              
    ## [206] "m4.2.9 := a4.2*b2.9"                                                                                                                                                              
    ## [207] "m4.2.10 := a4.2*b2.10"                                                                                                                                                            
    ## [208] "m4.2.11 := a4.2*b2.11"                                                                                                                                                            
    ## [209] "m4.2.12 := a4.2*b2.12"                                                                                                                                                            
    ## [210] "m4.2.13 := a4.2*b2.13"                                                                                                                                                            
    ## [211] "m4.2.14 := a4.2*b2.14"                                                                                                                                                            
    ## [212] "m4.2.15 := a4.2*b2.15"                                                                                                                                                            
    ## [213] "m4.2.16 := a4.2*b2.16"                                                                                                                                                            
    ## [214] "m1.3.1 := a1.3*b3.1"                                                                                                                                                              
    ## [215] "m1.3.2 := a1.3*b3.2"                                                                                                                                                              
    ## [216] "m1.3.3 := a1.3*b3.3"                                                                                                                                                              
    ## [217] "m1.3.4 := a1.3*b3.4"                                                                                                                                                              
    ## [218] "m1.3.5 := a1.3*b3.5"                                                                                                                                                              
    ## [219] "m1.3.6 := a1.3*b3.6"                                                                                                                                                              
    ## [220] "m1.3.7 := a1.3*b3.7"                                                                                                                                                              
    ## [221] "m1.3.8 := a1.3*b3.8"                                                                                                                                                              
    ## [222] "m1.3.9 := a1.3*b3.9"                                                                                                                                                              
    ## [223] "m1.3.10 := a1.3*b3.10"                                                                                                                                                            
    ## [224] "m1.3.11 := a1.3*b3.11"                                                                                                                                                            
    ## [225] "m1.3.12 := a1.3*b3.12"                                                                                                                                                            
    ## [226] "m1.3.13 := a1.3*b3.13"                                                                                                                                                            
    ## [227] "m1.3.14 := a1.3*b3.14"                                                                                                                                                            
    ## [228] "m1.3.15 := a1.3*b3.15"                                                                                                                                                            
    ## [229] "m1.3.16 := a1.3*b3.16"                                                                                                                                                            
    ## [230] "m2.3.1 := a2.3*b3.1"                                                                                                                                                              
    ## [231] "m2.3.2 := a2.3*b3.2"                                                                                                                                                              
    ## [232] "m2.3.3 := a2.3*b3.3"                                                                                                                                                              
    ## [233] "m2.3.4 := a2.3*b3.4"                                                                                                                                                              
    ## [234] "m2.3.5 := a2.3*b3.5"                                                                                                                                                              
    ## [235] "m2.3.6 := a2.3*b3.6"                                                                                                                                                              
    ## [236] "m2.3.7 := a2.3*b3.7"                                                                                                                                                              
    ## [237] "m2.3.8 := a2.3*b3.8"                                                                                                                                                              
    ## [238] "m2.3.9 := a2.3*b3.9"                                                                                                                                                              
    ## [239] "m2.3.10 := a2.3*b3.10"                                                                                                                                                            
    ## [240] "m2.3.11 := a2.3*b3.11"                                                                                                                                                            
    ## [241] "m2.3.12 := a2.3*b3.12"                                                                                                                                                            
    ## [242] "m2.3.13 := a2.3*b3.13"                                                                                                                                                            
    ## [243] "m2.3.14 := a2.3*b3.14"                                                                                                                                                            
    ## [244] "m2.3.15 := a2.3*b3.15"                                                                                                                                                            
    ## [245] "m2.3.16 := a2.3*b3.16"                                                                                                                                                            
    ## [246] "m3.3.1 := a3.3*b3.1"                                                                                                                                                              
    ## [247] "m3.3.2 := a3.3*b3.2"                                                                                                                                                              
    ## [248] "m3.3.3 := a3.3*b3.3"                                                                                                                                                              
    ## [249] "m3.3.4 := a3.3*b3.4"                                                                                                                                                              
    ## [250] "m3.3.5 := a3.3*b3.5"                                                                                                                                                              
    ## [251] "m3.3.6 := a3.3*b3.6"                                                                                                                                                              
    ## [252] "m3.3.7 := a3.3*b3.7"                                                                                                                                                              
    ## [253] "m3.3.8 := a3.3*b3.8"                                                                                                                                                              
    ## [254] "m3.3.9 := a3.3*b3.9"                                                                                                                                                              
    ## [255] "m3.3.10 := a3.3*b3.10"                                                                                                                                                            
    ## [256] "m3.3.11 := a3.3*b3.11"                                                                                                                                                            
    ## [257] "m3.3.12 := a3.3*b3.12"                                                                                                                                                            
    ## [258] "m3.3.13 := a3.3*b3.13"                                                                                                                                                            
    ## [259] "m3.3.14 := a3.3*b3.14"                                                                                                                                                            
    ## [260] "m3.3.15 := a3.3*b3.15"                                                                                                                                                            
    ## [261] "m3.3.16 := a3.3*b3.16"                                                                                                                                                            
    ## [262] "m4.3.1 := a4.3*b3.1"                                                                                                                                                              
    ## [263] "m4.3.2 := a4.3*b3.2"                                                                                                                                                              
    ## [264] "m4.3.3 := a4.3*b3.3"                                                                                                                                                              
    ## [265] "m4.3.4 := a4.3*b3.4"                                                                                                                                                              
    ## [266] "m4.3.5 := a4.3*b3.5"                                                                                                                                                              
    ## [267] "m4.3.6 := a4.3*b3.6"                                                                                                                                                              
    ## [268] "m4.3.7 := a4.3*b3.7"                                                                                                                                                              
    ## [269] "m4.3.8 := a4.3*b3.8"                                                                                                                                                              
    ## [270] "m4.3.9 := a4.3*b3.9"                                                                                                                                                              
    ## [271] "m4.3.10 := a4.3*b3.10"                                                                                                                                                            
    ## [272] "m4.3.11 := a4.3*b3.11"                                                                                                                                                            
    ## [273] "m4.3.12 := a4.3*b3.12"                                                                                                                                                            
    ## [274] "m4.3.13 := a4.3*b3.13"                                                                                                                                                            
    ## [275] "m4.3.14 := a4.3*b3.14"                                                                                                                                                            
    ## [276] "m4.3.15 := a4.3*b3.15"                                                                                                                                                            
    ## [277] "m4.3.16 := a4.3*b3.16"                                                                                                                                                            
    ## [278] "m1.4.1 := a1.4*b4.1"                                                                                                                                                              
    ## [279] "m1.4.2 := a1.4*b4.2"                                                                                                                                                              
    ## [280] "m1.4.3 := a1.4*b4.3"                                                                                                                                                              
    ## [281] "m1.4.4 := a1.4*b4.4"                                                                                                                                                              
    ## [282] "m1.4.5 := a1.4*b4.5"                                                                                                                                                              
    ## [283] "m1.4.6 := a1.4*b4.6"                                                                                                                                                              
    ## [284] "m1.4.7 := a1.4*b4.7"                                                                                                                                                              
    ## [285] "m1.4.8 := a1.4*b4.8"                                                                                                                                                              
    ## [286] "m1.4.9 := a1.4*b4.9"                                                                                                                                                              
    ## [287] "m1.4.10 := a1.4*b4.10"                                                                                                                                                            
    ## [288] "m1.4.11 := a1.4*b4.11"                                                                                                                                                            
    ## [289] "m1.4.12 := a1.4*b4.12"                                                                                                                                                            
    ## [290] "m1.4.13 := a1.4*b4.13"                                                                                                                                                            
    ## [291] "m1.4.14 := a1.4*b4.14"                                                                                                                                                            
    ## [292] "m1.4.15 := a1.4*b4.15"                                                                                                                                                            
    ## [293] "m1.4.16 := a1.4*b4.16"                                                                                                                                                            
    ## [294] "m2.4.1 := a2.4*b4.1"                                                                                                                                                              
    ## [295] "m2.4.2 := a2.4*b4.2"                                                                                                                                                              
    ## [296] "m2.4.3 := a2.4*b4.3"                                                                                                                                                              
    ## [297] "m2.4.4 := a2.4*b4.4"                                                                                                                                                              
    ## [298] "m2.4.5 := a2.4*b4.5"                                                                                                                                                              
    ## [299] "m2.4.6 := a2.4*b4.6"                                                                                                                                                              
    ## [300] "m2.4.7 := a2.4*b4.7"                                                                                                                                                              
    ## [301] "m2.4.8 := a2.4*b4.8"                                                                                                                                                              
    ## [302] "m2.4.9 := a2.4*b4.9"                                                                                                                                                              
    ## [303] "m2.4.10 := a2.4*b4.10"                                                                                                                                                            
    ## [304] "m2.4.11 := a2.4*b4.11"                                                                                                                                                            
    ## [305] "m2.4.12 := a2.4*b4.12"                                                                                                                                                            
    ## [306] "m2.4.13 := a2.4*b4.13"                                                                                                                                                            
    ## [307] "m2.4.14 := a2.4*b4.14"                                                                                                                                                            
    ## [308] "m2.4.15 := a2.4*b4.15"                                                                                                                                                            
    ## [309] "m2.4.16 := a2.4*b4.16"                                                                                                                                                            
    ## [310] "m3.4.1 := a3.4*b4.1"                                                                                                                                                              
    ## [311] "m3.4.2 := a3.4*b4.2"                                                                                                                                                              
    ## [312] "m3.4.3 := a3.4*b4.3"                                                                                                                                                              
    ## [313] "m3.4.4 := a3.4*b4.4"                                                                                                                                                              
    ## [314] "m3.4.5 := a3.4*b4.5"                                                                                                                                                              
    ## [315] "m3.4.6 := a3.4*b4.6"                                                                                                                                                              
    ## [316] "m3.4.7 := a3.4*b4.7"                                                                                                                                                              
    ## [317] "m3.4.8 := a3.4*b4.8"                                                                                                                                                              
    ## [318] "m3.4.9 := a3.4*b4.9"                                                                                                                                                              
    ## [319] "m3.4.10 := a3.4*b4.10"                                                                                                                                                            
    ## [320] "m3.4.11 := a3.4*b4.11"                                                                                                                                                            
    ## [321] "m3.4.12 := a3.4*b4.12"                                                                                                                                                            
    ## [322] "m3.4.13 := a3.4*b4.13"                                                                                                                                                            
    ## [323] "m3.4.14 := a3.4*b4.14"                                                                                                                                                            
    ## [324] "m3.4.15 := a3.4*b4.15"                                                                                                                                                            
    ## [325] "m3.4.16 := a3.4*b4.16"                                                                                                                                                            
    ## [326] "m4.4.1 := a4.4*b4.1"                                                                                                                                                              
    ## [327] "m4.4.2 := a4.4*b4.2"                                                                                                                                                              
    ## [328] "m4.4.3 := a4.4*b4.3"                                                                                                                                                              
    ## [329] "m4.4.4 := a4.4*b4.4"                                                                                                                                                              
    ## [330] "m4.4.5 := a4.4*b4.5"                                                                                                                                                              
    ## [331] "m4.4.6 := a4.4*b4.6"                                                                                                                                                              
    ## [332] "m4.4.7 := a4.4*b4.7"                                                                                                                                                              
    ## [333] "m4.4.8 := a4.4*b4.8"                                                                                                                                                              
    ## [334] "m4.4.9 := a4.4*b4.9"                                                                                                                                                              
    ## [335] "m4.4.10 := a4.4*b4.10"                                                                                                                                                            
    ## [336] "m4.4.11 := a4.4*b4.11"                                                                                                                                                            
    ## [337] "m4.4.12 := a4.4*b4.12"                                                                                                                                                            
    ## [338] "m4.4.13 := a4.4*b4.13"                                                                                                                                                            
    ## [339] "m4.4.14 := a4.4*b4.14"                                                                                                                                                            
    ## [340] "m4.4.15 := a4.4*b4.15"                                                                                                                                                            
    ## [341] "m4.4.16 := a4.4*b4.16"                                                                                                                                                            
    ## [342] "m1.5.1 := a1.5*b5.1"                                                                                                                                                              
    ## [343] "m1.5.2 := a1.5*b5.2"                                                                                                                                                              
    ## [344] "m1.5.3 := a1.5*b5.3"                                                                                                                                                              
    ## [345] "m1.5.4 := a1.5*b5.4"                                                                                                                                                              
    ## [346] "m1.5.5 := a1.5*b5.5"                                                                                                                                                              
    ## [347] "m1.5.6 := a1.5*b5.6"                                                                                                                                                              
    ## [348] "m1.5.7 := a1.5*b5.7"                                                                                                                                                              
    ## [349] "m1.5.8 := a1.5*b5.8"                                                                                                                                                              
    ## [350] "m1.5.9 := a1.5*b5.9"                                                                                                                                                              
    ## [351] "m1.5.10 := a1.5*b5.10"                                                                                                                                                            
    ## [352] "m1.5.11 := a1.5*b5.11"                                                                                                                                                            
    ## [353] "m1.5.12 := a1.5*b5.12"                                                                                                                                                            
    ## [354] "m1.5.13 := a1.5*b5.13"                                                                                                                                                            
    ## [355] "m1.5.14 := a1.5*b5.14"                                                                                                                                                            
    ## [356] "m1.5.15 := a1.5*b5.15"                                                                                                                                                            
    ## [357] "m1.5.16 := a1.5*b5.16"                                                                                                                                                            
    ## [358] "m2.5.1 := a2.5*b5.1"                                                                                                                                                              
    ## [359] "m2.5.2 := a2.5*b5.2"                                                                                                                                                              
    ## [360] "m2.5.3 := a2.5*b5.3"                                                                                                                                                              
    ## [361] "m2.5.4 := a2.5*b5.4"                                                                                                                                                              
    ## [362] "m2.5.5 := a2.5*b5.5"                                                                                                                                                              
    ## [363] "m2.5.6 := a2.5*b5.6"                                                                                                                                                              
    ## [364] "m2.5.7 := a2.5*b5.7"                                                                                                                                                              
    ## [365] "m2.5.8 := a2.5*b5.8"                                                                                                                                                              
    ## [366] "m2.5.9 := a2.5*b5.9"                                                                                                                                                              
    ## [367] "m2.5.10 := a2.5*b5.10"                                                                                                                                                            
    ## [368] "m2.5.11 := a2.5*b5.11"                                                                                                                                                            
    ## [369] "m2.5.12 := a2.5*b5.12"                                                                                                                                                            
    ## [370] "m2.5.13 := a2.5*b5.13"                                                                                                                                                            
    ## [371] "m2.5.14 := a2.5*b5.14"                                                                                                                                                            
    ## [372] "m2.5.15 := a2.5*b5.15"                                                                                                                                                            
    ## [373] "m2.5.16 := a2.5*b5.16"                                                                                                                                                            
    ## [374] "m3.5.1 := a3.5*b5.1"                                                                                                                                                              
    ## [375] "m3.5.2 := a3.5*b5.2"                                                                                                                                                              
    ## [376] "m3.5.3 := a3.5*b5.3"                                                                                                                                                              
    ## [377] "m3.5.4 := a3.5*b5.4"                                                                                                                                                              
    ## [378] "m3.5.5 := a3.5*b5.5"                                                                                                                                                              
    ## [379] "m3.5.6 := a3.5*b5.6"                                                                                                                                                              
    ## [380] "m3.5.7 := a3.5*b5.7"                                                                                                                                                              
    ## [381] "m3.5.8 := a3.5*b5.8"                                                                                                                                                              
    ## [382] "m3.5.9 := a3.5*b5.9"                                                                                                                                                              
    ## [383] "m3.5.10 := a3.5*b5.10"                                                                                                                                                            
    ## [384] "m3.5.11 := a3.5*b5.11"                                                                                                                                                            
    ## [385] "m3.5.12 := a3.5*b5.12"                                                                                                                                                            
    ## [386] "m3.5.13 := a3.5*b5.13"                                                                                                                                                            
    ## [387] "m3.5.14 := a3.5*b5.14"                                                                                                                                                            
    ## [388] "m3.5.15 := a3.5*b5.15"                                                                                                                                                            
    ## [389] "m3.5.16 := a3.5*b5.16"                                                                                                                                                            
    ## [390] "m4.5.1 := a4.5*b5.1"                                                                                                                                                              
    ## [391] "m4.5.2 := a4.5*b5.2"                                                                                                                                                              
    ## [392] "m4.5.3 := a4.5*b5.3"                                                                                                                                                              
    ## [393] "m4.5.4 := a4.5*b5.4"                                                                                                                                                              
    ## [394] "m4.5.5 := a4.5*b5.5"                                                                                                                                                              
    ## [395] "m4.5.6 := a4.5*b5.6"                                                                                                                                                              
    ## [396] "m4.5.7 := a4.5*b5.7"                                                                                                                                                              
    ## [397] "m4.5.8 := a4.5*b5.8"                                                                                                                                                              
    ## [398] "m4.5.9 := a4.5*b5.9"                                                                                                                                                              
    ## [399] "m4.5.10 := a4.5*b5.10"                                                                                                                                                            
    ## [400] "m4.5.11 := a4.5*b5.11"                                                                                                                                                            
    ## [401] "m4.5.12 := a4.5*b5.12"                                                                                                                                                            
    ## [402] "m4.5.13 := a4.5*b5.13"                                                                                                                                                            
    ## [403] "m4.5.14 := a4.5*b5.14"                                                                                                                                                            
    ## [404] "m4.5.15 := a4.5*b5.15"                                                                                                                                                            
    ## [405] "m4.5.16 := a4.5*b5.16"                                                                                                                                                            
    ## [406] "tot1.1 := a1.1*b1.1 + a1.2*b2.1 + a1.3*b3.1 + a1.4*b4.1 + a1.5*b5.1 + c1.1"                                                                                                       
    ## [407] "tot1.2 := a1.1*b1.2 + a1.2*b2.2 + a1.3*b3.2 + a1.4*b4.2 + a1.5*b5.2 + c1.2"                                                                                                       
    ## [408] "tot1.3 := a1.1*b1.3 + a1.2*b2.3 + a1.3*b3.3 + a1.4*b4.3 + a1.5*b5.3 + c1.3"                                                                                                       
    ## [409] "tot1.4 := a1.1*b1.4 + a1.2*b2.4 + a1.3*b3.4 + a1.4*b4.4 + a1.5*b5.4 + c1.4"                                                                                                       
    ## [410] "tot1.5 := a1.1*b1.5 + a1.2*b2.5 + a1.3*b3.5 + a1.4*b4.5 + a1.5*b5.5 + c1.5"                                                                                                       
    ## [411] "tot1.6 := a1.1*b1.6 + a1.2*b2.6 + a1.3*b3.6 + a1.4*b4.6 + a1.5*b5.6 + c1.6"                                                                                                       
    ## [412] "tot1.7 := a1.1*b1.7 + a1.2*b2.7 + a1.3*b3.7 + a1.4*b4.7 + a1.5*b5.7 + c1.7"                                                                                                       
    ## [413] "tot1.8 := a1.1*b1.8 + a1.2*b2.8 + a1.3*b3.8 + a1.4*b4.8 + a1.5*b5.8 + c1.8"                                                                                                       
    ## [414] "tot1.9 := a1.1*b1.9 + a1.2*b2.9 + a1.3*b3.9 + a1.4*b4.9 + a1.5*b5.9 + c1.9"                                                                                                       
    ## [415] "tot1.10 := a1.1*b1.10 + a1.2*b2.10 + a1.3*b3.10 + a1.4*b4.10 + a1.5*b5.10 + c1.10"                                                                                                
    ## [416] "tot1.11 := a1.1*b1.11 + a1.2*b2.11 + a1.3*b3.11 + a1.4*b4.11 + a1.5*b5.11 + c1.11"                                                                                                
    ## [417] "tot1.12 := a1.1*b1.12 + a1.2*b2.12 + a1.3*b3.12 + a1.4*b4.12 + a1.5*b5.12 + c1.12"                                                                                                
    ## [418] "tot1.13 := a1.1*b1.13 + a1.2*b2.13 + a1.3*b3.13 + a1.4*b4.13 + a1.5*b5.13 + c1.13"                                                                                                
    ## [419] "tot1.14 := a1.1*b1.14 + a1.2*b2.14 + a1.3*b3.14 + a1.4*b4.14 + a1.5*b5.14 + c1.14"                                                                                                
    ## [420] "tot1.15 := a1.1*b1.15 + a1.2*b2.15 + a1.3*b3.15 + a1.4*b4.15 + a1.5*b5.15 + c1.15"                                                                                                
    ## [421] "tot1.16 := a1.1*b1.16 + a1.2*b2.16 + a1.3*b3.16 + a1.4*b4.16 + a1.5*b5.16 + c1.16"                                                                                                
    ## [422] "tot2.1 := a2.1*b1.1 + a2.2*b2.1 + a2.3*b3.1 + a2.4*b4.1 + a2.5*b5.1 + c2.1"                                                                                                       
    ## [423] "tot2.2 := a2.1*b1.2 + a2.2*b2.2 + a2.3*b3.2 + a2.4*b4.2 + a2.5*b5.2 + c2.2"                                                                                                       
    ## [424] "tot2.3 := a2.1*b1.3 + a2.2*b2.3 + a2.3*b3.3 + a2.4*b4.3 + a2.5*b5.3 + c2.3"                                                                                                       
    ## [425] "tot2.4 := a2.1*b1.4 + a2.2*b2.4 + a2.3*b3.4 + a2.4*b4.4 + a2.5*b5.4 + c2.4"                                                                                                       
    ## [426] "tot2.5 := a2.1*b1.5 + a2.2*b2.5 + a2.3*b3.5 + a2.4*b4.5 + a2.5*b5.5 + c2.5"                                                                                                       
    ## [427] "tot2.6 := a2.1*b1.6 + a2.2*b2.6 + a2.3*b3.6 + a2.4*b4.6 + a2.5*b5.6 + c2.6"                                                                                                       
    ## [428] "tot2.7 := a2.1*b1.7 + a2.2*b2.7 + a2.3*b3.7 + a2.4*b4.7 + a2.5*b5.7 + c2.7"                                                                                                       
    ## [429] "tot2.8 := a2.1*b1.8 + a2.2*b2.8 + a2.3*b3.8 + a2.4*b4.8 + a2.5*b5.8 + c2.8"                                                                                                       
    ## [430] "tot2.9 := a2.1*b1.9 + a2.2*b2.9 + a2.3*b3.9 + a2.4*b4.9 + a2.5*b5.9 + c2.9"                                                                                                       
    ## [431] "tot2.10 := a2.1*b1.10 + a2.2*b2.10 + a2.3*b3.10 + a2.4*b4.10 + a2.5*b5.10 + c2.10"                                                                                                
    ## [432] "tot2.11 := a2.1*b1.11 + a2.2*b2.11 + a2.3*b3.11 + a2.4*b4.11 + a2.5*b5.11 + c2.11"                                                                                                
    ## [433] "tot2.12 := a2.1*b1.12 + a2.2*b2.12 + a2.3*b3.12 + a2.4*b4.12 + a2.5*b5.12 + c2.12"                                                                                                
    ## [434] "tot2.13 := a2.1*b1.13 + a2.2*b2.13 + a2.3*b3.13 + a2.4*b4.13 + a2.5*b5.13 + c2.13"                                                                                                
    ## [435] "tot2.14 := a2.1*b1.14 + a2.2*b2.14 + a2.3*b3.14 + a2.4*b4.14 + a2.5*b5.14 + c2.14"                                                                                                
    ## [436] "tot2.15 := a2.1*b1.15 + a2.2*b2.15 + a2.3*b3.15 + a2.4*b4.15 + a2.5*b5.15 + c2.15"                                                                                                
    ## [437] "tot2.16 := a2.1*b1.16 + a2.2*b2.16 + a2.3*b3.16 + a2.4*b4.16 + a2.5*b5.16 + c2.16"                                                                                                
    ## [438] "tot3.1 := a3.1*b1.1 + a3.2*b2.1 + a3.3*b3.1 + a3.4*b4.1 + a3.5*b5.1 + c3.1"                                                                                                       
    ## [439] "tot3.2 := a3.1*b1.2 + a3.2*b2.2 + a3.3*b3.2 + a3.4*b4.2 + a3.5*b5.2 + c3.2"                                                                                                       
    ## [440] "tot3.3 := a3.1*b1.3 + a3.2*b2.3 + a3.3*b3.3 + a3.4*b4.3 + a3.5*b5.3 + c3.3"                                                                                                       
    ## [441] "tot3.4 := a3.1*b1.4 + a3.2*b2.4 + a3.3*b3.4 + a3.4*b4.4 + a3.5*b5.4 + c3.4"                                                                                                       
    ## [442] "tot3.5 := a3.1*b1.5 + a3.2*b2.5 + a3.3*b3.5 + a3.4*b4.5 + a3.5*b5.5 + c3.5"                                                                                                       
    ## [443] "tot3.6 := a3.1*b1.6 + a3.2*b2.6 + a3.3*b3.6 + a3.4*b4.6 + a3.5*b5.6 + c3.6"                                                                                                       
    ## [444] "tot3.7 := a3.1*b1.7 + a3.2*b2.7 + a3.3*b3.7 + a3.4*b4.7 + a3.5*b5.7 + c3.7"                                                                                                       
    ## [445] "tot3.8 := a3.1*b1.8 + a3.2*b2.8 + a3.3*b3.8 + a3.4*b4.8 + a3.5*b5.8 + c3.8"                                                                                                       
    ## [446] "tot3.9 := a3.1*b1.9 + a3.2*b2.9 + a3.3*b3.9 + a3.4*b4.9 + a3.5*b5.9 + c3.9"                                                                                                       
    ## [447] "tot3.10 := a3.1*b1.10 + a3.2*b2.10 + a3.3*b3.10 + a3.4*b4.10 + a3.5*b5.10 + c3.10"                                                                                                
    ## [448] "tot3.11 := a3.1*b1.11 + a3.2*b2.11 + a3.3*b3.11 + a3.4*b4.11 + a3.5*b5.11 + c3.11"                                                                                                
    ## [449] "tot3.12 := a3.1*b1.12 + a3.2*b2.12 + a3.3*b3.12 + a3.4*b4.12 + a3.5*b5.12 + c3.12"                                                                                                
    ## [450] "tot3.13 := a3.1*b1.13 + a3.2*b2.13 + a3.3*b3.13 + a3.4*b4.13 + a3.5*b5.13 + c3.13"                                                                                                
    ## [451] "tot3.14 := a3.1*b1.14 + a3.2*b2.14 + a3.3*b3.14 + a3.4*b4.14 + a3.5*b5.14 + c3.14"                                                                                                
    ## [452] "tot3.15 := a3.1*b1.15 + a3.2*b2.15 + a3.3*b3.15 + a3.4*b4.15 + a3.5*b5.15 + c3.15"                                                                                                
    ## [453] "tot3.16 := a3.1*b1.16 + a3.2*b2.16 + a3.3*b3.16 + a3.4*b4.16 + a3.5*b5.16 + c3.16"                                                                                                
    ## [454] "tot4.1 := a4.1*b1.1 + a4.2*b2.1 + a4.3*b3.1 + a4.4*b4.1 + a4.5*b5.1 + c4.1"                                                                                                       
    ## [455] "tot4.2 := a4.1*b1.2 + a4.2*b2.2 + a4.3*b3.2 + a4.4*b4.2 + a4.5*b5.2 + c4.2"                                                                                                       
    ## [456] "tot4.3 := a4.1*b1.3 + a4.2*b2.3 + a4.3*b3.3 + a4.4*b4.3 + a4.5*b5.3 + c4.3"                                                                                                       
    ## [457] "tot4.4 := a4.1*b1.4 + a4.2*b2.4 + a4.3*b3.4 + a4.4*b4.4 + a4.5*b5.4 + c4.4"                                                                                                       
    ## [458] "tot4.5 := a4.1*b1.5 + a4.2*b2.5 + a4.3*b3.5 + a4.4*b4.5 + a4.5*b5.5 + c4.5"                                                                                                       
    ## [459] "tot4.6 := a4.1*b1.6 + a4.2*b2.6 + a4.3*b3.6 + a4.4*b4.6 + a4.5*b5.6 + c4.6"                                                                                                       
    ## [460] "tot4.7 := a4.1*b1.7 + a4.2*b2.7 + a4.3*b3.7 + a4.4*b4.7 + a4.5*b5.7 + c4.7"                                                                                                       
    ## [461] "tot4.8 := a4.1*b1.8 + a4.2*b2.8 + a4.3*b3.8 + a4.4*b4.8 + a4.5*b5.8 + c4.8"                                                                                                       
    ## [462] "tot4.9 := a4.1*b1.9 + a4.2*b2.9 + a4.3*b3.9 + a4.4*b4.9 + a4.5*b5.9 + c4.9"                                                                                                       
    ## [463] "tot4.10 := a4.1*b1.10 + a4.2*b2.10 + a4.3*b3.10 + a4.4*b4.10 + a4.5*b5.10 + c4.10"                                                                                                
    ## [464] "tot4.11 := a4.1*b1.11 + a4.2*b2.11 + a4.3*b3.11 + a4.4*b4.11 + a4.5*b5.11 + c4.11"                                                                                                
    ## [465] "tot4.12 := a4.1*b1.12 + a4.2*b2.12 + a4.3*b3.12 + a4.4*b4.12 + a4.5*b5.12 + c4.12"                                                                                                
    ## [466] "tot4.13 := a4.1*b1.13 + a4.2*b2.13 + a4.3*b3.13 + a4.4*b4.13 + a4.5*b5.13 + c4.13"                                                                                                
    ## [467] "tot4.14 := a4.1*b1.14 + a4.2*b2.14 + a4.3*b3.14 + a4.4*b4.14 + a4.5*b5.14 + c4.14"                                                                                                
    ## [468] "tot4.15 := a4.1*b1.15 + a4.2*b2.15 + a4.3*b3.15 + a4.4*b4.15 + a4.5*b5.15 + c4.15"                                                                                                
    ## [469] "tot4.16 := a4.1*b1.16 + a4.2*b2.16 + a4.3*b3.16 + a4.4*b4.16 + a4.5*b5.16 + c4.16"

### Apply model to imputed data sets

``` r
mult.pgs.fit <- sem.mi(model = mediation.model, data = data.imputed, ordered = c(paste0("O", 1:NO)), 
                        se = "robust", missing = "pairwise", parallel = "multicore", ncpus = 8)
```

### Extract and save parameters

``` r
mult.pgs.par = parameterEstimates.mi(mult.pgs.fit)
```

### Save standardised solution

``` r
mult.pgs.std = standardizedSolution.mi(mult.pgs.fit)
```

### Control for multiple comparisons

``` r
mult.pgs.df = adjust_multiple_df(mult.pgs.std)

sig.paths.std = mult.pgs.df %>% filter(pvalue.adj < 0.05)
```

## Single PGS model

### Depression

``` r
single.pgs.dep = single_model(data.lav, data.imputed, "dep.pgs", NO, NM, NCV)
```

### Schizophrenia

``` r
single.pgs.scz = single_model(data.lav, data.imputed, "scz.pgs", NO, NM, NCV)
```

### Bipolar Disoder

``` r
single.pgs.bd = single_model(data.lav, data.imputed, "bd.pgs", NO, NM, NCV)
```

### ADHD

``` r
single.pgs.adhd = single_model(data.lav, data.imputed, "adhd.pgs", NO, NM, NCV)
```

## Total scores model

### Create sum scores on imputed data

``` r
imp.tot <- data.frame(complete(data.imputed, include = TRUE, action = "long"))
imp.tot = imp.tot %>%
  rowwise() %>%
  mutate(tot.anx = sum(c(O1, O2, O16, O14, O15, O12, O13)),
         tot.dep = sum(c(O3, O4, O5, O6, O7, O8, O9, O10, O11)))
imp.tot <- as.mids(imp.tot)
```

### Demographics, means, sds on imputed data

``` r
imp.long = complete(imp.tot,action="long",include = FALSE)

columns_of_interest <- setdiff(names(imp.long), ".imp")  

pool_results <- list()

for (col_name in columns_of_interest) {
  pool_mean <- with(imp.long, by(imp.long, .imp, function(x) {
    c(mean = mean(x[[col_name]], na.rm = TRUE), sd = sd(x[[col_name]], na.rm = TRUE))
  }))
  
  pool_results[[col_name]] <- Reduce("+", pool_mean) / length(pool_mean)
}


pool_data <- list()


for (col_name in names(pool_results)) {
  
  pool_data[[col_name]] <- data.frame(
    Column = col_name, 
    Mean = pool_results[[col_name]][1],  # Mean
    SD = pool_results[[col_name]][2]     # Standard deviation
  )
}

pool_df <- do.call(rbind, pool_data)

pool_df_round = as.data.frame(lapply(pool_df[-1], round, digits = 2))
rownames(pool_df_round) = rownames(pool_df)
pool_df_round
model_names_df <- data.frame(Column = names(model_names), Description = model_names, stringsAsFactors = FALSE)
pool_df_round_names = merge(pool_df_round, model_names_df, by.x = "row.names", by.y = "Description") %>% select(Column, Mean, SD)

pool_df_round_names
```

### Create model

``` r
tot.model = make_multiple_model_cov(NO = 1, NE, NM, NCV)

tot.model.dep =  gsub("O1", "tot.dep", tot.model)
tot.model.anx =  gsub("O1", "tot.anx", tot.model)
```

### Apply model

``` r
tot.dep.fit <- sem.mi(model = tot.model.dep, data = imp.tot, se = "robust")

tot.anx.fit <- sem.mi(model = tot.model.anx, data = imp.tot, se = "robust")
```

### Save standardised solution

``` r
tot.dep.std = standardizedSolution.mi(tot.dep.fit)
tot.anx.std = standardizedSolution.mi(tot.anx.fit)
```

### Control for multiple comparisons

``` r
tot.dep.df = adjust_multiple_df(tot.dep.std)

tot.anx.df = adjust_multiple_df(tot.anx.std)
```

# Plot estimates

## Prepare results for plotting

``` r
# Create path for single and multiple pgs

multi.pgs = mult.pgs.df
multi.pgs$model = "multiple"
multi.pgs = create_path_label(multi.pgs, multi.pgs$label)


# Merge single PGS solutions
single.pgs.list = list("dep" = create_single_label(single.pgs.dep$adjusted.solution, single.pgs.dep$adjusted.solution$label, "dep.pgs"),
                     "adhd" = create_single_label(single.pgs.adhd$adjusted.solution, single.pgs.adhd$adjusted.solution$label, "adhd.pgs"),
                     "bd" = create_single_label(single.pgs.bd$adjusted.solution, single.pgs.bd$adjusted.solution$label, "bd.pgs"),
                     "scz" = create_single_label(single.pgs.scz$adjusted.solution, single.pgs.scz$adjusted.solution$label, "scz.pgs"))

df_list_named <- Map(function(df, name) {
  df$model <- name
  return(df)
}, single.pgs.list, names(single.pgs.list))

# Combine the dfs into one
single.pgs <- bind_rows(df_list_named)


# Merge single and multiple pgs

single.pgs = single.pgs %>% select(est.std, se, ci.lower, ci.upper, type, pvalue, pvalue.adj, type.label, model, prefix, suffix)
multi.pgs = multi.pgs %>% select(est.std, se, ci.lower, ci.upper, type, pvalue, pvalue.adj, type.label, model, prefix, suffix)

mediation.pgs = rbind(single.pgs, multi.pgs)

mediation.pgs$sig = ifelse(mediation.pgs$pvalue.adj < 0.05, 1, 0.2)
mediation.pgs$model = ifelse(mediation.pgs$model == "multiple", "multiple", "single")
```

## Plot results

### Total effects, single and multiple PGS models

``` r
(plot_ct = mediation.pgs %>%
  filter(type %in% c("t")) %>%
  ggplot(
    .,
    aes(
      y = factor(suffix),
      x = est.std,
      xmin = ci.lower,
      xmax = ci.upper,
      color = model,
      alpha = sig
    )
  ) +
  geom_point(size = 2, position = position_dodge(width = 0.8)) +
  geom_errorbarh(height = 0.6, position = position_dodge(width = 0.8)) +
  facet_wrap( ~ prefix, nrow = 1, labeller = labeller(
    prefix = c(
      `adhd.pgs` = "ADHD PGS",
      `bd.pgs` = "Bipolar PGS",
      `dep.pgs` = "Depression PGS",
      `scz.pgs` = "Schizophrenia PGS"
    )
  )) +
  scale_alpha_identity() +
  theme_classic() +
  ylab("") +
  xlab("Standardised estimate of total effect") +
  theme(legend.position = "bottom") +
  geom_vline(xintercept = 0, alpha = 0.2) +
  scale_color_manual(
    name = "Model",
    values = c("#E69F00", "black"),
    labels = c("Multiple PGS", "Single PGS")
  ) +
  scale_y_discrete(labels=c("anx" = "Anxious",
                            "wor" = "Can't control worry",
                            "psy" = "Psychomotor issues",
                            "dep" = "Depressed",
                            "ina" = "Feeling inadequate",
                            "tir" = "Tired",
                            "int" = "Lack of interest",
                            "app" = "Appetite problems",
                            "sui" = "Suicidal thoughts",
                            "con" = "Concentration issues",
                            "sle" = "Sleep issues",
                            "irr" = "Irritable",
                            "fore" = "Sense of foreboding",
                            "res" = "Restless",
                            "rel" = "Trouble relaxing",
                            "wor.t" = "Worrying too much")))
```

![](mediation_analysis_ukb17122024-114527_files/figure-gfm/unnamed-chunk-35-1.png)<!-- -->

### Total paths and unmediated paths comparisons

``` r
c_paths = mediation.pgs %>%
  filter(type == c("c") & model == "multiple") %>%
  select(est.std, type, type.label)

t_paths = mediation.pgs %>%
  filter(type == c("t") & model == "multiple")%>%
  select(est.std, type, type.label)


ct_paths = merge(c_paths, t_paths, by = "type.label") %>%
  select(type.label, est.std.x, est.std.y)
colnames(ct_paths) = c("type.label", "c", "t")

ct_paths$Difference = ct_paths$t - ct_paths$c


(
  deviation_plot = ggplot(ct_paths, aes(
    x = t, y = c, colour = Difference
  )) +
    geom_point(size = 2) +
    geom_abline() +
    ylab("Standardised estimate of unmediated paths") +
    xlab("Standardised estimate of total paths") +
    theme_minimal() +
    geom_label_repel(
      color = "black",
      size = 3,
      data = filter(ct_paths, Difference > 0.004),
      min.segment.length = 0,
      seed = 42,
      box.padding = 0.5 ,
      nudge_x = 0.009,
     aes(
        color = "black",
        label = c("ADHD PGS - Appetite",
                  "ADHD PGS - Lack of interest", 
                  "ADHD PGS - Psychomotor issues",
                  "ADHD PGS - Restlessness",
                  "ADHD PGS - Suicidal thoughts", 
                  "ADHD PGS - Tiredness")
      )
    )
)
```

![](mediation_analysis_ukb17122024-114527_files/figure-gfm/unnamed-chunk-36-1.png)<!-- -->

### PGS-mediator associations, single and multiple PGS models

``` r
(plot_a = mediation.pgs %>%
  filter(type %in% c("a")) %>%
  ggplot(
    .,
    aes(
      y = factor(suffix),
      x = est.std,
      xmin = ci.lower,
      xmax = ci.upper,
      color = model,
      alpha = sig
    )
  ) +
  geom_point(size = 2, position = position_dodge(width = 0.8)) +
  geom_errorbarh(height = 0.6, position = position_dodge(width = 0.8)) +
  facet_wrap( ~ prefix, nrow = 1, labeller = labeller(
    prefix = c(
      `adhd.pgs` = "ADHD PGS",
      `bd.pgs` = "Bipolar PGS",
      `dep.pgs` = "Depression PGS",
      `scz.pgs` = "Schizophrenia PGS"
    )
  )) +
  scale_alpha_identity() +
  theme_classic() +
   ylab("") +
  xlab("Standardised estimate") +
   scale_x_continuous(breaks = c(-0.075, -0.050, -0.025, 0, 0.025, 0.045)) +
  theme(legend.position = "bottom") +
  geom_vline(xintercept = 0, alpha = 0.2) +
  scale_color_manual(
    name = "Model",
    values = c("#3399CC", "black"),
    labels = c("Multiple PGS", "Single PGS")
  )  +
  scale_y_discrete(labels=c(`cing` = "Cingulate",
      `fusi` = "Fusiform",
      `hip` = "Hippocampus",
      `ins` = "Insula",
      `mofc` = "mOFC")) +
  theme(legend.position = "bottom", 
        axis.text.y=element_text(size=14), 
        strip.text.x = element_text(size = 14), 
        legend.title = element_text(size=14),
        legend.text = element_text(size = 12),
        axis.title.x = element_text(size = 14)))
```

![](mediation_analysis_ukb17122024-114527_files/figure-gfm/unnamed-chunk-37-1.png)<!-- -->

### Mediator-outcome associations

``` r
res.mi.cor <- micombine.cor(mi.res=data.imputed)
res.mi.pcor <-micombine.cor(mi.res=data.imputed, partial = ~CV1 + CV2)

cor.subset = res.mi.cor %>% filter(variable1 %in% c("O1", "O2", "O3", "O4", "O5", "O6", "O7", "O8", "O9", "O10", "O11", "O12", "O13", "O14", "O15", "O16")) %>%
  filter(variable2 %in% c("M1", "M2", 'M3', "M4", "M5")) %>%
  select(variable1, variable2, r, p, lower95, upper95)

cor.subset$pavlue.adj = p.adjust(cor.subset$p, method = "fdr")

pcor.subset = res.mi.pcor %>% filter(variable1 %in% c("O1", "O2", "O3", "O4", "O5", "O6", "O7", "O8", "O9", "O10", "O11", "O12", "O13", "O14", "O15", "O16")) %>%
  filter(variable2 %in% c("M1", "M2", 'M3', "M4", "M5")) %>%
  select(variable1, variable2, r, p, lower95, upper95)

pcor.subset$pavlue.adj = p.adjust(pcor.subset$p, method = "fdr")

cor.subset = cor.subset %>% mutate(variable1 = recode(variable1,
                                   "O1" = "anx",
                                   "O2" = "wor",
                                   "O3" = "psy",
                                   "O4" = "dep",
                                   "O5" = "ina",
                                   "O6" = "tir",
                                   "O7" = "int",
                                   "O8" = "app",
                                   "O9" = "sui",
                                   "O10" = "con",
                                   "O11" = "sle",
                                   "O12" = "irr",
                                   "O13" = "fore",
                                   "O14" = "res",
                                   "O15" = "rel",
                                   "O16" = "wor.t"),
                                   variable2 = recode(variable2,
                                    "M1" = "mofc",
                                    "M2" = "fusi",
                                    "M3" = "ins",
                                    "M4" = "hip",
                                    "M5" = "cing")) 

pcor.subset = pcor.subset %>% mutate(variable1 = recode(variable1,
                                   "O1" = "anx",
                                   "O2" = "wor",
                                   "O3" = "psy",
                                   "O4" = "dep",
                                   "O5" = "ina",
                                   "O6" = "tir",
                                   "O7" = "int",
                                   "O8" = "app",
                                   "O9" = "sui",
                                   "O10" = "con",
                                   "O11" = "sle",
                                   "O12" = "irr",
                                   "O13" = "fore",
                                   "O14" = "res",
                                   "O15" = "rel",
                                   "O16" = "wor.t"),
                                   variable2 = recode(variable2,
                                    "M1" = "mofc",
                                    "M2" = "fusi",
                                    "M3" = "ins",
                                    "M4" = "hip",
                                    "M5" = "cing")) %>% select(-p)



cor.subset$type.label = paste0(cor.subset$variable1, " - ", cor.subset$variable2)
cor.subset$model = "simple" 

pcor.subset$type.label = paste0(pcor.subset$variable1, " - ", pcor.subset$variable2)
pcor.subset$model = "partial (sex and age)" 


b_paths_multiple = mediation.pgs %>% filter(type == "b" & model == "multiple") %>% select(-c(type, se, sig))
colnames(b_paths_multiple) = c("r", "lower95", "upper95", "pvalue", "pavlue.adj", "type.label", "model", "variable2", "variable1")

b_paths_comp = rbind(cor.subset, pcor.subset, b_paths_multiple)
b_paths_comp$sig = ifelse(b_paths_comp$pavlue.adj < 0.05, 1, 0.3)
```

``` r
(b_comp_plot = b_paths_comp %>%
   filter(model != "partial (sex and age)") %>%
ggplot(.,
  aes(
    y = factor(variable1),
    x = r,
    xmin = lower95,
    xmax = upper95,
    color = model,
    alpha = sig
  )
) +
  facet_wrap( ~ variable2, nrow = 1, labeller = labeller(
    variable2 = c(
      `fusi` = "Fusiform",
      `ins` = "Insula",
      `hip` = "Hippocampus",
      `mofc` = "mOFC",
      `cing` = "Cingulate" 
    ))
  ) +
  geom_point(size = 2, position = position_dodge(width = 0.8)) +
  geom_errorbarh(height = 0.6, position = position_dodge(width = 0.8)) +
    theme_classic() +
  geom_vline(xintercept = 0, alpha = 0.2) +
  scale_y_discrete(labels=c("anx" = "Anxious",
                            "wor" = "Can't control worry",
                            "psy" = "Psychomotor issues",
                            "dep" = "Depressed",
                            "ina" = "Feeling inadequate",
                            "tir" = "Tired",
                            "int" = "Lack of interest",
                            "app" = "Appetite problems",
                            "sui" = "Suicidal thoughts",
                            "con" = "Concentration issues",
                            "sle" = "Sleep issues",
                            "irr" = "Irritable",
                            "fore" = "Sense of foreboding",
                            "res" = "Restless",
                            "rel" = "Trouble relaxing",
                            "wor.t" = "Worrying too much")) +
  scale_color_manual(
    name = "Estimate type",
    values = c("brown", "darkgreen"),
    labels = c("Multiple PGS mediation model", "Simple")
  ) +
  scale_alpha_identity() +
  xlab("Estimate") +
  ylab("") +
  theme(legend.position = "bottom", 
        axis.text.y=element_text(size=14), 
        strip.text.x = element_text(size = 14), 
        legend.title = element_text(size=14),
        legend.text = element_text(size = 12),
        axis.title.x = element_text(size = 14)))
```

![](mediation_analysis_ukb17122024-114527_files/figure-gfm/unnamed-chunk-39-1.png)<!-- -->

### Total mediation effects, single and multiple PGS models

``` r
(plot_mtot = mediation.pgs %>%
  filter(type == "m tot") %>%
  ggplot(
    .,
    aes(
      y = factor(suffix),
      x = est.std,
      xmin = ci.lower,
      xmax = ci.upper,
      color = model,
      alpha = sig
    )
  ) +
  geom_point(size = 2, position = position_dodge(width = 0.8)) +
  geom_errorbarh(height = 0.6, position = position_dodge(width = 0.8)) +
  facet_wrap( ~ prefix, nrow = 1, labeller = labeller(
    prefix = c(
      `adhd.pgs` = "ADHD PGS",
      `bd.pgs` = "Bipolar PGS",
      `dep.pgs` = "Depression PGS",
      `scz.pgs` = "Schizophrenia PGS"
    )
  )) +
  theme_classic() +
  xlab("Standardised estimate") +
   ylab("") +
  theme(legend.position = "bottom") +
  geom_vline(xintercept = 0, alpha = 0.2) +
  scale_color_manual(
    name = "Model",
    values = c("lightcoral", "black"),
    labels = c("Multiple PGS", "Single PGS")
  ) +
  scale_alpha_identity() +
  scale_y_discrete(labels=c("anx" = "Anxious",
                            "wor" = "Can't control worry",
                            "psy" = "Psychomotor issues",
                            "dep" = "Depressed",
                            "ina" = "Feeling inadequate",
                            "tir" = "Tired",
                            "int" = "Lack of interest",
                            "app" = "Appetite problems",
                            "sui" = "Suicidal thoughts",
                            "con" = "Concentration issues",
                            "sle" = "Sleep issues",
                            "irr" = "Irritable",
                            "fore" = "Sense of foreboding",
                            "res" = "Restless",
                            "rel" = "Trouble relaxing",
                            "wor.t" = "Worrying too much")))
```

![](mediation_analysis_ukb17122024-114527_files/figure-gfm/unnamed-chunk-40-1.png)<!-- -->

### Significant individual mediation effects, single and multiple PGS models

``` r
# No significant individual mediation effects

mediation.pgs %>%
  filter(type == "m" & pvalue.adj < 0.05) 
```

    ##  [1] est.std    se         ci.lower   ci.upper   type       pvalue.adj
    ##  [7] type.label model      prefix     suffix     sig       
    ## <0 rows> (or 0-length row.names)

### Sum scores results

``` r
tot.dep.pgs = create_path_label_tot(tot.dep.df, tot.dep.df$label)
tot.anx.pgs = create_path_label_tot(tot.anx.df, tot.anx.df$label)


tot.dep = tot.dep.pgs %>% filter(suffix == "total") %>% select(est.std, ci.lower, ci.upper, pvalue, pvalue.adj, prefix, suffix, type)
tot.anx = tot.anx.pgs %>% filter(suffix == "total") %>% select(est.std, ci.lower, ci.upper, pvalue, pvalue.adj, prefix, suffix, type)

tot.dep$suffix = "tot.dep"
tot.anx$suffix = "tot.anx"

tot.est = rbind(tot.dep, tot.anx)

tot.est = tot.est %>% mutate(prefix = recode(prefix,
                                   "mofc" = "mOFC",
                                   "fusi" = "Fusiform",
                                   "ins" = "Insula",
                                   "hip" = "Hippocampus",
                                   "cing" = "Cingulate",
                                   "dep.pgs" = "Depression PGS",
                                   "scz.pgs" = "Schizophrenia PGS",
                                   "bd.pgs" = "Bipolar disorder PGS",
                                   "adhd.pgs" = "ADHD PGS"),
                                   suffix = recode(suffix,
                                    "tot.dep" = "Total Depression",
                                    "tot.anx" = "Total Anxiety"))

tot.est$type.label = paste0(tot.est$prefix, " - ", tot.est$suffix)
tot.est$sig = ifelse(tot.est$pvalue.adj < 0.05, 1, 0.2)
```

``` r
tot_plot = tot.est %>%
  filter(type != "m tot") %>%
  mutate(type.label = fct_reorder(type.label, type)) %>% 
  ggplot(
    .,
    aes(
      y = type.label,
      x = est.std,
      xmin = ci.lower,
      xmax = ci.upper,
      color = factor(type),
      alpha = sig
    )
  ) +
  geom_point(size = 2, position = position_dodge(width = 0.8)) +
  geom_errorbarh(height = 0.6, position = position_dodge(width = 0.8)) +
  theme_classic() +
  xlab("Standardised effect estimate") +
  ylab("") +
  geom_vline(xintercept = 0, alpha = 0.5)+
  scale_color_manual(
    name = "Path type",
    values = c("brown", "orange", "darkgreen"),
    labels = c("PGS-Brain mediator", "Unmediated effect", "Total effect")
  ) +
  scale_alpha_identity()
tot_plot
```

![](mediation_analysis_ukb17122024-114527_files/figure-gfm/unnamed-chunk-43-1.png)<!-- -->

# Print estimates and package versions for draft

# Save files
