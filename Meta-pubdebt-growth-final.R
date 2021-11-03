rm(list = ls()) #clear list

#automatic installation of required packages
packages <- c("xlsx","calibrate","stargazer","sandwich","lmtest","getopt","CausalGAM","ggplot2","reshape2","xts",
              "lattice","gridExtra","gtable","plm","lfe","lmtest","car","tis","foreign","MASS","quantreg","ggrepel",
              "dplyr","stringr","datasets","rio","psych","systemfit","MatchIt","CRTgeeDR","eurostat","plyr","zoo","ggthemes",
              "robumeta","metafor","dplyr","clubSandwich","Hmisc","metafor","pracma","pkgs","broom","sjPlot", "here", "data.table", "pscore", "AER")
ipak(packages)

#load packages
library(xlsx) #Excel-Paket laden
library(calibrate) #Laden des Pakets, das f??r Datenbeschriftung n??tig ist
library (stargazer) #Laden des Pakets, mit dem R-Regressionsoutput in Latex-Tabellen ??bergef??hrt werden kann
library(sandwich)
library(lmtest)
library(getopt)
library(CausalGAM)
library(ggplot2)
library(reshape2)
library(xts)
library(lattice)
library(gridExtra)
library(gtable)
library(plm)
library(lfe)
library(lmtest)
library(car)
library(tis)
library(foreign)
library(MASS)
library(quantreg)
library(ggrepel)
library(dplyr)
library(stringr)
library(ggplot2)
library(datasets)
library(rio)
library(psych)
library(systemfit)
library(foreign)
library(MatchIt)
library(CRTgeeDR)
library(eurostat)
library(plyr)
library(zoo)
library(ggthemes)
library("robumeta")
library("metafor")
library("dplyr")
library(clubSandwich)
library(Hmisc)
library(metafor)
library(pracma)
library(broom)
library(sjPlot)
library(here)
library(data.table)
library(pscore)
library(AER)

#Load data
dat <- fread(here("public-debt-growth-coding-final.csv"))

#calculate the partial correlation coefficient
dat$PartialCorrelationCoefficient<- dat$TstatisticCorrected / (sqrt((dat$TstatisticCorrected^2)+dat$DegreesofFreedom))

dat$PartialCorrelationCoefficient <- dat$PartialCorrelationCoefficient*dat$Transform
dat$Tstatistic <- dat$TstatisticCorrected*dat$Transform

#calculate standard error of the partial correlation coefficient
dat$StandardErrorPartialCorrelation <- sqrt((1-(dat$PartialCorrelationCoefficient)^2)/dat$DegreesofFreedom)

#Precision of partial correlation
dat$PrecSE <- 1 / dat$StandardErrorPartialCorrelation
#Inverse standard error of partial correlation
dat$InverseSE <- 1 / dat$StandardErrorCorrected
#Inverse standard error of Standardised correlation coefficient
dat$InverseSECorrelationCoefficientCorrected <- 1 / dat$StandardErrorCorrected
#Variance of partial correlation
dat$Variance <- dat$StandardErrorPartialCorrelation^2
#Variance of Standardised correlation coefficient
dat$VarianceCorrelationCoefficientCorrected <- dat$StandardErrorCorrected^2
#PrecVariance partial correlation
dat$PrecVariance <- 1 / dat$Variance
#PrecVariance Standardised correlation coefficient
dat$PrecVarianceCorrelationCoefficientCorrected <- 1 / dat$VarianceCorrelationCoefficientCorrected
dat$PrecVarianceCorrelationCoefficientCorrected
#Proxy for standard error
dat$proxySE <- 1 / (sqrt(dat$Observations))

#adjust publication year
dat$PublicationYear <- dat$YearofPublication - mean(dat$YearofPublication)

#taking logs
#dat$YearofPublication <- log(dat$YearofPublication)
dat$NumberofCountries <- log(dat$NumberofCountries)
dat$Citations <- log(dat$Citations)

#variable type
dat$StartYear <- as.numeric(as.integer(dat$StartYear))
dat$EndYear <- as.numeric(as.integer(dat$EndYear))

#constructing additional variables
dat$LengthTimeSpan <- (dat$EndYear - dat$StartYear) + 1
dat$LengthTimeSpan <- log(dat$LengthTimeSpan)
dat$MeanYearData_c<- (dat$StartYear+dat$EndYear)/2
dat$MeanYearData<-dat$MeanYearData_c - mean(dat$MeanYearData_c)
dat$InstrumentSE <- 1 / (sqrt(dat$DegreesofFreedom))
dat$InstrumentVariance <- dat$InstrumentSE^2
dat$PrecInstrumentVariance <- 1 / dat$InstrumentVariance
#normalize journal impact factor
dat$JournalImpactFactor <- as.numeric(dat$JournalImpactFactor)
dat$MaxImpactFactor <- max(dat$JournalImpactFactor)
dat$MaxImpactFactor <- as.numeric(dat$MaxImpactFactor)
dat$NormalizedImpactFactor <- dat$JournalImpactFactor / max(dat$JournalImpactFactor)

#data subsample
dat_CorrelationCoefficientCorrected <- subset(dat, YNCorrected %in% c('1'))

#winsorising at the 2nd and 98th percentiles
dat_2p98p <- dat_CorrelationCoefficientCorrected

dat_2p98p$StandardErrorCorrected <- winsorizor(dat_CorrelationCoefficientCorrected$StandardErrorCorrected, c(0.02), na.rm=TRUE)
dat_2p98p$CorrelationCoefficientCorrected <- winsorizor(dat_CorrelationCoefficientCorrected$CorrelationCoefficientCorrected, c(0.02), na.rm=TRUE)
dat_2p98p$InstrumentSE <- winsorizor(dat_CorrelationCoefficientCorrected$InstrumentSE, c(0.02), na.rm=TRUE)

#calculate standard errors corrected after precision was sed
dat_2p98p$InverseSECorrected <- 1 / dat_2p98p$StandardErrorCorrected

#Variance of corrected correlation coefficient
dat_2p98p$VarianceSECorrected <- dat_2p98p$StandardErrorCorrected^2
#PrecVariance corrected
dat_2p98p$PrecVarianceSECorrected <- 1 / dat_2p98p$VarianceSECorrected

#preferred estimates only
dat_Preferred <- subset(dat_2p98p, Preferred %in% c('1'))

#Descriptive statistics
#all-set

#Standardised correlation coefficient
#unweighted average
uwa <- sum(dat_2p98p$CorrelationCoefficientCorrected, na.rm=TRUE) / 566
uwa
reguwa <- lm(CorrelationCoefficientCorrected~1, data=dat_2p98p)
confint(reguwa, level=0.95)
summary(reguwa)

coef_test(reguwa, vcov = "CR0", 
          cluster = dat_2p98p$paperid, test = "naive-t")

#confidence interval
confint(reguwa, level=0.95)

#descriptive stats
median(dat_2p98p$CorrelationCoefficientCorrected)
min(dat_2p98p$CorrelationCoefficientCorrected)
max(dat_2p98p$CorrelationCoefficientCorrected)
sd(dat_2p98p$CorrelationCoefficientCorrected)

#(precision-weighted) average
regwa <- lm(CorrelationCoefficientCorrected~1, data=dat_2p98p, weights=PrecVarianceCorrelationCoefficientCorrected)
summary(regwa)
#confidence interval
confint(regwa, level=0.95)

#Results on publication selection bias (table 1)

#FAT-PET test
#base with precision-weights
pubbias_1_CorrelationCoefficientCorrected <- lm(CorrelationCoefficientCorrected ~ StandardErrorCorrected, weights=PrecVarianceCorrelationCoefficientCorrected, data=dat_2p98p)
summary(pubbias_1_CorrelationCoefficientCorrected)

coef_test(pubbias_1_CorrelationCoefficientCorrected, vcov = "CR0", 
          cluster = dat_2p98p$paperid, test = "naive-t")


#FAT-PET test
#Median estimates only
pubbias_1_CorrelationCoefficientCorrected_median <- lm(CoeffMedian ~ StandardErrorCorrected, weights=PrecVarianceCorrelationCoefficientCorrected, data=dat_2p98p)
summary(pubbias_1_CorrelationCoefficientCorrected_median)

coef_test(pubbias_1_CorrelationCoefficientCorrected_median, vcov = "CR0", 
          cluster = dat_2p98p$paperid, test = "naive-t")

#Preferred estimates only
pubbias_1_CorrelationCoefficientCorrected_preferred <- lm(CorrelationCoefficientCorrected ~ StandardErrorCorrected, weights=PrecVarianceCorrelationCoefficientCorrected, data=dat_Preferred)
summary(pubbias_1_CorrelationCoefficientCorrected_preferred)

coef_test(pubbias_1_CorrelationCoefficientCorrected_preferred, vcov = "CR0", 
          cluster = dat_Preferred$paperid, test = "naive-t")

#IV
pubbias_1_CorrelationCoefficientCorrected_IV <- AER::ivreg(CorrelationCoefficientCorrected ~ StandardErrorCorrected | InstrumentSE, weights=1/InstrumentSE, data=dat_2p98p)
summary(pubbias_1_CorrelationCoefficientCorrected_IV)

coef_test(pubbias_1_CorrelationCoefficientCorrected_IV, vcov = "CR0", 
          cluster = dat_2p98p$paperid, test = "naive-t")

#partial correlation coefficient
pubbias_1_pcc <- lm(PartialCorrelationCoefficient ~ StandardErrorPartialCorrelation, weights=PrecVariance, data=dat_2p98p)
summary(pubbias_1_pcc)

coef_test(pubbias_1_pcc, vcov = "CR0", 
          cluster = dat_2p98p$paperid, test = "naive-t")

#Partial correlation
#unweighted average
uwa_partial <- sum(dat_2p98p$PartialCorrelationCoefficient, na.rm=TRUE) / 566
uwa_partial

regwa_partial <- lm(PartialCorrelationCoefficient~1, data=dat_2p98p, weights=PrecVariance)
summary(regwa_partial)
confint(regwa_partial, level=0.95)

#correct t-values and standard errors for stargazer table
ses_pubbias_1_CorrelationCoefficientCorrected <- list(coef_test(pubbias_1_CorrelationCoefficientCorrected, vcov = "CR0", 
                                                                    cluster = dat_2p98p$paperid, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals_pubbias_1_CorrelationCoefficientCorrected <- list(coef_test(pubbias_1_CorrelationCoefficientCorrected, vcov = "CR0", 
                                                                      cluster = dat_2p98p$paperid, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

pvals_pubbias_1_CorrelationCoefficientCorrected <- list(coef_test(pubbias_1_CorrelationCoefficientCorrected , vcov = "CR0", 
                                                                      cluster = dat_2p98p$paperid, test = "naive-t")[,4]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.


ses_pubbias_1_CorrelationCoefficientCorrected_median <- list(coef_test(pubbias_1_CorrelationCoefficientCorrected_median, vcov = "CR0", 
                                                                cluster = dat_2p98p$paperid, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals_pubbias_1_CorrelationCoefficientCorrected_median <- list(coef_test(pubbias_1_CorrelationCoefficientCorrected_median, vcov = "CR0", 
                                                                  cluster = dat_2p98p$paperid, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

pvals_pubbias_1_CorrelationCoefficientCorrected_median <- list(coef_test(pubbias_1_CorrelationCoefficientCorrected_median, vcov = "CR0", 
                                                                  cluster = dat_2p98p$paperid, test = "naive-t")[,4]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.


ses_pubbias_1_CorrelationCoefficientCorrected_preferred <- list(coef_test(pubbias_1_CorrelationCoefficientCorrected_preferred, vcov = "CR0", 
                                                                       cluster = dat_Preferred$paperid, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals_pubbias_1_CorrelationCoefficientCorrected_preferred <- list(coef_test(pubbias_1_CorrelationCoefficientCorrected_preferred, vcov = "CR0", 
                                                                         cluster = dat_Preferred$paperid, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

pvals_pubbias_1_CorrelationCoefficientCorrected_preferred <- list(coef_test(pubbias_1_CorrelationCoefficientCorrected_preferred, vcov = "CR0", 
                                                                         cluster = dat_Preferred$paperid, test = "naive-t")[,4]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.


ses_pubbias_1_CorrelationCoefficientCorrected_IV <- list(coef_test(pubbias_1_CorrelationCoefficientCorrected_IV, vcov = "CR0", 
                                                                          cluster = dat_2p98p$paperid, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals_pubbias_1_CorrelationCoefficientCorrected_IV <- list(coef_test(pubbias_1_CorrelationCoefficientCorrected_IV, vcov = "CR0", 
                                                                            cluster = dat_2p98p$paperid, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

pvals_pubbias_1_CorrelationCoefficientCorrected_IV <- list(coef_test(pubbias_1_CorrelationCoefficientCorrected_IV, vcov = "CR0", 
                                                                            cluster = dat_2p98p$paperid, test = "naive-t")[,4]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

ses_pubbias_1_CorrelationCoefficientCorrected_pcc <- list(coef_test(pubbias_1_pcc, vcov = "CR0", 
                                                                       cluster = dat_2p98p$paperid, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals_pubbias_1_CorrelationCoefficientCorrected_pcc <- list(coef_test(pubbias_1_pcc, vcov = "CR0", 
                                                                         cluster = dat_2p98p$paperid, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

pvals_pubbias_1_CorrelationCoefficientCorrected_pcc <- list(coef_test(pubbias_1_pcc, vcov = "CR0", 
                                                                         cluster = dat_2p98p$paperid, test = "naive-t")[,4]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

#stargazer table linear publication bias tests (table 1)
stargazer(pubbias_1_CorrelationCoefficientCorrected, pubbias_1_CorrelationCoefficientCorrected_median, pubbias_1_CorrelationCoefficientCorrected_preferred, pubbias_1_pcc, t=list(unlist(tvals_pubbias_1_CorrelationCoefficientCorrected), unlist(tvals_pubbias_1_CorrelationCoefficientCorrected_median), unlist(tvals_pubbias_1_CorrelationCoefficientCorrected_preferred), unlist(tvals_pubbias_1_CorrelationCoefficientCorrected_pcc)), se=list(unlist(ses_pubbias_1_CorrelationCoefficientCorrected), unlist(ses_pubbias_1_CorrelationCoefficientCorrected_median), unlist(ses_pubbias_1_CorrelationCoefficientCorrected_preferred), unlist(ses_pubbias_1_CorrelationCoefficientCorrected_pcc)), p=list(unlist(pvals_pubbias_1_CorrelationCoefficientCorrected), unlist(pvals_pubbias_1_CorrelationCoefficientCorrected_median), unlist(pvals_pubbias_1_CorrelationCoefficientCorrected_preferred), unlist(pvals_pubbias_1_CorrelationCoefficientCorrected_pcc)))

#Funnel plot (figure 3)

#Standardised coefficient
#Funnel plot
plot_funnel_CorrelationCoefficientCorrected_2p98p <- ggplot(data=dat_2p98p,
                                                            aes(x=CorrelationCoefficientCorrected, y=InverseSECorrected)) +
  geom_point(size=1.5, color="blue") +
  xlab("Standardised coefficient") +
  ylab("Inverse of standard error (precision)") +
  ggtitle("Standardised coefficients of public debt-growth estimates (n=566)")+
  theme(title=element_text(size=10, face='bold'))+
  theme(panel.background = element_rect(fill = "white")) +
  theme(legend.position="bottom")+
  theme(legend.title=element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.1),
        axis.line.y = element_line(color="black", size = 0.1)) +
  geom_vline(xintercept=-0, colour="black", linetype=2)+
  geom_vline(xintercept=-0.01420261, colour="black", linetype=1)+
  theme(legend.text = element_text(colour="black", size = 4))+
  annotate("text", x = -0.06, y = 800, size=3.5, label = "Unweighted mean: -0.014") +
  annotate("text", x = 0.02, y = 800, size=3.5, label = "Zero effect") +
  theme(axis.text.x=element_text(size=11))+
  theme(axis.title.x=element_text(size=11)) +
  theme(axis.text.y=element_text(size=11))+
  theme(axis.title.y=element_text(size=11))
plot_funnel_CorrelationCoefficientCorrected_2p98p

#Partial correlation
plot_funnel_PartialCorrelation <- ggplot(data=dat_2p98p,
                                         aes(x=PartialCorrelationCoefficient, y=PrecSE)) +
  geom_point(size=1.5, color="red") +
  xlab("Partial correlation coefficient") +
  ylab("Inverse of standard error (precision)") +
  ggtitle("Partial correlations of public debt-growth estimates")+
  theme(title=element_text(size=10, face='bold'))+
  theme(panel.background = element_rect(fill = "white")) +
  theme(legend.position="bottom")+
  theme(legend.title=element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.1),
        axis.line.y = element_line(color="black", size = 0.1)) +
  geom_vline(xintercept=-0.1375918, colour="black", linetype=1)+
  geom_vline(xintercept=0, colour="black", linetype=2)+
  annotate("text", x = 0.1, y = 55, size=3.5, label = "Zero effect") +
  theme(legend.text = element_text(colour="black", size = 4))+
  annotate("text", x = -0.4, y = 55, size=3.5, label = "Unweighted mean: -0.142") +
  theme(axis.text.x=element_text(size=11))+
  theme(axis.title.x=element_text(size=11)) +
  theme(axis.text.y=element_text(size=11))+
  theme(axis.title.y=element_text(size=11))
plot_funnel_PartialCorrelation

###

#Meta-regression results for table 2
#baseline
MRA_var_est_CorrelationCoefficientCorrected_1 <- lm(CorrelationCoefficientCorrected ~ StandardErrorCorrected + DevelopingCountriesOnly + MixofCountries + LongRunExplicit + ShortRunExplicit + TacklingEndogeneity + NormalizedImpactFactor + Citations, weights=PrecVarianceCorrelationCoefficientCorrected, data=dat_2p98p)
summary(MRA_var_est_CorrelationCoefficientCorrected_1)

coef_test(MRA_var_est_CorrelationCoefficientCorrected_1, vcov = "CR0", 
          cluster = dat_2p98p$paperid, test = "naive-t")

#+ Estimation + Data
MRA_var_est_CorrelationCoefficientCorrected_4 <- lm(CorrelationCoefficientCorrected ~ StandardErrorCorrected + DevelopingCountriesOnly + MixofCountries + LongRunExplicit + ShortRunExplicit + TacklingEndogeneity + NormalizedImpactFactor + Citations + GrowthPerCapita + LaggedPublicDebt + CrossSection + MeanYearData, weights=PrecVarianceCorrelationCoefficientCorrected, data=dat_2p98p)
summary(MRA_var_est_CorrelationCoefficientCorrected_4)

coef_test(MRA_var_est_CorrelationCoefficientCorrected_4, vcov = "CR0", 
          cluster = dat_2p98p$paperid, test = "naive-t")

#+ Publication characteristics
MRA_var_est_CorrelationCoefficientCorrected_5 <- lm(CorrelationCoefficientCorrected ~ StandardErrorCorrected + DevelopingCountriesOnly + MixofCountries + LongRunExplicit + ShortRunExplicit + TacklingEndogeneity + PublicationYear + NormalizedImpactFactor + Citations, weights=PrecVarianceCorrelationCoefficientCorrected, data=dat_2p98p)
summary(MRA_var_est_CorrelationCoefficientCorrected_5)

coef_test(MRA_var_est_CorrelationCoefficientCorrected_5, vcov = "CR0", 
          cluster = dat_2p98p$paperid, test = "naive-t")

#+ additional control variables in the underlying model specifications
MRA_var_est_CorrelationCoefficientCorrected_6 <- lm(CorrelationCoefficientCorrected ~ StandardErrorCorrected + DevelopingCountriesOnly + MixofCountries + LongRunExplicit + ShortRunExplicit + TacklingEndogeneity + NormalizedImpactFactor + Citations + Investment + Inflation, weights=PrecVarianceCorrelationCoefficientCorrected, data=dat_2p98p)
summary(MRA_var_est_CorrelationCoefficientCorrected_6)

coef_test(MRA_var_est_CorrelationCoefficientCorrected_6, vcov = "CR0", 
          cluster = dat_2p98p$paperid, test = "naive-t")

#Preferred
MRA_var_est_CorrelationCoefficientCorrected_Preferred <- lm(CorrelationCoefficientCorrected ~ StandardErrorCorrected + DevelopingCountriesOnly + MixofCountries + LongRunExplicit + ShortRunExplicit + TacklingEndogeneity + NormalizedImpactFactor + Citations, weights=PrecVarianceCorrelationCoefficientCorrected, data=dat_Preferred)
summary(MRA_var_est_CorrelationCoefficientCorrected_Preferred)

coef_test(MRA_var_est_CorrelationCoefficientCorrected_Preferred, vcov = "CR0", 
          cluster = dat_Preferred$paperid, test = "naive-t")

#partial correlation
MRA_var_est_PartialCorrelation <- lm(PartialCorrelationCoefficient ~ StandardErrorPartialCorrelation + DevelopingCountriesOnly + MixofCountries + LongRunExplicit + ShortRunExplicit + TacklingEndogeneity + NormalizedImpactFactor + Citations + GrowthPerCapita + LaggedPublicDebt + CrossSection + MeanYearData, weights=PrecVariance, data=dat_2p98p)
summary(MRA_var_est_PartialCorrelation)

coef_test(MRA_var_est_PartialCorrelation, vcov = "CR0", 
          cluster = dat_2p98p$paperid, test = "naive-t")

#correct t-values and standard errors for stargazer table
ses_MRA_var_est_CorrelationCoefficientCorrected_1 <- list(coef_test(MRA_var_est_CorrelationCoefficientCorrected_1, vcov = "CR0", 
                                                         cluster = dat_2p98p$paperid, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals_MRA_var_est_CorrelationCoefficientCorrected_1 <- list(coef_test(MRA_var_est_CorrelationCoefficientCorrected_1, vcov = "CR0", 
                                                           cluster = dat_2p98p$paperid, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

pvals_MRA_var_est_CorrelationCoefficientCorrected_1 <- list(coef_test(MRA_var_est_CorrelationCoefficientCorrected_1, vcov = "CR0", 
                                                           cluster = dat_2p98p$paperid, test = "naive-t")[,4]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

ses_MRA_var_est_CorrelationCoefficientCorrected_4 <- list(coef_test(MRA_var_est_CorrelationCoefficientCorrected_4, vcov = "CR0", 
                                                                        cluster = dat_2p98p$paperid, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals_MRA_var_est_CorrelationCoefficientCorrected_4 <- list(coef_test(MRA_var_est_CorrelationCoefficientCorrected_4, vcov = "CR0", 
                                                                          cluster = dat_2p98p$paperid, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

pvals_MRA_var_est_CorrelationCoefficientCorrected_4 <- list(coef_test(MRA_var_est_CorrelationCoefficientCorrected_4, vcov = "CR0", 
                                                                          cluster = dat_2p98p$paperid, test = "naive-t")[,4]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

ses_MRA_var_est_CorrelationCoefficientCorrected_5 <- list(coef_test(MRA_var_est_CorrelationCoefficientCorrected_5, vcov = "CR0", 
                                                                        cluster = dat_2p98p$paperid, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals_MRA_var_est_CorrelationCoefficientCorrected_5 <- list(coef_test(MRA_var_est_CorrelationCoefficientCorrected_5, vcov = "CR0", 
                                                                          cluster = dat_2p98p$paperid, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

pvals_MRA_var_est_CorrelationCoefficientCorrected_5 <- list(coef_test(MRA_var_est_CorrelationCoefficientCorrected_5, vcov = "CR0", 
                                                                          cluster = dat_2p98p$paperid, test = "naive-t")[,4]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

ses_MRA_var_est_CorrelationCoefficientCorrected_6 <- list(coef_test(MRA_var_est_CorrelationCoefficientCorrected_6, vcov = "CR0", 
                                                                        cluster = dat_2p98p$paperid, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals_MRA_var_est_CorrelationCoefficientCorrected_6 <- list(coef_test(MRA_var_est_CorrelationCoefficientCorrected_6, vcov = "CR0", 
                                                                          cluster = dat_2p98p$paperid, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

pvals_MRA_var_est_CorrelationCoefficientCorrected_6 <- list(coef_test(MRA_var_est_CorrelationCoefficientCorrected_6, vcov = "CR0", 
                                                                          cluster = dat_2p98p$paperid, test = "naive-t")[,4]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.


ses_MRA_var_est_CorrelationCoefficientCorrected_Preferred <- list(coef_test(MRA_var_est_CorrelationCoefficientCorrected_Preferred, vcov = "CR0", 
                                                                    cluster = dat_Preferred$paperid, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals_MRA_var_est_CorrelationCoefficientCorrected_Preferred <- list(coef_test(MRA_var_est_CorrelationCoefficientCorrected_Preferred, vcov = "CR0", 
                                                                      cluster = dat_Preferred$paperid, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

pvals_MRA_var_est_CorrelationCoefficientCorrected_Preferred <- list(coef_test(MRA_var_est_CorrelationCoefficientCorrected_Preferred, vcov = "CR0", 
                                                                      cluster = dat_Preferred$paperid, test = "naive-t")[,4]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.


ses_MRA_var_est_PartialCorrelation <- list(coef_test(MRA_var_est_PartialCorrelation, vcov = "CR0", 
                                                     cluster = dat_2p98p$paperid, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals_MRA_var_est_PartialCorrelation <- list(coef_test(MRA_var_est_PartialCorrelation, vcov = "CR0", 
                                                       cluster = dat_2p98p$paperid, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

pvals_MRA_var_est_PartialCorrelation <- list(coef_test(MRA_var_est_PartialCorrelation, vcov = "CR0", 
                                                           cluster = dat_2p98p$paperid, test = "naive-t")[,4]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

#-regression table (table 2)
stargazer(MRA_var_est_CorrelationCoefficientCorrected_1, MRA_var_est_CorrelationCoefficientCorrected_4, MRA_var_est_CorrelationCoefficientCorrected_6, MRA_var_est_CorrelationCoefficientCorrected_Preferred, MRA_var_est_PartialCorrelation, t=list(unlist(tvals_MRA_var_est_CorrelationCoefficientCorrected_1), unlist(tvals_MRA_var_est_CorrelationCoefficientCorrected_4), unlist(tvals_MRA_var_est_CorrelationCoefficientCorrected_6), unlist(tvals_MRA_var_est_CorrelationCoefficientCorrected_Preferred), unlist(tvals_MRA_var_est_PartialCorrelation)), se=list(unlist(ses_MRA_var_est_CorrelationCoefficientCorrected_1), unlist(ses_MRA_var_est_CorrelationCoefficientCorrected_4), unlist(ses_MRA_var_est_CorrelationCoefficientCorrected_6), unlist(ses_MRA_var_est_CorrelationCoefficientCorrected_Preferred), unlist(ses_MRA_var_est_PartialCorrelation)), p=list(unlist(pvals_MRA_var_est_CorrelationCoefficientCorrected_1), unlist(pvals_MRA_var_est_CorrelationCoefficientCorrected_4), unlist(pvals_MRA_var_est_CorrelationCoefficientCorrected_6), unlist(pvals_MRA_var_est_CorrelationCoefficientCorrected_Preferred), unlist(pvals_MRA_var_est_PartialCorrelation)))

#Threshold effects
dat_threshold <- fread(here("Threshold_pdebt_growth_final.csv"))
dat_threshold$InverseObservations <- 1/dat_threshold$Observations
dat_threshold$PublicationYear <- dat_threshold$YearofPublication - mean(dat_threshold$YearofPublication)
#normalize impact factor
dat_threshold$JournalImpactFactor <- as.numeric(dat_threshold$JournalImpactFactor)
dat_threshold$MaxImpactFactor <- max(dat_threshold$JournalImpactFactor)
dat_threshold$MaxImpactFactor <- as.numeric(dat_threshold$MaxImpactFactor)
dat_threshold$NormalizedImpactFactor <- dat_threshold$JournalImpactFactor / max(dat_threshold$JournalImpactFactor)
dat_threshold$Citations <- log(dat_threshold$Citations)
dat_threshold$MeanYearData_c<- (dat_threshold$StartYear+dat_threshold$EndYear)/2
dat_threshold$MeanYearData<-dat_threshold$MeanYearData_c - mean(dat_threshold$MeanYearData_c)

#winsorize
dat_threshold$Observations <- winsorizor(dat_threshold$Observations, c(0.02), na.rm=TRUE)
dat_threshold$Threshold <- winsorizor(dat_threshold$Threshold, c(0.02), na.rm=TRUE)

dat_threshold$SEProxy <- 1/sqrt(dat_threshold$Observations)
dat_threshold$PrecisionProxy <- sqrt(dat_threshold$Observations)

#descriptive statistics
median(dat_threshold$Threshold)
mean(dat_threshold$Threshold)
min(dat_threshold$Threshold)
max(dat_threshold$Threshold)
sd(dat_threshold$Threshold)

#preferred
dat_threshold_Preferred <- subset(dat_threshold, Preferred %in% c('1'))

#Standardised correlation coefficient
#unweighted average
uwa <- sum(dat_threshold$Threshold, na.rm=TRUE) / 260
uwa
reguwa <- lm(Threshold~1, data=dat_threshold)
summary(reguwa)

#(precision-weighted) average
numerator <- dat_threshold$Threshold*dat_threshold$PrecisionProxy
wa <- sum(numerator, na.rm=TRUE)/sum(dat_threshold$PrecisionProxy, na.rm=TRUE)
wa

regwa <- lm(Threshold~1, data=dat_threshold, weights=PrecisionProxy)
summary(regwa)
confint(regwa, level=0.95)

#FAT-PET test
pubbias_1_Threshold <- lm(Threshold ~ SEProxy, weights=PrecisionProxy, data=dat_threshold)
summary(pubbias_1_Threshold)

coef_test(pubbias_1_Threshold, vcov = "CR0", 
          cluster = dat_threshold$paperid, test = "naive-t")

confint(pubbias_1_Threshold, level=0.95)

#
plot_funnel_threshold <- ggplot(data=dat_threshold,
  aes(x=Threshold, y=PrecisionProxy)) +
  geom_point(size=1.5, color="blue") +
  xlab("Threshold (Public debt in % of GDP)") +
  ylab("Square root of sample size\n (as a proxy for precision)") +
  ggtitle("Funnel plot of public debt-growth threshold estimates (n=260)")+
  theme(title=element_text(size=10, face='bold'))+
  theme(panel.background = element_rect(fill = "white")) +
  theme(legend.position="bottom")+
  theme(legend.title=element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.1),
        axis.line.y = element_line(color="black", size = 0.1)) +
  geom_vline(xintercept=90.0, colour="black", linetype=2)+
  geom_vline(xintercept=59.79082, colour="black", linetype=1)+
  theme(legend.text = element_text(colour="black", size = 4))+
  theme(axis.text.x=element_text(size=13))+
  theme(axis.title.x=element_text(size=13))
plot_funnel_threshold

#save the Standardised correlation coefficient plot
#filename <- "figures/funnel_threshold.jpg"
#ggsave(filename, plot = plot_funnel_threshold, width = 10, height = 6)

#Meta-regressions (table 3)

MRA_threshold_1 <- lm(Threshold ~ SEProxy + DevelopingCountriesOnly + MixofCountries + SquaredTerm + TacklingEndogeneity + NormalizedImpactFactor + Citations, weights=PrecisionProxy, data=dat_threshold)
summary(MRA_threshold_1)

coef_test(MRA_threshold_1, vcov = "CR0", 
          cluster = dat_threshold$paperid, test = "naive-t")

MRA_threshold_2 <- lm(Threshold ~ SEProxy + DevelopingCountriesOnly + MixofCountries + SquaredTerm + TacklingEndogeneity + NormalizedImpactFactor + Citations + GrowthPerCapita + LaggedPublicDebt + MeanYearData, weights=PrecisionProxy, data=dat_threshold)
summary(MRA_threshold_2)

coef_test(MRA_threshold_2, vcov = "CR0", 
          cluster = dat_threshold$paperid, test = "naive-t")

MRA_threshold_3 <- lm(Threshold ~ SEProxy + DevelopingCountriesOnly + MixofCountries + SquaredTerm + TacklingEndogeneity + NormalizedImpactFactor + Citations + Investment + Inflation, weights=PrecisionProxy, data=dat_threshold)
summary(MRA_threshold_3)

coef_test(MRA_threshold_3, vcov = "CR0", 
          cluster = dat_threshold$paperid, test = "naive-t")

MRA_threshold_4 <- lm(Threshold ~ SEProxy + DevelopingCountriesOnly + MixofCountries + SquaredTerm + TacklingEndogeneity + NormalizedImpactFactor + Citations, weights=PrecisionProxy, data=dat_threshold_Preferred)
summary(MRA_threshold_4)

coef_test(MRA_threshold_4, vcov = "CR0", 
          cluster = dat_threshold_Preferred$paperid, test = "naive-t")

#correct t-values and standard errors for stargazer table
ses_MRA_threshold_1 <- list(coef_test(MRA_threshold_1, vcov = "CR0", 
                                                                    cluster = dat_threshold$paperid, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals_MRA_threshold_1 <- list(coef_test(MRA_threshold_1, vcov = "CR0", 
                                                                      cluster = dat_threshold$paperid, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

pvals_MRA_threshold_1 <- list(coef_test(MRA_threshold_1, vcov = "CR0", 
                                                                      cluster = dat_threshold$paperid, test = "naive-t")[,4]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

ses_MRA_threshold_2 <- list(coef_test(MRA_threshold_2, vcov = "CR0", 
                                      cluster = dat_threshold$paperid, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals_MRA_threshold_2 <- list(coef_test(MRA_threshold_2, vcov = "CR0", 
                                        cluster = dat_threshold$paperid, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

pvals_MRA_threshold_2 <- list(coef_test(MRA_threshold_2, vcov = "CR0", 
                                        cluster = dat_threshold$paperid, test = "naive-t")[,4]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

ses_MRA_threshold_3 <- list(coef_test(MRA_threshold_3, vcov = "CR0", 
                                      cluster = dat_threshold$paperid, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals_MRA_threshold_3 <- list(coef_test(MRA_threshold_3, vcov = "CR0", 
                                        cluster = dat_threshold$paperid, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

pvals_MRA_threshold_3 <- list(coef_test(MRA_threshold_3, vcov = "CR0", 
                                        cluster = dat_threshold$paperid, test = "naive-t")[,4]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

ses_MRA_threshold_4 <- list(coef_test(MRA_threshold_4, vcov = "CR0", 
                                      cluster = dat_threshold_Preferred$paperid, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals_MRA_threshold_4 <- list(coef_test(MRA_threshold_4, vcov = "CR0", 
                                        cluster = dat_threshold_Preferred$paperid, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

pvals_MRA_threshold_4 <- list(coef_test(MRA_threshold_4, vcov = "CR0", 
                                        cluster = dat_threshold_Preferred$paperid, test = "naive-t")[,4]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

#threshold regression results (table 3)
stargazer(MRA_threshold_1, MRA_threshold_2, MRA_threshold_3, MRA_threshold_4, t=list(unlist(tvals_MRA_threshold_1), unlist(tvals_MRA_threshold_2), unlist(tvals_MRA_threshold_3), unlist(tvals_MRA_threshold_4)), se=list(unlist(ses_MRA_threshold_1), unlist(ses_MRA_threshold_2), unlist(ses_MRA_threshold_3), unlist(ses_MRA_threshold_4)), p=list(unlist(pvals_MRA_threshold_1), unlist(pvals_MRA_threshold_2), unlist(pvals_MRA_threshold_3), unlist(pvals_MRA_threshold_4)))