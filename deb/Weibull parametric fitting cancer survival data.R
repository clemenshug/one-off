#Setting directory to R file location
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

#Required packages. Un-comment lines below if packages are not already installed
# install.packages("fitdistrplus")
# install.packages("ggplot2")
# install.packages("survival")
# install.packages("pracma")
# install.packages("tidyr")
# install.packages("dplyr")
# install.packages("pammtools")
# install.packages("utils")
# install.packages("stats")


#Loading required packages
library(fitdistrplus)
library(ggplot2)
library(survival)
library(pracma)
library(tidyr)
library(dplyr)
library(pammtools)
library(utils)
library(stats)


#Importing data. Can replace here with name of CSV file with trial IPD
TrialData <- read.csv("Checkmate057_1A.csv")
#TrialData <- read.csv("20050181_2A.csv")

#Extracting data corresponding to a single treatment arm of a clinical trial
#Can set "TreatmentArm" to name of therapy of interest from IPD csv
TreatmentArm<-unique(TrialData$Arm)[[1]]
IPDData<- TrialData[TrialData$Arm==TreatmentArm,]
CensoredIndex<-which(IPDData$Event==0)
DataInput= data.frame('right'=IPDData$Time, 'left'=IPDData$Time)
DataInput[CensoredIndex,1]=NA
Phase3Patients=length(IPDData$Time)

#Nonparametric fitting to trial data
KMfit <- survfit(Surv(Time,Event)
                 ~ 1, data = IPDData)
GroundTruthMedian<- as.numeric(summary(KMfit)$table[7])
GroundTruthYearSurvival= round(as.numeric(summary(KMfit, times = 12)$surv)*100, digits=3)
KMLabel= paste("Phase 3 data; n=", toString(Phase3Patients) )
KMDF= data.frame(X=KMfit$time, Y=KMfit$surv, category=replicate (length(KMfit$time), KMLabel))

#Parametric fitting to trial data
BestFit <- fitdistcens(DataInput, "weibull")
shape=as.numeric(BestFit$estimate[1])
scale=as.numeric(BestFit$estimate[2])
WeibullS= function(shape, scale,x) {
  return( exp(-(x/scale)^shape))
}
Time= seq(0,max(IPDData$Time),max(IPDData$Time)/(length(KMfit$time)-1))
WSurvival= WeibullS(shape, scale,Time)

#Assessing the quality of the Weibull fit using a Weibull plot and computing an R^2 value
CensoredTimes= IPDData[which(IPDData$Event==0),][,1]
ProgressionTimes= IPDData[which(IPDData$Event==1),][,1]
MiddleProgressionTimes= movavg(ProgressionTimes,2,type="s")[2:length(ProgressionTimes)]
SurvivalEstimate= as.numeric(summary(KMfit, times = MiddleProgressionTimes)$surv)
CorrectedS = log(-log(SurvivalEstimate))
rsq <- function (x, y) cor(x, y) ^ 2
R2=rsq(log(MiddleProgressionTimes), CorrectedS)
R2plot= round(R2, digits=3)
paste("R^2 value for Weibull fit: ", R2plot)

#Plotting the original and fit data
WDF= data.frame(X=Time, Y=WSurvival,category=replicate (length(KMfit$time), "Weibull Fit"))
WeibullLabel= expression(paste( R^2, "= ", R2plot))
DataPlot=rbind(WDF, KMDF)
DataPlotBest= ggplot(DataPlot, aes(x=X, y=100*Y, colour=category))+
  theme_bw()+
  labs(x="Time", y= "% Survival")+
  geom_step()+
  scale_color_manual(values=c('black','magenta'))+
  theme(text = element_text(size = 20))  +
  theme(legend.title=element_blank())+
  ylim(0,100)+
  ggtitle(bquote(""~R^2== .(R2plot)) )
print(DataPlotBest)

#Simulating Phase 2 trial
SimulatedPhase2=sample_n(IPDData,20, replace=FALSE)

#Computing 95% nonparametric CI:
KMfit <- survfit(Surv(Time,Event)
                 ~ 1, data = SimulatedPhase2)
Time= seq(0,max(SimulatedPhase2$Time),max(SimulatedPhase2$Time)/(length(SimulatedPhase2$Time)-1))
CILower= as.numeric(summary(KMfit, times = Time)$lower)
CIUpper= as.numeric(summary(KMfit, times = Time)$upper)
NP= as.numeric(summary(KMfit, times = Time)$surv)

#Computing 95% parametric CI:
CensoredIndex<-which(SimulatedPhase2$Event==0)
CensoredTimes= SimulatedPhase2[which(SimulatedPhase2$Event==0),][,1]
ProgressionTimes= SimulatedPhase2[which(SimulatedPhase2$Event==1),][,1]
DataInput= data.frame('right'=SimulatedPhase2$Time, 'left'=SimulatedPhase2$Time)
DataInput[CensoredIndex,1]=NA
BestFit <- fitdistcens(DataInput, "weibull")
shape=as.numeric(BestFit$estimate[1])
scale=as.numeric(BestFit$estimate[2])
WSurvival= WeibullS(shape, scale,Time)
WeibullLogCDF= function(shape, scale,x) {
  return( log(exp(-(x/scale)^shape)))
}
WeibullLogPDF= function(shape, scale,x) {
  return( log(((exp(-(x/scale)^shape))*shape*(x/scale)^(-1+shape))/scale))
}
alphas= seq(0.1,4,0.025)
betas= seq(1,scale*2,scale*2/150)
Likelihoods=c()
Parameters=c()
for (j in 1:length(alphas)) {
  for (k in 1:length(betas)) {
    Likelihood= sum(WeibullLogCDF(alphas[j],betas[k], CensoredTimes))+ sum(WeibullLogPDF(alphas[j],betas[k], ProgressionTimes))
    Likelihoods= c(Likelihoods,Likelihood )
    Label=paste(alphas[j], betas[k])
    Parameters= c(Parameters,Label )
  }}
SampleSize=1000
MaxLikelihood=max(Likelihoods)
AllLikelihoods= data.frame(parameters=Parameters, score=exp(-MaxLikelihood+Likelihoods))
ParameterSamples=sample(x=AllLikelihoods$parameters, size=SampleSize, replace=TRUE, prob=AllLikelihoods$score)
MinimumRelativeLikelihoodToSample = 1/SampleSize
AllLikelihoodsSelected= AllLikelihoods %>%
  filter(score >MinimumRelativeLikelihoodToSample)
ParameterSamples=sample(x=AllLikelihoodsSelected$parameters, size=SampleSize, replace=TRUE, prob=AllLikelihoodsSelected$score)
SampledAlphas=c()
SampledBetas=c()
for (j in 1:length(ParameterSamples)) {
SampleAlpha= as.numeric(strsplit(ParameterSamples[j], "\\s+")[[1]][1])
SampleBeta= as.numeric(strsplit(ParameterSamples[j], "\\s+")[[1]][2])
SampledAlphas= c(SampledAlphas,SampleAlpha)
SampledBetas= c(SampledBetas,SampleBeta)
}
CIValues= data.frame(SampledAlphas, SampledBetas)
ALPHA = 0.05
Wlower = c()
Wupper = c()
for(p in 1:length(Time)){
    TimeExp = Time[p]
    WCI= CIValues %>%
      mutate(Estimate=WeibullS(SampledAlphas, SampledBetas,TimeExp))
    Upper = quantile(WCI$Estimate, 1 - ALPHA/2)
    Lower= quantile(WCI$Estimate, ALPHA/2)
    Wlower= c(Wlower, as.numeric(Lower))
    Wupper= c(Wupper,as.numeric(Upper) )
}
CI_DF=data.frame(Time)

#Plotting subsampled data, with parametric and nonparametric confidence intervals
DataPlotCI_All= ggplot(CI_DF, aes(x=Time))+
  theme_bw()+
  labs(x="Time", y= "% Survival")+
  geom_step(aes(y = CILower*100, color="gray"))+
  geom_step(aes(y = CIUpper*100), color="gray") +
  geom_step(aes(y = NP*100, color="black"), size=1) +
  geom_stepribbon(aes(ymin=CILower*100,ymax=CIUpper*100), fill="gray", alpha=0.5) +
  geom_line(aes(y = Wlower*100, color="magenta"))+
  geom_line(aes(y = Wupper*100), color="magenta") +
  geom_ribbon(aes(ymin=Wlower*100,ymax=Wupper*100), fill="magenta", alpha=0.3) +
  scale_color_identity(name = "",
                       breaks = c("black", "gray", "magenta"),
                       labels = c("Simulated Phase 2 Data (n=20)", "Nonparametric CI", "Parametric CI"),
                       guide = "legend")+
  theme(text = element_text(size = 20))+
  scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(lim= c(0,100),expand = c(0, 0))
print(DataPlotCI_All)

#Comparing 12-month confidence interval
Endpoint= 12 #Change here for endpoint of interest
WCI= CIValues %>%
  mutate(Estimate=WeibullS(SampledAlphas, SampledBetas,Endpoint))
WUpper12 = quantile(WCI$Estimate, 1 - ALPHA/2)
WLower12= quantile(WCI$Estimate, ALPHA/2)
CILower12= as.numeric(summary(KMfit, times = Endpoint)$lower)
CIUpper12= as.numeric(summary(KMfit, times = Endpoint)$upper)

paste("Observed Phase 3 12-month survival (n=",Phase3Patients,"): ", GroundTruthYearSurvival, "%", sep = "")
paste(c("Parametric 95% CI at 12-months (n=20): ", round(as.numeric(WLower12), digits=3)*100, "% to ", round(as.numeric(WUpper12), digits=3)*100, "%"), collapse="")
paste(c("Nonparametric 95% CI at 12-months (n=20): ", round(CILower12, digits=3)*100, "% to ", round(CIUpper12, digits=3)*100, "%"), collapse="")

#Comparing median confidence intervals
Endpoint=0.5
WeibullInverseCDF= function(shape, scale,x){
  return(scale*(-log(1-x))^(1/shape))
}
MedianCI= CIValues %>%
  mutate(Estimate=WeibullInverseCDF(SampledAlphas, SampledBetas,Endpoint))
WLowerMedian=quantile(MedianCI$Estimate, 0.025)
WUpperMedian= quantile(MedianCI$Estimate, 0.975)
MedianLower=summary(KMfit)$table[8]
MedianUpper=summary(KMfit)$table[9]


paste("Observed Phase 3 median survival (n=",Phase3Patients,"): ", GroundTruthMedian, " months", sep = "")
paste(c("Parametric median 95% CI (n=20): ", round(WLowerMedian, digits=3), " months to ", round(WUpperMedian, digits=3), " months"), collapse="")
paste(c("Nonparametric median 95% CI (n=20): ", round(MedianLower, digits=3), " months to ", round(MedianUpper, digits=3), " months"), collapse="")
