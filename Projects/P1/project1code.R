##### set work directory and load dataset ##### 
setwd("/home/xmy/STA 101/Projects/P1")
HospFull<-read.csv("HospFull.csv", header = TRUE)
head(HospFull, n = 3)

##### load packages #####
library("ggplot2")
library("leaps")
library("MPV")

##### define functions #####
Partial.R2 = function(small.model,big.model){
  SSE1 = sum(small.model$residuals^2)
  SSE2 = sum(big.model$residuals^2)
  PR2 = (SSE1 - SSE2)/SSE1
  return(PR2)
}
All.Criteria = function(the.model){
  p = length(the.model$coefficients)
  n = length(the.model$residuals)
  the.BIC = BIC(the.model)
  the.LL = logLik(the.model)
  the.AIC = AIC(the.model)
  the.PRESS = PRESS(the.model)
  the.R2adj = summary(the.model)$adj.r.squared
  # the.CP = summary(the.model)$cp
  the.results = c(the.LL,p,n,the.AIC,the.BIC,the.PRESS,the.R2adj)
  names(the.results) = c("LL","p","n","AIC","BIC","PRESS","R2adj")
  return(the.results)
}


##### correlation #####
cor(HospFull$Length, HospFull$Infect)
cor(HospFull$Culture, HospFull$Infect)
cor(HospFull$Bed, HospFull$Infect)


##### Infect summary #####
summary(HospFull$Infect)
# grouped by MedSchool
aggregate(Infect ~ MedSchool, data = HospFull, summary)
# grouped by Region
aggregate(Infect ~ Region, data = HospFull, summary)

# plot(HospFull)
##### boxplots of Infect #####
require(ggplot2)
# boxplot grouped by MedSchool
ppi = 600
# Calculate the height and width (in pixels) for a 4x3-inch image at 600 ppi
png("group_boxplot_medschool.png", width=6*ppi, height=4*ppi, res=ppi)
ggplot(HospFull, aes(y=Infect, x = MedSchool))+ theme_gray() + geom_boxplot() + ylab("Probability of acquiring infection in hospital") + 
  xlab("category of MedSchool")+ coord_flip() 
#ggtitle("Boxplot of Infect grouped by Medchool") 
dev.off()

# boxplot grouped by Region
png("group_boxplot_region.png", width=6*ppi, height=4*ppi, res=ppi)
ggplot(HospFull, aes(y=Infect, x = Region))+ theme_gray() + geom_boxplot() + ylab("Probability of acquiring infection in hospital") + 
  xlab("category of Geographical region")+ coord_flip() 
#ggtitle("Boxplot of Infect grouped by Region") 
dev.off()


##### scatter plots of Infect #####
# scatter plot of Infect vs. Length
png("scatter_plot_length.png", width=6*ppi, height=4*ppi, res = ppi)
qplot(HospFull$Length, HospFull$Infect, data = HospFull) +xlab("length of stay") + ylab("probability of acquiring infection")
dev.off()

# scatter plot of Infect vs. Culture
png("scatter_plot_culture.png",  width=6*ppi, height=4*ppi, res = ppi)
qplot(HospFull$Culture, HospFull$Infect, data = HospFull) +xlab("culture/patients * 100") + ylab("probability of acquiring infection")
dev.off()

# scatter plot of Infect vs. Bed
png("scatter_plot_bed.png", width=6*ppi, height=4*ppi, res = ppi)
qplot(HospFull$Bed, HospFull$Infect, data = HospFull) +xlab("number of beds") + ylab("probability of acquiring infection")
dev.off()



##### remove outliers according to plots #####
# cover HospFull
the.original = HospFull
HospFull=HospFull[-which(HospFull$Length>15),]
HospFull=HospFull[-which(HospFull$Culture>60),]
HospFull=HospFull[-which(HospFull$MedSchool=="Y" & HospFull$Infect > 7),]
HospFull=HospFull[-which(HospFull$MedSchool=="N" & HospFull$Infect > 7),]
HospFull=HospFull[-which(HospFull$Region=="W" & HospFull$Infect < 3),]
HospFull=HospFull[-which(HospFull$Region=="NC" & HospFull$Infect < 2),]
length(the.original$Infect)
length(HospFull$Infect)
the.ratio = (length(the.original$Infect)-length(HospFull$Infect))/length(the.original$Infect)
the.ratio

##### subset models of Infect~. #####
# rename dataset for convenience
names(HospFull) = c("X1","Y","X2","X3","X4","X5")
full.model = lm(Y~ X1 + X2 + X3 + X4 + X5,data = HospFull)
round(full.model$coefficients,4)
bic.model = lm(Y~X1+X2+X5, data = HospFull)
round(bic.model$coefficients, 4)
all.models = c("Y~1","Y~X1","Y~X2","Y~X3","Y~X4","Y~X5",
               "Y~X1+X2","Y~X1+X3","Y~X1+X4","Y~X1+X5","Y~X2+X3","Y~X2+X4","Y~X2+X5","Y~X3+X4","Y~X3+X5","Y~X4+X5",
               "Y~X1+X2+X3","Y~X1+X2+X4","Y~X1+X2+X5","Y~X1+X3+X4","Y~X1+X3+X5","Y~X1+X4+X5","Y~X2+X3+X4","Y~X2+X3+X5","Y~X2+X4+X5","Y~X3+X4+X5",
               "Y~X1+X2+X3+X4","Y~X1+X2+X3+X5","Y~X1+X2+X4+X5","Y~X1+X3+X4+X5","Y~X2+X3+X4+X5",
               "Y~X1+X2+X3+X4+X5")
Infect.all.model.crit = t(sapply(all.models,function(M){
  current.model = lm(M,data = HospFull)
  All.Criteria(current.model)
}))
Infect.all.model.crit
Infect.all.model.crit = data.frame(Infect.all.model.crit)
# find the model with lowest BIC
Infect.all.model.crit[which(Infect.all.model.crit$BIC == min(Infect.all.model.crit[,5])),]
# find the model with lowest AIC
Infect.all.model.crit[which(Infect.all.model.crit$AIC == min(Infect.all.model.crit[,4])),]


##### anova analysis of X4 #####
summary(full.model)
summary(bic.model)
alpha = 0.05
the.CIs = confint(full.model,level = 1-alpha)
round(the.CIs, 4)
# drop X4
smaller.model = lm(Y~X1+X2+X3+X5, data = HospFull)
anova.small = anova(smaller.model)
larger.model = lm(Y~X1+X2+X3+X4+X5, data = HospFull)
anova.large = anova(larger.model)
anova(smaller.model,larger.model)


##### anova analysis of X3 #####
smaller.model = lm(Y~X1+X2+X5, data = HospFull)
anova.small = anova(smaller.model)
larger.model = lm(Y~X1+X2+X3+X5, data = HospFull)
anova.large = anova(larger.model)
anova(smaller.model,larger.model)
##### partial r2 of X3 #####
partial.R2=Partial.R2(smaller.model, larger.model)
partial.R2


##### considering interaction terms #####
# interaction term between X1 and X5
final.model = lm(Y~X1+X2+X5, data = HospFull)
final.model
X1.interation.model = lm(Y~X1+X2+X5+X1*X5, data = HospFull)
summary(X1.interation.model)
confint(X1.interation.model,level = 1-alpha)
anova(final.model,X1.interation.model)
partial.R2=Partial.R2(final.model,X1.interation.model)
partial.R2
# interaction term between X2 and X5
X2.interation.model = lm(Y~X1+X2+X5+X2*X5, data = HospFull)
X2.interation.model
summary(X2.interation.model)
confint(X2.interation.model,level = 1-alpha)
anova(final.model,X2.interation.model)
partial.R2=Partial.R2(final.model,X2.interation.model)
partial.R2


##### diagnose of model #####
final.model = lm(Y~X1+X2+X5, data = HospFull)
final.model
HospFull$ei = final.model$residuals
HospFull$yhat = final.model$fitted.values
## nomality
# qqplot
png("qqplot.png", width=6*ppi, height=4*ppi, res = ppi)
qqnorm(final.model$residuals)
qqline(final.model$residuals)
dev.off()
# S-W test
the.SWtest = shapiro.test(final.model$residuals)
the.SWtest

## constant variance
# ei-yi plot
png("scatter_plot_constant_variance.png", width=6*ppi, height=4*ppi, res = ppi)
qplot(yhat, ei, data = HospFull) +
  xlab("Fitted Values") + ylab("Errors") + geom_hline(yintercept = 0,col = "purple")
dev.off()
# F-K test
HospFull$ei = final.model$residuals
Group = rep("Lower",nrow(HospFull))
Group[HospFull$Y < median(HospFull$Y)] = "Upper"
Group = as.factor(Group)
HospFull$Group = Group
the.FKtest= fligner.test(HospFull$ei, HospFull$Group)
the.FKtest

## outliers
# cook's distance
cutoff = 0.10
CD = cooks.distance(final.model)
HospFull$CD = cooks.distance(final.model)
HospFull[which(HospFull$CD>cutoff),] 
# no outliers
png("cooks_distance.png", width=6*ppi, height=4*ppi, res = ppi)
plot(CD,ylab = "Cook's Distance")
abline(h = cutoff,color = "purple")
dev.off()

SR = stdres(final.model)
HospFull$SR = SR
cutoff= 3
png("standardized_error.png", width=6*ppi, height=4*ppi, res = ppi)
ggplot(HospFull,aes(x = SR))+geom_histogram(binwidth = 0.5,color = "black",fill = "white")+ xlab("standardized error")
dev.off()
SR[which(abs(SR) > cutoff)] 

##### final model #####
final.model
R2 = summary(final.model)$r.squared
R2

##### predict estimated values of Y #####
alpha = 0.05
x.star = data.frame(X1 = 8, X2 = 14, X5 = "W")
predict(final.model, x.star, interval = "confidence", level = 1-alpha)
predict(final.model, x.star, interval = "prediction", level = 1-alpha)


