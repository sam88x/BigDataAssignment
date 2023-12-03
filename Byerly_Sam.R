library(bootstrap)
library(car)
library(corrplot)
library(dplyr)
library(forcats)
library(ggplot2)
library(gvlma)
library(leaps)
library(multcomp)
library(mvtnorm)

fat = read.csv("bigData.csv")

head(fat)

fat[39,"Height"] = 69.5
rownames(fat) = fat[,1]
fat = fat[,-1]

pairs(fat,main = "Scatterplot Matrix")
cor(fat)
corrplot(cor(fat))

leaps <-regsubsets(Density ~ Age + Weight*Height + Neck + Chest + Abdomen + Hip + Thigh + Knee + Biceps + Forearm + Wrist, data=fat)
plot(leaps, scale="adjr2")

fit = lm(Density ~ Age + Weight*Height + Neck + Chest + Abdomen + Hip + Thigh + Knee + Biceps + Forearm + Wrist, data=fat)
summary(fit)

step(fit, direction="both")

fit2= lm(Density ~ Weight*Height + Neck + Abdomen + Hip + Thigh + Biceps, data=fat)
summary(fit2)

summary(powerTransform(fat$Density))
summary(powerTransform(fat$Weight))
summary(powerTransform(fat$Height))
summary(powerTransform(fat$Neck))
summary(powerTransform(fat$Abdomen))
summary(powerTransform(fat$Hip))
summary(powerTransform(fat$Thigh))
summary(powerTransform(fat$Biceps))

fit3 = lm(Density ~ Weight*Height + Neck + Abdomen + Hip + Thigh + Biceps,data=fat)
summary(fit3)
par(mfrow=c(2,2))
par(mar=c(2,2,2,2))
plot(fit3)

fat$sqrtDen = sqrt(fat$Density)
fat$logDen = log(fat$Density)
fit31 = lm(logDen ~ Weight*Height + Neck + Abdomen + Hip + Thigh + Biceps,data=fat)
summary(fit31)
plot(fit31)

std = rstandard(fit3)
hat = hatvalues(fit3)
keep = which(abs(std) < 3 & hat < length(fit3$coefficients)*3/235)
fatNoOut = fat[keep,]

gvlmaStat = gvlma(fit3)
summary(gvlmaStat)

library(car)
qqPlot(fit3, simulate=TRUE,
       id=list(method="identify"), main="Q-Q Plot")


durbinWatsonTest(fit3)
durbinWatsonTest(fit3, simulate=F)

ncvTest(fit3)
spreadLevelPlot(fit3)

vif(fit3) 
vif(fit3) > 10

influencePlot(fit, id="noteworthy", main="Influence Plot",
              sub="Circle size is proportional to Cook's distance")

summary(fat$Age)
for(i in 1:235){
  ifelse(fat$Age[i] < 30,fat$young[i] <- "<30",
         ifelse(fat$Age[i] < 40,fat$young[i] <- "30-39",
                ifelse(fat$Age[i] < 50,fat$young[i] <- "40-49",
                       ifelse(fat$Age[i] < 60,fat$young[i] <- "50-59",
                              ifelse(fat$Age[i] < 70,fat$young[i] <- "60-69",fat$young[i] <- ">70")))))
}

fat$young = factor(fat$young,levels=ordered)  
aov.fit <- aov(Density ~ young, data=fat)                                  
summary(aov.fit)

plotdata <- fat %>%
  group_by(young) %>%
  summarize(n = n(),
            mean = mean(Density),
            sd = sd(Density),
            ci = qt(0.975, df = n - 1) * sd / sqrt(n))
plotdata

library(forcats)
ordered = c("<30","30-39","40-49","50-59","60-69",">70")
ggplot(plotdata, 
       aes(x = young, y = mean, group = 1)) +
  geom_point(size = 3, color="red") +
  geom_line(linetype="dashed", color="darkgrey") +
  geom_errorbar(aes(ymin = mean - ci, 
                    ymax = mean + ci), 
                width = .2) +
  theme_bw() +
  xlab("Age") +
  ylab("Density")

library(multcomp)
tuk <- glht(aov.fit, linfct=mcp(young="Tukey")) 
summary(tuk)

labels1 <- cld(tuk, level=.05)$mcletters$Letters
labels2 <- paste(names(labels1), "\n", labels1)

ggplot(data=aov.fit$model, aes(x=young, y=Density)) +
  scale_x_discrete(breaks=names(labels1), labels=labels2) +
  geom_boxplot(fill="lightgrey") +
  theme_bw() +
  labs(x="Age",
       title="Distribution of Body Density by Age",
       subtitle="Groups without overlapping letters differ signifcantly (p < .05)")

aov.fit2 <- aov(Adiposity ~ young, data=fat)                                  
summary(aov.fit2)

plotdata <- fat %>%
  group_by(young) %>%
  summarize(n = n(),
            mean = mean(Adiposity),
            sd = sd(Adiposity),
            ci = qt(0.975, df = n - 1) * sd / sqrt(n))
plotdata

library(forcats)
ordered = c("<30","30-39","40-49","50-59","60-69",">70")
ggplot(plotdata, 
       aes(x = fct_inorder(ordered), y = mean, group = 1)) +
  geom_point(size = 3, color="red") +
  geom_line(linetype="dashed", color="darkgrey") +
  geom_errorbar(aes(ymin = mean - ci, 
                    ymax = mean + ci), 
                width = .2) +
  theme_bw() +
  xlab("Age") +
  ylab("Adiposity")

library(multcomp)
tuk <- glht(aov.fit2, linfct=mcp(young="Tukey")) 
summary(tuk)

labels1 <- cld(tuk, level=.05)$mcletters$Letters
labels2 <- paste(names(labels1), "\n", labels1)

ggplot(data=aov.fit2$model, aes(x=young, y=Adiposity)) +
  scale_x_discrete(breaks=names(labels1), labels=labels2) +
  geom_boxplot(fill="lightgrey") +
  theme_bw() +
  labs(x="Age",
       title="Distribution of Adiposity by Age",
       subtitle="Groups without overlapping letters differ signifcantly (p < .05)")

relweights <- function(fit,...){
  R <- cor(fit$model)
  nvar <- ncol(R)
  rxx <- R[2:nvar, 2:nvar]
  rxy <- R[2:nvar, 1]
  svd <- eigen(rxx)
  evec <- svd$vectors
  ev <- svd$values
  delta <- diag(sqrt(ev))
  lambda <- evec %*% delta %*% t(evec)
  lambdasq <- lambda ^ 2
  beta <- solve(lambda) %*% rxy
  rsquare <- colSums(beta ^ 2)
  rawwgt <- lambdasq %*% beta ^ 2
  import <- (rawwgt / rsquare) * 100
  import <- as.data.frame(import)
  row.names(import) <- names(fit$model[2:nvar])
  names(import) <- "Weights"
  import <- import[order(import),1, drop=FALSE]
  dotchart(import$Weights, labels=row.names(import),
           xlab="% of R-Square", pch=19,
           main="Relative Importance of Predictor Variables",
           sub=paste("Total R-Square=", round(rsquare, digits=3)),
           ...)
  return(import)
}


relweights(fit3, col="dodgerblue")

for(i in 1:235){
  ifelse(fat$Adiposity[i] < 18,fat$BMI[i] <- "Underweight",
         ifelse(fat$Adiposity[i] < 25,fat$BMI[i] <- "Normal Weight",
                ifelse(fat$Adiposity[i] < 30,fat$BMI[i] <- "Overweight",fat$BMI[i] <- "Obese")))
}

aov.fat3 = aov(Density ~ BMI*young,data=fat)
summary(aov.fat3)
fat$BMI <- factor(fat$BMI, levels = c("Normal Weight","Overweight","Obese"))

plotdata <- fat %>%
  group_by(young,BMI) %>%
  summarise(n=n(), mean=mean(Density), sd=sd(Density),
            ci = qt(0.975, df = n - 1) * sd / sqrt(n))
plotdata

pd <- position_dodge(0.2)
ggplot(plotdata, 
       aes(x = young, y = mean, 
           group=BMI, 
           color=BMI, 
           linetype=BMI)) +
  geom_point(size = 2, 
             position=pd) +
  geom_line(position=pd) +
  geom_errorbar(aes(ymin = mean - ci, ymax = mean + ci), 
                width = .1, 
                position=pd) +
  theme_bw() + 
  labs(x="Age",
       y="Mean Density",
       title="Mean Plot with 95% Confidence Interval") 

fat$ab.hips.ratio = fat$Abdomen/fat$Hip
plot(fat$Density,fat$ab.hips.ratio)
fit4 = lm(Density ~ ab.hips.ratio,data=fat)
summary(fit4)

ggplot(data = fat,
       mapping = aes(x = ab.hips.ratio, y = Density, 
                     color=young, shape=young, linetype=young)) +
  geom_point(alpha = .7, size = 3) +
  geom_smooth(method = "lm", se = FALSE, size = 1.5)



fit5 = lm(Density ~ Weight*Height + Neck + Abdomen + I(Hip^2) + ab.hips.ratio + Thigh + Biceps,data=fat)
relweights(fit5,color="dodgerblue")
summary(fit5)
par(mfrow=c(2,2))
par(mar=c(2,2,2,2))
plot(fit5)
gvlmaStat = gvlma(fit5)
summary(gvlmaStat)
ncvTest(fit5)
vif(fit5,type = "predict")
vif(fit5) > 10

leaps <-regsubsets(Density ~ Age + Weight*Height + Neck + Chest + Abdomen + Hip + Thigh + Knee + Biceps + Forearm + Wrist + ab.hips.ratio, data=fat)
plot(leaps, scale="adjr2")

fit = lm(Density ~ Age + Weight*Height + Neck + Chest + I(log(Abdomen)) + I(Hip^2) + Thigh + Knee + Biceps + Forearm + Wrist + ab.hips.ratio, data=fat)
summary(fit)

step(fit, direction="both")
