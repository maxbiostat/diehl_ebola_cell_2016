library(lme4)
library(ggplot2)
###########################################
source("diehl_aux.r")
##########################################
## Data preparation
raw <- read.csv("../data/GP82_raw_data.csv", header = TRUE) ## "all" data
EVD_GP82 <- read.csv("../data/EVD_case_fatality_GP82.csv", header = TRUE) ## clean data
coord_data <- read.csv("../data/coord_data_locations.csv", header = TRUE)
pop_level <- read.csv("../data/population_level.csv", header = TRUE)
EVD_GP82 <- merge(EVD_GP82, coord_data, by = "location")
EVD_GP82 <- merge(EVD_GP82, pop_level, by = "location")
EVD_GP82 <- subset(EVD_GP82, country.y == "GIN") ## restrict attention to Guinea
xtabs(~ Outcome + GP82, EVD_GP82)
nrow(EVD_GP82)
######
### Exploratory statistical analysis (ESA)

## Question 1: Do Ct values according to EBOV variant ?
# All data available
all_A <- raw$Assay.CT.Value[raw$GP82 == "A"]
sum(!is.na(all_A)) ## valid observations
all_V <- raw$Assay.CT.Value[raw$GP82 == "V"]
sum(!is.na(all_V))
t.test(all_A, all_V)

# Just the data we'll use for further analyses
filtered_A <- EVD_GP82$Assay.CT.Value[EVD_GP82$GP82 == "A"]
sum(!is.na(filtered_A))
filtered_V <-  EVD_GP82$Assay.CT.Value[EVD_GP82$GP82 == "V"]
sum(!is.na(filtered_V))
t.test(filtered_A, filtered_V)

#######
pdf("../plots/Ct_per_genotype.pdf")
ggplot(EVD_GP82) +
  geom_boxplot(aes(x = GP82, y = Assay.CT.Value, fill = GP82), colour = c("grey80", "black")) +
  scale_fill_manual(values = c("black", "red")) +
  theme_bw()
dev.off()
pdf("../plots/Ct_per_genotype_and_outcome.pdf")
ggplot(EVD_GP82) +
  geom_boxplot(aes(x = Outcome, y = Assay.CT.Value, fill = Outcome))+
  facet_grid(~GP82)
dev.off()

#### Now let's look at odds ratios
(CrossTab_gp82 <- xtabs(~ GP82 + Outcome , EVD_GP82) )
summary(CrossTab_gp82)

a <- CrossTab_gp82["V", "died"]
b <- CrossTab_gp82["V", "survived"]
c <- CrossTab_gp82["A", "died"]
d <- CrossTab_gp82["A", "survived"]  
SqrT <-  sqrt(1/a + 1/b + 1/c + 1/d)
(rawOR <- (a*d)/(b*c) )
alpha <- .95
( lwrOR <- exp(log(rawOR) + qnorm((1 - alpha)/2) * SqrT) ) 
( uprOR <- exp(log(rawOR) + qnorm((1 + alpha)/2) * SqrT) ) 

## Relative risk
m.RR <- log({a/(a+b)} / {c/(c + d)})
sdRR <- sqrt( (1/a + 1/c) -(1/(a + b) + 1/(c+d)))
Z <- qnorm(.975)
CI.RR <- c( exp(m.RR - Z*sdRR), exp(m.RR + Z*sdRR))
c(exp(m.RR), CI.RR)

### Ok, but what happens if we try to correct for other variables?
### First, some data preparation
forFit <- EVD_GP82
forFit$Assay.CT.Value <- -standz(EVD_GP82$Assay.CT.Value) ## invert sign because of lower CT = higher viral load
forFit$pdensMN <- standz(EVD_GP82$pdensMN)
forFit$geconMN <- standz(EVD_GP82$geconMN)
forFit$tt50kMN <- standz(EVD_GP82$tt50kMN)
forFit$cumCases <- standz(log(EVD_GP82$cumCases))
forFit$GP82 <- ifelse(EVD_GP82$GP82 == "V", 1, 0) ## GP82-A as baseline
forFit$Outcome <- ifelse(EVD_GP82$Outcome == "died", 1 , 0)  

simple_glm <- glm(Outcome ~ Assay.CT.Value + GP82, data = forFit, family = "binomial")
summary(simple_glm)
getOdds(simple_glm, hier = FALSE)

## (corrected)  relative risks
simple.pois.glm <- glm(Outcome ~ Assay.CT.Value + GP82, data = forFit, family = "poisson")
exp(cbind(coefficients(simple.pois.glm), confint(simple.pois.glm))) 

NewCt <- seq(min(forFit$Assay.CT.Value), max(forFit$Assay.CT.Value), length.out = 100)
mut <- c(0, 1)
newDat <- expand.grid(NewCt, mut)
names(newDat) <- c("Assay.CT.Value", "GP82")
Pred.logit <- predict(simple_glm, newdata = newDat, se.fit = TRUE, type = "response")
invlogit <- function(x) exp(x)/(1 + exp(x))
alpha <- 0.95
InvLogitPred <- data.frame(lwr = invlogit(Pred.logit$fit + qnorm((1-alpha)/2) * Pred.logit$se.fit) ,
                           mean = invlogit(Pred.logit$fit),
                           upr = invlogit(Pred.logit$fit + qnorm((1 + alpha)/2) * Pred.logit$se.fit))
toPlot_preds <- data.frame(newDat, InvLogitPred)
toPlot_preds$GP82 <- ifelse(toPlot_preds$GP82 == 1, "V", "A")

pdf("../plots/predicted_fatality_rates.pdf")
ggplot(toPlot_preds, aes(x = Assay.CT.Value, y = mean)) +
  # ggtitle("Predicted fatality rates") +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = GP82), alpha = .2) +
  scale_x_continuous("Transformed Ct", expand = c(0, 0)) +
  scale_y_continuous("Fatality rate", expand = c(0, 0)) +
  scale_fill_manual(values = c("black", "red")) +
  geom_line(aes(colour = GP82), size = 1) + 
  scale_color_manual(values = c("black", "red")) +
  guides(fill = FALSE) +
  theme_bw()
dev.off()
model_0f <- glmer(Outcome ~ Assay.CT.Value + GP82 + (cumCases|location),
                  data = forFit, family = "binomial")
summary(model_0f)
getOdds(model_0f)
RR <- ranef(model_0f, condVar = TRUE)
names(RR) <- " "
names(RR[[1]]) <- c("Location-specific intercept", "coefficient of cumulative number of cases")

pdf("../plots/varying_intercepts.pdf")
lattice::dotplot(RR)
dev.off()