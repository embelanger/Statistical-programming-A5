## Em Belanger: s2409525

## Estimating excess deaths in 2020

###############################################################################
## Function for simulating the predicted number deaths for a certain number of
## weeks based on the model:  
## Di = 0.9885*dj*qi*Ni (the 0.9885 is to account for using yearly, rather than 
## weekly, age data in an ageing population)
## Where:
## Di = the deaths for age group i in week j
## qi = 1-e^(-mi/52) (the expected proportion of Ni dying in a week)
## dj: alters qi to account for seasonal fluctuations in death rates
##
## The simulation ages the population for each iteration according to the following 
## model:
## Ni_star = Ni-Di
## Ni_plus = Ni_star*(51/52)+N(i-1)_star/52
##
## The simulation assumes a constant birth rate
##
## Inputs:
## female: a vector of the starting populations of females for age 0-100
## male: a vector of the starting populations of males for age groups 0-100
## mf: a vector of the corresponding mi values for the females
## mm: a vector of the corresponding mi values for the males
## d: a vector of dj values (not sex specific)
##
## Outputs:
## death_total: a vector of the predicted total number of deaths each week
###############################################################################

require(rjags)
require(coda)

death_sim <- function(female, male, mf, mm, d){
  ## Calculate qf and qm based on the supplied mf and mm
  qf <- 1-exp(-mf/52)
  qm <- 1-exp(-mm/52)
  
  ## initialize vectors to store death values, Df and Dm are to store that week's
  ## deaths for males and females. death_total is the sum of all deaths for that 
  ## week
  Df <- rep(0, length(female))
  Dm <- rep(0, length(male))
  death_total <- rep(0, length(d))
  
  ## Initialize vectors for updating populations
  
  Nf_star_prev <- rep(0, (length(female)+2))
  Nf_star <- rep(0, length(female))
  Nm_star_prev <- rep(0, (length(male)+2))
  Nm_star <- rep(0, length(male))
  
  ## function assumes constant birth rate, so N_0 should always be the original
  ## age 0 population
  N_0f <- female[1]
  N_0m <- male[1]
  
  for (j in 1:length(d)){
    ## female deaths
    Df <- 0.9885*d[j]*qf*female
    
    ## male deaths
    Dm <- 0.9885*d[j]*qm*male
    
    ## sum the deaths, store in the deaths vector
    death_total[j] <- sum(Df) + sum(Dm)
    
    ## population adjustment
    Nf_star_prev <- c(N_0f, (female - Df))
    Nf_star_prev <- head(Nf_star_prev, -1)
    Nf_star <- female - Df
    
    female <- Nf_star*(51/52) + Nf_star_prev/52
    
    Nm_star_prev <- c(N_0m, (male-Dm))
    Nm_star_prev <- head(Nm_star_prev, -1)
    Nm_star <- male - Dm
    
    male <- Nm_star*(51/52) + Nm_star_prev/52
    
  }
  return(death_total)
}

## Difference in total actual deaths and total predicted deaths from the start 
## of 2020 to the end of the data (excess deaths)

table1 <- read.table("lt1720uk.dat")
table2 <- read.table("death1722uk.dat")
female <- table1$fpop20
male <- table1$mpop20
mf <- table1$mf
mm <- table1$mm
d <- table2$d
d1 <- d[157:length(d)]
deaths <- table2$deaths

sum(deaths[157:length(deaths)] - death_sim(female, male, mf, mm, d1))

## Excess deaths for 2020 only

d2 <- d[157:208]
sum(deaths[157:208] - death_sim(female, male, mf, mm, d2))

## Plot observed deaths against week (make sure lower y-axis lim is zero!!!)
## overlay as a continuous curve the predicted deaths
## Report excess deaths in 2020 and overall in the plot title
plot(table2$week[157:length(deaths)], deaths[157:length(deaths)], 
     ylim = c(0, 25000), xlab = "Week", ylab = "Deaths", pch = 20, 
     main = "Excess deaths from 2020 onwards:  181, 248.4\nExcess deaths in 2020: 89, 602.38")
lines(table2$week[157:length(deaths)], death_sim(female, male, mf, mm, d1), col = "#AF1313", lwd = 3)
legend(x = "bottomright", pch = c(20, NA), lwd = c(NA, 3), 
       col = c("black", "#AF1313"), legend = c("Actual Deaths", "Predicted Deaths"))

## Compute vector of excess deaths each week and produce a plot of the 
## cumulative excess deaths by week (cumsum is useful) from 2020 onwards
excess_vec <- deaths[157:length(deaths)] - death_sim(female, male, mf, mm, d1)
excess_cum <- cumsum(excess_vec)
plot(table2$week[157:length(deaths)], excess_cum, xlab = "Week", ylab = "Deaths",
     main = "Cumulative Excess Deaths from 2020 Onwards", type = "l", lwd = 3)

## Use JAGS to estimate the parameters for the time series model for excess deaths in 2020:
## mu[1] <- alpha
## x[1] ~ dt(mu[1], tau, k)
## mu[i] <- (x[i-1] - alpha)*rho + alpha
## x[i] ~ dt(mu[i], tau, k)
##
## Use the following priors:
## tau ~ exp(1)
## rho ~ uniform(0, 0.9)
## alpha ~ normal(0, 0.0001)
## k ~ uniform(2, 100)
##
## Set weeks 51, 52, 53, 105, 106 to NA because of data collection issues (xmas/new years)
excess_vec[c(51, 52, 53, 105, 106)] <- NA
time_mod <- jags.model("model.jags", data=list(x=excess_vec, N=length(excess_vec)))
samps <- coda.samples(time_mod, c("mu", "rho", "k"), n.iter=10000)

## Produce trace plots and histograms for rho
plot(samps[[1]][,151], density = FALSE, main = "Trace Plot of Rho")
plot(samps[[1]][,151], trace = FALSE, main = "Density Plot of Rho")

## compute the posterior expected value vector for mu
E_mu_post <- rep(0, 150)

for (i in 2:150){
  E_mu_post[i] <- mean(samps[[1]][,i])
}
## (the first element didn't get filled with a mu because of the order that the parameters were estimated)
E_mu_post <- E_mu_post[-1]

## produce a plot showing every 50th sampled mu vector (against week) as a grey 
## curve, with the estimated expectation of mu overlaid in blue. Add the 
## observed excess deaths, xi, in black symbols, with the ones set to NA in red
plot(table2$week[157:length(deaths)], samps[[1]][50,2:150], type = "l", col = "#878787", 
     ylim = c(-7000, 15000), ylab = "Excess Deaths", xlab = "Week", 
     main = "Plotting Excess Death Estimates From 2020 Onwards")
 
for (i in 1:49){
  lines(table2$week[157:length(deaths)], samps[[1]][50*i,2:150], col = "#878787")
}

lines(table2$week[157:length(deaths)], E_mu_post, col = "#1C76D0")
points(table2$week[157:length(deaths)], excess_vec, pch = 4)

NA_vec <- rep(NA, 149)
excess_tot <- deaths[157:length(deaths)] - death_sim(female, male, mf, mm, d1)
NA_vec[c(51, 52, 53, 105, 106)] <- excess_tot[c(51, 52, 53, 105, 106)]

points(table2$week[157:length(deaths)], NA_vec, pch = 4, col = "#AF1313")

legend(x = "topright", legend = c("mu samples", "mean of mu samples", "excess estimates", "removed estimates"),
       lwd = c(1, 1, NA, NA), pch = c(NA, NA, 4, 4), col = c("#878787", "#1C76D0", "black", "#AF1313"))

## Plot the "residuals" from this model, xi - E(mui), against time.
time <- table2$week[157:length(deaths)]-157
resid <- excess_vec - E_mu_post
plot(time, resid, type = "p", main = "Residuals Plot", ylab = "xi - E(mui)", xlab = "Time" )







