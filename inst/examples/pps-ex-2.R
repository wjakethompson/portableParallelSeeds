library(portableParallelSeeds)

projSeeds <- seedCreator(50, 3, seed = 123456)
A1 <- projSeeds[[27]]
A1 ## shows states of 3 generators for one run

##' Receive a run-specific seed vector and conduct a simulation
##'
##' Rather than giving this function the whole seed warehouse that is used for
##' all the simulation runs, this one just receives one vector of seeds. Compare
##' to \code{\link{setSeeds}}.
##'
##' This function does not need to know what the run number is, since it accepts
##' a stream vector as its first argument.  Prehaps this is a simpler, more
##' elegant approach than in \code\link{{setSeeds}}.
##' @param streams A stream element selected from a "portableSeeds" object.
##'   Created by \code\link{{seedCreator}}
##' @param N Number of rows in simulated data set
##' @param m mean of simulated data
##' @param sd standard deviation of simulated data
##' @param beta regression coefficients
##' @return A list of estimated regressions fitted to the simulated values,
##'   including 3 linear models and one generalized linear model
##' @author Paul Johnson <pauljohn@@ku.edu>
runOneSimulation <- function(streams, N, m, sd, beta = c(0.3, 0.2, 0.1))
{
    setSeedVector(streams)
    
    datX <- rockchalk::mvrnorm(N,
                               mu = rep(m, 3),
                               Sigma = sd * diag(3))
    colnames(datX) <- c("X1", "X2", "X3")
    useStream(2)
    eta <- datX %*% beta
    SigmaE <- sqrt(mean(eta))
    Y1 <- eta + SigmaE * rnorm(N)
    datX <- data.frame(datX, Y1)
    lm1 <- lm(Y1 ~ X1 + X2 + X3, data = datX)
    ## lets see a log link, just for fun
    datX$Y2 <- rpois(N, lambda = exp(eta))
    lm2 <- lm(Y2 ~ X1 + X2 + X3, data = datX)
    lm3 <- lm(log(Y2+0.5) ~ X1 + X2 + X3, data = datX)
    glm3 <- glm(Y2 ~ X1 + X2 + X3,  data = datX, family = poisson(link = "log"))
    list("lm1" = lm1, "lm2" = lm2, "lm3" = lm3, "glm3" = glm3)
}


serial1 <- lapply(projSeeds[1:30], runOneSimulation, N = 1000, m = 2,
                  sd = 10, beta = c(0.5, 0.6, 0.1))

## Take the 27th simulated model
lapply(serial1[[27]], summary)

## Recall A1, the 27th seed vector, was created above.
## Repeat that run, explicitly:
serial1.27 <- runOneSimulation(A1, N = 1000, m = 2, sd = 10, beta =
c(0.5, 0.6, 0.1))

## Note they are same

sapply(serial1[[27]], coef)
sapply(serial1.27, coef)

## Check Again
identical(coef(serial1.27[[1]]), coef(serial1[[27]][[1]]))
identical(coef(serial1.27[[2]]), coef(serial1[[27]][[2]]))


