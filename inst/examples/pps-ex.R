library(portableParallelSeeds)

projSeeds <- seedCreator(2000, 3, seed = 123456, file = "fruits.rds")
A1 <- projSeeds[[787]]
A1 ## shows states of 3 generators for run 787


setSeeds(projSeeds, run = 787, verbose = TRUE)
.Random.seed
getCurrentStream()
runif(4)

## read from file, take run 787's seed
myFruitySeeds <- readRDS("fruits.rds")
B1 <- myFruitySeeds[[787]]

identical(A1, B1) # check


setSeeds("fruits.rds", run=787)
.Random.seed
runif(4)


runOneSimulation <- function(run, streamsource, N, m, sd){
    setSeeds(streamsource, run = run, verbose= FALSE)
    datX <- rnorm(N, mean = m, sd = sd)
    datXmean <- mean(datX)
    useStream(2)
    datY <- rpois(N, lambda = m)
    datYmean <- mean(datY)
    useStream(1)
    datXplusOne <- rnorm(1, mean = m, sd = sd)
    ## Should be N+1'th element from first stream
    c("datXmean" = datXmean, "datYmean" = datYmean, "datXplusOne" = datXplusOne)
}

\donttest{
## Give seed collection object to each simulation, let each pick desired seed
    serial1 <- lapply(1:1000, runOneSimulation, projSeeds, N=800, m = 14, sd = 10.1)
    
    
    ## First re-load the seed object, then give to simulations
    fruits2 <- readRDS("fruits.rds")
    serial2 <- lapply(1:1000, runOneSimulation, fruits2, N=800, m = 14, sd = 10.1)
    
    
    ## Re-load file separately in each run
    serial3 <- lapply(1:1000, runOneSimulation, "fruits.rds", N = 800, m = 14, sd = 10.1)
    identical(serial1, serial2)
    identical(serial1, serial3)
    ## The results form the 3 lapply statements are all identical.

    ## Lets check run 912, you'll see what I mean.
    serial1[[912]]
    serial2[[912]]
    
    
    ## Next question. Can we re-set the current interactive session's seed
    ## and re-draw that same value?

    ## Put the seeds from run 912 back into the current session
    setSeeds("fruits.rds", run = 912, verbose = FALSE)

    ## Draw 801 elements    
    X <- rnorm(801, m=14, sd = 10.1)
    ## The 801'th element should be same as saved in the third element
    ## of the serial1[[912]] object.
    X[801]
    
    if(!all.equal(unname(X[801] - serial1[[912]][3]), 0)){
        stop("The difference is not 0, something went wrong")
    } else {
        print("Celebrate. The value in X[801] equals the serial1[[912]][3] value")
    }
    
    ##Bingo. I'm right. Can draw a understandably replicatable streams of
    ## random numbers, whether we draw 800, switch to a different stream,
    ## and then change back to draw another, or if we just draw 801 in one
    ## series.
    
    unlink("fruits.rds") #delete file
}
