## Paul E. Johnson CRMDA <pauljohn@ku.edu>
## Portable Parallel Seeds Project.
## 2012-11-30


##' Selects among available random streams.
##'
##' R's global environment variable .Random.seed is re-set so that
##' random numbers generated after this call will follow from a
##' designated random number state. This will fail unless the
##' \code{setSeeds} function has already been called. This
##' function simply selects the n'th random stream for use,
##' presupposing that those stream states have already been set in the
##' global environment.
##'
##' This should handle three chores.  1) Copy the existing state of
##' the random generator (.Random.seed) into
##' currentStates[[currentStream]]. That allows us to "go back" to
##' that generator when we want to.  2) Change currentStream variable
##' to n, then 3) Re-set R's .Random.seed variable to the value from
##' currentStates[[currentStream]], so that successively drawn random
##' numbers follow the proper generator.
##' @export useStream
##' @param n An integer that selects which random stream should be
##'     used for the following work.
##' @param origin True or False. Should the stream be set at its
##'     original starting position, so as to re-generate the stream
##'     starting from the beginning? If FALSE, random numbers are
##'     drawn from the stream's current state.
##' @param verbose Requests detailed output for diagnostics.
##' @return Nothing. This function is run to change environment
##'     variables.
##' @author Paul E. Johnson <pauljohn@@ku.edu>
useStream <- function(n = NULL, origin = FALSE, verbose = FALSE){
    oldseed <-
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        else stop("in useStream, .Random.seed was NULL")
    ## get local copies of currentStream, currentStates
    curStream <- get("currentStream", envir = .pps, inherits = FALSE)
  curStates <- get("currentStates", envir = .pps, inherits = FALSE)
    
    if (n > length(curStates)) stop("requested stream does not exist")
    curStates[[curStream]] <- oldseed
    if (origin) {
        strtStates <- get("startStates", envir = .pps, inherits = FALSE)
        assign(".Random.seed", strtStates[[n]], envir = .GlobalEnv)
    } else {
        assign(".Random.seed", curStates[[n]], envir = .GlobalEnv)
    }
    ## put currentStream and currentStates back to .pps
    
    assign("currentStream", n, envir = .pps)
    assign("currentStates", curStates, envir = .pps)
    if (verbose){
        print("useStream useStream useStream useStream")
        print("CurrentStream CurrentStream CurrentStream")
        print( get("currentStream", envir = .pps, inherits = FALSE) )
        print("Current .Random.seed")
        print(.Random.seed)
    }
}


##' Return integer representing currently selected random stream
##'
##' The environment variable currentStates is a list of states
##' of possible random generators. The variable currentStream
##' indicates which of those generators is currently being used
##' by R.
##'
##' @export getCurrentStream
##' @return Integer index value of currently selected stream
##' @author Paul E. Johnson <pauljohn@@ku.edu>
##' @examples
##' mySeeds <- seedCreator(500, 5, file="mySeeds.rds", seed = 123123)
##' setSeeds(mySeeds, run = 17)
##' runif(2)
##' getCurrentStream()
##' useStream(2)
##' runif(2)
getCurrentStream <- function(){
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    else stop("in useStream, .Random.seed was NULL")
    ## get local copies of currentStream, currentStates
    curStream <- get("currentStream", envir = .pps, inherits = FALSE)
    curStream
}



##' Return one or all state vectors, either from original position or
##' current position.
##'
##' The environment variable currentStates is a list of states of
##' random generators. In most use cases, users will not need to use
##' this function because the details are handled by useStream.  Users
##' may need to explicitly access the stream states if they want
##' to extract that information and then pass it along to other
##' functions that do not use the built-in R random generator framework.
##'
##' @return Integer index value of currently selected stream
##' @author Paul E. Johnson <pauljohn@@ku.edu>
##' @param stream Optional. Selects a particular stream's state.  If
##'     no value is supplied, then the list including all of the
##'     stream states is supplied.
##' @param origin If FALSE (default), return the current, updated
##'     state of one or all random stream. If TRUE, return the initial
##'     stream states.
##' @export
##' @examples
##' mySeeds <- seedCreator(500, 5, file = "mySeeds.rds", seed = 123123)
##' setSeeds(mySeeds, run = 17)
##' runif(2)
##' getCurrentStream()
##' getState(origin=TRUE)
##' getState(stream = 2)
##' getState()
##' useStream(2)
##' runif(2)
getState <- function(stream, origin = FALSE){
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    else stop("in useStream, .Random.seed was NULL")
    ## get local copies of currentStream, currentStates
    if (origin == TRUE){
        if (exists("startStates", envir = .pps, inherits = FALSE))
            states <- get("startStates", envir = .pps, inherits = FALSE)
        else stop("startStates is not found in the environment.  Did you damage it?")
        if (missing(stream)) return(states)
        else if (length(states) >= stream) return (states[[stream]])
        else stop("Function getState: requested stream does not exist")
    } else {
        if (exists("currentStates", envir = .pps, inherits = FALSE))
            states <- get("currentStates", envir = .pps, inherits = FALSE)
        else stop("currentStates is not found in the environment.  Did you damage it?")

        if (missing(stream)) return(states)
        else if (length(states) >= stream) return (states[[stream]])
        else stop("Function getState: requested stream does not exist")
    }
    stop("Reached the end of getState without returning, some error occurred.")
}



##' Brings saved random streams back to life. Reads a portable
##' parallel seeds object (or file) and sets the seed collection in
##' the environment.
##'
##' The portable seeds object is created by the function
##' \link{seedCreator}. It is a list of lists. The list includes one
##' set of initializing states for each separate run of a simulation.
##' Within each of these sets, there will be enough information to
##' initialize one or more streams of random numbers.  These of
##' "initializing states" are the internal states of CMRG random
##' generators (see L'Ecuyer, 1999; L'Ecuyer, et al, 2002).
##'
##' This function scans the project's portable parallel seeds (either
##' an in-memory object or a named file), selects the desired run, and
##' then it writes 3 variables into the environment called ".pps". 1)
##' startStates is a collection of random generator states, one for
##' each random stream from which the user might wish to draw random
##' numbers. This is a fixed value which should not be altered during
##' the program. It can be used to reset the generators to their
##' initial positions. 2) currentStates is the collection of random
##' generator states that will be updated. When the program calls the
##' \link{useStream} function, the currentStates vector is updated.
##' 3) currentStream indicates which of the currentStates should be
##' used to draw the next random value.
##'
##' At the outset, startStates and currentStates are identical and
##' currentStream equals 1 (meaning the first element of currentStates
##' is taken as the state of the random generator).
##' @export setSeeds
##' @title setSeeds
##' @param projSeeds Required. Either an object of class portableSeeds
##' (created by \code{seedCreator}) or a text string giving the name
##' of an R saved file of the appropriate format (created by the
##' seedCreator function).
##' @param run Integer indicating which element from the portable seed
##' collection is to be selected
##' @param verbose Optional. Print out the state of the current
##' generator. Default = FALSE.
##' @return nothing is returned. This function is used for the side
##' effect of setting three objects in the global environment, the
##' startStates (list), currentStates (list), and currentStream (an
##' integer).
##' @author Paul E. Johnson <pauljohn@@ku.edu>
##' @seealso \code{seedCreator} to generate the input file for this
##' function and \code{useStream} to change from one stream to
##' another.
##' @example inst/examples/pps-ex-1.R
setSeeds <- function(projSeeds, run, verbose = FALSE){
    RNGkind("L'Ecuyer-CMRG")
    
    if (missing(projSeeds)) {
        stop(paste("setSeeds requires a seed object in order",
                   "to initialize the random streams"))
    } else if (is.character(projSeeds)){
        projSeeds <- readRDS(projSeeds)
    }
    if (class(projSeeds) != "portableSeeds"){
        stop(paste("Inappropriate project seed object supplied.",
                   "The projSeeds object must be created by the",
                   "seedCreator function, which would have set its",
                   "class as portableSeeds"))
    }
    
    
    if (missing(run)) stop(paste("run must be specified.",
                                 "Which replication is to be re-initialized?"))
    ## if (length(projSeeds) < run) stop(paste("The project seed object does not include enough elements to draw the one you are asking for. The seed object includes only ", length(projSeeds), " objects."))
    runSeeds <- projSeeds[[run]]
    assign("currentStream",  1L, envir = .pps)
    assign("startStates", runSeeds, envir = .pps)
    assign("currentStates", runSeeds, envir = .pps)
    assign(".Random.seed", runSeeds[[1L]],  envir = .GlobalEnv)
    if (verbose){
        print(paste("setSeeds, Run = ", run))
                               print(.Random.seed)
        print(paste("CurrentStream =", get("currentStream", envir = .pps, inherits = FALSE)))
                               print("All Current States")
        print(paste(get("currentStates", envir = .pps, inherits = FALSE)))
    }
}
NULL


##' Sets initial states for separate random generator
##' streams into the .pps environment.
##'
##' The main argument is the collection of seeds that is in the proper
##' format. It should be a "row" of a "portableSeeds" array object.
##' The currentStream argument determines which of these is set into
##' the R environment variable .Random.seed.
##'
##' This is different from \code{\link{setSeeds}} because setSeets accepts
##' the whole parallelSeeds object as input, whereas this function
##' wants only the input for one run of a simulation.
##'
##' In some specific use cases, particuarly the replication of individual
##' runs, it may be easier to use this function rather than setSeeds, but
##' both achieve the same purpose.
##'
##' @export setSeedVector
##' @param runSeeds Required. Seeds for the L'Ecuyer-CMRG random generator
##' @param currentStream Optional. Integer indicating which of the
##'     streams is to be used first. Default = 1.
##' @param verbose Optional. Default = FALSE.
##' @return nothing is returned. This function is used for the side
##'     effect of setting three objects in the environment, the
##'     startStates (list), currentStates (list), and currentStream
##'     (an integer).
##' @author Paul E. Johnson <pauljohn@@ku.edu>
##' @seealso \code{seedCreator}, \code{setSeeds}, and \code{useStream}.
##' @example inst/examples/pps-ex-2.R
setSeedVector <- function(runSeeds, currentStream = 1L, verbose = FALSE){
    RNGkind("L'Ecuyer-CMRG")

    if (missing(runSeeds)) {
        stop("setSeedVector requires a seed object in order to initialize the random streams")
    }

    ##TODO: find out what checks on the runSeeds elements are necessary.
    ## Should check each first element is 407? Make sure length of each is 7?
    ##
    assign("currentStream",  as.integer(currentStream), envir = .pps)
    assign("startStates", runSeeds, envir = .pps)
    assign("currentStates", runSeeds, envir = .pps)
    assign(".Random.seed", runSeeds[[1L]],  envir = .GlobalEnv)
    if (verbose){
        print(paste("setSeedVector"))
        print(.Random.seed)
        print(paste("CurrentStream =", get("currentStream", envir = .pps, inherits = FALSE)))
        print("All Current States")
        print(paste(get("currentStates", envir = .pps, inherits = FALSE)))
  }
}
NULL


##' The name previously used for the function
##' \code{\link{setSeeds}}. 
##'
##' This is replaced by \code{\link{setSeeds}}, which does exactly the
##' same thing, but with a name that is more immediatly understandable
##' to users.  This is now deprecated, will be removed.
##' @export
##' @param projSeeds Required. Either an object of class portableSeeds
##'     (created by \code{\link{seedCreator}}) or a text string giving the
##'     name of an R saved file of the appropriate format (created by
##'     the seedCreator function).
##' @param run Integer indicating which element from the portable seed
##'     collection is to be selected
##' @param verbose Optional. Print out the state of the current
##'     generator. Default = FALSE.
##' @return nothing is returned. This function is used for the side
##'     effect of setting three objects in the global environment, the
##'     startStates (list), currentStates (list), and currentStream
##'     (an integer).
##' @author Paul E. Johnson <pauljohn@@ku.edu>
##' @seealso \code{\link{setSeeds}} which is the replacement for this
##'     function.
initPortableStreams <- function(projSeeds, run, verbose = FALSE){
    setSeeds(projSeeds, run, verbose)
}
NULL
