## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

runlist <- function(...)
{
    object <- list(...)
    if (is.null(names(object)))
        names(object) <- rep("", length(object))
    unnamed <- names(object) == ""
    if (any(unnamed)) {
        dnames <- as.list(substitute(list(...)))[-1]
        for (i in seq_along(object)) {
            if (names(object)[i] == "") {
                names(object)[i] <- toString(deparse(dnames[[i]]), width = 12)
            }
        }
    }
    if (length(object) > 0)
        class(object) <- paste(class(object[[1]]), "runlist", sep = ".")
    class(object) <- unique(c(class(object), "runlist"))
    object
}

as.runlist <- function(x, ...)
    do.call("runlist", as.list(x))

"[.runlist" <- function (x, i, ...)
    structure(NextMethod("["), class = class(x))

print.runlist <-
    function(x, ...)
{
    cat("\nList of model runs:\n")
    print.default(x, ...)
    invisible(x)
}

print.hydromad.runlist <-
    function(x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("\nList of Hydromad model runs:\n")
    print.default(lapply(x, function(obj) {
        if (!is.null(obj$msg)) {
            list(call = obj$call, message = obj$msg)
        } else {
            obj$call
        }
    }), digits=digits, ...)
    invisible(x)
}

print.tf.runlist <-
    function(x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("\nList of transfer function model runs:\n")
    print.default(lapply(x, function(obj) {
        if (inherits(obj, "tf"))
            coef(obj) else obj
    }), digits=digits, ...)
    invisible(x)
}

errormasscurve.runlist <-
    function(x,
             coerce=byDays,
             superpose=TRUE,
             type=c("total", "average", "relative"),
             ...)
{
    type <- match.arg(type)
    if (type == "relative") {
        Q <- observed(x[[1]])
        Q[is.na(Q)] <- 0
        Qcum <- cumsum(Q)
        cumsumok <- function(Q) { Q[is.na(Q)] <- 0; cumsum(Q) / Qcum }
    } else if (type == "average") {
        cumsumok <- function(Q) { Q[is.na(Q)] <- 0; cumsum(Q) / seq_along(Q) }
    } else {
        cumsumok <- function(Q) { Q[is.na(Q)] <- 0; cumsum(Q) }
    }
    foo <- xyplot(x, residuals=TRUE,
                  coerce=coerce, trans=cumsumok,
                  superpose=superpose, ...)
    foo$call <- sys.call(sys.parent())
    foo
}

xyplot.runlist <-
    function(x, data,
             residuals=FALSE,
             coerce=byDays, trans=NULL,
             superpose=FALSE,
             panel=panel.superpose.cycle,
             ...)
{
    if (!missing(data) && !is.null(data))
        warning("'data' argument ignored.")

    Q <- if (!residuals) coerce(observed(x[[1]]))
    if (residuals) x <- lapply(x, residuals)
    else x <- lapply(x, fitted)
    x <- do.call(cbind, x)
    x <- coerce(x)

    transfn <- eval(trans)
    if (is.character(transfn)) transfn <- get(trans)
    if (!is.null(trans)) {
        x <- transfn(x)
        if (!residuals) Q <- transfn(Q)
    }

    if (residuals) {
        foo <- xyplot(x,
                      superpose = superpose,
                      panel=panel,
                      ...)

    } else {
        ## observed vs modelled
        foo <- xyplot(x,
                      superpose = superpose,
                      panel=panel,
                      ...)
        foo <- foo + layer(panel.lines(Q), data=list(Q=Q), style=if (superpose) (NCOL(x) + 1) else 2)
    }

                                        #if (FALSE &&             residuals == FALSE) {
                                        #	scales$y$limits <- extendrange(Q[is.finite(Q)])
                                        #	if (identical(trans, "log")) {
                                        #		# limit scales
                                        #		minQ <- min(Q[is.finite(Q)]) - 1 # log scale
                                        #		data[ (data[] < minQ) ] <- minQ
                                        #		Q[is.infinite(Q)] <- minQ
                                        #		scales$y$limits[1] <- minQ
                                        #		#c(minQ, max(c(Q[], data[]), na.rm=T) + 0.5)
                                        #	}
                                        #}

    foo$call <- sys.call(sys.parent())
    foo
}

qqmath.runlist <-
    function(x, data,
             residuals=FALSE,
             coerce=byDays, trans=NULL,
             superpose=TRUE,
             auto.key=TRUE,
             type = "b",
             pch=".", cex=2,
             ...)
{
    if (!missing(data) && !is.null(data))
        warning("'data' argument ignored.")
    Q <- if (!residuals) coerce(observed(x[[1]]))
    if (residuals) x <- lapply(x, stats::residuals)
    else x <- lapply(x, fitted)
                                        #x <- lapply(x, coerce)
    x <- do.call(cbind, x)
    x <- coerce(x)

    transfn <- eval(trans)
    if (is.character(transfn)) transfn <- get(trans)
    if (!is.null(trans)) {
        x <- transfn(x)
        if (!residuals) Q <- transfn(Q)
    }

    dat <- do.call(make.groups, as.data.frame(x))
    if (!residuals) dat <- rbind(make.groups(observed=Q), dat)

    if (!identical(auto.key, FALSE)) {
        if (isTRUE(auto.key)) auto.key <- list()
        if (is.null(auto.key$lines)) auto.key$lines <- TRUE
        if (is.null(auto.key$points)) auto.key$points <- FALSE
        #auto.key$type <- type[[1]]
        #auto.key$pch <- pch
        ## can not set cex here! need to convert to key
    }
                                        #	tsobj <- coerce(obsmod(x))
                                        #	transfn <- eval(trans)
                                        #	if (is.character(transfn)) transfn <- get(trans)
                                        #	if (!is.null(trans)) tsobj <- transfn(tsobj)
                                        #	dat <- make.groups(observed=Q, modelled=tsobj[,"mod"])
                                        # remove cases with any missing values
                                        #tsdat <- tsdat[complete.cases(tsdat),]
                                        #series <- factor(rep(colnames(tsdat), each=NROW(tsdat)))
                                        #tsdat <- as.vector(coredata(tsdat))
    foo <- qqmath(~ data, groups=which, data=dat, auto.key=auto.key,
                  type=type, pch=pch, cex=cex, ...)
    foo$call <- sys.call(sys.parent())
    foo
}

summary.hydromad.runlist <-
summary.tf.runlist <-
    function(object,
             pars=TRUE,
             flowstats=c("rel.bias","r.squared","r.sq.sqrt","r.sq.log","r.sq.monthly"),
             ...)
{
    stopifnot(is.list(object))
    object <- object[sapply(object, isValidModel)]
    if (length(object) == 0)
        return(NULL)
    ## if objects in runlist are different classes, set pars = FALSE
    classes <- sapply(object, function(x) class(x)[1])
    if (pars && any(classes != classes[1])) {
        pars <- FALSE
    }
    ## if any different parameter names, set pars = FALSE
    if (pars) {
        allpars <- lapply(object, coef)
        if (any(diff(sapply(allpars, length) != 0))) {
            ## number of parameters differs
            pars <- FALSE
        } else {
            ## paste parameter names together
            allparnames <- sapply(allpars, function(x) toString(names(x)))
            if (any(allparnames != allparnames[1])) {
                pars <- FALSE
            }
        }
    }

    summstats <- lapply(object, function(obj) {
        flow.summ <- summary(obj, which=flowstats, ...)
        c(if (pars) coef(obj),
          unlist(flow.summ[flowstats]))
    })

    ## pad out missing values with NAs
    ok.i <- which.max(sapply(summstats, length))
    varnames <- names(summstats[[ok.i]])
    summstats <- lapply(summstats, function(x) x[1:length(varnames)])
    summstats <- data.frame(do.call(rbind, summstats))
    names(summstats) <- varnames
    class(summstats) <- c("summary.runlist", class(summstats))
    summstats
}

print.summary.runlist <-
    function(x, digits = max(3, getOption("digits") - 3), ...)
{
    ## just simplify the printed output by rounding
    print.data.frame(x, digits=digits, ...)
    invisible(x)
}

