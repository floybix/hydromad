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
    class(object) <- unique(c(class(object), "runlist", "list"))
    object
}

as.runlist <- function(x, ...)
    do.call("runlist", as.list(x))

"[.runlist" <- function (x, i, ...)
    structure(NextMethod("["), class = class(x))

coef.runlist <-
    function(object, ...)
{
    summary(object, ..., FUN = coef)
}


summary.runlist <-
    function(object, ..., FUN = summary, items = NULL)
{
    stopifnot(is.list(object))
    if (length(object) == 0)
        return(NULL)
    ## extract elements from summary which are single numbers
    cc <- lapply(object, function(x) {
        tmp <- FUN(x, ...)
        if (is.null(items)) {
            tmp <- tmp[unlist(lapply(tmp, function(z) {
                is.numeric(z) && !is.matrix(z) &&
                (length(z) == 1)
            }))]
        } else {
            tmp <- tmp[items]
        }
        unlist(tmp)
    })
    ## pad out missing entries with NAs
    ## find the set of all names
    allnms <- unique(unlist(lapply(cc, names)))
    ans <- matrix(NA_real_,
                  nrow = length(object),
                  ncol = length(allnms),
                  dimnames = list(names(object), allnms))
    for (i in 1:NROW(ans))
        ans[i, names(cc[[i]])] <- cc[[i]]
    ans <- as.data.frame(ans)
    class(ans) <- c("summary.runlist", class(ans))
    ans
}

print.summary.runlist <-
    function(x, digits = max(3, getOption("digits") - 3), ...)
{
    ## just simplify the printed output by rounding
    print.data.frame(x, digits=digits, ...)
    invisible(x)
}

print.runlist <-
    function(x, ...)
{
    cat("\nList of model runs:\n")
    print.default(x, ...)
    invisible(x)
}

residuals.runlist <-
    function(object, ...)
{
    ans <- lapply(object, residuals, ...)
    bad <- sapply(ans, length) == 0
    if (any(bad))
        stop("residuals() returned nothing for items ",
             toString(names(ans)[bad]))
    do.call("cbind", ans)
}

fitted.runlist <-
    function(object, ...)
{
    ans <- lapply(object, fitted, ...)
    bad <- sapply(ans, length) == 0
    if (any(bad))
        stop("fitted() returned nothing for items ",
             toString(names(ans)[bad]))
    do.call("cbind", ans)
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

## TODO: xyplot.hydromad.runlist
## 

## TODO! clean this up...
xyplot.runlist <-
    function(x, data = NULL,
             residuals. = FALSE,
             coerce = byDays,
             superpose = TRUE,
             with.P = FALSE,
             ...)
{
    stopifnot(is.null(data))

    xdat <- lapply(x, if (residuals.) residuals else fitted)
    bad <- sapply(xdat, length) == 0
    if (any(bad))
        stop("residuals() or fitted() returned nothing for items ",
             toString(names(xdat)[bad]))
    xdat <- lapply(xdat, coerce)
    
    if (superpose) {
        xm <- do.call(cbind, xdat)
        if (!residuals.) {
            ## include observed series from item 1 (assuming all are the same!)
            xm <- cbind(obs = coerce(observed(x[[1]])), xm)
        }
        foo <- xyplot(xm, superpose = TRUE, ...)
    } else {
        ## juxtapose, not superpose: just call xyplot on each item
        foo <- xyplot.list(xdat, ...,
                           x.same = x.same, y.same = y.same, layout = layout)
        if (!residuals.) {
            foo <- foo +
                layer(panel.lines(obs), data = list(obs = coerce(observed(x[[1]]))),
                      style = 2)
        }
    }
    if (with.P) {
        foo <- c(foo, rainfall = xyplot(coerce(x[[1]]$data[,"P"]), ...),
                 x.same = x.same, y.same = NA, layout = layout)
    }
    #if (!residuals) {
    #    ## observed vs modelled: add 
    #    foo <- foo +
    #        layer(panel.lines(Q), data=list(Q=Q), style=if (superpose) (NCOL(x) + 1) else 2)
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
