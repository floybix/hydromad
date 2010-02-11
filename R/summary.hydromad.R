## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


## TODO: rework this
summary.hydromad <-
    function(object, breaks = NULL, coerce = byDays,
             which = c("rel.bias", "r.squared", "r.sq.sqrt", "r.sq.log", "r.sq.monthly"),
             na.action = na.exclude,
             ...)
{
    if (!isFullySpecified(object))
        stop("model parameters not fully specified")
    if (!isValidModel(object))
        stop("model is not valid")

    ans <- list()
    copyVars <- c("call")
    ans[copyVars] <- object[copyVars]

    ## TODO: should this be a stand-alone function?
    hydrostats <- function(chunk) {
        ans <- numeric()
        ans[["timesteps"]] <- NROW(chunk)
        ok <- complete.cases(chunk)
        ans[["missing"]] <- sum(!ok)
        ans[["mean.P"]] <- meanP <- mean(chunk[,"P"], na.rm=TRUE)
        ans[["mean.Q"]] <- meanQ <- mean(chunk[,"Q"], na.rm=TRUE)
        ans[["runoff"]] <- meanQ / meanP
        ans
    }

    DATA <- ts.intersect(P = observed(object, item = "P", all = TRUE),
                         Q = observed(object, all = TRUE),
                         U = fitted(object, U = TRUE, all = TRUE),
                         X = fitted(object, all = TRUE))

    if (!is.null(breaks)) {
        nms <- colnames(DATA)
        DATA <- coerce(na.trim(DATA))
        colnames(DATA) <- nms
        group <- cut(time(DATA), breaks = breaks)
        group <- factor(group)
        ## remove warmup
        DATA <- stripWarmup(DATA, object$warmup)
        if (object$warmup > 0)
            group <- group[-(1:object$warmup), drop = TRUE]
        ## compute stats
        basicstats <- by(DATA, group, hydrostats)
        basicstats <- do.call(rbind, basicstats)
        basicstats <- rbind(basicstats,
                            overall=hydrostats(DATA))

        chunkstats <- by(DATA, group, perfStats,
                         warmup=0, na.action=na.exclude, ...)
        ## TODO: ensure equal lengths of rows / intact names?
        chunkstats <- do.call(rbind, chunkstats)
        chunkstats <- rbind(chunkstats,
                            overall=perfStats(DATA, warmup=0, na.action=na.action, ...))
        ans <- cbind(basicstats, chunkstats)
        ans <- data.frame(ans)
        class(ans) <- c("summaryWithBreaks.hydromad", class(ans))
        return(ans)
    }

    ## TODO:
    #
    ans$used.rfit <- object$used.rfit

    ## ARPE
    if ("yic" %in% which) {
        which <- union(which, "arpe")
    }
    if ("arpe" %in% which) {
        if (is.null(vcov(object))) {
            ans$arpe <- NA
        } else {
            ans$coef.var <- diag(vcov(object))
            cc <- coef(object, "routing")
            cc2 <- tfParsConvert(cc, "a,b")
            cc[names(cc2)] <- cc2
            nms <- intersect(names(cc), names(ans$coef.var))
            ans$arpe <- mean(ans$coef.var[nms] / (cc[nms]^2))
            if ("yic" %in% which) {
                var.ratio <- (var(residuals(object), na.rm = TRUE) /
                              var(observed(object), na.rm = TRUE))
                ans$yic <- log(var.ratio) + log(ans$arpe)
            }
        }
    }

    ## call perfStats for the rest
    ans <- c(ans, perfStats(DATA, warmup = object$warmup, which = which, ...))

    class(ans) <- "summary.hydromad"
    ans
}

print.summaryWithBreaks.hydromad <-
    function(x, digits = max(3, getOption("digits") - 3), ...)
{
    ## just simplify the printed output by rounding
    print.data.frame(x, digits=digits, ...)
    invisible(x)
}

print.summary.hydromad <-
    function(x, digits = max(3, getOption("digits") - 3), ...)
{
    with(x, {
        cat("\nCall:\n")
        print(call)
        cat("\n")
        if (!is.null(x$yic))
            cat("YIC:", format(yic, digits=digits), "\n")
        if (!is.null(x$arpe))
            cat("ARPE (%):", format(arpe * 100, digits=digits), "\n")
        if (!is.null(x$r.squared))
            cat("R Squared:", format(r.squared, digits=digits), "\n")
        if (!is.null(x$r.sq.sqrt))
            cat("R Squared sqrt:", format(r.sq.sqrt, digits=digits), "\n")
        if (!is.null(x$r.sq.log))
            cat("R Squared log:", format(r.sq.log, digits=digits), "\n")
        if (!is.null(x$r.sq.monthly))
            cat("R Squared monthly:", format(r.sq.monthly, digits=digits), "\n")
        if (!is.null(x$bias))
            cat("Bias (units, total):", format(bias, digits=digits), "\n")
        if (!is.null(x$rel.bias))
            cat("Rel. Bias (%):", format(rel.bias * 100, digits=digits), "\n")
        if (!is.null(x$U1))
            cat("Error corr. with lag input (U1):", format(U1, digits=digits), "\n")
        if (!is.null(x$X1))
            cat("Error corr. with lag output (X1):", format(X1, digits=digits), "\n")
        if (!is.null(x$ssg))
            cat("Steady State Gain:", format(ssg, digits=digits), "\n")
    })
    invisible(x)
}
