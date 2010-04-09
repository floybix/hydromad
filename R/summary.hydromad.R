## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


summary.hydromad.runlist <-
    function(object, ...,
             stats = c("rel.bias", "r.squared", "r.sq.sqrt", "r.sq.log", "r.sq.monthly"),
             items = stats)
{
    summary.runlist(object, ...,
                    stats = stats,
                    items = items)
}

summary.hydromad <-
    function(object, breaks = NULL,
             coerce = byDays,
             stats = c("rel.bias", "r.squared", "r.sq.sqrt", "r.sq.log", "r.sq.monthly"),
             na.action = na.exclude,
             ...)
{
    if (!isFullySpecified(object))
        stop("model parameters not fully specified")
    if (!isValidModel(object))
        stop("model is not valid")

    hydrostats <- function(chunk) {
        ans <- numeric()
        ok <- complete.cases(chunk)
        meanP <- mean(chunk[ok, "P"])
        meanQ <- mean(chunk[ok, "Q"])
        ans <-
            c(timesteps = NROW(chunk),
              missing = sum(!ok),
              mean.P = meanP,
              mean.Q = meanQ,
              runoff = meanQ / meanP)
        ans
    }

    DATA <- ts.intersect(P = observed(object, item = "P", all = TRUE),
                         Q = observed(object, all = TRUE),
                         U = fitted(object, U = TRUE, all = TRUE),
                         X = fitted(object, all = TRUE))

    if (!is.null(breaks)) {
        nms <- colnames(DATA)
        DATA <- coerce(DATA)
        colnames(DATA) <- nms
        group <- cut(time(DATA), breaks = breaks)
        group <- factor(group)
        ## remove warmup
        DATA <- stripWarmup(DATA, object$warmup)
        if (object$warmup > 0)
            group <- group[-seq_len(object$warmup), drop = TRUE]

        ans <-
            eventapply(DATA, group, 
                       FUN = function(x, ...)
                         c(hydrostats(x), perfStats(x, ...)),
                       stats = stats, ...,
                       by.column = FALSE)
        ## copy the last entry with the final date, to mark the end of last period
        lastbit <- tail(ans, 1)
        time(lastbit) <- end(DATA)
        ans <- rbind(ans, lastbit)
        return(ans)
    }

    ans <- list(call = object$call)

    ## ARPE
    if ("yic" %in% stats) {
        stats <- union(stats, "arpe")
    }
    if ("arpe" %in% stats) {
        if (is.null(vcov(object))) {
            ans$arpe <- NA_real_
        } else {
            coef.var <- diag(vcov(object))
            cc <- coef(object, "routing")
            cc2 <- tfParsConvert(cc, "a,b")
            cc[names(cc2)] <- cc2
            nms <- intersect(names(cc), names(coef.var))
            ans$arpe <- mean(coef.var[nms] / (cc[nms]^2))
            if ("yic" %in% stats) {
                var.ratio <- (var(residuals(object), na.rm = TRUE) /
                              var(observed(object), na.rm = TRUE))
                ans$yic <- log(var.ratio) + log(ans$arpe)
            }
        }
    }

    ans <- c(ans, hydrostats(DATA))
    
    ## call perfStats for the rest
    ans <- c(ans, perfStats(DATA, warmup = object$warmup, stats = stats, ...))
    
    class(ans) <- "summary.hydromad"
    ans
}

print.summary.hydromad <-
    function(x, digits = max(3, getOption("digits") - 3), ...)
{
    with(x, {
        cat("\nCall:\n")
        print(call)
        cat("\n")
        if (!is.null(x$timesteps) && !is.null(x$missing))
            cat("Time steps: ", x$timesteps, " (", x$missing, " missing).\n")
        if (!is.null(x$mean.P) && !is.null(x$mean.Q))
            cat("Runoff ratio (Q/P): (",
                format(x$mean.Q, digits=digits), " / ",
                format(x$mean.P, digits=digits), ") = ",
                format(x$mean.Q / x$mean.P, digits=digits), "\n")
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

print.summaryWithBreaks.hydromad <-
    function(x, digits = max(3, getOption("digits") - 3), ...)
{
    ## just simplify the printed output by rounding
    NextMethod("print")
}
