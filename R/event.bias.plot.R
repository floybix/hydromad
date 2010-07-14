

event.bias.plot <- function(...)
    UseMethod("event.bias.plot")

## TODO: how to allow log scales (or sqrt scales etc) for particular covariates?

event.bias.plot.default <- 
    function(resids, events, ..., FUN.resids = sum,
             vars.sum = list(), vars.mean = list(), vars.max = list(), 
             vars.first = list(), vars.ante = list(),
             vars.custom = list(), with.duration = FALSE,
             prefix = TRUE, sep = ".",
             auto.key = NCOL(resids) > 1, layout = NULL,
             abline = list(h = 0),
             xlab = "covariate values",
             ylab = "event-aggregated residuals")
{
    stopifnot(length(resids) > 0)
    stopifnot(length(events) > 0)
    if (prefix) {
        donames <- function(x, prefix) {
            if (length(x) > 0) {
                if (is.null(names(x)))
                    names(x) <- seq_along(x)
                names(x) <- paste(prefix, make.names(names(x), unique = TRUE), sep = sep)
            }
            x
        }
        vars.sum <- donames(vars.sum, "sum")
        vars.mean <- donames(vars.mean, "mean")
        vars.max <- donames(vars.max, "max")
        vars.first <- donames(vars.first, "first")
        vars.ante <- donames(vars.ante, "ante")
    }
    first <- function(x) head(x, 1)
    allvars <- 
        c(list(residuals = eventapply(resids, events, FUN = FUN.resids)),
          if (with.duration) 
          list(duration = eventapply(obs, events, FUN = length)),
          lapply(vars.sum, function(v) eventapply(v, events, FUN = sum)),
          lapply(vars.mean, function(v) eventapply(v, events, FUN = mean)),
          lapply(vars.max, function(v) eventapply(v, events, FUN = max)),
          lapply(vars.first, function(v) eventapply(v, events, FUN = first)),
          lapply(vars.ante, function(v) eventapply(lag(v, -1), events, FUN = first)),
          lapply(vars.custom, function(v) eventapply(v$data, events, FUN = v$fun))
          )
    cvnames <- names(allvars)[-1]
    ## merge all series -- if zoo series they will be synchronised
    allvars <- do.call("cbind", allvars)
    covars <- allvars[,-(1:NCOL(resids))]
    resids <- allvars[, (1:NCOL(resids))]
    if (length(covars) == 0) stop("no covariates to plot")
    colnames(covars) <- cvnames
    residstack <- stack(as.data.frame(resids))
    foo <-
        xyplot.list(as.data.frame(covars),
                    data = residstack,
                    FUN = function(VAR, ...)
                    xyplot(values ~ rep(VAR, NCOL(resids)), groups = ind, ...),
                    default.scales = list(x = list(relation = "free")),
                    auto.key = auto.key, layout = layout,
                    abline = abline,
                    xlab = xlab, ylab = ylab, ...)
    foo$call <- sys.call(sys.parent())
    foo
}

event.bias.plot.hydromad.runlist <- 
event.bias.plot.hydromad <- 
    function(model, events, ..., boxcox = FALSE, start = NULL)
{
    foo <- 
        event.bias.plot(residuals(model, boxcox = boxcox, start = start), 
                    events = events, ...)
    foo$call <- sys.call(sys.parent())
    foo
}
