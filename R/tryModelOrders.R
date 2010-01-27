## ihacreslab: rainfall-runoff hydrology models and tools
##
## Copyright (c) 2008 Felix Andrews <felix@nfrac.org>
##


tryModelOrders <-
    function(expr, n = 0:3, m = 0:2,
             delay = ihacres.getOption("delay"),
             verbose = ihacres.getOption("trace"))
{
    expr <- substitute(expr)
    allOrders <- expand.grid(m = m, n = n, delay = delay)
    allOrders <- subset(allOrders, n >= m)
    results <- list()
    beSilent <- ihacres.getOption("quiet")
    ## restore previous setting when finished
    oopt <- ihacres.options("order", "delay", "quiet")
    on.exit(ihacres.options(oopt))
    ihacres.options(quiet = !verbose)
    for (i in 1:NROW(allOrders)) {
        orderSpec <- allOrders[i,]
        order <- c(n = orderSpec$n, m = orderSpec$m)
        d <- orderSpec$delay
        ihacres.options(order = order, delay = d)
        if (any(!is.na(delay))) {
            ihacres.options(delay = d)
            nm <- paste("(n=", order[1], ", m=", order[2],
                        ", d=", d, ")", sep = "")
        } else {
            nm <- paste("(n=", order[1], ", m=", order[2], ")", sep = "")
        }
        if (!isTRUE(beSilent))
            message("order: ", nm, "... ", appendLF = FALSE)
        ## execute the expression and store the result
        mod <- eval.parent(expr)
        results[[nm]] <- mod
        if (!isTRUE(beSilent)) {
            if (!isValidModel(mod)) {
                message(toString(mod))
            } else {
                modsumm <- summary(mod, which = c("yic", "r.squared"))
                YIC <- modsumm$yic
                R2 <- modsumm$r.squared
                message(" YIC = ", signif(YIC, 4),
                        ", R^2 = ", signif(R2, 4))
            }
        }
    }
    as.runlist(results)
}


ordersSummary <-
    function(x, flowstats = c("yic", "arpe", "r.squared", "r.squared.log", "ssg", "rel.bias"), ...)
{
    summary(x, pars = FALSE, flowstats = flowstats, ...)
}
