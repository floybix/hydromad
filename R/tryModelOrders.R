## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


tryModelOrders <-
    function(expr, n = 0:3, m = 0:2,
             delay = hydromad.getOption("delay"),
             verbose = hydromad.getOption("trace"))
{
    expr <- substitute(expr)
    allOrders <- expand.grid(m = m, n = n, delay = delay)
    allOrders <- subset(allOrders, n >= m)
    results <- list()
    beSilent <- hydromad.getOption("quiet")
    ## restore previous setting when finished
    oopt <- hydromad.options("order", "delay", "quiet")
    on.exit(hydromad.options(oopt))
    hydromad.options(quiet = !verbose)
    for (i in 1:NROW(allOrders)) {
        orderSpec <- allOrders[i,]
        order <- c(n = orderSpec$n, m = orderSpec$m)
        d <- orderSpec$delay
        hydromad.options(order = order, delay = d)
        if (any(!is.na(delay))) {
            hydromad.options(delay = d)
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
            if (!isFullySpecified(mod)) {
                message("model not fully specified!")
            } else if (!isValidModel(mod)) {
                message(toString(mod, width = 60))
            } else {
                modsumm <- summary(mod, which = c("yic", "r.squared"))
                YIC <- modsumm$yic
                R2 <- modsumm$r.squared
                message(" YIC = ", signif(YIC, 4),
                        ", R^2 = ", signif(R2, 4))
            }
        }
    }
    ans <- as.runlist(results)
    class(ans) <- c("tryModelOrders", class(ans))
    ans
}

summary.tryModelOrders <- 
    function(object,
             which = c("yic", "arpe", "r.squared", "r.squared.log", "ssg", "rel.bias"),
             ...)
{
    class(object) <- setdiff(class(object), "tryModelOrders")
    summary(object, which = which, ...)
}
