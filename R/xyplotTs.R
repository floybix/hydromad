##
## Copyright (c) 2008 Felix Andrews <felix@nfrac.org>
##


## for xyplot.zoo; type, col etc are lists
panel.superpose.cycle <- function(..., subscripts, groups, type, col, lty, lwd, pch) {
    inscreen <- seq_along(groups) %in% subscripts
    groups[!inscreen] <- NA
    groups <- groups[drop=TRUE]
    if (is.list(type)) type <- type[[1]]
    panel.superpose(..., subscripts=subscripts, groups=groups, type=type)
}
