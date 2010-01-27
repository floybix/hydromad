## ihacreslab: rainfall-runoff hydrology models and tools
##
## Copyright (c) 2008 Felix Andrews <felix@nfrac.org>
##

## local environment to store user options
.IhacresEnv <- new.env()
.IhacresEnv$options <- list()

.onLoad <- function(libname, pkgname)
{
    ihacres.options(.defaultIhacresOptions())
}
