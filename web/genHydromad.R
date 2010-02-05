## generate the web site

## setwd("X:/Packages/hydromad/web")
## Sys.setenv(R_GSCMD = "C:/Program Files/gs/gs8.63/bin/gswin32c.exe")
## source("genHydromad.R", echo = TRUE)

library(hydromad)

spec <- list()

spec[["modelling framework"]] <-
    list(
         "hydromad" = list(example = 1),
         "hydromad.object" = list(),
         #"predict.hydromad" = list(),
         #"summary.hydromad" = list(),
         "simulate.hydromad" = list(),
         "runlist" = list(),
         "hydromad.options" = list()
         )

spec[["soil moisture accounting"]] <-
    list(
         "scalar" = list(),
         "IHACRES.CWI.model" = list(),
         "IHACRES.CMD.model" = list(),
         "bucket" = list(),
         "sacramento" = list()
         )

spec[["calibration"]] <-
    list(
         "fitBySampling" = list(),
         "fitByOptim" = list(),
         "fitBySCE" = list(),
         "fitByDE" = list()
         )

spec[["routing"]] <-
    list(
         "expuh" = list(),
         "uh" = list()
         )

spec[["routing fitting"]] <-
    list(
         "uh.ls.fit" = list(),
         "uh.sriv.fit" = list(),
         "uh.inverse.fit" = list()
         )

spec[["event-based analysis"]] <-
    list(
         "event.clusters" = list(),
         "eventAttributes" = list()
         )

spec[["plot methods"]] <-
    list(
         "xyplot" = list(example = 1),
         "xyplot.runlist" = list(example = 1),
         "errormasscurve" = list(example = 1),
         "qqmath" = list(example = 1),
         "rollccf" = list(example = 1)
         )

spec[["utilities"]] <-
    list(
         "convertFlow" = list(),
         "estimateDelay" = list(),
         "tryModelOrders" = list(),
         "parameterSets" = list(),
         "fitStat" = list(),
         "objFunVal" = list(),
         "observed" = list()
         )

spec[["datasets"]] <-
    list(
         "BinghamTrib" = list(example = 1),
         "Canning" = list(example = 1),
         "Cotter" = list(example = 1),
         ##"Molonglo" = list(),
         "Murrindindi" = list(example = 1),
         "Queanbeyan" = list(example = 1),
         "SalmonBrook" = list(example = 1),
         "Wye" = list(example = 1)
         )

spec[["wetlands"]] <-
    list(
         "swimp" = list(example = 1),
         "poweroid" = list(example = 2)
         )

source("http://latticeextra.r-forge.r-project.org/generate.R")

## stop on errors
lattice.options(panel.error = NULL)

generateWebsite("hydromad", spec = spec, 
                man.src.dir = "../man/",
                imageSrcBase = "",
                themeNames = "custom_theme_2"
                do.examples = TRUE)
