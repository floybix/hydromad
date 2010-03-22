## generate the web site

## setwd("X:/Packages/hydromad/web")
## Sys.setenv(R_GSCMD = "C:/Program Files/gs/gs8.63/bin/gswin32c.exe")
## source("genHydromad.R", echo = TRUE)

library(hydromad)

spec <- list()

spec[["modelling framework"]] <-
    list(
         list("hydromad", example = 1),
         list("hydromad.object"),
         list("summary.hydromad"),
         list("simulate.hydromad"),
         list("runlist"),
         list("hydromad.options")
         )

spec[["soil moisture accounting"]] <-
    list(
         list("scalar"),
         list("IHACRES.CWI.model"),
         list("IHACRES.CMD.model"),
         list("bucket"),
         list("sacramento"),
         list("snow"),
         list("runoffratio"),
         list("dbm")
         )

spec[["calibration"]] <-
    list(
         list("objFunVal"),
         list("fitBySampling"),
         list("fitByOptim"),
         list("fitBySCE"),
         list("fitByDE"),
         #list("fitByDream")
         #list("mcmcByDream")
         list("tryModelOrders")
         )

spec[["routing"]] <-
    list(
         list("armax"),
         list("expuh"),
         list("lambda"),
         list("powuh")
         )

spec[["routing fitting"]] <-
    list(
         list("armax.ls.fit"),
         list("armax.sriv.fit"),
         list("armax.inverse.fit"),
         list("expuh.sriv.fit")
         #"lambda.inverse.fit"
         )

spec[["event-based analysis"]] <-
    list(
         list("eventseq"),
         list("eventAttributes")
         )

spec[["plot methods"]] <-
    list(
         list("xyplot", example = 1),
         list("xyplot.runlist", example = 1),
         list("errormasscurve", example = 1),
         list("qqmath", example = 1),
         list("rollccf", example = 1)
         )

spec[["utilities"]] <-
    list(
         list("convertFlow"),
         list("estimateDelay"),
         list("parameterSets"),
         list("fitStat"),
         list("observed")
         )

spec[["datasets"]] <-
    list(
         list("BinghamTrib", example = 1),
         list("Canning", example = 1),
         list("Cotter", example = 1),
         #list("Molonglo"),
         list("Murrindindi", example = 1),
         list("Queanbeyan", example = 1),
         list("SalmonBrook", example = 1),
         list("Wye", example = 1)
         )

spec[["wetlands"]] <-
    list(
         list("swimp", example = 1),
         list("poweroid", example = 2)
         )

source("http://latticeextra.r-forge.r-project.org/generate.R")

## stop on errors
lattice.options(panel.error = NULL)

generateWebsite("hydromad", spec = spec, 
                man.src.dir = "../man/",
                imageSrcBase = "",
                themeNames = "custom_theme_2"
                do.examples = TRUE)
