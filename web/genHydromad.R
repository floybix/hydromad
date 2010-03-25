## generate the web site

## setwd("X:/Packages/hydromad/web")
## Sys.setenv(R_GSCMD = "C:/Program Files/gs/gs8.63/bin/gswin32c.exe")
## source("genHydromad.R", echo = TRUE)

library(hydromad)

spec <- list()

spec[["modelling framework"]] <-
    list(
         list("hydromad"),
         list("hydromad.object", do.example = FALSE),
         list("summary", helpname = "summary.hydromad", do.example = FALSE),
         list("simulate.hydromad", helpname = "simulate.hydromad"),
         list("runlist", do.example = FALSE),
         list("hydromad.options", do.example = FALSE)
         )

spec[["soil moisture accounting"]] <-
    list(
         list("scalar"),
         list("IHACRES.CWI.model"),
         list("IHACRES.CMD.model"),
         list("bucket"),
         list("awbm"),
         list("sacramento"),
         list("snow"),
         list("runoffratio"),
         list("dbm")
         )

spec[["calibration"]] <-
    list(
         list("objFunVal", do.example = FALSE),
         list("fitBySampling", do.example = FALSE),
         list("fitByOptim", do.example = FALSE),
         list("fitBySCE", do.example = FALSE),
         list("fitByDE", do.example = FALSE),
         #list("fitByDream")
         #list("mcmcByDream")
         list("tryModelOrders", do.example = FALSE)
         )

spec[["routing"]] <-
    list(
         list("armax"),
         list("expuh"),
         list("lambda", do.example = FALSE),
         list("powuh", do.example = FALSE)
         )

spec[["routing fitting"]] <-
    list(
         list("armax.ls.fit"),
         list("armax.sriv.fit"),
         list("armax.inverse.fit"),
         list("expuh.sriv.fit", do.example = FALSE)
         #"lambda.inverse.fit"
         )

spec[["event-based analysis"]] <-
    list(
         list("eventseq"),
         list("eventinfo", helpname = "eventseq")
         )

spec[["plot methods"]] <-
    list(
         list("xyplot", helpname = "xyplot.hydromad"),
         list("xyplot.runlist"),
         list("errormasscurve", helpname = "xyplot.hydromad"),
         list("qqmath", helpname = "xyplot.hydromad"),
         list("rollccf")
         )

spec[["utilities"]] <-
    list(
         list("convertFlow", do.example = FALSE),
         list("estimateDelay", do.example = FALSE),
         list("parameterSets", do.example = FALSE),
         list("fitStat", do.example = FALSE),
         list("observed", do.example = FALSE)
         )

spec[["datasets"]] <-
    list(
         list("BinghamTrib"),
         list("Canning"),
         list("Cotter"),
         #list("Molonglo"),
         list("Murrindindi"),
         list("Queanbeyan"),
         list("SalmonBrook"),
         list("Wye"),
         list("HydroTestData")
         )

spec[["wetlands"]] <-
    list(
         list("swimp"),
         list("poweroid")
         )

source("http://latticeextra.r-forge.r-project.org/generate.R")

## stop on errors
lattice.options(panel.error = NULL)

generateWebsite("hydromad", spec = spec, 
                man.src.dir = "../man/",
                imageSrcBase = "",
                themeNames = "custom_theme_2",
                do.examples = TRUE)
