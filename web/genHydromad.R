## generate the web site

## setwd("X:/Packages/hydromad/web")
## Sys.setenv(R_GSCMD = "C:/Program Files/gs/gs8.63/bin/gswin32c.exe")
## source("genHydromad.R", echo = TRUE)

library(hydromad)

spec <- list()

spec[["modelling framework"]] <-
    list(
         list("hydromad"),
         list("methods...", helpname = "hydromad.object", do.example = FALSE,
              codefile = "coef.hydromad.R"),
         list("summary", helpname = "summary.hydromad", do.example = FALSE),
         list("simulate", helpname = "simulate.hydromad"),
         list("runlist", do.example = FALSE),
         list("hydromad.options", do.example = FALSE, codefile = "options.R")
         )

spec[["soil moisture accounting"]] <-
    list(
         list("scalar"),
         list("cwi", helpname = "IHACRES.CWI.model", codefile = "cwi.R"),
         list("cmd", helpname = "IHACRES.CMD.model", codefile = "cmd.R"),
         list("bucket"),
         list("awbm"),
         list("sacramento"),
         list("snow"),
         list("runoffratio"),
         list("dbm")
         )

spec[["calibration"]] <-
    list(
         list("objFunVal", codefile = "fitStat.R", do.example = FALSE),
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
         list("armax", helpname = "armax.sim"),
         list("expuh", helpname = "expuh.sim"),
         list("lambda", helpname = "lambda.sim", do.example = FALSE),
         list("powuh", helpname = "powuh.sim", do.example = FALSE)
         )

spec[["routing fitting"]] <-
    list(
         list("armax.ls.fit"),
         list("armax.sriv.fit"),
         list("armax.inverse.fit")
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
         list("BinghamTrib", codefile = "../data/BinghamTrib.R"),
         list("Canning", codefile = NA),
         list("Cotter", codefile = NA),
         #list("Molonglo"),
         list("Murrindindi", codefile = NA),
         list("Queanbeyan", codefile = NA),
         list("SalmonBrook", codefile = NA),
         list("Wye", codefile = NA),
         list("HydroTestData", codefile = "../data/HydroTestData.R")
         )

spec[["wetlands"]] <-
    list(
         list("swimp"),
         list("poweroid")
         )

source("http://latticeextra.r-forge.r-project.org/generate.R")
source("X:/Packages/latticeextra/www/generate.R")

## stop on errors
lattice.options(panel.error = NULL)

generateWebsite("hydromad", spec = spec, 
                man.src.dir = "../man/",
                imageSrcBase = "",
                codeSrcBase = "http://github.com/floybix/hydromad/tree/master/R/",
                themeNames = "custom_theme_2",
                do.examples = TRUE)
