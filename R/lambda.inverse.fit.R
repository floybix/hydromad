## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


lambda.inverse.fit <-
    function(DATA,
             ..., fit.method = "lambda") #with.lambda = TRUE)
{
    obj <- tf.inverse.fit(DATA, ..., fit.method = fit.method)#with.lambda = with.lambda)
    obj$call <- match.call()
    obj
}
