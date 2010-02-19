

## Time series comparison on corresponding times
ts_equals <- function(expected, ..., start = NULL, end = NULL, trim = FALSE) {
    function(actual) {
        windowts <- window(ts.intersect(expected, actual),
                           start = start, end = end)
        if (trim)
            windowts <- na.trim(windowts)
        equals(windowts[,1], ...)(windowts[,2])
    }
}
