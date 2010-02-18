

## Time series comparison on corresponding times
equals_in_window <- function(expected, ..., start = NULL, end = NULL) {
    function(actual) {
        windowts <- window(ts.intersect(expected, actual),
                           start = start, end = end)
        equals(windowts[,1], ...)(windowts[,2])
    }
}
