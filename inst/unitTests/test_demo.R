test_divideBy <- function() {
    checkEquals(4/ 2, 2)
    checkTrue(is.na(4/ 0))
    checkEqualsNumeric(4 /1.2345, 3.24, tolerance=1.0e-4)
}
