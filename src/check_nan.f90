subroutine check_nan
    use arrays
    implicit none

    if (any(isnan(cord))) stop 'cord becomes NaN!'
end subroutine check_nan