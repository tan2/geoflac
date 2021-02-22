subroutine check_nan
    use ieee_arithmetic
    use arrays
    implicit none

    !if (any(ieee_is_nan(cord))) stop 'cord becomes NaN!'
end subroutine check_nan
