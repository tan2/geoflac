module myrandom_mod

contains 

!
! A simple thread-safe pseudo random number generator
!

subroutine myrandom(iseed, r)
    !$ACC routine seq
    integer, intent(inout) :: iseed
    double precision, intent(out) :: r

    iseed = mod(iseed, 134456)
    iseed = mod(8121 * iseed + 28411, 134456 )
    r = (iseed * 1.0d0) / 134456

    return
end subroutine myrandom


end module myrandom_mod
