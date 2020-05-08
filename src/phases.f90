
module phases
    implicit none

    !! Phases the got special treatment in the code

    integer, parameter :: kocean0 = 1 ! basalt without dehydration
    integer, parameter :: kcont1 = 2
    integer, parameter :: kocean1 = 3
    integer, parameter :: kmant1 = 4
    integer, parameter :: kmetased = 5
    integer, parameter :: kcont2 = 6
    integer, parameter :: kocean2 = 7
    integer, parameter :: kmant2 = 8
    integer, parameter :: kserp = 9
    integer, parameter :: ksed1 = 10
    integer, parameter :: ksed2 = 11
    integer, parameter :: kweak = 12
    integer, parameter :: keclg = 13
    integer, parameter :: karc1 = 14
    integer, parameter :: kweakmc = 15
    integer, parameter :: khydmant = 16

    !$ACC declare copyin(kocean0, kcont1, kocean1, kmant1, kmetased, kcont2, kocean2, &
    !$ACC                kmant2, kserp, ksed1, ksed2, kweak, keclg, karc1, kweakmc, &
    !$ACC                khydmant)

end module phases
