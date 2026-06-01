subroutine fl_flexure
    use arrays
    use params
    implicit none

    double precision :: q(nx), w_flex(nx), Load(nx)
    double precision :: a_block(2, 2, nx), b_block(2, 2, nx), c_block(2, 2, nx)
    double precision :: r_block(2, nx), u_block(2, nx)
    double precision :: x(nx), dx(nx)
    double precision :: D(nx), E, nu, Te, rho_m, lambda, mu
    double precision :: cord_old_z, w_winkler_1, w_winkler_nx
    integer :: i, j, iph, bc_idx
    double precision :: a_i, c_i, b_i

    ! 1. Scan for the flexural boundary condition (nbc = 200) on side 2 (bottom)
    bc_idx = 0
    do i = 1, nofbc
        if (nbc(i) .eq. 200 .and. nofside(i) .eq. 2) then
            bc_idx = i
            exit
        endif
    enddo

    if (bc_idx .eq. 0) return  ! Flexural boundary condition not active

    Te = bca(bc_idx)       ! Elastic thickness in meters
    rho_m = bcb(bc_idx)    ! Asthenosphere density in kg/m3

    ! 2. Compute current vertical column weights (lithostatic loads)
    call get_column_weights(q)

    ! 3. If q_init is completely uninitialized (all zeros), initialize it with the starting column weights
    if (all(q_init .eq. 0.d0)) then
        q_init = q
    endif

    ! Compute excess load
    !$OMP parallel do
    do i = 1, nx
        Load(i) = q(i) - q_init(i)
    enddo

    ! 4. Set up grid coordinates and spacings in X
    !$OMP parallel do
    do i = 1, nx
        x(i) = cord(nz, i, 1)
    enddo
    !$OMP parallel do
    do i = 1, nx-1
        dx(i) = x(i+1) - x(i)
        if (dx(i) .le. 0.d0) dx(i) = 1.d0  ! Safeguard against zero/negative spacing
    enddo

    ! 5. Compute local flexural rigidity D(i) based on Lamé parameters of the bottom elements
    !$OMP parallel do private(iph, lambda, mu, E, nu)
    do i = 1, nx
        if (i .lt. nx) then
            iph = iphase(nz-1, i)
        else
            iph = iphase(nz-1, nx-1)
        endif

        ! Ensure phase is within valid range
        if (iph .lt. 1 .or. iph .gt. nphase) iph = 1

        lambda = rl(iph)
        mu = rm(iph)

        if (lambda + mu .gt. 0.d0) then
            E = mu * (3.0d0 * lambda + 2.0d0 * mu) / (lambda + mu)
            nu = lambda / (2.0d0 * (lambda + mu))
        else
            ! Fallback to standard mantle values if Lamé parameters are not set
            E = 1.0d11
            nu = 0.25d0
        endif

        if (nu .ge. 1.d0 .or. nu .le. -1.d0) nu = 0.25d0

        D(i) = (E * (Te**3)) / (12.0d0 * (1.0d0 - (nu**2)))
        if (D(i) .le. 0.d0) D(i) = 1.d15  ! Safeguard against zero/negative rigidity
    enddo

    ! 6. Construct the Block Tridiagonal System of Equations
    ! Boundary condition at left end: M_1 = 0, w_1 = q(1) / (rho_m * g)
    a_block(:, :, 1) = 0.d0
    c_block(:, :, 1) = 0.d0
    b_block(1, 1, 1) = 1.d0
    b_block(1, 2, 1) = 0.d0
    b_block(2, 1, 1) = 0.d0
    b_block(2, 2, 1) = 1.d0
    
    r_block(1, 1) = 0.d0
    if (rho_m * g .gt. 0.d0) then
        w_winkler_1 = Load(1) / (rho_m * g)
    else
        w_winkler_1 = 0.d0
    endif
    r_block(2, 1) = w_winkler_1

    ! Inner nodes using non-uniform 2nd-order central differences
    !$OMP parallel do private(a_i, c_i, b_i)
    do i = 2, nx-1
        a_i = 2.d0 / (dx(i-1) * (dx(i) + dx(i-1)))
        c_i = 2.d0 / (dx(i) * (dx(i) + dx(i-1)))
        b_i = -(a_i + c_i)

        ! Diagonal block elements
        a_block(1, 1, i) = a_i
        a_block(1, 2, i) = 0.d0
        a_block(2, 1, i) = 0.d0
        a_block(2, 2, i) = a_i

        c_block(1, 1, i) = c_i
        c_block(1, 2, i) = 0.d0
        c_block(2, 1, i) = 0.d0
        c_block(2, 2, i) = c_i

        b_block(1, 1, i) = b_i
        b_block(1, 2, i) = -rho_m * g
        b_block(2, 1, i) = 1.0d0 / D(i)
        b_block(2, 2, i) = b_i

        r_block(1, i) = -Load(i)
        r_block(2, i) = 0.d0
    enddo

    ! Boundary condition at right end: M_nx = 0, w_nx = q(nx) / (rho_m * g)
    a_block(:, :, nx) = 0.d0
    c_block(:, :, nx) = 0.d0
    b_block(1, 1, nx) = 1.d0
    b_block(1, 2, nx) = 0.d0
    b_block(2, 1, nx) = 0.d0
    b_block(2, 2, nx) = 1.d0
    
    r_block(1, nx) = 0.d0
    if (rho_m * g .gt. 0.d0) then
        w_winkler_nx = Load(nx) / (rho_m * g)
    else
        w_winkler_nx = 0.d0
    endif
    r_block(2, nx) = w_winkler_nx

    ! 7. Solve the block tridiagonal system directly
    call solve_block_tridiagonal(nx, a_block, b_block, c_block, r_block, u_block)

    ! 8. Apply new deflections, update bottom boundary coordinates & velocities
    w_flex = u_block(2, :)
    !$OMP parallel do private(cord_old_z)
    do i = 1, nx
        cord_old_z = cord(nz, i, 2)
        cord(nz, i, 2) = zoriginal(nz, i) - w_flex(i)
        vel(nz, i, 2) = (cord(nz, i, 2) - cord_old_z) / dt
    enddo

    !$ACC update device(cord, vel)

end subroutine fl_flexure


subroutine get_column_weights(q)
    use arrays
    use params
    implicit none
    double precision, intent(out) :: q(nx)
    double precision :: rho, dz
    integer :: i, j
    double precision, external :: Eff_dens

    !$OMP parallel do private(rho, dz, j)
    do i = 1, nx
        q(i) = 0.d0
        do j = 1, nz-1
            if (i .eq. 1) then
                rho = Eff_dens(j, 1)
            elseif (i .eq. nx) then
                rho = Eff_dens(j, nx-1)
            else
                rho = 0.5d0 * (Eff_dens(j, i-1) + Eff_dens(j, i))
            endif
            dz = cord(j, i, 2) - cord(j+1, i, 2)
            q(i) = q(i) + rho * g * dz
        enddo
    enddo
end subroutine get_column_weights


subroutine solve_block_tridiagonal(n, a, b, c, r, u)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: a(2, 2, n), b(2, 2, n), c(2, 2, n), r(2, n)
    double precision, intent(out) :: u(2, n)
    
    double precision :: cp(2, 2, n), rp(2, n)
    double precision :: denom(2, 2), inv_denom(2, 2), det
    double precision :: tmp_v(2)
    integer :: i
    
    ! For i = 1: B_1 is identity, so B_1^-1 is identity.
    cp(:, :, 1) = c(:, :, 1)
    rp(:, 1) = r(:, 1)
    
    ! Forward sweep:
    do i = 2, n-1
        ! denom = B_i - A_i * cp_{i-1}
        denom(1, 1) = b(1, 1, i) - a(1, 1, i) * cp(1, 1, i-1)
        denom(1, 2) = b(1, 2, i) - a(1, 1, i) * cp(1, 2, i-1)
        denom(2, 1) = b(2, 1, i) - a(2, 2, i) * cp(2, 1, i-1)
        denom(2, 2) = b(2, 2, i) - a(2, 2, i) * cp(2, 2, i-1)
        
        ! Invert denom
        det = denom(1, 1) * denom(2, 2) - denom(1, 2) * denom(2, 1)
        if (abs(det) .lt. 1.d-30) det = sign(1.d-30, det)
        
        inv_denom(1, 1) = denom(2, 2) / det
        inv_denom(1, 2) = -denom(1, 2) / det
        inv_denom(2, 1) = -denom(2, 1) / det
        inv_denom(2, 2) = denom(1, 1) / det
        
        ! cp_i = inv_denom * c_i
        cp(1, 1, i) = inv_denom(1, 1) * c(1, 1, i)
        cp(1, 2, i) = inv_denom(1, 2) * c(2, 2, i)
        cp(2, 1, i) = inv_denom(2, 1) * c(1, 1, i)
        cp(2, 2, i) = inv_denom(2, 2) * c(2, 2, i)
        
        ! rp_i = inv_denom * (r_i - A_i * rp_{i-1})
        tmp_v(1) = r(1, i) - a(1, 1, i) * rp(1, i-1)
        tmp_v(2) = r(2, i) - a(2, 2, i) * rp(2, i-1)
        
        rp(1, i) = inv_denom(1, 1) * tmp_v(1) + inv_denom(1, 2) * tmp_v(2)
        rp(2, i) = inv_denom(2, 1) * tmp_v(1) + inv_denom(2, 2) * tmp_v(2)
    enddo
    
    ! For i = n: A_n is 0, denom = B_n = identity.
    u(:, n) = r(:, n)
    
    ! Backward substitution:
    do i = n-1, 1, -1
        u(1, i) = rp(1, i) - (cp(1, 1, i) * u(1, i+1) + cp(1, 2, i) * u(2, i+1))
        u(2, i) = rp(2, i) - (cp(2, 1, i) * u(1, i+1) + cp(2, 2, i) * u(2, i+1))
    enddo
end subroutine solve_block_tridiagonal
