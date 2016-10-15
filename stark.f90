function eval(nmax,C,k,s0tau,shift) result(res)
    implicit none
    integer, parameter :: dp = kind(1.d0)
    integer, intent(in) :: nmax, k, shift
    integer :: n
    real(dp),intent(in) :: s0tau,C(0:nmax, *)

    real(dp) :: res, a

    res = 0.0_dp
    a = 1.0_dp
    do n = 0, nmax-1-shift
        res = res + C(n+shift, k) * a
        a = a * s0tau / (n + 1)
    end do
    res = res + C(nmax, k) * a
end function eval

subroutine evaldr0(rv0,vv0,pv,nmax,q,C,k,s0tau,alpha,shift,dadr0)
    implicit none
    integer, parameter :: dp = kind(1.d0)
    integer, intent(in) :: nmax, k, shift
    integer :: n
    real(dp),intent(in) :: rv0(3), vv0(3), pv(3), alpha
    real(dp),intent(in) :: q,s0tau,C(0:nmax, *)
    real(dp),intent(out) :: dadr0(3)

    real(dp) :: a
    real(dp) :: dqdz(3), dsdz(3), drdz(3)

    dqdz = -q**3 * rv0
    dsdz = vv0
    drdz = pv

    dadr0 = 0.0_dp
    a = 1.0_dp
    do n = 0, nmax-1-shift
        dadr0 = dadr0 + a * (&
            C(n+shift, 6*k-1) * dqdz + C(n+shift, 6*k) * dsdz + C(n+shift, 6*k+1) * drdz + &
            alpha * s0tau * q**2 * rv0 * C(n+shift+1, k))
        a = a * s0tau / (n + 1)
    end do
end subroutine evaldr0

subroutine evaldv0(rv0,vv0,pv,nmax,C,k,s0tau,shift,dadv0)
    implicit none
    integer, parameter :: dp = kind(1.d0)
    integer, intent(in) :: nmax, k, shift
    integer :: n
    real(dp),intent(in) :: rv0(3), vv0(3), pv(3)
    real(dp),intent(in) :: s0tau,C(0:nmax, *)
    real(dp),intent(out) :: dadv0(3)

    real(dp) :: a
    real(dp) :: dsdz(3), dwdz(3), dvdz(3)

    dsdz = rv0
    dwdz = pv
    dvdz = 2*vv0

    dadv0 = 0.0_dp
    a = 1.0_dp
    do n = 0, nmax-1-shift
        dadv0 = dadv0 + a * (C(n+shift, 6*k) * dsdz + C(n+shift, 6*k+2) * dwdz&
            + C(n+shift, 6*k+3) * dvdz)
        a = a * s0tau / (n + 1)
    end do
end subroutine evaldv0

subroutine evaldp(rv0,vv0,pv,nmax,C,k,s0tau,shift,dadp)
    implicit none
    integer, parameter :: dp = kind(1.d0)
    integer, intent(in) :: nmax, k, shift
    integer :: n
    real(dp),intent(in) :: rv0(3), vv0(3), pv(3)
    real(dp),intent(in) :: s0tau,C(0:nmax, *)
    real(dp),intent(out) :: dadp(3)

    real(dp) :: a
    real(dp) :: drdz(3), dwdz(3), dpdz(3)

    drdz = rv0
    dwdz = vv0
    dpdz = 2*pv

    dadp = 0.0_dp
    a = 1.0_dp
    do n = 0, nmax-1-shift
        dadp = dadp + a * (C(n+shift, 6*k+1) * drdz + C(n+shift, 6*k+2) * dwdz&
            + C(n+shift, 6*k+4) * dpdz)
        a = a * s0tau / (n + 1)
    end do
end subroutine evaldp

function det3(a, b, c)
    implicit none
    integer, parameter :: dp = kind(1.d0)
    real(dp), intent(in) :: a(3), b(3), c(3)
    real(dp) :: det3
    det3 = a(3)*(-(b(2)*c(1)) + b(1)*c(2)) + a(2)*(b(3)*c(1) - b(1)*c(3)) + &
    a(1)*(-(b(3)*c(2)) + b(2)*c(3))
end function det3

! Cost
! FG x 1
! eval x 7
! 7 x eval + 1 x coeff
subroutine FandG(rv0, vv0, p, tau, nmax, rv, vv, t, integr)
    implicit none
    integer, parameter :: dp = kind(1.d0)
    real(dp), parameter :: alpha = 1
    real(dp), intent(in) :: rv0(3), vv0(3), p(3), tau
    integer, intent(in) :: nmax
    real(dp) :: q,sigma,rho,omega,nu,theta,C(0:nmax, 4),r0,s0
    real(dp) :: s0tau, fr, gr, hr, fv, gv, hv, r, s, gam
    real(dp) :: eval
    real(dp), intent(out) :: rv(3), vv(3), t
    real(dp), intent(out) :: integr(3, 2)
    real(dp) :: ham0, ham, Lp0, Lp, Ap0, Ap, det3

    r0 = norm2(rv0)
    q = 1.0_dp / r0
    sigma = dot_product(rv0, vv0)
    rho = dot_product(rv0, p)
    omega = dot_product(vv0, p)
    nu = dot_product(vv0, vv0)
    theta = dot_product(p, p)

    ham0 = nu / 2 - q - rho
    Lp0 = det3(rv0, vv0, p)
    Ap0 = rho * (-nu + q + rho / 2) + omega * sigma - r0**2 * theta / 2

    s0 = r0**alpha
    s0tau = s0 * tau

    call FG(nmax,q,sigma,rho,omega,nu,theta,C)

    fr = eval(nmax,C,1,s0tau,0)
    gr = eval(nmax,C,2,s0tau,0)
    hr = eval(nmax,C,3,s0tau,0)
    t  = eval(nmax,C,4,s0tau,0)
    fv = eval(nmax,C,1,s0tau,1)
    gv = eval(nmax,C,2,s0tau,1)
    hv = eval(nmax,C,3,s0tau,1)

    rv = fr * rv0 + gr * vv0 + hr * p
    r = norm2(rv)
    s = r**alpha
    gam = s0 / s
    vv = gam * (fv * rv0 + gv * vv0 + hv * p)

    q = 1.0_dp / r
    sigma = dot_product(rv, vv)
    rho = dot_product(rv, p)
    omega = dot_product(vv, p)
    nu = dot_product(vv, vv)

    ham = nu / 2 - q - rho
    Lp = det3(rv, vv, p)
    Ap = rho * (-nu + q + rho / 2) + omega * sigma - r**2 * theta / 2

    integr(1, 1) = ham0
    integr(2, 1) = Lp0
    integr(3, 1) = Ap0
    integr(1, 2) = ham
    integr(2, 2) = Lp
    integr(3, 2) = Ap
end subroutine FandG

subroutine outer(a, b, c)
    implicit none
    integer, parameter :: dp = kind(1.d0)
    real(dp),intent(in) :: a(3), b(3)
    real(dp),intent(out) :: c(3, 3)
    integer :: i,j
    do i=1,3
        do j=1,3
            c(i, j) = a(i) * b(j)
        end do
    end do
end subroutine outer

! Cost
! FGex x 1
! eval x 7
! evaldr0 x 7
! evaldv0 x 7
! evaldp x 7
! 28 x eval + 1 x ex. coeff
!
! Approx. 4-5 times slower than FandG
subroutine FandGex(rv0, vv0, p, tau, nmax, &
        rv, vv, t, &
        drvdtau, dvvdtau, dtdtau, &
        drvdrv0, dvvdrv0, dtdrv0, &
        drvdvv0, dvvdvv0, dtdvv0, &
        drvdp, dvvdp, dtdp, integr)
    implicit none
    integer, parameter :: dp = kind(1.d0)
    real(dp), parameter :: alpha = 1
    real(dp), intent(in) :: rv0(3), vv0(3), p(3), tau
    integer, intent(in) :: nmax
    real(dp) :: q,sigma,rho,omega,nu,theta,C(0:nmax, 28),r0,s0
    real(dp) :: s0tau, fr, gr, hr, fv, gv, hv, r, s, gam
    real(dp) :: eval
    real(dp), intent(out) :: rv(3), vv(3), t
    real(dp), intent(out) :: drvdtau(3), dvvdtau(3), dtdtau
    real(dp), intent(out) :: drvdrv0(3, 3), dvvdrv0(3, 3), dtdrv0(3)
    real(dp), intent(out) :: drvdvv0(3, 3), dvvdvv0(3, 3), dtdvv0(3)
    real(dp), intent(out) :: drvdp(3, 3), dvvdp(3, 3), dtdp(3)
    real(dp), intent(out) :: integr(3, 2)
    real(dp) :: dtdz(3), dfrdz(3), dgrdz(3), dhrdz(3)
    real(dp) :: dfvdz(3), dgvdz(3), dhvdz(3)
    real(dp) :: dlngamdz(3), dr0dz(3), drdz(3)
    real(dp) :: tmp(3, 3)
    real(dp) :: ham0, ham, Lp0, Lp, Ap0, Ap, det3
    integer :: i

    r0 = norm2(rv0)
    q = 1.0_dp / r0
    sigma = dot_product(rv0, vv0)
    rho = dot_product(rv0, p)
    omega = dot_product(vv0, p)
    nu = dot_product(vv0, vv0)
    theta = dot_product(p, p)

    ham0 = nu / 2 - q - rho
    Lp0 = det3(rv0, vv0, p)
    Ap0 = rho * (-nu + q + rho / 2) + omega * sigma - r0**2 * theta / 2

    s0 = r0**alpha
    s0tau = s0 * tau

    call FGex(nmax,q,sigma,rho,omega,nu,theta,C)

    fr = eval(nmax,C,1,s0tau,0)
    gr = eval(nmax,C,2,s0tau,0)
    hr = eval(nmax,C,3,s0tau,0)
    t  = eval(nmax,C,4,s0tau,0)
    fv = eval(nmax,C,1,s0tau,1)
    gv = eval(nmax,C,2,s0tau,1)
    hv = eval(nmax,C,3,s0tau,1)

    rv = fr * rv0 + gr * vv0 + hr * p
    r = norm2(rv)
    s = r**alpha
    gam = s0 / s
    vv = gam * (fv * rv0 + gv * vv0 + hv * p)

    ! Derivatives wrt tau

    drvdtau = s * vv
    dvvdtau = s * (p - rv / r**3)
    dtdtau = s

    ! Derivatives wrt z = rv0

    call evaldr0(rv0,vv0,p,nmax,q,C,1,s0tau,alpha,0,dfrdz)
    call evaldr0(rv0,vv0,p,nmax,q,C,2,s0tau,alpha,0,dgrdz)
    call evaldr0(rv0,vv0,p,nmax,q,C,3,s0tau,alpha,0,dhrdz)
    call evaldr0(rv0,vv0,p,nmax,q,C,4,s0tau,alpha,0,dtdz)

    dtdrv0 = dtdz

    call outer(rv0, dfrdz, tmp)
    drvdrv0 = tmp
    call outer(vv0, dgrdz, tmp)
    drvdrv0 = drvdrv0 + tmp
    call outer(p, dhrdz, tmp)
    drvdrv0 = drvdrv0 + tmp
    do i = 1,3
        drvdrv0(i, i) = drvdrv0(i, i) + fr
    end do

    call evaldr0(rv0,vv0,p,nmax,q,C,1,s0tau,alpha,1,dfvdz)
    call evaldr0(rv0,vv0,p,nmax,q,C,2,s0tau,alpha,1,dgvdz)
    call evaldr0(rv0,vv0,p,nmax,q,C,3,s0tau,alpha,1,dhvdz)

    dr0dz = rv0 * q
    drdz = matmul(rv / r, drvdrv0)
    dlngamdz = alpha * (dr0dz * q - drdz / r)

    call outer(rv0, dfvdz, tmp)
    dvvdrv0 = gam * tmp
    call outer(vv0, dgvdz, tmp)
    dvvdrv0 = dvvdrv0 + gam * tmp
    call outer(p, dhvdz, tmp)
    dvvdrv0 = dvvdrv0 + gam * tmp
    call outer(vv, dlngamdz, tmp)
    dvvdrv0 = dvvdrv0 + tmp
    do i = 1,3
        dvvdrv0(i, i) = dvvdrv0(i, i) + gam * fv
    end do

    ! Derivatives wrt z = vv0

    call evaldv0(rv0,vv0,p,nmax,C,1,s0tau,0,dfrdz)
    call evaldv0(rv0,vv0,p,nmax,C,2,s0tau,0,dgrdz)
    call evaldv0(rv0,vv0,p,nmax,C,3,s0tau,0,dhrdz)
    call evaldv0(rv0,vv0,p,nmax,C,4,s0tau,0,dtdz)

    dtdvv0 = dtdz

    call outer(rv0, dfrdz, tmp)
    drvdvv0 = tmp
    call outer(vv0, dgrdz, tmp)
    drvdvv0 = drvdvv0 + tmp
    call outer(p, dhrdz, tmp)
    drvdvv0 = drvdvv0 + tmp
    do i = 1,3
        drvdvv0(i, i) = drvdvv0(i, i) + gr
    end do

    call evaldv0(rv0,vv0,p,nmax,C,1,s0tau,1,dfvdz)
    call evaldv0(rv0,vv0,p,nmax,C,2,s0tau,1,dgvdz)
    call evaldv0(rv0,vv0,p,nmax,C,3,s0tau,1,dhvdz)

    drdz = matmul(rv / r, drvdvv0)
    dlngamdz = alpha * (- drdz / r)

    call outer(rv0, dfvdz, tmp)
    dvvdvv0 = gam * tmp
    call outer(vv0, dgvdz, tmp)
    dvvdvv0 = dvvdvv0 + gam * tmp
    call outer(p, dhvdz, tmp)
    dvvdvv0 = dvvdvv0 + gam * tmp
    call outer(vv, dlngamdz, tmp)
    dvvdvv0 = dvvdvv0 + tmp
    do i = 1,3
        dvvdvv0(i, i) = dvvdvv0(i, i) + gam * gv
    end do

    ! Derivatives wrt z = p

    call evaldp(rv0,vv0,p,nmax,C,1,s0tau,0,dfrdz)
    call evaldp(rv0,vv0,p,nmax,C,2,s0tau,0,dgrdz)
    call evaldp(rv0,vv0,p,nmax,C,3,s0tau,0,dhrdz)
    call evaldp(rv0,vv0,p,nmax,C,4,s0tau,0,dtdz)

    dtdp = dtdz

    call outer(rv0, dfrdz, tmp)
    drvdp = tmp
    call outer(vv0, dgrdz, tmp)
    drvdp = drvdp + tmp
    call outer(p, dhrdz, tmp)
    drvdp = drvdp + tmp
    do i = 1,3
        drvdp(i, i) = drvdp(i, i) + hr
    end do

    call evaldp(rv0,vv0,p,nmax,C,1,s0tau,1,dfvdz)
    call evaldp(rv0,vv0,p,nmax,C,2,s0tau,1,dgvdz)
    call evaldp(rv0,vv0,p,nmax,C,3,s0tau,1,dhvdz)

    drdz = matmul(rv / r, drvdp)
    dlngamdz = alpha * (- drdz / r)

    call outer(rv0, dfvdz, tmp)
    dvvdp = gam * tmp
    call outer(vv0, dgvdz, tmp)
    dvvdp = dvvdp + gam * tmp
    call outer(p, dhvdz, tmp)
    dvvdp = dvvdp + gam * tmp
    call outer(vv, dlngamdz, tmp)
    dvvdp = dvvdp + tmp
    do i = 1,3
        dvvdp(i, i) = dvvdp(i, i) + gam * hv
    end do

    ! Check first integrals - energy, (L, p), (A, p)

    q = 1.0_dp / r
    sigma = dot_product(rv, vv)
    rho = dot_product(rv, p)
    omega = dot_product(vv, p)
    nu = dot_product(vv, vv)

    ham = nu / 2 - q - rho
    Lp = det3(rv, vv, p)
    Ap = rho * (-nu + q + rho / 2) + omega * sigma - r**2 * theta / 2

    integr(1, 1) = ham0
    integr(2, 1) = Lp0
    integr(3, 1) = Ap0
    integr(1, 2) = ham
    integr(2, 2) = Lp
    integr(3, 2) = Ap
end subroutine FandGex
