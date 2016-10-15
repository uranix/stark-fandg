program test
    implicit none
    integer, parameter :: dp = kind(1.d0)
    integer, parameter :: kmax = 15
    real(dp) :: rv0(3),vv0(3),p(3),tau
    real(dp) :: rv(3), vv(3), t
    real(dp) :: integr(3, 2)
    real(dp) :: rv1(3), vv1(3), t1, rv2(3), vv2(3), t2
    real(dp), parameter :: delta = 1e-6_dp
    real(dp), parameter :: dvec(3) = (/-0.470387_dp, -0.9156114_dp, 0.71380349_dp/)
    real(dp) :: drvdtau(3), dvvdtau(3), dtdtau
    real(dp) :: drvdrv0(3, 3), dvvdrv0(3, 3), dtdrv0(3)
    real(dp) :: drvdvv0(3, 3), dvvdvv0(3, 3), dtdvv0(3)
    real(dp) :: drvdp(3, 3), dvvdp(3, 3), dtdp(3)
    rv0 = (/ 1._dp, 2._dp, 3._dp/)
    vv0 = (/-1._dp, 1._dp, 2._dp/)
    p =   (/ 1._dp, 1._dp, 1._dp/)

    write(*,*) 'rv0 = ', rv0
    write(*,*) 'vv0 = ', vv0
    write(*,*) 'p = ', p

    tau = 0.1_dp
    call FandGex(rv0, vv0, p, tau, kmax, rv, vv, t, &
        drvdtau, dvvdtau, dtdtau, &
        drvdrv0, dvvdrv0, dtdrv0, &
        drvdvv0, dvvdvv0, dtdvv0, &
        drvdp, dvvdp, dtdp, integr)
    write(*,*) 'r = ', rv
    write(*,*) 'v = ', vv
    write(*,*) 't = ', t
    write(*,*) '|H0 - H| / |H0| = ', dabs(integr(1, 1) - integr(1,2)) / dabs(integr(1, 1))
    write(*,*) '============================== d/dtau ==============================='
    write(*,*) 'drvdtau =', drvdtau
    write(*,*) 'dvvdtau =', dvvdtau
    write(*,*) 'dtdtau  =', dtdtau
    write(*,*) '============================== d/drv0 ==============================='
    write(*,*) 'dtdrv0  =', dtdrv0
    write(*,*) '         ', drvdrv0(1, :)
    write(*,*) 'drvdrv0 =', drvdrv0(2, :)
    write(*,*) '         ', drvdrv0(3, :)
    write(*,*) '         ', dvvdrv0(1, :)
    write(*,*) 'dvvdrv0 =', dvvdrv0(2, :)
    write(*,*) '         ', dvvdrv0(3, :)
    write(*,*) '============================== d/dvv0 ==============================='
    write(*,*) 'dtdvv0  =', dtdvv0
    write(*,*) '         ', drvdvv0(1, :)
    write(*,*) 'drvdvv0 =', drvdvv0(2, :)
    write(*,*) '         ', drvdvv0(3, :)
    write(*,*) '         ', dvvdvv0(1, :)
    write(*,*) 'dvvdvv0 =', dvvdvv0(2, :)
    write(*,*) '         ', dvvdvv0(3, :)
    write(*,*) '=============================== d/dp ================================'
    write(*,*) 'dtdp    =', dtdp
    write(*,*) '         ', drvdp(1, :)
    write(*,*) 'drvdp   =', drvdp(2, :)
    write(*,*) '         ', drvdp(3, :)
    write(*,*) '         ', dvvdp(1, :)
    write(*,*) 'dvvdp   =', dvvdp(2, :)
    write(*,*) '         ', dvvdp(3, :)

    write(*,*) '====================== Numerical verification ======================='
    call FandG(rv0, vv0, p, tau * (1 - delta), kmax, rv1, vv1, t1, integr)
    call FandG(rv0, vv0, p, tau * (1 + delta), kmax, rv2, vv2, t2, integr)
    write(*,*) 'drvdtau ~ ', ((rv2-rv1) / (2*delta) - tau * drvdtau) / drvdtau
    write(*,*) 'dvvdtau ~ ', ((vv2-vv1) / (2*delta) - tau * dvvdtau) / dvvdtau
    write(*,*) 'dtdtau  ~ ', ((t2-t1) / (2*delta) - tau * dtdtau) / dtdtau

    call FandG(rv0 - delta * dvec, vv0, p, tau, kmax, rv1, vv1, t1, integr)
    call FandG(rv0 + delta * dvec, vv0, p, tau, kmax, rv2, vv2, t2, integr)
    write(*,*) 'drvdrv0 ~ ', ((rv2-rv1) / (2*delta) - matmul(drvdrv0, dvec)) / matmul(drvdrv0, dvec)
    write(*,*) 'dvvdrv0 ~ ', ((vv2-vv1) / (2*delta) - matmul(dvvdrv0, dvec)) / matmul(dvvdrv0, dvec)
    write(*,*) 'dtdrv0  ~ ', ((t2-t1) / (2*delta) - dot_product(dtdrv0, dvec)) / dot_product(dtdrv0, dvec)

    call FandG(rv0, vv0 - delta * dvec, p, tau, kmax, rv1, vv1, t1, integr)
    call FandG(rv0, vv0 + delta * dvec, p, tau, kmax, rv2, vv2, t2, integr)
    write(*,*) 'drvdvv0 ~ ', ((rv2-rv1) / (2*delta) - matmul(drvdvv0, dvec)) / matmul(drvdvv0, dvec)
    write(*,*) 'dvvdvv0 ~ ', ((vv2-vv1) / (2*delta) - matmul(dvvdvv0, dvec)) / matmul(dvvdvv0, dvec)
    write(*,*) 'dtdvv0  ~ ', ((t2-t1) / (2*delta) - dot_product(dtdvv0, dvec)) / dot_product(dtdvv0, dvec)

    call FandG(rv0, vv0, p - delta * dvec, tau, kmax, rv1, vv1, t1, integr)
    call FandG(rv0, vv0, p + delta * dvec, tau, kmax, rv2, vv2, t2, integr)
    write(*,*) 'drvdp   ~ ', ((rv2-rv1) / (2*delta) - matmul(drvdp, dvec)) / matmul(drvdp, dvec)
    write(*,*) 'dvvdp   ~ ', ((vv2-vv1) / (2*delta) - matmul(dvvdp, dvec)) / matmul(dvvdp, dvec)
    write(*,*) 'dtdp    ~ ', ((t2-t1) / (2*delta) - dot_product(dtdp, dvec)) / dot_product(dtdp, dvec)
end
