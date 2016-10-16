subroutine FG(alpha,nmax,q,s,r,w,v,p,C)
    integer, intent(in) :: nmax
    integer, parameter :: dp = kind(1.d0)
    real(dp), intent(in) :: alpha,q,s,r,w,v,p
    real(dp), intent(out) :: C(0:nmax, 4)
    if (alpha == 0._dp) then
        call fg_a0(nmax,q,s,r,w,v,p,C)
    else if (alpha == 1._dp) then
        call fg_a1(nmax,q,s,r,w,v,p,C)
    else if (alpha == 1.5_dp) then
        call fg_a3_2(nmax,q,s,r,w,v,p,C)
    else if (alpha == 2._dp) then
        call fg_a2(nmax,q,s,r,w,v,p,C)
    else
        stop "Bad alpha for FG"
    end if
end subroutine FG

subroutine FGex(alpha,nmax,q,s,r,w,v,p,C)
    integer, intent(in) :: nmax
    integer, parameter :: dp = kind(1.d0)
    real(dp), intent(in) :: alpha,q,s,r,w,v,p
    real(dp), intent(out) :: C(0:nmax, 4)
    if (alpha == 0._dp) then
        call fg_a0ex(nmax,q,s,r,w,v,p,C)
    else if (alpha == 1._dp) then
        call fg_a1ex(nmax,q,s,r,w,v,p,C)
    else if (alpha == 1.5_dp) then
        call fg_a3_2ex(nmax,q,s,r,w,v,p,C)
    else if (alpha == 2._dp) then
        call fg_a2ex(nmax,q,s,r,w,v,p,C)
    else
        stop "Bad alpha for FGex"
    end if
end subroutine FGex
