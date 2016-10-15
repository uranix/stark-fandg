subroutine FG(nmax,q,s,r,w,v,p,C)
    integer, intent(in) :: nmax
    integer, parameter :: dp = kind(1.d0)
    real(dp), intent(in) :: q,s,r,w,v,p
    real(dp), intent(out) :: C(0:nmax, 4)
    select case (nmax)
        case (5)
            call FG5(q,s,r,w,v,p,C)
        case (6)
            call FG6(q,s,r,w,v,p,C)
        case (7)
            call FG7(q,s,r,w,v,p,C)
        case (8)
            call FG8(q,s,r,w,v,p,C)
        case (9)
            call FG9(q,s,r,w,v,p,C)
        case (10)
            call FG10(q,s,r,w,v,p,C)
        case (11)
            call FG11(q,s,r,w,v,p,C)
        case (12)
            call FG12(q,s,r,w,v,p,C)
        case (13)
            call FG13(q,s,r,w,v,p,C)
        case (14)
            call FG14(q,s,r,w,v,p,C)
        case (15)
            call FG15(q,s,r,w,v,p,C)
!        case (16)
!            call FG16(q,s,r,w,v,p,C)
!        case (17)
!            call FG17(q,s,r,w,v,p,C)
!        case (18)
!            call FG18(q,s,r,w,v,p,C)
!        case (19)
!            call FG19(q,s,r,w,v,p,C)
!        case (20)
!            call FG20(q,s,r,w,v,p,C)
        case default
            stop 'Invalid nmax for FG'
    end select
end subroutine FG

subroutine FGex(nmax,q,s,r,w,v,p,C)
    integer, intent(in) :: nmax
    integer, parameter :: dp = kind(1.d0)
    real(dp), intent(in) :: q,s,r,w,v,p
    real(dp), intent(out) :: C(0:nmax, 28)
    select case (nmax)
        case (5)
            call FG5ex(q,s,r,w,v,p,C)
        case (6)
            call FG6ex(q,s,r,w,v,p,C)
        case (7)
            call FG7ex(q,s,r,w,v,p,C)
        case (8)
            call FG8ex(q,s,r,w,v,p,C)
        case (9)
            call FG9ex(q,s,r,w,v,p,C)
        case (10)
            call FG10ex(q,s,r,w,v,p,C)
        case (11)
            call FG11ex(q,s,r,w,v,p,C)
        case (12)
            call FG12ex(q,s,r,w,v,p,C)
        case (13)
            call FG13ex(q,s,r,w,v,p,C)
        case (14)
            call FG14ex(q,s,r,w,v,p,C)
        case (15)
            call FG15ex(q,s,r,w,v,p,C)
!        case (16)
!            call FG16ex(q,s,r,w,v,p,C)
!        case (17)
!            call FG17ex(q,s,r,w,v,p,C)
!        case (18)
!            call FG18ex(q,s,r,w,v,p,C)
!        case (19)
!            call FG19ex(q,s,r,w,v,p,C)
!        case (20)
!            call FG20ex(q,s,r,w,v,p,C)
        case default
            stop 'Invalid nmax for FG'
    end select
end subroutine FGex

