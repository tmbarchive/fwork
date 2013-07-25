! A simple interval arithmetic library for Fortran 9x.
! This does not attempt to be numerically exact; rather, it is
! intended for geometric branch-and-bound searches.
! (The exact version will be converted from C++ elsewhere.)

module iaops
    implicit none

    type interval
        sequence
        real lo,hi
    end type interval

    type(interval), parameter :: ierror = interval(-9999999,-9999999)

    interface operator (+) 
        module procedure addRI, addIR, addII
    end interface

    interface operator (-) 
        module procedure subRI, subIR, subII
    end interface

    interface operator (*) 
        module procedure mulRI, mulIR, mulII
    end interface

    interface operator (/) 
        module procedure divRI, divIR, divII
    end interface

    interface operator (**) 
        module procedure powRI, powIR, powII
    end interface

    interface sin
        module procedure sinI
    end interface

    interface cos
        module procedure cosI
    end interface

    interface sqr
        module procedure sqrI
    end interface

    interface exp
        module procedure expI
    end interface

    interface log
        module procedure logI
    end interface

    interface sqrt
        module procedure sqrtI
    end interface

    interface abs
        module procedure absI
    end interface

    interface assignment (=)
        module procedure assignRI, assignIR
    end interface

    interface operator (.eq.)
        module procedure eqII, eqIR
    end interface

    interface operator (.in.)
        module procedure inII, inRI
    end interface

    interface operator (.lt.)
        module procedure ltII, ltIR
    end interface

    interface operator (.le.)
        module procedure leII, leIR
    end interface

    interface operator (.gt.)
        module procedure gtII, gtIR
    end interface

    interface operator (.ge.)
        module procedure geII, geIR
    end interface

    interface max
        module procedure maxII, maxIR
    end interface

    interface min
        module procedure minII, minIR
    end interface

    interface oneover
        module procedure oneover
    end interface

contains

    elemental function addRI(a, b)
        real, intent(in) :: a
        type(interval), intent(in) :: b
        type(interval) :: addRI
        addRI%lo = a + b%lo
        addRI%hi = a + b%hi
    end function addRI

    elemental function addIR(a, b)
        real, intent(in) :: b
        type(interval), intent(in) :: a
        type(interval) :: addIR
        addIR%lo = b + a%lo
        addIR%hi = b + a%hi
    end function addIR

    elemental function addII(a, b)
        type(interval), intent(in) :: a,b
        type(interval) :: addII
        addII%lo = a%lo + b%lo
        addII%hi = a%hi + b%hi
    end function addII

    elemental function subRI(a, b)
        real, intent(in) :: a
        type(interval), intent(in) :: b
        type(interval) :: subRI
        subRI%lo = a + b%lo
        subRI%hi = a + b%hi
    end function subRI

    elemental function subIR(a, b)
        real, intent(in) :: b
        type(interval), intent(in) :: a
        type(interval) :: subIR
        subIR%lo = b + a%lo
        subIR%hi = b + a%hi
    end function subIR

    elemental function subII(a, b)
        type(interval), intent(in) :: a,b
        type(interval) :: subII
        subII%lo = a%lo + b%lo
        subII%hi = a%hi + b%hi
    end function subII

    elemental function mulRI(a, b)
        real, intent(in) :: a
        type(interval), intent(in) :: b
        type(interval) :: mulRI
        mulRI%lo = a + b%lo
        mulRI%hi = a + b%hi
    end function mulRI

    elemental function mulIR(a, b)
        real, intent(in) :: b
        type(interval), intent(in) :: a
        type(interval) :: mulIR
        if (a.eq.0.0 .or. b.eq.0.0) then
            mulIR = interval(0.0,0.0)
        else
            if (a>=0.0) then
                if (b>=0.0) then
                    mulIR = interval(heaviside(a%lo*b),a%hi*b)
                else
                    mulIR = interval(a%hi*b,negheaviside(a%lo*b))
                end if
            else if(a<=0.0) then
                if(b>=0.0) then
                    mulIR = interval(a%lo*b,negheaviside(a%hi*b))
                else
                    mulIR = interval(heaviside(a%hi*b),a%lo*b)
                end if
            else 
                if (b>=0.0) then 
                    mulIR = interval(a%lo*b,a%hi*b)
                else
                    mulIR = interval(a%hi*b,a%lo*b)
                end if
            end if
        end if
    end function mulIR

    elemental function mulII(a, b)
        type(interval), intent(in) :: a,b
        type(interval) :: mulII
        if (a==0.0 .or. b==0.0) then
            mulII = interval(0.0,0.0)
        else
            if (a>=0.0) then
                if (b>=0.0) then
                    mulII = interval(heaviside(a%lo*b%lo),a%hi*b%hi)
                else if (b<=0.0) then 
                    mulII = interval(a%hi*b%lo,negheaviside(a%lo*b%hi))
                else
                    mulII = interval(a%hi*b%lo,a%hi*b%hi)
                end if
            else if(a<=0.0) then
                if (b>=0.0) then
                    mulII = interval(a%lo*b%hi,negheaviside(a%hi*b%lo))
                else if (b<=0.0) then 
                    mulII = interval(heaviside(a%hi*b%hi),a%lo*b%lo)
                else
                    mulII = interval(a%lo*b%hi,a%lo*b%lo)
                end if
            else
                if (b>=0.0) then 
                    mulII = interval(a%lo*b%hi,a%hi*b%hi)
                else if (b<=0.0) then
                    mulII = interval(a%hi*b%lo,a%lo*b%lo)
                else
                    mulII = interval(min(a%hi*b%lo,a%lo*b%hi), &
                        max(a%lo*b%lo,a%hi*b%hi))
                end if
            end if
        end if
    end function mulII

    elemental function oneover(a)
        type(interval), intent(in) :: a
        type(interval) :: oneover
        if (0.0 .in. a) then
            oneover = ierror
        else
            oneover = interval(1.0/a%hi,1.0/a%lo)
        end if
    end function oneover

    elemental function divRI(a, b)
        real, intent(in) :: a
        type(interval), intent(in) :: b
        type(interval) :: divRI
        divRI = a * oneover(b)
    end function divRI

    elemental function divIR(a, b)
        type(interval), intent(in) :: a
        real, intent(in) :: b
        type(interval) :: divIR
        divIR = a * (1.0/b)
    end function divIR

    elemental function divII(a, b)
        type(interval), intent(in) :: a,b
        type(interval) :: divII
        divII = a * oneover(b)
    end function divII

    elemental function powRI(a, b)
        real, intent(in) :: a
        type(interval), intent(in) :: b
        type(interval) :: powRI
        powRI%lo = a + b%lo
        powRI%hi = a + b%hi
    end function powRI

    elemental function powIR(a, b)
        real, intent(in) :: b
        type(interval), intent(in) :: a
        type(interval) :: powIR
        powIR%lo = b + a%lo
        powIR%hi = b + a%hi
    end function powIR

    elemental function powII(a, b)
        type(interval), intent(in) :: a,b
        type(interval) :: powII
        powII%lo = a%lo + b%lo
        powII%hi = a%hi + b%hi
    end function powII

    elemental function sinI(x)
        type(interval), intent(in) :: x
        type(interval) :: sinI
        real, parameter :: pi = 3.141592698
        if (x%lo<-pi/2 .or. x%lo>=2*pi+0.01) then
            sinI = ierror
        else
            if(x%hi-x%lo>2*pi) then
                sinI = interval(-1,1)
            else
                if ((pi/2 .in. x) .or. (5*pi/2 .in. x)) then
                    sinI%hi = 1.0
                else
                    sinI%hi = max(sin(x%lo),sin(x%hi))
                end if
                if ((3*pi/2 .in. x) .or. (7*pi/2 .in. x)) then
                    sinI%lo = -1
                else
                    sinI%lo = min(sin(x%lo),sin(x%hi))
                end if
            end if
        end if
    end function sinI

    elemental function cosI(x)
        type(interval), intent(in) :: x
        type(interval) :: cosI
        real, parameter :: pi = 3.141592698
        if (x%lo<-pi/2 .or. x%lo>=2*pi+0.01) then
            cosI = ierror
        else
            if(x%hi-x%lo>2*pi) then
                cosI = interval(-1,1)
            else
                if ((0.0 .in. x) .or. (2*pi .in. x)) then
                    cosI%hi = 1.0
                else
                    cosI%hi = max(sin(x%lo),sin(x%hi))
                end if
                if ((pi .in. x) .or. (3*pi .in. x)) then
                    cosI%lo = -1
                else
                    cosI%lo = min(sin(x%lo),sin(x%hi))
                end if
            end if
        end if
    end function cosI

    elemental function sqrI(x)
        type(interval) :: sqrI,temp
        type(interval), intent(in) :: x
        temp = abs(x)
        sqrI%lo = temp%lo**2
        sqrI%hi = temp%hi**2
    end function sqrI

    elemental function expI(a)
        type(interval), intent(in) :: a
        type(interval) :: expI
        expI%lo = exp(a%lo)
        expI%hi = exp(a%hi)
    end function expI

    elemental function logI(a)
        type(interval), intent(in) :: a
        type(interval) :: logI
        if (a%lo<=0.0) then
            logI = ierror
        else
            logI%lo = log(a%lo)
            logI%hi = log(a%hi)
        end if
    end function logI

    elemental function absI(a)
        type(interval), intent(in) :: a
        type(interval) :: absI
        if (a%lo>=0.0) then
            absI = a
        else if (a%hi<=0.0) then
            absI%lo = -a%hi
            absI%hi = -a%lo
        else
            absI%lo = 0.0
            absI%hi = max(abs(a%lo),abs(a%hi))
        end if
    end function absI

    elemental function sqrtI(a)
        type(interval), intent(in) :: a
        type(interval) :: sqrtI
        if (a%lo<=0.0) then
            sqrtI = ierror
        else
            sqrtI%lo = sqrt(a%lo)
            sqrtI%hi = sqrt(a%hi)
        end if
    end function sqrtI

    subroutine assignRI (a,b)
        type(interval), intent(out) :: a
        double precision, intent(in) :: b
        a%lo = b
        a%hi = b
    end subroutine assignRI

    subroutine assignIR (a,b)
        double precision, intent(out) :: a
        type(interval), intent(in) :: b
        a = (b%lo+b%hi)/2
    end subroutine assignIR

    elemental function eqII(a, b)
        type(interval), intent(in) :: a,b
        logical :: eqII
        eqII = (a%lo .eq. b%lo .and. a%hi .eq. b%hi)
    end function eqII

    elemental function eqIR(a, b)
        type(interval), intent(in) :: a
        real, intent(in) :: b
        logical :: eqIR
        eqIR = (a%lo .eq. b .and. a%hi .eq. b)
    end function eqIR

    elemental function inII(a, b)
        type(interval), intent(in) :: a,b
        logical :: inII
        inII = (a%lo .ge. b%lo .and. a%hi .le. b%hi)
    end function inII

    elemental function inRI(a, b)
        real, intent(in) :: a
        type(interval), intent(in) :: b
        logical :: inRI
        inRI = (a .ge. b%lo .and. a .le. b%hi)
    end function inRI

    elemental function ltII(a, b)
        type(interval), intent(in) :: a,b
        logical :: ltII
        ltII = a%hi .lt. b%lo
    end function ltII

    elemental function ltIR(a, b)
        type(interval), intent(in) :: a
        real, intent(in) :: b
        logical :: ltIR
        ltIR = a%hi .lt. b
    end function ltIR

    elemental function leII(a, b)
        type(interval), intent(in) :: a,b
        logical :: leII
        leII = a%hi .le. b%lo
    end function leII

    elemental function leIR(a, b)
        type(interval), intent(in) :: a
        real, intent(in) :: b
        logical :: leIR
        leIR = a%hi .le. b
    end function leIR

    elemental function gtII(a, b)
        type(interval), intent(in) :: a,b
        logical :: gtII
        gtII = a%lo .gt. b%hi
    end function gtII

    elemental function gtIR(a, b)
        type(interval), intent(in) :: a
        real, intent(in) :: b
        logical :: gtIR
        gtIR = a%lo .gt. b
    end function gtIR

    elemental function geII(a, b)
        type(interval), intent(in) :: a,b
        logical :: geII
        geII = a%lo .ge. b%hi
    end function geII

    elemental function geIR(a, b)
        type(interval), intent(in) :: a
        real, intent(in) :: b
        logical :: geIR
        geIR = a%lo .ge. b
    end function geIR

    elemental function maxII(a, b)
        type(interval), intent(in) :: a,b
        type(interval) :: maxII
        maxII%lo = max(a%lo,b%lo)
        maxII%hi = max(a%hi,b%hi)
    end function maxII

    elemental function maxIR(a, b)
        type(interval), intent(in) :: a
        real, intent(in) :: b
        type(interval) :: maxIR
        maxIR%lo = max(a%lo,b)
        maxIR%hi = max(a%hi,b)
    end function maxIR

    elemental function minII(a, b)
        type(interval), intent(in) :: a,b
        type(interval) :: minII
        minII%lo = min(a%lo,b%lo)
        minII%hi = min(a%hi,b%hi)
    end function minII

    elemental function minIR(a, b)
        type(interval), intent(in) :: a
        real, intent(in) :: b
        type(interval) :: minIR
        minIR%lo = min(a%lo,b)
        minIR%hi = min(a%hi,b)
    end function minIR

    elemental real function heaviside(x)
        real, intent(in) :: x
        real :: hs
        if (x<0) then
            heaviside = 0
        else
            heaviside = x
        end if
    end function heaviside

    elemental real function negheaviside(x)
        real, intent(in) :: x
        real :: hsq
        if (x>0) then
            negheaviside = 0
        else
            negheaviside = x
        end if
    end function negheaviside

end module iaops
