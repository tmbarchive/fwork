module test_rast_mod
    use rast
    use jpegio
contains
    subroutine checkeq(a,b)
        real a,b
        if (max(a,b)>1e-6) then
            if (abs(a-b)/abs(max(a,b))>1e-6) then
                print *,"checkeq failed",a,b
                call abort()
            end if
        end if
    end subroutine checkeq

    subroutine checkle(a,b)
        real a,b
        if (a>b+abs(b*1e-6)) then
            print *,"checkle failed",a,b
            call abort()
        end if
    end subroutine checkle

    subroutine checkge(a,b)
        real a,b
        if (b>a+abs(a*1e-6)) then
            print *,"checkge failed",a,b
            call abort()
        end if
    end subroutine checkge

    subroutine test_segment_point_dist
        real c,s,a
        real x0,y0,x1,y1,x,y
        real d
        do i=0,5
            a = i*0.7
            c = cos(a)
            s = sin(a)
            x0 = c*1.0-s*1.0
            y0 = s*1.0+c*1.0
            x1 = c*4.0-s*1.0
            y1 = s*4.0+c*1.0
            x = c*2.0-s*2.0
            y = c*2.0+s*2.0
            d = segment_point_dist(x0,y0,x1,y1,x,y)
            call checkeq(d,1.0)
        end do
    end subroutine test_segment_point_dist

    function rand()
        real rand
        call random_number(rand)
    end function rand

    subroutine test_a_score
        real da,ra,aeps,score(2)
        integer i
        do i=1,1000
            da = rand()*10.0
            ra = rand()*0.1
            aeps = rand()*1.2
            score = a_score(da,ra,aeps)
            call checkle(score(1),score(2))
            score = a_score(-da,ra,aeps)
            call checkle(score(1),score(2))
        end do
        do i=1,1000
            da = rand()*0.6+0.9
            ra = rand()*0.1
            aeps = rand()*0.8
            score = a_score(da,ra,aeps)
            call checkle(score(2),1e-6)
            score = a_score(-da,ra,aeps)
            call checkle(score(2),1e-6)
        end do
        do i=1,1000
            da = rand()*2.5
            ra = rand()*0.1
            aeps = rand()*4.0+3.0
            score = a_score(da,ra,aeps)
            call checkge(score(1),1e-6)
            score = a_score(-da,ra,aeps)
            call checkge(score(1),1e-6)
        end do
    end subroutine test_a_score

    subroutine plot_test
        real image(512,512)
        type(segment) segments(3)
        type(point) points(1)
        integer(2) :: pairs(3,3)
        real :: tr(4,2)
        real score(2)
        real weights(3)

        points(1) = point(0.0,0.0,0.0)
        segments(1) = segment(70.0,50.0,310.0,50.0,0.0)
        segments(2) = segment(120.0,130.0,90.0,200.0,0.0)
        segments(3) = segment(400.0,410.0,310.0,300.0,0.0)
        pairs(:,1) = [1,1,1]
        pairs(:,2) = [1,2,3]
        weights = 1.0
        image = 0
        del = 0.0
        adel = 0.05
        sdel = 0.05
        eps = 5.0
        do i=1,size(image,1)
            do j=1,size(image,2)
                tr(:,1) = [i-del,j-del,pi/4-adel,1.0-sdel]
                tr(:,2) = [i+del,j+del,pi/4+adel,1.0+sdel]
                score = evaluate_ps(points,segments,pairs,weights,tr,eps,100.0)
                image(i,j) = score(2)
            end do
        end do
        print *,"writing bounds1.jpg",minval(image),maxval(image)
        image = 255 * image / maxval(image)
        call write_jpeg2("bounds1.jpg",image)
    end subroutine plot_test
        
end module test_rast_mod

program test_rast
    use test_rast_mod
    call test_segment_point_dist
    call test_a_score
    call plot_test
end program test_rast
