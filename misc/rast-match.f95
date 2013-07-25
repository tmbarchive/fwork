module rast_match
    use rast
    use jpegio
    use dataio
    use cedges
    use plutil
    use surfmod
    implicit none
contains

    subroutine read_points(points,fname)
        character(*) fname
        type(point), allocatable :: points(:)
        real, allocatable :: a(:,:)
        call read_table(a,fname)
        print *,"table",size(a,1),size(a,2)
        allocate(points(size(a,1)))
        points(:)%x = a(:,1)
        points(:)%y = a(:,2)
        points(:)%a = a(:,3)
    end subroutine read_points

    subroutine simple_problem(problem)
        character(len=1000) fname
        type(pp_problem) problem
        type(surfdesc), allocatable :: spoints(:),sipoints(:)
        integer nmodel,nimage,i,j,index
        fname = "points1.pts"
        call read_points(problem%model,fname)
        fname = "points2.pts"
        call read_points(problem%image,fname)
        nmodel = size(problem%model)
        nimage = size(problem%image)
        allocate(problem%pairs(nmodel*nimage,2))
        allocate(problem%weights(nmodel*nimage))
        do i=1,nmodel
            do j=1,nimage
                index = (i-1)*nimage+j
                problem%pairs(index,:) = [i,j]
                problem%weights(index) = 1
            end do
        end do
        print *,"model",nmodel
        print *,"image",nimage
    end subroutine simple_problem
    
    subroutine surf_problem(problem)
        character(len=1000) fname
        type(pp_problem) problem
        type(surfdesc), allocatable :: spoints(:),sipoints(:)
        call getarg(1,fname)
        fname = "lenna2.surf"
        call surf_read(spoints,fname)
        call getarg(2,fname)
        fname = "lenna.surf"
        call surf_read(sipoints,fname)
        call sconvert(problem%model,spoints)
        call sconvert(problem%image,sipoints)
        call initial_matchlist(problem%pairs,problem%weights,spoints,sipoints)
        problem%weights = 1.0
        print *,"model",size(spoints,1)
        print *,"image",size(sipoints,1)
    end subroutine surf_problem

    subroutine sconvert(points,spoints)
        type(point), allocatable :: points(:)
        type(surfdesc) spoints(:)
        integer i
        allocate(points(size(spoints)))
        do i=1,size(spoints)
            points(i)%x = spoints(i)%x
            points(i)%y = spoints(i)%y
            points(i)%a = 0
        end do
    end subroutine sconvert

    real function dist(u,v)
        real u(:),v(:)
        dist = sqrt(sum((u-v)**2))
    end function dist

    subroutine normalize(a)
        real a(:)
        a = (a-minval(a))/(maxval(a)-minval(a))
    end subroutine normalize

    subroutine initial_matchlist(pairs,weights,points,ipoints)
        integer(2), allocatable :: pairs(:,:)
        real, allocatable :: weights(:)
        type(surfdesc) :: points(:),ipoints(:)
        real mscore,thresh,d
        integer i,j,npairs,nmax,index

        nmax = size(points,1)*size(ipoints,1)
        allocate(weights(nmax))
        allocate(pairs(nmax,2))
        do i=1,size(points,1)
            do j=1,size(ipoints,1)
                index = (i-1)*size(ipoints,1)+j
                pairs(index,1) = i
                pairs(index,2) = j
                weights(index)= dist(points(i)%descriptor,ipoints(j)%descriptor)
            end do
        end do
        weights = max(0.0,weights/maxval(weights))
    end subroutine initial_matchlist

    subroutine plot_pp(problem,list)
        type(pp_problem) problem
        integer, allocatable :: list(:),l(:)
        real tr(4,2),score(2),image(500,500)
        integer i,j

        print *,"list",size(list)
        allocate(l(size(list)))
        do i=1,size(image,1)
            do j=1,size(image,2)
                tr(:,1) = [i,j,0,0]
                tr(:,2) = [i+1,j+1,0,0]
                l = list
                score = evaluate_pp(problem,l,tr)
                image(i,j) = score(2)
            end do
        end do
        print *,"range",minval(image),maxval(image)
        image = image/maxval(image)
        call write_jpeg2("map.jpg",image)
    end subroutine plot_pp

end module rast_match

program rast_match_main
    use rast
    use jpegio
    use cedges
    use plutil
    use surfmod
    use rast_match

    type(pp_problem) problem
    integer, allocatable :: list(:)
    integer i
    real tr(4,2)
    tr(:,1) = [0.0,0.0,0.0,1.0]
    tr(:,2) = [1000.0,1000.0,0.0,1.0]

    ! call surf_problem(problem)
    call simple_problem(problem)

    allocate(list(size(problem%weights)))
    forall (i=1:size(list)) list(i) = i
    problem%eps = 20.0
    problem%aeps = 0.1
    ! call search_pp(problem,list,tr,4.0)
    call plot_pp(problem,list)
end program rast_match_main
