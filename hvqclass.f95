module hvqclass

  implicit none
  
contains

  integer function randint(n)
    integer n
    real v
    call random_number(v)
    randint = modulo(floor(n*v),n)+1
  end function randint

  pure real function dist(u,v,n)
    integer, intent(in) :: n
    real, intent(in) :: u(n),v(n)
    dist = sqrt(sum((u-v)**2))
  end function dist

  subroutine argmindist2(i1,d1,i2,d2,v,means,k,d)
    integer, intent(out) :: i1,i2
    real, intent(out) :: d1,d2
    integer, intent(in) :: k,d
    real, intent(in) :: v(d),means(k,d)
    integer i
    real dists(k)
    ! FIXME add openmp comments here for parallelism 
    do i=1,k
       dists(i) = dist(v,means(i,:),d)
       if (dists(i)<0 .or. dists(i)>1e10) call abort()
    end do
    i1 = minloc(dists,1)
    d1 = dists(i1)
    dists(i1) = 1e30
    i2 = minloc(dists,1)
    d2 = dists(i2)
  end subroutine argmindist2

  integer function argmindist(v,means,k,d)
    integer, intent(in) :: k,d
    real, intent(in) :: v(d),means(k,d)
    integer i1,i2
    real d1,d2
    call argmindist2(i1,d1,i2,d2,v,means,k,d)
    argmindist = i1
  end function argmindist

  ! standard k-means algorithm; too slow for large
  ! datasets
  
  subroutine kmeans(means,data,maxiter,k,n,d)
    integer, intent(in) :: maxiter,k,n,d
    real, intent(out) :: means(k,d)
    real, intent(in) ::  data(n,d)
    integer assignment(n),last(n),counts(k)
    integer i,c,iter
    print *,"kmeans",k,n,d
    do i=1,k
       means(i,:) = data(i*(n/k),:)
    end do
    last(:) = 0
    do iter=1,maxiter
       print *,"iter",iter
       do i=1,n
          assignment(i) = argmindist(data(i,:),means,k,d)
       end do
       counts(:) = 0
       means(:,:) = 0
       do i=1,n
          c = assignment(i)
          counts(c) = counts(c) + 1
          means(c,:) = means(c,:) + data(i,:)
       end do
       if (all(assignment(:)==last(:))) exit
       do i=1,k
          if (counts(i)/=0) then
             means(i,:) = means(i,:) / real(counts(i))
          else
             means(i,:) = data(i,:) ! bad initializer
          end if
       end do
       last(:) = assignment(:)
    end do
    print *,"done"
  end subroutine kmeans

  ! incremental kmeans algorithm
  ! TODO: speed up with bounds

  subroutine incremental_kmeans(means,data,maxiter,rate,k,n,d)
    integer, intent(in) :: maxiter,k,n,d
    real, intent(out) :: means(k,d)
    real, intent(in) ::  data(n,d), rate
    integer c,i2,j,iter
    real d1,d2,dist,l
    do c=1,k
       j = randint(n)
       means(c,:) = data(j,:)
    end do
    dist = 0.0
    l = 1.0/1024.0
    do iter=1,maxiter
       if (modulo(iter,5000)==0) print *,"ikm",iter,dist
       j = randint(n)
       call argmindist2(c,d1,i2,d2,data(j,:),means,k,d)
       means(c,:) = means(c,:) + rate * (data(j,:) - means(c,:))
       dist = (1.0-l) * dist + l * d1
    end do
  end subroutine incremental_kmeans

  ! neural gas algorithm; basically the same as kmeans, but
  ! has a little bit of "regularization" in it (may lead
  ! to better balanced sets, fewer outliers)

  ! TODO: speed up with bounds

  subroutine neural_gas(means,data,maxiter,rate1,rate2,k,n,d)
    integer, intent(in) :: maxiter,k,n,d
    real, intent(out) :: means(k,d)
    real, intent(in) ::  data(n,d), rate1, rate2
    integer i1,i2,iter,i,j
    real d1,d2,dist,l
    do i=1,k
       i1 = randint(n)
       means(i,:) = data(i1,:)
    end do
    dist = 0.0
    l = 1.0/1024.0
    do iter=1,maxiter
       if (modulo(iter,10000)==0) print *,"ngq",iter,dist
       j = randint(n)
       call argmindist2(i1,d1,i2,d2,data(j,:),means,k,d)
       means(i1,:) = means(i1,:) + rate1 * (data(j,:) - means(i1,:))
       means(i2,:) = means(i2,:) + rate2 * (data(j,:) - means(i2,:))
       dist = (1.0-l) * dist + l * d1
    end do
  end subroutine neural_gas

end module
