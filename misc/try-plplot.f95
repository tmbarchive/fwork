program testpl
    use plplot
    real(plflt) :: xs(100),ys(100)
    call plinit
    call plenv(0.0_plflt,100.0_plflt,0.0_plflt,100.0_plflt,0,0)
    do i=1,100
        xs(i) = i
        ys(i) = 10.0*sqrt(i*1.0)
    end do
    call plcol0(4)
    call plpoin(xs,ys,100)
    call plline(xs,ys)
    call plend
end program testpl
