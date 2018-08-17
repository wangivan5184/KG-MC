	write(*,*)
    write(*,*) 'SCMC Simulation ', its
	write(*,12) 'vb0:', (vb0(i),i=1,m)
    write(22,*)
    write(22,*)
    write(22,*) 'SCMC Simulation ', its
	write(22,12) '                   vb0:', (vb0(i),i=1,m)


	close(22)	! Make the resfn accessible for the subroutine MC below
	
	call MC

    ! Output the simulation time, ensemble averages and errors
	open (22, file=resfn, access='append')
    write(22,*)
    write(22,2) runtime
    write(22,*)
    write(22,5) 
    write(22,6) 
    do i = 1, mA
		Asgm(i) = 3 * sqrt((A2bar(i)-Abar(i)**2)*(2*tau(i)+1)/Nsmp)
		dRe2(i) = Abar(i) - (i-1.0d0/3)/N
        write(22,3) i, Nc(i), 2*tau(i)+1, Abar(i), Asgm(i), dRe2(i), maxval(abs(Ps(:,i)-Pt(:,i)))
    end do
	write(*,12) '|dRe2|:', (abs(dRe2(i)), i=1,N-1)
	
	chisq = 0
	do i = 1, N-1
		chisq = chisq + (dRe2(i)/Asgm(i))**2
	end do
	write(22,4) chisq
	write(*,4) chisq
    close(22)

	good = .true.
    do i = 1, N-1
		if (abs(dRe2(i)) > Asgm(i))  then
			good = .false.
			exit
		end if
    end do
