	do imove = 1, Nmove
		rn = ran()
			
!****************************************
		if (rn < fhop)  then	! Hopping
!****************************************

			Shop = Shop + 1
			k = (N-1)*rn/fhop + 2				! Hopping segment k is in [2,N]

            !**************** Generate the hopping trial move ***************
            Rxk = Rx(k)                     	! Old position
            Ryk = Ry(k)
            Rzk = Rz(k)
            jb = 1
            je = k-1
            include 'ESeg.h'
            HCo = HCs
            jb = k+1
            je = N
            include 'ESeg.h'
            HCo = HCo + HCs

            Rxk = Rxk + (2*ran()-1) * dmax  	! New position
            Ryk = Ryk + (2*ran()-1) * dmax
            Rzk = Rzk + (2*ran()-1) * dmax
            jb = 1
            je = k-1
            include 'ESeg.h'
            HCn = HCs
            jb = k+1
            je = N
            include 'ESeg.h'
            HCn = HCn + HCs
            dE = HCn-HCo

            !**************** Decide on the acceptance ****************
            good = .true.
            if (dE > 0)  then
                Pacc = dexp(-dE)
                if (ran() > Pacc)  good = .false.
            end if

            if (good)  then                 	! Accept the hopping trial move
                Ahop = Ahop + 1
                Rx(k) = Rxk
                Ry(k) = Ryk
                Rz(k) = Rzk
                HC = HC + HCn-HCo
            end if

!**************************************
		else					! Pivot
!**************************************

			kp = (N-2)*(rn-fhop)/(1-fhop) + 2	! Pivot center segment k is in [2, N-1]
			if ((N/2*2 /= N) .and. (kp == N/2+1))  then
				brch = 1						! Rotate the right branch (i.e., segments kp+1 to N)
				dk = 1
			else
				brch = 2*kp/(N+1)				! Rotate the shorter branch (0 for left and 1 for right)
				dk = 2*brch - 1					! dk=-1 for left and 1 for right
			end if
			
			! Generate the 4 random numbers satisfying q0^2 + q1^2 + q2^2 + q3^2 = 1 using the method of Vesely,
			! and then the 3x3 rotation matrix T.
			p1 = 1
			do while (p1 >= 1)
				q0 = 2*ran() - 1
				q1 = 2*ran() - 1
				q02 = q0*q0
				q12 = q1*q1
				p1 = q02 + q12
			end do
			p2 = 1
			do while (p2 >= 1)
				q2 = 2*ran() - 1
				q3 = 2*ran() - 1
				p2 = q2*q2 + q3*q3
			end do
			q2 = q2*dsqrt((1-p1)/p2)
			q3 = q3*dsqrt((1-p1)/p2)
			q22 = q2*q2
			q32 = q3*q3
			t1 = q02 + q12 - q22 - q32
			t2 = 2 * (q1*q2 + q0*q3)
			t3 = 2 * (q1*q3 - q0*q2)
			t4 = 2 * (q1*q2 - q0*q3)
			t5 = q02 - q12 + q22 - q32
			t6 = 2 * (q2*q3 + q0*q1)
			t7 = 2 * (q1*q3 + q0*q2)
			t8 = 2 * (q2*q3 - q0*q1)
			t9 = q02 - q12 - q22 + q32
            
			HCo = 0
			HCn = 0
			jb = 1 + kp*(1-brch)				! jb=kp+1 for left and 1 for right
			je = N - (N-kp+1)*brch				! je=N for left and kp-1 for right; [jb,je] are unrotated segments (excluding the pivot center)
			do k = kp+dk, (N-1)*brch+1, dk		! do k=kp-1,1,-1 for left; and do k=kp+1,N,1 for right
				! Calculate the bonded energy of each rotating segment at its old position
				Rxk = Rx(k)
				Ryk = Ry(k)
				Rzk = Rz(k)
				include 'ESeg.h'
				HCo = HCo + HCs
				
				! Calculate the bonded energy of each rotating segment at its new position
				dRx = Rxk - Rx(kp)
				dRy = Ryk - Ry(kp)
				dRz = Rzk - Rz(kp)
				Rxk = dRx*t1 + dRy*t2 + dRz*t3 + Rx(kp)
				Ryk = dRx*t4 + dRy*t5 + dRz*t6 + Ry(kp)
				Rzk = dRx*t7 + dRy*t8 + dRz*t9 + Rz(kp)
				include 'ESeg.h'
				HCn = HCn + HCs

				! Store the new position in an array
				Rxp(k) = Rxk
				Ryp(k) = Ryk
				Rzp(k) = Rzk
			end do
            dE = HCn-HCo

            !**************** Decide on the acceptance ****************
            good = .true.
            if (dE > 0)  then
                Pacc = dexp(-dE)
                if (ran() > Pacc)  good = .false.
            end if

            if (good)  then                 	! Accept the pivot trial move
                Apvt = Apvt + 1
				HC = HC + dE
				do k = kp+dk, (N-1)*brch+1, dk	! do k=kp-1,1,-1 for left; and do k=kp+1,N,1 for right
					Rx(k) = Rxp(k)
					Ry(k) = Ryp(k)
					Rz(k) = Rzp(k)
				end do
				if (brch == 0)  then			! Shift the chain such that its first segment is at the origin
					Rxk = Rx(1)
					Ryk = Ry(1)
					Rzk = Rz(1)
					Rx = Rx - Rxk
					Ry = Ry - Ryk
					Rz = Rz - Rzk
				end if
            end if

		end if

	end do
