! Calculate the total bonded interaction energy HCt.

    HCt = 0
    do k = 2, N
        Rxk = Rx(k)
        Ryk = Ry(k)
        Rzk = Rz(k)

        ! Calculate the energy of the segment
        jb = 1
        je = k-1
        include "ESeg.h"
        HCt = HCt + HCs
    end do
