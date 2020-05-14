Program meshproducer
    implicit none

    integer :: i, j, k, loop1, loop2, loop3

    double precision :: startx, endx, starty, endy, startz, endz ! start and end of the fine domain
    double precision :: totalx, totaly, totalz ! size of the fine domain
    double precision :: dx, dy, dz ! size of the fine mesh
    double precision :: dd, ratio ! to calculate the size of the coarse mesh
    integer :: xnodect, ynodect, znodect ! size of the fine domain, later being overwritten as the size of the whole domain
    integer :: morexct, moreyct, morezct ! size of the coarse domain

    double precision, allocatable :: xnodecoord(:), ynodecoord(:), znodecoord(:) ! xyz nodes coordinate array
    integer :: allnodect, allcellct, allbcct, itype, ic, temp1, temp2 ! integers used to assign index to nodes and cells
    integer, allocatable :: nodescan(:,:,:) ! each node is assigned with an integer
    integer, allocatable :: cellscan(:,:,:) ! each cell is assigned with an integer

    integer :: NBSETS
    integer :: BCX1, BCX2, BCY1, BCY2, BCZ1, BCZ2

    real :: temp

    character(len = 40) :: fn

    !define geometry
    startx = -1.
    endx = 1.
    starty = -0.2
    endy = 0.5
    startz = -1.
    endz = 1.

    !define mesh size
    dx = 0.025
    dy = 0.025
    dz = 0.025

    !define boundary condition
    NBSETS = 1
    BCX1 = 1
    BCX2 = 1
    BCY1 = 1
    BCY2 = 1
    BCZ1 = 1
    BCZ2 = 1
    !these number >0 and <= NBSETS

    fn = "largemesh.neu"

    !calculate nodes number
    totalx = endx - startx
    totaly = endy - starty
    totalz = endz - startz

    xnodect = totalx/dx + 1.1
    ynodect = totaly/dy + 1.1
    znodect = totalz/dz + 1.1

    morexct = 10
    moreyct = 10
    morezct = 10

    ratio = 1.2

    print *, 'nodect', xnodect, ynodect, znodect

    allocate(xnodecoord(xnodect + 2 * morexct), ynodecoord(ynodect + moreyct), znodecoord(znodect + 2 * morezct))

    xnodecoord(morexct + 1) = startx
    do i = 2, xnodect
        xnodecoord(morexct + i) = startx + (i - 1) * dx
    end do

    dd = dx
    do i = 1, morexct
        dd = dd * ratio
        xnodecoord(morexct + 1 - i) = xnodecoord(morexct + 2 - i) - dd
        xnodecoord(morexct + xnodect + i) = xnodecoord(morexct + xnodect + i - 1) + dd
    end do

    do i = 1, size(xnodecoord)
        print *, 'xnode', i, xnodecoord(i)
    end do

    ynodecoord(moreyct + 1) = starty
    do i = 2, ynodect
        ynodecoord(moreyct + i) = starty + (i - 1) * dy
    end do

    dd = dy
    do i = 1, moreyct
        dd = dd * ratio
        ynodecoord(moreyct + 1 - i) = ynodecoord(moreyct + 2 - i) - dd
    end do

    do i = 1, size(ynodecoord)
        print *, 'ynode', i, ynodecoord(i)
    end do

    znodecoord(morezct + 1) = startz
    do i = 2, znodect
        znodecoord(morezct + i) = startz + (i - 1) * dz
    end do

    dd = dz
    do i = 1, morezct
        dd = dd * ratio
        znodecoord(morezct + 1 - i) = znodecoord(morezct + 2 - i) - dd
        znodecoord(morezct + znodect + i) = znodecoord(morezct + znodect - 1 + i) + dd
        print *, 'dd', i, dd
        print *, 'low z', morezct + 1 - i, znodecoord(morezct + 1 - i), znodecoord(morezct + 2 - i)
        print *, 'highz', morezct + znodect + i, znodecoord(morezct + znodect + i), znodecoord(morezct + znodect - 1 + i)
    end do

    do i = 1, size(znodecoord)
        print *, 'znode', i, znodecoord(i)
    end do

    xnodect = size(xnodecoord)
    ynodect = size(ynodecoord)
    znodect = size(znodecoord)

    open (1, file = fn)

    write (1, *) 'CONTROL INFO 2.4.6'
    write (1, *) '** GAMBIT NEUTRAL FILE'
    write (1, *) fn
    write (1, *) 'PROGRAM:                Gambit     VERSION:  2.4.6'
    write (1, *) '21 Nov 2012    14:27:59'
    write (1, *) 'NUMNP     NELEM     NGRPS    NBSETS     NDFCD     NDFVL'
    write (1, *) xnodect * ynodect * znodect, (xnodect - 1)*(ynodect - 1)*(znodect - 1), '       1', NBSETS, '       3', '       3'
    write (1, *) 'ENDOFSECTION'

    ! assign each node with an integer
    allocate(nodescan(xnodect, ynodect, znodect), cellscan(xnodect - 1, ynodect - 1, znodect - 1))

    nodescan = 0
    allnodect = 0
    ! write: node index, x, y, z in the sequence of edge, face, body
    write (1, *) 'NODAL COORDINATES 2.4.6'

    !write the first 4 edges along z coordinate
    do loop1 = 1, 4
        if (loop1 == 1) then
            i = 1
            j = ynodect
        else if (loop1 == 2) then
            i = 1
            j = 1
        else if (loop1 == 3) then
            i = xnodect
            j = 1
        else
            i = xnodect
            j = ynodect
        end if

        k = znodect
        if (nodescan(i, j, k) == 0) then
            allnodect = allnodect + 1
            write (1, *) allnodect, xnodecoord(i), ynodecoord(j), znodecoord(k) !i,j,k !
            nodescan(i, j, k) = allnodect
        end if

        k = 1
        if (nodescan(i, j, k) == 0) then
            allnodect = allnodect + 1
            write (1, *) allnodect, xnodecoord(i), ynodecoord(j), znodecoord(k) !i,j,k !
            nodescan(i, j, k) = allnodect
        end if

        do loop2 = 1, znodect - 2
            k = znodect - loop2
            if (nodescan(i, j, k) == 0) then
                allnodect = allnodect + 1
                write (1, *) allnodect, xnodecoord(i), ynodecoord(j), znodecoord(k) !i,j,k !
                nodescan(i, j, k) = allnodect
            end if
        end do
    end do

    print *, '4 edges', allnodect

    !the 5th edge
    do loop2 = 1, xnodect - 2
        i = xnodect - loop2
        j = ynodect
        k = znodect
        if (nodescan(i, j, k) == 0) then
            allnodect = allnodect + 1
            write (1, *) allnodect, xnodecoord(i), ynodecoord(j), znodecoord(k) !i,j,k !
            nodescan(i, j, k) = allnodect
        end if
    end do

    !the 6th edge
    do loop2 = 1, xnodect - 2
        i = 1 + loop2
        j = ynodect
        k = 1
        if (nodescan(i, j, k) == 0) then
            allnodect = allnodect + 1
            write (1, *) allnodect, xnodecoord(i), ynodecoord(j), znodecoord(k) !i,j,k !
            nodescan(i, j, k) = allnodect
        end if
    end do

    print *, '6 edges', allnodect

    !the 7th edge
    do loop2 = 1, ynodect - 2
        i = xnodect
        j = ynodect - loop2
        k = 1
        if (nodescan(i, j, k) == 0) then
            allnodect = allnodect + 1
            write (1, *) allnodect, xnodecoord(i), ynodecoord(j), znodecoord(k) !i,j,k !
            nodescan(i, j, k) = allnodect
        end if
    end do

    !the 8th edge
    do loop2 = 1, ynodect - 2
        i = xnodect
        j = 1 + loop2
        k = znodect
        if (nodescan(i, j, k) == 0) then
            allnodect = allnodect + 1
            write (1, *) allnodect, xnodecoord(i), ynodecoord(j), znodecoord(k) !i,j,k !
            nodescan(i, j, k) = allnodect
        end if
    end do

    print *, '8 edges', allnodect

    !the 9th edge
    do loop2 = 1, xnodect - 2
        i = 1 + loop2
        j = 1
        k = znodect
        if (nodescan(i, j, k) == 0) then
            allnodect = allnodect + 1
            write (1, *) allnodect, xnodecoord(i), ynodecoord(j), znodecoord(k) !i,j,k !
            nodescan(i, j, k) = allnodect
        end if
    end do

    !the 10h edge
    do loop2 = 1, xnodect - 2
        i = xnodect - loop2
        j = 1
        k = 1
        if (nodescan(i, j, k) == 0) then
            allnodect = allnodect + 1
            write (1, *) allnodect, xnodecoord(i), ynodecoord(j), znodecoord(k) !i,j,k !
            nodescan(i, j, k) = allnodect
        end if
    end do

    !the 11th edge
    do loop2 = 1, ynodect - 2
        i = 1
        j = 1 + loop2
        k = 1
        if (nodescan(i, j, k) == 0) then
            allnodect = allnodect + 1
            write (1, *) allnodect, xnodecoord(i), ynodecoord(j), znodecoord(k) !i,j,k !
            nodescan(i, j, k) = allnodect
        end if
    end do

    !the 12th edge
    do loop2 = 1, ynodect - 2
        i = 1
        j = ynodect - loop2
        k = znodect
        if (nodescan(i, j, k) == 0) then
            allnodect = allnodect + 1
            write (1, *) allnodect, xnodecoord(i), ynodecoord(j), znodecoord(k) !i,j,k !
            nodescan(i, j, k) = allnodect
        end if
    end do

    print *, '12 edges', allnodect

    !the 1st face at Zmax
    do loop1 = 1, xnodect
        do loop2 = 1, ynodect
            i = loop1
            j = loop2
            k = znodect
            if (nodescan(i, j, k) == 0) then
                allnodect = allnodect + 1
                write (1, *) allnodect, xnodecoord(i), ynodecoord(j), znodecoord(k) !i,j,k !
                nodescan(i, j, k) = allnodect
            end if
        end do
    end do

    print *, 'face1', allnodect

    !the 2nd face at Zmin
    do loop1 = 1, xnodect
        do loop2 = 1, ynodect
            i = loop1
            j = loop2
            k = 1
            if (nodescan(i, j, k) == 0) then
                allnodect = allnodect + 1
                write (1, *) allnodect, xnodecoord(i), ynodecoord(j), znodecoord(k) !i,j,k !
                nodescan(i, j, k) = allnodect
            end if
        end do
    end do

    print *, 'face2', allnodect

    !the 3rd face at Ymin
    do loop1 = 1, xnodect
        do loop2 = 1, znodect
            i = xnodect + 1 - loop1
            j = 1
            k = znodect + 1 - loop2
            if (nodescan(i, j, k) == 0) then
                allnodect = allnodect + 1
                write (1, *) allnodect, xnodecoord(i), ynodecoord(j), znodecoord(k) !i,j,k !
                nodescan(i, j, k) = allnodect
            end if
        end do
    end do

    print *, 'face3', allnodect

    !the 4th face at Xmin
    do loop1 = 1, ynodect
        do loop2 = 1, znodect
            i = 1
            j = loop1
            k = znodect + 1 - loop2
            if (nodescan(i, j, k) == 0) then
                allnodect = allnodect + 1
                write (1, *) allnodect, xnodecoord(i), ynodecoord(j), znodecoord(k) !i,j,k !
                nodescan(i, j, k) = allnodect
            end if
        end do
    end do

    print *, 'face4', allnodect

    !the 5th face at Ymax
    do loop1 = 1, xnodect
        do loop2 = 1, znodect
            i = loop1
            j = ynodect
            k = znodect + 1 - loop2
            if (nodescan(i, j, k) == 0) then
                allnodect = allnodect + 1
                write (1, *) allnodect, xnodecoord(i), ynodecoord(j), znodecoord(k) !i,j,k !
                nodescan(i, j, k) = allnodect
            end if
        end do
    end do

    print *, 'face5', allnodect

    !the 6th face at Xmax
    do loop1 = 1, ynodect
        do loop2 = 1, znodect
            i = xnodect
            j = ynodect + 1 - loop1
            k = znodect + 1 - loop2
            if (nodescan(i, j, k) == 0) then
                allnodect = allnodect + 1
                write (1, *) allnodect, xnodecoord(i), ynodecoord(j), znodecoord(k) !i,j,k !
                nodescan(i, j, k) = allnodect
            end if
        end do
    end do

    print *, 'face6', allnodect

    !body
    do loop1 = 1, ynodect
        do loop2 = 1, xnodect
            do loop3 = 1, znodect
                i = xnodect + 1 - loop2
                j = loop1
                k = znodect + 1 - loop3
                if (nodescan(i, j, k) == 0) then
                    allnodect = allnodect + 1
                    write (1, *) allnodect, xnodecoord(i), ynodecoord(j), znodecoord(k) !i,j,k !
                    nodescan(i, j, k) = allnodect
                end if
            end do
        end do
    end do

    print *, 'body', allnodect
    write(1, *) 'ENDOFSECTION'

    !cell connectivity
    ! write: cell index, 4, 8, and the index of the 8 nodes belonging to the cell
    write(1, *) 'ELEMENTS/CELLS 2.4.6'

    cellscan = 0
    allcellct = 0
    do loop1 = 1, ynodect - 1
        do loop2 = 1, xnodect - 1
            do loop3 = 1, znodect - 1
                i = xnodect - loop2
                j = loop1
                k = znodect - loop3
                allcellct = allcellct + 1
                cellscan(i, j, k) = allcellct
                write (1, *) allcellct, '4  8', &
                nodescan(i + 1, j, k + 1), nodescan(i + 1, j, k), &
                nodescan(i, j, k + 1), nodescan(i, j, k), &
                nodescan(i + 1, j + 1, k + 1), nodescan(i + 1, j + 1, k), &
                nodescan(i, j + 1, k + 1), nodescan(i, j + 1, k)
            end do
        end do
    end do

    if (allcellct /= (xnodect - 1)*(ynodect - 1)*(znodect - 1)) then
        print *, 'allcellct is wrong', allcellct, (xnodect - 1)*(ynodect - 1)*(znodect - 1)
        stop
    end if

    print *, 'connectivity', allcellct
    write(1, *) 'ENDOFSECTION'

    !element group
    ! write: from 1 to size(cellscan)
    itype = 2
    write (1, *) 'ELEMENT GROUP 2.4.6 GROUP: 1 ELEMENTS: X MATERIAL: X'
    write (1, '(28x,I10)') allcellct
    write (1, '(28x,I10)') itype
    write (1, *) 'NFLAGS: 1', '    2    0'

    do i = 1, allcellct
        write (1, *) i
    end do
    write(1, *) 'ENDOFSECTION'

    !boundary condition
    do loop1 = 1, NBSETS
        ! count how many cells are at the boundary
        allbcct = 0
        if (BCX1 == loop1) allbcct = allbcct + (ynodect - 1)*(znodect - 1)
        if (BCX2 == loop1) allbcct = allbcct + (ynodect - 1)*(znodect - 1)
        if (BCY1 == loop1) allbcct = allbcct + (znodect - 1)*(xnodect - 1)
        if (BCY2 == loop1) allbcct = allbcct + (znodect - 1)*(xnodect - 1)
        if (BCZ1 == loop1) allbcct = allbcct + (xnodect - 1)*(ynodect - 1)
        if (BCZ2 == loop1) allbcct = allbcct + (xnodect - 1)*(ynodect - 1)

        ic = 1
        temp1 = 0
        temp2 = 6
        print *, 'BC count', loop1, allbcct
        write (1, *) 'BOUNDARY CONDITIONS 2.4.6'
        write (1, *) loop1, ic, allbcct, temp1, temp2

        allbcct = 0

        ! write: the index of cell on boundary, 4, and a number from 1 to 6
        if (BCZ2 == loop1) then
            do loop2 = 1, ynodect - 1
                do loop3 = 1, xnodect - 1
                    i = loop3
                    j = loop2
                    k = znodect - 1
                    allbcct = allbcct + 1
                    write (1, *) cellscan(i, j, k), '4  4'
                end do
            end do
        end if

        print *, 'bz2', allbcct

        if (BCY2 == loop1) then
            do loop2 = 1, znodect - 1
                do loop3 = 1, xnodect - 1
                    i = loop3
                    j = ynodect - 1
                    k = znodect - loop2
                    allbcct = allbcct + 1
                    write (1, *) cellscan(i, j, k), '4  6'
                end do
            end do
        end if

        print *, 'by2', allbcct

        if (BCX2 == loop1) then
            do loop2 = 1, znodect - 1
                do loop3 = 1, ynodect - 1
                    i = xnodect - 1
                    j = ynodect - loop3
                    k = znodect - loop2
                    allbcct = allbcct + 1
                    write (1, *) cellscan(i, j, k), '4  1'
                end do
            end do
        end if

        print *, 'bx2', allbcct

        if (BCX1 == loop1) then
            do loop2 = 1, znodect - 1
                do loop3 = 1, ynodect - 1
                    i = 1
                    j = loop3
                    k = znodect - loop2
                    allbcct = allbcct + 1
                    write (1, *) cellscan(i, j, k), '4  3'
                end do
            end do
        end if

        print *, 'bx1', allbcct

        if (BCY1 == loop1) then
            do loop2 = 1, znodect - 1
                do loop3 = 1, xnodect - 1
                    i = xnodect - loop3
                    j = 1
                    k = znodect - loop2
                    allbcct = allbcct + 1
                    write (1, *) cellscan(i, j, k), '4  5'
                end do
            end do
        end if

        print *, 'by1', allbcct

        if (BCZ1 == loop1) then
            do loop2 = 1, xnodect - 1
                do loop3 = 1, ynodect - 1
                    i = loop2
                    j = loop3
                    k = 1
                    allbcct = allbcct + 1
                    write (1, *) cellscan(i, j, k), '4  2'
                end do
            end do
        end if

        print *, 'bz1', allbcct

        write(1, *) 'ENDOFSECTION'
        print *, 'real allbcct', allbcct

    end do

    deallocate(xnodecoord, ynodecoord, znodecoord)
    deallocate(nodescan, cellscan)

    close(1)

    print *, 'end of program'

end Program meshproducer
