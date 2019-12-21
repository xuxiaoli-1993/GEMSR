MODULE dfd_library
    use dfd_type
    interface operator(*)
        module procedure matrix_by_matrix, &
        matrixs_by_matrixs, &
        matrix_by_vector, &
        matrixs_by_vectors, &
        matrix_by_array, &
        matrix_by_scale, &
        matrixs_by_scale, &
        vector_by_array, &
        vector_by_scale, &
        vectors_by_vectors, &
        vectors_by_scale
    end interface

    interface data_product
        module procedure matrix_product, &
        vector_product
    end interface

    interface operator(+)
        module procedure vector_add_vector, &
        vectors_add_vectors, &
        matrix_add_matrix, &
        matrixs_add_matrixs, &
        matrixs_add_scale
    end interface

    interface operator(-)
        module procedure vector_min_vector, &
        vectors_min_vectors, &
        matrix_min_matrix, &
        matrixs_min_matrixs, &
        matrixs_min_scale
    end interface

contains

    function matrix_by_matrix(am, bm)
        type(matrix), intent(in) :: am(:), bm(:)
        type(matrix) :: matrix_by_matrix(size(am))
        integer :: i
        do i = 1, size(am)
            matrix_by_matrix(i) % e = matmul(am(i) % e, bm(i) % e)
        end do
    end function matrix_by_matrix

    function matrixs_by_matrixs(am, bm)
        type(matrix), intent(in) :: am, bm
        type(matrix) :: matrixs_by_matrixs
        matrixs_by_matrixs % e = matmul(am % e, bm % e)
    end function matrixs_by_matrixs


    function matrix_by_vector(am, bv)
        type(matrix), intent(in) :: am(:)
        type(vector), intent(in) :: bv(:)
        type(vector) :: matrix_by_vector(size(bv))
        integer :: i
        do i = 1, size(bv)
            matrix_by_vector(i) % v = matmul(am(i) % e, bv(i) % v)
        end do
    end function matrix_by_vector

    function matrixs_by_vectors(am, bv)
        type(matrix), intent(in) :: am
        type(vector), intent(in) :: bv
        type(vector) :: matrixs_by_vectors
        matrixs_by_vectors % v = matmul(am % e, bv % v)
    end function matrixs_by_vectors

    function matrix_by_array(am, a)
        type(matrix), intent(in) :: am(:)
        type(matrix) :: matrix_by_array(size(am))
        real(kind = rfp), intent(in) :: a(:)
        integer :: i
        do i = 1, size(am)
            matrix_by_array(i) % e = am(i) % e * a(i)
        end do
    end function matrix_by_array

    function matrix_by_scale(am, s)
        type(matrix), intent(in) :: am(:)
        type(matrix) :: matrix_by_scale(size(am))
        real(kind = rfp), intent(in) :: s
        integer :: i
        if (s == zero) then
            do i = 1, size(am)
                matrix_by_scale(i) % e = zero
            end do
        else
            do i = 1, size(am)
                matrix_by_scale(i) % e = am(i) % e * s
            end do
        end if
    end function matrix_by_scale

    function matrixs_by_scale(am, s)
        type(matrix), intent(in) :: am
        type(matrix) :: matrixs_by_scale
        real(kind = rfp), intent(in) :: s
        if (s == zero) then
            matrixs_by_scale % e = zero
        else
            matrixs_by_scale % e = am % e * s
        end if
    end function matrixs_by_scale

    function vectors_by_vectors(av, bv)
        type(vector), intent(in) :: av, bv
        type(vector) :: vectors_by_vectors
        vectors_by_vectors % v = av % v * bv % v
    end function vectors_by_vectors

    function vector_by_array(av, a)
        type(vector), intent(in) :: av(:)
        type(vector) :: vector_by_array(size(av))
        real(kind = rfp), intent(in) :: a(:)
        integer :: i
        do i = 1, size(av)
            if (a(i) == zero) then
                vector_by_array(i) % v = zero
            else
                vector_by_array(i) % v = av(i) % v * a(i)
            end if
        end do
    end function vector_by_array

    function vector_by_scale(av, s)
        type(vector), intent(in) :: av(:)
        type(vector) :: vector_by_scale(size(av))
        real(kind = rfp), intent(in) :: s
        integer :: i
        if (s == zero) then
            do i = 1, size(av)
                vector_by_scale(i) % v = zero
            end do
        else
            do i = 1, size(av)
                vector_by_scale(i) % v = av(i) % v * s
            end do
        end if
    end function vector_by_scale

    function vectors_by_scale(av, s)
        type(vector), intent(in) :: av
        type(vector) :: vectors_by_scale
        real(kind = rfp), intent(in) :: s
        if (s == zero) then
            vectors_by_scale % v = zero
        else
            vectors_by_scale % v = av % v * s
        end if
    end function vectors_by_scale

    function vector_add_vector(av, bv)
        type(vector), intent(in) :: av(:), bv(:)
        type(vector) :: vector_add_vector(size(av))
        integer :: i
        do i = 1, size(av)
            vector_add_vector(i) % v = av(i) % v + bv(i) % v
        end do
    end function vector_add_vector

    function vectors_add_vectors(av, bv)
        type(vector), intent(in) :: av, bv
        type(vector) :: vectors_add_vectors
        vectors_add_vectors % v = av % v + bv % v
    end function vectors_add_vectors

    function matrix_add_matrix(am, bm)
        type(matrix), intent(in) :: am(:), bm(:)
        type(matrix) :: matrix_add_matrix(size(am))
        integer :: i
        do i = 1, size(am)
            matrix_add_matrix(i) % e = am(i) % e + bm(i) % e
        end do
    end function matrix_add_matrix

    function matrixs_add_matrixs(am, bm)
        type(matrix), intent(in) :: am, bm
        type(matrix) :: matrixs_add_matrixs
        matrixs_add_matrixs % e = am % e + bm % e
    end function matrixs_add_matrixs

    function matrixs_add_scale(am, s)
        type(matrix), intent(in) :: am
        real(kind = rfp), intent(in) :: s
        type(matrix) :: matrixs_add_scale
        matrixs_add_scale % e = am % e + s
    end function matrixs_add_scale


    function vector_min_vector(av, bv)
        type(vector), intent(in) :: av(:), bv(:)
        type(vector) :: vector_min_vector(size(av))
        integer :: i
        do i = 1, size(av)
            vector_min_vector(i) % v = av(i) % v - bv(i) % v
        end do
    end function vector_min_vector

    function vectors_min_vectors(av, bv)
        type(vector), intent(in) :: av, bv
        type(vector) :: vectors_min_vectors
        vectors_min_vectors % v = av % v - bv % v
    end function vectors_min_vectors

    function matrix_min_matrix(am, bm)
        type(matrix), intent(in) :: am(:), bm(:)
        type(matrix) :: matrix_min_matrix(size(am))
        integer :: i
        do i = 1, size(am)
            matrix_min_matrix(i) % e = am(i) % e - bm(i) % e
        end do
    end function matrix_min_matrix

    function matrixs_min_matrixs(am, bm)
        type(matrix), intent(in) :: am, bm
        type(matrix) :: matrixs_min_matrixs
        matrixs_min_matrixs % e = am % e - bm % e
    end function matrixs_min_matrixs

    function matrixs_min_scale(am, s)
        type(matrix), intent(in) :: am
        real(kind = rfp), intent(in) :: s
        type(matrix) :: matrixs_min_scale
        matrixs_min_scale % e = am % e - s
    end function matrixs_min_scale

    function matrix_product(am, a)
        type(matrix), intent(in) :: am(:)
        type(matrix) :: matrix_product
        real(kind = rfp), intent(in) :: a(:)
        integer :: i
        matrix_product % e = zero
        do i = 1, size(a)
            matrix_product % e = matrix_product % e + am(i) % e * a(i)
        end do
    end function matrix_product

    function vector_product(av, bv)
        type(vector), intent(in) :: av(:), bv(:)
        type(vector) :: vector_product
        integer :: i
        vector_product % v = zero
        do i = 1, size(av)
            vector_product % v = vector_product % v + av(i) % v * bv(i) % v
        end do
    end function vector_product
    !
    !
    !Computing areas and volumes are important steps for any geometry.
    !The area of any polygon can be computed by decomposing
    !the polygon into a series of triangles.
    !The volume of a tetrahedron is also easily computed from the coordinates
    !of its vertices. 
    !
    ! volume of tetra
    function volume3d(pc, p1, p2, p3)
        real(rfp), intent(in) :: pc(ndim), p1(ndim), p2(ndim), p3(ndim)
        real(rfp) :: volume3d, x1, y1, z1, x2, y2, z2, x3, y3, z3

        x1 = p1(1) - pc(1)
        y1 = p1(2) - pc(2)
        z1 = p1(ndim) - pc(ndim)
        !
        x2 = p2(1) - pc(1)
        y2 = p2(2) - pc(2)
        z2 = p2(ndim) - pc(ndim)
        !
        x3 = p3(1) - pc(1)
        y3 = p3(2) - pc(2)
        z3 = p3(ndim) - pc(ndim)
        !
        !   volume3d =  x1 * y2 * z3 + x2 * y3 * z1 + x3 * z2 * y1   &
        !             - z1 * y2 * x3 - z2 * y3 * x1 - z3 * x2 * y1
        !
        volume3d = (x1 * z3 - z1 * x3) * y2 &
        +(x2 * z1 - z2 * x1) * y3 &
        +(x3 * z2 - z3 * x2) * y1


        volume3d = volume3d / 6.0_rfp

    end function volume3d

    ! area of triangle
    !
    function volume2d(p1, p2, p3)
        real(rfp), intent(in) :: p1(2), p2(2), p3(2)
        real(rfp) :: volume2d
        volume2d = p1(1) * p2(2) - p2(1) * p1(2) + &
        p2(1) * p3(2) - p3(1) * p2(2) + &
        p3(1) * p1(2) - p1(1) * p3(2)
        volume2d = volume2d * 0.5_rfp

    end function volume2d

    ! normal vector of 3d triangle
    !
    function area3d(p1, p2, p3)result(area)
        real(rfp), intent(in) :: p1(ndim), p2(ndim), p3(ndim)
        real(rfp) :: area(ndim), x1, y1, z1, x2, y2, z2, x3, y3, z3
        x1 = p1(1)
        y1 = p1(2)
        z1 = p1(ndim)
        !
        x2 = p2(1)
        y2 = p2(2)
        z2 = p2(ndim)
        !
        x3 = p3(1)
        y3 = p3(2)
        z3 = p3(ndim)

        area(1) = (y2 - y1) * (z3 - z1) - (z2 - z1) * (y3 - y1)
        area(2) = (z2 - z1) * (x3 - x1) - (x2 - x1) * (z3 - z1)
        area(ndim) = (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1)
        !
        area = area * 0.5_rfp
        !
    end function area3d
    !
    ! normal vector of 2d line
    !
    function area2d(p1, p2)
        real(rfp), intent(in) :: p1(2), p2(2)
        real(rfp) :: area2d(2)
        area2d(1) = p2(2) - p1(2)
        area2d(2) = p1(1) - p2(1)
    end function area2d

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! linear solution subroutines section
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !
    subroutine normalize(b, d, c)
        !
        ! Solve Bx = d or BX = c
        !
        implicit none
        real(kind = rfp), intent(inout) :: b(:,:)
        real(kind = rfp), intent(inout) :: d(:)
        real(kind = rfp), intent(inout), optional :: c(:,:)
        integer :: i, n, ipivot(size(d))
        n = size(b, 1)
        !
        !  Compute LU = b
        call ludecomp(b, ipivot)
        !
        ! d = B_inv * D
        call lusolve(b, d, ipivot)
        ! C = B_inv * C
        if (present(c))then
            do i = 1, n
                call lusolve(b, c(:, i), ipivot)
            end do
        end if
    end subroutine normalize
    !
    !****************************************************************
    !Using Crout's method to compute LU = A
    !****************************************************************
    subroutine ludecomp(ul, ips)
        implicit none
        real(kind = rfp), intent (inout) :: ul(:,:) ! Assumed-shape array arguments
        integer, intent (out) :: ips(:)
        !
        !Result LU will replace original A
        !
        integer :: i, n, iswap, k, ip, kp
        integer :: idxpiv
        real(kind = rfp) :: big, em, pivot
        !
        n = size(ul, 1)
        ! Perform factorization A = L * U with partial pivoting.
        !
        do i = 1, n; ips(i) = i; end do
            do k = 1, n - 1
                big = 0.0
                do i = k, n
                    ip = ips(i)
                    if (abs(ul(ip, k)) > big) then
                        big = abs(ul(ip, k))
                        idxpiv = i
                    end if
                end do
                if (idxpiv /= k) then
                    iswap = ips(k)
                    ips(k) = ips(idxpiv)
                    ips(idxpiv) = iswap
                end if
                kp = ips(k)
                if (abs(ul(kp, k)) <= mytiny) ul(kp, k) = mytiny
                pivot = one / ul(kp, k)
                do i = k + 1, n
                    ip = ips(i)
                    em = -ul(ip, k) * pivot
                    ul(ip, k) = -em
                    ul(ip, k + 1:n) = ul(ip, k + 1:n) + em * ul(kp, k + 1:n)
                end do
            end do
            return
        end subroutine ludecomp
        !
        !***************************************************************
        !Solve the linear system A x = b, using the LU decomposition of 
        ! A stored in LU. Results return to replace b
        !***************************************************************
        subroutine lusolve(ul, b, ips)
            implicit none
            real(kind = rfp), intent(in) :: ul(:,:)
            real(kind = rfp), intent(inout) :: b(:)
            real(kind = rfp) :: x(size(b))
            integer, intent(in) :: ips(:)
            integer :: ip, n, i
            n = size(ul, 1)
            ip = ips(1)
            x(1) = b(ip)
            do i = 2, n
                ip = ips(i)
                x(i) = b(ip) - dot_product(ul(ip, 1:i - 1), x(1:i - 1))
            end do
            ip = ips(n)
            b(n) = x(n)/ul(ip, n)
            do i = n - 1, 1, -1
                ip = ips(i)
                b(i) = (x(i) - dot_product(ul(ip, i + 1:n), b(i + 1:n)))/ul(ip, i)
            end do
            return
        end subroutine lusolve
        !
        subroutine matrixinv(am)
            real(rfp), intent(inout) :: am(:,:)
            real(rfp) :: bm(size(am, 1), size(am, 1))
            integer :: ipivot(size(am, 1)), n, i
            call ludecomp(am, ipivot)
            n = size(am, 1)
            bm = zero
            do i = 1, n
                bm(i, i) = one
                call lusolve(am, bm(:, i), ipivot)
            end do
            am = bm
        end subroutine matrixinv

        subroutine direct_inv(a)
            ! inverse for matrix 2X2 or 3X3
            real(rfp), intent(inout) :: a(:,:)
            real(rfp) :: det, am(size(a, 1), size(a, 2)), pivot(size(a, 1))
            integer :: i
            am = a
            pivot = maxval(abs(am), dim = 1)
            if (any(pivot == zero)) then
                print *, 'singular matrix', am
                return
            end if
            do i = 1, size(a, 1)
                am(:, i) = am(:, i) / pivot(i)
            end do
            select case(size(a, 1))
            case(2)
                det = am(1, 1) * am(2, 2) - am(1, 2) * am(2, 1)
                if (abs(det) < mytiny) then
                    a = zero
                    !  print *,'det = 0 in direct_inv'
                    return
                end if
                det = 1.0_rfp / det
                a(1, 1) = am(2, 2) * det
                a(1, 2) = -am(1, 2) * det
                a(2, 1) = -am(2, 1) * det
                a(2, 2) = am(1, 1) * det
            case(3)
                det = am(1, 1)*(am(2, 2) * am(3, 3) - am(2, 3) * am(3, 2)) - &
                am(2, 1)*(am(1, 2) * am(3, 3) - am(1, 3) * am(3, 2)) + &
                am(3, 1)*(am(1, 2) * am(2, 3) - am(2, 2) * am(1, 3))
                if (abs(det) < mytiny) then
                    a = zero
                    !  print *,'det = 0 in direct_inv'
                    return
                end if
                det = 1.0_rfp / det
                a(1, 1) = (am(2, 2) * am(3, 3) - am(2, 3) * am(3, 2)) * det
                a(1, 2) = -(am(1, 2) * am(3, 3) - am(1, 3) * am(3, 2)) * det
                a(1, 3) = (am(1, 2) * am(2, 3) - am(2, 2) * am(1, 3)) * det
                a(2, 1) = -(am(2, 1) * am(3, 3) - am(3, 1) * am(2, 3)) * det
                a(2, 2) = (am(1, 1) * am(3, 3) - am(1, 3) * am(3, 1)) * det
                a(2, 3) = -(am(1, 1) * am(2, 3) - am(2, 1) * am(1, 3)) * det
                a(3, 1) = (am(2, 1) * am(3, 2) - am(2, 2) * am(3, 1)) * det
                a(3, 2) = -(am(1, 1) * am(3, 2) - am(1, 2) * am(3, 1)) * det
                a(3, 3) = (am(1, 1) * am(2, 2) - am(1, 2) * am(2, 1)) * det
            end select

            do i = 1, size(a, 2)
                a(i,:) = a(i,:) / pivot(i)
            end do

        end subroutine direct_inv

        function cramer(a, b) result(x)
            implicit none
            real(rfp), intent(in) :: a(:,:), b(:)
            real(rfp) :: mya(size(b), size(b)), x(size(b)), det
            integer :: i
            !
            ! Cramer's rule is a determinant-based procedure 
            !  utilized to solve systems of equations 
            ! This sub is only desigened for 2X2 or 3X3 equations
            !
            det = determinant(a)
            !   if(abs(det) < mytiny) return 
            do i = 1, size(b)
                mya = a
                mya(:, i) = b
                x(i) = determinant(mya) / det
            end do
        end function cramer

        function determinant(a)result(det)
            real(rfp), intent(in) :: a(:,:)
            real(rfp) :: det
            select case(size(a, dim = 1))
            case(2)
                Det = a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1)
            case(3)
                Det = a(1, 1) *(a(2, 2) * a(3, 3) - a(2, 3) * a(3, 2)) &
                -a(1, 2) *(a(2, 1) * a(3, 3) - a(2, 3) * a(3, 1)) &
                +a(1, 3) *(a(2, 1) * a(3, 2) - a(2, 2) * a(3, 1))
            case default
                print *, 'error for determinant', size(a)
            end select
        end function determinant
        !
        !
        !*************************************************************
        ! Solve block tridiadinal linear equations backforward
        !*************************************************************
        subroutine block_tri_solve(rhs, am, bm, cm)
            type(vector), intent(inout) :: rhs(:)
            type(matrix) :: am(:), bm(:), cm(:)
            integer :: i, n, im1, ip1
            !
            n = size(rhs)
            call normalize(bm(1) % e, rhs(1) % v, cm(1) % e)
            do i = 2, n - 1
                im1 = i - 1
                bm(i) % e = bm(i) % e - matmul(am(i) % e, cm(im1) % e)
                rhs(i) % v = rhs(i) % v - matmul(am(i) % e, rhs(im1) % v)
                call normalize(bm(i) % e, rhs(i) % v, cm(i) % e)
            end do
            bm(n) % e = bm(n) % e - matmul(am(n) % e, cm(n - 1) % e)
            rhs(n) % v = rhs(n) % v - matmul(am(n) % e, rhs(n - 1) % v)
            call normalize(bm(n) % e, rhs(n) % v)
            do i = n - 1, 1, -1
                ip1 = i + 1
                rhs(i) % v = rhs(i) % v - matmul(cm(i) % e, rhs(ip1) % v)
            end do
        end subroutine block_tri_solve
        !
        !
        !******************************************************************************
        !  Given a matrix A, with dimensions M by N this
        !  routine computes its singular value decomposition, A = U.W.V'. The matrix U replaces 
        !  A on output. The diagonal matrix of singular values W is output as a vector W. The matrix 
        !  V (not the Transpose V') is output as V. M must be greater or equal to N; if it is smaller,
        !  then A should be filled up to square with zero rows.

        SUBROUTINE svdcmp_inv(a, ierr)
            implicit none
            INTEGER :: m, n, i, j, l, k, its, nm, jj, ierr
            REAL(rfp) :: g, sscale, anorm, s, f, h, x, z, y, c
            REAL(rfp), parameter :: EPS = 1.0e-10_rfp
            REAL(rfp) :: A(:,:), w(size(a, 2)), v(size(a, 2), size(a, 2)), rv1(size(a, 2))
            m = size(a, 1)
            n = size(a, 2)
            if (m .lt. n) then
                PRINT *, 'You must augment A with extra zero rows'
                call exit(10)
            ENDIF
            !Householder Reduction to bidiagonal form
            !(see Forsythe,Malcolm,Moler, "Computer Methods for Mathematical Computations"
            g = 0.0_rfp
            sscale = 0.0_rfp
            anorm = 0.0_rfp
            do i = 1, n
                l = i + 1
                rv1(i) = sscale * g
                g = 0.0_rfp
                s = 0.0_rfp
                sscale = 0.0_rfp
                if (i .le. m) then
                    sscale = sum(abs(a(i:m, i)))
                    !         if (sscale.ne.0.0_rfp) then
                    if (abs(sscale - 0.0_rfp) .gt. EPS) then
                        a(i:m, i) = a(i:m, i) / sscale
                        s = dot_product(a(i:m, i), a(i:m, i))
                        f = a(i, i)
                        g = -SIGN(SQRT(s), f)
                        h = f * g - s
                        a(i, i) = f - g
                        if (i .ne. n) then
                            do j = l, n
                                s = dot_product(a(i:m, i), a(i:m, j))
                                f = s / h
                                a(i:m, j) = a(i:m, j) + f * a(i:m, i)
                            end do ! j loop
                        end if
                        a(i:m, i) = sscale * a(i:m, i)
                    end if
                end if

                w(i) = sscale * g
                g = 0.0_rfp
                s = 0.0_rfp
                sscale = 0.0_rfp
                if ((i .le. m).AND.(i .ne. n)) then
                    sscale = sum(abs(a(i, l:n)))
                    !       if (sscale.ne.0.0_rfp) then
                    if (abs(sscale - 0.0_rfp) .gt. EPS) then
                        a(i, l:n) = a(i, l:n) /sscale
                        s = dot_product(a(i, l:n), a(i, l:n))
                        f = a(i, l)
                        g = -SIGN(SQRT(s), f)
                        h = f * g - s
                        a(i, l) = f - g
                        rv1(l:n) = a(i, l:n) / h
                        if (i .ne. m) then
                            do j = l, m
                                s = dot_product(a(j, l:n), a(i, l:n))
                                a(j, l:n) = a(j, l:n) + s * rv1(l:n)
                            end do ! j loop
                        end if
                        a(i, l:n) = sscale * a(i, l:n)
                    end if
                end if
                anorm = MAX(anorm, (abs(w(i)) + abs(rv1(i))))
            end do

            ! Accumulation of right-hand Transformations
            do i = n, 1, -1
                if (i .lt. n) then
                    !         if (g.ne.0.0_rfp) then
                    if (abs(g - 0.0_rfp) .gt. EPS) then
                        do j = l, n ! Double division to avoid possible overflow
                            v(j, i) = (a(i, j) / a(i, l)) / g
                        end do ! j loop
                        do j = l, n
                            s = dot_product(a(i, l:n), v(l:n, j))
                            v(l:n, j) = v(l:n, j) + s * v(l:n, i)
                        end do ! j loop
                    end if
                    v(i, l:n) = 0.0_rfp
                    v(l:n, i) = 0.0_rfp
                end if
                v(i, i) = 1.0_rfp
                g = rv1(i)
                l = i
            end do

            ! Accumulation of left-hand Transformations
            do i = n, 1, -1
                l = 1 + i
                g = w(i)
                if (i .lt. n) then
                    a(i, l:n) = 0.0_rfp
                end if
                !      if (g.ne.0.0_rfp) then
                if (abs(g - 0.0_rfp) .gt. EPS) then
                    g = 1.0_rfp / g
                    if (i .ne. n) then
                        do j = l, n
                            s = dot_product(a(l:m, i), a(l:m, j))
                            f = (s/a(i, i)) * g
                            a(i:m, j) = a(i:m, j) + f * a(i:m, i)
                        end do ! j loop
                    end if
                    a(i:m, i) = a(i:m, i) * g
                else
                    a(i:m, i) = 0.0_rfp
                end if
                a(i, i) = a(i, i) + 1.0_rfp
            end do ! i loop

            ! Diagonalization of the bidigonal form
            do k = n, 1, -1 !Loop over singular values
                do its = 1, 30 !Loop over allowed iterations
                    do l = k, 1, -1 !Test for splitting
                        nm = l - 1 ! Note that rv1(1) is always zero
                        !           if ( (abs(rv1(l))+anorm) .eq. anorm ) GO TO 2
                        !          	if ( (abs(w(nm))+anorm) .eq. anorm ) GO TO 1
                        if (abs((abs(rv1(l)) + anorm) - anorm) .lt. eps) GO TO 2
                        if (abs((abs(w(nm)) + anorm) - anorm) .lt. eps) GO TO 1
                    end do !  l loop

                    1 c = 0.0_rfp ! Cancellation of rv1(l), if l>1 :
                    s = 1.0_rfp
                    do i = l, k
                        f = s * rv1(i)
                        !            if ( (abs(f)+anorm) .ne. anorm ) then
                        if (abs((abs(f) + anorm) - anorm) .GT. eps) then

                            g = w(i)
                            h = SQRT(f * f + g * g)
                            w(i) = h
                            h = 1.0_rfp / h
                            c = g * h
                            s = -f * h
                            do j = 1, m
                                y = a(j, nm)
                                z = a(j, i)
                                a(j, nm) = (y * c) + (z * s)
                                a(j, i) = -(y * s) + (z * c)
                            end do ! j loop
                        end if
                    end do ! i loop
                    2 z = w(k)
                    if (l .eq. k) then ! convergence
                        if (z .lt. 0.0_rfp) then ! Singular value is made non-negative
                            w(k) = -z
                            v(1:n, k) = -v(1:n, k)
                        end if
                        GO TO 3
                    end if
                    if (its .eq. 30) then
                        !          PRINT*, 'No Convergence in 30 iterations'
                        ierr = ierr + 1
                        call exit(10)
                    ENDIF
                    x = w(l) ! Shift from bottom 2-by-2 minor
                    nm = k - 1
                    y = w(nm)
                    g = rv1(nm)
                    h = rv1(k)
                    f = ((y - z)*(y + z) + (g - h)*(g + h)) / (2.0_rfp * h * y)
                    g = SQRT(f * f + 1.0_rfp)
                    f = ((x - z)*(x + z) + h * ((y/(f + SIGN(g, f))) - h)) / x

                    ! Next   QR Transformation
                    c = 1.0_rfp
                    s = 1.0_rfp
                    do j = l, nm
                        i = j + 1
                        g = rv1(i)
                        y = w(i)
                        h = s * g
                        g = c * g
                        z = SQRT(f * f + h * h)
                        rv1(j) = z
                        c = f/z
                        s = h/z
                        f = (x * c) + (g * s)
                        g = -(x * s) + (g * c)
                        h = y * s
                        y = y * c
                        do jj = 1, n
                            x = v(jj, j)
                            z = v(jj, i)
                            v(jj, j) = (x * c) + (z * s)
                            v(jj, i) = -(x * s) + (z * c)
                        end do
                        z = SQRT(f * f + h * h)
                        w(j) = z
                        !            if (z.ne.0.0_rfp) then
                        if (abs(z - 0.0_rfp) .gt. eps) then
                            z = 1.0_rfp / z
                            c = f * z
                            s = h * z
                        end if
                        f = (g * c) + (y * s)
                        x = -(g * s) + (y * c)
                        do jj = 1, m
                            y = a(jj, j)
                            z = a(jj, i)
                            a(jj, j) = (y * c) + (z * s)
                            a(jj, i) = -(y * s) + (z * c)
                        end do
                    end do ! j loop
                    rv1(l) = 0.0_rfp
                    rv1(k) = f
                    w(k) = x
                end do ! its loop
                3 continue
            end do ! k loop
            ! generate 
            ![A.w.v^T]^-1 = [A.1/w.v^T]^T
            do i = 1, size(w)
                if (w(i) /= zero) then
                    v(:, i) = v(:, i) / w(i)
                else
                    v(:, i) = zero
                end if
            end do
            !  A (not the Transpose A')
            a = matmul(a, transpose(v))
            !
            !  x = matmul(b,a)
            return
        END SUBROUTINE svdcmp_inv

        !
        subroutine timera(icntrl, chr, time)
            ! this subroutine performs timing
            ! input: icntrl, chr 
            ! icntrl = (-1,0,1) = (initialize,ignore,read) clock
            ! clock should be initialized before it is read!
            ! chr = character variable for labeling timings
            ! time = elapsed time in seconds
            implicit none
            character*8 chr
            real*4 dtime, tarry, jclock, nclock, time
            dimension tarry(2)
            integer :: icntrl
            save jclock
            91 format (1x, a8, 1x, 'time = ', e14.7, 1x, 'sec')
            data jclock /0./
            if (icntrl .eq. 0) then
                time = dtime(tarry)
                return
            end if
            if (icntrl .eq. 1) go to 10
            ! initialize clock
            jclock = dtime(tarry)
            return
            ! read clock and write time difference from last clock initialization
            10 nclock = dtime(tarry)
            time = nclock
            write (6, 91) chr, time
            return
        end subroutine timera

        function polynomial(x, a)result(p)
            real(rfp), intent(in) :: x, a(:)
            real(rfp) :: p
            integer :: i, n
            n = size(a)
            if (n <= 0) then
                p = zero
            else
                p = a(n)
                do i = n - 1, 1, -1
                    p = x * p + a(i)
                end do
            end if

        end function polynomial

        subroutine sort(x, is)
            real(rfp) :: x(:), y(size(x)), s
            integer, optional :: is(:)
            integer :: i, j, n, k, ib
            n = size(x)
            !   is = (/(i,i=1,n)/)
            y = x
            do j = 2, n
                s = y(j)
                ib = is(j)
                do i = j - 1, 1, -1
                    if (y(i) < s) exit
                    y(i + 1) = y(i)
                    if (present(is)) is(i + 1) = is(i)
                end do
                y(i + 1) = s
                if (present(is)) is(i + 1) = ib
            end do
            x = y
        end subroutine sort

        function i2s(i)
            character(len = 4) :: i2s
            integer, intent(in) :: i
            i2s = ''
            select case(i)
            case(0:9)
                write(i2s, '(i1)') i
            case(10:99)
                write(i2s, '(i2)') i
            case(100:999)
                write(i2s, '(i3)') i
            case(1000:9999)
                write(i2s, '(i4)') i
            case default
                write(i2s, '(4h****)')
            end select
        end function i2s

        function i2s_pgi(i)result(i2s)
            character(len = 4) :: i2s, c(0:9)
            integer, intent(in) :: i
            integer :: j, k, m, n
            data c/'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'/
            i2s = ''
            select case(i)
            case(0:9)
                i2s = c(i)
            case(10:99)
                j = i / 10
                k = i - j * 10
                i2s = c(j) // c(k)
            case(100:999)
                j = i / 100
                k = (i - j * 100)/10
                m = i - j * 100 - k * 10
                i2s = c(j) // c(k) // c(m)
            case(1000:9999)
                j = i / 1000
                k = (i - j * 1000)/100
                m = (i - j * 1000 - k * 100)/10
                n = i - j * 1000 - k * 100 - m * 10
                i2s = c(j) // c(k) // c(m) // c(n)
            end select
        end function i2s_pgi

        function snrm2(x)
            real(rfp), intent(in) :: x(:)
            real(rfp) :: snrm2, scale, ssq
            integer :: n, i
            n = size(x)
            if (n < 1) then
                snrm2 = zero
            else if (n == 1) then
                snrm2 = abs(x(1))
            else
                scale = zero
                ssq = one
                do i = 1, n
                    absx = abs(x(i))
                    if (absx > mytiny) then
                        if (scale < absx) then
                            ssq = one + ssq * (scale / absx) **2
                            scale = absx
                        else
                            ssq = ssq + (absx / scale) **2
                        end if
                    end if
                end do
                snrm2 = scale * sqrt(ssq)
            end if
        end function snrm2

        subroutine put_diag(mat, diag, i1, i2)
            real(rfp) :: diag, mat(:,:)
            integer :: i, n, i1, i2
            ! n = size(diag)
            do i = i1, i2
                mat(i, i) = diag
            end do
        end subroutine put_diag

        subroutine diagadd(mat, diag)
            real(rfp) :: diag(:), mat(:,:)
            integer :: i, n
            n = size(diag)
            if (n /= size(mat, 1).or. n /= size(mat, 2)) return
            do i = 1, n
                mat(i, i) = mat(i, i) + diag(i)
            end do
        end subroutine diagadd

        function diagsum(mat)result(s)
            real(rfp) :: mat(:,:), s
            integer :: i, n
            n = min(size(mat, 1), size(mat, 2))
            s = zero
            do i = 1, n
                s = s + mat(i, i)
            end do

        end function diagsum

        function cross2mat(v)result(m)
            real(rfp) :: v(3), m(3, 3)
            m = zero
            m(1, 2) = -v(3); m(2, 1) = v(3)
            m(1, 3) = v(2); m(3, 1) = -v(2)
            m(2, 3) = -v(1); m(3, 2) = v(1)
        end function cross2mat

        function vec2mat(v1, v2)result(m)
            type(vector) :: v1, v2
            type(matrix) :: m
            integer :: i, j
            do i = 1, size(v2 % v)
                do j = 1, size(v1 % v)
                    m % e(i, j) = v1 % v(j) * v2 % v(i)
                end do
            end do
        end function vec2mat

        function matcrossvec(m, v)result(mo)
            real(rfp) :: v(3), m(3, 3), mo(3, 3)
            mo(:, 1) = acrossb(m(:, 1), v)
            mo(:, 2) = acrossb(m(:, 2), v)
            mo(:, 3) = acrossb(m(:, 3), v)
        end function matcrossvec

        function curl(a)result(c)
            real(rfp) :: a(3, 3), c(3)
            c(1) = a(3, 2) - a(2, 3)
            c(2) = a(1, 3) - a(3, 1)
            c(3) = a(2, 1) - a(1, 2)
        end function curl

        function acrossb(a, b)result(c)
            real(rfp), intent(in) :: a(3), b(3)
            real(rfp) :: c(3), am, bm
            c(1) = a(2) * b(3) - a(3) * b(2)
            c(2) = a(3) * b(1) - a(1) * b(3)
            c(3) = a(1) * b(2) - a(2) * b(1)
        end function acrossb

        function avto3(v)result(a)
            real(rfp) :: v(:), a(3)
            integer :: n
            if (n == 3) then
                a = v
            else
                a = zero; n = size(v)
                a(:n) = v
            end if
        end function avto3

        function tensor(a, b)result(m)
            real(rfp) :: a(:), b(:), m(size(a), size(b))
            integer :: i
            do i = 1, size(b)
                m(:, i) = a(:) * b(i)
            end do
        end function tensor

        function eval_spline(xf, x, y, b, c, d)result(yf)
            integer :: n, i
            real(rfp) :: x(:), y(:), b(:), c(:), d(:), xf, yf, t
            n = size(x)
            if (xf <= x(1)) then
                yf = y(1)
                return
            else if (xf >= x(n)) then
                yf = y(n)
                return
            end if
            do i = 2, n
                if (xf < x(i)) exit
            end do
            i = i - 1
            t = xf - x(i)
            yf = y(i) + (b(i) + (c(i) + d(i) * t) * t) * t
            return
        end function eval_spline

        subroutine spline(x, y, b, c, d)
            integer :: n
            real(rfp) :: x(:), y(:), b(:), c(:), d(:)
            !
            !  the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
            !  for a cubic interpolating spline
            !
            !    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
            !
            !    for  x(i) .le. x .le. x(i+1)
            !
            !  input..
            !
            !    n = the number of data points or knots (n.ge.2)
            !    x = the abscissas of the knots in strictly increasing order
            !    y = the ordinates of the knots
            !
            !  output..
            !
            !    b, c, d  = arrays of spline coefficients as defined above.
            !
            !  using  p  to denote differentiation,
            !
            !    y(i) = s(x(i))
            !    b(i) = sp(x(i))
            !    c(i) = spp(x(i))/2
            !    d(i) = sppp(x(i))/6  (derivative from the right)
            !
            !  the accompanying function subprogram  seval  can be used
            !  to evaluate the spline.
            !
            !
            integer :: nm1, ib, i
            real(rfp) :: t
            !
            n = size(x)
            nm1 = n - 1
            if (n .lt. 2) return
            if (n .lt. 3) go to 50
            !
            !  set up tridiagonal system
            !
            !  b = diagonal, d = offdiagonal, c = right hand side.
            !
            d(1) = x(2) - x(1)
            c(2) = (y(2) - y(1))/d(1)
            do 10 i = 2, nm1
                d(i) = x(i + 1) - x(i)
                b(i) = 2. * (d(i - 1) + d(i))
                c(i + 1) = (y(i + 1) - y(i))/d(i)
                c(i) = c(i + 1) - c(i)
                10 continue
                !
                !  end conditions.  third derivatives at  x(1)  and  x(n)
                !  obtained from divided differences
                !
                b(1) = -d(1)
                b(n) = -d(n - 1)
                c(1) = 0.
                c(n) = 0.
                if (n .eq. 3) go to 15
                c(1) = c(3)/(x(4) - x(2)) - c(2)/(x(3) - x(1))
                c(n) = c(n - 1)/(x(n) - x(n - 2)) - c(n - 2)/(x(n - 1) - x(n - 3))
                c(1) = c(1) * d(1)**2/(x(4) - x(1))
                c(n) = -c(n) * d(n - 1)**2/(x(n) - x(n - 3))
                !
                !  forward elimination
                !
                15 do 20 i = 2, n
                t = d(i - 1)/b(i - 1)
                b(i) = b(i) - t * d(i - 1)
                c(i) = c(i) - t * c(i - 1)
                20 continue
                !
                !  back substitution
                !
                c(n) = c(n)/b(n)
                do 30 ib = 1, nm1
                    i = n - ib
                    c(i) = (c(i) - d(i) * c(i + 1))/b(i)
                    30 continue
                    !
                    !  c(i) is now the sigma(i) of the text
                    !
                    !  compute polynomial coefficients
                    !
                    b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + 2. * c(n))
                    do 40 i = 1, nm1
                        b(i) = (y(i + 1) - y(i))/d(i) - d(i)*(c(i + 1) + 2. * c(i))
                        d(i) = (c(i + 1) - c(i))/d(i)
                        c(i) = 3. * c(i)
                        40 continue
                        c(n) = 3. * c(n)
                        d(n) = d(n - 1)
                        return
                        !
                        50 b(1) = (y(2) - y(1))/(x(2) - x(1))
                        c(1) = 0.
                        d(1) = 0.
                        b(2) = b(1)
                        c(2) = 0.
                        d(2) = 0.
                        return
                    end subroutine spline
                    !  Bessel function



                    SUBROUTINE JY01A(X, BJ0, DJ0, BJ1, DJ1, BY0, DY0, BY1, DY1)
                        !
                        !       =======================================================
                        !       Purpose: Compute Bessel functions J0(x), J1(x), Y0(x),
                        !                Y1(x), and their derivatives
                        !       Input :  x   --- Argument of Jn(x) & Yn(x) ( x ò 0 )
                        !       Output:  BJ0 --- J0(x)
                        !                DJ0 --- J0'(x)
                        !                BJ1 --- J1(x)
                        !                DJ1 --- J1'(x)
                        !                BY0 --- Y0(x)
                        !                DY0 --- Y0'(x)
                        !                BY1 --- Y1(x)
                        !                DY1 --- Y1'(x)
                        !       =======================================================
                        !
                        IMPLICIT DOUBLE PRECISION (A - H, O - Z)
                        DIMENSION A(12), B(12), A1(12), B1(12)
                        RP2 = 0.63661977236758D0
                        X2 = X * X
                        IF (X .EQ. 0.0D0) THEN
                            BJ0 = 1.0D0
                            BJ1 = 0.0D0
                            DJ0 = 0.0D0
                            DJ1 = 0.5D0
                            BY0 = -1.0D+300
                            BY1 = -1.0D+300
                            DY0 = 1.0D+300
                            DY1 = 1.0D+300
                            RETURN
                        ENDIF
                        IF (X .LE. 12.0D0) THEN
                            BJ0 = 1.0D0
                            R = 1.0D0
                            DO 5 K = 1, 30
                                R = -0.25D0 * R * X2/(K * K)
                                BJ0 = BJ0 + R
                                IF (DABS(R) .LT. DABS(BJ0) * 1.0D-15) GO TO 10
                                5 CONTINUE
                                10 BJ1 = 1.0D0
                                R = 1.0D0
                                DO 15 K = 1, 30
                                    R = -0.25D0 * R * X2/(K * (K + 1.0D0))
                                    BJ1 = BJ1 + R
                                    IF (DABS(R) .LT. DABS(BJ1) * 1.0D-15) GO TO 20
                                    15 CONTINUE
                                    20 BJ1 = 0.5D0 * X * BJ1
                                    EC = DLOG(X/2.0D0) + 0.5772156649015329D0
                                    CS0 = 0.0D0
                                    W0 = 0.0D0
                                    R0 = 1.0D0
                                    DO 25 K = 1, 30
                                        W0 = W0 + 1.0D0/K
                                        R0 = -0.25D0 * R0/(K * K) * X2
                                        R = R0 * W0
                                        CS0 = CS0 + R
                                        IF (DABS(R) .LT. DABS(CS0) * 1.0D-15) GO TO 30
                                        25 CONTINUE
                                        30 BY0 = RP2 * (EC * BJ0 - CS0)
                                        CS1 = 1.0D0
                                        W1 = 0.0D0
                                        R1 = 1.0D0
                                        DO 35 K = 1, 30
                                            W1 = W1 + 1.0D0/K
                                            R1 = -0.25D0 * R1/(K * (K + 1)) * X2
                                            R = R1 * (2.0D0 * W1 + 1.0D0/(K + 1.0D0))
                                            CS1 = CS1 + R
                                            IF (DABS(R) .LT. DABS(CS1) * 1.0D-15) GO TO 40
                                            35 CONTINUE
                                            40 BY1 = RP2 * (EC * BJ1 - 1.0D0/X - 0.25D0 * X * CS1)
                                            ELSE
                                                DATA A/-.7031250000000000D-01, .1121520996093750D+00, &
                                                -.5725014209747314D+00, .6074042001273483D+01, &
                                                -.1100171402692467D+03, .3038090510922384D+04, &
                                                -.1188384262567832D+06, .6252951493434797D+07, &
                                                -.4259392165047669D+09, .3646840080706556D+11, &
                                                -.3833534661393944D+13, .4854014686852901D+15/
                                                DATA B/ .7324218750000000D-01, -.2271080017089844D+00, &
                                                .1727727502584457D+01, -.2438052969955606D+02, &
                                                .5513358961220206D+03, -.1825775547429318D+05, &
                                                .8328593040162893D+06, -.5006958953198893D+08, &
                                                .3836255180230433D+10, -.3649010818849833D+12, &
                                                .4218971570284096D+14, -.5827244631566907D+16/
                                                DATA A1/.1171875000000000D+00, -.1441955566406250D+00, &
                                                .6765925884246826D+00, -.6883914268109947D+01, &
                                                .1215978918765359D+03, -.3302272294480852D+04, &
                                                .1276412726461746D+06, -.6656367718817688D+07, &
                                                .4502786003050393D+09, -.3833857520742790D+11, &
                                                .4011838599133198D+13, -.5060568503314727D+15/
                                                DATA B1/-.1025390625000000D+00, .2775764465332031D+00, &
                                                -.1993531733751297D+01, .2724882731126854D+02, &
                                                -.6038440767050702D+03, .1971837591223663D+05, &
                                                -.8902978767070678D+06, .5310411010968522D+08, &
                                                -.4043620325107754D+10, .3827011346598605D+12, &
                                                -.4406481417852278D+14, .6065091351222699D+16/
                                                K0 = 12
                                                IF (X .GE. 35.0) K0 = 10
                                                IF (X .GE. 50.0) K0 = 8
                                                T1 = X - 0.25D0 * PI
                                                P0 = 1.0D0
                                                Q0 = -0.125D0/X
                                                DO 45 K = 1, K0
                                                    P0 = P0 + A(K) * X**(-2 * K)
                                                    45 Q0 = Q0 + B(K) * X**(-2 * K - 1)
                                                    CU = DSQRT(RP2/X)
                                                    BJ0 = CU * (P0 * DCOS(T1) - Q0 * DSIN(T1))
                                                    BY0 = CU * (P0 * DSIN(T1) + Q0 * DCOS(T1))
                                                    T2 = X - 0.75D0 * PI
                                                    P1 = 1.0D0
                                                    Q1 = 0.375D0/X
                                                    DO 50 K = 1, K0
                                                        P1 = P1 + A1(K) * X**(-2 * K)
                                                        50 Q1 = Q1 + B1(K) * X**(-2 * K - 1)
                                                        CU = DSQRT(RP2/X)
                                                        BJ1 = CU * (P1 * DCOS(T2) - Q1 * DSIN(T2))
                                                        BY1 = CU * (P1 * DSIN(T2) + Q1 * DCOS(T2))
                                                    ENDIF
                                                    DJ0 = -BJ1
                                                    DJ1 = BJ0 - BJ1/X
                                                    DY0 = -BY1
                                                    DY1 = BY0 - BY1/X
                                                    RETURN
                                                END subroutine jy01a

                                            end MODULE dfd_library

