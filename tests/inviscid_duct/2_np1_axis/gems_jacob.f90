module gems_jacob
   use gems_state
   implicit none
contains
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! All Jacobian Matrix Module
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !***********************************************************************
   !  calculate Gamma jacobian or Gamma prima
   !***********************************************************************
   function gamma(is, qv, rhopp, cc, bspeed, cspeed)
      implicit none
      type(cell), optional :: cc
      real(rfp), optional :: cspeed, bspeed
      integer, intent(in) :: is
      type(vector) :: qv
      real(rfp) :: rhopp
      type(vector) :: rhox, hx
      real(rfp) :: rho, h0
      type(matrix) :: gamma
      integer :: k

      ! initial gamma
      gamma % e = zero
      if (neqf > 0) then

         rho = rhof(qv)
         h0 = h0f(qv)
         rhox = rhoxf(qv)
         hx = hxf(qv)

         if (is > 0) rhox % v(ico) = rhopp ! precondition
         !  if(is > 0) rhox%v(ien) = three * rhox%v(ien)
         !  
         ! continuity row
         if (nceq > 0) gamma % e(ico,:) = rhox % v
         ! momentum row
         if (nmeq > 0) then
            do k = imb, ime
               gamma % e(k,:) = rhox % v * qv % v(k)
               gamma % e(k, k) = rho
            end do
         end if
         !energy row
         if (neeq > 0) then
            gamma % e(ien,:) = rhox % v(:) * h0 + rho * hx % v(:)
            if (nceq > 0) gamma % e(ien, ico) = gamma % e(ien, ico) - one
         end if
         !
         !turbulence and multi-spice rows
         if(neqf > ien) then
            do k = ien + 1, neqf
               gamma % e(k,:) = rhox % v * qv % v(k)
               gamma % e(k, k) = gamma % e(k, k) + rho
            end do
         end if
         ! radiation 
         if (nrad > 0) then
            gamma % e(irad,:) = zero
            gamma % e(irad, irad) = one / s_cspeed
         end if
         !

      end if
      !

      ! Maxwell equation
      if (.not.present(cc)) return
      !
      if (nmag > 0) then
         rho = one / bspeed**2 !bspeed = one / bspeed**2
         h0 = one / cspeed**2 !espeed = one / espeed**2
         !if(is == 0) h0 = zero
         !  gamma%e(ibb:ife,:) = zero
         !  gamma%e(:,ibb:ife) = zero
         call put_diag(gamma % e, rho, ibb, ibe)
         call put_diag(gamma % e, h0, ieb, iee)
         call put_diag(gamma % e, h0, ifb, ifb)
         !  DivJ=0
         rho = rho / e_resistivity(qv, cc)
         call put_diag(gamma % e, rho, ife, ife)
      end if
   end function gamma
   !
   !***********************************************************************
   !  calculate flux A jacobian
   !***********************************************************************
   function ajacobian(vecn, qv, cc, r)
      implicit none
      type(cell) :: cc
      real(rfp), intent(in) :: vecn(:)
      type(vector) :: qv, rhox, h0x
      real(rfp) :: rho, h0, un, nv(3), r(ndim), unf
      type(matrix) :: ajacobian
      !  real(rfp)::rhopun
      integer :: i, im
      ajacobian % e = zero
      !
      if (neqf > 0) then
         rho = rhof(qv)
         h0 = h0f(qv)
         rhox = rhoxf(qv)
         h0x = hxf(qv)

         unf = sum(vecn * mesh_vel(r))
         un = dot_product(vecn, qv % v(imb:ime)) - unf
         rhox = rhox * un
         h0x = h0x * un
         !  rhopun = rhox%v(ico) * un
         ! continuity row(ico)
         if (nceq > 0) then
            ajacobian % e(ico,:) = rhox % v(:)
            ajacobian % e(ico, imb:ime) = rho * vecn(:)
         end if
         ! momentum rows(imb:ime)
         if (nmeq > 0) then
            i = 0
            do im = imb, ime
               i = i + 1
               ajacobian % e(im, ico) = rhox % v(ico) * qv % v(im) + vecn(i)
               ajacobian % e(im, imb:ime) = rho * qv % v(im) * vecn(:)
               ajacobian % e(im, im) = ajacobian % e(im, im) + rho * un
               ajacobian % e(im, ien:) = rhox % v(ien:) * qv % v(im)
            end do
         end if
         !energy row(ien)
         if (neeq > 0) then
            ajacobian % e(ien,:) = rhox % v * h0 + h0x % v * rho
            ajacobian % e(ien, ico) = ajacobian % e(ien, ico) + unf

            if (nmeq > 0) ajacobian % e(ien, imb:ime) = ajacobian % e(ien, imb:ime) + &
               rho * h0 * vecn(:)
         end if
         !
         !turbulence and multi-spice rows or auxiliary
         if(neqf > ien) then
            do i = ien + 1, neqf
               ajacobian % e(i,:) = rhox % v * qv % v(i)
               ajacobian % e(i, imb:ime) = rho * qv % v(i) * vecn(:)
               ajacobian % e(i, i) = ajacobian % e(i, i) + rho * un
            end do
         end if
         if (nrad > 0) ajacobian % e(irad,:) = zero ! radiation

      end if
      ! Maxwell
      if (nmag > 0) then
         nv = avto3(vecn)
         ajacobian % e(ibb:ife,:) = zero
         ajacobian % e(:, ibb:ife) = zero
         ajacobian % e(ibb:ibe, ieb:iee) = cross2mat(nv)
         ajacobian % e(ieb:iee, ibb:ibe) = -cross2mat(nv)
         ajacobian % e(ibb:ibe, ifb) = nv
         ! no divE=rho/epsi
         !  ajacobian%e(ieb:iee,ife    )     =  nv    !!!!!!!----
         ajacobian % e(ifb, ibb:ibe) = nv
         ajacobian % e(ife, ieb:iee) = nv !!!/ e_resistivity(qv,cc)  !!!!!
      end if
   end function ajacobian
   !
   !***********************************************************************
   !  calculate eigenvalues of A matrix
   !
   !  This subroution only used for hpp=hp,htp=ht,rhotp=rhot
   !
   !***********************************************************************
   subroutine eigenvalues(is, c1, c2, vecn, qv, rhopp, r)
      implicit none
      real(kind = rfp), intent(out) :: c1, c2
      type(vector) :: qv, rhox, h0x
      real(kind = rfp) :: rho, un, vecn(:), rhop, rhot, ht
      real(kind = rfp) :: rhopp, hp, r(ndim), unf
      real(kind = rfp) :: dpp, rhopmpp, ub
      integer, intent(in) :: is
      !
      if (neqf == 0) then
         c1 = one; c2 = one
         return
      end if
      unf = dot_product(vecn, mesh_vel(r))
      un = dot_product(vecn, qv % v(imb:ime)) - unf

      !no energy equation
      rhox = rhoxf(qv)
      h0x = hxf(qv)
      rho = rhof(qv)

      rhop = rhox % v(ico)
      rhot = rhox % v(ien)
      ht = h0x % v(ien)
      hp = h0x % v(ico)

      if(neeq == 0) then
         ht = htf(qv)
      end if

      if (is == 0) rhopp = rhop
      !
      if(rhot==0 .and. ht==0) then
         print*,'eigen value is infinity'
         stop      
      else
         dpp = rhopp + rhot * (one - rho * hp) / (rho * ht)
      end if

      if (abs(dpp) <= mytiny) then
         print *, rhopp, rhot, rho, ht, qv % v
         print *, 'error eigenvalue subroutine'
         stop
      end if
      dpp = one / dpp
      rhopmpp = (rhop - rhopp) * dpp
      ub = half * rhopmpp * un
      dpp = sqrt(abs(ub * ub + dpp))
      !  U + (rhop-rhopp)*U+ sqrt(rhop-rhopp)^2U^2 + 4N^2(rhopp+rhot'))
      !       ---------------------------------------------------------
      !                          2*(rhopp + rhot' * (1-rho*hp)
      !   rhot' = rhot / (rho*ht),  N^2 = nx^2 + ny^2 + nz^2 = 1
      !   if(dpp < 1.e-15_rfp) print *,' entropy fix'
      c1 = un + ub + dpp
      c2 = un + ub - dpp
   end subroutine eigenvalues

   subroutine eigenvector(c1, c2, zm, zminv, vecn, qv, rhopp, r)
      implicit none
      real(kind = rfp), intent(in) :: c1, c2, rhopp
      type(matrix) :: zm, zminv
      type(vector) :: qv, h0x, rhox
      real(kind = rfp), intent(in) :: vecn(:)
      real(kind = rfp) :: htinv, c1mun, c2mun, c2mc1inv, rho, un, ht, hp, d, rhop, rhot
      real(kind = rfp) :: rhoc2mc1inv, c1munoc2mc1, c2munoc2mc1, dnminv
      real(kind = rfp) :: mv(size(vecn)), nv(size(vecn)), sumvecn, unf, r(ndim)
      integer :: k

      rhox = rhoxf(qv)
      h0x = hxf(qv)
      rho = rhof(qv)
      !  rhop = rhox%v(ico)
      !  rhot = rhox%v(ien)
      ht = h0x % v(ien)
      hp = h0x % v(ico)

      ! no energy equation
      if(neeq == 0) then
         ht = htf(qv)
      end if
      !  d = one + (rhop - rhopp) / (rhopp + rhot * ( one - rho * hp)/(rho * ht))
      !  un = dot_product(vecn,qv%v(imb:ime)) * d
      unf = sum(vecn * mesh_vel(r))
      un = dot_product(vecn, qv % v(imb:ime)) - unf

      !
      if(ht==0) then
         print*,'eigen vector has infinity'
         stop
      end if

      htinv = (one - rho * hp) / ht
      c1mun = c1 - un
      c2mun = c2 - un
      !  c2mc1inv  = one / (c2 - c1 + mytiny)
      c2mc1inv = one / (c2 - c1)
      rhoc2mc1inv = c2mc1inv / rho
      c1munoc2mc1 = c1mun * c2mc1inv
      c2munoc2mc1 = c2mun * c2mc1inv
      !
      zm % e = zero ! initial eigenvector M    to zero
      zminv % e = zero ! initial eigenvector M^-1 to zero  
      select case(ndim)
      case(1) ! 1-dimension
         !  M  Matrix
         !  ico row
         zm % e(ico, imb) = c1mun * rho
         zm % e(ico, ien) = c2mun * rho
         !  imb row
         zm % e(imb, imb) = vecn(1)
         zm % e(imb, ien) = vecn(1)
         !  ien row
         zm % e(ien, ico) = one
         zm % e(ien, imb) = c1mun * htinv
         zm % e(ien, ien) = c2mun * htinv
         !  M^-1 Matrix
         ! ico row
         zminv % e(ico, ico) = -htinv / rho
         zminv % e(ico, ien) = one
         ! imb row
         zminv % e(imb, ico) = -rhoc2mc1inv
         zminv % e(imb, imb) = c2munoc2mc1 * vecn(1)
         ! ien row  
         zminv % e(ien, ico) = rhoc2mc1inv
         zminv % e(ien, imb) = -c1munoc2mc1 * vecn(1)
      case(2) ! 2-dimension
         mv(1) = vecn(2)
         mv(2) = -vecn(1)
         !  M  Matrix
         !  ico row
         zm % e(ico, ime) = c1mun * rho
         zm % e(ico, ien) = c2mun * rho
         !  imb:ime row
         zm % e(imb:ime, imb) = mv
         zm % e(imb:ime, ime) = vecn
         zm % e(imb:ime, ien) = vecn
         !  ien row
         zm % e(ien, ico) = one
         zm % e(ien, ime) = c1mun * htinv
         zm % e(ien, ien) = c2mun * htinv
         !  M^-1 Matrix
         ! ico row
         zminv % e(ico, ico) = -htinv / rho
         zminv % e(ico, ien) = one
         ! imb row
         zminv % e(imb, imb:ime) = mv
         ! ime row
         zminv % e(ime, ico) = -rhoc2mc1inv
         zminv % e(ime, imb:ime) = c2munoc2mc1 * vecn
         ! ien row  
         zminv % e(ien, ico) = rhoc2mc1inv
         zminv % e(ien, imb:ime) = -c1munoc2mc1 * vecn
      case(3)
         dnminv = one / sqrt(two * (one - vecn(1) * vecn(2) &
            -vecn(1) * vecn(3) &
            -vecn(2) * vecn(3)))
         mv(1) = (vecn(2) - vecn(3)) * dnminv
         mv(2) = (vecn(3) - vecn(1)) * dnminv
         mv(3) = (vecn(1) - vecn(2)) * dnminv
         !
         sumvecn = sum(vecn)
         nv(:) = (vecn(:) * sumvecn - one) * dnminv
         !  M  Matrix
         !  ico row
         zm % e(ico, ime) = c1mun * rho
         zm % e(ico, ien) = c2mun * rho
         !  imb:ime row
         zm % e(imb:ime, imb) = mv
         zm % e(imb:ime, imb + 1) = nv
         zm % e(imb:ime, ime) = vecn
         zm % e(imb:ime, ien) = vecn
         !  ien row
         zm % e(ien, ico) = one
         zm % e(ien, ime) = c1mun * htinv
         zm % e(ien, ien) = c2mun * htinv
         !  M^-1 Matrix
         ! ico row
         zminv % e(ico, ico) = -htinv / rho
         zminv % e(ico, ien) = one
         ! imb row
         zminv % e(imb, imb:ime) = mv(:)
         ! imb + 1 row
         zminv % e(imb + 1, imb:ime) = nv(:)
         ! ime row
         zminv % e(ime, ico) = -rhoc2mc1inv
         zminv % e(ime, imb:ime) = c2munoc2mc1 * vecn(:)
         ! ien row  
         zminv % e(ien, ico) = rhoc2mc1inv
         zminv % e(ien, imb:ime) = -c1munoc2mc1 * vecn(:)
      end select
      !
      if(neqf > ien) then
         do k = ien + 1, neqf
            zm % e(k, k) = one
            zminv % e(k, k) = one
         end do
      end if

   end subroutine eigenvector
   !
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !  End of Jacobian Section
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   !***********************************************************************
   !  flux terms E.N
   !***********************************************************************
   function fluxv(vecn, qv, cc, r)
      implicit none
      type(cell) :: cc
      real(rfp), intent(in) :: vecn(:)
      type(vector) :: fluxv, qv
      real(rfp) :: rho, p, h0, vel(ndim), un, rhoun, v(3), r(ndim), unf
      !
      if (neqf > 0) then
         rho = rhof(qv)
         p = qv % v(ico)
         h0 = h0f(qv)
         vel = qv % v(imb:ime)
         unf = sum(vecn * mesh_vel(r))
         un = dot_product(vecn, vel) - unf
         rhoun = rho * un
         fluxv % v(ico) = rhoun
         fluxv % v(imb:ime) = rhoun * vel(:) + p * vecn(:)
         fluxv % v(ien) = rhoun * h0 + unf * pf(qv)
         !fluxv % v(ien + 1:neqf) = rhoun * qv % v(ien + 1:neqf)
         if (neqf > ien) fluxv % v(ien + 1:neqf) = rhoun * qv % v(ien + 1:neqf)
         if (nrad > 0) fluxv % v(irad) = zero !radiation
      end if
      if (nmag > 0) then
         v = avto3(vecn)
         fluxv % v(ibb:ibe) = acrossb(v, qv % v(ieb:iee)) + v * qv % v(ifb)
         ! no divE=rho/epsi
         !  fluxv%v(ieb:iee) = -acrossb(v,qv%v(ibb:ibe))   + v * qv%v(ife)   !!!
         fluxv % v(ieb:iee) = -acrossb(v, qv % v(ibb:ibe))
         fluxv % v(ifb) = sum(v * qv % v(ibb:ibe))
         fluxv % v(ife) = sum(v * qv % v(ieb:iee)) !!/ e_resistivity(qv,cc)  !!!
      end if
      !
   end function fluxv
   !
   !***********************************************************************
   !  Conservation variables
   !***********************************************************************
   function conservative_qv(qv)result(q)
      implicit none
      type(vector) :: q, qv
      real(rfp) :: rho, p, h0, vel(ndim), un, rhoun
      ! 
      if (neqf > 0) then
         rho = rhof(qv)
         p = zero
         if (nceq > 0) q % v(ico) = rho
         if (nmeq > 0) then
            p = pf(qv)
            vel = qv % v(imb:ime)
            q % v(imb:ime) = rho * vel(:)
         end if
         if (neeq > 0) then
            h0 = h0f(qv)
            q % v(ien) = rho * h0 - p
         end if
         if (neqf > ien) q % v(ien + 1:neqf) = rho * qv % v(ien + 1:neqf)
         if (nrad > 0) q % v(irad) = qv % v(irad)
      end if
      if (neqm > 0) then
         q % v(ibb:ife) = qv % v(ibb:ife)
         q % v(ieb:iee) = qv % v(ieb:iee) / s_cspeed**2
         !   q%v(ieb:iee) = zero
         q % v(ifb) = qv % v(ifb) / s_cspeed**2
      end if
   end function conservative_qv
   !
   subroutine meigenmatrix(bspeed, espeed, v, zm)
      implicit none
      real(rfp) :: zm(:,:)
      real(rfp) :: dnminv, a, b, bspeed, espeed
      real(kind = rfp) :: v(:), vecn(3), mv(3), nv(3), sumvecn, w(3, 3)
      integer :: k
      !
      a = espeed / bspeed
      b = one / a !s_fspeed / s_bspeed
      zm = zero
      call put_diag(zm, a, ibb, ibe)
      call put_diag(zm, b, ieb, iee)
      call put_diag(zm, b, ifb, ifb)
      call put_diag(zm, a, ife, ife)
      return
      !

      vecn = avto3(v)
      dnminv = one / sqrt(two * (one - vecn(1) * vecn(2) &
         -vecn(1) * vecn(3) &
         -vecn(2) * vecn(3)))
      mv(1) = (vecn(2) - vecn(3)) * dnminv
      mv(2) = (vecn(3) - vecn(1)) * dnminv
      mv(3) = (vecn(1) - vecn(2)) * dnminv
      !
      sumvecn = sum(vecn)
      nv(:) = (vecn(:) * sumvecn - one) * dnminv
      zm(ibb:ife, ibb:ife) = zero
      w = tensor(mv, mv) + tensor(nv, nv) !  + tensor(vecn,vecn)
      zm(ibb:ibe, ibb:ibe) = (w + tensor(vecn, vecn)) * a
      !   zm(ibb:ibe,ibb:ibe) = w * a
      zm(ieb:iee, ieb:iee) = w * b
      zm(ifb, ifb) = b
      zm(ife, ife) = a
      !

      !
   end subroutine meigenmatrix

   subroutine aaws(bspeed, espeed, cc)
      implicit none
      real(kind = rfp) :: bspeed, espeed, fs, b(3), e(3), a, dl
      type(cell) :: cc
      a = g_mu0 * b_lenref / e_resistivity(cc % qv, cc) ! mu0*sigma*L
      a = g_mu0 * b_lenref * s_espeed / e_resistivity(cc % qv, cc)
      a = sqrt(a)
      ! a = max(one,a)
      !  a = max(1/s_cspeed,a)
      !  a = one
      !  a = s_espeed / s_dt_inv
      bspeed = a
      espeed = s_espeed / a
      !  bspeed = one   !a
      !  espeed = g_mu0 * b_lenref / e_resistivity(cc%qv,cc)    !s_cspeed / a
      !  espeed = min(one / espeed,s_espeed)
      !  espeed = cc%jv(iee) * 1.e3_rfp + one   !min(s_cspeed, two * 1.e7_rfp)

      !  bspeed = one
      !  espeed = one / a
      !  bs = s_bspeed; fs = s_fspeed; es = s_espeed;
   end subroutine aaws
   !
end module gems_jacob
