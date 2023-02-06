

module constants
 !defines constants to be used across program
 use precisn
 implicit none
 public hbarc, alpha, mred, pi, m1, m2, z1, z2
 real(wp), parameter :: hbarc = 197.327053 !MeV fm
 real(wp), parameter :: alpha = 1./137.035999679 
 real(wp), parameter :: m1 = 2. * 931.49406 !MeV
 real(wp), parameter :: m2 = 3. * 931.49406 !MeV
 real(wp), parameter :: mred = m1 * m2 / (m1 + m2) !MeV
 real(wp), parameter :: pi = acos(-1.)
 integer, parameter  :: z1 = 1
 integer, parameter  :: z2 = 1
end module constants

program mainprog

 use precisn
 use constants
 implicit none
 integer :: N, N_start, n_max, n_step,Nth, nume, i, j, B, lmax,l0, p, vswap, mm, plp, kk, pl,c_len
 real(wp) :: a, x1, x2, e, k, eta, ro, n_d, V, gamma, re_i, psiabs, h, plgndr, VijOR
 real(wp), allocatable, dimension(:, :) :: M, plplot, R
 real(wp), allocatable, dimension(:) :: x, w, Vij, VijDi
 real(wp), allocatable, dimension(:) ::  V_plot, V_l_lp
 real(wp), allocatable, dimension(:, :) :: C
 complex(wp), allocatable, dimension(:) :: in, out, inp, outp, ft
 complex(wp), allocatable, dimension(:, :) :: psiout, ulint, Z, S, ulint_temp, ulext
 character*1 strj, strl0
 character(len=:),allocatable::strn,stra,strk
 character*3 strth
 character*2 strn_d
 character*1 strnd



open (unit=1, file='input_Cexp_nd6_lmax6.dat')
 read (unit=1, fmt=10); 10 format(a) 
 read (unit=1, fmt=11); 11 format(a)
 read (unit=1, fmt=12) lmax; 12 format(i2) !orbital momentum
 read (unit=1, fmt=13); 13 format(a)
 read (unit=1, fmt=203) l0; 203 format(i2) !orbital momentum
 read (unit=1, fmt=204); 204 format(a)
 read (unit=1, fmt=14) a; 14 format(f6.2)!matching radius
 read (unit=1, fmt=15); 15 format(a)
 read (unit=1, fmt=16) n; 16 format(i3) !mesh pts
 read (unit=1, fmt=17); 17 format(a)
 read (unit=1, fmt=111) nth; 111 format(i3) !THETA mesh pts
 read (unit=1, fmt=112); 112 format(a)
 read (unit=1, fmt=200) e; 200 format(f5.2) !max energy
 read (unit=1, fmt=201); 201 format(a)
 read (unit=1, fmt=208) n_d; 208 format(f5.2) !n_d value
 read (unit=1, fmt=209); 209 format(a)
 read( unit=1, fmt=210) vswap; 210 format(i1) !selects a potential (0,couloumb)(1,vl)


do lmax = 28, 40, 2
a = 80

n = 250

 if (n.ge.100) then
	allocate(character(len=3) :: strn)
	write(strn,'(i3)'),n
 else
	allocate(character(len=2) :: strn)
	write(strn,'(i2)'),n
 end if

 if (a.ge.100) then
	allocate(character(len=3) :: stra)
	write(stra,'(i3)'), ceiling(a)
 else
	allocate(character(len=2) :: stra)
	write(stra,'(i2)'), ceiling(a)
 end if


 if (lmax.ge.10) then
	allocate(character(len=2) :: strk)
	write(strk,'(i2)'), lmax
 else
	allocate(character(len=1) :: strk)
	write(strk,'(i1)'), lmax
 end if

print *, strk


 p = lmax + 1;x1 = 0.; x2 = 1.
 c_len = p*N

 allocate(C(c_len,c_len), M(n, n), x(n), w(n))
 allocate(S(lmax+1,lmax+1))
 allocate(R(p,p), Z(p,p))

 k = sqrt(2. * mred * E) / hbarc; ro = k * a
 eta = alpha * z1 * z2 * mred / (hbarc * k)

call gauleg(N, x1, x2, x, w)
call cmatrix(N,Nth, lmax, a, x, w, C, n_d, vswap, E, c_len)


write(strl0,'(i1)') l0
write(strnd,'(i1)') ceiling(n_d)



open(unit=5, file='cmatrix_lmax'//strk//'_l0'//strl0//'_nd'//strnd//'_N'//strN//'_a'//stra//'.dat')
 write (unit=5, fmt=501) 'lmax, l0, a, n, nth, e_min, e_max, delta, n_d, vswap, C^(-1)ij'; 501 format(a) !brief desc.
 write (unit=5, fmt=502) lmax; 502 format(i2) !orbital momentum
 write (unit=5, fmt=503) l0; 503 format(i2) !orbital momentum
 write (unit=5, fmt=504) a; 504 format(f6.2)!matching radius
 write (unit=5, fmt=505) n; 505 format(i3) !mesh pts
 write (unit=5, fmt=506) nth; 506 format(i3) !THETA mesh pts
 write (unit=5, fmt=507) e; 507 format(f5.3) !min energy
 write (unit=5, fmt=510) n_d; 510 format(f5.2) !n_d value
 write (unit=5, fmt=511) vswap; 511 format(i1) !selects a potential (0,couloumb)(1,vl)


do i = 1,c_len
	do j = 1, c_len
		write(unit=5, fmt=500) C(i,j)
		500 format(ES23.16)
	end do
end do
close(unit=5)


 deallocate(C, M, x, w)
 deallocate(S)
 deallocate(R, Z, strn, stra,strk)


end do





end program mainprog





real(wp) function delta(i, j)
use precisn
implicit none
integer, intent(in) :: i, j
if(i == j) then
	delta = 1.
else 
	delta = 0.
end if
end function delta




subroutine cmatrix(N,Nth, lmax, a, x, w, C, n_d, vswap, E, c_len)
 !calculates the C-matrix
 use precisn
 use omp_lib
 implicit none
 integer, intent(in) :: N,Nth, lmax, vswap, c_len
 real(wp), intent(in) :: a,  n_d, E
 real(wp), dimension(1:N), intent(in) :: x, w
 real(wp), dimension(1:c_len,1:c_len), intent(out) :: C
 integer :: i, j,l,lp
 real(wp) :: Cij


 !$OMP PARALLEL DO PRIVATE(l,lp,i,j)
 do l = 0, lmax
	do lp = 0,lmax
 		do i = 1, N
 			do j = 1, N
				C(l*N + i,lp*N + j) = Cij(N,Nth, lmax, l, lp, a, i, j, w, x, x(i), x(j), n_d, vswap, E)
			end do
		 end do
	end do
 end do
 !$OMP END PARALLEL DO

 call inversg(C,c_len)

end subroutine cmatrix





real(wp) function Cij(N,Nth, lmax, l, lp, a, i, j, w, x, xi, xj, n_d, vswap, E)
 !calcultes the value of the C-matrix at i, j
 !C1 is the l = 0 contribution
 !C2 is the l > 0 contribution
 use precisn
 use constants
 implicit none
 integer, intent(in) :: N,Nth, lmax, l, lp, i, j, vswap
 real(wp), intent(in) :: a, xi, xj, n_d, E
 real(wp), dimension(1:N), intent(in) :: x, w
 real(wp), dimension(1:N) :: v_l_lp
 real(wp) :: Vij, V, C1, C2, C3




 if (i == j) then !diagonal points
    call Vij_matrixe(N,Nth, lmax, l, lp, x, w, n_d, vswap, a, V_l_lp)
    C1 = (4*N*(N+1)+3+(1-6*xi)/(xi*(1-xi)))/(6*a*a*xi*(1-xi))
    Vij = v_l_lp(i)
    C2 = 1.*l*(l+1)/(2.*a*a*xi*xi)
    C3 = -E
   
 else !off-diagonal points
    C1 = ((-1)**(i+j))*(N*(N+1)+1+(xi+xj-2*xi*xj)/((xi-xj)*(xi-xj))&
         -1/(1-xi)-1/(1-xj))/(2*a*a*sqrt(xi*xj*(1-xi)*(1-xj)))
    Vij = 0.
    C2 = 0.
    C3 = 0
 end if

 if (l.ne.lp) then
 C1 = 0
 C2 = 0 
 C3 = 0
 end if


 Cij = (hbarc * hbarc / mred) * (C1 + C2) + C3 + Vij
end function Cij





subroutine Din(rmax, rn, r_e, theta)
 !gives the limit of Din for a particular theta
 use precisn
 use constants
 implicit none
 real(wp), intent(in) :: rn, r_e, theta
 real(wp), intent(out) :: rmax
 real(wp) :: cut, thetar

 cut = atan(rn/r_e)
 thetar = theta

 !need theta in	range -pi, pi
 if (thetar.lt.-pi.or.thetar.gt.pi) then
    print *, 'theta is out of range'
    rmax = -1
 else if (thetar.lt.0) then 	 !fn is even but is written for	positive input 
    thetar = -thetar
end if



 !top hemisphere
 if (thetar.ge.0.and.thetar.lt.cut) then
    rmax = r_e * cos(thetar) + sqrt(r_e**2 * (cos(thetar)**2 - 1) + rn**2)
 !straight edge
 else if (thetar.ge.cut.and.thetar.lt.(pi - cut)) then
    rmax = rn / sin(thetar)
 !bottom hemisphere
 else if (thetar.ge.(pi-cut).and.thetar.le.pi) then
    rmax = -r_e * cos(thetar) + sqrt(r_e**2 * (cos(thetar)**2 - 1) + rn**2)
 end if

end subroutine Din






real(wp) function V(x,thetar,l,n_d,vswap)
 !calculates the potential at point x
 use precisn
 use constants
 implicit none
 real(wp) :: x, VN, VC, rn, U_0, eps_0, r_e, V_l
 real(wp) :: costh, p_l, plgndr, thetar, n_d, V_0, gamma
 real(wp) :: rmax
 integer :: l, i, j, vswap


 costh = cos(thetar)
 U_0 = 35. !MeV depth of nuclear well
 rn = 3.89 !fm width of nuclear well
 eps_0 = 0.05526349406 ! e^2 GeV^-1 fm^-1 permittivity of free space
 r_e = n_d * rn
 V_0 = 0.37
 V = 0
 
 !calculates the edge of the well
 call Din(rmax, rn, r_e, thetar)


 select case (vswap)

 case (0) !COULOMB AND SQUARE WELL POTENTIAL DT

 if (x.lt.rn) then
    V = -U_0
 else
    V = hbarc * alpha / (4 * pi * eps_0 * x)
 end if

 case (1) !V_EFF POTENTIAL

 

 if (x.lt.rmax) then
    V = -U_0 !eq9
 else if (mod(l,2).eq.0) then
    do j = 0, l, 2
       call V_l_out(V_l,j,x,rn,r_e,n_d)!eq10,11
       p_l = plgndr(j,0,costh)
       V = V + (V_0 / n_d) * p_l * V_l !eq9
    end do
	!V = 0 !no odd l V_eff terms 
 else 
	V = 0. !if odd l v = 0 outside of the well
 end if


 case (2) !JUST COULOMB

 V = hbarc * alpha / (4 * pi * eps_0 * x)


 end select

end function V





subroutine V_l_out(V_l,l,x,rn,r_e,n_d)
 !Calculates the V_l term of the potential
 use precisn
 use constants
 implicit none
 real(wp), intent(in) :: x, rn, r_e,  n_d
 real(wp), intent(out) :: V_l
 integer, intent(in) :: l
 real(wp) :: re_l,V_l_temp, sum, a, b, c, d, gamma
 integer :: i


 re_l = real(l)
 a = (re_l + 1.) / 2.
 b = re_l/2. + 1.
 V_l_temp = 0 


 if ((x / rn).le.n_d) then
    do i = 1, l/2
	c = i + 0.5
         d = real(i)
         V_l_temp = V_l_temp + ( -(x / r_e)**(2*i-l-2) + (x / r_e)**(-2*i+l) )&
           *sqrt(1 - (x/r_e)**2) * ( (gamma(a) * gamma(d)) / ( gamma(b) * gamma(c) ))
    end do
    V_l =   (1/pi) * V_l_temp + (1/pi)*( (x/r_e)**(-l-1) )*( 2*gamma(a) / (sqrt(pi) * gamma(b)) )*asin(x/r_e)&
            -(1/pi)*( (x/r_e)**l )*( 2*gamma(a) / (sqrt(pi) * gamma(b)) )*log(tan(asin(x/r_e)/2))
 else !x/rn gt n_d
    V_l = ((x / r_e)**(-l-1)) * ( gamma(a) / (sqrt(pi) * gamma(b)) ) !eq11
 end if

end subroutine V_l_out





subroutine yl0(y, l, costh)
!calculates the spherical harmonics for m = 0
 use precisn
 use constants
 implicit none
 integer, intent(in) :: l
 real(wp), intent(in) :: costh 
 real(wp) :: plgndr, p_l
 real(wp), intent(out) :: y

 p_l = plgndr(l,0,costh)

 y = sqrt((2*l + 1)/(4*pi)) * p_l


end subroutine yl0







subroutine Vij_matrixe(N, Nth, lmax, l, lp, x, w, n_d, vswap, a, Vij)
use precisn
use constants
implicit none
integer, intent(in):: l, lp, lmax, vswap, N, Nth
integer :: j, i, k
real(wp) :: y, yl, ylp, Vv, V, theta
real(wp), intent(in):: n_d, a
real(wp), dimension(1:Nth) :: ctheta0, w_ctheta
real(wp), dimension(1:N), intent(in) :: x, w
real(wp), dimension(1:N), intent(out) :: Vij

call gauleg(Nth, -1.d0, 1.d0, ctheta0, w_ctheta)


Vij = 0
do j = 1, N
Vv = 0
	do i = 1, Nth
		y = acos(ctheta0(i))
		call yl0(ylp, lp, ctheta0(i))
		call yl0(yl, l, ctheta0(i))
	   	Vv = Vv + yl * w_ctheta(i) *  V(a * x(j), y, 6, n_d, vswap) * ylp * 2 * pi
	end do
Vij(j) = Vv
end do

!theta = pi / 2

!do j = 1,N
!Vv = 0
!	call yl0(ylp, lp, cos(theta))
!	call yl0(yl, l, cos(theta))
!   	Vv =  yl *  V(a * x(j), theta, lmax, n_d, vswap) * ylp * 2 * pi
!	print *,'V', V(a * x(j), theta, lmax, n_d, vswap)
!print *,'Vv', Vv
!end do


end subroutine Vij_matrixe





