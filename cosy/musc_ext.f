	subroutine musc_ext(m2,p,rad_len,x_len,dth,dph,y,x)
!-------------------------------------------------------------------------
! MUSC - Simulate multiple scattering of any particle.
!        Used for extended scatterers.
!
! According to Particle Data Booklet, July 1994
!
!-------------------------------------------------------------------------

	implicit none

	real*8		Es, epsilon
	parameter	(Es = 13.6e-3) ! MeV
	parameter	(epsilon = 0.088)

	DOUBLE PRECISION rad_len, x_len, dth, dph, x, y
	DOUBLE PRECISION beta, g1, g2, theta_sigma
	DOUBLE PRECISION m2, p

	if (rad_len.eq.0.d0) return
	if (x_len.le.0.d0 .or. rad_len.lt.0.d0) then
	  write(6,*) 'x_len or rad_len < 0 in musc_ext.f'
	  write(6,*) 
     >    'This is bad.  Really bad.  Dont even ask how bad it is.'
	  write(6,*) 'Just fix it now.'
	  stop
	endif
	if (p.lt.50.d0) write(6,*) 
     >       'Momentum passed to musc.f should be in MeV'

	beta = p / sqrt(m2+p*p)
	theta_sigma=Es/p/beta*sqrt(rad_len)*(1+epsilon*log10(rad_len))

! Compute new trajectory angles and displacements (units are rad and cm)

	call gaussians(g1,g2)
	dth = dth + theta_sigma*g1
	x   = x   + theta_sigma*x_len*g2/sqrt(12.d0) + 
     >        theta_sigma*x_len*g1/2.d0

	call gaussians(g1,g2)
	dph = dph + theta_sigma*g1
	y   = y   + theta_sigma*x_len*g2/sqrt(12.d0) + 
     >        theta_sigma*x_len*g1/2.d0
	
	return
	end









