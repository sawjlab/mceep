	subroutine musc(m2,p,rad_len,dth,dph)
!--------------------------------------------------------------------
! MUSC - Simulate multiple scattering of any particle.
!
! ASSUMPTIONS: DTH and DPH given in milli-radians, RAD_LEN in radiation
!   lengths. The formula used is due to Rossi and Greisen (See the book
!   by Segre, NUCLEI AND PARTICLES, 1982, p. 48.) The formula assumes a
!   gaussian distribution for the scattering angle. It is further assumed
!   that the angles DTH and DPH are located at an angle with respect 
!   to the beam direction that is large enough so that DTH-DPH space can
!   be approximated as a cartesian space.
!
! D. Potterveld - Sept. 1985
!
! Add option for protons 
!   Precision not good enough for our thick targets. Latest values supplied:
!   Lynch and Dahl, NIM B58 (1991) 6.
!   Note, that there is a typo in Particle data booklet: eps = 0.088 instead
!   of 0.2! Precision: -5% deviation for H2, +8 for 238U at radiation 
!   lengths between ~ 0.01 and 1.
!
! H.J. Bulten - Aug. 1991
!--------------------------------------------------------------------

	implicit none

	DOUBLE PRECISION		Es, epsilon
	parameter (Es = 13.6)		!MeV
	parameter	(epsilon = 0.088)

	DOUBLE PRECISION		rad_len, dth, dph
	DOUBLE PRECISION		beta, g1, g2, theta_sigma
	DOUBLE PRECISION		m2, p

! Compute scattering angles, THETA_SCAT from a gaussian distribution,
! PHI_SCAT from uniform distribution.

	beta = p / sqrt(m2+p*p)
	theta_sigma= Es/p/beta*sqrt(rad_len)*(1+epsilon*
     >   log10(rad_len))

! Compute new trajectory angles (units are rad)

	call gaussians(g1,g2)
	dth = dth + theta_sigma * g1
	dph = dph + theta_sigma * g2

	return
	end


	subroutine gaussians(g1,g2)

!--------------------------------------------------------------------
! This subroutine generates a random number distributed about a Gaussian
! centered at zero with a standrd deviation of 1.
! Use as a function like: VAL' = VAL + WIDTH*X
!  
! Algorithm is from PHYSICS LETTERS B V.204, P.83 ''Review of particle
! properties", Particle Data Group.
!
!        Changes
!           Uses RANECU uniform random number generator,
!                         rather than VAX generator.  -WLH
!	
!--------------------------------------------------------------------

	implicit none

	INCLUDE 'ran.cmn'

	DOUBLE PRECISION u1,u2,v1,v2,s,g1,g2
        REAL             iran_vec(2)

1	CALL RANECU(IRAN_VEC,2)
	v1 = 2.d0*DBLE(IRAN_VEC(1))-1.d0
	v2 = 2.d0*DBLE(IRAN_VEC(2))-1.d0
	s  = v1**2+v2**2

	if (s.gt.1.d0 .or. s.eq.0.d0) goto 1

	g1 = v1*sqrt(-2.d0*log(s)/s)
	g2 = v2*sqrt(-2.d0*log(s)/s)

	return                                                           
	end                                                              
