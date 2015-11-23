C
C ------------------------------------------------------------------
C ------------------------------------------------------------------
C
C     Program   TOTERR:
C
C
C     NOTE:  You must do statistical averaging BEFORE
C            taking quadrature sums of errors here, as this sum
C            consists of positive definite terms and is therefore 
C            prone to statistics which won't get averaged out
C            properly.  So, before running this code, you need to
C            get the statistically averaged arrays first (by
C            executing the Fortran code, systerr, followed by the
C            PAW Kumac, systerr.kumac - see the README file). 
C
C     Calculate total error.  This is accomplished by reading
C     in vectors produced by running systerr.kumac.
C     The kumac scans the Ntuple produced by the Fortran code
C     systerr and produces vectors consisting of the average of 
C     the cross section weights vs. whatever kinematic variable.  
C     There are a total of ten vectors produced by systerr.kumac,
C     one for the unperturbed cross sections and nine for "standard"
C     variations of each of the nine kinematical quantities, in turn.
C     Note that the values of the kinematic quantities along the
C     independent axis are always the same - only the cross section
C     reflects the kinematic shifts.
C
C
C     INPUT (from err_file):
C
C            err(i,i):     These are the errors in units of 10^-3/1mr
C                            for each of the 9 variables in the order
C                            given in systerr.f (i.e. e1, ..., thetax)
C
C            i1,j1:        Index pair which implies that these variables
C                            are correlated.  This line and all lines
C                            below should not be present for an
C                            uncorrelated analysis.
C
C            i2,j2:        Next index pair
C              .
C              .
C              .
C            in,jn:        Last index pair
C
C
C      README (i.e. here's how the errors get added):
C
C      Suppose that the analysis involves shifts of variables 1, 4 and 7
C      in a correlated fashion, whereas all other variables are
C      uncorrelated.  (This is the case where, for instance, an elastic
C      measurement may force a correlation between the momenta of the
C      three particles.)  Then the fractional uncertainty in the 
C      cross section is:
C
C             del_sig = sqrt[(    dsig/dvar1*del1
C                              +  dsig/dvar4*del4
C                              +  dsig/dvar7*del7 ) ^2
C                              + (dsig/dvar2*del2)^2
C                              + (dsig/dvar3*del3)^2
C                              + (dsig/dvar5*del5)^2
C                              + (dsig/dvar6*del6)^2
C                              + (dsig/dvar8*del8)^2
C                              + (dsig/dvar9*del9)^2 ] / sig0
C
C      where dsig/dvari is the derivative of the cross section wrt
C      variable i, deli is the uncertainty of variable i and sig0
C      is the cross section evaluated at the nominal kinematics.
C      In this case, in addition to the nine uncertainties (deli)
C      given on a single line of the input file, the file should
C      contain 3 additional lines reading:
C                 1,4
C                 1,7
C                 4,7
C      If you expand the expression for del_sig you get
C      all cross terms involving those three indices (1,4,7)
C      That's why it's insufficient to have only two of the above
C      three index pairs in the input file!  The factor of two
C      arising from the cross terms is dealt with internally in the 
C      code by making the error matrix symmetric: 
C                  err(i,j) = err(j,i)
C
C ---------------------------------------------------------------------
C
C
      implicit none
c
      character*8      vec_file(0:9)
      character*80     err_file
      double precision err(9,9),vec_sig(0:9),diff(9),del_sig
      integer          i,j,k,lun
c
      data vec_file /'vec0.dat','vec1.dat','vec2.dat','vec3.dat',
     #               'vec4.dat','vec5.dat','vec6.dat','vec7.dat',
     #               'vec8.dat','vec9.dat'/
c
      write(6,*) ' Enter name of file containing error info > '
      read(5,*) err_file
c
c ---------------------------------------------------------------------
c     Produce the error matrix.
c ---------------------------------------------------------------------
c
      open(unit=1,file=err_file,status='old',form='formatted')
c
      do i=1,9                      ! initialize first
         do j=1,9
            err(i,j) = 0.d0
         enddo
      enddo
c
      read(1,*) (err(i,i),i=1,9)  ! diagonal elements (i.e. uncertainties)
      do i=1,100                  ! loop over correlated indices
         read(1,*,end=10) j,k     ! read in correlated index pairs
         err(j,k) = err(j,j)*err(k,k)
         err(k,j) = err(j,k)      ! symmetric matrix
      enddo
c
 10   close(unit=1)
c
      do i=1,9
         err(i,i) = err(i,i)**2
      enddo
c
c ---------------------------------------------------------------------
c     Open the vector files produced by systerr.kumac.
c     Also open the outfile which will contain the total errors.
c ---------------------------------------------------------------------
c
      do i=0,9                        ! open the vector files
         lun=i+10
         open(unit=lun,file=vec_file(i),status='old',
     #        form='formatted')
      enddo
      open(unit=20,name='toterr.dat',status='unknown',
     #     form='formatted')
c
c ---------------------------------------------------------------------
c     Loop over the elements of each vector, computing the total
c     (quadrature sum) error as we go.
c ---------------------------------------------------------------------
c
      do k=1,1000000 
         do i=0,9
            lun=i+10
            read(lun,*,end=200) vec_sig(i)
         enddo
         do i=1,9
            diff(i) = (vec_sig(i)-vec_sig(0))  ! dsig for 10^-3/1mr shifts
         enddo
         del_sig = 0.d0       ! initialize
         do i=1,9
            do j=1,9
               del_sig = del_sig + err(i,j)*diff(i)*diff(j)
            enddo
         enddo
         if(vec_sig(0) .ne. 0.d0) then
            del_sig = sqrt(del_sig)/vec_sig(0)  ! fractional error
         else
            del_sig = 0.d0
         endif
         write(20,190) del_sig
 190     format(2x,f12.5)
      enddo
c
 200  do k=0,9                               ! close the files
         lun=k+10
         close(unit=lun)
      enddo
      close(unit=20)
c
      stop
      end
