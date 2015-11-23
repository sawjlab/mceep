c  MAD spectrometer transfer functions for MAD in its 35 degree
c  minimum angle configuration. Formulated 4/11/02 -JJL
c 
c general description:
c functions take trajectories from plane to plane in the spectrometer.
c The planes are numbered 0 through 9
c plane 0: the target plane                                 
c plane 1: end of 1st drift                                  
c plane 2: entrance to magnet 1                              
c plane 3: middle of magnet 1                                
c plane 4: exit of magnet 1                                  
c plane 5: entrance to magnet 2                              
c plane 6: halfway into magnet 2 (good place for an x-collimator) 
c plane 7: middle of magnet 2                               
c plane 8: middle of magnet 2 + 0.5 m                        
c plane 9: middle of magnet 2 + 1.0 m                        
c plane 10: middle of magnet 2 + 1.5 m                       
c plane 11: middle of magnet 2 + 2.0 m also exit magnet 2   
c plane 12: end of significant magnetic field for magnet 2  
c plane 13: 1st drift chamber                               
c
c useful endplane information
c all functioons give the locations and directions of a trjectory relative to the central one
c x = x location in the local (magnet coordinate system) of the central trajectory 
c cx = x-direction cosine of the central trajectory (approx. theta)
c cy = z-direction cosine of the central trajectory
c cz = y-direction cosine of the central trajectory (approx. phi)
c l = central trajectory path length in mm
c ep	         x	       cx		cz		cy		l
c 1          0.000000	     0.000000	     1.000000	     0.000000	  1200.000000
c 2	   -70.027237	     0.086628	     0.996241	     0.000000	  1360.037598
c 3	    42.759571       -0.005988	     0.999982	     0.000002	  3364.012207
c 4	  -120.611023	    -0.096570	     0.995326	     0.000006	  5670.668457
c 5	  -120.390579	     0.144034	     0.989573	     0.000006	  5874.040527
c 6	   -26.365952	     0.118558	     0.992947	     0.000004	  6580.346680
c 7	    63.262596	     0.057633	     0.998338	     0.000000	  7584.516113
c 8	    84.443748	     0.027914	     0.999610	    -0.000002	  8084.983398
c 9	    92.592392	     0.006456	     0.999979	    -0.000004	  8585.059570
c 10	    92.738144	    -0.004253	     0.999991	    -0.000006	  9085.061523
c 11	    89.496307	    -0.008001	     0.999968	    -0.000006	  9585.072266
c 12	    85.170967	    -0.009050        0.999959	    -0.000006	 10085.08984
c 13	   401.250183	     0.078137	     0.996943	    -0.000006	 11345.88086


c
c for various pairs of planes there are 5 functions:
c  x_m35_i_j(x,m): gives x position at plane j as a function of trajectory
c              parameters at plane i.
c  t_m35_i_j(x,m): gives theta at plane j
c  y_m35_i_j(x,m): gives y at plane j
c  p_m35_i_j(x,m): gives phi at plane j
c  l_m35_i_j(x,m): gives difference in length from the "central trajectory
c              between planes i and j
c
c input: m=5
c        x is a 5 element array (REAL)
c        x(1)= x position of trajectory at plane i
c        x(2)= theta of trajectory at plane i
c        x(3)= y position of trajectory at plane i
c        x(4)= phi of trajectory at plane i
c        x(1)= delta for the trajectory 
c             delta=(p-p0)/p0 as usual 
c                   p= momentum of trajectory
c                   p0= momentum of the central trajectory
c Units: so called "natural units" are used. (Meters, tangent angle,
c        fractional deltas). DON'T BE FOOLED BY THE mm's ABOVE!
c 
c Coordinates are always relative to the central trajectory
c  for the central trajectory start at plane 0 with
c       x0 = 0., 0., 0., 0., 0. 
c  at plane 2
c      x2(1) = x_m02(x0,5)
c      x2(2) = t_m02(x0,5)
c      x2(3) = y_m02(x0,5)
c      x2(4) = p_m02(x0,5)
c      x2(5) = x0(5)
c then subsequently at plane 3
c      x3(1) = x_m23(x2,5)
c      x3(2) = t_m23(x2,5)
c      x3(3) = y_m23(x2,5)
c      x3(4) = p_m23(x2,5)
c      x3(5) = x2(5)
c etc... 
c for the central trajectory you'll get zeros all the way through
c  (that's why it's the central trajectory!)
c
c For the functions to be of real use you'll want to know the
c parameters of the central trajectory at each plane in the 
c coordinate system centered on the magnet axis. That is tabulated
c above: (see useful endplane information)
c The magnets are cylinders 1.2 m in diameter extending from the
c entrance to the exit
c y and phi are always zero for the central trajectory
c +x is down.
c x^ x y^ = z^ (z^ points downstream along the magnet axis)
c     
c  l_m35_i_j(x,5) gives the difference between the ray in question and
c             the central ray going from plane i to plane j
c 


      function x_m35_0_1(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 3)
      data ncoeff/  2/
      data avdat/  0.5417693E-03/
      data xmin/
     1 -0.99979E-04,-0.21093E+00,-0.59989E-01,-0.36658E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.21092E+00, 0.59989E-01, 0.36658E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.55125792E-03, 0.25311120E+00,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
c          set up monomials   functions
      x21 = x2
c
c                  function
c
      x_m35_0_1=avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
c
      return
      end
      function t_m35_0_1(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 3)
      data ncoeff/  2/
      data avdat/  0.4523665E-03/
      data xmin/
     1 -0.99979E-04,-0.21093E+00,-0.59989E-01,-0.36658E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.21092E+00, 0.59989E-01, 0.36658E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.45937650E-03, 0.21092500E+00,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
c          set up monomials   functions
      x21 = x2
c
c                  function
c
      t_m35_0_1=avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
c
      return
      end
      function y_m35_0_1(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 3)
      data ncoeff/  2/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.99979E-04,-0.21093E+00,-0.59989E-01,-0.36658E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.21092E+00, 0.59989E-01, 0.36658E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.59988703E-01, 0.43989960E-01,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
c          set up monomials   functions
      x31 = x3
      x41 = x4
c
c                  function
c
      y_m35_0_1=avdat
     1  +coeff( 1)        *x31        
     2  +coeff( 2)            *x41    
c
      return
      end
      function p_m35_0_1(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 3)
      data ncoeff/  2/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.99979E-04,-0.21093E+00,-0.59989E-01,-0.36658E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.21092E+00, 0.59989E-01, 0.36658E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.59971670E-11, 0.36658300E-01,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
c          set up monomials   functions
      x31 = x3
      x41 = x4
c
c                  function
c
      p_m35_0_1=avdat
     1  +coeff( 1)        *x31        
     2  +coeff( 2)            *x41    
c
      return
      end
      function l_m35_0_1(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 4)
      data ncoeff/  3/
      data avdat/ -0.8823982E-02/
      data xmin/
     1 -0.99979E-04,-0.21093E+00,-0.59989E-01,-0.36658E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.21092E+00, 0.59989E-01, 0.36658E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.87984860E-02,-0.26439381E-01,-0.80139520E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x41 = x4
      x42 = x41*x4
c
c                  function
c
      l_m35_0_1=avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x22            
     3  +coeff( 3)            *x42    
c
      return
      end
      function x_m35_1_2   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 4)
      data ncoeff/  3/
      data avdat/  0.2328079E-02/
      data xmin/
     1 -0.25316E+00,-0.21093E+00,-0.10185E+00,-0.36658E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.25308E+00, 0.21092E+00, 0.10185E+00, 0.36658E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.23387390E-02, 0.28796663E+00, 0.53015244E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
c
c                  function
c
      x_m35_1_2   =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)    *x22            
c
      return
      end
      function t_m35_1_2   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 4)
      data ncoeff/  3/
      data avdat/  0.4287185E-03/
      data xmin/
     1 -0.25316E+00,-0.21093E+00,-0.10185E+00,-0.36658E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.25308E+00, 0.21092E+00, 0.10185E+00, 0.36658E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.45980210E-03, 0.21025604E+00, 0.87730110E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x21 = x2
      x51 = x5
c
c                  function
c
      t_m35_1_2   =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)                *x51
c
      return
      end
      function y_m35_1_2   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 5)
      data ncoeff/  4/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.25316E+00,-0.21093E+00,-0.10185E+00,-0.36658E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.25308E+00, 0.21092E+00, 0.10185E+00, 0.36658E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.10185780E+00, 0.58723804E-02,-0.27154390E-05, 0.91975260E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
c          set up monomials   functions
      x21 = x2
      x31 = x3
      x41 = x4
c
c                  function
c
      y_m35_1_2   =avdat
     1  +coeff( 1)        *x31        
     2  +coeff( 2)            *x41    
     3  +coeff( 3)    *x21*x31        
     4  +coeff( 4)    *x21    *x41    
c
      return
      end
      function p_m35_1_2   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 4)
      data ncoeff/  3/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.25316E+00,-0.21093E+00,-0.10185E+00,-0.36658E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.25308E+00, 0.21092E+00, 0.10185E+00, 0.36658E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.16711393E-03, 0.36800650E-01, 0.49313570E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
c          set up monomials   functions
      x21 = x2
      x31 = x3
      x41 = x4
c
c                  function
c
      p_m35_1_2   =avdat
     1  +coeff( 1)        *x31        
     2  +coeff( 2)            *x41    
     3  +coeff( 3)    *x21*x31        
c
      return
      end
      function l_m35_1_2   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 4)
      data ncoeff/  3/
      data avdat/ -0.1380061E-02/
      data xmin/
     1 -0.25316E+00,-0.21093E+00,-0.10185E+00,-0.36658E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.25308E+00, 0.21092E+00, 0.10185E+00, 0.36658E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.13452793E-02,-0.25422750E-01,-0.39954222E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
c
c                  function
c
      l_m35_1_2   =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)    *x22            
c
      return
      end
      function x_m35_2_3   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(15)
      data ncoeff/ 14/
      data avdat/ -0.1945024E-01/
      data xmin/
     1 -0.28276E+00,-0.21036E+00,-0.10799E+00,-0.37129E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.27707E+00, 0.19855E+00, 0.10799E+00, 0.37129E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.26607224E-02, 0.58154410E+00, 0.83938700E-02, 0.15762573E-01,
     +  0.62103620E-02,-0.54027810E-02,-0.12141060E-02,-0.22815340E-02,
     +  0.16326220E-02, 0.27823564E-03,-0.74605940E-04, 0.62688061E-03,
     + -0.14668462E-02, 0.21776260E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_m35_2_3   =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)                *x51
     4  +coeff( 4)    *x21        *x51
     5  +coeff( 5)    *x22            
     6  +coeff( 6)    *x23            
     7  +coeff( 7)                *x52
     8  +coeff( 8)    *x21        *x52
      x_m35_2_3   =x_m35_2_3   
     9  +coeff( 9)    *x24            
     1  +coeff(10)    *x22        *x51
     2  +coeff(11)    *x22*x32        
     3  +coeff(12)    *x23        *x53
     4  +coeff(13)*x12*x22*x32        
     5  +coeff(14)        *x32        
c
      return
      end
      function t_m35_2_3   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(24)
      data ncoeff/ 23/
      data avdat/ -0.3954356E-02/
      data xmin/
     1 -0.28276E+00,-0.21036E+00,-0.10799E+00,-0.37129E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.27707E+00, 0.19855E+00, 0.10799E+00, 0.37129E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.80038200E-03, 0.13730690E+00, 0.61635584E-02,-0.83649990E-01,
     +  0.87644670E-01,-0.27104580E-03, 0.35558690E-03,-0.95237774E-04,
     +  0.35368820E+00,-0.19125760E-02,-0.87962850E-01,-0.33747693E-02,
     + -0.33315460E+00, 0.69379370E-02, 0.99279200E-04,-0.29729662E-02,
     +  0.29959084E-03,-0.33990580E-02, 0.18886610E-02,-0.40882182E-03,
     + -0.17820000E-02,-0.46216353E-03, 0.22484414E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_m35_2_3   =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)                *x51
     4  +coeff( 4)*x11                
     5  +coeff( 5)    *x22            
     6  +coeff( 6)        *x32        
     7  +coeff( 7)        *x31*x41    
     8  +coeff( 8)            *x42    
      t_m35_2_3   =t_m35_2_3   
     9  +coeff( 9)    *x21        *x51
     1  +coeff(10)                *x52
     2  +coeff(11)*x11*x21            
     3  +coeff(12)    *x23            
     4  +coeff(13)*x11            *x51
     5  +coeff(14)    *x22        *x51
     6  +coeff(15)        *x32    *x51
     7  +coeff(16)    *x21        *x52
     8  +coeff(17)                *x53
      t_m35_2_3   =t_m35_2_3   
     9  +coeff(18)*x12                
     1  +coeff(19)    *x24            
     2  +coeff(20)    *x22*x32        
     3  +coeff(21)    *x22*x31*x41    
     4  +coeff(22)*x11*x21        *x51
     5  +coeff(23)    *x23        *x51
c
      return
      end
      function y_m35_2_3   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(21)
      data ncoeff/ 20/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.28276E+00,-0.21036E+00,-0.10799E+00,-0.37129E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.27707E+00, 0.19855E+00, 0.10799E+00, 0.37129E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.13488982E+00, 0.81760090E-01,-0.93235330E-02, 0.28595801E-01,
     + -0.40178570E-02,-0.12731614E-02, 0.60420060E-02,-0.29529750E-01,
     + -0.13846340E-01, 0.67005030E-03,-0.25196100E-03,-0.74096950E-03,
     +  0.19741860E-03, 0.15838842E-03, 0.14876290E-03, 0.28650190E-03,
     +  0.14256443E-01, 0.31041892E-03,-0.13823701E-02,-0.14651660E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      y_m35_2_3   =avdat
     1  +coeff( 1)        *x31        
     2  +coeff( 2)            *x41    
     3  +coeff( 3)    *x21*x31        
     4  +coeff( 4)    *x21    *x41    
     5  +coeff( 5)        *x31    *x51
     6  +coeff( 6)            *x41*x51
     7  +coeff( 7)*x11    *x31        
     8  +coeff( 8)*x11        *x41    
      y_m35_2_3   =y_m35_2_3   
     9  +coeff( 9)    *x21*x31    *x51
     1  +coeff(10)        *x31    *x52
     2  +coeff(11)    *x23*x31        
     3  +coeff(12)    *x23    *x41    
     4  +coeff(13)            *x41*x52
     5  +coeff(14)        *x32*x41    
     6  +coeff(15)    *x21    *x41*x51
     7  +coeff(16)    *x21*x33        
     8  +coeff(17)*x11    *x31    *x51
      y_m35_2_3   =y_m35_2_3   
     9  +coeff(18)    *x23*x31    *x51
     1  +coeff(19)*x11*x24*x31        
     2  +coeff(20)    *x22*x31        
c
      return
      end
      function p_m35_2_3   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(21)
      data ncoeff/ 20/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.28276E+00,-0.21036E+00,-0.10799E+00,-0.37129E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.27707E+00, 0.19855E+00, 0.10799E+00, 0.37129E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.35352830E-01, 0.50429612E-01,-0.15396301E-01, 0.35914611E-01,
     + -0.57900464E-02,-0.22779284E-02, 0.13619700E-01,-0.14334600E-01,
     +  0.14208092E-03,-0.36409430E-01, 0.43288452E-03, 0.68901324E-04,
     +  0.93795632E-03, 0.37465780E-03, 0.13323782E-01,-0.48935570E-04,
     + -0.14026943E-02, 0.24494982E-03, 0.12804040E-02,-0.30469424E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      p_m35_2_3   =avdat
     1  +coeff( 1)        *x31        
     2  +coeff( 2)            *x41    
     3  +coeff( 3)    *x21*x31        
     4  +coeff( 4)    *x21    *x41    
     5  +coeff( 5)        *x31    *x51
     6  +coeff( 6)            *x41*x51
     7  +coeff( 7)*x11    *x31        
     8  +coeff( 8)    *x22*x31        
      p_m35_2_3   =p_m35_2_3   
     9  +coeff( 9)        *x33        
     1  +coeff(10)*x11        *x41    
     2  +coeff(11)    *x21*x31    *x51
     3  +coeff(12)    *x21    *x41*x51
     4  +coeff(13)        *x31    *x52
     5  +coeff(14)            *x41*x52
     6  +coeff(15)*x11*x21*x31        
     7  +coeff(16)    *x21*x33        
     8  +coeff(17)    *x23    *x41    
      p_m35_2_3   =p_m35_2_3   
     9  +coeff(18)        *x31*x44    
     1  +coeff(19)    *x23*x32*x41    
     2  +coeff(20)*x11*x24*x31        
c
      return
      end
      function l_m35_2_3   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(16)
      data ncoeff/ 15/
      data avdat/ -0.7802763E-02/
      data xmin/
     1 -0.28276E+00,-0.21036E+00,-0.10799E+00,-0.37129E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.27707E+00, 0.19855E+00, 0.10799E+00, 0.37129E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.84143860E-02,-0.22708991E-01, 0.44797240E-02,-0.19714140E-02,
     + -0.19951600E-01,-0.12163131E-02,-0.16990393E-02, 0.92726360E-01,
     + -0.94147540E-01, 0.29206570E-03, 0.21986314E-04, 0.75342511E-03,
     + -0.29283520E-03, 0.30549513E-03,-0.50610340E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
c
c                  function
c
      l_m35_2_3   =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)*x11                
     4  +coeff( 4)                *x51
     5  +coeff( 5)    *x22            
     6  +coeff( 6)        *x31*x41    
     7  +coeff( 7)            *x42    
     8  +coeff( 8)    *x21        *x51
      l_m35_2_3   =l_m35_2_3   
     9  +coeff( 9)*x11            *x51
     1  +coeff(10)    *x21*x31*x41    
     2  +coeff(11)        *x32    *x51
     3  +coeff(12)    *x24            
     4  +coeff(13)        *x32        
     5  +coeff(14)        *x31*x41*x51
     6  +coeff(15)*x11*x21            
c
      return
      end
      function x_m35_3_4   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/ -0.2748773E-02/
      data xmin/
     1 -0.56383E+00,-0.56436E-01,-0.21189E+00,-0.88687E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.55229E+00, 0.67502E-01, 0.21189E+00, 0.88687E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.26717310E-01, 0.28711100E+00,-0.62263463E-02,-0.22009630E-01,
     +  0.53661644E-01,-0.61304920E-02, 0.17800760E+00,-0.47716131E-03,
     +  0.98305102E-02,-0.63949720E-02, 0.16883714E-01, 0.20692000E-01,
     + -0.82453554E-02,-0.12854990E-01,-0.10877380E-01, 0.45435640E-02,
     +  0.31138623E-01,-0.93650180E-02,-0.14520640E-01, 0.41304621E-03,
     +  0.21962501E-01,-0.41290380E-03,-0.27831792E-01,-0.20197140E-01,
     +  0.12884982E-02,-0.10442494E-01,-0.32376144E-01, 0.73284213E-02,
     + -0.22224960E-01, 0.37551400E-01, 0.13690570E-02,-0.39721320E-01,
     + -0.22297690E-02,-0.80587472E-02, 0.17563360E+00,-0.58322292E-02,
     + -0.11383330E+00,-0.13799681E-01, 0.21190490E-02, 0.11418340E-02,
     + -0.55351834E-01, 0.44659632E-02, 0.27261620E-02,-0.33781144E-02,
     +  0.18336564E-01,-0.14354163E-01, 0.29003501E-02, 0.95144660E-02,
     +  0.47327661E-02,-0.47593220E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x25 = x24*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x35 = x34*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      x_m35_3_4   =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)                *x51
     4  +coeff( 4)    *x22            
     5  +coeff( 5)    *x21        *x51
     6  +coeff( 6)                *x52
     7  +coeff( 7)*x11                
     8  +coeff( 8)        *x32        
      x_m35_3_4   =x_m35_3_4   
     9  +coeff( 9)*x11    *x32        
     1  +coeff(10)        *x31*x41    
     2  +coeff(11)    *x22        *x51
     3  +coeff(12)    *x23*x32        
     4  +coeff(13)*x12    *x32        
     5  +coeff(14)    *x23            
     6  +coeff(15)    *x21*x32        
     7  +coeff(16)    *x21*x31*x41    
     8  +coeff(17)*x11*x21            
      x_m35_3_4   =x_m35_3_4   
     9  +coeff(18)*x11*x21        *x51
     1  +coeff(19)*x11    *x31*x41    
     2  +coeff(20)    *x25            
     3  +coeff(21)    *x24*x32        
     4  +coeff(22)    *x22*x34        
     5  +coeff(23)    *x23*x34        
     6  +coeff(24)*x11*x23*x31*x41    
     7  +coeff(25)    *x24        *x51
     8  +coeff(26)    *x22*x32    *x51
      x_m35_3_4   =x_m35_3_4   
     9  +coeff(27)    *x23*x33*x41    
     1  +coeff(28)    *x25    *x42    
     2  +coeff(29)    *x25*x34        
     3  +coeff(30)    *x23*x35*x41    
     4  +coeff(31)*x12    *x32*x42    
     5  +coeff(32)*x12*x23*x31*x41    
     6  +coeff(33)*x12*x22*x32    *x51
     7  +coeff(34)    *x25*x34    *x51
     8  +coeff(35)*x11*x24*x33*x41    
      x_m35_3_4   =x_m35_3_4   
     9  +coeff(36)*x12*x24*x32        
     1  +coeff(37)    *x25*x34*x42    
     2  +coeff(38)*x12*x25*x32        
     3  +coeff(39)        *x32    *x51
     4  +coeff(40)    *x21        *x52
     5  +coeff(41)*x11            *x51
     6  +coeff(42)    *x22*x32        
     7  +coeff(43)    *x23        *x51
     8  +coeff(44)    *x21*x32    *x51
      x_m35_3_4   =x_m35_3_4   
     9  +coeff(45)*x11*x22            
     1  +coeff(46)*x11            *x52
     2  +coeff(47)    *x21*x34        
     3  +coeff(48)    *x22*x31*x41*x51
     4  +coeff(49)*x11    *x31*x41*x51
     5  +coeff(50)*x11*x21        *x52
c
      return
      end
      function t_m35_3_4   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/ -0.4158733E-02/
      data xmin/
     1 -0.56383E+00,-0.56436E-01,-0.21189E+00,-0.88687E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.55229E+00, 0.67502E-01, 0.21189E+00, 0.88687E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.12140954E-01, 0.79759450E-01,-0.46244713E-04,-0.17528050E+00,
     +  0.11001740E+00,-0.12141420E-02,-0.24261230E-02,-0.97808130E-01,
     +  0.13515993E-01,-0.83825410E-01, 0.43201081E-01,-0.11935153E-01,
     +  0.18632850E-01, 0.82290962E-01, 0.11187090E-02,-0.17137514E-02,
     + -0.11591220E-01,-0.16764310E-02,-0.12181860E-01, 0.69412840E-02,
     + -0.55115630E-03,-0.19420964E-01, 0.34829311E-01, 0.20241610E-01,
     +  0.15676610E-01,-0.12786212E-01,-0.24496370E-01,-0.26385001E-01,
     + -0.61707783E-01, 0.42206370E-01, 0.42352914E-01,-0.26758540E-02,
     +  0.54056950E-02,-0.65358910E-02,-0.48410800E-02,-0.45821940E-04,
     + -0.88627412E-02,-0.58328583E-02,-0.52571830E-01,-0.99212950E-02,
     +  0.26408690E-02, 0.20718682E-01, 0.99729122E-02,-0.45092360E-02,
     + -0.31496380E-02, 0.53744590E-02, 0.97120400E-03, 0.83862461E-04,
     +  0.97393310E-03,-0.93529870E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_m35_3_4   =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)                *x51
     4  +coeff( 4)*x11                
     5  +coeff( 5)    *x22            
     6  +coeff( 6)        *x32        
     7  +coeff( 7)        *x31*x41    
     8  +coeff( 8)    *x21        *x51
      t_m35_3_4   =t_m35_3_4   
     9  +coeff( 9)                *x52
     1  +coeff(10)*x11*x21            
     2  +coeff(11)    *x23            
     3  +coeff(12)    *x21*x32        
     4  +coeff(13)    *x21*x31*x41    
     5  +coeff(14)*x11            *x51
     6  +coeff(15)        *x32    *x51
     7  +coeff(16)        *x31*x41*x51
     8  +coeff(17)    *x21        *x52
      t_m35_3_4   =t_m35_3_4   
     9  +coeff(18)*x12                
     1  +coeff(19)*x11*x22            
     2  +coeff(20)*x11    *x32        
     3  +coeff(21)    *x22*x32        
     4  +coeff(22)*x11    *x31*x41    
     5  +coeff(23)*x11            *x52
     6  +coeff(24)    *x23*x32        
     7  +coeff(25)    *x21*x34        
     8  +coeff(26)    *x21*x33*x41    
      t_m35_3_4   =t_m35_3_4   
     9  +coeff(27)*x11*x22*x31*x41    
     1  +coeff(28)*x12*x21        *x51
     2  +coeff(29)    *x23*x34        
     3  +coeff(30)    *x23*x32*x42    
     4  +coeff(31)*x11*x22*x34        
     5  +coeff(32)    *x22*x31*x41    
     6  +coeff(33)    *x21*x34    *x51
     7  +coeff(34)*x12*x21*x32        
     8  +coeff(35)    *x21*x34*x42    
      t_m35_3_4   =t_m35_3_4   
     9  +coeff(36)*x11*x24*x32        
     1  +coeff(37)    *x21*x34*x42*x51
     2  +coeff(38)    *x24            
     3  +coeff(39)*x11*x21        *x51
     4  +coeff(40)*x12*x21            
     5  +coeff(41)*x11*x21*x32        
     6  +coeff(42)*x12*x22            
     7  +coeff(43)*x11*x24            
     8  +coeff(44)    *x24*x32        
      t_m35_3_4   =t_m35_3_4   
     9  +coeff(45)    *x23*x32    *x51
     1  +coeff(46)*x12            *x52
     2  +coeff(47)*x11*x22        *x52
     3  +coeff(48)    *x24        *x52
     4  +coeff(49)    *x21*x32    *x53
     5  +coeff(50)*x12*x23            
c
      return
      end
      function y_m35_3_4   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.56383E+00,-0.56436E-01,-0.21189E+00,-0.88687E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.55229E+00, 0.67502E-01, 0.21189E+00, 0.88687E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.30300652E+00, 0.23122082E+00,-0.11410540E+00,-0.18770954E-03,
     +  0.11748520E-01,-0.46640341E-02, 0.93420043E-01,-0.14532811E+00,
     + -0.26228420E-02, 0.11294714E+00,-0.14881643E-01, 0.20249854E+00,
     + -0.11851283E+00,-0.39183041E-02, 0.13146131E-01,-0.54111440E-01,
     + -0.67468903E-01, 0.23729590E+00, 0.28468130E-02, 0.69255850E-02,
     + -0.12986910E+00, 0.50593394E-01, 0.55572181E-02, 0.80804610E-01,
     + -0.99949724E-02,-0.85395844E-02, 0.25540221E-01,-0.56779150E-01,
     + -0.16211740E-01,-0.23905150E-01,-0.97017700E-02, 0.52964771E-02,
     +  0.56588900E-02,-0.56280391E-02,-0.56489560E-02, 0.52592091E-01,
     + -0.55623050E-01,-0.11767683E-01,-0.21076410E-01, 0.10307560E-02,
     + -0.49346190E-02,-0.53900804E-01,-0.27875721E-01, 0.45195524E-02,
     +  0.73143653E-03,-0.18428083E-01, 0.36404814E-01,-0.94494810E-02,
     + -0.23892760E-03, 0.26493990E-01,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      y_m35_3_4   =avdat
     1  +coeff( 1)        *x31        
     2  +coeff( 2)            *x41    
     3  +coeff( 3)    *x21*x31        
     4  +coeff( 4)    *x21    *x41    
     5  +coeff( 5)        *x31    *x51
     6  +coeff( 6)            *x41*x51
     7  +coeff( 7)*x11    *x31        
     8  +coeff( 8)    *x22*x31        
      y_m35_3_4   =y_m35_3_4   
     9  +coeff( 9)        *x33        
     1  +coeff(10)    *x21*x31    *x51
     2  +coeff(11)        *x31    *x52
     3  +coeff(12)*x11*x21*x31        
     4  +coeff(13)    *x23*x31        
     5  +coeff(14)    *x21*x33        
     6  +coeff(15)    *x23    *x41    
     7  +coeff(16)*x11    *x31    *x51
     8  +coeff(17)*x12    *x31        
      y_m35_3_4   =y_m35_3_4   
     9  +coeff(18)*x11*x22*x31        
     1  +coeff(19)*x11    *x33        
     2  +coeff(20)    *x22*x33        
     3  +coeff(21)*x12*x21*x31        
     4  +coeff(22)*x11*x24*x31        
     5  +coeff(23)        *x32*x41    
     6  +coeff(24)    *x22*x31    *x51
     7  +coeff(25)    *x21*x31    *x52
     8  +coeff(26)*x11*x22*x33        
      y_m35_3_4   =y_m35_3_4   
     9  +coeff(27)    *x22*x34*x41    
     1  +coeff(28)*x12*x23*x31        
     2  +coeff(29)*x11*x22    *x41    
     3  +coeff(30)    *x22*x32*x41    
     4  +coeff(31)    *x22*x31    *x52
     5  +coeff(32)*x12        *x43    
     6  +coeff(33)    *x24*x31    *x53
     7  +coeff(34)*x12    *x34*x41    
     8  +coeff(35)    *x22*x34*x43    
      y_m35_3_4   =y_m35_3_4   
     9  +coeff(36)*x12*x24*x33        
     1  +coeff(37)*x12*x22*x34*x41    
     2  +coeff(38)*x12*x23*x33    *x51
     3  +coeff(39)*x12*x24*x33    *x51
     4  +coeff(40)        *x31    *x53
     5  +coeff(41)    *x24*x31        
     6  +coeff(42)*x11*x21*x31    *x51
     7  +coeff(43)*x11    *x31    *x52
     8  +coeff(44)    *x21*x31    *x53
      y_m35_3_4   =y_m35_3_4   
     9  +coeff(45)*x11*x21*x33        
     1  +coeff(46)    *x24*x33        
     2  +coeff(47)    *x24*x32*x41    
     3  +coeff(48)    *x24*x31    *x52
     4  +coeff(49)    *x21*x32*x41*x53
     5  +coeff(50)*x12*x21*x33        
c
      return
      end
      function p_m35_3_4   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.56383E+00,-0.56436E-01,-0.21189E+00,-0.88687E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.55229E+00, 0.67502E-01, 0.21189E+00, 0.88687E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.61282310E-01, 0.11165000E+00,-0.95070153E-01,-0.18888990E-02,
     +  0.68378932E-02,-0.32410330E-02, 0.76457160E-01,-0.66664911E-01,
     + -0.71995640E-02, 0.10385982E-01,-0.30508960E-02, 0.21882490E-01,
     + -0.18292370E-02, 0.92362950E-01,-0.26755520E-01,-0.53330743E-02,
     + -0.76838531E-02,-0.24061670E-02, 0.11892213E-01, 0.10269430E-02,
     +  0.77310280E-01, 0.28744204E-02, 0.58064941E-01,-0.11449231E+00,
     +  0.12354673E-02,-0.63093320E-02,-0.55928380E-01, 0.56008980E-02,
     +  0.35671282E-01, 0.64584960E-03, 0.43418630E-01, 0.43813490E-02,
     + -0.47900222E-01,-0.26544800E-01, 0.58440933E-02,-0.63771504E-03,
     + -0.37181600E-01,-0.51253191E-02,-0.63592512E-02,-0.37389550E-01,
     +  0.20250120E-01, 0.39149392E-02, 0.31850520E-01, 0.19977733E-02,
     +  0.41370873E-03,-0.49137940E-03, 0.10121590E-01,-0.90108142E-03,
     + -0.23791230E-02, 0.32878450E-01,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_m35_3_4   =avdat
     1  +coeff( 1)        *x31        
     2  +coeff( 2)            *x41    
     3  +coeff( 3)    *x21*x31        
     4  +coeff( 4)    *x21    *x41    
     5  +coeff( 5)        *x31    *x51
     6  +coeff( 6)            *x41*x51
     7  +coeff( 7)*x11    *x31        
     8  +coeff( 8)    *x22*x31        
      p_m35_3_4   =p_m35_3_4   
     9  +coeff( 9)        *x33        
     1  +coeff(10)        *x32*x41    
     2  +coeff(11)        *x31*x42    
     3  +coeff(12)    *x21*x31    *x51
     4  +coeff(13)        *x31    *x52
     5  +coeff(14)*x11*x21*x31        
     6  +coeff(15)    *x23*x31        
     7  +coeff(16)    *x21*x33        
     8  +coeff(17)*x11*x21    *x41    
      p_m35_3_4   =p_m35_3_4   
     9  +coeff(18)    *x23    *x41    
     1  +coeff(19)*x11    *x31    *x51
     2  +coeff(20)*x11        *x41*x51
     3  +coeff(21)*x11*x22*x31        
     4  +coeff(22)*x11    *x33        
     5  +coeff(23)    *x22*x33        
     6  +coeff(24)    *x22*x32*x41    
     7  +coeff(25)        *x34*x41    
     8  +coeff(26)*x11    *x31    *x52
      p_m35_3_4   =p_m35_3_4   
     9  +coeff(27)*x12*x21*x31        
     1  +coeff(28)    *x23*x33        
     2  +coeff(29)*x11*x21*x32*x41    
     3  +coeff(30)        *x33    *x53
     4  +coeff(31)*x11*x24*x31        
     5  +coeff(32)    *x22*x34*x41    
     6  +coeff(33)*x12*x23*x31        
     7  +coeff(34)*x11*x21*x34*x41    
     8  +coeff(35)    *x22    *x41    
      p_m35_3_4   =p_m35_3_4   
     9  +coeff(36)        *x31    *x53
     1  +coeff(37)*x12    *x31        
     2  +coeff(38)*x11*x21*x31    *x51
     3  +coeff(39)    *x23*x31    *x51
     4  +coeff(40)    *x24*x33        
     5  +coeff(41)    *x24*x32*x41    
     6  +coeff(42)    *x21*x34*x41*x51
     7  +coeff(43)*x11*x23*x32*x41    
     8  +coeff(44)        *x34*x43*x51
      p_m35_3_4   =p_m35_3_4   
     9  +coeff(45)            *x41*x52
     1  +coeff(46)    *x21*x31*x42    
     2  +coeff(47)    *x22*x31    *x51
     3  +coeff(48)        *x31*x42*x51
     4  +coeff(49)    *x24    *x41    
     5  +coeff(50)    *x22*x31*x42    
c
      return
      end
      function l_m35_3_4   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(31)
      data ncoeff/ 30/
      data avdat/ -0.8340841E-02/
      data xmin/
     1 -0.56383E+00,-0.56436E-01,-0.21189E+00,-0.88687E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.55229E+00, 0.67502E-01, 0.21189E+00, 0.88687E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.10070934E-01, 0.15610230E-01,-0.74050230E-04,-0.24444210E-01,
     + -0.16155421E-02,-0.28104300E-01,-0.12777860E-01, 0.33729710E-01,
     + -0.99659900E-02,-0.47449690E-03,-0.86682802E-03,-0.44736373E-02,
     + -0.10875750E-01, 0.35568724E-02, 0.32487810E-02,-0.94556080E-03,
     + -0.23465224E-02, 0.25951040E-02, 0.16290361E-02,-0.49022640E-03,
     + -0.26396430E-03,-0.35058911E-02, 0.15896240E-02, 0.43321060E-02,
     + -0.22093241E-02,-0.34440451E-03, 0.17977904E-02,-0.79755480E-03,
     +  0.17632710E-03,-0.10675490E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_m35_3_4   =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)                *x51
     4  +coeff( 4)*x11                
     5  +coeff( 5)        *x32        
     6  +coeff( 6)*x12                
     7  +coeff( 7)        *x31*x41    
     8  +coeff( 8)*x11*x21            
      l_m35_3_4   =l_m35_3_4   
     9  +coeff( 9)    *x22            
     1  +coeff(10)        *x32    *x51
     2  +coeff(11)    *x23            
     3  +coeff(12)*x11    *x32        
     4  +coeff(13)            *x42    
     5  +coeff(14)    *x21*x32        
     6  +coeff(15)        *x31*x41*x51
     7  +coeff(16)*x11*x22            
     8  +coeff(17)    *x22*x32        
      l_m35_3_4   =l_m35_3_4   
     9  +coeff(18)*x12*x21            
     1  +coeff(19)    *x22*x32    *x51
     2  +coeff(20)    *x22*x32    *x52
     3  +coeff(21)    *x21*x32    *x51
     4  +coeff(22)    *x23*x32        
     5  +coeff(23)*x11*x21*x31*x41    
     6  +coeff(24)*x11*x22*x32        
     7  +coeff(25)*x12        *x42    
     8  +coeff(26)        *x34*x42    
      l_m35_3_4   =l_m35_3_4   
     9  +coeff(27)    *x22*x31*x43    
     1  +coeff(28)            *x44*x52
     2  +coeff(29)*x11*x22        *x52
     3  +coeff(30)*x11*x23*x32        
c
      return
      end
      function x_m35_4_5   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(32)
      data ncoeff/ 31/
      data avdat/ -0.5592477E-02/
      data xmin/
     1 -0.42252E+00,-0.19647E+00,-0.54023E+00,-0.18368E+00,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.48443E+00, 0.14310E+00, 0.54023E+00, 0.18368E+00, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.32194584E-01, 0.25381210E-01, 0.30093844E-02, 0.12537691E-02,
     +  0.56353650E-04,-0.34660933E-03, 0.10122662E-03,-0.54142843E-02,
     +  0.49275491E-03, 0.44553470E+00, 0.19805550E-01, 0.13066762E-03,
     + -0.12207764E-04, 0.95017731E-03,-0.48099573E-04, 0.18080910E-03,
     + -0.29877884E-03,-0.61022740E-03,-0.45669361E-03,-0.10305200E-02,
     +  0.25705780E-03, 0.74450740E-03, 0.57231360E-03,-0.38496470E-02,
     +  0.27391673E-02, 0.18268240E-03,-0.41624920E-03, 0.11260690E-02,
     + -0.10550324E-02, 0.59883830E-03,-0.19343360E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x25 = x24*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x35 = x34*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      x_m35_4_5   =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)                *x51
     4  +coeff( 4)    *x22            
     5  +coeff( 5)        *x32        
     6  +coeff( 6)        *x31*x41    
     7  +coeff( 7)            *x42    
     8  +coeff( 8)    *x21        *x51
      x_m35_4_5   =x_m35_4_5   
     9  +coeff( 9)                *x52
     1  +coeff(10)*x11                
     2  +coeff(11)*x11*x21            
     3  +coeff(12)    *x23            
     4  +coeff(13)*x11*x21*x32        
     5  +coeff(14)*x11    *x34        
     6  +coeff(15)        *x32    *x51
     7  +coeff(16)        *x34        
     8  +coeff(17)    *x23        *x51
      x_m35_4_5   =x_m35_4_5   
     9  +coeff(18)    *x21        *x52
     1  +coeff(19)    *x22*x32        
     2  +coeff(20)*x11    *x32        
     3  +coeff(21)*x11            *x52
     4  +coeff(22)    *x25            
     5  +coeff(23)    *x25        *x51
     6  +coeff(24)*x12    *x32        
     7  +coeff(25)*x12    *x31*x41    
     8  +coeff(26)*x12*x21        *x51
      x_m35_4_5   =x_m35_4_5   
     9  +coeff(27)*x11*x25            
     1  +coeff(28)*x12*x21*x32        
     2  +coeff(29)*x11    *x35*x41    
     3  +coeff(30)*x12*x23        *x51
     4  +coeff(31)*x12*x23*x33*x41    
c
      return
      end
      function t_m35_4_5   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/ -0.3624947E-02/
      data xmin/
     1 -0.42252E+00,-0.19647E+00,-0.54023E+00,-0.18368E+00,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.48443E+00, 0.14310E+00, 0.54023E+00, 0.18368E+00, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.22521601E-01, 0.16742821E+00, 0.14428500E-02,-0.15683330E-03,
     + -0.12337120E-01,-0.24468931E-02, 0.13787790E-01,-0.67724003E-02,
     + -0.10940582E-01,-0.25241590E-01, 0.67795560E-03, 0.57897474E-02,
     +  0.15075932E-02, 0.72804060E-02,-0.10755310E-02,-0.62030390E-04,
     +  0.11406080E-02, 0.15040942E-02,-0.89618563E-02,-0.14533240E-01,
     +  0.81469900E-02, 0.16949821E-02,-0.85058150E-03,-0.82417484E-02,
     + -0.31139723E-01,-0.56035274E-04,-0.39721000E-02, 0.22848570E-02,
     +  0.63485262E-03, 0.64356060E-04,-0.20591642E-01, 0.19233200E-02,
     + -0.18346462E-02, 0.88264903E-03,-0.29757884E-02,-0.10520251E-02,
     +  0.20579544E-02, 0.18777954E-02, 0.25197050E-02,-0.44173550E-02,
     + -0.14983410E-02, 0.24737843E-02, 0.22669970E-02,-0.23060930E-02,
     + -0.47371350E-02, 0.59830340E-01,-0.27945460E-02, 0.35015100E-03,
     +  0.83656840E-04, 0.20829271E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_m35_4_5   =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)                *x51
     4  +coeff( 4)*x11                
     5  +coeff( 5)    *x22            
     6  +coeff( 6)        *x32        
     7  +coeff( 7)        *x31*x41    
     8  +coeff( 8)            *x42    
      t_m35_4_5   =t_m35_4_5   
     9  +coeff( 9)*x11*x21            
     1  +coeff(10)    *x21*x32        
     2  +coeff(11)    *x21*x31*x41    
     3  +coeff(12)        *x32    *x51
     4  +coeff(13)    *x24            
     5  +coeff(14)    *x22*x32        
     6  +coeff(15)        *x34        
     7  +coeff(16)*x11    *x31*x41    
     8  +coeff(17)        *x33*x41    
      t_m35_4_5   =t_m35_4_5   
     9  +coeff(18)*x11*x21        *x51
     1  +coeff(19)    *x21*x32    *x51
     2  +coeff(20)*x11*x21*x32        
     3  +coeff(21)    *x21*x34        
     4  +coeff(22)    *x21        *x51
     5  +coeff(23)*x12*x21            
     6  +coeff(24)    *x21*x33*x41    
     7  +coeff(25)*x11*x21    *x42    
     8  +coeff(26)        *x34    *x51
      t_m35_4_5   =t_m35_4_5   
     9  +coeff(27)*x12*x21*x32        
     1  +coeff(28)    *x23            
     2  +coeff(29)*x11            *x51
     3  +coeff(30)    *x21        *x52
     4  +coeff(31)*x11    *x32        
     5  +coeff(32)    *x23*x32        
     6  +coeff(33)    *x24*x32        
     7  +coeff(34)        *x34*x42    
     8  +coeff(35)*x12    *x34        
      t_m35_4_5   =t_m35_4_5   
     9  +coeff(36)*x11    *x34*x42    
     1  +coeff(37)*x11*x23            
     2  +coeff(38)*x11    *x32    *x51
     3  +coeff(39)*x12    *x32        
     4  +coeff(40)*x11    *x32*x42    
     5  +coeff(41)    *x22*x32*x42    
     6  +coeff(42)*x12*x22*x32        
     7  +coeff(43)    *x24*x32*x42    
     8  +coeff(44)    *x22*x34*x42    
      t_m35_4_5   =t_m35_4_5   
     9  +coeff(45)    *x22        *x51
     1  +coeff(46)*x11*x21*x31*x41    
     2  +coeff(47)*x11*x22        *x51
     3  +coeff(48)*x11*x21        *x52
     4  +coeff(49)        *x32    *x53
     5  +coeff(50)*x12*x22            
c
      return
      end
      function y_m35_4_5   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 9)
      data ncoeff/  8/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.42252E+00,-0.19647E+00,-0.54023E+00,-0.18368E+00,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.48443E+00, 0.14310E+00, 0.54023E+00, 0.18368E+00, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.53834891E+00, 0.39821570E-01, 0.23170381E-01,-0.22210782E-01,
     + -0.55369432E-02, 0.56402771E-02, 0.19202502E-01,-0.68928620E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x41 = x4
      x51 = x5
c
c                  function
c
      y_m35_4_5   =avdat
     1  +coeff( 1)        *x31        
     2  +coeff( 2)            *x41    
     3  +coeff( 3)    *x21*x31        
     4  +coeff( 4)    *x21    *x41    
     5  +coeff( 5)        *x31    *x51
     6  +coeff( 6)            *x41*x51
     7  +coeff( 7)*x11    *x31        
     8  +coeff( 8)    *x22*x31        
c
      return
      end
      function p_m35_4_5   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.42252E+00,-0.19647E+00,-0.54023E+00,-0.18368E+00,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.48443E+00, 0.14310E+00, 0.54023E+00, 0.18368E+00, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.48974600E-02, 0.18331514E+00, 0.26618111E-02, 0.68364220E-02,
     + -0.41171760E-03,-0.16050183E-02, 0.62874943E-03, 0.20677690E-02,
     +  0.13001980E-02,-0.16214620E-02,-0.55982120E-03,-0.66750533E-02,
     + -0.40815430E-02, 0.34480480E-02,-0.93938022E-01,-0.63643880E-01,
     + -0.23439973E-01, 0.11599513E-02, 0.54499260E-01,-0.30695550E-01,
     +  0.14553910E-02,-0.48237343E-03, 0.55481223E-02,-0.10464143E-03,
     +  0.16169004E-01,-0.63605900E-02,-0.29152804E-01, 0.25224090E-02,
     +  0.98700060E-02, 0.10691091E+00,-0.26571940E-02, 0.96801680E+00,
     +  0.81389494E-01,-0.19790230E+01,-0.43415240E-02,-0.57873114E-01,
     +  0.39449851E-02, 0.46030632E-03,-0.18263700E-02, 0.66174594E-02,
     + -0.14554200E-01, 0.36941400E-02,-0.18161244E+00,-0.83288200E-03,
     +  0.10157410E+01, 0.88435600E-03, 0.19937093E+00,-0.13529682E+00,
     +  0.82235224E-02,-0.77345862E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      p_m35_4_5   =avdat
     1  +coeff( 1)        *x31        
     2  +coeff( 2)            *x41    
     3  +coeff( 3)    *x21*x31        
     4  +coeff( 4)    *x21    *x41    
     5  +coeff( 5)*x11    *x31        
     6  +coeff( 6)    *x22*x31        
     7  +coeff( 7)        *x33        
     8  +coeff( 8)        *x32*x41    
      p_m35_4_5   =p_m35_4_5   
     9  +coeff( 9)    *x21*x31    *x51
     1  +coeff(10)*x11*x21*x31        
     2  +coeff(11)    *x23*x31        
     3  +coeff(12)    *x21*x33        
     4  +coeff(13)    *x21*x32*x41    
     5  +coeff(14)*x11*x22*x31        
     6  +coeff(15)    *x22*x33        
     7  +coeff(16)*x11*x21*x33        
     8  +coeff(17)        *x34*x41    
      p_m35_4_5   =p_m35_4_5   
     9  +coeff(18)*x12    *x31        
     1  +coeff(19)        *x33*x42    
     2  +coeff(20)        *x32*x43    
     3  +coeff(21)*x12    *x33        
     4  +coeff(22)*x11        *x41    
     5  +coeff(23)        *x33    *x51
     6  +coeff(24)    *x21*x31    *x52
     7  +coeff(25)    *x21*x33    *x51
     8  +coeff(26)    *x21*x34*x41    
      p_m35_4_5   =p_m35_4_5   
     9  +coeff(27)    *x22*x33    *x51
     1  +coeff(28)*x12*x22*x31        
     2  +coeff(29)*x11*x22*x33        
     3  +coeff(30)    *x24*x33        
     4  +coeff(31)*x11*x24    *x41    
     5  +coeff(32)    *x22*x34*x41    
     6  +coeff(33)    *x24*x31*x42    
     7  +coeff(34)    *x22*x33*x42    
     8  +coeff(35)*x12*x21*x31    *x51
      p_m35_4_5   =p_m35_4_5   
     9  +coeff(36)*x11*x21*x33    *x51
     1  +coeff(37)*x11*x23    *x41*x51
     2  +coeff(38)*x12    *x31    *x52
     3  +coeff(39)*x12*x23*x31        
     4  +coeff(40)*x12*x22*x31    *x51
     5  +coeff(41)*x11    *x33        
     6  +coeff(42)*x12*x21*x31        
     7  +coeff(43)    *x24*x32*x41    
     8  +coeff(44)*x11    *x32*x43    
      p_m35_4_5   =p_m35_4_5   
     9  +coeff(45)    *x22*x32*x43    
     1  +coeff(46)*x12*x21    *x41*x51
     2  +coeff(47)*x11*x21*x32*x41*x51
     3  +coeff(48)*x11*x21*x31*x42*x51
     4  +coeff(49)*x11*x21    *x43*x51
     5  +coeff(50)    *x24*x31    *x52
c
      return
      end
      function l_m35_4_5   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(25)
      data ncoeff/ 24/
      data avdat/ -0.1691846E-03/
      data xmin/
     1 -0.42252E+00,-0.19647E+00,-0.54023E+00,-0.18368E+00,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.48443E+00, 0.14310E+00, 0.54023E+00, 0.18368E+00, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.66838753E-02,-0.39660152E-02,-0.90027370E-04,-0.11028360E+00,
     + -0.78147030E-03,-0.21140110E-04,-0.31458370E-02, 0.34865222E-02,
     + -0.15326624E-02, 0.14551061E-02,-0.22927524E-02, 0.29177062E-04,
     +  0.20785040E-02,-0.66723110E-02,-0.34406752E-04, 0.47905490E-04,
     +  0.29628840E-02,-0.37781540E-02,-0.46963370E-04,-0.52607490E-04,
     +  0.47850291E-04,-0.22296090E-03,-0.41144452E-02, 0.10189050E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
c
c                  function
c
      l_m35_4_5   =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)                *x51
     4  +coeff( 4)*x11                
     5  +coeff( 5)        *x32        
     6  +coeff( 6)    *x22        *x51
     7  +coeff( 7)    *x22            
     8  +coeff( 8)    *x21*x32        
      l_m35_4_5   =l_m35_4_5   
     9  +coeff( 9)*x11*x22            
     1  +coeff(10)        *x31*x41    
     2  +coeff(11)*x11*x21            
     3  +coeff(12)        *x32    *x51
     4  +coeff(13)*x11    *x32        
     5  +coeff(14)    *x21*x31*x41    
     6  +coeff(15)        *x31*x41*x51
     7  +coeff(16)    *x22*x32        
     8  +coeff(17)    *x21    *x42    
      l_m35_4_5   =l_m35_4_5   
     9  +coeff(18)*x11    *x31*x41    
     1  +coeff(19)    *x23            
     2  +coeff(20)*x11            *x51
     3  +coeff(21)    *x24            
     4  +coeff(22)*x11*x21*x32        
     5  +coeff(23)            *x42    
     6  +coeff(24)    *x21        *x51
c
      return
      end
      function x_m35_5_6   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(30)
      data ncoeff/ 29/
      data avdat/ -0.1230867E-01/
      data xmin/
     1 -0.40971E+00,-0.18752E+00,-0.50997E+00,-0.16300E+00,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.42437E+00, 0.14282E+00, 0.50997E+00, 0.16300E+00, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.37782622E-02, 0.11793090E+00, 0.10002343E-02, 0.31054810E-02,
     +  0.44435053E-02,-0.70236100E-02, 0.37296260E-02, 0.12602803E-02,
     + -0.61310880E-03, 0.42712150E+00, 0.36974910E-04,-0.12432590E-02,
     + -0.94016583E-03, 0.11087450E-02, 0.32630053E-03,-0.12893552E-02,
     +  0.45073880E-03,-0.15877500E-02, 0.50130200E-03,-0.91537960E-03,
     + -0.10912801E-02, 0.69595410E-03, 0.12177751E-02, 0.29116250E-02,
     + -0.13368703E-01, 0.17533523E-01, 0.32425951E-03, 0.58458790E-02,
     + -0.28440854E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x35 = x34*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_m35_5_6   =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)                *x51
     4  +coeff( 4)    *x22            
     5  +coeff( 5)        *x32        
     6  +coeff( 6)        *x31*x41    
     7  +coeff( 7)            *x42    
     8  +coeff( 8)    *x21        *x51
      x_m35_5_6   =x_m35_5_6   
     9  +coeff( 9)                *x52
     1  +coeff(10)*x11                
     2  +coeff(11)    *x23            
     3  +coeff(12)    *x22*x32        
     4  +coeff(13)    *x21*x32        
     5  +coeff(14)    *x22        *x51
     6  +coeff(15)    *x21        *x52
     7  +coeff(16)    *x21*x35*x41*x52
     8  +coeff(17)    *x22        *x52
      x_m35_5_6   =x_m35_5_6   
     9  +coeff(18)*x11    *x32        
     1  +coeff(19)*x11*x21*x32        
     2  +coeff(20)    *x24*x32        
     3  +coeff(21)    *x23*x32    *x51
     4  +coeff(22)    *x21*x31*x41*x53
     5  +coeff(23)*x11    *x34        
     6  +coeff(24)*x11*x22*x31*x41    
     7  +coeff(25)    *x23*x34        
     8  +coeff(26)    *x23*x33*x41    
      x_m35_5_6   =x_m35_5_6   
     9  +coeff(27)*x11*x21*x32*x42    
     1  +coeff(28)    *x24*x33*x41    
     2  +coeff(29)    *x22*x35*x41    
c
      return
      end
      function t_m35_5_6   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(44)
      data ncoeff/ 43/
      data avdat/ -0.3432297E-02/
      data xmin/
     1 -0.40971E+00,-0.18752E+00,-0.50997E+00,-0.16300E+00,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.42437E+00, 0.14282E+00, 0.50997E+00, 0.16300E+00, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.24610092E-01, 0.21657900E+00,-0.89151244E-02, 0.67921320E-01,
     + -0.50230000E-02,-0.36116451E-01, 0.83935880E-01, 0.16617660E-01,
     + -0.58566970E-03,-0.10807590E-01,-0.79304020E-02, 0.95907750E+00,
     + -0.19206881E+01, 0.61520403E+00,-0.96226744E-02,-0.84027390E-01,
     +  0.18361980E+00, 0.58119813E-03,-0.15040750E-01, 0.24163742E+00,
     + -0.43301840E+00,-0.25267544E+00,-0.10772440E-01,-0.14117410E-02,
     + -0.72151840E-01, 0.87090684E-02,-0.99902313E+00, 0.14347040E-01,
     + -0.10137061E+00, 0.16379511E+01,-0.59335730E+00,-0.15478564E+00,
     +  0.16358490E+00, 0.34212800E+00,-0.31152340E+00,-0.14516203E+00,
     +  0.45966690E-01,-0.90049500E-04, 0.23654513E-03,-0.50890883E-02,
     +  0.10286940E-02, 0.15091180E-02, 0.17233180E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_m35_5_6   =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)                *x51
     4  +coeff( 4)*x11                
     5  +coeff( 5)    *x22            
     6  +coeff( 6)        *x32        
     7  +coeff( 7)        *x31*x41    
     8  +coeff( 8)    *x21        *x51
      t_m35_5_6   =t_m35_5_6   
     9  +coeff( 9)                *x52
     1  +coeff(10)*x11*x21            
     2  +coeff(11)    *x23            
     3  +coeff(12)    *x21*x32        
     4  +coeff(13)    *x21*x31*x41    
     5  +coeff(14)    *x21    *x42    
     6  +coeff(15)*x11            *x51
     7  +coeff(16)        *x32    *x51
     8  +coeff(17)        *x31*x41*x51
      t_m35_5_6   =t_m35_5_6   
     9  +coeff(18)                *x53
     1  +coeff(19)*x11    *x32        
     2  +coeff(20)    *x22*x32        
     3  +coeff(21)    *x22*x31*x41    
     4  +coeff(22)*x11        *x42    
     5  +coeff(23)            *x44    
     6  +coeff(24)    *x23        *x51
     7  +coeff(25)    *x21*x32    *x51
     8  +coeff(26)*x11*x21*x32        
      t_m35_5_6   =t_m35_5_6   
     9  +coeff(27)    *x23*x32        
     1  +coeff(28)    *x21*x34        
     2  +coeff(29)*x11*x21*x31*x41    
     3  +coeff(30)    *x23*x31*x41    
     4  +coeff(31)    *x23    *x42    
     5  +coeff(32)*x11    *x32    *x51
     6  +coeff(33)    *x22*x32    *x51
     7  +coeff(34)*x11    *x31*x41*x51
     8  +coeff(35)    *x22*x31*x41*x51
      t_m35_5_6   =t_m35_5_6   
     9  +coeff(36)*x11        *x42*x51
     1  +coeff(37)    *x22    *x42*x51
     2  +coeff(38)        *x32*x42*x51
     3  +coeff(39)        *x31*x43*x51
     4  +coeff(40)            *x44*x51
     5  +coeff(41)*x11*x21        *x52
     6  +coeff(42)    *x21*x32    *x52
     7  +coeff(43)    *x21*x31*x41*x52
c
      return
      end
      function y_m35_5_6   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(21)
      data ncoeff/ 20/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.40971E+00,-0.18752E+00,-0.50997E+00,-0.16300E+00,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.42437E+00, 0.14282E+00, 0.50997E+00, 0.16300E+00, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.49752631E+00, 0.11268080E+00, 0.24609440E-02,-0.19888660E-02,
     +  0.60458760E-03,-0.53795100E-03,-0.40797020E-03,-0.30152772E-02,
     + -0.19604750E-03, 0.23504190E-02, 0.97266080E-04,-0.14882230E-03,
     + -0.91060560E-02,-0.30522371E-03, 0.15837810E-01,-0.48567893E-03,
     + -0.90358760E-03, 0.17213100E-02,-0.31115210E-02, 0.80826120E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      y_m35_5_6   =avdat
     1  +coeff( 1)        *x31        
     2  +coeff( 2)            *x41    
     3  +coeff( 3)        *x31    *x51
     4  +coeff( 4)*x11    *x31        
     5  +coeff( 5)*x11        *x41    
     6  +coeff( 6)        *x32*x41    
     7  +coeff( 7)        *x31*x42    
     8  +coeff( 8)    *x21*x31        
      y_m35_5_6   =y_m35_5_6   
     9  +coeff( 9)    *x22*x31        
     1  +coeff(10)*x11*x21*x31        
     2  +coeff(11)*x11    *x32*x41    
     3  +coeff(12)        *x33    *x52
     4  +coeff(13)*x12*x21*x33        
     5  +coeff(14)        *x31    *x52
     6  +coeff(15)*x12*x21*x34*x41    
     7  +coeff(16)        *x34*x41    
     8  +coeff(17)    *x21*x33    *x51
      y_m35_5_6   =y_m35_5_6   
     9  +coeff(18)*x12*x21*x31        
     1  +coeff(19)    *x21*x34*x41    
     2  +coeff(20)        *x34*x41*x51
c
      return
      end
      function p_m35_5_6   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.40971E+00,-0.18752E+00,-0.50997E+00,-0.16300E+00,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.42437E+00, 0.14282E+00, 0.50997E+00, 0.16300E+00, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.43930480E-01, 0.15634800E+00, 0.11868300E-01,-0.27435524E-02,
     +  0.86771021E-03, 0.18465630E-02, 0.13052402E-01, 0.15100404E-01,
     +  0.52517582E-02, 0.41805062E-03, 0.54158220E-02,-0.27936420E-01,
     +  0.21385742E-01,-0.10269040E-02, 0.24367660E-04, 0.76813073E-02,
     +  0.11952530E-01,-0.65612164E-02,-0.14078770E-02, 0.84483850E-04,
     + -0.97996370E-02,-0.15107020E-01, 0.29901000E-01, 0.10068550E+00,
     + -0.41305440E-01,-0.41070662E-01, 0.11747024E-01,-0.21229050E-02,
     + -0.19223140E-01, 0.10840260E-01,-0.74704892E-01, 0.72660870E-01,
     + -0.19648450E-02, 0.10746400E-01,-0.10887760E-02,-0.26625333E-01,
     + -0.10515640E-01,-0.49579740E-01,-0.63540791E-02, 0.27059360E-01,
     + -0.45071100E-01,-0.11846041E+00,-0.31079740E-01, 0.10122640E-01,
     +  0.26990260E-01,-0.46314573E-02, 0.11427750E-01, 0.70804240E-01,
     + -0.28207791E-02, 0.37677770E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_m35_5_6   =avdat
     1  +coeff( 1)        *x31        
     2  +coeff( 2)            *x41    
     3  +coeff( 3)    *x21*x31        
     4  +coeff( 4)    *x21    *x41    
     5  +coeff( 5)        *x31    *x51
     6  +coeff( 6)            *x41*x51
     7  +coeff( 7)*x11    *x31        
     8  +coeff( 8)    *x22*x31        
      p_m35_5_6   =p_m35_5_6   
     9  +coeff( 9)        *x33        
     1  +coeff(10)*x11        *x41    
     2  +coeff(11)    *x22    *x41    
     3  +coeff(12)        *x32*x41    
     4  +coeff(13)        *x31*x42    
     5  +coeff(14)        *x31    *x52
     6  +coeff(15)            *x41*x52
     7  +coeff(16)*x11*x21*x31        
     8  +coeff(17)    *x22*x31    *x51
      p_m35_5_6   =p_m35_5_6   
     9  +coeff(18)        *x33    *x51
     1  +coeff(19)    *x22    *x41*x51
     2  +coeff(20)        *x31    *x53
     3  +coeff(21)*x12    *x31        
     4  +coeff(22)*x11*x22*x31        
     5  +coeff(23)*x11    *x33        
     6  +coeff(24)    *x22*x33        
     7  +coeff(25)*x11    *x32*x41    
     8  +coeff(26)        *x34*x41    
      p_m35_5_6   =p_m35_5_6   
     9  +coeff(27)*x11*x21*x31*x42    
     1  +coeff(28)    *x21*x33    *x52
     2  +coeff(29)*x11*x24*x31        
     3  +coeff(30)*x12    *x33        
     4  +coeff(31)    *x22*x34*x41    
     5  +coeff(32)    *x22*x33*x42    
     6  +coeff(33)        *x34*x43    
     7  +coeff(34)*x11*x21    *x43*x51
     8  +coeff(35)        *x34*x41*x52
      p_m35_5_6   =p_m35_5_6   
     9  +coeff(36)*x12*x23*x31        
     1  +coeff(37)    *x21*x34*x43    
     2  +coeff(38)*x11*x21*x31*x44    
     3  +coeff(39)    *x23*x31        
     4  +coeff(40)    *x21*x33        
     5  +coeff(41)    *x21*x32*x41    
     6  +coeff(42)    *x22*x32*x41    
     7  +coeff(43)        *x31*x44    
     8  +coeff(44)    *x21*x34*x41    
      p_m35_5_6   =p_m35_5_6   
     9  +coeff(45)*x11*x21*x34*x41    
     1  +coeff(46)*x11*x21*x33    *x52
     2  +coeff(47)        *x32*x41*x51
     3  +coeff(48)        *x33*x42    
     4  +coeff(49)*x11*x21*x31    *x51
     5  +coeff(50)    *x23*x31    *x51
c
      return
      end
      function l_m35_5_6   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(25)
      data ncoeff/ 24/
      data avdat/ -0.2720841E-02/
      data xmin/
     1 -0.40971E+00,-0.18752E+00,-0.50997E+00,-0.16300E+00,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.42437E+00, 0.14282E+00, 0.50997E+00, 0.16300E+00, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.44793770E-02,-0.11636650E-01,-0.10155690E-01,-0.21493640E-02,
     + -0.52384741E-03, 0.60570822E-02,-0.32284634E-03, 0.43380644E-03,
     + -0.66076332E-04, 0.10140810E-03,-0.66812941E-03, 0.20204132E-03,
     +  0.65944580E-05,-0.19804070E-02, 0.66253310E-03,-0.20450504E-03,
     +  0.33384151E-03,-0.11186672E-03,-0.11101300E-01,-0.89432810E-03,
     + -0.38625910E-03,-0.67369511E-03, 0.79025980E-03, 0.84058910E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_m35_5_6   =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)    *x22            
     4  +coeff( 4)        *x32        
     5  +coeff( 5)                *x51
     6  +coeff( 6)        *x31*x41    
     7  +coeff( 7)    *x22        *x51
     8  +coeff( 8)    *x21        *x51
      l_m35_5_6   =l_m35_5_6   
     9  +coeff( 9)    *x23            
     1  +coeff(10)    *x21        *x52
     2  +coeff(11)*x11*x21*x32        
     3  +coeff(12)        *x34*x42    
     4  +coeff(13)*x11                
     5  +coeff(14)*x11*x21            
     6  +coeff(15)        *x32    *x51
     7  +coeff(16)*x12                
     8  +coeff(17)*x11*x22            
      l_m35_5_6   =l_m35_5_6   
     9  +coeff(18)    *x21*x34        
     1  +coeff(19)            *x42    
     2  +coeff(20)        *x31*x41*x51
     3  +coeff(21)*x11    *x34*x42    
     4  +coeff(22)    *x21*x31*x41    
     5  +coeff(23)    *x21    *x42    
     6  +coeff(24)*x11    *x32    *x51
c
      return
      end
      function x_m35_6_7   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(13)
      data ncoeff/ 12/
      data avdat/ -0.2131009E-01/
      data xmin/
     1 -0.33254E+00,-0.17988E+00,-0.55661E+00,-0.96502E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.37600E+00, 0.10672E+00, 0.55661E+00, 0.96502E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.80211674E-02, 0.14450281E+00, 0.46764262E-02,-0.16771400E-02,
     +  0.40303650E-02,-0.59081630E-02,-0.27164920E-03, 0.55385170E-02,
     + -0.39415191E-02, 0.38185720E+00, 0.42379400E-02,-0.18305823E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      x_m35_6_7   =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)                *x51
     4  +coeff( 4)    *x22            
     5  +coeff( 5)        *x32        
     6  +coeff( 6)        *x31*x41    
     7  +coeff( 7)            *x42    
     8  +coeff( 8)    *x21        *x51
      x_m35_6_7   =x_m35_6_7   
     9  +coeff( 9)                *x52
     1  +coeff(10)*x11                
     2  +coeff(11)    *x21        *x52
     3  +coeff(12)    *x22        *x51
c
      return
      end
      function t_m35_6_7   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/ -0.8645021E-02/
      data xmin/
     1 -0.33254E+00,-0.17988E+00,-0.55661E+00,-0.96502E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.37600E+00, 0.10672E+00, 0.55661E+00, 0.96502E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.20324183E-01, 0.12105912E+00, 0.21921563E-01, 0.38440864E-01,
     +  0.18151500E-01, 0.20754700E-01,-0.31152073E-01, 0.75528281E-02,
     + -0.11717983E-01,-0.83109131E-02, 0.17080760E-01,-0.28801440E-01,
     + -0.70569610E+00, 0.13942571E+01,-0.69402420E+00, 0.56968410E-02,
     +  0.24867640E-01, 0.66276760E-02, 0.54364614E-02,-0.90537802E-03,
     + -0.23129410E-01,-0.18825860E-01,-0.22754720E+00, 0.20781510E-02,
     + -0.47363300E-02, 0.16163311E+00, 0.91470871E-02, 0.40114562E-01,
     + -0.14130310E-02,-0.35600871E-02, 0.18464000E-01,-0.51228240E-02,
     + -0.63106352E-02,-0.99815850E-03, 0.73008470E-02, 0.40441720E-02,
     + -0.85618561E-02, 0.30984152E-02, 0.52996650E-02, 0.24194410E-02,
     +  0.19884330E-02,-0.14193534E-02, 0.65831122E-02, 0.28207513E-02,
     +  0.15552822E-02, 0.21717820E-02,-0.15187150E-03, 0.11902000E-02,
     + -0.23541213E-02,-0.11640771E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_m35_6_7   =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)                *x51
     4  +coeff( 4)*x11                
     5  +coeff( 5)    *x22            
     6  +coeff( 6)        *x32        
     7  +coeff( 7)        *x31*x41    
     8  +coeff( 8)            *x42    
      t_m35_6_7   =t_m35_6_7   
     9  +coeff( 9)    *x21        *x51
     1  +coeff(10)                *x52
     2  +coeff(11)*x11*x21            
     3  +coeff(12)    *x23            
     4  +coeff(13)    *x21*x32        
     5  +coeff(14)    *x21*x31*x41    
     6  +coeff(15)    *x21    *x42    
     7  +coeff(16)*x11            *x51
     8  +coeff(17)    *x22        *x51
      t_m35_6_7   =t_m35_6_7   
     9  +coeff(18)        *x32    *x51
     1  +coeff(19)    *x21        *x52
     2  +coeff(20)                *x53
     3  +coeff(21)*x11*x22            
     4  +coeff(22)*x11    *x32        
     5  +coeff(23)    *x22*x32        
     6  +coeff(24)    *x22*x32*x42    
     7  +coeff(25)    *x24            
     8  +coeff(26)    *x22*x31*x41    
      t_m35_6_7   =t_m35_6_7   
     9  +coeff(27)*x11        *x42    
     1  +coeff(28)    *x22    *x42    
     2  +coeff(29)    *x21        *x53
     3  +coeff(30)    *x23    *x42    
     4  +coeff(31)*x11*x22        *x51
     5  +coeff(32)    *x22*x32    *x51
     6  +coeff(33)    *x22*x31*x41*x51
     7  +coeff(34)*x11*x21        *x52
     8  +coeff(35)    *x23        *x52
      t_m35_6_7   =t_m35_6_7   
     9  +coeff(36)    *x22        *x53
     1  +coeff(37)*x11*x24            
     2  +coeff(38)    *x24*x32        
     3  +coeff(39)    *x22*x34        
     4  +coeff(40)*x11*x22*x31*x41    
     5  +coeff(41)    *x24*x31*x41    
     6  +coeff(42)*x11    *x33*x41    
     7  +coeff(43)    *x22*x33*x41    
     8  +coeff(44)*x11*x22    *x42    
      t_m35_6_7   =t_m35_6_7   
     9  +coeff(45)    *x24    *x42    
     1  +coeff(46)*x11    *x32*x42    
     2  +coeff(47)        *x34*x42    
     3  +coeff(48)*x11    *x31*x43    
     4  +coeff(49)    *x22*x31*x43    
     5  +coeff(50)        *x33*x43    
c
      return
      end
      function y_m35_6_7   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.33254E+00,-0.17988E+00,-0.55661E+00,-0.96502E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.37600E+00, 0.10672E+00, 0.55661E+00, 0.96502E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.51533430E+00, 0.93349680E-01,-0.25948450E-01,-0.83877830E-04,
     +  0.18526641E-01,-0.15062254E-01,-0.70985460E-03,-0.19287444E-02,
     + -0.18634623E-02, 0.90807834E-02,-0.14197263E+00,-0.14554694E-01,
     + -0.10728340E-01, 0.26152380E-01, 0.36181311E-02, 0.38511880E-02,
     +  0.19124560E-02,-0.21350143E-03, 0.79890720E-04, 0.44137440E-02,
     +  0.30405740E-02, 0.28418860E+00,-0.13624450E+00, 0.22039194E-02,
     + -0.39486330E-03, 0.73584610E-03, 0.13931320E-02,-0.37069520E-02,
     +  0.16018480E-01,-0.14758041E-01, 0.15929561E-02,-0.36196850E-03,
     +  0.15426764E-03, 0.40613594E-02,-0.37376730E-03, 0.15375473E-02,
     +  0.55943572E+00,-0.21715180E-01,-0.10318770E+01, 0.37122420E+00,
     +  0.97566634E-01, 0.19796350E-01, 0.18958362E-02,-0.19599360E-01,
     +  0.22149700E-02, 0.21298040E-01,-0.10937454E-01,-0.34626610E-03,
     + -0.43540261E-02, 0.47304780E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      y_m35_6_7   =avdat
     1  +coeff( 1)        *x31        
     2  +coeff( 2)            *x41    
     3  +coeff( 3)    *x21*x31        
     4  +coeff( 4)    *x21    *x41    
     5  +coeff( 5)        *x31    *x51
     6  +coeff( 6)*x11    *x31        
     7  +coeff( 7)        *x31    *x52
     8  +coeff( 8)        *x33        
      y_m35_6_7   =y_m35_6_7   
     9  +coeff( 9)    *x22    *x41    
     1  +coeff(10)*x11*x21*x31        
     2  +coeff(11)*x11    *x33        
     3  +coeff(12)    *x21*x31    *x51
     4  +coeff(13)*x11    *x34*x41    
     5  +coeff(14)*x11    *x34*x43    
     6  +coeff(15)    *x23*x31        
     7  +coeff(16)    *x21*x33        
     8  +coeff(17)    *x22    *x41*x51
      y_m35_6_7   =y_m35_6_7   
     9  +coeff(18)        *x31    *x53
     1  +coeff(19)            *x41*x51
     2  +coeff(20)    *x22*x31        
     3  +coeff(21)*x11        *x41*x51
     4  +coeff(22)*x11    *x32*x41    
     5  +coeff(23)*x11    *x31*x42    
     6  +coeff(24)    *x23*x31    *x51
     7  +coeff(25)    *x22*x31    *x52
     8  +coeff(26)    *x21*x31    *x53
      y_m35_6_7   =y_m35_6_7   
     9  +coeff(27)*x12*x21*x31        
     1  +coeff(28)    *x21*x34*x41    
     2  +coeff(29)    *x23*x31*x42    
     3  +coeff(30)    *x23    *x43    
     4  +coeff(31)*x11*x22*x31    *x51
     5  +coeff(32)    *x23    *x41*x52
     6  +coeff(33)        *x33*x44    
     7  +coeff(34)*x11*x23*x33        
     8  +coeff(35)    *x23*x34*x41    
      y_m35_6_7   =y_m35_6_7   
     9  +coeff(36)    *x21*x34*x43    
     1  +coeff(37)*x12    *x33    *x51
     2  +coeff(38)*x11*x22*x33    *x51
     3  +coeff(39)*x12    *x32*x41*x51
     4  +coeff(40)*x12    *x31*x42*x51
     5  +coeff(41)*x12        *x43*x51
     6  +coeff(42)*x11*x22    *x43*x51
     7  +coeff(43)*x12*x24    *x41    
     8  +coeff(44)*x11    *x33*x44    
      y_m35_6_7   =y_m35_6_7   
     9  +coeff(45)*x12*x21*x33    *x51
     1  +coeff(46)*x12*x21*x32*x41*x51
     2  +coeff(47)*x12*x21*x31*x42*x51
     3  +coeff(48)*x11*x23*x31*x42*x51
     4  +coeff(49)*x12*x21    *x43*x51
     5  +coeff(50)*x11*x23    *x43*x51
c
      return
      end
      function p_m35_6_7   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.33254E+00,-0.17988E+00,-0.55661E+00,-0.96502E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.37600E+00, 0.10672E+00, 0.55661E+00, 0.96502E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.94551790E-01, 0.92417010E-01,-0.97746170E-02,-0.60271512E-03,
     +  0.14771200E-01, 0.29904404E-02,-0.69160540E-02, 0.42911432E-02,
     +  0.10165740E+01, 0.73794582E-02, 0.31846810E-02,-0.32889070E+01,
     +  0.35473372E+01,-0.12768970E+01,-0.84566710E-02, 0.46669380E-03,
     + -0.11212914E-02, 0.40663131E-02,-0.31876452E-02, 0.43239650E+00,
     + -0.97373150E+00, 0.56481081E+00,-0.23978840E-01,-0.11223141E-01,
     +  0.50907813E-01, 0.17431540E-02, 0.15525740E-01,-0.66436360E-01,
     +  0.14899890E-01, 0.15347291E-02,-0.82894082E-03,-0.83688413E-03,
     +  0.65462650E-02,-0.52145490E-02, 0.58484752E-02, 0.55005750E-02,
     +  0.14657890E-02,-0.98518021E-02,-0.12057550E-01,-0.20303730E-02,
     + -0.10360962E-02, 0.25939830E-02,-0.37767680E-02,-0.74415053E-02,
     +  0.28697340E-02,-0.19557750E-02,-0.10967440E-04, 0.32197963E-03,
     + -0.47415420E-02,-0.22938852E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_m35_6_7   =avdat
     1  +coeff( 1)        *x31        
     2  +coeff( 2)            *x41    
     3  +coeff( 3)    *x21*x31        
     4  +coeff( 4)    *x21    *x41    
     5  +coeff( 5)        *x31    *x51
     6  +coeff( 6)            *x41*x51
     7  +coeff( 7)*x11    *x31        
     8  +coeff( 8)    *x22*x31        
      p_m35_6_7   =p_m35_6_7   
     9  +coeff( 9)        *x33        
     1  +coeff(10)*x11        *x41    
     2  +coeff(11)    *x22    *x41    
     3  +coeff(12)        *x32*x41    
     4  +coeff(13)        *x31*x42    
     5  +coeff(14)            *x43    
     6  +coeff(15)    *x21*x31    *x51
     7  +coeff(16)    *x21    *x41*x51
     8  +coeff(17)        *x31    *x52
      p_m35_6_7   =p_m35_6_7   
     9  +coeff(18)*x11*x21*x31        
     1  +coeff(19)    *x23*x31        
     2  +coeff(20)    *x21*x33        
     3  +coeff(21)    *x21*x32*x41    
     4  +coeff(22)    *x21*x31*x42    
     5  +coeff(23)    *x21    *x43    
     6  +coeff(24)    *x22*x31    *x51
     7  +coeff(25)        *x33    *x51
     8  +coeff(26)*x11        *x41*x51
      p_m35_6_7   =p_m35_6_7   
     9  +coeff(27)    *x22    *x41*x51
     1  +coeff(28)        *x32*x41*x51
     2  +coeff(29)        *x31*x42*x51
     3  +coeff(30)            *x43*x51
     4  +coeff(31)    *x21*x31    *x52
     5  +coeff(32)    *x21    *x41*x52
     6  +coeff(33)        *x31    *x53
     7  +coeff(34)            *x41*x53
     8  +coeff(35)*x12    *x31        
      p_m35_6_7   =p_m35_6_7   
     9  +coeff(36)*x11*x22*x31        
     1  +coeff(37)    *x22*x33        
     2  +coeff(38)*x12        *x41    
     3  +coeff(39)*x11*x22    *x41    
     4  +coeff(40)    *x24    *x41    
     5  +coeff(41)*x11    *x32*x41    
     6  +coeff(42)    *x22*x32*x41    
     7  +coeff(43)*x11    *x31*x42    
     8  +coeff(44)    *x22*x31*x42    
      p_m35_6_7   =p_m35_6_7   
     9  +coeff(45)*x11        *x43    
     1  +coeff(46)    *x22    *x43    
     2  +coeff(47)        *x32*x43    
     3  +coeff(48)        *x31*x44    
     4  +coeff(49)*x11*x21*x31    *x51
     5  +coeff(50)    *x23*x31    *x51
c
      return
      end
      function l_m35_6_7   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(25)
      data ncoeff/ 24/
      data avdat/ -0.1216648E-02/
      data xmin/
     1 -0.33254E+00,-0.17988E+00,-0.55661E+00,-0.96502E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.37600E+00, 0.10672E+00, 0.55661E+00, 0.96502E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.37421080E-02,-0.77913770E-02,-0.10287132E-01,-0.31850060E-03,
     + -0.38388310E-02,-0.66390890E-03,-0.12187020E-02,-0.63794070E-03,
     + -0.39563880E-03, 0.29250140E-02, 0.38387250E-03, 0.37655961E-03,
     +  0.11294234E-03, 0.11974840E-03, 0.19306810E-03, 0.58876770E-04,
     + -0.46168183E-03,-0.55413973E-04,-0.32657210E-03,-0.37115920E-02,
     + -0.59074984E-03, 0.48582430E-03,-0.74015200E-04,-0.19153800E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      l_m35_6_7   =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)    *x22            
     4  +coeff( 4)                *x51
     5  +coeff( 5)*x11*x21            
     6  +coeff( 6)        *x32        
     7  +coeff( 7)*x11                
     8  +coeff( 8)    *x21        *x51
      l_m35_6_7   =l_m35_6_7   
     9  +coeff( 9)*x12                
     1  +coeff(10)        *x31*x41    
     2  +coeff(11)    *x21*x32        
     3  +coeff(12)        *x32    *x51
     4  +coeff(13)    *x21        *x52
     5  +coeff(14)                *x53
     6  +coeff(15)*x11*x21        *x51
     7  +coeff(16)*x12*x21            
     8  +coeff(17)    *x22*x32        
      l_m35_6_7   =l_m35_6_7   
     9  +coeff(18)*x11            *x52
     1  +coeff(19)*x11*x21*x32        
     2  +coeff(20)            *x42    
     3  +coeff(21)        *x31*x41*x51
     4  +coeff(22)    *x22    *x42    
     5  +coeff(23)    *x22        *x52
     6  +coeff(24)    *x21        *x53
c
      return
      end
      function x_m35_7_8   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/ -0.2741403E-01/
      data xmin/
     1 -0.27559E+00,-0.19762E+00,-0.59695E+00,-0.27135E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.37594E+00, 0.49014E-01, 0.59695E+00, 0.27135E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.41537303E-01, 0.61444044E-01, 0.16322651E-02,-0.28744550E-03,
     + -0.80649540E-03,-0.17658350E-04,-0.85364280E-04,-0.56961120E-03,
     + -0.26763960E-04, 0.33253820E+00,-0.61075790E-03,-0.93120600E-03,
     + -0.31383980E-04,-0.49931580E-03, 0.23472492E-03, 0.17505184E-03,
     + -0.98589160E-03, 0.15491551E-03, 0.32878310E-05, 0.69386942E-03,
     +  0.15013120E-03, 0.59450574E-03,-0.25204950E-03,-0.30143640E-03,
     + -0.96109492E-04,-0.52499293E-04, 0.21225680E-02,-0.20953624E-02,
     +  0.51298300E-03,-0.36456490E-03, 0.45208414E-03,-0.90280680E-03,
     +  0.10976551E-02,-0.79742560E-03, 0.70442631E-03, 0.65060163E-03,
     +  0.41483121E-03,-0.19371422E-03,-0.31197200E-04, 0.41986364E-03,
     +  0.45751011E-03,-0.30256414E-03,-0.57377410E-03,-0.15986261E-04,
     +  0.45027583E-03,-0.33368682E-03, 0.10464010E-02, 0.18469060E-03,
     + -0.28634060E-03, 0.45227643E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x35 = x34*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      x_m35_7_8   =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)                *x51
     4  +coeff( 4)    *x22            
     5  +coeff( 5)        *x32        
     6  +coeff( 6)        *x31*x41    
     7  +coeff( 7)            *x42    
     8  +coeff( 8)    *x21        *x51
      x_m35_7_8   =x_m35_7_8   
     9  +coeff( 9)                *x52
     1  +coeff(10)*x11                
     2  +coeff(11)*x11            *x51
     3  +coeff(12)    *x22        *x51
     4  +coeff(13)    *x21        *x52
     5  +coeff(14)*x12                
     6  +coeff(15)        *x32    *x51
     7  +coeff(16)*x11            *x52
     8  +coeff(17)*x11*x22        *x51
      x_m35_7_8   =x_m35_7_8   
     9  +coeff(18)    *x21    *x42    
     1  +coeff(19)    *x21*x34        
     2  +coeff(20)*x11    *x32    *x51
     3  +coeff(21)*x11*x21        *x52
     4  +coeff(22)*x11*x24            
     5  +coeff(23)*x11    *x34        
     6  +coeff(24)*x11    *x32    *x52
     7  +coeff(25)*x11*x21        *x53
     8  +coeff(26)*x12    *x31*x41    
      x_m35_7_8   =x_m35_7_8   
     9  +coeff(27)*x12    *x32    *x51
     1  +coeff(28)*x12*x21*x32    *x51
     2  +coeff(29)*x12*x21*x31*x41*x51
     3  +coeff(30)    *x23    *x42*x54
     4  +coeff(31)*x11*x21*x34    *x52
     5  +coeff(32)*x12*x21*x34        
     6  +coeff(33)*x12    *x34    *x51
     7  +coeff(34)*x12*x21*x32    *x52
     8  +coeff(35)*x12    *x34    *x52
      x_m35_7_8   =x_m35_7_8   
     9  +coeff(36)    *x23            
     1  +coeff(37)    *x24            
     2  +coeff(38)        *x34        
     3  +coeff(39)        *x33*x41    
     4  +coeff(40)        *x31*x43    
     5  +coeff(41)    *x23        *x51
     6  +coeff(42)    *x21*x31*x41*x51
     7  +coeff(43)    *x21*x32*x42    
     8  +coeff(44)*x11*x21*x32        
      x_m35_7_8   =x_m35_7_8   
     9  +coeff(45)*x12*x21            
     1  +coeff(46)        *x35*x41    
     2  +coeff(47)    *x21*x33*x41*x51
     3  +coeff(48)            *x42*x54
     4  +coeff(49)*x11*x22*x32        
     5  +coeff(50)*x11    *x33*x41    
c
      return
      end
      function t_m35_7_8   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/ -0.1253666E-01/
      data xmin/
     1 -0.27559E+00,-0.19762E+00,-0.59695E+00,-0.27135E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.37594E+00, 0.49014E-01, 0.59695E+00, 0.27135E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.60131520E-01, 0.12791700E+00, 0.21699702E-03, 0.28494230E-01,
     +  0.75288041E-03,-0.28230700E-02, 0.11162620E-03, 0.19117271E-03,
     +  0.37922482E-02,-0.38558750E-02,-0.13748153E-02,-0.40392822E-03,
     + -0.47593251E-02, 0.19801010E-02, 0.83202990E-03, 0.24674390E-02,
     + -0.70243160E-04, 0.56406721E-03,-0.49321730E-03,-0.14690133E-02,
     +  0.72895420E-02, 0.38529622E-04, 0.47674362E-04,-0.13232670E-01,
     +  0.12296260E-02,-0.57568342E-03,-0.91711041E-03,-0.34468674E-04,
     + -0.17254860E-02,-0.52752410E-04, 0.12939120E-03, 0.77355520E-03,
     +  0.19785380E-02,-0.64605171E-02, 0.62779393E-02, 0.71936440E-02,
     +  0.13436980E-02, 0.72888780E-04,-0.32129262E-02, 0.41506720E-02,
     + -0.12690582E-02, 0.73385174E-03,-0.59659061E-02,-0.29127963E-02,
     +  0.13900981E-02, 0.12427860E-02, 0.10888240E-02,-0.18812700E-03,
     +  0.20318490E-02, 0.20199434E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_m35_7_8   =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)                *x51
     4  +coeff( 4)*x11                
     5  +coeff( 5)    *x22            
     6  +coeff( 6)        *x32        
     7  +coeff( 7)        *x31*x41    
     8  +coeff( 8)            *x42    
      t_m35_7_8   =t_m35_7_8   
     9  +coeff( 9)    *x21        *x51
     1  +coeff(10)                *x52
     2  +coeff(11)    *x23            
     3  +coeff(12)    *x21*x31*x41    
     4  +coeff(13)*x11            *x51
     5  +coeff(14)    *x22        *x51
     6  +coeff(15)        *x32    *x51
     7  +coeff(16)    *x21        *x52
     8  +coeff(17)*x12                
      t_m35_7_8   =t_m35_7_8   
     9  +coeff(18)*x11    *x32        
     1  +coeff(19)        *x34        
     2  +coeff(20)    *x23        *x51
     3  +coeff(21)*x12    *x32        
     4  +coeff(22)    *x22*x34        
     5  +coeff(23)*x12        *x42    
     6  +coeff(24)*x12*x21*x32        
     7  +coeff(25)    *x21*x31*x41*x51
     8  +coeff(26)    *x21*x34        
      t_m35_7_8   =t_m35_7_8   
     9  +coeff(27)*x12            *x51
     1  +coeff(28)*x12*x22            
     2  +coeff(29)*x11*x22    *x42    
     3  +coeff(30)    *x22*x32*x42    
     4  +coeff(31)    *x21*x34    *x51
     5  +coeff(32)    *x21*x33*x41*x51
     6  +coeff(33)*x11    *x32    *x52
     7  +coeff(34)*x11*x23*x32        
     8  +coeff(35)*x12    *x32    *x51
      t_m35_7_8   =t_m35_7_8   
     9  +coeff(36)*x11*x22*x32    *x51
     1  +coeff(37)*x11    *x34    *x51
     2  +coeff(38)    *x22*x34    *x51
     3  +coeff(39)    *x24*x31*x41*x51
     4  +coeff(40)*x11    *x33*x41*x51
     5  +coeff(41)        *x34*x42*x51
     6  +coeff(42)*x12*x21        *x52
     7  +coeff(43)*x11*x21*x32    *x52
     8  +coeff(44)    *x23*x32    *x52
      t_m35_7_8   =t_m35_7_8   
     9  +coeff(45)    *x21*x34    *x52
     1  +coeff(46)    *x23*x31*x41*x52
     2  +coeff(47)*x11*x21    *x42*x52
     3  +coeff(48)*x12            *x53
     4  +coeff(49)*x11    *x32    *x53
     5  +coeff(50)    *x22*x32    *x53
c
      return
      end
      function y_m35_7_8   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(12)
      data ncoeff/ 11/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.27559E+00,-0.19762E+00,-0.59695E+00,-0.27135E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.37594E+00, 0.49014E-01, 0.59695E+00, 0.27135E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.58441690E+00, 0.13537170E-01, 0.19556980E-02, 0.45284780E-03,
     + -0.12319160E-03,-0.41914560E-03, 0.18437560E-03, 0.17430630E-03,
     +  0.12591771E-02,-0.15269000E-03, 0.70030990E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      y_m35_7_8   =avdat
     1  +coeff( 1)        *x31        
     2  +coeff( 2)            *x41    
     3  +coeff( 3)        *x31    *x51
     4  +coeff( 4)*x11    *x31        
     5  +coeff( 5)    *x21    *x41    
     6  +coeff( 6)        *x31    *x52
     7  +coeff( 7)        *x31*x42*x51
     8  +coeff( 8)    *x21*x31    *x52
      y_m35_7_8   =y_m35_7_8   
     9  +coeff( 9)*x12    *x34*x41    
     1  +coeff(10)        *x33    *x51
     2  +coeff(11)*x12    *x31*x42    
c
      return
      end
      function p_m35_7_8   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.27559E+00,-0.19762E+00,-0.59695E+00,-0.27135E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.37594E+00, 0.49014E-01, 0.59695E+00, 0.27135E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.54784782E-01, 0.26324311E-01, 0.86593320E-02, 0.39589434E-03,
     + -0.57811250E-03, 0.35254261E-03, 0.10898100E-01, 0.40027634E-02,
     + -0.45330751E-04, 0.27486260E-02, 0.51639470E-02, 0.37666410E-02,
     +  0.44217531E-03,-0.73815160E-02, 0.10177291E-01, 0.69912990E-02,
     + -0.30762160E-02, 0.81403480E-03, 0.60830560E-02,-0.38376471E-02,
     + -0.53458503E-03,-0.32671610E-02, 0.13861680E-02, 0.70750461E-02,
     + -0.63219210E-02,-0.23493340E-02, 0.51977970E-02, 0.34682143E-02,
     +  0.88303030E-03,-0.81241613E-03,-0.18417850E-02,-0.34279810E-02,
     + -0.20679990E-02,-0.36521123E-02, 0.26328460E-02,-0.55338830E-02,
     +  0.42527700E-02,-0.16335530E-03,-0.43235940E-03,-0.11919641E-02,
     +  0.18646433E-02, 0.67080610E-03,-0.30989340E-02,-0.41783060E-02,
     + -0.10732510E-01,-0.34347882E-02,-0.71639631E-03,-0.46349610E-03,
     + -0.16385090E-02, 0.11377941E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_m35_7_8   =avdat
     1  +coeff( 1)        *x31        
     2  +coeff( 2)            *x41    
     3  +coeff( 3)    *x21*x31        
     4  +coeff( 4)    *x21    *x41    
     5  +coeff( 5)        *x31    *x51
     6  +coeff( 6)            *x41*x51
     7  +coeff( 7)*x11    *x31        
     8  +coeff( 8)    *x22*x31        
      p_m35_7_8   =p_m35_7_8   
     9  +coeff( 9)        *x33        
     1  +coeff(10)            *x43    
     2  +coeff(11)    *x21*x31    *x51
     3  +coeff(12)        *x31    *x52
     4  +coeff(13)            *x41*x52
     5  +coeff(14)*x11*x21*x31    *x51
     6  +coeff(15)*x11    *x31    *x52
     7  +coeff(16)*x11    *x33*x42    
     8  +coeff(17)*x11*x22*x31    *x52
      p_m35_7_8   =p_m35_7_8   
     9  +coeff(18)*x12*x21*x33        
     1  +coeff(19)*x12    *x32*x41*x51
     2  +coeff(20)        *x33    *x51
     3  +coeff(21)    *x22    *x41*x51
     4  +coeff(22)    *x24*x31        
     5  +coeff(23)*x11    *x33        
     6  +coeff(24)*x11        *x43    
     7  +coeff(25)        *x32*x43    
     8  +coeff(26)    *x22*x31    *x52
      p_m35_7_8   =p_m35_7_8   
     9  +coeff(27)        *x32*x41*x52
     1  +coeff(28)    *x23*x33        
     2  +coeff(29)    *x23*x31*x42    
     3  +coeff(30)*x11    *x33    *x51
     4  +coeff(31)    *x22    *x43*x51
     5  +coeff(32)*x12*x22*x31        
     6  +coeff(33)*x11*x22*x31*x42    
     7  +coeff(34)        *x33*x44    
     8  +coeff(35)*x11*x21*x32*x41*x51
      p_m35_7_8   =p_m35_7_8   
     9  +coeff(36)    *x23*x32*x41*x51
     1  +coeff(37)*x12    *x31    *x52
     2  +coeff(38)*x12        *x41*x52
     3  +coeff(39)    *x23*x31    *x53
     4  +coeff(40)*x12*x22*x31    *x51
     5  +coeff(41)*x11*x24*x31    *x51
     6  +coeff(42)*x12    *x33    *x51
     7  +coeff(43)*x11*x22*x33    *x51
     8  +coeff(44)*x11*x21*x31*x42*x52
      p_m35_7_8   =p_m35_7_8   
     9  +coeff(45)*x11*x21*x31        
     1  +coeff(46)    *x23*x31        
     2  +coeff(47)    *x21*x33        
     3  +coeff(48)    *x23    *x41    
     4  +coeff(49)    *x21*x31*x42    
     5  +coeff(50)    *x22*x33        
c
      return
      end
      function l_m35_7_8   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(25)
      data ncoeff/ 24/
      data avdat/ -0.5040885E-03/
      data xmin/
     1 -0.27559E+00,-0.19762E+00,-0.59695E+00,-0.27135E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.37594E+00, 0.49014E-01, 0.59695E+00, 0.27135E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.72592881E-03, 0.36915920E-04,-0.37642670E-02, 0.45140733E-03,
     + -0.99326240E-04,-0.18979130E-03,-0.81366960E-03, 0.54113920E-04,
     +  0.18484300E-02,-0.20555172E-03, 0.50416554E-04, 0.23106380E-03,
     + -0.11559393E-03, 0.83325292E-04,-0.40755843E-04,-0.19246972E-03,
     +  0.55107140E-04,-0.15148571E-04,-0.77167132E-04, 0.16519330E-03,
     + -0.21028730E-04,-0.33622010E-04,-0.88306071E-04, 0.37320330E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_m35_7_8   =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)                *x51
     3  +coeff( 3)    *x22            
     4  +coeff( 4)        *x31*x41    
     5  +coeff( 5)            *x42    
     6  +coeff( 6)    *x21        *x51
     7  +coeff( 7)*x11*x21            
     8  +coeff( 8)    *x22        *x51
      l_m35_7_8   =l_m35_7_8   
     9  +coeff( 9)    *x21            
     1  +coeff(10)        *x32        
     2  +coeff(11)    *x21*x32        
     3  +coeff(12)*x11                
     4  +coeff(13)*x11            *x51
     5  +coeff(14)        *x32    *x51
     6  +coeff(15)*x12                
     7  +coeff(16)    *x21*x31*x41    
     8  +coeff(17)*x11    *x32    *x52
      l_m35_7_8   =l_m35_7_8   
     9  +coeff(18)    *x22*x32        
     1  +coeff(19)*x11    *x31*x41    
     2  +coeff(20)*x11*x21        *x51
     3  +coeff(21)*x11*x21        *x52
     4  +coeff(22)    *x23            
     5  +coeff(23)    *x21    *x42    
     6  +coeff(24)            *x42*x51
c
      return
      end
      function x_m35_8_9   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(13)
      data ncoeff/ 12/
      data avdat/ -0.3512337E-01/
      data xmin/
     1 -0.26329E+00,-0.22254E+00,-0.58060E+00,-0.86625E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.40506E+00, 0.80612E-01, 0.58060E+00, 0.86625E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.72223380E-01, 0.75367413E-01, 0.20055550E-02,-0.45287400E-03,
     + -0.10616034E-02,-0.21644450E-02,-0.16397942E-02,-0.74083820E-03,
     + -0.17575700E-03, 0.34074800E+00,-0.66156960E-03, 0.57017010E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      x_m35_8_9   =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)                *x51
     4  +coeff( 4)    *x22            
     5  +coeff( 5)        *x32        
     6  +coeff( 6)        *x31*x41    
     7  +coeff( 7)            *x42    
     8  +coeff( 8)    *x21        *x51
      x_m35_8_9   =x_m35_8_9   
     9  +coeff( 9)                *x52
     1  +coeff(10)*x11                
     2  +coeff(11)*x12                
     3  +coeff(12)*x12*x21            
c
      return
      end
      function t_m35_8_9   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/ -0.1642303E-01/
      data xmin/
     1 -0.26329E+00,-0.22254E+00,-0.58060E+00,-0.86625E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.40506E+00, 0.80612E-01, 0.58060E+00, 0.86625E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.49812780E-01, 0.15560334E+00, 0.30256790E-02, 0.27424840E-01,
     + -0.41340910E-02,-0.22207540E-02, 0.12557851E-02, 0.17342483E-02,
     +  0.33324323E-02,-0.19298100E-02,-0.15083520E-02,-0.72411350E-04,
     + -0.29493034E-02,-0.51252022E-02,-0.39024630E-02, 0.46079380E-02,
     +  0.13372474E-02,-0.17828030E-02,-0.55990261E-02,-0.24743950E-02,
     +  0.17548202E-02,-0.14549372E-02, 0.53959270E-02, 0.10478470E-02,
     + -0.15408460E-02,-0.18571832E-02, 0.42523513E-02, 0.11060650E-02,
     +  0.68263510E-04, 0.33081060E-02,-0.16122360E-02,-0.11054780E-02,
     + -0.27666762E-02,-0.11932770E-02, 0.21690920E-02, 0.19816390E-02,
     +  0.23501543E-02, 0.15350744E-02, 0.68252190E-03,-0.44416300E-02,
     +  0.31026680E-02, 0.12981613E-02, 0.27780702E-02,-0.16144780E-02,
     + -0.91727220E-03,-0.14821810E-02, 0.24729373E-02,-0.16300202E-02,
     + -0.34369190E-02,-0.16622380E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_m35_8_9   =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)                *x51
     4  +coeff( 4)*x11                
     5  +coeff( 5)    *x22            
     6  +coeff( 6)        *x32        
     7  +coeff( 7)        *x31*x41    
     8  +coeff( 8)            *x42    
      t_m35_8_9   =t_m35_8_9   
     9  +coeff( 9)    *x21        *x51
     1  +coeff(10)                *x52
     2  +coeff(11)*x11*x21            
     3  +coeff(12)    *x23            
     4  +coeff(13)*x11            *x51
     5  +coeff(14)    *x22        *x51
     6  +coeff(15)        *x32    *x51
     7  +coeff(16)    *x22*x32        
     8  +coeff(17)        *x34        
      t_m35_8_9   =t_m35_8_9   
     9  +coeff(18)*x11    *x31*x41    
     1  +coeff(19)*x11*x21*x32        
     2  +coeff(20)*x12    *x32        
     3  +coeff(21)*x11    *x34        
     4  +coeff(22)*x12    *x31*x41    
     5  +coeff(23)*x12    *x34        
     6  +coeff(24)    *x21        *x52
     7  +coeff(25)        *x31*x41*x51
     8  +coeff(26)*x12                
      t_m35_8_9   =t_m35_8_9   
     9  +coeff(27)*x11*x22            
     1  +coeff(28)*x11    *x32        
     2  +coeff(29)    *x22        *x52
     3  +coeff(30)*x12*x21            
     4  +coeff(31)        *x34    *x51
     5  +coeff(32)*x12            *x52
     6  +coeff(33)*x12*x21*x32        
     7  +coeff(34)    *x23*x34        
     8  +coeff(35)*x12    *x32    *x51
      t_m35_8_9   =t_m35_8_9   
     9  +coeff(36)    *x24*x32    *x51
     1  +coeff(37)    *x22*x34    *x51
     2  +coeff(38)*x11    *x33*x41*x51
     3  +coeff(39)    *x22*x31*x41*x53
     4  +coeff(40)*x12*x22*x32        
     5  +coeff(41)    *x24*x34        
     6  +coeff(42)    *x21*x32        
     7  +coeff(43)    *x24            
     8  +coeff(44)*x11        *x42    
      t_m35_8_9   =t_m35_8_9   
     9  +coeff(45)    *x21*x31*x41*x51
     1  +coeff(46)*x11            *x52
     2  +coeff(47)    *x24        *x51
     3  +coeff(48)*x11*x24            
     4  +coeff(49)*x11    *x33*x41    
     5  +coeff(50)*x11*x21*x32    *x51
c
      return
      end
      function y_m35_8_9   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(31)
      data ncoeff/ 30/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.26329E+00,-0.22254E+00,-0.58060E+00,-0.86625E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.40506E+00, 0.80612E-01, 0.58060E+00, 0.86625E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.56843030E+00, 0.42910170E-01, 0.20196272E-02, 0.52553263E-03,
     + -0.13259540E-03, 0.90995920E-04, 0.16679870E-03, 0.26383690E-03,
     + -0.51896670E-03, 0.20759490E-02,-0.20196870E-03, 0.15431181E-04,
     +  0.38012261E-04, 0.11677380E-02, 0.56810950E-03,-0.16039780E-03,
     + -0.30527210E-03,-0.19909381E-03, 0.41157842E-03, 0.14974932E-03,
     + -0.28858330E-05,-0.37511520E-03,-0.13678020E-03, 0.56646680E-04,
     + -0.27829280E-03, 0.27173131E-02, 0.63394260E-03, 0.58414070E-03,
     +  0.16173394E-02, 0.19425670E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      y_m35_8_9   =avdat
     1  +coeff( 1)        *x31        
     2  +coeff( 2)            *x41    
     3  +coeff( 3)        *x31    *x51
     4  +coeff( 4)    *x21*x31        
     5  +coeff( 5)            *x41*x51
     6  +coeff( 6)*x11    *x31        
     7  +coeff( 7)    *x21    *x41    
     8  +coeff( 8)        *x33        
      y_m35_8_9   =y_m35_8_9   
     9  +coeff( 9)    *x21*x31    *x51
     1  +coeff(10)*x11    *x33        
     2  +coeff(11)    *x22*x31        
     3  +coeff(12)*x11        *x41    
     4  +coeff(13)    *x23*x33        
     5  +coeff(14)*x12    *x33        
     6  +coeff(15)*x12*x22    *x41    
     7  +coeff(16)*x12    *x31    *x52
     8  +coeff(17)*x11*x23*x33        
      y_m35_8_9   =y_m35_8_9   
     9  +coeff(18)    *x22    *x41    
     1  +coeff(19)        *x32*x41    
     2  +coeff(20)            *x41*x52
     3  +coeff(21)*x11    *x31    *x51
     4  +coeff(22)*x12    *x31        
     5  +coeff(23)    *x22*x33        
     6  +coeff(24)    *x24    *x41    
     7  +coeff(25)*x12    *x32*x41    
     8  +coeff(26)*x11    *x34*x41    
      y_m35_8_9   =y_m35_8_9   
     9  +coeff(27)    *x22*x34*x41    
     1  +coeff(28)        *x34*x43    
     2  +coeff(29)*x12*x22*x33        
     3  +coeff(30)*x12    *x34*x41    
c
      return
      end
      function p_m35_8_9   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.26329E+00,-0.22254E+00,-0.58060E+00,-0.86625E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.40506E+00, 0.80612E-01, 0.58060E+00, 0.86625E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.48797953E-01, 0.84622204E-01, 0.24217392E-02, 0.61685554E-02,
     +  0.62407704E-03, 0.20209022E-02, 0.64890530E-03,-0.81866100E-03,
     + -0.32349160E-02, 0.10802331E-02,-0.13738710E-02, 0.71667260E-03,
     + -0.56303130E-02, 0.15942520E-02, 0.45614461E-02,-0.12029841E-02,
     + -0.77565400E-02,-0.69056130E-02,-0.10765750E-01,-0.69933631E-02,
     + -0.69907624E-02, 0.74934180E-02, 0.12665610E-02, 0.66576483E-02,
     +  0.69887763E-02, 0.43611220E-02,-0.10802800E-02, 0.44672080E-03,
     +  0.31461700E-01, 0.21775880E-01, 0.66868710E-02,-0.40617492E-02,
     +  0.53356760E-03,-0.61331200E-02, 0.20784370E-02, 0.84283994E-02,
     +  0.68917000E-02,-0.13827800E-04, 0.33556594E-03, 0.27435870E-02,
     +  0.82037951E-02, 0.15330250E-02,-0.34427590E-03,-0.22756401E-02,
     + -0.68302272E-03,-0.75933882E-02,-0.92389900E-02,-0.22689034E-02,
     + -0.93113400E-02, 0.21751022E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      p_m35_8_9   =avdat
     1  +coeff( 1)        *x31        
     2  +coeff( 2)            *x41    
     3  +coeff( 3)    *x21*x31        
     4  +coeff( 4)        *x31    *x51
     5  +coeff( 5)            *x41*x51
     6  +coeff( 6)*x11    *x31        
     7  +coeff( 7)    *x22*x31        
     8  +coeff( 8)        *x33        
      p_m35_8_9   =p_m35_8_9   
     9  +coeff( 9)        *x32*x41    
     1  +coeff(10)    *x21*x31    *x51
     2  +coeff(11)        *x31    *x52
     3  +coeff(12)            *x41*x52
     4  +coeff(13)*x11*x21*x31        
     5  +coeff(14)*x12    *x31        
     6  +coeff(15)*x11    *x33        
     7  +coeff(16)    *x22*x33        
     8  +coeff(17)*x11    *x32*x41    
      p_m35_8_9   =p_m35_8_9   
     9  +coeff(18)*x11    *x31*x42    
     1  +coeff(19)*x12    *x33        
     2  +coeff(20)*x12    *x32*x41    
     3  +coeff(21)*x12*x21*x31        
     4  +coeff(22)*x11    *x34*x41    
     5  +coeff(23)*x11    *x33*x42    
     6  +coeff(24)        *x34*x43    
     7  +coeff(25)        *x33*x44    
     8  +coeff(26)            *x43    
      p_m35_8_9   =p_m35_8_9   
     9  +coeff(27)    *x21    *x41*x52
     1  +coeff(28)        *x34*x41    
     2  +coeff(29)*x12*x21*x33        
     3  +coeff(30)*x12*x21*x32*x41    
     4  +coeff(31)    *x23*x34*x41    
     5  +coeff(32)*x11*x23*x31*x42    
     6  +coeff(33)*x11*x24*x31    *x51
     7  +coeff(34)*x11    *x33*x42*x51
     8  +coeff(35)*x11*x22*x31        
      p_m35_8_9   =p_m35_8_9   
     9  +coeff(36)    *x22*x31*x42    
     1  +coeff(37)    *x22    *x43    
     2  +coeff(38)    *x23*x31    *x51
     3  +coeff(39)*x11*x23*x31        
     4  +coeff(40)*x12    *x31    *x51
     5  +coeff(41)*x11    *x33    *x51
     6  +coeff(42)    *x22*x33    *x51
     7  +coeff(43)*x11*x21*x31    *x52
     8  +coeff(44)*x11*x22*x33        
      p_m35_8_9   =p_m35_8_9   
     9  +coeff(45)    *x24*x33        
     1  +coeff(46)    *x22*x34*x41    
     2  +coeff(47)    *x22*x33*x42    
     3  +coeff(48)*x12*x21*x31    *x51
     4  +coeff(49)*x11*x21*x33    *x51
     5  +coeff(50)    *x21*x34*x41*x51
c
      return
      end
      function l_m35_8_9   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(25)
      data ncoeff/ 24/
      data avdat/ -0.1173386E-02/
      data xmin/
     1 -0.26329E+00,-0.22254E+00,-0.58060E+00,-0.86625E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.40506E+00, 0.80612E-01, 0.58060E+00, 0.86625E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.57832300E-03, 0.39354264E-02, 0.84604780E-04,-0.57426550E-02,
     + -0.11906154E-03, 0.13977190E-02,-0.17982801E-03, 0.35133384E-03,
     + -0.99847970E-03,-0.15664211E-02,-0.12033850E-03, 0.24356722E-03,
     +  0.91250774E-04,-0.98740290E-04, 0.14807120E-04,-0.63189080E-04,
     +  0.14051783E-03, 0.63323372E-04, 0.91766502E-04, 0.98686490E-04,
     +  0.21888150E-03,-0.87696484E-04,-0.81324940E-04,-0.12523933E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
c
c                  function
c
      l_m35_8_9   =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)                *x51
     4  +coeff( 4)    *x22            
     5  +coeff( 5)        *x32        
     6  +coeff( 6)        *x31*x41    
     7  +coeff( 7)    *x21        *x51
     8  +coeff( 8)*x11                
      l_m35_8_9   =l_m35_8_9   
     9  +coeff( 9)*x11*x21            
     1  +coeff(10)            *x42    
     2  +coeff(11)    *x23            
     3  +coeff(12)            *x42*x51
     4  +coeff(13)    *x22        *x51
     5  +coeff(14)*x12                
     6  +coeff(15)    *x21    *x42*x51
     7  +coeff(16)*x11*x21*x32        
     8  +coeff(17)    *x21*x32        
      l_m35_8_9   =l_m35_8_9   
     9  +coeff(18)*x12*x21            
     1  +coeff(19)*x12    *x32    *x51
     2  +coeff(20)*x11    *x32        
     3  +coeff(21)*x11*x21        *x51
     4  +coeff(22)*x11*x23            
     5  +coeff(23)*x11            *x51
     6  +coeff(24)        *x34        
c
      return
      end
      function x_m35_9_10  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(43)
      data ncoeff/ 42/
      data avdat/ -0.4732611E-01/
      data xmin/
     1 -0.35053E+00,-0.25790E+00,-0.54116E+00,-0.14564E+00,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.43230E+00, 0.10910E+00, 0.54116E+00, 0.14564E+00, 0.14986E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.51580842E-01, 0.92487812E-01, 0.44168144E-03,-0.48987080E-03,
     + -0.11210520E-02,-0.21669660E-02,-0.11459560E-02, 0.43181480E-03,
     + -0.22470110E-03, 0.39859200E+00,-0.69981743E-03,-0.45614050E-03,
     +  0.25691662E-03,-0.37886070E-03, 0.18169052E-02, 0.13674090E-02,
     + -0.57499580E-03,-0.54529344E-03,-0.27610560E-03, 0.23118150E-02,
     +  0.16238171E-02,-0.15356243E-02, 0.19872950E-02, 0.16121482E-02,
     +  0.35186850E-04,-0.17964250E-03, 0.67699530E-03,-0.12123870E-02,
     +  0.46484220E-03,-0.99985243E-03, 0.21560940E-02,-0.30383851E-04,
     + -0.26655721E-02,-0.24864070E-02,-0.38946542E-03, 0.12116710E-02,
     +  0.19137850E-02, 0.34813422E-02,-0.15441540E-02,-0.35383721E-03,
     +  0.31874190E-03, 0.69270160E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x35 = x34*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_m35_9_10  =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)                *x51
     4  +coeff( 4)    *x22            
     5  +coeff( 5)        *x32        
     6  +coeff( 6)        *x31*x41    
     7  +coeff( 7)            *x42    
     8  +coeff( 8)    *x21        *x51
      x_m35_9_10  =x_m35_9_10  
     9  +coeff( 9)                *x52
     1  +coeff(10)*x11                
     2  +coeff(11)    *x23            
     3  +coeff(12)*x12    *x32        
     4  +coeff(13)*x12*x21        *x51
     5  +coeff(14)*x12                
     6  +coeff(15)    *x21*x34        
     7  +coeff(16)*x11*x21*x35*x41    
     8  +coeff(17)        *x32    *x51
      x_m35_9_10  =x_m35_9_10  
     9  +coeff(18)*x11            *x51
     1  +coeff(19)*x11    *x32        
     2  +coeff(20)    *x21*x33*x41    
     3  +coeff(21)*x12*x21            
     4  +coeff(22)*x11*x22*x32        
     5  +coeff(23)*x11    *x34        
     6  +coeff(24)*x11*x22*x34        
     7  +coeff(25)                *x53
     8  +coeff(26)    *x22*x32        
      x_m35_9_10  =x_m35_9_10  
     9  +coeff(27)        *x34        
     1  +coeff(28)    *x22*x31*x41    
     2  +coeff(29)    *x21*x32    *x51
     3  +coeff(30)*x12    *x31*x41    
     4  +coeff(31)    *x23*x33*x41    
     5  +coeff(32)*x11    *x34    *x51
     6  +coeff(33)*x12*x21*x32        
     7  +coeff(34)*x12*x22*x32        
     8  +coeff(35)    *x23    *x44*x52
      x_m35_9_10  =x_m35_9_10  
     9  +coeff(36)*x11*x23*x32    *x52
     1  +coeff(37)    *x21*x31*x41*x51
     2  +coeff(38)    *x23*x32        
     3  +coeff(39)*x11    *x32    *x51
     4  +coeff(40)*x11*x23        *x52
     5  +coeff(41)*x12*x24            
     6  +coeff(42)*x12    *x34        
c
      return
      end
      function t_m35_9_10  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/ -0.2032666E-01/
      data xmin/
     1 -0.35053E+00,-0.25790E+00,-0.54116E+00,-0.14564E+00,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.43230E+00, 0.10910E+00, 0.54116E+00, 0.14564E+00, 0.14986E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.52016284E-01, 0.18574150E+00, 0.25723050E-02, 0.25261650E-01,
     +  0.66619773E-03,-0.25064761E-02,-0.36270800E-02,-0.17068370E-02,
     + -0.33035580E-02, 0.56099402E-02,-0.34513980E-02,-0.35090630E-02,
     +  0.77222590E-02,-0.46178590E-03,-0.10313532E-02,-0.30070080E-02,
     + -0.72754050E-04, 0.49849064E-02, 0.25240820E-02,-0.79900963E-03,
     +  0.26116220E-02,-0.26700980E-03,-0.10955920E-02,-0.10191630E-02,
     +  0.10424300E-02,-0.81355750E-02, 0.22021604E-01, 0.72795260E-03,
     + -0.68254051E-02, 0.71392051E-03, 0.11071600E-01, 0.70457281E-02,
     +  0.69205110E-02,-0.35604812E-01, 0.27036890E-02,-0.13514970E-02,
     + -0.65126543E-03,-0.12399702E-02, 0.40770634E-02,-0.38801140E-02,
     +  0.38657370E-02,-0.36363301E-02,-0.10055610E-02, 0.44992120E-02,
     +  0.28979191E-02,-0.20812781E-01, 0.25943530E-03,-0.85538364E-02,
     +  0.12018701E-01, 0.24446650E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      t_m35_9_10  =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)                *x51
     4  +coeff( 4)*x11                
     5  +coeff( 5)    *x22            
     6  +coeff( 6)        *x32        
     7  +coeff( 7)        *x31*x41    
     8  +coeff( 8)            *x42    
      t_m35_9_10  =t_m35_9_10  
     9  +coeff( 9)    *x21        *x51
     1  +coeff(10)*x11*x21            
     2  +coeff(11)    *x23            
     3  +coeff(12)*x11            *x51
     4  +coeff(13)    *x22        *x51
     5  +coeff(14)*x12                
     6  +coeff(15)*x11*x22            
     7  +coeff(16)    *x24            
     8  +coeff(17)*x11    *x32        
      t_m35_9_10  =t_m35_9_10  
     9  +coeff(18)    *x22*x32        
     1  +coeff(19)        *x34        
     2  +coeff(20)*x11    *x31*x41    
     3  +coeff(21)    *x22*x31*x41    
     4  +coeff(22)        *x33*x41    
     5  +coeff(23)*x11            *x52
     6  +coeff(24)        *x32    *x52
     7  +coeff(25)            *x42*x52
     8  +coeff(26)*x11*x21*x32        
      t_m35_9_10  =t_m35_9_10  
     9  +coeff(27)    *x23*x32        
     1  +coeff(28)*x11*x21*x31*x41    
     2  +coeff(29)*x11    *x32    *x51
     3  +coeff(30)        *x34    *x51
     4  +coeff(31)*x12*x22            
     5  +coeff(32)*x11*x22*x32        
     6  +coeff(33)*x11    *x34        
     7  +coeff(34)*x12*x21*x32        
     8  +coeff(35)*x12    *x34        
      t_m35_9_10  =t_m35_9_10  
     9  +coeff(36)                *x52
     1  +coeff(37)    *x23        *x51
     2  +coeff(38)    *x21*x32    *x51
     3  +coeff(39)*x12    *x32        
     4  +coeff(40)        *x32    *x51
     5  +coeff(41)*x11*x21        *x51
     6  +coeff(42)    *x22*x32    *x51
     7  +coeff(43)    *x22*x34        
     8  +coeff(44)    *x24    *x42    
      t_m35_9_10  =t_m35_9_10  
     9  +coeff(45)*x11*x21*x31*x41*x51
     1  +coeff(46)*x12*x21*x31*x41    
     2  +coeff(47)    *x21*x32        
     3  +coeff(48)*x11*x23            
     4  +coeff(49)    *x23*x31*x41    
     5  +coeff(50)    *x21*x32    *x52
c
      return
      end
      function y_m35_9_10  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(17)
      data ncoeff/ 16/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.35053E+00,-0.25790E+00,-0.54116E+00,-0.14564E+00,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.43230E+00, 0.10910E+00, 0.54116E+00, 0.14564E+00, 0.14986E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.53144651E+00, 0.72608970E-01, 0.52274460E-03, 0.74341881E-03,
     +  0.13037350E-02,-0.25032880E-02,-0.14294820E-02, 0.18847350E-02,
     +  0.26777074E-02, 0.26859910E-02,-0.10395170E-02, 0.94309150E-03,
     +  0.18319860E-02,-0.48027070E-03,-0.31660814E-03,-0.48364864E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x51 = x5
c
c                  function
c
      y_m35_9_10  =avdat
     1  +coeff( 1)        *x31        
     2  +coeff( 2)            *x41    
     3  +coeff( 3)    *x22*x31        
     4  +coeff( 4)        *x31    *x51
     5  +coeff( 5)        *x33        
     6  +coeff( 6)*x11*x21*x31        
     7  +coeff( 7)*x12    *x31        
     8  +coeff( 8)*x11    *x33        
      y_m35_9_10  =y_m35_9_10  
     9  +coeff( 9)*x12    *x33        
     1  +coeff(10)        *x32*x41    
     2  +coeff(11)            *x41*x51
     3  +coeff(12)*x11    *x31        
     4  +coeff(13)*x11        *x41    
     5  +coeff(14)        *x33    *x51
     6  +coeff(15)*x12*x21*x31        
     7  +coeff(16)*x11*x23*x31        
c
      return
      end
      function p_m35_9_10  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(48)
      data ncoeff/ 47/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.35053E+00,-0.25790E+00,-0.54116E+00,-0.14564E+00,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.43230E+00, 0.10910E+00, 0.54116E+00, 0.14564E+00, 0.14986E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.42673200E-01, 0.13897450E+00, 0.74383750E-02, 0.50974491E-03,
     + -0.14526802E-01,-0.17705682E-01, 0.10972280E-01, 0.20598812E-01,
     +  0.35662350E+00, 0.15015100E-01, 0.31654490E-01, 0.10833944E+01,
     +  0.10467902E+01, 0.30627760E+00, 0.14094191E-01, 0.10166463E-01,
     + -0.13312504E-02,-0.14431044E-02,-0.24481060E-01, 0.62461660E-02,
     + -0.40174490E+00,-0.83981650E+00,-0.30358153E+00,-0.13645670E-01,
     +  0.27653590E-01,-0.96948854E-02,-0.37206740E-02,-0.34972820E-01,
     + -0.30137581E-01, 0.15542800E-01, 0.14346950E+00, 0.14676640E-01,
     +  0.12139590E-02, 0.15758630E-01, 0.26606990E-01, 0.23716210E+00,
     +  0.25067713E+00,-0.31181620E-01, 0.10977000E-01,-0.24585781E-01,
     +  0.14958040E-01,-0.24404290E-02,-0.29552740E-01, 0.14453900E+00,
     +  0.20616200E-01,-0.47824060E-02, 0.82147420E-01,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      p_m35_9_10  =avdat
     1  +coeff( 1)        *x31        
     2  +coeff( 2)            *x41    
     3  +coeff( 3)    *x21*x31        
     4  +coeff( 4)    *x21    *x41    
     5  +coeff( 5)        *x31    *x51
     6  +coeff( 6)            *x41*x51
     7  +coeff( 7)*x11    *x31        
     8  +coeff( 8)    *x22*x31        
      p_m35_9_10  =p_m35_9_10  
     9  +coeff( 9)        *x33        
     1  +coeff(10)*x11        *x41    
     2  +coeff(11)    *x22    *x41    
     3  +coeff(12)        *x32*x41    
     4  +coeff(13)        *x31*x42    
     5  +coeff(14)            *x43    
     6  +coeff(15)    *x21*x31    *x51
     7  +coeff(16)    *x21    *x41*x51
     8  +coeff(17)        *x31    *x52
      p_m35_9_10  =p_m35_9_10  
     9  +coeff(18)            *x41*x52
     1  +coeff(19)*x11*x21*x31        
     2  +coeff(20)    *x23*x31        
     3  +coeff(21)    *x21*x33        
     4  +coeff(22)    *x21*x32*x41    
     5  +coeff(23)    *x21*x31*x42    
     6  +coeff(24)*x11    *x31    *x51
     7  +coeff(25)        *x33    *x51
     8  +coeff(26)*x11        *x41*x51
      p_m35_9_10  =p_m35_9_10  
     9  +coeff(27)*x12    *x31        
     1  +coeff(28)*x11*x22*x31        
     2  +coeff(29)*x11    *x33        
     3  +coeff(30)        *x34*x41    
     4  +coeff(31)*x11    *x31*x42    
     5  +coeff(32)        *x33*x42    
     6  +coeff(33)*x11*x21*x31    *x51
     7  +coeff(34)*x12*x21*x31        
     8  +coeff(35)*x11*x23*x31        
      p_m35_9_10  =p_m35_9_10  
     9  +coeff(36)*x11*x21*x33        
     1  +coeff(37)*x11*x21*x32*x41    
     2  +coeff(38)*x12*x22*x31        
     3  +coeff(39)*x11*x24*x31        
     4  +coeff(40)*x12    *x33    *x51
     5  +coeff(41)    *x24*x32*x41*x51
     6  +coeff(42)*x12*x21*x31    *x52
     7  +coeff(43)*x11*x21    *x41    
     8  +coeff(44)    *x21    *x43    
      p_m35_9_10  =p_m35_9_10  
     9  +coeff(45)        *x32*x41*x51
     1  +coeff(46)    *x22*x33        
     2  +coeff(47)*x11    *x32*x41    
c
      return
      end
      function l_m35_9_10  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(25)
      data ncoeff/ 24/
      data avdat/ -0.2190780E-02/
      data xmin/
     1 -0.35053E+00,-0.25790E+00,-0.54116E+00,-0.14564E+00,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.43230E+00, 0.10910E+00, 0.54116E+00, 0.14564E+00, 0.14986E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.83867344E-03, 0.67456970E-02,-0.32754070E-04, 0.53498500E-03,
     + -0.83882580E-02,-0.29138140E-03, 0.94317760E-03,-0.13512832E-02,
     + -0.54652960E-02,-0.16501630E-05,-0.69270274E-04,-0.21417211E-03,
     +  0.27534000E-03,-0.61148421E-04, 0.43386942E-03,-0.16380460E-03,
     +  0.92876490E-05,-0.17203430E-03,-0.16093101E-03, 0.25238180E-03,
     + -0.15905870E-03,-0.72931703E-04, 0.19451601E-03,-0.35399110E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      l_m35_9_10  =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)                *x51
     4  +coeff( 4)*x11                
     5  +coeff( 5)    *x22            
     6  +coeff( 6)        *x32        
     7  +coeff( 7)        *x31*x41    
     8  +coeff( 8)*x11*x21            
      l_m35_9_10  =l_m35_9_10  
     9  +coeff( 9)            *x42    
     1  +coeff(10)                *x53
     2  +coeff(11)    *x23            
     3  +coeff(12)*x11            *x51
     4  +coeff(13)        *x32    *x51
     5  +coeff(14)        *x34        
     6  +coeff(15)*x11*x21        *x51
     7  +coeff(16)    *x21*x32    *x51
     8  +coeff(17)*x11            *x52
      l_m35_9_10  =l_m35_9_10  
     9  +coeff(18)*x12*x21            
     1  +coeff(19)*x11*x22        *x51
     2  +coeff(20)    *x22*x31*x41    
     3  +coeff(21)            *x44    
     4  +coeff(22)            *x42*x52
     5  +coeff(23)    *x24*x32        
     6  +coeff(24)*x12    *x32        
c
      return
      end
      function x_m35_10_11 (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(13)
      data ncoeff/ 12/
      data avdat/ -0.6386141E-01/
      data xmin/
     1 -0.48171E+00,-0.29529E+00,-0.48998E+00,-0.19107E+00,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.44277E+00, 0.12999E+00, 0.48998E+00, 0.19107E+00, 0.14986E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.29797520E-02, 0.10642910E+00, 0.87458010E-03,-0.74742780E-03,
     + -0.33814244E-03, 0.10978501E-02,-0.10900464E-02, 0.33445090E-04,
     +  0.46562750E+00, 0.10895530E-02,-0.17225982E-02, 0.15597400E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      x_m35_10_11 =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)                *x51
     4  +coeff( 4)    *x22            
     5  +coeff( 5)        *x32        
     6  +coeff( 6)            *x42    
     7  +coeff( 7)    *x21        *x51
     8  +coeff( 8)                *x52
      x_m35_10_11 =x_m35_10_11 
     9  +coeff( 9)*x11                
     1  +coeff(10)*x11*x21            
     2  +coeff(11)*x11    *x34        
     3  +coeff(12)    *x23*x32        
c
      return
      end
      function t_m35_10_11 (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/ -0.2280975E-01/
      data xmin/
     1 -0.48171E+00,-0.29529E+00,-0.48998E+00,-0.19107E+00,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.44277E+00, 0.12999E+00, 0.48998E+00, 0.19107E+00, 0.14986E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.60964170E-01, 0.21394044E+00, 0.13920511E-01, 0.15334630E-02,
     +  0.16461210E-01, 0.32223740E-01, 0.17791220E-01, 0.14942940E-02,
     + -0.46218044E-03, 0.40119811E-02,-0.31901694E-01,-0.56455910E-01,
     +  0.96608340E-03, 0.29431912E-02, 0.24250620E-02, 0.25393420E-02,
     +  0.16406720E-01, 0.82601350E-02, 0.53773240E-01, 0.90232550E-02,
     + -0.37299364E-02,-0.59138020E-03, 0.15207702E-01,-0.28670351E-02,
     +  0.31831180E-01, 0.97644510E-02,-0.32479780E-03,-0.20799821E-01,
     + -0.43820664E-02, 0.20005843E-04,-0.32217204E-02,-0.52419914E-02,
     + -0.21809710E-03,-0.23619912E-01, 0.24568694E-02,-0.19328610E-01,
     + -0.58388580E-02,-0.65707373E-02,-0.15425000E-01, 0.28744670E-01,
     +  0.12455880E-01,-0.42688734E-02, 0.16586080E-01, 0.23287180E-04,
     +  0.34755894E-02,-0.20431370E-01,-0.24363812E-01, 0.37385781E-02,
     +  0.23181301E-02, 0.31378120E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      t_m35_10_11 =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)*x11                
     4  +coeff( 4)    *x22            
     5  +coeff( 5)        *x32        
     6  +coeff( 6)        *x31*x41    
     7  +coeff( 7)            *x42    
     8  +coeff( 8)    *x21        *x51
      t_m35_10_11 =t_m35_10_11 
     9  +coeff( 9)                *x52
     1  +coeff(10)    *x23            
     2  +coeff(11)    *x21*x32        
     3  +coeff(12)    *x21*x31*x41    
     4  +coeff(13)*x11            *x51
     5  +coeff(14)        *x32    *x51
     6  +coeff(15)        *x31*x41*x51
     7  +coeff(16)*x12                
     8  +coeff(17)*x11    *x32        
      t_m35_10_11 =t_m35_10_11 
     9  +coeff(18)    *x22*x32        
     1  +coeff(19)*x11    *x31*x41    
     2  +coeff(20)    *x22*x31*x41    
     3  +coeff(21)*x11*x21*x32        
     4  +coeff(22)    *x21*x34        
     5  +coeff(23)*x11*x24            
     6  +coeff(24)*x11    *x34        
     7  +coeff(25)*x11*x23        *x51
     8  +coeff(26)*x12*x21*x32        
      t_m35_10_11 =t_m35_10_11 
     9  +coeff(27)*x12*x21            
     1  +coeff(28)*x12*x23            
     2  +coeff(29)*x12    *x32    *x51
     3  +coeff(30)                *x51
     4  +coeff(31)*x12    *x32        
     5  +coeff(32)*x11*x22        *x52
     6  +coeff(33)*x11    *x31*x41*x52
     7  +coeff(34)    *x21    *x42    
     8  +coeff(35)    *x22        *x52
      t_m35_10_11 =t_m35_10_11 
     9  +coeff(36)*x11*x21*x31*x41    
     1  +coeff(37)*x11*x21            
     2  +coeff(38)*x11*x22            
     3  +coeff(39)    *x24            
     4  +coeff(40)*x11        *x42    
     5  +coeff(41)    *x23        *x51
     6  +coeff(42)    *x21*x32    *x51
     7  +coeff(43)*x11*x23            
     8  +coeff(44)    *x23    *x42    
      t_m35_10_11 =t_m35_10_11 
     9  +coeff(45)*x12            *x51
     1  +coeff(46)*x11*x22        *x51
     2  +coeff(47)    *x24        *x51
     3  +coeff(48)    *x22*x32    *x51
     4  +coeff(49)*x11    *x31*x41*x51
     5  +coeff(50)*x11*x21        *x52
c
      return
      end
      function y_m35_10_11 (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 9)
      data ncoeff/  8/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.48171E+00,-0.29529E+00,-0.48998E+00,-0.19107E+00,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.44277E+00, 0.12999E+00, 0.48998E+00, 0.19107E+00, 0.14986E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.48586330E+00, 0.95410970E-01,-0.13799710E-02, 0.15187413E-02,
     +  0.47098910E-03,-0.49039070E-03,-0.23098290E-02, 0.33715190E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x51 = x5
c
c                  function
c
      y_m35_10_11 =avdat
     1  +coeff( 1)        *x31        
     2  +coeff( 2)            *x41    
     3  +coeff( 3)*x11*x21*x31        
     4  +coeff( 4)*x12*x22*x31        
     5  +coeff( 5)        *x31    *x51
     6  +coeff( 6)    *x22*x32*x41    
     7  +coeff( 7)*x11*x21*x33        
     8  +coeff( 8)    *x21*x31        
c
      return
      end
      function p_m35_10_11 (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(36)
      data ncoeff/ 35/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.48171E+00,-0.29529E+00,-0.48998E+00,-0.19107E+00,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.44277E+00, 0.12999E+00, 0.48998E+00, 0.19107E+00, 0.14986E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.13375610E-01, 0.19018980E+00,-0.19249940E-02, 0.77767483E-03,
     + -0.11694413E-02,-0.96572330E-03, 0.76292290E-03,-0.48633080E-02,
     + -0.34178074E-01,-0.77787674E-01, 0.10775790E-02, 0.12052563E-01,
     +  0.34210470E-01,-0.30658680E-02,-0.17179680E-01,-0.62510440E-02,
     +  0.36972950E-02,-0.80192312E-02, 0.25211230E-02,-0.27497363E-03,
     +  0.51006340E-02, 0.38413531E-02,-0.26248754E-02,-0.44403230E-01,
     +  0.31975440E-02, 0.25354780E-02,-0.13359953E-01,-0.25923990E-01,
     +  0.21169700E-02,-0.10926620E-01, 0.13978893E-01,-0.44839790E-02,
     +  0.17361663E-01, 0.38669850E-02, 0.39147300E-01,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
c
c                  function
c
      p_m35_10_11 =avdat
     1  +coeff( 1)        *x31        
     2  +coeff( 2)            *x41    
     3  +coeff( 3)    *x21    *x41    
     4  +coeff( 4)        *x31    *x51
     5  +coeff( 5)            *x41*x51
     6  +coeff( 6)*x11    *x31        
     7  +coeff( 7)        *x33        
     8  +coeff( 8)*x11*x21*x31        
      p_m35_10_11 =p_m35_10_11 
     9  +coeff( 9)    *x21*x33        
     1  +coeff(10)    *x21*x32*x41    
     2  +coeff(11)        *x33    *x51
     3  +coeff(12)    *x24*x31        
     4  +coeff(13)    *x22*x33        
     5  +coeff(14)*x12        *x41    
     6  +coeff(15)*x11*x23*x31        
     7  +coeff(16)*x11*x21*x33        
     8  +coeff(17)*x12*x22*x31        
      p_m35_10_11 =p_m35_10_11 
     9  +coeff(18)*x12    *x31*x42    
     1  +coeff(19)    *x22*x31        
     2  +coeff(20)        *x34*x41    
     3  +coeff(21)    *x22    *x41    
     4  +coeff(22)    *x23*x31        
     5  +coeff(23)*x11*x21    *x41    
     6  +coeff(24)    *x21*x31*x42    
     7  +coeff(25)    *x22    *x41*x51
     8  +coeff(26)*x11    *x33        
      p_m35_10_11 =p_m35_10_11 
     9  +coeff(27)    *x23*x33        
     1  +coeff(28)*x11*x23    *x41    
     2  +coeff(29)*x11    *x33    *x51
     3  +coeff(30)*x11*x24*x31        
     4  +coeff(31)*x12*x23*x31        
     5  +coeff(32)*x11*x21*x34*x41    
     6  +coeff(33)    *x24    *x41    
     7  +coeff(34)*x11    *x32*x41    
     8  +coeff(35)    *x22*x32*x41    
c
      return
      end
      function l_m35_10_11 (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(25)
      data ncoeff/ 24/
      data avdat/ -0.2929458E-02/
      data xmin/
     1 -0.48171E+00,-0.29529E+00,-0.48998E+00,-0.19107E+00,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.44277E+00, 0.12999E+00, 0.48998E+00, 0.19107E+00, 0.14986E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.94140800E-03, 0.94947680E-02, 0.39863950E-04, 0.36309272E-03,
     + -0.11246910E-01, 0.55268420E-03, 0.22494494E-02,-0.81033590E-02,
     + -0.11987720E-03,-0.87410461E-03,-0.76930600E-04,-0.16154490E-03,
     +  0.22346720E-03, 0.18690683E-03, 0.23959940E-03,-0.16845972E-03,
     + -0.57907462E-04,-0.33519854E-05, 0.31047844E-03, 0.22545842E-03,
     + -0.24447593E-03,-0.56174241E-04, 0.59705780E-04,-0.13457771E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_m35_10_11 =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)                *x51
     4  +coeff( 4)*x11                
     5  +coeff( 5)    *x22            
     6  +coeff( 6)        *x32        
     7  +coeff( 7)        *x31*x41    
     8  +coeff( 8)            *x42    
      l_m35_10_11 =l_m35_10_11 
     9  +coeff( 9)    *x21        *x51
     1  +coeff(10)*x11*x21            
     2  +coeff(11)    *x21*x32        
     3  +coeff(12)        *x32    *x51
     4  +coeff(13)    *x24            
     5  +coeff(14)*x11    *x31*x41    
     6  +coeff(15)    *x22    *x42    
     7  +coeff(16)    *x23            
     8  +coeff(17)        *x34        
      l_m35_10_11 =l_m35_10_11 
     9  +coeff(18)                *x52
     1  +coeff(19)            *x42*x51
     2  +coeff(20)    *x21*x32    *x51
     3  +coeff(21)    *x24*x32        
     4  +coeff(22)*x11    *x34        
     5  +coeff(23)    *x22        *x51
     6  +coeff(24)        *x33*x41    
c
      return
      end
      function x_m35_11_12 (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(19)
      data ncoeff/ 18/
      data avdat/ -0.7622619E-01/
      data xmin/
     1 -0.51878E+00,-0.25621E+00,-0.43034E+00,-0.19074E+00,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.45836E+00, 0.12542E+00, 0.43034E+00, 0.19074E+00, 0.14986E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.13270232E-01, 0.95472812E-01, 0.15283600E-03,-0.15417890E-03,
     + -0.13332880E-03,-0.58117700E-03,-0.14402700E-03,-0.26863910E-03,
     +  0.16451762E-04, 0.48908470E+00, 0.11675924E-03, 0.25084860E-03,
     +  0.18090760E-03,-0.74171100E-03, 0.14511471E-03, 0.23142911E-03,
     +  0.25774820E-03,-0.68677123E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      x_m35_11_12 =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)                *x51
     4  +coeff( 4)    *x22            
     5  +coeff( 5)        *x32        
     6  +coeff( 6)        *x31*x41    
     7  +coeff( 7)            *x42    
     8  +coeff( 8)    *x21        *x51
      x_m35_11_12 =x_m35_11_12 
     9  +coeff( 9)                *x52
     1  +coeff(10)*x11                
     2  +coeff(11)    *x21*x32        
     3  +coeff(12)    *x23            
     4  +coeff(13)    *x21*x31*x41    
     5  +coeff(14)*x11    *x32        
     6  +coeff(15)*x11*x22            
     7  +coeff(16)*x12                
     8  +coeff(17)    *x23*x32        
      x_m35_11_12 =x_m35_11_12 
     9  +coeff(18)*x12*x21            
c
      return
      end
      function t_m35_11_12 (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(38)
      data ncoeff/ 37/
      data avdat/ -0.2211770E-01/
      data xmin/
     1 -0.51878E+00,-0.25621E+00,-0.43034E+00,-0.19074E+00,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.45836E+00, 0.12542E+00, 0.43034E+00, 0.19074E+00, 0.14986E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.43519520E-01, 0.19104294E+00,-0.17012020E-03, 0.17368280E-02,
     +  0.68152631E-03,-0.27583510E-03,-0.13819360E-02, 0.21927393E-03,
     + -0.63406570E-03,-0.20909523E-02, 0.54212124E-03,-0.11287600E-03,
     +  0.96331794E-04,-0.13996512E-04,-0.32351044E-02,-0.50836950E-03,
     + -0.11443100E-02, 0.19620930E-03,-0.40103800E-03,-0.78330370E-04,
     +  0.70769240E-03, 0.13693690E-02,-0.21908571E-02,-0.19005523E-04,
     +  0.70196161E-04,-0.22532534E-03,-0.13008800E-02, 0.18750383E-02,
     +  0.10411804E-02,-0.10308972E-03,-0.17682440E-03, 0.58387231E-05,
     +  0.36459500E-02,-0.12561454E-03, 0.66092514E-04,-0.18323350E-02,
     +  0.46518960E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_m35_11_12 =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)                *x51
     4  +coeff( 4)*x11                
     5  +coeff( 5)        *x32        
     6  +coeff( 6)        *x31*x41    
     7  +coeff( 7)    *x23            
     8  +coeff( 8)    *x21*x32        
      t_m35_11_12 =t_m35_11_12 
     9  +coeff( 9)*x11            *x51
     1  +coeff(10)*x11    *x32        
     2  +coeff(11)    *x21*x31*x41    
     3  +coeff(12)        *x32    *x51
     4  +coeff(13)    *x21        *x52
     5  +coeff(14)    *x22        *x52
     6  +coeff(15)*x12*x21            
     7  +coeff(16)*x11*x21*x32        
     8  +coeff(17)*x11*x21    *x42    
      t_m35_11_12 =t_m35_11_12 
     9  +coeff(18)*x11*x24            
     1  +coeff(19)*x12    *x32        
     2  +coeff(20)*x12*x23            
     3  +coeff(21)    *x21        *x51
     4  +coeff(22)*x12                
     5  +coeff(23)*x11*x21            
     6  +coeff(24)            *x42*x51
     7  +coeff(25)    *x24            
     8  +coeff(26)    *x21*x31*x41*x51
      t_m35_11_12 =t_m35_11_12 
     9  +coeff(27)*x11*x22*x32        
     1  +coeff(28)*x12*x21*x32        
     2  +coeff(29)    *x22            
     3  +coeff(30)                *x52
     4  +coeff(31)    *x22        *x51
     5  +coeff(32)                *x53
     6  +coeff(33)*x11*x22            
     7  +coeff(34)    *x21*x32    *x51
     8  +coeff(35)    *x21*x34        
      t_m35_11_12 =t_m35_11_12 
     9  +coeff(36)*x11*x21*x31*x41    
     1  +coeff(37)*x11    *x32    *x51
c
      return
      end
      function y_m35_11_12 (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 4)
      data ncoeff/  3/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.51878E+00,-0.25621E+00,-0.43034E+00,-0.19074E+00,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.45836E+00, 0.12542E+00, 0.43034E+00, 0.19074E+00, 0.14986E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.43006211E+00, 0.95511481E-01, 0.52564820E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x31 = x3
      x41 = x4
c
c                  function
c
      y_m35_11_12 =avdat
     1  +coeff( 1)        *x31        
     2  +coeff( 2)            *x41    
     3  +coeff( 3)*x12    *x31        
c
      return
      end
      function p_m35_11_12 (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 8)
      data ncoeff/  7/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.51878E+00,-0.25621E+00,-0.43034E+00,-0.19074E+00,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.45836E+00, 0.12542E+00, 0.43034E+00, 0.19074E+00, 0.14986E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.11124330E-02, 0.19094340E+00, 0.22435003E-03,-0.38293810E-03,
     +  0.18640770E-02,-0.23382970E-03, 0.40564604E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x51 = x5
c
c                  function
c
      p_m35_11_12 =avdat
     1  +coeff( 1)        *x31        
     2  +coeff( 2)            *x41    
     3  +coeff( 3)        *x31    *x51
     4  +coeff( 4)    *x22*x31        
     5  +coeff( 5)*x12    *x31        
     6  +coeff( 6)*x11    *x31        
     7  +coeff( 7)        *x33        
c
      return
      end
      function l_m35_11_12 (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(25)
      data ncoeff/ 24/
      data avdat/ -0.2934002E-02/
      data xmin/
     1 -0.51878E+00,-0.25621E+00,-0.43034E+00,-0.19074E+00,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.45836E+00, 0.12542E+00, 0.43034E+00, 0.19074E+00, 0.14986E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.15762450E-02, 0.70823760E-02, 0.32658652E-04,-0.90998750E-02,
     +  0.37791283E-04, 0.17787500E-03,-0.89992910E-02,-0.37102432E-06,
     + -0.10776850E-03,-0.33861360E-04, 0.48017614E-05, 0.98710840E-04,
     +  0.13893700E-03,-0.30595380E-04,-0.12129480E-03, 0.13355830E-03,
     + -0.45927972E-04,-0.62013550E-04, 0.24199620E-03,-0.81334540E-04,
     + -0.27609261E-04,-0.22028700E-03,-0.36508140E-04,-0.77617813E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_m35_11_12 =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)*x11                
     4  +coeff( 4)    *x22            
     5  +coeff( 5)        *x32        
     6  +coeff( 6)        *x31*x41    
     7  +coeff( 7)            *x42    
     8  +coeff( 8)                *x52
      l_m35_11_12 =l_m35_11_12 
     9  +coeff( 9)    *x23            
     1  +coeff(10)    *x24    *x42    
     2  +coeff(11)                *x51
     3  +coeff(12)    *x21*x31*x41    
     4  +coeff(13)    *x24            
     5  +coeff(14)    *x22*x32        
     6  +coeff(15)        *x34        
     7  +coeff(16)*x12    *x32        
     8  +coeff(17)*x11        *x44    
      l_m35_11_12 =l_m35_11_12 
     9  +coeff(18)    *x21    *x42    
     1  +coeff(19)    *x22    *x42    
     2  +coeff(20)*x11*x21            
     3  +coeff(21)    *x22        *x51
     4  +coeff(22)        *x33*x41    
     5  +coeff(23)*x11        *x42    
     6  +coeff(24)    *x21*x31*x41*x51
c
      return
      end
      function x_m35_12_13 (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 5)
      data ncoeff/  4/
      data avdat/ -0.1024518E+00/
      data xmin/
     1 -0.64408E+00,-0.25697E+00,-0.36868E+00,-0.19085E+00,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.51160E+00, 0.12609E+00, 0.36868E+00, 0.19085E+00, 0.14986E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.45504450E-01, 0.23514580E+00, 0.57901954E+00, 0.12618530E-01,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
c
c                  function
c
      x_m35_12_13 =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)*x11                
     4  +coeff( 4)    *x22            
c
      return
      end
      function t_m35_12_13 (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 3)
      data ncoeff/  2/
      data avdat/ -0.2211773E-01/
      data xmin/
     1 -0.64408E+00,-0.25697E+00,-0.36868E+00,-0.19085E+00,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.51160E+00, 0.12609E+00, 0.36868E+00, 0.19085E+00, 0.14986E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.43321803E-01, 0.19153350E+00,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
c          set up monomials   functions
      x21 = x2
c
c                  function
c
      t_m35_12_13 =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
c
      return
      end
      function y_m35_12_13 (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 9)
      data ncoeff/  8/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.64408E+00,-0.25697E+00,-0.36868E+00,-0.19085E+00,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.51160E+00, 0.12609E+00, 0.36868E+00, 0.19085E+00, 0.14986E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.36886160E+00, 0.23833300E+00,-0.57643290E-03, 0.34042010E-02,
     +  0.11397642E-03, 0.12943561E-03, 0.76488320E-04, 0.97006320E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x31 = x3
      x41 = x4
      x51 = x5
c
c                  function
c
      y_m35_12_13 =avdat
     1  +coeff( 1)        *x31        
     2  +coeff( 2)            *x41    
     3  +coeff( 3)    *x21*x31        
     4  +coeff( 4)    *x21    *x41    
     5  +coeff( 5)        *x31    *x51
     6  +coeff( 6)            *x41*x51
     7  +coeff( 7)*x11    *x31        
     8  +coeff( 8)*x11        *x41    
c
      return
      end
      function p_m35_12_13 (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 5)
      data ncoeff/  4/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.64408E+00,-0.25697E+00,-0.36868E+00,-0.19085E+00,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.51160E+00, 0.12609E+00, 0.36868E+00, 0.19085E+00, 0.14986E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.17648391E-04, 0.19035580E+00,-0.80411970E-04, 0.31011670E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
c          set up monomials   functions
      x21 = x2
      x31 = x3
      x41 = x4
c
c                  function
c
      p_m35_12_13 =avdat
     1  +coeff( 1)        *x31        
     2  +coeff( 2)            *x41    
     3  +coeff( 3)    *x21*x31        
     4  +coeff( 4)    *x21    *x41    
c
      return
      end
      function l_m35_12_13 (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(25)
      data ncoeff/ 24/
      data avdat/  0.1650415E-02/
      data xmin/
     1 -0.64408E+00,-0.25697E+00,-0.36868E+00,-0.19085E+00,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.51160E+00, 0.12609E+00, 0.36868E+00, 0.19085E+00, 0.14986E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.78749530E-02,-0.30815140E-02,-0.50369560E-01,-0.22786070E-01,
     +  0.34468621E-03, 0.65353553E-03,-0.22357990E-01,-0.58734960E-03,
     + -0.11302924E-02,-0.21244670E-02,-0.12499504E-04,-0.34747610E-05,
     +  0.87205050E-03, 0.32429574E-03,-0.16660994E-02,-0.92748630E-03,
     +  0.16231942E-03,-0.31434330E-03, 0.77234480E-04, 0.17742493E-02,
     +  0.28531360E-05, 0.99029700E-04,-0.12394963E-03,-0.76905023E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
c
c                  function
c
      l_m35_12_13 =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)*x11                
     4  +coeff( 4)    *x22            
     5  +coeff( 5)        *x32        
     6  +coeff( 6)        *x31*x41    
     7  +coeff( 7)            *x42    
     8  +coeff( 8)    *x23            
      l_m35_12_13 =l_m35_12_13 
     9  +coeff( 9)    *x21*x32        
     1  +coeff(10)    *x21*x31*x41    
     2  +coeff(11)    *x21        *x51
     3  +coeff(12)*x12                
     4  +coeff(13)*x11    *x32        
     5  +coeff(14)    *x22    *x42    
     6  +coeff(15)    *x21    *x42    
     7  +coeff(16)*x11*x22            
     8  +coeff(17)    *x24            
      l_m35_12_13 =l_m35_12_13 
     9  +coeff(18)        *x31*x43    
     1  +coeff(19)        *x32    *x51
     2  +coeff(20)*x11    *x31*x41    
     3  +coeff(21)                *x51
     4  +coeff(22)        *x33*x41    
     5  +coeff(23)*x11*x21            
     6  +coeff(24)            *x42*x51
c
      return
      end
      function x_m35_0_2   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 5)
      data ncoeff/  4/
      data avdat/  0.2328079E-02/
      data xmin/
     1 -0.99979E-04,-0.21093E+00,-0.59989E-01,-0.36658E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.21092E+00, 0.59989E-01, 0.36658E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.23377714E-02, 0.28796550E+00, 0.53018550E-02, 0.10059720E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
c
c                  function
c
      x_m35_0_2   =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)    *x22            
     4  +coeff( 4)*x11                
c
      return
      end
      function t_m35_0_2   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 7)
      data ncoeff/  6/
      data avdat/  0.4287185E-03/
      data xmin/
     1 -0.99979E-04,-0.21093E+00,-0.59989E-01,-0.36658E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.21092E+00, 0.59989E-01, 0.36658E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.43397722E-03, 0.21022264E+00, 0.84519110E-04, 0.10192570E-03,
     + -0.79149550E-04, 0.56310030E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x51 = x5
c
c                  function
c
      t_m35_0_2   =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)                *x51
     4  +coeff( 4)    *x21        *x51
     5  +coeff( 5)    *x22            
     6  +coeff( 6)    *x23            
c
      return
      end
      function y_m35_0_2   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 4)
      data ncoeff/  3/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.99979E-04,-0.21093E+00,-0.59989E-01,-0.36658E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.21092E+00, 0.59989E-01, 0.36658E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.59995610E-01, 0.49867402E-01, 0.91860984E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
c          set up monomials   functions
      x21 = x2
      x31 = x3
      x41 = x4
c
c                  function
c
      y_m35_0_2   =avdat
     1  +coeff( 1)        *x31        
     2  +coeff( 2)            *x41    
     3  +coeff( 3)    *x21    *x41    
c
      return
      end
      function p_m35_0_2   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 7)
      data ncoeff/  6/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.99979E-04,-0.21093E+00,-0.59989E-01,-0.36658E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.21092E+00, 0.59989E-01, 0.36658E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.13374642E-03, 0.36901760E-01,-0.30019050E-04, 0.64972800E-03,
     + -0.11144390E-03,-0.80885482E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x41 = x4
c
c                  function
c
      p_m35_0_2   =avdat
     1  +coeff( 1)        *x31        
     2  +coeff( 2)            *x41    
     3  +coeff( 3)    *x21*x31        
     4  +coeff( 4)    *x21    *x41    
     5  +coeff( 5)    *x22*x31        
     6  +coeff( 6)    *x22    *x41    
c
      return
      end
      function l_m35_0_2   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(10)
      data ncoeff/  9/
      data avdat/ -0.1020404E-01/
      data xmin/
     1 -0.99979E-04,-0.21093E+00,-0.59989E-01,-0.36658E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.21092E+00, 0.59989E-01, 0.36658E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.10204930E-01,-0.25083744E-01,-0.30700450E-01,-0.91381470E-03,
     + -0.55757040E-03, 0.31566910E-03,-0.88356360E-05,-0.16951080E-04,
     +  0.19350370E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x41 = x4
      x42 = x41*x4
c
c                  function
c
      l_m35_0_2   =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)    *x22            
     4  +coeff( 4)            *x42    
     5  +coeff( 5)    *x23            
     6  +coeff( 6)    *x24            
     7  +coeff( 7)*x11                
     8  +coeff( 8)    *x21    *x42    
      l_m35_0_2   =l_m35_0_2   
     9  +coeff( 9)    *x22    *x42    
c
      return
      end
      function x_m35_0_3   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(19)
      data ncoeff/ 18/
      data avdat/ -0.1945024E-01/
      data xmin/
     1 -0.99979E-04,-0.21093E+00,-0.59989E-01,-0.36593E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36593E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.32089150E-02, 0.58161770E+00, 0.85149390E-02, 0.15790862E-01,
     +  0.59185420E-02,-0.52470224E-02,-0.12505464E-02,-0.23152264E-02,
     +  0.14599210E-02, 0.21945290E-03, 0.32382944E-03, 0.77063014E-04,
     +  0.20080210E-03, 0.18667412E-03,-0.49542973E-03, 0.24381684E-03,
     +  0.46447594E-03, 0.34352870E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_m35_0_3   =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)                *x51
     4  +coeff( 4)    *x21        *x51
     5  +coeff( 5)    *x22            
     6  +coeff( 6)    *x23            
     7  +coeff( 7)                *x52
     8  +coeff( 8)    *x21        *x52
      x_m35_0_3   =x_m35_0_3   
     9  +coeff( 9)    *x24            
     1  +coeff(10)    *x21*x31*x41    
     2  +coeff(11)    *x22        *x51
     3  +coeff(12)*x11                
     4  +coeff(13)    *x21    *x42    
     5  +coeff(14)                *x53
     6  +coeff(15)    *x22*x31*x41    
     7  +coeff(16)        *x33*x41    
     8  +coeff(17)    *x23        *x51
      x_m35_0_3   =x_m35_0_3   
     9  +coeff(18)    *x21        *x53
c
      return
      end
      function t_m35_0_3   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/ -0.3954356E-02/
      data xmin/
     1 -0.99979E-04,-0.21093E+00,-0.59989E-01,-0.36593E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36593E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.24156081E-02, 0.55546530E-01, 0.12513840E-01,-0.54387082E-02,
     + -0.56628050E-04,-0.57822660E-04, 0.20623683E-01,-0.17656790E-02,
     + -0.51365992E-02,-0.74103770E-04, 0.37505113E-04,-0.48479540E-06,
     +  0.40419102E-03, 0.31780850E-04,-0.28330530E-02, 0.25268550E-03,
     +  0.22464094E-02,-0.10284410E-03,-0.13570832E-03,-0.17405940E-03,
     +  0.41216890E-03,-0.13358390E-02,-0.10925034E-02,-0.31165320E-04,
     +  0.18416640E-03, 0.33184200E-03, 0.25793022E-04,-0.12415370E-03,
     +  0.68825010E-04,-0.28669800E-04, 0.13311110E-03,-0.37314850E-04,
     +  0.33832010E-04,-0.75918700E-05,-0.30568020E-04, 0.43206810E-04,
     + -0.50942230E-04,-0.11248980E-03, 0.69808651E-04, 0.70237710E-04,
     + -0.92725120E-04,-0.56229290E-04,-0.24820300E-03, 0.98115153E-04,
     +  0.42812630E-04, 0.59520640E-04,-0.64500720E-04,-0.23059690E-04,
     +  0.84247651E-04, 0.18605050E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      t_m35_0_3   =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)                *x51
     4  +coeff( 4)    *x22            
     5  +coeff( 5)        *x31*x41    
     6  +coeff( 6)            *x42    
     7  +coeff( 7)    *x21        *x51
     8  +coeff( 8)                *x52
      t_m35_0_3   =t_m35_0_3   
     9  +coeff( 9)    *x23            
     1  +coeff(10)    *x21*x32        
     2  +coeff(11)    *x21*x31*x41    
     3  +coeff(12)    *x21    *x42    
     4  +coeff(13)    *x22        *x51
     5  +coeff(14)        *x31*x41*x51
     6  +coeff(15)    *x21        *x52
     7  +coeff(16)                *x53
     8  +coeff(17)    *x24            
      t_m35_0_3   =t_m35_0_3   
     9  +coeff(18)    *x22*x32        
     1  +coeff(19)    *x22*x31*x41    
     2  +coeff(20)    *x22    *x42    
     3  +coeff(21)    *x21        *x53
     4  +coeff(22)    *x24*x31*x41    
     5  +coeff(23)    *x24    *x42    
     6  +coeff(24)*x11                
     7  +coeff(25)    *x23        *x51
     8  +coeff(26)    *x24*x31*x43    
      t_m35_0_3   =t_m35_0_3   
     9  +coeff(27)            *x42*x51
     1  +coeff(28)    *x24        *x51
     2  +coeff(29)    *x22*x32    *x51
     3  +coeff(30)    *x22*x31*x41*x51
     4  +coeff(31)    *x22    *x44    
     5  +coeff(32)        *x32        
     6  +coeff(33)*x11*x22            
     7  +coeff(34)        *x32*x42    
     8  +coeff(35)                *x54
      t_m35_0_3   =t_m35_0_3   
     9  +coeff(36)    *x21*x34        
     1  +coeff(37)    *x23*x31*x41    
     2  +coeff(38)    *x23    *x42    
     3  +coeff(39)    *x22    *x42*x51
     4  +coeff(40)    *x23        *x52
     5  +coeff(41)    *x21        *x54
     6  +coeff(42)*x11*x24            
     7  +coeff(43)    *x24*x32        
     8  +coeff(44)    *x22*x32*x42    
      t_m35_0_3   =t_m35_0_3   
     9  +coeff(45)    *x21*x32*x42*x51
     1  +coeff(46)    *x24        *x52
     2  +coeff(47)    *x22*x32    *x52
     3  +coeff(48)*x11*x23*x32        
     4  +coeff(49)    *x23    *x44    
     5  +coeff(50)    *x24*x31*x41*x51
c
      return
      end
      function y_m35_0_3   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.99979E-04,-0.21093E+00,-0.59989E-01,-0.36593E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36593E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.75167514E-01, 0.14370781E+00,-0.16366160E-02, 0.29503132E-03,
     + -0.23845320E-02,-0.32172140E-02,-0.85199680E-03,-0.14010670E-02,
     +  0.22234210E-03, 0.24621770E-03, 0.38733530E-03, 0.51398010E-03,
     + -0.10619333E-02,-0.16755100E-02, 0.92476180E-04, 0.31974710E-05,
     +  0.92665861E-04, 0.14491472E-03, 0.71556314E-04, 0.29347040E-04,
     + -0.12930754E-04,-0.88296360E-04, 0.16616220E-03, 0.23989300E-03,
     +  0.30705073E-03, 0.46082834E-04,-0.46249670E-04,-0.43566160E-04,
     + -0.59523780E-04,-0.39869600E-03, 0.19597592E-03, 0.46793020E-03,
     + -0.12178630E-04,-0.16173991E-03,-0.28042550E-04,-0.57684270E-04,
     +  0.31893840E-04, 0.13127810E-03, 0.11083603E-03, 0.11526804E-03,
     +  0.40356582E-04, 0.13133050E-03, 0.64427120E-04,-0.84558094E-04,
     + -0.30349874E-04, 0.15857450E-03, 0.57329880E-04,-0.42245243E-03,
     +  0.17212270E-03,-0.11123304E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      y_m35_0_3   =avdat
     1  +coeff( 1)        *x31        
     2  +coeff( 2)            *x41    
     3  +coeff( 3)    *x21*x31        
     4  +coeff( 4)    *x21    *x41    
     5  +coeff( 5)        *x31    *x51
     6  +coeff( 6)            *x41*x51
     7  +coeff( 7)    *x22*x31        
     8  +coeff( 8)    *x22    *x41    
      y_m35_0_3   =y_m35_0_3   
     9  +coeff( 9)    *x21*x31    *x51
     1  +coeff(10)    *x21    *x41*x51
     2  +coeff(11)        *x31    *x52
     3  +coeff(12)            *x41*x52
     4  +coeff(13)    *x23*x31        
     5  +coeff(14)    *x23    *x41    
     6  +coeff(15)        *x31*x42    
     7  +coeff(16)    *x21*x31*x42    
     8  +coeff(17)    *x22*x31    *x51
      y_m35_0_3   =y_m35_0_3   
     9  +coeff(18)    *x22    *x41*x51
     1  +coeff(19)            *x43    
     2  +coeff(20)    *x21*x32*x41    
     3  +coeff(21)    *x21    *x43    
     4  +coeff(22)            *x41*x53
     5  +coeff(23)    *x23*x31    *x51
     6  +coeff(24)    *x23    *x41*x51
     7  +coeff(25)    *x23*x31*x42    
     8  +coeff(26)    *x24*x32*x41    
      y_m35_0_3   =y_m35_0_3   
     9  +coeff(27)    *x21*x31    *x52
     1  +coeff(28)    *x21    *x41*x52
     2  +coeff(29)        *x31    *x53
     3  +coeff(30)    *x22*x31*x42    
     4  +coeff(31)    *x23    *x43    
     5  +coeff(32)    *x24*x31*x42    
     6  +coeff(33)        *x31*x42*x51
     7  +coeff(34)    *x24*x31        
     8  +coeff(35)    *x22    *x43    
      y_m35_0_3   =y_m35_0_3   
     9  +coeff(36)    *x22*x31    *x52
     1  +coeff(37)    *x21*x31    *x53
     2  +coeff(38)    *x23*x32*x41    
     3  +coeff(39)    *x24*x31    *x51
     4  +coeff(40)    *x24    *x41*x51
     5  +coeff(41)    *x24    *x43    
     6  +coeff(42)    *x22*x31*x44    
     7  +coeff(43)    *x24*x31*x42*x51
     8  +coeff(44)    *x24    *x43*x51
      y_m35_0_3   =y_m35_0_3   
     9  +coeff(45)*x11*x24*x31*x42    
     1  +coeff(46)    *x24*x32*x43    
     2  +coeff(47)    *x22*x33    *x54
     3  +coeff(48)    *x24*x31*x44*x53
     4  +coeff(49)    *x22*x33*x44*x53
     5  +coeff(50)*x11*x23*x31*x42*x54
c
      return
      end
      function p_m35_0_3   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.99979E-04,-0.21093E+00,-0.59989E-01,-0.36593E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36593E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.19684903E-01, 0.66846140E-01,-0.80129210E-03, 0.41101410E-04,
     + -0.32214804E-02,-0.48829033E-02,-0.54079282E-03,-0.12264580E-02,
     +  0.21051142E-03, 0.12295870E-03, 0.74963892E-04,-0.10595560E-04,
     +  0.52218050E-03, 0.81694324E-03,-0.17037880E-02,-0.28050090E-02,
     + -0.66748521E-04,-0.14065533E-04, 0.53628550E-04, 0.14429932E-03,
     + -0.11530743E-03, 0.89240034E-04, 0.31298593E-03, 0.37597670E-03,
     +  0.81232143E-03, 0.41346703E-03, 0.66060700E-03, 0.15544261E-03,
     + -0.84360333E-04, 0.71612080E-04,-0.27145710E-03, 0.39576030E-04,
     + -0.37917170E-04,-0.70609571E-03,-0.47187593E-04, 0.16530281E-03,
     +  0.23603460E-03, 0.53875774E-04, 0.26317130E-03, 0.14602413E-03,
     +  0.14316922E-03,-0.10205684E-03,-0.12831603E-03,-0.38910450E-04,
     +  0.16851600E-03, 0.39532760E-03,-0.10245330E-03,-0.24187384E-03,
     +  0.65279244E-04,-0.33344150E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_m35_0_3   =avdat
     1  +coeff( 1)        *x31        
     2  +coeff( 2)            *x41    
     3  +coeff( 3)    *x21*x31        
     4  +coeff( 4)    *x21    *x41    
     5  +coeff( 5)        *x31    *x51
     6  +coeff( 6)            *x41*x51
     7  +coeff( 7)    *x22*x31        
     8  +coeff( 8)    *x22    *x41    
      p_m35_0_3   =p_m35_0_3   
     9  +coeff( 9)        *x31*x42    
     1  +coeff(10)            *x43    
     2  +coeff(11)    *x21*x31    *x51
     3  +coeff(12)    *x21    *x41*x51
     4  +coeff(13)        *x31    *x52
     5  +coeff(14)            *x41*x52
     6  +coeff(15)    *x23*x31        
     7  +coeff(16)    *x23    *x41    
     8  +coeff(17)    *x21*x31*x42    
      p_m35_0_3   =p_m35_0_3   
     9  +coeff(18)    *x21*x32*x43    
     1  +coeff(19)    *x22*x31    *x51
     2  +coeff(20)    *x22    *x41*x51
     3  +coeff(21)            *x41*x53
     4  +coeff(22)    *x22*x32*x41    
     5  +coeff(23)    *x23*x31    *x51
     6  +coeff(24)    *x23    *x41*x51
     7  +coeff(25)    *x23*x31*x42    
     8  +coeff(26)    *x23    *x43    
      p_m35_0_3   =p_m35_0_3   
     9  +coeff(27)    *x24*x31*x42    
     1  +coeff(28)    *x23*x34*x41    
     2  +coeff(29)        *x31    *x53
     3  +coeff(30)    *x24    *x43    
     4  +coeff(31)    *x22*x31*x42*x53
     5  +coeff(32)    *x21*x33        
     6  +coeff(33)    *x21*x31    *x52
     7  +coeff(34)    *x22*x31*x42    
     8  +coeff(35)        *x31*x44    
      p_m35_0_3   =p_m35_0_3   
     9  +coeff(36)    *x21*x32*x41*x51
     1  +coeff(37)    *x22*x31    *x52
     2  +coeff(38)    *x21*x31    *x53
     3  +coeff(39)    *x23*x32*x41    
     4  +coeff(40)    *x24*x31    *x51
     5  +coeff(41)    *x24    *x41*x51
     6  +coeff(42)    *x22    *x43*x51
     7  +coeff(43)    *x22    *x41*x53
     8  +coeff(44)*x11*x22*x31*x42    
      p_m35_0_3   =p_m35_0_3   
     9  +coeff(45)    *x22*x32*x43    
     1  +coeff(46)    *x22*x31*x44    
     2  +coeff(47)    *x23*x33    *x51
     3  +coeff(48)    *x21*x34*x41*x51
     4  +coeff(49)*x11*x21*x31*x42*x51
     5  +coeff(50)    *x24*x31    *x52
c
      return
      end
      function l_m35_0_3   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(25)
      data ncoeff/ 24/
      data avdat/ -0.1632585E-01/
      data xmin/
     1 -0.99979E-04,-0.21093E+00,-0.59989E-01,-0.36593E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36593E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.17520073E-01,-0.40949933E-01,-0.53999174E-01,-0.32164970E-02,
     + -0.23677870E-03,-0.13938002E-02,-0.83230790E-03,-0.17472110E-02,
     + -0.68545562E-03, 0.10298063E-02,-0.95936011E-04, 0.35315533E-04,
     +  0.17096303E-03, 0.17955331E-03,-0.99807102E-05,-0.16845910E-04,
     +  0.77940410E-04, 0.13605270E-03,-0.84434570E-05, 0.30333870E-04,
     +  0.85760290E-04, 0.14690641E-03, 0.14956200E-03, 0.22800362E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_m35_0_3   =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)    *x22            
     4  +coeff( 4)            *x42    
     5  +coeff( 5)                *x51
     6  +coeff( 6)    *x21        *x51
     7  +coeff( 7)        *x31*x41    
     8  +coeff( 8)    *x22        *x51
      l_m35_0_3   =l_m35_0_3   
     9  +coeff( 9)    *x23            
     1  +coeff(10)    *x24            
     2  +coeff(11)        *x32        
     3  +coeff(12)    *x21*x31*x41    
     4  +coeff(13)        *x31*x41*x51
     5  +coeff(14)            *x42*x51
     6  +coeff(15)                *x52
     7  +coeff(16)    *x21    *x42    
     8  +coeff(17)    *x21        *x52
      l_m35_0_3   =l_m35_0_3   
     9  +coeff(18)    *x22        *x52
     1  +coeff(19)*x11                
     2  +coeff(20)        *x32    *x51
     3  +coeff(21)    *x22    *x42    
     4  +coeff(22)    *x23*x31*x41    
     5  +coeff(23)    *x23    *x42    
     6  +coeff(24)    *x21*x32        
c
      return
      end
      function x_m35_0_4   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/ -0.2748773E-02/
      data xmin/
     1 -0.99979E-04,-0.20540E+00,-0.59989E-01,-0.36418E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36418E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.36571384E-02, 0.43636652E+00, 0.51350671E-01,-0.13373111E-01,
     + -0.28049710E-02,-0.33500820E-02, 0.88340021E-01,-0.63099912E-02,
     + -0.16189620E-01,-0.11167920E-01,-0.27804840E-02,-0.65089440E-02,
     +  0.10282340E-01,-0.64117010E-03,-0.12730371E-02,-0.21732244E-02,
     + -0.21775392E-02,-0.13082720E-02,-0.67007791E-03, 0.38952352E-02,
     +  0.10745890E-02, 0.10599510E-02, 0.82997710E-03, 0.14976040E-02,
     + -0.20039654E-02,-0.71599744E-02,-0.76389424E-02,-0.80723980E-03,
     + -0.20067144E-02, 0.28321621E-03, 0.14592292E-02, 0.78668540E-03,
     +  0.13750044E-02,-0.28114201E-03,-0.34777622E-03, 0.11315030E-01,
     +  0.13031650E-02,-0.15761010E-01, 0.73756584E-02,-0.97577571E-02,
     + -0.98985880E-03,-0.25738380E-02,-0.54743431E-03, 0.15300290E-01,
     +  0.23316473E-02,-0.78976643E-03, 0.68600650E-03,-0.26912250E-02,
     +  0.68320570E-03,-0.46858732E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x25 = x24*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_m35_0_4   =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)                *x51
     4  +coeff( 4)    *x22            
     5  +coeff( 5)        *x31*x41    
     6  +coeff( 6)            *x42    
     7  +coeff( 7)    *x21        *x51
     8  +coeff( 8)                *x52
      x_m35_0_4   =x_m35_0_4   
     9  +coeff( 9)    *x23            
     1  +coeff(10)    *x21        *x52
     2  +coeff(11)    *x21*x31*x41    
     3  +coeff(12)    *x21    *x42    
     4  +coeff(13)    *x24            
     5  +coeff(14)    *x24        *x51
     6  +coeff(15)    *x21*x32        
     7  +coeff(16)    *x22*x31*x41    
     8  +coeff(17)    *x22    *x42    
      x_m35_0_4   =x_m35_0_4   
     9  +coeff(18)    *x24*x32        
     1  +coeff(19)        *x32        
     2  +coeff(20)    *x22        *x51
     3  +coeff(21)        *x31*x41*x51
     4  +coeff(22)            *x42*x51
     5  +coeff(23)                *x53
     6  +coeff(24)    *x21        *x53
     7  +coeff(25)    *x25            
     8  +coeff(26)    *x24*x31*x41    
      x_m35_0_4   =x_m35_0_4   
     9  +coeff(27)    *x24    *x42    
     1  +coeff(28)    *x22        *x52
     2  +coeff(29)    *x21*x31*x43    
     3  +coeff(30)        *x32    *x51
     4  +coeff(31)    *x23        *x51
     5  +coeff(32)    *x21*x31*x41*x51
     6  +coeff(33)    *x21    *x42*x51
     7  +coeff(34)        *x31*x41*x52
     8  +coeff(35)            *x42*x52
      x_m35_0_4   =x_m35_0_4   
     9  +coeff(36)    *x23    *x42    
     1  +coeff(37)    *x21    *x44    
     2  +coeff(38)    *x25    *x42    
     3  +coeff(39)    *x23*x31*x43    
     4  +coeff(40)    *x23    *x44    
     5  +coeff(41)    *x23*x31*x41*x52
     6  +coeff(42)    *x25    *x42*x51
     7  +coeff(43)    *x25*x33*x41    
     8  +coeff(44)    *x25    *x44    
      x_m35_0_4   =x_m35_0_4   
     9  +coeff(45)    *x25*x32*x42*x51
     1  +coeff(46)    *x22*x32        
     2  +coeff(47)    *x23*x32        
     3  +coeff(48)    *x23*x31*x41    
     4  +coeff(49)    *x21*x32*x42    
     5  +coeff(50)    *x21    *x42*x52
c
      return
      end
      function t_m35_0_4   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/ -0.4158733E-02/
      data xmin/
     1 -0.99979E-04,-0.20540E+00,-0.59989E-01,-0.36418E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36418E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.59655490E-02,-0.11905990E+00, 0.19141280E-01,-0.67132404E-02,
     + -0.40298710E-03,-0.19298591E-02,-0.24410790E-02, 0.35301472E-01,
     + -0.19243551E-02,-0.13032370E-04,-0.52474443E-02,-0.85302790E-03,
     + -0.16178060E-02,-0.18165591E-02, 0.15880831E-02, 0.48419460E-03,
     +  0.71516470E-03,-0.38121123E-02, 0.45187813E-02,-0.16032170E-02,
     +  0.46309660E-03, 0.15235590E-02, 0.10946774E-02, 0.14364353E-02,
     +  0.38804090E-03,-0.57228184E-02,-0.85908530E-02,-0.21404924E-02,
     + -0.19876380E-02, 0.34254320E-03,-0.77778050E-03, 0.31378932E-03,
     +  0.37965192E-02, 0.12174180E-01, 0.92142050E-02,-0.93418400E-02,
     + -0.34717580E-02,-0.13224522E-02,-0.44292660E-04, 0.17271464E-03,
     + -0.30740870E-03,-0.15103740E-02, 0.37052313E-03, 0.50789040E-03,
     + -0.21784940E-02,-0.15439400E-02, 0.23665523E-02,-0.13076522E-02,
     + -0.19776180E-02, 0.14766490E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_m35_0_4   =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)                *x51
     4  +coeff( 4)    *x22            
     5  +coeff( 5)        *x32        
     6  +coeff( 6)        *x31*x41    
     7  +coeff( 7)            *x42    
     8  +coeff( 8)    *x21        *x51
      t_m35_0_4   =t_m35_0_4   
     9  +coeff( 9)                *x52
     1  +coeff(10)*x11*x21            
     2  +coeff(11)    *x23            
     3  +coeff(12)    *x21*x32        
     4  +coeff(13)    *x21*x31*x41    
     5  +coeff(14)    *x21    *x42    
     6  +coeff(15)    *x22        *x51
     7  +coeff(16)        *x31*x41*x51
     8  +coeff(17)            *x42*x51
      t_m35_0_4   =t_m35_0_4   
     9  +coeff(18)    *x21        *x52
     1  +coeff(19)    *x24            
     2  +coeff(20)    *x22*x31*x41    
     3  +coeff(21)    *x22    *x42    
     4  +coeff(22)    *x23        *x51
     5  +coeff(23)    *x21*x31*x41*x51
     6  +coeff(24)    *x21    *x42*x51
     7  +coeff(25)    *x21        *x53
     8  +coeff(26)    *x23*x31*x41    
      t_m35_0_4   =t_m35_0_4   
     9  +coeff(27)    *x23    *x42    
     1  +coeff(28)    *x21*x31*x43    
     2  +coeff(29)    *x21    *x44    
     3  +coeff(30)    *x22*x32    *x51
     4  +coeff(31)    *x24*x32        
     5  +coeff(32)    *x24*x31*x41    
     6  +coeff(33)    *x23*x32*x42    
     7  +coeff(34)    *x23*x31*x43    
     8  +coeff(35)    *x23    *x44    
      t_m35_0_4   =t_m35_0_4   
     9  +coeff(36)    *x24*x31*x43    
     1  +coeff(37)    *x24    *x44    
     2  +coeff(38)    *x24    *x42*x52
     3  +coeff(39)*x11                
     4  +coeff(40)                *x53
     5  +coeff(41)        *x31*x41*x52
     6  +coeff(42)    *x22    *x42*x51
     7  +coeff(43)        *x32*x42*x51
     8  +coeff(44)        *x31*x43*x51
      t_m35_0_4   =t_m35_0_4   
     9  +coeff(45)    *x24    *x42    
     1  +coeff(46)    *x22*x32*x42    
     2  +coeff(47)    *x22*x31*x43    
     3  +coeff(48)    *x23*x31*x41*x51
     4  +coeff(49)    *x23    *x42*x51
     5  +coeff(50)    *x22    *x44*x51
c
      return
      end
      function y_m35_0_4   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.99979E-04,-0.20540E+00,-0.59989E-01,-0.36418E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36418E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.16222750E+00, 0.38445773E+00,-0.33426871E-02, 0.27475870E-02,
     + -0.18557474E-01,-0.31495761E-01,-0.12204864E-02,-0.17191000E-02,
     +  0.90538104E-03, 0.15795210E-02, 0.17119860E-02, 0.59015920E-03,
     +  0.19906443E-03, 0.35862550E-02, 0.58231372E-02,-0.12419330E-01,
     + -0.23211814E-01,-0.26533260E-02,-0.11656040E-02, 0.72072534E-02,
     +  0.12721230E-01, 0.78718671E-02,-0.11091960E-02, 0.16568813E-02,
     +  0.21733670E-02,-0.12326522E-02,-0.53086474E-04,-0.37424480E-03,
     + -0.42977430E-03,-0.69018220E-03, 0.56706930E-03,-0.46687810E-02,
     + -0.53022970E-02,-0.60604010E-02, 0.15218873E-03,-0.17345551E-03,
     + -0.22018060E-02,-0.73185320E-03,-0.30103651E-03,-0.82585064E-03,
     + -0.82185060E-03, 0.36766610E-02, 0.46894961E-03, 0.71179070E-03,
     +  0.31482134E-03, 0.58876930E-03,-0.33256530E-03, 0.13169083E-02,
     +  0.90366040E-03, 0.11208160E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      y_m35_0_4   =avdat
     1  +coeff( 1)        *x31        
     2  +coeff( 2)            *x41    
     3  +coeff( 3)    *x21*x31        
     4  +coeff( 4)    *x21    *x41    
     5  +coeff( 5)        *x31    *x51
     6  +coeff( 6)            *x41*x51
     7  +coeff( 7)    *x22*x31        
     8  +coeff( 8)    *x22    *x41    
      y_m35_0_4   =y_m35_0_4   
     9  +coeff( 9)        *x32*x41    
     1  +coeff(10)        *x31*x42    
     2  +coeff(11)            *x43    
     3  +coeff(12)    *x21*x31    *x51
     4  +coeff(13)    *x21    *x41*x51
     5  +coeff(14)        *x31    *x52
     6  +coeff(15)            *x41*x52
     7  +coeff(16)    *x23*x31        
     8  +coeff(17)    *x23    *x41    
      y_m35_0_4   =y_m35_0_4   
     9  +coeff(18)    *x21*x31*x42    
     1  +coeff(19)    *x21    *x43    
     2  +coeff(20)    *x23*x32*x41    
     3  +coeff(21)    *x23*x31*x42    
     4  +coeff(22)    *x23    *x43    
     5  +coeff(23)            *x41*x53
     6  +coeff(24)    *x23*x31    *x51
     7  +coeff(25)    *x23    *x41*x51
     8  +coeff(26)    *x24    *x41*x51
      y_m35_0_4   =y_m35_0_4   
     9  +coeff(27)        *x31*x42*x53
     1  +coeff(28)        *x32*x41*x51
     2  +coeff(29)            *x43*x51
     3  +coeff(30)        *x31    *x53
     4  +coeff(31)    *x21*x32*x43    
     5  +coeff(32)    *x22*x31*x44    
     6  +coeff(33)    *x24*x32*x43    
     7  +coeff(34)    *x24*x33*x44    
     8  +coeff(35)        *x33        
      y_m35_0_4   =y_m35_0_4   
     9  +coeff(36)    *x21*x33        
     1  +coeff(37)    *x21*x32*x41    
     2  +coeff(38)        *x31*x42*x51
     3  +coeff(39)    *x21*x31    *x52
     4  +coeff(40)    *x24*x31        
     5  +coeff(41)    *x24    *x41    
     6  +coeff(42)    *x22*x31*x42    
     7  +coeff(43)    *x22    *x43    
     8  +coeff(44)        *x31*x44    
      y_m35_0_4   =y_m35_0_4   
     9  +coeff(45)    *x21*x31*x42*x51
     1  +coeff(46)    *x22    *x41*x52
     2  +coeff(47)    *x21    *x41*x53
     3  +coeff(48)    *x23*x33        
     4  +coeff(49)    *x21*x34*x41    
     5  +coeff(50)    *x21*x31*x44    
c
      return
      end
      function p_m35_0_4   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.99979E-04,-0.20540E+00,-0.59989E-01,-0.36418E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36418E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.49701390E-01, 0.13121970E+00,-0.25051943E-02,-0.19678400E-02,
     + -0.98041370E-02,-0.17838220E-01,-0.23018294E-02,-0.44354690E-02,
     +  0.67002001E-04, 0.10190700E-03, 0.48609280E-03, 0.11362050E-03,
     + -0.16906764E-03, 0.20511024E-02, 0.34858242E-02,-0.53271000E-02,
     + -0.11454722E-01,-0.56838232E-03, 0.34891450E-04,-0.55072881E-03,
     + -0.20329924E-04,-0.68947710E-03, 0.77648861E-02, 0.30857373E-02,
     +  0.87204633E-03, 0.26893322E-02, 0.88073620E-03,-0.78085414E-02,
     + -0.42733260E-03, 0.28766444E-02, 0.45732173E-03, 0.10297273E-02,
     +  0.57136170E-03, 0.94345830E-03, 0.43746480E-03, 0.30992180E-02,
     + -0.14437664E-02,-0.43667260E-03,-0.47461101E-02, 0.11943780E-02,
     +  0.53488850E-02,-0.23365570E-03,-0.19989401E-03,-0.19508430E-03,
     + -0.94656721E-03, 0.56354061E-03,-0.65562642E-04,-0.26014490E-03,
     + -0.13097580E-02,-0.18114640E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_m35_0_4   =avdat
     1  +coeff( 1)        *x31        
     2  +coeff( 2)            *x41    
     3  +coeff( 3)    *x21*x31        
     4  +coeff( 4)    *x21    *x41    
     5  +coeff( 5)        *x31    *x51
     6  +coeff( 6)            *x41*x51
     7  +coeff( 7)    *x22*x31        
     8  +coeff( 8)    *x22    *x41    
      p_m35_0_4   =p_m35_0_4   
     9  +coeff( 9)        *x32*x41    
     1  +coeff(10)        *x31*x42    
     2  +coeff(11)            *x43    
     3  +coeff(12)    *x21*x31    *x51
     4  +coeff(13)    *x21    *x41*x51
     5  +coeff(14)        *x31    *x52
     6  +coeff(15)            *x41*x52
     7  +coeff(16)    *x23*x31        
     8  +coeff(17)    *x23    *x41    
      p_m35_0_4   =p_m35_0_4   
     9  +coeff(18)    *x21*x32*x41    
     1  +coeff(19)    *x21*x31*x42    
     2  +coeff(20)    *x21    *x43    
     3  +coeff(21)    *x22    *x41*x51
     4  +coeff(22)            *x41*x53
     5  +coeff(23)    *x22*x31*x42    
     6  +coeff(24)    *x22    *x43    
     7  +coeff(25)        *x31*x44    
     8  +coeff(26)    *x23*x32*x41    
      p_m35_0_4   =p_m35_0_4   
     9  +coeff(27)    *x23*x31*x42    
     1  +coeff(28)    *x22*x31*x44    
     2  +coeff(29)        *x31    *x53
     3  +coeff(30)    *x22*x32*x41    
     4  +coeff(31)        *x32*x43    
     5  +coeff(32)    *x23*x31    *x51
     6  +coeff(33)    *x23    *x41*x51
     7  +coeff(34)    *x21*x31*x42*x51
     8  +coeff(35)    *x22    *x41*x52
      p_m35_0_4   =p_m35_0_4   
     9  +coeff(36)    *x23    *x43    
     1  +coeff(37)    *x21*x31*x44    
     2  +coeff(38)        *x31*x44*x51
     3  +coeff(39)    *x22*x32*x43    
     4  +coeff(40)    *x23    *x43*x51
     5  +coeff(41)    *x23*x31*x44    
     6  +coeff(42)        *x32*x41*x51
     7  +coeff(43)            *x43*x51
     8  +coeff(44)    *x21*x31    *x52
      p_m35_0_4   =p_m35_0_4   
     9  +coeff(45)    *x24    *x41    
     1  +coeff(46)    *x21*x32*x41*x51
     2  +coeff(47)    *x24*x31    *x51
     3  +coeff(48)    *x22*x33    *x51
     4  +coeff(49)    *x24    *x41*x51
     5  +coeff(50)*x12*x22    *x41    
c
      return
      end
      function l_m35_0_4   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(25)
      data ncoeff/ 24/
      data avdat/ -0.2520346E-01/
      data xmin/
     1 -0.99979E-04,-0.20540E+00,-0.59989E-01,-0.36418E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36418E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.26071822E-01,-0.54943930E-01,-0.61380490E-01,-0.16770830E-01,
     +  0.26255060E-02,-0.10823960E-01, 0.70819750E-02,-0.20663540E-02,
     + -0.78579894E-03, 0.44113532E-02, 0.35468910E-02, 0.38009130E-02,
     + -0.24169940E-02,-0.18840220E-02, 0.80531530E-03,-0.19568814E-02,
     + -0.31184890E-03,-0.69984170E-03,-0.89286523E-03, 0.81159680E-03,
     + -0.93917360E-03,-0.97038540E-03, 0.34040960E-03, 0.42507532E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      l_m35_0_4   =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)    *x22            
     4  +coeff( 4)            *x42    
     5  +coeff( 5)                *x51
     6  +coeff( 6)        *x31*x41    
     7  +coeff( 7)    *x21        *x51
     8  +coeff( 8)        *x32        
      l_m35_0_4   =l_m35_0_4   
     9  +coeff( 9)                *x52
     1  +coeff(10)    *x22        *x51
     2  +coeff(11)        *x31*x41*x51
     3  +coeff(12)            *x42*x51
     4  +coeff(13)    *x21        *x52
     5  +coeff(14)    *x23            
     6  +coeff(15)        *x32    *x51
     7  +coeff(16)    *x22        *x52
     8  +coeff(17)    *x21    *x42    
      l_m35_0_4   =l_m35_0_4   
     9  +coeff(18)    *x22*x31*x41    
     1  +coeff(19)    *x22    *x42    
     2  +coeff(20)    *x23        *x51
     3  +coeff(21)        *x31*x41*x52
     4  +coeff(22)            *x42*x52
     5  +coeff(23)                *x53
     6  +coeff(24)    *x24            
c
      return
      end
      function x_m35_0_5   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/ -0.5592477E-02/
      data xmin/
     1 -0.99979E-04,-0.20540E+00,-0.59989E-01,-0.36418E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36418E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.45616601E-03, 0.41491851E+00, 0.55664390E-01,-0.27311360E-01,
     + -0.31292370E-02,-0.37991370E-02, 0.96617333E-01,-0.65262424E-02,
     + -0.16982914E-01,-0.74428622E-02,-0.11129180E-01,-0.40639261E-02,
     +  0.51497210E-02, 0.11695130E-01, 0.88119360E-02,-0.25174252E-02,
     + -0.20067270E-02, 0.32486700E-03,-0.27076851E-02,-0.80013570E-03,
     + -0.35007400E-02, 0.10108240E-02, 0.11782824E-02, 0.74922584E-03,
     + -0.19907124E-02,-0.84417090E-02,-0.93375680E-02,-0.12239180E-02,
     +  0.13711380E-02,-0.25620390E-02, 0.14796982E-02, 0.10138020E-01,
     +  0.26893460E-03, 0.14618592E-02,-0.79197390E-04, 0.19114020E-02,
     +  0.12092354E-01,-0.88376490E-03,-0.72449434E-03,-0.32472580E-02,
     + -0.49639301E-03, 0.15848590E-02,-0.20515720E-01,-0.11115990E-01,
     + -0.10262430E-02, 0.20757742E-01, 0.51218020E-02,-0.78198453E-02,
     + -0.36694790E-02,-0.32481130E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x25 = x24*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x45 = x44*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_m35_0_5   =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)                *x51
     4  +coeff( 4)    *x22            
     5  +coeff( 5)        *x31*x41    
     6  +coeff( 6)            *x42    
     7  +coeff( 7)    *x21        *x51
     8  +coeff( 8)                *x52
      x_m35_0_5   =x_m35_0_5   
     9  +coeff( 9)    *x23            
     1  +coeff(10)    *x21    *x42    
     2  +coeff(11)    *x21        *x52
     3  +coeff(12)    *x23*x31*x41    
     4  +coeff(13)    *x22        *x51
     5  +coeff(14)    *x24            
     6  +coeff(15)    *x25*x32*x42    
     7  +coeff(16)    *x22*x31*x41    
     8  +coeff(17)    *x22    *x42    
      x_m35_0_5   =x_m35_0_5   
     9  +coeff(18)    *x21*x33*x41    
     1  +coeff(19)    *x24*x32        
     2  +coeff(20)        *x32        
     3  +coeff(21)    *x21*x31*x41    
     4  +coeff(22)        *x31*x41*x51
     5  +coeff(23)            *x42*x51
     6  +coeff(24)                *x53
     7  +coeff(25)    *x25            
     8  +coeff(26)    *x24*x31*x41    
      x_m35_0_5   =x_m35_0_5   
     9  +coeff(27)    *x24    *x42    
     1  +coeff(28)    *x21*x32        
     2  +coeff(29)    *x21        *x53
     3  +coeff(30)    *x21*x31*x43    
     4  +coeff(31)    *x21    *x44    
     5  +coeff(32)    *x23*x31*x43    
     6  +coeff(33)        *x32    *x51
     7  +coeff(34)    *x23        *x51
     8  +coeff(35)    *x21*x31*x41*x51
      x_m35_0_5   =x_m35_0_5   
     9  +coeff(36)    *x21    *x42*x51
     1  +coeff(37)    *x23    *x42    
     2  +coeff(38)    *x24        *x51
     3  +coeff(39)    *x22    *x42*x51
     4  +coeff(40)    *x23    *x42*x51
     5  +coeff(41)    *x21*x32*x42*x51
     6  +coeff(42)    *x21*x31*x43*x51
     7  +coeff(43)    *x25    *x42    
     8  +coeff(44)    *x23    *x44    
      x_m35_0_5   =x_m35_0_5   
     9  +coeff(45)    *x22    *x44*x52
     1  +coeff(46)    *x25    *x44    
     2  +coeff(47)    *x23*x32*x44*x51
     3  +coeff(48)    *x25*x32*x44    
     4  +coeff(49)    *x23*x33*x45*x52
     5  +coeff(50)*x11                
c
      return
      end
      function t_m35_0_5   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/ -0.3624947E-02/
      data xmin/
     1 -0.99979E-04,-0.20540E+00,-0.59989E-01,-0.36418E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36418E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.53393380E-02,-0.11724640E+00, 0.19848933E-01,-0.43075140E-02,
     + -0.14559510E-03,-0.10658133E-02,-0.11042590E-02, 0.36047562E-01,
     + -0.18926821E-02,-0.57384741E-05,-0.59564411E-02,-0.89100270E-03,
     + -0.41511822E-02,-0.30143803E-02, 0.20743950E-02,-0.39293081E-02,
     +  0.41562390E-02,-0.15458890E-02,-0.31330380E-02, 0.42852322E-03,
     +  0.10697890E-02, 0.48693473E-03, 0.10081650E-02, 0.54189060E-03,
     +  0.45318960E-03,-0.38439160E-03,-0.51002380E-02,-0.51181230E-03,
     + -0.47711082E-03, 0.59719260E-02,-0.45843833E-04,-0.91596530E-04,
     +  0.13890830E-03,-0.82941940E-03, 0.39254364E-03,-0.18446080E-02,
     + -0.27407620E-02,-0.16682592E-02,-0.23277131E-02, 0.57807060E-03,
     +  0.35220150E-02,-0.67325640E-03, 0.46916920E-03,-0.22203363E-02,
     +  0.14627083E-02,-0.85779593E-03,-0.16822204E-02,-0.89315042E-04,
     +  0.14170860E-03, 0.18540074E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_m35_0_5   =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)                *x51
     4  +coeff( 4)    *x22            
     5  +coeff( 5)        *x32        
     6  +coeff( 6)        *x31*x41    
     7  +coeff( 7)            *x42    
     8  +coeff( 8)    *x21        *x51
      t_m35_0_5   =t_m35_0_5   
     9  +coeff( 9)                *x52
     1  +coeff(10)*x11*x21            
     2  +coeff(11)    *x23            
     3  +coeff(12)    *x21*x32        
     4  +coeff(13)    *x21*x31*x41    
     5  +coeff(14)    *x21    *x42    
     6  +coeff(15)    *x22        *x51
     7  +coeff(16)    *x21        *x52
     8  +coeff(17)    *x24            
      t_m35_0_5   =t_m35_0_5   
     9  +coeff(18)    *x22*x31*x41    
     1  +coeff(19)    *x22    *x42    
     2  +coeff(20)        *x32*x42    
     3  +coeff(21)        *x31*x43    
     4  +coeff(22)            *x44    
     5  +coeff(23)    *x23        *x51
     6  +coeff(24)    *x21*x31*x41*x51
     7  +coeff(25)    *x21    *x42*x51
     8  +coeff(26)    *x22        *x52
      t_m35_0_5   =t_m35_0_5   
     9  +coeff(27)    *x23    *x42    
     1  +coeff(28)    *x21*x31*x43    
     2  +coeff(29)    *x24*x32        
     3  +coeff(30)    *x23    *x44    
     4  +coeff(31)*x11                
     5  +coeff(32)            *x42*x51
     6  +coeff(33)                *x53
     7  +coeff(34)    *x22*x32        
     8  +coeff(35)    *x21        *x53
      t_m35_0_5   =t_m35_0_5   
     9  +coeff(36)    *x21    *x44    
     1  +coeff(37)    *x24*x31*x41    
     2  +coeff(38)    *x24    *x42    
     3  +coeff(39)    *x22*x31*x43    
     4  +coeff(40)    *x23*x34        
     5  +coeff(41)    *x23*x31*x43    
     6  +coeff(42)    *x24    *x42*x51
     7  +coeff(43)    *x22*x31*x43*x51
     8  +coeff(44)    *x23*x31*x41*x52
      t_m35_0_5   =t_m35_0_5   
     9  +coeff(45)    *x21*x31*x43*x52
     1  +coeff(46)    *x22*x31*x41*x53
     2  +coeff(47)    *x24    *x44    
     3  +coeff(48)        *x31*x41*x51
     4  +coeff(49)        *x33*x41    
     5  +coeff(50)    *x21*x32    *x51
c
      return
      end
      function y_m35_0_5   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.99979E-04,-0.20540E+00,-0.59989E-01,-0.36418E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36418E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.17210990E+00, 0.41041290E+00, 0.28814040E-03, 0.14751180E-01,
     + -0.19834440E-01,-0.33476360E-01,-0.19496240E-02,-0.31359542E-02,
     +  0.81396044E-03, 0.17066943E-02, 0.20041820E-02, 0.85571920E-03,
     +  0.99496640E-03, 0.37647870E-02, 0.60143033E-02,-0.12800562E-01,
     + -0.26007963E-01,-0.95193093E-03,-0.86638530E-03, 0.57900130E-02,
     +  0.12760830E-01, 0.85466310E-02,-0.11805080E-02, 0.16373570E-02,
     +  0.22745244E-02,-0.14541460E-02,-0.92952800E-03,-0.72863174E-03,
     + -0.20544221E-02, 0.59969770E-03,-0.10953490E-02, 0.87345822E-03,
     + -0.87548760E-03,-0.87543192E-03,-0.10706521E-01, 0.13993531E-02,
     +  0.10105250E-03, 0.70957641E-03,-0.35528652E-03,-0.16095380E-02,
     +  0.18928461E-02, 0.27474912E-03, 0.47275430E-03, 0.52327814E-03,
     +  0.80111413E-03,-0.15215724E-02,-0.22376930E-03,-0.14416063E-02,
     +  0.66777660E-02,-0.25900410E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      y_m35_0_5   =avdat
     1  +coeff( 1)        *x31        
     2  +coeff( 2)            *x41    
     3  +coeff( 3)    *x21*x31        
     4  +coeff( 4)    *x21    *x41    
     5  +coeff( 5)        *x31    *x51
     6  +coeff( 6)            *x41*x51
     7  +coeff( 7)    *x22*x31        
     8  +coeff( 8)    *x22    *x41    
      y_m35_0_5   =y_m35_0_5   
     9  +coeff( 9)        *x32*x41    
     1  +coeff(10)        *x31*x42    
     2  +coeff(11)            *x43    
     3  +coeff(12)    *x21*x31    *x51
     4  +coeff(13)    *x21    *x41*x51
     5  +coeff(14)        *x31    *x52
     6  +coeff(15)            *x41*x52
     7  +coeff(16)    *x23*x31        
     8  +coeff(17)    *x23    *x41    
      y_m35_0_5   =y_m35_0_5   
     9  +coeff(18)    *x21*x31*x42    
     1  +coeff(19)    *x21    *x43    
     2  +coeff(20)    *x23*x32*x41    
     3  +coeff(21)    *x23*x31*x42    
     4  +coeff(22)    *x23    *x43    
     5  +coeff(23)            *x41*x53
     6  +coeff(24)    *x23*x31    *x51
     7  +coeff(25)    *x23    *x41*x51
     8  +coeff(26)    *x24    *x41*x51
      y_m35_0_5   =y_m35_0_5   
     9  +coeff(27)        *x31*x42*x51
     1  +coeff(28)        *x31    *x53
     2  +coeff(29)    *x24    *x41    
     3  +coeff(30)        *x31*x44    
     4  +coeff(31)    *x21*x32*x43    
     5  +coeff(32)    *x22*x31*x42*x51
     6  +coeff(33)        *x32*x43*x51
     7  +coeff(34)    *x23*x31    *x52
     8  +coeff(35)    *x24*x31*x44    
      y_m35_0_5   =y_m35_0_5   
     9  +coeff(36)    *x24*x32*x43*x51
     1  +coeff(37)        *x33        
     2  +coeff(38)    *x21*x33        
     3  +coeff(39)    *x21    *x41*x52
     4  +coeff(40)    *x24*x31        
     5  +coeff(41)    *x22*x32*x41    
     6  +coeff(42)    *x22*x31*x42    
     7  +coeff(43)    *x21*x31*x42*x51
     8  +coeff(44)    *x21    *x43*x51
      y_m35_0_5   =y_m35_0_5   
     9  +coeff(45)    *x22    *x41*x52
     1  +coeff(46)    *x21*x33*x42    
     2  +coeff(47)*x12*x22    *x41    
     3  +coeff(48)    *x24*x32*x41    
     4  +coeff(49)    *x24*x31*x42    
     5  +coeff(50)    *x22*x32*x43    
c
      return
      end
      function p_m35_0_5   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.99979E-04,-0.20540E+00,-0.59989E-01,-0.36418E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36418E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.48568762E-01, 0.12837864E+00,-0.46048741E-02,-0.79172990E-02,
     + -0.95431630E-02,-0.16823610E-01,-0.13530943E-02,-0.16845570E-02,
     +  0.52993630E-03, 0.82366191E-03, 0.11579200E-02, 0.93798390E-03,
     +  0.13815731E-02, 0.17494470E-02, 0.32893014E-02,-0.56763454E-02,
     + -0.11391520E-01, 0.29082620E-03, 0.20325490E-03, 0.82684872E-03,
     + -0.37116944E-03,-0.68209622E-03,-0.19138203E-02, 0.12549700E-02,
     +  0.23714600E-02, 0.32462861E-02,-0.67069380E-03,-0.39916863E-04,
     + -0.26381180E-03, 0.97556610E-03, 0.62605901E-03, 0.19769612E-02,
     + -0.25480020E-03, 0.23114380E-02, 0.36136650E-02,-0.46566390E-02,
     +  0.82272880E-03,-0.11552303E-02,-0.22722550E-03,-0.20109260E-03,
     + -0.16613201E-02, 0.26234440E-02,-0.46173582E-03, 0.42875774E-03,
     +  0.32718510E-03, 0.20218810E-03, 0.66229264E-03,-0.74392324E-03,
     +  0.66123770E-03,-0.45202210E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_m35_0_5   =avdat
     1  +coeff( 1)        *x31        
     2  +coeff( 2)            *x41    
     3  +coeff( 3)    *x21*x31        
     4  +coeff( 4)    *x21    *x41    
     5  +coeff( 5)        *x31    *x51
     6  +coeff( 6)            *x41*x51
     7  +coeff( 7)    *x22*x31        
     8  +coeff( 8)    *x22    *x41    
      p_m35_0_5   =p_m35_0_5   
     9  +coeff( 9)        *x32*x41    
     1  +coeff(10)        *x31*x42    
     2  +coeff(11)            *x43    
     3  +coeff(12)    *x21*x31    *x51
     4  +coeff(13)    *x21    *x41*x51
     5  +coeff(14)        *x31    *x52
     6  +coeff(15)            *x41*x52
     7  +coeff(16)    *x23*x31        
     8  +coeff(17)    *x23    *x41    
      p_m35_0_5   =p_m35_0_5   
     9  +coeff(18)    *x21*x32*x41    
     1  +coeff(19)    *x21*x31*x42    
     2  +coeff(20)    *x21    *x43    
     3  +coeff(21)    *x22    *x41*x51
     4  +coeff(22)            *x41*x53
     5  +coeff(23)    *x22    *x43    
     6  +coeff(24)    *x23    *x41*x51
     7  +coeff(25)    *x23*x32*x41    
     8  +coeff(26)    *x23*x31*x42    
      p_m35_0_5   =p_m35_0_5   
     9  +coeff(27)        *x31*x42*x53
     1  +coeff(28)        *x32*x41*x51
     2  +coeff(29)            *x43*x51
     3  +coeff(30)    *x23*x31    *x51
     4  +coeff(31)    *x22    *x41*x52
     5  +coeff(32)    *x23    *x43    
     6  +coeff(33)        *x33    *x53
     7  +coeff(34)    *x24*x31*x42    
     8  +coeff(35)    *x24    *x43    
      p_m35_0_5   =p_m35_0_5   
     9  +coeff(36)    *x22*x31*x44    
     1  +coeff(37)        *x33*x44    
     2  +coeff(38)    *x23*x32*x43    
     3  +coeff(39)    *x21*x31    *x52
     4  +coeff(40)    *x24*x31        
     5  +coeff(41)    *x24    *x41    
     6  +coeff(42)    *x22*x31*x42    
     7  +coeff(43)        *x33*x42    
     8  +coeff(44)        *x31*x44    
      p_m35_0_5   =p_m35_0_5   
     9  +coeff(45)    *x22*x31    *x52
     1  +coeff(46)        *x33    *x52
     2  +coeff(47)    *x21*x31*x44    
     3  +coeff(48)    *x24    *x41*x51
     4  +coeff(49)    *x22    *x43*x51
     5  +coeff(50)        *x32*x43*x51
c
      return
      end
      function l_m35_0_5   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(25)
      data ncoeff/ 24/
      data avdat/ -0.2544361E-01/
      data xmin/
     1 -0.99979E-04,-0.20540E+00,-0.59989E-01,-0.36418E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36418E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.27813260E-01,-0.15805742E+00,-0.57700920E-01,-0.10425490E-01,
     + -0.17694810E-01,-0.15431933E-01,-0.11503500E-01,-0.20202080E-02,
     +  0.45082913E-02, 0.35928250E-02, 0.38519084E-02, 0.18220360E-02,
     + -0.19018110E-02,-0.20270671E-02, 0.71238481E-03, 0.69485820E-03,
     + -0.23823330E-03, 0.76441770E-03, 0.11560890E-02, 0.97707090E-03,
     +  0.77944440E-03,-0.86980100E-03,-0.87907930E-03, 0.12347510E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_m35_0_5   =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)    *x22            
     4  +coeff( 4)                *x51
     5  +coeff( 5)            *x42    
     6  +coeff( 6)    *x21        *x51
     7  +coeff( 7)        *x31*x41    
     8  +coeff( 8)        *x32        
      l_m35_0_5   =l_m35_0_5   
     9  +coeff( 9)    *x22        *x51
     1  +coeff(10)        *x31*x41*x51
     2  +coeff(11)            *x42*x51
     3  +coeff(12)    *x23            
     4  +coeff(13)    *x24            
     5  +coeff(14)    *x22        *x52
     6  +coeff(15)                *x52
     7  +coeff(16)    *x21*x31*x41    
     8  +coeff(17)    *x21    *x42    
      l_m35_0_5   =l_m35_0_5   
     9  +coeff(18)        *x32    *x51
     1  +coeff(19)    *x22*x31*x41    
     2  +coeff(20)    *x22    *x42    
     3  +coeff(21)    *x23        *x51
     4  +coeff(22)        *x31*x41*x52
     5  +coeff(23)            *x42*x52
     6  +coeff(24)    *x23    *x42    
c
      return
      end
      function x_m35_0_6   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/ -0.1230867E-01/
      data xmin/
     1 -0.99979E-04,-0.20540E+00,-0.59989E-01,-0.36288E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36288E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.73342323E-02, 0.34039330E+00, 0.72528580E-01,-0.29194210E-01,
     + -0.32651550E-02,-0.39760721E-02, 0.12309404E+00,-0.84018353E-02,
     + -0.19735602E-01,-0.14088972E-01,-0.54513450E-02,-0.10082660E-01,
     +  0.14514231E-01, 0.56251394E-02,-0.67519330E-02,-0.64130090E-02,
     + -0.17240560E-02,-0.23075030E-02,-0.40978654E-02,-0.79027570E-03,
     +  0.85687282E-03, 0.61866210E-03, 0.91993550E-03, 0.16262040E-02,
     + -0.49345370E-03,-0.25404340E-02,-0.72937180E-02,-0.85489080E-02,
     +  0.83298441E-02, 0.28910000E-03, 0.21114830E-02, 0.10706390E-02,
     +  0.19530071E-02, 0.55318833E-02,-0.73842732E-02, 0.11365820E-01,
     + -0.50804340E-02,-0.40798230E-02, 0.12143204E-02,-0.66695341E-04,
     + -0.15429840E-02, 0.62367730E-03, 0.39165630E-03, 0.43262424E-03,
     + -0.54186992E-02,-0.21147940E-03, 0.23851641E-02, 0.45114670E-03,
     + -0.24585620E-02,-0.73465623E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x25 = x24*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_m35_0_6   =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)                *x51
     4  +coeff( 4)    *x22            
     5  +coeff( 5)        *x31*x41    
     6  +coeff( 6)            *x42    
     7  +coeff( 7)    *x21        *x51
     8  +coeff( 8)                *x52
      x_m35_0_6   =x_m35_0_6   
     9  +coeff( 9)    *x23            
     1  +coeff(10)    *x21        *x52
     2  +coeff(11)    *x21*x31*x41    
     3  +coeff(12)    *x21    *x42    
     4  +coeff(13)    *x24            
     5  +coeff(14)    *x22        *x51
     6  +coeff(15)    *x22*x31*x41    
     7  +coeff(16)    *x22    *x42    
     8  +coeff(17)    *x24*x32        
      x_m35_0_6   =x_m35_0_6   
     9  +coeff(18)    *x21*x32        
     1  +coeff(19)    *x25            
     2  +coeff(20)        *x32        
     3  +coeff(21)        *x31*x41*x51
     4  +coeff(22)            *x42*x51
     5  +coeff(23)                *x53
     6  +coeff(24)    *x21        *x53
     7  +coeff(25)    *x23    *x42    
     8  +coeff(26)    *x21*x31*x43    
      x_m35_0_6   =x_m35_0_6   
     9  +coeff(27)    *x24*x31*x41    
     1  +coeff(28)    *x24    *x42    
     2  +coeff(29)    *x23    *x44    
     3  +coeff(30)        *x32    *x51
     4  +coeff(31)    *x23        *x51
     5  +coeff(32)    *x21*x31*x41*x51
     6  +coeff(33)    *x21    *x42*x51
     7  +coeff(34)    *x21*x32*x42    
     8  +coeff(35)    *x25    *x42    
      x_m35_0_6   =x_m35_0_6   
     9  +coeff(36)    *x23*x31*x43    
     1  +coeff(37)    *x21*x32*x44    
     2  +coeff(38)    *x24    *x42*x51
     3  +coeff(39)    *x25*x34        
     4  +coeff(40)*x11                
     5  +coeff(41)    *x22*x32        
     6  +coeff(42)        *x31*x43    
     7  +coeff(43)            *x44    
     8  +coeff(44)    *x21*x32    *x51
      x_m35_0_6   =x_m35_0_6   
     9  +coeff(45)    *x23*x31*x41    
     1  +coeff(46)    *x24        *x51
     2  +coeff(47)    *x22    *x42*x51
     3  +coeff(48)        *x34*x42    
     4  +coeff(49)    *x23    *x42*x51
     5  +coeff(50)    *x22    *x42*x52
c
      return
      end
      function t_m35_0_6   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/ -0.3432297E-02/
      data xmin/
     1 -0.99979E-04,-0.20540E+00,-0.59989E-01,-0.36288E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36288E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.46534854E-02,-0.87302620E-01, 0.28828712E-01,-0.52328963E-04,
     + -0.53056760E-02,-0.54649710E-04,-0.40629322E-03,-0.33813534E-03,
     +  0.40798970E-01,-0.38559141E-02, 0.58073801E-04,-0.69589214E-02,
     + -0.83198462E-03,-0.23289170E-02,-0.26838120E-02, 0.15085454E-04,
     +  0.33429381E-02, 0.29429040E-04, 0.13321350E-03,-0.55572640E-02,
     +  0.54238911E-03, 0.43875183E-02,-0.56801182E-02,-0.68717682E-02,
     +  0.36578334E-03, 0.22153210E-03, 0.15234230E-02,-0.11813520E-03,
     + -0.29256060E-03,-0.67380780E-03, 0.76345034E-03,-0.14360790E-03,
     + -0.43900140E-02, 0.62259690E-03,-0.64205713E-02, 0.30701440E-02,
     +  0.23688410E-02, 0.83349272E-03, 0.33339940E-03, 0.36994280E-03,
     +  0.53949132E-02, 0.46413570E-02,-0.21036210E-02,-0.16394570E-02,
     +  0.97322610E-03, 0.39945500E-03, 0.55626200E-03,-0.79817030E-03,
     + -0.21069720E-03, 0.34358900E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_m35_0_6   =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)                *x51
     4  +coeff( 4)*x11                
     5  +coeff( 5)    *x22            
     6  +coeff( 6)        *x32        
     7  +coeff( 7)        *x31*x41    
     8  +coeff( 8)            *x42    
      t_m35_0_6   =t_m35_0_6   
     9  +coeff( 9)    *x21        *x51
     1  +coeff(10)                *x52
     2  +coeff(11)*x11*x21            
     3  +coeff(12)    *x23            
     4  +coeff(13)    *x21*x32        
     5  +coeff(14)    *x21*x31*x41    
     6  +coeff(15)    *x21    *x42    
     7  +coeff(16)*x11            *x51
     8  +coeff(17)    *x22        *x51
      t_m35_0_6   =t_m35_0_6   
     9  +coeff(18)        *x31*x41*x51
     1  +coeff(19)            *x42*x51
     2  +coeff(20)    *x21        *x52
     3  +coeff(21)                *x53
     4  +coeff(22)    *x24            
     5  +coeff(23)    *x22*x31*x41    
     6  +coeff(24)    *x22    *x42    
     7  +coeff(25)        *x31*x43    
     8  +coeff(26)            *x44    
      t_m35_0_6   =t_m35_0_6   
     9  +coeff(27)    *x23        *x51
     1  +coeff(28)    *x21*x31*x41*x51
     2  +coeff(29)    *x21    *x42*x51
     3  +coeff(30)    *x22        *x52
     4  +coeff(31)    *x21        *x53
     5  +coeff(32)*x11*x23            
     6  +coeff(33)    *x23*x31*x41    
     7  +coeff(34)    *x21*x33*x41    
     8  +coeff(35)    *x23    *x42    
      t_m35_0_6   =t_m35_0_6   
     9  +coeff(36)    *x21*x32*x42    
     1  +coeff(37)    *x21*x31*x43    
     2  +coeff(38)    *x21    *x44    
     3  +coeff(39)    *x22*x32    *x51
     4  +coeff(40)    *x24*x32        
     5  +coeff(41)    *x23*x31*x43    
     6  +coeff(42)    *x23    *x44    
     7  +coeff(43)    *x24    *x42*x51
     8  +coeff(44)    *x22*x32        
      t_m35_0_6   =t_m35_0_6   
     9  +coeff(45)    *x22    *x42*x51
     1  +coeff(46)        *x32*x42*x51
     2  +coeff(47)        *x31*x43*x51
     3  +coeff(48)    *x22*x31*x43    
     4  +coeff(49)            *x44*x52
     5  +coeff(50)*x12*x21*x31*x41    
c
      return
      end
      function y_m35_0_6   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.99979E-04,-0.20540E+00,-0.59989E-01,-0.36288E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36288E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.20146684E+00, 0.48725512E+00,-0.26307380E-02, 0.97484740E-02,
     + -0.25281464E-01,-0.42490770E-01,-0.39457152E-02,-0.53160991E-02,
     +  0.16352750E-02, 0.11160300E-02, 0.13580560E-02, 0.45977532E-02,
     +  0.77951923E-02,-0.15847420E-01,-0.32926954E-01,-0.12710143E-02,
     +  0.33762780E-03, 0.93246320E-02, 0.71400240E-02, 0.12465020E-01,
     +  0.91241070E-02, 0.97955640E-03,-0.13738990E-02,-0.14975764E-02,
     + -0.20021741E-03, 0.29766990E-02, 0.31554163E-02,-0.11644810E-01,
     + -0.86456350E-03,-0.41103560E-02,-0.88785640E-03,-0.94292000E-03,
     +  0.25514570E-02,-0.69179292E-02, 0.35316400E-02, 0.87527593E-03,
     +  0.48021334E-03,-0.10302223E-02, 0.23081580E-02, 0.10344092E-02,
     + -0.11682043E-05, 0.35900410E-02, 0.96458810E-03,-0.22810960E-02,
     + -0.81960790E-03, 0.11642590E-02,-0.29005570E-02, 0.42652813E-02,
     + -0.33222900E-02, 0.10824190E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      y_m35_0_6   =avdat
     1  +coeff( 1)        *x31        
     2  +coeff( 2)            *x41    
     3  +coeff( 3)    *x21*x31        
     4  +coeff( 4)    *x21    *x41    
     5  +coeff( 5)        *x31    *x51
     6  +coeff( 6)            *x41*x51
     7  +coeff( 7)    *x22*x31        
     8  +coeff( 8)    *x22    *x41    
      y_m35_0_6   =y_m35_0_6   
     9  +coeff( 9)            *x43    
     1  +coeff(10)    *x21*x31    *x51
     2  +coeff(11)    *x21    *x41*x51
     3  +coeff(12)        *x31    *x52
     4  +coeff(13)            *x41*x52
     5  +coeff(14)    *x23*x31        
     6  +coeff(15)    *x23    *x41    
     7  +coeff(16)    *x21*x31*x42    
     8  +coeff(17)    *x21    *x43    
      y_m35_0_6   =y_m35_0_6   
     9  +coeff(18)    *x22*x31*x42    
     1  +coeff(19)    *x23*x32*x41    
     2  +coeff(20)    *x23*x31*x42    
     3  +coeff(21)    *x23    *x43    
     4  +coeff(22)        *x31*x42    
     5  +coeff(23)    *x22    *x41*x51
     6  +coeff(24)            *x41*x53
     7  +coeff(25)        *x34*x41    
     8  +coeff(26)    *x23*x31    *x51
      y_m35_0_6   =y_m35_0_6   
     9  +coeff(27)    *x23    *x41*x51
     1  +coeff(28)    *x22*x31*x44    
     2  +coeff(29)        *x31    *x53
     3  +coeff(30)    *x24    *x41    
     4  +coeff(31)    *x22    *x41*x52
     5  +coeff(32)    *x23*x31    *x52
     6  +coeff(33)    *x24    *x43    
     7  +coeff(34)    *x22*x32*x43    
     8  +coeff(35)    *x24    *x41*x52
      y_m35_0_6   =y_m35_0_6   
     9  +coeff(36)        *x32*x41    
     1  +coeff(37)    *x21*x33        
     2  +coeff(38)    *x24*x31        
     3  +coeff(39)    *x22*x32*x41    
     4  +coeff(40)        *x31*x44    
     5  +coeff(41)    *x22*x31    *x52
     6  +coeff(42)    *x21*x33*x42    
     7  +coeff(43)    *x21*x31*x44    
     8  +coeff(44)    *x22*x33*x42    
      y_m35_0_6   =y_m35_0_6   
     9  +coeff(45)    *x23*x33    *x51
     1  +coeff(46)    *x22*x33    *x52
     2  +coeff(47)    *x23*x33*x42    
     3  +coeff(48)    *x23*x31*x44    
     4  +coeff(49)    *x21*x33*x44    
     5  +coeff(50)    *x22    *x43*x53
c
      return
      end
      function p_m35_0_6   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.99979E-04,-0.20540E+00,-0.59989E-01,-0.36288E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36288E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.32073691E-01, 0.88781811E-01,-0.36397650E-02,-0.57799024E-02,
     + -0.51655611E-02,-0.82644550E-02, 0.18449631E-04,-0.28372240E-02,
     + -0.60918820E-02,-0.79409911E-03,-0.18822120E-02,-0.19536040E-02,
     + -0.17479382E-03,-0.85752680E-03, 0.83589123E-03, 0.13535680E-02,
     + -0.38438520E-02,-0.79933630E-02, 0.13574490E-02, 0.39902003E-02,
     +  0.21222670E-02,-0.66476640E-03, 0.20234654E-02, 0.15513480E-02,
     +  0.32605012E-03,-0.13805080E-02, 0.37827450E-02, 0.20526430E-02,
     + -0.15506354E-02,-0.12176480E-02, 0.72045833E-03, 0.16049850E-03,
     + -0.34297250E-04,-0.28714642E-03,-0.79705710E-03, 0.26639233E-03,
     +  0.87777330E-03, 0.23753651E-03,-0.16028310E-03,-0.21651780E-03,
     +  0.16888723E-02, 0.41414101E-03, 0.48671753E-03, 0.66561350E-03,
     +  0.44086530E-03,-0.70864440E-03,-0.66940040E-03,-0.12153490E-02,
     +  0.69956380E-03, 0.15555090E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_m35_0_6   =avdat
     1  +coeff( 1)        *x31        
     2  +coeff( 2)            *x41    
     3  +coeff( 3)    *x21*x31        
     4  +coeff( 4)    *x21    *x41    
     5  +coeff( 5)        *x31    *x51
     6  +coeff( 6)            *x41*x51
     7  +coeff( 7)*x11    *x31        
     8  +coeff( 8)    *x22*x31        
      p_m35_0_6   =p_m35_0_6   
     9  +coeff( 9)    *x22    *x41    
     1  +coeff(10)        *x32*x41    
     2  +coeff(11)        *x31*x42    
     3  +coeff(12)            *x43    
     4  +coeff(13)    *x21*x31    *x51
     5  +coeff(14)    *x21    *x41*x51
     6  +coeff(15)        *x31    *x52
     7  +coeff(16)            *x41*x52
     8  +coeff(17)    *x23*x31        
      p_m35_0_6   =p_m35_0_6   
     9  +coeff(18)    *x23    *x41    
     1  +coeff(19)    *x21*x32*x41    
     2  +coeff(20)    *x21*x31*x42    
     3  +coeff(21)    *x21    *x43    
     4  +coeff(22)    *x22    *x41*x51
     5  +coeff(23)        *x31*x42*x51
     6  +coeff(24)            *x43*x51
     7  +coeff(25)    *x21    *x41*x52
     8  +coeff(26)    *x24    *x41    
      p_m35_0_6   =p_m35_0_6   
     9  +coeff(27)    *x22*x31*x42    
     1  +coeff(28)    *x22    *x43    
     2  +coeff(29)        *x32*x43    
     3  +coeff(30)        *x31*x44    
     4  +coeff(31)    *x22    *x41*x52
     5  +coeff(32)    *x22*x32*x41*x51
     6  +coeff(33)        *x34*x41*x51
     7  +coeff(34)    *x24*x32*x41    
     8  +coeff(35)        *x33*x44    
      p_m35_0_6   =p_m35_0_6   
     9  +coeff(36)    *x21*x33        
     1  +coeff(37)        *x32*x41*x51
     2  +coeff(38)    *x21*x31    *x52
     3  +coeff(39)        *x31    *x53
     4  +coeff(40)            *x41*x53
     5  +coeff(41)    *x22*x32*x41    
     6  +coeff(42)    *x23*x31    *x51
     7  +coeff(43)    *x21*x31*x42*x51
     8  +coeff(44)    *x21    *x43*x51
      p_m35_0_6   =p_m35_0_6   
     9  +coeff(45)    *x22*x31    *x52
     1  +coeff(46)        *x31*x42*x52
     2  +coeff(47)            *x43*x52
     3  +coeff(48)    *x21*x31*x44    
     4  +coeff(49)    *x23    *x41*x52
     5  +coeff(50)    *x24    *x43    
c
      return
      end
      function l_m35_0_6   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(25)
      data ncoeff/ 24/
      data avdat/ -0.2555332E-01/
      data xmin/
     1 -0.99979E-04,-0.20540E+00,-0.59989E-01,-0.36288E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36288E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.27870243E-01,-0.14834704E+00,-0.61151050E-01,-0.12624730E-01,
     + -0.22274880E-01,-0.17227630E-01,-0.76805840E-04,-0.14993270E-01,
     +  0.73863110E-02,-0.26106880E-02, 0.44878004E-02, 0.48357130E-02,
     +  0.29503491E-02,-0.29613631E-02, 0.14691520E-02, 0.13551084E-02,
     +  0.93018230E-03,-0.29924700E-02, 0.41098974E-03, 0.33014520E-03,
     + -0.45570910E-03, 0.15472230E-02, 0.14525891E-02, 0.94318593E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_m35_0_6   =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)    *x22            
     4  +coeff( 4)                *x51
     5  +coeff( 5)            *x42    
     6  +coeff( 6)    *x21        *x51
     7  +coeff( 7)        *x33*x41    
     8  +coeff( 8)        *x31*x41    
      l_m35_0_6   =l_m35_0_6   
     9  +coeff( 9)    *x22        *x51
     1  +coeff(10)        *x32        
     2  +coeff(11)        *x31*x41*x51
     3  +coeff(12)            *x42*x51
     4  +coeff(13)    *x23            
     5  +coeff(14)    *x22        *x52
     6  +coeff(15)    *x21*x31*x41    
     7  +coeff(16)    *x21    *x42    
     8  +coeff(17)        *x32    *x51
      l_m35_0_6   =l_m35_0_6   
     9  +coeff(18)    *x24            
     1  +coeff(19)                *x52
     2  +coeff(20)    *x21*x32        
     3  +coeff(21)    *x21        *x52
     4  +coeff(22)    *x22*x31*x41    
     5  +coeff(23)    *x22    *x42    
     6  +coeff(24)    *x23        *x51
c
      return
      end
      function x_m35_0_7   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/ -0.2131009E-01/
      data xmin/
     1 -0.99979E-04,-0.20540E+00,-0.59989E-01,-0.36252E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36252E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.17189390E-01, 0.27831801E+00, 0.11327281E+00,-0.37559310E-01,
     + -0.54086381E-02,-0.64385524E-02, 0.17192330E+00,-0.14728322E-01,
     + -0.25737820E-01,-0.85898500E-02,-0.22038350E-01,-0.11992672E-01,
     +  0.87897930E-02, 0.20592690E-01,-0.12630410E-01,-0.12764701E-01,
     +  0.17690814E-03,-0.61505693E-02,-0.61876340E-02,-0.11033844E-02,
     + -0.16956890E-02, 0.17342624E-02, 0.15079344E-02, 0.18011650E-02,
     +  0.27533100E-02, 0.30578750E-02,-0.10539340E-01,-0.12497720E-03,
     +  0.38710250E-03,-0.46794750E-02, 0.14111940E-02, 0.10742350E-02,
     + -0.18395513E-02,-0.28252131E-02,-0.94725810E-02,-0.11528160E-01,
     + -0.93599362E-02, 0.83292120E-02, 0.21341683E-01, 0.17559350E-01,
     + -0.83868430E-02, 0.18069151E-02,-0.82495020E-03,-0.68157480E-03,
     + -0.11695382E-02, 0.27815540E-02, 0.72285014E-03,-0.20197401E-02,
     +  0.13845664E-02, 0.23972683E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x25 = x24*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_m35_0_7   =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)                *x51
     4  +coeff( 4)    *x22            
     5  +coeff( 5)        *x31*x41    
     6  +coeff( 6)            *x42    
     7  +coeff( 7)    *x21        *x51
     8  +coeff( 8)                *x52
      x_m35_0_7   =x_m35_0_7   
     9  +coeff( 9)    *x23            
     1  +coeff(10)    *x21    *x42    
     2  +coeff(11)    *x21        *x52
     3  +coeff(12)    *x23*x31*x41    
     4  +coeff(13)    *x22        *x51
     5  +coeff(14)    *x24            
     6  +coeff(15)    *x22*x31*x41    
     7  +coeff(16)    *x22    *x42    
     8  +coeff(17)    *x22*x34        
      x_m35_0_7   =x_m35_0_7   
     9  +coeff(18)    *x21*x31*x41    
     1  +coeff(19)    *x25            
     2  +coeff(20)        *x32        
     3  +coeff(21)    *x21*x32        
     4  +coeff(22)        *x31*x41*x51
     5  +coeff(23)            *x42*x51
     6  +coeff(24)                *x53
     7  +coeff(25)    *x23        *x51
     8  +coeff(26)    *x21        *x53
      x_m35_0_7   =x_m35_0_7   
     9  +coeff(27)    *x23    *x42    
     1  +coeff(28)*x11                
     2  +coeff(29)        *x32    *x51
     3  +coeff(30)    *x22*x32        
     4  +coeff(31)        *x31*x43    
     5  +coeff(32)            *x44    
     6  +coeff(33)    *x21*x31*x43    
     7  +coeff(34)    *x21    *x44    
     8  +coeff(35)    *x24*x31*x41    
      x_m35_0_7   =x_m35_0_7   
     9  +coeff(36)    *x24    *x42    
     1  +coeff(37)    *x25    *x42    
     2  +coeff(38)    *x23*x32*x42    
     3  +coeff(39)    *x23*x31*x43    
     4  +coeff(40)    *x23    *x44    
     5  +coeff(41)    *x24    *x42*x51
     6  +coeff(42)    *x21*x31*x41*x51
     7  +coeff(43)    *x22        *x52
     8  +coeff(44)            *x42*x52
      x_m35_0_7   =x_m35_0_7   
     9  +coeff(45)    *x23*x32        
     1  +coeff(46)    *x22    *x42*x51
     2  +coeff(47)        *x34*x42    
     3  +coeff(48)    *x21*x33*x41*x51
     4  +coeff(49)    *x24*x32    *x51
     5  +coeff(50)    *x22    *x42*x53
c
      return
      end
      function t_m35_0_7   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/ -0.8645021E-02/
      data xmin/
     1 -0.99979E-04,-0.20540E+00,-0.59989E-01,-0.36252E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36252E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.91628963E-02,-0.35104081E-01, 0.53203370E-01,-0.70935440E-04,
     + -0.12980430E-01,-0.10101763E-02,-0.42849420E-02,-0.52871392E-02,
     +  0.57023120E-01,-0.91448444E-02, 0.78549172E-04,-0.11299841E-01,
     + -0.70899750E-03,-0.38641500E-02,-0.49005961E-02, 0.10235171E-05,
     +  0.52533960E-02, 0.18481860E-02, 0.19574640E-02,-0.11436080E-01,
     +  0.19441040E-02, 0.83015663E-02,-0.62191002E-02,-0.53022640E-04,
     + -0.71804760E-02, 0.86780320E-03, 0.81066082E-03, 0.26401460E-02,
     +  0.29410900E-04,-0.88861954E-04,-0.14533742E-02,-0.63020730E-03,
     + -0.77137374E-03, 0.24986700E-02,-0.56262690E-03,-0.24716630E-03,
     + -0.14838993E-02, 0.46197540E-03,-0.35087810E-02, 0.28656134E-02,
     +  0.56370520E-02, 0.38332182E-02,-0.94916723E-03, 0.66941062E-03,
     +  0.38626490E-03, 0.52659750E-03,-0.56894750E-03,-0.24493781E-02,
     + -0.36140393E-02,-0.32432291E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      t_m35_0_7   =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)                *x51
     4  +coeff( 4)*x11                
     5  +coeff( 5)    *x22            
     6  +coeff( 6)        *x32        
     7  +coeff( 7)        *x31*x41    
     8  +coeff( 8)            *x42    
      t_m35_0_7   =t_m35_0_7   
     9  +coeff( 9)    *x21        *x51
     1  +coeff(10)                *x52
     2  +coeff(11)*x11*x21            
     3  +coeff(12)    *x23            
     4  +coeff(13)    *x21*x32        
     5  +coeff(14)    *x21*x31*x41    
     6  +coeff(15)    *x21    *x42    
     7  +coeff(16)*x11            *x51
     8  +coeff(17)    *x22        *x51
      t_m35_0_7   =t_m35_0_7   
     9  +coeff(18)        *x31*x41*x51
     1  +coeff(19)            *x42*x51
     2  +coeff(20)    *x21        *x52
     3  +coeff(21)                *x53
     4  +coeff(22)    *x24            
     5  +coeff(23)    *x22*x31*x41    
     6  +coeff(24)        *x33*x41    
     7  +coeff(25)    *x22    *x42    
     8  +coeff(26)        *x31*x43    
      t_m35_0_7   =t_m35_0_7   
     9  +coeff(27)            *x44    
     1  +coeff(28)    *x23        *x51
     2  +coeff(29)    *x21*x31*x41*x51
     3  +coeff(30)    *x21    *x42*x51
     4  +coeff(31)    *x22        *x52
     5  +coeff(32)        *x31*x41*x52
     6  +coeff(33)            *x42*x52
     7  +coeff(34)    *x21        *x53
     8  +coeff(35)                *x54
      t_m35_0_7   =t_m35_0_7   
     9  +coeff(36)*x11*x23            
     1  +coeff(37)    *x23*x31*x41    
     2  +coeff(38)    *x21*x33*x41    
     3  +coeff(39)    *x23    *x42    
     4  +coeff(40)    *x21*x32*x42    
     5  +coeff(41)    *x21*x31*x43    
     6  +coeff(42)    *x21    *x44    
     7  +coeff(43)    *x24        *x51
     8  +coeff(44)    *x22*x32    *x51
      t_m35_0_7   =t_m35_0_7   
     9  +coeff(45)        *x34    *x51
     1  +coeff(46)    *x22*x31*x41*x51
     2  +coeff(47)    *x23        *x52
     3  +coeff(48)    *x24*x32        
     4  +coeff(49)    *x24*x31*x41    
     5  +coeff(50)    *x24    *x42    
c
      return
      end
      function y_m35_0_7   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.99979E-04,-0.20540E+00,-0.59989E-01,-0.36252E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36252E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.21552420E+00, 0.53170890E+00,-0.41993190E-02, 0.71050850E-02,
     + -0.25364460E-01,-0.39822320E-01,-0.73044110E-02,-0.12148461E-01,
     + -0.21901520E-02,-0.14648640E-02, 0.86113432E-03,-0.46131100E-03,
     +  0.43284874E-02, 0.64716320E-02,-0.19535824E-01,-0.39656393E-01,
     +  0.13208293E-02, 0.15935230E-02,-0.40322930E-02, 0.10533920E-01,
     +  0.17112250E-02, 0.97683570E-02, 0.18364811E-01,-0.10186320E-02,
     +  0.18675570E-02, 0.92830002E-03,-0.12516050E-02, 0.15359490E-02,
     +  0.11435850E-01,-0.83110290E-02, 0.79205930E-03, 0.54041390E-03,
     + -0.35654310E-04,-0.73775544E-03,-0.48759672E-02, 0.42002680E-02,
     + -0.75603760E-03, 0.24496412E-02, 0.15105942E-02, 0.16462820E-02,
     + -0.24167231E-02, 0.43240990E-02, 0.16599972E-02, 0.45597600E-02,
     + -0.55178470E-02,-0.11700020E-02, 0.22623540E-02, 0.40698840E-02,
     +  0.99292760E-02,-0.13653040E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      y_m35_0_7   =avdat
     1  +coeff( 1)        *x31        
     2  +coeff( 2)            *x41    
     3  +coeff( 3)    *x21*x31        
     4  +coeff( 4)    *x21    *x41    
     5  +coeff( 5)        *x31    *x51
     6  +coeff( 6)            *x41*x51
     7  +coeff( 7)    *x22*x31        
     8  +coeff( 8)    *x22    *x41    
      y_m35_0_7   =y_m35_0_7   
     9  +coeff( 9)        *x31*x42    
     1  +coeff(10)            *x43    
     2  +coeff(11)    *x21*x31    *x51
     3  +coeff(12)    *x21    *x41*x51
     4  +coeff(13)        *x31    *x52
     5  +coeff(14)            *x41*x52
     6  +coeff(15)    *x23*x31        
     7  +coeff(16)    *x23    *x41    
     8  +coeff(17)    *x21*x31*x42    
      y_m35_0_7   =y_m35_0_7   
     9  +coeff(18)    *x21    *x43    
     1  +coeff(19)    *x22    *x41*x51
     2  +coeff(20)    *x22*x31*x42    
     3  +coeff(21)    *x22    *x43    
     4  +coeff(22)    *x23*x32*x41    
     5  +coeff(23)    *x23*x31*x42    
     6  +coeff(24)        *x32*x41    
     7  +coeff(25)        *x31*x42*x51
     8  +coeff(26)            *x43*x51
      y_m35_0_7   =y_m35_0_7   
     9  +coeff(27)            *x41*x53
     1  +coeff(28)    *x23*x31    *x51
     2  +coeff(29)    *x23    *x43    
     3  +coeff(30)    *x22*x31*x44    
     4  +coeff(31)    *x21*x33        
     5  +coeff(32)        *x32*x41*x51
     6  +coeff(33)    *x21    *x41*x52
     7  +coeff(34)        *x31    *x53
     8  +coeff(35)    *x24    *x41    
      y_m35_0_7   =y_m35_0_7   
     9  +coeff(36)    *x22*x32*x41    
     1  +coeff(37)        *x32*x43    
     2  +coeff(38)    *x23    *x41*x51
     3  +coeff(39)    *x21    *x43*x51
     4  +coeff(40)    *x22    *x41*x52
     5  +coeff(41)    *x22*x31*x42*x51
     6  +coeff(42)    *x22    *x43*x51
     7  +coeff(43)    *x23    *x41*x52
     8  +coeff(44)    *x24    *x43    
      y_m35_0_7   =y_m35_0_7   
     9  +coeff(45)    *x22*x32*x43    
     1  +coeff(46)        *x33*x44    
     2  +coeff(47)    *x23*x31*x42*x51
     3  +coeff(48)    *x24*x32*x41*x51
     4  +coeff(49)    *x24*x31*x42*x51
     5  +coeff(50)    *x24*x31    *x53
c
      return
      end
      function p_m35_0_7   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.99979E-04,-0.20540E+00,-0.59989E-01,-0.36252E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36252E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.42707770E-02,-0.31849232E-03,-0.63399760E-03,-0.17845261E-02,
     +  0.48446953E-02, 0.12613480E-01,-0.37515630E-04,-0.25705213E-02,
     + -0.13517263E-03,-0.33966083E-04,-0.57179164E-02,-0.13375670E-02,
     + -0.32217580E-02,-0.32611633E-02,-0.80693790E-03,-0.10762213E-02,
     + -0.17004060E-02,-0.35737962E-02,-0.80196290E-03,-0.45074262E-04,
     + -0.19857122E-02, 0.10511451E-02, 0.25420970E-02, 0.19938591E-02,
     + -0.12119970E-04, 0.17004720E-03,-0.16736890E-03, 0.11840040E-02,
     +  0.25621860E-02, 0.18277450E-02, 0.38728570E-03, 0.55716524E-03,
     +  0.38011780E-03, 0.77649630E-03, 0.79432393E-04, 0.20063164E-03,
     + -0.90001000E-03, 0.13585530E-02,-0.61348090E-05, 0.99800934E-04,
     +  0.30582810E-02,-0.73204620E-03, 0.84073650E-04, 0.23015490E-02,
     + -0.15215313E-02,-0.13741280E-02,-0.10753994E-03, 0.13256070E-03,
     + -0.12720460E-03,-0.10651362E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_m35_0_7   =avdat
     1  +coeff( 1)        *x31        
     2  +coeff( 2)            *x41    
     3  +coeff( 3)    *x21*x31        
     4  +coeff( 4)    *x21    *x41    
     5  +coeff( 5)        *x31    *x51
     6  +coeff( 6)            *x41*x51
     7  +coeff( 7)*x11    *x31        
     8  +coeff( 8)    *x22*x31        
      p_m35_0_7   =p_m35_0_7   
     9  +coeff( 9)        *x33        
     1  +coeff(10)*x11        *x41    
     2  +coeff(11)    *x22    *x41    
     3  +coeff(12)        *x32*x41    
     4  +coeff(13)        *x31*x42    
     5  +coeff(14)            *x43    
     6  +coeff(15)    *x21*x31    *x51
     7  +coeff(16)    *x21    *x41*x51
     8  +coeff(17)        *x31    *x52
      p_m35_0_7   =p_m35_0_7   
     9  +coeff(18)            *x41*x52
     1  +coeff(19)    *x23*x31        
     2  +coeff(20)*x11*x21    *x41    
     3  +coeff(21)    *x23    *x41    
     4  +coeff(22)    *x21*x32*x41    
     5  +coeff(23)    *x21*x31*x42    
     6  +coeff(24)    *x21    *x43    
     7  +coeff(25)*x11    *x31    *x51
     8  +coeff(26)    *x22*x31    *x51
      p_m35_0_7   =p_m35_0_7   
     9  +coeff(27)    *x22    *x41*x51
     1  +coeff(28)        *x32*x41*x51
     2  +coeff(29)        *x31*x42*x51
     3  +coeff(30)            *x43*x51
     4  +coeff(31)    *x21*x31    *x52
     5  +coeff(32)    *x21    *x41*x52
     6  +coeff(33)        *x31    *x53
     7  +coeff(34)            *x41*x53
     8  +coeff(35)*x11*x22*x31        
      p_m35_0_7   =p_m35_0_7   
     9  +coeff(36)    *x22*x33        
     1  +coeff(37)    *x24    *x41    
     2  +coeff(38)    *x22*x32*x41    
     3  +coeff(39)        *x34*x41    
     4  +coeff(40)*x11    *x31*x42    
     5  +coeff(41)    *x22*x31*x42    
     6  +coeff(42)        *x33*x42    
     7  +coeff(43)*x11        *x43    
     8  +coeff(44)    *x22    *x43    
      p_m35_0_7   =p_m35_0_7   
     9  +coeff(45)        *x32*x43    
     1  +coeff(46)        *x31*x44    
     2  +coeff(47)    *x23*x31    *x51
     3  +coeff(48)    *x21*x33    *x51
     4  +coeff(49)*x11*x21    *x41*x51
     5  +coeff(50)    *x23    *x41*x51
c
      return
      end
      function l_m35_0_7   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(25)
      data ncoeff/ 24/
      data avdat/ -0.2470803E-01/
      data xmin/
     1 -0.99979E-04,-0.20540E+00,-0.59989E-01,-0.36252E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36252E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.26926650E-01,-0.14294171E+00,-0.61999384E-01,-0.16444830E-01,
     + -0.23206511E-01,-0.17811010E-01,-0.74367133E-04,-0.15575040E-01,
     +  0.11006460E-01,-0.27113100E-02, 0.48222513E-02, 0.51425350E-02,
     + -0.50278200E-02,-0.14494110E-02, 0.40012961E-02,-0.44908402E-02,
     +  0.15141594E-02, 0.14479880E-02, 0.10531274E-02,-0.26915180E-02,
     +  0.34197910E-03, 0.50266430E-03, 0.16695734E-02, 0.16594720E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      l_m35_0_7   =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)    *x22            
     4  +coeff( 4)                *x51
     5  +coeff( 5)            *x42    
     6  +coeff( 6)    *x21        *x51
     7  +coeff( 7)        *x33*x41    
     8  +coeff( 8)        *x31*x41    
      l_m35_0_7   =l_m35_0_7   
     9  +coeff( 9)    *x22        *x51
     1  +coeff(10)        *x32        
     2  +coeff(11)        *x31*x41*x51
     3  +coeff(12)            *x42*x51
     4  +coeff(13)    *x22        *x52
     5  +coeff(14)    *x21        *x54
     6  +coeff(15)    *x23            
     7  +coeff(16)    *x24            
     8  +coeff(17)    *x21*x31*x41    
      l_m35_0_7   =l_m35_0_7   
     9  +coeff(18)    *x21    *x42    
     1  +coeff(19)        *x32    *x51
     2  +coeff(20)    *x23        *x52
     3  +coeff(21)    *x21*x32        
     4  +coeff(22)                *x53
     5  +coeff(23)    *x22*x31*x41    
     6  +coeff(24)    *x22    *x42    
c
      return
      end
      function x_m35_0_8   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/ -0.2741403E-01/
      data xmin/
     1 -0.99979E-04,-0.20540E+00,-0.59989E-01,-0.36252E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36252E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.23430680E-01, 0.26711540E+00, 0.14356391E+00,-0.44661190E-01,
     + -0.80202960E-02,-0.99495090E-02, 0.20345390E+00,-0.19850013E-01,
     + -0.33991582E-01,-0.18681144E-01, 0.12218211E-01,-0.28803190E-01,
     +  0.24807980E-01,-0.16914382E-01,-0.15124703E-01,-0.95202150E-02,
     + -0.17845500E-02,-0.27688810E-02,-0.70575461E-02, 0.29933534E-02,
     + -0.54331421E-02,-0.17183203E-01, 0.39767120E-02, 0.44599250E-02,
     + -0.43010320E-01, 0.30570950E-02, 0.30000593E-02,-0.47344560E-02,
     +  0.51603460E-03, 0.68287090E-02,-0.16667930E-03, 0.62349421E-03,
     +  0.20329730E-02, 0.13893691E-02,-0.16835862E-02, 0.54970980E-02,
     + -0.85876030E-03,-0.12724080E-01, 0.22199450E-01,-0.22861110E-01,
     + -0.63822143E-02, 0.40651012E-01, 0.52503620E-02, 0.10097450E-02,
     +  0.57304860E-03,-0.97635621E-03,-0.12506300E-02,-0.66909042E-03,
     +  0.23179983E-01, 0.11141704E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x25 = x24*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      x_m35_0_8   =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)                *x51
     4  +coeff( 4)    *x22            
     5  +coeff( 5)        *x31*x41    
     6  +coeff( 6)            *x42    
     7  +coeff( 7)    *x21        *x51
     8  +coeff( 8)                *x52
      x_m35_0_8   =x_m35_0_8   
     9  +coeff( 9)    *x23            
     1  +coeff(10)    *x21    *x42    
     2  +coeff(11)    *x22        *x51
     3  +coeff(12)    *x21        *x52
     4  +coeff(13)    *x24            
     5  +coeff(14)    *x22    *x42    
     6  +coeff(15)    *x23*x31*x41    
     7  +coeff(16)    *x24*x31*x41    
     8  +coeff(17)        *x32        
      x_m35_0_8   =x_m35_0_8   
     9  +coeff(18)    *x21*x32        
     1  +coeff(19)    *x21*x31*x41    
     2  +coeff(20)                *x53
     3  +coeff(21)    *x22*x32        
     4  +coeff(22)    *x22*x31*x41    
     5  +coeff(23)    *x23        *x51
     6  +coeff(24)    *x21        *x53
     7  +coeff(25)    *x25    *x42    
     8  +coeff(26)        *x31*x41*x51
      x_m35_0_8   =x_m35_0_8   
     9  +coeff(27)            *x42*x51
     1  +coeff(28)    *x25            
     2  +coeff(29)    *x21*x31*x43    
     3  +coeff(30)    *x21    *x44    
     4  +coeff(31)*x11                
     5  +coeff(32)        *x32    *x51
     6  +coeff(33)        *x31*x43    
     7  +coeff(34)            *x44    
     8  +coeff(35)    *x22        *x52
      x_m35_0_8   =x_m35_0_8   
     9  +coeff(36)    *x21*x32*x42    
     1  +coeff(37)    *x24        *x51
     2  +coeff(38)    *x24    *x42    
     3  +coeff(39)    *x23*x31*x43    
     4  +coeff(40)    *x23    *x44    
     5  +coeff(41)    *x24    *x42*x51
     6  +coeff(42)    *x25    *x44    
     7  +coeff(43)    *x22    *x44*x53
     8  +coeff(44)        *x32*x42    
      x_m35_0_8   =x_m35_0_8   
     9  +coeff(45)    *x21*x31*x41*x51
     1  +coeff(46)        *x31*x41*x52
     2  +coeff(47)            *x42*x52
     3  +coeff(48)                *x54
     4  +coeff(49)    *x23    *x42    
     5  +coeff(50)    *x22*x32    *x51
c
      return
      end
      function t_m35_0_8   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/ -0.1253666E-01/
      data xmin/
     1 -0.99979E-04,-0.20540E+00,-0.59989E-01,-0.36252E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36252E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.12692860E-01,-0.12001230E-01, 0.68517480E-01,-0.91890820E-04,
     + -0.16326540E-01,-0.16878163E-02,-0.61642450E-02,-0.80911554E-02,
     +  0.68801510E-01,-0.12888224E-01, 0.56691350E-04,-0.14167810E-01,
     + -0.80656370E-03,-0.42938120E-02,-0.57375933E-02,-0.61143090E-05,
     +  0.60823110E-02, 0.30171081E-02, 0.27676611E-02,-0.14728790E-01,
     +  0.27857131E-02, 0.95377163E-02,-0.10529072E-01,-0.86639454E-04,
     + -0.11253120E-01, 0.56866350E-03, 0.36854024E-02, 0.60522270E-03,
     +  0.66031020E-03,-0.21582732E-02,-0.97306900E-03,-0.12297820E-02,
     +  0.40045820E-02,-0.89517390E-03,-0.25461983E-03,-0.21144910E-02,
     +  0.14122170E-03,-0.43933084E-02, 0.23455102E-02, 0.58178552E-02,
     +  0.42131943E-02,-0.15030030E-02, 0.85031922E-03, 0.63963112E-03,
     +  0.12038920E-02, 0.78185380E-03,-0.11212400E-02, 0.11932690E-02,
     + -0.13333570E-02,-0.27341763E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      t_m35_0_8   =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)                *x51
     4  +coeff( 4)*x11                
     5  +coeff( 5)    *x22            
     6  +coeff( 6)        *x32        
     7  +coeff( 7)        *x31*x41    
     8  +coeff( 8)            *x42    
      t_m35_0_8   =t_m35_0_8   
     9  +coeff( 9)    *x21        *x51
     1  +coeff(10)                *x52
     2  +coeff(11)*x11*x21            
     3  +coeff(12)    *x23            
     4  +coeff(13)    *x21*x32        
     5  +coeff(14)    *x21*x31*x41    
     6  +coeff(15)    *x21    *x42    
     7  +coeff(16)*x11            *x51
     8  +coeff(17)    *x22        *x51
      t_m35_0_8   =t_m35_0_8   
     9  +coeff(18)        *x31*x41*x51
     1  +coeff(19)            *x42*x51
     2  +coeff(20)    *x21        *x52
     3  +coeff(21)                *x53
     4  +coeff(22)    *x24            
     5  +coeff(23)    *x22*x31*x41    
     6  +coeff(24)        *x33*x41    
     7  +coeff(25)    *x22    *x42    
     8  +coeff(26)            *x44    
      t_m35_0_8   =t_m35_0_8   
     9  +coeff(27)    *x23        *x51
     1  +coeff(28)    *x21*x31*x41*x51
     2  +coeff(29)    *x21    *x42*x51
     3  +coeff(30)    *x22        *x52
     4  +coeff(31)        *x31*x41*x52
     5  +coeff(32)            *x42*x52
     6  +coeff(33)    *x21        *x53
     7  +coeff(34)                *x54
     8  +coeff(35)*x11*x23            
      t_m35_0_8   =t_m35_0_8   
     9  +coeff(36)    *x23*x31*x41    
     1  +coeff(37)    *x21*x33*x41    
     2  +coeff(38)    *x23    *x42    
     3  +coeff(39)    *x21*x32*x42    
     4  +coeff(40)    *x21*x31*x43    
     5  +coeff(41)    *x21    *x44    
     6  +coeff(42)    *x24        *x51
     7  +coeff(43)    *x22*x32    *x51
     8  +coeff(44)        *x34    *x51
      t_m35_0_8   =t_m35_0_8   
     9  +coeff(45)    *x22*x31*x41*x51
     1  +coeff(46)            *x44*x51
     2  +coeff(47)    *x23        *x52
     3  +coeff(48)    *x22        *x53
     4  +coeff(49)    *x21        *x54
     5  +coeff(50)    *x24*x32        
c
      return
      end
      function y_m35_0_8   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.99979E-04,-0.20540E+00,-0.59989E-01,-0.36252E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36252E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.20886793E+00, 0.51994080E+00,-0.44678290E-02, 0.59610740E-02,
     + -0.21745374E-01,-0.31022112E-01,-0.79306170E-02,-0.12277813E-01,
     + -0.46468190E-02,-0.27929281E-02, 0.31982120E-04,-0.87805690E-03,
     +  0.29125550E-02, 0.41429540E-02,-0.19282110E-01,-0.39118964E-01,
     +  0.31492140E-02, 0.26967330E-02,-0.29728484E-02,-0.11760332E-02,
     +  0.14869390E-01, 0.75966250E-03, 0.84567650E-02, 0.17141684E-01,
     +  0.30915292E-02, 0.20783030E-02,-0.73901550E-02, 0.22249461E-02,
     +  0.11165080E-01, 0.17456813E-02, 0.26833564E-02,-0.10888590E-01,
     +  0.27624080E-02,-0.26288460E-03, 0.91898320E-03, 0.13356820E-02,
     +  0.83265034E-03,-0.49796560E-03,-0.57167460E-03,-0.11054133E-02,
     + -0.13179960E-02, 0.16814321E-02, 0.18061681E-02, 0.20094981E-02,
     +  0.11523762E-02, 0.22421721E-02, 0.31511073E-02,-0.46522810E-02,
     +  0.62662820E-02,-0.52522951E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      y_m35_0_8   =avdat
     1  +coeff( 1)        *x31        
     2  +coeff( 2)            *x41    
     3  +coeff( 3)    *x21*x31        
     4  +coeff( 4)    *x21    *x41    
     5  +coeff( 5)        *x31    *x51
     6  +coeff( 6)            *x41*x51
     7  +coeff( 7)    *x22*x31        
     8  +coeff( 8)    *x22    *x41    
      y_m35_0_8   =y_m35_0_8   
     9  +coeff( 9)        *x31*x42    
     1  +coeff(10)            *x43    
     2  +coeff(11)    *x21*x31    *x51
     3  +coeff(12)    *x21    *x41*x51
     4  +coeff(13)        *x31    *x52
     5  +coeff(14)            *x41*x52
     6  +coeff(15)    *x23*x31        
     7  +coeff(16)    *x23    *x41    
     8  +coeff(17)    *x21*x31*x42    
      y_m35_0_8   =y_m35_0_8   
     9  +coeff(18)    *x21    *x43    
     1  +coeff(19)    *x22    *x41*x51
     2  +coeff(20)        *x34*x41    
     3  +coeff(21)    *x22*x31*x42    
     4  +coeff(22)    *x22    *x43    
     5  +coeff(23)    *x23*x32*x41    
     6  +coeff(24)    *x23*x31*x42    
     7  +coeff(25)        *x31*x42*x51
     8  +coeff(26)            *x43*x51
      y_m35_0_8   =y_m35_0_8   
     9  +coeff(27)    *x24    *x41    
     1  +coeff(28)    *x23*x31    *x51
     2  +coeff(29)    *x23    *x43    
     3  +coeff(30)        *x34*x41*x51
     4  +coeff(31)    *x22*x34*x41    
     5  +coeff(32)    *x22*x31*x44    
     6  +coeff(33)    *x24    *x41*x52
     7  +coeff(34)        *x32*x41    
     8  +coeff(35)    *x21*x33        
      y_m35_0_8   =y_m35_0_8   
     9  +coeff(36)    *x21*x32*x41    
     1  +coeff(37)    *x21    *x41*x52
     2  +coeff(38)        *x31    *x53
     3  +coeff(39)            *x41*x53
     4  +coeff(40)    *x24*x31        
     5  +coeff(41)        *x32*x43    
     6  +coeff(42)    *x23    *x41*x51
     7  +coeff(43)    *x21*x31*x42*x51
     8  +coeff(44)    *x21    *x43*x51
      y_m35_0_8   =y_m35_0_8   
     9  +coeff(45)    *x22*x31    *x52
     1  +coeff(46)    *x22*x31*x42*x51
     2  +coeff(47)    *x22    *x43*x51
     3  +coeff(48)    *x22*x33*x42    
     4  +coeff(49)    *x24    *x43    
     5  +coeff(50)    *x22*x32*x43    
c
      return
      end
      function p_m35_0_8   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.99979E-04,-0.20540E+00,-0.59989E-01,-0.36252E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36252E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.22145980E-01,-0.44903144E-01, 0.27810350E-03,-0.37901730E-03,
     +  0.98290010E-02, 0.23369380E-01,-0.19577720E-02,-0.33960462E-03,
     + -0.13704030E-04,-0.44795100E-02,-0.13277394E-02,-0.35363930E-02,
     + -0.33036174E-02,-0.44583250E-03,-0.29053850E-02,-0.58607770E-02,
     +  0.61286822E-03, 0.47133540E-03, 0.30776034E-02, 0.11127284E-02,
     +  0.19689110E-02, 0.13396691E-02,-0.11681712E-03, 0.87478774E-03,
     +  0.16188661E-02,-0.62094430E-03, 0.11572852E-03, 0.19880414E-02,
     +  0.16437980E-02,-0.13550273E-02,-0.12539520E-02,-0.88121750E-03,
     + -0.24607460E-03,-0.20479390E-03,-0.14004460E-02,-0.30918931E-03,
     +  0.68511580E-03,-0.57506054E-03, 0.69981184E-03,-0.37992330E-03,
     +  0.49581090E-03,-0.80470710E-03,-0.15162080E-02,-0.26390381E-03,
     +  0.62712802E-04, 0.58244290E-03, 0.11068022E-02, 0.95670380E-03,
     +  0.85580230E-03, 0.13187492E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      p_m35_0_8   =avdat
     1  +coeff( 1)        *x31        
     2  +coeff( 2)            *x41    
     3  +coeff( 3)    *x21*x31        
     4  +coeff( 4)    *x21    *x41    
     5  +coeff( 5)        *x31    *x51
     6  +coeff( 6)            *x41*x51
     7  +coeff( 7)    *x22*x31        
     8  +coeff( 8)        *x33        
      p_m35_0_8   =p_m35_0_8   
     9  +coeff( 9)*x11        *x41    
     1  +coeff(10)    *x22    *x41    
     2  +coeff(11)        *x32*x41    
     3  +coeff(12)        *x31*x42    
     4  +coeff(13)            *x43    
     5  +coeff(14)    *x21    *x41*x51
     6  +coeff(15)        *x31    *x52
     7  +coeff(16)            *x41*x52
     8  +coeff(17)    *x23*x31        
      p_m35_0_8   =p_m35_0_8   
     9  +coeff(18)    *x23    *x41    
     1  +coeff(19)    *x21*x31*x42    
     2  +coeff(20)    *x21    *x43    
     3  +coeff(21)        *x31*x42*x51
     4  +coeff(22)            *x43*x51
     5  +coeff(23)    *x21    *x41*x52
     6  +coeff(24)        *x31    *x53
     7  +coeff(25)            *x41*x53
     8  +coeff(26)    *x24    *x41    
      p_m35_0_8   =p_m35_0_8   
     9  +coeff(27)        *x34*x41    
     1  +coeff(28)    *x22*x31*x42    
     2  +coeff(29)    *x22    *x43    
     3  +coeff(30)        *x32*x43    
     4  +coeff(31)        *x31*x44    
     5  +coeff(32)    *x23*x31    *x51
     6  +coeff(33)    *x21*x33    *x51
     7  +coeff(34)*x11*x21    *x41*x51
     8  +coeff(35)    *x23    *x41*x51
      p_m35_0_8   =p_m35_0_8   
     9  +coeff(36)    *x21*x32*x41*x51
     1  +coeff(37)    *x21    *x43*x51
     2  +coeff(38)        *x31*x42*x52
     3  +coeff(39)    *x21    *x41*x53
     4  +coeff(40)            *x41*x54
     5  +coeff(41)    *x21*x34*x41    
     6  +coeff(42)    *x23*x31*x42    
     7  +coeff(43)    *x21*x31*x44    
     8  +coeff(44)*x11*x22    *x41*x51
      p_m35_0_8   =p_m35_0_8   
     9  +coeff(45)    *x22*x32*x41*x51
     1  +coeff(46)        *x34*x41*x51
     2  +coeff(47)    *x22*x31*x42*x51
     3  +coeff(48)        *x33*x42*x51
     4  +coeff(49)    *x22    *x43*x51
     5  +coeff(50)        *x32*x43*x51
c
      return
      end
      function l_m35_0_8   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(25)
      data ncoeff/ 24/
      data avdat/ -0.2486142E-01/
      data xmin/
     1 -0.99979E-04,-0.20540E+00,-0.59989E-01,-0.36252E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36252E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.26858860E-01,-0.14178721E+00,-0.61113201E-01,-0.17918832E-01,
     + -0.22627500E-01,-0.19722264E-01,-0.86244240E-04,-0.15410703E-01,
     +  0.12776924E-01,-0.47718930E-02,-0.69110090E-02,-0.28036820E-02,
     +  0.52592963E-02, 0.56765930E-02, 0.34556760E-02,-0.52247940E-02,
     + -0.35977642E-03, 0.12720110E-02, 0.13142150E-02, 0.11461120E-02,
     +  0.10682002E-02, 0.22944472E-02,-0.12614160E-02, 0.20234890E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      l_m35_0_8   =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)    *x22            
     4  +coeff( 4)                *x51
     5  +coeff( 5)            *x42    
     6  +coeff( 6)    *x21        *x51
     7  +coeff( 7)        *x33*x41    
     8  +coeff( 8)        *x31*x41    
      l_m35_0_8   =l_m35_0_8   
     9  +coeff( 9)    *x22        *x51
     1  +coeff(10)    *x21        *x52
     2  +coeff(11)    *x22        *x52
     3  +coeff(12)        *x32        
     4  +coeff(13)        *x31*x41*x51
     5  +coeff(14)            *x42*x51
     6  +coeff(15)    *x23            
     7  +coeff(16)    *x24            
     8  +coeff(17)                *x52
      l_m35_0_8   =l_m35_0_8   
     9  +coeff(18)    *x21*x31*x41    
     1  +coeff(19)    *x21    *x42    
     2  +coeff(20)        *x32    *x51
     3  +coeff(21)                *x53
     4  +coeff(22)    *x23        *x51
     5  +coeff(23)            *x42*x52
     6  +coeff(24)    *x21        *x53
c
      return
      end
      function x_m35_0_9   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/ -0.3512337E-01/
      data xmin/
     1 -0.99979E-04,-0.20540E+00,-0.59989E-01,-0.36252E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36252E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.31113151E-01, 0.26680591E+00, 0.18163812E+00,-0.53950560E-01,
     + -0.11330091E-01,-0.14073071E-01, 0.24321170E+00,-0.28400810E-01,
     + -0.41015412E-01,-0.23245320E-01, 0.16697050E-01,-0.38240663E-01,
     +  0.30209954E-01,-0.22360913E-01,-0.22190770E-01,-0.20941440E-01,
     + -0.22921164E-02,-0.66051413E-02, 0.48176134E-02, 0.29940570E-02,
     + -0.40985392E-02, 0.49166833E-02, 0.50535820E-02,-0.71227103E-02,
     + -0.55169750E-02, 0.22913331E-01, 0.13706110E-02,-0.29456482E-02,
     +  0.10093690E-01,-0.33653840E-02, 0.32833450E-01,-0.23223450E-01,
     +  0.40042940E-02,-0.20801572E-03, 0.16934581E-02, 0.15995813E-02,
     + -0.19238920E-02,-0.21740023E-02, 0.54707061E-02, 0.13462483E-01,
     + -0.12537112E-02,-0.11354740E-01,-0.14918280E-01, 0.38895623E-02,
     + -0.49368701E-01,-0.11032100E-01,-0.16575813E-01, 0.48058730E-01,
     +  0.11157574E-01, 0.71634720E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x25 = x24*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      x_m35_0_9   =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)                *x51
     4  +coeff( 4)    *x22            
     5  +coeff( 5)        *x31*x41    
     6  +coeff( 6)            *x42    
     7  +coeff( 7)    *x21        *x51
     8  +coeff( 8)                *x52
      x_m35_0_9   =x_m35_0_9   
     9  +coeff( 9)    *x23            
     1  +coeff(10)    *x21    *x42    
     2  +coeff(11)    *x22        *x51
     3  +coeff(12)    *x21        *x52
     4  +coeff(13)    *x24            
     5  +coeff(14)    *x22*x31*x41    
     6  +coeff(15)    *x22    *x42    
     7  +coeff(16)    *x23*x31*x41    
     8  +coeff(17)        *x32        
      x_m35_0_9   =x_m35_0_9   
     9  +coeff(18)    *x21*x31*x41    
     1  +coeff(19)                *x53
     2  +coeff(20)    *x23        *x53
     3  +coeff(21)    *x21*x32        
     4  +coeff(22)        *x31*x41*x51
     5  +coeff(23)            *x42*x51
     6  +coeff(24)    *x22*x32        
     7  +coeff(25)    *x25            
     8  +coeff(26)    *x23    *x42    
      x_m35_0_9   =x_m35_0_9   
     9  +coeff(27)        *x32    *x51
     1  +coeff(28)    *x21*x31*x43    
     2  +coeff(29)    *x21    *x44    
     3  +coeff(30)    *x22        *x54
     4  +coeff(31)    *x23*x31*x43    
     5  +coeff(32)    *x23    *x44    
     6  +coeff(33)    *x23*x34*x42    
     7  +coeff(34)*x11                
     8  +coeff(35)        *x31*x43    
      x_m35_0_9   =x_m35_0_9   
     9  +coeff(36)            *x44    
     1  +coeff(37)        *x31*x41*x52
     2  +coeff(38)            *x42*x52
     3  +coeff(39)    *x21        *x53
     4  +coeff(40)    *x21*x32*x42    
     5  +coeff(41)    *x24        *x51
     6  +coeff(42)    *x24*x31*x41    
     7  +coeff(43)    *x24    *x42    
     8  +coeff(44)    *x25        *x51
      x_m35_0_9   =x_m35_0_9   
     9  +coeff(45)    *x25    *x42    
     1  +coeff(46)    *x21*x32*x44    
     2  +coeff(47)    *x24    *x42*x51
     3  +coeff(48)    *x25    *x44    
     4  +coeff(49)    *x24    *x44*x51
     5  +coeff(50)    *x24    *x42*x53
c
      return
      end
      function t_m35_0_9   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/ -0.1642303E-01/
      data xmin/
     1 -0.99979E-04,-0.20540E+00,-0.59989E-01,-0.36252E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36252E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.16341503E-01, 0.10837570E-01, 0.84667760E-01,-0.91979964E-04,
     + -0.20499700E-01,-0.24835920E-02,-0.80312650E-02,-0.98465440E-02,
     +  0.82764800E-01,-0.17114480E-01, 0.18528371E-03,-0.18195260E-01,
     + -0.12128000E-02,-0.53369310E-02,-0.76925350E-02, 0.13495890E-04,
     +  0.82635660E-02, 0.38878650E-02, 0.44661043E-02,-0.22476714E-01,
     +  0.48902831E-02, 0.27840200E-04,-0.51671200E-04, 0.11044540E-01,
     +  0.17080082E-03,-0.12335511E-01, 0.10115540E-03,-0.14163640E-01,
     +  0.16655900E-02, 0.15996860E-02,-0.12274271E-03, 0.55930000E-02,
     +  0.10124050E-02, 0.11166610E-02,-0.34355720E-02,-0.18592784E-02,
     + -0.24651280E-02, 0.69612553E-02,-0.14761710E-02,-0.52851420E-03,
     + -0.31873660E-02,-0.64295360E-05,-0.55537670E-02, 0.31356604E-02,
     +  0.76085664E-02, 0.57606190E-02,-0.17431090E-02, 0.27009651E-02,
     +  0.68445754E-03, 0.22679620E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      t_m35_0_9   =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)                *x51
     4  +coeff( 4)*x11                
     5  +coeff( 5)    *x22            
     6  +coeff( 6)        *x32        
     7  +coeff( 7)        *x31*x41    
     8  +coeff( 8)            *x42    
      t_m35_0_9   =t_m35_0_9   
     9  +coeff( 9)    *x21        *x51
     1  +coeff(10)                *x52
     2  +coeff(11)*x11*x21            
     3  +coeff(12)    *x23            
     4  +coeff(13)    *x21*x32        
     5  +coeff(14)    *x21*x31*x41    
     6  +coeff(15)    *x21    *x42    
     7  +coeff(16)*x11            *x51
     8  +coeff(17)    *x22        *x51
      t_m35_0_9   =t_m35_0_9   
     9  +coeff(18)        *x31*x41*x51
     1  +coeff(19)            *x42*x51
     2  +coeff(20)    *x21        *x52
     3  +coeff(21)                *x53
     4  +coeff(22)*x12                
     5  +coeff(23)*x11*x22            
     6  +coeff(24)    *x24            
     7  +coeff(25)        *x34        
     8  +coeff(26)    *x22*x31*x41    
      t_m35_0_9   =t_m35_0_9   
     9  +coeff(27)        *x33*x41    
     1  +coeff(28)    *x22    *x42    
     2  +coeff(29)        *x31*x43    
     3  +coeff(30)            *x44    
     4  +coeff(31)*x11*x21        *x51
     5  +coeff(32)    *x23        *x51
     6  +coeff(33)    *x21*x31*x41*x51
     7  +coeff(34)    *x21    *x42*x51
     8  +coeff(35)    *x22        *x52
      t_m35_0_9   =t_m35_0_9   
     9  +coeff(36)        *x31*x41*x52
     1  +coeff(37)            *x42*x52
     2  +coeff(38)    *x21        *x53
     3  +coeff(39)                *x54
     4  +coeff(40)*x11*x23            
     5  +coeff(41)    *x23*x31*x41    
     6  +coeff(42)    *x21*x33*x41    
     7  +coeff(43)    *x23    *x42    
     8  +coeff(44)    *x21*x32*x42    
      t_m35_0_9   =t_m35_0_9   
     9  +coeff(45)    *x21*x31*x43    
     1  +coeff(46)    *x21    *x44    
     2  +coeff(47)    *x24        *x51
     3  +coeff(48)    *x22*x32    *x51
     4  +coeff(49)        *x34    *x51
     5  +coeff(50)    *x22*x31*x41*x51
c
      return
      end
      function y_m35_0_9   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.99979E-04,-0.20540E+00,-0.59989E-01,-0.36252E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36252E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.19366754E+00, 0.48726970E+00,-0.40995553E-02, 0.59216083E-02,
     + -0.15604973E-01,-0.16927540E-01,-0.98049800E-02,-0.16285400E-01,
     + -0.71340752E-02,-0.49699590E-02,-0.40029230E-03,-0.59722930E-03,
     +  0.10171330E-02,-0.18670170E-01,-0.38276970E-01, 0.41555361E-02,
     +  0.35797260E-02,-0.30698070E-02, 0.48898120E-02, 0.31554230E-02,
     +  0.56634680E-02,-0.30657490E-03, 0.17790094E-01, 0.38875692E-02,
     +  0.83268871E-02, 0.16466584E-01,-0.26172590E-02, 0.42201334E-03,
     +  0.17552820E-02, 0.22902920E-02,-0.65639550E-02, 0.23676050E-02,
     +  0.32304213E-02, 0.31816333E-02, 0.10867750E-01,-0.12064794E-01,
     +  0.23608903E-02, 0.10290670E-01,-0.57607120E-02,-0.13216730E-03,
     +  0.93860324E-03,-0.57522110E-03,-0.35810410E-03, 0.13152942E-02,
     + -0.19151330E-02, 0.39414321E-02,-0.52056922E-02, 0.49733040E-02,
     + -0.93816211E-02, 0.28023002E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      y_m35_0_9   =avdat
     1  +coeff( 1)        *x31        
     2  +coeff( 2)            *x41    
     3  +coeff( 3)    *x21*x31        
     4  +coeff( 4)    *x21    *x41    
     5  +coeff( 5)        *x31    *x51
     6  +coeff( 6)            *x41*x51
     7  +coeff( 7)    *x22*x31        
     8  +coeff( 8)    *x22    *x41    
      y_m35_0_9   =y_m35_0_9   
     9  +coeff( 9)        *x31*x42    
     1  +coeff(10)            *x43    
     2  +coeff(11)    *x21*x31    *x51
     3  +coeff(12)    *x21    *x41*x51
     4  +coeff(13)        *x31    *x52
     5  +coeff(14)    *x23*x31        
     6  +coeff(15)    *x23    *x41    
     7  +coeff(16)    *x21*x31*x42    
     8  +coeff(17)    *x21    *x43    
      y_m35_0_9   =y_m35_0_9   
     9  +coeff(18)    *x22    *x41*x51
     1  +coeff(19)        *x31*x42*x51
     2  +coeff(20)            *x43*x51
     3  +coeff(21)    *x22*x32*x41    
     4  +coeff(22)        *x34*x41    
     5  +coeff(23)    *x22*x31*x42    
     6  +coeff(24)    *x22    *x43    
     7  +coeff(25)    *x23*x32*x41    
     8  +coeff(26)    *x23*x31*x42    
      y_m35_0_9   =y_m35_0_9   
     9  +coeff(27)        *x32*x41    
     1  +coeff(28)            *x41*x52
     2  +coeff(29)    *x21*x32*x41    
     3  +coeff(30)        *x32*x41*x51
     4  +coeff(31)    *x24    *x41    
     5  +coeff(32)    *x23*x31    *x51
     6  +coeff(33)    *x21*x31*x42*x51
     7  +coeff(34)    *x21    *x43*x51
     8  +coeff(35)    *x23    *x43    
      y_m35_0_9   =y_m35_0_9   
     9  +coeff(36)    *x22*x31*x44    
     1  +coeff(37)    *x24    *x41*x52
     2  +coeff(38)    *x24*x31*x42*x51
     3  +coeff(39)    *x24*x31*x42*x53
     4  +coeff(40)        *x33        
     5  +coeff(41)    *x21*x33        
     6  +coeff(42)        *x32*x43    
     7  +coeff(43)    *x21*x32*x41*x51
     8  +coeff(44)    *x22*x31    *x52
      y_m35_0_9   =y_m35_0_9   
     9  +coeff(45)    *x22*x31*x42*x51
     1  +coeff(46)    *x22    *x43*x51
     2  +coeff(47)    *x22*x33*x42    
     3  +coeff(48)    *x24    *x43    
     4  +coeff(49)    *x22*x32*x43    
     5  +coeff(50)    *x23*x32*x41*x51
c
      return
      end
      function p_m35_0_9   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.99979E-04,-0.20540E+00,-0.59989E-01,-0.36252E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36252E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.38917070E-01,-0.86292680E-01, 0.81804784E-03,-0.19205510E-03,
     +  0.14177630E-01, 0.32588541E-01,-0.14806950E-02,-0.18529663E-04,
     + -0.46299510E-02,-0.23748120E-02,-0.47699413E-02,-0.42676710E-02,
     + -0.79715833E-03,-0.75429904E-03,-0.44610170E-02,-0.82633644E-02,
     +  0.21042940E-02, 0.31490220E-02, 0.12171980E-02, 0.43049210E-02,
     +  0.12975290E-02, 0.42027090E-02, 0.26563010E-02,-0.55189660E-03,
     + -0.70889951E-03, 0.13752390E-02, 0.25183460E-02,-0.44179413E-03,
     +  0.81568720E-03, 0.28260440E-02, 0.22187240E-02,-0.80309453E-03,
     + -0.11621392E-02, 0.34939980E-02, 0.33255080E-02,-0.71770090E-03,
     +  0.15387440E-02,-0.79380320E-03,-0.18190644E-02,-0.20574440E-02,
     + -0.22574840E-02, 0.13696720E-03, 0.15187443E-02, 0.19924920E-02,
     +  0.10862353E-02, 0.12166091E-02, 0.71200804E-03,-0.13977670E-02,
     + -0.15443110E-02, 0.98264280E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      p_m35_0_9   =avdat
     1  +coeff( 1)        *x31        
     2  +coeff( 2)            *x41    
     3  +coeff( 3)    *x21*x31        
     4  +coeff( 4)    *x21    *x41    
     5  +coeff( 5)        *x31    *x51
     6  +coeff( 6)            *x41*x51
     7  +coeff( 7)    *x22*x31        
     8  +coeff( 8)*x11        *x41    
      p_m35_0_9   =p_m35_0_9   
     9  +coeff( 9)    *x22    *x41    
     1  +coeff(10)        *x32*x41    
     2  +coeff(11)        *x31*x42    
     3  +coeff(12)            *x43    
     4  +coeff(13)    *x21*x31    *x51
     5  +coeff(14)    *x21    *x41*x51
     6  +coeff(15)        *x31    *x52
     7  +coeff(16)            *x41*x52
     8  +coeff(17)    *x23*x31        
      p_m35_0_9   =p_m35_0_9   
     9  +coeff(18)    *x23    *x41    
     1  +coeff(19)    *x21*x32*x41    
     2  +coeff(20)    *x21*x31*x42    
     3  +coeff(21)    *x21    *x43    
     4  +coeff(22)        *x31*x42*x51
     5  +coeff(23)            *x43*x51
     6  +coeff(24)    *x21*x31    *x52
     7  +coeff(25)    *x21    *x41*x52
     8  +coeff(26)        *x31    *x53
      p_m35_0_9   =p_m35_0_9   
     9  +coeff(27)            *x41*x53
     1  +coeff(28)    *x22*x33        
     2  +coeff(29)    *x22*x32*x41    
     3  +coeff(30)    *x22*x31*x42    
     4  +coeff(31)    *x22    *x43    
     5  +coeff(32)        *x31*x44    
     6  +coeff(33)    *x23    *x41*x51
     7  +coeff(34)    *x21*x31*x42*x51
     8  +coeff(35)    *x21    *x43*x51
      p_m35_0_9   =p_m35_0_9   
     9  +coeff(36)    *x22    *x41*x52
     1  +coeff(37)    *x21    *x41*x53
     2  +coeff(38)            *x41*x54
     3  +coeff(39)    *x23*x32*x41    
     4  +coeff(40)    *x23*x31*x42    
     5  +coeff(41)    *x21*x31*x44    
     6  +coeff(42)    *x22*x32*x41*x51
     7  +coeff(43)        *x34*x41*x51
     8  +coeff(44)    *x22*x31*x42*x51
      p_m35_0_9   =p_m35_0_9   
     9  +coeff(45)        *x33*x42*x51
     1  +coeff(46)    *x22    *x43*x51
     2  +coeff(47)        *x32*x43*x51
     3  +coeff(48)        *x31*x44*x51
     4  +coeff(49)    *x21    *x41*x54
     5  +coeff(50)    *x21*x34*x41*x51
c
      return
      end
      function l_m35_0_9   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(25)
      data ncoeff/ 24/
      data avdat/ -0.2586983E-01/
      data xmin/
     1 -0.99979E-04,-0.20540E+00,-0.59989E-01,-0.36252E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36252E-01, 0.14989E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.28190880E-01,-0.14116263E+00,-0.61636500E-01,-0.18753992E-01,
     + -0.24513600E-01,-0.18858522E-01,-0.26445934E-05,-0.16249060E-01,
     + -0.27185690E-02, 0.14437670E-01,-0.31249872E-02,-0.84128260E-02,
     + -0.30059780E-02, 0.67230470E-02, 0.71909860E-02, 0.29641870E-02,
     +  0.14114272E-02, 0.19929580E-02, 0.94124470E-03, 0.51131390E-02,
     + -0.43944320E-02, 0.17458660E-02,-0.73202303E-02,-0.18165570E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      l_m35_0_9   =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)    *x22            
     4  +coeff( 4)                *x51
     5  +coeff( 5)            *x42    
     6  +coeff( 6)    *x21        *x51
     7  +coeff( 7)        *x33*x41    
     8  +coeff( 8)        *x31*x41    
      l_m35_0_9   =l_m35_0_9   
     9  +coeff( 9)                *x52
     1  +coeff(10)    *x22        *x51
     2  +coeff(11)        *x32        
     3  +coeff(12)    *x21        *x52
     4  +coeff(13)    *x24        *x52
     5  +coeff(14)        *x31*x41*x51
     6  +coeff(15)            *x42*x51
     7  +coeff(16)    *x23            
     8  +coeff(17)        *x32    *x51
      l_m35_0_9   =l_m35_0_9   
     9  +coeff(18)                *x53
     1  +coeff(19)    *x21        *x53
     2  +coeff(20)    *x23        *x53
     3  +coeff(21)    *x24            
     4  +coeff(22)    *x21    *x42*x51
     5  +coeff(23)    *x22        *x52
     6  +coeff(24)        *x31*x41*x52
c
      return
      end
      function x_m35_0_10  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/ -0.4732611E-01/
      data xmin/
     1 -0.99979E-04,-0.20540E+00,-0.59989E-01,-0.36252E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36252E-01, 0.14986E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.42973354E-01, 0.27593740E+00, 0.22829510E+00,-0.63522660E-01,
     + -0.15889660E-01,-0.18349554E-01, 0.28719040E+00,-0.37274860E-01,
     + -0.47981900E-01,-0.22217791E-01, 0.20270320E-01,-0.48645270E-01,
     +  0.35333490E-01,-0.28313403E-01,-0.37471370E-01,-0.58069900E-02,
     + -0.33842143E-02,-0.16127470E-01, 0.71746581E-02,-0.90261020E-02,
     +  0.11320720E-01, 0.17880181E-02,-0.37664144E-02, 0.72720423E-02,
     +  0.70792160E-02, 0.89925560E-02,-0.40441830E-01, 0.19194071E-02,
     + -0.40781010E-02,-0.33619361E-02,-0.37567054E-02, 0.15361314E-01,
     +  0.57211504E-02, 0.31521510E-02, 0.23355122E-01,-0.25101600E-01,
     +  0.38278251E-01,-0.23941460E-03, 0.39832284E-02, 0.16368370E-02,
     +  0.74410820E-02,-0.18458960E-02,-0.79305041E-02, 0.82530984E-02,
     + -0.43119480E-02,-0.27900533E-02,-0.12368840E-01,-0.10360650E-01,
     +  0.19297700E-02, 0.59925830E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x25 = x24*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      x_m35_0_10  =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)                *x51
     4  +coeff( 4)    *x22            
     5  +coeff( 5)        *x31*x41    
     6  +coeff( 6)            *x42    
     7  +coeff( 7)    *x21        *x51
     8  +coeff( 8)                *x52
      x_m35_0_10  =x_m35_0_10  
     9  +coeff( 9)    *x23            
     1  +coeff(10)    *x21    *x42    
     2  +coeff(11)    *x22        *x51
     3  +coeff(12)    *x21        *x52
     4  +coeff(13)    *x24            
     5  +coeff(14)    *x22*x31*x41    
     6  +coeff(15)    *x22    *x42    
     7  +coeff(16)    *x23*x31*x41    
     8  +coeff(17)        *x32        
      x_m35_0_10  =x_m35_0_10  
     9  +coeff(18)    *x21*x31*x41    
     1  +coeff(19)                *x53
     2  +coeff(20)    *x22*x32        
     3  +coeff(21)    *x21        *x53
     4  +coeff(22)    *x25        *x51
     5  +coeff(23)    *x21*x32        
     6  +coeff(24)        *x31*x41*x51
     7  +coeff(25)            *x42*x51
     8  +coeff(26)    *x23    *x42    
      x_m35_0_10  =x_m35_0_10  
     9  +coeff(27)    *x25    *x42    
     1  +coeff(28)        *x32    *x51
     2  +coeff(29)    *x22        *x52
     3  +coeff(30)        *x31*x41*x52
     4  +coeff(31)            *x42*x52
     5  +coeff(32)    *x21*x31*x43    
     6  +coeff(33)    *x21    *x44    
     7  +coeff(34)    *x22        *x53
     8  +coeff(35)    *x25    *x44    
      x_m35_0_10  =x_m35_0_10  
     9  +coeff(36)    *x25*x31*x41*x52
     1  +coeff(37)    *x25*x31*x43*x52
     2  +coeff(38)*x11                
     3  +coeff(39)        *x31*x43    
     4  +coeff(40)            *x44    
     5  +coeff(41)    *x23        *x51
     6  +coeff(42)                *x54
     7  +coeff(43)    *x25            
     8  +coeff(44)    *x21*x32*x42    
      x_m35_0_10  =x_m35_0_10  
     9  +coeff(45)    *x24        *x51
     1  +coeff(46)    *x21        *x54
     2  +coeff(47)    *x24*x31*x41    
     3  +coeff(48)    *x24    *x42    
     4  +coeff(49)        *x34*x42    
     5  +coeff(50)    *x22    *x44    
c
      return
      end
      function t_m35_0_10  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/ -0.2032666E-01/
      data xmin/
     1 -0.99979E-04,-0.20540E+00,-0.59989E-01,-0.36252E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36252E-01, 0.14986E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.19914900E-01, 0.29430424E-01, 0.99435791E-01,-0.85993771E-04,
     + -0.23650584E-01,-0.25935983E-02,-0.79332600E-02,-0.98471320E-02,
     +  0.97896374E-01,-0.20830821E-01, 0.18565430E-03,-0.21579233E-01,
     + -0.15702990E-02,-0.73512024E-02,-0.10508941E-01, 0.45891032E-04,
     +  0.11622152E-01, 0.83618913E-03, 0.53140511E-02, 0.57135700E-02,
     + -0.27345202E-01, 0.69222990E-02, 0.55376243E-04, 0.12990484E-01,
     +  0.65222244E-04,-0.16418810E-01, 0.30518100E-04,-0.10252623E-03,
     + -0.18923142E-01, 0.27262372E-02, 0.24574224E-02,-0.32168121E-04,
     +  0.76292030E-02, 0.82887813E-03, 0.71368104E-03,-0.32045390E-02,
     + -0.35925370E-02,-0.45105300E-02, 0.10079240E-01,-0.21531893E-02,
     + -0.50782710E-03,-0.41438923E-02, 0.23319860E-03,-0.68292760E-02,
     +  0.51281373E-02, 0.11921033E-01, 0.87115243E-02,-0.21659340E-02,
     +  0.34502870E-02, 0.19265110E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      t_m35_0_10  =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)                *x51
     4  +coeff( 4)*x11                
     5  +coeff( 5)    *x22            
     6  +coeff( 6)        *x32        
     7  +coeff( 7)        *x31*x41    
     8  +coeff( 8)            *x42    
      t_m35_0_10  =t_m35_0_10  
     9  +coeff( 9)    *x21        *x51
     1  +coeff(10)                *x52
     2  +coeff(11)*x11*x21            
     3  +coeff(12)    *x23            
     4  +coeff(13)    *x21*x32        
     5  +coeff(14)    *x21*x31*x41    
     6  +coeff(15)    *x21    *x42    
     7  +coeff(16)*x11            *x51
     8  +coeff(17)    *x22        *x51
      t_m35_0_10  =t_m35_0_10  
     9  +coeff(18)        *x32    *x51
     1  +coeff(19)        *x31*x41*x51
     2  +coeff(20)            *x42*x51
     3  +coeff(21)    *x21        *x52
     4  +coeff(22)                *x53
     5  +coeff(23)*x12                
     6  +coeff(24)    *x24            
     7  +coeff(25)        *x34        
     8  +coeff(26)    *x22*x31*x41    
      t_m35_0_10  =t_m35_0_10  
     9  +coeff(27)        *x33*x41    
     1  +coeff(28)*x11        *x42    
     2  +coeff(29)    *x22    *x42    
     3  +coeff(30)        *x31*x43    
     4  +coeff(31)            *x44    
     5  +coeff(32)*x11*x21        *x51
     6  +coeff(33)    *x23        *x51
     7  +coeff(34)    *x21*x31*x41*x51
     8  +coeff(35)    *x21    *x42*x51
      t_m35_0_10  =t_m35_0_10  
     9  +coeff(36)    *x22        *x52
     1  +coeff(37)        *x31*x41*x52
     2  +coeff(38)            *x42*x52
     3  +coeff(39)    *x21        *x53
     4  +coeff(40)                *x54
     5  +coeff(41)*x11*x23            
     6  +coeff(42)    *x23*x31*x41    
     7  +coeff(43)    *x21*x33*x41    
     8  +coeff(44)    *x23    *x42    
      t_m35_0_10  =t_m35_0_10  
     9  +coeff(45)    *x21*x32*x42    
     1  +coeff(46)    *x21*x31*x43    
     2  +coeff(47)    *x21    *x44    
     3  +coeff(48)    *x24        *x51
     4  +coeff(49)    *x22*x32    *x51
     5  +coeff(50)    *x22*x31*x41*x51
c
      return
      end
      function y_m35_0_10  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.99979E-04,-0.20540E+00,-0.59989E-01,-0.36252E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36252E-01, 0.14986E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.17070731E+00, 0.43631452E+00,-0.36260252E-02, 0.54245300E-02,
     + -0.73129902E-02, 0.32454550E-03,-0.10350991E-01,-0.47591902E-05,
     + -0.22666132E-01,-0.47609400E-02,-0.95952810E-02,-0.83849320E-02,
     + -0.19710281E-03, 0.98921114E-03,-0.49666820E-02,-0.17471690E-01,
     + -0.36387074E-01, 0.63378014E-02, 0.50052960E-02,-0.25074820E-02,
     +  0.74817612E-02, 0.59357741E-02,-0.25936100E-02,-0.21251000E-02,
     +  0.16067180E-01, 0.10612730E-01, 0.29419563E-03,-0.23834450E-03,
     +  0.72051220E-02, 0.14878100E-01,-0.60552184E-03, 0.10303093E-03,
     +  0.14514300E-01, 0.12205664E-01,-0.12872910E-02, 0.99087122E-03,
     +  0.27069500E-02, 0.32258080E-02,-0.11069780E-02, 0.14926850E-02,
     +  0.86225681E-02,-0.98256650E-03,-0.61772312E-02,-0.11223940E-02,
     +  0.10408370E-01,-0.13890410E-01,-0.68051160E-02,-0.37485530E-03,
     +  0.79537530E-02,-0.89461630E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      y_m35_0_10  =avdat
     1  +coeff( 1)        *x31        
     2  +coeff( 2)            *x41    
     3  +coeff( 3)    *x21*x31        
     4  +coeff( 4)    *x21    *x41    
     5  +coeff( 5)        *x31    *x51
     6  +coeff( 6)            *x41*x51
     7  +coeff( 7)    *x22*x31        
     8  +coeff( 8)*x11        *x41    
      y_m35_0_10  =y_m35_0_10  
     9  +coeff( 9)    *x22    *x41    
     1  +coeff(10)        *x32*x41    
     2  +coeff(11)        *x31*x42    
     3  +coeff(12)            *x43    
     4  +coeff(13)    *x21*x31    *x51
     5  +coeff(14)    *x21    *x41*x51
     6  +coeff(15)            *x41*x52
     7  +coeff(16)    *x23*x31        
     8  +coeff(17)    *x23    *x41    
      y_m35_0_10  =y_m35_0_10  
     9  +coeff(18)    *x21*x31*x42    
     1  +coeff(19)    *x21    *x43    
     2  +coeff(20)    *x22    *x41*x51
     3  +coeff(21)        *x31*x42*x51
     4  +coeff(22)            *x43*x51
     5  +coeff(23)    *x21    *x41*x52
     6  +coeff(24)    *x24    *x41    
     7  +coeff(25)    *x22*x31*x42    
     8  +coeff(26)    *x22    *x43    
      y_m35_0_10  =y_m35_0_10  
     9  +coeff(27)    *x22*x31    *x52
     1  +coeff(28)        *x33    *x52
     2  +coeff(29)    *x23*x32*x41    
     3  +coeff(30)    *x23*x31*x42    
     4  +coeff(31)    *x22*x32*x41*x51
     5  +coeff(32)        *x34*x41*x51
     6  +coeff(33)    *x23*x31*x42*x51
     7  +coeff(34)    *x23    *x43*x51
     8  +coeff(35)        *x31    *x52
      y_m35_0_10  =y_m35_0_10  
     9  +coeff(36)    *x21*x33        
     1  +coeff(37)    *x21*x32*x41    
     2  +coeff(38)        *x32*x41*x51
     3  +coeff(39)    *x21*x31    *x52
     4  +coeff(40)            *x41*x53
     5  +coeff(41)    *x22*x32*x41    
     6  +coeff(42)        *x31*x44    
     7  +coeff(43)    *x23    *x41*x51
     8  +coeff(44)    *x21*x32*x41*x51
      y_m35_0_10  =y_m35_0_10  
     9  +coeff(45)    *x23    *x43    
     1  +coeff(46)    *x22*x32*x43    
     2  +coeff(47)    *x22*x31*x44    
     3  +coeff(48)        *x33*x44    
     4  +coeff(49)    *x23*x32*x41*x51
     5  +coeff(50)    *x22*x33*x44    
c
      return
      end
      function p_m35_0_10  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.99979E-04,-0.20540E+00,-0.59989E-01,-0.36252E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36252E-01, 0.14986E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.51149410E-01,-0.11756840E+00, 0.72553450E-03,-0.21032420E-02,
     +  0.15924240E-01, 0.37732694E-01,-0.27064912E-03,-0.24049780E-05,
     + -0.35933300E-02,-0.17244330E-02,-0.46169070E-02,-0.51028733E-02,
     + -0.44628250E-02,-0.62480112E-02,-0.10595730E-01, 0.26303610E-02,
     +  0.48611220E-02, 0.11978110E-02, 0.56535240E-02, 0.33312754E-02,
     +  0.84300270E-03,-0.63109200E-02, 0.60242512E-02, 0.38085663E-02,
     + -0.29354020E-02,-0.77187502E-02, 0.19409040E-02, 0.31087961E-02,
     + -0.11958234E-02,-0.10343170E-02,-0.10348630E-02, 0.15928370E-02,
     + -0.26268681E-02,-0.30577233E-02,-0.29826160E-02,-0.27449650E-02,
     +  0.17219950E-02, 0.62049250E-02, 0.64865471E-02,-0.57979091E-02,
     + -0.67272860E-03, 0.21221270E-02,-0.11513190E-02,-0.23577841E-02,
     + -0.22224440E-02, 0.20756123E-02, 0.19974552E-02, 0.70711360E-02,
     +  0.70350230E-02, 0.15759953E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      p_m35_0_10  =avdat
     1  +coeff( 1)        *x31        
     2  +coeff( 2)            *x41    
     3  +coeff( 3)    *x21*x31        
     4  +coeff( 4)    *x21    *x41    
     5  +coeff( 5)        *x31    *x51
     6  +coeff( 6)            *x41*x51
     7  +coeff( 7)        *x33        
     8  +coeff( 8)*x11        *x41    
      p_m35_0_10  =p_m35_0_10  
     9  +coeff( 9)    *x22    *x41    
     1  +coeff(10)        *x32*x41    
     2  +coeff(11)        *x31*x42    
     3  +coeff(12)            *x43    
     4  +coeff(13)    *x21    *x41*x51
     5  +coeff(14)        *x31    *x52
     6  +coeff(15)            *x41*x52
     7  +coeff(16)    *x23*x31        
     8  +coeff(17)    *x23    *x41    
      p_m35_0_10  =p_m35_0_10  
     9  +coeff(18)    *x21*x32*x41    
     1  +coeff(19)    *x21*x31*x42    
     2  +coeff(20)    *x21    *x43    
     3  +coeff(21)        *x33    *x51
     4  +coeff(22)    *x22    *x41*x51
     5  +coeff(23)        *x31*x42*x51
     6  +coeff(24)            *x43*x51
     7  +coeff(25)    *x21*x31    *x52
     8  +coeff(26)    *x21    *x41*x52
      p_m35_0_10  =p_m35_0_10  
     9  +coeff(27)        *x31    *x53
     1  +coeff(28)            *x41*x53
     2  +coeff(29)    *x24*x31        
     3  +coeff(30)    *x22*x33        
     4  +coeff(31)    *x22*x32*x41    
     5  +coeff(32)    *x22    *x43    
     6  +coeff(33)        *x32*x43    
     7  +coeff(34)        *x31*x44    
     8  +coeff(35)    *x23*x31    *x51
      p_m35_0_10  =p_m35_0_10  
     9  +coeff(36)    *x23    *x41*x51
     1  +coeff(37)    *x21*x32*x41*x51
     2  +coeff(38)    *x21*x31*x42*x51
     3  +coeff(39)    *x21    *x43*x51
     4  +coeff(40)    *x22    *x41*x52
     5  +coeff(41)    *x21*x31    *x53
     6  +coeff(42)    *x21    *x41*x53
     7  +coeff(43)            *x41*x54
     8  +coeff(44)    *x21*x31*x44    
      p_m35_0_10  =p_m35_0_10  
     9  +coeff(45)    *x24*x31    *x51
     1  +coeff(46)    *x22*x32*x41*x51
     2  +coeff(47)        *x34*x41*x51
     3  +coeff(48)    *x22*x31*x42*x51
     4  +coeff(49)    *x22    *x43*x51
     5  +coeff(50)        *x32*x43*x51
c
      return
      end
      function l_m35_0_10  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(25)
      data ncoeff/ 24/
      data avdat/ -0.2724607E-01/
      data xmin/
     1 -0.99979E-04,-0.20540E+00,-0.59989E-01,-0.36252E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36252E-01, 0.14986E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.29804040E-01,-0.14093640E+00,-0.61801460E-01,-0.18402861E-01,
     + -0.27979884E-01,-0.13807413E-03, 0.56261580E-02,-0.19909810E-01,
     + -0.22745041E-01,-0.55025182E-02,-0.13729810E-01,-0.38297662E-02,
     +  0.94875684E-02, 0.10414832E-01,-0.11556380E-01, 0.29580080E-02,
     +  0.73676891E-02, 0.28740374E-02, 0.11000082E-01, 0.20267170E-02,
     +  0.23197373E-02,-0.46861222E-02, 0.59862690E-02, 0.32053280E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      l_m35_0_10  =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)    *x22            
     4  +coeff( 4)                *x51
     5  +coeff( 5)            *x42    
     6  +coeff( 6)        *x33*x41    
     7  +coeff( 7)    *x22        *x53
     8  +coeff( 8)        *x31*x41    
      l_m35_0_10  =l_m35_0_10  
     9  +coeff( 9)    *x21        *x51
     1  +coeff(10)                *x52
     2  +coeff(11)    *x21        *x52
     3  +coeff(12)        *x32        
     4  +coeff(13)        *x31*x41*x51
     5  +coeff(14)            *x42*x51
     6  +coeff(15)    *x22        *x52
     7  +coeff(16)    *x24        *x51
     8  +coeff(17)    *x21        *x53
      l_m35_0_10  =l_m35_0_10  
     9  +coeff(18)    *x23            
     1  +coeff(19)    *x22        *x51
     2  +coeff(20)        *x32    *x51
     3  +coeff(21)                *x53
     4  +coeff(22)    *x24            
     5  +coeff(23)    *x23        *x51
     6  +coeff(24)    *x21    *x42*x51
c
      return
      end
      function x_m35_0_11  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/ -0.6386141E-01/
      data xmin/
     1 -0.99979E-04,-0.20540E+00,-0.59989E-01,-0.36252E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36252E-01, 0.14986E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.59184280E-01, 0.29143690E+00, 0.28047730E+00,-0.76136420E-01,
     + -0.19182390E-01,-0.23723904E-01, 0.33822230E+00,-0.48948430E-01,
     + -0.53472190E-01,-0.20411520E-01, 0.27026330E-01, 0.93113800E-02,
     +  0.94351120E-02,-0.61015270E-01, 0.42779393E-01,-0.36151293E-01,
     + -0.39589970E-01,-0.31707484E-01, 0.13102130E-02,-0.82348290E-02,
     +  0.10063280E-01, 0.12818383E-01, 0.16173130E-01,-0.37641380E-02,
     + -0.48624030E-02,-0.13669130E-01,-0.20813540E-01, 0.48400111E-01,
     +  0.37133000E-01, 0.26129570E-02,-0.12367622E-01,-0.47829691E-02,
     + -0.46572600E-02,-0.54578610E-02,-0.24183450E-02, 0.11277833E-01,
     + -0.56345434E-02, 0.54790121E-02,-0.32258740E-03, 0.39116680E-02,
     +  0.35275504E-02, 0.13180970E-02,-0.44306251E-02,-0.50924230E-02,
     + -0.35719020E-02,-0.14232030E-01,-0.15051412E-01, 0.36917950E-02,
     +  0.29304320E-02,-0.19067890E-01,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x25 = x24*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
      x55 = x54*x5
c
c                  function
c
      x_m35_0_11  =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)                *x51
     4  +coeff( 4)    *x22            
     5  +coeff( 5)        *x31*x41    
     6  +coeff( 6)            *x42    
     7  +coeff( 7)    *x21        *x51
     8  +coeff( 8)                *x52
      x_m35_0_11  =x_m35_0_11  
     9  +coeff( 9)    *x23            
     1  +coeff(10)    *x21    *x42    
     2  +coeff(11)    *x22        *x51
     3  +coeff(12)        *x31*x41*x51
     4  +coeff(13)            *x42*x51
     5  +coeff(14)    *x21        *x52
     6  +coeff(15)    *x24            
     7  +coeff(16)    *x22*x31*x41    
     8  +coeff(17)    *x22    *x42    
      x_m35_0_11  =x_m35_0_11  
     9  +coeff(18)    *x23*x31*x41    
     1  +coeff(19)    *x22*x34        
     2  +coeff(20)    *x21*x31*x41    
     3  +coeff(21)                *x53
     4  +coeff(22)    *x23        *x51
     5  +coeff(23)    *x21        *x53
     6  +coeff(24)        *x32        
     7  +coeff(25)    *x21*x32        
     8  +coeff(26)    *x25            
      x_m35_0_11  =x_m35_0_11  
     9  +coeff(27)    *x23    *x42    
     1  +coeff(28)    *x23*x31*x43    
     2  +coeff(29)    *x23    *x44    
     3  +coeff(30)        *x32    *x51
     4  +coeff(31)    *x22*x32        
     5  +coeff(32)    *x22        *x52
     6  +coeff(33)        *x31*x41*x52
     7  +coeff(34)            *x42*x52
     8  +coeff(35)                *x54
      x_m35_0_11  =x_m35_0_11  
     9  +coeff(36)    *x21*x32*x42    
     1  +coeff(37)    *x24        *x51
     2  +coeff(38)    *x22        *x55
     3  +coeff(39)*x11                
     4  +coeff(40)        *x31*x43    
     5  +coeff(41)            *x44    
     6  +coeff(42)    *x22*x31*x41*x51
     7  +coeff(43)    *x21*x31*x41*x52
     8  +coeff(44)    *x21    *x42*x52
      x_m35_0_11  =x_m35_0_11  
     9  +coeff(45)    *x21        *x54
     1  +coeff(46)    *x24*x31*x41    
     2  +coeff(47)    *x24    *x42    
     3  +coeff(48)    *x21*x31*x41*x53
     4  +coeff(49)    *x21    *x42*x53
     5  +coeff(50)    *x25    *x42    
c
      return
      end
      function t_m35_0_11  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/ -0.2280975E-01/
      data xmin/
     1 -0.99979E-04,-0.20540E+00,-0.59989E-01,-0.36252E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36252E-01, 0.14986E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.22282732E-01, 0.36660600E-01, 0.10672840E+00,-0.12265780E-03,
     + -0.25163590E-01,-0.24388774E-02,-0.75152690E-02,-0.93502630E-02,
     +  0.10555400E+00,-0.22736081E-01, 0.20088170E-03,-0.23576230E-01,
     + -0.18461500E-02,-0.82012000E-02,-0.11588343E-01, 0.16676183E-04,
     +  0.12534320E-01, 0.32794170E-02, 0.37005630E-02,-0.29585080E-01,
     +  0.68606960E-02, 0.33677522E-04, 0.13705314E-01, 0.14164533E-03,
     + -0.15829853E-01,-0.22911683E-04,-0.14721210E-04,-0.18202190E-01,
     +  0.27072730E-02, 0.24505530E-02,-0.17557933E-03, 0.72260110E-02,
     +  0.30714162E-03, 0.17907492E-03,-0.39553040E-02,-0.22094971E-02,
     + -0.28584792E-02, 0.96255750E-02,-0.22758700E-02,-0.60916511E-03,
     + -0.39207160E-02, 0.55742402E-04,-0.68947034E-02, 0.46371673E-02,
     +  0.10151470E-01, 0.76907450E-02,-0.25865300E-02, 0.34650870E-02,
     +  0.53774470E-03, 0.27889400E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      t_m35_0_11  =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)                *x51
     4  +coeff( 4)*x11                
     5  +coeff( 5)    *x22            
     6  +coeff( 6)        *x32        
     7  +coeff( 7)        *x31*x41    
     8  +coeff( 8)            *x42    
      t_m35_0_11  =t_m35_0_11  
     9  +coeff( 9)    *x21        *x51
     1  +coeff(10)                *x52
     2  +coeff(11)*x11*x21            
     3  +coeff(12)    *x23            
     4  +coeff(13)    *x21*x32        
     5  +coeff(14)    *x21*x31*x41    
     6  +coeff(15)    *x21    *x42    
     7  +coeff(16)*x11            *x51
     8  +coeff(17)    *x22        *x51
      t_m35_0_11  =t_m35_0_11  
     9  +coeff(18)        *x31*x41*x51
     1  +coeff(19)            *x42*x51
     2  +coeff(20)    *x21        *x52
     3  +coeff(21)                *x53
     4  +coeff(22)*x12                
     5  +coeff(23)    *x24            
     6  +coeff(24)        *x34        
     7  +coeff(25)    *x22*x31*x41    
     8  +coeff(26)        *x33*x41    
      t_m35_0_11  =t_m35_0_11  
     9  +coeff(27)*x11        *x42    
     1  +coeff(28)    *x22    *x42    
     2  +coeff(29)        *x31*x43    
     3  +coeff(30)            *x44    
     4  +coeff(31)*x11*x21        *x51
     5  +coeff(32)    *x23        *x51
     6  +coeff(33)    *x21*x31*x41*x51
     7  +coeff(34)    *x21    *x42*x51
     8  +coeff(35)    *x22        *x52
      t_m35_0_11  =t_m35_0_11  
     9  +coeff(36)        *x31*x41*x52
     1  +coeff(37)            *x42*x52
     2  +coeff(38)    *x21        *x53
     3  +coeff(39)                *x54
     4  +coeff(40)*x11*x23            
     5  +coeff(41)    *x23*x31*x41    
     6  +coeff(42)    *x21*x33*x41    
     7  +coeff(43)    *x23    *x42    
     8  +coeff(44)    *x21*x32*x42    
      t_m35_0_11  =t_m35_0_11  
     9  +coeff(45)    *x21*x31*x43    
     1  +coeff(46)    *x21    *x44    
     2  +coeff(47)    *x24        *x51
     3  +coeff(48)    *x22*x32    *x51
     4  +coeff(49)        *x34    *x51
     5  +coeff(50)    *x22*x31*x41*x51
c
      return
      end
      function y_m35_0_11  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.99979E-04,-0.20540E+00,-0.59989E-01,-0.36252E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36252E-01, 0.14986E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.14420001E+00, 0.37385690E+00,-0.32502922E-02, 0.31102540E-02,
     +  0.59269141E-03, 0.18280852E-01,-0.10560850E-01, 0.12109683E-05,
     + -0.23481650E-01,-0.45487060E-02,-0.11251790E-01,-0.10349914E-01,
     + -0.10473260E-02,-0.56124630E-02,-0.45732710E-02,-0.11196553E-01,
     + -0.15539973E-01, 0.89126850E-03,-0.31165720E-01, 0.88249920E-02,
     +  0.73352763E-02,-0.29047080E-02, 0.11482050E-01, 0.90163000E-02,
     + -0.26823780E-02,-0.57690422E-02, 0.37300830E-02,-0.21422480E-02,
     +  0.10575450E-01, 0.91315340E-02,-0.43303220E-02, 0.87924273E-02,
     +  0.82050110E-02, 0.30584884E-02, 0.39390280E-02,-0.41330582E-03,
     +  0.10277980E-01, 0.68440041E-02,-0.77761670E-03,-0.67425760E-04,
     +  0.90636924E-03,-0.11693960E-03,-0.85338053E-03, 0.44520970E-02,
     +  0.52033271E-02, 0.10884184E-02, 0.37694800E-02,-0.38635030E-02,
     +  0.30868190E-02,-0.14892682E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      y_m35_0_11  =avdat
     1  +coeff( 1)        *x31        
     2  +coeff( 2)            *x41    
     3  +coeff( 3)    *x21*x31        
     4  +coeff( 4)    *x21    *x41    
     5  +coeff( 5)        *x31    *x51
     6  +coeff( 6)            *x41*x51
     7  +coeff( 7)    *x22*x31        
     8  +coeff( 8)*x11        *x41    
      y_m35_0_11  =y_m35_0_11  
     9  +coeff( 9)    *x22    *x41    
     1  +coeff(10)        *x32*x41    
     2  +coeff(11)        *x31*x42    
     3  +coeff(12)            *x43    
     4  +coeff(13)    *x21*x31    *x51
     5  +coeff(14)    *x21    *x41*x51
     6  +coeff(15)        *x31    *x52
     7  +coeff(16)            *x41*x52
     8  +coeff(17)    *x23*x31        
      y_m35_0_11  =y_m35_0_11  
     9  +coeff(18)    *x21*x33        
     1  +coeff(19)    *x23    *x41    
     2  +coeff(20)    *x21*x31*x42    
     3  +coeff(21)    *x21    *x43    
     4  +coeff(22)    *x22    *x41*x51
     5  +coeff(23)        *x31*x42*x51
     6  +coeff(24)            *x43*x51
     7  +coeff(25)    *x21*x31    *x52
     8  +coeff(26)    *x21    *x41*x52
      y_m35_0_11  =y_m35_0_11  
     9  +coeff(27)            *x41*x53
     1  +coeff(28)    *x24    *x41    
     2  +coeff(29)    *x22*x31*x42    
     3  +coeff(30)    *x22    *x43    
     4  +coeff(31)        *x31*x44    
     5  +coeff(32)    *x21*x31*x42*x51
     6  +coeff(33)    *x21    *x43*x51
     7  +coeff(34)    *x21    *x41*x53
     8  +coeff(35)    *x23*x32*x41    
      y_m35_0_11  =y_m35_0_11  
     9  +coeff(36)    *x21*x34*x41    
     1  +coeff(37)    *x23*x31*x42    
     2  +coeff(38)    *x23    *x43    
     3  +coeff(39)    *x22*x32*x41*x51
     4  +coeff(40)        *x34*x41*x51
     5  +coeff(41)        *x33    *x53
     6  +coeff(42)    *x24*x33        
     7  +coeff(43)        *x33        
     8  +coeff(44)    *x21*x32*x41    
      y_m35_0_11  =y_m35_0_11  
     9  +coeff(45)        *x32*x41*x51
     1  +coeff(46)        *x31    *x53
     2  +coeff(47)    *x22*x32*x41    
     3  +coeff(48)        *x32*x43    
     4  +coeff(49)    *x21*x32*x41*x51
     5  +coeff(50)    *x22    *x41*x52
c
      return
      end
      function p_m35_0_11  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.99979E-04,-0.20540E+00,-0.59989E-01,-0.36252E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36252E-01, 0.14986E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.55323640E-01,-0.12825980E+00, 0.59525070E-03,-0.28170880E-02,
     +  0.16931570E-01, 0.37635460E-01,-0.47655020E-03,-0.28678290E-04,
     + -0.18810560E-02,-0.26099991E-02,-0.45587252E-02,-0.48356773E-02,
     + -0.35338230E-02,-0.56840730E-02,-0.11492230E-01, 0.32969801E-02,
     +  0.62769560E-02,-0.19816250E-03, 0.14864292E-02, 0.47291520E-02,
     +  0.44043692E-02,-0.20737061E-02,-0.64186034E-02, 0.22554420E-02,
     +  0.41994890E-02,-0.59229700E-03,-0.68878463E-03,-0.18250452E-02,
     + -0.45329201E-02,-0.28177220E-02,-0.34925251E-03, 0.41870280E-03,
     + -0.32424610E-02,-0.27801794E-02,-0.11245621E-02, 0.43865470E-02,
     +  0.36329020E-02,-0.25595930E-02, 0.16596390E-02, 0.29593540E-02,
     +  0.28832670E-02, 0.20243311E-02, 0.84974850E-02, 0.60307010E-02,
     +  0.16051920E-02,-0.21927200E-02, 0.28493623E-02, 0.25830050E-02,
     +  0.37347720E-02,-0.34815530E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_m35_0_11  =avdat
     1  +coeff( 1)        *x31        
     2  +coeff( 2)            *x41    
     3  +coeff( 3)    *x21*x31        
     4  +coeff( 4)    *x21    *x41    
     5  +coeff( 5)        *x31    *x51
     6  +coeff( 6)            *x41*x51
     7  +coeff( 7)    *x22*x31        
     8  +coeff( 8)*x11        *x41    
      p_m35_0_11  =p_m35_0_11  
     9  +coeff( 9)    *x22    *x41    
     1  +coeff(10)        *x32*x41    
     2  +coeff(11)        *x31*x42    
     3  +coeff(12)            *x43    
     4  +coeff(13)    *x21    *x41*x51
     5  +coeff(14)        *x31    *x52
     6  +coeff(15)            *x41*x52
     7  +coeff(16)    *x23*x31        
     8  +coeff(17)    *x23    *x41    
      p_m35_0_11  =p_m35_0_11  
     9  +coeff(18)    *x21*x32*x41    
     1  +coeff(19)    *x21    *x43    
     2  +coeff(20)        *x31*x42*x51
     3  +coeff(21)            *x43*x51
     4  +coeff(22)    *x21*x31    *x52
     5  +coeff(23)    *x21    *x41*x52
     6  +coeff(24)        *x31    *x53
     7  +coeff(25)            *x41*x53
     8  +coeff(26)    *x22*x33        
      p_m35_0_11  =p_m35_0_11  
     9  +coeff(27)    *x22*x32*x41    
     1  +coeff(28)        *x31*x44    
     2  +coeff(29)    *x23*x31    *x51
     3  +coeff(30)    *x23    *x41*x51
     4  +coeff(31)    *x21*x31*x42*x51
     5  +coeff(32)    *x21    *x43*x51
     6  +coeff(33)    *x22    *x41*x52
     7  +coeff(34)        *x31*x42*x52
     8  +coeff(35)            *x43*x52
      p_m35_0_11  =p_m35_0_11  
     9  +coeff(36)    *x21    *x41*x53
     1  +coeff(37)    *x21*x33*x42    
     2  +coeff(38)    *x24    *x41*x51
     3  +coeff(39)    *x22*x32*x41*x51
     4  +coeff(40)        *x34*x41*x51
     5  +coeff(41)    *x22*x31*x42*x51
     6  +coeff(42)        *x33*x42*x51
     7  +coeff(43)    *x23*x31*x42*x51
     8  +coeff(44)    *x23    *x43*x51
      p_m35_0_11  =p_m35_0_11  
     9  +coeff(45)    *x21*x32*x43*x51
     1  +coeff(46)        *x32*x43*x52
     2  +coeff(47)    *x23*x31    *x53
     3  +coeff(48)    *x21*x34*x43    
     4  +coeff(49)    *x24    *x43*x51
     5  +coeff(50)        *x33        
c
      return
      end
      function l_m35_0_11  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(25)
      data ncoeff/ 24/
      data avdat/ -0.2850040E-01/
      data xmin/
     1 -0.99979E-04,-0.20540E+00,-0.59989E-01,-0.36252E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36252E-01, 0.14986E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.30730854E-01,-0.13906392E+00,-0.61425500E-01,-0.30943540E-01,
     + -0.22464234E-01,-0.73864054E-02,-0.18725110E-01,-0.22726520E-01,
     +  0.13645311E-01,-0.21007720E-01,-0.47361692E-02, 0.12966840E-01,
     +  0.14519363E-01, 0.41470870E-02,-0.16719270E-01, 0.28627160E-02,
     +  0.10564601E-01, 0.86677140E-02,-0.50365480E-02, 0.53079350E-02,
     +  0.34917511E-02,-0.53989280E-02,-0.61457972E-02, 0.74568320E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      l_m35_0_11  =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)    *x22            
     4  +coeff( 4)            *x42    
     5  +coeff( 5)        *x31*x41    
     6  +coeff( 6)                *x52
     7  +coeff( 7)                *x51
     8  +coeff( 8)    *x21        *x51
      l_m35_0_11  =l_m35_0_11  
     9  +coeff( 9)    *x22        *x51
     1  +coeff(10)    *x21        *x52
     2  +coeff(11)        *x32        
     3  +coeff(12)        *x31*x41*x51
     4  +coeff(13)            *x42*x51
     5  +coeff(14)                *x53
     6  +coeff(15)    *x22        *x52
     7  +coeff(16)        *x32    *x51
     8  +coeff(17)    *x21        *x53
      l_m35_0_11  =l_m35_0_11  
     9  +coeff(18)    *x23    *x42*x51
     1  +coeff(19)    *x24            
     2  +coeff(20)    *x23        *x51
     3  +coeff(21)    *x21*x31*x41*x51
     4  +coeff(22)        *x31*x41*x52
     5  +coeff(23)            *x42*x52
     6  +coeff(24)    *x22        *x53
c
      return
      end
      function x_m35_0_12  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(101)
      data ncoeff/100/
      data avdat/ -0.7622619E-01/
      data xmin/
     1 -0.99979E-04,-0.20540E+00,-0.59989E-01,-0.36252E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36252E-01, 0.14986E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.71474604E-01, 0.31038010E+00, 0.33520130E+00,-0.36285540E-03,
     + -0.91502440E-01, 0.39288282E+00,-0.22900060E-01,-0.30541530E-01,
     + -0.60821663E-01,-0.79993430E-04,-0.81953153E-01, 0.25955440E-01,
     + -0.30673373E-01,-0.31788690E-01,-0.70335500E-01, 0.54516904E-01,
     +  0.11333893E-01, 0.10819791E-01, 0.12019030E-01, 0.16986021E-01,
     + -0.11039870E-01,-0.36061330E-01,-0.28190870E-01, 0.18040940E-01,
     +  0.17034120E+00,-0.56450571E-02,-0.40755490E-02, 0.26505020E-01,
     +  0.56068141E-01, 0.66718380E-02, 0.37267781E-01, 0.16896490E-01,
     +  0.24717680E-02, 0.49283500E-02, 0.68338760E-02,-0.61346711E-02,
     + -0.74707311E-02, 0.22007330E-01,-0.33311380E-01,-0.53973120E-01,
     + -0.21600890E-01, 0.64948224E-02,-0.77124270E-02,-0.63177384E-01,
     + -0.11916594E+00, 0.11583930E-02,-0.29210880E-01, 0.33744750E-02,
     +  0.67262310E-02, 0.15776280E-01,-0.15590162E-02,-0.42507810E-01,
     + -0.76786614E-02,-0.22084521E-01,-0.13677510E-01,-0.10035353E+00,
     + -0.64420122E-02,-0.62383571E-02,-0.10265723E-01,-0.72223720E-02,
     + -0.13074700E-01,-0.50874054E-02, 0.43054310E-01, 0.29519790E-02,
     + -0.22582190E+00, 0.20601330E+00,-0.81591784E-01,-0.11696234E-03,
     +  0.24323570E-01, 0.10520482E-01, 0.46826012E-01, 0.89983050E-01,
     + -0.37240570E-01,-0.57111834E-02,-0.55963110E-01,-0.98185580E-02,
     +  0.26009150E-02, 0.13245400E-01, 0.20413500E-01,-0.38927093E-01,
     +  0.82042840E-03, 0.44398230E-01, 0.48203870E-01, 0.16847260E+00,
     +  0.33046350E-01, 0.11204670E-01, 0.63575893E-01, 0.35892270E-01,
     + -0.11363570E-01, 0.34038912E-01, 0.13106603E-01, 0.75773842E-03,
     + -0.23247881E-02,-0.50272050E-02,-0.11013090E-01, 0.28778391E-02,
     + -0.68687082E-03,-0.33090484E-02,-0.26461214E-02,-0.36960111E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x25 = x24*x2
      x26 = x25*x2
      x27 = x26*x2
      x28 = x27*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x35 = x34*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x45 = x44*x4
      x46 = x45*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
      x55 = x54*x5
c
c                  function
c
      x_m35_0_12  =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)                *x51
     4  +coeff( 4)*x11                
     5  +coeff( 5)    *x22            
     6  +coeff( 6)    *x21        *x51
     7  +coeff( 7)        *x31*x41    
     8  +coeff( 8)            *x42    
      x_m35_0_12  =x_m35_0_12  
     9  +coeff( 9)                *x52
     1  +coeff(10)*x11*x21            
     2  +coeff(11)    *x23            
     3  +coeff(12)    *x22        *x51
     4  +coeff(13)    *x21*x31*x41    
     5  +coeff(14)    *x21    *x42    
     6  +coeff(15)    *x21        *x52
     7  +coeff(16)    *x24            
     8  +coeff(17)        *x31*x41*x51
      x_m35_0_12  =x_m35_0_12  
     9  +coeff(18)            *x42*x51
     1  +coeff(19)                *x53
     2  +coeff(20)    *x23        *x51
     3  +coeff(21)    *x22*x32        
     4  +coeff(22)    *x22*x31*x41    
     5  +coeff(23)    *x22    *x42    
     6  +coeff(24)    *x21        *x53
     7  +coeff(25)    *x25    *x42    
     8  +coeff(26)        *x32        
      x_m35_0_12  =x_m35_0_12  
     9  +coeff(27)    *x21*x32        
     1  +coeff(28)    *x25            
     2  +coeff(29)    *x23*x31*x41    
     3  +coeff(30)    *x22        *x53
     4  +coeff(31)    *x21*x31*x43    
     5  +coeff(32)    *x21    *x44    
     6  +coeff(33)        *x32    *x51
     7  +coeff(34)        *x31*x43    
     8  +coeff(35)            *x44    
      x_m35_0_12  =x_m35_0_12  
     9  +coeff(36)        *x31*x41*x52
     1  +coeff(37)            *x42*x52
     2  +coeff(38)    *x21*x32*x42    
     3  +coeff(39)    *x24*x31*x41    
     4  +coeff(40)    *x24    *x42    
     5  +coeff(41)    *x27            
     6  +coeff(42)    *x23*x31*x41*x51
     7  +coeff(43)    *x22        *x54
     8  +coeff(44)    *x24    *x42*x51
      x_m35_0_12  =x_m35_0_12  
     9  +coeff(45)    *x23*x31*x43    
     1  +coeff(46)    *x23    *x44    
     2  +coeff(47)    *x23    *x42*x52
     3  +coeff(48)    *x25*x31*x41*x52
     4  +coeff(49)        *x32*x42    
     5  +coeff(50)    *x24        *x51
     6  +coeff(51)        *x32    *x52
     7  +coeff(52)    *x23    *x42    
     8  +coeff(53)    *x21        *x54
      x_m35_0_12  =x_m35_0_12  
     9  +coeff(54)    *x22    *x44    
     1  +coeff(55)    *x26        *x51
     2  +coeff(56)    *x25*x31*x41    
     3  +coeff(57)    *x25        *x52
     4  +coeff(58)        *x32*x44    
     5  +coeff(59)    *x23*x32    *x52
     6  +coeff(60)    *x22    *x44*x51
     7  +coeff(61)    *x21*x32*x44    
     8  +coeff(62)    *x21*x31*x45    
      x_m35_0_12  =x_m35_0_12  
     9  +coeff(63)    *x24    *x44    
     1  +coeff(64)            *x42*x55
     2  +coeff(65)    *x27    *x42    
     3  +coeff(66)    *x25*x31*x43    
     4  +coeff(67)    *x25    *x44    
     5  +coeff(68)    *x21    *x46*x51
     6  +coeff(69)    *x25    *x42*x52
     7  +coeff(70)    *x28*x31*x41    
     8  +coeff(71)    *x24*x31*x43*x51
      x_m35_0_12  =x_m35_0_12  
     9  +coeff(72)    *x24    *x44*x51
     1  +coeff(73)    *x24*x31*x41*x53
     2  +coeff(74)    *x23*x34*x42    
     3  +coeff(75)    *x27    *x42*x51
     4  +coeff(76)    *x23*x31*x41*x54
     5  +coeff(77)    *x22*x35*x41*x51
     6  +coeff(78)    *x22*x34*x42*x51
     7  +coeff(79)    *x22*x33*x43*x51
     8  +coeff(80)    *x25*x33*x41*x51
      x_m35_0_12  =x_m35_0_12  
     9  +coeff(81)    *x25*x32*x42*x51
     1  +coeff(82)    *x25    *x44*x51
     2  +coeff(83)    *x27*x32*x42    
     3  +coeff(84)    *x27    *x44    
     4  +coeff(85)    *x23*x33*x43*x51
     5  +coeff(86)    *x23    *x42*x55
     6  +coeff(87)    *x26*x32*x42*x51
     7  +coeff(88)    *x26*x31*x43*x51
     8  +coeff(89)    *x26*x32    *x53
      x_m35_0_12  =x_m35_0_12  
     9  +coeff(90)    *x24*x31*x41*x55
     1  +coeff(91)    *x23*x34    *x54
     2  +coeff(92)        *x34        
     3  +coeff(93)                *x54
     4  +coeff(94)    *x23*x32        
     5  +coeff(95)    *x22*x31*x41*x51
     6  +coeff(96)    *x21*x33*x41    
     7  +coeff(97)*x11*x23        *x51
     8  +coeff(98)    *x21*x31*x41*x52
      x_m35_0_12  =x_m35_0_12  
     9  +coeff(99)    *x24*x32        
     1  +coeff(100)    *x24        *x52
c
      return
      end
      function t_m35_0_12  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(101)
      data ncoeff/100/
      data avdat/ -0.2211770E-01/
      data xmin/
     1 -0.99979E-04,-0.20540E+00,-0.59989E-01,-0.36252E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36252E-01, 0.14986E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.21461034E-01, 0.35432342E-01, 0.10848430E+00,-0.25532713E-01,
     + -0.99079180E-02, 0.10815870E+00,-0.23623954E-01,-0.77098560E-02,
     + -0.15469120E-02,-0.19542270E-01,-0.14024680E-01, 0.12777260E-01,
     +  0.22920714E-02,-0.24537030E-01, 0.53088492E-02,-0.87082032E-02,
     +  0.18730781E-02,-0.15047010E-02, 0.57050280E-03, 0.15518323E-01,
     + -0.88549982E-02, 0.32571430E-02, 0.45727253E-02,-0.19567254E-03,
     + -0.14450460E-02,-0.23946670E-02, 0.75307302E-02,-0.13722840E-02,
     + -0.10251214E-01, 0.69831321E-02, 0.28530010E-03,-0.23869270E-02,
     + -0.97947550E-05,-0.42592590E-02, 0.43904050E-02,-0.58657500E-03,
     + -0.17797880E-02, 0.28373163E-01, 0.89098350E-02,-0.23963940E-02,
     + -0.13525770E-01,-0.11164430E-02,-0.34060380E-02,-0.62462482E-02,
     +  0.17380523E-02, 0.20627050E-02,-0.34204791E-02, 0.14152783E-01,
     +  0.12328060E-01,-0.65466412E-03,-0.36133830E-02, 0.18298723E-02,
     + -0.17438001E-02, 0.99973184E-02,-0.16672950E-01,-0.81373723E-02,
     + -0.28678620E-02,-0.17739923E-02,-0.26975840E-02,-0.10701421E-01,
     + -0.73749600E-02,-0.32096060E-02,-0.49828620E-01,-0.35430662E-01,
     +  0.17062130E-01,-0.27723320E-01,-0.38552330E-01, 0.17843460E-01,
     + -0.57608741E-02,-0.14027804E-01,-0.60011470E-02,-0.23537040E-03,
     +  0.76184570E-02, 0.16284111E-01,-0.16052500E-01, 0.41791661E-02,
     +  0.12710770E-01, 0.19818183E-01, 0.59821230E-02,-0.81616680E-02,
     +  0.55980172E-01, 0.63701380E-01, 0.50879451E-02, 0.22376670E-01,
     +  0.21734900E-01,-0.10369630E-03,-0.37476043E-02, 0.18842840E-03,
     +  0.22151000E-02,-0.41393013E-02,-0.81815421E-02,-0.55104433E-02,
     + -0.48945430E-03, 0.37153670E-03,-0.43244304E-04,-0.12944660E-02,
     +  0.23312843E-02, 0.48834942E-02,-0.29117332E-02, 0.46408553E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x25 = x24*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x45 = x44*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
      x55 = x54*x5
c
c                  function
c
      t_m35_0_12  =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)                *x51
     4  +coeff( 4)    *x22            
     5  +coeff( 5)            *x42    
     6  +coeff( 6)    *x21        *x51
     7  +coeff( 7)                *x52
     8  +coeff( 8)        *x31*x41    
      t_m35_0_12  =t_m35_0_12  
     9  +coeff( 9)        *x32        
     1  +coeff(10)    *x23            
     2  +coeff(11)    *x21    *x42    
     3  +coeff(12)    *x22        *x51
     4  +coeff(13)            *x42*x51
     5  +coeff(14)    *x21        *x52
     6  +coeff(15)                *x53
     7  +coeff(16)    *x21*x31*x41    
     8  +coeff(17)        *x31*x41*x51
      t_m35_0_12  =t_m35_0_12  
     9  +coeff(18)    *x21*x32        
     1  +coeff(19)        *x32    *x51
     2  +coeff(20)    *x24            
     3  +coeff(21)    *x22    *x42    
     4  +coeff(22)            *x44    
     5  +coeff(23)    *x23        *x51
     6  +coeff(24)    *x21    *x42*x51
     7  +coeff(25)    *x22        *x52
     8  +coeff(26)            *x42*x52
      t_m35_0_12  =t_m35_0_12  
     9  +coeff(27)    *x21        *x53
     1  +coeff(28)                *x54
     2  +coeff(29)    *x22*x31*x41    
     3  +coeff(30)        *x31*x43    
     4  +coeff(31)    *x21*x31*x41*x51
     5  +coeff(32)        *x31*x41*x52
     6  +coeff(33)*x11*x22            
     7  +coeff(34)    *x22*x32        
     8  +coeff(35)        *x32*x42    
      t_m35_0_12  =t_m35_0_12  
     9  +coeff(36)        *x32    *x52
     1  +coeff(37)    *x25            
     2  +coeff(38)    *x23    *x42    
     3  +coeff(39)    *x21    *x44    
     4  +coeff(40)    *x24        *x51
     5  +coeff(41)    *x22    *x42*x51
     6  +coeff(42)            *x44*x51
     7  +coeff(43)    *x23        *x52
     8  +coeff(44)    *x21    *x42*x52
      t_m35_0_12  =t_m35_0_12  
     9  +coeff(45)    *x22        *x53
     1  +coeff(46)            *x42*x53
     2  +coeff(47)    *x21        *x54
     3  +coeff(48)    *x23*x31*x41    
     4  +coeff(49)    *x21*x31*x43    
     5  +coeff(50)        *x31*x43*x51
     6  +coeff(51)    *x21*x31*x41*x52
     7  +coeff(52)        *x31*x41*x53
     8  +coeff(53)    *x23*x32        
      t_m35_0_12  =t_m35_0_12  
     9  +coeff(54)    *x21*x32*x42    
     1  +coeff(55)    *x24    *x42    
     2  +coeff(56)    *x22    *x44    
     3  +coeff(57)    *x21    *x44*x51
     4  +coeff(58)    *x22    *x42*x52
     5  +coeff(59)    *x22        *x54
     6  +coeff(60)    *x24*x31*x41    
     7  +coeff(61)    *x22*x31*x43    
     8  +coeff(62)    *x23*x31*x41*x51
      t_m35_0_12  =t_m35_0_12  
     9  +coeff(63)    *x25    *x42    
     1  +coeff(64)    *x23    *x44    
     2  +coeff(65)    *x22    *x44*x51
     3  +coeff(66)    *x25*x31*x41    
     4  +coeff(67)    *x23*x31*x43    
     5  +coeff(68)    *x22*x31*x43*x51
     6  +coeff(69)    *x22*x31*x41*x53
     7  +coeff(70)    *x23*x32*x42    
     8  +coeff(71)    *x21*x32*x44    
      t_m35_0_12  =t_m35_0_12  
     9  +coeff(72)    *x24*x32    *x51
     1  +coeff(73)    *x22*x32*x42*x51
     2  +coeff(74)    *x24    *x44    
     3  +coeff(75)    *x25    *x42*x51
     4  +coeff(76)    *x21    *x44*x53
     5  +coeff(77)    *x24*x31*x43    
     6  +coeff(78)    *x23*x31*x43*x51
     7  +coeff(79)    *x21*x31*x43*x53
     8  +coeff(80)    *x23*x33*x41*x51
      t_m35_0_12  =t_m35_0_12  
     9  +coeff(81)    *x25    *x44    
     1  +coeff(82)    *x25*x31*x43    
     2  +coeff(83)    *x22*x31*x41*x55
     3  +coeff(84)    *x25*x32*x42    
     4  +coeff(85)    *x25    *x44*x51
     5  +coeff(86)*x11                
     6  +coeff(87)    *x22*x31*x41*x51
     7  +coeff(88)    *x21*x33*x41    
     8  +coeff(89)    *x25        *x51
      t_m35_0_12  =t_m35_0_12  
     9  +coeff(90)        *x31*x45    
     1  +coeff(91)    *x21*x31*x43*x51
     2  +coeff(92)        *x32*x44    
     3  +coeff(93)*x11*x23        *x51
     4  +coeff(94)*x11*x21        *x53
     5  +coeff(95)    *x22*x33*x41    
     6  +coeff(96)        *x33*x43    
     7  +coeff(97)    *x21*x33*x41*x51
     8  +coeff(98)    *x21    *x44*x52
      t_m35_0_12  =t_m35_0_12  
     9  +coeff(99)    *x21*x31*x45    
     1  +coeff(100)    *x21*x31*x43*x52
c
      return
      end
      function y_m35_0_12  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(101)
      data ncoeff/100/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.99979E-04,-0.20540E+00,-0.59989E-01,-0.36252E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36252E-01, 0.14986E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.11603351E+00, 0.30825942E+00,-0.28874021E-02, 0.21506760E-02,
     +  0.90564191E-02, 0.37996420E-01,-0.11540853E-01,-0.16411470E-05,
     + -0.19427830E-01,-0.45428834E-02,-0.11435901E-01,-0.10950940E-01,
     + -0.31741950E-02,-0.45860330E-02,-0.61150430E-02,-0.13206601E-01,
     + -0.15256540E-01, 0.91032433E-03, 0.18889600E-03,-0.28959821E-01,
     +  0.77026202E-02, 0.69725080E-02,-0.34288240E-02, 0.14053530E-01,
     +  0.10420394E-01,-0.41828202E-02,-0.10037980E-01, 0.24801234E-02,
     +  0.54979900E-02,-0.90234093E-02, 0.15743010E-01, 0.10761320E-02,
     + -0.57741561E-02, 0.99104480E-03, 0.82605700E-02,-0.10249660E-03,
     + -0.15411633E-01,-0.63464330E-02,-0.52120303E-02, 0.71548302E-02,
     +  0.14484890E-02,-0.22565361E-02, 0.24747373E-01, 0.79856542E-02,
     +  0.21506862E-02, 0.19276393E-03, 0.63540083E-02,-0.37572350E-02,
     + -0.28381021E-02, 0.26422000E-01, 0.16464831E-01,-0.10453154E-02,
     +  0.10265330E-01,-0.68204480E-03,-0.15491060E-01,-0.11401722E-02,
     +  0.46171220E-02, 0.50225200E-03, 0.95604202E-03, 0.56096450E-02,
     +  0.69829123E-02,-0.59812860E-02,-0.23304540E-02,-0.25508841E-02,
     +  0.41165770E-02, 0.29425870E-02,-0.19204802E-01,-0.84464633E-02,
     +  0.24849840E-02,-0.13029860E-01, 0.10374700E-01,-0.56888074E-02,
     + -0.51866653E-02,-0.14831490E-02, 0.77607990E-02,-0.16508510E-01,
     +  0.16108792E-01, 0.99704451E-02,-0.14491181E-02,-0.18698643E-01,
     +  0.45328740E-01, 0.25561350E-01, 0.77952053E-02, 0.48790900E-02,
     + -0.77702663E-02, 0.22807192E-01, 0.11872880E-01,-0.10869100E-01,
     + -0.65644500E-02, 0.13651980E-01, 0.16836560E-01, 0.42971721E-02,
     + -0.39635004E-03,-0.56505524E-02,-0.27425610E-02,-0.52001234E-03,
     + -0.58621651E-03, 0.29779672E-02,-0.42802700E-02,-0.77879300E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      y_m35_0_12  =avdat
     1  +coeff( 1)        *x31        
     2  +coeff( 2)            *x41    
     3  +coeff( 3)    *x21*x31        
     4  +coeff( 4)    *x21    *x41    
     5  +coeff( 5)        *x31    *x51
     6  +coeff( 6)            *x41*x51
     7  +coeff( 7)    *x22*x31        
     8  +coeff( 8)*x11        *x41    
      y_m35_0_12  =y_m35_0_12  
     9  +coeff( 9)    *x22    *x41    
     1  +coeff(10)        *x32*x41    
     2  +coeff(11)        *x31*x42    
     3  +coeff(12)            *x43    
     4  +coeff(13)    *x21*x31    *x51
     5  +coeff(14)    *x21    *x41*x51
     6  +coeff(15)        *x31    *x52
     7  +coeff(16)            *x41*x52
     8  +coeff(17)    *x23*x31        
      y_m35_0_12  =y_m35_0_12  
     9  +coeff(18)    *x21*x33        
     1  +coeff(19)*x11*x21    *x41    
     2  +coeff(20)    *x23    *x41    
     3  +coeff(21)    *x21*x31*x42    
     4  +coeff(22)    *x21    *x43    
     5  +coeff(23)    *x22    *x41*x51
     6  +coeff(24)        *x31*x42*x51
     7  +coeff(25)            *x43*x51
     8  +coeff(26)    *x21*x31    *x52
      y_m35_0_12  =y_m35_0_12  
     9  +coeff(27)    *x21    *x41*x52
     1  +coeff(28)        *x31    *x53
     2  +coeff(29)            *x41*x53
     3  +coeff(30)    *x24    *x41    
     4  +coeff(31)    *x22*x31*x42    
     5  +coeff(32)    *x22    *x43    
     6  +coeff(33)        *x31*x44    
     7  +coeff(34)    *x21*x32*x41*x51
     8  +coeff(35)    *x21*x31*x42*x51
      y_m35_0_12  =y_m35_0_12  
     9  +coeff(36)    *x21    *x43*x51
     1  +coeff(37)    *x22    *x41*x52
     2  +coeff(38)        *x31*x42*x52
     3  +coeff(39)            *x43*x52
     4  +coeff(40)    *x21    *x41*x53
     5  +coeff(41)    *x23*x32*x41    
     6  +coeff(42)    *x21*x34*x41    
     7  +coeff(43)    *x23*x31*x42    
     8  +coeff(44)    *x23    *x43    
      y_m35_0_12  =y_m35_0_12  
     9  +coeff(45)    *x22*x32*x41*x51
     1  +coeff(46)        *x34*x41*x51
     2  +coeff(47)    *x23    *x41*x52
     3  +coeff(48)    *x24*x33        
     4  +coeff(49)    *x22*x34*x41    
     5  +coeff(50)    *x23*x31*x42*x51
     6  +coeff(51)    *x23    *x43*x51
     7  +coeff(52)        *x32*x43*x52
     8  +coeff(53)    *x23*x34*x41    
      y_m35_0_12  =y_m35_0_12  
     9  +coeff(54)    *x24    *x43*x53
     1  +coeff(55)    *x24*x33*x44*x51
     2  +coeff(56)        *x33        
     3  +coeff(57)    *x21*x32*x41    
     4  +coeff(58)    *x22*x31    *x51
     5  +coeff(59)        *x33    *x51
     6  +coeff(60)        *x32*x41*x51
     7  +coeff(61)    *x22*x32*x41    
     8  +coeff(62)        *x32*x43    
      y_m35_0_12  =y_m35_0_12  
     9  +coeff(63)    *x23*x31    *x51
     1  +coeff(64)        *x32*x41*x52
     2  +coeff(65)    *x21*x31    *x53
     3  +coeff(66)    *x21*x31*x44    
     4  +coeff(67)    *x22*x31*x42*x51
     5  +coeff(68)    *x22    *x43*x51
     6  +coeff(69)    *x23*x31    *x52
     7  +coeff(70)    *x22*x33*x42    
     8  +coeff(71)    *x24    *x43    
      y_m35_0_12  =y_m35_0_12  
     9  +coeff(72)    *x22*x32*x43    
     1  +coeff(73)    *x22*x31*x44    
     2  +coeff(74)        *x33*x44    
     3  +coeff(75)    *x23*x32*x41*x51
     4  +coeff(76)    *x21*x31*x44*x51
     5  +coeff(77)    *x24    *x41*x52
     6  +coeff(78)    *x22    *x43*x52
     7  +coeff(79)    *x23*x31    *x53
     8  +coeff(80)    *x23*x31*x44    
      y_m35_0_12  =y_m35_0_12  
     9  +coeff(81)    *x24*x31*x42*x51
     1  +coeff(82)    *x24    *x43*x51
     2  +coeff(83)    *x24    *x41*x53
     3  +coeff(84)    *x24*x33*x42    
     4  +coeff(85)    *x21*x34*x43*x51
     5  +coeff(86)    *x24*x31*x42*x52
     6  +coeff(87)    *x22*x32*x43*x52
     7  +coeff(88)    *x23*x34*x43    
     8  +coeff(89)    *x23*x34*x41*x52
      y_m35_0_12  =y_m35_0_12  
     9  +coeff(90)    *x22*x33*x42*x53
     1  +coeff(91)    *x24*x33*x42*x52
     2  +coeff(92)    *x22*x33        
     3  +coeff(93)*x11*x21    *x41*x51
     4  +coeff(94)    *x23    *x41*x51
     5  +coeff(95)    *x22*x31    *x52
     6  +coeff(96)            *x41*x54
     7  +coeff(97)*x11*x23    *x41    
     8  +coeff(98)    *x21*x32*x43    
      y_m35_0_12  =y_m35_0_12  
     9  +coeff(99)    *x24*x31    *x51
     1  +coeff(100)    *x24    *x41*x51
c
      return
      end
      function p_m35_0_12  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(101)
      data ncoeff/100/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.99979E-04,-0.20540E+00,-0.59989E-01,-0.36252E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36252E-01, 0.14986E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.55773694E-01,-0.12957051E+00, 0.31217320E-03,-0.36463000E-02,
     +  0.17228862E-01, 0.38720190E-01,-0.34064280E-04,-0.78973510E-03,
     + -0.25967734E-04,-0.19619041E-02,-0.15908510E-02,-0.36923240E-02,
     + -0.40968232E-02,-0.25153502E-02,-0.51927730E-02,-0.10113963E-01,
     +  0.33792983E-02, 0.74678140E-02, 0.12993780E-02, 0.22966360E-02,
     +  0.66142482E-02, 0.29209290E-02,-0.10485331E-02,-0.56114932E-02,
     +  0.14428860E-02, 0.24792940E-02,-0.27828600E-02,-0.15669280E-02,
     + -0.35498414E-02, 0.14093540E-02,-0.81797420E-03,-0.38296790E-02,
     + -0.31632133E-02,-0.24617740E-02, 0.32561973E-02, 0.15183773E-02,
     +  0.27809040E-02, 0.11790590E-03, 0.31380774E-03, 0.92027010E-02,
     +  0.80741383E-02,-0.20731550E-03,-0.64166010E-03,-0.32911580E-03,
     +  0.20839920E-02, 0.43956890E-03, 0.23381580E-02,-0.23630620E-03,
     + -0.75083144E-03,-0.60482794E-03,-0.25632532E-02, 0.36613603E-05,
     + -0.31117623E-03,-0.92644774E-03, 0.16608570E-02,-0.40757813E-03,
     + -0.13930860E-02, 0.28657060E-02,-0.16442431E-02, 0.49435342E-02,
     + -0.24589514E-02, 0.67505780E-03, 0.13606140E-02, 0.17889371E-02,
     +  0.26787940E-02, 0.48878170E-02,-0.59079710E-03,-0.14652081E-02,
     + -0.28021503E-02,-0.90464204E-02, 0.26559580E-02,-0.10926280E-02,
     +  0.25133993E-02,-0.19331800E-02,-0.11833970E-02,-0.54042004E-02,
     + -0.11607690E-03, 0.30784610E-03,-0.51314354E-03,-0.26831690E-03,
     +  0.23180264E-02, 0.19361910E-02, 0.40252931E-03, 0.88617760E-03,
     +  0.12677310E-02, 0.34345814E-02, 0.20336430E-02,-0.89750811E-03,
     +  0.25308064E-02,-0.46309482E-03,-0.27596114E-02, 0.44744390E-02,
     +  0.55643680E-03, 0.21119760E-02, 0.25089550E-02,-0.30658970E-03,
     + -0.67233160E-03, 0.71727973E-03,-0.26641192E-02,-0.58560470E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      p_m35_0_12  =avdat
     1  +coeff( 1)        *x31        
     2  +coeff( 2)            *x41    
     3  +coeff( 3)    *x21*x31        
     4  +coeff( 4)    *x21    *x41    
     5  +coeff( 5)        *x31    *x51
     6  +coeff( 6)            *x41*x51
     7  +coeff( 7)*x11    *x31        
     8  +coeff( 8)    *x22*x31        
      p_m35_0_12  =p_m35_0_12  
     9  +coeff( 9)*x11        *x41    
     1  +coeff(10)    *x22    *x41    
     2  +coeff(11)        *x32*x41    
     3  +coeff(12)        *x31*x42    
     4  +coeff(13)            *x43    
     5  +coeff(14)    *x21    *x41*x51
     6  +coeff(15)        *x31    *x52
     7  +coeff(16)            *x41*x52
     8  +coeff(17)    *x23*x31        
      p_m35_0_12  =p_m35_0_12  
     9  +coeff(18)    *x23    *x41    
     1  +coeff(19)    *x21*x32*x41    
     2  +coeff(20)    *x21    *x43    
     3  +coeff(21)        *x31*x42*x51
     4  +coeff(22)            *x43*x51
     5  +coeff(23)    *x21*x31    *x52
     6  +coeff(24)    *x21    *x41*x52
     7  +coeff(25)        *x31    *x53
     8  +coeff(26)            *x41*x53
      p_m35_0_12  =p_m35_0_12  
     9  +coeff(27)        *x31*x44    
     1  +coeff(28)    *x23*x31    *x51
     2  +coeff(29)    *x23    *x41*x51
     3  +coeff(30)    *x21*x31*x42*x51
     4  +coeff(31)    *x21    *x43*x51
     5  +coeff(32)    *x22    *x41*x52
     6  +coeff(33)        *x31*x42*x52
     7  +coeff(34)            *x43*x52
     8  +coeff(35)    *x21    *x41*x53
      p_m35_0_12  =p_m35_0_12  
     9  +coeff(36)    *x21*x31*x44    
     1  +coeff(37)    *x22*x32*x41*x51
     2  +coeff(38)        *x34*x41*x51
     3  +coeff(39)    *x24*x33        
     4  +coeff(40)    *x23*x31*x42*x51
     5  +coeff(41)    *x23    *x43*x51
     6  +coeff(42)        *x32*x43*x52
     7  +coeff(43)    *x22*x34*x41*x51
     8  +coeff(44)        *x33        
      p_m35_0_12  =p_m35_0_12  
     9  +coeff(45)    *x21*x31*x42    
     1  +coeff(46)        *x33    *x51
     2  +coeff(47)        *x32*x41*x51
     3  +coeff(48)    *x22*x32*x41    
     4  +coeff(49)        *x33*x42    
     5  +coeff(50)    *x22    *x43    
     6  +coeff(51)        *x32*x43    
     7  +coeff(52)    *x21*x32*x41*x51
     8  +coeff(53)    *x22*x31    *x52
      p_m35_0_12  =p_m35_0_12  
     9  +coeff(54)        *x32*x41*x52
     1  +coeff(55)    *x21*x31    *x53
     2  +coeff(56)            *x41*x54
     3  +coeff(57)    *x23*x32*x41    
     4  +coeff(58)    *x23*x31*x42    
     5  +coeff(59)    *x23    *x43    
     6  +coeff(60)    *x22    *x43*x51
     7  +coeff(61)        *x31*x44*x51
     8  +coeff(62)    *x24*x31*x42    
      p_m35_0_12  =p_m35_0_12  
     9  +coeff(63)    *x24    *x43    
     1  +coeff(64)    *x22*x31*x44    
     2  +coeff(65)    *x23*x32*x41*x51
     3  +coeff(66)    *x22    *x43*x52
     4  +coeff(67)    *x23*x31    *x53
     5  +coeff(68)    *x22*x31    *x54
     6  +coeff(69)    *x23*x32*x43    
     7  +coeff(70)    *x23*x31*x44    
     8  +coeff(71)    *x24*x31*x42*x51
      p_m35_0_12  =p_m35_0_12  
     9  +coeff(72)    *x22    *x43*x53
     1  +coeff(73)        *x31*x44*x53
     2  +coeff(74)    *x21*x31    *x51
     3  +coeff(75)    *x22*x31    *x51
     4  +coeff(76)    *x22    *x41*x51
     5  +coeff(77)*x11*x21*x31    *x51
     6  +coeff(78)    *x21*x33    *x51
     7  +coeff(79)*x11*x21    *x41*x51
     8  +coeff(80)    *x24*x31    *x51
      p_m35_0_12  =p_m35_0_12  
     9  +coeff(81)    *x24    *x41*x51
     1  +coeff(82)    *x22*x31*x42*x51
     2  +coeff(83)    *x23    *x41*x52
     3  +coeff(84)    *x21    *x43*x52
     4  +coeff(85)    *x22*x31    *x53
     5  +coeff(86)    *x22    *x41*x53
     6  +coeff(87)            *x43*x53
     7  +coeff(88)    *x21*x31    *x54
     8  +coeff(89)    *x22*x32*x43    
      p_m35_0_12  =p_m35_0_12  
     9  +coeff(90)*x11*x21*x31*x42*x51
     1  +coeff(91)    *x21*x31*x44*x51
     2  +coeff(92)    *x22*x31*x42*x52
     3  +coeff(93)*x11*x21    *x41*x53
     4  +coeff(94)    *x21*x31*x42*x53
     5  +coeff(95)    *x21    *x43*x53
     6  +coeff(96)*x12*x23*x31        
     7  +coeff(97)*x12*x23    *x41    
     8  +coeff(98)*x12*x21*x32*x41    
      p_m35_0_12  =p_m35_0_12  
     9  +coeff(99)    *x24*x32*x41*x51
     1  +coeff(100)*x11*x22*x31*x42*x51
c
      return
      end
      function l_m35_0_12  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/ -0.2846769E-01/
      data xmin/
     1 -0.99979E-04,-0.20540E+00,-0.59989E-01,-0.36252E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36252E-01, 0.14986E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.30444782E-01,-0.14145540E+00,-0.62395480E-01,-0.34545960E-01,
     + -0.25502042E-01,-0.85625910E-02,-0.50501222E-02,-0.24302350E-01,
     +  0.13741100E-01,-0.21384291E-01,-0.20079340E-01,-0.30002702E-04,
     + -0.17350083E-01, 0.83058430E-02, 0.15139614E-01, 0.37645340E-02,
     +  0.49780640E-02, 0.73690190E-02,-0.69519490E-02,-0.73802020E-02,
     +  0.12139470E-01, 0.31491671E-02,-0.44298810E-02, 0.31005380E-02,
     +  0.10051120E-01, 0.67807040E-02, 0.20365650E-02,-0.15788944E-02,
     +  0.47656890E-02, 0.59707360E-02, 0.76687640E-02, 0.75435540E-02,
     +  0.26308380E-02,-0.16595300E-02,-0.20046050E-02,-0.12717650E-02,
     +  0.17740600E-02,-0.24769161E-02, 0.26830192E-02,-0.44070680E-02,
     + -0.39639901E-02, 0.26245801E-02, 0.26704310E-02,-0.38395242E-02,
     + -0.11684792E-02, 0.18471493E-02,-0.25348960E-03,-0.41336940E-02,
     +  0.40685190E-03, 0.13392910E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      l_m35_0_12  =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)    *x22            
     4  +coeff( 4)            *x42    
     5  +coeff( 5)        *x31*x41    
     6  +coeff( 6)                *x52
     7  +coeff( 7)        *x32        
     8  +coeff( 8)    *x21        *x51
      l_m35_0_12  =l_m35_0_12  
     9  +coeff( 9)        *x31*x41*x51
     1  +coeff(10)    *x21        *x52
     2  +coeff(11)    *x22        *x52
     3  +coeff(12)*x12            *x51
     4  +coeff(13)                *x51
     5  +coeff(14)    *x22        *x51
     6  +coeff(15)            *x42*x51
     7  +coeff(16)        *x32    *x51
     8  +coeff(17)                *x53
      l_m35_0_12  =l_m35_0_12  
     9  +coeff(18)    *x23        *x51
     1  +coeff(19)        *x31*x41*x52
     2  +coeff(20)            *x42*x52
     3  +coeff(21)    *x21        *x53
     4  +coeff(22)    *x23            
     5  +coeff(23)    *x24            
     6  +coeff(24)    *x21    *x42*x51
     7  +coeff(25)    *x22        *x53
     8  +coeff(26)    *x23*x31*x41*x51
      l_m35_0_12  =l_m35_0_12  
     9  +coeff(27)    *x21*x31*x41    
     1  +coeff(28)        *x32    *x52
     2  +coeff(29)    *x24        *x51
     3  +coeff(30)    *x22*x31*x41*x51
     4  +coeff(31)    *x22    *x42*x51
     5  +coeff(32)    *x23    *x42*x51
     6  +coeff(33)    *x23        *x53
     7  +coeff(34)    *x23    *x44    
     8  +coeff(35)        *x31*x43    
      l_m35_0_12  =l_m35_0_12  
     9  +coeff(36)            *x44    
     1  +coeff(37)    *x21*x31*x41*x51
     2  +coeff(38)                *x54
     3  +coeff(39)    *x23    *x42    
     4  +coeff(40)    *x23        *x52
     5  +coeff(41)    *x21    *x42*x52
     6  +coeff(42)        *x31*x41*x53
     7  +coeff(43)            *x42*x53
     8  +coeff(44)    *x21        *x54
      l_m35_0_12  =l_m35_0_12  
     9  +coeff(45)        *x34*x42    
     1  +coeff(46)    *x23*x32    *x51
     2  +coeff(47)    *x21*x32*x44    
     3  +coeff(48)    *x21*x31*x43*x52
     4  +coeff(49)    *x21*x32        
     5  +coeff(50)    *x21    *x42    
c
      return
      end
      function x_m35_0_13  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(101)
      data ncoeff/100/
      data avdat/ -0.1024518E+00/
      data xmin/
     1 -0.99979E-04,-0.20540E+00,-0.59989E-01,-0.36252E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36252E-01, 0.14986E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.97269760E-01, 0.35912560E+00, 0.47238670E+00,-0.51799503E-03,
     + -0.12616990E+00, 0.53204894E+00,-0.60467012E-02,-0.34143444E-01,
     + -0.43269433E-01,-0.89897260E-01,-0.11681490E-03,-0.11929820E+00,
     + -0.48237050E-04, 0.44560710E-01,-0.94424900E-02,-0.31396770E-01,
     + -0.59186550E-01,-0.10261040E+00, 0.13563400E-01, 0.13382510E-01,
     +  0.18373781E-01, 0.77720650E-01, 0.26569630E-01,-0.43597880E-01,
     + -0.18828872E-01,-0.51522062E-03, 0.26858190E-01, 0.46428374E-01,
     +  0.11338281E-01,-0.51730801E-02, 0.23379843E-01, 0.36993860E-01,
     + -0.12500091E-01,-0.66759360E-01, 0.10642614E-01,-0.13040280E+00,
     +  0.37671062E-02,-0.10275911E-01, 0.12853940E-01, 0.10736150E-01,
     + -0.80774873E-02,-0.91564124E-02, 0.86695070E-02, 0.97438400E-02,
     +  0.49917120E-01,-0.58658370E-01, 0.12956260E-01, 0.98494172E-01,
     + -0.27103951E-01, 0.15913003E-02,-0.13194434E+00,-0.86680750E-01,
     + -0.27324543E-02,-0.37801290E-02, 0.14158660E-02,-0.60238622E-01,
     + -0.10722580E+00, 0.48947160E-01,-0.97564980E-02,-0.36652430E-01,
     + -0.69078100E-02,-0.45462090E-01,-0.44320341E-01,-0.10543190E-01,
     +  0.62133190E-02,-0.15986270E-02,-0.40708862E-01, 0.87570592E-01,
     +  0.11377163E+00,-0.14297644E-01,-0.20391073E-01, 0.84838610E-02,
     + -0.46233691E-01, 0.10346882E-01,-0.32739180E-01,-0.47011312E-01,
     +  0.25352213E-02,-0.59276884E-02, 0.16294622E+00, 0.30735701E-01,
     + -0.84620742E-02,-0.10190250E+00,-0.59420280E-01,-0.10965570E-01,
     +  0.11742960E+00,-0.59393510E-01,-0.82578934E-01, 0.84558870E-01,
     +  0.29010504E+00, 0.15921250E+00,-0.11971110E+00, 0.28424162E-01,
     +  0.74712120E-01, 0.45108474E-01, 0.50047540E-01,-0.16876344E+00,
     +  0.20292030E+00,-0.60075080E-01,-0.84206730E-01,-0.15014713E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x25 = x24*x2
      x26 = x25*x2
      x27 = x26*x2
      x28 = x27*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x35 = x34*x3
      x36 = x35*x3
      x37 = x36*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x45 = x44*x4
      x46 = x45*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_m35_0_13  =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)                *x51
     4  +coeff( 4)*x11                
     5  +coeff( 5)    *x22            
     6  +coeff( 6)    *x21        *x51
     7  +coeff( 7)        *x32        
     8  +coeff( 8)        *x31*x41    
      x_m35_0_13  =x_m35_0_13  
     9  +coeff( 9)            *x42    
     1  +coeff(10)                *x52
     2  +coeff(11)*x11*x21            
     3  +coeff(12)    *x23            
     4  +coeff(13)*x11            *x51
     5  +coeff(14)    *x22        *x51
     6  +coeff(15)    *x21*x32        
     7  +coeff(16)    *x21*x31*x41    
     8  +coeff(17)    *x21    *x42    
      x_m35_0_13  =x_m35_0_13  
     9  +coeff(18)    *x21        *x52
     1  +coeff(19)        *x31*x41*x51
     2  +coeff(20)            *x42*x51
     3  +coeff(21)                *x53
     4  +coeff(22)    *x24            
     5  +coeff(23)    *x23        *x51
     6  +coeff(24)    *x22*x31*x41    
     7  +coeff(25)    *x22    *x42    
     8  +coeff(26)    *x22        *x52
      x_m35_0_13  =x_m35_0_13  
     9  +coeff(27)    *x21        *x53
     1  +coeff(28)    *x25            
     2  +coeff(29)    *x23    *x42    
     3  +coeff(30)    *x22*x32    *x51
     4  +coeff(31)    *x21*x31*x43    
     5  +coeff(32)    *x21    *x44    
     6  +coeff(33)    *x24*x32        
     7  +coeff(34)    *x25*x31*x41    
     8  +coeff(35)    *x25    *x42    
      x_m35_0_13  =x_m35_0_13  
     9  +coeff(36)    *x27    *x42    
     1  +coeff(37)        *x32    *x51
     2  +coeff(38)    *x22*x32        
     3  +coeff(39)        *x31*x43    
     4  +coeff(40)            *x44    
     5  +coeff(41)        *x31*x41*x52
     6  +coeff(42)            *x42*x52
     7  +coeff(43)    *x24        *x51
     8  +coeff(44)    *x22        *x53
      x_m35_0_13  =x_m35_0_13  
     9  +coeff(45)    *x21*x32*x42    
     1  +coeff(46)    *x22    *x44    
     2  +coeff(47)        *x34*x42    
     3  +coeff(48)    *x23*x31*x43    
     4  +coeff(49)    *x23    *x44    
     5  +coeff(50)    *x22*x31*x43*x51
     6  +coeff(51)    *x25    *x42*x51
     7  +coeff(52)    *x23    *x44*x52
     8  +coeff(53)        *x32    *x52
      x_m35_0_13  =x_m35_0_13  
     9  +coeff(54)    *x21*x31*x41*x52
     1  +coeff(55)            *x44*x51
     2  +coeff(56)    *x24*x31*x41    
     3  +coeff(57)    *x24    *x42    
     4  +coeff(58)    *x23    *x42*x51
     5  +coeff(59)    *x21*x32*x42*x51
     6  +coeff(60)    *x27            
     7  +coeff(61)    *x26        *x51
     8  +coeff(62)    *x24*x31*x41*x51
      x_m35_0_13  =x_m35_0_13  
     9  +coeff(63)    *x24    *x42*x51
     1  +coeff(64)    *x23*x32    *x52
     2  +coeff(65)    *x22*x33*x41*x51
     3  +coeff(66)    *x22    *x44*x51
     4  +coeff(67)    *x21*x32*x44    
     5  +coeff(68)    *x24*x31*x43    
     6  +coeff(69)    *x24    *x44    
     7  +coeff(70)    *x24    *x42*x52
     8  +coeff(71)    *x23*x33*x41*x51
      x_m35_0_13  =x_m35_0_13  
     9  +coeff(72)    *x23*x31*x43*x51
     1  +coeff(73)    *x23    *x44*x51
     2  +coeff(74)    *x23    *x42*x53
     3  +coeff(75)    *x22*x33*x43    
     4  +coeff(76)    *x22*x31*x45    
     5  +coeff(77)        *x37*x41    
     6  +coeff(78)        *x34*x44    
     7  +coeff(79)    *x25    *x44    
     8  +coeff(80)    *x24*x32*x42*x51
      x_m35_0_13  =x_m35_0_13  
     9  +coeff(81)    *x24    *x44*x51
     1  +coeff(82)    *x23*x31*x45    
     2  +coeff(83)    *x23    *x46    
     3  +coeff(84)    *x23*x33*x41*x52
     4  +coeff(85)    *x25    *x44*x51
     5  +coeff(86)    *x22*x34*x44    
     6  +coeff(87)    *x28    *x42*x51
     7  +coeff(88)    *x27*x31*x43    
     8  +coeff(89)    *x26*x31*x43*x51
      x_m35_0_13  =x_m35_0_13  
     9  +coeff(90)    *x26    *x44*x51
     1  +coeff(91)    *x25*x34*x42    
     2  +coeff(92)    *x23*x37*x41    
     3  +coeff(93)    *x23    *x46*x52
     4  +coeff(94)    *x26*x34*x42    
     5  +coeff(95)    *x28*x32*x42*x51
     6  +coeff(96)    *x28*x31*x43*x51
     7  +coeff(97)    *x27*x34*x42    
     8  +coeff(98)    *x23*x37*x43    
      x_m35_0_13  =x_m35_0_13  
     9  +coeff(99)    *x27*x34*x44    
     1  +coeff(100)        *x34        
c
      return
      end
      function t_m35_0_13  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(101)
      data ncoeff/100/
      data avdat/ -0.2211773E-01/
      data xmin/
     1 -0.99979E-04,-0.20540E+00,-0.59989E-01,-0.36252E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36252E-01, 0.14986E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.21461021E-01, 0.35432294E-01, 0.10848430E+00,-0.25532694E-01,
     + -0.99078850E-02, 0.10815870E+00,-0.23623934E-01,-0.77098673E-02,
     + -0.15469061E-02,-0.19542090E-01,-0.14024463E-01, 0.12777180E-01,
     +  0.22919440E-02,-0.24536940E-01, 0.53088800E-02,-0.87080430E-02,
     +  0.18731203E-02,-0.15046800E-02, 0.57050430E-03, 0.15518300E-01,
     + -0.88550610E-02, 0.32571121E-02, 0.45727020E-02,-0.19566350E-03,
     + -0.14449480E-02,-0.23947132E-02, 0.75307810E-02,-0.13722951E-02,
     + -0.10251201E-01, 0.69831530E-02, 0.28523740E-03,-0.23869430E-02,
     + -0.98050780E-05,-0.42592692E-02, 0.43904162E-02,-0.58660283E-03,
     + -0.17799631E-02, 0.28372550E-01, 0.89096000E-02,-0.23962290E-02,
     + -0.13525620E-01,-0.11163001E-02,-0.34060883E-02,-0.62463884E-02,
     +  0.17380220E-02, 0.20627374E-02,-0.34204872E-02, 0.14152063E-01,
     +  0.12327713E-01,-0.65479864E-03,-0.36131271E-02, 0.18299540E-02,
     + -0.17438032E-02, 0.99972672E-02,-0.16672830E-01,-0.81372242E-02,
     + -0.28677682E-02,-0.17739400E-02,-0.26976340E-02,-0.10701333E-01,
     + -0.73750503E-02,-0.32096680E-02,-0.49828194E-01,-0.35430114E-01,
     +  0.17061830E-01,-0.27722760E-01,-0.38551073E-01, 0.17843861E-01,
     + -0.57613270E-02,-0.14028080E-01,-0.60010930E-02,-0.23544044E-03,
     +  0.76185920E-02, 0.16283812E-01,-0.16052523E-01, 0.41789990E-02,
     +  0.12710673E-01, 0.19817454E-01, 0.59821573E-02,-0.81611640E-02,
     +  0.55979941E-01, 0.63700470E-01, 0.50881640E-02, 0.22376973E-01,
     +  0.21734813E-01,-0.10369650E-03,-0.37476932E-02, 0.18838760E-03,
     +  0.22151570E-02,-0.41392480E-02,-0.81810980E-02,-0.55104480E-02,
     + -0.48938102E-03, 0.37149800E-03,-0.43221003E-04,-0.12945250E-02,
     +  0.23310190E-02, 0.48836960E-02,-0.29116590E-02, 0.46405150E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x25 = x24*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x45 = x44*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
      x55 = x54*x5
c
c                  function
c
      t_m35_0_13  =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)                *x51
     4  +coeff( 4)    *x22            
     5  +coeff( 5)            *x42    
     6  +coeff( 6)    *x21        *x51
     7  +coeff( 7)                *x52
     8  +coeff( 8)        *x31*x41    
      t_m35_0_13  =t_m35_0_13  
     9  +coeff( 9)        *x32        
     1  +coeff(10)    *x23            
     2  +coeff(11)    *x21    *x42    
     3  +coeff(12)    *x22        *x51
     4  +coeff(13)            *x42*x51
     5  +coeff(14)    *x21        *x52
     6  +coeff(15)                *x53
     7  +coeff(16)    *x21*x31*x41    
     8  +coeff(17)        *x31*x41*x51
      t_m35_0_13  =t_m35_0_13  
     9  +coeff(18)    *x21*x32        
     1  +coeff(19)        *x32    *x51
     2  +coeff(20)    *x24            
     3  +coeff(21)    *x22    *x42    
     4  +coeff(22)            *x44    
     5  +coeff(23)    *x23        *x51
     6  +coeff(24)    *x21    *x42*x51
     7  +coeff(25)    *x22        *x52
     8  +coeff(26)            *x42*x52
      t_m35_0_13  =t_m35_0_13  
     9  +coeff(27)    *x21        *x53
     1  +coeff(28)                *x54
     2  +coeff(29)    *x22*x31*x41    
     3  +coeff(30)        *x31*x43    
     4  +coeff(31)    *x21*x31*x41*x51
     5  +coeff(32)        *x31*x41*x52
     6  +coeff(33)*x11*x22            
     7  +coeff(34)    *x22*x32        
     8  +coeff(35)        *x32*x42    
      t_m35_0_13  =t_m35_0_13  
     9  +coeff(36)        *x32    *x52
     1  +coeff(37)    *x25            
     2  +coeff(38)    *x23    *x42    
     3  +coeff(39)    *x21    *x44    
     4  +coeff(40)    *x24        *x51
     5  +coeff(41)    *x22    *x42*x51
     6  +coeff(42)            *x44*x51
     7  +coeff(43)    *x23        *x52
     8  +coeff(44)    *x21    *x42*x52
      t_m35_0_13  =t_m35_0_13  
     9  +coeff(45)    *x22        *x53
     1  +coeff(46)            *x42*x53
     2  +coeff(47)    *x21        *x54
     3  +coeff(48)    *x23*x31*x41    
     4  +coeff(49)    *x21*x31*x43    
     5  +coeff(50)        *x31*x43*x51
     6  +coeff(51)    *x21*x31*x41*x52
     7  +coeff(52)        *x31*x41*x53
     8  +coeff(53)    *x23*x32        
      t_m35_0_13  =t_m35_0_13  
     9  +coeff(54)    *x21*x32*x42    
     1  +coeff(55)    *x24    *x42    
     2  +coeff(56)    *x22    *x44    
     3  +coeff(57)    *x21    *x44*x51
     4  +coeff(58)    *x22    *x42*x52
     5  +coeff(59)    *x22        *x54
     6  +coeff(60)    *x24*x31*x41    
     7  +coeff(61)    *x22*x31*x43    
     8  +coeff(62)    *x23*x31*x41*x51
      t_m35_0_13  =t_m35_0_13  
     9  +coeff(63)    *x25    *x42    
     1  +coeff(64)    *x23    *x44    
     2  +coeff(65)    *x22    *x44*x51
     3  +coeff(66)    *x25*x31*x41    
     4  +coeff(67)    *x23*x31*x43    
     5  +coeff(68)    *x22*x31*x43*x51
     6  +coeff(69)    *x22*x31*x41*x53
     7  +coeff(70)    *x23*x32*x42    
     8  +coeff(71)    *x21*x32*x44    
      t_m35_0_13  =t_m35_0_13  
     9  +coeff(72)    *x24*x32    *x51
     1  +coeff(73)    *x22*x32*x42*x51
     2  +coeff(74)    *x24    *x44    
     3  +coeff(75)    *x25    *x42*x51
     4  +coeff(76)    *x21    *x44*x53
     5  +coeff(77)    *x24*x31*x43    
     6  +coeff(78)    *x23*x31*x43*x51
     7  +coeff(79)    *x21*x31*x43*x53
     8  +coeff(80)    *x23*x33*x41*x51
      t_m35_0_13  =t_m35_0_13  
     9  +coeff(81)    *x25    *x44    
     1  +coeff(82)    *x25*x31*x43    
     2  +coeff(83)    *x22*x31*x41*x55
     3  +coeff(84)    *x25*x32*x42    
     4  +coeff(85)    *x25    *x44*x51
     5  +coeff(86)*x11                
     6  +coeff(87)    *x22*x31*x41*x51
     7  +coeff(88)    *x21*x33*x41    
     8  +coeff(89)    *x25        *x51
      t_m35_0_13  =t_m35_0_13  
     9  +coeff(90)        *x31*x45    
     1  +coeff(91)    *x21*x31*x43*x51
     2  +coeff(92)        *x32*x44    
     3  +coeff(93)*x11*x23        *x51
     4  +coeff(94)*x11*x21        *x53
     5  +coeff(95)    *x22*x33*x41    
     6  +coeff(96)        *x33*x43    
     7  +coeff(97)    *x21*x33*x41*x51
     8  +coeff(98)    *x21    *x44*x52
      t_m35_0_13  =t_m35_0_13  
     9  +coeff(99)    *x21*x31*x45    
     1  +coeff(100)    *x21*x31*x43*x52
c
      return
      end
      function y_m35_0_13  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(101)
      data ncoeff/100/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.99979E-04,-0.20540E+00,-0.59989E-01,-0.36252E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36252E-01, 0.14986E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.45858480E-01, 0.14567620E+00,-0.49571162E-02,-0.80802110E-02,
     +  0.27405820E-01, 0.81861780E-01,-0.81345870E-04,-0.12091951E-01,
     + -0.15037494E-02, 0.11955181E-04,-0.26575240E-01,-0.70648230E-02,
     + -0.16409862E-01,-0.16764421E-01,-0.93699563E-02,-0.15648760E-01,
     + -0.11749424E-01,-0.24923880E-01,-0.83596054E-02, 0.40182954E-03,
     +  0.75781270E-03,-0.15041670E-01, 0.73683774E-02, 0.23605093E-01,
     +  0.13186760E-01,-0.25045953E-03, 0.16766844E-02, 0.10629751E-03,
     + -0.17359530E-01, 0.93304370E-02, 0.21040763E-01, 0.13291724E-01,
     + -0.57915253E-02,-0.13649304E-01, 0.53521832E-02, 0.68591440E-02,
     +  0.15486980E-02, 0.31796592E-03,-0.18966993E-02, 0.57600470E-02,
     +  0.11816080E-03, 0.15223820E-01, 0.11172080E-01,-0.88795780E-02,
     + -0.86748800E-02, 0.17144550E-03,-0.50672311E-02, 0.42020720E-02,
     +  0.99482690E-02, 0.41043470E-02,-0.12238740E-01,-0.92504192E-02,
     + -0.66203274E-02, 0.82410583E-02, 0.10809770E-01,-0.57437750E-03,
     + -0.20457570E-02,-0.10134490E-02,-0.10365900E-01,-0.30690370E-02,
     +  0.85714794E-02,-0.50306650E-02, 0.98158720E-02, 0.40240241E-02,
     + -0.20061244E-02, 0.64463303E-02,-0.18092760E-02, 0.83462580E-02,
     +  0.22343290E-02,-0.27754481E-02,-0.30170411E-02,-0.21890732E-02,
     + -0.40375664E-02, 0.24161500E-01, 0.12078721E-02,-0.17381810E-02,
     +  0.17880750E-01,-0.10002530E-01,-0.33650572E-02, 0.81326460E-02,
     +  0.11857111E-02,-0.31541062E-02, 0.12296510E-01, 0.98378940E-02,
     + -0.65885963E-02,-0.26685101E-02,-0.19444824E-02, 0.20262560E-02,
     + -0.13558533E-02,-0.57205860E-02, 0.31088062E-02, 0.23244920E-02,
     +  0.55541500E-02,-0.36235593E-02, 0.35148840E-02,-0.18375450E-02,
     + -0.14636300E-02, 0.44668220E-02, 0.18289500E-01, 0.13581720E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      y_m35_0_13  =avdat
     1  +coeff( 1)        *x31        
     2  +coeff( 2)            *x41    
     3  +coeff( 3)    *x21*x31        
     4  +coeff( 4)    *x21    *x41    
     5  +coeff( 5)        *x31    *x51
     6  +coeff( 6)            *x41*x51
     7  +coeff( 7)*x11    *x31        
     8  +coeff( 8)    *x22*x31        
      y_m35_0_13  =y_m35_0_13  
     9  +coeff( 9)        *x33        
     1  +coeff(10)*x11        *x41    
     2  +coeff(11)    *x22    *x41    
     3  +coeff(12)        *x32*x41    
     4  +coeff(13)        *x31*x42    
     5  +coeff(14)            *x43    
     6  +coeff(15)    *x21*x31    *x51
     7  +coeff(16)    *x21    *x41*x51
     8  +coeff(17)        *x31    *x52
      y_m35_0_13  =y_m35_0_13  
     9  +coeff(18)            *x41*x52
     1  +coeff(19)    *x23*x31        
     2  +coeff(20)    *x21*x33        
     3  +coeff(21)*x11*x21    *x41    
     4  +coeff(22)    *x23    *x41    
     5  +coeff(23)    *x21*x32*x41    
     6  +coeff(24)    *x21*x31*x42    
     7  +coeff(25)    *x21    *x43    
     8  +coeff(26)*x11    *x31    *x51
      y_m35_0_13  =y_m35_0_13  
     9  +coeff(27)        *x33    *x51
     1  +coeff(28)*x11        *x41*x51
     2  +coeff(29)    *x22    *x41*x51
     3  +coeff(30)        *x32*x41*x51
     4  +coeff(31)        *x31*x42*x51
     5  +coeff(32)            *x43*x51
     6  +coeff(33)    *x21*x31    *x52
     7  +coeff(34)    *x21    *x41*x52
     8  +coeff(35)        *x31    *x53
      y_m35_0_13  =y_m35_0_13  
     9  +coeff(36)            *x41*x53
     1  +coeff(37)    *x22*x33        
     2  +coeff(38)*x11*x22    *x41    
     3  +coeff(39)    *x24    *x41    
     4  +coeff(40)    *x22*x32*x41    
     5  +coeff(41)        *x34*x41    
     6  +coeff(42)    *x22*x31*x42    
     7  +coeff(43)    *x22    *x43    
     8  +coeff(44)        *x32*x43    
      y_m35_0_13  =y_m35_0_13  
     9  +coeff(45)        *x31*x44    
     1  +coeff(46)*x11*x21    *x41*x51
     2  +coeff(47)    *x23    *x41*x51
     3  +coeff(48)    *x21*x32*x41*x51
     4  +coeff(49)    *x21*x31*x42*x51
     5  +coeff(50)    *x21    *x43*x51
     6  +coeff(51)    *x22    *x41*x52
     7  +coeff(52)        *x31*x42*x52
     8  +coeff(53)            *x43*x52
      y_m35_0_13  =y_m35_0_13  
     9  +coeff(54)    *x21*x31    *x53
     1  +coeff(55)    *x21    *x41*x53
     2  +coeff(56)            *x41*x54
     3  +coeff(57)*x11*x23    *x41    
     4  +coeff(58)    *x23*x32*x41    
     5  +coeff(59)    *x21*x31*x44    
     6  +coeff(60)*x11*x22    *x41*x51
     7  +coeff(61)    *x24    *x41*x51
     8  +coeff(62)    *x22*x31*x42*x51
      y_m35_0_13  =y_m35_0_13  
     9  +coeff(63)    *x22    *x43*x51
     1  +coeff(64)    *x23*x31    *x52
     2  +coeff(65)*x11*x21    *x41*x52
     3  +coeff(66)    *x23    *x41*x52
     4  +coeff(67)    *x22*x31    *x53
     5  +coeff(68)    *x22    *x41*x53
     6  +coeff(69)            *x43*x53
     7  +coeff(70)    *x22*x33*x42    
     8  +coeff(71)        *x33*x44    
      y_m35_0_13  =y_m35_0_13  
     9  +coeff(72)*x11*x23    *x41*x51
     1  +coeff(73)*x11*x21*x31*x42*x51
     2  +coeff(74)    *x23*x31*x42*x51
     3  +coeff(75)    *x21*x33*x42*x51
     4  +coeff(76)*x11*x21    *x43*x51
     5  +coeff(77)    *x23    *x43*x51
     6  +coeff(78)    *x21*x31*x44*x51
     7  +coeff(79)*x11*x22    *x41*x52
     8  +coeff(80)    *x24    *x41*x52
      y_m35_0_13  =y_m35_0_13  
     9  +coeff(81)    *x22*x32*x41*x52
     1  +coeff(82)        *x34*x41*x52
     2  +coeff(83)    *x22*x31*x42*x52
     3  +coeff(84)    *x22    *x43*x52
     4  +coeff(85)    *x23*x31    *x53
     5  +coeff(86)    *x22*x31    *x54
     6  +coeff(87)*x12*x23    *x41    
     7  +coeff(88)*x12*x21*x32*x41    
     8  +coeff(89)*x11*x23*x32*x41    
      y_m35_0_13  =y_m35_0_13  
     9  +coeff(90)*x11*x23*x31*x42    
     1  +coeff(91)*x11*x21*x33*x42    
     2  +coeff(92)    *x23*x33*x42    
     3  +coeff(93)*x11*x21*x32*x43    
     4  +coeff(94)    *x21*x34*x43    
     5  +coeff(95)*x11*x21*x31*x44    
     6  +coeff(96)*x12*x22*x31    *x51
     7  +coeff(97)*x12*x22    *x41*x51
     8  +coeff(98)*x11*x22*x32*x41*x51
      y_m35_0_13  =y_m35_0_13  
     9  +coeff(99)    *x24*x31*x42*x51
     1  +coeff(100)*x11    *x33*x42*x51
c
      return
      end
      function p_m35_0_13  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(101)
      data ncoeff/100/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.99979E-04,-0.20540E+00,-0.59989E-01,-0.36252E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36252E-01, 0.14986E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.56032080E-01,-0.12979604E+00, 0.28954060E-03,-0.43045431E-02,
     +  0.16751473E-01, 0.37502920E-01,-0.24156850E-04,-0.46947100E-03,
     + -0.16492134E-04,-0.27877781E-02,-0.18689650E-02,-0.35025973E-02,
     + -0.42791020E-02,-0.43863850E-02,-0.47788410E-02,-0.97025893E-02,
     +  0.32407504E-02, 0.79918680E-02, 0.18138790E-02, 0.28121531E-02,
     +  0.63530812E-02, 0.29575900E-02,-0.44711241E-02, 0.13458120E-02,
     +  0.23110744E-02,-0.25894502E-02,-0.16655942E-02, 0.15275221E-03,
     + -0.31579530E-02, 0.35154744E-02, 0.91464420E-03,-0.23346930E-02,
     + -0.24105021E-02,-0.20727682E-02, 0.39430670E-02, 0.10242510E-02,
     +  0.12945330E-02, 0.35164301E-03, 0.34124772E-02, 0.12892280E-02,
     + -0.36557484E-03, 0.24381473E-03, 0.20314290E-03, 0.79597340E-02,
     +  0.68131300E-02,-0.16206410E-02,-0.22101020E-03,-0.31914174E-03,
     + -0.14667931E-02,-0.22284460E-03,-0.26358790E-02, 0.23455070E-02,
     +  0.41431561E-03, 0.22604404E-02,-0.16165730E-02, 0.87868830E-03,
     +  0.16646971E-03,-0.93536150E-03, 0.13364810E-02,-0.18635680E-02,
     +  0.94481714E-03, 0.20272054E-02,-0.23978310E-02, 0.20628730E-02,
     + -0.21623280E-02, 0.47726940E-02,-0.23366541E-02,-0.44610300E-02,
     + -0.69118570E-02,-0.22966591E-02, 0.21295691E-02,-0.13608214E-02,
     + -0.48179990E-02, 0.81010541E-03,-0.56031212E-03,-0.56547910E-03,
     + -0.11570700E-02,-0.36854530E-03,-0.31649090E-03,-0.10130610E-02,
     + -0.24034030E-03, 0.13655603E-02,-0.27696690E-03, 0.63116830E-03,
     +  0.14693310E-02, 0.37771263E-02, 0.17871860E-02, 0.15635304E-02,
     +  0.23169370E-02,-0.23910333E-02,-0.31166351E-03,-0.70421870E-03,
     +  0.38910120E-02, 0.28072590E-02,-0.13154462E-02, 0.72262343E-03,
     +  0.49177493E-03,-0.33093942E-03,-0.62166330E-03, 0.55868680E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x( 1)-xmax( 1))*scale( 1)
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      p_m35_0_13  =avdat
     1  +coeff( 1)        *x31        
     2  +coeff( 2)            *x41    
     3  +coeff( 3)    *x21*x31        
     4  +coeff( 4)    *x21    *x41    
     5  +coeff( 5)        *x31    *x51
     6  +coeff( 6)            *x41*x51
     7  +coeff( 7)*x11    *x31        
     8  +coeff( 8)    *x22*x31        
      p_m35_0_13  =p_m35_0_13  
     9  +coeff( 9)*x11        *x41    
     1  +coeff(10)    *x22    *x41    
     2  +coeff(11)        *x32*x41    
     3  +coeff(12)        *x31*x42    
     4  +coeff(13)            *x43    
     5  +coeff(14)    *x21    *x41*x51
     6  +coeff(15)        *x31    *x52
     7  +coeff(16)            *x41*x52
     8  +coeff(17)    *x23*x31        
      p_m35_0_13  =p_m35_0_13  
     9  +coeff(18)    *x23    *x41    
     1  +coeff(19)    *x21*x32*x41    
     2  +coeff(20)    *x21    *x43    
     3  +coeff(21)        *x31*x42*x51
     4  +coeff(22)            *x43*x51
     5  +coeff(23)    *x21    *x41*x52
     6  +coeff(24)        *x31    *x53
     7  +coeff(25)            *x41*x53
     8  +coeff(26)        *x31*x44    
      p_m35_0_13  =p_m35_0_13  
     9  +coeff(27)    *x23*x31    *x51
     1  +coeff(28)    *x21*x33    *x51
     2  +coeff(29)    *x23    *x41*x51
     3  +coeff(30)    *x21*x31*x42*x51
     4  +coeff(31)    *x21    *x43*x51
     5  +coeff(32)    *x22    *x41*x52
     6  +coeff(33)        *x31*x42*x52
     7  +coeff(34)            *x43*x52
     8  +coeff(35)    *x21    *x41*x53
      p_m35_0_13  =p_m35_0_13  
     9  +coeff(36)    *x21*x31*x44    
     1  +coeff(37)    *x22*x32*x41*x51
     2  +coeff(38)        *x34*x41*x51
     3  +coeff(39)    *x22*x31*x42*x51
     4  +coeff(40)    *x23*x31    *x52
     5  +coeff(41)    *x21*x31    *x54
     6  +coeff(42)    *x24*x33        
     7  +coeff(43)    *x23*x33    *x51
     8  +coeff(44)    *x23*x31*x42*x51
      p_m35_0_13  =p_m35_0_13  
     9  +coeff(45)    *x23    *x43*x51
     1  +coeff(46)        *x32*x43*x52
     2  +coeff(47)    *x21*x33    *x53
     3  +coeff(48)        *x33    *x54
     4  +coeff(49)    *x22*x34*x41*x51
     5  +coeff(50)        *x33        
     6  +coeff(51)    *x21*x31    *x51
     7  +coeff(52)    *x21*x31*x42    
     8  +coeff(53)        *x33    *x51
      p_m35_0_13  =p_m35_0_13  
     9  +coeff(54)        *x32*x41*x51
     1  +coeff(55)    *x21*x31    *x52
     2  +coeff(56)    *x22*x32*x41    
     3  +coeff(57)    *x22*x31*x42    
     4  +coeff(58)        *x33*x42    
     5  +coeff(59)    *x22    *x43    
     6  +coeff(60)        *x32*x43    
     7  +coeff(61)    *x21*x32*x41*x51
     8  +coeff(62)    *x21*x31    *x53
      p_m35_0_13  =p_m35_0_13  
     9  +coeff(63)    *x23*x32*x41    
     1  +coeff(64)    *x23*x31*x42    
     2  +coeff(65)    *x23    *x43    
     3  +coeff(66)    *x22    *x43*x51
     4  +coeff(67)        *x31*x44*x51
     5  +coeff(68)    *x21*x31*x44*x51
     6  +coeff(69)    *x23*x31*x44    
     7  +coeff(70)    *x22    *x43*x53
     8  +coeff(71)        *x31*x44*x53
      p_m35_0_13  =p_m35_0_13  
     9  +coeff(72)    *x22*x31    *x51
     1  +coeff(73)    *x22    *x41*x51
     2  +coeff(74)    *x24    *x41    
     3  +coeff(75)*x11*x21*x31    *x51
     4  +coeff(76)*x11*x21    *x41*x51
     5  +coeff(77)    *x22*x31    *x52
     6  +coeff(78)        *x32*x41*x52
     7  +coeff(79)            *x41*x54
     8  +coeff(80)    *x21*x32*x43    
      p_m35_0_13  =p_m35_0_13  
     9  +coeff(81)*x11*x22    *x41*x51
     1  +coeff(82)    *x24    *x41*x51
     2  +coeff(83)*x11*x21    *x41*x52
     3  +coeff(84)    *x23    *x41*x52
     4  +coeff(85)    *x22*x31    *x53
     5  +coeff(86)    *x22    *x41*x53
     6  +coeff(87)            *x43*x53
     7  +coeff(88)    *x22*x31*x44    
     8  +coeff(89)    *x23*x32*x41*x51
      p_m35_0_13  =p_m35_0_13  
     9  +coeff(90)    *x21*x32*x43*x51
     1  +coeff(91)*x11*x22    *x41*x52
     2  +coeff(92)    *x24    *x41*x52
     3  +coeff(93)    *x22*x31*x42*x52
     4  +coeff(94)    *x22    *x43*x52
     5  +coeff(95)        *x31*x44*x52
     6  +coeff(96)*x11*x21*x31    *x53
     7  +coeff(97)*x11*x21    *x41*x53
     8  +coeff(98)*x12*x23*x31        
      p_m35_0_13  =p_m35_0_13  
     9  +coeff(99)*x12*x23    *x41    
     1  +coeff(100)*x12*x21*x32*x41    
c
      return
      end
      function l_m35_0_13  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(51)
      data ncoeff/ 50/
      data avdat/ -0.2681724E-01/
      data xmin/
     1 -0.99979E-04,-0.20540E+00,-0.59989E-01,-0.36252E-01,-0.14999E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99965E-04, 0.19942E+00, 0.59989E-01, 0.36252E-01, 0.14986E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.29241730E-01,-0.17180700E+00,-0.57514742E-01,-0.52889253E-01,
     + -0.41936870E-01,-0.75321540E-01,-0.31078900E-01,-0.80847000E-02,
     + -0.27597364E-01, 0.22864912E-02, 0.19530460E-01, 0.21689280E-01,
     + -0.30352070E-01,-0.65359454E-02, 0.12508413E-01, 0.69807580E-02,
     +  0.18806580E-01, 0.53723030E-02,-0.93748880E-02, 0.12585462E-01,
     +  0.58994102E-02,-0.97474921E-02,-0.10240390E-01, 0.13342940E-01,
     +  0.46502630E-02, 0.29811380E-02, 0.68329772E-02, 0.12954073E-01,
     + -0.49642210E-03, 0.10343702E-02, 0.50840293E-02,-0.47477353E-02,
     + -0.23821180E-02,-0.35695670E-02, 0.47251800E-02, 0.76940003E-02,
     +  0.98671531E-02,-0.67523554E-02,-0.45914160E-02, 0.34700390E-02,
     +  0.37348382E-02,-0.52343462E-02,-0.32418884E-02, 0.85464082E-02,
     + -0.65661421E-02, 0.65401450E-03, 0.11830851E-02,-0.99670900E-03,
     +  0.63468310E-02,-0.20519710E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x2 =1.+(x( 2)-xmax( 2))*scale( 2)
      x3 =1.+(x( 3)-xmax( 3))*scale( 3)
      x4 =1.+(x( 4)-xmax( 4))*scale( 4)
      x5 =1.+(x( 5)-xmax( 5))*scale( 5)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      l_m35_0_13  =avdat
     1  +coeff( 1)                    
     2  +coeff( 2)    *x21            
     3  +coeff( 3)                *x51
     4  +coeff( 4)    *x22            
     5  +coeff( 5)            *x42    
     6  +coeff( 6)    *x21        *x51
     7  +coeff( 7)        *x31*x41    
     8  +coeff( 8)                *x52
      l_m35_0_13  =l_m35_0_13  
     9  +coeff( 9)    *x21        *x52
     1  +coeff(10)    *x22        *x51
     2  +coeff(11)        *x31*x41*x51
     3  +coeff(12)            *x42*x51
     4  +coeff(13)    *x22        *x52
     5  +coeff(14)        *x32        
     6  +coeff(15)    *x23            
     7  +coeff(16)                *x53
     8  +coeff(17)    *x21        *x53
      l_m35_0_13  =l_m35_0_13  
     9  +coeff(18)        *x32    *x51
     1  +coeff(19)    *x24            
     2  +coeff(20)    *x23        *x51
     3  +coeff(21)    *x21    *x42*x51
     4  +coeff(22)        *x31*x41*x52
     5  +coeff(23)            *x42*x52
     6  +coeff(24)    *x22        *x53
     7  +coeff(25)    *x21*x31*x41    
     8  +coeff(26)    *x21    *x42    
      l_m35_0_13  =l_m35_0_13  
     9  +coeff(27)    *x21*x31*x41*x51
     1  +coeff(28)    *x22    *x42*x51
     2  +coeff(29)    *x24    *x42    
     3  +coeff(30)    *x24*x33*x41    
     4  +coeff(31)    *x22*x31*x41    
     5  +coeff(32)        *x31*x43    
     6  +coeff(33)        *x32    *x52
     7  +coeff(34)                *x54
     8  +coeff(35)    *x23    *x42    
      l_m35_0_13  =l_m35_0_13  
     9  +coeff(36)    *x24        *x51
     1  +coeff(37)    *x22*x31*x41*x51
     2  +coeff(38)    *x23        *x52
     3  +coeff(39)    *x21    *x42*x52
     4  +coeff(40)        *x31*x41*x53
     5  +coeff(41)            *x42*x53
     6  +coeff(42)    *x21        *x54
     7  +coeff(43)        *x32*x44    
     8  +coeff(44)    *x23    *x42*x51
      l_m35_0_13  =l_m35_0_13  
     9  +coeff(45)    *x21*x31*x43*x52
     1  +coeff(46)    *x21*x32        
     2  +coeff(47)    *x22*x32        
     3  +coeff(48)        *x33*x41    
     4  +coeff(49)    *x22    *x42    
     5  +coeff(50)            *x44    
c
      return
      end


