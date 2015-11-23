c  MAD spectrometer forward and reverse transfer functions for MAD in its 
c  no quadrupoles 12 degree configuration
c  see line 6117 and following for instructions on the reverse functions 
c  minimum angle configuration. Formulated 12/13/04 -JJL  
c   
c forward function general description:  
c functions take trajectories from target to various planes in the spectrometer.  
c The planes are numbered 0 through 15 (0,2,3,4,5,7,9,15)  
c plane 0: the target plane                                   
c plane 2: entrance to magnet 1                                
c plane 3: middle of magnet 1                                  
c plane 4: exit of magnet 1                                    
c plane 5: entrance to magnet 2                                
c plane 7: middle of magnet 2                                 
c plane 9: exit of magnet 2                          
c plane 11: 1st drift chamber
c plane 15: Calorimeter (an acceptance defining aperture)                                 
c  
c useful endplane information  
c all functions give the locations and directions of a trjectory relative to the central one  
c x = x location in the local (magnet coordinate system) of the central trajectory (+x is down)  
c cx = x-direction cosine of the central trajectory (approx. theta)  
c cz = z-direction cosine of the central trajectory (~1) 
c cy = y-direction cosine of the central trajectory (approx. phi)(zero)   
c l = central trajectory path length in mm
C ep#	     x	          cx	          cz            cy		l
C  2	-0.671598	0.085681	0.996323	0             6999.552246
C  3	112.408371	0.001202	0.999999	0             9003.407227
C  4	-21.088808	-0.084298	0.99644	        0            11308.055664
C  5	71.681625	0.169168	0.985587	0            11628.569336
C  7	243.803665	0.007286	0.999973	0            13339.333008
C  9	48.602173	-0.161692	0.986841	0            15351.593750
C  11	-88.020126	0.009318	0.999956	0            16303.638672
C  15	-50.74712	0.009318	0.999956	0            20303.812500
C  
c  
c for various planes there are 5 functions:  
c  x_mnq_0_j(x,m): gives x position at plane j as a function of trajectory  
c              parameters at plane 0.  
c  t_mnq_0_j(x,m): gives theta at plane j  
c  y_mnq_0_j(x,m): gives y at plane j  
c  p_mnq_0_j(x,m): gives phi at plane j  
c  pl_mnq_0_j(x,m): gives difference in length from the "central" trajectory  
c              between planes 0 and j  
c  
c input: m=5  
c        x is a 5 element array (REAL)  
c        x(1)= x position of trajectory at plane 0  
c        x(2)= theta of trajectory at plane 0  
c        x(3)= y position of trajectory at plane 0 
c        x(4)= phi of trajectory at plane 0  
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
c      x2(1) = x_mnq _0_2(x0,5)  
c      x2(2) = t_mnq _0_2(x0,5)  
c      x2(3) = y_mnq _0_2(x0,5)  
c      x2(4) = p_mnq _0_2(x0,5)  
c      x2(5) = x0(5)  
c  at plane 3 
c      x3(1) = x_mnq _0_3(x0,5)  etc...   
c for the central trajectory you'll get zeros all the way through  
c  (that's why it's the central trajectory!)  
c  
c For the functions to be of real use you'll want to know the  
c parameters of the central trajectory at each plane in the   
c coordinate system centered on the magnet axis. That is tabulated  
c above: (see useful endplane information)  
c The magnets are cylinders 1.2 m in diameter extending from the  
c entrance to the exit  
c z and phi are always zero for the central trajectory  
c +x is down.   
c       
c  pl_mnq _0_j(x,5) gives the difference between the ray in question and  
c             the central ray going from plane 0 to plane j  
c 
      function x_mnq_0_2   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(  9)
      data ncoeff/  8/
      data avdat/ -0.3709274E-02/
      data xmin/
     1 -0.49981E-02,-0.50031E-01,-0.10395E+00,-0.50044E-01,-0.49981E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49938E-02, 0.49989E-01, 0.10395E+00, 0.50044E-01, 0.49994E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.35651389E-02, 0.35137123E+00, 0.15987114E-02, 0.50151343E-02,
     +  0.16099235E-03, 0.18641409E-03, 0.10842620E-03,-0.87868451E-04,
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
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      x_mnq_0_2   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x22            
     4  +coeff(  4)*x11                
     5  +coeff(  5)                *x51
     6  +coeff(  6)            *x42    
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)                *x52
c
      return
      end
      function t_mnq_0_2   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 41)
      data ncoeff/ 40/
      data avdat/ -0.3555324E-03/
      data xmin/
     1 -0.49981E-02,-0.50031E-01,-0.10395E+00,-0.50044E-01,-0.49981E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49938E-02, 0.49989E-01, 0.10395E+00, 0.50044E-01, 0.49994E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.35229107E-03, 0.49956195E-01, 0.72011928E-03, 0.29133732E-03,
     +  0.63068270E-04, 0.43107491E-03, 0.72809734E-03,-0.46880313E-03,
     + -0.14397399E-03,-0.20611264E-03,-0.36057684E-03, 0.25123436E-03,
     +  0.23469745E-03,-0.72942734E-04, 0.92281880E-04, 0.13581192E-03,
     +  0.27344489E-04, 0.23048719E-04,-0.26440701E-04,-0.44002070E-04,
     + -0.69673675E-04,-0.11991502E-03, 0.13388746E-04,-0.11500539E-04,
     +  0.82803153E-05, 0.44215521E-04,-0.42371634E-04, 0.11508417E-04,
     + -0.10603075E-04,-0.19839759E-04,-0.72456883E-05,-0.11500134E-04,
     +  0.17737637E-04,-0.61476421E-05, 0.24110104E-04,-0.60230300E-05,
     + -0.16675563E-05,-0.14984626E-04,-0.13210097E-04,-0.24185907E-04,
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
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_mnq_0_2   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22            
     5  +coeff(  5)        *x32        
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)                *x52
      t_mnq_0_2   =t_mnq_0_2   
     9  +coeff(  9)    *x22        *x51
     1  +coeff( 10)        *x31*x41*x51
     2  +coeff( 11)            *x42*x51
     3  +coeff( 12)                *x53
     4  +coeff( 13)            *x42*x52
     5  +coeff( 14)    *x22    *x42    
     6  +coeff( 15)    *x22        *x52
     7  +coeff( 16)        *x31*x41*x52
     8  +coeff( 17)    *x21        *x51
      t_mnq_0_2   =t_mnq_0_2   
     9  +coeff( 18)    *x21    *x42    
     1  +coeff( 19)        *x32    *x51
     2  +coeff( 20)    *x22*x31*x41    
     3  +coeff( 21)        *x31*x41*x53
     4  +coeff( 22)            *x42*x53
     5  +coeff( 23)    *x21*x31*x41    
     6  +coeff( 24)    *x21        *x52
     7  +coeff( 25)*x11*x21            
     8  +coeff( 26)    *x22    *x42*x51
      t_mnq_0_2   =t_mnq_0_2   
     9  +coeff( 27)    *x22        *x53
     1  +coeff( 28)    *x23            
     2  +coeff( 29)        *x32*x42    
     3  +coeff( 30)        *x31*x43    
     4  +coeff( 31)    *x21*x31*x41*x51
     5  +coeff( 32)    *x21    *x42*x51
     6  +coeff( 33)        *x32    *x52
     7  +coeff( 34)*x11*x21        *x51
     8  +coeff( 35)    *x22*x31*x41*x51
      t_mnq_0_2   =t_mnq_0_2   
     9  +coeff( 36)    *x22*x32        
     1  +coeff( 37)*x11    *x32        
     2  +coeff( 38)        *x32    *x53
     3  +coeff( 39)        *x33*x41*x52
     4  +coeff( 40)    *x22    *x42*x52
c
      return
      end
      function y_mnq_0_2   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(  5)
      data ncoeff/  4/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49981E-02,-0.50031E-01,-0.10395E+00,-0.50044E-01,-0.49981E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49938E-02, 0.49989E-01, 0.10395E+00, 0.50044E-01, 0.49994E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.10393745E+00, 0.35026276E+00, 0.13838187E-02,-0.45106328E-04,
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
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
c          set up monomials   functions
      x21 = x2
      x31 = x3
      x41 = x4
c
c                  function
c
      y_mnq_0_2   =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)    *x21*x31        
c
      return
      end
      function p_mnq_0_2   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 35)
      data ncoeff/ 34/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49981E-02,-0.50031E-01,-0.10395E+00,-0.50044E-01,-0.49981E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49938E-02, 0.49989E-01, 0.10395E+00, 0.50044E-01, 0.49994E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.38603175E-04, 0.50103568E-01,-0.22089550E-03,-0.49269194E-03,
     +  0.58820551E-04, 0.11016111E-03, 0.35503888E-03,-0.21179968E-03,
     +  0.10511072E-03, 0.21861657E-04,-0.36484638E-04, 0.46569046E-04,
     +  0.39105929E-04,-0.58663700E-04, 0.10771229E-03, 0.12456304E-04,
     +  0.14763868E-04, 0.16201198E-04,-0.10448987E-04, 0.36091609E-04,
     + -0.82983988E-05,-0.59082846E-04, 0.15649872E-04, 0.19748417E-04,
     + -0.18045261E-04,-0.25613990E-04,-0.23431592E-04, 0.31155549E-04,
     +  0.40193154E-04, 0.45225538E-05, 0.46494974E-05,-0.24594092E-05,
     + -0.57609823E-05,-0.84917274E-05,
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
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_mnq_0_2   =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21*x31        
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)            *x41*x51
     6  +coeff(  6)    *x21*x31    *x51
     7  +coeff(  7)    *x21    *x41*x51
     8  +coeff(  8)    *x21    *x41*x52
      p_mnq_0_2   =p_mnq_0_2   
     9  +coeff(  9)    *x23    *x41    
     1  +coeff( 10)        *x31    *x51
     2  +coeff( 11)            *x41*x52
     3  +coeff( 12)    *x21*x31*x42    
     4  +coeff( 13)    *x21    *x43    
     5  +coeff( 14)    *x21*x31    *x52
     6  +coeff( 15)    *x21    *x41*x53
     7  +coeff( 16)    *x22    *x41    
     8  +coeff( 17)        *x31*x42    
      p_mnq_0_2   =p_mnq_0_2   
     9  +coeff( 18)            *x43    
     1  +coeff( 19)        *x31    *x52
     2  +coeff( 20)    *x23*x31        
     3  +coeff( 21)*x11        *x41    
     4  +coeff( 22)    *x23    *x41*x51
     5  +coeff( 23)    *x21*x32*x41    
     6  +coeff( 24)            *x41*x53
     7  +coeff( 25)    *x23*x31    *x51
     8  +coeff( 26)    *x21*x31*x42*x51
      p_mnq_0_2   =p_mnq_0_2   
     9  +coeff( 27)    *x21    *x43*x51
     1  +coeff( 28)    *x21*x31    *x53
     2  +coeff( 29)    *x23    *x43*x52
     3  +coeff( 30)    *x22*x31        
     4  +coeff( 31)        *x32*x41    
     5  +coeff( 32)*x11    *x31        
     6  +coeff( 33)    *x22    *x41*x51
     7  +coeff( 34)        *x31*x42*x51
c
      return
      end
      function pl_mnq_0_2  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 27)
      data ncoeff/ 26/
      data avdat/ -0.5540815E-02/
      data xmin/
     1 -0.49981E-02,-0.50031E-01,-0.10395E+00,-0.50044E-01,-0.49981E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49938E-02, 0.49989E-01, 0.10395E+00, 0.50044E-01, 0.49994E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.55536358E-02,-0.30598179E-01, 0.12555889E-04,-0.88921487E-02,
     + -0.87799402E-02,-0.43744073E-03,-0.15381695E-04,-0.42736065E-04,
     + -0.41036732E-04,-0.82473316E-05,-0.82440720E-05, 0.91523916E-05,
     + -0.37454006E-05, 0.48977399E-05, 0.79959427E-05, 0.52563419E-05,
     + -0.13323373E-05,-0.20096384E-05, 0.22679355E-05,-0.40639216E-05,
     +  0.56940053E-05, 0.65673498E-05,-0.44360709E-05, 0.20540551E-05,
     + -0.25472107E-05,-0.27847821E-05,
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
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
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
      pl_mnq_0_2  =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x22    *x42    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)            *x42    
     6  +coeff(  6)*x11                
     7  +coeff(  7)                *x51
     8  +coeff(  8)    *x23            
      pl_mnq_0_2  =pl_mnq_0_2  
     9  +coeff(  9)    *x21    *x42    
     1  +coeff( 10)    *x21        *x51
     2  +coeff( 11)        *x31*x41    
     3  +coeff( 12)                *x52
     4  +coeff( 13)    *x21*x31*x41    
     5  +coeff( 14)        *x31*x41*x51
     6  +coeff( 15)            *x42*x51
     7  +coeff( 16)    *x21        *x52
     8  +coeff( 17)        *x32        
      pl_mnq_0_2  =pl_mnq_0_2  
     9  +coeff( 18)*x11*x21            
     1  +coeff( 19)    *x22        *x51
     2  +coeff( 20)                *x53
     3  +coeff( 21)    *x24            
     4  +coeff( 22)            *x44    
     5  +coeff( 23)            *x42*x52
     6  +coeff( 24)    *x21*x31*x41*x51
     7  +coeff( 25)        *x31*x41*x52
     8  +coeff( 26)    *x21        *x53
c
      return
      end
      function x_mnq_0_3   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.2649106E-01/
      data xmin/
     1 -0.49981E-02,-0.50031E-01,-0.10395E+00,-0.50044E-01,-0.49981E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49997E-02, 0.49989E-01, 0.10395E+00, 0.50044E-01, 0.49994E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.26648853E-01, 0.45137087E+00, 0.30537514E-01, 0.49963314E-02,
     +  0.14859164E-02, 0.33157149E-02, 0.40663872E-03,-0.19754304E-01,
     +  0.50546578E-02,-0.17125515E-02,-0.81136968E-03,-0.19795257E-02,
     +  0.10395232E-01,-0.41257306E-02,-0.27652029E-02, 0.20899664E-03,
     +  0.79546793E-03,-0.43592337E-03,-0.32760913E-02, 0.32723436E-03,
     +  0.11279294E-02,-0.10940555E-02, 0.16684696E-02,-0.44947296E-05,
     +  0.10141134E-02, 0.67825511E-03, 0.16577486E-02,-0.38095331E-03,
     +  0.26168686E-02, 0.14072129E-02, 0.27470249E-02, 0.48491602E-04,
     + -0.16652227E-03,-0.17545861E-03, 0.15433900E-03,-0.71756508E-05,
     +  0.33203381E-04,-0.35730776E-03, 0.77388610E-03, 0.76595519E-03,
     + -0.70496800E-03, 0.72844594E-03,-0.10364681E-03,-0.14626735E-03,
     + -0.16666354E-02,-0.36558154E-03,-0.14883558E-02,-0.85687975E-03,
     + -0.56426070E-03, 0.71631098E-03,-0.23503262E-03,-0.81054401E-04,
     + -0.41776628E-03,-0.37048908E-03, 0.26369864E-04,-0.12757420E-03,
     + -0.80824815E-04,-0.26985058E-05, 0.43218755E-03, 0.22150540E-03,
     +  0.29545627E-03, 0.42920728E-03,-0.89022203E-03,-0.10222216E-02,
     +  0.66112974E-04,-0.18748503E-03,-0.32050473E-04,-0.89425390E-04,
     +  0.63745218E-03,-0.24020643E-03,-0.24240025E-03, 0.21791665E-02,
     +  0.44294886E-03,-0.11820476E-03, 0.22666738E-03,-0.41059079E-03,
     +  0.19347734E-02, 0.89491799E-03, 0.89485443E-03, 0.17714786E-02,
     + -0.11599481E-02, 0.61029493E-03,-0.22258842E-02, 0.22368187E-02,
     + -0.31926783E-03,-0.22731064E-03,-0.93567767E-04, 0.18035083E-03,
     +  0.17104165E-03, 0.42222138E-03,
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
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_mnq_0_3   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22            
     5  +coeff(  5)        *x31*x41    
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)                *x52
      x_mnq_0_3   =x_mnq_0_3   
     9  +coeff(  9)*x11                
     1  +coeff( 10)    *x22        *x51
     2  +coeff( 11)        *x31*x41*x51
     3  +coeff( 12)            *x42*x51
     4  +coeff( 13)                *x53
     5  +coeff( 14)    *x22    *x42    
     6  +coeff( 15)    *x23    *x42    
     7  +coeff( 16)        *x32        
     8  +coeff( 17)    *x23            
      x_mnq_0_3   =x_mnq_0_3   
     9  +coeff( 18)    *x21*x31*x41    
     1  +coeff( 19)    *x22*x31*x41    
     2  +coeff( 20)        *x31*x43    
     3  +coeff( 21)            *x42*x52
     4  +coeff( 22)    *x23*x31*x41    
     5  +coeff( 23)    *x22    *x42*x51
     6  +coeff( 24)    *x21    *x42*x51
     7  +coeff( 25)    *x22        *x52
     8  +coeff( 26)        *x31*x41*x52
      x_mnq_0_3   =x_mnq_0_3   
     9  +coeff( 27)    *x22*x31*x41*x51
     1  +coeff( 28)        *x31*x43*x51
     2  +coeff( 29)    *x22*x31*x43    
     3  +coeff( 30)    *x21*x31*x41*x53
     4  +coeff( 31)    *x23    *x42*x53
     5  +coeff( 32)    *x21*x32        
     6  +coeff( 33)        *x32    *x51
     7  +coeff( 34)    *x21        *x52
     8  +coeff( 35)*x11*x21            
      x_mnq_0_3   =x_mnq_0_3   
     9  +coeff( 36)        *x33*x41    
     1  +coeff( 37)        *x32*x42    
     2  +coeff( 38)    *x21        *x53
     3  +coeff( 39)    *x21*x31*x43    
     4  +coeff( 40)    *x23        *x52
     5  +coeff( 41)    *x21*x31*x43*x51
     6  +coeff( 42)    *x22*x31*x41*x52
     7  +coeff( 43)    *x23        *x53
     8  +coeff( 44)*x11    *x32*x42    
      x_mnq_0_3   =x_mnq_0_3   
     9  +coeff( 45)    *x22*x31*x43*x51
     1  +coeff( 46)    *x23*x31*x41*x52
     2  +coeff( 47)    *x23    *x42*x52
     3  +coeff( 48)        *x32*x42*x53
     4  +coeff( 49)*x11*x22*x32    *x51
     5  +coeff( 50)*x11*x21*x33*x43*x51
     6  +coeff( 51)*x11*x23*x31*x41*x53
     7  +coeff( 52)*x11            *x51
     8  +coeff( 53)    *x23        *x51
      x_mnq_0_3   =x_mnq_0_3   
     9  +coeff( 54)    *x23*x32        
     1  +coeff( 55)*x11    *x31*x41    
     2  +coeff( 56)*x11*x21        *x51
     3  +coeff( 57)*x11*x21*x32        
     4  +coeff( 58)*x11*x21    *x42    
     5  +coeff( 59)    *x22*x32*x42    
     6  +coeff( 60)*x11*x22        *x51
     7  +coeff( 61)    *x23*x32    *x51
     8  +coeff( 62)    *x23*x31*x41*x51
      x_mnq_0_3   =x_mnq_0_3   
     9  +coeff( 63)        *x33*x41*x52
     1  +coeff( 64)    *x22    *x42*x52
     2  +coeff( 65)        *x32*x42*x52
     3  +coeff( 66)*x11*x22    *x42    
     4  +coeff( 67)    *x23*x32*x42    
     5  +coeff( 68)*x11    *x32    *x52
     6  +coeff( 69)    *x22    *x42*x53
     7  +coeff( 70)*x11*x21*x33*x41    
     8  +coeff( 71)*x11*x23    *x42    
      x_mnq_0_3   =x_mnq_0_3   
     9  +coeff( 72)    *x22*x32*x42*x52
     1  +coeff( 73)        *x33*x43*x52
     2  +coeff( 74)*x11*x22        *x53
     3  +coeff( 75)*x11*x23    *x42*x51
     4  +coeff( 76)    *x23*x33*x41*x52
     5  +coeff( 77)    *x23*x32*x42*x52
     6  +coeff( 78)*x11*x22*x32*x42*x51
     7  +coeff( 79)*x11*x22*x32    *x53
     8  +coeff( 80)*x11*x23*x32*x42*x51
      x_mnq_0_3   =x_mnq_0_3   
     9  +coeff( 81)*x11*x21*x32*x42*x53
     1  +coeff( 82)*x11*x22*x33*x41*x53
     2  +coeff( 83)*x11*x23*x31*x43*x53
     3  +coeff( 84)*x11*x21*x33*x43*x53
     4  +coeff( 85)    *x22*x32        
     5  +coeff( 86)    *x21*x31*x41*x51
     6  +coeff( 87)        *x32    *x52
     7  +coeff( 88)    *x22*x32    *x51
     8  +coeff( 89)        *x33*x41*x51
      x_mnq_0_3   =x_mnq_0_3   
     9  +coeff( 90)        *x32*x42*x51
c
      return
      end
      function t_mnq_0_3   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.1169026E-01/
      data xmin/
     1 -0.49981E-02,-0.50031E-01,-0.10395E+00,-0.50044E-01,-0.49981E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49997E-02, 0.49989E-01, 0.10395E+00, 0.50044E-01, 0.49994E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.12177239E-01, 0.48906814E-01, 0.42104021E-01, 0.24947078E-04,
     +  0.13214283E-03, 0.33373432E-03, 0.99943543E-03,-0.27469482E-01,
     +  0.10337107E-03, 0.96826896E-03,-0.13855133E-03,-0.50857285E-03,
     +  0.64029469E-03, 0.36047344E-03,-0.35623985E-03, 0.14800804E-01,
     + -0.39790844E-03,-0.47231670E-02,-0.67394561E-04,-0.52494667E-02,
     +  0.11181395E-04,-0.16452405E-02, 0.16595716E-03,-0.13400642E-02,
     + -0.15086900E-03,-0.19360395E-02,-0.49962499E-02, 0.17554957E-02,
     +  0.16356743E-03, 0.23902166E-02,-0.10006307E-02, 0.70189708E-03,
     + -0.29651387E-03,-0.37630642E-03,-0.19617651E-03, 0.47290386E-02,
     +  0.32218965E-03, 0.14377094E-02, 0.30645810E-02,-0.47686789E-03,
     +  0.93186783E-04, 0.39638815E-03, 0.56955492E-03, 0.24634252E-04,
     +  0.18915996E-02, 0.17371729E-02, 0.89706416E-03,-0.55562885E-03,
     +  0.66067732E-04,-0.27811862E-03,-0.16527253E-03,-0.11368766E-02,
     + -0.28949461E-03,-0.10784114E-03,-0.22228622E-03,-0.15249148E-03,
     +  0.82923181E-03, 0.29258701E-03,-0.49795894E-04, 0.51949418E-03,
     + -0.33234214E-03,-0.14268469E-03,-0.17444903E-03, 0.18009061E-02,
     +  0.12020746E-03, 0.15651491E-02,-0.73427602E-03, 0.80267544E-03,
     +  0.90264081E-03, 0.12646658E-03, 0.97900062E-04,-0.96539740E-03,
     + -0.21212141E-03,-0.32118033E-03, 0.46444937E-04,-0.60375361E-04,
     +  0.64889499E-03, 0.14415946E-03, 0.15100269E-03, 0.26207251E-03,
     +  0.33177415E-03, 0.35665624E-03, 0.39033868E-03,-0.62670457E-04,
     +  0.36647829E-03,-0.31656251E-04,-0.20106739E-03,-0.41180178E-04,
     + -0.51389798E-04, 0.12959992E-04,
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
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_mnq_0_3   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)        *x32        
     5  +coeff(  5)        *x31*x41    
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)                *x52
      t_mnq_0_3   =t_mnq_0_3   
     9  +coeff(  9)*x11                
     1  +coeff( 10)    *x23            
     2  +coeff( 11)    *x21*x32        
     3  +coeff( 12)    *x21*x31*x41    
     4  +coeff( 13)    *x21    *x42    
     5  +coeff( 14)        *x31*x41*x51
     6  +coeff( 15)            *x42*x51
     7  +coeff( 16)                *x53
     8  +coeff( 17)    *x22*x32        
      t_mnq_0_3   =t_mnq_0_3   
     9  +coeff( 18)    *x22*x31*x41    
     1  +coeff( 19)        *x33*x41    
     2  +coeff( 20)    *x22    *x42    
     3  +coeff( 21)        *x31*x43    
     4  +coeff( 22)    *x21    *x42*x51
     5  +coeff( 23)    *x22        *x52
     6  +coeff( 24)    *x21        *x53
     7  +coeff( 25)*x11    *x32        
     8  +coeff( 26)    *x23*x31*x41    
      t_mnq_0_3   =t_mnq_0_3   
     9  +coeff( 27)    *x23    *x42    
     1  +coeff( 28)    *x21*x31*x43    
     2  +coeff( 29)    *x22*x32    *x51
     3  +coeff( 30)    *x22    *x42*x51
     4  +coeff( 31)        *x31*x43*x51
     5  +coeff( 32)    *x23        *x52
     6  +coeff( 33)    *x21*x31*x41*x52
     7  +coeff( 34)    *x21    *x42*x52
     8  +coeff( 35)*x11*x21*x31*x41    
      t_mnq_0_3   =t_mnq_0_3   
     9  +coeff( 36)    *x22*x31*x43    
     1  +coeff( 37)*x11*x22        *x51
     2  +coeff( 38)    *x23*x31*x41*x51
     3  +coeff( 39)    *x23    *x42*x51
     4  +coeff( 40)    *x21*x31*x43*x51
     5  +coeff( 41)    *x22*x32    *x52
     6  +coeff( 42)    *x22    *x42*x52
     7  +coeff( 43)        *x32*x42*x52
     8  +coeff( 44)*x11            *x53
      t_mnq_0_3   =t_mnq_0_3   
     9  +coeff( 45)    *x21*x31*x41*x53
     1  +coeff( 46)    *x21    *x42*x53
     2  +coeff( 47)    *x22            
     3  +coeff( 48)    *x22        *x51
     4  +coeff( 49)    *x21        *x52
     5  +coeff( 50)        *x32*x42    
     6  +coeff( 51)*x11            *x51
     7  +coeff( 52)    *x23        *x51
     8  +coeff( 53)    *x21*x32    *x51
      t_mnq_0_3   =t_mnq_0_3   
     9  +coeff( 54)*x11*x22            
     1  +coeff( 55)    *x23*x32        
     2  +coeff( 56)*x11        *x42    
     3  +coeff( 57)    *x22*x31*x41*x51
     4  +coeff( 58)        *x33*x41*x51
     5  +coeff( 59)        *x32*x42*x51
     6  +coeff( 60)    *x22        *x53
     7  +coeff( 61)        *x32    *x53
     8  +coeff( 62)*x11*x21*x32        
      t_mnq_0_3   =t_mnq_0_3   
     9  +coeff( 63)*x11*x21    *x42    
     1  +coeff( 64)    *x22*x32*x42    
     2  +coeff( 65)*x11        *x42*x51
     3  +coeff( 66)    *x22*x31*x41*x52
     4  +coeff( 67)        *x33*x41*x52
     5  +coeff( 68)        *x31*x43*x52
     6  +coeff( 69)    *x23        *x53
     7  +coeff( 70)        *x32    *x51
     8  +coeff( 71)*x11*x21            
      t_mnq_0_3   =t_mnq_0_3   
     9  +coeff( 72)    *x21*x31*x41*x51
     1  +coeff( 73)        *x32    *x52
     2  +coeff( 74)        *x31*x41*x52
     3  +coeff( 75)*x11    *x31*x41    
     4  +coeff( 76)    *x21*x33*x41    
     5  +coeff( 77)    *x21*x32*x42    
     6  +coeff( 78)    *x21*x32    *x52
     7  +coeff( 79)        *x31*x41*x53
     8  +coeff( 80)            *x42*x53
      t_mnq_0_3   =t_mnq_0_3   
     9  +coeff( 81)    *x22*x33*x41    
     1  +coeff( 82)    *x23*x32    *x51
     2  +coeff( 83)    *x21*x32*x42*x51
     3  +coeff( 84)*x11*x21        *x52
     4  +coeff( 85)    *x21*x32    *x53
     5  +coeff( 86)*x11*x21        *x51
     6  +coeff( 87)        *x33*x43    
     7  +coeff( 88)*x11    *x31*x41*x51
     8  +coeff( 89)            *x42*x52
      t_mnq_0_3   =t_mnq_0_3   
     9  +coeff( 90)*x11            *x52
c
      return
      end
      function y_mnq_0_3   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 77)
      data ncoeff/ 76/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49981E-02,-0.50031E-01,-0.10395E+00,-0.50044E-01,-0.49981E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49997E-02, 0.49989E-01, 0.10395E+00, 0.50044E-01, 0.49994E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.10333942E+00, 0.44841018E+00,-0.17961873E-02,-0.46305386E-02,
     +  0.93766960E-03,-0.26960480E-02, 0.34186172E-02,-0.42139706E-02,
     +  0.69058855E-03,-0.21102170E-02, 0.65041339E-03, 0.26437207E-03,
     + -0.68436493E-03, 0.83496643E-03,-0.35791140E-03,-0.11151793E-02,
     +  0.28845418E-03, 0.36765437E-02, 0.51937071E-02,-0.12416182E-03,
     +  0.13436865E-02, 0.12292276E-02,-0.74640167E-03, 0.48377173E-04,
     +  0.27773252E-02,-0.22335142E-04,-0.40255938E-03, 0.40403958E-02,
     + -0.20248634E-03,-0.10641183E-03,-0.72991490E-04, 0.30486551E-03,
     +  0.24719053E-03, 0.53869502E-03,-0.36632468E-03, 0.33196621E-03,
     +  0.90297451E-03, 0.26774721E-03,-0.16279948E-03, 0.89304624E-04,
     +  0.96248339E-04,-0.27349667E-03,-0.15350962E-02,-0.93721010E-03,
     +  0.53913740E-03,-0.21228872E-02,-0.19279612E-02, 0.36832897E-03,
     + -0.14057124E-03,-0.93192493E-05, 0.32514404E-03,-0.75350149E-03,
     + -0.74088509E-03, 0.62824826E-03,-0.37463786E-04, 0.17546128E-04,
     +  0.38582832E-04, 0.14235823E-02,-0.17835850E-03,-0.14150853E-03,
     + -0.63109910E-04, 0.76517079E-03,-0.85054067E-04, 0.49530685E-04,
     +  0.51653387E-04,-0.32460186E-03,-0.11663100E-02, 0.34560275E-04,
     +  0.71226124E-04, 0.21891325E-03, 0.14861702E-03,-0.98087701E-04,
     +  0.85052970E-03, 0.27971450E-03,-0.24726504E-03, 0.76059334E-03,
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
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      y_mnq_0_3   =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21*x31        
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)            *x41*x51
     6  +coeff(  6)    *x22    *x41    
     7  +coeff(  7)    *x21    *x41*x51
     8  +coeff(  8)    *x23    *x41    
      y_mnq_0_3   =y_mnq_0_3   
     9  +coeff(  9)    *x21    *x43    
     1  +coeff( 10)    *x21    *x41*x52
     2  +coeff( 11)    *x23*x31    *x51
     3  +coeff( 12)        *x31    *x51
     4  +coeff( 13)    *x22*x31        
     5  +coeff( 14)    *x21*x31    *x51
     6  +coeff( 15)            *x41*x52
     7  +coeff( 16)    *x23*x31        
     8  +coeff( 17)    *x21*x31*x42    
      y_mnq_0_3   =y_mnq_0_3   
     9  +coeff( 18)    *x22    *x43    
     1  +coeff( 19)    *x23    *x43    
     2  +coeff( 20)    *x22*x33*x42    
     3  +coeff( 21)    *x23    *x41*x53
     4  +coeff( 22)    *x22    *x41*x51
     5  +coeff( 23)    *x21*x31    *x52
     6  +coeff( 24)    *x22*x33        
     7  +coeff( 25)    *x22*x31*x42    
     8  +coeff( 26)    *x21*x31*x42*x51
      y_mnq_0_3   =y_mnq_0_3   
     9  +coeff( 27)    *x21    *x43*x51
     1  +coeff( 28)    *x23*x31*x42    
     2  +coeff( 29)            *x43    
     3  +coeff( 30)        *x31    *x52
     4  +coeff( 31)*x11        *x41    
     5  +coeff( 32)    *x21*x32*x41    
     6  +coeff( 33)    *x22*x31    *x51
     7  +coeff( 34)    *x22*x32*x41    
     8  +coeff( 35)    *x22    *x41*x52
      y_mnq_0_3   =y_mnq_0_3   
     9  +coeff( 36)    *x21*x31    *x53
     1  +coeff( 37)    *x23*x32*x41    
     2  +coeff( 38)    *x21    *x43*x52
     3  +coeff( 39)        *x31*x42    
     4  +coeff( 40)            *x43*x51
     5  +coeff( 41)            *x41*x53
     6  +coeff( 42)    *x22*x32*x41*x51
     7  +coeff( 43)    *x22    *x43*x51
     8  +coeff( 44)    *x23    *x41*x52
      y_mnq_0_3   =y_mnq_0_3   
     9  +coeff( 45)    *x21*x31*x42*x52
     1  +coeff( 46)    *x23*x31*x42*x51
     2  +coeff( 47)    *x23    *x43*x51
     3  +coeff( 48)    *x22*x33*x42*x51
     4  +coeff( 49)    *x22*x33    *x53
     5  +coeff( 50)    *x22*x32*x41*x53
     6  +coeff( 51)    *x22    *x43*x53
     7  +coeff( 52)    *x23*x32*x41*x53
     8  +coeff( 53)    *x23    *x43*x53
      y_mnq_0_3   =y_mnq_0_3   
     9  +coeff( 54)    *x21*x32*x43*x53
     1  +coeff( 55)*x11    *x31        
     2  +coeff( 56)*x11*x21    *x41    
     3  +coeff( 57)*x11        *x41*x51
     4  +coeff( 58)    *x23    *x41*x51
     5  +coeff( 59)    *x21*x32*x41*x51
     6  +coeff( 60)    *x22*x31    *x52
     7  +coeff( 61)        *x32*x41*x52
     8  +coeff( 62)    *x21    *x41*x53
      y_mnq_0_3   =y_mnq_0_3   
     9  +coeff( 63)*x11*x22    *x41    
     1  +coeff( 64)*x11    *x31*x42    
     2  +coeff( 65)*x11        *x43    
     3  +coeff( 66)    *x21*x32*x43    
     4  +coeff( 67)    *x22*x31*x42*x51
     5  +coeff( 68)        *x33*x42*x51
     6  +coeff( 69)        *x32*x43*x51
     7  +coeff( 70)    *x22*x31    *x53
     8  +coeff( 71)        *x31*x42*x53
      y_mnq_0_3   =y_mnq_0_3   
     9  +coeff( 72)*x11*x23    *x41    
     1  +coeff( 73)    *x23    *x43*x52
     2  +coeff( 74)*x11*x23    *x43*x51
     3  +coeff( 75)*x11*x21*x32*x43*x51
     4  +coeff( 76)    *x22*x32*x43*x53
c
      return
      end
      function p_mnq_0_3   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49981E-02,-0.50031E-01,-0.10395E+00,-0.50044E-01,-0.49981E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49997E-02, 0.49989E-01, 0.10395E+00, 0.50044E-01, 0.49994E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.25515541E-03, 0.49199332E-01,-0.47371152E-03,-0.11916212E-02,
     +  0.62061932E-04, 0.35211595E-03,-0.88525383E-03,-0.39190999E-02,
     + -0.34836448E-04,-0.40999800E-03,-0.53739670E-03, 0.14964747E-03,
     +  0.72175398E-03,-0.12192049E-02,-0.49489363E-05,-0.56376457E-02,
     +  0.36256583E-03,-0.20806587E-03,-0.45061635E-04, 0.13565775E-02,
     + -0.25337904E-05,-0.58078440E-03,-0.16465760E-03, 0.46647140E-02,
     +  0.67964504E-02, 0.18073347E-02, 0.36911402E-03, 0.32758512E-03,
     + -0.15882862E-03,-0.36195721E-03,-0.25276601E-03, 0.58718179E-02,
     +  0.88500902E-02,-0.67338464E-03,-0.15750020E-03,-0.20722654E-02,
     + -0.58359074E-05,-0.29825598E-02, 0.16919209E-02, 0.85845758E-03,
     +  0.93998889E-04,-0.20146408E-03, 0.67977307E-04, 0.12751516E-02,
     + -0.45980756E-04, 0.12268611E-02, 0.24878524E-04,-0.20611091E-02,
     +  0.44005096E-03,-0.10020818E-02,-0.54318475E-03,-0.34964527E-04,
     + -0.24683645E-02,-0.56628400E-03, 0.32436437E-03, 0.29505984E-03,
     +  0.33257171E-03, 0.18053380E-02, 0.13185672E-02, 0.12566478E-03,
     +  0.38799347E-03,-0.37503825E-03,-0.15264235E-02,-0.16474628E-02,
     + -0.18037249E-02, 0.79040835E-03,-0.64242567E-03, 0.81475812E-03,
     + -0.85382781E-03,-0.52527688E-03,-0.75106177E-03,-0.98193344E-03,
     +  0.39769761E-05, 0.49748016E-03, 0.17408263E-03,-0.21054282E-04,
     +  0.18247785E-04, 0.62807789E-03, 0.22442164E-03, 0.67468442E-04,
     + -0.66089247E-04,-0.55943336E-03, 0.98540593E-04,-0.55605290E-03,
     +  0.10786574E-03,-0.24129407E-03,-0.82542858E-04, 0.14522111E-03,
     +  0.56572870E-04,-0.14545914E-03,
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
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_mnq_0_3   =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21*x31        
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)        *x31    *x51
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)    *x22*x31        
     8  +coeff(  8)    *x22    *x41    
      p_mnq_0_3   =p_mnq_0_3   
     9  +coeff(  9)        *x32*x41    
     1  +coeff( 10)        *x31*x42    
     2  +coeff( 11)            *x43    
     3  +coeff( 12)    *x21*x31    *x51
     4  +coeff( 13)    *x21    *x41*x51
     5  +coeff( 14)    *x23*x31        
     6  +coeff( 15)    *x21*x33        
     7  +coeff( 16)    *x23    *x41    
     8  +coeff( 17)    *x21*x32*x41    
      p_mnq_0_3   =p_mnq_0_3   
     9  +coeff( 18)    *x21*x31*x42    
     1  +coeff( 19)    *x21    *x43    
     2  +coeff( 20)    *x22    *x41*x51
     3  +coeff( 21)    *x21*x31    *x52
     4  +coeff( 22)    *x21    *x41*x52
     5  +coeff( 23)            *x41*x53
     6  +coeff( 24)    *x22*x31*x42    
     7  +coeff( 25)    *x22    *x43    
     8  +coeff( 26)    *x23    *x41*x51
      p_mnq_0_3   =p_mnq_0_3   
     9  +coeff( 27)    *x21*x31*x42*x51
     1  +coeff( 28)    *x21    *x43*x51
     2  +coeff( 29)    *x22*x31    *x52
     3  +coeff( 30)    *x22    *x41*x52
     4  +coeff( 31)*x11*x22    *x41    
     5  +coeff( 32)    *x23*x31*x42    
     6  +coeff( 33)    *x23    *x43    
     7  +coeff( 34)    *x21*x32*x43    
     8  +coeff( 35)    *x22*x33    *x51
      p_mnq_0_3   =p_mnq_0_3   
     9  +coeff( 36)    *x22    *x43*x51
     1  +coeff( 37)    *x21    *x43*x52
     2  +coeff( 38)    *x23    *x43*x51
     3  +coeff( 39)    *x23    *x41*x53
     4  +coeff( 40)    *x22*x33*x42*x51
     5  +coeff( 41)            *x41*x52
     6  +coeff( 42)        *x32*x41*x51
     7  +coeff( 43)            *x43*x51
     8  +coeff( 44)    *x22*x32*x41    
      p_mnq_0_3   =p_mnq_0_3   
     9  +coeff( 45)        *x32*x43    
     1  +coeff( 46)    *x23*x32*x41    
     2  +coeff( 47)*x11    *x31*x42    
     3  +coeff( 48)    *x22*x31*x42*x51
     4  +coeff( 49)        *x32*x43*x51
     5  +coeff( 50)    *x23    *x41*x52
     6  +coeff( 51)    *x21*x31*x42*x52
     7  +coeff( 52)*x11*x23    *x41    
     8  +coeff( 53)    *x23*x31*x42*x51
      p_mnq_0_3   =p_mnq_0_3   
     9  +coeff( 54)    *x22*x32*x41*x52
     1  +coeff( 55)    *x23*x31    *x53
     2  +coeff( 56)    *x21*x31*x42*x53
     3  +coeff( 57)*x11*x22    *x43    
     4  +coeff( 58)    *x23*x31*x42*x52
     5  +coeff( 59)    *x23    *x43*x52
     6  +coeff( 60)        *x33*x42*x53
     7  +coeff( 61)    *x22    *x43*x53
     8  +coeff( 62)*x11*x22*x31*x42*x51
      p_mnq_0_3   =p_mnq_0_3   
     9  +coeff( 63)    *x23*x32*x41*x53
     1  +coeff( 64)    *x23*x31*x42*x53
     2  +coeff( 65)    *x23    *x43*x53
     3  +coeff( 66)    *x21*x32*x43*x53
     4  +coeff( 67)*x11*x21*x33*x42*x51
     5  +coeff( 68)*x11*x22*x31*x42*x52
     6  +coeff( 69)*x11*x22*x32*x43*x51
     7  +coeff( 70)*x11*x23    *x43*x52
     8  +coeff( 71)*x11*x22*x31*x42*x53
      p_mnq_0_3   =p_mnq_0_3   
     9  +coeff( 72)*x11*x21*x32*x43*x53
     1  +coeff( 73)*x11    *x31        
     2  +coeff( 74)    *x22*x31    *x51
     3  +coeff( 75)        *x31*x42*x51
     4  +coeff( 76)*x11*x21*x31        
     5  +coeff( 77)*x11    *x31    *x51
     6  +coeff( 78)    *x23*x31    *x51
     7  +coeff( 79)    *x21*x33*x42    
     8  +coeff( 80)*x11*x21*x31    *x51
      p_mnq_0_3   =p_mnq_0_3   
     9  +coeff( 81)*x11    *x31    *x52
     1  +coeff( 82)    *x23*x31    *x52
     2  +coeff( 83)*x11*x21    *x43    
     3  +coeff( 84)    *x22*x32*x43    
     4  +coeff( 85)*x11*x22*x31    *x51
     5  +coeff( 86)    *x22*x31*x42*x52
     6  +coeff( 87)*x11*x22*x33        
     7  +coeff( 88)*x11*x21    *x43*x51
     8  +coeff( 89)*x11    *x32*x41*x52
      p_mnq_0_3   =p_mnq_0_3   
     9  +coeff( 90)*x11*x23*x32*x41    
c
      return
      end
      function pl_mnq_0_3  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 41)
      data ncoeff/ 40/
      data avdat/ -0.5441641E-02/
      data xmin/
     1 -0.49981E-02,-0.50031E-01,-0.10395E+00,-0.50044E-01,-0.49981E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49997E-02, 0.49989E-01, 0.10395E+00, 0.50044E-01, 0.49994E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.54544699E-02,-0.36242612E-01, 0.35619861E-03,-0.11556229E-01,
     + -0.11367097E-01,-0.11528927E-02,-0.14987387E-02,-0.43795284E-03,
     +  0.71774167E-03, 0.36294447E-03,-0.25789623E-03, 0.10185998E-03,
     + -0.51960908E-03,-0.36353851E-03,-0.78090641E-04, 0.22087894E-03,
     + -0.16447731E-03, 0.62394206E-03,-0.38275317E-04, 0.28158500E-03,
     + -0.13156201E-04, 0.42077700E-05, 0.39627819E-04, 0.25142050E-04,
     +  0.18674780E-03, 0.16609923E-04, 0.51696265E-04, 0.25569284E-03,
     + -0.20568354E-04, 0.11375590E-03,-0.26196335E-03,-0.35389551E-03,
     + -0.37328416E-03,-0.49005792E-03,-0.52957794E-05,-0.11637693E-04,
     +  0.13794614E-04,-0.50097773E-04, 0.39743307E-04, 0.83520761E-04,
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
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
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
      x54 = x53*x5
c
c                  function
c
      pl_mnq_0_3  =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x22    *x42    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)            *x42    
     6  +coeff(  6)                *x51
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)*x11                
      pl_mnq_0_3  =pl_mnq_0_3  
     9  +coeff(  9)    *x21        *x52
     1  +coeff( 10)                *x53
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x21    *x42    
     4  +coeff( 13)    *x21        *x53
     5  +coeff( 14)                *x54
     6  +coeff( 15)        *x31*x41    
     7  +coeff( 16)                *x52
     8  +coeff( 17)    *x21    *x42*x51
      pl_mnq_0_3  =pl_mnq_0_3  
     9  +coeff( 18)    *x23    *x42    
     1  +coeff( 19)    *x21    *x44    
     2  +coeff( 20)    *x21        *x54
     3  +coeff( 21)        *x32        
     4  +coeff( 22)    *x21*x31*x41    
     5  +coeff( 23)    *x22        *x51
     6  +coeff( 24)        *x31*x41*x51
     7  +coeff( 25)    *x22*x31*x41    
     8  +coeff( 26)            *x44    
      pl_mnq_0_3  =pl_mnq_0_3  
     9  +coeff( 27)    *x23        *x51
     1  +coeff( 28)    *x23*x31*x41    
     2  +coeff( 29)    *x21*x31*x43    
     3  +coeff( 30)    *x24    *x42    
     4  +coeff( 31)    *x22*x31*x43    
     5  +coeff( 32)    *x22    *x44    
     6  +coeff( 33)    *x23*x31*x43    
     7  +coeff( 34)    *x23    *x44    
     8  +coeff( 35)*x11*x21            
      pl_mnq_0_3  =pl_mnq_0_3  
     9  +coeff( 36)    *x21*x32        
     1  +coeff( 37)            *x42*x51
     2  +coeff( 38)    *x24            
     3  +coeff( 39)    *x21    *x42*x52
     4  +coeff( 40)    *x21    *x44*x51
c
      return
      end
      function x_mnq_0_4   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.1682478E-02/
      data xmin/
     1 -0.49978E-02,-0.50007E-01,-0.10395E+00,-0.50032E-01,-0.49981E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49997E-02, 0.49989E-01, 0.10395E+00, 0.50032E-01, 0.49998E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.41791131E-02, 0.56317019E+00, 0.19598725E+00, 0.51008100E-02,
     + -0.22109134E-03,-0.52489853E-03,-0.32811619E-02,-0.25418864E-02,
     + -0.13034700E+00, 0.48841666E-02, 0.79882229E-02, 0.33701528E-03,
     +  0.75384732E-02,-0.36813953E-03, 0.66230184E-03, 0.27458728E-02,
     + -0.62692701E-02, 0.69154210E-01,-0.26262634E-01,-0.36907863E-01,
     + -0.85373910E-03,-0.14128906E-01,-0.13338838E-01, 0.96153775E-02,
     + -0.22918297E-01,-0.42062983E-01, 0.38650357E-02,-0.41835578E-02,
     + -0.71526133E-02,-0.21554792E-01, 0.19192331E-01,-0.56348201E-02,
     + -0.45849726E-03,-0.10513120E-02, 0.16202202E-01,-0.40291692E-03,
     +  0.22353901E-01, 0.41097715E-02,-0.42381627E-02, 0.82951583E-01,
     + -0.81501175E-02,-0.26311824E-03, 0.38785546E-02,-0.17552576E-02,
     + -0.19532156E-02, 0.16496340E-01, 0.25699899E-03,-0.21026812E-02,
     +  0.14631096E-02, 0.50258734E-02,-0.18584801E-01,-0.70805154E-02,
     +  0.89250365E-02, 0.25598228E-01,-0.78319717E-03, 0.21135253E-02,
     +  0.31120123E-01,-0.66352696E-02, 0.59262015E-01,-0.50081369E-02,
     +  0.46478296E-02, 0.84313024E-02, 0.29735785E-01, 0.18269233E-01,
     + -0.16074395E-01,-0.17961324E-02, 0.18536302E-03, 0.52868407E-02,
     + -0.36455330E-01,-0.21586761E-01,-0.84292023E-02,-0.18662501E-01,
     +  0.79084830E-02, 0.15798112E-03,-0.23332266E-01,-0.62831980E-02,
     + -0.42689778E-02,-0.43374533E-02, 0.30126187E-02, 0.43357229E-02,
     +  0.23358289E-01, 0.15905887E-01,-0.31801019E-03, 0.94114320E-03,
     +  0.40128926E-03, 0.38853651E-02, 0.14135040E-02,-0.73852614E-02,
     +  0.40654739E-03,-0.88957166E-02,
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
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_mnq_0_4   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22            
     5  +coeff(  5)        *x32        
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21        *x51
      x_mnq_0_4   =x_mnq_0_4   
     9  +coeff(  9)                *x52
     1  +coeff( 10)*x11                
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x21*x32        
     4  +coeff( 13)    *x21    *x42    
     5  +coeff( 14)        *x32    *x51
     6  +coeff( 15)        *x31*x41*x51
     7  +coeff( 16)            *x42*x51
     8  +coeff( 17)    *x21        *x52
      x_mnq_0_4   =x_mnq_0_4   
     9  +coeff( 18)                *x53
     1  +coeff( 19)    *x22*x31*x41    
     2  +coeff( 20)    *x22    *x42    
     3  +coeff( 21)        *x31*x43    
     4  +coeff( 22)    *x21*x31*x41*x51
     5  +coeff( 23)    *x21    *x42*x51
     6  +coeff( 24)    *x21        *x53
     7  +coeff( 25)    *x23*x31*x41    
     8  +coeff( 26)    *x23    *x42    
      x_mnq_0_4   =x_mnq_0_4   
     9  +coeff( 27)    *x21*x31*x43    
     1  +coeff( 28)        *x31*x43*x51
     2  +coeff( 29)    *x23        *x52
     3  +coeff( 30)    *x22        *x53
     4  +coeff( 31)    *x22*x31*x43    
     5  +coeff( 32)    *x21*x31*x43*x51
     6  +coeff( 33)    *x22*x32    *x52
     7  +coeff( 34)    *x21*x32    *x53
     8  +coeff( 35)    *x21    *x42*x53
      x_mnq_0_4   =x_mnq_0_4   
     9  +coeff( 36)*x11    *x32*x42    
     1  +coeff( 37)    *x22*x32*x42*x52
     2  +coeff( 38)    *x23*x32    *x53
     3  +coeff( 39)    *x21*x33*x41*x53
     4  +coeff( 40)    *x23    *x42*x53
     5  +coeff( 41)*x11*x23    *x42*x52
     6  +coeff( 42)    *x23*x32*x42*x53
     7  +coeff( 43)    *x22        *x51
     8  +coeff( 44)    *x22*x32        
      x_mnq_0_4   =x_mnq_0_4   
     9  +coeff( 45)    *x23        *x51
     1  +coeff( 46)    *x22        *x52
     2  +coeff( 47)        *x32    *x52
     3  +coeff( 48)    *x23*x32        
     4  +coeff( 49)    *x21*x32*x42    
     5  +coeff( 50)        *x33*x41*x51
     6  +coeff( 51)    *x22    *x42*x51
     7  +coeff( 52)    *x21    *x42*x52
     8  +coeff( 53)    *x22*x33*x41    
      x_mnq_0_4   =x_mnq_0_4   
     9  +coeff( 54)    *x23*x31*x41*x51
     1  +coeff( 55)    *x21*x33*x41*x51
     2  +coeff( 56)*x11*x21        *x52
     3  +coeff( 57)    *x22*x31*x41*x52
     4  +coeff( 58)        *x33*x41*x52
     5  +coeff( 59)    *x22    *x42*x52
     6  +coeff( 60)        *x32*x42*x52
     7  +coeff( 61)        *x31*x43*x52
     8  +coeff( 62)    *x23        *x53
      x_mnq_0_4   =x_mnq_0_4   
     9  +coeff( 63)    *x21*x31*x41*x53
     1  +coeff( 64)    *x23*x31*x43    
     2  +coeff( 65)    *x22*x33*x41*x51
     3  +coeff( 66)*x11    *x32    *x52
     4  +coeff( 67)    *x23*x31*x41*x52
     5  +coeff( 68)    *x21*x33*x41*x52
     6  +coeff( 69)    *x23    *x42*x52
     7  +coeff( 70)    *x22    *x42*x53
     8  +coeff( 71)    *x22*x33*x43    
      x_mnq_0_4   =x_mnq_0_4   
     9  +coeff( 72)    *x23*x33*x43    
     1  +coeff( 73)*x11*x23    *x42*x51
     2  +coeff( 74)*x11*x21*x31*x43*x51
     3  +coeff( 75)    *x23*x33*x41*x52
     4  +coeff( 76)*x11*x21*x32    *x53
     5  +coeff( 77)*x11*x21    *x42*x53
     6  +coeff( 78)*x11*x21*x33*x41*x52
     7  +coeff( 79)*x11*x22*x32    *x53
     8  +coeff( 80)*x11*x23*x32    *x53
      x_mnq_0_4   =x_mnq_0_4   
     9  +coeff( 81)*x11*x23*x32*x42*x53
     1  +coeff( 82)*x11*x21*x33*x43*x53
     2  +coeff( 83)*x11*x21            
     3  +coeff( 84)        *x32*x42    
     4  +coeff( 85)*x11*x22            
     5  +coeff( 86)    *x21*x33*x41    
     6  +coeff( 87)    *x22*x32    *x51
     7  +coeff( 88)    *x22*x31*x41*x51
     8  +coeff( 89)*x11            *x52
      x_mnq_0_4   =x_mnq_0_4   
     9  +coeff( 90)    *x21*x31*x41*x52
c
      return
      end
      function t_mnq_0_4   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.9347861E-02/
      data xmin/
     1 -0.49978E-02,-0.50007E-01,-0.10395E+00,-0.50032E-01,-0.49981E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49997E-02, 0.49989E-01, 0.10395E+00, 0.50032E-01, 0.49998E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.10307515E-01, 0.48118085E-01, 0.84617369E-01, 0.33732567E-02,
     +  0.12242235E-03,-0.73899620E-03,-0.16870657E-02, 0.31995706E-02,
     + -0.55441454E-01, 0.15581639E-03, 0.48922030E-02, 0.44227432E-04,
     +  0.22630396E-02, 0.55968044E-02, 0.18746480E-02,-0.29497972E-03,
     +  0.28085392E-02,-0.35712309E-02, 0.29716816E-01,-0.10676944E-01,
     +  0.65755501E-03,-0.17876854E-01,-0.15217491E-03,-0.76711859E-03,
     + -0.71328417E-02,-0.16991582E-01,-0.25850153E-03, 0.16158246E-02,
     + -0.24761824E-03,-0.13263137E-01,-0.27200913E-01, 0.17659473E-02,
     +  0.45589326E-03, 0.13531268E-02,-0.15307083E-01,-0.25152299E-02,
     + -0.56701833E-02,-0.93042972E-02,-0.10077735E-01, 0.36436290E-03,
     + -0.33332620E-03,-0.50443341E-03, 0.81262141E-02,-0.84412162E-03,
     +  0.14746597E-02, 0.14240976E-01, 0.26315117E-01, 0.66760631E-03,
     +  0.35895026E-03, 0.23966176E-02, 0.17110646E-01,-0.28375459E-02,
     +  0.28241372E-01, 0.11128059E-01, 0.13043026E-01, 0.22196496E-01,
     + -0.73682728E-04,-0.16695014E-02,-0.62564813E-03,-0.59187110E-02,
     +  0.64604282E-02,-0.61767339E-03,-0.56453224E-03,-0.86294668E-03,
     + -0.27648977E-03, 0.90684305E-03,-0.97477073E-02, 0.25743622E-03,
     + -0.12500719E-03,-0.66916985E-02,-0.31133119E-02,-0.17546245E-03,
     +  0.10177830E-02, 0.51797531E-02, 0.10138329E-02, 0.13649735E-04,
     + -0.20343296E-02, 0.56952744E-03, 0.38383896E-02,-0.33200231E-05,
     +  0.21620523E-02,-0.25028648E-03,-0.22898830E-03, 0.47777290E-03,
     + -0.18966007E-03,-0.11144805E-03,-0.19246654E-02, 0.33867441E-03,
     +  0.82078547E-03, 0.71928033E-03,
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
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_mnq_0_4   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22            
     5  +coeff(  5)        *x32        
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21        *x51
      t_mnq_0_4   =t_mnq_0_4   
     9  +coeff(  9)                *x52
     1  +coeff( 10)*x11                
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x21*x32        
     4  +coeff( 13)    *x21*x31*x41    
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)    *x22        *x51
     7  +coeff( 16)        *x32    *x51
     8  +coeff( 17)            *x42*x51
      t_mnq_0_4   =t_mnq_0_4   
     9  +coeff( 18)    *x21        *x52
     1  +coeff( 19)                *x53
     2  +coeff( 20)    *x22*x31*x41    
     3  +coeff( 21)        *x33*x41    
     4  +coeff( 22)    *x22    *x42    
     5  +coeff( 23)        *x31*x43    
     6  +coeff( 24)    *x21*x32    *x51
     7  +coeff( 25)    *x21*x31*x41*x51
     8  +coeff( 26)    *x21    *x42*x51
      t_mnq_0_4   =t_mnq_0_4   
     9  +coeff( 27)            *x42*x52
     1  +coeff( 28)    *x21        *x53
     2  +coeff( 29)*x11    *x32        
     3  +coeff( 30)    *x23*x31*x41    
     4  +coeff( 31)    *x23    *x42    
     5  +coeff( 32)    *x21*x31*x43    
     6  +coeff( 33)    *x22*x32    *x51
     7  +coeff( 34)        *x33*x41*x51
     8  +coeff( 35)    *x22    *x42*x51
      t_mnq_0_4   =t_mnq_0_4   
     9  +coeff( 36)        *x31*x43*x51
     1  +coeff( 37)    *x23        *x52
     2  +coeff( 38)    *x21    *x42*x52
     3  +coeff( 39)    *x22        *x53
     4  +coeff( 40)*x11*x23            
     5  +coeff( 41)*x11*x21*x31*x41    
     6  +coeff( 42)*x11*x21    *x42    
     7  +coeff( 43)    *x22*x31*x43    
     8  +coeff( 44)        *x33*x43    
      t_mnq_0_4   =t_mnq_0_4   
     9  +coeff( 45)    *x23*x32    *x51
     1  +coeff( 46)    *x23*x31*x41*x51
     2  +coeff( 47)    *x23    *x42*x51
     3  +coeff( 48)    *x21*x32*x42*x51
     4  +coeff( 49)    *x21*x31*x43*x51
     5  +coeff( 50)    *x22*x32    *x52
     6  +coeff( 51)    *x22*x31*x41*x52
     7  +coeff( 52)        *x33*x41*x52
     8  +coeff( 53)    *x22    *x42*x52
      t_mnq_0_4   =t_mnq_0_4   
     9  +coeff( 54)    *x23        *x53
     1  +coeff( 55)    *x21*x31*x41*x53
     2  +coeff( 56)    *x21    *x42*x53
     3  +coeff( 57)*x11*x21            
     4  +coeff( 58)    *x22*x32        
     5  +coeff( 59)        *x32*x42    
     6  +coeff( 60)    *x23        *x51
     7  +coeff( 61)    *x22        *x52
     8  +coeff( 62)        *x32    *x52
      t_mnq_0_4   =t_mnq_0_4   
     9  +coeff( 63)        *x31*x41*x52
     1  +coeff( 64)    *x23*x32        
     2  +coeff( 65)*x11        *x42    
     3  +coeff( 66)    *x21*x32*x42    
     4  +coeff( 67)    *x22*x31*x41*x51
     5  +coeff( 68)        *x32*x42*x51
     6  +coeff( 69)    *x21*x32    *x52
     7  +coeff( 70)    *x21*x31*x41*x52
     8  +coeff( 71)            *x42*x53
      t_mnq_0_4   =t_mnq_0_4   
     9  +coeff( 72)*x11*x21*x32        
     1  +coeff( 73)    *x22*x33*x41    
     2  +coeff( 74)    *x22*x32*x42    
     3  +coeff( 75)*x11*x22        *x51
     4  +coeff( 76)*x11    *x32    *x51
     5  +coeff( 77)    *x21*x33*x41*x51
     6  +coeff( 78)*x11*x21        *x52
     7  +coeff( 79)        *x31*x43*x52
     8  +coeff( 80)*x11            *x53
      t_mnq_0_4   =t_mnq_0_4   
     9  +coeff( 81)        *x31*x41*x51
     1  +coeff( 82)*x11            *x51
     2  +coeff( 83)*x11*x22            
     3  +coeff( 84)    *x21*x33*x41    
     4  +coeff( 85)*x11*x21        *x51
     5  +coeff( 86)*x11            *x52
     6  +coeff( 87)        *x31*x41*x53
     7  +coeff( 88)*x11        *x42*x51
     8  +coeff( 89)        *x32*x42*x52
      t_mnq_0_4   =t_mnq_0_4   
     9  +coeff( 90)    *x21*x32    *x53
c
      return
      end
      function y_mnq_0_4   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49978E-02,-0.50007E-01,-0.10395E+00,-0.50032E-01,-0.49981E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49997E-02, 0.49989E-01, 0.10395E+00, 0.50032E-01, 0.49998E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.10272294E+00, 0.56159526E+00,-0.10761094E-02, 0.14677473E-02,
     +  0.11547239E-02, 0.46103760E-02,-0.20744160E-01,-0.13182492E-02,
     + -0.21147260E-02,-0.12531849E-02,-0.67182761E-02,-0.57461760E-02,
     + -0.99104065E-02,-0.41617662E-01,-0.45036948E-02,-0.45080911E-02,
     +  0.26541902E-02, 0.22637514E-01, 0.36647744E-01, 0.82671316E-02,
     +  0.30217594E-01, 0.40557110E-02, 0.56676380E-01, 0.75427227E-01,
     + -0.17891162E-02, 0.51038042E-02,-0.44705169E-02,-0.22770149E-03,
     + -0.14606504E-02,-0.27558370E-03,-0.51305396E-03,-0.39667925E-02,
     +  0.79436094E-03, 0.47094086E-02, 0.22883369E-01, 0.70449682E-02,
     +  0.10842499E-01, 0.17601559E-01,-0.81059374E-02,-0.72624860E-02,
     + -0.40927198E-03,-0.12633983E-02, 0.10835516E-01, 0.15197633E-01,
     + -0.20275335E-02, 0.15818447E-01,-0.42387009E-01,-0.53744167E-01,
     + -0.21080276E-01,-0.35108760E-01, 0.40996708E-02,-0.91234286E-03,
     +  0.12334277E-01, 0.16222035E-02,-0.31421576E-02, 0.81389106E-03,
     + -0.20098176E-02, 0.54505328E-03, 0.17846907E-02, 0.25172876E-02,
     +  0.75227971E-03,-0.49897720E-03,-0.88162618E-04,-0.14528410E-01,
     +  0.19070994E-02,-0.10663007E-01,-0.61549086E-04,-0.33690230E-03,
     +  0.15436898E-02,-0.34612841E-02, 0.32422880E-02,-0.62268795E-04,
     +  0.57920176E-02, 0.26353290E-02,-0.27927474E-03,-0.25928352E-03,
     + -0.32993318E-02, 0.11352406E-01, 0.20738907E-01,-0.17298590E-02,
     +  0.15450003E-01, 0.10091391E-02, 0.57410228E-03, 0.18746376E-02,
     + -0.29534849E-02,-0.22262740E-02, 0.17823509E-02,-0.12206252E-01,
     + -0.36751132E-02, 0.26092180E-02,
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
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      y_mnq_0_4   =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21*x31        
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)        *x31    *x51
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)    *x22    *x41    
     8  +coeff(  8)        *x31*x42    
      y_mnq_0_4   =y_mnq_0_4   
     9  +coeff(  9)            *x43    
     1  +coeff( 10)    *x21*x31    *x51
     2  +coeff( 11)    *x21    *x41*x51
     3  +coeff( 12)            *x41*x52
     4  +coeff( 13)    *x23*x31        
     5  +coeff( 14)    *x23    *x41    
     6  +coeff( 15)    *x21*x31*x42    
     7  +coeff( 16)    *x21    *x43    
     8  +coeff( 17)            *x41*x53
      y_mnq_0_4   =y_mnq_0_4   
     9  +coeff( 18)    *x22*x31*x42    
     1  +coeff( 19)    *x22    *x43    
     2  +coeff( 20)    *x23*x31    *x51
     3  +coeff( 21)    *x23    *x41*x51
     4  +coeff( 22)    *x22*x31    *x52
     5  +coeff( 23)    *x23*x31*x42    
     6  +coeff( 24)    *x23    *x43    
     7  +coeff( 25)    *x21*x32*x43    
     8  +coeff( 26)    *x23*x32*x43    
      y_mnq_0_4   =y_mnq_0_4   
     9  +coeff( 27)    *x22*x31        
     1  +coeff( 28)        *x32*x41    
     2  +coeff( 29)        *x31    *x52
     3  +coeff( 30)    *x21*x32*x41    
     4  +coeff( 31)            *x43*x51
     5  +coeff( 32)    *x21    *x41*x52
     6  +coeff( 33)        *x31    *x53
     7  +coeff( 34)    *x22*x32*x41    
     8  +coeff( 35)    *x22    *x41*x52
      y_mnq_0_4   =y_mnq_0_4   
     9  +coeff( 36)    *x21    *x41*x53
     1  +coeff( 37)    *x23*x32*x41    
     2  +coeff( 38)    *x22*x31*x42*x51
     3  +coeff( 39)    *x21    *x43*x53
     4  +coeff( 40)    *x22    *x41*x51
     5  +coeff( 41)        *x31*x42*x51
     6  +coeff( 42)    *x21*x31    *x52
     7  +coeff( 43)    *x21*x31*x42*x51
     8  +coeff( 44)    *x21    *x43*x51
      y_mnq_0_4   =y_mnq_0_4   
     9  +coeff( 45)*x11*x22    *x41    
     1  +coeff( 46)    *x22    *x43*x51
     2  +coeff( 47)    *x23*x31*x42*x51
     3  +coeff( 48)    *x23    *x43*x51
     4  +coeff( 49)    *x22*x31*x42*x52
     5  +coeff( 50)    *x22    *x43*x52
     6  +coeff( 51)*x11*x22    *x43    
     7  +coeff( 52)*x11*x21*x31*x42*x51
     8  +coeff( 53)    *x22*x32*x43*x51
      y_mnq_0_4   =y_mnq_0_4   
     9  +coeff( 54)*x11*x21*x31*x42*x52
     1  +coeff( 55)*x11*x21*x32*x43*x51
     2  +coeff( 56)    *x21*x33        
     3  +coeff( 57)    *x22*x31    *x51
     4  +coeff( 58)    *x22*x33        
     5  +coeff( 59)        *x31*x42*x52
     6  +coeff( 60)            *x43*x52
     7  +coeff( 61)    *x21*x31    *x53
     8  +coeff( 62)*x11*x22*x31        
      y_mnq_0_4   =y_mnq_0_4   
     9  +coeff( 63)*x11*x21*x31    *x51
     1  +coeff( 64)    *x23    *x41*x52
     2  +coeff( 65)    *x21    *x43*x52
     3  +coeff( 66)    *x22    *x41*x53
     4  +coeff( 67)*x11*x21    *x43    
     5  +coeff( 68)*x11*x22*x31    *x51
     6  +coeff( 69)*x11*x22    *x41*x51
     7  +coeff( 70)    *x23*x32*x41*x51
     8  +coeff( 71)    *x21*x32*x43*x51
      y_mnq_0_4   =y_mnq_0_4   
     9  +coeff( 72)*x11*x21*x31    *x52
     1  +coeff( 73)    *x23    *x41*x53
     2  +coeff( 74)*x11*x22*x31*x42    
     3  +coeff( 75)*x11*x21    *x43*x51
     4  +coeff( 76)*x11    *x33    *x52
     5  +coeff( 77)    *x23*x33    *x52
     6  +coeff( 78)    *x23*x31*x42*x52
     7  +coeff( 79)    *x23    *x43*x52
     8  +coeff( 80)    *x22*x33    *x53
      y_mnq_0_4   =y_mnq_0_4   
     9  +coeff( 81)    *x22    *x43*x53
     1  +coeff( 82)*x11*x23*x31*x42    
     2  +coeff( 83)*x11*x21*x32*x43    
     3  +coeff( 84)*x11*x22*x33    *x51
     4  +coeff( 85)*x11*x22    *x43*x51
     5  +coeff( 86)    *x23*x32*x41*x53
     6  +coeff( 87)*x11*x21    *x43*x53
     7  +coeff( 88)    *x22*x32*x43*x53
     8  +coeff( 89)*x11*x22*x33*x42*x51
      y_mnq_0_4   =y_mnq_0_4   
     9  +coeff( 90)*x11*x21*x32*x43*x52
c
      return
      end
      function p_mnq_0_4   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49978E-02,-0.50007E-01,-0.10395E+00,-0.50032E-01,-0.49981E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49997E-02, 0.49989E-01, 0.10395E+00, 0.50032E-01, 0.49998E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.50087098E-03, 0.48063427E-01,-0.10840533E-03,-0.69584476E-03,
     +  0.53277018E-03, 0.23185920E-02,-0.20629824E-02,-0.97954031E-02,
     + -0.77124067E-04,-0.66280138E-03,-0.98452554E-03,-0.57013554E-03,
     + -0.32558418E-02,-0.70762204E-03,-0.25316104E-02,-0.47447504E-02,
     + -0.15818816E-04,-0.21880355E-01,-0.33680993E-03,-0.24850611E-02,
     + -0.23450847E-02,-0.68801302E-02, 0.18090088E-03, 0.56873243E-04,
     + -0.37119110E-03,-0.50764537E-03,-0.29303473E-02, 0.41972488E-03,
     +  0.92124351E-03,-0.24625851E-03, 0.20213155E-02, 0.10560267E-01,
     +  0.48718484E-04, 0.17394945E-01, 0.13969429E-01, 0.49334611E-02,
     +  0.65960847E-02, 0.23149266E-02, 0.11146022E-01, 0.57707075E-03,
     +  0.65089599E-03, 0.26528235E-03, 0.43432419E-02,-0.78787090E-03,
     +  0.58300500E-02, 0.28007083E-01, 0.58588671E-03, 0.94407633E-05,
     +  0.39820399E-01, 0.85270097E-02, 0.13374600E-01,-0.77820075E-03,
     +  0.64778799E-03, 0.21619599E-02,-0.31836212E-02,-0.19433293E-03,
     +  0.67846687E-03, 0.79973548E-03,-0.27947894E-02,-0.17527860E-01,
     + -0.21402724E-01, 0.74489886E-03,-0.93470877E-02,-0.15139341E-01,
     +  0.50360558E-03,-0.30244617E-02, 0.13858232E-02,-0.85178704E-03,
     +  0.11417510E-02, 0.77332149E-03,-0.43129822E-03,-0.77888842E-04,
     + -0.14708658E-02, 0.16835335E-03,-0.63101492E-04, 0.10683864E-02,
     + -0.13687043E-02,-0.88401226E-04, 0.19538036E-03, 0.32203926E-02,
     +  0.91877556E-03, 0.12641934E-03, 0.22067490E-03, 0.14559422E-02,
     + -0.98200275E-04,-0.33624666E-02,-0.88264933E-03,-0.63725497E-03,
     + -0.37187850E-03, 0.96926338E-03,
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
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_mnq_0_4   =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21*x31        
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)        *x31    *x51
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)    *x22*x31        
     8  +coeff(  8)    *x22    *x41    
      p_mnq_0_4   =p_mnq_0_4   
     9  +coeff(  9)        *x32*x41    
     1  +coeff( 10)        *x31*x42    
     2  +coeff( 11)            *x43    
     3  +coeff( 12)    *x21*x31    *x51
     4  +coeff( 13)    *x21    *x41*x51
     5  +coeff( 14)        *x31    *x52
     6  +coeff( 15)            *x41*x52
     7  +coeff( 16)    *x23*x31        
     8  +coeff( 17)*x11        *x41    
      p_mnq_0_4   =p_mnq_0_4   
     9  +coeff( 18)    *x23    *x41    
     1  +coeff( 19)    *x21*x32*x41    
     2  +coeff( 20)    *x21*x31*x42    
     3  +coeff( 21)    *x21    *x43    
     4  +coeff( 22)    *x22    *x41*x51
     5  +coeff( 23)        *x32*x41*x51
     6  +coeff( 24)        *x31*x42*x51
     7  +coeff( 25)            *x43*x51
     8  +coeff( 26)    *x21*x31    *x52
      p_mnq_0_4   =p_mnq_0_4   
     9  +coeff( 27)    *x21    *x41*x52
     1  +coeff( 28)        *x31    *x53
     2  +coeff( 29)            *x41*x53
     3  +coeff( 30)*x11*x21    *x41    
     4  +coeff( 31)    *x22*x32*x41    
     5  +coeff( 32)    *x22*x31*x42    
     6  +coeff( 33)        *x33*x42    
     7  +coeff( 34)    *x22    *x43    
     8  +coeff( 35)    *x23    *x41*x51
      p_mnq_0_4   =p_mnq_0_4   
     9  +coeff( 36)    *x21*x31*x42*x51
     1  +coeff( 37)    *x21    *x43*x51
     2  +coeff( 38)    *x22*x31    *x52
     3  +coeff( 39)    *x22    *x41*x52
     4  +coeff( 40)        *x31*x42*x52
     5  +coeff( 41)            *x43*x52
     6  +coeff( 42)    *x21*x31    *x53
     7  +coeff( 43)    *x21    *x41*x53
     8  +coeff( 44)*x11*x22    *x41    
      p_mnq_0_4   =p_mnq_0_4   
     9  +coeff( 45)    *x23*x32*x41    
     1  +coeff( 46)    *x23*x31*x42    
     2  +coeff( 47)    *x21*x33*x42    
     3  +coeff( 48)*x11        *x43    
     4  +coeff( 49)    *x23    *x43    
     5  +coeff( 50)    *x22*x31*x42*x51
     6  +coeff( 51)    *x22    *x43*x51
     7  +coeff( 52)    *x23*x31    *x52
     8  +coeff( 53)    *x21*x31*x42*x52
      p_mnq_0_4   =p_mnq_0_4   
     9  +coeff( 54)    *x21    *x43*x52
     1  +coeff( 55)    *x22    *x41*x53
     2  +coeff( 56)*x11*x23    *x41    
     3  +coeff( 57)*x11*x21*x31*x42    
     4  +coeff( 58)*x11*x21    *x43    
     5  +coeff( 59)    *x23*x32*x41*x51
     6  +coeff( 60)    *x23*x31*x42*x51
     7  +coeff( 61)    *x23    *x43*x51
     8  +coeff( 62)    *x21*x32*x43*x51
      p_mnq_0_4   =p_mnq_0_4   
     9  +coeff( 63)    *x22*x31*x42*x52
     1  +coeff( 64)    *x22    *x43*x52
     2  +coeff( 65)    *x23*x31    *x53
     3  +coeff( 66)    *x21    *x43*x53
     4  +coeff( 67)*x11*x22    *x43    
     5  +coeff( 68)*x11*x21*x31*x42*x51
     6  +coeff( 69)    *x23*x32*x41*x52
     7  +coeff( 70)    *x23    *x43*x52
     8  +coeff( 71)    *x22*x33    *x53
      p_mnq_0_4   =p_mnq_0_4   
     9  +coeff( 72)*x11*x23*x33        
     1  +coeff( 73)*x11*x21*x32*x43*x51
     2  +coeff( 74)*x11*x23*x31    *x53
     3  +coeff( 75)    *x22*x33*x42*x53
     4  +coeff( 76)*x11*x23    *x43*x53
     5  +coeff( 77)    *x22*x31    *x51
     6  +coeff( 78)*x11*x21*x31        
     7  +coeff( 79)    *x22*x33        
     8  +coeff( 80)    *x23*x31    *x51
      p_mnq_0_4   =p_mnq_0_4   
     9  +coeff( 81)    *x21*x32*x41*x51
     1  +coeff( 82)*x11    *x31*x42    
     2  +coeff( 83)    *x22*x33    *x51
     3  +coeff( 84)    *x22*x32*x41*x51
     4  +coeff( 85)*x11    *x31    *x52
     5  +coeff( 86)    *x23    *x41*x52
     6  +coeff( 87)    *x21*x32*x41*x52
     7  +coeff( 88)    *x22*x31    *x53
     8  +coeff( 89)        *x32*x41*x53
      p_mnq_0_4   =p_mnq_0_4   
     9  +coeff( 90)            *x43*x53
c
      return
      end
      function pl_mnq_0_4  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 41)
      data ncoeff/ 40/
      data avdat/ -0.1246754E-01/
      data xmin/
     1 -0.49978E-02,-0.50007E-01,-0.10395E+00,-0.50032E-01,-0.49981E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49997E-02, 0.49989E-01, 0.10395E+00, 0.50032E-01, 0.49998E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.12388290E-01,-0.29844437E-01, 0.79958681E-02,-0.10097589E-01,
     + -0.14240578E-01,-0.14526688E-01,-0.10359921E-01, 0.16018908E-01,
     +  0.68469294E-02,-0.11094577E-01,-0.44458604E-03, 0.76416315E-03,
     + -0.12274605E-02, 0.19647793E-02,-0.35799772E-03, 0.83692663E-03,
     + -0.12524308E-03, 0.11352124E-02, 0.66085230E-02, 0.57242683E-03,
     + -0.16480290E-02, 0.25447021E-03,-0.47965805E-03,-0.81097445E-03,
     +  0.43334221E-03,-0.11324702E-02,-0.65600249E-03,-0.23190291E-02,
     + -0.16156518E-02,-0.12718051E-02, 0.11103577E-02, 0.57438965E-03,
     + -0.34030522E-02, 0.23531877E-02, 0.42795902E-02,-0.91977231E-03,
     + -0.56851394E-02, 0.49861060E-02, 0.42076223E-02,-0.16839878E-03,
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
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
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
      pl_mnq_0_4  =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)                *x52
     5  +coeff(  5)    *x22            
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)                *x53
      pl_mnq_0_4  =pl_mnq_0_4  
     9  +coeff(  9)    *x21        *x52
     1  +coeff( 10)                *x54
     2  +coeff( 11)*x11                
     3  +coeff( 12)            *x42*x51
     4  +coeff( 13)    *x21        *x53
     5  +coeff( 14)    *x23    *x42*x51
     6  +coeff( 15)        *x31*x41    
     7  +coeff( 16)        *x31*x41*x51
     8  +coeff( 17)    *x21*x31*x41*x51
      pl_mnq_0_4  =pl_mnq_0_4  
     9  +coeff( 18)    *x23    *x42    
     1  +coeff( 19)    *x24    *x42*x51
     2  +coeff( 20)    *x22*x31*x41    
     3  +coeff( 21)    *x22    *x42    
     4  +coeff( 22)    *x23        *x51
     5  +coeff( 23)    *x21    *x42*x51
     6  +coeff( 24)            *x42*x52
     7  +coeff( 25)    *x23*x31*x41    
     8  +coeff( 26)    *x22*x31*x41*x51
      pl_mnq_0_4  =pl_mnq_0_4  
     9  +coeff( 27)        *x31*x43*x51
     1  +coeff( 28)    *x22        *x53
     2  +coeff( 29)    *x21        *x54
     3  +coeff( 30)    *x24*x31*x41    
     4  +coeff( 31)    *x24    *x42    
     5  +coeff( 32)    *x23*x31*x41*x51
     6  +coeff( 33)    *x22*x31*x41*x52
     7  +coeff( 34)    *x22        *x54
     8  +coeff( 35)    *x24*x31*x41*x51
      pl_mnq_0_4  =pl_mnq_0_4  
     9  +coeff( 36)    *x23        *x54
     1  +coeff( 37)    *x24    *x42*x52
     2  +coeff( 38)    *x22*x31*x43*x52
     3  +coeff( 39)    *x22    *x44*x52
     4  +coeff( 40)    *x23            
c
      return
      end
      function x_mnq_0_5   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.4886806E-01/
      data xmin/
     1 -0.49978E-02,-0.50007E-01,-0.10395E+00,-0.50032E-01,-0.49981E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49997E-02, 0.49989E-01, 0.10395E+00, 0.50032E-01, 0.49998E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.51591255E-01, 0.58590961E+00, 0.22722895E+00, 0.12594271E-01,
     + -0.79686142E-04, 0.25896056E-03,-0.38017631E-02, 0.21901133E-01,
     + -0.14052455E+00, 0.51440326E-02, 0.80422303E-02, 0.14055386E-02,
     +  0.62080339E-03, 0.51832716E-02,-0.17071217E-02,-0.75768476E-03,
     +  0.25006668E-02,-0.17119745E-01, 0.70208810E-01,-0.26808986E-02,
     + -0.20985581E-01,-0.49224935E-01,-0.19125254E-02,-0.56213055E-01,
     +  0.15055812E-01,-0.11868570E-02, 0.42131268E-02,-0.20530969E-01,
     + -0.59628900E-01, 0.20094444E-02,-0.34687579E-01, 0.36201214E-02,
     +  0.10052244E+00, 0.26347358E-02, 0.11539402E+00, 0.38631078E-01,
     +  0.78894235E-01, 0.81677974E-03,-0.15029057E-01, 0.45253988E-01,
     +  0.45397684E-01,-0.59425924E-02, 0.63240066E-01,-0.34278247E-02,
     + -0.26085760E-01, 0.29117111E-03,-0.81165535E-02,-0.42570726E-03,
     + -0.46292269E-02,-0.22652993E-03,-0.67732618E-02, 0.37535168E-02,
     +  0.92899408E-02, 0.38611505E-02,-0.10484250E-03, 0.31809360E-01,
     + -0.12876674E-01,-0.97469315E-02,-0.59819507E-03, 0.10000409E-01,
     + -0.73430716E-03,-0.40182441E-01, 0.26459333E-02,-0.86412048E-02,
     + -0.23521250E-01,-0.91884080E-02, 0.70439056E-01, 0.28816866E-01,
     + -0.81982939E-02,-0.35806664E-02,-0.25943499E-01,-0.48007178E-02,
     + -0.84854066E-02,-0.63401666E-02, 0.39927796E-01, 0.94787180E-02,
     +  0.15611955E-02, 0.23235958E-01, 0.12264060E-01,-0.12268126E-01,
     +  0.17495859E-01, 0.24391305E-01, 0.23577120E-01, 0.79961820E-03,
     + -0.94477477E-03, 0.34823039E-03,-0.12836922E-02, 0.56359195E-03,
     +  0.18928079E-02, 0.15644344E-02,
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
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_mnq_0_5   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22            
     5  +coeff(  5)        *x32        
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21        *x51
      x_mnq_0_5   =x_mnq_0_5   
     9  +coeff(  9)                *x52
     1  +coeff( 10)*x11                
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x21*x32        
     4  +coeff( 13)    *x21*x31*x41    
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)    *x22        *x51
     7  +coeff( 16)        *x32    *x51
     8  +coeff( 17)            *x42*x51
      x_mnq_0_5   =x_mnq_0_5   
     9  +coeff( 18)    *x21        *x52
     1  +coeff( 19)                *x53
     2  +coeff( 20)    *x22*x32        
     3  +coeff( 21)    *x22*x31*x41    
     4  +coeff( 22)    *x22    *x42    
     5  +coeff( 23)        *x31*x43    
     6  +coeff( 24)    *x21    *x42*x51
     7  +coeff( 25)    *x22        *x52
     8  +coeff( 26)            *x42*x52
      x_mnq_0_5   =x_mnq_0_5   
     9  +coeff( 27)    *x21        *x53
     1  +coeff( 28)    *x23*x31*x41    
     2  +coeff( 29)    *x23    *x42    
     3  +coeff( 30)    *x21*x31*x43    
     4  +coeff( 31)    *x22    *x42*x51
     5  +coeff( 32)    *x22*x31*x43    
     6  +coeff( 33)    *x23    *x42*x51
     7  +coeff( 34)    *x21*x32*x42*x51
     8  +coeff( 35)    *x22    *x42*x52
      x_mnq_0_5   =x_mnq_0_5   
     9  +coeff( 36)    *x23        *x53
     1  +coeff( 37)    *x21    *x42*x53
     2  +coeff( 38)    *x23*x33*x41*x51
     3  +coeff( 39)    *x21*x33*x43*x51
     4  +coeff( 40)    *x22*x33*x41*x52
     5  +coeff( 41)    *x21*x33*x41*x53
     6  +coeff( 42)    *x23    *x42*x53
     7  +coeff( 43)    *x23*x33*x43*x51
     8  +coeff( 44)        *x33*x41    
      x_mnq_0_5   =x_mnq_0_5   
     9  +coeff( 45)    *x23        *x51
     1  +coeff( 46)    *x21*x32    *x51
     2  +coeff( 47)    *x21*x31*x41*x51
     3  +coeff( 48)*x11    *x32        
     4  +coeff( 49)    *x23*x32        
     5  +coeff( 50)    *x21    *x42*x52
     6  +coeff( 51)    *x22        *x53
     7  +coeff( 52)    *x22*x33*x41    
     8  +coeff( 53)    *x22*x32*x42    
      x_mnq_0_5   =x_mnq_0_5   
     9  +coeff( 54)        *x33*x43    
     1  +coeff( 55)*x11*x22        *x51
     2  +coeff( 56)    *x23*x31*x41*x51
     3  +coeff( 57)    *x21*x33*x41*x51
     4  +coeff( 58)    *x21*x31*x43*x51
     5  +coeff( 59)*x11*x22    *x42    
     6  +coeff( 60)    *x23*x32*x42    
     7  +coeff( 61)    *x23*x31*x43    
     8  +coeff( 62)    *x22*x31*x43*x51
      x_mnq_0_5   =x_mnq_0_5   
     9  +coeff( 63)        *x33*x43*x51
     1  +coeff( 64)    *x23*x31*x41*x52
     2  +coeff( 65)    *x23    *x42*x52
     3  +coeff( 66)*x11*x22    *x42*x51
     4  +coeff( 67)    *x22*x31*x43*x52
     5  +coeff( 68)    *x21*x31*x43*x53
     6  +coeff( 69)*x11*x23    *x42*x51
     7  +coeff( 70)*x11    *x32*x42*x52
     8  +coeff( 71)    *x23*x31*x43*x52
      x_mnq_0_5   =x_mnq_0_5   
     9  +coeff( 72)*x11*x21*x32    *x53
     1  +coeff( 73)*x11*x21    *x42*x53
     2  +coeff( 74)*x11*x23*x33*x41    
     3  +coeff( 75)*x11*x22*x32*x42*x51
     4  +coeff( 76)*x11*x23*x31*x41*x52
     5  +coeff( 77)*x11*x21*x33*x41*x52
     6  +coeff( 78)*x11*x21*x32*x42*x52
     7  +coeff( 79)*x11*x22    *x42*x53
     8  +coeff( 80)*x11*x22*x33*x43    
      x_mnq_0_5   =x_mnq_0_5   
     9  +coeff( 81)*x11*x23*x32*x42*x51
     1  +coeff( 82)*x11*x23    *x42*x53
     2  +coeff( 83)*x11*x22*x33*x43*x51
     3  +coeff( 84)        *x31*x41*x51
     4  +coeff( 85)        *x31*x41*x52
     5  +coeff( 86)*x11    *x31*x41    
     6  +coeff( 87)    *x21*x32*x42    
     7  +coeff( 88)*x11*x21        *x51
     8  +coeff( 89)    *x22*x32    *x51
      x_mnq_0_5   =x_mnq_0_5   
     9  +coeff( 90)        *x32*x42*x51
c
      return
      end
      function t_mnq_0_5   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.1385349E-01/
      data xmin/
     1 -0.49978E-02,-0.50007E-01,-0.10395E+00,-0.50032E-01,-0.49981E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49997E-02, 0.49989E-01, 0.10395E+00, 0.50032E-01, 0.49998E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.14801719E-01, 0.47283743E-01, 0.88438489E-01, 0.65008062E-02,
     +  0.10460950E-03, 0.15044415E-02, 0.37368152E-02, 0.77651469E-02,
     + -0.55742893E-01, 0.12392699E-03, 0.40503656E-02, 0.17651387E-03,
     +  0.14676850E-02, 0.49947146E-02,-0.32073897E-02,-0.36962054E-03,
     + -0.97693480E-03,-0.81762876E-02, 0.27880093E-01,-0.11771136E-02,
     + -0.11504877E-01,-0.58594049E-03,-0.21788886E-01,-0.35077486E-04,
     + -0.44170304E-03,-0.26039319E-01, 0.51831603E-02,-0.18434416E-03,
     +  0.29378373E-02,-0.22554562E-03,-0.10335612E-02,-0.11287630E-01,
     + -0.26066350E-01, 0.21371325E-02,-0.12129402E-01,-0.20756666E-02,
     + -0.10127476E-02, 0.78559788E-04,-0.19140299E-03, 0.47645760E-02,
     +  0.10347868E-02, 0.28380703E-02, 0.41562904E-01, 0.19508387E-02,
     + -0.20701528E-03, 0.27501473E-01, 0.48870143E-01, 0.20037209E-03,
     +  0.15825156E-01, 0.19145364E-01, 0.32564513E-01, 0.37118827E-03,
     + -0.20234330E-03,-0.12092355E-01,-0.63105527E-03,-0.13193847E-01,
     + -0.20823362E-02, 0.69243315E-03, 0.17549563E-02, 0.94415416E-03,
     + -0.62605795E-02, 0.13477776E-02, 0.75550831E-03,-0.15474553E-03,
     + -0.31864722E-02, 0.33943243E-02, 0.14777883E-03, 0.56167557E-02,
     +  0.21969518E-01, 0.20486838E-02, 0.43021765E-03, 0.74527366E-03,
     +  0.23890035E-02,-0.89733745E-03, 0.45816749E-02,-0.30321104E-03,
     + -0.24501226E-03,-0.15440852E-03,-0.20964880E-03,-0.69193321E-03,
     +  0.13909943E-03,-0.21193202E-02,-0.56889089E-03, 0.36345617E-03,
     +  0.97222102E-03, 0.17328515E-03, 0.29581774E-03,-0.16623126E-03,
     + -0.40400053E-04,-0.10271765E-03,
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
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_mnq_0_5   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22            
     5  +coeff(  5)        *x32        
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21        *x51
      t_mnq_0_5   =t_mnq_0_5   
     9  +coeff(  9)                *x52
     1  +coeff( 10)*x11                
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x21*x32        
     4  +coeff( 13)    *x21*x31*x41    
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)    *x22        *x51
     7  +coeff( 16)        *x32    *x51
     8  +coeff( 17)            *x42*x51
      t_mnq_0_5   =t_mnq_0_5   
     9  +coeff( 18)    *x21        *x52
     1  +coeff( 19)                *x53
     2  +coeff( 20)    *x22*x32        
     3  +coeff( 21)    *x22*x31*x41    
     4  +coeff( 22)        *x33*x41    
     5  +coeff( 23)    *x22    *x42    
     6  +coeff( 24)        *x32*x42    
     7  +coeff( 25)        *x31*x43    
     8  +coeff( 26)    *x21    *x42*x51
      t_mnq_0_5   =t_mnq_0_5   
     9  +coeff( 27)    *x22        *x52
     1  +coeff( 28)            *x42*x52
     2  +coeff( 29)    *x21        *x53
     3  +coeff( 30)*x11    *x32        
     4  +coeff( 31)    *x23*x32        
     5  +coeff( 32)    *x23*x31*x41    
     6  +coeff( 33)    *x23    *x42    
     7  +coeff( 34)    *x21*x31*x43    
     8  +coeff( 35)    *x22    *x42*x51
      t_mnq_0_5   =t_mnq_0_5   
     9  +coeff( 36)        *x31*x43*x51
     1  +coeff( 37)    *x22        *x53
     2  +coeff( 38)            *x42*x53
     3  +coeff( 39)*x11*x21*x31*x41    
     4  +coeff( 40)    *x22*x31*x43    
     5  +coeff( 41)*x11*x22        *x51
     6  +coeff( 42)    *x23*x32    *x51
     7  +coeff( 43)    *x23    *x42*x51
     8  +coeff( 44)    *x21*x32*x42*x51
      t_mnq_0_5   =t_mnq_0_5   
     9  +coeff( 45)    *x21*x31*x43*x51
     1  +coeff( 46)    *x22*x31*x41*x52
     2  +coeff( 47)    *x22    *x42*x52
     3  +coeff( 48)        *x32*x42*x52
     4  +coeff( 49)    *x23        *x53
     5  +coeff( 50)    *x21*x31*x41*x53
     6  +coeff( 51)    *x21    *x42*x53
     7  +coeff( 52)        *x31*x41*x51
     8  +coeff( 53)*x11*x21            
      t_mnq_0_5   =t_mnq_0_5   
     9  +coeff( 54)    *x23        *x51
     1  +coeff( 55)    *x21*x32    *x51
     2  +coeff( 56)    *x21*x31*x41*x51
     3  +coeff( 57)        *x31*x41*x52
     4  +coeff( 58)    *x21*x33*x41    
     5  +coeff( 59)    *x21*x32*x42    
     6  +coeff( 60)    *x22*x32    *x51
     7  +coeff( 61)    *x22*x31*x41*x51
     8  +coeff( 62)        *x33*x41*x51
      t_mnq_0_5   =t_mnq_0_5   
     9  +coeff( 63)        *x32*x42*x51
     1  +coeff( 64)*x11            *x52
     2  +coeff( 65)    *x21*x31*x41*x52
     3  +coeff( 66)    *x22*x33*x41    
     4  +coeff( 67)*x11*x21    *x42    
     5  +coeff( 68)    *x22*x32*x42    
     6  +coeff( 69)    *x23*x31*x41*x51
     7  +coeff( 70)    *x21*x33*x41*x51
     8  +coeff( 71)*x11        *x42*x51
      t_mnq_0_5   =t_mnq_0_5   
     9  +coeff( 72)*x11*x21        *x52
     1  +coeff( 73)    *x22*x32    *x52
     2  +coeff( 74)        *x33*x41*x52
     3  +coeff( 75)        *x31*x43*x52
     4  +coeff( 76)*x11            *x51
     5  +coeff( 77)*x11*x22            
     6  +coeff( 78)*x11        *x42    
     7  +coeff( 79)*x11*x21        *x51
     8  +coeff( 80)    *x23        *x52
      t_mnq_0_5   =t_mnq_0_5   
     9  +coeff( 81)    *x21*x32    *x52
     1  +coeff( 82)    *x21    *x42*x52
     2  +coeff( 83)        *x31*x41*x53
     3  +coeff( 84)*x11*x23            
     4  +coeff( 85)        *x33*x43    
     5  +coeff( 86)*x11    *x31*x41*x51
     6  +coeff( 87)*x11            *x53
     7  +coeff( 88)        *x32    *x52
     8  +coeff( 89)*x11    *x31*x41    
      t_mnq_0_5   =t_mnq_0_5   
     9  +coeff( 90)        *x32    *x53
c
      return
      end
      function y_mnq_0_5   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49978E-02,-0.50007E-01,-0.10395E+00,-0.50032E-01,-0.49981E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49997E-02, 0.49989E-01, 0.10395E+00, 0.50032E-01, 0.49998E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.10244253E+00, 0.57605666E+00,-0.17090220E-02, 0.38065112E-02,
     +  0.65415157E-02,-0.56291916E-02,-0.23931490E-01,-0.19123659E-02,
     + -0.11540441E-01,-0.74655777E-02,-0.10835983E-01,-0.42546481E-01,
     + -0.16404891E-02,-0.27256517E-03,-0.10359217E-02, 0.14532974E-03,
     +  0.69491804E-03,-0.17204027E-02, 0.39316625E-02, 0.29662762E-01,
     +  0.46387456E-01, 0.11132294E-01, 0.46983805E-01, 0.19309981E-01,
     +  0.22900533E-02, 0.62741727E-01, 0.87616913E-01,-0.27772379E-02,
     + -0.21354067E-02,-0.60335055E-01, 0.35343855E-02, 0.10050565E-02,
     +  0.10575056E-02, 0.11278954E-02,-0.10844517E-02, 0.58063203E-02,
     +  0.24551596E-01, 0.28618881E-01, 0.99143852E-02,-0.75364351E-01,
     + -0.51003419E-01,-0.16645456E-02, 0.60206320E-03, 0.98592862E-02,
     +  0.12271386E-01, 0.12314156E-01,-0.13812848E-01,-0.32961377E-03,
     + -0.14953990E-01,-0.89822942E-02, 0.30204933E-03,-0.38380254E-01,
     + -0.11932080E-02,-0.24214613E-02, 0.27948632E-03, 0.32386356E-02,
     +  0.76243076E-02, 0.31504536E-02, 0.18921643E-02,-0.18674156E-02,
     + -0.19354853E-03,-0.51317114E-04,-0.15646077E-02,-0.11953664E-02,
     + -0.88920438E-03, 0.85373892E-03, 0.13557913E-03, 0.80536981E-03,
     +  0.12261821E-02,-0.82484316E-02,-0.11252233E-01, 0.11122139E-03,
     +  0.24679455E-02, 0.23124563E-02, 0.48538411E-03,-0.23247767E-02,
     +  0.10024324E-02, 0.87248534E-02, 0.89735218E-03, 0.21116654E-02,
     + -0.48231883E-02, 0.54626330E-02,-0.37201268E-02,-0.23947046E-02,
     +  0.56418432E-02,-0.12742016E-03, 0.40311387E-03, 0.13155605E-03,
     +  0.31457545E-03, 0.11755099E-03,
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
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      y_mnq_0_5   =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21*x31        
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)            *x41*x51
     6  +coeff(  6)    *x22*x31        
     7  +coeff(  7)    *x22    *x41    
     8  +coeff(  8)            *x43    
      y_mnq_0_5   =y_mnq_0_5   
     9  +coeff(  9)    *x21    *x41*x51
     1  +coeff( 10)            *x41*x52
     2  +coeff( 11)    *x23*x31        
     3  +coeff( 12)    *x23    *x41    
     4  +coeff( 13)    *x21*x31*x42    
     5  +coeff( 14)    *x21    *x43    
     6  +coeff( 15)    *x22*x31    *x51
     7  +coeff( 16)        *x33    *x51
     8  +coeff( 17)        *x31*x42*x51
      y_mnq_0_5   =y_mnq_0_5   
     9  +coeff( 18)    *x21    *x41*x52
     1  +coeff( 19)            *x41*x53
     2  +coeff( 20)    *x22*x31*x42    
     3  +coeff( 21)    *x22    *x43    
     4  +coeff( 22)    *x23*x31    *x51
     5  +coeff( 23)    *x23    *x41*x51
     6  +coeff( 24)    *x21*x31*x42*x51
     7  +coeff( 25)        *x31*x42*x52
     8  +coeff( 26)    *x23*x31*x42    
      y_mnq_0_5   =y_mnq_0_5   
     9  +coeff( 27)    *x23    *x43    
     1  +coeff( 28)    *x21*x32*x43    
     2  +coeff( 29)    *x22*x32*x43    
     3  +coeff( 30)    *x23*x31*x42*x51
     4  +coeff( 31)    *x23*x32*x43    
     5  +coeff( 32)        *x31    *x51
     6  +coeff( 33)    *x21*x32*x41    
     7  +coeff( 34)            *x43*x51
     8  +coeff( 35)    *x21*x31    *x52
      y_mnq_0_5   =y_mnq_0_5   
     9  +coeff( 36)    *x22*x32*x41    
     1  +coeff( 37)    *x21    *x43*x51
     2  +coeff( 38)    *x22    *x41*x52
     3  +coeff( 39)    *x21    *x41*x53
     4  +coeff( 40)    *x23    *x43*x51
     5  +coeff( 41)    *x22    *x43*x52
     6  +coeff( 42)        *x31    *x52
     7  +coeff( 43)        *x31    *x53
     8  +coeff( 44)    *x23*x32*x41    
      y_mnq_0_5   =y_mnq_0_5   
     9  +coeff( 45)    *x22*x31*x42*x51
     1  +coeff( 46)    *x22    *x43*x51
     2  +coeff( 47)    *x23    *x41*x52
     3  +coeff( 48)    *x21    *x43*x52
     4  +coeff( 49)    *x22    *x41*x53
     5  +coeff( 50)    *x23*x32*x41*x51
     6  +coeff( 51)    *x21*x32*x43*x51
     7  +coeff( 52)    *x22*x31*x42*x52
     8  +coeff( 53)        *x31*x42    
      y_mnq_0_5   =y_mnq_0_5   
     9  +coeff( 54)    *x21*x31    *x51
     1  +coeff( 55)        *x32*x41*x51
     2  +coeff( 56)    *x21*x32*x41*x51
     3  +coeff( 57)    *x22*x31    *x52
     4  +coeff( 58)            *x43*x52
     5  +coeff( 59)    *x21*x31    *x53
     6  +coeff( 60)*x11*x22    *x41    
     7  +coeff( 61)*x11    *x32*x41    
     8  +coeff( 62)*x11        *x43    
      y_mnq_0_5   =y_mnq_0_5   
     9  +coeff( 63)    *x21*x31*x42*x52
     1  +coeff( 64)*x11*x23    *x41    
     2  +coeff( 65)*x11*x21*x32*x41    
     3  +coeff( 66)*x11*x21*x31*x42    
     4  +coeff( 67)*x11*x21    *x43    
     5  +coeff( 68)*x11*x22*x31    *x51
     6  +coeff( 69)    *x23*x31    *x53
     7  +coeff( 70)    *x21*x31*x42*x53
     8  +coeff( 71)    *x21    *x43*x53
      y_mnq_0_5   =y_mnq_0_5   
     9  +coeff( 72)*x11*x22*x32*x41    
     1  +coeff( 73)*x11*x22*x31*x42    
     2  +coeff( 74)*x11*x22    *x43    
     3  +coeff( 75)*x11*x23*x31    *x51
     4  +coeff( 76)*x11*x22*x31    *x52
     5  +coeff( 77)    *x23*x31*x42*x52
     6  +coeff( 78)    *x22    *x43*x53
     7  +coeff( 79)*x11*x23    *x43    
     8  +coeff( 80)*x11*x21*x32*x43    
      y_mnq_0_5   =y_mnq_0_5   
     9  +coeff( 81)*x11*x22*x31*x42*x51
     1  +coeff( 82)*x11*x22*x32*x43    
     2  +coeff( 83)*x11*x23*x31*x42*x51
     3  +coeff( 84)*x11*x23*x31    *x53
     4  +coeff( 85)*x11*x23*x32*x43    
     5  +coeff( 86)        *x32*x41    
     6  +coeff( 87)    *x21*x33        
     7  +coeff( 88)*x11        *x41    
     8  +coeff( 89)    *x22*x33        
      y_mnq_0_5   =y_mnq_0_5   
     9  +coeff( 90)*x11    *x31    *x51
c
      return
      end
      function p_mnq_0_5   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49978E-02,-0.50007E-01,-0.10395E+00,-0.50032E-01,-0.49981E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49997E-02, 0.49989E-01, 0.10395E+00, 0.50032E-01, 0.49998E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.10557618E-02, 0.45534100E-01,-0.18197171E-02,-0.91905193E-02,
     +  0.64203865E-03,-0.21552076E-02,-0.10501133E-01,-0.11073837E-03,
     + -0.53989887E-03,-0.89882198E-03,-0.10977680E-03,-0.80949959E-03,
     +  0.82564290E-03,-0.42408537E-02,-0.60775612E-04,-0.18520251E-01,
     + -0.74147811E-03,-0.78988017E-03,-0.86371432E-03, 0.13111853E-03,
     + -0.24579749E-02,-0.28229484E-04, 0.10256266E-03, 0.21463186E-03,
     + -0.32389082E-03,-0.87118259E-03,-0.98792394E-03, 0.10640053E-03,
     +  0.63150497E-06, 0.28816937E-02, 0.10190853E-01,-0.35069446E-04,
     +  0.16310969E-01, 0.15456439E-01, 0.53079790E-02, 0.70175459E-02,
     + -0.42409578E-04, 0.97724162E-02, 0.69490995E-03, 0.10313023E-02,
     + -0.11182953E-03, 0.45673558E-02, 0.23927625E-01, 0.34803465E-01,
     +  0.73899250E-02, 0.68687787E-02,-0.42317226E-02,-0.11022280E-02,
     +  0.39409005E-03,-0.19743644E-01,-0.25492633E-01,-0.20545369E-02,
     + -0.10908818E-01,-0.14417303E-01,-0.15752499E-03, 0.14326827E-02,
     +  0.20618443E-03, 0.35120891E-02,-0.24304734E-03, 0.36243761E-02,
     + -0.32615603E-03, 0.10861581E-02,-0.14009811E-02, 0.24127112E-03,
     + -0.42761001E-02, 0.22909644E-02, 0.50996910E-02, 0.10993242E-04,
     + -0.18656156E-02, 0.20490236E-05, 0.47364854E-04,-0.50424169E-04,
     +  0.18464215E-03,-0.58385449E-04, 0.33956559E-02,-0.15208429E-03,
     +  0.53146785E-04, 0.65514864E-03, 0.20323754E-02, 0.14533187E-02,
     + -0.38650320E-03,-0.15058876E-03, 0.83654822E-03,-0.36785237E-02,
     + -0.99936698E-03,-0.13166573E-02,-0.17374675E-03,-0.21804649E-03,
     + -0.36682165E-03, 0.30240789E-03,
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
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_mnq_0_5   =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21*x31        
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)            *x41*x51
     6  +coeff(  6)    *x22*x31        
     7  +coeff(  7)    *x22    *x41    
     8  +coeff(  8)        *x32*x41    
      p_mnq_0_5   =p_mnq_0_5   
     9  +coeff(  9)        *x31*x42    
     1  +coeff( 10)            *x43    
     2  +coeff( 11)    *x21*x31    *x51
     3  +coeff( 12)    *x21    *x41*x51
     4  +coeff( 13)            *x41*x52
     5  +coeff( 14)    *x23*x31        
     6  +coeff( 15)*x11        *x41    
     7  +coeff( 16)    *x23    *x41    
     8  +coeff( 17)    *x21*x31*x42    
      p_mnq_0_5   =p_mnq_0_5   
     9  +coeff( 18)    *x21    *x43    
     1  +coeff( 19)    *x22*x31    *x51
     2  +coeff( 20)        *x33    *x51
     3  +coeff( 21)    *x22    *x41*x51
     4  +coeff( 22)        *x32*x41*x51
     5  +coeff( 23)        *x31*x42*x51
     6  +coeff( 24)            *x43*x51
     7  +coeff( 25)    *x21*x31    *x52
     8  +coeff( 26)    *x21    *x41*x52
      p_mnq_0_5   =p_mnq_0_5   
     9  +coeff( 27)            *x41*x53
     1  +coeff( 28)    *x22*x33        
     2  +coeff( 29)*x11*x21    *x41    
     3  +coeff( 30)    *x22*x32*x41    
     4  +coeff( 31)    *x22*x31*x42    
     5  +coeff( 32)        *x33*x42    
     6  +coeff( 33)    *x22    *x43    
     7  +coeff( 34)    *x23    *x41*x51
     8  +coeff( 35)    *x21*x31*x42*x51
      p_mnq_0_5   =p_mnq_0_5   
     9  +coeff( 36)    *x21    *x43*x51
     1  +coeff( 37)        *x33    *x52
     2  +coeff( 38)    *x22    *x41*x52
     3  +coeff( 39)        *x31*x42*x52
     4  +coeff( 40)            *x43*x52
     5  +coeff( 41)    *x23*x33        
     6  +coeff( 42)    *x23*x32*x41    
     7  +coeff( 43)    *x23*x31*x42    
     8  +coeff( 44)    *x23    *x43    
      p_mnq_0_5   =p_mnq_0_5   
     9  +coeff( 45)    *x22*x31*x42*x51
     1  +coeff( 46)    *x22    *x43*x51
     2  +coeff( 47)    *x22    *x41*x53
     3  +coeff( 48)    *x22*x32*x43    
     4  +coeff( 49)    *x23*x33    *x51
     5  +coeff( 50)    *x23*x31*x42*x51
     6  +coeff( 51)    *x23    *x43*x51
     7  +coeff( 52)    *x22*x32*x41*x52
     8  +coeff( 53)    *x22*x31*x42*x52
      p_mnq_0_5   =p_mnq_0_5   
     9  +coeff( 54)    *x22    *x43*x52
     1  +coeff( 55)    *x21*x33    *x53
     2  +coeff( 56)*x11*x22*x31*x42    
     3  +coeff( 57)*x11    *x32*x43    
     4  +coeff( 58)    *x22*x32*x43*x51
     5  +coeff( 59)*x11*x22*x31    *x52
     6  +coeff( 60)    *x23*x31*x42*x52
     7  +coeff( 61)    *x23    *x43*x52
     8  +coeff( 62)*x11*x21*x32*x43    
      p_mnq_0_5   =p_mnq_0_5   
     9  +coeff( 63)*x11*x22*x31*x42*x51
     1  +coeff( 64)*x11    *x33    *x53
     2  +coeff( 65)    *x23*x32*x41*x53
     3  +coeff( 66)*x11*x22*x32*x43    
     4  +coeff( 67)    *x23*x32*x43*x52
     5  +coeff( 68)*x11*x23    *x41*x53
     6  +coeff( 69)*x11*x22*x32*x43*x51
     7  +coeff( 70)        *x33        
     8  +coeff( 71)        *x31    *x52
      p_mnq_0_5   =p_mnq_0_5   
     9  +coeff( 72)*x11    *x31        
     1  +coeff( 73)    *x21*x33        
     2  +coeff( 74)        *x31    *x53
     3  +coeff( 75)    *x23*x31    *x51
     4  +coeff( 76)    *x21*x33    *x51
     5  +coeff( 77)*x11        *x41*x51
     6  +coeff( 78)    *x21*x32*x41*x51
     7  +coeff( 79)    *x22*x31    *x52
     8  +coeff( 80)    *x21    *x41*x53
      p_mnq_0_5   =p_mnq_0_5   
     9  +coeff( 81)*x11*x22    *x41    
     1  +coeff( 82)*x11    *x32*x41    
     2  +coeff( 83)    *x22*x32*x41*x51
     3  +coeff( 84)    *x23    *x41*x52
     4  +coeff( 85)    *x21*x31*x42*x52
     5  +coeff( 86)    *x21    *x43*x52
     6  +coeff( 87)        *x33    *x53
     7  +coeff( 88)*x11*x23    *x41    
     8  +coeff( 89)*x11*x21*x32*x41    
      p_mnq_0_5   =p_mnq_0_5   
     9  +coeff( 90)*x11*x21*x31*x42    
c
      return
      end
      function pl_mnq_0_5  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 41)
      data ncoeff/ 40/
      data avdat/  0.1616230E-02/
      data xmin/
     1 -0.49978E-02,-0.50007E-01,-0.10395E+00,-0.50032E-01,-0.49981E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49997E-02, 0.49989E-01, 0.10395E+00, 0.50032E-01, 0.49998E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.15343524E-02,-0.18047170E+00,-0.48322819E-01, 0.15751665E-01,
     + -0.13716107E-01,-0.17454838E-02,-0.17282393E-01,-0.14441638E-01,
     +  0.72444165E-02, 0.77057863E-04, 0.68063205E-02, 0.27117747E-03,
     + -0.13664544E-02, 0.33985469E-02, 0.15293190E-02,-0.20602141E-02,
     +  0.14007235E-01,-0.85716468E-03,-0.93261953E-02,-0.14646087E-02,
     + -0.21784303E-02, 0.15103164E-02,-0.42618951E-03,-0.16749422E-02,
     + -0.23751582E-03,-0.20391054E-02, 0.26351623E-02,-0.33316873E-02,
     +  0.35289451E-02,-0.86485213E-02,-0.21051102E-03,-0.39114512E-03,
     +  0.28246508E-02,-0.10220166E-02,-0.16488610E-02,-0.10705828E-02,
     +  0.11136318E-01,-0.65371613E-02, 0.64834249E-02,-0.17506349E-02,
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
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
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
      x54 = x53*x5
c
c                  function
c
      pl_mnq_0_5  =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)                *x52
     5  +coeff(  5)            *x42    
     6  +coeff(  6)*x11                
     7  +coeff(  7)    *x22            
     8  +coeff(  8)    *x21        *x51
      pl_mnq_0_5  =pl_mnq_0_5  
     9  +coeff(  9)    *x21        *x52
     1  +coeff( 10)    *x21    *x42    
     2  +coeff( 11)    *x22    *x42    
     3  +coeff( 12)        *x31*x41    
     4  +coeff( 13)                *x53
     5  +coeff( 14)    *x22*x31*x41    
     6  +coeff( 15)    *x22        *x52
     7  +coeff( 16)                *x54
     8  +coeff( 17)    *x23    *x42    
      pl_mnq_0_5  =pl_mnq_0_5  
     9  +coeff( 18)    *x21    *x44    
     1  +coeff( 19)    *x22    *x44    
     2  +coeff( 20)    *x24*x31*x43    
     3  +coeff( 21)    *x23            
     4  +coeff( 22)    *x21*x31*x41    
     5  +coeff( 23)    *x22        *x51
     6  +coeff( 24)    *x24            
     7  +coeff( 25)        *x31*x43    
     8  +coeff( 26)    *x21        *x53
      pl_mnq_0_5  =pl_mnq_0_5  
     9  +coeff( 27)    *x23*x31*x41    
     1  +coeff( 28)    *x21*x31*x43    
     2  +coeff( 29)    *x24*x31*x41    
     3  +coeff( 30)    *x23    *x44    
     4  +coeff( 31)        *x32*x42    
     5  +coeff( 32)    *x21    *x42*x51
     6  +coeff( 33)    *x22    *x42*x51
     7  +coeff( 34)    *x21*x31*x41*x52
     8  +coeff( 35)    *x22        *x53
      pl_mnq_0_5  =pl_mnq_0_5  
     9  +coeff( 36)    *x21        *x54
     1  +coeff( 37)    *x24    *x42    
     2  +coeff( 38)    *x22*x31*x43    
     3  +coeff( 39)    *x23    *x42*x51
     4  +coeff( 40)    *x22*x31*x41*x52
c
      return
      end
      function x_mnq_0_7   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.1913084E+00/
      data xmin/
     1 -0.49978E-02,-0.50007E-01,-0.10393E+00,-0.49955E-01,-0.49981E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49997E-02, 0.49989E-01, 0.10393E+00, 0.49955E-01, 0.49998E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.19557291E+00, 0.66731340E+00, 0.43393472E+00, 0.24114610E-02,
     +  0.25504382E-03,-0.76770247E-02,-0.15125794E-01,-0.13282268E-01,
     + -0.26207703E+00, 0.51568258E-02,-0.17351294E-01,-0.27481902E-02,
     + -0.16609093E-01,-0.49370389E-01,-0.32684695E-01,-0.18000522E-02,
     + -0.17697983E-02,-0.25359470E-01,-0.11062209E-01, 0.13884832E+00,
     + -0.80235712E-02,-0.12537502E+00, 0.43685124E-02,-0.89011062E-03,
     + -0.28918257E-01,-0.10287879E+00, 0.14350176E+00,-0.31012387E-02,
     +  0.11061830E-01, 0.31002017E-01, 0.14044847E+00,-0.57555025E-03,
     + -0.21344153E-02,-0.29464532E-01,-0.40053714E-01, 0.60023386E-02,
     +  0.34337621E-01, 0.96663214E-01, 0.57449313E-02, 0.10268173E+00,
     +  0.55691255E-02, 0.66749252E-01, 0.22890751E+00, 0.11440796E+00,
     +  0.80471314E-01, 0.71784239E-02, 0.25489544E-01,-0.21535954E-02,
     +  0.91482028E-01,-0.68210252E-02, 0.25468281E+00,-0.10510105E-01,
     +  0.12994696E+00, 0.41666135E+00, 0.49517494E-01,-0.80779921E-02,
     +  0.25121045E+00, 0.84075006E-02, 0.47704627E-03,-0.64498946E-01,
     +  0.33677917E-01, 0.29501412E-01, 0.29101519E-01, 0.26250806E-01,
     +  0.28282328E-01, 0.13306387E-01, 0.42357665E-01, 0.79418652E-01,
     +  0.25506860E-01, 0.16816631E-01, 0.28600816E-01, 0.98841358E-02,
     + -0.25529790E+00,-0.17133161E+00,-0.46366893E-01,-0.39239731E-02,
     +  0.88068126E-02,-0.30473219E-02,-0.69044009E-02,-0.14605013E-02,
     +  0.10978600E-02, 0.13898327E-02, 0.73604472E-01, 0.30244039E-01,
     +  0.15187950E-01, 0.52101943E-02,-0.39118133E-01, 0.32647778E-02,
     + -0.12036609E-01,-0.14677010E-02,
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
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_mnq_0_7   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22            
     5  +coeff(  5)        *x32        
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21        *x51
      x_mnq_0_7   =x_mnq_0_7   
     9  +coeff(  9)                *x52
     1  +coeff( 10)*x11                
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x21*x32        
     4  +coeff( 13)    *x21*x31*x41    
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)    *x22        *x51
     7  +coeff( 16)        *x32    *x51
     8  +coeff( 17)        *x31*x41*x51
      x_mnq_0_7   =x_mnq_0_7   
     9  +coeff( 18)            *x42*x51
     1  +coeff( 19)    *x21        *x52
     2  +coeff( 20)                *x53
     3  +coeff( 21)    *x22*x32        
     4  +coeff( 22)    *x22    *x42    
     5  +coeff( 23)        *x31*x43    
     6  +coeff( 24)    *x23        *x51
     7  +coeff( 25)    *x21*x31*x41*x51
     8  +coeff( 26)    *x21    *x42*x51
      x_mnq_0_7   =x_mnq_0_7   
     9  +coeff( 27)    *x22        *x52
     1  +coeff( 28)        *x32    *x52
     2  +coeff( 29)        *x31*x41*x52
     3  +coeff( 30)            *x42*x52
     4  +coeff( 31)    *x21        *x53
     5  +coeff( 32)*x11    *x32        
     6  +coeff( 33)    *x23*x32        
     7  +coeff( 34)    *x23*x31*x41    
     8  +coeff( 35)    *x23    *x42    
      x_mnq_0_7   =x_mnq_0_7   
     9  +coeff( 36)    *x21*x31*x43    
     1  +coeff( 37)    *x22*x31*x41*x51
     2  +coeff( 38)    *x22    *x42*x51
     3  +coeff( 39)        *x32*x42*x51
     4  +coeff( 40)    *x23        *x52
     5  +coeff( 41)    *x21*x32    *x52
     6  +coeff( 42)    *x21*x31*x41*x52
     7  +coeff( 43)    *x21    *x42*x52
     8  +coeff( 44)    *x22        *x53
      x_mnq_0_7   =x_mnq_0_7   
     9  +coeff( 45)            *x42*x53
     1  +coeff( 46)    *x22*x33*x41    
     2  +coeff( 47)    *x22*x32*x42    
     3  +coeff( 48)        *x33*x43    
     4  +coeff( 49)    *x23*x31*x41*x51
     5  +coeff( 50)    *x21*x33*x41*x51
     6  +coeff( 51)    *x23    *x42*x51
     7  +coeff( 52)    *x21*x31*x43*x51
     8  +coeff( 53)    *x22*x31*x41*x52
      x_mnq_0_7   =x_mnq_0_7   
     9  +coeff( 54)    *x22    *x42*x52
     1  +coeff( 55)    *x23        *x53
     2  +coeff( 56)    *x21*x32    *x53
     3  +coeff( 57)    *x21    *x42*x53
     4  +coeff( 58)    *x23*x32*x42    
     5  +coeff( 59)    *x23*x31*x41*x52
     6  +coeff( 60)    *x23    *x42*x52
     7  +coeff( 61)    *x22*x31*x41*x53
     8  +coeff( 62)        *x33*x41*x53
      x_mnq_0_7   =x_mnq_0_7   
     9  +coeff( 63)        *x31*x43*x53
     1  +coeff( 64)    *x22*x33*x43    
     2  +coeff( 65)    *x21*x33*x41*x53
     3  +coeff( 66)*x11*x22*x31*x43    
     4  +coeff( 67)    *x22*x33*x43*x51
     5  +coeff( 68)    *x21*x33*x43*x52
     6  +coeff( 69)*x11*x23*x31*x43    
     7  +coeff( 70)*x11*x22*x33*x41*x51
     8  +coeff( 71)    *x23*x32*x42*x53
      x_mnq_0_7   =x_mnq_0_7   
     9  +coeff( 72)*x11*x23*x33*x41*x51
     1  +coeff( 73)    *x23*x33*x43*x52
     2  +coeff( 74)    *x22*x33*x43*x53
     3  +coeff( 75)    *x22*x31*x41    
     4  +coeff( 76)    *x21*x33*x41    
     5  +coeff( 77)    *x21*x32*x42    
     6  +coeff( 78)*x11*x21        *x51
     7  +coeff( 79)        *x33*x41*x51
     8  +coeff( 80)*x11            *x52
      x_mnq_0_7   =x_mnq_0_7   
     9  +coeff( 81)*x11*x21    *x42    
     1  +coeff( 82)        *x33*x41*x52
     2  +coeff( 83)    *x21*x31*x41*x53
     3  +coeff( 84)    *x23*x33*x41    
     4  +coeff( 85)    *x23*x31*x43    
     5  +coeff( 86)*x11*x23        *x51
     6  +coeff( 87)    *x22*x31*x43*x51
     7  +coeff( 88)*x11*x22        *x52
     8  +coeff( 89)    *x23*x32    *x52
      x_mnq_0_7   =x_mnq_0_7   
     9  +coeff( 90)*x11    *x31*x41*x52
c
      return
      end
      function t_mnq_0_7   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.4824348E-01/
      data xmin/
     1 -0.49978E-02,-0.50007E-01,-0.10393E+00,-0.49955E-01,-0.49981E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49997E-02, 0.49989E-01, 0.10393E+00, 0.49955E-01, 0.49998E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.50859969E-01, 0.43815468E-01, 0.16338752E+00,-0.10991500E-01,
     + -0.45600723E-03,-0.44252309E-02,-0.11492127E-01,-0.21434672E-01,
     + -0.11018211E+00, 0.84992615E-04,-0.11098385E-01,-0.13971358E-02,
     + -0.11876676E-01,-0.25409568E-01,-0.15403210E-01,-0.18622122E-02,
     + -0.72546569E-02,-0.14779290E-01, 0.10253688E-01, 0.66202708E-01,
     + -0.87947963E-03,-0.61938446E-01, 0.19104347E-02,-0.25604860E-03,
     +  0.49769916E-02,-0.18579870E-01,-0.60616355E-01, 0.85807227E-01,
     + -0.11243614E-02, 0.59385695E-02, 0.14220056E-01, 0.72114006E-01,
     + -0.43411460E-03, 0.18483492E-02, 0.10380502E-02,-0.19833036E-01,
     +  0.70403204E-02, 0.31992134E-01, 0.35173293E-01, 0.51045127E-01,
     +  0.43685913E-01, 0.91791287E-01, 0.62956132E-01, 0.40949702E-01,
     +  0.80433553E-02, 0.38875684E-02, 0.45705765E-01, 0.87447397E-01,
     +  0.77702895E-01, 0.17512675E+00, 0.28722642E-01, 0.16691374E-02,
     +  0.12364018E+00,-0.32917544E-03,-0.18496536E-01,-0.11389379E-02,
     +  0.83947191E-02,-0.12724563E-02, 0.40807254E-02, 0.56221201E-02,
     +  0.16953766E-01, 0.11913074E-01, 0.38301472E-02, 0.13373146E-01,
     +  0.14054030E-01,-0.50834141E-03, 0.10614480E-02, 0.21880087E-02,
     +  0.49312081E-01,-0.15387210E-02, 0.45807189E-02, 0.27725247E-02,
     + -0.50738931E-03,-0.26653666E-03, 0.48921898E-03,-0.12365002E-03,
     +  0.12442450E-01, 0.20977217E-03, 0.78115482E-02,-0.48353049E-03,
     +  0.59775898E-03,-0.34404182E-03,-0.64739399E-03, 0.15175999E-03,
     +  0.10249957E-01, 0.69414149E-02, 0.13790645E-02,-0.28323784E-03,
     + -0.51620364E-03, 0.37936030E-02,
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
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_mnq_0_7   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22            
     5  +coeff(  5)        *x32        
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21        *x51
      t_mnq_0_7   =t_mnq_0_7   
     9  +coeff(  9)                *x52
     1  +coeff( 10)*x11                
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x21*x32        
     4  +coeff( 13)    *x21*x31*x41    
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)    *x22        *x51
     7  +coeff( 16)        *x32    *x51
     8  +coeff( 17)        *x31*x41*x51
      t_mnq_0_7   =t_mnq_0_7   
     9  +coeff( 18)            *x42*x51
     1  +coeff( 19)    *x21        *x52
     2  +coeff( 20)                *x53
     3  +coeff( 21)    *x22*x32        
     4  +coeff( 22)    *x22    *x42    
     5  +coeff( 23)        *x32*x42    
     6  +coeff( 24)        *x31*x43    
     7  +coeff( 25)    *x23        *x51
     8  +coeff( 26)    *x21*x31*x41*x51
      t_mnq_0_7   =t_mnq_0_7   
     9  +coeff( 27)    *x21    *x42*x51
     1  +coeff( 28)    *x22        *x52
     2  +coeff( 29)        *x32    *x52
     3  +coeff( 30)        *x31*x41*x52
     4  +coeff( 31)            *x42*x52
     5  +coeff( 32)    *x21        *x53
     6  +coeff( 33)*x11    *x32        
     7  +coeff( 34)    *x23*x32        
     8  +coeff( 35)    *x23*x31*x41    
      t_mnq_0_7   =t_mnq_0_7   
     9  +coeff( 36)    *x23    *x42    
     1  +coeff( 37)    *x21*x31*x43    
     2  +coeff( 38)    *x22*x31*x41*x51
     3  +coeff( 39)    *x22    *x42*x51
     4  +coeff( 40)    *x23        *x52
     5  +coeff( 41)    *x21*x31*x41*x52
     6  +coeff( 42)    *x21    *x42*x52
     7  +coeff( 43)    *x22        *x53
     8  +coeff( 44)            *x42*x53
      t_mnq_0_7   =t_mnq_0_7   
     9  +coeff( 45)    *x22*x33*x41    
     1  +coeff( 46)        *x33*x43    
     2  +coeff( 47)    *x23*x31*x41*x51
     3  +coeff( 48)    *x23    *x42*x51
     4  +coeff( 49)    *x22*x31*x41*x52
     5  +coeff( 50)    *x22    *x42*x52
     6  +coeff( 51)    *x23        *x53
     7  +coeff( 52)    *x21*x32    *x53
     8  +coeff( 53)    *x21    *x42*x53
      t_mnq_0_7   =t_mnq_0_7   
     9  +coeff( 54)*x11*x21            
     1  +coeff( 55)    *x22*x31*x41    
     2  +coeff( 56)        *x33*x41    
     3  +coeff( 57)    *x21*x32*x42    
     4  +coeff( 58)*x11*x21        *x51
     5  +coeff( 59)        *x33*x41*x51
     6  +coeff( 60)        *x32*x42*x51
     7  +coeff( 61)        *x31*x41*x53
     8  +coeff( 62)    *x22*x32*x42    
      t_mnq_0_7   =t_mnq_0_7   
     9  +coeff( 63)    *x22*x31*x43    
     1  +coeff( 64)    *x21*x32*x42*x51
     2  +coeff( 65)    *x22*x32    *x52
     3  +coeff( 66)        *x33*x41*x52
     4  +coeff( 67)        *x32*x42*x52
     5  +coeff( 68)        *x31*x43*x52
     6  +coeff( 69)    *x21*x31*x41*x53
     7  +coeff( 70)    *x21*x32    *x51
     8  +coeff( 71)    *x21*x33*x41    
      t_mnq_0_7   =t_mnq_0_7   
     9  +coeff( 72)        *x31*x43*x51
     1  +coeff( 73)*x11            *x52
     2  +coeff( 74)*x11*x21*x32        
     3  +coeff( 75)*x11*x21*x31*x41    
     4  +coeff( 76)*x11    *x32    *x51
     5  +coeff( 77)    *x23*x32    *x51
     6  +coeff( 78)*x11    *x31*x41*x51
     7  +coeff( 79)    *x21*x33*x41*x51
     8  +coeff( 80)*x11        *x42*x51
      t_mnq_0_7   =t_mnq_0_7   
     9  +coeff( 81)*x11            *x53
     1  +coeff( 82)*x11            *x51
     2  +coeff( 83)*x11*x22            
     3  +coeff( 84)*x11    *x31*x41    
     4  +coeff( 85)    *x22*x32    *x51
     5  +coeff( 86)    *x21*x32    *x52
     6  +coeff( 87)        *x32    *x53
     7  +coeff( 88)*x11*x21    *x42    
     8  +coeff( 89)*x11*x22        *x51
      t_mnq_0_7   =t_mnq_0_7   
     9  +coeff( 90)    *x21*x31*x43*x51
c
      return
      end
      function y_mnq_0_7   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49978E-02,-0.50007E-01,-0.10393E+00,-0.49955E-01,-0.49981E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49997E-02, 0.49989E-01, 0.10393E+00, 0.49955E-01, 0.49998E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.99210210E-01, 0.64383346E+00,-0.75065969E-02,-0.21167666E-01,
     +  0.33368676E-02,-0.13090344E-01,-0.66539444E-01,-0.98342796E-04,
     +  0.67592831E-03,-0.18776996E-01,-0.95075302E-01, 0.21999937E-01,
     +  0.37231125E-01,-0.72446223E-02,-0.12711580E-03,-0.34508880E-01,
     +  0.16878017E-02, 0.74724113E-02, 0.13901552E-01,-0.39742137E-02,
     +  0.94003551E-01, 0.15181567E+00, 0.16628476E-01, 0.86107910E-01,
     +  0.16393738E-01, 0.99075466E-01, 0.86662225E-01,-0.19318791E-03,
     +  0.27466847E-02, 0.21129327E-02, 0.27181285E-01, 0.13140975E+00,
     +  0.20906422E+00, 0.14149250E-01, 0.13517876E+00, 0.28917615E-02,
     + -0.23457237E-01,-0.74617460E-01,-0.73577175E-02,-0.37194140E-01,
     + -0.10224234E-02, 0.21066152E-01, 0.45962889E-01,-0.12664346E-02,
     + -0.23972904E-01, 0.27120067E-03,-0.99745125E-01,-0.37487295E-01,
     + -0.11288488E+00,-0.13764982E+00, 0.48431032E-02,-0.77116005E-01,
     +  0.76612621E-03,-0.56468011E-02, 0.95470394E-04, 0.12703069E-02,
     + -0.54063235E-03, 0.57069678E-01, 0.17162476E-01,-0.65496488E-03,
     +  0.81451677E-01,-0.51403143E-02, 0.29243808E-01, 0.62553869E-02,
     + -0.52021956E-03,-0.25111646E-02,-0.22011368E-01,-0.43988526E-01,
     + -0.28657515E-01,-0.92390366E-02,-0.17117576E-02, 0.71174378E-03,
     +  0.47631152E-02, 0.54109758E-02,-0.27825573E-03,-0.93280161E-02,
     + -0.52889688E-02, 0.80349663E-03,-0.13554237E-02, 0.22534998E-02,
     + -0.22453098E-01, 0.66381167E-02, 0.19630454E-01,-0.12145649E-01,
     + -0.28502464E-02,-0.87070893E-02,-0.21915345E-02, 0.14064233E-02,
     + -0.28058223E-02, 0.24357426E-02,
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
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      y_mnq_0_7   =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21*x31        
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)            *x41*x51
     6  +coeff(  6)    *x22*x31        
     7  +coeff(  7)    *x22    *x41    
     8  +coeff(  8)        *x31*x42    
      y_mnq_0_7   =y_mnq_0_7   
     9  +coeff(  9)            *x43    
     1  +coeff( 10)    *x23*x31        
     2  +coeff( 11)    *x23    *x41    
     3  +coeff( 12)    *x21*x31*x42    
     4  +coeff( 13)    *x21    *x43    
     5  +coeff( 14)    *x22*x31    *x51
     6  +coeff( 15)        *x33    *x51
     7  +coeff( 16)    *x22    *x41*x51
     8  +coeff( 17)        *x32*x41*x51
      y_mnq_0_7   =y_mnq_0_7   
     9  +coeff( 18)        *x31*x42*x51
     1  +coeff( 19)            *x43*x51
     2  +coeff( 20)    *x21    *x41*x52
     3  +coeff( 21)    *x22*x31*x42    
     4  +coeff( 22)    *x22    *x43    
     5  +coeff( 23)    *x23*x31    *x51
     6  +coeff( 24)    *x23    *x41*x51
     7  +coeff( 25)    *x21*x32*x41*x51
     8  +coeff( 26)    *x21    *x43*x51
      y_mnq_0_7   =y_mnq_0_7   
     9  +coeff( 27)    *x22    *x41*x52
     1  +coeff( 28)        *x32*x41*x52
     2  +coeff( 29)        *x31*x42*x52
     3  +coeff( 30)            *x43*x52
     4  +coeff( 31)    *x23*x32*x41    
     5  +coeff( 32)    *x23*x31*x42    
     6  +coeff( 33)    *x23    *x43    
     7  +coeff( 34)    *x22*x32*x41*x51
     8  +coeff( 35)    *x22    *x43*x51
      y_mnq_0_7   =y_mnq_0_7   
     9  +coeff( 36)    *x21*x31*x42*x52
     1  +coeff( 37)    *x22    *x41*x53
     2  +coeff( 38)    *x23*x31*x42*x51
     3  +coeff( 39)    *x23    *x41*x53
     4  +coeff( 40)    *x21    *x41*x51
     5  +coeff( 41)        *x31    *x52
     6  +coeff( 42)    *x22*x32*x41    
     7  +coeff( 43)    *x21    *x41*x53
     8  +coeff( 44)*x11*x22    *x41    
      y_mnq_0_7   =y_mnq_0_7   
     9  +coeff( 45)    *x23    *x41*x52
     1  +coeff( 46)        *x31*x42*x53
     2  +coeff( 47)    *x23    *x43*x51
     3  +coeff( 48)    *x22*x32*x41*x52
     4  +coeff( 49)    *x22*x31*x42*x52
     5  +coeff( 50)    *x22    *x43*x52
     6  +coeff( 51)    *x23*x31    *x53
     7  +coeff( 52)    *x21    *x43*x53
     8  +coeff( 53)        *x31    *x51
      y_mnq_0_7   =y_mnq_0_7   
     9  +coeff( 54)    *x21*x31    *x51
     1  +coeff( 55)    *x21*x31    *x52
     2  +coeff( 56)            *x41*x53
     3  +coeff( 57)*x11*x21    *x41    
     4  +coeff( 58)    *x21*x31*x42*x51
     5  +coeff( 59)    *x22*x31    *x52
     6  +coeff( 60)*x11*x22*x31        
     7  +coeff( 61)    *x22*x31*x42*x51
     8  +coeff( 62)    *x21*x32*x41*x52
      y_mnq_0_7   =y_mnq_0_7   
     9  +coeff( 63)    *x21    *x43*x52
     1  +coeff( 64)            *x43*x53
     2  +coeff( 65)*x11*x23*x31        
     3  +coeff( 66)    *x22*x32*x43    
     4  +coeff( 67)    *x23*x32*x41*x51
     5  +coeff( 68)    *x21*x31*x42*x53
     6  +coeff( 69)    *x23*x31*x42*x53
     7  +coeff( 70)*x11*x23*x32*x43*x51
     8  +coeff( 71)            *x41*x52
      y_mnq_0_7   =y_mnq_0_7   
     9  +coeff( 72)    *x21*x33        
     1  +coeff( 73)    *x21*x32*x41    
     2  +coeff( 74)    *x21*x31    *x53
     3  +coeff( 75)*x11    *x32*x41    
     4  +coeff( 76)    *x23*x31    *x52
     5  +coeff( 77)    *x22*x31    *x53
     6  +coeff( 78)        *x33    *x53
     7  +coeff( 79)    *x23*x33    *x51
     8  +coeff( 80)    *x21*x33*x42*x51
      y_mnq_0_7   =y_mnq_0_7   
     9  +coeff( 81)    *x21*x32*x41*x53
     1  +coeff( 82)    *x22*x32*x43*x51
     2  +coeff( 83)    *x23*x31*x42*x52
     3  +coeff( 84)    *x22*x32*x41*x53
     4  +coeff( 85)        *x32*x43*x53
     5  +coeff( 86)    *x23*x32*x41*x53
     6  +coeff( 87)*x11*x23*x31*x42*x51
     7  +coeff( 88)*x11*x23    *x41*x53
     8  +coeff( 89)*x11*x23    *x43*x52
      y_mnq_0_7   =y_mnq_0_7   
     9  +coeff( 90)*x11*x21*x32*x43*x52
c
      return
      end
      function p_mnq_0_7   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49978E-02,-0.50007E-01,-0.10393E+00,-0.49955E-01,-0.49981E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49997E-02, 0.49989E-01, 0.10393E+00, 0.49955E-01, 0.49998E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.15253415E-02, 0.41766960E-01,-0.15951194E-02,-0.69361064E-02,
     + -0.11085926E-03,-0.52934390E-03,-0.76612541E-02,-0.42654861E-04,
     + -0.41638192E-01, 0.18508477E-04, 0.77638996E-03, 0.72369265E-03,
     + -0.44789370E-02,-0.31546369E-01,-0.66433544E-03,-0.21097267E-02,
     + -0.98370528E-02,-0.35558725E-04,-0.53863369E-01, 0.15608181E-02,
     +  0.13616685E-01, 0.25086084E-01,-0.10237482E-01,-0.53901553E-01,
     +  0.14512027E-02, 0.60114292E-02, 0.10134039E-01,-0.31766857E-02,
     + -0.16903304E-01, 0.34352345E-03,-0.19180779E-02,-0.51214074E-03,
     +  0.69074719E-02, 0.51698528E-01,-0.43549473E-03, 0.10660470E+00,
     +  0.11244658E-01, 0.40024091E-02, 0.39433204E-01, 0.88351615E-01,
     +  0.46734754E-02, 0.32609373E-01,-0.65586652E-03, 0.16594322E-02,
     +  0.75759073E-02, 0.22031525E-02, 0.20900510E-01,-0.23776503E-04,
     + -0.11307328E-02, 0.10274679E-01, 0.67112222E-01,-0.61891228E-03,
     +  0.13328832E+00, 0.11626208E-02, 0.70014234E-04,-0.39224638E-03,
     +  0.80776932E-02, 0.82447223E-01, 0.18589337E+00,-0.23599290E-02,
     +  0.12182831E-02, 0.35828082E-02, 0.23035951E-01, 0.71483895E-01,
     + -0.24031717E-02, 0.14717272E-03,-0.28826827E-02, 0.90393797E-02,
     +  0.52965904E-03, 0.20558001E-03,-0.15657356E-02,-0.52484931E-04,
     +  0.26839197E-01,-0.26906538E-03, 0.17563890E-02,-0.23977526E-01,
     + -0.19236090E-01, 0.14585737E-02,-0.14343185E-03, 0.24257686E-02,
     + -0.58491383E-03,-0.27029144E-03,-0.22988778E-01,-0.37040353E-01,
     +  0.74067764E-03, 0.70983375E-03,-0.46120398E-02, 0.75294324E-02,
     +  0.95814392E-02, 0.14756820E-02,
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
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_mnq_0_7   =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21*x31        
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)        *x31    *x51
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)    *x22*x31        
     8  +coeff(  8)        *x33        
      p_mnq_0_7   =p_mnq_0_7   
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)        *x32*x41    
     2  +coeff( 11)        *x31*x42    
     3  +coeff( 12)            *x43    
     4  +coeff( 13)    *x21*x31    *x51
     5  +coeff( 14)    *x21    *x41*x51
     6  +coeff( 15)        *x31    *x52
     7  +coeff( 16)            *x41*x52
     8  +coeff( 17)    *x23*x31        
      p_mnq_0_7   =p_mnq_0_7   
     9  +coeff( 18)*x11        *x41    
     1  +coeff( 19)    *x23    *x41    
     2  +coeff( 20)    *x21*x32*x41    
     3  +coeff( 21)    *x21*x31*x42    
     4  +coeff( 22)    *x21    *x43    
     5  +coeff( 23)    *x22*x31    *x51
     6  +coeff( 24)    *x22    *x41*x51
     7  +coeff( 25)        *x32*x41*x51
     8  +coeff( 26)        *x31*x42*x51
      p_mnq_0_7   =p_mnq_0_7   
     9  +coeff( 27)            *x43*x51
     1  +coeff( 28)    *x21*x31    *x52
     2  +coeff( 29)    *x21    *x41*x52
     3  +coeff( 30)        *x31    *x53
     4  +coeff( 31)            *x41*x53
     5  +coeff( 32)*x11*x21    *x41    
     6  +coeff( 33)    *x22*x32*x41    
     7  +coeff( 34)    *x22*x31*x42    
     8  +coeff( 35)        *x33*x42    
      p_mnq_0_7   =p_mnq_0_7   
     9  +coeff( 36)    *x22    *x43    
     1  +coeff( 37)    *x23    *x41*x51
     2  +coeff( 38)    *x21*x32*x41*x51
     3  +coeff( 39)    *x21*x31*x42*x51
     4  +coeff( 40)    *x21    *x43*x51
     5  +coeff( 41)    *x22*x31    *x52
     6  +coeff( 42)    *x22    *x41*x52
     7  +coeff( 43)        *x32*x41*x52
     8  +coeff( 44)        *x31*x42*x52
      p_mnq_0_7   =p_mnq_0_7   
     9  +coeff( 45)            *x43*x52
     1  +coeff( 46)    *x21*x31    *x53
     2  +coeff( 47)    *x21    *x41*x53
     3  +coeff( 48)    *x23*x33        
     4  +coeff( 49)*x11*x22    *x41    
     5  +coeff( 50)    *x23*x32*x41    
     6  +coeff( 51)    *x23*x31*x42    
     7  +coeff( 52)    *x21*x33*x42    
     8  +coeff( 53)    *x23    *x43    
      p_mnq_0_7   =p_mnq_0_7   
     9  +coeff( 54)    *x21*x32*x43    
     1  +coeff( 55)*x11*x21*x31    *x51
     2  +coeff( 56)*x11*x21    *x41*x51
     3  +coeff( 57)    *x22*x32*x41*x51
     4  +coeff( 58)    *x22*x31*x42*x51
     5  +coeff( 59)    *x22    *x43*x51
     6  +coeff( 60)    *x23*x31    *x52
     7  +coeff( 61)    *x21*x33    *x52
     8  +coeff( 62)    *x21*x32*x41*x52
      p_mnq_0_7   =p_mnq_0_7   
     9  +coeff( 63)    *x21*x31*x42*x52
     1  +coeff( 64)    *x21    *x43*x52
     2  +coeff( 65)    *x22*x31    *x53
     3  +coeff( 66)        *x33    *x53
     4  +coeff( 67)    *x22    *x41*x53
     5  +coeff( 68)            *x43*x53
     6  +coeff( 69)*x11*x21*x32*x41    
     7  +coeff( 70)*x11*x21*x31*x42    
     8  +coeff( 71)    *x22*x32*x43    
      p_mnq_0_7   =p_mnq_0_7   
     9  +coeff( 72)*x11        *x43*x51
     1  +coeff( 73)    *x23    *x43*x51
     2  +coeff( 74)*x11*x21*x31    *x52
     3  +coeff( 75)    *x22*x32*x41*x52
     4  +coeff( 76)    *x22*x31*x42*x52
     5  +coeff( 77)    *x22    *x43*x52
     6  +coeff( 78)        *x32*x43*x52
     7  +coeff( 79)*x11    *x31    *x53
     8  +coeff( 80)    *x23*x31    *x53
      p_mnq_0_7   =p_mnq_0_7   
     9  +coeff( 81)    *x21*x33    *x53
     1  +coeff( 82)*x11        *x41*x53
     2  +coeff( 83)    *x21*x31*x42*x53
     3  +coeff( 84)    *x21    *x43*x53
     4  +coeff( 85)*x11*x22*x32*x41    
     5  +coeff( 86)*x11*x22    *x43    
     6  +coeff( 87)    *x23*x32*x43    
     7  +coeff( 88)    *x23*x31*x42*x52
     8  +coeff( 89)    *x23    *x43*x52
      p_mnq_0_7   =p_mnq_0_7   
     9  +coeff( 90)    *x22*x33    *x53
c
      return
      end
      function pl_mnq_0_7  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 41)
      data ncoeff/ 40/
      data avdat/  0.2556555E-01/
      data xmin/
     1 -0.49978E-02,-0.50007E-01,-0.10393E+00,-0.49955E-01,-0.49981E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49997E-02, 0.49989E-01, 0.10393E+00, 0.49955E-01, 0.49998E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.25597509E-01,-0.18951036E+00,-0.68876930E-01,-0.15313345E-01,
     + -0.23686649E-01,-0.16032000E-02,-0.17327195E-02,-0.21123473E-01,
     +  0.14644511E-01,-0.77259936E-02, 0.84858462E-02,-0.64011216E-02,
     +  0.15265765E-01, 0.41793892E-02, 0.20716771E-01, 0.75353645E-02,
     + -0.18370744E-01,-0.12320419E-01, 0.64224415E-02, 0.19786802E-03,
     + -0.25476848E-02, 0.29235512E-01, 0.45592440E-02,-0.29836021E-01,
     + -0.13401114E-01, 0.16391000E-02,-0.87263643E-05, 0.10441703E-01,
     +  0.13853201E-01,-0.33266463E-02, 0.18051172E-01, 0.62577045E-02,
     + -0.30324379E-01,-0.80030123E-02, 0.27875910E-02,-0.52751177E-02,
     +  0.51447279E-02,-0.22655188E-02, 0.30087100E-02, 0.97198542E-02,
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
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
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
      pl_mnq_0_7  =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)            *x42    
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)    *x24            
     7  +coeff(  7)*x11                
     8  +coeff(  8)    *x22            
      pl_mnq_0_7  =pl_mnq_0_7  
     9  +coeff(  9)                *x53
     1  +coeff( 10)    *x24        *x51
     2  +coeff( 11)    *x22        *x51
     3  +coeff( 12)    *x21        *x53
     4  +coeff( 13)                *x52
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)    *x22    *x42    
     7  +coeff( 16)    *x22        *x52
     8  +coeff( 17)                *x54
      pl_mnq_0_7  =pl_mnq_0_7  
     9  +coeff( 18)    *x23        *x52
     1  +coeff( 19)    *x23        *x54
     2  +coeff( 20)        *x31*x41    
     3  +coeff( 21)    *x23            
     4  +coeff( 22)    *x21        *x52
     5  +coeff( 23)    *x22*x31*x41    
     6  +coeff( 24)    *x21        *x54
     7  +coeff( 25)    *x22    *x44    
     8  +coeff( 26)            *x42*x51
      pl_mnq_0_7  =pl_mnq_0_7  
     9  +coeff( 27)        *x31*x43    
     1  +coeff( 28)    *x21    *x42*x51
     2  +coeff( 29)    *x23    *x42    
     3  +coeff( 30)    *x21    *x44    
     4  +coeff( 31)    *x22    *x42*x51
     5  +coeff( 32)    *x21    *x42*x52
     6  +coeff( 33)    *x22        *x53
     7  +coeff( 34)    *x22*x31*x43    
     8  +coeff( 35)            *x44*x52
      pl_mnq_0_7  =pl_mnq_0_7  
     9  +coeff( 36)    *x22        *x54
     1  +coeff( 37)    *x23*x33*x41    
     2  +coeff( 38)    *x21*x33*x43    
     3  +coeff( 39)    *x22*x33*x41*x51
     4  +coeff( 40)    *x24        *x53
c
      return
      end
      function x_mnq_0_9   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.7938691E-01/
      data xmin/
     1 -0.49978E-02,-0.50007E-01,-0.10393E+00,-0.46630E-01,-0.42140E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49966E-02, 0.49989E-01, 0.10393E+00, 0.46630E-01, 0.49976E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.14652123E+00, 0.75314659E+00, 0.75889868E+00, 0.18798964E-01,
     + -0.15112488E-02,-0.13128483E-01,-0.37018321E-01, 0.62534581E-02,
     + -0.36694697E+00, 0.44992659E-02, 0.77216686E-02, 0.25422242E-02,
     + -0.75331405E-01,-0.54104552E-01, 0.21695502E-02,-0.15265602E-01,
     + -0.48549928E-01,-0.57839602E-01, 0.14406417E+00,-0.11381907E-02,
     + -0.45489292E-02,-0.15823531E+00,-0.35959715E+00,-0.55705854E-02,
     + -0.52719016E-03,-0.83399683E-01,-0.86583253E-02,-0.17447324E+00,
     + -0.38115588E+00, 0.46419583E-01,-0.67297919E-02,-0.33429366E-01,
     + -0.69605656E-01, 0.15778928E+00,-0.11650857E-03,-0.10246017E-01,
     + -0.12517112E+00, 0.45670732E-02, 0.21028951E-03,-0.32995790E+00,
     + -0.54085599E-02, 0.37735982E-02, 0.75975019E-02,-0.19235471E+00,
     + -0.33915765E-02, 0.13873336E+00, 0.45730677E-01, 0.24542400E+00,
     +  0.18749240E+00, 0.20403078E+00, 0.28662292E-01,-0.16787385E-02,
     +  0.13996227E+00, 0.40036369E-01, 0.59075278E+00,-0.30923884E-02,
     +  0.20913944E+00, 0.10679877E+01, 0.41144967E-01, 0.43276094E-01,
     +  0.13349664E+00, 0.56912643E+00, 0.27605737E-01,-0.16045759E-01,
     + -0.14045230E-01, 0.51870111E-01, 0.71900912E-01, 0.14525181E-01,
     + -0.29896980E-01, 0.29432154E-02,-0.35455027E-02,-0.10260959E+00,
     + -0.52185259E-02, 0.48723500E-01, 0.40788320E-02, 0.31987843E-02,
     +  0.65020330E-01, 0.56932000E-02, 0.48676711E-01, 0.84575579E-01,
     +  0.42380080E-01, 0.91771729E-03,-0.45582810E-02,-0.37989961E-02,
     +  0.39237563E-01, 0.24283614E-02,-0.30945661E-01,-0.26688389E-01,
     +  0.58696298E-02,-0.31358924E-01,
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
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_mnq_0_9   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22            
     5  +coeff(  5)        *x32        
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21        *x51
      x_mnq_0_9   =x_mnq_0_9   
     9  +coeff(  9)                *x52
     1  +coeff( 10)*x11                
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x21*x32        
     4  +coeff( 13)    *x21    *x42    
     5  +coeff( 14)    *x22        *x51
     6  +coeff( 15)        *x32    *x51
     7  +coeff( 16)        *x31*x41*x51
     8  +coeff( 17)            *x42*x51
      x_mnq_0_9   =x_mnq_0_9   
     9  +coeff( 18)    *x21        *x52
     1  +coeff( 19)                *x53
     2  +coeff( 20)*x11*x21            
     3  +coeff( 21)    *x22*x32        
     4  +coeff( 22)    *x22*x31*x41    
     5  +coeff( 23)    *x22    *x42    
     6  +coeff( 24)        *x31*x43    
     7  +coeff( 25)*x11            *x51
     8  +coeff( 26)    *x23        *x51
      x_mnq_0_9   =x_mnq_0_9   
     9  +coeff( 27)    *x21*x32    *x51
     1  +coeff( 28)    *x21*x31*x41*x51
     2  +coeff( 29)    *x21    *x42*x51
     3  +coeff( 30)    *x22        *x52
     4  +coeff( 31)        *x32    *x52
     5  +coeff( 32)        *x31*x41*x52
     6  +coeff( 33)            *x42*x52
     7  +coeff( 34)    *x21        *x53
     8  +coeff( 35)*x11*x22            
      x_mnq_0_9   =x_mnq_0_9   
     9  +coeff( 36)    *x23*x32        
     1  +coeff( 37)    *x23*x31*x41    
     2  +coeff( 38)    *x21*x33*x41    
     3  +coeff( 39)*x11        *x42    
     4  +coeff( 40)    *x23    *x42    
     5  +coeff( 41)    *x21*x32*x42    
     6  +coeff( 42)    *x21*x31*x43    
     7  +coeff( 43)        *x33*x41*x51
     8  +coeff( 44)    *x22    *x42*x51
      x_mnq_0_9   =x_mnq_0_9   
     9  +coeff( 45)        *x31*x43*x51
     1  +coeff( 46)    *x23        *x52
     2  +coeff( 47)    *x21*x31*x41*x52
     3  +coeff( 48)    *x21    *x42*x52
     4  +coeff( 49)    *x22        *x53
     5  +coeff( 50)            *x42*x53
     6  +coeff( 51)    *x22*x33*x41    
     7  +coeff( 52)    *x22*x32*x42    
     8  +coeff( 53)    *x23*x31*x41*x51
      x_mnq_0_9   =x_mnq_0_9   
     9  +coeff( 54)    *x21*x33*x41*x51
     1  +coeff( 55)    *x23    *x42*x51
     2  +coeff( 56)*x11*x21        *x52
     3  +coeff( 57)    *x22*x31*x41*x52
     4  +coeff( 58)    *x22    *x42*x52
     5  +coeff( 59)        *x31*x43*x52
     6  +coeff( 60)    *x23        *x53
     7  +coeff( 61)    *x21*x31*x41*x53
     8  +coeff( 62)    *x21    *x42*x53
      x_mnq_0_9   =x_mnq_0_9   
     9  +coeff( 63)    *x23*x33*x41    
     1  +coeff( 64)        *x33*x41*x53
     2  +coeff( 65)*x11*x23    *x42    
     3  +coeff( 66)    *x23*x33*x41*x51
     4  +coeff( 67)    *x22*x33*x41*x52
     5  +coeff( 68)*x11*x23    *x42*x51
     6  +coeff( 69)    *x21*x31*x41    
     7  +coeff( 70)        *x32*x42    
     8  +coeff( 71)*x11*x21        *x51
      x_mnq_0_9   =x_mnq_0_9   
     9  +coeff( 72)    *x22*x31*x41*x51
     1  +coeff( 73)        *x32*x42*x51
     2  +coeff( 74)        *x31*x41*x53
     3  +coeff( 75)*x11*x23            
     4  +coeff( 76)*x11*x21    *x42    
     5  +coeff( 77)    *x22*x31*x43    
     6  +coeff( 78)        *x33*x43    
     7  +coeff( 79)    *x23*x32    *x51
     8  +coeff( 80)    *x21*x31*x43*x51
      x_mnq_0_9   =x_mnq_0_9   
     9  +coeff( 81)    *x22*x32    *x52
     1  +coeff( 82)*x11*x22*x32        
     2  +coeff( 83)*x11*x22*x31*x41    
     3  +coeff( 84)*x11*x21*x31*x41*x51
     4  +coeff( 85)    *x22*x33*x41*x51
     5  +coeff( 86)*x11*x22        *x52
     6  +coeff( 87)    *x23*x32    *x52
     7  +coeff( 88)    *x23*x31*x41*x52
     8  +coeff( 89)*x11*x21        *x53
      x_mnq_0_9   =x_mnq_0_9   
     9  +coeff( 90)    *x22*x32    *x53
c
      return
      end
      function t_mnq_0_9   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.5018460E-02/
      data xmin/
     1 -0.49978E-02,-0.50007E-01,-0.10393E+00,-0.46630E-01,-0.42140E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49966E-02, 0.49989E-01, 0.10393E+00, 0.46630E-01, 0.49976E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.13896543E-01, 0.43563645E-01, 0.21872087E+00, 0.98897452E-02,
     + -0.43831041E-03,-0.44103060E-02,-0.15588412E-01, 0.12829172E-01,
     + -0.99437833E-01, 0.26917836E-04, 0.95780818E-02, 0.11782527E-02,
     + -0.41056763E-01, 0.35063177E-02, 0.76037872E-03,-0.84083993E-02,
     + -0.27999571E-01,-0.81669940E-02, 0.40391751E-01,-0.10788791E-03,
     + -0.23063386E-02,-0.17287149E+00, 0.81947143E-03,-0.36325287E-02,
     + -0.14687581E-01,-0.40873019E-02,-0.75124942E-01,-0.21252948E+00,
     +  0.19563848E-01,-0.22655784E-02,-0.18612027E-01,-0.52130610E-01,
     +  0.46037629E-01,-0.11916151E-02,-0.51463516E-02,-0.47002919E-01,
     +  0.15672725E-02,-0.31502952E-03,-0.15220161E+00,-0.27430074E-02,
     + -0.60657423E-01, 0.12198174E-02,-0.20073059E+00,-0.20830946E-02,
     +  0.38414482E-01,-0.29427859E-02,-0.37874819E-02, 0.16369317E-01,
     +  0.55217624E-01,-0.13745152E-02, 0.15502636E-01, 0.72394803E-01,
     +  0.37167794E-02,-0.42589610E-02, 0.17833952E-01, 0.32668295E-02,
     +  0.85144741E-02, 0.38365863E-01, 0.30239872E-02, 0.80278583E-01,
     + -0.24214685E-02, 0.14625229E-01,-0.71486941E-03, 0.48806842E-01,
     +  0.19356775E+00,-0.12550072E-01,-0.60269311E-01,-0.40719933E-02,
     + -0.18633378E-02,-0.17543766E-02, 0.12509026E+00, 0.31030861E+00,
     + -0.26996341E-03,-0.29298197E-03,-0.43920521E-02,-0.83106325E-03,
     + -0.10916427E-02,-0.42062008E-03, 0.32277692E-01,-0.13932279E-03,
     +  0.70405793E-02, 0.16465982E-01, 0.84818504E-03,-0.44294837E-03,
     + -0.82628925E-04,-0.80185197E-03,-0.31079605E-03,-0.17038418E-03,
     + -0.13965775E-02,-0.17463972E-03,
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
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_mnq_0_9   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22            
     5  +coeff(  5)        *x32        
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21        *x51
      t_mnq_0_9   =t_mnq_0_9   
     9  +coeff(  9)                *x52
     1  +coeff( 10)*x11                
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x21*x32        
     4  +coeff( 13)    *x21    *x42    
     5  +coeff( 14)    *x22        *x51
     6  +coeff( 15)        *x32    *x51
     7  +coeff( 16)        *x31*x41*x51
     8  +coeff( 17)            *x42*x51
      t_mnq_0_9   =t_mnq_0_9   
     9  +coeff( 18)    *x21        *x52
     1  +coeff( 19)                *x53
     2  +coeff( 20)*x11*x21            
     3  +coeff( 21)    *x22*x32        
     4  +coeff( 22)    *x22    *x42    
     5  +coeff( 23)        *x32*x42    
     6  +coeff( 24)        *x31*x43    
     7  +coeff( 25)    *x23        *x51
     8  +coeff( 26)    *x21*x32    *x51
      t_mnq_0_9   =t_mnq_0_9   
     9  +coeff( 27)    *x21*x31*x41*x51
     1  +coeff( 28)    *x21    *x42*x51
     2  +coeff( 29)    *x22        *x52
     3  +coeff( 30)        *x32    *x52
     4  +coeff( 31)        *x31*x41*x52
     5  +coeff( 32)            *x42*x52
     6  +coeff( 33)    *x21        *x53
     7  +coeff( 34)*x11*x22            
     8  +coeff( 35)    *x23*x32        
      t_mnq_0_9   =t_mnq_0_9   
     9  +coeff( 36)    *x23*x31*x41    
     1  +coeff( 37)    *x21*x33*x41    
     2  +coeff( 38)*x11        *x42    
     3  +coeff( 39)    *x23    *x42    
     4  +coeff( 40)    *x21*x31*x43    
     5  +coeff( 41)    *x22*x31*x41*x51
     6  +coeff( 42)        *x33*x41*x51
     7  +coeff( 43)    *x22    *x42*x51
     8  +coeff( 44)        *x32*x42*x51
      t_mnq_0_9   =t_mnq_0_9   
     9  +coeff( 45)    *x23        *x52
     1  +coeff( 46)    *x21*x32    *x52
     2  +coeff( 47)    *x21*x31*x41*x52
     3  +coeff( 48)    *x21    *x42*x52
     4  +coeff( 49)    *x22        *x53
     5  +coeff( 50)        *x32    *x53
     6  +coeff( 51)        *x31*x41*x53
     7  +coeff( 52)            *x42*x53
     8  +coeff( 53)    *x22*x33*x41    
      t_mnq_0_9   =t_mnq_0_9   
     9  +coeff( 54)    *x22*x32*x42    
     1  +coeff( 55)    *x22*x31*x43    
     2  +coeff( 56)        *x33*x43    
     3  +coeff( 57)    *x23*x32    *x51
     4  +coeff( 58)    *x23*x31*x41*x51
     5  +coeff( 59)    *x21*x33*x41*x51
     6  +coeff( 60)    *x22*x31*x41*x52
     7  +coeff( 61)        *x33*x41*x52
     8  +coeff( 62)    *x23        *x53
      t_mnq_0_9   =t_mnq_0_9   
     9  +coeff( 63)    *x21*x32    *x53
     1  +coeff( 64)    *x21*x31*x41*x53
     2  +coeff( 65)    *x21    *x42*x53
     3  +coeff( 66)    *x21*x31*x41    
     4  +coeff( 67)    *x22*x31*x41    
     5  +coeff( 68)    *x21*x32*x42    
     6  +coeff( 69)*x11*x21        *x51
     7  +coeff( 70)        *x31*x43*x51
     8  +coeff( 71)    *x23    *x42*x51
      t_mnq_0_9   =t_mnq_0_9   
     9  +coeff( 72)    *x22    *x42*x52
     1  +coeff( 73)*x11            *x51
     2  +coeff( 74)*x11    *x31*x41    
     3  +coeff( 75)    *x22*x32    *x51
     4  +coeff( 76)*x11            *x52
     5  +coeff( 77)*x11*x21    *x42    
     6  +coeff( 78)*x11        *x42*x51
     7  +coeff( 79)    *x21*x31*x43*x51
     8  +coeff( 80)*x11*x21        *x52
      t_mnq_0_9   =t_mnq_0_9   
     9  +coeff( 81)    *x22*x32    *x52
     1  +coeff( 82)        *x31*x43*x52
     2  +coeff( 83)*x11            *x53
     3  +coeff( 84)        *x33*x41    
     4  +coeff( 85)*x11    *x32        
     5  +coeff( 86)*x11*x23            
     6  +coeff( 87)*x11*x21*x32        
     7  +coeff( 88)*x11*x21*x31*x41    
     8  +coeff( 89)*x11*x22        *x51
      t_mnq_0_9   =t_mnq_0_9   
     9  +coeff( 90)*x11    *x32    *x51
c
      return
      end
      function y_mnq_0_9   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49978E-02,-0.50007E-01,-0.10393E+00,-0.46630E-01,-0.42140E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49966E-02, 0.49989E-01, 0.10393E+00, 0.46630E-01, 0.49976E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.96610606E-01, 0.68422359E+00,-0.97147953E-02,-0.32781180E-01,
     +  0.19178600E-02, 0.92713293E-02,-0.42187676E-01,-0.24942885E+00,
     + -0.57337811E-03,-0.16164324E-02,-0.31161405E-01,-0.22002220E+00,
     + -0.10709154E-01,-0.67243159E-01,-0.56275833E-01,-0.32266048E+00,
     +  0.54646190E-02, 0.67141339E-01, 0.13960643E+00,-0.36040357E+00,
     +  0.75550303E-01,-0.99491604E-01, 0.17825951E-02,-0.10844850E-02,
     +  0.60624961E-01, 0.35152701E+00,-0.23563574E-02, 0.71041226E+00,
     +  0.30972978E-01, 0.17304000E+00, 0.33383021E+00, 0.70723903E+00,
     +  0.33677357E+00, 0.64829350E-01, 0.14985597E+00, 0.14544247E+00,
     + -0.18064936E-02, 0.48595202E+00, 0.91480696E+00, 0.94653063E-01,
     + -0.56320713E-02, 0.12908082E+01, 0.28851242E-02,-0.80759294E-01,
     +  0.22703362E+00, 0.48523378E+00,-0.10393612E-01, 0.70836605E-03,
     + -0.56428693E-01, 0.16387791E-02,-0.88243090E-01,-0.20113476E+00,
     + -0.12190867E-01,-0.23872162E+00,-0.39477208E+00, 0.89202570E-02,
     +  0.48089992E-01,-0.48634339E-01,-0.20811889E-01,-0.51163277E-02,
     +  0.12566728E-02,-0.20624808E-03,-0.50651561E-01,-0.36994476E-03,
     +  0.38790874E-01, 0.22443483E-01, 0.50860487E-01, 0.10421725E-01,
     +  0.19200854E-01, 0.99810548E-01, 0.60321164E+00,-0.42297882E+00,
     + -0.27727058E-01,-0.28174698E-01, 0.67059780E-02,-0.12147587E-01,
     +  0.42690644E-02,-0.82098186E-03,-0.44768551E-03, 0.50602932E-01,
     + -0.97243412E-03, 0.24549343E-03, 0.32902875E-02,-0.15620287E-01,
     +  0.23056623E-01,-0.57978305E-03,-0.10549921E-02,-0.15332982E+00,
     + -0.68901861E+00,-0.21418786E+00,
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
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      y_mnq_0_9   =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21*x31        
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)        *x31    *x51
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)    *x22*x31        
     8  +coeff(  8)    *x22    *x41    
      y_mnq_0_9   =y_mnq_0_9   
     9  +coeff(  9)        *x32*x41    
     1  +coeff( 10)            *x43    
     2  +coeff( 11)    *x21*x31    *x51
     3  +coeff( 12)    *x21    *x41*x51
     4  +coeff( 13)        *x31    *x52
     5  +coeff( 14)            *x41*x52
     6  +coeff( 15)    *x23*x31        
     7  +coeff( 16)    *x23    *x41    
     8  +coeff( 17)    *x21*x32*x41    
      y_mnq_0_9   =y_mnq_0_9   
     9  +coeff( 18)    *x21*x31*x42    
     1  +coeff( 19)    *x21    *x43    
     2  +coeff( 20)    *x22    *x41*x51
     3  +coeff( 21)            *x43*x51
     4  +coeff( 22)    *x21    *x41*x52
     5  +coeff( 23)    *x22*x33        
     6  +coeff( 24)*x11*x21    *x41    
     7  +coeff( 25)    *x22*x32*x41    
     8  +coeff( 26)    *x22*x31*x42    
      y_mnq_0_9   =y_mnq_0_9   
     9  +coeff( 27)        *x33*x42    
     1  +coeff( 28)    *x22    *x43    
     2  +coeff( 29)    *x23*x31    *x51
     3  +coeff( 30)    *x23    *x41*x51
     4  +coeff( 31)    *x21*x31*x42*x51
     5  +coeff( 32)    *x21    *x43*x51
     6  +coeff( 33)    *x22    *x41*x52
     7  +coeff( 34)        *x31*x42*x52
     8  +coeff( 35)            *x43*x52
      y_mnq_0_9   =y_mnq_0_9   
     9  +coeff( 36)    *x21    *x41*x53
     1  +coeff( 37)*x11*x22    *x41    
     2  +coeff( 38)    *x23*x31*x42    
     3  +coeff( 39)    *x23    *x43    
     4  +coeff( 40)    *x22*x32*x41*x51
     5  +coeff( 41)        *x33*x42*x51
     6  +coeff( 42)    *x22    *x43*x51
     7  +coeff( 43)        *x32*x43*x51
     8  +coeff( 44)    *x23    *x41*x52
      y_mnq_0_9   =y_mnq_0_9   
     9  +coeff( 45)    *x21*x31*x42*x52
     1  +coeff( 46)    *x21    *x43*x52
     2  +coeff( 47)    *x22*x31    *x53
     3  +coeff( 48)        *x33    *x53
     4  +coeff( 49)    *x22    *x41*x53
     5  +coeff( 50)        *x32*x41*x53
     6  +coeff( 51)    *x23*x32*x41*x51
     7  +coeff( 52)    *x23*x31*x42*x51
     8  +coeff( 53)    *x21*x33*x42*x51
      y_mnq_0_9   =y_mnq_0_9   
     9  +coeff( 54)    *x23    *x43*x51
     1  +coeff( 55)    *x22*x31*x42*x52
     2  +coeff( 56)    *x23*x31    *x53
     3  +coeff( 57)    *x23    *x41*x53
     4  +coeff( 58)    *x21*x32*x41*x53
     5  +coeff( 59)    *x21*x33*x42*x52
     6  +coeff( 60)    *x23*x33*x42*x52
     7  +coeff( 61)        *x31*x42    
     8  +coeff( 62)    *x21*x33        
      y_mnq_0_9   =y_mnq_0_9   
     9  +coeff( 63)    *x22*x31    *x51
     1  +coeff( 64)        *x33    *x51
     2  +coeff( 65)        *x31*x42*x51
     3  +coeff( 66)            *x41*x53
     4  +coeff( 67)    *x22*x31    *x52
     5  +coeff( 68)        *x32*x41*x52
     6  +coeff( 69)    *x21*x31    *x53
     7  +coeff( 70)    *x23*x32*x41    
     8  +coeff( 71)    *x22*x31*x42*x51
      y_mnq_0_9   =y_mnq_0_9   
     9  +coeff( 72)    *x21    *x43*x53
     1  +coeff( 73)    *x22*x33*x42*x52
     2  +coeff( 74)    *x23*x32*x41*x53
     3  +coeff( 75)    *x23    *x43*x53
     4  +coeff( 76)    *x21*x31    *x52
     5  +coeff( 77)        *x31    *x53
     6  +coeff( 78)    *x21*x33    *x51
     7  +coeff( 79)*x11        *x41*x51
     8  +coeff( 80)    *x21*x32*x41*x51
      y_mnq_0_9   =y_mnq_0_9   
     9  +coeff( 81)*x11*x22*x31        
     1  +coeff( 82)*x11    *x33        
     2  +coeff( 83)    *x23*x33        
     3  +coeff( 84)    *x23*x31    *x52
     4  +coeff( 85)    *x21*x32*x41*x52
     5  +coeff( 86)*x11*x23*x31        
     6  +coeff( 87)*x11*x21*x32*x41    
     7  +coeff( 88)    *x22*x32*x41*x52
     8  +coeff( 89)    *x22    *x43*x52
      y_mnq_0_9   =y_mnq_0_9   
     9  +coeff( 90)    *x21*x31*x42*x53
c
      return
      end
      function p_mnq_0_9   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49978E-02,-0.50007E-01,-0.10393E+00,-0.46630E-01,-0.42140E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49966E-02, 0.49989E-01, 0.10393E+00, 0.46630E-01, 0.49976E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.19688995E-02, 0.33939544E-01,-0.32814196E-02,-0.17340003E-01,
     + -0.99362840E-03,-0.50633843E-03,-0.11548369E-01,-0.59460342E-03,
     + -0.90850011E-01, 0.14128631E-02, 0.24131364E-02,-0.10079984E-01,
     + -0.10986417E+00,-0.53276881E-02,-0.35716940E-01,-0.19476755E-03,
     + -0.10703862E-01,-0.57261408E-03,-0.34186873E-03,-0.10686316E+00,
     +  0.57779825E-02, 0.37536893E-01, 0.72331928E-01, 0.60236093E-03,
     + -0.15464872E+00, 0.47714505E-02, 0.27480057E-01, 0.51793627E-01,
     +  0.90902010E-02,-0.60899094E-01, 0.63565020E-02, 0.63085193E-02,
     + -0.54333272E-03,-0.79650083E-03, 0.21928800E-02, 0.12113445E+00,
     +  0.18292114E-02, 0.26241803E+00,-0.63841138E-03,-0.21453491E-03,
     +  0.37393908E-02,-0.34131948E-02, 0.28951000E-01, 0.13338690E+00,
     +  0.32170048E+00, 0.39693327E-02,-0.70374081E-03, 0.46815153E-01,
     + -0.54098666E-02, 0.29857142E-01, 0.77185079E-01, 0.32690294E-01,
     + -0.18474324E-03, 0.32120905E-03,-0.83194289E-03, 0.11882890E-03,
     +  0.39384764E-03, 0.12866546E+00,-0.16483251E-02, 0.28575066E-03,
     +  0.30516300E+00,-0.22474579E-02,-0.10525758E-01,-0.10491646E-02,
     + -0.21114895E-01, 0.18386911E+00, 0.49590179E+00,-0.62589103E-03,
     +  0.19219851E-02,-0.10888753E-01, 0.49138756E-03,-0.30244463E-02,
     + -0.18017968E-01, 0.63848399E-01, 0.21094376E+00,-0.33459249E-02,
     + -0.15959035E-02,-0.34391757E-02,-0.14582603E-02, 0.19056065E-01,
     +  0.73433638E-03,-0.16835381E-02, 0.14705681E-02,-0.59912410E-02,
     +  0.48520134E-03, 0.21567384E-02, 0.42684548E-03, 0.11083988E-02,
     + -0.36381150E-02,-0.73423390E-02,
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
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_mnq_0_9   =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21*x31        
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)        *x31    *x51
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)    *x22*x31        
     8  +coeff(  8)        *x33        
      p_mnq_0_9   =p_mnq_0_9   
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)        *x32*x41    
     2  +coeff( 11)            *x43    
     3  +coeff( 12)    *x21*x31    *x51
     4  +coeff( 13)    *x21    *x41*x51
     5  +coeff( 14)        *x31    *x52
     6  +coeff( 15)            *x41*x52
     7  +coeff( 16)*x11    *x31        
     8  +coeff( 17)    *x23*x31        
      p_mnq_0_9   =p_mnq_0_9   
     9  +coeff( 18)    *x21*x33        
     1  +coeff( 19)*x11        *x41    
     2  +coeff( 20)    *x23    *x41    
     3  +coeff( 21)    *x21*x32*x41    
     4  +coeff( 22)    *x21*x31*x42    
     5  +coeff( 23)    *x21    *x43    
     6  +coeff( 24)        *x33    *x51
     7  +coeff( 25)    *x22    *x41*x51
     8  +coeff( 26)        *x32*x41*x51
      p_mnq_0_9   =p_mnq_0_9   
     9  +coeff( 27)        *x31*x42*x51
     1  +coeff( 28)            *x43*x51
     2  +coeff( 29)    *x21*x31    *x52
     3  +coeff( 30)    *x21    *x41*x52
     4  +coeff( 31)        *x31    *x53
     5  +coeff( 32)            *x41*x53
     6  +coeff( 33)*x11*x21*x31        
     7  +coeff( 34)*x11*x21    *x41    
     8  +coeff( 35)    *x22*x32*x41    
      p_mnq_0_9   =p_mnq_0_9   
     9  +coeff( 36)    *x22*x31*x42    
     1  +coeff( 37)        *x33*x42    
     2  +coeff( 38)    *x22    *x43    
     3  +coeff( 39)        *x32*x43    
     4  +coeff( 40)*x11    *x31    *x51
     5  +coeff( 41)    *x23*x31    *x51
     6  +coeff( 42)    *x21*x33    *x51
     7  +coeff( 43)    *x23    *x41*x51
     8  +coeff( 44)    *x21*x31*x42*x51
      p_mnq_0_9   =p_mnq_0_9   
     9  +coeff( 45)    *x21    *x43*x51
     1  +coeff( 46)    *x22*x31    *x52
     2  +coeff( 47)        *x33    *x52
     3  +coeff( 48)    *x22    *x41*x52
     4  +coeff( 49)        *x32*x41*x52
     5  +coeff( 50)        *x31*x42*x52
     6  +coeff( 51)            *x43*x52
     7  +coeff( 52)    *x21    *x41*x53
     8  +coeff( 53)*x11*x22*x31        
      p_mnq_0_9   =p_mnq_0_9   
     9  +coeff( 54)*x11    *x33        
     1  +coeff( 55)*x11*x22    *x41    
     2  +coeff( 56)*x11    *x32*x41    
     3  +coeff( 57)*x11    *x31*x42    
     4  +coeff( 58)    *x23*x31*x42    
     5  +coeff( 59)    *x21*x33*x42    
     6  +coeff( 60)*x11        *x43    
     7  +coeff( 61)    *x23    *x43    
     8  +coeff( 62)    *x21*x32*x43    
      p_mnq_0_9   =p_mnq_0_9   
     9  +coeff( 63)    *x22*x33    *x51
     1  +coeff( 64)*x11*x21    *x41*x51
     2  +coeff( 65)    *x22*x32*x41*x51
     3  +coeff( 66)    *x22*x31*x42*x51
     4  +coeff( 67)    *x22    *x43*x51
     5  +coeff( 68)        *x32*x43*x51
     6  +coeff( 69)    *x23*x31    *x52
     7  +coeff( 70)    *x21*x33    *x52
     8  +coeff( 71)*x11        *x41*x52
      p_mnq_0_9   =p_mnq_0_9   
     9  +coeff( 72)    *x23    *x41*x52
     1  +coeff( 73)    *x21*x32*x41*x52
     2  +coeff( 74)    *x21*x31*x42*x52
     3  +coeff( 75)    *x21    *x43*x52
     4  +coeff( 76)        *x33    *x53
     5  +coeff( 77)    *x22    *x41*x53
     6  +coeff( 78)        *x32*x41*x53
     7  +coeff( 79)        *x31*x42*x53
     8  +coeff( 80)            *x43*x53
      p_mnq_0_9   =p_mnq_0_9   
     9  +coeff( 81)*x11*x21*x33        
     1  +coeff( 82)*x11*x23    *x41    
     2  +coeff( 83)*x11*x21*x31*x42    
     3  +coeff( 84)    *x22*x33*x42    
     4  +coeff( 85)*x11*x21    *x43    
     5  +coeff( 86)    *x22*x32*x43    
     6  +coeff( 87)*x11    *x33    *x51
     7  +coeff( 88)    *x23*x33    *x51
     8  +coeff( 89)*x11*x22    *x41*x51
      p_mnq_0_9   =p_mnq_0_9   
     9  +coeff( 90)    *x23*x32*x41*x51
c
      return
      end
      function pl_mnq_0_9  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 41)
      data ncoeff/ 40/
      data avdat/  0.1220282E-01/
      data xmin/
     1 -0.49978E-02,-0.50007E-01,-0.10393E+00,-0.46630E-01,-0.42140E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49966E-02, 0.49989E-01, 0.10393E+00, 0.46630E-01, 0.49976E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.14203655E-01,-0.18534650E+00,-0.39892480E-01,-0.41589539E-01,
     +  0.49634796E-01,-0.28783191E-01,-0.17074071E-01, 0.57012338E-01,
     + -0.34378931E-01,-0.17735203E-02,-0.23037575E-01, 0.12376651E-02,
     + -0.78377128E-01, 0.10636464E-02, 0.69160522E-02, 0.43915599E-02,
     + -0.23640406E-02, 0.74573938E-04,-0.85657761E-02,-0.58240453E-02,
     +  0.14802325E-01,-0.63318230E-01, 0.87806897E-03,-0.31358234E-03,
     +  0.34855762E-02, 0.18331381E-01, 0.50209172E-03, 0.29554211E-02,
     + -0.53687178E-03, 0.25365522E-01,-0.42226417E-02, 0.48563894E-01,
     + -0.15260311E-02,-0.39826590E-02,-0.90747736E-02, 0.21543691E-03,
     + -0.14970075E-01,-0.26445176E-01,-0.78511588E-01,-0.70984408E-01,
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
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
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
      x54 = x53*x5
c
c                  function
c
      pl_mnq_0_9  =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21        *x51
     4  +coeff(  4)                *x52
     5  +coeff(  5)    *x21        *x52
     6  +coeff(  6)                *x51
     7  +coeff(  7)            *x42    
     8  +coeff(  8)                *x53
      pl_mnq_0_9  =pl_mnq_0_9  
     9  +coeff(  9)                *x54
     1  +coeff( 10)*x11                
     2  +coeff( 11)    *x22            
     3  +coeff( 12)    *x21    *x42    
     4  +coeff( 13)    *x21        *x54
     5  +coeff( 14)    *x21*x31*x41    
     6  +coeff( 15)            *x42*x51
     7  +coeff( 16)    *x24            
     8  +coeff( 17)    *x22*x31*x41    
      pl_mnq_0_9  =pl_mnq_0_9  
     9  +coeff( 18)        *x33*x41    
     1  +coeff( 19)    *x22    *x42    
     2  +coeff( 20)            *x42*x52
     3  +coeff( 21)    *x21        *x53
     4  +coeff( 22)    *x22        *x53
     5  +coeff( 23)    *x24    *x42*x51
     6  +coeff( 24)        *x31*x41    
     7  +coeff( 25)    *x23            
     8  +coeff( 26)    *x22        *x51
      pl_mnq_0_9  =pl_mnq_0_9  
     9  +coeff( 27)        *x32    *x51
     1  +coeff( 28)        *x31*x41*x51
     2  +coeff( 29)    *x22*x32        
     3  +coeff( 30)    *x23        *x51
     4  +coeff( 31)    *x21    *x42*x51
     5  +coeff( 32)    *x22        *x52
     6  +coeff( 33)        *x31*x41*x52
     7  +coeff( 34)    *x23*x31*x41    
     8  +coeff( 35)    *x23    *x42    
      pl_mnq_0_9  =pl_mnq_0_9  
     9  +coeff( 36)*x11*x22        *x51
     1  +coeff( 37)    *x23        *x52
     2  +coeff( 38)    *x24        *x52
     3  +coeff( 39)    *x23        *x53
     4  +coeff( 40)    *x22        *x54
c
      return
      end
      function x_mnq_0_11  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.4183537E-01/
      data xmin/
     1 -0.49978E-02,-0.50007E-01,-0.10393E+00,-0.45852E-01,-0.39669E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49966E-02, 0.49912E-01, 0.10393E+00, 0.45852E-01, 0.49976E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.14961095E+00, 0.78833199E+00, 0.91763389E+00, 0.53142037E-01,
     + -0.13868996E-02,-0.16187083E-01,-0.44279363E-01, 0.84765092E-01,
     + -0.37376538E+00, 0.44949083E-02, 0.39499983E-01, 0.25176541E-02,
     + -0.41725922E-01,-0.12333430E+00,-0.17878823E-01, 0.43569587E-03,
     + -0.25056146E-01,-0.82426324E-01,-0.76197036E-01, 0.12273461E+00,
     + -0.12617026E-01,-0.19284575E+00,-0.93099443E-04,-0.55212289E+00,
     +  0.29337301E-02,-0.62593366E-02,-0.11976163E+00,-0.60345316E+00,
     + -0.38642664E-01,-0.85219964E-02,-0.37654329E-01,-0.13039590E+00,
     +  0.10322398E+00,-0.17961126E-01,-0.16941096E+00,-0.24338986E-03,
     + -0.50347984E+00,-0.42819660E-02,-0.15184173E+00,-0.37703019E+00,
     +  0.11940251E+00, 0.42778393E-02, 0.46532154E-01, 0.27484557E+00,
     +  0.16389447E+00, 0.47553470E-03, 0.58249552E-01, 0.26505825E+00,
     +  0.47478261E-02, 0.48284144E-02, 0.44389024E-01, 0.20206779E+00,
     +  0.10106649E-01, 0.77486199E+00,-0.23512472E-02, 0.36332464E+00,
     + -0.92346063E-02, 0.14478958E+01, 0.39635021E-01, 0.40220367E-02,
     +  0.14446504E+00, 0.74466908E+00, 0.14078918E-01,-0.50976174E-02,
     + -0.52997854E-01,-0.61879924E-03, 0.19463206E-03,-0.20362858E-01,
     + -0.58557233E-02,-0.35315342E-01,-0.93893800E-02,-0.95832953E-02,
     + -0.39873734E-01, 0.69113970E-01,-0.23846335E-02, 0.32892404E-02,
     + -0.21393606E-01,-0.34352589E-01,-0.22086060E-01,-0.28106252E-01,
     + -0.16608480E-01, 0.33314146E-01,-0.17840944E-01,-0.18543400E+00,
     +  0.11090611E-02, 0.76277978E-02, 0.35586074E-01, 0.90492470E-02,
     + -0.23000264E-02, 0.38254503E-01,
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
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_mnq_0_11  =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22            
     5  +coeff(  5)        *x32        
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21        *x51
      x_mnq_0_11  =x_mnq_0_11  
     9  +coeff(  9)                *x52
     1  +coeff( 10)*x11                
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x21*x32        
     4  +coeff( 13)    *x21*x31*x41    
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)    *x22        *x51
     7  +coeff( 16)        *x32    *x51
     8  +coeff( 17)        *x31*x41*x51
      x_mnq_0_11  =x_mnq_0_11  
     9  +coeff( 18)            *x42*x51
     1  +coeff( 19)    *x21        *x52
     2  +coeff( 20)                *x53
     3  +coeff( 21)    *x22*x32        
     4  +coeff( 22)    *x22*x31*x41    
     5  +coeff( 23)        *x33*x41    
     6  +coeff( 24)    *x22    *x42    
     7  +coeff( 25)        *x32*x42    
     8  +coeff( 26)        *x31*x43    
      x_mnq_0_11  =x_mnq_0_11  
     9  +coeff( 27)    *x23        *x51
     1  +coeff( 28)    *x21    *x42*x51
     2  +coeff( 29)    *x22        *x52
     3  +coeff( 30)        *x32    *x52
     4  +coeff( 31)        *x31*x41*x52
     5  +coeff( 32)            *x42*x52
     6  +coeff( 33)    *x21        *x53
     7  +coeff( 34)    *x23*x32        
     8  +coeff( 35)    *x23*x31*x41    
      x_mnq_0_11  =x_mnq_0_11  
     9  +coeff( 36)*x11        *x42    
     1  +coeff( 37)    *x23    *x42    
     2  +coeff( 38)    *x21*x32*x42    
     3  +coeff( 39)    *x22*x31*x41*x51
     4  +coeff( 40)    *x22    *x42*x51
     5  +coeff( 41)    *x23        *x52
     6  +coeff( 42)    *x21*x32    *x52
     7  +coeff( 43)    *x21*x31*x41*x52
     8  +coeff( 44)    *x21    *x42*x52
      x_mnq_0_11  =x_mnq_0_11  
     9  +coeff( 45)    *x22        *x53
     1  +coeff( 46)        *x32    *x53
     2  +coeff( 47)        *x31*x41*x53
     3  +coeff( 48)            *x42*x53
     4  +coeff( 49)*x11*x23            
     5  +coeff( 50)*x11*x22        *x51
     6  +coeff( 51)    *x23*x32    *x51
     7  +coeff( 52)    *x23*x31*x41*x51
     8  +coeff( 53)    *x21*x33*x41*x51
      x_mnq_0_11  =x_mnq_0_11  
     9  +coeff( 54)    *x23    *x42*x51
     1  +coeff( 55)    *x21*x31*x43*x51
     2  +coeff( 56)    *x22*x31*x41*x52
     3  +coeff( 57)        *x33*x41*x52
     4  +coeff( 58)    *x22    *x42*x52
     5  +coeff( 59)    *x23        *x53
     6  +coeff( 60)    *x21*x32    *x53
     7  +coeff( 61)    *x21*x31*x41*x53
     8  +coeff( 62)    *x21    *x42*x53
      x_mnq_0_11  =x_mnq_0_11  
     9  +coeff( 63)    *x23*x33*x41    
     1  +coeff( 64)    *x22*x31*x43*x51
     2  +coeff( 65)    *x23*x31*x41*x52
     3  +coeff( 66)    *x21*x33*x41*x52
     4  +coeff( 67)*x11        *x42*x52
     5  +coeff( 68)    *x22*x31*x41*x53
     6  +coeff( 69)        *x33*x41*x53
     7  +coeff( 70)    *x22    *x42*x53
     8  +coeff( 71)*x11*x23    *x42    
      x_mnq_0_11  =x_mnq_0_11  
     9  +coeff( 72)    *x23*x33*x41*x51
     1  +coeff( 73)    *x22*x31*x43*x52
     2  +coeff( 74)    *x23*x31*x41*x53
     3  +coeff( 75)*x11*x23    *x42*x51
     4  +coeff( 76)*x11*x22*x32    *x52
     5  +coeff( 77)    *x23*x33*x43*x51
     6  +coeff( 78)*x11*x21*x33*x41*x52
     7  +coeff( 79)*x11*x22*x31*x41*x53
     8  +coeff( 80)    *x23*x33*x41*x53
      x_mnq_0_11  =x_mnq_0_11  
     9  +coeff( 81)*x11*x22*x33*x43    
     1  +coeff( 82)*x11*x23    *x42*x53
     2  +coeff( 83)    *x21*x32    *x51
     3  +coeff( 84)    *x21*x31*x41*x51
     4  +coeff( 85)*x11*x22            
     5  +coeff( 86)        *x33*x41*x51
     6  +coeff( 87)    *x22*x31*x43    
     7  +coeff( 88)        *x33*x43    
     8  +coeff( 89)*x11    *x31*x41*x51
      x_mnq_0_11  =x_mnq_0_11  
     9  +coeff( 90)    *x22*x32    *x52
c
      return
      end
      function t_mnq_0_11  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.1426885E-01/
      data xmin/
     1 -0.49978E-02,-0.50007E-01,-0.10393E+00,-0.45852E-01,-0.39669E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49966E-02, 0.49912E-01, 0.10393E+00, 0.45852E-01, 0.49976E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.10272787E-01, 0.44507608E-01, 0.21012844E+00, 0.16099250E-01,
     + -0.26141846E-03,-0.27570073E-02,-0.10349214E-01, 0.23636589E-01,
     + -0.85604645E-01, 0.14517169E-01, 0.11343396E-02,-0.14364465E-01,
     + -0.44276033E-01, 0.83774263E-02, 0.30223877E-03,-0.97328015E-02,
     + -0.30673599E-01,-0.13896993E-01, 0.29224180E-01,-0.20038665E-02,
     + -0.76292275E-03,-0.18260425E+00, 0.35509744E-03,-0.44324193E-02,
     + -0.23127770E-01,-0.76227650E-01,-0.21709785E+00,-0.66832229E-02,
     + -0.28612469E-02,-0.18099831E-01,-0.51225662E-01, 0.26446585E-01,
     + -0.14172134E-02,-0.84629057E-04,-0.39752475E-02,-0.16378522E+00,
     + -0.62497724E-01,-0.20007561E+00,-0.30894813E-02, 0.30346610E-01,
     +  0.12443452E-02, 0.22986748E-02, 0.30298209E-01, 0.42453427E-01,
     +  0.17739842E-01, 0.72248451E-01,-0.56373206E-03, 0.41156914E-02,
     + -0.59670494E-02, 0.14333377E-01, 0.50414298E-01, 0.21147050E-02,
     +  0.27066311E-01, 0.10122740E+00,-0.30751836E-02, 0.34748337E+00,
     +  0.10326186E-01, 0.56677584E-01, 0.20151642E+00, 0.70590249E-04,
     + -0.11415604E-03,-0.63508995E-01,-0.47548665E-02,-0.52247480E-01,
     + -0.55447337E-02,-0.21163383E-02,-0.10314385E-03,-0.90902526E-03,
     +  0.46434769E-03,-0.14310905E-02, 0.41279085E-02, 0.10090765E-01,
     +  0.14997648E+00, 0.80509614E-02, 0.15922330E-01,-0.32138100E-03,
     +  0.25293047E-02,-0.42738431E-03,-0.32192809E-02, 0.16874050E-02,
     + -0.73788641E-03, 0.34329016E-02, 0.63227350E-03,-0.23574347E-03,
     + -0.20217423E-02,-0.11367488E-03,-0.58302854E-03,-0.87217864E-03,
     +  0.10244570E-03,-0.34404616E-03,
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
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_mnq_0_11  =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22            
     5  +coeff(  5)        *x32        
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21        *x51
      t_mnq_0_11  =t_mnq_0_11  
     9  +coeff(  9)                *x52
     1  +coeff( 10)    *x23            
     2  +coeff( 11)    *x21*x32        
     3  +coeff( 12)    *x21*x31*x41    
     4  +coeff( 13)    *x21    *x42    
     5  +coeff( 14)    *x22        *x51
     6  +coeff( 15)        *x32    *x51
     7  +coeff( 16)        *x31*x41*x51
     8  +coeff( 17)            *x42*x51
      t_mnq_0_11  =t_mnq_0_11  
     9  +coeff( 18)    *x21        *x52
     1  +coeff( 19)                *x53
     2  +coeff( 20)    *x22*x32        
     3  +coeff( 21)        *x33*x41    
     4  +coeff( 22)    *x22    *x42    
     5  +coeff( 23)        *x32*x42    
     6  +coeff( 24)        *x31*x43    
     7  +coeff( 25)    *x23        *x51
     8  +coeff( 26)    *x21*x31*x41*x51
      t_mnq_0_11  =t_mnq_0_11  
     9  +coeff( 27)    *x21    *x42*x51
     1  +coeff( 28)    *x22        *x52
     2  +coeff( 29)        *x32    *x52
     3  +coeff( 30)        *x31*x41*x52
     4  +coeff( 31)            *x42*x52
     5  +coeff( 32)    *x21        *x53
     6  +coeff( 33)*x11*x22            
     7  +coeff( 34)*x11    *x32        
     8  +coeff( 35)    *x23*x32        
      t_mnq_0_11  =t_mnq_0_11  
     9  +coeff( 36)    *x23    *x42    
     1  +coeff( 37)    *x22*x31*x41*x51
     2  +coeff( 38)    *x22    *x42*x51
     3  +coeff( 39)        *x32*x42*x51
     4  +coeff( 40)    *x23        *x52
     5  +coeff( 41)    *x21*x32    *x52
     6  +coeff( 42)    *x21*x31*x41*x52
     7  +coeff( 43)    *x21    *x42*x52
     8  +coeff( 44)    *x22        *x53
      t_mnq_0_11  =t_mnq_0_11  
     9  +coeff( 45)        *x31*x41*x53
     1  +coeff( 46)            *x42*x53
     2  +coeff( 47)*x11*x23            
     3  +coeff( 48)    *x22*x33*x41    
     4  +coeff( 49)    *x22*x32*x42    
     5  +coeff( 50)    *x22*x31*x43    
     6  +coeff( 51)    *x23*x31*x41*x51
     7  +coeff( 52)    *x21*x33*x41*x51
     8  +coeff( 53)    *x21*x31*x43*x51
      t_mnq_0_11  =t_mnq_0_11  
     9  +coeff( 54)    *x22*x31*x41*x52
     1  +coeff( 55)        *x33*x41*x52
     2  +coeff( 56)    *x22    *x42*x52
     3  +coeff( 57)    *x23        *x53
     4  +coeff( 58)    *x21*x31*x41*x53
     5  +coeff( 59)    *x21    *x42*x53
     6  +coeff( 60)*x11                
     7  +coeff( 61)*x11*x21            
     8  +coeff( 62)    *x22*x31*x41    
      t_mnq_0_11  =t_mnq_0_11  
     9  +coeff( 63)    *x21*x32    *x51
     1  +coeff( 64)    *x23*x31*x41    
     2  +coeff( 65)    *x21*x32*x42    
     3  +coeff( 66)*x11*x21        *x51
     4  +coeff( 67)    *x22*x32    *x51
     5  +coeff( 68)*x11            *x52
     6  +coeff( 69)        *x32    *x53
     7  +coeff( 70)*x11*x21    *x42    
     8  +coeff( 71)        *x33*x43    
      t_mnq_0_11  =t_mnq_0_11  
     9  +coeff( 72)    *x23*x32    *x51
     1  +coeff( 73)    *x23    *x42*x51
     2  +coeff( 74)    *x22*x32    *x52
     3  +coeff( 75)        *x31*x43*x52
     4  +coeff( 76)*x11    *x31*x41    
     5  +coeff( 77)    *x21*x33*x41    
     6  +coeff( 78)*x11        *x42    
     7  +coeff( 79)    *x21*x31*x43    
     8  +coeff( 80)        *x33*x41*x51
      t_mnq_0_11  =t_mnq_0_11  
     9  +coeff( 81)*x11        *x42*x51
     1  +coeff( 82)        *x32*x42*x52
     2  +coeff( 83)*x11            *x53
     3  +coeff( 84)*x11            *x51
     4  +coeff( 85)        *x31*x43*x51
     5  +coeff( 86)*x11*x21*x32        
     6  +coeff( 87)*x11*x21*x31*x41    
     7  +coeff( 88)*x11*x22        *x51
     8  +coeff( 89)*x11    *x32    *x51
      t_mnq_0_11  =t_mnq_0_11  
     9  +coeff( 90)*x11    *x31*x41*x51
c
      return
      end
      function y_mnq_0_11  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49978E-02,-0.50007E-01,-0.10393E+00,-0.45852E-01,-0.39669E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49966E-02, 0.49912E-01, 0.10393E+00, 0.45852E-01, 0.49976E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.94563365E-01, 0.70498377E+00,-0.14068826E-01,-0.60974102E-01,
     + -0.92763628E-04,-0.41377228E-02,-0.61706752E-01,-0.37227127E+00,
     +  0.79574238E-03, 0.14595757E-02,-0.34555715E+00,-0.17967219E-01,
     + -0.10681180E+00,-0.78405298E-01,-0.44655025E+00, 0.12094828E-01,
     +  0.12096473E+00, 0.24769419E+00,-0.70166804E-01,-0.44498643E+00,
     +  0.16436031E+00,-0.27752377E-01, 0.13828060E-01,-0.26055200E-02,
     +  0.40228438E-01, 0.50747174E+00,-0.43416554E-02, 0.10994262E+01,
     +  0.50878756E-01, 0.32410529E+00, 0.48697343E+00, 0.11302730E+01,
     +  0.10466554E+00, 0.64894730E+00, 0.10505918E+00, 0.24349017E+00,
     +  0.41163109E-01, 0.27609840E+00,-0.69712300E-03, 0.63960755E+00,
     +  0.13111901E+01,-0.31311274E-01,-0.83064465E-02, 0.15619972E+01,
     + -0.16072340E-02, 0.10633042E-01,-0.79688005E-01, 0.31562693E-01,
     +  0.19195046E+00,-0.80385514E-01,-0.18674729E-01,-0.12525223E+00,
     + -0.23872225E+00, 0.18861124E-01,-0.12832084E+00,-0.45445547E+00,
     + -0.76316828E+00,-0.92366701E+00,-0.18017597E+01, 0.66484608E-01,
     + -0.29060721E-01,-0.41299215E+00,-0.95528996E+00,-0.61798305E-02,
     +  0.16324112E-01,-0.35877272E-01, 0.43081790E-02,-0.48697133E-01,
     + -0.25627303E-02,-0.89516671E-03, 0.73640207E-02, 0.79685070E-01,
     + -0.89164926E-02, 0.91973610E-01,-0.48726117E-02,-0.90812817E-02,
     +  0.63115917E-01, 0.66962099E+00,-0.74031830E-01,-0.53304071E-02,
     + -0.23637167E+00, 0.17718907E-01,-0.79039320E-01, 0.76512108E-02,
     +  0.25947807E-01, 0.47485858E-01,-0.82383007E-01,-0.63052736E-01,
     + -0.73418417E-02,-0.35129848E-03,
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
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      y_mnq_0_11  =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21*x31        
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)        *x31    *x51
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)    *x22*x31        
     8  +coeff(  8)    *x22    *x41    
      y_mnq_0_11  =y_mnq_0_11  
     9  +coeff(  9)        *x32*x41    
     1  +coeff( 10)            *x43    
     2  +coeff( 11)    *x21    *x41*x51
     3  +coeff( 12)        *x31    *x52
     4  +coeff( 13)            *x41*x52
     5  +coeff( 14)    *x23*x31        
     6  +coeff( 15)    *x23    *x41    
     7  +coeff( 16)    *x21*x32*x41    
     8  +coeff( 17)    *x21*x31*x42    
      y_mnq_0_11  =y_mnq_0_11  
     9  +coeff( 18)    *x21    *x43    
     1  +coeff( 19)    *x22*x31    *x51
     2  +coeff( 20)    *x22    *x41*x51
     3  +coeff( 21)            *x43*x51
     4  +coeff( 22)    *x21    *x41*x52
     5  +coeff( 23)        *x31    *x53
     6  +coeff( 24)*x11*x21    *x41    
     7  +coeff( 25)    *x22*x32*x41    
     8  +coeff( 26)    *x22*x31*x42    
      y_mnq_0_11  =y_mnq_0_11  
     9  +coeff( 27)        *x33*x42    
     1  +coeff( 28)    *x22    *x43    
     2  +coeff( 29)    *x23*x31    *x51
     3  +coeff( 30)    *x23    *x41*x51
     4  +coeff( 31)    *x21*x31*x42*x51
     5  +coeff( 32)    *x21    *x43*x51
     6  +coeff( 33)    *x22*x31    *x52
     7  +coeff( 34)    *x22    *x41*x52
     8  +coeff( 35)        *x31*x42*x52
      y_mnq_0_11  =y_mnq_0_11  
     9  +coeff( 36)            *x43*x52
     1  +coeff( 37)    *x21*x31    *x53
     2  +coeff( 38)    *x21    *x41*x53
     3  +coeff( 39)*x11*x22    *x41    
     4  +coeff( 40)    *x23*x31*x42    
     5  +coeff( 41)    *x23    *x43    
     6  +coeff( 42)    *x22*x32*x41*x51
     7  +coeff( 43)        *x33*x42*x51
     8  +coeff( 44)    *x22    *x43*x51
      y_mnq_0_11  =y_mnq_0_11  
     9  +coeff( 45)        *x32*x43*x51
     1  +coeff( 46)    *x21*x33    *x52
     2  +coeff( 47)    *x23    *x41*x52
     3  +coeff( 48)    *x21*x31*x42*x52
     4  +coeff( 49)    *x21    *x43*x52
     5  +coeff( 50)    *x22    *x41*x53
     6  +coeff( 51)        *x32*x41*x53
     7  +coeff( 52)        *x31*x42*x53
     8  +coeff( 53)            *x43*x53
      y_mnq_0_11  =y_mnq_0_11  
     9  +coeff( 54)    *x22*x33*x42    
     1  +coeff( 55)    *x23*x32*x41*x51
     2  +coeff( 56)    *x23*x31*x42*x51
     3  +coeff( 57)    *x23    *x43*x51
     4  +coeff( 58)    *x22*x31*x42*x52
     5  +coeff( 59)    *x22    *x43*x52
     6  +coeff( 60)    *x23    *x41*x53
     7  +coeff( 61)    *x21*x32*x41*x53
     8  +coeff( 62)    *x21*x31*x42*x53
      y_mnq_0_11  =y_mnq_0_11  
     9  +coeff( 63)    *x21    *x43*x53
     1  +coeff( 64)    *x23*x33    *x52
     2  +coeff( 65)    *x21*x33*x42*x52
     3  +coeff( 66)    *x23*x33*x42*x52
     4  +coeff( 67)        *x31*x42    
     5  +coeff( 68)    *x21*x31    *x51
     6  +coeff( 69)    *x21*x33        
     7  +coeff( 70)        *x33    *x51
     8  +coeff( 71)        *x32*x41*x51
      y_mnq_0_11  =y_mnq_0_11  
     9  +coeff( 72)        *x31*x42*x51
     1  +coeff( 73)    *x21*x31    *x52
     2  +coeff( 74)            *x41*x53
     3  +coeff( 75)    *x22*x33        
     4  +coeff( 76)    *x21*x33    *x51
     5  +coeff( 77)    *x23*x32*x41    
     6  +coeff( 78)    *x22*x31*x42*x51
     7  +coeff( 79)    *x21*x32*x41*x52
     8  +coeff( 80)*x11*x22    *x41*x51
      y_mnq_0_11  =y_mnq_0_11  
     9  +coeff( 81)    *x22*x32*x41*x52
     1  +coeff( 82)    *x23*x31    *x53
     2  +coeff( 83)    *x23    *x43*x52
     3  +coeff( 84)    *x22*x32*x41*x53
     4  +coeff( 85)        *x33*x42*x53
     5  +coeff( 86)    *x23*x33*x42*x51
     6  +coeff( 87)    *x23*x32*x41*x53
     7  +coeff( 88)    *x23*x31*x42*x53
     8  +coeff( 89)*x11*x23*x32*x41*x52
      y_mnq_0_11  =y_mnq_0_11  
     9  +coeff( 90)*x11        *x41    
c
      return
      end
      function p_mnq_0_11  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49978E-02,-0.50007E-01,-0.10393E+00,-0.45852E-01,-0.39669E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49966E-02, 0.49912E-01, 0.10393E+00, 0.45852E-01, 0.49976E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.16644489E-02, 0.32341130E-01,-0.23773074E-01, 0.20506042E-02,
     + -0.37811871E-02,-0.15106282E-01,-0.11793287E-02,-0.96225649E-01,
     + -0.12841923E-03, 0.29868556E-02,-0.14712608E-01,-0.11229252E+00,
     + -0.62493519E-02,-0.35611421E-01,-0.81248058E-04,-0.19352542E-01,
     + -0.64445278E-02,-0.19533167E-03,-0.10167070E+00, 0.33048277E-02,
     +  0.35984967E-01, 0.76526068E-01,-0.26667958E-01,-0.45991205E-02,
     + -0.13949256E+00, 0.14316881E-02, 0.24610559E-01, 0.52899517E-01,
     + -0.15319862E-01,-0.41841600E-01,-0.12635755E-02, 0.10586669E-01,
     + -0.88429492E-03, 0.27879467E-02,-0.55600930E-03, 0.70427526E-02,
     +  0.91300733E-01, 0.40213028E-02, 0.24751662E+00, 0.16296060E-02,
     + -0.19587245E-03, 0.14979450E-02, 0.16182859E-02, 0.22260848E-01,
     +  0.40793787E-02, 0.90187095E-01, 0.29636148E+00, 0.36573627E-02,
     +  0.14604999E-02, 0.34891471E-01, 0.18182082E-01, 0.72525546E-01,
     +  0.20013943E-02, 0.20744642E-01,-0.16960481E-02, 0.31605083E-03,
     +  0.13384633E-01, 0.56826574E-03, 0.85088002E-04, 0.48352577E-01,
     +  0.73263841E-02, 0.24162515E+00, 0.42707454E-02,-0.29117821E-02,
     +  0.28421568E-01, 0.49527462E-04,-0.72146626E-02, 0.80211749E-02,
     +  0.33626422E+00, 0.51393919E-02,-0.12038630E-02,-0.40216981E-02,
     +  0.27855974E-01,-0.84125903E-02,-0.73129371E-01, 0.76514795E-01,
     + -0.45944527E-02, 0.90406677E-02,-0.61474121E-02, 0.20121622E-02,
     + -0.35766233E-01,-0.23279749E-01,-0.51476358E-03, 0.11019282E-02,
     +  0.96927362E-03, 0.82318566E-03,-0.87493649E-02, 0.51735691E-03,
     + -0.12380993E-02, 0.53154887E-03,
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
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_mnq_0_11  =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)        *x31    *x51
     5  +coeff(  5)            *x41*x51
     6  +coeff(  6)    *x22*x31        
     7  +coeff(  7)        *x33        
     8  +coeff(  8)    *x22    *x41    
      p_mnq_0_11  =p_mnq_0_11  
     9  +coeff(  9)        *x32*x41    
     1  +coeff( 10)            *x43    
     2  +coeff( 11)    *x21*x31    *x51
     3  +coeff( 12)    *x21    *x41*x51
     4  +coeff( 13)        *x31    *x52
     5  +coeff( 14)            *x41*x52
     6  +coeff( 15)*x11    *x31        
     7  +coeff( 16)    *x23*x31        
     8  +coeff( 17)    *x21*x33        
      p_mnq_0_11  =p_mnq_0_11  
     9  +coeff( 18)*x11        *x41    
     1  +coeff( 19)    *x23    *x41    
     2  +coeff( 20)    *x21*x32*x41    
     3  +coeff( 21)    *x21*x31*x42    
     4  +coeff( 22)    *x21    *x43    
     5  +coeff( 23)    *x22*x31    *x51
     6  +coeff( 24)        *x33    *x51
     7  +coeff( 25)    *x22    *x41*x51
     8  +coeff( 26)        *x32*x41*x51
      p_mnq_0_11  =p_mnq_0_11  
     9  +coeff( 27)        *x31*x42*x51
     1  +coeff( 28)            *x43*x51
     2  +coeff( 29)    *x21*x31    *x52
     3  +coeff( 30)    *x21    *x41*x52
     4  +coeff( 31)        *x31    *x53
     5  +coeff( 32)            *x41*x53
     6  +coeff( 33)*x11*x21*x31        
     7  +coeff( 34)    *x22*x33        
     8  +coeff( 35)*x11*x21    *x41    
      p_mnq_0_11  =p_mnq_0_11  
     9  +coeff( 36)    *x22*x32*x41    
     1  +coeff( 37)    *x22*x31*x42    
     2  +coeff( 38)        *x33*x42    
     3  +coeff( 39)    *x22    *x43    
     4  +coeff( 40)        *x32*x43    
     5  +coeff( 41)*x11    *x31    *x51
     6  +coeff( 42)    *x23*x31    *x51
     7  +coeff( 43)    *x21*x33    *x51
     8  +coeff( 44)    *x23    *x41*x51
      p_mnq_0_11  =p_mnq_0_11  
     9  +coeff( 45)    *x21*x32*x41*x51
     1  +coeff( 46)    *x21*x31*x42*x51
     2  +coeff( 47)    *x21    *x43*x51
     3  +coeff( 48)    *x22*x31    *x52
     4  +coeff( 49)        *x33    *x52
     5  +coeff( 50)    *x22    *x41*x52
     6  +coeff( 51)        *x31*x42*x52
     7  +coeff( 52)            *x43*x52
     8  +coeff( 53)    *x21*x31    *x53
      p_mnq_0_11  =p_mnq_0_11  
     9  +coeff( 54)    *x21    *x41*x53
     1  +coeff( 55)*x11*x22*x31        
     2  +coeff( 56)*x11    *x33        
     3  +coeff( 57)    *x23*x33        
     4  +coeff( 58)*x11*x22    *x41    
     5  +coeff( 59)*x11    *x31*x42    
     6  +coeff( 60)    *x23*x31*x42    
     7  +coeff( 61)    *x21*x33*x42    
     8  +coeff( 62)    *x23    *x43    
      p_mnq_0_11  =p_mnq_0_11  
     9  +coeff( 63)    *x21*x32*x43    
     1  +coeff( 64)*x11*x21*x31    *x51
     2  +coeff( 65)    *x22*x33    *x51
     3  +coeff( 66)*x11*x21    *x41*x51
     4  +coeff( 67)    *x22*x32*x41*x51
     5  +coeff( 68)        *x33*x42*x51
     6  +coeff( 69)    *x22    *x43*x51
     7  +coeff( 70)        *x32*x43*x51
     8  +coeff( 71)*x11    *x31    *x52
      p_mnq_0_11  =p_mnq_0_11  
     9  +coeff( 72)    *x23*x31    *x52
     1  +coeff( 73)    *x21*x33    *x52
     2  +coeff( 74)    *x23    *x41*x52
     3  +coeff( 75)    *x21*x31*x42*x52
     4  +coeff( 76)    *x21    *x43*x52
     5  +coeff( 77)    *x22*x31    *x53
     6  +coeff( 78)        *x33    *x53
     7  +coeff( 79)    *x22    *x41*x53
     8  +coeff( 80)        *x32*x41*x53
      p_mnq_0_11  =p_mnq_0_11  
     9  +coeff( 81)        *x31*x42*x53
     1  +coeff( 82)            *x43*x53
     2  +coeff( 83)*x11*x23*x31        
     3  +coeff( 84)*x11*x21*x33        
     4  +coeff( 85)*x11*x23    *x41    
     5  +coeff( 86)*x11*x21*x31*x42    
     6  +coeff( 87)    *x22*x33*x42    
     7  +coeff( 88)*x11*x21    *x43    
     8  +coeff( 89)*x11*x22*x31    *x51
      p_mnq_0_11  =p_mnq_0_11  
     9  +coeff( 90)*x11    *x33    *x51
c
      return
      end
      function pl_mnq_0_11 (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 41)
      data ncoeff/ 40/
      data avdat/  0.2026542E-01/
      data xmin/
     1 -0.49978E-02,-0.50007E-01,-0.10393E+00,-0.45852E-01,-0.39669E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49966E-02, 0.49912E-01, 0.10393E+00, 0.45852E-01, 0.49976E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.38108803E-01,-0.31537014E+00,-0.15952028E+00,-0.19160416E-01,
     +  0.16974121E-01,-0.25913308E-02, 0.29444103E-02,-0.74186292E-02,
     + -0.49232036E-01,-0.15476746E-01, 0.82637963E-03,-0.27771246E-01,
     + -0.21681713E-02, 0.33204682E-01,-0.28291652E-01, 0.80740619E-02,
     + -0.26836997E-01,-0.17176092E-01, 0.13244323E-01,-0.20061506E-01,
     +  0.63783387E-02,-0.62779733E-02,-0.49224189E-02,-0.58462038E-02,
     +  0.56134914E-02,-0.37192954E-02,-0.43931664E-02, 0.33589549E-01,
     + -0.26093975E-02,-0.83231945E-02, 0.19873118E-01, 0.33763055E-01,
     +  0.37593443E-01,-0.19936545E-01, 0.19324829E-02, 0.13560189E-01,
     +  0.40442315E-02, 0.15203436E-03, 0.82269334E-03, 0.45185382E-02,
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
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
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
      pl_mnq_0_11 =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22        *x51
     5  +coeff(  5)    *x21        *x52
     6  +coeff(  6)*x11                
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)            *x42    
      pl_mnq_0_11 =pl_mnq_0_11 
     9  +coeff(  9)    *x21        *x51
     1  +coeff( 10)    *x24            
     2  +coeff( 11)    *x22        *x52
     3  +coeff( 12)    *x22            
     4  +coeff( 13)                *x52
     5  +coeff( 14)                *x53
     6  +coeff( 15)    *x23        *x51
     7  +coeff( 16)    *x21        *x53
     8  +coeff( 17)    *x21        *x54
      pl_mnq_0_11 =pl_mnq_0_11 
     9  +coeff( 18)    *x23            
     1  +coeff( 19)    *x22    *x42    
     2  +coeff( 20)                *x54
     3  +coeff( 21)    *x22*x31*x41    
     4  +coeff( 22)            *x44    
     5  +coeff( 23)    *x22*x31*x43    
     6  +coeff( 24)    *x23        *x53
     7  +coeff( 25)            *x42*x51
     8  +coeff( 26)        *x31*x43    
      pl_mnq_0_11 =pl_mnq_0_11 
     9  +coeff( 27)            *x42*x52
     1  +coeff( 28)    *x23    *x42    
     2  +coeff( 29)    *x21    *x44    
     3  +coeff( 30)    *x24        *x51
     4  +coeff( 31)    *x22    *x42*x51
     5  +coeff( 32)    *x24    *x42    
     6  +coeff( 33)    *x23    *x42*x51
     7  +coeff( 34)    *x22        *x54
     8  +coeff( 35)    *x23*x33*x41    
      pl_mnq_0_11 =pl_mnq_0_11 
     9  +coeff( 36)    *x24    *x42*x51
     1  +coeff( 37)        *x33*x43*x51
     2  +coeff( 38)        *x32        
     3  +coeff( 39)    *x21*x31*x41    
     4  +coeff( 40)    *x21    *x42    
c
      return
      end
      function x_mnq_0_15  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.1558709E-01/
      data xmin/
     1 -0.49978E-02,-0.50007E-01,-0.10393E+00,-0.45852E-01,-0.39669E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49966E-02, 0.49912E-01, 0.10393E+00, 0.45852E-01, 0.49976E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.19046794E+00, 0.96728855E+00, 0.17587003E+01, 0.11623566E+00,
     + -0.27237956E-02,-0.28785808E-01,-0.87247200E-01, 0.17495342E+00,
     + -0.71572924E+00, 0.38714851E-02, 0.92975333E-01, 0.43791286E-02,
     + -0.10965247E+00,-0.30595604E+00, 0.15065094E-01,-0.19214209E-03,
     + -0.61688662E-01,-0.20436625E+00,-0.12921025E+00, 0.24041191E+00,
     + -0.17533334E-01,-0.83589827E-03,-0.12728466E+01, 0.87005710E-02,
     + -0.22219734E-01,-0.19515269E+00,-0.52365988E+00,-0.14541826E+01,
     + -0.38958594E-01,-0.18665660E-01,-0.11925712E+00,-0.32394877E+00,
     +  0.22378139E+00,-0.21365348E-01, 0.86169655E-03,-0.11291790E+01,
     + -0.44576347E+00,-0.11687914E+01, 0.68533130E-03, 0.23593675E+00,
     + -0.44059721E-02, 0.36096737E-01, 0.39487028E+00, 0.32439435E+00,
     +  0.11167950E+00, 0.54288667E+00, 0.59710224E-02,-0.15136286E-01,
     +  0.18316574E-01,-0.25713826E-01, 0.66178702E-01, 0.46077531E-01,
     +  0.39918482E+00, 0.75136200E-02, 0.13523684E+01, 0.14038195E+00,
     +  0.66217917E+00,-0.16525237E-01, 0.27947018E+01, 0.99101566E-01,
     +  0.77040702E-01,-0.35782307E-01, 0.35870454E+00, 0.15283810E+01,
     + -0.62875129E-01,-0.84439248E-01,-0.23380937E-01,-0.14733538E-01,
     +  0.14313908E+00,-0.40614102E-01, 0.56766335E-01, 0.12917416E+00,
     +  0.72141767E-01,-0.14349109E-01, 0.38730055E-01,-0.29223306E-01,
     + -0.36224523E-02,-0.43962830E+00,-0.27985388E-01,-0.20445313E-02,
     + -0.37546355E+00, 0.23901029E-01,-0.15782706E-02,-0.54314104E-02,
     +  0.16662423E-01, 0.64536324E-02, 0.21677636E-01, 0.19529020E-01,
     + -0.87329566E-01,-0.25509965E-01,
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
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_mnq_0_15  =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22            
     5  +coeff(  5)        *x32        
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21        *x51
      x_mnq_0_15  =x_mnq_0_15  
     9  +coeff(  9)                *x52
     1  +coeff( 10)*x11                
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x21*x32        
     4  +coeff( 13)    *x21*x31*x41    
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)    *x22        *x51
     7  +coeff( 16)        *x32    *x51
     8  +coeff( 17)        *x31*x41*x51
      x_mnq_0_15  =x_mnq_0_15  
     9  +coeff( 18)            *x42*x51
     1  +coeff( 19)    *x21        *x52
     2  +coeff( 20)                *x53
     3  +coeff( 21)    *x22*x32        
     4  +coeff( 22)        *x33*x41    
     5  +coeff( 23)    *x22    *x42    
     6  +coeff( 24)        *x32*x42    
     7  +coeff( 25)        *x31*x43    
     8  +coeff( 26)    *x23        *x51
      x_mnq_0_15  =x_mnq_0_15  
     9  +coeff( 27)    *x21*x31*x41*x51
     1  +coeff( 28)    *x21    *x42*x51
     2  +coeff( 29)    *x22        *x52
     3  +coeff( 30)        *x32    *x52
     4  +coeff( 31)        *x31*x41*x52
     5  +coeff( 32)            *x42*x52
     6  +coeff( 33)    *x21        *x53
     7  +coeff( 34)    *x23*x32        
     8  +coeff( 35)*x11        *x42    
      x_mnq_0_15  =x_mnq_0_15  
     9  +coeff( 36)    *x23    *x42    
     1  +coeff( 37)    *x22*x31*x41*x51
     2  +coeff( 38)    *x22    *x42*x51
     3  +coeff( 39)        *x32*x42*x51
     4  +coeff( 40)    *x23        *x52
     5  +coeff( 41)    *x21*x32    *x52
     6  +coeff( 42)    *x21*x31*x41*x52
     7  +coeff( 43)    *x21    *x42*x52
     8  +coeff( 44)    *x22        *x53
      x_mnq_0_15  =x_mnq_0_15  
     9  +coeff( 45)        *x31*x41*x53
     1  +coeff( 46)            *x42*x53
     2  +coeff( 47)*x11*x23            
     3  +coeff( 48)    *x22*x33*x41    
     4  +coeff( 49)*x11*x21    *x42    
     5  +coeff( 50)    *x22*x32*x42    
     6  +coeff( 51)    *x22*x31*x43    
     7  +coeff( 52)    *x23*x32    *x51
     8  +coeff( 53)    *x23*x31*x41*x51
      x_mnq_0_15  =x_mnq_0_15  
     9  +coeff( 54)    *x21*x33*x41*x51
     1  +coeff( 55)    *x23    *x42*x51
     2  +coeff( 56)    *x21*x31*x43*x51
     3  +coeff( 57)    *x22*x31*x41*x52
     4  +coeff( 58)        *x33*x41*x52
     5  +coeff( 59)    *x22    *x42*x52
     6  +coeff( 60)        *x31*x43*x52
     7  +coeff( 61)    *x23        *x53
     8  +coeff( 62)    *x21*x32    *x53
      x_mnq_0_15  =x_mnq_0_15  
     9  +coeff( 63)    *x21*x31*x41*x53
     1  +coeff( 64)    *x21    *x42*x53
     2  +coeff( 65)    *x22*x32*x42*x51
     3  +coeff( 66)    *x23*x31*x41*x52
     4  +coeff( 67)*x11        *x42*x52
     5  +coeff( 68)    *x22*x31*x41*x53
     6  +coeff( 69)    *x22*x33*x43    
     7  +coeff( 70)    *x23*x33*x41*x51
     8  +coeff( 71)*x11*x22    *x42*x51
      x_mnq_0_15  =x_mnq_0_15  
     9  +coeff( 72)    *x21*x33*x43*x51
     1  +coeff( 73)*x11*x23    *x42*x51
     2  +coeff( 74)    *x22*x32*x42*x53
     3  +coeff( 75)*x11*x22*x33*x41*x51
     4  +coeff( 76)*x11*x22*x33*x43    
     5  +coeff( 77)*x11*x21            
     6  +coeff( 78)    *x22*x31*x41    
     7  +coeff( 79)    *x21*x32    *x51
     8  +coeff( 80)*x11*x22            
      x_mnq_0_15  =x_mnq_0_15  
     9  +coeff( 81)    *x23*x31*x41    
     1  +coeff( 82)    *x21*x33*x41    
     2  +coeff( 83)    *x21*x32*x42    
     3  +coeff( 84)*x11*x21        *x51
     4  +coeff( 85)        *x33*x41*x51
     5  +coeff( 86)*x11*x21*x31*x41    
     6  +coeff( 87)        *x33*x43    
     7  +coeff( 88)*x11*x22    *x42    
     8  +coeff( 89)    *x23*x32*x42    
      x_mnq_0_15  =x_mnq_0_15  
     9  +coeff( 90)*x11*x21    *x42*x51
c
      return
      end
      function t_mnq_0_15  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.1426885E-01/
      data xmin/
     1 -0.49978E-02,-0.50007E-01,-0.10393E+00,-0.45852E-01,-0.39669E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49966E-02, 0.49912E-01, 0.10393E+00, 0.45852E-01, 0.49976E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.10272787E-01, 0.44507608E-01, 0.21012844E+00, 0.16099250E-01,
     + -0.26141846E-03,-0.27570073E-02,-0.10349214E-01, 0.23636589E-01,
     + -0.85604645E-01, 0.14517169E-01, 0.11343396E-02,-0.14364465E-01,
     + -0.44276033E-01, 0.83774263E-02, 0.30223877E-03,-0.97328015E-02,
     + -0.30673599E-01,-0.13896993E-01, 0.29224180E-01,-0.20038665E-02,
     + -0.76292275E-03,-0.18260425E+00, 0.35509744E-03,-0.44324193E-02,
     + -0.23127770E-01,-0.76227650E-01,-0.21709785E+00,-0.66832229E-02,
     + -0.28612469E-02,-0.18099831E-01,-0.51225662E-01, 0.26446585E-01,
     + -0.14172134E-02,-0.84629057E-04,-0.39752475E-02,-0.16378522E+00,
     + -0.62497724E-01,-0.20007561E+00,-0.30894813E-02, 0.30346610E-01,
     +  0.12443452E-02, 0.22986748E-02, 0.30298209E-01, 0.42453427E-01,
     +  0.17739842E-01, 0.72248451E-01,-0.56373206E-03, 0.41156914E-02,
     + -0.59670494E-02, 0.14333377E-01, 0.50414298E-01, 0.21147050E-02,
     +  0.27066311E-01, 0.10122740E+00,-0.30751836E-02, 0.34748337E+00,
     +  0.10326186E-01, 0.56677584E-01, 0.20151642E+00, 0.70590249E-04,
     + -0.11415604E-03,-0.63508995E-01,-0.47548665E-02,-0.52247480E-01,
     + -0.55447337E-02,-0.21163383E-02,-0.10314385E-03,-0.90902526E-03,
     +  0.46434769E-03,-0.14310905E-02, 0.41279085E-02, 0.10090765E-01,
     +  0.14997648E+00, 0.80509614E-02, 0.15922330E-01,-0.32138100E-03,
     +  0.25293047E-02,-0.42738431E-03,-0.32192809E-02, 0.16874050E-02,
     + -0.73788641E-03, 0.34329016E-02, 0.63227350E-03,-0.23574347E-03,
     + -0.20217423E-02,-0.11367488E-03,-0.58302854E-03,-0.87217864E-03,
     +  0.10244570E-03,-0.34404616E-03,
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
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_mnq_0_15  =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22            
     5  +coeff(  5)        *x32        
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21        *x51
      t_mnq_0_15  =t_mnq_0_15  
     9  +coeff(  9)                *x52
     1  +coeff( 10)    *x23            
     2  +coeff( 11)    *x21*x32        
     3  +coeff( 12)    *x21*x31*x41    
     4  +coeff( 13)    *x21    *x42    
     5  +coeff( 14)    *x22        *x51
     6  +coeff( 15)        *x32    *x51
     7  +coeff( 16)        *x31*x41*x51
     8  +coeff( 17)            *x42*x51
      t_mnq_0_15  =t_mnq_0_15  
     9  +coeff( 18)    *x21        *x52
     1  +coeff( 19)                *x53
     2  +coeff( 20)    *x22*x32        
     3  +coeff( 21)        *x33*x41    
     4  +coeff( 22)    *x22    *x42    
     5  +coeff( 23)        *x32*x42    
     6  +coeff( 24)        *x31*x43    
     7  +coeff( 25)    *x23        *x51
     8  +coeff( 26)    *x21*x31*x41*x51
      t_mnq_0_15  =t_mnq_0_15  
     9  +coeff( 27)    *x21    *x42*x51
     1  +coeff( 28)    *x22        *x52
     2  +coeff( 29)        *x32    *x52
     3  +coeff( 30)        *x31*x41*x52
     4  +coeff( 31)            *x42*x52
     5  +coeff( 32)    *x21        *x53
     6  +coeff( 33)*x11*x22            
     7  +coeff( 34)*x11    *x32        
     8  +coeff( 35)    *x23*x32        
      t_mnq_0_15  =t_mnq_0_15  
     9  +coeff( 36)    *x23    *x42    
     1  +coeff( 37)    *x22*x31*x41*x51
     2  +coeff( 38)    *x22    *x42*x51
     3  +coeff( 39)        *x32*x42*x51
     4  +coeff( 40)    *x23        *x52
     5  +coeff( 41)    *x21*x32    *x52
     6  +coeff( 42)    *x21*x31*x41*x52
     7  +coeff( 43)    *x21    *x42*x52
     8  +coeff( 44)    *x22        *x53
      t_mnq_0_15  =t_mnq_0_15  
     9  +coeff( 45)        *x31*x41*x53
     1  +coeff( 46)            *x42*x53
     2  +coeff( 47)*x11*x23            
     3  +coeff( 48)    *x22*x33*x41    
     4  +coeff( 49)    *x22*x32*x42    
     5  +coeff( 50)    *x22*x31*x43    
     6  +coeff( 51)    *x23*x31*x41*x51
     7  +coeff( 52)    *x21*x33*x41*x51
     8  +coeff( 53)    *x21*x31*x43*x51
      t_mnq_0_15  =t_mnq_0_15  
     9  +coeff( 54)    *x22*x31*x41*x52
     1  +coeff( 55)        *x33*x41*x52
     2  +coeff( 56)    *x22    *x42*x52
     3  +coeff( 57)    *x23        *x53
     4  +coeff( 58)    *x21*x31*x41*x53
     5  +coeff( 59)    *x21    *x42*x53
     6  +coeff( 60)*x11                
     7  +coeff( 61)*x11*x21            
     8  +coeff( 62)    *x22*x31*x41    
      t_mnq_0_15  =t_mnq_0_15  
     9  +coeff( 63)    *x21*x32    *x51
     1  +coeff( 64)    *x23*x31*x41    
     2  +coeff( 65)    *x21*x32*x42    
     3  +coeff( 66)*x11*x21        *x51
     4  +coeff( 67)    *x22*x32    *x51
     5  +coeff( 68)*x11            *x52
     6  +coeff( 69)        *x32    *x53
     7  +coeff( 70)*x11*x21    *x42    
     8  +coeff( 71)        *x33*x43    
      t_mnq_0_15  =t_mnq_0_15  
     9  +coeff( 72)    *x23*x32    *x51
     1  +coeff( 73)    *x23    *x42*x51
     2  +coeff( 74)    *x22*x32    *x52
     3  +coeff( 75)        *x31*x43*x52
     4  +coeff( 76)*x11    *x31*x41    
     5  +coeff( 77)    *x21*x33*x41    
     6  +coeff( 78)*x11        *x42    
     7  +coeff( 79)    *x21*x31*x43    
     8  +coeff( 80)        *x33*x41*x51
      t_mnq_0_15  =t_mnq_0_15  
     9  +coeff( 81)*x11        *x42*x51
     1  +coeff( 82)        *x32*x42*x52
     2  +coeff( 83)*x11            *x53
     3  +coeff( 84)*x11            *x51
     4  +coeff( 85)        *x31*x43*x51
     5  +coeff( 86)*x11*x21*x32        
     6  +coeff( 87)*x11*x21*x31*x41    
     7  +coeff( 88)*x11*x22        *x51
     8  +coeff( 89)*x11    *x32    *x51
      t_mnq_0_15  =t_mnq_0_15  
     9  +coeff( 90)*x11    *x31*x41*x51
c
      return
      end
      function y_mnq_0_15  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49978E-02,-0.50007E-01,-0.10393E+00,-0.45852E-01,-0.39669E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49966E-02, 0.49912E-01, 0.10393E+00, 0.45852E-01, 0.49976E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.83141223E-01, 0.83603787E+00,-0.26393831E-01,-0.14886977E+00,
     +  0.10733387E-01, 0.12179423E-01,-0.60958754E-01,-0.72252268E+00,
     +  0.30399725E-03, 0.73222308E-02,-0.70325655E+00,-0.12553912E-01,
     + -0.22376524E+00,-0.61757181E-01,-0.81236553E+00, 0.18790280E-03,
     +  0.18401366E+00, 0.51126373E+00,-0.59229000E-02,-0.11095666E+01,
     +  0.27868956E+00, 0.16488489E-01,-0.40309218E+00,-0.12154061E-01,
     + -0.38402283E-03,-0.62485440E-02,-0.27416114E-01, 0.51321900E+00,
     +  0.49690478E-02, 0.20008104E+01, 0.51060035E-02, 0.82557350E-02,
     +  0.42713568E-01,-0.73128410E-01, 0.19101618E+00, 0.20003183E+01,
     + -0.50235377E-02, 0.22028873E-01, 0.42745063E+00,-0.93991041E-01,
     + -0.13443456E-01, 0.54511744E+00,-0.15057600E-02, 0.22428069E+01,
     + -0.86619975E-02,-0.19566019E+00, 0.11238361E+00, 0.30295892E+01,
     + -0.54033078E-01, 0.17077476E-02,-0.20816651E+00,-0.53814662E+00,
     +  0.85565716E+00,-0.45242731E-01,-0.14741497E+00, 0.47874186E-01,
     +  0.10944167E+00,-0.38147084E-01,-0.18410647E-01,-0.11346550E-01,
     + -0.16440284E+00,-0.43241170E+00,-0.17028673E+00,-0.39825347E+00,
     + -0.54310572E+00,-0.14341553E-01, 0.75588569E-01,-0.59919227E-01,
     +  0.15659679E+00, 0.47969496E+00,-0.30451411E+00,-0.67996155E-02,
     +  0.16140284E-01, 0.56631707E-01, 0.20097193E+00, 0.45523497E+00,
     +  0.52532412E-01, 0.10900761E+00,-0.74510761E-01,-0.21462968E-01,
     +  0.93824401E-01,-0.30043229E+00,-0.28486171E+00, 0.22815156E+00,
     + -0.25641998E-01,-0.55323583E+00, 0.23738839E-01, 0.23592299E+00,
     + -0.15999389E-04, 0.10783287E-01,
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
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      y_mnq_0_15  =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21*x31        
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)        *x31    *x51
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)    *x22*x31        
     8  +coeff(  8)    *x22    *x41    
      y_mnq_0_15  =y_mnq_0_15  
     9  +coeff(  9)        *x32*x41    
     1  +coeff( 10)            *x43    
     2  +coeff( 11)    *x21    *x41*x51
     3  +coeff( 12)        *x31    *x52
     4  +coeff( 13)            *x41*x52
     5  +coeff( 14)    *x23*x31        
     6  +coeff( 15)    *x23    *x41    
     7  +coeff( 16)    *x21*x32*x41    
     8  +coeff( 17)    *x21*x31*x42    
      y_mnq_0_15  =y_mnq_0_15  
     9  +coeff( 18)    *x21    *x43    
     1  +coeff( 19)    *x22*x31    *x51
     2  +coeff( 20)    *x22    *x41*x51
     3  +coeff( 21)            *x43*x51
     4  +coeff( 22)    *x21*x31    *x52
     5  +coeff( 23)    *x21    *x41*x52
     6  +coeff( 24)        *x31    *x53
     7  +coeff( 25)*x11*x21*x31        
     8  +coeff( 26)*x11*x21    *x41    
      y_mnq_0_15  =y_mnq_0_15  
     9  +coeff( 27)    *x22*x32*x41    
     1  +coeff( 28)    *x22*x31*x42    
     2  +coeff( 29)        *x33*x42    
     3  +coeff( 30)    *x22    *x43    
     4  +coeff( 31)        *x32*x43    
     5  +coeff( 32)    *x23*x31    *x51
     6  +coeff( 33)    *x23    *x41*x51
     7  +coeff( 34)    *x21*x32*x41*x51
     8  +coeff( 35)    *x21*x31*x42*x51
      y_mnq_0_15  =y_mnq_0_15  
     9  +coeff( 36)    *x21    *x43*x51
     1  +coeff( 37)        *x33    *x52
     2  +coeff( 38)        *x31*x42*x52
     3  +coeff( 39)            *x43*x52
     4  +coeff( 40)    *x21*x31    *x53
     5  +coeff( 41)*x11*x22    *x41    
     6  +coeff( 42)    *x23*x31*x42    
     7  +coeff( 43)*x11        *x43    
     8  +coeff( 44)    *x23    *x43    
      y_mnq_0_15  =y_mnq_0_15  
     9  +coeff( 45)*x11*x21    *x41*x51
     1  +coeff( 46)    *x22*x32*x41*x51
     2  +coeff( 47)        *x33*x42*x51
     3  +coeff( 48)    *x22    *x43*x51
     4  +coeff( 49)    *x23*x31    *x52
     5  +coeff( 50)    *x21*x33    *x52
     6  +coeff( 51)    *x23    *x41*x52
     7  +coeff( 52)    *x21*x31*x42*x52
     8  +coeff( 53)    *x21    *x43*x52
      y_mnq_0_15  =y_mnq_0_15  
     9  +coeff( 54)    *x22*x31    *x53
     1  +coeff( 55)    *x22    *x41*x53
     2  +coeff( 56)        *x32*x41*x53
     3  +coeff( 57)        *x31*x42*x53
     4  +coeff( 58)    *x22*x33*x42    
     5  +coeff( 59)    *x23*x33    *x51
     6  +coeff( 60)*x11*x22    *x41*x51
     7  +coeff( 61)    *x23*x32*x41*x51
     8  +coeff( 62)    *x23*x31*x42*x51
      y_mnq_0_15  =y_mnq_0_15  
     9  +coeff( 63)    *x23    *x43*x51
     1  +coeff( 64)    *x22*x31*x42*x52
     2  +coeff( 65)    *x22    *x43*x52
     3  +coeff( 66)        *x32*x43*x52
     4  +coeff( 67)    *x23*x31    *x53
     5  +coeff( 68)    *x23    *x41*x53
     6  +coeff( 69)    *x21*x32*x41*x53
     7  +coeff( 70)    *x21*x31*x42*x53
     8  +coeff( 71)    *x21    *x43*x53
      y_mnq_0_15  =y_mnq_0_15  
     9  +coeff( 72)*x11*x23    *x41*x51
     1  +coeff( 73)*x11*x21*x32*x41*x51
     2  +coeff( 74)    *x23*x33    *x52
     3  +coeff( 75)    *x23*x31*x42*x52
     4  +coeff( 76)    *x21*x33*x42*x52
     5  +coeff( 77)    *x22*x33    *x53
     6  +coeff( 78)    *x22*x32*x41*x53
     7  +coeff( 79)        *x33*x42*x53
     8  +coeff( 80)*x11*x23*x32*x41    
      y_mnq_0_15  =y_mnq_0_15  
     9  +coeff( 81)    *x23*x33*x42*x51
     1  +coeff( 82)    *x23*x31*x42*x53
     2  +coeff( 83)    *x21*x33*x42*x53
     3  +coeff( 84)    *x23    *x43*x53
     4  +coeff( 85)*x11*x23    *x43*x51
     5  +coeff( 86)    *x23*x33*x42*x52
     6  +coeff( 87)*x11*x23    *x43*x52
     7  +coeff( 88)    *x23*x33*x42*x53
     8  +coeff( 89)        *x33        
      y_mnq_0_15  =y_mnq_0_15  
     9  +coeff( 90)        *x31*x42    
c
      return
      end
      function p_mnq_0_15  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49978E-02,-0.50007E-01,-0.10393E+00,-0.45852E-01,-0.39669E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49966E-02, 0.49912E-01, 0.10393E+00, 0.45852E-01, 0.49976E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.16644489E-02, 0.32341130E-01,-0.23773074E-01, 0.20506042E-02,
     + -0.37811871E-02,-0.15106282E-01,-0.11793287E-02,-0.96225649E-01,
     + -0.12841923E-03, 0.29868556E-02,-0.14712608E-01,-0.11229252E+00,
     + -0.62493519E-02,-0.35611421E-01,-0.81248058E-04,-0.19352542E-01,
     + -0.64445278E-02,-0.19533167E-03,-0.10167070E+00, 0.33048277E-02,
     +  0.35984967E-01, 0.76526068E-01,-0.26667958E-01,-0.45991205E-02,
     + -0.13949256E+00, 0.14316881E-02, 0.24610559E-01, 0.52899517E-01,
     + -0.15319862E-01,-0.41841600E-01,-0.12635755E-02, 0.10586669E-01,
     + -0.88429492E-03, 0.27879467E-02,-0.55600930E-03, 0.70427526E-02,
     +  0.91300733E-01, 0.40213028E-02, 0.24751662E+00, 0.16296060E-02,
     + -0.19587245E-03, 0.14979450E-02, 0.16182859E-02, 0.22260848E-01,
     +  0.40793787E-02, 0.90187095E-01, 0.29636148E+00, 0.36573627E-02,
     +  0.14604999E-02, 0.34891471E-01, 0.18182082E-01, 0.72525546E-01,
     +  0.20013943E-02, 0.20744642E-01,-0.16960481E-02, 0.31605083E-03,
     +  0.13384633E-01, 0.56826574E-03, 0.85088002E-04, 0.48352577E-01,
     +  0.73263841E-02, 0.24162515E+00, 0.42707454E-02,-0.29117821E-02,
     +  0.28421568E-01, 0.49527462E-04,-0.72146626E-02, 0.80211749E-02,
     +  0.33626422E+00, 0.51393919E-02,-0.12038630E-02,-0.40216981E-02,
     +  0.27855974E-01,-0.84125903E-02,-0.73129371E-01, 0.76514795E-01,
     + -0.45944527E-02, 0.90406677E-02,-0.61474121E-02, 0.20121622E-02,
     + -0.35766233E-01,-0.23279749E-01,-0.51476358E-03, 0.11019282E-02,
     +  0.96927362E-03, 0.82318566E-03,-0.87493649E-02, 0.51735691E-03,
     + -0.12380993E-02, 0.53154887E-03,
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
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_mnq_0_15  =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)        *x31    *x51
     5  +coeff(  5)            *x41*x51
     6  +coeff(  6)    *x22*x31        
     7  +coeff(  7)        *x33        
     8  +coeff(  8)    *x22    *x41    
      p_mnq_0_15  =p_mnq_0_15  
     9  +coeff(  9)        *x32*x41    
     1  +coeff( 10)            *x43    
     2  +coeff( 11)    *x21*x31    *x51
     3  +coeff( 12)    *x21    *x41*x51
     4  +coeff( 13)        *x31    *x52
     5  +coeff( 14)            *x41*x52
     6  +coeff( 15)*x11    *x31        
     7  +coeff( 16)    *x23*x31        
     8  +coeff( 17)    *x21*x33        
      p_mnq_0_15  =p_mnq_0_15  
     9  +coeff( 18)*x11        *x41    
     1  +coeff( 19)    *x23    *x41    
     2  +coeff( 20)    *x21*x32*x41    
     3  +coeff( 21)    *x21*x31*x42    
     4  +coeff( 22)    *x21    *x43    
     5  +coeff( 23)    *x22*x31    *x51
     6  +coeff( 24)        *x33    *x51
     7  +coeff( 25)    *x22    *x41*x51
     8  +coeff( 26)        *x32*x41*x51
      p_mnq_0_15  =p_mnq_0_15  
     9  +coeff( 27)        *x31*x42*x51
     1  +coeff( 28)            *x43*x51
     2  +coeff( 29)    *x21*x31    *x52
     3  +coeff( 30)    *x21    *x41*x52
     4  +coeff( 31)        *x31    *x53
     5  +coeff( 32)            *x41*x53
     6  +coeff( 33)*x11*x21*x31        
     7  +coeff( 34)    *x22*x33        
     8  +coeff( 35)*x11*x21    *x41    
      p_mnq_0_15  =p_mnq_0_15  
     9  +coeff( 36)    *x22*x32*x41    
     1  +coeff( 37)    *x22*x31*x42    
     2  +coeff( 38)        *x33*x42    
     3  +coeff( 39)    *x22    *x43    
     4  +coeff( 40)        *x32*x43    
     5  +coeff( 41)*x11    *x31    *x51
     6  +coeff( 42)    *x23*x31    *x51
     7  +coeff( 43)    *x21*x33    *x51
     8  +coeff( 44)    *x23    *x41*x51
      p_mnq_0_15  =p_mnq_0_15  
     9  +coeff( 45)    *x21*x32*x41*x51
     1  +coeff( 46)    *x21*x31*x42*x51
     2  +coeff( 47)    *x21    *x43*x51
     3  +coeff( 48)    *x22*x31    *x52
     4  +coeff( 49)        *x33    *x52
     5  +coeff( 50)    *x22    *x41*x52
     6  +coeff( 51)        *x31*x42*x52
     7  +coeff( 52)            *x43*x52
     8  +coeff( 53)    *x21*x31    *x53
      p_mnq_0_15  =p_mnq_0_15  
     9  +coeff( 54)    *x21    *x41*x53
     1  +coeff( 55)*x11*x22*x31        
     2  +coeff( 56)*x11    *x33        
     3  +coeff( 57)    *x23*x33        
     4  +coeff( 58)*x11*x22    *x41    
     5  +coeff( 59)*x11    *x31*x42    
     6  +coeff( 60)    *x23*x31*x42    
     7  +coeff( 61)    *x21*x33*x42    
     8  +coeff( 62)    *x23    *x43    
      p_mnq_0_15  =p_mnq_0_15  
     9  +coeff( 63)    *x21*x32*x43    
     1  +coeff( 64)*x11*x21*x31    *x51
     2  +coeff( 65)    *x22*x33    *x51
     3  +coeff( 66)*x11*x21    *x41*x51
     4  +coeff( 67)    *x22*x32*x41*x51
     5  +coeff( 68)        *x33*x42*x51
     6  +coeff( 69)    *x22    *x43*x51
     7  +coeff( 70)        *x32*x43*x51
     8  +coeff( 71)*x11    *x31    *x52
      p_mnq_0_15  =p_mnq_0_15  
     9  +coeff( 72)    *x23*x31    *x52
     1  +coeff( 73)    *x21*x33    *x52
     2  +coeff( 74)    *x23    *x41*x52
     3  +coeff( 75)    *x21*x31*x42*x52
     4  +coeff( 76)    *x21    *x43*x52
     5  +coeff( 77)    *x22*x31    *x53
     6  +coeff( 78)        *x33    *x53
     7  +coeff( 79)    *x22    *x41*x53
     8  +coeff( 80)        *x32*x41*x53
      p_mnq_0_15  =p_mnq_0_15  
     9  +coeff( 81)        *x31*x42*x53
     1  +coeff( 82)            *x43*x53
     2  +coeff( 83)*x11*x23*x31        
     3  +coeff( 84)*x11*x21*x33        
     4  +coeff( 85)*x11*x23    *x41    
     5  +coeff( 86)*x11*x21*x31*x42    
     6  +coeff( 87)    *x22*x33*x42    
     7  +coeff( 88)*x11*x21    *x43    
     8  +coeff( 89)*x11*x22*x31    *x51
      p_mnq_0_15  =p_mnq_0_15  
     9  +coeff( 90)*x11    *x33    *x51
c
      return
      end
      function pl_mnq_0_15 (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 41)
      data ncoeff/ 40/
      data avdat/  0.1021436E-02/
      data xmin/
     1 -0.49978E-02,-0.50007E-01,-0.10393E+00,-0.45852E-01,-0.39669E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49966E-02, 0.49912E-01, 0.10393E+00, 0.45852E-01, 0.49976E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.20697575E-01,-0.32348025E+00,-0.18951699E+00,-0.90061918E-01,
     + -0.78563258E-01,-0.63795671E-02, 0.51746331E-01, 0.39885771E-01,
     + -0.99898977E-02, 0.10884161E+00,-0.26055055E-02,-0.36324315E-01,
     + -0.58916062E-01,-0.12379142E+00, 0.30269283E-02, 0.23257975E-01,
     +  0.79863090E-02, 0.31957980E-01,-0.43594092E-02, 0.36023591E-01,
     +  0.58758419E-01,-0.14291805E-01,-0.64858705E-01,-0.14676929E+00,
     + -0.13311657E+00,-0.10054222E-01, 0.10855325E-01, 0.35060797E-01,
     +  0.11950491E-01,-0.11611378E-01, 0.15820126E-02,-0.19729491E-01,
     + -0.55758607E-01, 0.15689231E-01, 0.11032579E-01,-0.18581783E-02,
     + -0.77256272E-02,-0.88194488E-02,-0.19246707E-01,-0.70577851E-02,
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
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
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
      x54 = x53*x5
c
c                  function
c
      pl_mnq_0_15 =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x21        *x51
     5  +coeff(  5)                *x52
     6  +coeff(  6)    *x22        *x51
     7  +coeff(  7)    *x21        *x52
     8  +coeff(  8)    *x21        *x53
      pl_mnq_0_15 =pl_mnq_0_15 
     9  +coeff(  9)            *x42    
     1  +coeff( 10)                *x53
     2  +coeff( 11)*x11                
     3  +coeff( 12)    *x22            
     4  +coeff( 13)                *x54
     5  +coeff( 14)    *x21        *x54
     6  +coeff( 15)        *x31*x41    
     7  +coeff( 16)    *x21    *x42    
     8  +coeff( 17)    *x21*x31*x41    
      pl_mnq_0_15 =pl_mnq_0_15 
     9  +coeff( 18)            *x42*x51
     1  +coeff( 19)    *x24            
     2  +coeff( 20)    *x21    *x42*x51
     3  +coeff( 21)    *x22        *x52
     4  +coeff( 22)    *x23    *x42    
     5  +coeff( 23)    *x22        *x53
     6  +coeff( 24)    *x23        *x53
     7  +coeff( 25)    *x22        *x54
     8  +coeff( 26)    *x23            
      pl_mnq_0_15 =pl_mnq_0_15 
     9  +coeff( 27)        *x31*x41*x51
     1  +coeff( 28)    *x22    *x42    
     2  +coeff( 29)    *x21*x31*x41*x51
     3  +coeff( 30)    *x23*x31*x41    
     4  +coeff( 31)    *x22*x32    *x51
     5  +coeff( 32)    *x24    *x42    
     6  +coeff( 33)    *x24        *x52
     7  +coeff( 34)    *x24    *x42*x51
     8  +coeff( 35)    *x22*x31*x41    
      pl_mnq_0_15 =pl_mnq_0_15 
     9  +coeff( 36)            *x44    
     1  +coeff( 37)    *x21    *x44    
     2  +coeff( 38)        *x31*x43*x51
     3  +coeff( 39)            *x44*x51
     4  +coeff( 40)    *x24*x31*x41    
c
      return
      end

c***********************************************************************************   
c***********************************************************************************       
c Reverse functions for MAD with 35 degree bend angle and    
c 12 degree minimum scattering angle quads off    
c                                              -JJL 12/13/04    
c    
c re usage:    
c There are 4 fortran functions:  (TXFIT is not needed in the no quads tune)  
c    
c delta_nq(x,m) for determining delta   
c theta_nq(x,m)  for determining theta_0 (vertical angle at the target) 
c phi_nq(x,m) For determining phi_0 (horizontal angle at the target)   
c y00_nq(x,m) For determining y_0 (transverse/horizontal position perpendicular to the optic axis at the target)   
    
c Calling convention is as follows. x is an array in which x(1->5) are xf, thetaf, yf, phif    
c and x0 (the detector coordinates in the transport coordinate system and the vertical beam    
c position at the target). m=5.    
c    
c Procedure:    
c I think it's self-explanitory:    
c    
c delta0 = delta_nq(x,5)    
c theta0= theta_nq(x,5)    
c phi0=phi_nq(x,5)    
c y0=y00_nq(x,5)

      function    delta_nq (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 96)
      data ncoeff/ 95/
      data avdat/  0.9984455E-01/
      data xmin/
     1 -0.71654E+00,-0.30874E+00,-0.62118E+00,-0.39822E-01,-0.49978E-02,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.68856E+00, 0.16818E+00, 0.62118E+00, 0.39822E-01, 0.49966E-02,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.24973683E+00,-0.90310745E-01, 0.45132911E+00,-0.51132999E-01,
     +  0.23940196E-01, 0.18808697E+00, 0.42354099E-01,-0.22699319E-01,
     +  0.65827444E-02, 0.64186798E-03,-0.15874464E-01, 0.53739525E-01,
     + -0.10100927E+00, 0.12750185E+00, 0.23409103E+00,-0.19987890E+00,
     + -0.55673119E-01, 0.72636724E-01,-0.25488378E-01, 0.99350465E-03,
     +  0.46198678E-03,-0.12459370E-03, 0.96522402E-02,-0.47403537E-01,
     +  0.81462435E-01,-0.13734142E+00, 0.11406457E+00, 0.62736905E+00,
     + -0.97129059E+00, 0.43361372E+00,-0.19772579E+00, 0.30860886E+00,
     + -0.12063400E+00,-0.90623572E-02, 0.14540930E+00, 0.85961536E-01,
     +  0.42256672E-01, 0.56386924E+00,-0.93460417E+00, 0.62971354E+00,
     + -0.20587285E+00,-0.32746622E+00, 0.49089850E-02, 0.19685830E+00,
     + -0.27692392E+00,-0.17092253E+00, 0.78649569E+00, 0.13064294E+00,
     + -0.23923454E+00, 0.32149887E+00,-0.24005215E+00,-0.10069960E+00,
     +  0.98086055E-03,-0.23722136E+00,-0.21816803E-01,-0.49910862E-02,
     +  0.25104927E-01,-0.50647402E+00, 0.66649369E-02, 0.38177233E-01,
     + -0.29712611E+00, 0.35072437E-02,-0.24833722E+00, 0.23689678E+00,
     +  0.83451547E-01, 0.41190851E+00,-0.70220846E+00, 0.97765468E-01,
     + -0.27155709E-01,-0.12406183E-01,-0.11105172E+00, 0.27028292E-01,
     +  0.58516290E-01, 0.44382017E-03,-0.88283662E-02, 0.38411492E-02,
     + -0.75448960E-01, 0.86627253E-04, 0.19308396E-02, 0.30233299E-01,
     + -0.37556019E-01, 0.45163465E+00,-0.65084167E-01, 0.17691420E-02,
     +  0.10872652E+00,-0.22533554E+00, 0.31229785E+00,-0.15964693E+00,
     + -0.29291697E-01, 0.21730296E-01, 0.52990213E-01,-0.54513026E-01,
     +  0.12300364E+00, 0.54219848E+00,-0.29376468E+00,
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
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x13 = x12*x1
      x14 = x13*x1
      x15 = x14*x1
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
c
c                  function
c
         delta_nq =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)*x11                
     3  +coeff(  3)    *x21            
     4  +coeff(  4)*x12                
     5  +coeff(  5)*x11*x21            
     6  +coeff(  6)    *x22            
     7  +coeff(  7)        *x32        
     8  +coeff(  8)        *x31*x41    
         delta_nq =   delta_nq 
     9  +coeff(  9)            *x42    
     1  +coeff( 10)                *x51
     2  +coeff( 11)*x13                
     3  +coeff( 12)*x12*x21            
     4  +coeff( 13)*x11*x22            
     5  +coeff( 14)    *x23            
     6  +coeff( 15)*x11    *x32        
     7  +coeff( 16)    *x21*x32        
     8  +coeff( 17)*x11    *x31*x41    
         delta_nq =   delta_nq 
     9  +coeff( 18)    *x21*x31*x41    
     1  +coeff( 19)*x11        *x42    
     2  +coeff( 20)    *x21    *x42    
     3  +coeff( 21)*x11            *x51
     4  +coeff( 22)    *x21        *x51
     5  +coeff( 23)*x14                
     6  +coeff( 24)*x13*x21            
     7  +coeff( 25)*x12*x22            
     8  +coeff( 26)*x11*x23            
         delta_nq =   delta_nq 
     9  +coeff( 27)    *x24            
     1  +coeff( 28)*x12    *x32        
     2  +coeff( 29)*x11*x21*x32        
     3  +coeff( 30)    *x22*x32        
     4  +coeff( 31)*x12    *x31*x41    
     5  +coeff( 32)*x11*x21*x31*x41    
     6  +coeff( 33)    *x22*x31*x41    
     7  +coeff( 34)*x15                
     8  +coeff( 35)*x14*x21            
         delta_nq =   delta_nq 
     9  +coeff( 36)*x11*x24            
     1  +coeff( 37)    *x25            
     2  +coeff( 38)*x13    *x32        
     3  +coeff( 39)*x12*x21*x32        
     4  +coeff( 40)*x11*x22*x32        
     5  +coeff( 41)    *x23*x32        
     6  +coeff( 42)*x13    *x31*x41    
     7  +coeff( 43)*x12*x21        *x51
     8  +coeff( 44)*x14    *x32        
         delta_nq =   delta_nq 
     9  +coeff( 45)*x13*x21*x32        
     1  +coeff( 46)*x14    *x31*x41    
     2  +coeff( 47)*x13*x21*x31*x41    
     3  +coeff( 48)*x12    *x31*x43    
     4  +coeff( 49)*x14*x22*x32        
     5  +coeff( 50)*x12*x24*x32        
     6  +coeff( 51)*x13*x23*x31*x41    
     7  +coeff( 52)*x13*x21*x33*x41    
     8  +coeff( 53)*x12            *x51
         delta_nq =   delta_nq 
     9  +coeff( 54)*x13*x22            
     1  +coeff( 55)*x13        *x42    
     2  +coeff( 56)*x11*x22        *x51
     3  +coeff( 57)*x11    *x34        
     4  +coeff( 58)*x12*x22*x31*x41    
     5  +coeff( 59)    *x24    *x42    
     6  +coeff( 60)*x15*x22            
     7  +coeff( 61)*x14*x23            
     8  +coeff( 62)*x11*x21*x32    *x51
         delta_nq =   delta_nq 
     9  +coeff( 63)*x14*x21*x32        
     1  +coeff( 64)*x12*x23*x32        
     2  +coeff( 65)*x15    *x31*x41    
     3  +coeff( 66)*x14*x21*x31*x41    
     4  +coeff( 67)*x13*x22*x31*x41    
     5  +coeff( 68)*x14*x21    *x42    
     6  +coeff( 69)*x14*x24            
     7  +coeff( 70)*x12*x21*x34        
     8  +coeff( 71)*x11*x22*x32*x42    
         delta_nq =   delta_nq 
     9  +coeff( 72)*x13*x23    *x42    
     1  +coeff( 73)    *x25*x31*x43    
     2  +coeff( 74)        *x34*x44    
     3  +coeff( 75)    *x25    *x42*x51
     4  +coeff( 76)    *x23*x34    *x51
     5  +coeff( 77)*x12        *x42    
     6  +coeff( 78)    *x22        *x51
     7  +coeff( 79)        *x34        
     8  +coeff( 80)        *x33*x41    
         delta_nq =   delta_nq 
     9  +coeff( 81)        *x32*x42    
     1  +coeff( 82)*x12*x21*x31*x41    
     2  +coeff( 83)*x11*x22*x31*x41    
     3  +coeff( 84)    *x23        *x51
     4  +coeff( 85)*x15*x21            
     5  +coeff( 86)*x14*x22            
     6  +coeff( 87)*x12*x24            
     7  +coeff( 88)*x11*x25            
     8  +coeff( 89)    *x21*x34        
         delta_nq =   delta_nq 
     9  +coeff( 90)    *x21*x33*x41    
     1  +coeff( 91)*x11    *x32*x42    
     2  +coeff( 92)*x11    *x31*x43    
     3  +coeff( 93)*x12*x22    *x42    
     4  +coeff( 94)*x13*x24            
     5  +coeff( 95)*x12*x25            
c
      return
      end
      function    theta_nq (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 96)
      data ncoeff/ 95/
      data avdat/ -0.8959717E-02/
      data xmin/
     1 -0.71654E+00,-0.30874E+00,-0.62118E+00,-0.39822E-01,-0.49978E-02,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.68856E+00, 0.16818E+00, 0.62118E+00, 0.39822E-01, 0.49966E-02,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.33237368E-01, 0.58816697E-01,-0.85582443E-01,-0.12384327E-03,
     +  0.49463417E-02,-0.69615282E-02, 0.14350672E-02,-0.38600997E-02,
     +  0.17227845E-02,-0.42350401E-03,-0.29934535E-02, 0.13918588E-01,
     + -0.22770699E-01, 0.12429207E-01,-0.85917348E-02, 0.13899537E-02,
     +  0.40890677E-02,-0.19569052E-02,-0.18677308E-02,-0.44150930E-02,
     +  0.20539783E-01,-0.34266345E-01, 0.26508830E-01,-0.27036671E-01,
     +  0.21404466E-01, 0.18210032E-02, 0.49794186E-02,-0.23728954E-02,
     + -0.96519572E-04, 0.89834780E-02, 0.26773436E-01, 0.46676959E-03,
     + -0.13521871E-01,-0.12555291E-01, 0.52383570E-02,-0.57483464E-02,
     +  0.29326396E-02,-0.41116737E-02,-0.27221482E-01, 0.16092778E-02,
     + -0.84380833E-02,-0.49733283E-03,-0.49873791E-02, 0.10407210E-02,
     + -0.13779785E-01, 0.74374741E-02,-0.46897572E-02, 0.21769732E-01,
     +  0.74133752E-02,-0.11420131E-03,-0.80722673E-02, 0.89898901E-02,
     + -0.49323426E-02, 0.81970713E-04,-0.43649483E-03,-0.11345570E-02,
     + -0.11739567E-02, 0.78760674E-02,-0.15421468E-01, 0.13072611E-01,
     +  0.12590858E-02, 0.18610405E-02,-0.15570632E-02, 0.94386153E-02,
     + -0.12528563E-03, 0.25810514E-03, 0.45136795E-02, 0.14361700E-02,
     + -0.59358403E-03,-0.10009510E-01,-0.12290332E-02, 0.15322524E-01,
     + -0.62519987E-02, 0.22437006E-03,-0.21702710E-02, 0.49780305E-02,
     +  0.15468887E-02,-0.99281957E-02, 0.51485226E-02, 0.25714082E-02,
     + -0.49606833E-03,-0.41391253E-02, 0.24071780E-02,-0.33482458E-03,
     +  0.28274441E-03,-0.20879910E-02, 0.86177379E-03,-0.36448389E-03,
     + -0.94672729E-03,-0.33198434E-03, 0.81619452E-03,-0.52353443E-03,
     + -0.37354121E-02, 0.49330207E-03,-0.34950854E-03,
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
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x13 = x12*x1
      x14 = x13*x1
      x15 = x14*x1
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
c
c                  function
c
         theta_nq =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)*x11                
     3  +coeff(  3)    *x21            
     4  +coeff(  4)*x12                
     5  +coeff(  5)*x11*x21            
     6  +coeff(  6)    *x22            
     7  +coeff(  7)        *x32        
     8  +coeff(  8)        *x31*x41    
         theta_nq =   theta_nq 
     9  +coeff(  9)            *x42    
     1  +coeff( 10)                *x51
     2  +coeff( 11)*x13                
     3  +coeff( 12)*x12*x21            
     4  +coeff( 13)*x11*x22            
     5  +coeff( 14)    *x23            
     6  +coeff( 15)*x11    *x32        
     7  +coeff( 16)*x11    *x31*x41    
     8  +coeff( 17)    *x21*x31*x41    
         theta_nq =   theta_nq 
     9  +coeff( 18)*x11        *x42    
     1  +coeff( 19)    *x21    *x42    
     2  +coeff( 20)*x14                
     3  +coeff( 21)*x13*x21            
     4  +coeff( 22)*x12*x22            
     5  +coeff( 23)*x11*x23            
     6  +coeff( 24)*x12    *x32        
     7  +coeff( 25)*x11*x21*x32        
     8  +coeff( 26)    *x22*x32        
         theta_nq =   theta_nq 
     9  +coeff( 27)*x12    *x31*x41    
     1  +coeff( 28)*x11*x21*x31*x41    
     2  +coeff( 29)    *x22*x31*x41    
     3  +coeff( 30)*x11*x21    *x42    
     4  +coeff( 31)*x12*x21*x32        
     5  +coeff( 32)*x11*x22*x32        
     6  +coeff( 33)    *x23*x32        
     7  +coeff( 34)*x13    *x31*x41    
     8  +coeff( 35)*x12*x21*x31*x41    
         theta_nq =   theta_nq 
     9  +coeff( 36)*x11*x22    *x42    
     1  +coeff( 37)    *x23    *x42    
     2  +coeff( 38)*x13*x21    *x42    
     3  +coeff( 39)*x11*x24*x32        
     4  +coeff( 40)    *x21*x32        
     5  +coeff( 41)    *x24            
     6  +coeff( 42)*x12            *x51
     7  +coeff( 43)*x11*x24            
     8  +coeff( 44)    *x25            
         theta_nq =   theta_nq 
     9  +coeff( 45)*x13    *x32        
     1  +coeff( 46)    *x24*x32        
     2  +coeff( 47)*x11*x23*x31*x41    
     3  +coeff( 48)*x12*x23*x32        
     4  +coeff( 49)    *x25*x32        
     5  +coeff( 50)*x14*x21        *x51
     6  +coeff( 51)*x14    *x31*x43    
     7  +coeff( 52)*x13*x21*x31*x43    
     8  +coeff( 53)*x15*x22*x32        
         theta_nq =   theta_nq 
     9  +coeff( 54)    *x25    *x42*x51
     1  +coeff( 55)*x12*x21*x34    *x51
     2  +coeff( 56)*x12        *x42    
     3  +coeff( 57)*x15                
     4  +coeff( 58)*x14*x21            
     5  +coeff( 59)*x13*x22            
     6  +coeff( 60)*x12*x23            
     7  +coeff( 61)        *x32*x42    
     8  +coeff( 62)        *x31*x43    
         theta_nq =   theta_nq 
     9  +coeff( 63)            *x44    
     1  +coeff( 64)*x13        *x42    
     2  +coeff( 65)*x11    *x34        
     3  +coeff( 66)    *x21*x32    *x51
     4  +coeff( 67)*x12    *x34        
     5  +coeff( 68)*x12    *x33*x41    
     6  +coeff( 69)    *x22*x33*x41    
     7  +coeff( 70)*x12    *x32*x42    
     8  +coeff( 71)    *x22    *x42*x51
         theta_nq =   theta_nq 
     9  +coeff( 72)*x12*x23*x31*x41    
     1  +coeff( 73)*x13*x22    *x42    
     2  +coeff( 74)*x11*x22*x32    *x51
     3  +coeff( 75)    *x23*x32    *x51
     4  +coeff( 76)*x15*x21*x32        
     5  +coeff( 77)*x13*x23*x32        
     6  +coeff( 78)*x12*x24*x32        
     7  +coeff( 79)    *x24*x31*x41*x51
     8  +coeff( 80)*x15    *x34        
         theta_nq =   theta_nq 
     9  +coeff( 81)*x15*x21*x32    *x51
     1  +coeff( 82)*x11*x25*x31*x41*x51
     2  +coeff( 83)*x12*x24    *x42*x51
     3  +coeff( 84)*x11            *x51
     4  +coeff( 85)    *x21        *x51
     5  +coeff( 86)    *x22    *x42    
     6  +coeff( 87)*x11*x21        *x51
     7  +coeff( 88)    *x22        *x51
     8  +coeff( 89)        *x34        
         theta_nq =   theta_nq 
     9  +coeff( 90)        *x32    *x51
     1  +coeff( 91)        *x31*x41*x51
     2  +coeff( 92)            *x42*x51
     3  +coeff( 93)    *x23*x31*x41    
     4  +coeff( 94)*x12*x21        *x51
     5  +coeff( 95)*x11*x22        *x51
c
      return
      end
      function     y00_nq  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 96)
      data ncoeff/ 95/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.71654E+00,-0.30874E+00,-0.62118E+00,-0.39822E-01,-0.49978E-02,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.68856E+00, 0.16818E+00, 0.62118E+00, 0.39822E-01, 0.49966E-02,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.19317921E+00,-0.54024762E+00,-0.12692316E+01, 0.16116366E+01,
     + -0.46763337E+00, 0.27814679E+01,-0.53437054E+00, 0.26960301E+00,
     +  0.43400931E+00,-0.36668834E+00, 0.10296767E+00,-0.26381147E+01,
     + -0.30178246E+01, 0.20242734E+01, 0.73788002E-01,-0.35551705E+01,
     + -0.17727075E+01,-0.83380966E+01, 0.51152563E+01,-0.12477097E+02,
     +  0.26642940E+01,-0.59688205E+00, 0.11408877E+00, 0.93315041E+00,
     + -0.12892734E+01, 0.25117818E+00, 0.78723945E-01,-0.12169954E+01,
     +  0.21630566E+00,-0.37267822E+00, 0.18930998E+02,-0.30455062E+00,
     +  0.61485741E-01,-0.12736297E+01, 0.24236388E+01,-0.61053452E+01,
     +  0.41489344E+01,-0.27645178E+01, 0.33430018E+01,-0.19407346E+01,
     + -0.45730063E+00, 0.59655447E+01, 0.22293082E+00, 0.77249950E+00,
     +  0.70505452E+00, 0.62771797E+01, 0.36887851E+01, 0.68950486E+01,
     + -0.11438419E+01,-0.10242534E+01,-0.33269474E+00,-0.23051620E+01,
     +  0.29939144E+01, 0.18705090E+01, 0.93773171E-01,-0.14307291E+02,
     + -0.85102425E+01, 0.10035162E+02,-0.10203300E+02,-0.10408479E+01,
     + -0.61317915E+00, 0.60540628E+00, 0.54979282E-02,-0.59470093E+00,
     +  0.20110095E-01,-0.10071170E+00, 0.40835466E+01,-0.33869693E+00,
     +  0.20312366E+01, 0.20513511E+01,-0.14749131E+01, 0.20232497E-01,
     + -0.80499631E+00,-0.15104734E-01,-0.32375777E+00, 0.33878615E+00,
     +  0.14401078E+00, 0.25890374E+00, 0.33050352E+00,-0.59879196E+00,
     + -0.11212237E+01, 0.12771906E-02,-0.20216027E-01,-0.21623006E-01,
     + -0.20523916E-03, 0.11401569E-01, 0.70692796E+00, 0.26242018E+00,
     + -0.17313005E+01, 0.11422200E+01,-0.33655623E+00, 0.10503547E-01,
     +  0.11479689E+01, 0.19419877E-01,-0.27382545E-01,
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
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x13 = x12*x1
      x14 = x13*x1
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
c
c                  function
c
          y00_nq  =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)*x11    *x31        
     4  +coeff(  4)    *x21*x31        
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)*x11    *x33        
     7  +coeff(  7)*x11    *x32*x41    
     8  +coeff(  8)*x11        *x41    
          y00_nq  =    y00_nq  
     9  +coeff(  9)        *x33        
     1  +coeff( 10)        *x32*x41    
     2  +coeff( 11)            *x43    
     3  +coeff( 12)*x12    *x31        
     4  +coeff( 13)    *x21*x33        
     5  +coeff( 14)    *x21*x32*x41    
     6  +coeff( 15)        *x35        
     7  +coeff( 16)    *x22*x31        
     8  +coeff( 17)*x13    *x31        
          y00_nq  =    y00_nq  
     9  +coeff( 18)*x11*x22*x31        
     1  +coeff( 19)*x12    *x33        
     2  +coeff( 20)*x11*x21*x33        
     3  +coeff( 21)*x11*x21*x32*x41    
     4  +coeff( 22)*x12    *x31*x42    
     5  +coeff( 23)    *x22    *x43    
     6  +coeff( 24)*x13*x21*x31        
     7  +coeff( 25)    *x24*x31        
     8  +coeff( 26)*x14        *x41    
          y00_nq  =    y00_nq  
     9  +coeff( 27)*x13*x21    *x41    
     1  +coeff( 28)*x12*x22    *x41    
     2  +coeff( 29)*x11    *x35        
     3  +coeff( 30)    *x21*x35        
     4  +coeff( 31)*x11*x22*x33        
     5  +coeff( 32)*x13    *x32*x41    
     6  +coeff( 33)*x12*x21*x32*x41    
     7  +coeff( 34)*x12*x21*x31*x42    
     8  +coeff( 35)*x14    *x33        
          y00_nq  =    y00_nq  
     9  +coeff( 36)*x13*x21*x33        
     1  +coeff( 37)    *x24*x33        
     2  +coeff( 38)*x14    *x32*x41    
     3  +coeff( 39)*x13*x21*x32*x41    
     4  +coeff( 40)*x12*x22*x32*x41    
     5  +coeff( 41)    *x23*x35        
     6  +coeff( 42)*x11*x21*x31        
     7  +coeff( 43)*x12        *x41    
     8  +coeff( 44)    *x22    *x41    
          y00_nq  =    y00_nq  
     9  +coeff( 45)*x11        *x43    
     1  +coeff( 46)*x12*x21*x31        
     2  +coeff( 47)    *x23*x31        
     3  +coeff( 48)    *x22*x33        
     4  +coeff( 49)    *x22*x32*x41    
     5  +coeff( 50)    *x22*x31*x42    
     6  +coeff( 51)*x14    *x31        
     7  +coeff( 52)*x12*x22*x31        
     8  +coeff( 53)*x11*x23*x31        
          y00_nq  =    y00_nq  
     9  +coeff( 54)*x11*x23    *x41    
     1  +coeff( 55)*x11    *x34*x41    
     2  +coeff( 56)*x12*x21*x33        
     3  +coeff( 57)    *x23*x33        
     4  +coeff( 58)*x12*x22*x33        
     5  +coeff( 59)*x11*x23*x33        
     6  +coeff( 60)*x11*x21    *x41    
     7  +coeff( 61)    *x21*x31*x42    
     8  +coeff( 62)*x12    *x32*x41    
          y00_nq  =    y00_nq  
     9  +coeff( 63)        *x33    *x51
     1  +coeff( 64)    *x24    *x41    
     2  +coeff( 65)*x12        *x41*x51
     3  +coeff( 66)*x11    *x31*x44    
     4  +coeff( 67)*x13    *x33        
     5  +coeff( 68)*x11*x22*x31*x42    
     6  +coeff( 69)    *x23*x31*x42    
     7  +coeff( 70)*x12*x21    *x43    
     8  +coeff( 71)*x11*x22    *x43    
          y00_nq  =    y00_nq  
     9  +coeff( 72)*x14*x21*x31        
     1  +coeff( 73)    *x24*x31*x42    
     2  +coeff( 74)*x12*x24*x31        
     3  +coeff( 75)*x12*x24    *x41    
     4  +coeff( 76)*x11*x22*x35        
     5  +coeff( 77)*x14*x23*x31        
     6  +coeff( 78)*x12*x21*x35*x42    
     7  +coeff( 79)*x14*x24    *x43    
     8  +coeff( 80)*x14*x23*x35        
          y00_nq  =    y00_nq  
     9  +coeff( 81)*x11    *x31*x42    
     1  +coeff( 82)    *x21    *x41*x51
     2  +coeff( 83)        *x32*x41*x51
     3  +coeff( 84)*x11*x21    *x41*x51
     4  +coeff( 85)    *x21*x32*x41*x51
     5  +coeff( 86)    *x23    *x41*x51
     6  +coeff( 87)    *x22*x35        
     7  +coeff( 88)*x12    *x34*x41    
     8  +coeff( 89)*x11*x21*x34*x41    
          y00_nq  =    y00_nq  
     9  +coeff( 90)*x11*x21*x33*x42    
     1  +coeff( 91)*x12    *x32*x43    
     2  +coeff( 92)        *x33*x42*x51
     3  +coeff( 93)*x13*x21*x31*x42    
     4  +coeff( 94)*x11*x21*x33    *x51
     5  +coeff( 95)*x12    *x32*x41*x51
c
      return
      end
      function     phi_nq  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 96)
      data ncoeff/ 95/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.71654E+00,-0.30874E+00,-0.62118E+00,-0.39822E-01,-0.49978E-02,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.68856E+00, 0.16818E+00, 0.62118E+00, 0.39822E-01, 0.49966E-02,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.34925394E-01, 0.29465601E-01, 0.10617533E+00,-0.13527393E+00,
     +  0.45352455E-01,-0.28634382E-01,-0.25578803E+00,-0.37582450E-01,
     +  0.29617386E-01,-0.33964494E-02, 0.21928546E+00, 0.26685300E+00,
     + -0.18573001E+00,-0.14727386E-01,-0.33932105E-02,-0.27751001E-02,
     +  0.30020145E+00,-0.89651756E-01, 0.12645844E+00,-0.56508714E-02,
     +  0.14205028E+00, 0.70434803E+00, 0.61546043E-01,-0.41963091E+00,
     +  0.99008971E+00,-0.20851016E+00, 0.33217806E-01, 0.52758932E-01,
     + -0.64567156E-01, 0.10789455E+00,-0.18862395E-01,-0.35574932E-01,
     +  0.35488419E-01, 0.55790097E-02,-0.13111914E+01,-0.83576195E-01,
     + -0.21762791E+00, 0.63707963E-01,-0.92873998E-01, 0.65527134E-01,
     + -0.49697703E+00, 0.12408296E+00, 0.10050838E+00,-0.52130175E+00,
     + -0.31232673E+00, 0.11089904E+00,-0.14815505E+00,-0.57896399E+00,
     +  0.13635403E+00,-0.16914412E-01, 0.85864952E-02, 0.13626559E-01,
     +  0.19618250E+00,-0.25927055E+00, 0.50461110E-01,-0.18769780E+00,
     +  0.21387479E+00,-0.25806138E+00, 0.92069286E+00, 0.63635945E+00,
     +  0.58124151E-01, 0.11951684E-01,-0.14117996E+00,-0.29021913E+00,
     +  0.52504903E+00,-0.41753784E-01, 0.12695915E+00,-0.16237740E+00,
     +  0.10865751E+00, 0.10302030E+00,-0.45842361E-02,-0.37683826E-03,
     +  0.43455344E-01, 0.12779565E+00,-0.31889670E-01,-0.14129721E+00,
     +  0.74026398E-01,-0.54250132E-01,-0.46150096E-01,-0.27902586E-01,
     + -0.31360364E-03, 0.10069006E-02, 0.47111772E-02,-0.40989104E-03,
     + -0.41511604E+00, 0.22823921E+00, 0.12734393E-01,-0.68288662E-01,
     + -0.58834709E-03,-0.75453827E-02, 0.71083372E-02, 0.40267952E-03,
     +  0.11051905E-01,-0.28844725E-01, 0.10746124E+00,
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
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x13 = x12*x1
      x14 = x13*x1
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
      x45 = x44*x4
      x51 = x5
c
c                  function
c
          phi_nq  =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)*x11    *x31        
     4  +coeff(  4)    *x21*x31        
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)*x11        *x41    
     7  +coeff(  7)*x11    *x33        
     8  +coeff(  8)        *x33        
          phi_nq  =    phi_nq  
     9  +coeff(  9)        *x32*x41    
     1  +coeff( 10)        *x31*x42    
     2  +coeff( 11)*x12    *x31        
     3  +coeff( 12)    *x21*x33        
     4  +coeff( 13)    *x21*x32*x41    
     5  +coeff( 14)*x11    *x31*x42    
     6  +coeff( 15)        *x35        
     7  +coeff( 16)            *x43    
     8  +coeff( 17)    *x22*x31        
          phi_nq  =    phi_nq  
     9  +coeff( 18)    *x22    *x41    
     1  +coeff( 19)*x11    *x32*x41    
     2  +coeff( 20)*x11        *x43    
     3  +coeff( 21)*x13    *x31        
     4  +coeff( 22)*x11*x22*x31        
     5  +coeff( 23)    *x23    *x41    
     6  +coeff( 24)*x12    *x33        
     7  +coeff( 25)*x11*x21*x33        
     8  +coeff( 26)*x11*x21*x32*x41    
          phi_nq  =    phi_nq  
     9  +coeff( 27)*x12    *x31*x42    
     1  +coeff( 28)    *x22*x31*x42    
     2  +coeff( 29)*x13*x21*x31        
     3  +coeff( 30)    *x24*x31        
     4  +coeff( 31)*x14        *x41    
     5  +coeff( 32)*x12*x22    *x41    
     6  +coeff( 33)*x11    *x35        
     7  +coeff( 34)    *x21*x35        
     8  +coeff( 35)*x11*x22*x33        
          phi_nq  =    phi_nq  
     9  +coeff( 36)*x14    *x33        
     1  +coeff( 37)    *x24*x33        
     2  +coeff( 38)*x14    *x32*x41    
     3  +coeff( 39)*x12*x22*x32*x41    
     4  +coeff( 40)*x13    *x35        
     5  +coeff( 41)*x11*x21*x31        
     6  +coeff( 42)*x11*x21    *x41    
     7  +coeff( 43)    *x21*x31*x42    
     8  +coeff( 44)*x12*x21*x31        
          phi_nq  =    phi_nq  
     9  +coeff( 45)    *x23*x31        
     1  +coeff( 46)*x12*x21    *x41    
     2  +coeff( 47)*x11*x22    *x41    
     3  +coeff( 48)    *x22*x33        
     4  +coeff( 49)    *x22*x32*x41    
     5  +coeff( 50)*x12        *x43    
     6  +coeff( 51)*x11*x21    *x43    
     7  +coeff( 52)*x14    *x31        
     8  +coeff( 53)*x12*x22*x31        
          phi_nq  =    phi_nq  
     9  +coeff( 54)*x11*x23*x31        
     1  +coeff( 55)*x13*x21    *x41    
     2  +coeff( 56)*x11    *x34*x41    
     3  +coeff( 57)*x11    *x33*x42    
     4  +coeff( 58)*x13    *x33        
     5  +coeff( 59)*x12*x21*x33        
     6  +coeff( 60)    *x23*x33        
     7  +coeff( 61)*x13    *x32*x41    
     8  +coeff( 62)*x12*x21*x32*x41    
          phi_nq  =    phi_nq  
     9  +coeff( 63)*x13    *x31*x42    
     1  +coeff( 64)*x12*x22*x33        
     2  +coeff( 65)*x11*x23*x33        
     3  +coeff( 66)*x12        *x41    
     4  +coeff( 67)*x11    *x32*x43    
     5  +coeff( 68)    *x23*x32*x41    
     6  +coeff( 69)*x12*x21*x31*x42    
     7  +coeff( 70)*x11*x22*x31*x42    
     8  +coeff( 71)*x14*x21*x31        
          phi_nq  =    phi_nq  
     9  +coeff( 72)        *x35    *x51
     1  +coeff( 73)*x11*x24*x32*x41    
     2  +coeff( 74)*x14    *x35        
     3  +coeff( 75)*x13*x21*x35        
     4  +coeff( 76)*x12*x22*x35        
     5  +coeff( 77)*x11*x23*x34*x41    
     6  +coeff( 78)*x14    *x33*x42    
     7  +coeff( 79)    *x21    *x43    
     8  +coeff( 80)*x13        *x41    
          phi_nq  =    phi_nq  
     9  +coeff( 81)    *x21    *x41*x51
     1  +coeff( 82)        *x32*x41*x51
     2  +coeff( 83)    *x24    *x41    
     3  +coeff( 84)*x12    *x31    *x51
     4  +coeff( 85)*x11    *x31*x44    
     5  +coeff( 86)*x11        *x45    
     6  +coeff( 87)*x13        *x43    
     7  +coeff( 88)*x12*x21    *x43    
     8  +coeff( 89)    *x21    *x43*x51
          phi_nq  =    phi_nq  
     9  +coeff( 90)*x12*x23*x31        
     1  +coeff( 91)*x12*x23    *x41    
     2  +coeff( 92)*x11*x22*x31    *x51
     3  +coeff( 93)*x12    *x35        
     4  +coeff( 94)*x11*x21    *x45    
     5  +coeff( 95)*x13*x21*x33        
c
      return
      end
