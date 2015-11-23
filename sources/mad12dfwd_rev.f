c  MAD spectrometer forward and reverse transfer functions for MAD in its 12 degree
c  see line 5402 and following for instructions on the reverse functions
c  minimum angle configuration. Formulated 12/07/04 -JJL 
c  
c forward function general description: 
c functions take trajectories from target to various planes in the spectrometer. 
c The planes are numbered 0 through 11 (0,2,3,4,5,7,9,11) 
c plane 0: the target plane                                  
c plane 2: entrance to magnet 1                               
c plane 3: middle of magnet 1                                 
c plane 4: exit of magnet 1                                   
c plane 5: entrance to magnet 2                               
c plane 7: middle of magnet 2                                
c plane 9: exit of magnet 2                         
c plane 11: 1st drift chamber                                
c 
c useful endplane information 
c all functions give the locations and directions of a trjectory relative to the central one 
c x = x location in the local (magnet coordinate system) of the central trajectory (+x is down) 
c cx = x-direction cosine of the central trajectory (approx. theta) 
c cz = z-direction cosine of the central trajectory (~1)
c cy = y-direction cosine of the central trajectory (approx. phi)(zero)  
c l = central trajectory path length in mm 
c ep	         x	   cx		cz		cy		l 
c   2	    -0.666043	0.085701	0.996321	0		6999.552246   
c   3	   104.630898  -0.010819	0.999942	0.000001	9003.19043   
c   4	   -73.731644  -0.105757	0.994393	0.000002	11311.05469   
c   5	    11.757298	0.148056	0.988979	0.000002	11616.64648   
c   7	   156.909363	0.001901	0.999999	0		13324.50391   
c   9	   -26.31488   -0.150983	0.988537	-0.000002	15335.125   
c  11	  -151.757401	0.020135	0.999797	-0.000003	16274.31348   

 
c 
c for various planes there are 5 functions: 
c  x_m12_0_j(x,m): gives x position at plane j as a function of trajectory 
c              parameters at plane 0. 
c  t_m12_0_j(x,m): gives theta at plane j 
c  y_m12_0_j(x,m): gives y at plane j 
c  p_m12_0_j(x,m): gives phi at plane j 
c  l_m12_0_j(x,m): gives difference in length from the "central trajectory 
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
c      x2(1) = x_m12_0_2(x0,5) 
c      x2(2) = t_m12_0_2(x0,5) 
c      x2(3) = y_m12_0_2(x0,5) 
c      x2(4) = p_m12_0_2(x0,5) 
c      x2(5) = x0(5) 
c  at plane 3
c      x3(1) = x_m12_0_3(x0,5)  etc...  
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
c  l_m12_0_j(x,5) gives the difference between the ray in question and 
c             the central ray going from plane 0 to plane j 
c 
      function x_m12_0_2   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)   !,xmean(10)
      dimension coeff( 13)
      data ncoeff/ 12/
      data avdat/  0.8178358E-03/
      data xmin/
     1 -0.49998E-02,-0.79804E-01,-0.99991E-01,-0.37074E-01,-0.42491E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49964E-02, 0.81140E-01, 0.99991E-01, 0.37074E-01, 0.42491E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.38782253E-02, 0.56532139E+00, 0.41188491E-02, 0.50164228E-02,
     +  0.13115193E-03, 0.74945921E-04, 0.97564960E-04, 0.44290588E-04,
     + -0.59440044E-04, 0.11734436E-03, 0.61505838E-04, 0.37713449E-04,
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
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      x_m12_0_2   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x22            
     4  +coeff(  4)*x11                
     5  +coeff(  5)                *x51
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21        *x51
      x_m12_0_2   =x_m12_0_2   
     9  +coeff(  9)                *x52
     1  +coeff( 10)    *x23            
     2  +coeff( 11)    *x21    *x42    
     3  +coeff( 12)*x11*x21            
c
      return
      end
      function t_m12_0_2   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)   !,xmean(10)
      dimension coeff( 53)
      data ncoeff/ 52/
      data avdat/  0.2053692E-03/
      data xmin/
     1 -0.49998E-02,-0.79804E-01,-0.99991E-01,-0.37074E-01,-0.42491E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49964E-02, 0.81140E-01, 0.99991E-01, 0.37074E-01, 0.42491E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.47260657E-03, 0.79461306E-01, 0.61160483E-03, 0.64902869E-03,
     +  0.30211059E-03, 0.39966949E-03, 0.42847535E-03,-0.31148322E-03,
     +  0.65587374E-03, 0.34732220E-03, 0.49127761E-03,-0.26625095E-03,
     + -0.11628612E-03,-0.16920306E-03, 0.54646349E-04,-0.20592879E-03,
     +  0.14115633E-03,-0.27780130E-03,-0.15862711E-03,-0.20985543E-03,
     +  0.71317321E-04, 0.13490376E-03, 0.64866435E-04, 0.82271035E-04,
     + -0.17397149E-03,-0.26579448E-04, 0.11419757E-04,-0.43503278E-04,
     + -0.70611197E-04, 0.93086055E-04,-0.11687205E-03, 0.12761570E-03,
     +  0.65631422E-04, 0.87286928E-04,-0.76149340E-05,-0.26891854E-04,
     +  0.15130616E-04,-0.66912624E-04,-0.36924663E-04,-0.37645172E-04,
     +  0.76433826E-04,-0.91068641E-05,-0.10228237E-04, 0.12043294E-04,
     + -0.19378822E-04, 0.26286645E-04,-0.42718816E-05, 0.29554780E-04,
     +  0.55808596E-04,-0.58142483E-04,-0.31156102E-04,-0.40461637E-04,
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
      t_m12_0_2   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22            
     5  +coeff(  5)        *x31*x41    
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)                *x52
      t_m12_0_2   =t_m12_0_2   
     9  +coeff(  9)    *x23            
     1  +coeff( 10)    *x21*x31*x41    
     2  +coeff( 11)    *x21    *x42    
     3  +coeff( 12)    *x22        *x51
     4  +coeff( 13)        *x31*x41*x51
     5  +coeff( 14)            *x42*x51
     6  +coeff( 15)        *x32        
     7  +coeff( 16)    *x21        *x52
     8  +coeff( 17)                *x53
      t_m12_0_2   =t_m12_0_2   
     9  +coeff( 18)    *x23        *x51
     1  +coeff( 19)    *x21*x31*x41*x51
     2  +coeff( 20)    *x21    *x42*x51
     3  +coeff( 21)    *x21*x32        
     4  +coeff( 22)    *x22        *x52
     5  +coeff( 23)        *x31*x41*x52
     6  +coeff( 24)            *x42*x52
     7  +coeff( 25)    *x23    *x42    
     8  +coeff( 26)        *x32    *x51
      t_m12_0_2   =t_m12_0_2   
     9  +coeff( 27)*x11*x21            
     1  +coeff( 28)    *x22*x31*x41    
     2  +coeff( 29)    *x22    *x42    
     3  +coeff( 30)    *x21        *x53
     4  +coeff( 31)    *x23*x31*x41    
     5  +coeff( 32)    *x23        *x52
     6  +coeff( 33)    *x21*x31*x41*x52
     7  +coeff( 34)    *x21    *x42*x52
     8  +coeff( 35)*x11                
      t_m12_0_2   =t_m12_0_2   
     9  +coeff( 36)    *x21*x32    *x51
     1  +coeff( 37)*x11*x22            
     2  +coeff( 38)    *x22        *x53
     3  +coeff( 39)        *x31*x41*x53
     4  +coeff( 40)            *x42*x53
     5  +coeff( 41)    *x23    *x42*x51
     6  +coeff( 42)    *x22*x32        
     7  +coeff( 43)        *x31*x43    
     8  +coeff( 44)        *x32    *x52
      t_m12_0_2   =t_m12_0_2   
     9  +coeff( 45)    *x23*x32        
     1  +coeff( 46)    *x21*x31*x43    
     2  +coeff( 47)*x11*x21        *x51
     3  +coeff( 48)    *x22    *x42*x51
     4  +coeff( 49)    *x23*x31*x41*x51
     5  +coeff( 50)    *x23        *x53
     6  +coeff( 51)    *x21*x31*x41*x53
     7  +coeff( 52)    *x21    *x42*x53
c
      return
      end
      function y_m12_0_2   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)   !,xmean(10)
      dimension coeff(  6)
      data ncoeff/  5/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49998E-02,-0.79804E-01,-0.99991E-01,-0.37074E-01,-0.42491E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49964E-02, 0.81140E-01, 0.99991E-01, 0.37074E-01, 0.42491E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.99986814E-01, 0.25955087E+00, 0.16585028E-02,-0.67097237E-04,
     + -0.10274105E-03,
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
      x22 = x21*x2
      x31 = x3
      x41 = x4
c
c                  function
c
      y_m12_0_2   =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)    *x21*x31        
     5  +coeff(  5)    *x22    *x41    
c
      return
      end
      function p_m12_0_2   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)   !,xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49998E-02,-0.79804E-01,-0.99991E-01,-0.37074E-01,-0.42491E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49964E-02, 0.81140E-01, 0.99991E-01, 0.37074E-01, 0.42491E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.12171332E-03, 0.37523359E-01,-0.36336467E-03,-0.66587894E-03,
     + -0.32186674E-03,-0.82327862E-03, 0.16994818E-03, 0.37981087E-03,
     +  0.29901622E-03, 0.34832137E-03, 0.15504191E-04, 0.29373607E-04,
     + -0.20096160E-03,-0.13603065E-03,-0.53495311E-04,-0.42940963E-04,
     +  0.11085768E-03, 0.11512142E-03, 0.10606045E-05,-0.85995285E-04,
     + -0.11774975E-03,-0.51326995E-04,-0.24184441E-04, 0.26582244E-04,
     + -0.18197570E-03, 0.99332635E-06,-0.14479983E-04, 0.94263625E-04,
     +  0.73421470E-04, 0.12511157E-04, 0.24290437E-04, 0.23056880E-04,
     + -0.23167104E-04,-0.69563270E-04,-0.76509205E-04,-0.33287970E-05,
     +  0.20958489E-04, 0.74231677E-04,-0.27602841E-05, 0.31280462E-04,
     + -0.41914882E-05,-0.83029281E-05, 0.23228333E-04,-0.24367944E-04,
     +  0.12589803E-04, 0.54324526E-04, 0.47945399E-04, 0.88588582E-04,
     + -0.23268155E-05,-0.11138703E-04, 0.87672915E-05, 0.39703722E-04,
     + -0.64077614E-06, 0.94941715E-05, 0.52971754E-04,-0.10318915E-04,
     +  0.36929901E-04,-0.41228726E-04,-0.18697353E-03,-0.44929315E-04,
     +  0.76195924E-04,-0.10702814E-04,-0.24447327E-04,-0.94390607E-05,
     +  0.46458567E-05,-0.11522658E-04, 0.12339031E-04, 0.27250095E-04,
     +  0.13116833E-04, 0.73606966E-05, 0.10868825E-04, 0.49690613E-06,
     + -0.12393856E-04, 0.35354307E-04,-0.58444070E-05,-0.23524231E-04,
     + -0.12759690E-04,-0.19663439E-04,-0.32771877E-04,-0.13248093E-04,
     +  0.18190141E-04,-0.21707287E-04, 0.13780039E-04,-0.23799721E-04,
     + -0.13525508E-04,-0.35801037E-04, 0.54758948E-04,-0.27412170E-04,
     +  0.20309591E-04, 0.50640534E-04,
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
      p_m12_0_2   =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21*x31        
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)    *x22*x31        
     6  +coeff(  6)    *x22    *x41    
     7  +coeff(  7)    *x21*x31    *x51
     8  +coeff(  8)    *x21    *x41*x51
      p_m12_0_2   =p_m12_0_2   
     9  +coeff(  9)    *x23    *x41    
     1  +coeff( 10)    *x22    *x41*x51
     2  +coeff( 11)        *x32*x41*x51
     3  +coeff( 12)            *x43*x51
     4  +coeff( 13)    *x21    *x41*x52
     5  +coeff( 14)            *x41*x51
     6  +coeff( 15)        *x31*x42    
     7  +coeff( 16)            *x43    
     8  +coeff( 17)    *x23*x31        
      p_m12_0_2   =p_m12_0_2   
     9  +coeff( 18)    *x22*x31    *x51
     1  +coeff( 19)        *x33    *x51
     2  +coeff( 20)    *x21*x31    *x52
     3  +coeff( 21)    *x23    *x41*x51
     4  +coeff( 22)        *x31    *x51
     5  +coeff( 23)        *x32*x41    
     6  +coeff( 24)    *x21*x31*x42    
     7  +coeff( 25)    *x22    *x41*x52
     8  +coeff( 26)        *x32*x41*x52
      p_m12_0_2   =p_m12_0_2   
     9  +coeff( 27)            *x43*x52
     1  +coeff( 28)    *x21    *x41*x53
     2  +coeff( 29)            *x41*x52
     3  +coeff( 30)    *x21*x32*x41    
     4  +coeff( 31)    *x21    *x43    
     5  +coeff( 32)        *x31*x42*x51
     6  +coeff( 33)*x11*x21    *x41    
     7  +coeff( 34)    *x23*x31    *x51
     8  +coeff( 35)    *x22*x31    *x52
      p_m12_0_2   =p_m12_0_2   
     9  +coeff( 36)        *x33    *x52
     1  +coeff( 37)    *x21*x31    *x53
     2  +coeff( 38)    *x23    *x41*x52
     3  +coeff( 39)        *x33        
     4  +coeff( 40)        *x31    *x52
     5  +coeff( 41)*x11    *x31        
     6  +coeff( 42)*x11        *x41    
     7  +coeff( 43)    *x22*x31*x42    
     8  +coeff( 44)        *x31*x42*x52
      p_m12_0_2   =p_m12_0_2   
     9  +coeff( 45)*x11*x22    *x41    
     1  +coeff( 46)    *x23*x31    *x52
     2  +coeff( 47)    *x22*x31    *x53
     3  +coeff( 48)    *x22    *x41*x53
     4  +coeff( 49)        *x32*x41*x53
     5  +coeff( 50)            *x43*x53
     6  +coeff( 51)*x11*x23    *x41    
     7  +coeff( 52)    *x22*x32*x43    
     8  +coeff( 53)*x11*x21*x32*x41*x51
      p_m12_0_2   =p_m12_0_2   
     9  +coeff( 54)*x11*x21*x32*x43    
     1  +coeff( 55)*x11*x23*x31    *x52
     2  +coeff( 56)*x11*x21*x33    *x52
     3  +coeff( 57)*x11*x23    *x41*x52
     4  +coeff( 58)*x11*x23    *x41*x53
     5  +coeff( 59)*x11*x21*x32*x43*x52
     6  +coeff( 60)*x11    *x32*x43*x53
     7  +coeff( 61)*x11*x21*x32*x43*x53
     8  +coeff( 62)        *x31    *x53
      p_m12_0_2   =p_m12_0_2   
     9  +coeff( 63)            *x41*x53
     1  +coeff( 64)*x11*x21*x31        
     2  +coeff( 65)*x11        *x41*x51
     3  +coeff( 66)    *x21*x31*x42*x51
     4  +coeff( 67)*x11*x22*x31        
     5  +coeff( 68)    *x21*x33*x42    
     6  +coeff( 69)    *x22*x33    *x51
     7  +coeff( 70)*x11*x21    *x41*x51
     8  +coeff( 71)    *x22*x31*x42*x51
      p_m12_0_2   =p_m12_0_2   
     9  +coeff( 72)*x11    *x31    *x52
     1  +coeff( 73)    *x21*x33    *x52
     2  +coeff( 74)*x11*x21*x31*x42    
     3  +coeff( 75)    *x23*x31*x42*x51
     4  +coeff( 76)    *x21*x33*x42*x51
     5  +coeff( 77)    *x21*x32*x43*x51
     6  +coeff( 78)        *x32*x43*x52
     7  +coeff( 79)    *x23    *x41*x53
     8  +coeff( 80)*x11*x22*x31*x42    
      p_m12_0_2   =p_m12_0_2   
     9  +coeff( 81)    *x23*x32*x43    
     1  +coeff( 82)*x11*x23*x31    *x51
     2  +coeff( 83)*x11*x21*x33    *x51
     3  +coeff( 84)    *x22*x32*x43*x51
     4  +coeff( 85)*x11*x22*x31    *x52
     5  +coeff( 86)*x11*x21*x33*x42    
     6  +coeff( 87)*x11*x21*x32*x41*x52
     7  +coeff( 88)*x11*x21*x31*x42*x52
     8  +coeff( 89)*x11*x21    *x43*x52
      p_m12_0_2   =p_m12_0_2   
     9  +coeff( 90)    *x23*x31*x42*x53
c
      return
      end
      real function l_m12_0_2   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)   !,xmean(10)
      dimension coeff( 14)
      data ncoeff/ 13/
      data avdat/ -0.8915826E-02/
      data xmin/
     1 -0.49998E-02,-0.79804E-01,-0.99991E-01,-0.37074E-01,-0.42491E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49964E-02, 0.81140E-01, 0.99991E-01, 0.37074E-01, 0.42491E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.85028987E-02,-0.49617205E-01,-0.23008579E-01,-0.48145135E-02,
     + -0.43779719E-03,-0.18072620E-03,-0.12609211E-04,-0.40670697E-04,
     + -0.78836892E-05,-0.13832213E-04, 0.57028151E-05,-0.80536229E-05,
     +  0.28135304E-04,
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
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_m12_0_2   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x22            
     4  +coeff(  4)            *x42    
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x23            
     7  +coeff(  7)                *x51
     8  +coeff(  8)    *x21    *x42    
      l_m12_0_2   =l_m12_0_2   
     9  +coeff(  9)        *x31*x41    
     1  +coeff( 10)    *x21        *x51
     2  +coeff( 11)                *x52
     3  +coeff( 12)    *x21*x31*x41    
     4  +coeff( 13)    *x24            
c
      return
      end
      function x_m12_0_3   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)   !,xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.9422942E-01/
      data xmin/
     1 -0.49998E-02,-0.79804E-01,-0.99991E-01,-0.37038E-01,-0.42491E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49981E-02, 0.64043E-01, 0.99991E-01, 0.37038E-01, 0.42498E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.30563785E-01, 0.58000904E+00, 0.25823209E-01, 0.94654309E-02,
     +  0.23002439E-03, 0.12279121E-02, 0.19923707E-02, 0.28580286E-01,
     + -0.12631257E-01, 0.44209510E-02,-0.13895431E-01, 0.54479465E-02,
     + -0.54809493E-02,-0.64985552E-02, 0.60125529E-02,-0.50180196E-03,
     +  0.51634535E-04,-0.27363866E-02, 0.29284388E-03, 0.46733175E-02,
     +  0.11929334E-02, 0.14951058E-02,-0.83546696E-03,-0.92767784E-03,
     + -0.32334050E-03, 0.10207526E-02,-0.27984125E-02,-0.11141909E-02,
     +  0.25981998E-02, 0.30220104E-02, 0.78137935E-03, 0.12503630E-03,
     + -0.10672101E-02, 0.29218910E-03,-0.33316417E-02, 0.31726824E-02,
     + -0.22528727E-03, 0.22339420E-02, 0.52245092E-02,-0.34782002E-04,
     + -0.25500190E-02, 0.47330213E-04, 0.18058573E-03,-0.41733615E-03,
     + -0.81024948E-03,-0.73293777E-03,-0.54057874E-03, 0.53165207E-03,
     + -0.15742781E-03,-0.73574614E-04, 0.74070337E-03, 0.16984331E-02,
     + -0.85163105E-03,-0.47207915E-03, 0.12526262E-02,-0.38312361E-03,
     + -0.14789986E-03, 0.31103095E-03,-0.22106725E-02,-0.54423162E-03,
     +  0.95785473E-03,-0.96803717E-03,-0.40414848E-03,-0.12291590E-03,
     +  0.24076428E-03,-0.90696973E-04, 0.14750054E-03,-0.33425727E-04,
     +  0.11228614E-03,-0.18288688E-03,-0.10104124E-03,-0.72087400E-03,
     + -0.24941866E-03, 0.91736074E-04,-0.33765883E-03,-0.10854466E-02,
     +  0.43712993E-03, 0.17761385E-03, 0.37249242E-03, 0.47403251E-03,
     +  0.22871202E-03, 0.49783121E-03, 0.15455409E-03, 0.83431609E-04,
     + -0.40133356E-03,-0.11651024E-03, 0.16617100E-03, 0.19067027E-03,
     + -0.19273392E-03, 0.79953141E-03,
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
      x_m12_0_3   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22            
     5  +coeff(  5)        *x32        
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21        *x51
      x_m12_0_3   =x_m12_0_3   
     9  +coeff(  9)                *x52
     1  +coeff( 10)*x11                
     2  +coeff( 11)    *x21        *x52
     3  +coeff( 12)                *x53
     4  +coeff( 13)    *x22*x31*x41    
     5  +coeff( 14)    *x22    *x42    
     6  +coeff( 15)    *x21        *x53
     7  +coeff( 16)    *x22*x32    *x51
     8  +coeff( 17)    *x22        *x53
      x_m12_0_3   =x_m12_0_3   
     9  +coeff( 18)    *x22        *x51
     1  +coeff( 19)*x11            *x51
     2  +coeff( 20)    *x22*x31*x43    
     3  +coeff( 21)    *x21*x31*x41    
     4  +coeff( 22)    *x21    *x42    
     5  +coeff( 23)        *x31*x41*x51
     6  +coeff( 24)            *x42*x51
     7  +coeff( 25)        *x32*x42    
     8  +coeff( 26)    *x22        *x52
      x_m12_0_3   =x_m12_0_3   
     9  +coeff( 27)    *x23    *x42    
     1  +coeff( 28)    *x21*x31*x43    
     2  +coeff( 29)    *x22*x31*x41*x51
     3  +coeff( 30)    *x22    *x42*x51
     4  +coeff( 31)    *x23            
     5  +coeff( 32)*x11*x21            
     6  +coeff( 33)    *x22*x32        
     7  +coeff( 34)            *x42*x52
     8  +coeff( 35)    *x23*x31*x41    
      x_m12_0_3   =x_m12_0_3   
     9  +coeff( 36)    *x22*x32*x42    
     1  +coeff( 37)    *x21*x31*x43*x51
     2  +coeff( 38)    *x23*x32*x42    
     3  +coeff( 39)    *x23*x31*x43    
     4  +coeff( 40)    *x22*x32*x42*x51
     5  +coeff( 41)    *x22*x31*x43*x51
     6  +coeff( 42)    *x21*x32        
     7  +coeff( 43)        *x33*x41    
     8  +coeff( 44)        *x31*x43    
      x_m12_0_3   =x_m12_0_3   
     9  +coeff( 45)    *x23        *x51
     1  +coeff( 46)    *x21    *x42*x51
     2  +coeff( 47)    *x23*x32        
     3  +coeff( 48)        *x31*x43*x51
     4  +coeff( 49)*x11            *x52
     5  +coeff( 50)            *x42*x53
     6  +coeff( 51)    *x23*x31*x41*x51
     7  +coeff( 52)    *x23    *x42*x51
     8  +coeff( 53)    *x22    *x42*x52
      x_m12_0_3   =x_m12_0_3   
     9  +coeff( 54)    *x21    *x42*x53
     1  +coeff( 55)    *x22*x32    *x53
     2  +coeff( 56)*x11*x21*x31*x43    
     3  +coeff( 57)*x11    *x32*x42*x51
     4  +coeff( 58)*x11*x22*x32    *x52
     5  +coeff( 59)    *x22*x32*x42*x53
     6  +coeff( 60)        *x33*x43*x53
     7  +coeff( 61)*x11    *x32*x42*x53
     8  +coeff( 62)*x11*x21*x31*x43*x53
      x_m12_0_3   =x_m12_0_3   
     9  +coeff( 63)    *x21*x31*x41*x51
     1  +coeff( 64)        *x32    *x52
     2  +coeff( 65)        *x31*x41*x52
     3  +coeff( 66)*x11*x22            
     4  +coeff( 67)    *x21*x33*x41    
     5  +coeff( 68)*x11*x21        *x51
     6  +coeff( 69)        *x32*x42*x51
     7  +coeff( 70)*x11*x21    *x42    
     8  +coeff( 71)*x11    *x32    *x51
      x_m12_0_3   =x_m12_0_3   
     9  +coeff( 72)    *x23*x32    *x51
     1  +coeff( 73)    *x21*x32*x42*x51
     2  +coeff( 74)*x11*x21        *x52
     3  +coeff( 75)    *x22*x32    *x52
     4  +coeff( 76)    *x22*x31*x41*x52
     5  +coeff( 77)        *x32*x42*x52
     6  +coeff( 78)        *x31*x43*x52
     7  +coeff( 79)    *x23        *x53
     8  +coeff( 80)    *x21*x32    *x53
      x_m12_0_3   =x_m12_0_3   
     9  +coeff( 81)*x11*x21*x31*x41*x51
     1  +coeff( 82)    *x22*x33*x41*x51
     2  +coeff( 83)*x11*x22        *x52
     3  +coeff( 84)    *x23*x31*x41*x52
     4  +coeff( 85)    *x21*x31*x43*x52
     5  +coeff( 86)*x11*x23*x32        
     6  +coeff( 87)*x11*x21*x33*x41    
     7  +coeff( 88)*x11*x23    *x42    
     8  +coeff( 89)*x11*x22    *x42*x51
      x_m12_0_3   =x_m12_0_3   
     9  +coeff( 90)    *x23*x32*x42*x51
c
      return
      end
      function t_m12_0_3   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)   !,xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.5044550E-02/
      data xmin/
     1 -0.49998E-02,-0.79804E-01,-0.99991E-01,-0.37038E-01,-0.42491E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49981E-02, 0.64043E-01, 0.99991E-01, 0.37038E-01, 0.42498E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.65747476E-02,-0.13865067E-01, 0.36103491E-01, 0.19576524E-02,
     +  0.89356974E-04, 0.32928865E-03, 0.45340133E-03, 0.34892634E-01,
     + -0.17350912E-01,-0.75714348E-03,-0.54667122E-03,-0.21279848E-03,
     + -0.64321508E-03,-0.30791896E-03,-0.10509423E-02, 0.86444215E-05,
     +  0.29649906E-03, 0.15400404E-03,-0.16575634E-01, 0.74653486E-02,
     + -0.82524968E-02,-0.86824344E-02,-0.57534839E-03,-0.75346534E-03,
     +  0.20694263E-03,-0.70694520E-03,-0.22682834E-03,-0.49398473E-03,
     + -0.31781840E-03,-0.19465620E-03,-0.27250306E-03,-0.18328891E-03,
     + -0.34654423E-03, 0.68691135E-02,-0.68690494E-03,-0.19901434E-02,
     +  0.11363842E-02, 0.17707573E-02,-0.14957905E-03, 0.29836339E-02,
     + -0.19275928E-03,-0.43358843E-03,-0.34675072E-03,-0.42568831E-03,
     +  0.88052603E-03,-0.19245490E-03, 0.51549841E-02, 0.81238849E-02,
     + -0.45451252E-04, 0.20702062E-02,-0.26403653E-03,-0.11780348E-02,
     +  0.98574557E-03,-0.14391661E-02, 0.13544444E-03, 0.10098125E-02,
     + -0.82998558E-04, 0.25356060E-03,-0.11342084E-04,-0.15040780E-03,
     +  0.64801273E-03,-0.47258087E-03, 0.76972827E-03, 0.95129602E-04,
     +  0.79743179E-04,-0.45017365E-03, 0.59535418E-03, 0.16781401E-03,
     +  0.69292367E-03,-0.44099175E-03,-0.15656310E-03,-0.11891742E-03,
     + -0.21986171E-03, 0.30558393E-03, 0.32317729E-03, 0.29657243E-03,
     +  0.26977903E-03,-0.22700406E-04,-0.17916964E-04,-0.90383161E-04,
     + -0.27250373E-04, 0.27580843E-04, 0.88252375E-04, 0.15842932E-04,
     +  0.11051339E-04,-0.98249075E-05,-0.32237913E-04,-0.29489624E-04,
     +  0.14423899E-04, 0.72998992E-04,
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
      t_m12_0_3   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22            
     5  +coeff(  5)        *x32        
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21        *x51
      t_m12_0_3   =t_m12_0_3   
     9  +coeff(  9)                *x52
     1  +coeff( 10)*x11                
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x21*x32        
     4  +coeff( 13)    *x21*x31*x41    
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)    *x22        *x51
     7  +coeff( 16)        *x32    *x51
     8  +coeff( 17)        *x31*x41*x51
      t_m12_0_3   =t_m12_0_3   
     9  +coeff( 18)            *x42*x51
     1  +coeff( 19)    *x21        *x52
     2  +coeff( 20)                *x53
     3  +coeff( 21)    *x22*x31*x41    
     4  +coeff( 22)    *x22    *x42    
     5  +coeff( 23)        *x32*x42    
     6  +coeff( 24)        *x31*x43    
     7  +coeff( 25)*x11            *x51
     8  +coeff( 26)    *x23        *x51
      t_m12_0_3   =t_m12_0_3   
     9  +coeff( 27)    *x21*x32    *x51
     1  +coeff( 28)    *x21*x31*x41*x51
     2  +coeff( 29)    *x21    *x42*x51
     3  +coeff( 30)    *x22        *x52
     4  +coeff( 31)        *x32    *x52
     5  +coeff( 32)        *x31*x41*x52
     6  +coeff( 33)            *x42*x52
     7  +coeff( 34)    *x21        *x53
     8  +coeff( 35)    *x23*x31*x41    
      t_m12_0_3   =t_m12_0_3   
     9  +coeff( 36)    *x23    *x42    
     1  +coeff( 37)    *x21*x32*x42    
     2  +coeff( 38)    *x22*x31*x41*x51
     3  +coeff( 39)        *x33*x41*x51
     4  +coeff( 40)    *x22    *x42*x51
     5  +coeff( 41)        *x32*x42*x51
     6  +coeff( 42)        *x31*x43*x51
     7  +coeff( 43)    *x21*x32    *x52
     8  +coeff( 44)    *x21*x31*x41*x52
      t_m12_0_3   =t_m12_0_3   
     9  +coeff( 45)    *x22        *x53
     1  +coeff( 46)*x11*x21    *x42    
     2  +coeff( 47)    *x22*x32*x42    
     3  +coeff( 48)    *x22*x31*x43    
     4  +coeff( 49)        *x33*x43    
     5  +coeff( 50)    *x23    *x42*x51
     6  +coeff( 51)    *x22*x32    *x52
     7  +coeff( 52)    *x22*x31*x41*x52
     8  +coeff( 53)    *x23        *x53
      t_m12_0_3   =t_m12_0_3   
     9  +coeff( 54)    *x22*x32        
     1  +coeff( 55)        *x33*x41    
     2  +coeff( 56)    *x21*x31*x43    
     3  +coeff( 57)*x11            *x52
     4  +coeff( 58)        *x32    *x53
     5  +coeff( 59)*x11*x23            
     6  +coeff( 60)*x11*x21*x31*x41    
     7  +coeff( 61)    *x22*x33*x41    
     8  +coeff( 62)    *x23*x32    *x51
      t_m12_0_3   =t_m12_0_3   
     9  +coeff( 63)    *x23*x31*x41*x51
     1  +coeff( 64)*x11        *x42*x51
     2  +coeff( 65)*x11*x21        *x52
     3  +coeff( 66)    *x22    *x42*x52
     4  +coeff( 67)        *x32*x42*x52
     5  +coeff( 68)*x11            *x53
     6  +coeff( 69)    *x21*x32    *x53
     7  +coeff( 70)    *x21    *x42*x53
     8  +coeff( 71)    *x23*x32        
      t_m12_0_3   =t_m12_0_3   
     9  +coeff( 72)    *x22*x32    *x51
     1  +coeff( 73)    *x23        *x52
     2  +coeff( 74)    *x21*x32*x42*x51
     3  +coeff( 75)    *x21*x31*x43*x51
     4  +coeff( 76)        *x31*x43*x52
     5  +coeff( 77)    *x21*x31*x41*x53
     6  +coeff( 78)*x11*x22            
     7  +coeff( 79)*x11    *x31*x41    
     8  +coeff( 80)            *x42*x53
      t_m12_0_3   =t_m12_0_3   
     9  +coeff( 81)*x11*x22        *x51
     1  +coeff( 82)*x11    *x32    *x51
     2  +coeff( 83)        *x33*x41*x52
     3  +coeff( 84)*x11*x21            
     4  +coeff( 85)*x11    *x32        
     5  +coeff( 86)*x11        *x42    
     6  +coeff( 87)    *x21    *x42*x52
     7  +coeff( 88)        *x31*x41*x53
     8  +coeff( 89)*x11*x21*x32        
      t_m12_0_3   =t_m12_0_3   
     9  +coeff( 90)    *x21*x33*x41*x51
c
      return
      end
      function y_m12_0_3   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)   !,xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49998E-02,-0.79804E-01,-0.99991E-01,-0.37038E-01,-0.42491E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49981E-02, 0.64043E-01, 0.99991E-01, 0.37038E-01, 0.42498E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.11215901E+00, 0.36877632E+00,-0.26523543E-02,-0.47448864E-02,
     + -0.52749342E-02,-0.15447399E-01,-0.19398802E-02,-0.52673924E-02,
     +  0.22362184E-02, 0.28352651E-02, 0.81062671E-02,-0.30510018E-02,
     + -0.88912370E-02, 0.55380344E-04,-0.10093149E-03,-0.13476883E-02,
     + -0.37587530E-02, 0.18773332E-02,-0.30424385E-03, 0.85088989E-03,
     +  0.25084647E-02,-0.13105881E-02, 0.27418637E-02, 0.33450408E-02,
     +  0.48060855E-02, 0.62946244E-02, 0.65123434E-02,-0.13006441E-02,
     +  0.27069575E-03,-0.24676527E-03, 0.40439071E-03, 0.94565522E-03,
     + -0.17602727E-03, 0.63113036E-03, 0.67351910E-03, 0.13044019E-02,
     + -0.74305040E-04, 0.14910758E-02,-0.63975551E-03, 0.69883390E-03,
     + -0.15875861E-02,-0.15601158E-02,-0.84039907E-03,-0.16409567E-02,
     + -0.39242548E-02,-0.40936009E-02, 0.22587542E-03,-0.37703667E-04,
     + -0.10226020E-02, 0.29674760E-03,-0.22607970E-03,-0.42261137E-04,
     + -0.46832269E-04, 0.16092301E-03, 0.23245635E-03, 0.16518755E-03,
     +  0.65678760E-03,-0.42019298E-03, 0.15778196E-04,-0.12747828E-02,
     + -0.11497163E-02,-0.67681342E-03,-0.14105080E-03,-0.12299359E-03,
     +  0.12235019E-03, 0.21279918E-02, 0.27821565E-02,-0.13621726E-03,
     +  0.12086879E-02,-0.14254630E-02, 0.13166101E-02,-0.62487525E-03,
     +  0.62128709E-03, 0.68727188E-03,-0.63703024E-04, 0.11966147E-03,
     +  0.18400919E-03, 0.67641263E-05,-0.16426471E-03,-0.81283943E-04,
     +  0.91060654E-04,-0.71771748E-04, 0.64407373E-04, 0.82987688E-04,
     +  0.12084757E-03, 0.19098823E-03,-0.19635515E-03, 0.57126152E-04,
     + -0.47908892E-04, 0.10013494E-03,
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
      y_m12_0_3   =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21*x31        
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)        *x31    *x51
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)    *x22*x31        
     8  +coeff(  8)    *x22    *x41    
      y_m12_0_3   =y_m12_0_3   
     9  +coeff(  9)    *x21    *x41*x51
     1  +coeff( 10)        *x31    *x52
     2  +coeff( 11)            *x41*x52
     3  +coeff( 12)    *x23*x31        
     4  +coeff( 13)    *x23    *x41    
     5  +coeff( 14)    *x21*x31*x42    
     6  +coeff( 15)    *x21    *x43    
     7  +coeff( 16)    *x21    *x41*x52
     8  +coeff( 17)            *x41*x53
      y_m12_0_3   =y_m12_0_3   
     9  +coeff( 18)    *x23*x31    *x51
     1  +coeff( 19)            *x43    
     2  +coeff( 20)    *x21*x31    *x51
     3  +coeff( 21)    *x22    *x41*x51
     4  +coeff( 22)        *x31    *x53
     5  +coeff( 23)    *x22*x31*x42    
     6  +coeff( 24)    *x22    *x43    
     7  +coeff( 25)    *x23    *x41*x51
     8  +coeff( 26)    *x23*x31*x42    
      y_m12_0_3   =y_m12_0_3   
     9  +coeff( 27)    *x23    *x43    
     1  +coeff( 28)    *x21*x32*x43    
     2  +coeff( 29)    *x22*x33    *x51
     3  +coeff( 30)        *x31*x42    
     4  +coeff( 31)    *x21*x32*x41    
     5  +coeff( 32)    *x22*x31    *x51
     6  +coeff( 33)    *x21*x31    *x52
     7  +coeff( 34)    *x21*x31*x42*x51
     8  +coeff( 35)    *x21    *x43*x51
      y_m12_0_3   =y_m12_0_3   
     9  +coeff( 36)    *x23*x32*x41    
     1  +coeff( 37)*x11        *x41    
     2  +coeff( 38)    *x22*x32*x41    
     3  +coeff( 39)    *x22    *x41*x52
     4  +coeff( 40)    *x21    *x41*x53
     5  +coeff( 41)    *x22    *x43*x51
     6  +coeff( 42)    *x23    *x41*x52
     7  +coeff( 43)    *x21*x31*x42*x52
     8  +coeff( 44)    *x22*x32*x43    
      y_m12_0_3   =y_m12_0_3   
     9  +coeff( 45)    *x23*x31*x42*x51
     1  +coeff( 46)    *x23    *x43*x51
     2  +coeff( 47)    *x21*x32*x41*x53
     3  +coeff( 48)    *x22*x33*x42*x51
     4  +coeff( 49)*x11*x22    *x41*x52
     5  +coeff( 50)        *x33*x42*x53
     6  +coeff( 51)        *x32*x41    
     7  +coeff( 52)*x11    *x31        
     8  +coeff( 53)    *x21*x33        
      y_m12_0_3   =y_m12_0_3   
     9  +coeff( 54)            *x43*x51
     1  +coeff( 55)        *x32*x43    
     2  +coeff( 56)    *x21*x31    *x53
     3  +coeff( 57)    *x23*x33        
     4  +coeff( 58)    *x21*x33*x42    
     5  +coeff( 59)*x11        *x43    
     6  +coeff( 60)    *x22*x31*x42*x51
     7  +coeff( 61)    *x23*x31    *x52
     8  +coeff( 62)    *x22*x33*x42    
      y_m12_0_3   =y_m12_0_3   
     9  +coeff( 63)    *x23*x32*x41*x51
     1  +coeff( 64)*x11*x21    *x41*x52
     2  +coeff( 65)*x11    *x33*x42    
     3  +coeff( 66)    *x23*x32*x43    
     4  +coeff( 67)    *x23*x31*x42*x52
     5  +coeff( 68)*x11        *x43*x52
     6  +coeff( 69)    *x23    *x43*x52
     7  +coeff( 70)    *x23*x32*x41*x53
     8  +coeff( 71)*x11*x22    *x43*x52
      y_m12_0_3   =y_m12_0_3   
     9  +coeff( 72)*x11*x23*x33    *x52
     1  +coeff( 73)*x11*x21*x33*x42*x52
     2  +coeff( 74)*x11*x21*x32*x43*x52
     3  +coeff( 75)        *x33    *x51
     4  +coeff( 76)        *x31*x42*x51
     5  +coeff( 77)    *x22*x33        
     6  +coeff( 78)*x11        *x41*x51
     7  +coeff( 79)    *x22*x31    *x52
     8  +coeff( 80)        *x33    *x52
      y_m12_0_3   =y_m12_0_3   
     9  +coeff( 81)            *x43*x52
     1  +coeff( 82)*x11*x22    *x41    
     2  +coeff( 83)*x11    *x32*x41    
     3  +coeff( 84)*x11*x21    *x41*x51
     4  +coeff( 85)*x11        *x41*x52
     5  +coeff( 86)    *x21*x32*x41*x52
     6  +coeff( 87)    *x22*x31    *x53
     7  +coeff( 88)*x11*x23*x31        
     8  +coeff( 89)*x11*x21    *x43    
      y_m12_0_3   =y_m12_0_3   
     9  +coeff( 90)*x11*x22    *x41*x51
c
      return
      end
      function p_m12_0_3   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)   !,xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49998E-02,-0.79804E-01,-0.99991E-01,-0.37038E-01,-0.42491E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49981E-02, 0.64043E-01, 0.99991E-01, 0.37038E-01, 0.42498E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.15446544E-01, 0.83851010E-01,-0.36660096E-03,-0.61424536E-03,
     + -0.69020879E-02,-0.20724582E-01,-0.13537668E-02,-0.34176232E-02,
     +  0.69690359E-05, 0.74300398E-04,-0.41418918E-03, 0.37345944E-02,
     +  0.11189468E-01,-0.46965415E-02,-0.12881731E-01,-0.14223122E-02,
     + -0.11647926E-02, 0.52416653E-04,-0.13232361E-03,-0.17494163E-02,
     + -0.53411978E-02, 0.36021308E-02, 0.34555276E-02, 0.22261688E-02,
     +  0.68807444E-02, 0.10571890E-02, 0.11251972E-02, 0.42049909E-04,
     +  0.40153544E-02, 0.12502870E-01, 0.12905304E-01,-0.11227299E-02,
     + -0.32796018E-03,-0.12117581E-02,-0.24925310E-02,-0.38950978E-03,
     + -0.13992330E-02,-0.68615591E-02,-0.66910293E-02,-0.66143075E-05,
     + -0.35993377E-04, 0.87681640E-03, 0.18416821E-02,-0.33262803E-03,
     +  0.90742944E-03,-0.17280737E-02,-0.43511778E-03,-0.12413355E-02,
     +  0.15406894E-03, 0.35809926E-02,-0.90336296E-04,-0.49142557E-03,
     + -0.70885631E-04, 0.34171692E-03, 0.71427065E-04,-0.14044883E-03,
     +  0.30741011E-03, 0.31181358E-03, 0.84444636E-03,-0.26632985E-03,
     +  0.68185698E-04,-0.33745315E-03, 0.10022546E-03, 0.12541690E-03,
     + -0.12934423E-03, 0.44518001E-04,-0.86136139E-03,-0.66483993E-03,
     +  0.65259228E-04, 0.18696710E-03, 0.70805388E-03,-0.23900505E-03,
     + -0.11652465E-02, 0.29810271E-03, 0.31226859E-03,-0.57679415E-03,
     +  0.45296259E-03, 0.53966866E-03, 0.22867694E-02, 0.15164295E-02,
     + -0.62522019E-03, 0.58815046E-03,-0.13973789E-02, 0.10116976E-02,
     +  0.47939343E-03,-0.50544413E-03, 0.56683767E-03,-0.75687625E-03,
     + -0.28303137E-04, 0.34028795E-03,
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
      p_m12_0_3   =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21*x31        
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)        *x31    *x51
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)    *x22*x31        
     8  +coeff(  8)    *x22    *x41    
      p_m12_0_3   =p_m12_0_3   
     9  +coeff(  9)            *x43    
     1  +coeff( 10)    *x21*x31    *x51
     2  +coeff( 11)    *x21    *x41*x51
     3  +coeff( 12)        *x31    *x52
     4  +coeff( 13)            *x41*x52
     5  +coeff( 14)    *x23*x31        
     6  +coeff( 15)    *x23    *x41    
     7  +coeff( 16)    *x21*x31*x42    
     8  +coeff( 17)    *x21    *x43    
      p_m12_0_3   =p_m12_0_3   
     9  +coeff( 18)        *x31*x42*x51
     1  +coeff( 19)            *x43*x51
     2  +coeff( 20)        *x31    *x53
     3  +coeff( 21)            *x41*x53
     4  +coeff( 22)    *x22*x31*x42    
     5  +coeff( 23)    *x22    *x43    
     6  +coeff( 24)    *x23*x31    *x51
     7  +coeff( 25)    *x23    *x41*x51
     8  +coeff( 26)    *x21*x31*x42*x51
      p_m12_0_3   =p_m12_0_3   
     9  +coeff( 27)    *x21    *x43*x51
     1  +coeff( 28)        *x31*x42*x52
     2  +coeff( 29)    *x23*x32*x41    
     3  +coeff( 30)    *x23*x31*x42    
     4  +coeff( 31)    *x23    *x43    
     5  +coeff( 32)    *x21*x32*x43    
     6  +coeff( 33)    *x22*x32*x41*x51
     7  +coeff( 34)    *x22    *x43*x51
     8  +coeff( 35)    *x23    *x41*x52
      p_m12_0_3   =p_m12_0_3   
     9  +coeff( 36)    *x22*x31    *x53
     1  +coeff( 37)    *x22*x33*x42    
     2  +coeff( 38)    *x23*x31*x42*x51
     3  +coeff( 39)    *x23    *x43*x51
     4  +coeff( 40)        *x32*x41    
     5  +coeff( 41)*x11        *x41    
     6  +coeff( 42)    *x22*x31    *x51
     7  +coeff( 43)    *x22    *x41*x51
     8  +coeff( 44)    *x21*x31    *x52
      p_m12_0_3   =p_m12_0_3   
     9  +coeff( 45)    *x22*x32*x41    
     1  +coeff( 46)    *x22*x31*x42*x51
     2  +coeff( 47)    *x22*x32*x43    
     3  +coeff( 48)    *x23*x32*x41*x51
     4  +coeff( 49)*x11    *x32*x43    
     5  +coeff( 50)    *x23    *x43*x52
     6  +coeff( 51)*x11    *x32*x43*x51
     7  +coeff( 52)*x11    *x33*x42*x53
     8  +coeff( 53)    *x21*x33        
      p_m12_0_3   =p_m12_0_3   
     9  +coeff( 54)    *x22*x33        
     1  +coeff( 55)*x11        *x41*x51
     2  +coeff( 56)    *x21*x32*x41*x51
     3  +coeff( 57)            *x43*x52
     4  +coeff( 58)    *x21    *x41*x53
     5  +coeff( 59)    *x23*x33        
     6  +coeff( 60)*x11*x22    *x41    
     7  +coeff( 61)*x11    *x31*x42    
     8  +coeff( 62)    *x21*x33*x42    
      p_m12_0_3   =p_m12_0_3   
     9  +coeff( 63)    *x22*x33    *x51
     1  +coeff( 64)*x11*x21    *x41*x51
     2  +coeff( 65)    *x21*x33    *x52
     3  +coeff( 66)    *x21*x31*x42*x52
     4  +coeff( 67)    *x21    *x43*x52
     5  +coeff( 68)    *x22    *x41*x53
     6  +coeff( 69)*x11    *x33    *x51
     7  +coeff( 70)    *x21*x33*x42*x51
     8  +coeff( 71)    *x21*x32*x43*x51
      p_m12_0_3   =p_m12_0_3   
     9  +coeff( 72)    *x21*x33    *x53
     1  +coeff( 73)    *x23    *x41*x53
     2  +coeff( 74)*x11*x22    *x43    
     3  +coeff( 75)*x11*x23    *x41*x51
     4  +coeff( 76)*x11*x21*x32*x41*x51
     5  +coeff( 77)    *x22*x33*x42*x51
     6  +coeff( 78)    *x22*x32*x43*x51
     7  +coeff( 79)    *x23*x31*x42*x52
     8  +coeff( 80)    *x23*x33*x42*x51
      p_m12_0_3   =p_m12_0_3   
     9  +coeff( 81)    *x22*x32*x43*x52
     1  +coeff( 82)*x11*x22*x31*x42*x52
     2  +coeff( 83)    *x23*x33*x42*x52
     3  +coeff( 84)    *x22*x33*x42*x53
     4  +coeff( 85)*x11*x21*x33*x42*x52
     5  +coeff( 86)*x11    *x32*x43*x53
     6  +coeff( 87)*x11*x23*x33    *x53
     7  +coeff( 88)*x11*x21*x33*x42*x53
     8  +coeff( 89)*x11    *x31        
      p_m12_0_3   =p_m12_0_3   
     9  +coeff( 90)    *x21    *x41*x52
c
      return
      end
      real function l_m12_0_3   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)   !,xmean(10)
      dimension coeff( 41)
      data ncoeff/ 40/
      data avdat/ -0.1421042E-02/
      data xmin/
     1 -0.49998E-02,-0.79804E-01,-0.99991E-01,-0.37038E-01,-0.42491E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49981E-02, 0.64043E-01, 0.99991E-01, 0.37038E-01, 0.42498E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.66023613E-02,-0.45043606E-01,-0.20679636E-01,-0.80873929E-02,
     + -0.82873832E-03,-0.12523599E-02,-0.42215668E-03,-0.56393910E-03,
     + -0.71025366E-03, 0.51278458E-03, 0.10319541E-02, 0.21457032E-03,
     +  0.72494709E-04,-0.38050115E-03, 0.17271418E-03,-0.45643945E-03,
     + -0.34162062E-04,-0.39588902E-03, 0.45996136E-03,-0.80729849E-04,
     +  0.28137304E-03,-0.20636304E-03,-0.27122407E-03,-0.23068873E-03,
     +  0.40777944E-03, 0.71000541E-04,-0.80504338E-04, 0.16922546E-03,
     +  0.62130491E-03,-0.25651918E-03, 0.29775151E-03, 0.22658559E-03,
     + -0.34023353E-03,-0.41299217E-03, 0.43694381E-04,-0.80721729E-04,
     +  0.10190562E-03,-0.19736162E-03, 0.73306851E-05,-0.27326463E-04,
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
      l_m12_0_3   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x22            
     4  +coeff(  4)            *x42    
     5  +coeff(  5)        *x31*x41    
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)*x11                
     8  +coeff(  8)            *x42*x52
      l_m12_0_3   =l_m12_0_3   
     9  +coeff(  9)                *x51
     1  +coeff( 10)        *x31*x41*x51
     2  +coeff( 11)            *x42*x51
     3  +coeff( 12)                *x53
     4  +coeff( 13)                *x52
     5  +coeff( 14)    *x23            
     6  +coeff( 15)    *x21    *x42    
     7  +coeff( 16)    *x22        *x51
     8  +coeff( 17)    *x21        *x52
      l_m12_0_3   =l_m12_0_3   
     9  +coeff( 18)        *x31*x41*x52
     1  +coeff( 19)    *x21        *x53
     2  +coeff( 20)        *x32        
     3  +coeff( 21)    *x21*x31*x41    
     4  +coeff( 22)    *x21    *x42*x51
     5  +coeff( 23)    *x22        *x52
     6  +coeff( 24)                *x54
     7  +coeff( 25)            *x42*x53
     8  +coeff( 26)        *x32    *x51
      l_m12_0_3   =l_m12_0_3   
     9  +coeff( 27)    *x24            
     1  +coeff( 28)    *x22    *x42    
     2  +coeff( 29)    *x23    *x42    
     3  +coeff( 30)    *x21*x31*x43    
     4  +coeff( 31)    *x22        *x53
     5  +coeff( 32)        *x31*x41*x53
     6  +coeff( 33)    *x21        *x54
     7  +coeff( 34)    *x23    *x44    
     8  +coeff( 35)    *x22*x31*x41    
      l_m12_0_3   =l_m12_0_3   
     9  +coeff( 36)    *x21*x31*x41*x51
     1  +coeff( 37)    *x21    *x42*x52
     2  +coeff( 38)            *x42*x54
     3  +coeff( 39)*x11            *x51
     4  +coeff( 40)        *x31*x43    
c
      return
      end
      function x_m12_0_4   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)   !,xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.5920973E-01/
      data xmin/
     1 -0.49998E-02,-0.79858E-01,-0.99988E-01,-0.37038E-01,-0.42461E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49981E-02, 0.64228E-01, 0.99988E-01, 0.37038E-01, 0.42498E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.13447052E-01, 0.42227197E+00, 0.15636854E+00, 0.99194432E-02,
     + -0.49386208E-03,-0.38924341E-02,-0.69175926E-02, 0.15452670E+00,
     + -0.69968358E-01, 0.17833301E-02,-0.20439464E-02, 0.47711952E-03,
     + -0.24212405E-03,-0.25165984E-02, 0.14854626E-03, 0.49570082E-02,
     +  0.95266318E-02,-0.71560316E-01, 0.28047198E-01,-0.28012522E-01,
     + -0.43049321E-01, 0.12193258E-02,-0.22165852E-02,-0.44882377E-02,
     + -0.26389165E-03, 0.26004091E-01,-0.11102702E-01,-0.13944603E-01,
     +  0.13709958E-01, 0.64442592E-03, 0.55428693E-03, 0.52156048E-02,
     +  0.79644797E-03, 0.39282804E-02,-0.23915423E-02,-0.12403495E-01,
     +  0.26316932E-03,-0.15122270E-01, 0.72684865E-02,-0.21931685E-03,
     +  0.70436425E-02,-0.19773778E-02,-0.43157646E-02,-0.32818425E-02,
     +  0.12227091E-01,-0.43650783E-03,-0.64217504E-02, 0.11605751E-01,
     +  0.11739761E-01, 0.16273337E-02,-0.13918323E-02,-0.19872061E-02,
     + -0.40579680E-02,-0.24164448E-03,-0.29738999E-02,-0.31277090E-02,
     +  0.83800908E-02,-0.35926928E-02,-0.29703937E-02, 0.13036701E-02,
     +  0.39889677E-02, 0.61768258E-03, 0.44242102E-02, 0.29920540E-02,
     +  0.12485331E-02, 0.15286504E-01, 0.15644668E-03, 0.27378669E-02,
     + -0.76484616E-03,-0.18032128E-01, 0.71982387E-02, 0.11453943E-01,
     + -0.51681104E-03,-0.13538455E-02, 0.24977780E-03,-0.22107665E-02,
     + -0.25686517E-02, 0.63220540E-03, 0.11876763E-01,-0.12484459E-01,
     +  0.78672118E-03,-0.11808830E-02, 0.30908843E-04,-0.18119803E-01,
     + -0.21483162E-01,-0.11812475E-01, 0.43206871E-02,-0.18415066E-02,
     + -0.36202881E-02,-0.22573934E-02,
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
      x_m12_0_4   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22            
     5  +coeff(  5)        *x32        
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21        *x51
      x_m12_0_4   =x_m12_0_4   
     9  +coeff(  9)                *x52
     1  +coeff( 10)*x11                
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x21*x32        
     4  +coeff( 13)    *x21*x31*x41    
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)    *x22        *x51
     7  +coeff( 16)        *x31*x41*x51
     8  +coeff( 17)            *x42*x51
      x_m12_0_4   =x_m12_0_4   
     9  +coeff( 18)    *x21        *x52
     1  +coeff( 19)                *x53
     2  +coeff( 20)    *x22*x31*x41    
     3  +coeff( 21)    *x22    *x42    
     4  +coeff( 22)*x11            *x51
     5  +coeff( 23)    *x23        *x51
     6  +coeff( 24)    *x22        *x52
     7  +coeff( 25)            *x42*x52
     8  +coeff( 26)    *x21        *x53
      x_m12_0_4   =x_m12_0_4   
     9  +coeff( 27)    *x23*x31*x41    
     1  +coeff( 28)    *x23    *x42    
     2  +coeff( 29)    *x22    *x42*x51
     3  +coeff( 30)    *x23        *x52
     4  +coeff( 31)    *x21*x31*x41*x52
     5  +coeff( 32)    *x22        *x53
     6  +coeff( 33)        *x32    *x53
     7  +coeff( 34)    *x22*x32*x42    
     8  +coeff( 35)    *x22*x32    *x52
      x_m12_0_4   =x_m12_0_4   
     9  +coeff( 36)    *x22*x31*x41*x52
     1  +coeff( 37)        *x33*x41*x52
     2  +coeff( 38)    *x22    *x42*x52
     3  +coeff( 39)    *x23        *x53
     4  +coeff( 40)    *x22*x33*x41*x51
     5  +coeff( 41)    *x22*x31*x41*x53
     6  +coeff( 42)    *x22*x33*x41*x53
     7  +coeff( 43)    *x22*x32        
     8  +coeff( 44)    *x23*x32        
      x_m12_0_4   =x_m12_0_4   
     9  +coeff( 45)    *x22*x31*x41*x51
     1  +coeff( 46)*x11            *x52
     2  +coeff( 47)            *x42*x53
     3  +coeff( 48)    *x22*x31*x43    
     4  +coeff( 49)    *x23*x32*x42    
     5  +coeff( 50)*x11        *x42*x53
     6  +coeff( 51)    *x21*x32    *x51
     7  +coeff( 52)    *x21*x31*x41*x51
     8  +coeff( 53)    *x21    *x42*x51
      x_m12_0_4   =x_m12_0_4   
     9  +coeff( 54)        *x32    *x52
     1  +coeff( 55)    *x21*x32*x42    
     2  +coeff( 56)    *x22*x32    *x51
     3  +coeff( 57)    *x21    *x42*x52
     4  +coeff( 58)        *x31*x41*x53
     5  +coeff( 59)    *x23*x32    *x51
     6  +coeff( 60)    *x23*x31*x41*x51
     7  +coeff( 61)    *x21*x32*x42*x51
     8  +coeff( 62)*x11*x21        *x52
      x_m12_0_4   =x_m12_0_4   
     9  +coeff( 63)    *x21*x32    *x53
     1  +coeff( 64)    *x23*x33*x41    
     2  +coeff( 65)*x11*x21    *x42*x51
     3  +coeff( 66)    *x22*x32*x42*x51
     4  +coeff( 67)*x11    *x32    *x52
     5  +coeff( 68)    *x23*x31*x41*x52
     6  +coeff( 69)*x11        *x42*x52
     7  +coeff( 70)    *x23    *x42*x52
     8  +coeff( 71)    *x22*x32    *x53
      x_m12_0_4   =x_m12_0_4   
     9  +coeff( 72)    *x22    *x42*x53
     1  +coeff( 73)*x11*x23*x32        
     2  +coeff( 74)*x11*x22*x32    *x51
     3  +coeff( 75)*x11*x21*x32    *x52
     4  +coeff( 76)*x11*x21*x31*x41*x52
     5  +coeff( 77)*x11*x21    *x42*x52
     6  +coeff( 78)*x11*x22        *x53
     7  +coeff( 79)    *x23    *x42*x53
     8  +coeff( 80)    *x21*x32*x42*x53
      x_m12_0_4   =x_m12_0_4   
     9  +coeff( 81)    *x21*x31*x43*x53
     1  +coeff( 82)*x11*x23*x32    *x51
     2  +coeff( 83)    *x23*x32*x42*x52
     3  +coeff( 84)    *x23*x31*x43*x52
     4  +coeff( 85)    *x22*x32*x42*x53
     5  +coeff( 86)    *x22*x31*x43*x53
     6  +coeff( 87)*x11*x23*x33*x41    
     7  +coeff( 88)*x11*x23*x32*x42    
     8  +coeff( 89)*x11*x23*x31*x43    
      x_m12_0_4   =x_m12_0_4   
     9  +coeff( 90)*x11*x21*x33*x43    
c
      return
      end
      function t_m12_0_4   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)   !,xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.1402796E-03/
      data xmin/
     1 -0.49998E-02,-0.79858E-01,-0.99988E-01,-0.37038E-01,-0.42461E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49981E-02, 0.64228E-01, 0.99988E-01, 0.37038E-01, 0.42498E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.10252186E-01,-0.97150937E-01, 0.62849030E-01, 0.22149526E-02,
     + -0.19667739E-03,-0.20579242E-02,-0.43295044E-02, 0.62001124E-01,
     + -0.24372758E-01,-0.13178685E-02,-0.15549712E-02,-0.15475608E-03,
     + -0.32925647E-03,-0.16502542E-02, 0.14902814E-02, 0.28071934E-03,
     +  0.26941593E-02, 0.54516671E-02,-0.25088785E-01, 0.83770594E-02,
     + -0.12190904E-01,-0.18193701E-01,-0.52530516E-03, 0.37668229E-03,
     +  0.96895051E-03, 0.68959553E-03,-0.29653900E-02,-0.16809910E-02,
     + -0.31274266E-02, 0.81103332E-02,-0.49477285E-02, 0.32213712E-03,
     +  0.32079578E-02, 0.52472851E-02,-0.13770221E-02,-0.18822885E-02,
     +  0.33410289E-02, 0.33477654E-02,-0.70706842E-03,-0.40926822E-02,
     + -0.47360668E-02, 0.33072091E-02,-0.19119020E-02,-0.40946470E-03,
     + -0.35884455E-02, 0.81734860E-03,-0.17394903E-03,-0.11231698E-02,
     + -0.20056857E-03, 0.51104361E-02,-0.14724494E-02, 0.32691497E-03,
     +  0.26409631E-02,-0.87785121E-03,-0.64025674E-03, 0.36594796E-03,
     +  0.70541115E-04,-0.33920286E-02,-0.36840304E-02,-0.13095851E-02,
     + -0.44984970E-03, 0.12491047E-03, 0.15926339E-03, 0.12907074E-02,
     + -0.35061277E-03, 0.24073409E-03, 0.56007400E-03,-0.51401901E-04,
     +  0.47079775E-04,-0.50609247E-04, 0.18194717E-03, 0.10527499E-02,
     +  0.11393878E-02,-0.75780539E-04, 0.11513623E-02, 0.83595514E-04,
     + -0.96766962E-04,-0.75317826E-03, 0.58241090E-03, 0.83520281E-04,
     + -0.28467061E-05, 0.28355367E-04, 0.29838420E-03, 0.24821138E-03,
     + -0.11470712E-03,-0.13066626E-03,-0.10535546E-03,-0.18939099E-03,
     + -0.26909742E-03,-0.27053140E-03,
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
      t_m12_0_4   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22            
     5  +coeff(  5)        *x32        
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21        *x51
      t_m12_0_4   =t_m12_0_4   
     9  +coeff(  9)                *x52
     1  +coeff( 10)*x11                
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x21*x32        
     4  +coeff( 13)    *x21*x31*x41    
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)    *x22        *x51
     7  +coeff( 16)        *x32    *x51
     8  +coeff( 17)        *x31*x41*x51
      t_m12_0_4   =t_m12_0_4   
     9  +coeff( 18)            *x42*x51
     1  +coeff( 19)    *x21        *x52
     2  +coeff( 20)                *x53
     3  +coeff( 21)    *x22*x31*x41    
     4  +coeff( 22)    *x22    *x42    
     5  +coeff( 23)        *x31*x43    
     6  +coeff( 24)*x11            *x51
     7  +coeff( 25)    *x23        *x51
     8  +coeff( 26)    *x21    *x42*x51
      t_m12_0_4   =t_m12_0_4   
     9  +coeff( 27)    *x22        *x52
     1  +coeff( 28)        *x31*x41*x52
     2  +coeff( 29)            *x42*x52
     3  +coeff( 30)    *x21        *x53
     4  +coeff( 31)    *x23    *x42    
     5  +coeff( 32)    *x21*x31*x43    
     6  +coeff( 33)    *x22*x31*x41*x51
     7  +coeff( 34)    *x22    *x42*x51
     8  +coeff( 35)    *x23        *x52
      t_m12_0_4   =t_m12_0_4   
     9  +coeff( 36)    *x21    *x42*x52
     1  +coeff( 37)    *x22        *x53
     2  +coeff( 38)    *x22*x32*x42    
     3  +coeff( 39)    *x22*x32    *x52
     4  +coeff( 40)    *x22*x31*x41*x52
     5  +coeff( 41)    *x22    *x42*x52
     6  +coeff( 42)    *x23        *x53
     7  +coeff( 43)    *x22*x32        
     8  +coeff( 44)        *x32    *x52
      t_m12_0_4   =t_m12_0_4   
     9  +coeff( 45)    *x23*x31*x41    
     1  +coeff( 46)    *x21*x32*x42    
     2  +coeff( 47)*x11            *x52
     3  +coeff( 48)    *x21*x31*x41*x52
     4  +coeff( 49)*x11*x23            
     5  +coeff( 50)    *x22*x31*x43    
     6  +coeff( 51)    *x23*x32    *x51
     7  +coeff( 52)*x11        *x42*x51
     8  +coeff( 53)    *x21    *x42*x53
      t_m12_0_4   =t_m12_0_4   
     9  +coeff( 54)        *x32*x42    
     1  +coeff( 55)    *x23*x32        
     2  +coeff( 56)        *x32    *x53
     3  +coeff( 57)*x11    *x31*x41*x51
     4  +coeff( 58)    *x23*x31*x41*x51
     5  +coeff( 59)    *x23    *x42*x51
     6  +coeff( 60)    *x21*x32*x42*x51
     7  +coeff( 61)    *x21*x31*x43*x51
     8  +coeff( 62)*x11*x21        *x52
      t_m12_0_4   =t_m12_0_4   
     9  +coeff( 63)*x11            *x53
     1  +coeff( 64)    *x21*x31*x41*x53
     2  +coeff( 65)        *x33*x41    
     3  +coeff( 66)    *x21*x32    *x51
     4  +coeff( 67)    *x21*x31*x41*x51
     5  +coeff( 68)*x11*x22            
     6  +coeff( 69)*x11    *x32        
     7  +coeff( 70)*x11        *x42    
     8  +coeff( 71)        *x33*x41*x51
      t_m12_0_4   =t_m12_0_4   
     9  +coeff( 72)        *x32*x42*x51
     1  +coeff( 73)        *x31*x43*x51
     2  +coeff( 74)    *x21*x32    *x52
     3  +coeff( 75)    *x22*x33*x41    
     4  +coeff( 76)*x11*x22        *x51
     5  +coeff( 77)*x11    *x32    *x51
     6  +coeff( 78)        *x31*x43*x52
     7  +coeff( 79)    *x21*x32    *x53
     8  +coeff( 80)*x11*x21            
      t_m12_0_4   =t_m12_0_4   
     9  +coeff( 81)*x11    *x31*x41    
     1  +coeff( 82)*x11*x21        *x51
     2  +coeff( 83)        *x31*x41*x53
     3  +coeff( 84)            *x42*x53
     4  +coeff( 85)*x11*x21*x32        
     5  +coeff( 86)*x11*x21*x31*x41    
     6  +coeff( 87)*x11*x21    *x42    
     7  +coeff( 88)        *x33*x43    
     8  +coeff( 89)    *x21*x33*x41*x51
      t_m12_0_4   =t_m12_0_4   
     9  +coeff( 90)        *x32*x42*x52
c
      return
      end
      function y_m12_0_4   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)   !,xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49998E-02,-0.79858E-01,-0.99988E-01,-0.37038E-01,-0.42461E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49981E-02, 0.64228E-01, 0.99988E-01, 0.37038E-01, 0.42498E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.17526855E+00, 0.65721214E+00,-0.33807030E-02,-0.20537742E-02,
     + -0.34383308E-01,-0.11021338E+00,-0.53093433E-02,-0.14500580E-01,
     +  0.38681086E-02, 0.13531087E-02, 0.19502181E-02, 0.21465410E-01,
     +  0.67323364E-01,-0.20702083E-01,-0.65233335E-01,-0.72239577E-02,
     + -0.72425823E-02,-0.11802320E-01,-0.34764014E-01, 0.16448895E-01,
     +  0.17163001E-01, 0.72137424E-03, 0.29081276E-01, 0.57363599E-02,
     +  0.69874828E-02, 0.21503543E-01, 0.70932657E-01, 0.72386019E-01,
     + -0.15055229E-01, 0.43017310E-02,-0.88304169E-02,-0.45503498E-03,
     +  0.71112101E-03,-0.70162229E-02, 0.99420745E-03, 0.59050000E-02,
     +  0.12327392E-01,-0.58502360E-03,-0.74402755E-02,-0.82559595E-02,
     +  0.21119516E-01,-0.37516407E-02,-0.23189118E-01,-0.18326543E-01,
     +  0.17818015E-01,-0.42755753E-04, 0.59121312E-02, 0.17346522E-03,
     +  0.34821483E-02,-0.10555463E-02, 0.76137790E-02, 0.20061110E-02,
     +  0.23646142E-01,-0.93419029E-03,-0.63599413E-03, 0.60448807E-03,
     +  0.13004312E-01,-0.27536066E-01,-0.12224112E-01, 0.28523493E-02,
     +  0.28613445E-02,-0.28521875E-02,-0.36288728E-03,-0.18558973E-02,
     + -0.73477640E-02,-0.95841521E-03,-0.18869283E-02,-0.14463738E-01,
     + -0.14386834E-01, 0.12102599E-02,-0.74898393E-03,-0.20031871E-02,
     + -0.12170257E-01, 0.56172197E-03,-0.29288320E-03, 0.49350024E-02,
     +  0.84146094E-02, 0.58971774E-02, 0.37119230E-02, 0.17705874E-03,
     +  0.49590552E-03, 0.16884431E-01, 0.19721435E-01,-0.29177032E-02,
     + -0.11868814E-01,-0.15363569E-01,-0.78530813E-03, 0.26337772E-02,
     + -0.64767017E-02, 0.19834803E-02,
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
      y_m12_0_4   =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21*x31        
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)        *x31    *x51
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)    *x22*x31        
     8  +coeff(  8)    *x22    *x41    
      y_m12_0_4   =y_m12_0_4   
     9  +coeff(  9)            *x43    
     1  +coeff( 10)    *x21*x31    *x51
     2  +coeff( 11)    *x21    *x41*x51
     3  +coeff( 12)        *x31    *x52
     4  +coeff( 13)            *x41*x52
     5  +coeff( 14)    *x23*x31        
     6  +coeff( 15)    *x23    *x41    
     7  +coeff( 16)    *x21*x31*x42    
     8  +coeff( 17)    *x21    *x43    
      y_m12_0_4   =y_m12_0_4   
     9  +coeff( 18)        *x31    *x53
     1  +coeff( 19)            *x41*x53
     2  +coeff( 20)    *x22*x31*x42    
     3  +coeff( 21)    *x22    *x43    
     4  +coeff( 22)        *x32*x43    
     5  +coeff( 23)    *x23    *x41*x51
     6  +coeff( 24)    *x21*x31*x42*x51
     7  +coeff( 25)    *x21    *x43*x51
     8  +coeff( 26)    *x23*x32*x41    
      y_m12_0_4   =y_m12_0_4   
     9  +coeff( 27)    *x23*x31*x42    
     1  +coeff( 28)    *x23    *x43    
     2  +coeff( 29)    *x23    *x41*x52
     3  +coeff( 30)    *x22*x32*x43    
     4  +coeff( 31)    *x23*x32*x41*x51
     5  +coeff( 32)    *x21*x33        
     6  +coeff( 33)    *x22*x31    *x51
     7  +coeff( 34)            *x43*x51
     8  +coeff( 35)    *x21*x31    *x52
      y_m12_0_4   =y_m12_0_4   
     9  +coeff( 36)    *x22*x32*x41    
     1  +coeff( 37)    *x23*x31    *x51
     2  +coeff( 38)    *x22    *x41*x52
     3  +coeff( 39)    *x23*x31    *x52
     4  +coeff( 40)    *x22    *x41*x53
     5  +coeff( 41)            *x43*x53
     6  +coeff( 42)    *x23*x33    *x51
     7  +coeff( 43)    *x23*x31*x42*x51
     8  +coeff( 44)    *x23    *x43*x51
      y_m12_0_4   =y_m12_0_4   
     9  +coeff( 45)    *x23*x33*x42*x52
     1  +coeff( 46)*x11        *x41    
     2  +coeff( 47)    *x22    *x41*x51
     3  +coeff( 48)    *x21    *x41*x53
     4  +coeff( 49)    *x23*x33        
     5  +coeff( 50)*x11*x22    *x41    
     6  +coeff( 51)    *x21    *x43*x52
     7  +coeff( 52)    *x22*x31    *x53
     8  +coeff( 53)        *x31*x42*x53
      y_m12_0_4   =y_m12_0_4   
     9  +coeff( 54)    *x22*x33    *x52
     1  +coeff( 55)    *x21*x32*x41*x53
     2  +coeff( 56)*x11*x21*x33    *x51
     3  +coeff( 57)    *x23    *x43*x52
     4  +coeff( 58)    *x22*x31*x42*x53
     5  +coeff( 59)    *x22*x32*x43*x53
     6  +coeff( 60)*x11*x23*x33    *x53
     7  +coeff( 61)        *x31*x42    
     8  +coeff( 62)    *x21*x32*x41    
      y_m12_0_4   =y_m12_0_4   
     9  +coeff( 63)        *x33    *x51
     1  +coeff( 64)        *x32*x41*x51
     2  +coeff( 65)        *x31*x42*x51
     3  +coeff( 66)    *x21    *x41*x52
     4  +coeff( 67)        *x32*x41*x52
     5  +coeff( 68)        *x31*x42*x52
     6  +coeff( 69)            *x43*x52
     7  +coeff( 70)    *x21*x33*x42    
     8  +coeff( 71)*x11*x21*x31    *x51
      y_m12_0_4   =y_m12_0_4   
     9  +coeff( 72)    *x22*x32*x41*x51
     1  +coeff( 73)    *x22    *x43*x51
     2  +coeff( 74)        *x32*x43*x51
     3  +coeff( 75)    *x21*x33    *x52
     4  +coeff( 76)    *x21*x32*x41*x52
     5  +coeff( 77)    *x21*x31*x42*x52
     6  +coeff( 78)        *x32*x41*x53
     7  +coeff( 79)    *x22*x33*x42    
     8  +coeff( 80)*x11    *x33    *x51
      y_m12_0_4   =y_m12_0_4   
     9  +coeff( 81)*x11*x22    *x41*x51
     1  +coeff( 82)    *x22*x31*x42*x52
     2  +coeff( 83)    *x22    *x43*x52
     3  +coeff( 84)        *x32*x43*x52
     4  +coeff( 85)    *x21*x31*x42*x53
     5  +coeff( 86)    *x21    *x43*x53
     6  +coeff( 87)*x11*x22*x33        
     7  +coeff( 88)*x11*x22*x31*x42    
     8  +coeff( 89)    *x23*x33*x42    
      y_m12_0_4   =y_m12_0_4   
     9  +coeff( 90)*x11*x22    *x43    
c
      return
      end
      function p_m12_0_4   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)   !,xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49998E-02,-0.79858E-01,-0.99988E-01,-0.37038E-01,-0.42461E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49981E-02, 0.64228E-01, 0.99988E-01, 0.37038E-01, 0.42498E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.34747835E-01, 0.15086095E+00,-0.42764922E-02,-0.16882230E-01,
     + -0.55107806E-01,-0.25162508E-03,-0.65643531E-02,-0.28103128E-04,
     + -0.70379674E-03, 0.86910080E-03, 0.74486202E-03, 0.15740438E-02,
     +  0.10245937E-01, 0.36024060E-01,-0.87646320E-02, 0.26605945E-03,
     + -0.30626496E-01,-0.43948302E-02, 0.97169854E-04,-0.21524262E-02,
     + -0.36179267E-02, 0.15947153E-02,-0.59625418E-02,-0.19122254E-01,
     +  0.96077072E-02, 0.85965917E-02, 0.45538247E-02, 0.12669190E-01,
     +  0.91449666E-03,-0.12431690E-02, 0.51233178E-03,-0.39997217E-03,
     +  0.31258069E-01, 0.36794137E-01, 0.10509846E-02,-0.18704036E-02,
     + -0.45742458E-02,-0.14085856E-02,-0.76502203E-02, 0.81020891E-03,
     + -0.59085973E-02,-0.43053529E-02, 0.94432784E-02, 0.29237922E-02,
     + -0.11125061E-01,-0.59552132E-02,-0.10332727E-02, 0.15142772E-02,
     + -0.21921538E-02, 0.43764194E-02,-0.53015852E-03, 0.27660790E-02,
     +  0.45521655E-02, 0.15014264E-01, 0.26173971E-04,-0.82230130E-02,
     + -0.64897933E-04,-0.91014756E-03,-0.22901241E-02,-0.22046592E-02,
     + -0.10337920E-02,-0.24752463E-02, 0.42067879E-03, 0.18158627E-02,
     +  0.39656674E-02, 0.53339638E-03,-0.52437386E-02, 0.90764957E-02,
     + -0.89205854E-03,-0.51734247E-02,-0.33963381E-03, 0.69438391E-02,
     + -0.70010358E-02,-0.56729433E-02,-0.15460222E-03,-0.18258528E-03,
     + -0.61135609E-02, 0.11491596E-01, 0.43830606E-02, 0.75850799E-02,
     +  0.46701400E-03,-0.83743008E-02,-0.84220627E-02, 0.12780910E-02,
     +  0.42142422E-03, 0.13569137E-02, 0.29196721E-02, 0.49567218E-02,
     +  0.11690778E-02,-0.96502498E-03,
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
      p_m12_0_4   =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)        *x31    *x51
     5  +coeff(  5)            *x41*x51
     6  +coeff(  6)        *x33        
     7  +coeff(  7)    *x22    *x41    
     8  +coeff(  8)        *x32*x41    
      p_m12_0_4   =p_m12_0_4   
     9  +coeff(  9)        *x31*x42    
     1  +coeff( 10)            *x43    
     2  +coeff( 11)    *x21*x31    *x51
     3  +coeff( 12)    *x21    *x41*x51
     4  +coeff( 13)        *x31    *x52
     5  +coeff( 14)            *x41*x52
     6  +coeff( 15)    *x23*x31        
     7  +coeff( 16)    *x21*x33        
     8  +coeff( 17)    *x23    *x41    
      p_m12_0_4   =p_m12_0_4   
     9  +coeff( 18)    *x21    *x43    
     1  +coeff( 19)    *x22    *x41*x51
     2  +coeff( 20)        *x31*x42*x51
     3  +coeff( 21)            *x43*x51
     4  +coeff( 22)    *x21*x31    *x52
     5  +coeff( 23)        *x31    *x53
     6  +coeff( 24)            *x41*x53
     7  +coeff( 25)    *x22*x31*x42    
     8  +coeff( 26)    *x22    *x43    
      p_m12_0_4   =p_m12_0_4   
     9  +coeff( 27)    *x23*x31    *x51
     1  +coeff( 28)    *x23    *x41*x51
     2  +coeff( 29)    *x22*x31    *x52
     3  +coeff( 30)    *x21    *x41*x53
     4  +coeff( 31)    *x23*x33        
     5  +coeff( 32)*x11*x22    *x41    
     6  +coeff( 33)    *x23*x31*x42    
     7  +coeff( 34)    *x23    *x43    
     8  +coeff( 35)    *x21*x32*x43    
      p_m12_0_4   =p_m12_0_4   
     9  +coeff( 36)    *x22*x32*x41*x51
     1  +coeff( 37)    *x22*x31*x42*x51
     2  +coeff( 38)    *x21*x33    *x52
     3  +coeff( 39)    *x23    *x41*x52
     4  +coeff( 40)    *x21*x32*x41*x52
     5  +coeff( 41)    *x21*x31*x42*x52
     6  +coeff( 42)    *x22    *x41*x53
     7  +coeff( 43)            *x43*x53
     8  +coeff( 44)    *x22*x32*x43    
      p_m12_0_4   =p_m12_0_4   
     9  +coeff( 45)    *x23*x31*x42*x51
     1  +coeff( 46)    *x23    *x43*x51
     2  +coeff( 47)    *x22*x33    *x52
     3  +coeff( 48)    *x22*x32*x41*x52
     4  +coeff( 49)    *x22*x31*x42*x52
     5  +coeff( 50)    *x22    *x43*x52
     6  +coeff( 51)    *x23*x32*x43    
     7  +coeff( 52)    *x23*x33    *x52
     8  +coeff( 53)    *x23*x32*x41*x52
      p_m12_0_4   =p_m12_0_4   
     9  +coeff( 54)    *x23*x31*x42*x52
     1  +coeff( 55)    *x21*x33*x42*x52
     2  +coeff( 56)    *x21*x32*x43*x52
     3  +coeff( 57)    *x22*x33    *x53
     4  +coeff( 58)    *x23*x33    *x53
     5  +coeff( 59)    *x21*x31        
     6  +coeff( 60)    *x22*x31        
     7  +coeff( 61)    *x21*x32*x41    
     8  +coeff( 62)    *x21*x31*x42    
      p_m12_0_4   =p_m12_0_4   
     9  +coeff( 63)    *x22*x31    *x51
     1  +coeff( 64)    *x22*x32*x41    
     2  +coeff( 65)    *x21    *x43*x51
     3  +coeff( 66)    *x22    *x41*x52
     4  +coeff( 67)            *x43*x52
     5  +coeff( 68)    *x23*x32*x41    
     6  +coeff( 69)    *x22    *x43*x51
     7  +coeff( 70)    *x23*x31    *x52
     8  +coeff( 71)    *x22*x31    *x53
      p_m12_0_4   =p_m12_0_4   
     9  +coeff( 72)        *x31*x42*x53
     1  +coeff( 73)    *x23*x32*x41*x51
     2  +coeff( 74)        *x33*x42*x52
     3  +coeff( 75)*x11*x22*x33        
     4  +coeff( 76)*x11*x21*x33    *x51
     5  +coeff( 77)    *x22*x32*x43*x51
     6  +coeff( 78)    *x23    *x43*x52
     7  +coeff( 79)    *x22*x32*x41*x53
     8  +coeff( 80)    *x22*x33*x42*x52
      p_m12_0_4   =p_m12_0_4   
     9  +coeff( 81)*x11    *x31*x42*x53
     1  +coeff( 82)    *x23    *x43*x53
     2  +coeff( 83)    *x23*x33*x42*x53
     3  +coeff( 84)*x11*x23*x33    *x53
     4  +coeff( 85)        *x32*x41*x51
     5  +coeff( 86)        *x33*x42    
     6  +coeff( 87)    *x21*x32*x41*x51
     7  +coeff( 88)    *x21*x31*x42*x51
     8  +coeff( 89)        *x33    *x52
      p_m12_0_4   =p_m12_0_4   
     9  +coeff( 90)        *x32*x41*x52
c
      return
      end
      real function l_m12_0_4   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)   !,xmean(10)
      dimension coeff( 41)
      data ncoeff/ 40/
      data avdat/ -0.1244915E-01/
      data xmin/
     1 -0.49998E-02,-0.79858E-01,-0.99988E-01,-0.37038E-01,-0.42461E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49981E-02, 0.64228E-01, 0.99988E-01, 0.37038E-01, 0.42498E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.18996032E-01,-0.58241498E-01, 0.84788147E-02,-0.27389666E-01,
     +  0.16658423E-01,-0.27238764E-01,-0.63956059E-02,-0.12742088E-01,
     + -0.93264347E-02, 0.83712237E-02, 0.14759088E-01, 0.58366284E-02,
     + -0.68569053E-02, 0.11920200E-01, 0.16927219E-03, 0.87195709E-02,
     + -0.13799022E-01,-0.67690230E-03,-0.11363636E-02,-0.81090266E-02,
     + -0.37669884E-02,-0.53358668E-04, 0.74272030E-02,-0.70320126E-02,
     + -0.46208203E-02, 0.90648857E-03, 0.60959170E-02, 0.79616848E-02,
     + -0.16606657E-02, 0.51986752E-02,-0.25705170E-03, 0.30257727E-03,
     +  0.16086003E-02,-0.21278749E-02,-0.11419294E-02,-0.15364368E-02,
     +  0.32720317E-02,-0.44250651E-02,-0.49318383E-02,-0.39145566E-03,
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
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      l_m12_0_4   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)            *x42    
     7  +coeff(  7)                *x52
     8  +coeff(  8)    *x21        *x52
      l_m12_0_4   =l_m12_0_4   
     9  +coeff(  9)        *x31*x41    
     1  +coeff( 10)    *x22        *x51
     2  +coeff( 11)            *x42*x51
     3  +coeff( 12)                *x53
     4  +coeff( 13)    *x22        *x52
     5  +coeff( 14)    *x21        *x53
     6  +coeff( 15)        *x33*x41*x51
     7  +coeff( 16)        *x31*x41*x51
     8  +coeff( 17)            *x42*x52
      l_m12_0_4   =l_m12_0_4   
     9  +coeff( 18)*x11                
     1  +coeff( 19)        *x32        
     2  +coeff( 20)        *x31*x41*x52
     3  +coeff( 21)                *x54
     4  +coeff( 22)        *x34    *x51
     5  +coeff( 23)    *x22        *x53
     6  +coeff( 24)    *x21        *x54
     7  +coeff( 25)    *x22    *x42    
     8  +coeff( 26)    *x23        *x51
      l_m12_0_4   =l_m12_0_4   
     9  +coeff( 27)    *x22    *x42*x51
     1  +coeff( 28)            *x42*x53
     2  +coeff( 29)    *x24*x31*x41    
     3  +coeff( 30)    *x22*x31*x41*x53
     4  +coeff( 31)*x11*x21            
     5  +coeff( 32)*x11            *x51
     6  +coeff( 33)        *x32    *x51
     7  +coeff( 34)    *x22*x31*x41    
     8  +coeff( 35)        *x32    *x52
      l_m12_0_4   =l_m12_0_4   
     9  +coeff( 36)    *x23        *x52
     1  +coeff( 37)        *x31*x41*x53
     2  +coeff( 38)    *x22        *x54
     3  +coeff( 39)    *x24    *x42*x52
     4  +coeff( 40)    *x21    *x42    
c
      return
      end
      function x_m12_0_5   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)   !,xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.6877057E-01/
      data xmin/
     1 -0.49998E-02,-0.79858E-01,-0.99988E-01,-0.37028E-01,-0.42461E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49981E-02, 0.64228E-01, 0.99988E-01, 0.37028E-01, 0.42498E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.25744922E-01, 0.39696348E+00, 0.17689392E+00, 0.29372363E-03,
     + -0.55592315E-03,-0.32807051E-02,-0.66929883E-02, 0.17698520E+00,
     + -0.73842861E-01, 0.14075277E-02,-0.27265870E-02, 0.22183127E-03,
     + -0.19992343E-02,-0.54136016E-02, 0.27931449E-02, 0.51728995E-02,
     +  0.92078652E-02,-0.72454683E-01, 0.28053202E-01,-0.34305722E-01,
     + -0.52379835E-01, 0.14304443E-02,-0.12697461E-02,-0.19976052E-02,
     + -0.83811581E-03, 0.24294194E-01,-0.13890869E-01, 0.99413013E-02,
     +  0.15022238E-01, 0.16686744E-03, 0.40485831E-02, 0.19888347E-03,
     +  0.44371518E-02,-0.27718947E-02,-0.11498498E-01,-0.12782885E-01,
     +  0.64221374E-02,-0.10531803E-03,-0.46701315E-02,-0.86330064E-02,
     + -0.28563691E-02,-0.22018660E-03,-0.62648067E-02, 0.97883102E-02,
     +  0.28271712E-02, 0.14111228E-01, 0.30247071E-02,-0.10260072E-02,
     + -0.44527357E-02,-0.57287496E-02,-0.23073082E-03,-0.37058324E-03,
     + -0.30678662E-02, 0.36268789E-03,-0.79826750E-04, 0.69663003E-02,
     + -0.40935213E-02,-0.46980203E-03, 0.19896573E-02,-0.41525182E-03,
     + -0.17361883E-02, 0.11615799E-02, 0.77063503E-03, 0.25454685E-02,
     +  0.66219983E-02,-0.96098514E-03, 0.10053809E-01, 0.68897258E-04,
     + -0.23192051E-02,-0.37169482E-02, 0.51642121E-02,-0.35514850E-02,
     + -0.54229321E-02, 0.12110804E-01,-0.12186612E-02,-0.69002965E-02,
     +  0.14494224E-01,-0.12638820E-02,-0.11536038E-01, 0.25356568E-02,
     + -0.10608708E-02,-0.15211819E-01,-0.20062276E-02,-0.13495788E-01,
     +  0.14391719E-01, 0.26182071E-02, 0.20327521E-01, 0.46899807E-03,
     + -0.81345922E-03,-0.22133645E-03,
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
      x_m12_0_5   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22            
     5  +coeff(  5)        *x32        
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21        *x51
      x_m12_0_5   =x_m12_0_5   
     9  +coeff(  9)                *x52
     1  +coeff( 10)*x11                
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x21*x32        
     4  +coeff( 13)    *x21*x31*x41    
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)    *x22        *x51
     7  +coeff( 16)        *x31*x41*x51
     8  +coeff( 17)            *x42*x51
      x_m12_0_5   =x_m12_0_5   
     9  +coeff( 18)    *x21        *x52
     1  +coeff( 19)                *x53
     2  +coeff( 20)    *x22*x31*x41    
     3  +coeff( 21)    *x22    *x42    
     4  +coeff( 22)*x11            *x51
     5  +coeff( 23)    *x23        *x51
     6  +coeff( 24)    *x22        *x52
     7  +coeff( 25)            *x42*x52
     8  +coeff( 26)    *x21        *x53
      x_m12_0_5   =x_m12_0_5   
     9  +coeff( 27)    *x23    *x42    
     1  +coeff( 28)    *x22*x31*x41*x51
     2  +coeff( 29)    *x22    *x42*x51
     3  +coeff( 30)    *x23        *x52
     4  +coeff( 31)    *x22        *x53
     5  +coeff( 32)        *x32    *x53
     6  +coeff( 33)    *x22*x32*x42    
     7  +coeff( 34)    *x22*x32    *x52
     8  +coeff( 35)    *x22*x31*x41*x52
      x_m12_0_5   =x_m12_0_5   
     9  +coeff( 36)    *x22    *x42*x52
     1  +coeff( 37)    *x23        *x53
     2  +coeff( 38)*x11*x21            
     3  +coeff( 39)    *x22*x32        
     4  +coeff( 40)    *x23*x31*x41    
     5  +coeff( 41)    *x21*x32*x42    
     6  +coeff( 42)*x11            *x52
     7  +coeff( 43)            *x42*x53
     8  +coeff( 44)    *x22*x31*x43    
      x_m12_0_5   =x_m12_0_5   
     9  +coeff( 45)    *x22*x32    *x53
     1  +coeff( 46)    *x22    *x42*x53
     2  +coeff( 47)*x11        *x42*x53
     3  +coeff( 48)    *x21*x32    *x51
     4  +coeff( 49)    *x21*x31*x41*x51
     5  +coeff( 50)    *x21    *x42*x51
     6  +coeff( 51)        *x32    *x52
     7  +coeff( 52)        *x31*x41*x52
     8  +coeff( 53)    *x23*x32        
      x_m12_0_5   =x_m12_0_5   
     9  +coeff( 54)*x11        *x42    
     1  +coeff( 55)        *x33*x41*x51
     2  +coeff( 56)    *x21    *x42*x52
     3  +coeff( 57)        *x31*x41*x53
     4  +coeff( 58)*x11*x23            
     5  +coeff( 59)    *x22*x33*x41    
     6  +coeff( 60)*x11    *x32    *x51
     7  +coeff( 61)    *x23*x32    *x51
     8  +coeff( 62)*x11    *x31*x41*x51
      x_m12_0_5   =x_m12_0_5   
     9  +coeff( 63)*x11*x21        *x52
     1  +coeff( 64)    *x21*x32    *x53
     2  +coeff( 65)    *x21*x31*x41*x53
     3  +coeff( 66)    *x21    *x42*x53
     4  +coeff( 67)    *x23*x32*x42    
     5  +coeff( 68)*x11    *x31*x43    
     6  +coeff( 69)    *x23*x31*x43    
     7  +coeff( 70)    *x22*x33*x41*x51
     8  +coeff( 71)    *x23*x31*x41*x52
      x_m12_0_5   =x_m12_0_5   
     9  +coeff( 72)*x11        *x42*x52
     1  +coeff( 73)    *x23    *x42*x52
     2  +coeff( 74)    *x22*x31*x41*x53
     3  +coeff( 75)    *x22*x32*x42*x52
     4  +coeff( 76)    *x23*x31*x41*x53
     5  +coeff( 77)    *x23    *x42*x53
     6  +coeff( 78)*x11*x23*x32    *x51
     7  +coeff( 79)*x11*x22*x31*x41*x52
     8  +coeff( 80)*x11    *x32*x42*x52
      x_m12_0_5   =x_m12_0_5   
     9  +coeff( 81)*x11*x21*x31*x41*x53
     1  +coeff( 82)    *x23*x32*x42*x53
     2  +coeff( 83)*x11    *x31*x43*x53
     3  +coeff( 84)    *x21*x33*x43*x53
     4  +coeff( 85)*x11*x22*x33*x41*x52
     5  +coeff( 86)*x11*x23    *x42*x53
     6  +coeff( 87)    *x23*x33*x43*x53
     7  +coeff( 88)        *x32    *x51
     8  +coeff( 89)        *x31*x43    
      x_m12_0_5   =x_m12_0_5   
     9  +coeff( 90)*x11*x22            
c
      return
      end
      function t_m12_0_5   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)   !,xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.3190489E-02/
      data xmin/
     1 -0.49998E-02,-0.79858E-01,-0.99988E-01,-0.37028E-01,-0.42461E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49981E-02, 0.64228E-01, 0.99988E-01, 0.37028E-01, 0.42498E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.69848271E-02,-0.97335577E-01, 0.66031821E-01, 0.43924018E-02,
     +  0.49954874E-03, 0.31321496E-02, 0.55890540E-02, 0.63914061E-01,
     + -0.24970712E-01,-0.13195669E-02,-0.18662683E-02,-0.43351937E-03,
     + -0.26082883E-02,-0.55388599E-02, 0.23106653E-02,-0.42686489E-03,
     + -0.26257036E-02,-0.44325781E-02,-0.25487980E-01, 0.80360211E-02,
     + -0.13565922E-01,-0.21351893E-01, 0.34043918E-03, 0.95102109E-03,
     + -0.52987953E-03,-0.40115560E-02, 0.52939621E-02, 0.78465128E-02,
     + -0.34858289E-02,-0.16181704E-02, 0.35964048E-02,-0.38339752E-02,
     +  0.37153410E-02,-0.31144434E-03,-0.20579556E-02, 0.15919025E-03,
     +  0.25329200E-03, 0.31214969E-02,-0.21623510E-02,-0.48385895E-03,
     + -0.10048938E-02, 0.27960902E-02,-0.24564436E-02, 0.94891456E-03,
     +  0.53028669E-02,-0.13834187E-03,-0.20869384E-02, 0.34643270E-02,
     + -0.65284135E-03, 0.13593155E-03,-0.79251084E-04,-0.51836570E-03,
     +  0.66384753E-04, 0.36666574E-04, 0.69813621E-04, 0.34246584E-02,
     + -0.66918947E-04, 0.11269693E-02,-0.19460121E-03, 0.16020262E-02,
     + -0.12010286E-02,-0.17450457E-02, 0.15789302E-03,-0.21847489E-02,
     +  0.19443098E-03, 0.66067063E-03, 0.76424556E-04, 0.23519270E-03,
     +  0.40895109E-04, 0.86969452E-03,-0.13294682E-03, 0.14947669E-03,
     + -0.59515925E-04, 0.26985069E-03,-0.11244027E-02,-0.18469850E-02,
     + -0.16926079E-02, 0.11215348E-02, 0.15787259E-02,-0.16968135E-04,
     + -0.72448187E-04, 0.12124984E-03,-0.56997535E-04,-0.83214996E-04,
     + -0.33282395E-04, 0.24760779E-03, 0.97724864E-04,-0.38725440E-03,
     +  0.20626748E-03, 0.98456992E-04,
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
      t_m12_0_5   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22            
     5  +coeff(  5)        *x32        
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21        *x51
      t_m12_0_5   =t_m12_0_5   
     9  +coeff(  9)                *x52
     1  +coeff( 10)*x11                
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x21*x32        
     4  +coeff( 13)    *x21*x31*x41    
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)    *x22        *x51
     7  +coeff( 16)        *x32    *x51
     8  +coeff( 17)        *x31*x41*x51
      t_m12_0_5   =t_m12_0_5   
     9  +coeff( 18)            *x42*x51
     1  +coeff( 19)    *x21        *x52
     2  +coeff( 20)                *x53
     3  +coeff( 21)    *x22*x31*x41    
     4  +coeff( 22)    *x22    *x42    
     5  +coeff( 23)*x11            *x51
     6  +coeff( 24)    *x23        *x51
     7  +coeff( 25)    *x21    *x42*x51
     8  +coeff( 26)    *x22        *x52
      t_m12_0_5   =t_m12_0_5   
     9  +coeff( 27)            *x42*x52
     1  +coeff( 28)    *x21        *x53
     2  +coeff( 29)    *x23    *x42    
     3  +coeff( 30)    *x23        *x52
     4  +coeff( 31)    *x22        *x53
     5  +coeff( 32)            *x42*x53
     6  +coeff( 33)    *x22*x32*x42    
     7  +coeff( 34)    *x22*x32    *x52
     8  +coeff( 35)    *x22*x31*x41*x52
      t_m12_0_5   =t_m12_0_5   
     9  +coeff( 36)        *x33*x41*x52
     1  +coeff( 37)        *x32*x42*x52
     2  +coeff( 38)    *x23        *x53
     3  +coeff( 39)    *x22*x32        
     4  +coeff( 40)        *x33*x41    
     5  +coeff( 41)        *x32*x42    
     6  +coeff( 42)        *x31*x41*x52
     7  +coeff( 43)    *x23*x31*x41    
     8  +coeff( 44)    *x21*x32*x42    
      t_m12_0_5   =t_m12_0_5   
     9  +coeff( 45)    *x22    *x42*x51
     1  +coeff( 46)*x11            *x52
     2  +coeff( 47)        *x31*x41*x53
     3  +coeff( 48)    *x22*x31*x43    
     4  +coeff( 49)        *x31*x43    
     5  +coeff( 50)    *x21*x32    *x51
     6  +coeff( 51)*x11*x22            
     7  +coeff( 52)    *x23*x32        
     8  +coeff( 53)    *x21*x31*x43    
      t_m12_0_5   =t_m12_0_5   
     9  +coeff( 54)*x11*x21        *x51
     1  +coeff( 55)    *x22*x32    *x51
     2  +coeff( 56)    *x22*x31*x41*x51
     3  +coeff( 57)        *x33*x41*x51
     4  +coeff( 58)        *x31*x43*x51
     5  +coeff( 59)*x11*x23            
     6  +coeff( 60)    *x22*x33*x41    
     7  +coeff( 61)    *x23*x32    *x51
     8  +coeff( 62)    *x23*x31*x41*x51
      t_m12_0_5   =t_m12_0_5   
     9  +coeff( 63)*x11*x21        *x52
     1  +coeff( 64)    *x22    *x42*x52
     2  +coeff( 65)*x11            *x53
     3  +coeff( 66)    *x21*x32    *x53
     4  +coeff( 67)*x11*x21            
     5  +coeff( 68)        *x32    *x52
     6  +coeff( 69)*x11    *x32        
     7  +coeff( 70)        *x32*x42*x51
     8  +coeff( 71)*x11*x21*x32        
      t_m12_0_5   =t_m12_0_5   
     9  +coeff( 72)*x11*x22        *x51
     1  +coeff( 73)*x11    *x32    *x51
     2  +coeff( 74)*x11        *x42*x51
     3  +coeff( 75)    *x23    *x42*x51
     4  +coeff( 76)    *x21*x32*x42*x51
     5  +coeff( 77)    *x21*x31*x43*x51
     6  +coeff( 78)    *x21*x31*x41*x53
     7  +coeff( 79)    *x21    *x42*x53
     8  +coeff( 80)*x11    *x31*x41    
      t_m12_0_5   =t_m12_0_5   
     9  +coeff( 81)*x11        *x42    
     1  +coeff( 82)    *x21*x31*x41*x52
     2  +coeff( 83)        *x32    *x53
     3  +coeff( 84)*x11*x21*x31*x41    
     4  +coeff( 85)*x11*x21    *x42    
     5  +coeff( 86)        *x33*x43    
     6  +coeff( 87)*x11    *x31*x41*x51
     7  +coeff( 88)    *x21*x33*x41*x51
     8  +coeff( 89)    *x21*x31*x41*x51
      t_m12_0_5   =t_m12_0_5   
     9  +coeff( 90)    *x21*x33*x41    
c
      return
      end
      function y_m12_0_5   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)   !,xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49998E-02,-0.79858E-01,-0.99988E-01,-0.37028E-01,-0.42461E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49981E-02, 0.64228E-01, 0.99988E-01, 0.37028E-01, 0.42498E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.18539825E+00, 0.70100892E+00,-0.59795281E-03, 0.11215420E-01,
     + -0.38225625E-01,-0.12029304E+00,-0.70314147E-02,-0.18928112E-01,
     +  0.36681211E-03, 0.32115937E-02, 0.22482952E-01, 0.71171261E-01,
     + -0.23154559E-01,-0.74435540E-01,-0.74463682E-02,-0.53661414E-02,
     + -0.35619508E-01, 0.24111573E-01, 0.22857012E-01, 0.12823098E-01,
     +  0.31771727E-01, 0.78953825E-01, 0.82265891E-01,-0.15453338E-02,
     + -0.16603502E-01, 0.44426611E-02, 0.64791255E-02,-0.24131048E-03,
     +  0.10942570E-02, 0.35458780E-02, 0.10268630E-01, 0.47513316E-02,
     + -0.25538485E-02,-0.22674452E-02, 0.93585895E-02, 0.15992368E-01,
     +  0.27691668E-02,-0.23680222E-02, 0.30072783E-02,-0.78009791E-03,
     +  0.48827380E-04,-0.11934044E-01, 0.52852244E-02,-0.29762110E-02,
     + -0.10183830E-02, 0.23777718E-01,-0.75952932E-02,-0.77009080E-02,
     +  0.19323468E-01,-0.31663396E-02, 0.23462038E-01, 0.15155886E-02,
     +  0.10862404E-02,-0.31836209E-03,-0.74693570E-02,-0.30163038E-02,
     + -0.14120308E-01,-0.32629012E-03, 0.26869669E-02,-0.22043584E-01,
     + -0.17940043E-01,-0.61661325E-03, 0.32223540E-02, 0.43461714E-02,
     +  0.17470315E-01, 0.17529385E-01,-0.15427055E-02, 0.15235232E-02,
     + -0.58439290E-02, 0.19051017E-01, 0.66107712E-02,-0.30944977E-01,
     + -0.32241482E-01, 0.27729166E-02, 0.96272503E-04, 0.59844535E-02,
     + -0.31133157E-02, 0.28568128E-03, 0.63552620E-03, 0.15857054E-02,
     +  0.10056166E-03,-0.39562373E-02,-0.12193805E-01, 0.82022999E-03,
     +  0.81019179E-03,-0.51835779E-03,-0.38931012E-03,-0.59822742E-02,
     +  0.40149101E-03,-0.62373499E-02,
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
      y_m12_0_5   =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21*x31        
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)        *x31    *x51
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)    *x22*x31        
     8  +coeff(  8)    *x22    *x41    
      y_m12_0_5   =y_m12_0_5   
     9  +coeff(  9)        *x32*x41    
     1  +coeff( 10)            *x43    
     2  +coeff( 11)        *x31    *x52
     3  +coeff( 12)            *x41*x52
     4  +coeff( 13)    *x23*x31        
     5  +coeff( 14)    *x23    *x41    
     6  +coeff( 15)    *x21    *x43    
     7  +coeff( 16)            *x43*x51
     8  +coeff( 17)            *x41*x53
      y_m12_0_5   =y_m12_0_5   
     9  +coeff( 18)    *x22*x31*x42    
     1  +coeff( 19)    *x22    *x43    
     2  +coeff( 20)    *x23*x31    *x51
     3  +coeff( 21)    *x23    *x41*x51
     4  +coeff( 22)    *x23*x31*x42    
     5  +coeff( 23)    *x23    *x43    
     6  +coeff( 24)    *x21*x32*x43    
     7  +coeff( 25)    *x23    *x41*x52
     8  +coeff( 26)    *x21*x32*x41*x52
      y_m12_0_5   =y_m12_0_5   
     9  +coeff( 27)    *x21*x31*x42*x52
     1  +coeff( 28)    *x22*x31    *x53
     2  +coeff( 29)        *x33    *x53
     3  +coeff( 30)    *x22*x32*x43    
     4  +coeff( 31)    *x22*x32*x41*x52
     5  +coeff( 32)    *x23*x32*x43    
     6  +coeff( 33)    *x23*x32*x41*x52
     7  +coeff( 34)    *x22*x33    *x53
     8  +coeff( 35)    *x22*x33*x42*x52
      y_m12_0_5   =y_m12_0_5   
     9  +coeff( 36)    *x23*x33*x42*x52
     1  +coeff( 37)    *x21    *x41*x51
     2  +coeff( 38)    *x21*x32*x41    
     3  +coeff( 39)    *x22*x31    *x51
     4  +coeff( 40)        *x33    *x51
     5  +coeff( 41)    *x21*x31    *x52
     6  +coeff( 42)        *x31    *x53
     7  +coeff( 43)    *x22*x32*x41    
     8  +coeff( 44)    *x22    *x41*x52
      y_m12_0_5   =y_m12_0_5   
     9  +coeff( 45)*x11*x22    *x41    
     1  +coeff( 46)    *x23*x32*x41    
     2  +coeff( 47)    *x23*x31    *x52
     3  +coeff( 48)    *x22    *x41*x53
     4  +coeff( 49)            *x43*x53
     5  +coeff( 50)    *x23*x33    *x51
     6  +coeff( 51)    *x22    *x43*x52
     7  +coeff( 52)        *x31*x42    
     8  +coeff( 53)    *x21*x31    *x51
      y_m12_0_5   =y_m12_0_5   
     9  +coeff( 54)    *x21*x33        
     1  +coeff( 55)    *x21*x31*x42    
     2  +coeff( 56)    *x21    *x41*x52
     3  +coeff( 57)            *x43*x52
     4  +coeff( 58)*x11*x22*x31        
     5  +coeff( 59)    *x23*x33        
     6  +coeff( 60)    *x22*x31*x42*x51
     7  +coeff( 61)    *x22    *x43*x51
     8  +coeff( 62)        *x32*x43*x51
      y_m12_0_5   =y_m12_0_5   
     9  +coeff( 63)    *x21    *x43*x52
     1  +coeff( 64)        *x32*x41*x53
     2  +coeff( 65)        *x31*x42*x53
     3  +coeff( 66)    *x22*x31*x42*x52
     4  +coeff( 67)        *x33*x42*x52
     5  +coeff( 68)    *x22*x33*x42*x51
     6  +coeff( 69)    *x22*x32*x43*x51
     7  +coeff( 70)    *x23    *x43*x52
     8  +coeff( 71)    *x22    *x43*x53
      y_m12_0_5   =y_m12_0_5   
     9  +coeff( 72)    *x23*x31*x42*x53
     1  +coeff( 73)    *x23    *x43*x53
     2  +coeff( 74)*x11*x23*x33    *x53
     3  +coeff( 75)        *x33        
     4  +coeff( 76)    *x22    *x41*x51
     5  +coeff( 77)        *x31*x42*x51
     6  +coeff( 78)        *x33*x42    
     7  +coeff( 79)    *x21*x32*x41*x51
     8  +coeff( 80)    *x21    *x43*x51
      y_m12_0_5   =y_m12_0_5   
     9  +coeff( 81)        *x33    *x52
     1  +coeff( 82)        *x32*x41*x52
     2  +coeff( 83)        *x31*x42*x52
     3  +coeff( 84)    *x21*x31    *x53
     4  +coeff( 85)    *x21    *x41*x53
     5  +coeff( 86)*x11*x21*x31    *x51
     6  +coeff( 87)*x11*x21    *x41*x51
     7  +coeff( 88)    *x22*x32*x41*x51
     8  +coeff( 89)*x11*x22    *x41*x51
      y_m12_0_5   =y_m12_0_5   
     9  +coeff( 90)    *x23*x32*x41*x51
c
      return
      end
      function p_m12_0_5   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)   !,xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49998E-02,-0.79858E-01,-0.99988E-01,-0.37028E-01,-0.42461E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49981E-02, 0.64228E-01, 0.99988E-01, 0.37028E-01, 0.42498E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.33478376E-01, 0.14701878E+00,-0.38252557E-02,-0.11393061E-01,
     + -0.16103728E-01,-0.52246708E-01,-0.12928060E-02, 0.46291421E-04,
     + -0.37793221E-02, 0.57377951E-03, 0.29875040E-02, 0.41686986E-02,
     +  0.45975004E-02, 0.11599871E-01, 0.37405774E-01,-0.64944339E-04,
     + -0.95550111E-02,-0.32008041E-01,-0.19449636E-02,-0.57851961E-02,
     + -0.69223437E-02,-0.27788966E-03,-0.25336284E-02,-0.14936682E-02,
     + -0.70599564E-02,-0.88050077E-02,-0.13563334E-02,-0.67836526E-02,
     + -0.22402190E-01, 0.66850157E-02, 0.54894970E-02, 0.10143383E-02,
     +  0.13044340E-01, 0.51015355E-02, 0.59703439E-02, 0.10907952E-02,
     +  0.68262621E-03, 0.46091774E-03, 0.10426855E-02,-0.44675698E-03,
     +  0.36054540E-01, 0.38808428E-01,-0.30063184E-02,-0.69556129E-02,
     + -0.81396865E-03, 0.12001183E-01, 0.11000737E-02, 0.13577649E-02,
     + -0.22751691E-02,-0.30837581E-02, 0.13994124E-02, 0.44684019E-03,
     + -0.59211144E-03, 0.27282226E-02, 0.46828066E-03, 0.49352511E-02,
     +  0.15405579E-01, 0.21174127E-01,-0.39302083E-03, 0.69740601E-03,
     + -0.98706698E-02, 0.10744936E-02, 0.43261731E-04, 0.66039674E-02,
     +  0.23114961E-02, 0.10337013E-01, 0.44262293E-03,-0.47670619E-03,
     +  0.14681900E-02, 0.10463438E-02, 0.23902492E-02, 0.99823065E-02,
     + -0.54423901E-03,-0.17002604E-02,-0.38407096E-02,-0.70930514E-02,
     +  0.77284448E-03,-0.31879488E-01,-0.28459575E-01,-0.14140758E-02,
     +  0.15052905E-02,-0.70337388E-04,-0.42892387E-03, 0.16486809E-02,
     +  0.29471438E-03,-0.22850161E-04,-0.32269722E-02,-0.39944313E-02,
     + -0.26215380E-03,-0.12298864E-02,
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
      p_m12_0_5   =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21*x31        
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)        *x31    *x51
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)    *x22*x31        
     8  +coeff(  8)        *x33        
      p_m12_0_5   =p_m12_0_5   
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)        *x32*x41    
     2  +coeff( 11)        *x31*x42    
     3  +coeff( 12)            *x43    
     4  +coeff( 13)    *x21    *x41*x51
     5  +coeff( 14)        *x31    *x52
     6  +coeff( 15)            *x41*x52
     7  +coeff( 16)*x11    *x31        
     8  +coeff( 17)    *x23*x31        
      p_m12_0_5   =p_m12_0_5   
     9  +coeff( 18)    *x23    *x41    
     1  +coeff( 19)    *x21*x32*x41    
     2  +coeff( 20)    *x21*x31*x42    
     3  +coeff( 21)    *x21    *x43    
     4  +coeff( 22)        *x33    *x51
     5  +coeff( 23)    *x22    *x41*x51
     6  +coeff( 24)        *x32*x41*x51
     7  +coeff( 25)        *x31*x42*x51
     8  +coeff( 26)            *x43*x51
      p_m12_0_5   =p_m12_0_5   
     9  +coeff( 27)    *x21    *x41*x52
     1  +coeff( 28)        *x31    *x53
     2  +coeff( 29)            *x41*x53
     3  +coeff( 30)    *x22    *x43    
     4  +coeff( 31)    *x23*x31    *x51
     5  +coeff( 32)    *x21*x33    *x51
     6  +coeff( 33)    *x23    *x41*x51
     7  +coeff( 34)    *x21*x31*x42*x51
     8  +coeff( 35)    *x21    *x43*x51
      p_m12_0_5   =p_m12_0_5   
     9  +coeff( 36)    *x22    *x41*x52
     1  +coeff( 37)    *x21*x31    *x53
     2  +coeff( 38)    *x21    *x41*x53
     3  +coeff( 39)    *x23*x33        
     4  +coeff( 40)*x11*x22    *x41    
     5  +coeff( 41)    *x23*x31*x42    
     6  +coeff( 42)    *x23    *x43    
     7  +coeff( 43)    *x23*x31    *x52
     8  +coeff( 44)    *x23    *x41*x52
      p_m12_0_5   =p_m12_0_5   
     9  +coeff( 45)    *x22    *x41*x53
     1  +coeff( 46)            *x43*x53
     2  +coeff( 47)    *x22*x33*x42    
     3  +coeff( 48)    *x22*x32*x43    
     4  +coeff( 49)    *x23*x33    *x51
     5  +coeff( 50)    *x23*x32*x41*x51
     6  +coeff( 51)    *x23*x31*x42*x51
     7  +coeff( 52)    *x23    *x43*x51
     8  +coeff( 53)    *x21*x33    *x53
      p_m12_0_5   =p_m12_0_5   
     9  +coeff( 54)    *x23*x32*x43    
     1  +coeff( 55)*x11*x21*x33    *x51
     2  +coeff( 56)    *x23*x32*x41*x52
     3  +coeff( 57)    *x23*x31*x42*x52
     4  +coeff( 58)    *x23    *x43*x52
     5  +coeff( 59)    *x22*x33    *x53
     6  +coeff( 60)        *x33*x42*x53
     7  +coeff( 61)    *x23*x32*x43*x52
     8  +coeff( 62)    *x21*x31    *x51
      p_m12_0_5   =p_m12_0_5   
     9  +coeff( 63)    *x22*x33        
     1  +coeff( 64)    *x22*x31*x42    
     2  +coeff( 65)    *x21*x32*x41*x51
     3  +coeff( 66)    *x23*x32*x41    
     4  +coeff( 67)    *x22*x33    *x51
     5  +coeff( 68)        *x33*x42*x51
     6  +coeff( 69)    *x22    *x43*x51
     7  +coeff( 70)    *x21*x32*x41*x52
     8  +coeff( 71)        *x32*x41*x53
      p_m12_0_5   =p_m12_0_5   
     9  +coeff( 72)        *x31*x42*x53
     1  +coeff( 73)*x11        *x41*x53
     2  +coeff( 74)    *x21*x31*x42*x53
     3  +coeff( 75)    *x21    *x43*x53
     4  +coeff( 76)    *x23*x32*x41*x53
     5  +coeff( 77)*x11    *x31*x42*x53
     6  +coeff( 78)    *x23*x31*x42*x53
     7  +coeff( 79)    *x23    *x43*x53
     8  +coeff( 80)*x11*x21    *x43*x53
      p_m12_0_5   =p_m12_0_5   
     9  +coeff( 81)*x11*x22*x32*x41*x53
     1  +coeff( 82)*x11        *x41    
     2  +coeff( 83)    *x22*x31    *x51
     3  +coeff( 84)    *x22*x32*x41    
     4  +coeff( 85)*x11        *x41*x51
     5  +coeff( 86)    *x22*x31    *x52
     6  +coeff( 87)        *x31*x42*x52
     7  +coeff( 88)            *x43*x52
     8  +coeff( 89)*x11*x21*x31    *x51
      p_m12_0_5   =p_m12_0_5   
     9  +coeff( 90)        *x32*x43*x51
c
      return
      end
      real function l_m12_0_5   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)   !,xmean(10)
      dimension coeff( 41)
      data ncoeff/ 40/
      data avdat/  0.5128684E-02/
      data xmin/
     1 -0.49998E-02,-0.79858E-01,-0.99988E-01,-0.37028E-01,-0.42461E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49981E-02, 0.64228E-01, 0.99988E-01, 0.37028E-01, 0.42498E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.12806107E-01,-0.16379207E+00,-0.35672721E-01,-0.26426058E-01,
     + -0.28135156E-01,-0.24550335E-01,-0.95198387E-02, 0.97866161E-02,
     +  0.10251816E-01, 0.15322209E-01,-0.23384858E-02,-0.10444614E-02,
     +  0.81954077E-02, 0.32661634E-03,-0.97811373E-03,-0.13564615E-01,
     +  0.51126159E-02, 0.61737485E-02, 0.12562558E-01, 0.22588300E-04,
     + -0.33063459E-03, 0.14098983E-02,-0.95824404E-02,-0.76684519E-02,
     +  0.76863845E-02,-0.37711873E-02, 0.76652667E-02,-0.26355282E-03,
     + -0.39049205E-02, 0.31034759E-03, 0.16505737E-02,-0.10860686E-02,
     + -0.16552921E-02, 0.42898972E-02, 0.10669131E-02, 0.17605163E-03,
     + -0.94263401E-03, 0.78459270E-03, 0.12338293E-02, 0.16148799E-02,
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
      l_m12_0_5   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22            
     5  +coeff(  5)            *x42    
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)                *x52
      l_m12_0_5   =l_m12_0_5   
     9  +coeff(  9)    *x22        *x51
     1  +coeff( 10)            *x42*x51
     2  +coeff( 11)    *x24        *x54
     3  +coeff( 12)*x11                
     4  +coeff( 13)        *x31*x41*x51
     5  +coeff( 14)    *x21        *x52
     6  +coeff( 15)        *x32        
     7  +coeff( 16)            *x42*x52
     8  +coeff( 17)    *x21        *x53
      l_m12_0_5   =l_m12_0_5   
     9  +coeff( 18)    *x22        *x53
     1  +coeff( 19)    *x24    *x42    
     2  +coeff( 20)        *x33*x41*x52
     3  +coeff( 21)    *x24*x33*x41    
     4  +coeff( 22)        *x32    *x51
     5  +coeff( 23)    *x22        *x52
     6  +coeff( 24)        *x31*x41*x52
     7  +coeff( 25)            *x42*x53
     8  +coeff( 26)    *x21        *x54
      l_m12_0_5   =l_m12_0_5   
     9  +coeff( 27)    *x24*x31*x41    
     1  +coeff( 28)*x11*x21            
     2  +coeff( 29)    *x24            
     3  +coeff( 30)*x11*x21        *x51
     4  +coeff( 31)    *x21    *x42*x51
     5  +coeff( 32)        *x32    *x52
     6  +coeff( 33)                *x54
     7  +coeff( 34)        *x31*x41*x53
     8  +coeff( 35)    *x24*x32        
      l_m12_0_5   =l_m12_0_5   
     9  +coeff( 36)    *x21    *x42    
     1  +coeff( 37)                *x53
     2  +coeff( 38)    *x21*x31*x41*x51
     3  +coeff( 39)    *x23*x31*x41    
     4  +coeff( 40)    *x23    *x42    
c
      return
      end
      function x_m12_0_7   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)   !,xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.1070197E+00/
      data xmin/
     1 -0.49998E-02,-0.79858E-01,-0.99988E-01,-0.35976E-01,-0.42461E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49981E-02, 0.64054E-01, 0.99988E-01, 0.35976E-01, 0.42498E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.78712434E-01, 0.26452401E+00, 0.35841990E+00, 0.89025348E-02,
     + -0.58338256E-03,-0.32828841E-02,-0.63162018E-02, 0.29126796E+00,
     + -0.16513805E+00,-0.68775757E-03,-0.85048461E-02,-0.33389758E-02,
     + -0.16434098E-01,-0.31627797E-01, 0.17186308E-02, 0.36661953E-03,
     +  0.24835046E-02, 0.83181309E-02,-0.12859783E+00, 0.68229824E-01,
     +  0.86650572E-04,-0.71367003E-01,-0.12237111E+00, 0.19033302E-02,
     + -0.36778485E-02,-0.35805060E-02,-0.32837652E-01,-0.23792011E-02,
     +  0.79040881E-02, 0.51183883E-01, 0.19374426E-02, 0.43135315E-01,
     +  0.33648741E-01,-0.15268828E-02, 0.16118194E-02, 0.39508052E-01,
     +  0.25610490E-01,-0.17028505E-01, 0.78872601E-02, 0.19352260E-03,
     +  0.27243715E-01,-0.28046297E-02, 0.36836429E-02, 0.62651313E-02,
     +  0.25599017E-02, 0.22851897E-01, 0.20030490E-02, 0.17691730E-02,
     + -0.11427308E-01,-0.26707320E-01, 0.58839861E-02,-0.86778710E-02,
     +  0.15944945E-01,-0.10684007E-01,-0.25634008E-01,-0.29540100E-03,
     +  0.15532142E-01, 0.40580649E-01,-0.88547971E-02,-0.11935204E-03,
     + -0.11439207E-02,-0.22936461E-02,-0.20705671E-02, 0.14260432E-02,
     +  0.68625221E-02,-0.71907835E-03, 0.14980695E-01, 0.94871093E-02,
     + -0.32419838E-01, 0.42886022E-02, 0.68636646E-03,-0.19430041E-01,
     + -0.29099526E-01,-0.73972745E-02, 0.13647407E-01,-0.15911579E-01,
     + -0.26479231E-01, 0.40938020E-01,-0.95952470E-02,-0.24986200E-02,
     +  0.11234570E-01,-0.44894028E-02,-0.11266388E-01,-0.61362344E-02,
     + -0.42815786E-01, 0.42415978E-02, 0.35407525E-01, 0.88617140E-02,
     + -0.13269567E-02, 0.63326862E-02,
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
      x_m12_0_7   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22            
     5  +coeff(  5)        *x32        
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21        *x51
      x_m12_0_7   =x_m12_0_7   
     9  +coeff(  9)                *x52
     1  +coeff( 10)*x11                
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x21*x32        
     4  +coeff( 13)    *x21*x31*x41    
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)    *x22        *x51
     7  +coeff( 16)        *x32    *x51
     8  +coeff( 17)        *x31*x41*x51
      x_m12_0_7   =x_m12_0_7   
     9  +coeff( 18)            *x42*x51
     1  +coeff( 19)    *x21        *x52
     2  +coeff( 20)                *x53
     3  +coeff( 21)*x11*x21            
     4  +coeff( 22)    *x22*x31*x41    
     5  +coeff( 23)    *x22    *x42    
     6  +coeff( 24)*x11            *x51
     7  +coeff( 25)    *x23        *x51
     8  +coeff( 26)    *x21*x32    *x51
      x_m12_0_7   =x_m12_0_7   
     9  +coeff( 27)    *x21    *x42*x51
     1  +coeff( 28)    *x22        *x52
     2  +coeff( 29)            *x42*x52
     3  +coeff( 30)    *x21        *x53
     4  +coeff( 31)    *x21*x33*x41    
     5  +coeff( 32)    *x22*x31*x41*x51
     6  +coeff( 33)    *x22    *x42*x51
     7  +coeff( 34)*x11            *x52
     8  +coeff( 35)    *x23        *x52
      x_m12_0_7   =x_m12_0_7   
     9  +coeff( 36)    *x21    *x42*x52
     1  +coeff( 37)    *x22        *x53
     2  +coeff( 38)            *x42*x53
     3  +coeff( 39)    *x22*x32*x42    
     4  +coeff( 40)        *x33*x43    
     5  +coeff( 41)    *x23*x31*x41*x51
     6  +coeff( 42)    *x22*x32    *x52
     7  +coeff( 43)    *x22*x31*x41*x52
     8  +coeff( 44)        *x32*x42*x52
      x_m12_0_7   =x_m12_0_7   
     9  +coeff( 45)        *x31*x43*x52
     1  +coeff( 46)    *x23        *x53
     2  +coeff( 47)    *x21*x31*x41*x53
     3  +coeff( 48)    *x23*x32    *x52
     4  +coeff( 49)    *x23*x31*x41*x52
     5  +coeff( 50)    *x21*x31*x43*x52
     6  +coeff( 51)    *x22*x32    *x53
     7  +coeff( 52)    *x22*x31*x41*x53
     8  +coeff( 53)    *x21*x33*x41*x53
      x_m12_0_7   =x_m12_0_7   
     9  +coeff( 54)    *x22*x32        
     1  +coeff( 55)    *x21*x31*x41*x51
     2  +coeff( 56)    *x23*x31*x41    
     3  +coeff( 57)    *x21*x32*x42    
     4  +coeff( 58)    *x21*x31*x41*x52
     5  +coeff( 59)        *x31*x41*x53
     6  +coeff( 60)*x11*x21    *x42    
     7  +coeff( 61)*x11    *x32    *x51
     8  +coeff( 62)    *x23*x32    *x51
      x_m12_0_7   =x_m12_0_7   
     9  +coeff( 63)    *x23    *x42*x51
     1  +coeff( 64)*x11            *x53
     2  +coeff( 65)    *x21*x32    *x53
     3  +coeff( 66)*x11*x22    *x42    
     4  +coeff( 67)    *x21*x33*x43    
     5  +coeff( 68)    *x22*x32*x42*x51
     6  +coeff( 69)    *x22*x31*x43*x51
     7  +coeff( 70)        *x33*x43*x51
     8  +coeff( 71)*x11    *x32    *x52
      x_m12_0_7   =x_m12_0_7   
     9  +coeff( 72)    *x21*x33*x41*x52
     1  +coeff( 73)    *x23    *x42*x52
     2  +coeff( 74)*x11*x23    *x42    
     3  +coeff( 75)    *x22*x33*x43    
     4  +coeff( 76)    *x22*x33*x41*x52
     5  +coeff( 77)    *x23*x31*x41*x53
     6  +coeff( 78)    *x23    *x42*x53
     7  +coeff( 79)    *x21*x32*x42*x53
     8  +coeff( 80)*x11*x23*x32    *x51
      x_m12_0_7   =x_m12_0_7   
     9  +coeff( 81)*x11*x23    *x42*x51
     1  +coeff( 82)*x11*x22*x31*x41*x52
     2  +coeff( 83)*x11*x23*x31*x43    
     3  +coeff( 84)*x11    *x33*x43*x51
     4  +coeff( 85)    *x23*x33*x43*x51
     5  +coeff( 86)*x11    *x33*x41*x53
     6  +coeff( 87)    *x23*x31*x43*x53
     7  +coeff( 88)*x11*x23*x31*x41*x53
     8  +coeff( 89)        *x31*x43    
      x_m12_0_7   =x_m12_0_7   
     9  +coeff( 90)        *x31*x41*x52
c
      return
      end
      function t_m12_0_7   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)   !,xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.1761815E-01/
      data xmin/
     1 -0.49998E-02,-0.79858E-01,-0.99988E-01,-0.35976E-01,-0.42461E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49981E-02, 0.64054E-01, 0.99988E-01, 0.35976E-01, 0.42498E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.25654830E-01,-0.54169219E-01, 0.15951107E+00, 0.25297569E-02,
     + -0.79858984E-03,-0.53495285E-02,-0.10239994E-01, 0.72846428E-01,
     + -0.10273220E+00,-0.96742337E-03,-0.29308919E-02, 0.49939874E-04,
     + -0.27104076E-02,-0.10398205E-01,-0.29931080E-02, 0.15346035E-02,
     +  0.15191221E-01,-0.60735468E-01, 0.54770298E-01, 0.35105966E-03,
     + -0.32871019E-02,-0.25806701E-01, 0.15205295E-03,-0.44228625E-01,
     +  0.10781111E-02,-0.17446266E-02, 0.32203583E-03,-0.31483823E-02,
     + -0.12436938E-01,-0.17777335E-01, 0.50920718E-02, 0.11689775E-01,
     +  0.40157437E-01,-0.12269196E-02,-0.31710090E-02, 0.17552413E-01,
     +  0.24753367E-02, 0.24673095E-01,-0.11713512E-02, 0.11276123E-01,
     +  0.22040935E-01, 0.41887891E-01, 0.34463651E-01,-0.21645389E-01,
     + -0.54432130E-02, 0.23164342E-02,-0.34774288E-02, 0.12516381E-02,
     + -0.11483261E-02, 0.95569179E-03, 0.27373577E-01,-0.46285517E-02,
     +  0.78590102E-02, 0.11716039E-01,-0.14717086E-02,-0.54603204E-03,
     + -0.28536032E-03, 0.49711880E-03, 0.12607625E-02,-0.14611368E-01,
     +  0.46934004E-03,-0.63623823E-02, 0.20080125E-02,-0.10489206E-02,
     + -0.43099359E-03, 0.41807466E-02, 0.20950676E-02, 0.89289090E-02,
     + -0.10772302E-01,-0.79895984E-02, 0.76229539E-03, 0.50409068E-04,
     + -0.11693992E-02, 0.12169700E-02,-0.13674144E-02,-0.56532223E-03,
     + -0.69313450E-03, 0.14225985E-02,-0.46268763E-03,-0.77642649E-02,
     + -0.64379289E-02, 0.85167045E-03,-0.30993912E-03,-0.58564445E-03,
     + -0.13801988E-02, 0.11069872E-03,-0.93083654E-03,-0.10748178E-02,
     + -0.59205331E-02,-0.44699237E-02,
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
      t_m12_0_7   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22            
     5  +coeff(  5)        *x32        
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21        *x51
      t_m12_0_7   =t_m12_0_7   
     9  +coeff(  9)                *x52
     1  +coeff( 10)*x11                
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x21*x32        
     4  +coeff( 13)    *x21*x31*x41    
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)    *x22        *x51
     7  +coeff( 16)        *x32    *x51
     8  +coeff( 17)            *x42*x51
      t_m12_0_7   =t_m12_0_7   
     9  +coeff( 18)    *x21        *x52
     1  +coeff( 19)                *x53
     2  +coeff( 20)*x11*x21            
     3  +coeff( 21)    *x22*x32        
     4  +coeff( 22)    *x22*x31*x41    
     5  +coeff( 23)        *x33*x41    
     6  +coeff( 24)    *x22    *x42    
     7  +coeff( 25)        *x32*x42    
     8  +coeff( 26)        *x31*x43    
      t_m12_0_7   =t_m12_0_7   
     9  +coeff( 27)*x11            *x51
     1  +coeff( 28)    *x23        *x51
     2  +coeff( 29)    *x21*x31*x41*x51
     3  +coeff( 30)    *x21    *x42*x51
     4  +coeff( 31)    *x22        *x52
     5  +coeff( 32)            *x42*x52
     6  +coeff( 33)    *x21        *x53
     7  +coeff( 34)    *x23*x31*x41    
     8  +coeff( 35)    *x21*x33*x41    
      t_m12_0_7   =t_m12_0_7   
     9  +coeff( 36)    *x22*x31*x41*x51
     1  +coeff( 37)        *x33*x41*x51
     2  +coeff( 38)    *x22    *x42*x51
     3  +coeff( 39)*x11            *x52
     4  +coeff( 40)    *x23        *x52
     5  +coeff( 41)    *x21*x31*x41*x52
     6  +coeff( 42)    *x21    *x42*x52
     7  +coeff( 43)    *x22        *x53
     8  +coeff( 44)            *x42*x53
      t_m12_0_7   =t_m12_0_7   
     9  +coeff( 45)    *x23*x32    *x51
     1  +coeff( 46)    *x21*x33*x41*x51
     2  +coeff( 47)    *x22*x32    *x52
     3  +coeff( 48)    *x22    *x42*x52
     4  +coeff( 49)        *x31*x43*x52
     5  +coeff( 50)*x11            *x53
     6  +coeff( 51)    *x23        *x53
     7  +coeff( 52)    *x21*x31*x41*x53
     8  +coeff( 53)        *x31*x41*x51
      t_m12_0_7   =t_m12_0_7   
     9  +coeff( 54)        *x31*x41*x52
     1  +coeff( 55)    *x23*x32        
     2  +coeff( 56)*x11*x21        *x51
     3  +coeff( 57)    *x22*x32    *x51
     4  +coeff( 58)        *x31*x43*x51
     5  +coeff( 59)    *x21*x32    *x52
     6  +coeff( 60)        *x31*x41*x53
     7  +coeff( 61)*x11*x21*x31*x41    
     8  +coeff( 62)    *x22*x32*x42    
      t_m12_0_7   =t_m12_0_7   
     9  +coeff( 63)        *x33*x43    
     1  +coeff( 64)*x11*x22        *x51
     2  +coeff( 65)*x11    *x32    *x51
     3  +coeff( 66)    *x23*x31*x41*x51
     4  +coeff( 67)*x11        *x42*x51
     5  +coeff( 68)    *x23    *x42*x51
     6  +coeff( 69)    *x21*x32*x42*x51
     7  +coeff( 70)        *x33*x41*x52
     8  +coeff( 71)    *x21*x32    *x51
      t_m12_0_7   =t_m12_0_7   
     9  +coeff( 72)*x11    *x32        
     1  +coeff( 73)*x11        *x42    
     2  +coeff( 74)        *x32*x42*x51
     3  +coeff( 75)        *x32    *x53
     4  +coeff( 76)*x11*x23            
     5  +coeff( 77)*x11*x21    *x42    
     6  +coeff( 78)*x11    *x31*x41*x51
     7  +coeff( 79)*x11*x21        *x52
     8  +coeff( 80)    *x22*x31*x41*x52
      t_m12_0_7   =t_m12_0_7   
     9  +coeff( 81)    *x21    *x42*x53
     1  +coeff( 82)        *x32    *x52
     2  +coeff( 83)*x11*x22            
     3  +coeff( 84)*x11    *x31*x41    
     4  +coeff( 85)    *x23    *x42    
     5  +coeff( 86)*x11*x21*x32        
     6  +coeff( 87)    *x22*x33*x41    
     7  +coeff( 88)    *x22*x31*x43    
     8  +coeff( 89)    *x21*x31*x43*x51
      t_m12_0_7   =t_m12_0_7   
     9  +coeff( 90)        *x32*x42*x52
c
      return
      end
      function y_m12_0_7   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)   !,xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49998E-02,-0.79858E-01,-0.99988E-01,-0.35976E-01,-0.42461E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49981E-02, 0.64054E-01, 0.99988E-01, 0.35976E-01, 0.42498E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.22026159E+00, 0.84122157E+00,-0.47919094E-02, 0.35989655E-02,
     + -0.52449413E-01,-0.15851477E+00,-0.10671340E-01,-0.23174302E-03,
     + -0.30686693E-01,-0.16345944E-02,-0.34721574E-03, 0.20894832E-02,
     +  0.31141490E-01, 0.93917139E-01,-0.35104766E-01,-0.12000224E+00,
     +  0.10151427E-01, 0.72764587E-02,-0.35853020E-03,-0.34861201E-02,
     +  0.69086752E-02,-0.16790938E-02,-0.65530660E-02,-0.48572827E-01,
     +  0.11647854E-01, 0.42599503E-01,-0.46211787E-03, 0.36949903E-01,
     +  0.13636750E+00, 0.16050577E+00, 0.13358379E-01,-0.26995773E-03,
     +  0.29881818E-01,-0.24610420E-03, 0.16146125E-02, 0.32580143E-03,
     +  0.33899709E-02,-0.10416504E-02, 0.28347294E-02,-0.82262431E-03,
     +  0.69130398E-03, 0.58943806E-02,-0.16581455E-01, 0.32005813E-01,
     +  0.30077903E-01,-0.55886987E-02,-0.40511172E-02,-0.11359112E-01,
     + -0.54206207E-01, 0.13038474E-02,-0.71513760E-02,-0.13060689E-01,
     +  0.63488379E-01, 0.31698331E-01,-0.17883059E-02,-0.71946373E-02,
     + -0.40313799E-01,-0.42300666E-03,-0.81081346E-01,-0.70402148E-03,
     + -0.62077194E-02, 0.77723133E-04, 0.50581708E-02, 0.13229229E-01,
     +  0.14011345E-02, 0.32791682E-02,-0.28800806E-02,-0.37015995E-03,
     + -0.95926011E-02,-0.47143426E-01,-0.86444098E-03,-0.16464436E-02,
     +  0.10688042E-01, 0.14065020E-01,-0.10090588E-01, 0.13470232E-01,
     +  0.52333593E-01, 0.12986372E-02, 0.49073324E-01,-0.14556461E-01,
     + -0.78807950E-01,-0.30752982E-02,-0.44705991E-01, 0.86400900E-02,
     + -0.84205903E-02,-0.18949510E-03, 0.10810330E-02, 0.65452829E-02,
     + -0.42621186E-02, 0.35490326E-02,
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
      y_m12_0_7   =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21*x31        
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)        *x31    *x51
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)    *x22*x31        
     8  +coeff(  8)        *x33        
      y_m12_0_7   =y_m12_0_7   
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)        *x32*x41    
     2  +coeff( 11)            *x43    
     3  +coeff( 12)    *x21*x31    *x51
     4  +coeff( 13)        *x31    *x52
     5  +coeff( 14)            *x41*x52
     6  +coeff( 15)    *x23*x31        
     7  +coeff( 16)    *x23    *x41    
     8  +coeff( 17)    *x21*x31*x42    
      y_m12_0_7   =y_m12_0_7   
     9  +coeff( 18)    *x21    *x43    
     1  +coeff( 19)    *x22*x31    *x51
     2  +coeff( 20)    *x22    *x41*x51
     3  +coeff( 21)            *x43*x51
     4  +coeff( 22)    *x21*x31    *x52
     5  +coeff( 23)    *x21    *x41*x52
     6  +coeff( 24)            *x41*x53
     7  +coeff( 25)    *x23*x31    *x51
     8  +coeff( 26)    *x23    *x41*x51
      y_m12_0_7   =y_m12_0_7   
     9  +coeff( 27)    *x21*x31    *x53
     1  +coeff( 28)    *x23*x32*x41    
     2  +coeff( 29)    *x23*x31*x42    
     3  +coeff( 30)    *x23    *x43    
     4  +coeff( 31)    *x22*x31*x42*x51
     5  +coeff( 32)        *x33*x42*x51
     6  +coeff( 33)    *x22    *x43*x51
     7  +coeff( 34)    *x22*x31    *x53
     8  +coeff( 35)        *x33    *x53
      y_m12_0_7   =y_m12_0_7   
     9  +coeff( 36)    *x21*x33    *x53
     1  +coeff( 37)    *x22*x32*x43*x51
     2  +coeff( 38)    *x22*x33    *x53
     3  +coeff( 39)    *x21*x33        
     4  +coeff( 40)        *x33    *x51
     5  +coeff( 41)        *x32*x41*x51
     6  +coeff( 42)        *x31*x42*x51
     7  +coeff( 43)        *x31    *x53
     8  +coeff( 44)    *x22*x31*x42    
      y_m12_0_7   =y_m12_0_7   
     9  +coeff( 45)    *x22    *x43    
     1  +coeff( 46)    *x21*x32*x41*x51
     2  +coeff( 47)    *x21*x31*x42*x51
     3  +coeff( 48)    *x22    *x41*x52
     4  +coeff( 49)            *x43*x52
     5  +coeff( 50)    *x21*x32*x43    
     6  +coeff( 51)    *x23*x31    *x52
     7  +coeff( 52)    *x23    *x41*x52
     8  +coeff( 53)            *x43*x53
      y_m12_0_7   =y_m12_0_7   
     9  +coeff( 54)    *x22*x31*x42*x52
     1  +coeff( 55)        *x33*x42*x52
     2  +coeff( 56)        *x32*x43*x52
     3  +coeff( 57)    *x22*x31*x42*x53
     4  +coeff( 58)        *x33*x42*x53
     5  +coeff( 59)    *x23    *x43*x53
     6  +coeff( 60)        *x31*x42    
     7  +coeff( 61)    *x21    *x41*x51
     8  +coeff( 62)*x11        *x41    
      y_m12_0_7   =y_m12_0_7   
     9  +coeff( 63)    *x21*x32*x41    
     1  +coeff( 64)    *x22*x32*x41    
     2  +coeff( 65)        *x33*x42    
     3  +coeff( 66)        *x32*x43    
     4  +coeff( 67)    *x21*x33    *x51
     5  +coeff( 68)        *x33    *x52
     6  +coeff( 69)        *x32*x41*x52
     7  +coeff( 70)        *x31*x42*x52
     8  +coeff( 71)*x11*x22*x31        
      y_m12_0_7   =y_m12_0_7   
     9  +coeff( 72)*x11*x22    *x41    
     1  +coeff( 73)    *x21*x31*x42*x52
     2  +coeff( 74)    *x21    *x43*x52
     3  +coeff( 75)    *x22    *x41*x53
     4  +coeff( 76)        *x32*x41*x53
     5  +coeff( 77)        *x31*x42*x53
     6  +coeff( 78)    *x23*x31*x42*x51
     7  +coeff( 79)    *x22    *x43*x52
     8  +coeff( 80)    *x23    *x41*x53
      y_m12_0_7   =y_m12_0_7   
     9  +coeff( 81)    *x22    *x43*x53
     1  +coeff( 82)*x11*x22*x31*x42*x51
     2  +coeff( 83)    *x23*x31*x42*x53
     3  +coeff( 84)*x11*x22*x31*x42*x52
     4  +coeff( 85)*x11*x23    *x43*x52
     5  +coeff( 86)*x11*x21*x31        
     6  +coeff( 87)    *x22*x33        
     7  +coeff( 88)    *x21    *x43*x51
     8  +coeff( 89)    *x22*x31    *x52
      y_m12_0_7   =y_m12_0_7   
     9  +coeff( 90)    *x21    *x41*x53
c
      return
      end
      function p_m12_0_7   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)   !,xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49998E-02,-0.79858E-01,-0.99988E-01,-0.35976E-01,-0.42461E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49981E-02, 0.64054E-01, 0.99988E-01, 0.35976E-01, 0.42498E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.45115035E-02, 0.34028720E-01, 0.19086774E-03, 0.20900774E-02,
     +  0.42399717E-02, 0.14960937E-01,-0.17664677E-02, 0.47151483E-04,
     + -0.74487920E-02,-0.13797578E-02,-0.44083879E-02,-0.68171383E-02,
     + -0.16369181E-02,-0.69754049E-02,-0.50212792E-02,-0.18460186E-01,
     +  0.17097336E-04,-0.54625557E-02, 0.64838707E-03, 0.23216417E-04,
     + -0.19554479E-01, 0.69079562E-02, 0.15710428E-01, 0.13048436E-01,
     + -0.39893724E-02,-0.22466884E-02,-0.14351311E-01, 0.14168646E-01,
     +  0.19120628E-01,-0.13771061E-03,-0.23128230E-02, 0.10144055E-01,
     + -0.51261182E-03,-0.10251141E-02,-0.12908660E-04,-0.76474986E-04,
     +  0.13300202E-03,-0.69870852E-03,-0.19229319E-02,-0.10601781E-01,
     +  0.32475829E-03, 0.96703442E-02, 0.81165014E-02,-0.16367751E-02,
     + -0.87454682E-02,-0.79816621E-03,-0.34605355E-02,-0.52868323E-02,
     + -0.25355380E-02,-0.24441052E-02, 0.68905836E-04, 0.88647001E-04,
     +  0.21053109E-01,-0.29146874E-02, 0.35128228E-01,-0.17242940E-02,
     +  0.50152151E-03, 0.20728777E-03, 0.66854972E-02, 0.29282829E-01,
     +  0.17326578E-02, 0.47359444E-01, 0.33842903E-02,-0.13338677E-03,
     +  0.28286502E-03, 0.42771711E-03,-0.24733873E-03, 0.14185940E-02,
     +  0.13439409E-01, 0.28664229E-01, 0.17382260E-02, 0.46488857E-02,
     +  0.80505881E-04, 0.58763702E-02, 0.65605552E-03,-0.29102778E-02,
     + -0.40958214E-03, 0.46391888E-02,-0.54086151E-03, 0.85102338E-02,
     + -0.17437787E-03, 0.24935193E-01, 0.19198752E-02, 0.21422810E-03,
     + -0.42381450E-02, 0.12809684E-01,-0.66539636E-02, 0.40646296E-01,
     + -0.92079220E-02, 0.38462116E-02,
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
      p_m12_0_7   =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21*x31        
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)        *x31    *x51
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)    *x22*x31        
     8  +coeff(  8)        *x33        
      p_m12_0_7   =p_m12_0_7   
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)        *x32*x41    
     2  +coeff( 11)        *x31*x42    
     3  +coeff( 12)            *x43    
     4  +coeff( 13)    *x21*x31    *x51
     5  +coeff( 14)    *x21    *x41*x51
     6  +coeff( 15)        *x31    *x52
     7  +coeff( 16)            *x41*x52
     8  +coeff( 17)*x11    *x31        
      p_m12_0_7   =p_m12_0_7   
     9  +coeff( 18)    *x23*x31        
     1  +coeff( 19)    *x21*x33        
     2  +coeff( 20)*x11        *x41    
     3  +coeff( 21)    *x23    *x41    
     4  +coeff( 22)    *x21*x32*x41    
     5  +coeff( 23)    *x21*x31*x42    
     6  +coeff( 24)    *x21    *x43    
     7  +coeff( 25)    *x22*x31    *x51
     8  +coeff( 26)        *x33    *x51
      p_m12_0_7   =p_m12_0_7   
     9  +coeff( 27)    *x22    *x41*x51
     1  +coeff( 28)        *x31*x42*x51
     2  +coeff( 29)            *x43*x51
     3  +coeff( 30)    *x21*x31    *x52
     4  +coeff( 31)    *x21    *x41*x52
     5  +coeff( 32)            *x41*x53
     6  +coeff( 33)    *x22*x33        
     7  +coeff( 34)        *x33*x42    
     8  +coeff( 35)    *x22    *x43    
      p_m12_0_7   =p_m12_0_7   
     9  +coeff( 36)        *x32*x43    
     1  +coeff( 37)*x11    *x31    *x51
     2  +coeff( 38)    *x23*x31    *x51
     3  +coeff( 39)    *x21*x33    *x51
     4  +coeff( 40)    *x23    *x41*x51
     5  +coeff( 41)    *x21*x32*x41*x51
     6  +coeff( 42)    *x21*x31*x42*x51
     7  +coeff( 43)    *x21    *x43*x51
     8  +coeff( 44)    *x22*x31    *x52
      p_m12_0_7   =p_m12_0_7   
     9  +coeff( 45)    *x22    *x41*x52
     1  +coeff( 46)        *x32*x41*x52
     2  +coeff( 47)        *x31*x42*x52
     3  +coeff( 48)            *x43*x52
     4  +coeff( 49)    *x21*x31    *x53
     5  +coeff( 50)    *x21    *x41*x53
     6  +coeff( 51)*x11*x22*x31        
     7  +coeff( 52)*x11    *x32*x41    
     8  +coeff( 53)    *x23*x31*x42    
      p_m12_0_7   =p_m12_0_7   
     9  +coeff( 54)    *x21*x33*x42    
     1  +coeff( 55)    *x23    *x43    
     2  +coeff( 56)    *x21*x32*x43    
     3  +coeff( 57)    *x22*x33    *x51
     4  +coeff( 58)*x11*x21    *x41*x51
     5  +coeff( 59)    *x22*x32*x41*x51
     6  +coeff( 60)    *x22*x31*x42*x51
     7  +coeff( 61)        *x33*x42*x51
     8  +coeff( 62)    *x22    *x43*x51
      p_m12_0_7   =p_m12_0_7   
     9  +coeff( 63)        *x32*x43*x51
     1  +coeff( 64)*x11    *x31    *x52
     2  +coeff( 65)    *x23*x31    *x52
     3  +coeff( 66)    *x21*x33    *x52
     4  +coeff( 67)*x11        *x41*x52
     5  +coeff( 68)    *x23    *x41*x52
     6  +coeff( 69)    *x21*x31*x42*x52
     7  +coeff( 70)    *x21    *x43*x52
     8  +coeff( 71)    *x22*x31    *x53
      p_m12_0_7   =p_m12_0_7   
     9  +coeff( 72)        *x33    *x53
     1  +coeff( 73)    *x22    *x41*x53
     2  +coeff( 74)        *x32*x41*x53
     3  +coeff( 75)        *x31*x42*x53
     4  +coeff( 76)            *x43*x53
     5  +coeff( 77)*x11*x21*x31*x42    
     6  +coeff( 78)    *x22*x33*x42    
     7  +coeff( 79)*x11*x21    *x43    
     8  +coeff( 80)    *x22*x32*x43    
      p_m12_0_7   =p_m12_0_7   
     9  +coeff( 81)*x11*x22*x31    *x51
     1  +coeff( 82)    *x23    *x43*x51
     2  +coeff( 83)    *x21*x32*x43*x51
     3  +coeff( 84)*x11*x21*x31    *x52
     4  +coeff( 85)    *x22*x32*x41*x52
     5  +coeff( 86)    *x22*x31*x42*x52
     6  +coeff( 87)        *x33*x42*x52
     7  +coeff( 88)    *x22    *x43*x52
     8  +coeff( 89)        *x32*x43*x52
      p_m12_0_7   =p_m12_0_7   
     9  +coeff( 90)    *x21*x33    *x53
c
      return
      end
      real function l_m12_0_7   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)   !,xmean(10)
      dimension coeff( 41)
      data ncoeff/ 40/
      data avdat/  0.1548579E-01/
      data xmin/
     1 -0.49998E-02,-0.79858E-01,-0.99988E-01,-0.35976E-01,-0.42461E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49981E-02, 0.64054E-01, 0.99988E-01, 0.35976E-01, 0.42498E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.12553940E-02,-0.14958182E+00,-0.51714685E-01,-0.31805634E-01,
     + -0.34831684E-01, 0.18825127E-01,-0.12143780E-01,-0.15589804E-01,
     +  0.22689413E-01,-0.12824934E-01,-0.21912191E-01, 0.11026592E-01,
     + -0.87730965E-03, 0.58026137E-02, 0.22312367E-01,-0.93749622E-02,
     +  0.62547141E-03,-0.20085075E-02, 0.12727027E-01,-0.15363492E-02,
     +  0.14064342E-01,-0.53930520E-02,-0.18415796E-01,-0.17561488E-01,
     +  0.12211588E-02,-0.46606272E-03, 0.19678720E-02, 0.87935460E-03,
     + -0.11526139E-01, 0.95478278E-02, 0.14301832E-01,-0.96614640E-02,
     +  0.13529055E-02, 0.12498542E-02, 0.60501724E-03,-0.14388530E-02,
     +  0.12457328E-02, 0.32409088E-03, 0.65704589E-02, 0.44839741E-02,
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
      l_m12_0_7   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22            
     5  +coeff(  5)            *x42    
     6  +coeff(  6)    *x22        *x51
     7  +coeff(  7)    *x21        *x52
     8  +coeff(  8)    *x22        *x52
      l_m12_0_7   =l_m12_0_7   
     9  +coeff(  9)    *x21        *x53
     1  +coeff( 10)        *x31*x41    
     2  +coeff( 11)    *x21        *x51
     3  +coeff( 12)                *x53
     4  +coeff( 13)*x11                
     5  +coeff( 14)                *x52
     6  +coeff( 15)            *x42*x51
     7  +coeff( 16)                *x54
     8  +coeff( 17)    *x22*x31*x41*x51
      l_m12_0_7   =l_m12_0_7   
     9  +coeff( 18)        *x33*x41*x51
     1  +coeff( 19)    *x22        *x53
     2  +coeff( 20)        *x32        
     3  +coeff( 21)        *x31*x41*x51
     4  +coeff( 22)    *x24            
     5  +coeff( 23)            *x42*x52
     6  +coeff( 24)    *x21        *x54
     7  +coeff( 25)        *x31*x41*x54
     8  +coeff( 26)*x11*x21            
      l_m12_0_7   =l_m12_0_7   
     9  +coeff( 27)        *x32    *x51
     1  +coeff( 28)    *x22    *x42    
     2  +coeff( 29)        *x31*x41*x52
     3  +coeff( 30)    *x24*x31*x41    
     4  +coeff( 31)    *x24    *x42    
     5  +coeff( 32)    *x22        *x54
     6  +coeff( 33)    *x21    *x42    
     7  +coeff( 34)    *x22*x32        
     8  +coeff( 35)*x11*x21        *x51
      l_m12_0_7   =l_m12_0_7   
     9  +coeff( 36)        *x32    *x52
     1  +coeff( 37)    *x23*x31*x41    
     2  +coeff( 38)*x11            *x53
     3  +coeff( 39)            *x42*x53
     4  +coeff( 40)        *x33*x41*x53
c
      return
      end
      function x_m12_0_9   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)   !,xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.3910198E-01/
      data xmin/
     1 -0.49995E-02,-0.79858E-01,-0.99988E-01,-0.35144E-01,-0.25572E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49981E-02, 0.61732E-01, 0.99988E-01, 0.35144E-01, 0.42498E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.87848544E-01, 0.26952633E+00, 0.56594253E+00, 0.22386258E-01,
     + -0.23259779E-02,-0.16194554E-01,-0.30843070E-01, 0.32342881E+00,
     + -0.22024716E+00,-0.29228579E-02,-0.46840301E-02,-0.55911973E-01,
     +  0.10959051E-02,-0.26156567E-03, 0.46980273E-03, 0.38990262E-03,
     + -0.10339246E+00, 0.78624822E-01, 0.23006706E-03,-0.13240907E+00,
     + -0.22528540E+00, 0.27139350E-02,-0.65026805E-02,-0.53269684E-01,
     + -0.94717935E-01, 0.11845183E-02,-0.88480562E-02,-0.21877168E-01,
     +  0.39446898E-01, 0.36227496E-03,-0.13448308E-01,-0.10016907E-01,
     + -0.47541666E-02,-0.10301435E-01,-0.16987344E-01, 0.62766722E-02,
     + -0.55671044E-01, 0.17607365E-01,-0.28968861E-01,-0.69906586E-02,
     + -0.11730463E+00, 0.33640355E-01, 0.44367900E-02, 0.27638169E-01,
     + -0.16988741E-01, 0.28517402E-02,-0.29291278E-01,-0.14410096E-01,
     +  0.23666890E-01, 0.39484468E-02,-0.19993817E-01,-0.62320633E-02,
     + -0.17559461E+00,-0.15252965E-01,-0.14568388E-01, 0.42886548E-01,
     +  0.56759279E-01,-0.44276832E-04, 0.47185514E-01, 0.33914743E-02,
     + -0.31455252E-01, 0.70841260E-01, 0.74985236E-01, 0.46229500E-01,
     +  0.11574425E-01, 0.54075297E-01,-0.57665273E-02,-0.43941177E-01,
     +  0.66044880E-02,-0.12135608E+00, 0.17911132E-01,-0.71548447E-02,
     +  0.76636836E-01, 0.37473451E-01, 0.39408123E-02,-0.53237891E-02,
     + -0.10727276E+00, 0.50460417E-01,-0.48905239E-01, 0.46938069E-01,
     + -0.21484151E-01,-0.22600582E-01, 0.57666665E-02,-0.57680784E-02,
     +  0.31575758E-01, 0.11946619E-01,-0.35126727E-01,-0.11678522E-01,
     + -0.18962807E-02,-0.91199100E-01,
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
      x_m12_0_9   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22            
     5  +coeff(  5)        *x32        
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21        *x51
      x_m12_0_9   =x_m12_0_9   
     9  +coeff(  9)                *x52
     1  +coeff( 10)*x11                
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x21    *x42    
     4  +coeff( 13)    *x22        *x51
     5  +coeff( 14)        *x32    *x51
     6  +coeff( 15)        *x31*x41*x51
     7  +coeff( 16)            *x42*x51
     8  +coeff( 17)    *x21        *x52
      x_m12_0_9   =x_m12_0_9   
     9  +coeff( 18)                *x53
     1  +coeff( 19)*x11*x21            
     2  +coeff( 20)    *x22*x31*x41    
     3  +coeff( 21)    *x22    *x42    
     4  +coeff( 22)*x11            *x51
     5  +coeff( 23)    *x23        *x51
     6  +coeff( 24)    *x21*x31*x41*x51
     7  +coeff( 25)    *x21    *x42*x51
     8  +coeff( 26)    *x22        *x52
      x_m12_0_9   =x_m12_0_9   
     9  +coeff( 27)        *x31*x41*x52
     1  +coeff( 28)            *x42*x52
     2  +coeff( 29)    *x21        *x53
     3  +coeff( 30)*x11*x22            
     4  +coeff( 31)    *x23*x32        
     5  +coeff( 32)    *x23*x31*x41    
     6  +coeff( 33)    *x21*x33*x41    
     7  +coeff( 34)    *x23    *x42    
     8  +coeff( 35)    *x21*x31*x43    
      x_m12_0_9   =x_m12_0_9   
     9  +coeff( 36)        *x33*x41*x51
     1  +coeff( 37)    *x22    *x42*x51
     2  +coeff( 38)        *x32*x42*x51
     3  +coeff( 39)    *x23        *x52
     4  +coeff( 40)    *x21*x32    *x52
     5  +coeff( 41)    *x21    *x42*x52
     6  +coeff( 42)    *x22        *x53
     7  +coeff( 43)    *x22*x33*x41    
     8  +coeff( 44)    *x22*x32*x42    
      x_m12_0_9   =x_m12_0_9   
     9  +coeff( 45)    *x22*x31*x43    
     1  +coeff( 46)        *x33*x43    
     2  +coeff( 47)    *x23*x31*x41*x51
     3  +coeff( 48)    *x21*x32*x42*x51
     4  +coeff( 49)    *x21*x31*x43*x51
     5  +coeff( 50)*x11*x21        *x52
     6  +coeff( 51)    *x22*x32    *x52
     7  +coeff( 52)        *x33*x41*x52
     8  +coeff( 53)    *x22    *x42*x52
      x_m12_0_9   =x_m12_0_9   
     9  +coeff( 54)        *x32*x42*x52
     1  +coeff( 55)        *x31*x43*x52
     2  +coeff( 56)    *x23        *x53
     3  +coeff( 57)    *x21    *x42*x53
     4  +coeff( 58)    *x23*x33*x41    
     5  +coeff( 59)    *x21*x33*x43    
     6  +coeff( 60)*x11*x21*x32    *x51
     7  +coeff( 61)    *x22*x33*x41*x51
     8  +coeff( 62)    *x23*x31*x41*x52
      x_m12_0_9   =x_m12_0_9   
     9  +coeff( 63)    *x22*x31*x41*x53
     1  +coeff( 64)    *x22    *x42*x53
     2  +coeff( 65)        *x31*x43*x53
     3  +coeff( 66)    *x22*x33*x43    
     4  +coeff( 67)*x11*x22*x32    *x51
     5  +coeff( 68)    *x21*x33*x43*x51
     6  +coeff( 69)*x11*x21*x31*x41*x52
     7  +coeff( 70)    *x22*x33*x41*x52
     8  +coeff( 71)        *x33*x43*x52
      x_m12_0_9   =x_m12_0_9   
     9  +coeff( 72)    *x23*x32    *x53
     1  +coeff( 73)    *x23*x31*x41*x53
     2  +coeff( 74)    *x21*x31*x43*x53
     3  +coeff( 75)*x11*x22*x33*x41    
     4  +coeff( 76)*x11*x23*x32    *x51
     5  +coeff( 77)    *x23*x33*x41*x52
     6  +coeff( 78)*x11*x22*x32    *x53
     7  +coeff( 79)    *x23*x33*x41*x53
     8  +coeff( 80)*x11*x23*x32    *x53
      x_m12_0_9   =x_m12_0_9   
     9  +coeff( 81)    *x21*x31*x41    
     1  +coeff( 82)    *x22*x32        
     2  +coeff( 83)        *x32*x42    
     3  +coeff( 84)    *x21*x32    *x51
     4  +coeff( 85)    *x21*x32*x42    
     5  +coeff( 86)        *x31*x43*x51
     6  +coeff( 87)    *x21*x31*x41*x52
     7  +coeff( 88)        *x31*x41*x53
     8  +coeff( 89)*x11*x23            
      x_m12_0_9   =x_m12_0_9   
     9  +coeff( 90)    *x23    *x42*x51
c
      return
      end
      function t_m12_0_9   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)   !,xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.3956924E-01/
      data xmin/
     1 -0.49995E-02,-0.79858E-01,-0.99988E-01,-0.35144E-01,-0.25572E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49981E-02, 0.61732E-01, 0.99988E-01, 0.35144E-01, 0.42498E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.16598593E-01,-0.65053990E-02, 0.19631791E+00, 0.75900150E-02,
     + -0.69083890E-03,-0.53700744E-02,-0.11991108E-01, 0.70299067E-01,
     + -0.86301327E-01,-0.14854609E-02, 0.25828660E-03,-0.17799189E-01,
     +  0.52847872E-02,-0.60788775E-03,-0.47841016E-02,-0.95232464E-02,
     + -0.26997041E-01, 0.34338787E-01, 0.22347971E-04,-0.44363558E-01,
     + -0.87172515E-03,-0.73100232E-01, 0.60675718E-03, 0.61579584E-02,
     + -0.25434023E-01,-0.45062002E-01, 0.88599129E-02,-0.51515335E-02,
     + -0.12516170E-01, 0.22828871E-01,-0.34703703E-02,-0.93985163E-02,
     +  0.29772434E-02, 0.38394155E-02,-0.44317436E-02, 0.37988401E-02,
     +  0.94334660E-02, 0.59284889E-02,-0.43010856E-02,-0.30606680E-02,
     + -0.86147673E-01, 0.16970668E-01,-0.38260194E-02, 0.46671494E-02,
     +  0.37318880E-02, 0.68175956E-03, 0.34700630E-02,-0.34659673E-01,
     +  0.35224209E-03,-0.31814747E-02,-0.66256295E-02,-0.66788249E-01,
     + -0.39982172E-02,-0.12324227E+00,-0.22849310E-02, 0.10527954E-01,
     + -0.25820392E-02,-0.13173495E-01,-0.81680417E-02,-0.90642059E-02,
     + -0.74456064E-02, 0.23438223E-02,-0.13775098E-02,-0.22102860E-02,
     + -0.11168640E-02, 0.70631347E-03, 0.72577520E-03,-0.13532690E-01,
     +  0.96670873E-02, 0.88428549E-03,-0.29131632E-01,-0.60065012E-01,
     +  0.29395358E-04,-0.40835649E-01, 0.29804971E-03,-0.57558999E-02,
     + -0.85077627E-03,-0.98868192E-03,-0.74621998E-01, 0.86852763E-03,
     + -0.70786616E-02, 0.58003515E-03,-0.10243051E-02, 0.33369471E-03,
     + -0.60483992E-04,-0.52653211E-02,-0.12416794E-03,-0.17814532E-03,
     +  0.64880500E-03, 0.68220391E-03,
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
      t_m12_0_9   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22            
     5  +coeff(  5)        *x32        
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21        *x51
      t_m12_0_9   =t_m12_0_9   
     9  +coeff(  9)                *x52
     1  +coeff( 10)*x11                
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x21    *x42    
     4  +coeff( 13)    *x22        *x51
     5  +coeff( 14)        *x32    *x51
     6  +coeff( 15)        *x31*x41*x51
     7  +coeff( 16)            *x42*x51
     8  +coeff( 17)    *x21        *x52
      t_m12_0_9   =t_m12_0_9   
     9  +coeff( 18)                *x53
     1  +coeff( 19)*x11*x21            
     2  +coeff( 20)    *x22*x31*x41    
     3  +coeff( 21)        *x33*x41    
     4  +coeff( 22)    *x22    *x42    
     5  +coeff( 23)*x11            *x51
     6  +coeff( 24)    *x23        *x51
     7  +coeff( 25)    *x21*x31*x41*x51
     8  +coeff( 26)    *x21    *x42*x51
      t_m12_0_9   =t_m12_0_9   
     9  +coeff( 27)    *x22        *x52
     1  +coeff( 28)        *x31*x41*x52
     2  +coeff( 29)            *x42*x52
     3  +coeff( 30)    *x21        *x53
     4  +coeff( 31)    *x23*x32        
     5  +coeff( 32)    *x23*x31*x41    
     6  +coeff( 33)    *x21*x33*x41    
     7  +coeff( 34)    *x21*x31*x43    
     8  +coeff( 35)    *x22*x32    *x51
      t_m12_0_9   =t_m12_0_9   
     9  +coeff( 36)        *x33*x41*x51
     1  +coeff( 37)        *x32*x42*x51
     2  +coeff( 38)        *x31*x43*x51
     3  +coeff( 39)    *x23        *x52
     4  +coeff( 40)    *x21*x32    *x52
     5  +coeff( 41)    *x21    *x42*x52
     6  +coeff( 42)    *x22        *x53
     7  +coeff( 43)            *x42*x53
     8  +coeff( 44)    *x22*x33*x41    
      t_m12_0_9   =t_m12_0_9   
     9  +coeff( 45)    *x22*x32*x42    
     1  +coeff( 46)    *x22*x31*x43    
     2  +coeff( 47)        *x33*x43    
     3  +coeff( 48)    *x23*x31*x41*x51
     4  +coeff( 49)    *x21*x33*x41*x51
     5  +coeff( 50)    *x21*x32*x42*x51
     6  +coeff( 51)    *x22*x32    *x52
     7  +coeff( 52)    *x22*x31*x41*x52
     8  +coeff( 53)        *x33*x41*x52
      t_m12_0_9   =t_m12_0_9   
     9  +coeff( 54)    *x22    *x42*x52
     1  +coeff( 55)        *x31*x43*x52
     2  +coeff( 56)    *x23        *x53
     3  +coeff( 57)    *x21*x32    *x53
     4  +coeff( 58)    *x21*x31*x41*x53
     5  +coeff( 59)    *x21    *x42*x53
     6  +coeff( 60)    *x21*x31*x41    
     7  +coeff( 61)    *x22*x32        
     8  +coeff( 62)        *x32*x42    
      t_m12_0_9   =t_m12_0_9   
     9  +coeff( 63)        *x31*x43    
     1  +coeff( 64)    *x21*x32    *x51
     2  +coeff( 65)        *x32    *x52
     3  +coeff( 66)*x11    *x31*x41    
     4  +coeff( 67)*x11        *x42    
     5  +coeff( 68)    *x23    *x42    
     6  +coeff( 69)    *x21*x32*x42    
     7  +coeff( 70)*x11*x21        *x51
     8  +coeff( 71)    *x22*x31*x41*x51
      t_m12_0_9   =t_m12_0_9   
     9  +coeff( 72)    *x22    *x42*x51
     1  +coeff( 73)*x11            *x52
     2  +coeff( 74)    *x21*x31*x41*x52
     3  +coeff( 75)*x11*x21*x31*x41    
     4  +coeff( 76)    *x23*x32    *x51
     5  +coeff( 77)*x11    *x31*x41*x51
     6  +coeff( 78)*x11        *x42*x51
     7  +coeff( 79)    *x23    *x42*x51
     8  +coeff( 80)    *x21*x31*x43*x51
      t_m12_0_9   =t_m12_0_9   
     9  +coeff( 81)        *x32*x42*x52
     1  +coeff( 82)*x11            *x53
     2  +coeff( 83)    *x21*x32        
     3  +coeff( 84)*x11*x22            
     4  +coeff( 85)        *x32    *x53
     5  +coeff( 86)        *x31*x41*x53
     6  +coeff( 87)*x11*x23            
     7  +coeff( 88)*x11*x21*x32        
     8  +coeff( 89)*x11*x22        *x51
      t_m12_0_9   =t_m12_0_9   
     9  +coeff( 90)*x11*x21        *x52
c
      return
      end
      function y_m12_0_9   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)   !,xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49995E-02,-0.79858E-01,-0.99988E-01,-0.35144E-01,-0.25572E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49981E-02, 0.61732E-01, 0.99988E-01, 0.35144E-01, 0.42498E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.19038931E+00, 0.74670756E+00, 0.65437605E-03, 0.12269432E-01,
     + -0.14100486E-01,-0.23934845E-01,-0.17555127E-01,-0.50544642E-01,
     + -0.34486854E-02,-0.15820669E-01,-0.74832309E-02,-0.28405549E-01,
     + -0.69351508E-02,-0.31496108E-01, 0.11547618E-03,-0.47761332E-01,
     + -0.14427456E+00, 0.49566817E-01, 0.35408877E-01,-0.83388701E-01,
     +  0.11302146E-01, 0.50551422E-01,-0.55620592E-01, 0.17624054E-01,
     +  0.38211542E-03, 0.15626846E-01, 0.10170913E+00,-0.34202664E-03,
     + -0.81334921E-03,-0.35159830E-01, 0.20192135E-01, 0.89503385E-01,
     + -0.64002544E-01, 0.14684326E-01, 0.60029067E-01, 0.22993010E+00,
     +  0.28408217E+00, 0.14967092E-02, 0.28774717E+00, 0.68869034E-03,
     +  0.41109377E+00,-0.22113169E-01, 0.82425671E-02,-0.11750922E-01,
     + -0.98627293E-02, 0.24855895E+00, 0.21222189E-01,-0.14918759E-01,
     + -0.98173460E-03,-0.13251824E-01, 0.10841373E-01,-0.23523676E-02,
     +  0.37902944E-01,-0.11385813E-01, 0.29255734E-02, 0.10530447E+00,
     + -0.10769754E-02, 0.12241744E+00,-0.45516086E-02,-0.32820784E-01,
     +  0.45885675E-01,-0.17045634E-01, 0.23226145E+00,-0.14738275E-01,
     + -0.18251242E-01,-0.28855260E-01, 0.14589806E+00,-0.62768185E-03,
     +  0.97629078E-01, 0.27227929E+00, 0.12478407E-01, 0.15028772E-01,
     + -0.19864524E-01,-0.37444576E-02, 0.16359954E+00, 0.15922907E+00,
     +  0.60993787E-01, 0.10749206E-01, 0.22905510E-01,-0.18820760E+00,
     + -0.86398161E-03,-0.18215513E-01, 0.29145626E-02, 0.65155770E-03,
     + -0.67403931E-02,-0.68768896E-02,-0.90847380E-03, 0.10441930E-01,
     + -0.49255195E-03,-0.16990418E-01,
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
      y_m12_0_9   =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21*x31        
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)        *x31    *x51
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)    *x22*x31        
     8  +coeff(  8)    *x22    *x41    
      y_m12_0_9   =y_m12_0_9   
     9  +coeff(  9)        *x32*x41    
     1  +coeff( 10)            *x43    
     2  +coeff( 11)    *x21*x31    *x51
     3  +coeff( 12)    *x21    *x41*x51
     4  +coeff( 13)        *x31    *x52
     5  +coeff( 14)            *x41*x52
     6  +coeff( 15)*x11    *x31        
     7  +coeff( 16)    *x23*x31        
     8  +coeff( 17)    *x23    *x41    
      y_m12_0_9   =y_m12_0_9   
     9  +coeff( 18)    *x21*x31*x42    
     1  +coeff( 19)    *x21    *x43    
     2  +coeff( 20)    *x22    *x41*x51
     3  +coeff( 21)        *x32*x41*x51
     4  +coeff( 22)            *x43*x51
     5  +coeff( 23)    *x21    *x41*x52
     6  +coeff( 24)            *x41*x53
     7  +coeff( 25)*x11*x21    *x41    
     8  +coeff( 26)    *x22*x32*x41    
      y_m12_0_9   =y_m12_0_9   
     9  +coeff( 27)    *x22*x31*x42    
     1  +coeff( 28)        *x33*x42    
     2  +coeff( 29)    *x23*x31    *x51
     3  +coeff( 30)    *x23    *x41*x51
     4  +coeff( 31)    *x21*x32*x41*x51
     5  +coeff( 32)    *x21*x31*x42*x51
     6  +coeff( 33)    *x22    *x41*x52
     7  +coeff( 34)            *x43*x52
     8  +coeff( 35)    *x23*x32*x41    
      y_m12_0_9   =y_m12_0_9   
     9  +coeff( 36)    *x23*x31*x42    
     1  +coeff( 37)    *x23    *x43    
     2  +coeff( 38)    *x22*x33    *x51
     3  +coeff( 39)    *x22*x31*x42*x51
     4  +coeff( 40)        *x33*x42*x51
     5  +coeff( 41)    *x22    *x43*x51
     6  +coeff( 42)    *x23    *x41*x52
     7  +coeff( 43)        *x33    *x53
     8  +coeff( 44)    *x23*x33    *x51
      y_m12_0_9   =y_m12_0_9   
     9  +coeff( 45)    *x23*x32*x41*x51
     1  +coeff( 46)    *x23    *x43*x51
     2  +coeff( 47)        *x33*x42*x52
     3  +coeff( 48)    *x22*x33*x42*x51
     4  +coeff( 49)        *x33        
     5  +coeff( 50)        *x31*x42    
     6  +coeff( 51)    *x21*x32*x41    
     7  +coeff( 52)        *x33    *x51
     8  +coeff( 53)        *x31*x42*x51
      y_m12_0_9   =y_m12_0_9   
     9  +coeff( 54)    *x21*x31    *x52
     1  +coeff( 55)    *x22*x33        
     2  +coeff( 56)    *x22    *x43    
     3  +coeff( 57)        *x32*x43    
     4  +coeff( 58)    *x21    *x43*x51
     5  +coeff( 59)        *x31*x42*x52
     6  +coeff( 60)    *x21*x33*x42    
     7  +coeff( 61)    *x22*x32*x41*x51
     8  +coeff( 62)        *x32*x43*x51
      y_m12_0_9   =y_m12_0_9   
     9  +coeff( 63)    *x21    *x43*x52
     1  +coeff( 64)    *x22    *x41*x53
     2  +coeff( 65)        *x31*x42*x53
     3  +coeff( 66)    *x22*x33*x42    
     4  +coeff( 67)    *x23*x31*x42*x51
     5  +coeff( 68)    *x22*x33    *x52
     6  +coeff( 69)    *x22*x31*x42*x52
     7  +coeff( 70)    *x22    *x43*x52
     8  +coeff( 71)    *x21*x33    *x53
      y_m12_0_9   =y_m12_0_9   
     9  +coeff( 72)    *x21*x32*x41*x53
     1  +coeff( 73)    *x21*x31*x42*x53
     2  +coeff( 74)*x11*x22    *x41*x52
     3  +coeff( 75)    *x23*x31*x42*x52
     4  +coeff( 76)    *x21*x33*x42*x52
     5  +coeff( 77)    *x23    *x43*x52
     6  +coeff( 78)        *x32*x43*x53
     7  +coeff( 79)    *x22*x33*x42*x52
     8  +coeff( 80)    *x23*x33*x42*x52
      y_m12_0_9   =y_m12_0_9   
     9  +coeff( 81)    *x21*x33        
     1  +coeff( 82)    *x22*x31    *x51
     2  +coeff( 83)        *x31    *x53
     3  +coeff( 84)*x11        *x41*x51
     4  +coeff( 85)    *x22*x31    *x52
     5  +coeff( 86)    *x21    *x41*x53
     6  +coeff( 87)*x11*x22*x31        
     7  +coeff( 88)    *x23*x33        
     8  +coeff( 89)*x11        *x43    
      y_m12_0_9   =y_m12_0_9   
     9  +coeff( 90)    *x21*x32*x43    
c
      return
      end
      function p_m12_0_9   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)   !,xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49995E-02,-0.79858E-01,-0.99988E-01,-0.35144E-01,-0.25572E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49981E-02, 0.61732E-01, 0.99988E-01, 0.35144E-01, 0.42498E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.21361386E-01,-0.65962948E-01, 0.23793953E-02, 0.28477242E-03,
     +  0.15598949E-01, 0.62234905E-01,-0.15776058E-02,-0.35211988E-03,
     + -0.82437340E-02,-0.21856695E-02,-0.43185204E-02,-0.51339129E-02,
     + -0.45529669E-02,-0.11947159E-01,-0.15151583E-01,-0.52889857E-01,
     +  0.30345062E-04,-0.31288937E-02, 0.14269886E-03, 0.14722104E-01,
     + -0.72478724E-03,-0.29807793E-01, 0.82615716E-02, 0.35264309E-01,
     +  0.32742299E-01,-0.64070597E-02,-0.39426159E-01, 0.78877416E-02,
     +  0.17573847E-01, 0.41363764E-03,-0.16897173E-02, 0.13525239E-02,
     + -0.18728442E-02,-0.15214951E-02, 0.18016397E-02, 0.52215140E-02,
     +  0.32250630E-02,-0.54800248E-03, 0.47013781E-03,-0.29305998E-01,
     +  0.42704973E-01, 0.45131512E-01,-0.35639308E-01,-0.13512205E-01,
     + -0.11133980E-01, 0.24545018E-02,-0.68819057E-03,-0.29327968E-03,
     +  0.17692416E-02, 0.12259594E-02, 0.17658173E-03, 0.18773122E-02,
     +  0.22270158E-01,-0.52638579E-03, 0.16468950E-01, 0.60517463E-03,
     + -0.38279595E-02, 0.28455886E-02,-0.46105129E-02, 0.12027806E-02,
     +  0.63204087E-01, 0.48339696E-03,-0.37259089E-02,-0.19219976E-01,
     +  0.13831013E-01, 0.48821107E-01, 0.90787575E-01,-0.30632520E-02,
     + -0.10761647E-01, 0.26605383E-02, 0.30330345E-01,-0.27421961E-03,
     +  0.17223716E-01, 0.22017831E-01, 0.65690913E-03, 0.18388206E-03,
     + -0.59649386E-02, 0.38791180E-02,-0.25884295E-01, 0.14164157E-01,
     +  0.13558345E-01,-0.56971400E-02, 0.72356185E-03, 0.56702336E-02,
     +  0.56822640E-02, 0.47047567E-01, 0.10478655E-02,-0.10293466E-01,
     + -0.35192762E-01, 0.57473395E-03,
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
      p_m12_0_9   =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21*x31        
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)        *x31    *x51
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)    *x22*x31        
     8  +coeff(  8)        *x33        
      p_m12_0_9   =p_m12_0_9   
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)        *x32*x41    
     2  +coeff( 11)        *x31*x42    
     3  +coeff( 12)            *x43    
     4  +coeff( 13)    *x21*x31    *x51
     5  +coeff( 14)    *x21    *x41*x51
     6  +coeff( 15)        *x31    *x52
     7  +coeff( 16)            *x41*x52
     8  +coeff( 17)*x11    *x31        
      p_m12_0_9   =p_m12_0_9   
     9  +coeff( 18)    *x21*x33        
     1  +coeff( 19)*x11        *x41    
     2  +coeff( 20)    *x21    *x43    
     3  +coeff( 21)    *x22*x31    *x51
     4  +coeff( 22)    *x22    *x41*x51
     5  +coeff( 23)        *x32*x41*x51
     6  +coeff( 24)        *x31*x42*x51
     7  +coeff( 25)            *x43*x51
     8  +coeff( 26)    *x21*x31    *x52
      p_m12_0_9   =p_m12_0_9   
     9  +coeff( 27)    *x21    *x41*x52
     1  +coeff( 28)        *x31    *x53
     2  +coeff( 29)            *x41*x53
     3  +coeff( 30)*x11*x21*x31        
     4  +coeff( 31)    *x22*x33        
     5  +coeff( 32)*x11*x21    *x41    
     6  +coeff( 33)    *x22*x32*x41    
     7  +coeff( 34)    *x22*x31*x42    
     8  +coeff( 35)        *x33*x42    
      p_m12_0_9   =p_m12_0_9   
     9  +coeff( 36)    *x22    *x43    
     1  +coeff( 37)        *x32*x43    
     2  +coeff( 38)    *x21*x33    *x51
     3  +coeff( 39)*x11        *x41*x51
     4  +coeff( 40)    *x23    *x41*x51
     5  +coeff( 41)    *x21*x31*x42*x51
     6  +coeff( 42)    *x21    *x43*x51
     7  +coeff( 43)    *x22    *x41*x52
     8  +coeff( 44)        *x31*x42*x52
      p_m12_0_9   =p_m12_0_9   
     9  +coeff( 45)            *x43*x52
     1  +coeff( 46)    *x21*x31    *x53
     2  +coeff( 47)    *x21    *x41*x53
     3  +coeff( 48)*x11*x22*x31        
     4  +coeff( 49)*x11*x22    *x41    
     5  +coeff( 50)    *x23*x32*x41    
     6  +coeff( 51)*x11    *x31*x42    
     7  +coeff( 52)    *x23*x31*x42    
     8  +coeff( 53)    *x21*x33*x42    
      p_m12_0_9   =p_m12_0_9   
     9  +coeff( 54)*x11        *x43    
     1  +coeff( 55)    *x21*x32*x43    
     2  +coeff( 56)*x11*x21*x31    *x51
     3  +coeff( 57)    *x22*x33    *x51
     4  +coeff( 58)*x11*x21    *x41*x51
     5  +coeff( 59)    *x22*x32*x41*x51
     6  +coeff( 60)        *x33*x42*x51
     7  +coeff( 61)    *x22    *x43*x51
     8  +coeff( 62)        *x32*x43*x51
      p_m12_0_9   =p_m12_0_9   
     9  +coeff( 63)    *x23*x31    *x52
     1  +coeff( 64)    *x23    *x41*x52
     2  +coeff( 65)    *x21*x32*x41*x52
     3  +coeff( 66)    *x21*x31*x42*x52
     4  +coeff( 67)    *x21    *x43*x52
     5  +coeff( 68)    *x22*x31    *x53
     6  +coeff( 69)    *x22    *x41*x53
     7  +coeff( 70)        *x31*x42*x53
     8  +coeff( 71)            *x43*x53
      p_m12_0_9   =p_m12_0_9   
     9  +coeff( 72)*x11*x23*x31        
     1  +coeff( 73)    *x22*x33*x42    
     2  +coeff( 74)    *x22*x32*x43    
     3  +coeff( 75)*x11*x22*x31    *x51
     4  +coeff( 76)*x11    *x33    *x51
     5  +coeff( 77)    *x23*x33    *x51
     6  +coeff( 78)*x11*x22    *x41*x51
     7  +coeff( 79)    *x23*x31*x42*x51
     8  +coeff( 80)    *x23    *x43*x51
      p_m12_0_9   =p_m12_0_9   
     9  +coeff( 81)    *x21*x32*x43*x51
     1  +coeff( 82)    *x22*x33    *x52
     2  +coeff( 83)*x11*x21    *x41*x52
     3  +coeff( 84)    *x22*x31*x42*x52
     4  +coeff( 85)        *x33*x42*x52
     5  +coeff( 86)    *x22    *x43*x52
     6  +coeff( 87)    *x21*x33    *x53
     7  +coeff( 88)    *x21*x32*x41*x53
     8  +coeff( 89)    *x21*x31*x42*x53
      p_m12_0_9   =p_m12_0_9   
     9  +coeff( 90)*x11*x22*x33        
c
      return
      end
      real function l_m12_0_9   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)   !,xmean(10)
      dimension coeff( 41)
      data ncoeff/ 40/
      data avdat/  0.1695417E-01/
      data xmin/
     1 -0.49995E-02,-0.79858E-01,-0.99988E-01,-0.35144E-01,-0.25572E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49981E-02, 0.61732E-01, 0.99988E-01, 0.35144E-01, 0.42498E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.43124871E-03,-0.15167284E+00,-0.28707093E-01,-0.27323334E-01,
     + -0.33492561E-01,-0.30230159E-01,-0.12486684E-01,-0.24140501E-01,
     + -0.93537159E-02, 0.14286153E-01, 0.34362171E-01,-0.16626935E-01,
     +  0.32660630E-01, 0.16008282E-01,-0.93529886E-03, 0.13358057E-01,
     +  0.21868730E-01,-0.19644005E-01,-0.13959949E-02, 0.25921427E-02,
     + -0.32546390E-02,-0.24482150E-01, 0.10026784E-01,-0.16581712E-01,
     +  0.53815061E-03, 0.22445053E-02,-0.46505022E-03, 0.61388372E-03,
     +  0.55572845E-03,-0.16224446E-01,-0.29502690E-02, 0.72555915E-02,
     +  0.15230780E-01, 0.48667556E-02, 0.46129432E-03, 0.18098515E-02,
     +  0.28439043E-02,-0.51588420E-03,-0.19774262E-02, 0.60427412E-02,
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
      l_m12_0_9   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x22            
     4  +coeff(  4)                *x52
     5  +coeff(  5)            *x42    
     6  +coeff(  6)    *x21        *x52
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)                *x51
      l_m12_0_9   =l_m12_0_9   
     9  +coeff(  9)    *x21        *x51
     1  +coeff( 10)    *x22        *x51
     2  +coeff( 11)                *x53
     3  +coeff( 12)    *x22        *x52
     4  +coeff( 13)    *x21        *x53
     5  +coeff( 14)    *x22    *x42*x51
     6  +coeff( 15)*x11                
     7  +coeff( 16)        *x31*x41*x51
     8  +coeff( 17)            *x42*x51
      l_m12_0_9   =l_m12_0_9   
     9  +coeff( 18)                *x54
     1  +coeff( 19)        *x32        
     2  +coeff( 20)        *x32    *x51
     3  +coeff( 21)    *x24            
     4  +coeff( 22)            *x42*x52
     5  +coeff( 23)    *x22        *x53
     6  +coeff( 24)    *x21        *x54
     7  +coeff( 25)        *x33*x41*x52
     8  +coeff( 26)        *x31*x41*x54
      l_m12_0_9   =l_m12_0_9   
     9  +coeff( 27)*x11*x21            
     1  +coeff( 28)*x11            *x51
     2  +coeff( 29)*x11*x21        *x51
     3  +coeff( 30)        *x31*x41*x52
     4  +coeff( 31)    *x24        *x51
     5  +coeff( 32)    *x22*x31*x41*x51
     6  +coeff( 33)            *x42*x53
     7  +coeff( 34)    *x21    *x42*x53
     8  +coeff( 35)    *x23            
      l_m12_0_9   =l_m12_0_9   
     9  +coeff( 36)    *x22    *x42    
     1  +coeff( 37)    *x21    *x42*x51
     2  +coeff( 38)*x11            *x52
     3  +coeff( 39)        *x32    *x52
     4  +coeff( 40)        *x31*x41*x53
c
      return
      end
      function x_m12_0_11  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)   !,xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.1048844E+00/
      data xmin/
     1 -0.49995E-02,-0.79804E-01,-0.99988E-01,-0.35144E-01,-0.21223E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49981E-02, 0.61732E-01, 0.99988E-01, 0.35144E-01, 0.42498E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.12106877E+00, 0.28752497E+00, 0.67640460E+00, 0.32142211E-01,
     + -0.14848100E-02,-0.19035213E-01,-0.35838023E-01, 0.36431351E+00,
     + -0.22150937E+00,-0.40544546E-02,-0.58656174E-03,-0.13100653E-02,
     + -0.73669389E-01, 0.17259203E-01,-0.20133040E-02,-0.87939957E-02,
     + -0.11182774E-01,-0.85021153E-01, 0.74540250E-01, 0.26486591E-02,
     + -0.34314699E-01,-0.17547739E+00, 0.32905838E-04,-0.31064335E+00,
     +  0.30874026E-02,-0.29438213E-02,-0.99348770E-02,-0.94770424E-01,
     + -0.16624799E+00,-0.30726746E-02,-0.25940364E-01,-0.41564193E-01,
     +  0.27702453E-01, 0.79888274E-03,-0.22189254E-01,-0.36713764E-01,
     + -0.13041997E-01,-0.72841972E-01,-0.11080801E-01, 0.91890171E-02,
     +  0.16390786E-01, 0.17666321E-01,-0.48854891E-01,-0.19612793E-01,
     + -0.22364780E+00,-0.16529748E-01, 0.63255660E-01,-0.20074131E+00,
     + -0.31194773E+00, 0.20051403E-01,-0.32331410E-02, 0.27011923E-01,
     +  0.29714989E-01, 0.66506468E-01, 0.58567233E-01, 0.19865859E-01,
     +  0.28444810E-01, 0.14922996E-01,-0.15778234E-01, 0.95709451E-01,
     + -0.37894135E-02, 0.70009761E-01, 0.14631560E-01,-0.50909162E-01,
     +  0.20168997E-01,-0.31630840E-01, 0.16014299E-02,-0.62899212E-02,
     +  0.30410266E-01, 0.21367171E-02,-0.10709734E+00,-0.23064864E+00,
     + -0.96823186E-01,-0.38668588E-02,-0.36491645E-02,-0.10620375E-01,
     + -0.26033486E-02,-0.78900844E-01, 0.99054761E-02, 0.11105212E-01,
     + -0.32414016E-03,-0.13858901E+00,-0.26146155E-02, 0.13130269E-02,
     +  0.45774425E-02, 0.33725102E-01,-0.33720855E-02, 0.49054265E-01,
     + -0.23895148E-02, 0.15148030E-01,
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
      x_m12_0_11  =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22            
     5  +coeff(  5)        *x32        
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21        *x51
      x_m12_0_11  =x_m12_0_11  
     9  +coeff(  9)                *x52
     1  +coeff( 10)*x11                
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x21*x32        
     4  +coeff( 13)    *x21    *x42    
     5  +coeff( 14)    *x22        *x51
     6  +coeff( 15)        *x32    *x51
     7  +coeff( 16)        *x31*x41*x51
     8  +coeff( 17)            *x42*x51
      x_m12_0_11  =x_m12_0_11  
     9  +coeff( 18)    *x21        *x52
     1  +coeff( 19)                *x53
     2  +coeff( 20)*x11*x21            
     3  +coeff( 21)    *x22*x32        
     4  +coeff( 22)    *x22*x31*x41    
     5  +coeff( 23)        *x33*x41    
     6  +coeff( 24)    *x22    *x42    
     7  +coeff( 25)*x11            *x51
     8  +coeff( 26)    *x23        *x51
      x_m12_0_11  =x_m12_0_11  
     9  +coeff( 27)    *x21*x32    *x51
     1  +coeff( 28)    *x21*x31*x41*x51
     2  +coeff( 29)    *x21    *x42*x51
     3  +coeff( 30)    *x22        *x52
     4  +coeff( 31)        *x31*x41*x52
     5  +coeff( 32)            *x42*x52
     6  +coeff( 33)    *x21        *x53
     7  +coeff( 34)*x11*x22            
     8  +coeff( 35)    *x23*x32        
      x_m12_0_11  =x_m12_0_11  
     9  +coeff( 36)    *x23*x31*x41    
     1  +coeff( 37)    *x21*x33*x41    
     2  +coeff( 38)    *x23    *x42    
     3  +coeff( 39)    *x21*x31*x43    
     4  +coeff( 40)        *x33*x41*x51
     5  +coeff( 41)        *x32*x42*x51
     6  +coeff( 42)        *x31*x43*x51
     7  +coeff( 43)    *x23        *x52
     8  +coeff( 44)    *x21*x32    *x52
      x_m12_0_11  =x_m12_0_11  
     9  +coeff( 45)    *x21    *x42*x52
     1  +coeff( 46)            *x42*x53
     2  +coeff( 47)    *x22*x32*x42    
     3  +coeff( 48)    *x23    *x42*x51
     4  +coeff( 49)    *x22    *x42*x52
     5  +coeff( 50)    *x23        *x53
     6  +coeff( 51)    *x21*x31*x41*x53
     7  +coeff( 52)    *x21    *x42*x53
     8  +coeff( 53)    *x23*x33*x41    
      x_m12_0_11  =x_m12_0_11  
     9  +coeff( 54)    *x23*x32*x42    
     1  +coeff( 55)    *x21*x33*x43    
     2  +coeff( 56)    *x22*x33*x41*x51
     3  +coeff( 57)    *x22*x32*x42*x51
     4  +coeff( 58)    *x23    *x42*x52
     5  +coeff( 59)        *x33*x41*x53
     6  +coeff( 60)    *x22    *x42*x53
     7  +coeff( 61)*x11*x23*x32        
     8  +coeff( 62)    *x22*x33*x43    
      x_m12_0_11  =x_m12_0_11  
     9  +coeff( 63)    *x23*x32*x42*x51
     1  +coeff( 64)    *x22*x31*x43*x52
     2  +coeff( 65)        *x33*x43*x52
     3  +coeff( 66)    *x21*x31*x41    
     4  +coeff( 67)        *x32*x42    
     5  +coeff( 68)        *x32    *x52
     6  +coeff( 69)    *x21*x32*x42    
     7  +coeff( 70)*x11*x21        *x51
     8  +coeff( 71)    *x22*x31*x41*x51
      x_m12_0_11  =x_m12_0_11  
     9  +coeff( 72)    *x22    *x42*x51
     1  +coeff( 73)    *x21*x31*x41*x52
     2  +coeff( 74)        *x31*x41*x53
     3  +coeff( 75)*x11*x23            
     4  +coeff( 76)*x11*x21    *x42    
     5  +coeff( 77)    *x23*x32    *x51
     6  +coeff( 78)    *x23*x31*x41*x51
     7  +coeff( 79)    *x21*x33*x41*x51
     8  +coeff( 80)    *x21*x31*x43*x51
      x_m12_0_11  =x_m12_0_11  
     9  +coeff( 81)    *x22*x32    *x52
     1  +coeff( 82)    *x22*x31*x41*x52
     2  +coeff( 83)*x11*x22*x32        
     3  +coeff( 84)*x11*x21*x31*x41*x51
     4  +coeff( 85)*x11*x21    *x42*x51
     5  +coeff( 86)    *x23*x32    *x52
     6  +coeff( 87)*x11*x21        *x53
     7  +coeff( 88)    *x22*x31*x41*x53
     8  +coeff( 89)*x11*x21*x33*x41    
      x_m12_0_11  =x_m12_0_11  
     9  +coeff( 90)*x11*x23    *x42    
c
      return
      end
      function t_m12_0_11  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)   !,xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.4781470E-01/
      data xmin/
     1 -0.49995E-02,-0.79804E-01,-0.99988E-01,-0.35144E-01,-0.21223E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49981E-02, 0.61732E-01, 0.99988E-01, 0.35144E-01, 0.42498E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.20845897E-01,-0.16522192E-02, 0.17616484E+00, 0.84403148E-02,
     + -0.36520645E-03,-0.27504738E-02,-0.69455910E-02, 0.64830482E-01,
     + -0.67857310E-01,-0.14313029E-02, 0.87941153E-03,-0.11231332E-02,
     + -0.20946227E-01, 0.76167849E-02,-0.11025799E-02,-0.63464139E-02,
     + -0.13482980E-01,-0.18861314E-01, 0.25350396E-01, 0.23536611E-03,
     + -0.77558164E-02,-0.46264105E-01,-0.80115825E-03,-0.78238375E-01,
     +  0.80818398E-03, 0.50335387E-02,-0.28863789E-02,-0.28371869E-01,
     + -0.52102048E-01, 0.79138903E-02,-0.77629169E-02,-0.13530451E-01,
     +  0.15832836E-01,-0.39902166E-02,-0.13230030E-01, 0.29098399E-02,
     +  0.30424555E-02, 0.20038879E-02, 0.59271622E-02, 0.36105667E-02,
     + -0.53709657E-02,-0.43348181E-02,-0.87630779E-01, 0.86411191E-02,
     + -0.43301350E-02, 0.14517157E-02,-0.77706650E-01,-0.12078549E+00,
     +  0.65455502E-02,-0.13215495E-01,-0.16479502E-01,-0.10717214E-01,
     +  0.32359816E-02,-0.12049330E-02,-0.21412892E-01, 0.89561446E-02,
     +  0.50310435E-03,-0.51630959E-02,-0.41024480E-01,-0.81260085E-01,
     + -0.43511961E-01,-0.40101386E-02,-0.32003951E-03,-0.25460686E-03,
     +  0.37390394E-02,-0.20268069E-03, 0.37074557E-02,-0.37926059E-01,
     +  0.24079843E-03,-0.27965968E-02,-0.67743473E-02,-0.63573994E-01,
     + -0.54913647E-02,-0.15876686E-02, 0.26515021E-03,-0.33931417E-05,
     +  0.58116775E-03, 0.54832606E-03, 0.36910686E-03, 0.10733300E-03,
     + -0.52012657E-02,-0.64657378E-03,-0.70415233E-03,-0.62250057E-02,
     + -0.12644884E-02,-0.15581039E-02,-0.71946735E-03, 0.13175198E-03,
     + -0.25680459E-02, 0.27672574E-03,
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
      t_m12_0_11  =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22            
     5  +coeff(  5)        *x32        
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21        *x51
      t_m12_0_11  =t_m12_0_11  
     9  +coeff(  9)                *x52
     1  +coeff( 10)*x11                
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x21*x32        
     4  +coeff( 13)    *x21    *x42    
     5  +coeff( 14)    *x22        *x51
     6  +coeff( 15)        *x32    *x51
     7  +coeff( 16)        *x31*x41*x51
     8  +coeff( 17)            *x42*x51
      t_m12_0_11  =t_m12_0_11  
     9  +coeff( 18)    *x21        *x52
     1  +coeff( 19)                *x53
     2  +coeff( 20)*x11*x21            
     3  +coeff( 21)    *x22*x32        
     4  +coeff( 22)    *x22*x31*x41    
     5  +coeff( 23)        *x33*x41    
     6  +coeff( 24)    *x22    *x42    
     7  +coeff( 25)*x11            *x51
     8  +coeff( 26)    *x23        *x51
      t_m12_0_11  =t_m12_0_11  
     9  +coeff( 27)    *x21*x32    *x51
     1  +coeff( 28)    *x21*x31*x41*x51
     2  +coeff( 29)    *x21    *x42*x51
     3  +coeff( 30)    *x22        *x52
     4  +coeff( 31)        *x31*x41*x52
     5  +coeff( 32)            *x42*x52
     6  +coeff( 33)    *x21        *x53
     7  +coeff( 34)    *x23*x32        
     8  +coeff( 35)    *x23*x31*x41    
      t_m12_0_11  =t_m12_0_11  
     9  +coeff( 36)    *x21*x33*x41    
     1  +coeff( 37)    *x21*x31*x43    
     2  +coeff( 38)        *x33*x41*x51
     3  +coeff( 39)        *x32*x42*x51
     4  +coeff( 40)        *x31*x43*x51
     5  +coeff( 41)    *x23        *x52
     6  +coeff( 42)    *x21*x32    *x52
     7  +coeff( 43)    *x21    *x42*x52
     8  +coeff( 44)    *x22        *x53
      t_m12_0_11  =t_m12_0_11  
     9  +coeff( 45)            *x42*x53
     1  +coeff( 46)    *x22*x32*x42    
     2  +coeff( 47)    *x23    *x42*x51
     3  +coeff( 48)    *x22    *x42*x52
     4  +coeff( 49)    *x23        *x53
     5  +coeff( 50)    *x21*x31*x41*x53
     6  +coeff( 51)    *x21    *x42*x53
     7  +coeff( 52)    *x21*x31*x41    
     8  +coeff( 53)        *x32*x42    
      t_m12_0_11  =t_m12_0_11  
     9  +coeff( 54)        *x32    *x52
     1  +coeff( 55)    *x23    *x42    
     2  +coeff( 56)    *x21*x32*x42    
     3  +coeff( 57)*x11*x21        *x51
     4  +coeff( 58)    *x22*x32    *x51
     5  +coeff( 59)    *x22*x31*x41*x51
     6  +coeff( 60)    *x22    *x42*x51
     7  +coeff( 61)    *x21*x31*x41*x52
     8  +coeff( 62)        *x31*x41*x53
      t_m12_0_11  =t_m12_0_11  
     9  +coeff( 63)*x11*x23            
     1  +coeff( 64)*x11*x21*x32        
     2  +coeff( 65)    *x22*x33*x41    
     3  +coeff( 66)*x11*x21    *x42    
     4  +coeff( 67)        *x33*x43    
     5  +coeff( 68)    *x23*x31*x41*x51
     6  +coeff( 69)    *x21*x33*x41*x51
     7  +coeff( 70)    *x21*x31*x43*x51
     8  +coeff( 71)    *x22*x32    *x52
      t_m12_0_11  =t_m12_0_11  
     9  +coeff( 72)    *x22*x31*x41*x52
     1  +coeff( 73)        *x32*x42*x52
     2  +coeff( 74)    *x21*x32    *x53
     3  +coeff( 75)*x11*x22            
     4  +coeff( 76)*x11    *x32        
     5  +coeff( 77)*x11    *x31*x41    
     6  +coeff( 78)*x11        *x42    
     7  +coeff( 79)        *x32    *x53
     8  +coeff( 80)*x11*x21*x31*x41    
      t_m12_0_11  =t_m12_0_11  
     9  +coeff( 81)    *x23*x32    *x51
     1  +coeff( 82)*x11    *x31*x41*x51
     2  +coeff( 83)*x11        *x42*x51
     3  +coeff( 84)    *x21*x32*x42*x51
     4  +coeff( 85)        *x33*x41*x52
     5  +coeff( 86)        *x31*x43*x52
     6  +coeff( 87)        *x31*x43    
     7  +coeff( 88)*x11            *x52
     8  +coeff( 89)    *x22*x31*x43    
      t_m12_0_11  =t_m12_0_11  
     9  +coeff( 90)*x11*x21        *x52
c
      return
      end
      function y_m12_0_11  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)   !,xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49995E-02,-0.79804E-01,-0.99988E-01,-0.35144E-01,-0.21223E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49981E-02, 0.61732E-01, 0.99988E-01, 0.35144E-01, 0.42498E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.16872309E+00, 0.68396485E+00,-0.11620645E-03, 0.10271362E-01,
     + -0.36421365E-02, 0.17726770E-01,-0.59171749E-03,-0.76652549E-01,
     + -0.44612391E-02,-0.17518282E-01,-0.12003219E-01,-0.58438301E-01,
     + -0.52807622E-01,-0.79140642E-04,-0.49813915E-01,-0.38191152E-03,
     + -0.16991495E+00, 0.37831198E-01, 0.41959766E-01,-0.19971855E-01,
     + -0.14561857E+00, 0.13236955E-01, 0.72783820E-01, 0.84002636E-01,
     + -0.17731439E-01,-0.89388162E-01, 0.20028014E-01, 0.77354372E-03,
     +  0.66904438E-03, 0.30683914E-01, 0.13424851E+00, 0.19084552E+00,
     + -0.23304450E-02, 0.26699126E-01, 0.16166772E+00, 0.25596812E+00,
     +  0.16279147E-02,-0.33737604E-04,-0.91272928E-01,-0.47262893E-02,
     +  0.81616919E-02, 0.24630908E-01, 0.43401821E-02, 0.64967363E-02,
     +  0.74869029E-01, 0.29965821E+00,-0.32610210E-03, 0.39109540E+00,
     +  0.68701625E-01, 0.46823096E+00,-0.83599836E-02, 0.18233601E-01,
     + -0.32091031E-02,-0.13940545E-01, 0.11405788E+00, 0.33233517E+00,
     + -0.17061399E-01, 0.14674156E-01,-0.18830860E+00, 0.15095100E+00,
     + -0.14011085E-01,-0.54001528E+00,-0.28564206E+00,-0.20265909E+00,
     +  0.18847764E-02, 0.33070537E+00,-0.10244620E-01,-0.75407869E+00,
     +  0.90657875E-01,-0.39324716E+00,-0.22402222E+00,-0.17394597E-01,
     + -0.23014493E+00, 0.19294646E-01, 0.37327451E+00, 0.50791707E-01,
     +  0.58537108E+00, 0.67754090E-02, 0.17929125E+00, 0.56935179E+00,
     +  0.16162586E+00,-0.20475307E-01,-0.10200376E-01,-0.12851685E-01,
     +  0.13008821E-01,-0.14197762E-02, 0.65814024E-02,-0.63905108E-03,
     +  0.45659272E-02,-0.74837811E-01,
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
      y_m12_0_11  =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21*x31        
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)        *x31    *x51
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)        *x33        
     8  +coeff(  8)    *x22    *x41    
      y_m12_0_11  =y_m12_0_11  
     9  +coeff(  9)        *x32*x41    
     1  +coeff( 10)            *x43    
     2  +coeff( 11)    *x21*x31    *x51
     3  +coeff( 12)    *x21    *x41*x51
     4  +coeff( 13)            *x41*x52
     5  +coeff( 14)*x11    *x31        
     6  +coeff( 15)    *x23*x31        
     7  +coeff( 16)*x11        *x41    
     8  +coeff( 17)    *x23    *x41    
      y_m12_0_11  =y_m12_0_11  
     9  +coeff( 18)    *x21*x31*x42    
     1  +coeff( 19)    *x21    *x43    
     2  +coeff( 20)    *x22*x31    *x51
     3  +coeff( 21)    *x22    *x41*x51
     4  +coeff( 22)        *x32*x41*x51
     5  +coeff( 23)        *x31*x42*x51
     6  +coeff( 24)            *x43*x51
     7  +coeff( 25)    *x21*x31    *x52
     8  +coeff( 26)    *x21    *x41*x52
      y_m12_0_11  =y_m12_0_11  
     9  +coeff( 27)            *x41*x53
     1  +coeff( 28)    *x22*x33        
     2  +coeff( 29)*x11*x21    *x41    
     3  +coeff( 30)    *x22*x32*x41    
     4  +coeff( 31)    *x22*x31*x42    
     5  +coeff( 32)    *x22    *x43    
     6  +coeff( 33)    *x23*x31    *x51
     7  +coeff( 34)    *x21*x32*x41*x51
     8  +coeff( 35)    *x21*x31*x42*x51
      y_m12_0_11  =y_m12_0_11  
     9  +coeff( 36)    *x21    *x43*x51
     1  +coeff( 37)    *x22*x31    *x52
     2  +coeff( 38)        *x33    *x52
     3  +coeff( 39)    *x22    *x41*x52
     4  +coeff( 40)        *x32*x41*x52
     5  +coeff( 41)        *x31*x42*x52
     6  +coeff( 42)            *x43*x52
     7  +coeff( 43)    *x21*x31    *x53
     8  +coeff( 44)    *x23*x33        
      y_m12_0_11  =y_m12_0_11  
     9  +coeff( 45)    *x23*x32*x41    
     1  +coeff( 46)    *x23*x31*x42    
     2  +coeff( 47)    *x21*x33*x42    
     3  +coeff( 48)    *x23    *x43    
     4  +coeff( 49)    *x22*x32*x41*x51
     5  +coeff( 50)    *x22    *x43*x51
     6  +coeff( 51)        *x32*x43*x51
     7  +coeff( 52)    *x23*x31    *x52
     8  +coeff( 53)    *x21*x33    *x52
      y_m12_0_11  =y_m12_0_11  
     9  +coeff( 54)    *x23    *x41*x52
     1  +coeff( 55)    *x21*x31*x42*x52
     2  +coeff( 56)    *x21    *x43*x52
     3  +coeff( 57)    *x23*x33    *x51
     4  +coeff( 58)    *x23*x32*x41*x51
     5  +coeff( 59)    *x23*x31*x42*x51
     6  +coeff( 60)    *x23    *x43*x51
     7  +coeff( 61)    *x22*x33    *x52
     8  +coeff( 62)    *x22*x31*x42*x52
      y_m12_0_11  =y_m12_0_11  
     9  +coeff( 63)    *x21*x31*x42*x53
     1  +coeff( 64)    *x21    *x43*x53
     2  +coeff( 65)    *x23*x33*x42    
     3  +coeff( 66)    *x22*x33*x42*x51
     4  +coeff( 67)    *x23*x33    *x52
     5  +coeff( 68)    *x23*x31*x42*x52
     6  +coeff( 69)    *x21*x33*x42*x52
     7  +coeff( 70)    *x23    *x43*x52
     8  +coeff( 71)    *x22*x31*x42*x53
      y_m12_0_11  =y_m12_0_11  
     9  +coeff( 72)        *x33*x42*x53
     1  +coeff( 73)    *x22    *x43*x53
     2  +coeff( 74)        *x32*x43*x53
     3  +coeff( 75)    *x23*x33*x42*x51
     4  +coeff( 76)    *x23*x32*x43*x51
     5  +coeff( 77)    *x22*x33*x42*x52
     6  +coeff( 78)    *x23*x33    *x53
     7  +coeff( 79)    *x21*x33*x42*x53
     8  +coeff( 80)    *x23*x33*x42*x52
      y_m12_0_11  =y_m12_0_11  
     9  +coeff( 81)    *x22*x33*x42*x53
     1  +coeff( 82)    *x22*x31        
     2  +coeff( 83)        *x31*x42    
     3  +coeff( 84)        *x31    *x52
     4  +coeff( 85)    *x21*x32*x41    
     5  +coeff( 86)        *x33    *x51
     6  +coeff( 87)        *x31    *x53
     7  +coeff( 88)        *x33*x42    
     8  +coeff( 89)        *x32*x43    
      y_m12_0_11  =y_m12_0_11  
     9  +coeff( 90)    *x23    *x41*x51
c
      return
      end
      function p_m12_0_11  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)   !,xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49995E-02,-0.79804E-01,-0.99988E-01,-0.35144E-01,-0.21223E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49981E-02, 0.61732E-01, 0.99988E-01, 0.35144E-01, 0.42498E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.21217234E-01,-0.64921446E-01,-0.89027663E-03,-0.56505157E-02,
     +  0.15275474E-01, 0.53110074E-01,-0.26885720E-02, 0.25902092E-03,
     + -0.11783465E-01,-0.18686557E-03,-0.25847827E-02,-0.21224648E-02,
     + -0.45458754E-02,-0.20256799E-01,-0.98129921E-02,-0.37170567E-01,
     +  0.38513294E-05, 0.13091408E-02, 0.47834284E-04, 0.17228113E-01,
     +  0.26468124E-01,-0.10320991E-01,-0.43801501E-01,-0.30462479E-02,
     +  0.10808182E-01,-0.73532700E-02,-0.31789411E-01, 0.57110139E-02,
     + -0.17361168E-02, 0.94580319E-03, 0.42405251E-01, 0.59151050E-03,
     +  0.47468700E-01, 0.26533973E-04, 0.30295021E-03,-0.69330824E-02,
     + -0.11760045E-02, 0.45302513E-03,-0.32362364E-01, 0.56069456E-01,
     +  0.81591889E-01,-0.70168888E-02,-0.84378914E-03,-0.26943605E-01,
     + -0.93682390E-03, 0.15922459E-02, 0.71764742E-02,-0.47460175E-03,
     +  0.49247714E-02,-0.96674594E-04,-0.35403369E-02, 0.14909821E-02,
     +  0.45680585E-02, 0.44563536E-01,-0.64745918E-02,-0.62388653E-03,
     +  0.42815294E-01, 0.31219607E-02,-0.66068661E-02, 0.17973440E-02,
     +  0.26849523E-01, 0.16929358E+00, 0.52475516E-04, 0.20957385E+00,
     + -0.49312720E-02,-0.30893078E-02,-0.58706146E-05,-0.15263114E-01,
     +  0.18933535E-01, 0.12022507E+00, 0.15660332E+00, 0.11415058E-02,
     + -0.20586085E-02, 0.14729339E-01, 0.46190310E-01, 0.49437325E-01,
     + -0.37444173E-02, 0.14115550E-01,-0.15224442E-03,-0.74667744E-02,
     +  0.15609006E-02, 0.18847838E-01, 0.10444856E+00, 0.47375308E-02,
     +  0.11911366E+00, 0.17655483E-01,-0.67503401E-02, 0.81022456E-01,
     +  0.10797419E+00,-0.49091864E-03,
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
      p_m12_0_11  =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21*x31        
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)        *x31    *x51
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)    *x22*x31        
     8  +coeff(  8)        *x33        
      p_m12_0_11  =p_m12_0_11  
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)        *x32*x41    
     2  +coeff( 11)        *x31*x42    
     3  +coeff( 12)            *x43    
     4  +coeff( 13)    *x21*x31    *x51
     5  +coeff( 14)    *x21    *x41*x51
     6  +coeff( 15)        *x31    *x52
     7  +coeff( 16)            *x41*x52
     8  +coeff( 17)*x11    *x31        
      p_m12_0_11  =p_m12_0_11  
     9  +coeff( 18)    *x21*x33        
     1  +coeff( 19)*x11        *x41    
     2  +coeff( 20)    *x21*x31*x42    
     3  +coeff( 21)    *x21    *x43    
     4  +coeff( 22)    *x22*x31    *x51
     5  +coeff( 23)    *x22    *x41*x51
     6  +coeff( 24)        *x32*x41*x51
     7  +coeff( 25)            *x43*x51
     8  +coeff( 26)    *x21*x31    *x52
      p_m12_0_11  =p_m12_0_11  
     9  +coeff( 27)    *x21    *x41*x52
     1  +coeff( 28)            *x41*x53
     2  +coeff( 29)    *x22*x33        
     3  +coeff( 30)*x11*x21    *x41    
     4  +coeff( 31)    *x22*x31*x42    
     5  +coeff( 32)        *x33*x42    
     6  +coeff( 33)    *x22    *x43    
     7  +coeff( 34)        *x32*x43    
     8  +coeff( 35)*x11    *x31    *x51
      p_m12_0_11  =p_m12_0_11  
     9  +coeff( 36)    *x23*x31    *x51
     1  +coeff( 37)    *x21*x33    *x51
     2  +coeff( 38)*x11        *x41*x51
     3  +coeff( 39)    *x23    *x41*x51
     4  +coeff( 40)    *x21*x31*x42*x51
     5  +coeff( 41)    *x21    *x43*x51
     6  +coeff( 42)    *x22*x31    *x52
     7  +coeff( 43)        *x33    *x52
     8  +coeff( 44)    *x22    *x41*x52
      p_m12_0_11  =p_m12_0_11  
     9  +coeff( 45)        *x32*x41*x52
     1  +coeff( 46)        *x31*x42*x52
     2  +coeff( 47)            *x43*x52
     3  +coeff( 48)    *x21*x31    *x53
     4  +coeff( 49)    *x21    *x41*x53
     5  +coeff( 50)*x11*x22*x31        
     6  +coeff( 51)    *x23*x33        
     7  +coeff( 52)*x11*x22    *x41    
     8  +coeff( 53)    *x23*x32*x41    
      p_m12_0_11  =p_m12_0_11  
     9  +coeff( 54)    *x23*x31*x42    
     1  +coeff( 55)    *x21*x33*x42    
     2  +coeff( 56)*x11        *x43    
     3  +coeff( 57)    *x23    *x43    
     4  +coeff( 58)    *x21*x32*x43    
     5  +coeff( 59)    *x22*x33    *x51
     6  +coeff( 60)*x11*x21    *x41*x51
     7  +coeff( 61)    *x22*x32*x41*x51
     8  +coeff( 62)    *x22*x31*x42*x51
      p_m12_0_11  =p_m12_0_11  
     9  +coeff( 63)        *x33*x42*x51
     1  +coeff( 64)    *x22    *x43*x51
     2  +coeff( 65)    *x23*x31    *x52
     3  +coeff( 66)    *x21*x33    *x52
     4  +coeff( 67)*x11        *x41*x52
     5  +coeff( 68)    *x23    *x41*x52
     6  +coeff( 69)    *x21*x32*x41*x52
     7  +coeff( 70)    *x21*x31*x42*x52
     8  +coeff( 71)    *x21    *x43*x52
      p_m12_0_11  =p_m12_0_11  
     9  +coeff( 72)        *x33    *x53
     1  +coeff( 73)    *x22    *x41*x53
     2  +coeff( 74)        *x32*x41*x53
     3  +coeff( 75)        *x31*x42*x53
     4  +coeff( 76)            *x43*x53
     5  +coeff( 77)    *x22*x33*x42    
     6  +coeff( 78)    *x22*x32*x43    
     7  +coeff( 79)*x11*x22*x31    *x51
     8  +coeff( 80)    *x23*x33    *x51
      p_m12_0_11  =p_m12_0_11  
     9  +coeff( 81)*x11*x22    *x41*x51
     1  +coeff( 82)    *x23*x32*x41*x51
     2  +coeff( 83)    *x23*x31*x42*x51
     3  +coeff( 84)    *x21*x33*x42*x51
     4  +coeff( 85)    *x23    *x43*x51
     5  +coeff( 86)    *x21*x32*x43*x51
     6  +coeff( 87)    *x22*x33    *x52
     7  +coeff( 88)    *x22*x31*x42*x52
     8  +coeff( 89)    *x22    *x43*x52
      p_m12_0_11  =p_m12_0_11  
     9  +coeff( 90)*x11    *x31    *x53
c
      return
      end
      real function l_m12_0_11  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)   !,xmean(10)
      dimension coeff( 41)
      data ncoeff/ 40/
      data avdat/  0.1597054E-02/
      data xmin/
     1 -0.49995E-02,-0.79804E-01,-0.99988E-01,-0.35144E-01,-0.21223E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49981E-02, 0.61732E-01, 0.99988E-01, 0.35144E-01, 0.42498E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.18533411E-01,-0.20176940E+00,-0.12831977E+00,-0.67979604E-01,
     + -0.23786521E-01,-0.26702764E-01,-0.67269332E-02,-0.23972165E-01,
     + -0.94070267E-02, 0.21917744E-01,-0.12563294E-01, 0.27927954E-01,
     + -0.16208321E-01, 0.82551436E-02, 0.17731810E-01, 0.26553128E-01,
     + -0.71367142E-02,-0.11770256E-02, 0.13204712E-01, 0.39631374E-01,
     +  0.12126390E-01, 0.17372094E-01, 0.35460457E-01,-0.87985611E-02,
     +  0.38597070E-01, 0.56224499E-01,-0.50109439E-03,-0.14663636E-02,
     +  0.10263203E-02, 0.64913691E-02, 0.76674325E-02, 0.22200574E-02,
     +  0.41207708E-02,-0.40333918E-02,-0.39191982E-02,-0.84301119E-03,
     +  0.21461757E-01, 0.13334877E-01, 0.23121158E-01, 0.70610452E-02,
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
      l_m12_0_11  =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x21        *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)            *x42    
     7  +coeff(  7)                *x52
     8  +coeff(  8)    *x21        *x52
      l_m12_0_11  =l_m12_0_11  
     9  +coeff(  9)        *x31*x41    
     1  +coeff( 10)                *x53
     2  +coeff( 11)    *x22        *x52
     3  +coeff( 12)            *x42*x51
     4  +coeff( 13)    *x24            
     5  +coeff( 14)    *x22    *x42    
     6  +coeff( 15)    *x21*x31*x41*x51
     7  +coeff( 16)    *x21        *x53
     8  +coeff( 17)    *x24        *x51
      l_m12_0_11  =l_m12_0_11  
     9  +coeff( 18)        *x33*x41*x51
     1  +coeff( 19)    *x22        *x53
     2  +coeff( 20)    *x24*x31*x41    
     3  +coeff( 21)    *x21    *x42    
     4  +coeff( 22)        *x31*x41*x51
     5  +coeff( 23)    *x21    *x42*x51
     6  +coeff( 24)                *x54
     7  +coeff( 25)    *x22    *x42*x51
     8  +coeff( 26)    *x24    *x42    
      l_m12_0_11  =l_m12_0_11  
     9  +coeff( 27)*x11                
     1  +coeff( 28)        *x32        
     2  +coeff( 29)    *x21*x32        
     3  +coeff( 30)    *x21*x31*x41    
     4  +coeff( 31)    *x22        *x51
     5  +coeff( 32)        *x32    *x51
     6  +coeff( 33)    *x22*x32        
     7  +coeff( 34)    *x23        *x51
     8  +coeff( 35)            *x42*x52
      l_m12_0_11  =l_m12_0_11  
     9  +coeff( 36)*x11*x23            
     1  +coeff( 37)    *x22*x31*x41*x51
     2  +coeff( 38)    *x21*x31*x41*x52
     3  +coeff( 39)    *x21    *x42*x52
     4  +coeff( 40)            *x42*x53
c
      return
      end
c***********************************************************************************
c***********************************************************************************      
c Reverse functions for MAD with 35 degree bend angle and   
c 12 degree minimum scattering angle    
c                                              -JJL 12/7/04   
c   
c re usage:   
c There are 5 fortran functions:   
c   
c TXFIT_12d(x,m)   
c delta_12(x,m)   
c theta_12(x,m)   
c phi_12(x,m)   
c y00_12(x,m)   
   
c Calling convention is as follows. x is an array in which x(1->5) are xf, thetaf, yf, phif   
c and x0 (the detector coordinates in the transport coordinate system and the vertical beam   
c position at the target). For TXFIT m=1 for all others m=5.   
c   
c Procedure:   
c Convert theta to a new orthogonalized theta' as follows:   
c   
c theta'=theta-TXFIT_12(x,1)   
c   
c Put the new theta' into your x array (or x(2)=x(2)-TXFIT_12(x,1) should work)   
c   
c Then I think it's self-explanitory:   
c   
c delta0 = delta_12(x,5)   
c theta0= theta_12(x,5)   
c phi0=phi_12(x,5)   
c y0=y00_12(x,5)
      function     txfit_12d   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)
      dimension coeff(  5)
      data ncoeff/  4/
      data avdat/ -0.2169641E-01/
      data xmin/
     1 -0.94556E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.51091E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.37064292E-01, 0.19944827E+00,-0.86163562E-02, 0.75834477E-03,
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
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x13 = x12*x1
c
c                  function
c
          txfit_12d   =avdat
     1  +coeff(  1)    
     2  +coeff(  2)*x11
     3  +coeff(  3)*x12
     4  +coeff(  4)*x13
c
      return
      end
      function    delta_12 (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 96)
      data ncoeff/ 95/
      data avdat/  0.1045559E+00/
      data xmin/
     1 -0.56442E+00,-0.61299E-01,-0.51957E+00,-0.89836E-01,-0.49995E-02,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.81280E+00, 0.10645E+00, 0.51957E+00, 0.89836E-01, 0.49981E-02,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.24586268E-01, 0.33275667E+00, 0.13156047E+00, 0.12517072E+00,
     +  0.26518315E+00, 0.69176130E-01, 0.24281401E-01, 0.62989064E-01,
     +  0.59394259E-01, 0.26033872E-02, 0.51660787E-01, 0.15277949E+00,
     +  0.13470389E+00, 0.40945482E-01,-0.54485973E-01, 0.88406093E-02,
     + -0.15463565E+00,-0.89213654E-01, 0.39638869E-01, 0.15665422E-02,
     +  0.15014862E-02,-0.21842318E-01, 0.11139997E-01, 0.56103524E-01,
     +  0.15335392E-01, 0.34196216E-02, 0.61987907E-01,-0.38702003E-03,
     + -0.10811804E+00, 0.31151703E-01,-0.36849687E+00,-0.97785413E-01,
     +  0.11400052E-01,-0.11068811E+00,-0.25410170E-01,-0.10664346E+00,
     + -0.10938760E-01,-0.12062149E+00,-0.21989513E-01,-0.15334544E-01,
     +  0.34719855E-01, 0.18525918E+00, 0.33938250E+00, 0.27060273E+00,
     +  0.84247336E-01,-0.11057190E-01,-0.15676405E+00, 0.18130745E+00,
     +  0.31219419E-01,-0.11399706E+00, 0.16210362E+00,-0.81142433E-01,
     + -0.13607357E+00,-0.21597424E-02, 0.94299763E-03,-0.95851463E-03,
     +  0.14185937E-02,-0.27995312E+00, 0.84809154E-01,-0.11014417E+00,
     + -0.17460691E-01,-0.26247771E-01,-0.98115928E-03,-0.34344971E-01,
     +  0.54868206E-01,-0.19493084E+00, 0.15917148E+00, 0.94514489E-01,
     +  0.36450726E+00,-0.20341215E+00,-0.61310176E-02,-0.51045464E-02,
     +  0.32895394E-01,-0.36246464E+00, 0.33149701E+00, 0.69789207E+00,
     +  0.30678049E+00,-0.60312585E-02, 0.53077757E-01,-0.26165012E+00,
     + -0.89516856E-01,-0.77802032E-01,-0.59587955E+00, 0.68758838E-02,
     + -0.52525974E-02,-0.18267448E+00,-0.73022902E-01,-0.16975284E+00,
     + -0.11510786E+00,-0.25261549E-01, 0.12535074E+00,-0.13764973E+00,
     + -0.16547281E+00,-0.34962371E+00,-0.14314376E+00,
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
         delta_12 =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)*x11                
     3  +coeff(  3)    *x21            
     4  +coeff(  4)*x12                
     5  +coeff(  5)*x11*x21            
     6  +coeff(  6)    *x22            
     7  +coeff(  7)        *x32        
     8  +coeff(  8)        *x31*x41    
         delta_12 =   delta_12 
     9  +coeff(  9)            *x42    
     1  +coeff( 10)                *x51
     2  +coeff( 11)*x13                
     3  +coeff( 12)*x12*x21            
     4  +coeff( 13)*x11*x22            
     5  +coeff( 14)    *x23            
     6  +coeff( 15)*x11    *x32        
     7  +coeff( 16)    *x21*x32        
     8  +coeff( 17)*x11    *x31*x41    
         delta_12 =   delta_12 
     9  +coeff( 18)*x11        *x42    
     1  +coeff( 19)    *x21    *x42    
     2  +coeff( 20)*x11            *x51
     3  +coeff( 21)    *x21        *x51
     4  +coeff( 22)*x14                
     5  +coeff( 23)*x13*x21            
     6  +coeff( 24)*x12*x22            
     7  +coeff( 25)*x11*x23            
     8  +coeff( 26)    *x24            
         delta_12 =   delta_12 
     9  +coeff( 27)*x12    *x32        
     1  +coeff( 28)*x11*x21*x32        
     2  +coeff( 29)    *x22*x32        
     3  +coeff( 30)*x12    *x31*x41    
     4  +coeff( 31)    *x22*x31*x41    
     5  +coeff( 32)*x12        *x42    
     6  +coeff( 33)*x11*x21    *x42    
     7  +coeff( 34)    *x22    *x42    
     8  +coeff( 35)*x15                
         delta_12 =   delta_12 
     9  +coeff( 36)*x14*x21            
     1  +coeff( 37)*x13*x22            
     2  +coeff( 38)*x12*x23            
     3  +coeff( 39)*x11*x24            
     4  +coeff( 40)    *x25            
     5  +coeff( 41)        *x34        
     6  +coeff( 42)        *x33*x41    
     7  +coeff( 43)        *x32*x42    
     8  +coeff( 44)        *x31*x43    
         delta_12 =   delta_12 
     9  +coeff( 45)            *x44    
     1  +coeff( 46)*x13    *x32        
     2  +coeff( 47)*x12*x21*x32        
     3  +coeff( 48)*x11*x22*x32        
     4  +coeff( 49)    *x23*x32        
     5  +coeff( 50)*x12*x21*x31*x41    
     6  +coeff( 51)*x11*x22*x31*x41    
     7  +coeff( 52)    *x23*x31*x41    
     8  +coeff( 53)    *x23    *x42    
         delta_12 =   delta_12 
     9  +coeff( 54)*x13            *x51
     1  +coeff( 55)*x12*x21        *x51
     2  +coeff( 56)    *x23        *x51
     3  +coeff( 57)*x15*x21            
     4  +coeff( 58)*x14*x22            
     5  +coeff( 59)*x13*x23            
     6  +coeff( 60)*x12*x24            
     7  +coeff( 61)    *x21*x34        
     8  +coeff( 62)    *x21*x33*x41    
         delta_12 =   delta_12 
     9  +coeff( 63)    *x21*x32*x42    
     1  +coeff( 64)*x14    *x32        
     2  +coeff( 65)    *x24*x32        
     3  +coeff( 66)*x14    *x31*x41    
     4  +coeff( 67)*x12*x22*x31*x41    
     5  +coeff( 68)    *x24*x31*x41    
     6  +coeff( 69)*x12*x22    *x42    
     7  +coeff( 70)*x11*x23    *x42    
     8  +coeff( 71)    *x24    *x42    
         delta_12 =   delta_12 
     9  +coeff( 72)*x13*x21        *x51
     1  +coeff( 73)*x15*x22            
     2  +coeff( 74)*x14*x23            
     3  +coeff( 75)    *x22*x32*x42    
     4  +coeff( 76)    *x22*x31*x43    
     5  +coeff( 77)    *x22    *x44    
     6  +coeff( 78)*x11*x21*x31*x41*x51
     7  +coeff( 79)*x14*x21*x32        
     8  +coeff( 80)*x13*x22*x32        
         delta_12 =   delta_12 
     9  +coeff( 81)*x12*x23*x32        
     1  +coeff( 82)*x15    *x31*x41    
     2  +coeff( 83)*x12*x23*x31*x41    
     3  +coeff( 84)*x11*x24    *x42    
     4  +coeff( 85)*x13*x22        *x51
     5  +coeff( 86)*x14*x24            
     6  +coeff( 87)*x13*x25            
     7  +coeff( 88)*x12*x21*x32*x42    
     8  +coeff( 89)    *x23*x32*x42    
         delta_12 =   delta_12 
     9  +coeff( 90)*x13    *x31*x43    
     1  +coeff( 91)*x11*x22    *x44    
     2  +coeff( 92)*x15*x21*x32        
     3  +coeff( 93)*x13*x23*x32        
     4  +coeff( 94)*x12*x24*x31*x41    
     5  +coeff( 95)*x11*x25*x31*x41    
c
      return
      end
      function    theta_12 (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 88)
      data ncoeff/ 87/
      data avdat/ -0.1949095E-01/
      data xmin/
     1 -0.56442E+00,-0.61299E-01,-0.51957E+00,-0.89836E-01,-0.49995E-02,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.81280E+00, 0.10645E+00, 0.51957E+00, 0.89836E-01, 0.49981E-02,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.23861702E-02, 0.49947593E-02,-0.76850615E-01,-0.54880008E-02,
     +  0.29143855E-01, 0.11518307E-01, 0.81482558E-02, 0.20488529E-01,
     +  0.18362874E-01,-0.24638412E-03, 0.45473393E-03,-0.15176488E-01,
     +  0.95629562E-02, 0.33070255E-03,-0.13556735E-01,-0.97072115E-02,
     + -0.14666531E-01,-0.21186268E-01,-0.81761340E-02,-0.19517954E-03,
     +  0.24621696E-02,-0.61357906E-02,-0.62964605E-02,-0.54896581E-02,
     + -0.15387335E-02,-0.80755614E-02,-0.38562187E-02,-0.29452166E-01,
     + -0.21857761E-02,-0.21407073E-02,-0.50718978E-03,-0.49148157E-03,
     + -0.32866795E-01, 0.82487054E-02,-0.32521486E-02, 0.11786449E-02,
     + -0.31281670E-02,-0.45581548E-02, 0.14145626E-01, 0.22989286E-01,
     +  0.66966265E-02, 0.21117464E-01, 0.25679085E-01, 0.24015648E-01,
     +  0.24964865E-02, 0.16151682E-01,-0.18711239E-01,-0.78508276E-02,
     +  0.81252521E-02, 0.39427322E-02,-0.16175896E-01,-0.14178577E-03,
     + -0.16254855E-02, 0.45330040E-02, 0.54369103E-01, 0.19282395E-01,
     +  0.15084945E-01, 0.38963608E-04, 0.54203026E-01, 0.54661229E-01,
     +  0.22946909E-01, 0.71976549E-03, 0.17757748E-02, 0.10564500E+00,
     + -0.41456591E-01, 0.75280339E-01,-0.91677625E-02, 0.69127558E-02,
     + -0.42695016E-01,-0.83999813E-01,-0.69443695E-03, 0.86770475E-01,
     + -0.66304244E-02, 0.30526433E-02, 0.23162854E+00,-0.28378705E-01,
     +  0.49740586E-01, 0.18806273E+00,-0.27878698E-01, 0.64504594E-01,
     + -0.17649952E-01,-0.16553177E-01,-0.54108329E-01, 0.15433392E+00,
     + -0.11277263E+00, 0.10738762E+00, 0.16928419E-02,
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
         theta_12 =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)*x11                
     3  +coeff(  3)    *x21            
     4  +coeff(  4)*x12                
     5  +coeff(  5)*x11*x21            
     6  +coeff(  6)    *x22            
     7  +coeff(  7)        *x32        
     8  +coeff(  8)        *x31*x41    
         theta_12 =   theta_12 
     9  +coeff(  9)            *x42    
     1  +coeff( 10)                *x51
     2  +coeff( 11)*x13                
     3  +coeff( 12)*x12*x21            
     4  +coeff( 13)*x11*x22            
     5  +coeff( 14)    *x23            
     6  +coeff( 15)*x11    *x32        
     7  +coeff( 16)    *x21*x32        
     8  +coeff( 17)*x11    *x31*x41    
         theta_12 =   theta_12 
     9  +coeff( 18)    *x21*x31*x41    
     1  +coeff( 19)    *x21    *x42    
     2  +coeff( 20)*x11            *x51
     3  +coeff( 21)*x14                
     4  +coeff( 22)*x13*x21            
     5  +coeff( 23)*x12*x22            
     6  +coeff( 24)*x11*x23            
     7  +coeff( 25)*x12    *x32        
     8  +coeff( 26)    *x22*x32        
         theta_12 =   theta_12 
     9  +coeff( 27)*x12    *x31*x41    
     1  +coeff( 28)    *x22*x31*x41    
     2  +coeff( 29)*x12        *x42    
     3  +coeff( 30)    *x22    *x42    
     4  +coeff( 31)*x11*x21        *x51
     5  +coeff( 32)*x14*x21            
     6  +coeff( 33)*x13*x22            
     7  +coeff( 34)*x12*x23            
     8  +coeff( 35)*x11*x24            
         theta_12 =   theta_12 
     9  +coeff( 36)    *x25            
     1  +coeff( 37)        *x34        
     2  +coeff( 38)        *x33*x41    
     3  +coeff( 39)*x12*x21*x32        
     4  +coeff( 40)*x11*x22*x32        
     5  +coeff( 41)    *x23*x32        
     6  +coeff( 42)*x12*x21*x31*x41    
     7  +coeff( 43)*x11*x22*x31*x41    
     8  +coeff( 44)*x12*x21    *x42    
         theta_12 =   theta_12 
     9  +coeff( 45)*x12*x21        *x51
     1  +coeff( 46)*x15*x21            
     2  +coeff( 47)*x14*x22            
     3  +coeff( 48)*x13*x23            
     4  +coeff( 49)    *x21*x32*x42    
     5  +coeff( 50)*x11        *x44    
     6  +coeff( 51)    *x21    *x44    
     7  +coeff( 52)*x11        *x42*x51
     8  +coeff( 53)    *x21    *x42*x51
         theta_12 =   theta_12 
     9  +coeff( 54)*x14    *x32        
     1  +coeff( 55)*x13*x21*x32        
     2  +coeff( 56)*x12*x22*x32        
     3  +coeff( 57)*x11*x23*x32        
     4  +coeff( 58)*x14    *x31*x41    
     5  +coeff( 59)*x13*x21*x31*x41    
     6  +coeff( 60)*x12*x22*x31*x41    
     7  +coeff( 61)*x11*x23*x31*x41    
     8  +coeff( 62)*x14            *x51
         theta_12 =   theta_12 
     9  +coeff( 63)*x12*x22        *x51
     1  +coeff( 64)*x15*x22            
     2  +coeff( 65)*x14*x23            
     3  +coeff( 66)*x13*x24            
     4  +coeff( 67)*x12*x25            
     5  +coeff( 68)*x12    *x34        
     6  +coeff( 69)    *x22*x33*x41    
     7  +coeff( 70)    *x22*x32*x42    
     8  +coeff( 71)*x11*x21*x31*x41*x51
         theta_12 =   theta_12 
     9  +coeff( 72)*x13*x22*x32        
     1  +coeff( 73)*x12*x23*x32        
     2  +coeff( 74)*x11*x24*x32        
     3  +coeff( 75)*x13*x22*x31*x41    
     4  +coeff( 76)*x12*x23*x31*x41    
     5  +coeff( 77)*x11*x24*x31*x41    
     6  +coeff( 78)*x15*x23            
     7  +coeff( 79)*x14*x24            
     8  +coeff( 80)*x13*x25            
         theta_12 =   theta_12 
     9  +coeff( 81)*x11*x22*x34        
     1  +coeff( 82)    *x23*x32*x42    
     2  +coeff( 83)*x12*x24*x32        
     3  +coeff( 84)*x13*x23*x31*x41    
     4  +coeff( 85)*x12*x24*x31*x41    
     5  +coeff( 86)*x15*x24            
     6  +coeff( 87)    *x25*x32    *x51
c
      return
      end
      function     y00_12  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 58)
      data ncoeff/ 57/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.56442E+00,-0.61299E-01,-0.51957E+00,-0.89836E-01,-0.49995E-02,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.81280E+00, 0.10645E+00, 0.51957E+00, 0.89836E-01, 0.49981E-02,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.10162787E+01,-0.17358043E+01, 0.91096514E+00, 0.40611893E+00,
     + -0.97879693E-01, 0.94740845E-01,-0.32691962E+00, 0.36146399E+00,
     +  0.12600388E+00,-0.33428073E+00, 0.91888577E+00, 0.11264547E+01,
     + -0.36496103E-01,-0.29634383E+00, 0.11709670E+01,-0.88792913E-01,
     + -0.11312437E+00, 0.14116300E+01, 0.36023841E+01, 0.92666912E+00,
     +  0.60180652E+00,-0.23728045E-01, 0.77989566E+00,-0.62857580E+00,
     + -0.96877702E-01, 0.67287651E-02, 0.40692153E+01,-0.79342216E-01,
     + -0.10796915E+01, 0.19241075E+00, 0.41805158E+01, 0.51975787E+00,
     + -0.14094219E+00,-0.21525349E+00,-0.52222586E+00,-0.72348765E-02,
     +  0.20832242E+00, 0.22318900E+01, 0.32411590E-01, 0.14296093E+01,
     +  0.29416997E+01,-0.36808440E+00,-0.85581811E-02, 0.12733635E+00,
     +  0.21771505E+00, 0.25092273E-02,-0.16232964E-01, 0.18439721E-01,
     +  0.16106263E+00, 0.32601702E+01,-0.57054205E+01,-0.19498261E+01,
     + -0.11879611E+02,-0.70895761E+00,-0.61978493E+01,-0.97314566E+00,
     +  0.12332698E+01,
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
      x51 = x5
c
c                  function
c
          y00_12  =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)*x11    *x31        
     4  +coeff(  4)    *x21*x31        
     5  +coeff(  5)*x11        *x41    
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)*x12    *x31        
     8  +coeff(  8)*x11*x21*x31        
          y00_12  =    y00_12  
     9  +coeff(  9)*x11*x21    *x41    
     1  +coeff( 10)    *x21*x33        
     2  +coeff( 11)*x11    *x32*x41    
     3  +coeff( 12)    *x22*x33        
     4  +coeff( 13)    *x23*x31        
     5  +coeff( 14)*x12*x22*x31        
     6  +coeff( 15)*x11*x23*x33        
     7  +coeff( 16)        *x33        
     8  +coeff( 17)        *x32*x41    
          y00_12  =    y00_12  
     9  +coeff( 18)    *x21*x31*x42    
     1  +coeff( 19)*x11        *x43    
     2  +coeff( 20)    *x21    *x43    
     3  +coeff( 21)*x12    *x33        
     4  +coeff( 22)*x11    *x35        
     5  +coeff( 23)*x12*x21*x33        
     6  +coeff( 24)*x12*x23*x31        
     7  +coeff( 25)        *x31*x42    
     8  +coeff( 26)        *x31    *x51
          y00_12  =    y00_12  
     9  +coeff( 27)*x11    *x31*x42    
     1  +coeff( 28)*x13    *x31        
     2  +coeff( 29)    *x23*x33        
     3  +coeff( 30)*x14    *x33        
     4  +coeff( 31)*x11*x23*x32*x41    
     5  +coeff( 32)*x14*x23*x31        
     6  +coeff( 33)*x13*x24*x31        
     7  +coeff( 34)*x11*x22*x31        
     8  +coeff( 35)    *x23    *x41    
          y00_12  =    y00_12  
     9  +coeff( 36)    *x21*x31    *x51
     1  +coeff( 37)*x12    *x32*x41    
     2  +coeff( 38)    *x22*x32*x41    
     3  +coeff( 39)*x11*x21*x31*x42    
     4  +coeff( 40)    *x22*x31*x42    
     5  +coeff( 41)*x12        *x43    
     6  +coeff( 42)*x11*x21    *x43    
     7  +coeff( 43)        *x33    *x51
     8  +coeff( 44)*x13*x21*x31        
          y00_12  =    y00_12  
     9  +coeff( 45)*x11*x23*x31        
     1  +coeff( 46)    *x24*x31        
     2  +coeff( 47)*x11*x21*x31    *x51
     3  +coeff( 48)    *x22*x31    *x51
     4  +coeff( 49)*x13    *x33        
     5  +coeff( 50)*x12*x21*x32*x41    
     6  +coeff( 51)    *x23*x32*x41    
     7  +coeff( 52)*x13    *x31*x42    
     8  +coeff( 53)    *x23*x31*x42    
          y00_12  =    y00_12  
     9  +coeff( 54)*x12*x21    *x43    
     1  +coeff( 55)    *x23    *x43    
     2  +coeff( 56)*x13*x21*x33        
     3  +coeff( 57)*x14*x24*x31*x42    
c
      return
      end
      function     phi_12  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 53)
      data ncoeff/ 52/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.56442E+00,-0.61299E-01,-0.51957E+00,-0.89836E-01,-0.49995E-02,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.81280E+00, 0.10645E+00, 0.51957E+00, 0.89836E-01, 0.49981E-02,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.11390240E+00, 0.14874293E+00,-0.83443716E-01,-0.35192788E-01,
     +  0.11123690E-02,-0.83252816E-02, 0.48874788E-01, 0.54901987E-01,
     +  0.37320618E-01,-0.41921142E-01, 0.24584227E-02,-0.26820971E-01,
     +  0.26954839E-02,-0.51758446E-01,-0.24943782E+00, 0.49331099E-01,
     +  0.11886330E+00,-0.20919034E-01, 0.78152589E-01, 0.17844118E-01,
     + -0.64815185E-02,-0.67603675E-03,-0.26516430E-03, 0.58701267E-02,
     + -0.26369607E+00, 0.63009560E-01, 0.28105574E-01, 0.24539918E-01,
     +  0.14309634E-01, 0.18829811E-01,-0.35199702E+00,-0.23454648E-01,
     +  0.55851322E-02,-0.48750412E-03,-0.63773222E-01,-0.90385295E-01,
     +  0.45590144E-01,-0.17168494E-01,-0.19447497E-02,-0.38128269E+00,
     + -0.64957693E-01, 0.40769065E-02,-0.24281444E-02, 0.67877921E-03,
     +  0.91292540E-03,-0.14203310E-02,-0.62406957E-01, 0.33628768E+00,
     + -0.76849908E-01, 0.82115769E+00, 0.43210882E+00,-0.10671967E+00,
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
      x51 = x5
c
c                  function
c
          phi_12  =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)*x11    *x31        
     4  +coeff(  4)    *x21*x31        
     5  +coeff(  5)*x11        *x41    
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)        *x31*x42    
     8  +coeff(  8)            *x43    
          phi_12  =    phi_12  
     9  +coeff(  9)*x12    *x31        
     1  +coeff( 10)*x11*x21*x31        
     2  +coeff( 11)    *x22*x31        
     3  +coeff( 12)*x11*x21    *x41    
     4  +coeff( 13)    *x22    *x41    
     5  +coeff( 14)    *x21*x33        
     6  +coeff( 15)    *x21*x32*x41    
     7  +coeff( 16)    *x23    *x41    
     8  +coeff( 17)*x11*x21*x33        
          phi_12  =    phi_12  
     9  +coeff( 18)*x14    *x35        
     1  +coeff( 19)    *x22*x32*x41    
     2  +coeff( 20)*x11*x23    *x41    
     3  +coeff( 21)*x12*x23*x31        
     4  +coeff( 22)        *x31    *x51
     5  +coeff( 23)*x11    *x33        
     6  +coeff( 24)*x11    *x31*x42    
     7  +coeff( 25)    *x21*x31*x42    
     8  +coeff( 26)*x11        *x43    
          phi_12  =    phi_12  
     9  +coeff( 27)*x12*x21*x31        
     1  +coeff( 28)*x11*x22*x31        
     2  +coeff( 29)*x12*x21    *x41    
     3  +coeff( 30)    *x22*x33        
     4  +coeff( 31)*x11*x21*x31*x42    
     5  +coeff( 32)*x13*x21*x31        
     6  +coeff( 33)*x11*x23*x31        
     7  +coeff( 34)    *x21*x35        
     8  +coeff( 35)*x13    *x33        
          phi_12  =    phi_12  
     9  +coeff( 36)*x11*x22*x33        
     1  +coeff( 37)    *x23*x33        
     2  +coeff( 38)*x14*x21*x31        
     3  +coeff( 39)    *x24*x33        
     4  +coeff( 40)*x11*x23*x32*x41    
     5  +coeff( 41)*x11*x23*x31*x42    
     6  +coeff( 42)*x12        *x41    
     7  +coeff( 43)*x13        *x41    
     8  +coeff( 44)    *x21*x31    *x51
          phi_12  =    phi_12  
     9  +coeff( 45)        *x33    *x51
     1  +coeff( 46)    *x22*x31    *x51
     2  +coeff( 47)*x13    *x32*x41    
     3  +coeff( 48)    *x23*x32*x41    
     4  +coeff( 49)*x12*x21*x31*x42    
     5  +coeff( 50)    *x23*x31*x42    
     6  +coeff( 51)    *x23    *x43    
     7  +coeff( 52)*x11*x23*x33        
c
      return
      end
