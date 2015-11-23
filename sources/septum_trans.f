C Forward transfer functions for left hrs with septum based on Leftseptum_dir.dat
c                     -JJL 4/4/04


c typical call: answer = function(x,5)
c INPUTS: x = 5 or more element array 
c              x(1)=x0  (meters)
c              x(2)=theta0 (really tan(theta0))
c              x(3)=y0   (meters)
c              x(4)=phi0 (really tan(phi0))
c              x(5)=delta (fractional value NOT percent)
c         M=5
c
c OUTPUT: units are the same as inputs
c 
c NOMENCLATURE: function name = prefix + _sr_ +suffix
c           prefixes:     x means xfinal
c                         t means thetafinal
c                         y means yfinal
c                         p means phifinal
c                         l means path length
c     
c           suffixes:     fp means target to focus
c                         q1ex means target to Q1 exit
c                         dent means target to Dipole entrance
c                         dext means target to dipole exit
c                         q3en means target to Q3 entrance
c                         q3ex means target to Q3 exit
c                         ep3  means septum entrance
c                         ep4  means one quarter the way into the septum
c                         ep5  means halfway through the septum
c                         ep6  means three quarters the way through the septum
c                         ep7  means septum exit
c
c          _sl_ is for left hrs with septum
c
c APERTURES:
c    Coordinate systems are different in the spectrometer with septum model 
c    compared to the spectrometer without septum. So some apertures even in 
c    the body of the spectrometer appear to be different. Numbers given here 
c    supercede any numbers given in regard to transfer functions for the 
c    spectrometers alone.
c
c     ep3: -0.1486233 < x < -0.08869672
c          -0.110 < y < 0.110
c     ep4: -0.1792231 < x < -0.1089169
c          -0.110 < y < 0.110
c     ep5: -0.2209211 < x < -0.1353789
c          -0.110 < y < 0.110
c     ep6: -0.2763536 < x < -0.1697464
c          -0.110 < y < 0.110
c     ep7: -0.3485396 < x < -0.2156404
c          -0.110 < y < 0.110
c     q1ex is a circle of radius 0.1492 m
c     dent is a trapazoid:
c                                   -5.22008 < x < -4.98099
c             -(-0.192436784*x -0.192436784) < y < -0.192436784*x -0.192436784  
c     dext is also a trapazoid: 
c                                   -0.46188 < x < 0.46188
c                   -(-0.01610808*x + 0.125) < y < -0.01610808*x + 0.125
c     q3en is a circle of radius 0.3 m
c     q3ex is a circle of radius 0.3 m
c

      function x_sl_ep3    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff( 11)
      data ncoeff/ 10/
      data avdat/ -0.1159293E+00/
      data xmin/
     1 -0.19937E-02,-0.54847E-01,-0.19935E-01,-0.29994E-01,-0.49883E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19971E-02, 0.54952E-01, 0.19987E-01, 0.29928E-01, 0.49858E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.20842899E-02,-0.19864492E-01,-0.32466602E-01,-0.48664642E-04,
     +  0.75161574E-05, 0.16019378E-05,-0.31904840E-05,-0.90302291E-06,
     + -0.10595068E-04, 0.11970922E-05,
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
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
c
c                  function
c
      x_sl_ep3    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x22    *x41    
     5  +coeff(  5)    *x22            
     6  +coeff(  6)                *x51
     7  +coeff(  7)            *x42    
     8  +coeff(  8)        *x31*x42    
      x_sl_ep3    =x_sl_ep3    
     9  +coeff(  9)            *x43    
     1  +coeff( 10)        *x32        
c
      return
      end
      function t_sl_ep3    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff(  6)
      data ncoeff/  5/
      data avdat/ -0.1067509E+00/
      data xmin/
     1 -0.19937E-02,-0.54847E-01,-0.19935E-01,-0.29994E-01,-0.49883E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19971E-02, 0.54952E-01, 0.19987E-01, 0.29928E-01, 0.49858E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.10833655E-02,-0.15684851E-04,-0.29986490E-01, 0.38844041E-03,
     +  0.46003821E-04,
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
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x41 = x4
      x51 = x5
c
c                  function
c
      t_sl_ep3    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)                *x51
c
      return
      end
      function y_sl_ep3    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff(  4)
      data ncoeff/  3/
      data avdat/ -0.1313174E-02/
      data xmin/
     1 -0.19937E-02,-0.54847E-01,-0.19935E-01,-0.29994E-01,-0.49883E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19971E-02, 0.54952E-01, 0.19987E-01, 0.29928E-01, 0.49858E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.13721157E-02, 0.59822123E-01, 0.19945537E-02,
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
c          set up monomials   functions
      x11 = x1
      x21 = x2
c
c                  function
c
      y_sl_ep3    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
c
      return
      end
      function p_sl_ep3    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff(  4)
      data ncoeff/  3/
      data avdat/ -0.1186335E-02/
      data xmin/
     1 -0.19937E-02,-0.54847E-01,-0.19935E-01,-0.29994E-01,-0.49883E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19971E-02, 0.54952E-01, 0.19987E-01, 0.29928E-01, 0.49858E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.12386889E-02, 0.54867808E-01, 0.12650639E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
c          set up monomials   functions
      x21 = x2
      x41 = x4
c
c                  function
c
      p_sl_ep3    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
c
      return
      end
      function l_sl_ep3    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff( 15)
      data ncoeff/ 14/
      data avdat/  0.1089571E+01/
      data xmin/
     1 -0.19937E-02,-0.54847E-01,-0.19935E-01,-0.29994E-01,-0.49883E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19971E-02, 0.54952E-01, 0.19987E-01, 0.29928E-01, 0.49858E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.67956746E-03, 0.19444669E-03,-0.79449045E-07,-0.35980692E-08,
     +  0.16431817E-02,-0.10784547E-08, 0.50341379E-06, 0.40624084E-06,
     +  0.49014221E-03,-0.11419254E-05,-0.45034591E-06, 0.67545096E-06,
     +  0.10892550E-05, 0.11918529E-03,
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
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
c
c                  function
c
      l_sl_ep3    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)        *x32        
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)        *x31*x41    
      l_sl_ep3    =l_sl_ep3    
     9  +coeff(  9)            *x42    
     1  +coeff( 10)    *x22*x31        
     2  +coeff( 11)        *x33        
     3  +coeff( 12)*x12    *x31        
     4  +coeff( 13)    *x22*x33        
     5  +coeff( 14)        *x31        
c
      return
      end
      function x_sl_ep4    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.1407679E+00/
      data xmin/
     1 -0.19937E-02,-0.54847E-01,-0.19935E-01,-0.29994E-01,-0.49883E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19971E-02, 0.54952E-01, 0.19987E-01, 0.29928E-01, 0.49858E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.23048890E-02,-0.19808570E-01,-0.38886923E-01, 0.86243446E-04,
     +  0.99380180E-04,-0.62419494E-04, 0.11646053E-03, 0.25984780E-04,
     +  0.98425699E-04,-0.11290584E-03, 0.36705635E-04, 0.14792397E-04,
     +  0.83239713E-04,-0.40876844E-05,-0.10424437E-05, 0.46539740E-05,
     + -0.18002879E-05, 0.44225933E-04,-0.73397273E-05,-0.36079349E-04,
     + -0.14804231E-03,-0.57919686E-04,-0.31065585E-04,-0.26219173E-04,
     + -0.29425868E-04, 0.61970240E-04, 0.12600183E-04, 0.32080316E-04,
     +  0.15864221E-04,-0.22726072E-05,-0.27196664E-04,-0.13712169E-04,
     + -0.20181758E-04,-0.80934831E-06,-0.18357950E-05, 0.76169031E-05,
     +  0.25296147E-05, 0.33598000E-04, 0.54135862E-05,-0.31302447E-04,
     +  0.67250988E-06,-0.13564144E-03, 0.72464445E-05, 0.74099726E-05,
     +  0.53812760E-05,-0.23936907E-05,-0.28978715E-04,-0.33200413E-04,
     + -0.22608658E-06,-0.34770765E-05,-0.12329814E-04,-0.15536189E-04,
     + -0.21788789E-03,-0.54696857E-05,-0.24093165E-05,-0.81475889E-06,
     +  0.23831453E-05,-0.77292070E-05, 0.25877382E-04,-0.67455994E-05,
     +  0.29991870E-04,-0.14235714E-03,-0.10665875E-03, 0.31033753E-04,
     +  0.21422535E-04,-0.23231558E-04,-0.16674489E-04,-0.38218527E-05,
     + -0.34305507E-04, 0.18774157E-04, 0.98996370E-05, 0.49271148E-04,
     + -0.67299297E-05,-0.12112937E-04,-0.21946889E-05, 0.19997755E-04,
     + -0.35281540E-03,-0.36491503E-03,-0.39636416E-03,-0.47885542E-03,
     + -0.19331777E-03,-0.18706340E-03,-0.38303762E-04,-0.15574695E-04,
     +  0.65634763E-05,-0.60203388E-05, 0.67240162E-05, 0.11569564E-04,
     +  0.13287942E-04, 0.11369504E-04,
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
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      x_sl_ep4    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x22    *x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)        *x32        
      x_sl_ep4    =x_sl_ep4    
     9  +coeff(  9)        *x31*x41    
     1  +coeff( 10)    *x22    *x42    
     2  +coeff( 11)    *x24*x32        
     3  +coeff( 12)    *x24*x31*x41    
     4  +coeff( 13)    *x24    *x42    
     5  +coeff( 14)                *x52
     6  +coeff( 15)            *x43    
     7  +coeff( 16)*x11*x21            
     8  +coeff( 17)    *x21        *x51
      x_sl_ep4    =x_sl_ep4    
     9  +coeff( 18)        *x31*x42    
     1  +coeff( 19)    *x24            
     2  +coeff( 20)    *x22*x31*x41    
     3  +coeff( 21)        *x31*x43    
     4  +coeff( 22)            *x44    
     5  +coeff( 23)        *x34*x41    
     6  +coeff( 24)    *x22*x32    *x51
     7  +coeff( 25)    *x22*x34        
     8  +coeff( 26)        *x34*x42    
      x_sl_ep4    =x_sl_ep4    
     9  +coeff( 27)    *x23    *x41*x52
     1  +coeff( 28)    *x24    *x43    
     2  +coeff( 29)*x11*x22*x31    *x51
     3  +coeff( 30)*x11*x22*x33        
     4  +coeff( 31)*x11*x24*x31    *x51
     5  +coeff( 32)*x11*x24*x32    *x51
     6  +coeff( 33)*x11*x24*x31*x41*x51
     7  +coeff( 34)    *x21*x31        
     8  +coeff( 35)    *x21    *x41    
      x_sl_ep4    =x_sl_ep4    
     9  +coeff( 36)    *x22*x31        
     1  +coeff( 37)        *x33        
     2  +coeff( 38)        *x32*x41    
     3  +coeff( 39)    *x21    *x41*x51
     4  +coeff( 40)    *x22*x32        
     5  +coeff( 41)        *x34        
     6  +coeff( 42)        *x32*x42    
     7  +coeff( 43)    *x23        *x51
     8  +coeff( 44)    *x22*x31    *x51
      x_sl_ep4    =x_sl_ep4    
     9  +coeff( 45)        *x31*x41*x52
     1  +coeff( 46)            *x42*x52
     2  +coeff( 47)    *x22*x31*x42    
     3  +coeff( 48)        *x33*x42    
     4  +coeff( 49)        *x34    *x51
     5  +coeff( 50)        *x33*x41*x51
     6  +coeff( 51)    *x21    *x43*x51
     7  +coeff( 52)    *x22        *x53
     8  +coeff( 53)    *x22*x33*x41    
      x_sl_ep4    =x_sl_ep4    
     9  +coeff( 54)    *x24*x31    *x51
     1  +coeff( 55)    *x23*x32    *x51
     2  +coeff( 56)*x11            *x52
     3  +coeff( 57)        *x32    *x54
     4  +coeff( 58)*x11*x21*x31*x41    
     5  +coeff( 59)    *x22*x34*x41    
     6  +coeff( 60)*x11*x21    *x42    
     7  +coeff( 61)    *x24*x31*x42    
     8  +coeff( 62)    *x22*x32*x43    
      x_sl_ep4    =x_sl_ep4    
     9  +coeff( 63)    *x22*x31*x44    
     1  +coeff( 64)    *x22*x34    *x51
     2  +coeff( 65)    *x22*x33*x41*x51
     3  +coeff( 66)    *x24    *x42*x51
     4  +coeff( 67)    *x23*x31*x42*x51
     5  +coeff( 68)*x11        *x41*x52
     6  +coeff( 69)    *x22    *x43*x52
     7  +coeff( 70)    *x24        *x53
     8  +coeff( 71)    *x22*x31*x41*x53
      x_sl_ep4    =x_sl_ep4    
     9  +coeff( 72)    *x22    *x42*x53
     1  +coeff( 73)        *x31*x43*x53
     2  +coeff( 74)            *x44*x53
     3  +coeff( 75)    *x23        *x54
     4  +coeff( 76)    *x22    *x41*x54
     5  +coeff( 77)    *x24*x32*x42    
     6  +coeff( 78)    *x22*x34*x42    
     7  +coeff( 79)    *x24*x31*x43    
     8  +coeff( 80)    *x22*x33*x43    
      x_sl_ep4    =x_sl_ep4    
     9  +coeff( 81)    *x24    *x44    
     1  +coeff( 82)    *x22*x32*x44    
     2  +coeff( 83)        *x34*x44    
     3  +coeff( 84)    *x24*x33    *x51
     4  +coeff( 85)*x11*x22    *x41*x51
     5  +coeff( 86)    *x24        *x54
     6  +coeff( 87)    *x23*x31    *x54
     7  +coeff( 88)    *x22*x32    *x54
     8  +coeff( 89)*x11*x23*x31*x41    
      x_sl_ep4    =x_sl_ep4    
     9  +coeff( 90)*x11*x22*x32*x41    
c
      return
      end
      function t_sl_ep4    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff( 11)
      data ncoeff/ 10/
      data avdat/ -0.1264541E+00/
      data xmin/
     1 -0.19937E-02,-0.54847E-01,-0.19935E-01,-0.29994E-01,-0.49883E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19971E-02, 0.54952E-01, 0.19987E-01, 0.29928E-01, 0.49858E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.13134300E-02,-0.15654045E-03,-0.30272529E-01, 0.10292840E-02,
     + -0.13824236E-02,-0.19035802E-04,-0.12305843E-02,-0.93163784E-04,
     +  0.29923141E-03, 0.14061246E-03,
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
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
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
c
c                  function
c
      t_sl_ep4    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22    *x42    
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)    *x22*x31*x41    
     8  +coeff(  8)        *x32        
      t_sl_ep4    =t_sl_ep4    
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)    *x24*x31        
c
      return
      end
      function y_sl_ep4    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff(  4)
      data ncoeff/  3/
      data avdat/ -0.1571737E-02/
      data xmin/
     1 -0.19937E-02,-0.54847E-01,-0.19935E-01,-0.29994E-01,-0.49883E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19971E-02, 0.54952E-01, 0.19987E-01, 0.29928E-01, 0.49858E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.16418438E-02, 0.71774013E-01, 0.19889101E-02,
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
c          set up monomials   functions
      x11 = x1
      x21 = x2
c
c                  function
c
      y_sl_ep4    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
c
      return
      end
      function p_sl_ep4    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff(  5)
      data ncoeff/  4/
      data avdat/ -0.1169962E-02/
      data xmin/
     1 -0.19937E-02,-0.54847E-01,-0.19935E-01,-0.29994E-01,-0.49883E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19971E-02, 0.54952E-01, 0.19987E-01, 0.29928E-01, 0.49858E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.12198533E-02, 0.53944774E-01, 0.11733044E-02, 0.53302600E-03,
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
      p_sl_ep4    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)    *x21*x31        
c
      return
      end
      function l_sl_ep4    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff( 30)
      data ncoeff/ 29/
      data avdat/  0.1308969E+01/
      data xmin/
     1 -0.19937E-02,-0.54847E-01,-0.19935E-01,-0.29994E-01,-0.49883E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19971E-02, 0.54952E-01, 0.19987E-01, 0.29928E-01, 0.49858E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.69465500E-03,-0.11466041E-02, 0.36030692E-05,-0.16312405E-05,
     +  0.19662094E-02,-0.48486118E-05,-0.27855318E-04,-0.49246452E-03,
     + -0.39090887E-05,-0.79210062E-03,-0.79258916E-03,-0.88679401E-04,
     + -0.21552497E-03,-0.86679257E-03,-0.54738688E-03,-0.99079436E-04,
     +  0.86360716E-03, 0.45866014E-04, 0.19043771E-02, 0.19795897E-02,
     +  0.10001712E-02,-0.84780688E-04, 0.24951165E-03, 0.64245018E-03,
     + -0.47890440E-03, 0.38720533E-03, 0.82442659E-03, 0.28975311E-03,
     +  0.56008040E-03,
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
      l_sl_ep4    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21*x31        
     7  +coeff(  7)        *x32        
     8  +coeff(  8)        *x31*x41    
      l_sl_ep4    =l_sl_ep4    
     9  +coeff(  9)*x11*x21            
     1  +coeff( 10)    *x22*x31        
     2  +coeff( 11)        *x33        
     3  +coeff( 12)    *x22    *x41    
     4  +coeff( 13)        *x32*x41    
     5  +coeff( 14)        *x31*x42    
     6  +coeff( 15)        *x31    *x52
     7  +coeff( 16)            *x41*x52
     8  +coeff( 17)        *x33*x41    
      l_sl_ep4    =l_sl_ep4    
     9  +coeff( 18)    *x22    *x42    
     1  +coeff( 19)        *x32*x42    
     2  +coeff( 20)        *x31*x43    
     3  +coeff( 21)            *x44    
     4  +coeff( 22)*x12    *x31        
     5  +coeff( 23)    *x24*x31        
     6  +coeff( 24)    *x22*x33        
     7  +coeff( 25)        *x34*x41    
     8  +coeff( 26)        *x32*x43    
      l_sl_ep4    =l_sl_ep4    
     9  +coeff( 27)        *x31*x44    
     1  +coeff( 28)    *x22*x31    *x52
     2  +coeff( 29)        *x33    *x52
c
      return
      end
      function x_sl_ep5    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.1711894E+00/
      data xmin/
     1 -0.19937E-02,-0.54847E-01,-0.19935E-01,-0.29994E-01,-0.49883E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19971E-02, 0.54952E-01, 0.19987E-01, 0.29928E-01, 0.49858E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.26540679E-02,-0.19680865E-01,-0.44834692E-01, 0.48633415E-03,
     +  0.34484551E-05, 0.95697033E-04,-0.18354525E-02,-0.61026146E-03,
     +  0.20898509E-03, 0.54698129E-03, 0.97366705E-04, 0.20312634E-03,
     + -0.22448285E-04, 0.64974178E-04, 0.41810439E-04,-0.75239846E-05,
     + -0.28819422E-03, 0.36006726E-03, 0.38509006E-04, 0.17145518E-03,
     + -0.18414657E-04,-0.32548076E-05,-0.43964692E-05,-0.22631001E-02,
     + -0.19899765E-02,-0.57774241E-03,-0.18353766E-02,-0.56294311E-03,
     +  0.14356535E-03, 0.20926348E-03, 0.30050151E-04,-0.27846306E-05,
     + -0.96276890E-05,-0.84874910E-05, 0.69832304E-05,-0.13548440E-05,
     + -0.50499995E-06,-0.26395710E-05, 0.39005604E-05,-0.86373266E-05,
     +  0.46822282E-04,-0.11991718E-04, 0.21664515E-04, 0.10429694E-03,
     + -0.50434355E-04, 0.85367879E-03,-0.54158325E-04, 0.12042032E-02,
     + -0.81464881E-04, 0.67305961E-03,-0.50150196E-04, 0.86062297E-04,
     + -0.15182939E-04,-0.37045935E-02,-0.45199278E-04,-0.27993663E-04,
     + -0.18722260E-03,-0.21678487E-04, 0.29723172E-03, 0.40734975E-03,
     +  0.13013529E-03, 0.14225300E-03, 0.16979009E-03, 0.90648275E-04,
     +  0.59403006E-04, 0.40101339E-04, 0.52150872E-05,-0.13081029E-03,
     + -0.90351641E-04,-0.61284984E-04, 0.65854954E-04, 0.58470239E-04,
     +  0.80161881E-04, 0.77891047E-04,-0.23812859E-03,-0.14194509E-03,
     +  0.68185240E-04, 0.38211172E-04,-0.69707537E-04,-0.67661153E-04,
     + -0.42983543E-05,-0.51429961E-05, 0.12599425E-04,-0.17940813E-04,
     + -0.23099898E-04, 0.31245875E-05, 0.13115990E-04,-0.38526387E-05,
     + -0.26938436E-03,-0.56544936E-05,
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
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      x_sl_ep5    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x24    *x42    
     7  +coeff(  7)    *x24*x33*x41    
     8  +coeff(  8)    *x22*x31*x41    
      x_sl_ep5    =x_sl_ep5    
     9  +coeff(  9)    *x22*x34        
     1  +coeff( 10)    *x24*x31*x41    
     2  +coeff( 11)    *x24    *x43    
     3  +coeff( 12)            *x42    
     4  +coeff( 13)                *x52
     5  +coeff( 14)    *x22*x31        
     6  +coeff( 15)    *x22    *x41    
     7  +coeff( 16)    *x24            
     8  +coeff( 17)    *x22    *x42    
      x_sl_ep5    =x_sl_ep5    
     9  +coeff( 18)    *x24*x32        
     1  +coeff( 19)        *x32        
     2  +coeff( 20)        *x31*x41    
     3  +coeff( 21)    *x21        *x51
     4  +coeff( 22)    *x22        *x51
     5  +coeff( 23)*x11*x21            
     6  +coeff( 24)    *x22*x32*x42    
     7  +coeff( 25)    *x22    *x44    
     8  +coeff( 26)    *x24*x34        
      x_sl_ep5    =x_sl_ep5    
     9  +coeff( 27)    *x24*x32*x42    
     1  +coeff( 28)    *x24*x31*x43    
     2  +coeff( 29)    *x23*x32*x43    
     3  +coeff( 30)    *x24    *x44    
     4  +coeff( 31)    *x24*x33*x42    
     5  +coeff( 32)    *x21            
     6  +coeff( 33)            *x41*x51
     7  +coeff( 34)        *x33        
     8  +coeff( 35)        *x31*x42    
      x_sl_ep5    =x_sl_ep5    
     9  +coeff( 36)        *x31*x41*x51
     1  +coeff( 37)            *x41*x52
     2  +coeff( 38)*x11                
     3  +coeff( 39)    *x23*x31        
     4  +coeff( 40)        *x31*x43    
     5  +coeff( 41)    *x23        *x51
     6  +coeff( 42)            *x42*x52
     7  +coeff( 43)    *x21        *x53
     8  +coeff( 44)    *x22*x33        
      x_sl_ep5    =x_sl_ep5    
     9  +coeff( 45)    *x24    *x41    
     1  +coeff( 46)    *x22*x32*x41    
     2  +coeff( 47)        *x34*x41    
     3  +coeff( 48)    *x22*x31*x42    
     4  +coeff( 49)        *x33*x42    
     5  +coeff( 50)    *x22    *x43    
     6  +coeff( 51)    *x22    *x42*x51
     7  +coeff( 52)    *x21    *x42*x52
     8  +coeff( 53)            *x43*x52
      x_sl_ep5    =x_sl_ep5    
     9  +coeff( 54)    *x22*x31*x43    
     1  +coeff( 55)    *x23        *x53
     2  +coeff( 56)    *x23*x34        
     3  +coeff( 57)    *x22*x34*x41    
     4  +coeff( 58)*x11*x21    *x42    
     5  +coeff( 59)    *x23*x32*x42    
     6  +coeff( 60)    *x23*x31*x43    
     7  +coeff( 61)    *x23    *x44    
     8  +coeff( 62)    *x24*x31*x41*x51
      x_sl_ep5    =x_sl_ep5    
     9  +coeff( 63)    *x24    *x42*x51
     1  +coeff( 64)    *x23*x31*x44    
     2  +coeff( 65)*x11*x22    *x41*x51
     3  +coeff( 66)*x11*x22*x33        
     4  +coeff( 67)    *x24*x34    *x51
     5  +coeff( 68)    *x21*x32*x42*x54
     6  +coeff( 69)*x11*x22*x32*x41*x51
     7  +coeff( 70)*x11*x24*x33        
     8  +coeff( 71)*x11*x24*x32*x41    
      x_sl_ep5    =x_sl_ep5    
     9  +coeff( 72)*x11*x24*x31    *x52
     1  +coeff( 73)*x11*x23    *x42*x52
     2  +coeff( 74)    *x24*x34    *x53
     3  +coeff( 75)*x11*x22*x31*x41*x53
     4  +coeff( 76)*x11*x22    *x42*x53
     5  +coeff( 77)*x11*x22*x32*x44    
     6  +coeff( 78)*x11*x23*x33    *x52
     7  +coeff( 79)*x11*x23*x32*x44    
     8  +coeff( 80)*x11*x24*x34    *x51
      x_sl_ep5    =x_sl_ep5    
     9  +coeff( 81)        *x31    *x51
     1  +coeff( 82)    *x23            
     2  +coeff( 83)    *x21*x32        
     3  +coeff( 84)    *x21    *x42    
     4  +coeff( 85)            *x43    
     5  +coeff( 86)    *x21*x31    *x51
     6  +coeff( 87)    *x21    *x41*x51
     7  +coeff( 88)                *x53
     8  +coeff( 89)    *x22*x32        
      x_sl_ep5    =x_sl_ep5    
     9  +coeff( 90)    *x23    *x41    
c
      return
      end
      function t_sl_ep5    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff( 21)
      data ncoeff/ 20/
      data avdat/ -0.1630271E+00/
      data xmin/
     1 -0.19937E-02,-0.54847E-01,-0.19935E-01,-0.29994E-01,-0.49883E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19971E-02, 0.54952E-01, 0.19987E-01, 0.29928E-01, 0.49858E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.16984198E-02,-0.68164816E-04,-0.30001730E-01, 0.28301000E-02,
     + -0.29515405E-03, 0.14339187E-03, 0.20397073E-02,-0.57917032E-02,
     + -0.33616429E-03, 0.13714866E-03, 0.87304378E-03,-0.11540094E-02,
     + -0.51537380E-02,-0.14182026E-03, 0.66660234E-03,-0.87494584E-04,
     +  0.18602690E-04, 0.18678203E-02, 0.10628838E-03,-0.16942548E-02,
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
c
c                  function
c
      t_sl_ep5    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x22    *x41    
     8  +coeff(  8)    *x22    *x42    
      t_sl_ep5    =t_sl_ep5    
     9  +coeff(  9)    *x24*x31*x41    
     1  +coeff( 10)        *x31*x41    
     2  +coeff( 11)    *x22*x31        
     3  +coeff( 12)    *x22*x32        
     4  +coeff( 13)    *x22*x31*x41    
     5  +coeff( 14)                *x52
     6  +coeff( 15)    *x22    *x43    
     7  +coeff( 16)*x11*x21            
     8  +coeff( 17)        *x31*x42    
      t_sl_ep5    =t_sl_ep5    
     9  +coeff( 18)    *x22*x31*x42    
     1  +coeff( 19)    *x24        *x51
     2  +coeff( 20)    *x22*x31*x44    
c
      return
      end
      function y_sl_ep5    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff(  5)
      data ncoeff/  4/
      data avdat/ -0.1840557E-02/
      data xmin/
     1 -0.19937E-02,-0.54847E-01,-0.19935E-01,-0.29994E-01,-0.49883E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19971E-02, 0.54952E-01, 0.19987E-01, 0.29928E-01, 0.49858E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.19206532E-02, 0.83128229E-01, 0.19799103E-02, 0.20162747E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x41 = x4
c
c                  function
c
      y_sl_ep5    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
     4  +coeff(  4)    *x21    *x41    
c
      return
      end
      function p_sl_ep5    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff( 10)
      data ncoeff/  9/
      data avdat/ -0.1181728E-02/
      data xmin/
     1 -0.19937E-02,-0.54847E-01,-0.19935E-01,-0.29994E-01,-0.49883E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19971E-02, 0.54952E-01, 0.19987E-01, 0.29928E-01, 0.49858E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.12311526E-02, 0.53784274E-01, 0.12060371E-02, 0.56832336E-03,
     + -0.14726954E-02,-0.12551694E-05,-0.10337839E-02, 0.86084683E-03,
     +  0.80989097E-03,
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
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
c
c                  function
c
      p_sl_ep5    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)    *x21*x31        
     5  +coeff(  5)    *x21    *x42    
     6  +coeff(  6)    *x23*x31*x41    
     7  +coeff(  7)    *x21*x31*x41    
     8  +coeff(  8)    *x21    *x43    
      p_sl_ep5    =p_sl_ep5    
     9  +coeff(  9)    *x21*x33*x42    
c
      return
      end
      function l_sl_ep5    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff(  9)
      data ncoeff/  8/
      data avdat/  0.1519837E+01/
      data xmin/
     1 -0.19937E-02,-0.54847E-01,-0.19935E-01,-0.29994E-01,-0.49883E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19971E-02, 0.54952E-01, 0.19987E-01, 0.29928E-01, 0.49858E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.70498954E-03,-0.17008985E-02,-0.35857258E-02, 0.23313241E-04,
     +  0.22741158E-02,-0.31789350E-05,-0.72018988E-05, 0.67814591E-03,
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
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
c
c                  function
c
      l_sl_ep5    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)        *x32        
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)            *x42    
c
      return
      end
      function x_sl_ep6    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.2095920E+00/
      data xmin/
     1 -0.19937E-02,-0.54847E-01,-0.19935E-01,-0.29994E-01,-0.49883E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19971E-02, 0.54952E-01, 0.19987E-01, 0.29928E-01, 0.49858E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.30572843E-02,-0.19532645E-01,-0.50611146E-01, 0.12554757E-02,
     + -0.21279462E-03, 0.33782053E-03, 0.25291249E-03,-0.11414835E-02,
     +  0.58911566E-03, 0.19303017E-03,-0.87577973E-04,-0.14723986E-02,
     +  0.27521010E-03,-0.63334781E-04,-0.23827776E-03,-0.12011078E-04,
     +  0.25321946E-02,-0.11653418E-03, 0.34861943E-04,-0.68117392E-05,
     +  0.35689495E-03, 0.27743389E-02, 0.40650428E-02,-0.45211168E-03,
     +  0.16361726E-04,-0.70828209E-02,-0.76402030E-02,-0.39209686E-02,
     + -0.12947983E-02,-0.23629731E-02,-0.20828668E-02,-0.32664993E-03,
     + -0.32220915E-03, 0.26636559E-03,-0.22397587E-05,-0.10649434E-04,
     + -0.32105680E-04,-0.14070009E-04,-0.23598695E-03, 0.10278772E-04,
     + -0.11917518E-05, 0.46279656E-04, 0.56118013E-04,-0.49075083E-03,
     +  0.71501679E-04,-0.21878830E-02,-0.25666226E-04, 0.33648121E-04,
     + -0.21367929E-04,-0.65078167E-03,-0.10407310E-03, 0.54217363E-03,
     +  0.34952700E-05,-0.15738327E-02, 0.29386465E-04,-0.43075616E-02,
     + -0.19643963E-02,-0.70082519E-04, 0.43161665E-03, 0.17551619E-03,
     +  0.29309056E-03, 0.22183267E-03, 0.21638493E-03, 0.10418983E-03,
     + -0.25805290E-03,-0.13183578E-03,-0.15284825E-03, 0.15014407E-03,
     + -0.12059968E-04,-0.13521266E-04,-0.27603985E-05, 0.28962178E-04,
     +  0.51785843E-04, 0.21366717E-04, 0.40943569E-05,-0.43851429E-04,
     +  0.19574753E-03, 0.26990616E-03, 0.10596341E-03, 0.53466913E-04,
     + -0.26741071E-04, 0.24951058E-04,-0.12222685E-04, 0.79131787E-05,
     + -0.54271564E-04, 0.58154204E-04,-0.15489063E-03,-0.41533283E-04,
     + -0.36668938E-04,-0.69510796E-04,
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
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      x_sl_ep6    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x22    *x41    
     8  +coeff(  8)    *x22    *x42    
      x_sl_ep6    =x_sl_ep6    
     9  +coeff(  9)    *x24*x31*x41    
     1  +coeff( 10)    *x22*x31        
     2  +coeff( 11)    *x22*x32        
     3  +coeff( 12)    *x22*x31*x41    
     4  +coeff( 13)        *x31*x41    
     5  +coeff( 14)                *x52
     6  +coeff( 15)            *x43    
     7  +coeff( 16)*x11*x21            
     8  +coeff( 17)    *x22    *x43    
      x_sl_ep6    =x_sl_ep6    
     9  +coeff( 18)    *x22*x33*x42    
     1  +coeff( 19)        *x32        
     2  +coeff( 20)            *x44    
     3  +coeff( 21)    *x22*x33        
     4  +coeff( 22)    *x22*x32*x41    
     5  +coeff( 23)    *x22*x31*x42    
     6  +coeff( 24)        *x31*x44    
     7  +coeff( 25)    *x24        *x51
     8  +coeff( 26)    *x22*x32*x42    
      x_sl_ep6    =x_sl_ep6    
     9  +coeff( 27)    *x22*x31*x43    
     1  +coeff( 28)    *x22    *x44    
     2  +coeff( 29)    *x24*x33*x41    
     3  +coeff( 30)    *x24*x32*x42    
     4  +coeff( 31)    *x24*x31*x43    
     5  +coeff( 32)    *x22    *x44*x52
     6  +coeff( 33)*x11*x21*x31*x44    
     7  +coeff( 34)*x11*x22    *x41*x53
     8  +coeff( 35)    *x21            
      x_sl_ep6    =x_sl_ep6    
     9  +coeff( 36)    *x21        *x51
     1  +coeff( 37)            *x41*x51
     2  +coeff( 38)    *x23            
     3  +coeff( 39)        *x32*x41    
     4  +coeff( 40)    *x21        *x52
     5  +coeff( 41)*x11                
     6  +coeff( 42)    *x24            
     7  +coeff( 43)    *x23*x33        
     8  +coeff( 44)    *x22*x34        
      x_sl_ep6    =x_sl_ep6    
     9  +coeff( 45)    *x23*x32*x41    
     1  +coeff( 46)    *x22*x33*x41    
     2  +coeff( 47)    *x22*x31    *x53
     3  +coeff( 48)            *x43*x53
     4  +coeff( 49)*x11*x21*x31*x41    
     5  +coeff( 50)    *x22*x34*x41    
     6  +coeff( 51)*x11*x21    *x42    
     7  +coeff( 52)        *x33*x44    
     8  +coeff( 53)*x11*x21*x33        
      x_sl_ep6    =x_sl_ep6    
     9  +coeff( 54)    *x22*x34*x42    
     1  +coeff( 55)*x11*x21    *x43    
     2  +coeff( 56)    *x22*x33*x43    
     3  +coeff( 57)    *x22*x32*x44    
     4  +coeff( 58)*x11*x22    *x41*x51
     5  +coeff( 59)    *x22*x33*x42*x51
     6  +coeff( 60)    *x22    *x42*x54
     7  +coeff( 61)*x11*x21*x33*x42    
     8  +coeff( 62)*x11*x22*x31*x42*x51
      x_sl_ep6    =x_sl_ep6    
     9  +coeff( 63)    *x23*x33*x42*x52
     1  +coeff( 64)*x11*x24*x32*x41    
     2  +coeff( 65)*x11*x24*x31*x41*x51
     3  +coeff( 66)*x11*x24*x32    *x52
     4  +coeff( 67)*x11*x24*x34    *x51
     5  +coeff( 68)*x11*x24*x34*x42    
     6  +coeff( 69)    *x21*x31        
     7  +coeff( 70)    *x21*x31    *x51
     8  +coeff( 71)        *x31*x41*x51
      x_sl_ep6    =x_sl_ep6    
     9  +coeff( 72)        *x34        
     1  +coeff( 73)    *x23        *x51
     2  +coeff( 74)    *x21        *x53
     3  +coeff( 75)*x11    *x31        
     4  +coeff( 76)    *x24    *x41    
     5  +coeff( 77)        *x34*x41    
     6  +coeff( 78)        *x32*x43    
     7  +coeff( 79)    *x22*x31*x41*x51
     8  +coeff( 80)    *x22    *x42*x51
      x_sl_ep6    =x_sl_ep6    
     9  +coeff( 81)    *x21    *x43*x51
     1  +coeff( 82)    *x21    *x42*x52
     2  +coeff( 83)*x11*x22            
     3  +coeff( 84)*x11*x21        *x51
     4  +coeff( 85)    *x22*x33    *x51
     5  +coeff( 86)    *x23    *x42*x51
     6  +coeff( 87)    *x22*x31*x42*x51
     7  +coeff( 88)    *x21*x32*x42*x51
     8  +coeff( 89)        *x33*x42*x51
      x_sl_ep6    =x_sl_ep6    
     9  +coeff( 90)    *x23        *x53
c
      return
      end
      function t_sl_ep6    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff( 66)
      data ncoeff/ 65/
      data avdat/ -0.1999042E+00/
      data xmin/
     1 -0.19937E-02,-0.54847E-01,-0.19935E-01,-0.29994E-01,-0.49883E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19971E-02, 0.54952E-01, 0.19987E-01, 0.29928E-01, 0.49858E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.21643052E-02,-0.35641835E-04,-0.29628497E-01, 0.46498869E-02,
     + -0.11801572E-02, 0.34375702E-04, 0.18031120E-02,-0.13154118E-03,
     + -0.45762910E-02, 0.26340010E-02, 0.14490891E-02,-0.70277830E-02,
     +  0.10286656E-01, 0.18912512E-03,-0.20712422E-03, 0.14816165E-03,
     + -0.10327647E-02, 0.14657086E-01, 0.44372797E-03, 0.15690271E-03,
     + -0.13865622E-03, 0.11453238E-03,-0.76520484E-03, 0.59739007E-02,
     +  0.61915384E-03,-0.19043137E-03,-0.11563759E-01,-0.19617251E-02,
     +  0.76079258E-03,-0.40975995E-02,-0.14583202E-01, 0.15523431E-02,
     + -0.11614764E-02,-0.94325356E-02,-0.40584418E-03,-0.14107416E-01,
     + -0.10912966E-01, 0.33179764E-02, 0.40627266E-02, 0.44187286E-03,
     + -0.18375880E-02,-0.69662486E-02,-0.20890608E-02, 0.34116025E-03,
     +  0.16766477E-03, 0.90876705E-03,-0.30449656E-03, 0.67935861E-03,
     + -0.60064095E-03, 0.78069320E-03,-0.13120404E-01,-0.25459514E-02,
     + -0.30145359E-02,-0.45146695E-02,-0.84158854E-03,-0.19021002E-02,
     + -0.95317466E-02,-0.15343317E-01, 0.10737813E-02, 0.13759787E-02,
     + -0.16301398E-02, 0.71337406E-03, 0.32955112E-02, 0.13093791E-02,
     + -0.16203118E-02,
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
      t_sl_ep6    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x22    *x41    
     8  +coeff(  8)        *x33*x41    
      t_sl_ep6    =t_sl_ep6    
     9  +coeff(  9)    *x22    *x42    
     1  +coeff( 10)    *x24*x31*x41    
     2  +coeff( 11)    *x22*x31        
     3  +coeff( 12)    *x22*x31*x41    
     4  +coeff( 13)    *x22    *x43    
     5  +coeff( 14)    *x24*x32        
     6  +coeff( 15)                *x52
     7  +coeff( 16)        *x32*x41    
     8  +coeff( 17)    *x22*x32        
      t_sl_ep6    =t_sl_ep6    
     9  +coeff( 18)    *x22*x31*x42    
     1  +coeff( 19)    *x22    *x43*x51
     2  +coeff( 20)    *x22    *x44*x51
     3  +coeff( 21)*x11*x21            
     4  +coeff( 22)            *x44    
     5  +coeff( 23)    *x22    *x41*x51
     6  +coeff( 24)    *x22*x32*x41    
     7  +coeff( 25)        *x34*x41    
     8  +coeff( 26)*x11*x22    *x42    
      t_sl_ep6    =t_sl_ep6    
     9  +coeff( 27)    *x22*x32*x42    
     1  +coeff( 28)        *x34*x42    
     2  +coeff( 29)*x11*x21    *x43    
     3  +coeff( 30)        *x33*x43    
     4  +coeff( 31)    *x22    *x44    
     5  +coeff( 32)    *x22*x32*x41*x51
     6  +coeff( 33)    *x22*x33*x42    
     7  +coeff( 34)    *x22*x31*x44    
     8  +coeff( 35)    *x23*x34*x41    
      t_sl_ep6    =t_sl_ep6    
     9  +coeff( 36)    *x24*x31*x43    
     1  +coeff( 37)    *x22*x33*x43    
     2  +coeff( 38)    *x22*x33*x42*x51
     3  +coeff( 39)    *x22*x33*x44    
     4  +coeff( 40)*x12*x22    *x41*x52
     5  +coeff( 41)    *x22    *x44*x53
     6  +coeff( 42)    *x24*x34*x42    
     7  +coeff( 43)    *x22*x32*x44*x54
     8  +coeff( 44)        *x33        
      t_sl_ep6    =t_sl_ep6    
     9  +coeff( 45)*x11*x21*x31        
     1  +coeff( 46)    *x24    *x41    
     2  +coeff( 47)*x11*x21*x31*x41    
     3  +coeff( 48)*x11*x21    *x42    
     4  +coeff( 49)    *x22*x34        
     5  +coeff( 50)    *x24    *x42    
     6  +coeff( 51)    *x22*x31*x43    
     7  +coeff( 52)        *x32*x44    
     8  +coeff( 53)    *x24    *x43    
      t_sl_ep6    =t_sl_ep6    
     9  +coeff( 54)    *x22*x32*x43    
     1  +coeff( 55)        *x34*x43    
     2  +coeff( 56)*x11*x21    *x44    
     3  +coeff( 57)    *x24*x33*x41    
     4  +coeff( 58)    *x24*x32*x42    
     5  +coeff( 59)    *x24    *x43*x51
     6  +coeff( 60)    *x22    *x43*x53
     7  +coeff( 61)*x11*x21*x34*x42    
     8  +coeff( 62)*x12*x24    *x41*x51
      t_sl_ep6    =t_sl_ep6    
     9  +coeff( 63)*x11*x21*x34*x44    
     1  +coeff( 64)    *x22*x34*x44*x51
     2  +coeff( 65)*x11*x21*x32*x44*x52
c
      return
      end
      function y_sl_ep6    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff(  5)
      data ncoeff/  4/
      data avdat/ -0.2079399E-02/
      data xmin/
     1 -0.19937E-02,-0.54847E-01,-0.19935E-01,-0.29994E-01,-0.49883E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19971E-02, 0.54952E-01, 0.19987E-01, 0.29928E-01, 0.49858E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.21680333E-02, 0.94517022E-01, 0.19662313E-02, 0.35358741E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x41 = x4
c
c                  function
c
      y_sl_ep6    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
     4  +coeff(  4)    *x21    *x41    
c
      return
      end
      function p_sl_ep6    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff( 15)
      data ncoeff/ 14/
      data avdat/ -0.1180994E-02/
      data xmin/
     1 -0.19937E-02,-0.54847E-01,-0.19935E-01,-0.29994E-01,-0.49883E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19971E-02, 0.54952E-01, 0.19987E-01, 0.29928E-01, 0.49858E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.12295432E-02, 0.54067232E-01, 0.33274808E-03,-0.22634193E-02,
     +  0.40505352E-03,-0.51554287E-03,-0.28629757E-02, 0.45278529E-02,
     +  0.28097141E-02, 0.32735360E-02, 0.31040795E-02,-0.27902434E-02,
     + -0.41475762E-02,-0.60472819E-02,
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
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
c
c                  function
c
      p_sl_ep6    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)    *x21    *x42    
     5  +coeff(  5)    *x21*x31        
     6  +coeff(  6)    *x21*x32        
     7  +coeff(  7)    *x21*x31*x41    
     8  +coeff(  8)    *x21    *x43    
      p_sl_ep6    =p_sl_ep6    
     9  +coeff(  9)    *x21*x34*x41    
     1  +coeff( 10)    *x21*x33*x42    
     2  +coeff( 11)    *x21*x31*x42    
     3  +coeff( 12)    *x21    *x44    
     4  +coeff( 13)    *x21*x34*x42    
     5  +coeff( 14)    *x21*x33*x43    
c
      return
      end
      function l_sl_ep6    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff(  8)
      data ncoeff/  7/
      data avdat/  0.1732405E+01/
      data xmin/
     1 -0.19937E-02,-0.54847E-01,-0.19935E-01,-0.29994E-01,-0.49883E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19971E-02, 0.54952E-01, 0.19987E-01, 0.29928E-01, 0.49858E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.68762410E-03,-0.25485421E-02,-0.58466783E-02, 0.83711297E-04,
     +  0.25553477E-02,-0.27040685E-04, 0.73924253E-03,
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
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
c
c                  function
c
      l_sl_ep6    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
c
      return
      end
      function x_sl_ep7    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.2586834E+00/
      data xmin/
     1 -0.19937E-02,-0.54847E-01,-0.19935E-01,-0.29994E-01,-0.49883E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19971E-02, 0.54952E-01, 0.19987E-01, 0.29928E-01, 0.49858E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.37696422E-02,-0.19085709E-01,-0.56077320E-01, 0.23908974E-02,
     + -0.78972452E-03, 0.45386094E-03, 0.85465296E-03,-0.16831243E-03,
     + -0.21078542E-02, 0.13368718E-02, 0.25814283E-03,-0.23713533E-03,
     + -0.33440948E-02, 0.48195771E-02, 0.41686458E-03,-0.11418525E-03,
     + -0.16421509E-03, 0.91327475E-02,-0.24759420E-03,-0.14674204E-02,
     + -0.57745645E-04,-0.40382478E-04, 0.42601284E-02, 0.58075684E-05,
     + -0.13936278E-03,-0.12615506E-01,-0.78033637E-02,-0.45650024E-04,
     + -0.12158271E-02,-0.32770203E-02, 0.14336308E-03, 0.33445216E-04,
     + -0.41327467E-02,-0.49188938E-02,-0.47525112E-02,-0.23753986E-03,
     + -0.11739997E-02,-0.24836510E-02,-0.28095548E-02,-0.87223307E-03,
     +  0.89120446E-03, 0.12875847E-04,-0.78807325E-05, 0.96655400E-04,
     + -0.26419156E-03, 0.82500355E-03, 0.67181769E-04,-0.82799431E-03,
     +  0.36951173E-04,-0.11942947E-01,-0.24154224E-04,-0.91462360E-04,
     +  0.33370455E-03, 0.51263349E-04, 0.23407930E-03,-0.52673724E-02,
     + -0.25969022E-03,-0.16269082E-02, 0.32213305E-02,-0.45158694E-03,
     +  0.37429892E-03, 0.54574467E-03,-0.97570795E-03,-0.12307393E-02,
     +  0.63035644E-04,-0.63269275E-04,-0.32202213E-04,-0.15062284E-03,
     + -0.97713591E-05,-0.39027491E-04, 0.17806712E-03, 0.17855158E-04,
     + -0.43349096E-03,-0.60488481E-04, 0.29835331E-06, 0.67888184E-04,
     +  0.92798880E-04, 0.86899781E-04, 0.11336972E-03,-0.11296054E-03,
     + -0.58357700E-05, 0.81435937E-04,-0.60169878E-04,-0.23275586E-04,
     + -0.26962850E-02,-0.11535339E-02,-0.18184728E-04, 0.92030612E-04,
     +  0.28102691E-03, 0.89589150E-04,
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
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      x_sl_ep7    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x22    *x41    
     8  +coeff(  8)        *x33*x41    
      x_sl_ep7    =x_sl_ep7    
     9  +coeff(  9)    *x22    *x42    
     1  +coeff( 10)    *x24*x31*x41    
     2  +coeff( 11)    *x22*x31        
     3  +coeff( 12)    *x22*x32        
     4  +coeff( 13)    *x22*x31*x41    
     5  +coeff( 14)    *x22    *x43    
     6  +coeff( 15)        *x31*x41    
     7  +coeff( 16)                *x52
     8  +coeff( 17)        *x31*x42    
      x_sl_ep7    =x_sl_ep7    
     9  +coeff( 18)    *x22*x31*x42    
     1  +coeff( 19)    *x22    *x44*x52
     2  +coeff( 20)    *x24*x31*x44*x53
     3  +coeff( 21)        *x32*x41    
     4  +coeff( 22)*x11*x21            
     5  +coeff( 23)    *x22*x32*x41    
     6  +coeff( 24)        *x34*x41    
     7  +coeff( 25)            *x42*x53
     8  +coeff( 26)    *x22*x32*x42    
      x_sl_ep7    =x_sl_ep7    
     9  +coeff( 27)    *x22    *x44    
     1  +coeff( 28)        *x32*x44    
     2  +coeff( 29)    *x22*x33*x42    
     3  +coeff( 30)    *x22*x31*x44    
     4  +coeff( 31)    *x23    *x43*x51
     5  +coeff( 32)*x11*x22*x31*x41    
     6  +coeff( 33)    *x24*x33*x41    
     7  +coeff( 34)    *x24*x31*x43    
     8  +coeff( 35)    *x22*x33*x43    
      x_sl_ep7    =x_sl_ep7    
     9  +coeff( 36)    *x24    *x44    
     1  +coeff( 37)        *x34*x44    
     2  +coeff( 38)    *x22*x32*x43*x51
     3  +coeff( 39)    *x24*x34*x42    
     4  +coeff( 40)    *x24*x33*x43    
     5  +coeff( 41)*x11*x24*x31*x41*x54
     6  +coeff( 42)    *x21            
     7  +coeff( 43)        *x33        
     8  +coeff( 44)    *x24            
      x_sl_ep7    =x_sl_ep7    
     9  +coeff( 45)        *x31*x43    
     1  +coeff( 46)    *x22*x33        
     2  +coeff( 47)    *x22*x31*x41*x51
     3  +coeff( 48)    *x22*x34        
     4  +coeff( 49)*x11*x21    *x41    
     5  +coeff( 50)    *x22*x31*x43    
     6  +coeff( 51)    *x23        *x53
     7  +coeff( 52)*x11*x21    *x42    
     8  +coeff( 53)    *x24*x31*x41*x51
      x_sl_ep7    =x_sl_ep7    
     9  +coeff( 54)*x11*x23*x31        
     1  +coeff( 55)*x11*x21*x31*x42    
     2  +coeff( 56)    *x24*x32*x42    
     3  +coeff( 57)*x11*x22    *x41*x51
     4  +coeff( 58)    *x22*x31*x44*x51
     5  +coeff( 59)    *x22*x33*x44    
     6  +coeff( 60)*x11*x21*x31*x44    
     7  +coeff( 61)    *x24*x34*x41*x51
     8  +coeff( 62)*x11*x22    *x41*x53
      x_sl_ep7    =x_sl_ep7    
     9  +coeff( 63)    *x23*x31*x44*x53
     1  +coeff( 64)    *x24*x34*x44*x52
     2  +coeff( 65)        *x32        
     3  +coeff( 66)            *x41*x51
     4  +coeff( 67)    *x23            
     5  +coeff( 68)            *x43    
     6  +coeff( 69)    *x21*x31    *x51
     7  +coeff( 70)    *x21    *x41*x51
     8  +coeff( 71)            *x42*x51
      x_sl_ep7    =x_sl_ep7    
     9  +coeff( 72)    *x21*x33        
     1  +coeff( 73)        *x32*x42    
     2  +coeff( 74)    *x22*x31    *x51
     3  +coeff( 75)        *x33    *x51
     4  +coeff( 76)        *x32*x41*x51
     5  +coeff( 77)    *x24*x31        
     6  +coeff( 78)    *x22*x32    *x51
     7  +coeff( 79)    *x21*x31*x42*x51
     8  +coeff( 80)            *x44*x51
      x_sl_ep7    =x_sl_ep7    
     9  +coeff( 81)    *x21*x32    *x52
     1  +coeff( 82)    *x21*x31*x41*x52
     2  +coeff( 83)            *x43*x52
     3  +coeff( 84)*x11*x22            
     4  +coeff( 85)    *x22*x33*x41    
     5  +coeff( 86)        *x33*x43    
     6  +coeff( 87)    *x23    *x42*x51
     7  +coeff( 88)        *x33*x41*x52
     8  +coeff( 89)    *x22    *x42*x52
      x_sl_ep7    =x_sl_ep7    
     9  +coeff( 90)    *x21*x31*x41*x53
c
      return
      end
      function t_sl_ep7    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.2210230E+00/
      data xmin/
     1 -0.19937E-02,-0.54847E-01,-0.19935E-01,-0.29994E-01,-0.49883E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19971E-02, 0.54952E-01, 0.19987E-01, 0.29928E-01, 0.49858E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.26494064E-02,-0.35552242E-04,-0.29884681E-01, 0.56848023E-02,
     + -0.17578003E-02, 0.67309811E-04, 0.36483724E-02, 0.68843749E-03,
     + -0.10688721E-01, 0.92408393E-03, 0.17123900E-02,-0.10798687E-01,
     +  0.79726567E-02, 0.10294260E-02,-0.17458433E-01, 0.52472576E-02,
     + -0.12791139E-01,-0.49373969E-02,-0.11815183E-01,-0.29085183E-02,
     +  0.12327428E-03, 0.20752932E-03,-0.29353268E-03,-0.43308045E-03,
     +  0.15520173E-03,-0.37618603E-02, 0.25621694E-03,-0.33646362E-03,
     +  0.16140414E-01, 0.48695397E-03, 0.13299522E-05, 0.58422848E-02,
     +  0.77172235E-03, 0.21069467E-02, 0.16183474E-02,-0.18666329E-01,
     + -0.45179608E-02, 0.80084224E-03,-0.37880046E-02,-0.10771444E-01,
     +  0.10933101E-02,-0.19263494E-02,-0.26583229E-02,-0.32279757E-02,
     + -0.81860322E-04,-0.67595589E-04,-0.15355652E-03,-0.40516745E-04,
     + -0.28130706E-03,-0.38488040E-03,-0.24891814E-03,-0.18383679E-03,
     + -0.92246810E-04,-0.11461529E-03, 0.10263058E-03, 0.57327914E-04,
     +  0.29680907E-03, 0.63933018E-02, 0.15555845E-02, 0.38855735E-03,
     +  0.50865981E-03, 0.16582500E-02, 0.42337496E-03, 0.68575836E-03,
     + -0.71785913E-03,-0.40725223E-02,-0.12721326E-02, 0.10591863E-02,
     +  0.29149934E-03,-0.59262896E-03, 0.77737286E-03,-0.10000979E-01,
     + -0.40192605E-03,-0.95905755E-02,-0.18829175E-02,-0.48871071E-03,
     + -0.12598147E-02, 0.98342879E-03, 0.84518781E-03,-0.12974433E-02,
     + -0.32254632E-03, 0.42775451E-03,-0.26958800E-03, 0.26443019E-03,
     + -0.32867995E-05, 0.56321645E-03,-0.13679297E-01,-0.99102012E-03,
     +  0.34599062E-02,-0.85987244E-03,
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
      t_sl_ep7    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x22    *x41    
     8  +coeff(  8)        *x33*x41    
      t_sl_ep7    =t_sl_ep7    
     9  +coeff(  9)    *x22    *x42    
     1  +coeff( 10)    *x24*x31*x41    
     2  +coeff( 11)    *x22*x31        
     3  +coeff( 12)    *x22*x31*x41    
     4  +coeff( 13)    *x22    *x43    
     5  +coeff( 14)    *x24*x32        
     6  +coeff( 15)    *x22*x33*x42    
     7  +coeff( 16)    *x22*x31*x43*x52
     8  +coeff( 17)    *x22*x31*x44*x53
      t_sl_ep7    =t_sl_ep7    
     9  +coeff( 18)    *x22*x34*x41    
     1  +coeff( 19)    *x22*x32*x43*x53
     2  +coeff( 20)*x12*x22    *x42*x53
     3  +coeff( 21)        *x32        
     4  +coeff( 22)    *x21        *x51
     5  +coeff( 23)                *x52
     6  +coeff( 24)    *x21    *x41*x51
     7  +coeff( 25)    *x24            
     8  +coeff( 26)    *x22*x32        
      t_sl_ep7    =t_sl_ep7    
     9  +coeff( 27)    *x21    *x43    
     1  +coeff( 28)    *x21    *x41*x52
     2  +coeff( 29)    *x22*x31*x42    
     3  +coeff( 30)    *x23    *x41*x51
     4  +coeff( 31)*x11*x22*x32        
     5  +coeff( 32)    *x22    *x42*x52
     6  +coeff( 33)*x12*x22    *x41    
     7  +coeff( 34)    *x22*x31*x43*x51
     8  +coeff( 35)*x11*x22*x31*x43    
      t_sl_ep7    =t_sl_ep7    
     9  +coeff( 36)    *x24*x31*x43    
     1  +coeff( 37)    *x24    *x44    
     2  +coeff( 38)*x12*x22    *x41*x51
     3  +coeff( 39)    *x22    *x42*x54
     4  +coeff( 40)    *x22*x34*x43    
     5  +coeff( 41)*x12*x22    *x42*x51
     6  +coeff( 42)    *x24*x34*x42    
     7  +coeff( 43)*x12*x23*x31*x43    
     8  +coeff( 44)*x12*x23    *x42*x52
      t_sl_ep7    =t_sl_ep7    
     9  +coeff( 45)            *x41*x51
     1  +coeff( 46)*x11*x21            
     2  +coeff( 47)    *x21*x31    *x51
     3  +coeff( 48)        *x34        
     4  +coeff( 49)*x11*x21    *x41    
     5  +coeff( 50)    *x23        *x51
     6  +coeff( 51)    *x21*x31*x41*x51
     7  +coeff( 52)    *x21*x31    *x52
     8  +coeff( 53)    *x21        *x53
      t_sl_ep7    =t_sl_ep7    
     9  +coeff( 54)*x11*x23            
     1  +coeff( 55)*x12    *x31        
     2  +coeff( 56)*x11*x22    *x41    
     3  +coeff( 57)*x11*x21*x31*x41    
     4  +coeff( 58)    *x22*x32*x41    
     5  +coeff( 59)    *x22*x32    *x51
     6  +coeff( 60)        *x32*x42*x51
     7  +coeff( 61)    *x21    *x42*x52
     8  +coeff( 62)    *x22*x34        
      t_sl_ep7    =t_sl_ep7    
     9  +coeff( 63)*x11*x23    *x41    
     1  +coeff( 64)    *x21*x33*x42    
     2  +coeff( 65)        *x33*x43    
     3  +coeff( 66)    *x22    *x44    
     4  +coeff( 67)    *x21*x31*x42*x52
     5  +coeff( 68)    *x21*x31*x41*x53
     6  +coeff( 69)*x12*x22*x31        
     7  +coeff( 70)*x11*x23*x31*x41    
     8  +coeff( 71)*x11*x23    *x42    
      t_sl_ep7    =t_sl_ep7    
     9  +coeff( 72)    *x22*x32*x43    
     1  +coeff( 73)*x11    *x31*x44    
     2  +coeff( 74)    *x22*x31*x44    
     3  +coeff( 75)    *x24*x32    *x51
     4  +coeff( 76)*x11*x23    *x41*x51
     5  +coeff( 77)*x11*x22*x31*x41*x51
     6  +coeff( 78)    *x23    *x43*x51
     7  +coeff( 79)*x11*x21*x31*x41*x52
     8  +coeff( 80)    *x21    *x44*x52
      t_sl_ep7    =t_sl_ep7    
     9  +coeff( 81)*x11*x22        *x53
     1  +coeff( 82)    *x23*x31    *x53
     2  +coeff( 83)        *x34    *x53
     3  +coeff( 84)    *x22    *x41*x54
     4  +coeff( 85)*x11*x23*x32*x41    
     5  +coeff( 86)*x12*x22    *x42    
     6  +coeff( 87)    *x24*x32*x42    
     7  +coeff( 88)*x11*x21*x33*x42    
     8  +coeff( 89)    *x22*x34*x42    
      t_sl_ep7    =t_sl_ep7    
     9  +coeff( 90)    *x23*x32*x43    
c
      return
      end
      function y_sl_ep7    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff(  7)
      data ncoeff/  6/
      data avdat/ -0.2361445E-02/
      data xmin/
     1 -0.19937E-02,-0.54847E-01,-0.19935E-01,-0.29994E-01,-0.49883E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19971E-02, 0.54952E-01, 0.19987E-01, 0.29928E-01, 0.49858E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.24625417E-02, 0.10707038E+00, 0.19613858E-02,-0.98968600E-03,
     +  0.53056318E-03,-0.90006803E-03,
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
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
c
c                  function
c
      y_sl_ep7    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
     4  +coeff(  4)    *x21    *x42    
     5  +coeff(  5)    *x21    *x43    
     6  +coeff(  6)    *x21*x33*x41    
c
      return
      end
      function p_sl_ep7    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff( 16)
      data ncoeff/ 15/
      data avdat/ -0.1173859E-02/
      data xmin/
     1 -0.19937E-02,-0.54847E-01,-0.19935E-01,-0.29994E-01,-0.49883E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19971E-02, 0.54952E-01, 0.19987E-01, 0.29928E-01, 0.49858E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.12228359E-02, 0.53733408E-01, 0.57194341E-03,-0.49765687E-02,
     +  0.45103702E-03,-0.40581268E-02,-0.82721008E-03, 0.47583580E-02,
     +  0.18395043E-03, 0.22737966E-02, 0.10280502E-03, 0.21465728E-03,
     +  0.40483638E-02, 0.26179797E-02, 0.52534195E-03,
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
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
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
      x51 = x5
c
c                  function
c
      p_sl_ep7    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)    *x21    *x42    
     5  +coeff(  5)    *x21*x31        
     6  +coeff(  6)    *x21*x31*x41    
     7  +coeff(  7)    *x21*x32        
     8  +coeff(  8)    *x21    *x43    
      p_sl_ep7    =p_sl_ep7    
     9  +coeff(  9)    *x23*x32*x41    
     1  +coeff( 10)    *x21*x33*x42    
     2  +coeff( 11)    *x21        *x51
     3  +coeff( 12)    *x21*x33        
     4  +coeff( 13)    *x21*x31*x42    
     5  +coeff( 14)    *x21*x34*x41    
     6  +coeff( 15)    *x23    *x43*x51
c
      return
      end
      function l_sl_ep7    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff( 10)
      data ncoeff/  9/
      data avdat/  0.1963088E+01/
      data xmin/
     1 -0.19937E-02,-0.54847E-01,-0.19935E-01,-0.29994E-01,-0.49883E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19971E-02, 0.54952E-01, 0.19987E-01, 0.29928E-01, 0.49858E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.42563619E-03,-0.44119740E-02,-0.11551324E-01, 0.34379959E-03,
     +  0.28227067E-02,-0.12977577E-03, 0.82502782E-03, 0.20832254E-03,
     + -0.34278538E-03,
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
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
c
c                  function
c
      l_sl_ep7    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x22    *x41    
      l_sl_ep7    =l_sl_ep7    
     9  +coeff(  9)    *x22    *x42    
c
      return
      end
      function x_sl_q1ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.3039816E-02/
      data xmin/
     1 -0.19937E-02,-0.54458E-01,-0.19935E-01,-0.29994E-01,-0.49883E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19971E-02, 0.54459E-01, 0.19987E-01, 0.29928E-01, 0.49858E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.30510749E-02, 0.14510249E+00, 0.31367748E-02, 0.24070508E-03,
     + -0.57038888E-02, 0.11009210E-02, 0.31529888E-03,-0.45846715E-02,
     + -0.86786848E-03, 0.12469088E-01,-0.71875268E-03,-0.47001438E-02,
     + -0.22703824E-03, 0.73859608E-02, 0.16859731E-01, 0.36230084E-03,
     + -0.20996406E-03, 0.11492996E-02,-0.88413870E-02,-0.28670658E-03,
     +  0.54347690E-03,-0.55362633E-02,-0.65504638E-02,-0.35037319E-03,
     +  0.14743416E-03,-0.13811247E-01, 0.48686095E-04,-0.28542473E-03,
     +  0.56453969E-05,-0.72030896E-04,-0.11707053E-01,-0.13665970E-01,
     + -0.58447514E-02,-0.23910715E-02, 0.17896922E-02, 0.14029285E-03,
     + -0.27205137E-03, 0.97167548E-02, 0.52391165E-02, 0.45158959E-03,
     +  0.14431935E-04,-0.25306176E-03,-0.30685768E-02, 0.17052885E-02,
     + -0.43413095E-03,-0.56686287E-03,-0.43399301E-02, 0.10432590E-02,
     + -0.87636411E-02, 0.84205094E-03,-0.19303948E-02, 0.50280977E-03,
     +  0.66711902E-04, 0.23372278E-02, 0.25616330E-03,-0.95577394E-04,
     +  0.13903929E-03,-0.34895007E-03,-0.52062878E-02, 0.45346475E-03,
     +  0.14337383E-02, 0.51304549E-02, 0.63863197E-04, 0.85363956E-03,
     + -0.22135812E-02, 0.15705773E-02,-0.20621091E-03,-0.87477807E-02,
     +  0.82190549E-02,-0.13061179E-02, 0.15262783E-01, 0.36138117E-05,
     + -0.50666433E-04,-0.52737720E-04, 0.27495484E-04, 0.58351165E-04,
     +  0.54376804E-04, 0.43943903E-04, 0.31273335E-04, 0.49944036E-03,
     + -0.21080628E-03,-0.82083087E-03, 0.49253268E-03,-0.67556633E-04,
     +  0.29939885E-04,-0.33117238E-04,-0.20001081E-03,-0.62101460E-02,
     +  0.13585857E-03,-0.26012536E-02,
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
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      x_sl_q1ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21        *x51
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)    *x21    *x42    
     6  +coeff(  6)*x11                
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)    *x21*x31*x41    
      x_sl_q1ex   =x_sl_q1ex   
     9  +coeff(  9)    *x21*x32        
     1  +coeff( 10)    *x21    *x43    
     2  +coeff( 11)    *x21*x34*x41    
     3  +coeff( 12)    *x21*x33*x42    
     4  +coeff( 13)    *x23            
     5  +coeff( 14)    *x21*x32*x41    
     6  +coeff( 15)    *x21*x31*x42    
     7  +coeff( 16)    *x23    *x41*x51
     8  +coeff( 17)    *x21        *x52
      x_sl_q1ex   =x_sl_q1ex   
     9  +coeff( 18)    *x21*x33        
     1  +coeff( 19)    *x21    *x44    
     2  +coeff( 20)    *x23    *x42*x51
     3  +coeff( 21)    *x22    *x42*x52
     4  +coeff( 22)    *x21*x34*x42    
     5  +coeff( 23)    *x21*x33*x43    
     6  +coeff( 24)    *x23*x31*x41*x52
     7  +coeff( 25)    *x22*x31*x41    
     8  +coeff( 26)    *x21*x31*x43    
      x_sl_q1ex   =x_sl_q1ex   
     9  +coeff( 27)*x11            *x51
     1  +coeff( 28)    *x21    *x43*x51
     2  +coeff( 29)    *x22*x32*x42    
     3  +coeff( 30)    *x22    *x41*x53
     4  +coeff( 31)    *x23*x32*x42    
     5  +coeff( 32)    *x23*x31*x43    
     6  +coeff( 33)    *x23    *x44    
     7  +coeff( 34)    *x21*x32*x44    
     8  +coeff( 35)    *x23    *x43*x51
      x_sl_q1ex   =x_sl_q1ex   
     9  +coeff( 36)*x11        *x42*x52
     1  +coeff( 37)*x11*x23    *x42    
     2  +coeff( 38)    *x23*x34*x42    
     3  +coeff( 39)    *x23*x33*x43    
     4  +coeff( 40)    *x23*x34*x41*x51
     5  +coeff( 41)*x11        *x43*x52
     6  +coeff( 42)*x11        *x42*x54
     7  +coeff( 43)    *x24    *x44*x54
     8  +coeff( 44)    *x24*x34*x43*x52
      x_sl_q1ex   =x_sl_q1ex   
     9  +coeff( 45)    *x21    *x41*x52
     1  +coeff( 46)    *x21*x34        
     2  +coeff( 47)    *x21*x33*x41    
     3  +coeff( 48)    *x23    *x42    
     4  +coeff( 49)    *x21*x32*x42    
     5  +coeff( 50)    *x21    *x42*x52
     6  +coeff( 51)    *x23    *x43    
     7  +coeff( 52)    *x23    *x41*x52
     8  +coeff( 53)*x11*x22*x31        
      x_sl_q1ex   =x_sl_q1ex   
     9  +coeff( 54)    *x23*x31*x42*x51
     1  +coeff( 55)    *x23*x32    *x52
     2  +coeff( 56)    *x22*x32*x41*x52
     3  +coeff( 57)    *x21*x33*x41*x52
     4  +coeff( 58)    *x22    *x42*x53
     5  +coeff( 59)    *x23*x31*x44    
     6  +coeff( 60)    *x21*x32*x41*x54
     7  +coeff( 61)    *x23*x31*x43*x52
     8  +coeff( 62)    *x23*x33*x44    
      x_sl_q1ex   =x_sl_q1ex   
     9  +coeff( 63)*x11*x22*x34*x41    
     1  +coeff( 64)*x11*x22*x34*x42    
     2  +coeff( 65)*x11*x24*x31*x43    
     3  +coeff( 66)*x11*x24    *x44    
     4  +coeff( 67)*x11    *x33*x41*x54
     5  +coeff( 68)    *x24*x31*x44*x54
     6  +coeff( 69)*x11*x24*x33*x43    
     7  +coeff( 70)*x11*x22    *x44*x54
     8  +coeff( 71)    *x24*x33*x44*x54
      x_sl_q1ex   =x_sl_q1ex   
     9  +coeff( 72)                *x51
     1  +coeff( 73)            *x41*x51
     2  +coeff( 74)                *x52
     3  +coeff( 75)    *x21*x31    *x51
     4  +coeff( 76)            *x42*x51
     5  +coeff( 77)            *x41*x53
     6  +coeff( 78)                *x54
     7  +coeff( 79)*x11        *x41    
     8  +coeff( 80)    *x23*x31*x41    
      x_sl_q1ex   =x_sl_q1ex   
     9  +coeff( 81)    *x21*x32*x41*x51
     1  +coeff( 82)    *x21*x31*x42*x51
     2  +coeff( 83)    *x21*x31*x41*x52
     3  +coeff( 84)            *x42*x53
     4  +coeff( 85)*x11*x21    *x41    
     5  +coeff( 86)*x11        *x42    
     6  +coeff( 87)    *x23*x31*x42    
     7  +coeff( 88)    *x21*x32*x43    
     8  +coeff( 89)    *x22    *x44    
      x_sl_q1ex   =x_sl_q1ex   
     9  +coeff( 90)    *x21*x31*x44    
c
      return
      end
      function t_sl_q1ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff( 31)
      data ncoeff/ 30/
      data avdat/  0.2312069E-03/
      data xmin/
     1 -0.19937E-02,-0.54458E-01,-0.19935E-01,-0.29994E-01,-0.49883E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19971E-02, 0.54459E-01, 0.19987E-01, 0.29928E-01, 0.49858E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.23191998E-03,-0.92335697E-02,-0.82722859E-03, 0.11232383E-03,
     +  0.29868258E-02,-0.20311507E-02, 0.10035558E-03,-0.37047407E-03,
     + -0.19159602E-02,-0.93478803E-03,-0.24122167E-03, 0.17156181E-03,
     + -0.15081986E-03, 0.25659776E-02, 0.24888411E-03,-0.27594829E-03,
     +  0.36327721E-04, 0.16068388E-03, 0.29148939E-02,-0.17363932E-02,
     +  0.39338827E-03, 0.81609141E-04,-0.11778154E-03, 0.77650271E-03,
     +  0.71554699E-04,-0.27210548E-03, 0.10582149E-02, 0.27519098E-03,
     + -0.12873583E-02, 0.27884549E-03,
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
      t_sl_q1ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)    *x21    *x42    
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)    *x21*x32        
      t_sl_q1ex   =t_sl_q1ex   
     9  +coeff(  9)    *x21*x31*x41    
     1  +coeff( 10)    *x23    *x42    
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x21    *x41*x51
     4  +coeff( 13)    *x21        *x52
     5  +coeff( 14)    *x21    *x43    
     6  +coeff( 15)    *x23*x31*x42    
     7  +coeff( 16)    *x21*x32*x43    
     8  +coeff( 17)*x11            *x51
      t_sl_q1ex   =t_sl_q1ex   
     9  +coeff( 18)    *x21*x33        
     1  +coeff( 19)    *x21*x31*x42    
     2  +coeff( 20)    *x21    *x44    
     3  +coeff( 21)    *x23*x32*x41    
     4  +coeff( 22)    *x21*x34*x41    
     5  +coeff( 23)    *x21*x33*x43    
     6  +coeff( 24)    *x23    *x43*x51
     7  +coeff( 25)    *x21*x31    *x51
     8  +coeff( 26)    *x23    *x41    
      t_sl_q1ex   =t_sl_q1ex   
     9  +coeff( 27)    *x21*x32*x41    
     1  +coeff( 28)    *x21*x33*x41    
     2  +coeff( 29)    *x21*x31*x43    
     3  +coeff( 30)*x12*x23    *x41*x51
c
      return
      end
      function y_sl_q1ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff( 39)
      data ncoeff/ 38/
      data avdat/  0.9249179E-02/
      data xmin/
     1 -0.19937E-02,-0.54458E-01,-0.19935E-01,-0.29994E-01,-0.49883E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19971E-02, 0.54459E-01, 0.19987E-01, 0.29928E-01, 0.49858E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.71000657E-02, 0.29154683E-01, 0.15293354E+00,-0.16561029E-01,
     + -0.75657747E-03,-0.17208076E-02, 0.48638652E-02,-0.11397893E-01,
     + -0.17629287E-04, 0.35866618E-01,-0.35999727E-02,-0.21764375E-02,
     +  0.10287268E-02,-0.67834837E-04,-0.69727155E-03, 0.86113467E-03,
     +  0.15676366E-03,-0.28114580E-02, 0.84385229E-02, 0.32456819E-01,
     +  0.39928500E-01,-0.14898805E+00,-0.15419769E-01, 0.19524235E-01,
     + -0.34219949E-02, 0.24864037E-03,-0.14452631E-01,-0.50192349E-01,
     + -0.95905527E-03,-0.22331211E-02,-0.35205367E-03, 0.11604355E-02,
     +  0.52138411E-01,-0.16860457E-02,-0.15305842E-02,-0.10912696E-01,
     + -0.64020775E-01,-0.51972526E-02,
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
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      y_sl_q1ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)            *x42    
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)    *x22            
     8  +coeff(  8)    *x22    *x41    
      y_sl_q1ex   =y_sl_q1ex   
     9  +coeff(  9)        *x33*x41    
     1  +coeff( 10)    *x22    *x42    
     2  +coeff( 11)    *x22*x33        
     3  +coeff( 12)    *x22*x33*x41    
     4  +coeff( 13)    *x24*x31*x41    
     5  +coeff( 14)        *x32        
     6  +coeff( 15)        *x31*x41    
     7  +coeff( 16)                *x52
     8  +coeff( 17)        *x33        
      y_sl_q1ex   =y_sl_q1ex   
     9  +coeff( 18)    *x22*x31        
     1  +coeff( 19)    *x22*x32        
     2  +coeff( 20)    *x22*x31*x41    
     3  +coeff( 21)*x12*x22    *x43*x51
     4  +coeff( 22)*x12*x24    *x45*x51
     5  +coeff( 23)    *x22    *x43    
     6  +coeff( 24)    *x22*x33*x42    
     7  +coeff( 25)    *x22*x35    *x53
     8  +coeff( 26)*x11*x21            
      y_sl_q1ex   =y_sl_q1ex   
     9  +coeff( 27)    *x22*x32*x41    
     1  +coeff( 28)    *x22*x31*x42    
     2  +coeff( 29)    *x21*x31*x44    
     3  +coeff( 30)    *x22*x34        
     4  +coeff( 31)*x11    *x31*x43    
     5  +coeff( 32)*x11*x22    *x42    
     6  +coeff( 33)    *x22*x31*x44    
     7  +coeff( 34)*x12*x22*x31        
     8  +coeff( 35)*x12*x22    *x41    
      y_sl_q1ex   =y_sl_q1ex   
     9  +coeff( 36)    *x22*x34*x43    
     1  +coeff( 37)    *x22*x33*x44    
     2  +coeff( 38)*x12*x22    *x41*x51
c
      return
      end
      function p_sl_q1ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff( 57)
      data ncoeff/ 56/
      data avdat/  0.4267549E-02/
      data xmin/
     1 -0.19937E-02,-0.54458E-01,-0.19935E-01,-0.29994E-01,-0.49883E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19971E-02, 0.54459E-01, 0.19987E-01, 0.29928E-01, 0.49858E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.34473275E-02,-0.99775416E-05, 0.93665598E-02, 0.69449626E-01,
     + -0.92471782E-02, 0.28652630E-02,-0.26699525E-03,-0.20499381E-02,
     + -0.62980363E-02,-0.20141857E-03,-0.15989132E-02, 0.20741168E-01,
     + -0.16423305E-03,-0.10873976E-03,-0.46758426E-03, 0.66156418E-03,
     +  0.19266633E-01,-0.16633065E-02,-0.42936037E-03, 0.67457566E-02,
     +  0.23630023E-01, 0.51682116E-02,-0.26409399E-01,-0.10613680E-01,
     + -0.50203071E-03, 0.25507335E-02,-0.28522419E-01,-0.45704059E-03,
     + -0.40300557E-03, 0.19609049E-03, 0.28921105E-03,-0.20700560E-02,
     +  0.96613547E-03,-0.88016903E-02, 0.86790015E-03, 0.68041123E-03,
     + -0.49226014E-02, 0.39621433E-02, 0.26361451E-02,-0.12175565E-02,
     + -0.13432923E-02,-0.15915890E-02,-0.19384342E-02, 0.27017435E-01,
     + -0.13374622E-02,-0.72581728E-03, 0.44428135E-03,-0.53734398E-02,
     + -0.70034922E-03, 0.46743080E-02,-0.18280810E-01,-0.38057152E-01,
     + -0.17203037E-02, 0.30990066E-02, 0.72670439E-02,-0.17331593E-01,
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
      p_sl_q1ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)            *x42    
     8  +coeff(  8)            *x41*x51
      p_sl_q1ex   =p_sl_q1ex   
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)        *x31*x41    
     2  +coeff( 11)    *x22*x31        
     3  +coeff( 12)    *x22    *x42    
     4  +coeff( 13)    *x24*x31*x41    
     5  +coeff( 14)        *x32        
     6  +coeff( 15)        *x31    *x51
     7  +coeff( 16)                *x52
     8  +coeff( 17)    *x22*x31*x41    
      p_sl_q1ex   =p_sl_q1ex   
     9  +coeff( 18)    *x24*x32        
     1  +coeff( 19)    *x24    *x43    
     2  +coeff( 20)*x12*x22    *x43*x51
     3  +coeff( 21)*x12*x24    *x44*x52
     4  +coeff( 22)    *x22*x32        
     5  +coeff( 23)    *x22*x31*x42    
     6  +coeff( 24)    *x22    *x43    
     7  +coeff( 25)    *x24*x32*x41    
     8  +coeff( 26)*x11*x21    *x42*x54
      p_sl_q1ex   =p_sl_q1ex   
     9  +coeff( 27)*x12*x24    *x43*x51
     1  +coeff( 28)        *x32*x41    
     2  +coeff( 29)        *x31*x42    
     3  +coeff( 30)            *x41*x52
     4  +coeff( 31)    *x21    *x43    
     5  +coeff( 32)    *x22*x33        
     6  +coeff( 33)*x11*x21*x31*x41    
     7  +coeff( 34)    *x22*x32*x41    
     8  +coeff( 35)*x11*x22    *x42    
      p_sl_q1ex   =p_sl_q1ex   
     9  +coeff( 36)    *x24    *x41*x51
     1  +coeff( 37)    *x22    *x42*x52
     2  +coeff( 38)        *x32*x42*x52
     3  +coeff( 39)        *x31*x43*x52
     4  +coeff( 40)*x12*x22*x31        
     5  +coeff( 41)*x12*x22    *x41    
     6  +coeff( 42)*x11*x21*x33*x41    
     7  +coeff( 43)*x11*x22    *x43    
     8  +coeff( 44)    *x22*x31*x44    
      p_sl_q1ex   =p_sl_q1ex   
     9  +coeff( 45)    *x21*x31*x44*x51
     1  +coeff( 46)*x12*x23    *x41    
     2  +coeff( 47)*x12*x23        *x51
     3  +coeff( 48)    *x24*x31*x41*x52
     4  +coeff( 49)        *x34    *x54
     5  +coeff( 50)*x12*x22    *x43    
     6  +coeff( 51)    *x22*x34*x43    
     7  +coeff( 52)    *x22*x33*x44    
     8  +coeff( 53)*x12*x22*x33*x41    
      p_sl_q1ex   =p_sl_q1ex   
     9  +coeff( 54)    *x23*x34*x42*x51
     1  +coeff( 55)*x11*x24    *x43*x52
     2  +coeff( 56)*x12*x24*x32*x43    
c
      return
      end
      function l_sl_q1ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff( 22)
      data ncoeff/ 21/
      data avdat/  0.3845210E+01/
      data xmin/
     1 -0.19937E-02,-0.54458E-01,-0.19935E-01,-0.29994E-01,-0.49883E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19971E-02, 0.54459E-01, 0.19987E-01, 0.29928E-01, 0.49858E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.15571980E-02,-0.22431212E-02,-0.51103658E-02, 0.49752663E-04,
     + -0.31424506E-05, 0.38860501E-02, 0.13545714E-04, 0.53542992E-03,
     +  0.34607046E-02,-0.68075307E-04,-0.80641324E-03, 0.30224133E-03,
     + -0.97276573E-03,-0.48435590E-03, 0.18070909E-02,-0.73502318E-03,
     +  0.23529154E-04, 0.60881710E-04, 0.70786940E-04, 0.12864695E-02,
     + -0.35279014E-03,
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
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_sl_q1ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)        *x32        
     8  +coeff(  8)        *x31*x41    
      l_sl_q1ex   =l_sl_q1ex   
     9  +coeff(  9)            *x42    
     1  +coeff( 10)        *x31    *x51
     2  +coeff( 11)            *x41*x51
     3  +coeff( 12)    *x22    *x41    
     4  +coeff( 13)    *x22    *x42    
     5  +coeff( 14)    *x22*x31*x41    
     6  +coeff( 15)    *x22    *x43    
     7  +coeff( 16)    *x22    *x42*x51
     8  +coeff( 17)    *x24*x31*x42    
      l_sl_q1ex   =l_sl_q1ex   
     9  +coeff( 18)                *x52
     1  +coeff( 19)            *x41*x52
     2  +coeff( 20)    *x22*x31*x42    
     3  +coeff( 21)    *x22*x31*x41*x51
c
      return
      end
      function x_sl_dent   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.5098072E+01/
      data xmin/
     1 -0.19937E-02,-0.52993E-01,-0.19935E-01,-0.29994E-01,-0.49883E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19971E-02, 0.54222E-01, 0.19987E-01, 0.29928E-01, 0.49858E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.33852882E-02,-0.12283091E+00,-0.70358696E-02, 0.12478214E-01,
     +  0.34995170E-02,-0.52676926E-03, 0.88614915E-02, 0.14539872E-02,
     + -0.69570215E-03,-0.21989251E-01,-0.55363483E-03, 0.10595961E-02,
     +  0.18473391E-02,-0.61133946E-03,-0.18920865E-02,-0.27615098E-01,
     + -0.13900088E-02, 0.41740181E-03, 0.64538041E-03,-0.11657003E-01,
     +  0.16535744E-01, 0.14346086E-01,-0.11891678E-03, 0.11460640E-02,
     +  0.20457584E-03, 0.85210884E-02, 0.25151217E-01,-0.18307340E-03,
     +  0.24097508E-01,-0.18764210E-02,-0.13559039E-01,-0.22616887E-02,
     + -0.52857613E-02, 0.34116732E-04, 0.97533812E-04,-0.15834483E-03,
     +  0.47488688E-03,-0.11523976E-03,-0.45875035E-03, 0.15131678E-01,
     +  0.32879502E-03,-0.24097445E-02,-0.61995478E-03, 0.12100966E-01,
     +  0.70104445E-02, 0.46007009E-03,-0.95464679E-03,-0.22122382E-04,
     +  0.11916432E-03, 0.13044686E-03, 0.27759323E-04,-0.86620312E-04,
     +  0.34596739E-03,-0.26722639E-03,-0.69992820E-03, 0.83747657E-03,
     +  0.11175343E-04,-0.14310514E-02,-0.14887923E-03,-0.11164482E-03,
     +  0.54846061E-02,-0.89336716E-03, 0.39457993E-03,-0.22412132E-03,
     + -0.19627718E-03, 0.16898735E-03, 0.10528603E-01, 0.23652839E-01,
     + -0.25787388E-02, 0.94541250E-03,-0.14419837E-04, 0.11328339E-02,
     + -0.26969030E-02, 0.29997879E-02, 0.85181855E-02,-0.19527759E-02,
     + -0.68202215E-04,-0.50765218E-03, 0.21937673E-02,-0.25798257E-02,
     + -0.10869855E-02,-0.20890094E-02, 0.75817364E-02, 0.86619006E-02,
     + -0.27070224E-03,-0.28624979E-02, 0.10468995E-02, 0.19306445E-02,
     +  0.21161560E-03,-0.14044269E-02,
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
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      x_sl_dent   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21        *x51
     4  +coeff(  4)    *x21    *x42    
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x21*x31*x41    
     8  +coeff(  8)*x11                
      x_sl_dent   =x_sl_dent   
     9  +coeff(  9)    *x21*x31        
     1  +coeff( 10)    *x21    *x43    
     2  +coeff( 11)    *x21*x33*x42    
     3  +coeff( 12)    *x23            
     4  +coeff( 13)    *x21*x32        
     5  +coeff( 14)    *x21    *x41*x51
     6  +coeff( 15)    *x21*x33        
     7  +coeff( 16)    *x21*x31*x42    
     8  +coeff( 17)    *x21*x34*x41    
      x_sl_dent   =x_sl_dent   
     9  +coeff( 18)    *x22        *x51
     1  +coeff( 19)    *x21        *x52
     2  +coeff( 20)    *x21*x32*x41    
     3  +coeff( 21)    *x21    *x44    
     4  +coeff( 22)    *x21*x33*x43    
     5  +coeff( 23)                *x51
     6  +coeff( 24)    *x21*x34        
     7  +coeff( 25)    *x23*x31*x41    
     8  +coeff( 26)    *x21*x33*x41    
      x_sl_dent   =x_sl_dent   
     9  +coeff( 27)    *x21*x31*x43    
     1  +coeff( 28)*x11            *x51
     2  +coeff( 29)    *x23*x32*x42    
     3  +coeff( 30)    *x23    *x43*x51
     4  +coeff( 31)    *x23*x34*x42    
     5  +coeff( 32)    *x23*x33*x43    
     6  +coeff( 33)*x11*x22*x31*x43*x54
     7  +coeff( 34)            *x42    
     8  +coeff( 35)    *x22    *x41    
      x_sl_dent   =x_sl_dent   
     9  +coeff( 36)    *x21*x31    *x51
     1  +coeff( 37)    *x23    *x41    
     2  +coeff( 38)*x11*x21            
     3  +coeff( 39)    *x23*x32        
     4  +coeff( 40)    *x21*x32*x42    
     5  +coeff( 41)    *x22    *x43    
     6  +coeff( 42)    *x22    *x44    
     7  +coeff( 43)    *x22*x31*x41*x52
     8  +coeff( 44)    *x23    *x44    
      x_sl_dent   =x_sl_dent   
     9  +coeff( 45)    *x21*x32*x44    
     1  +coeff( 46)*x11*x23    *x42    
     2  +coeff( 47)*x11*x22*x32*x41*x51
     3  +coeff( 48)        *x31        
     4  +coeff( 49)            *x41*x51
     5  +coeff( 50)    *x22*x31        
     6  +coeff( 51)        *x31*x42    
     7  +coeff( 52)            *x43    
     8  +coeff( 53)    *x23*x31        
      x_sl_dent   =x_sl_dent   
     9  +coeff( 54)    *x22*x31*x41    
     1  +coeff( 55)    *x21    *x42*x51
     2  +coeff( 56)    *x21    *x41*x52
     3  +coeff( 57)        *x31    *x53
     4  +coeff( 58)    *x23    *x41*x51
     5  +coeff( 59)    *x21*x32    *x52
     6  +coeff( 60)    *x24*x32        
     7  +coeff( 61)    *x23    *x43    
     8  +coeff( 62)    *x23    *x41*x52
      x_sl_dent   =x_sl_dent   
     9  +coeff( 63)    *x21*x32*x41*x52
     1  +coeff( 64)    *x23        *x53
     2  +coeff( 65)            *x43*x53
     3  +coeff( 66)    *x21    *x41*x54
     4  +coeff( 67)    *x21*x34*x42    
     5  +coeff( 68)    *x23*x31*x43    
     6  +coeff( 69)    *x23*x31*x42*x51
     7  +coeff( 70)    *x21*x31*x44*x51
     8  +coeff( 71)    *x22*x31*x41*x53
      x_sl_dent   =x_sl_dent   
     9  +coeff( 72)*x11*x22*x31*x41    
     1  +coeff( 73)    *x24*x31*x43    
     2  +coeff( 74)    *x23*x32*x43    
     3  +coeff( 75)    *x23*x31*x44    
     4  +coeff( 76)    *x21*x33*x44    
     5  +coeff( 77)    *x23*x34    *x51
     6  +coeff( 78)    *x23*x33    *x52
     7  +coeff( 79)    *x22*x33*x41*x52
     8  +coeff( 80)    *x22*x31*x43*x52
      x_sl_dent   =x_sl_dent   
     9  +coeff( 81)    *x21*x33*x41*x53
     1  +coeff( 82)    *x21*x32*x41*x54
     2  +coeff( 83)    *x23*x32*x43*x51
     3  +coeff( 84)    *x23*x31*x44*x51
     4  +coeff( 85)*x11*x22*x31    *x52
     5  +coeff( 86)    *x22*x33*x42*x52
     6  +coeff( 87)    *x21*x34*x42*x52
     7  +coeff( 88)    *x22*x31*x42*x54
     8  +coeff( 89)*x11*x24*x32        
      x_sl_dent   =x_sl_dent   
     9  +coeff( 90)*x11*x22*x33*x41    
c
      return
      end
      function t_sl_dent   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff( 40)
      data ncoeff/ 39/
      data avdat/  0.1298843E+01/
      data xmin/
     1 -0.19937E-02,-0.52993E-01,-0.19935E-01,-0.29994E-01,-0.49883E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19971E-02, 0.54222E-01, 0.19987E-01, 0.29928E-01, 0.49858E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.14487928E-03, 0.55103965E-01,-0.92326518E-03, 0.42479425E-02,
     +  0.31344031E-02,-0.92247697E-02, 0.66114980E-03, 0.10725661E-02,
     + -0.73830551E-02,-0.92197512E-03, 0.90285984E-03,-0.15562925E-02,
     +  0.71934040E-03,-0.13807388E-02,-0.69344713E-03, 0.12515081E-03,
     + -0.39176352E-03,-0.16291604E-02, 0.66802329E-02, 0.61936053E-02,
     + -0.17623348E-02, 0.92951199E-02, 0.23529206E-02, 0.42690593E-02,
     +  0.22466351E-04,-0.32010273E-03,-0.11800952E-03,-0.67526969E-03,
     +  0.23958943E-03, 0.72042830E-03,-0.56290784E-03, 0.54865270E-02,
     +  0.53422531E-03,-0.25826355E-02, 0.13189856E-02, 0.25492439E-02,
     + -0.84137387E-03,-0.15297454E-02, 0.36744301E-02,
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
      t_sl_dent   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)    *x21    *x42    
     7  +coeff(  7)                *x51
     8  +coeff(  8)    *x21    *x41    
      t_sl_dent   =t_sl_dent   
     9  +coeff(  9)    *x21*x31*x41    
     1  +coeff( 10)    *x23*x31        
     2  +coeff( 11)    *x21*x31        
     3  +coeff( 12)    *x21*x32        
     4  +coeff( 13)    *x21    *x41*x51
     5  +coeff( 14)    *x23    *x42    
     6  +coeff( 15)        *x31*x41    
     7  +coeff( 16)    *x22    *x41    
     8  +coeff( 17)    *x21        *x52
      t_sl_dent   =t_sl_dent   
     9  +coeff( 18)    *x23    *x41    
     1  +coeff( 19)    *x21    *x43    
     2  +coeff( 20)    *x23*x32*x41    
     3  +coeff( 21)    *x24    *x42    
     4  +coeff( 22)    *x23*x31*x42    
     5  +coeff( 23)    *x21*x33*x42    
     6  +coeff( 24)    *x23    *x43    
     7  +coeff( 25)        *x31    *x51
     8  +coeff( 26)            *x41*x51
      t_sl_dent   =t_sl_dent   
     9  +coeff( 27)*x11*x21            
     1  +coeff( 28)    *x23            
     2  +coeff( 29)        *x31*x41*x51
     3  +coeff( 30)            *x42*x51
     4  +coeff( 31)    *x24            
     5  +coeff( 32)    *x21*x31*x42    
     6  +coeff( 33)    *x22    *x41*x51
     7  +coeff( 34)    *x21    *x44    
     8  +coeff( 35)    *x23*x33        
      t_sl_dent   =t_sl_dent   
     9  +coeff( 36)    *x21*x34*x41    
     1  +coeff( 37)    *x22    *x42*x52
     2  +coeff( 38)    *x21*x33*x43    
     3  +coeff( 39)*x12*x23    *x43*x51
c
      return
      end
      function y_sl_dent   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff( 79)
      data ncoeff/ 78/
      data avdat/  0.6207331E-02/
      data xmin/
     1 -0.19937E-02,-0.52993E-01,-0.19935E-01,-0.29994E-01,-0.49883E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19971E-02, 0.54222E-01, 0.19987E-01, 0.29928E-01, 0.49858E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.50725290E-02, 0.62354365E-02, 0.11155653E+00,-0.18051835E-01,
     + -0.10814343E-02, 0.14273601E-02, 0.56155524E-02, 0.24179493E-02,
     +  0.17261973E-01, 0.59649111E-02,-0.15267462E-02,-0.18619681E-01,
     +  0.11546015E-02, 0.46772597E-03, 0.29951965E-01, 0.27431636E-02,
     + -0.39394698E-02,-0.94572146E-03,-0.27866188E-01, 0.42333258E-02,
     +  0.17147033E-02,-0.54809090E-03, 0.18547975E-01,-0.35327582E-02,
     + -0.41410946E-02,-0.39842149E-03, 0.27255872E-01,-0.12130269E-02,
     + -0.36585387E-01, 0.52078497E-02,-0.44662273E-03,-0.11707095E-01,
     + -0.62040586E-01, 0.34223506E-03,-0.28794032E-03,-0.83190965E-03,
     +  0.45141336E-03, 0.58599217E-02, 0.97022485E-03,-0.22060115E-01,
     + -0.42179629E-01, 0.26443405E-02,-0.15919370E-02, 0.25847820E-02,
     +  0.30362867E-02, 0.13762055E-01,-0.11222024E-02,-0.52558124E-03,
     +  0.96371950E-03, 0.45945673E-03, 0.25535440E-02,-0.36944717E-03,
     +  0.34340352E-03, 0.38302090E-01, 0.17839551E-01,-0.48358816E-04,
     + -0.44560698E-02,-0.82855746E-02, 0.22449382E-02,-0.39791125E-02,
     + -0.13983664E-02,-0.12275262E-02,-0.15363558E-02,-0.20083464E-02,
     +  0.63284084E-01, 0.18567578E+00, 0.25846890E+00, 0.12903234E+00,
     + -0.32589615E-02,-0.13960409E-02, 0.93439355E-03, 0.59382808E-02,
     + -0.38751515E-02, 0.79018399E-02, 0.24064798E-02,-0.20806745E-02,
     +  0.26198060E-02,-0.88207546E-03,
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
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      y_sl_dent   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)            *x42    
     6  +coeff(  6)    *x21*x31        
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)        *x31    *x51
      y_sl_dent   =y_sl_dent   
     9  +coeff(  9)            *x41*x51
     1  +coeff( 10)    *x22            
     2  +coeff( 11)                *x52
     3  +coeff( 12)    *x22    *x41    
     4  +coeff( 13)    *x22        *x51
     5  +coeff( 14)        *x33*x41    
     6  +coeff( 15)    *x22    *x42    
     7  +coeff( 16)    *x22    *x41*x51
     8  +coeff( 17)    *x22*x33        
      y_sl_dent   =y_sl_dent   
     9  +coeff( 18)    *x24*x31        
     1  +coeff( 19)    *x22*x33*x41    
     2  +coeff( 20)    *x24*x31*x41    
     3  +coeff( 21)    *x24    *x42    
     4  +coeff( 22)    *x24*x34        
     5  +coeff( 23)*x12*x23    *x42*x52
     6  +coeff( 24)            *x43    
     7  +coeff( 25)    *x22*x31        
     8  +coeff( 26)    *x21    *x41*x51
      y_sl_dent   =y_sl_dent   
     9  +coeff( 27)    *x22*x31*x41    
     1  +coeff( 28)    *x24*x32        
     2  +coeff( 29)    *x22    *x45    
     3  +coeff( 30)    *x24*x31*x42    
     4  +coeff( 31)*x12*x21    *x41*x51
     5  +coeff( 32)    *x22*x33*x44    
     6  +coeff( 33)*x12*x24    *x45*x51
     7  +coeff( 34)    *x21            
     8  +coeff( 35)        *x32        
      y_sl_dent   =y_sl_dent   
     9  +coeff( 36)            *x41*x52
     1  +coeff( 37)*x11*x21            
     2  +coeff( 38)    *x22*x32        
     3  +coeff( 39)            *x44*x51
     4  +coeff( 40)    *x22*x32*x41    
     5  +coeff( 41)    *x22*x31*x42    
     6  +coeff( 42)    *x22*x33    *x51
     7  +coeff( 43)*x11*x22    *x42    
     8  +coeff( 44)    *x22*x32*x41*x53
      y_sl_dent   =y_sl_dent   
     9  +coeff( 45)*x12*x21*x31*x41*x52
     1  +coeff( 46)*x11*x24    *x44    
     2  +coeff( 47)        *x31*x41    
     3  +coeff( 48)    *x21        *x51
     4  +coeff( 49)    *x21    *x42    
     5  +coeff( 50)    *x21*x31*x41*x51
     6  +coeff( 51)            *x45    
     7  +coeff( 52)    *x21*x31*x44    
     8  +coeff( 53)*x11*x22        *x51
      y_sl_dent   =y_sl_dent   
     9  +coeff( 54)    *x22*x31*x43    
     1  +coeff( 55)    *x22    *x44    
     2  +coeff( 56)        *x32*x41*x53
     3  +coeff( 57)*x11*x22*x31    *x51
     4  +coeff( 58)    *x24    *x41*x51
     5  +coeff( 59)*x12*x21*x31*x41    
     6  +coeff( 60)*x11*x23    *x42    
     7  +coeff( 61)*x11*x21*x31*x41*x52
     8  +coeff( 62)*x12*x22*x31        
      y_sl_dent   =y_sl_dent   
     9  +coeff( 63)*x12*x22    *x41    
     1  +coeff( 64)*x11*x24    *x41    
     2  +coeff( 65)    *x22*x35*x41    
     3  +coeff( 66)    *x22*x34*x42    
     4  +coeff( 67)    *x22*x33*x43    
     5  +coeff( 68)    *x22*x32*x44    
     6  +coeff( 69)    *x21*x35*x41*x51
     7  +coeff( 70)        *x35*x41*x52
     8  +coeff( 71)*x12*x22        *x51
      y_sl_dent   =y_sl_dent   
     9  +coeff( 72)*x11*x22*x33    *x51
     1  +coeff( 73)    *x24*x33    *x51
     2  +coeff( 74)    *x24*x32*x41*x51
     3  +coeff( 75)*x12*x22*x32        
     4  +coeff( 76)*x12*x23    *x41    
     5  +coeff( 77)*x11*x21*x34*x42    
     6  +coeff( 78)*x12*x24            
c
      return
      end
      function p_sl_dent   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.1674092E-02/
      data xmin/
     1 -0.19937E-02,-0.52993E-01,-0.19935E-01,-0.29994E-01,-0.49883E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19971E-02, 0.54222E-01, 0.19987E-01, 0.29928E-01, 0.49858E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.12783958E-02,-0.71728537E-04,-0.67677535E-02,-0.24278637E-01,
     +  0.17420255E-02,-0.57612685E-03,-0.19193348E-03,-0.35208194E-02,
     +  0.74495608E-03, 0.46113455E-02,-0.62937994E-03, 0.60595310E-03,
     +  0.39905266E-03, 0.13165596E-02,-0.40156567E-02,-0.33129894E-03,
     + -0.11551188E-02, 0.61766128E-03,-0.59589063E-03,-0.33153179E-02,
     + -0.36849218E-04,-0.48697902E-05, 0.33237180E-03, 0.47176669E-03,
     + -0.38416992E-03, 0.43550399E-03,-0.64189947E-03,-0.10310147E-02,
     + -0.70940241E-05, 0.40337504E-03, 0.28628316E-02,-0.86431438E-03,
     + -0.78663824E-03, 0.59481109E-02,-0.13752303E-03, 0.49526094E-04,
     + -0.81151857E-05, 0.26832597E-03, 0.32410887E-03,-0.90935675E-04,
     +  0.23820247E-03,-0.71634079E-03, 0.27212081E-03,-0.21703290E-03,
     + -0.10143255E-02, 0.70333271E-03, 0.20901829E-02, 0.42972367E-03,
     +  0.12230652E-03,-0.23249497E-03,
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
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      p_sl_dent   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)    *x21    *x41    
      p_sl_dent   =p_sl_dent   
     9  +coeff(  9)        *x31    *x51
     1  +coeff( 10)            *x41*x51
     2  +coeff( 11)                *x52
     3  +coeff( 12)    *x21        *x51
     4  +coeff( 13)    *x22*x31        
     5  +coeff( 14)    *x23    *x41    
     6  +coeff( 15)    *x22    *x42    
     7  +coeff( 16)    *x24*x31*x41    
     8  +coeff( 17)    *x22    *x41    
      p_sl_dent   =p_sl_dent   
     9  +coeff( 18)    *x22        *x51
     1  +coeff( 19)    *x21    *x41*x51
     2  +coeff( 20)    *x22*x31*x41    
     3  +coeff( 21)    *x24*x32        
     4  +coeff( 22)        *x32        
     5  +coeff( 23)        *x31*x42    
     6  +coeff( 24)            *x42*x51
     7  +coeff( 25)            *x41*x52
     8  +coeff( 26)    *x21    *x43    
      p_sl_dent   =p_sl_dent   
     9  +coeff( 27)    *x22    *x41*x51
     1  +coeff( 28)    *x23    *x42    
     2  +coeff( 29)    *x21*x32*x42    
     3  +coeff( 30)    *x23*x33*x41    
     4  +coeff( 31)    *x24    *x43    
     5  +coeff( 32)*x12*x21    *x43*x52
     6  +coeff( 33)*x12*x21    *x42*x53
     7  +coeff( 34)*x12*x23    *x43*x54
     8  +coeff( 35)    *x23            
      p_sl_dent   =p_sl_dent   
     9  +coeff( 36)*x11        *x41    
     1  +coeff( 37)    *x21*x31*x41    
     2  +coeff( 38)        *x32*x41    
     3  +coeff( 39)        *x31*x41*x51
     4  +coeff( 40)        *x31    *x52
     5  +coeff( 41)    *x23*x31        
     6  +coeff( 42)    *x22*x32        
     7  +coeff( 43)    *x21*x31*x42    
     8  +coeff( 44)    *x23        *x51
      p_sl_dent   =p_sl_dent   
     9  +coeff( 45)    *x23*x31*x41    
     1  +coeff( 46)    *x22*x32*x41    
     2  +coeff( 47)    *x22*x31*x42    
     3  +coeff( 48)    *x23    *x41*x51
     4  +coeff( 49)        *x32    *x53
     5  +coeff( 50)    *x24*x31    *x51
c
      return
      end
      function l_sl_dent   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff( 14)
      data ncoeff/ 13/
      data avdat/  0.1076003E+02/
      data xmin/
     1 -0.19937E-02,-0.52993E-01,-0.19935E-01,-0.29994E-01,-0.49883E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19971E-02, 0.54222E-01, 0.19987E-01, 0.29928E-01, 0.49858E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.10449916E-01,-0.23525193E+00, 0.16623868E-01,-0.13967809E-01,
     +  0.40876430E-01,-0.41196244E-02, 0.29445291E-02, 0.54975222E-02,
     +  0.32346409E-01,-0.19543075E-02,-0.67250589E-02, 0.12829898E-01,
     + -0.55691735E-02,
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
c
c                  function
c
      l_sl_dent   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x22            
     4  +coeff(  4)    *x21        *x51
     5  +coeff(  5)    *x23    *x42    
     6  +coeff(  6)            *x41    
     7  +coeff(  7)*x11                
     8  +coeff(  8)            *x42    
      l_sl_dent   =l_sl_dent   
     9  +coeff(  9)    *x23*x31*x41    
     1  +coeff( 10)        *x31        
     2  +coeff( 11)    *x21    *x41    
     3  +coeff( 12)    *x21    *x42    
     4  +coeff( 13)    *x23*x31        
c
      return
      end
      function x_sl_dext   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.6484990E-02/
      data xmin/
     1 -0.19937E-02,-0.52993E-01,-0.19935E-01,-0.29994E-01,-0.49883E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19971E-02, 0.54222E-01, 0.19987E-01, 0.29928E-01, 0.49858E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.10782188E-01, 0.33667040E+00, 0.13345107E+00, 0.44310771E-01,
     + -0.43818362E-01,-0.56163147E-02, 0.12042374E-01, 0.15065764E-02,
     + -0.32480523E-01, 0.20788040E-02,-0.35366470E-02,-0.61770491E-02,
     +  0.33798404E-02,-0.41119945E-02,-0.78934181E-03, 0.67015536E-01,
     + -0.89093028E-02,-0.10469859E-01, 0.35320584E-01, 0.61521644E-03,
     +  0.98457956E-03, 0.87662943E-01,-0.51938217E-01,-0.23381410E-02,
     + -0.48934970E-01,-0.15754901E-01,-0.12009065E-02,-0.15678704E-02,
     +  0.56974166E-02, 0.43494716E-01,-0.61419536E-03, 0.15080574E-02,
     + -0.21322684E-01,-0.65984376E-01,-0.63397473E-03,-0.31154177E-02,
     + -0.70577592E-01, 0.16259057E-01, 0.63383386E-01, 0.17280350E-01,
     +  0.22234404E-03,-0.16139030E-02,-0.19765196E-02, 0.94179128E-03,
     + -0.19036452E-03,-0.34787129E-02,-0.38573068E-01, 0.24198560E-03,
     +  0.23247686E-03,-0.37729554E-01,-0.16172437E-02, 0.99864472E-02,
     + -0.19133873E-02, 0.73156430E-03, 0.99461155E-04,-0.24837983E-03,
     + -0.15459272E-03,-0.53044374E-03,-0.31324266E-02, 0.23886934E-03,
     + -0.23204547E-02, 0.24374724E-03, 0.30290319E-02, 0.86533034E-03,
     +  0.49588969E-02, 0.55986652E-02, 0.10814396E-02,-0.22694678E-03,
     + -0.22017928E-02, 0.19843914E-02,-0.16438188E-01,-0.30829690E-01,
     + -0.15926940E-01, 0.28119194E-02, 0.60882448E-03,-0.84887750E-01,
     + -0.29944981E-01,-0.16861457E-02, 0.29415761E-02,-0.34445697E-02,
     +  0.28549882E-02, 0.99722098E-03,-0.17517297E-01,-0.47773138E-01,
     + -0.47234867E-01, 0.28563140E-01, 0.36329968E-03, 0.60558133E-03,
     + -0.36420999E-02, 0.49112983E-01,
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
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      x_sl_dext   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x21        *x51
     5  +coeff(  5)    *x21    *x42    
     6  +coeff(  6)*x11                
     7  +coeff(  7)    *x22            
     8  +coeff(  8)    *x21    *x41    
      x_sl_dext   =x_sl_dext   
     9  +coeff(  9)    *x21*x31*x41    
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)                *x52
     3  +coeff( 12)    *x21*x32        
     4  +coeff( 13)    *x21    *x41*x51
     5  +coeff( 14)            *x42    
     6  +coeff( 15)    *x23            
     7  +coeff( 16)    *x21    *x43    
     8  +coeff( 17)    *x21*x34*x41    
      x_sl_dext   =x_sl_dext   
     9  +coeff( 18)    *x21*x33*x42    
     1  +coeff( 19)    *x21*x34*x43    
     2  +coeff( 20)            *x41*x51
     3  +coeff( 21)    *x21*x31    *x51
     4  +coeff( 22)    *x21*x31*x42    
     5  +coeff( 23)    *x21    *x44    
     6  +coeff( 24)    *x23    *x42*x51
     7  +coeff( 25)    *x21*x33*x43    
     8  +coeff( 26)    *x21*x32*x44    
      x_sl_dext   =x_sl_dext   
     9  +coeff( 27)        *x31*x41    
     1  +coeff( 28)    *x21        *x52
     2  +coeff( 29)    *x21*x33        
     3  +coeff( 30)    *x21*x32*x41    
     4  +coeff( 31)    *x21*x31*x41*x51
     5  +coeff( 32)    *x23*x31*x41    
     6  +coeff( 33)    *x21*x33*x41    
     7  +coeff( 34)    *x21*x31*x43    
     8  +coeff( 35)            *x42*x53
      x_sl_dext   =x_sl_dext   
     9  +coeff( 36)    *x22    *x42*x52
     1  +coeff( 37)    *x23*x32*x42    
     2  +coeff( 38)    *x23    *x43*x51
     3  +coeff( 39)    *x23*x34*x42    
     4  +coeff( 40)*x11*x22*x31*x43*x54
     5  +coeff( 41)    *x22*x31        
     6  +coeff( 42)    *x24            
     7  +coeff( 43)    *x23    *x41    
     8  +coeff( 44)    *x23        *x51
      x_sl_dext   =x_sl_dext   
     9  +coeff( 45)*x11*x21            
     1  +coeff( 46)    *x23    *x42    
     2  +coeff( 47)    *x21*x32*x42    
     3  +coeff( 48)*x11            *x51
     4  +coeff( 49)    *x22*x31*x41*x51
     5  +coeff( 50)    *x21*x34*x42    
     6  +coeff( 51)    *x23*x32*x41*x51
     7  +coeff( 52)    *x23*x31*x42*x51
     8  +coeff( 53)*x11*x23    *x42    
      x_sl_dext   =x_sl_dext   
     9  +coeff( 54)    *x22*x32*x42*x54
     1  +coeff( 55)            *x41    
     2  +coeff( 56)        *x32        
     3  +coeff( 57)        *x32*x41    
     4  +coeff( 58)    *x23*x31        
     5  +coeff( 59)    *x21    *x41*x52
     6  +coeff( 60)            *x41*x53
     7  +coeff( 61)    *x21*x34        
     8  +coeff( 62)    *x22*x31*x42    
      x_sl_dext   =x_sl_dext   
     9  +coeff( 63)    *x23    *x41*x51
     1  +coeff( 64)    *x21*x32*x41*x51
     2  +coeff( 65)    *x21*x31*x41*x52
     3  +coeff( 66)    *x21    *x42*x52
     4  +coeff( 67)            *x43*x52
     5  +coeff( 68)            *x41*x54
     6  +coeff( 69)    *x24    *x42    
     7  +coeff( 70)    *x22*x32*x42    
     8  +coeff( 71)    *x23    *x43    
      x_sl_dext   =x_sl_dext   
     9  +coeff( 72)    *x21*x32*x43    
     1  +coeff( 73)    *x21*x31*x44    
     2  +coeff( 74)    *x23    *x41*x52
     3  +coeff( 75)            *x43*x53
     4  +coeff( 76)    *x23*x31*x43    
     5  +coeff( 77)    *x23    *x44    
     6  +coeff( 78)    *x21*x31*x44*x51
     7  +coeff( 79)    *x23*x32    *x52
     8  +coeff( 80)    *x23    *x42*x52
      x_sl_dext   =x_sl_dext   
     9  +coeff( 81)    *x22*x31*x42*x52
     1  +coeff( 82)*x11*x22    *x42    
     2  +coeff( 83)    *x23*x33*x42    
     3  +coeff( 84)    *x23*x32*x43    
     4  +coeff( 85)    *x23*x31*x44    
     5  +coeff( 86)    *x21*x33*x44    
     6  +coeff( 87)*x11    *x32    *x52
     7  +coeff( 88)    *x21*x33*x42*x52
     8  +coeff( 89)    *x21*x31*x42*x54
      x_sl_dext   =x_sl_dext   
     9  +coeff( 90)    *x23*x33*x43    
c
      return
      end
      function t_sl_dext   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff( 26)
      data ncoeff/ 25/
      data avdat/  0.5334208E+00/
      data xmin/
     1 -0.19937E-02,-0.52993E-01,-0.19935E-01,-0.29994E-01,-0.49883E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19971E-02, 0.54222E-01, 0.19987E-01, 0.29928E-01, 0.49858E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.24697569E-03,-0.77090755E-01, 0.22669706E-01,-0.33494285E-02,
     + -0.40802327E-02, 0.11348788E-01, 0.92344568E-03,-0.87511970E-03,
     + -0.12146849E-02, 0.82615055E-02,-0.10310973E-02, 0.16475643E-02,
     + -0.85803209E-03, 0.52955630E-03,-0.10092754E-01,-0.57463879E-02,
     + -0.54166168E-02, 0.25373622E-03,-0.23041094E-03,-0.15105777E-03,
     +  0.17667658E-03, 0.66835259E-03,-0.85113924E-02, 0.17077340E-02,
     +  0.35958563E-02,
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
      t_sl_dext   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)    *x21    *x42    
     7  +coeff(  7)*x11                
     8  +coeff(  8)    *x21    *x41    
      t_sl_dext   =t_sl_dext   
     9  +coeff(  9)                *x52
     1  +coeff( 10)    *x21*x31*x41    
     2  +coeff( 11)    *x21*x31        
     3  +coeff( 12)    *x21*x32        
     4  +coeff( 13)    *x21    *x41*x51
     5  +coeff( 14)    *x21        *x52
     6  +coeff( 15)    *x21    *x43    
     7  +coeff( 16)    *x21*x34*x41    
     8  +coeff( 17)    *x21*x33*x42    
      t_sl_dext   =t_sl_dext   
     9  +coeff( 18)    *x21*x34*x43    
     1  +coeff( 19)        *x32        
     2  +coeff( 20)        *x31    *x51
     3  +coeff( 21)*x11*x21            
     4  +coeff( 22)    *x22    *x42    
     5  +coeff( 23)    *x21*x31*x42    
     6  +coeff( 24)    *x23    *x42*x51
     7  +coeff( 25)    *x23    *x44    
c
      return
      end
      function y_sl_dext   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff( 96)
      data ncoeff/ 95/
      data avdat/ -0.1470632E-02/
      data xmin/
     1 -0.19937E-02,-0.52993E-01,-0.19935E-01,-0.29994E-01,-0.49883E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19971E-02, 0.54222E-01, 0.19987E-01, 0.29928E-01, 0.49858E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.22568041E-03,-0.37907962E-01, 0.86745024E-02, 0.29422922E-03,
     + -0.17228726E-01,-0.96255091E-04,-0.71781926E-03,-0.56196726E-02,
     + -0.51183023E-01, 0.86596645E-02, 0.56449048E-01, 0.64197555E-02,
     +  0.64837718E-02,-0.63181883E-02,-0.64563397E-02,-0.83111636E-02,
     +  0.15419806E-03, 0.25936430E-02,-0.79300469E-02,-0.42999309E-01,
     + -0.44573567E-03,-0.55812127E-02,-0.37145088E-02,-0.34041705E-02,
     +  0.64959330E-02,-0.51185250E-03, 0.87444484E-02, 0.32559615E-01,
     +  0.77380999E-02,-0.17483546E-02,-0.69171605E-02,-0.58724708E-02,
     + -0.28929373E-02,-0.70257825E-02, 0.46792785E-02,-0.28080307E-01,
     + -0.13546729E-02, 0.10030115E-01, 0.12609963E-01, 0.62525866E-03,
     + -0.12465767E-02, 0.74963312E-03, 0.70386147E-03, 0.90690274E-02,
     +  0.33080534E-02, 0.33380356E-01, 0.27090386E-02,-0.16838772E-02,
     + -0.21809382E-01, 0.62721595E-02, 0.42585758E-02,-0.16450109E-01,
     +  0.42440128E-02,-0.39130524E-02,-0.10806850E-01, 0.63387100E-02,
     + -0.16944384E-02,-0.44265068E-02,-0.34345747E-04, 0.56823855E-03,
     + -0.59496530E-03, 0.85197733E-03, 0.39605792E-02,-0.41879662E-02,
     + -0.56972411E-02, 0.93307375E-03, 0.18072856E-02, 0.29935595E-02,
     +  0.22441696E-01, 0.21987129E-01, 0.20491506E-02,-0.41450127E-02,
     + -0.18789265E-01, 0.72615286E-02,-0.23509547E-03,-0.55097202E-02,
     + -0.17328311E-01,-0.94365152E-02,-0.10303336E-02,-0.10024737E-01,
     + -0.19644483E-02, 0.14485984E-01, 0.44177944E-03,-0.96088886E-03,
     +  0.14639533E-01,-0.24833700E-01, 0.21456642E-01, 0.12005156E-01,
     + -0.66784066E-02, 0.73202229E-02,-0.71514193E-02,-0.42151228E-01,
     + -0.15878243E-01, 0.26890755E-03,-0.10394468E-03,
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
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      y_sl_dext   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x21            
     5  +coeff(  5)                *x51
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21*x31        
      y_sl_dext   =y_sl_dext   
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)        *x31    *x51
     2  +coeff( 11)            *x41*x51
     3  +coeff( 12)    *x22            
     4  +coeff( 13)    *x21        *x51
     5  +coeff( 14)                *x52
     6  +coeff( 15)    *x21*x31*x41    
     7  +coeff( 16)    *x21    *x42    
     8  +coeff( 17)        *x32    *x51
      y_sl_dext   =y_sl_dext   
     9  +coeff( 18)            *x42*x51
     1  +coeff( 19)    *x22*x31        
     2  +coeff( 20)    *x22    *x41    
     3  +coeff( 21)    *x21*x31    *x51
     4  +coeff( 22)    *x21    *x41*x51
     5  +coeff( 23)            *x41*x52
     6  +coeff( 24)    *x23            
     7  +coeff( 25)    *x22        *x51
     8  +coeff( 26)        *x34        
      y_sl_dext   =y_sl_dext   
     9  +coeff( 27)    *x21    *x43    
     1  +coeff( 28)    *x22    *x42    
     2  +coeff( 29)    *x23    *x41    
     3  +coeff( 30)    *x22*x31    *x51
     4  +coeff( 31)    *x22    *x41*x51
     5  +coeff( 32)        *x31*x44    
     6  +coeff( 33)            *x45    
     7  +coeff( 34)    *x22    *x43    
     8  +coeff( 35)    *x22    *x42*x51
      y_sl_dext   =y_sl_dext   
     9  +coeff( 36)    *x22*x33*x41    
     1  +coeff( 37)    *x24*x31*x41    
     2  +coeff( 38)    *x24    *x44*x51
     3  +coeff( 39)    *x24*x34*x41*x51
     4  +coeff( 40)        *x32        
     5  +coeff( 41)    *x21*x32        
     6  +coeff( 42)        *x31*x41*x51
     7  +coeff( 43)                *x53
     8  +coeff( 44)    *x21*x31*x42    
      y_sl_dext   =y_sl_dext   
     9  +coeff( 45)    *x22*x32        
     1  +coeff( 46)    *x22*x31*x41    
     2  +coeff( 47)    *x23*x31        
     3  +coeff( 48)    *x24            
     4  +coeff( 49)    *x22*x31*x42    
     5  +coeff( 50)    *x22*x31*x41*x51
     6  +coeff( 51)    *x21    *x42*x52
     7  +coeff( 52)    *x22    *x45    
     8  +coeff( 53)*x11*x21    *x42*x52
      y_sl_dext   =y_sl_dext   
     9  +coeff( 54)    *x23*x32*x43    
     1  +coeff( 55)*x11*x23*x31*x43    
     2  +coeff( 56)*x12*x22    *x43*x51
     3  +coeff( 57)        *x32*x41    
     4  +coeff( 58)            *x43    
     5  +coeff( 59)*x11    *x31        
     6  +coeff( 60)*x11        *x41    
     7  +coeff( 61)        *x31    *x52
     8  +coeff( 62)    *x21*x32    *x51
      y_sl_dext   =y_sl_dext   
     9  +coeff( 63)    *x21    *x41*x52
     1  +coeff( 64)    *x21    *x44*x51
     2  +coeff( 65)    *x21    *x43*x52
     3  +coeff( 66)    *x23*x31    *x52
     4  +coeff( 67)*x11*x21*x33*x41    
     5  +coeff( 68)    *x21*x33*x41*x52
     6  +coeff( 69)    *x22*x35*x41    
     7  +coeff( 70)    *x23    *x45    
     8  +coeff( 71)    *x21*x34*x41*x52
      y_sl_dext   =y_sl_dext   
     9  +coeff( 72)*x11*x22    *x44    
     1  +coeff( 73)    *x24    *x44    
     2  +coeff( 74)    *x21    *x44*x53
     3  +coeff( 75)*x12        *x43*x51
     4  +coeff( 76)*x12*x22*x31*x41    
     5  +coeff( 77)    *x24    *x42*x52
     6  +coeff( 78)    *x22*x31*x41*x54
     7  +coeff( 79)*x12*x23    *x41    
     8  +coeff( 80)    *x23*x33*x43    
      y_sl_dext   =y_sl_dext   
     9  +coeff( 81)*x12*x21    *x41*x52
     1  +coeff( 82)    *x24*x31*x44    
     2  +coeff( 83)*x12*x23        *x51
     3  +coeff( 84)*x11*x23        *x53
     4  +coeff( 85)*x11*x23*x32*x42    
     5  +coeff( 86)*x11*x23    *x44    
     6  +coeff( 87)*x11*x24    *x44    
     7  +coeff( 88)*x12*x22    *x44*x51
     8  +coeff( 89)*x11*x23*x34    *x52
      y_sl_dext   =y_sl_dext   
     9  +coeff( 90)*x11*x23*x32*x44*x51
     1  +coeff( 91)*x12*x23*x31    *x53
     2  +coeff( 92)*x12*x24    *x43*x51
     3  +coeff( 93)*x12*x23*x32*x41*x53
     4  +coeff( 94)*x11*x21            
     5  +coeff( 95)*x11            *x51
c
      return
      end
      function p_sl_dext   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff( 75)
      data ncoeff/ 74/
      data avdat/ -0.1266671E-02/
      data xmin/
     1 -0.19937E-02,-0.52993E-01,-0.19935E-01,-0.29994E-01,-0.49883E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19971E-02, 0.54222E-01, 0.19987E-01, 0.29928E-01, 0.49858E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.83914382E-03, 0.68947412E-06,-0.92731258E-02,-0.13882218E-01,
     + -0.12287310E-02,-0.42704160E-05, 0.68654225E-03,-0.11905688E-02,
     + -0.10466918E-01, 0.40355801E-04, 0.14939836E-02, 0.13826580E-02,
     +  0.10284891E-01,-0.12689906E-02,-0.53304667E-03,-0.64235921E-02,
     + -0.92192367E-03, 0.13811761E-02,-0.13201935E-02,-0.18833440E-02,
     + -0.67356479E-03,-0.12065597E-02, 0.43887438E-03, 0.21289401E-02,
     +  0.24073746E-02,-0.16885644E-02,-0.66345627E-03, 0.15903334E-03,
     +  0.17239292E-03,-0.37842756E-03,-0.43475899E-03, 0.60015405E-03,
     +  0.18153118E-02,-0.32511593E-02,-0.26724942E-03, 0.10802624E-03,
     +  0.13499848E-03, 0.10017233E-03, 0.55912556E-03, 0.21852241E-02,
     + -0.28637471E-03,-0.35555387E-03,-0.45018567E-03, 0.13084275E-02,
     +  0.15957933E-03, 0.72455371E-03,-0.13550764E-03, 0.19201898E-02,
     +  0.26196642E-02, 0.10740687E-02, 0.75768403E-04, 0.63017615E-04,
     +  0.32305140E-05,-0.22878647E-04,-0.87796916E-05,-0.30457179E-03,
     + -0.16687611E-03, 0.46233770E-04,-0.11690089E-02,-0.87038294E-03,
     +  0.59849944E-03,-0.12442489E-03, 0.12285236E-02,-0.15559005E-03,
     + -0.37407374E-03, 0.59774169E-03, 0.20892306E-02,-0.41927039E-03,
     +  0.26906759E-03,-0.11739894E-02, 0.74470881E-03, 0.92290645E-03,
     +  0.18480035E-02,-0.71337621E-03,
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
      p_sl_dext   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)*x11                
     7  +coeff(  7)    *x22            
     8  +coeff(  8)    *x21*x31        
      p_sl_dext   =p_sl_dext   
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)            *x42    
     2  +coeff( 11)    *x21        *x51
     3  +coeff( 12)        *x31    *x51
     4  +coeff( 13)            *x41*x51
     5  +coeff( 14)                *x52
     6  +coeff( 15)    *x22*x31        
     7  +coeff( 16)    *x22    *x41    
     8  +coeff( 17)            *x43    
      p_sl_dext   =p_sl_dext   
     9  +coeff( 18)    *x22        *x51
     1  +coeff( 19)    *x21    *x41*x51
     2  +coeff( 20)    *x22    *x41*x51
     3  +coeff( 21)    *x23            
     4  +coeff( 22)    *x21    *x42    
     5  +coeff( 23)            *x42*x51
     6  +coeff( 24)    *x23    *x41    
     7  +coeff( 25)    *x21    *x43    
     8  +coeff( 26)    *x23*x31*x41    
      p_sl_dext   =p_sl_dext   
     9  +coeff( 27)    *x21*x31*x41    
     1  +coeff( 28)        *x31*x41*x51
     2  +coeff( 29)    *x21        *x52
     3  +coeff( 30)            *x41*x52
     4  +coeff( 31)    *x24            
     5  +coeff( 32)    *x21*x32*x41    
     6  +coeff( 33)    *x22    *x42*x51
     7  +coeff( 34)    *x23    *x44    
     8  +coeff( 35)    *x21*x32        
      p_sl_dext   =p_sl_dext   
     9  +coeff( 36)*x11        *x41    
     1  +coeff( 37)        *x32*x41    
     2  +coeff( 38)                *x53
     3  +coeff( 39)    *x23*x31        
     4  +coeff( 40)    *x21*x31*x42    
     5  +coeff( 41)    *x23        *x51
     6  +coeff( 42)    *x22*x31    *x51
     7  +coeff( 43)    *x21    *x42*x51
     8  +coeff( 44)    *x22*x31*x41*x51
      p_sl_dext   =p_sl_dext   
     9  +coeff( 45)    *x21*x31    *x53
     1  +coeff( 46)    *x24*x31*x41    
     2  +coeff( 47)    *x22    *x43*x51
     3  +coeff( 48)*x11*x22*x31*x44    
     4  +coeff( 49)    *x23    *x44*x52
     5  +coeff( 50)    *x23*x31*x41*x54
     6  +coeff( 51)        *x33        
     7  +coeff( 52)        *x32    *x51
     8  +coeff( 53)*x11*x21    *x41    
      p_sl_dext   =p_sl_dext   
     9  +coeff( 54)*x11*x21        *x51
     1  +coeff( 55)    *x21*x32    *x51
     2  +coeff( 56)    *x21*x31*x41*x51
     3  +coeff( 57)            *x42*x52
     4  +coeff( 58)*x12*x21            
     5  +coeff( 59)    *x22*x31*x42    
     6  +coeff( 60)    *x21*x31*x43    
     7  +coeff( 61)    *x23    *x41*x51
     8  +coeff( 62)*x12*x21    *x41    
      p_sl_dext   =p_sl_dext   
     9  +coeff( 63)    *x22*x32*x42    
     1  +coeff( 64)    *x22*x33    *x51
     2  +coeff( 65)    *x21    *x43*x52
     3  +coeff( 66)    *x23*x33*x41    
     4  +coeff( 67)    *x22*x31*x44    
     5  +coeff( 68)        *x33*x44    
     6  +coeff( 69)*x11*x24    *x42    
     7  +coeff( 70)    *x22*x34*x42    
     8  +coeff( 71)*x11*x23    *x43    
      p_sl_dext   =p_sl_dext   
     9  +coeff( 72)    *x24*x32*x41*x51
     1  +coeff( 73)*x11*x22*x32*x43    
     2  +coeff( 74)*x12*x24    *x41*x53
c
      return
      end
      function l_sl_dext   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff( 14)
      data ncoeff/ 13/
      data avdat/  0.1734474E+02/
      data xmin/
     1 -0.19937E-02,-0.52993E-01,-0.19935E-01,-0.29994E-01,-0.49883E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19971E-02, 0.54222E-01, 0.19987E-01, 0.29928E-01, 0.49858E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.57046423E-02, 0.47976759E+00, 0.99055469E-01,-0.70931041E-02,
     +  0.30653503E-01, 0.45108709E-01,-0.67118756E-01, 0.18129779E-01,
     + -0.56265317E-01, 0.10304668E-01,-0.10142681E-01,-0.78506902E-01,
     +  0.42579759E-01,
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
      x44 = x43*x4
      x51 = x5
c
c                  function
c
      l_sl_dext   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)    *x21    *x42    
     8  +coeff(  8)    *x21    *x41    
      l_sl_dext   =l_sl_dext   
     9  +coeff(  9)    *x21*x31*x41    
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)    *x21*x32        
     3  +coeff( 12)    *x23    *x44    
     4  +coeff( 13)    *x23    *x43*x51
c
      return
      end
      function x_sl_q3en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.8498874E-02/
      data xmin/
     1 -0.19937E-02,-0.52993E-01,-0.19935E-01,-0.29994E-01,-0.49883E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19971E-02, 0.54222E-01, 0.19987E-01, 0.29928E-01, 0.49858E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.12528533E-01, 0.21310301E+00, 0.13863331E+00, 0.18989410E-01,
     +  0.34921251E-01,-0.27756438E-01,-0.59906091E-02,-0.19021343E-01,
     + -0.38925454E-02, 0.46702445E-01, 0.11397181E-02,-0.36125679E-02,
     +  0.61328754E-01,-0.47122361E-02, 0.32745691E-02, 0.27032971E-01,
     + -0.70448930E-03,-0.89959690E-03, 0.84559171E-03, 0.39453362E-02,
     + -0.41888250E-03, 0.30665139E-02, 0.78580109E-03,-0.16639654E-02,
     + -0.18107881E-02,-0.51760586E-03,-0.52548945E-02,-0.41586466E-01,
     + -0.75549167E-03,-0.12185035E-02,-0.21157391E-01,-0.32807935E-01,
     +  0.13934145E-01,-0.13140525E-03, 0.29373434E-03, 0.30576787E-03,
     +  0.41142240E-03,-0.16557607E-02,-0.27946080E-02,-0.17337770E-02,
     + -0.82034282E-02,-0.64807959E-01, 0.93638984E-04,-0.26182231E-03,
     +  0.18802745E-03,-0.16120194E-01,-0.95399478E-02,-0.21889750E-01,
     + -0.89899143E-02,-0.22106681E-02,-0.19131722E-01, 0.15101366E-01,
     +  0.38998867E-02, 0.19927446E-01,-0.12918313E-01,-0.44352976E-02,
     +  0.64946078E-04, 0.45410224E-03, 0.65510429E-03,-0.29849841E-05,
     + -0.57043624E-03, 0.77380508E-03, 0.13085338E-03,-0.61864121E-03,
     +  0.34131075E-03,-0.12879570E-01,-0.33544876E-01,-0.10316843E-02,
     +  0.71342353E-03,-0.26744721E-03,-0.40960168E-02,-0.30501818E-02,
     + -0.74396622E-02,-0.28465025E-02,-0.23695778E-02, 0.10496109E-02,
     +  0.14888379E-02, 0.11238386E-03,-0.14158176E-02,-0.16345001E-03,
     + -0.19886201E-01, 0.60112979E-02,-0.43250252E-02,-0.59706951E-02,
     +  0.27365759E-02, 0.19458554E-02, 0.14609244E-02, 0.40541936E-03,
     + -0.30586987E-02,-0.23292257E-03,
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
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      x_sl_q3en   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)    *x21    *x42    
     7  +coeff(  7)                *x52
     8  +coeff(  8)    *x21*x31*x41    
      x_sl_q3en   =x_sl_q3en   
     9  +coeff(  9)*x11                
     1  +coeff( 10)    *x21    *x43    
     2  +coeff( 11)    *x21*x31        
     3  +coeff( 12)            *x42    
     4  +coeff( 13)    *x21*x31*x42    
     5  +coeff( 14)    *x21*x32        
     6  +coeff( 15)    *x21    *x41*x51
     7  +coeff( 16)    *x21*x32*x41    
     8  +coeff( 17)    *x24        *x51
      x_sl_q3en   =x_sl_q3en   
     9  +coeff( 18)        *x31*x41    
     1  +coeff( 19)            *x41*x51
     2  +coeff( 20)    *x21*x33        
     3  +coeff( 21)*x11*x21            
     4  +coeff( 22)    *x22        *x51
     5  +coeff( 23)    *x21*x31    *x51
     6  +coeff( 24)    *x24            
     7  +coeff( 25)    *x22    *x42    
     8  +coeff( 26)    *x21    *x42*x51
      x_sl_q3en   =x_sl_q3en   
     9  +coeff( 27)    *x23*x31*x41    
     1  +coeff( 28)    *x21    *x44    
     2  +coeff( 29)    *x21        *x54
     3  +coeff( 30)    *x24*x31*x41    
     4  +coeff( 31)    *x21*x34*x42    
     5  +coeff( 32)    *x21*x33*x43    
     6  +coeff( 33)    *x23    *x43*x51
     7  +coeff( 34)        *x32        
     8  +coeff( 35)    *x22*x31        
      x_sl_q3en   =x_sl_q3en   
     9  +coeff( 36)    *x22    *x41    
     1  +coeff( 37)                *x53
     2  +coeff( 38)    *x23    *x41    
     3  +coeff( 39)    *x22*x31*x41    
     4  +coeff( 40)    *x21*x31*x41*x51
     5  +coeff( 41)    *x23    *x42    
     6  +coeff( 42)    *x21*x31*x43    
     7  +coeff( 43)*x11            *x51
     8  +coeff( 44)        *x33*x41*x51
      x_sl_q3en   =x_sl_q3en   
     9  +coeff( 45)    *x21*x32    *x52
     1  +coeff( 46)    *x21*x33*x42    
     2  +coeff( 47)    *x23    *x43    
     3  +coeff( 48)    *x21*x32*x43    
     4  +coeff( 49)    *x21*x31*x44    
     5  +coeff( 50)    *x21    *x43*x52
     6  +coeff( 51)    *x23*x32*x42    
     7  +coeff( 52)    *x23*x31*x42*x51
     8  +coeff( 53)    *x23*x31*x42*x52
      x_sl_q3en   =x_sl_q3en   
     9  +coeff( 54)    *x23*x34*x42    
     1  +coeff( 55)    *x23*x31*x43*x53
     2  +coeff( 56)*x11*x22    *x44*x54
     3  +coeff( 57)            *x41    
     4  +coeff( 58)    *x21    *x41    
     5  +coeff( 59)    *x23            
     6  +coeff( 60)    *x23*x31        
     7  +coeff( 61)    *x22*x32        
     8  +coeff( 62)    *x23        *x51
      x_sl_q3en   =x_sl_q3en   
     9  +coeff( 63)    *x22*x31    *x51
     1  +coeff( 64)    *x21*x32    *x51
     2  +coeff( 65)    *x22*x33        
     3  +coeff( 66)    *x21*x33*x41    
     4  +coeff( 67)    *x21*x32*x42    
     5  +coeff( 68)        *x32*x42*x51
     6  +coeff( 69)            *x43*x52
     7  +coeff( 70)            *x41*x54
     8  +coeff( 71)    *x21*x34*x41    
      x_sl_q3en   =x_sl_q3en   
     9  +coeff( 72)    *x24    *x42    
     1  +coeff( 73)    *x23*x31*x42    
     2  +coeff( 74)    *x23    *x42*x51
     3  +coeff( 75)    *x21*x32*x42*x51
     4  +coeff( 76)    *x22*x31*x41*x52
     5  +coeff( 77)    *x21*x32*x41*x52
     6  +coeff( 78)            *x43*x53
     7  +coeff( 79)    *x23*x34        
     8  +coeff( 80)*x11*x21    *x42    
      x_sl_q3en   =x_sl_q3en   
     9  +coeff( 81)    *x21*x32*x44    
     1  +coeff( 82)    *x23*x32*x41*x51
     2  +coeff( 83)    *x21*x32*x43*x51
     3  +coeff( 84)    *x21*x31*x44*x51
     4  +coeff( 85)    *x21*x33*x41*x52
     5  +coeff( 86)    *x22*x31*x42*x52
     6  +coeff( 87)    *x21    *x42*x54
     7  +coeff( 88)*x11*x23    *x41    
     8  +coeff( 89)    *x24    *x44    
      x_sl_q3en   =x_sl_q3en   
     9  +coeff( 90)*x11*x23        *x51
c
      return
      end
      function t_sl_q3en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff( 26)
      data ncoeff/ 25/
      data avdat/ -0.1245997E-02/
      data xmin/
     1 -0.19937E-02,-0.52993E-01,-0.19935E-01,-0.29994E-01,-0.49883E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19971E-02, 0.54222E-01, 0.19987E-01, 0.29928E-01, 0.49858E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.24844197E-03,-0.77568568E-01, 0.22737240E-01,-0.47574057E-02,
     + -0.52620126E-02, 0.92496991E-03,-0.67676668E-03,-0.13958720E-02,
     +  0.12484854E-01, 0.12176114E-03,-0.94116590E-03, 0.88893203E-02,
     +  0.17109242E-02,-0.10730088E-01, 0.81195455E-03,-0.57167299E-02,
     + -0.58307406E-02,-0.75666088E-03,-0.13783836E-03, 0.22537365E-03,
     + -0.72192395E-03, 0.20010537E-02,-0.88094855E-02,-0.45372662E-03,
     +  0.13626391E-02,
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
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      t_sl_q3en   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)*x11                
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)                *x52
      t_sl_q3en   =t_sl_q3en   
     9  +coeff(  9)    *x21    *x42    
     1  +coeff( 10)    *x21*x33*x41    
     2  +coeff( 11)    *x21*x31        
     3  +coeff( 12)    *x21*x31*x41    
     4  +coeff( 13)    *x21*x32        
     5  +coeff( 14)    *x21    *x43    
     6  +coeff( 15)    *x23        *x52
     7  +coeff( 16)    *x21*x34*x41    
     8  +coeff( 17)    *x21*x33*x42    
      t_sl_q3en   =t_sl_q3en   
     9  +coeff( 18)    *x21*x34*x43    
     1  +coeff( 19)            *x42    
     2  +coeff( 20)*x11*x21            
     3  +coeff( 21)    *x21    *x41*x51
     4  +coeff( 22)    *x22    *x42    
     5  +coeff( 23)    *x21*x31*x42    
     6  +coeff( 24)    *x24        *x51
     7  +coeff( 25)    *x24*x31*x41    
c
      return
      end
      function y_sl_q3en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff( 96)
      data ncoeff/ 95/
      data avdat/ -0.2633865E-02/
      data xmin/
     1 -0.19937E-02,-0.52993E-01,-0.19935E-01,-0.29994E-01,-0.49883E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19971E-02, 0.54222E-01, 0.19987E-01, 0.29928E-01, 0.49858E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.10038784E-02,-0.48838757E-01,-0.53431843E-02, 0.18930522E-03,
     + -0.18967500E-01,-0.63207513E-03,-0.95709943E-03,-0.57505555E-02,
     + -0.61528206E-01, 0.11174984E-01, 0.69766916E-01, 0.70063937E-02,
     +  0.88256244E-02,-0.77441046E-02,-0.67816819E-02,-0.68572531E-02,
     + -0.10100804E-01,-0.51137365E-01,-0.95617649E-03,-0.79342173E-02,
     + -0.47054091E-02,-0.45575444E-02, 0.71859444E-02, 0.13139450E-01,
     +  0.35423778E-01, 0.11673777E-01,-0.23047146E-02,-0.10572344E-01,
     + -0.13686701E-02, 0.70932191E-02, 0.20820221E-02, 0.44151731E-02,
     + -0.26249059E-02, 0.49566790E-01,-0.15768144E-01,-0.20867086E-02,
     + -0.88375346E-04,-0.26001786E-02,-0.87556959E-03,-0.66177542E-02,
     +  0.83008071E-03, 0.95116172E-03, 0.21043315E-02, 0.13544585E-01,
     +  0.33980706E-02, 0.35282698E-01, 0.30090557E-02, 0.41951751E-02,
     + -0.20007100E-02, 0.15397589E-01, 0.83888695E-03, 0.86639328E-02,
     +  0.11979281E-03,-0.22319127E-01,-0.12899696E-01, 0.54373634E-02,
     +  0.12716913E-01, 0.94310259E-02,-0.66154485E-03, 0.50943840E-03,
     +  0.48285765E-02, 0.16818603E-02, 0.41996482E-04,-0.12575183E-02,
     +  0.12148869E-02,-0.19636292E-01,-0.78882463E-02, 0.12420091E-02,
     + -0.10830882E-02,-0.25100701E-02,-0.70676985E-02, 0.74057246E-03,
     +  0.12152456E-02,-0.99501852E-02,-0.11854392E-01, 0.95682871E-02,
     + -0.11469017E-02,-0.81647076E-02,-0.43036086E-02,-0.57343883E-02,
     +  0.39660181E-02, 0.43026231E-01,-0.13477630E-02,-0.10389551E-01,
     +  0.96352780E-02, 0.76911017E-01, 0.34582242E-02, 0.24893822E-02,
     + -0.16172091E-01, 0.18211734E-02, 0.33005155E-02, 0.15458383E-01,
     +  0.40414450E-02, 0.85076308E-02,-0.59425429E-01,
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
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      y_sl_q3en   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x21            
     5  +coeff(  5)                *x51
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21*x31        
      y_sl_q3en   =y_sl_q3en   
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)        *x31    *x51
     2  +coeff( 11)            *x41*x51
     3  +coeff( 12)    *x22            
     4  +coeff( 13)    *x21        *x51
     5  +coeff( 14)                *x52
     6  +coeff( 15)            *x43    
     7  +coeff( 16)    *x21    *x42    
     8  +coeff( 17)    *x22*x31        
      y_sl_q3en   =y_sl_q3en   
     9  +coeff( 18)    *x22    *x41    
     1  +coeff( 19)    *x21*x31    *x51
     2  +coeff( 20)    *x21    *x41*x51
     3  +coeff( 21)            *x41*x52
     4  +coeff( 22)    *x23            
     5  +coeff( 23)    *x22        *x51
     6  +coeff( 24)    *x21    *x43    
     7  +coeff( 25)    *x22    *x42    
     8  +coeff( 26)    *x23    *x41    
      y_sl_q3en   =y_sl_q3en   
     9  +coeff( 27)    *x22*x31    *x51
     1  +coeff( 28)    *x22    *x41*x51
     2  +coeff( 29)        *x34*x41    
     3  +coeff( 30)    *x21*x33*x41    
     4  +coeff( 31)        *x31*x43*x51
     5  +coeff( 32)            *x44*x51
     6  +coeff( 33)    *x22*x33*x41    
     7  +coeff( 34)    *x22*x32*x42    
     8  +coeff( 35)    *x24*x31*x41    
      y_sl_q3en   =y_sl_q3en   
     9  +coeff( 36)    *x24*x34        
     1  +coeff( 37)        *x32        
     2  +coeff( 38)        *x31*x42    
     3  +coeff( 39)    *x21*x32        
     4  +coeff( 40)    *x21*x31*x41    
     5  +coeff( 41)    *x21        *x52
     6  +coeff( 42)                *x53
     7  +coeff( 43)    *x21*x32*x41    
     8  +coeff( 44)    *x21*x31*x42    
      y_sl_q3en   =y_sl_q3en   
     9  +coeff( 45)            *x43*x51
     1  +coeff( 46)    *x22*x31*x41    
     2  +coeff( 47)    *x21    *x42*x51
     3  +coeff( 48)    *x23*x31        
     4  +coeff( 49)    *x24            
     5  +coeff( 50)    *x22    *x42*x51
     6  +coeff( 51)*x11        *x41*x52
     7  +coeff( 52)    *x22*x32*x41*x51
     8  +coeff( 53)    *x22*x31*x42*x51
      y_sl_q3en   =y_sl_q3en   
     9  +coeff( 54)    *x22    *x45    
     1  +coeff( 55)    *x23    *x44    
     2  +coeff( 56)    *x23    *x43*x51
     3  +coeff( 57)    *x24*x31*x41*x51
     4  +coeff( 58)    *x22*x35*x41    
     5  +coeff( 59)    *x22*x35*x42    
     6  +coeff( 60)*x11*x21            
     7  +coeff( 61)    *x22*x32        
     8  +coeff( 62)    *x21*x32    *x51
      y_sl_q3en   =y_sl_q3en   
     9  +coeff( 63)            *x42*x52
     1  +coeff( 64)    *x23        *x51
     2  +coeff( 65)    *x22*x33        
     3  +coeff( 66)    *x22*x31*x42    
     4  +coeff( 67)    *x23*x31*x41    
     5  +coeff( 68)    *x24*x31        
     6  +coeff( 69)    *x22*x31    *x52
     7  +coeff( 70)            *x45*x51
     8  +coeff( 71)    *x21    *x44*x51
      y_sl_q3en   =y_sl_q3en   
     9  +coeff( 72)        *x33*x41*x52
     1  +coeff( 73)*x11*x21    *x43    
     2  +coeff( 74)    *x24    *x42    
     3  +coeff( 75)    *x22    *x42*x52
     4  +coeff( 76)    *x21*x34*x42    
     5  +coeff( 77)    *x21*x34*x41*x51
     6  +coeff( 78)    *x22*x32*x42*x51
     7  +coeff( 79)    *x23*x32*x41*x51
     8  +coeff( 80)*x11*x23    *x42    
      y_sl_q3en   =y_sl_q3en   
     9  +coeff( 81)    *x24*x32    *x51
     1  +coeff( 82)    *x22*x33*x43    
     2  +coeff( 83)*x11*x21*x34*x41    
     3  +coeff( 84)    *x23*x31*x44    
     4  +coeff( 85)    *x22*x31*x44*x51
     5  +coeff( 86)    *x24*x31*x43    
     6  +coeff( 87)*x11*x22    *x44    
     7  +coeff( 88)        *x34*x45    
     8  +coeff( 89)    *x22*x33*x41*x52
      y_sl_q3en   =y_sl_q3en   
     9  +coeff( 90)*x12*x21*x31    *x52
     1  +coeff( 91)*x11*x22*x34*x41    
     2  +coeff( 92)    *x24*x31*x44    
     3  +coeff( 93)*x12*x21*x33*x41    
     4  +coeff( 94)*x11*x23    *x42*x52
     5  +coeff( 95)    *x24*x31*x45    
c
      return
      end
      function p_sl_q3en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.1376353E-02/
      data xmin/
     1 -0.19937E-02,-0.52993E-01,-0.19935E-01,-0.29994E-01,-0.49883E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19971E-02, 0.54222E-01, 0.19987E-01, 0.29928E-01, 0.49858E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.88024570E-03,-0.98516997E-04,-0.11132711E-01,-0.13507579E-01,
     + -0.21021485E-02, 0.20355469E-05, 0.95088518E-03,-0.14517449E-02,
     + -0.13036741E-01,-0.62723491E-04, 0.17289488E-02, 0.16041868E-02,
     +  0.13164841E-01,-0.16930889E-02,-0.12220475E-02,-0.90644881E-02,
     + -0.15482890E-02, 0.17014102E-02,-0.16677394E-02,-0.20975452E-02,
     + -0.24845222E-04,-0.58919151E-03, 0.50670200E-03, 0.73868595E-02,
     + -0.36865121E-03, 0.33460750E-04, 0.10208774E-03, 0.12294942E-03,
     + -0.26162789E-03,-0.48137541E-03, 0.20853344E-02,-0.94097917E-03,
     +  0.18230703E-03,-0.24711704E-02,-0.24434358E-02, 0.25544346E-02,
     +  0.13925256E-03,-0.10567618E-03,-0.52654490E-04, 0.15643635E-03,
     + -0.54125205E-04,-0.96205244E-03, 0.12775540E-03, 0.33879708E-02,
     +  0.49229083E-03,-0.29791155E-03,-0.31597976E-03,-0.32344472E-02,
     + -0.34846287E-03, 0.20791858E-02,-0.64383593E-03, 0.16820363E-02,
     + -0.39296490E-03,-0.58914768E-03,-0.24967673E-02,-0.75300311E-03,
     +  0.75510115E-03,-0.12198380E-01,-0.27707096E-02, 0.22513231E-02,
     +  0.41094157E-02,-0.26720983E-02, 0.73416886E-04, 0.49237856E-05,
     +  0.21642577E-02, 0.26151759E-02, 0.17026655E-03,-0.47292162E-03,
     +  0.69647998E-03,-0.19568154E-02,-0.47436072E-02, 0.37874037E-03,
     +  0.62845298E-03, 0.31038860E-03, 0.96080272E-04, 0.74720464E-03,
     + -0.84973726E-04, 0.18570243E-03,-0.24650712E-02, 0.55808108E-03,
     +  0.31100318E-02,-0.97057266E-04,-0.10515341E-02, 0.26914233E-03,
     +  0.27928693E-03,-0.33488477E-03,-0.41011447E-03, 0.19636140E-03,
     +  0.55316638E-03, 0.35807971E-03,
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
      p_sl_q3en   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)*x11                
     7  +coeff(  7)    *x22            
     8  +coeff(  8)    *x21*x31        
      p_sl_q3en   =p_sl_q3en   
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)            *x42    
     2  +coeff( 11)    *x21        *x51
     3  +coeff( 12)        *x31    *x51
     4  +coeff( 13)            *x41*x51
     5  +coeff( 14)                *x52
     6  +coeff( 15)    *x22*x31        
     7  +coeff( 16)    *x22    *x41    
     8  +coeff( 17)            *x43    
      p_sl_q3en   =p_sl_q3en   
     9  +coeff( 18)    *x22        *x51
     1  +coeff( 19)    *x21    *x41*x51
     2  +coeff( 20)    *x22    *x41*x51
     3  +coeff( 21)    *x23    *x43    
     4  +coeff( 22)    *x23            
     5  +coeff( 23)            *x42*x51
     6  +coeff( 24)    *x23*x31*x42    
     7  +coeff( 25)        *x32*x41    
     8  +coeff( 26)    *x21    *x42    
      p_sl_q3en   =p_sl_q3en   
     9  +coeff( 27)        *x31*x41*x51
     1  +coeff( 28)    *x21        *x52
     2  +coeff( 29)            *x41*x52
     3  +coeff( 30)    *x24            
     4  +coeff( 31)    *x23    *x41    
     5  +coeff( 32)    *x22*x31    *x51
     6  +coeff( 33)    *x21*x34        
     7  +coeff( 34)    *x23*x31*x41    
     8  +coeff( 35)    *x23    *x42    
      p_sl_q3en   =p_sl_q3en   
     9  +coeff( 36)    *x22    *x42*x51
     1  +coeff( 37)        *x31*x41    
     2  +coeff( 38)    *x21*x32        
     3  +coeff( 39)        *x33        
     4  +coeff( 40)*x11        *x41    
     5  +coeff( 41)    *x21*x31*x41    
     6  +coeff( 42)        *x31*x42    
     7  +coeff( 43)                *x53
     8  +coeff( 44)    *x22    *x42    
      p_sl_q3en   =p_sl_q3en   
     9  +coeff( 45)        *x32*x42    
     1  +coeff( 46)    *x23        *x51
     2  +coeff( 47)    *x21    *x42*x51
     3  +coeff( 48)    *x21    *x44    
     4  +coeff( 49)    *x23*x31    *x51
     5  +coeff( 50)    *x22*x31*x41*x51
     6  +coeff( 51)    *x24*x31*x41    
     7  +coeff( 52)    *x23*x32*x41    
     8  +coeff( 53)        *x34*x42    
      p_sl_q3en   =p_sl_q3en   
     9  +coeff( 54)    *x22    *x43*x51
     1  +coeff( 55)    *x23*x32*x41*x51
     2  +coeff( 56)    *x23*x31*x42*x51
     3  +coeff( 57)*x11*x23    *x43    
     4  +coeff( 58)    *x23*x31*x44    
     5  +coeff( 59)    *x22    *x44*x52
     6  +coeff( 60)    *x24*x31*x44    
     7  +coeff( 61)    *x23*x33*x44    
     8  +coeff( 62)    *x24*x33*x41*x52
      p_sl_q3en   =p_sl_q3en   
     9  +coeff( 63)        *x32        
     1  +coeff( 64)        *x32    *x51
     2  +coeff( 65)    *x22*x31*x41    
     3  +coeff( 66)    *x21    *x43    
     4  +coeff( 67)    *x21*x32    *x51
     5  +coeff( 68)    *x23*x32        
     6  +coeff( 69)    *x22*x32*x41    
     7  +coeff( 70)    *x21*x32*x42    
     8  +coeff( 71)    *x21*x31*x43    
      p_sl_q3en   =p_sl_q3en   
     9  +coeff( 72)    *x22*x32    *x51
     1  +coeff( 73)    *x23    *x41*x51
     2  +coeff( 74)    *x21*x32*x41*x51
     3  +coeff( 75)*x11*x21        *x52
     4  +coeff( 76)    *x21    *x42*x52
     5  +coeff( 77)*x12*x21*x31        
     6  +coeff( 78)    *x24*x32        
     7  +coeff( 79)    *x24    *x42    
     8  +coeff( 80)    *x21*x32*x43    
      p_sl_q3en   =p_sl_q3en   
     9  +coeff( 81)    *x21*x31*x44    
     1  +coeff( 82)*x11*x23        *x51
     2  +coeff( 83)    *x24    *x41*x51
     3  +coeff( 84)    *x24*x33        
     4  +coeff( 85)*x12*x21*x31*x41    
     5  +coeff( 86)*x11*x23    *x42    
     6  +coeff( 87)*x11*x22    *x43    
     7  +coeff( 88)*x11*x23*x31    *x51
     8  +coeff( 89)    *x21*x33*x41*x52
      p_sl_q3en   =p_sl_q3en   
     9  +coeff( 90)*x12*x23*x31        
c
      return
      end
      function l_sl_q3en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff( 14)
      data ncoeff/ 13/
      data avdat/  0.1836096E+02/
      data xmin/
     1 -0.19937E-02,-0.52993E-01,-0.19935E-01,-0.29994E-01,-0.49883E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19971E-02, 0.54222E-01, 0.19987E-01, 0.29928E-01, 0.49858E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.49658613E-02, 0.31058845E+00, 0.32453444E-01,-0.43560178E-02,
     +  0.28199537E-01, 0.20619366E-01,-0.46545744E-01, 0.11977985E-01,
     + -0.26181718E-01,-0.22353781E-02, 0.62221107E-02,-0.59273312E-04,
     +  0.50448664E-02,
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
c
c                  function
c
      l_sl_q3en   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)    *x21    *x42    
     8  +coeff(  8)    *x21    *x41    
      l_sl_q3en   =l_sl_q3en   
     9  +coeff(  9)    *x21*x31*x41    
     1  +coeff( 10)            *x41    
     2  +coeff( 11)    *x21*x31        
     3  +coeff( 12)        *x31*x41    
     4  +coeff( 13)            *x42    
c
      return
      end
      function x_sl_q3ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.4573884E-03/
      data xmin/
     1 -0.19937E-02,-0.52993E-01,-0.19935E-01,-0.29861E-01,-0.45102E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19971E-02, 0.54222E-01, 0.19987E-01, 0.29928E-01, 0.49800E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.17187519E-01, 0.71122140E-01, 0.31320220E+00, 0.34826059E-01,
     + -0.20356311E-01,-0.42601973E-02,-0.20016376E-01,-0.32314558E-02,
     +  0.15264207E-01,-0.64669591E-02,-0.12531599E-01, 0.99504126E-04,
     +  0.32450301E-02, 0.29596023E-03, 0.53907157E-03,-0.77261335E-04,
     + -0.17571774E-02, 0.14811328E-02,-0.34751040E-02, 0.35078353E-02,
     + -0.47382887E-03, 0.23122285E-02, 0.25665624E-01,-0.14021257E-01,
     + -0.38978444E-04, 0.34427866E-01,-0.94526901E-03,-0.44621089E-02,
     + -0.82451748E-02, 0.23577485E-03, 0.57792489E-03,-0.38809326E-03,
     +  0.44250567E-03,-0.33679730E-03,-0.51000511E-03, 0.14680927E-01,
     + -0.21423597E-01,-0.22108087E-01,-0.29100361E-02, 0.59244743E-04,
     + -0.23891278E-03, 0.89211455E-04,-0.53393550E-03,-0.18036256E-01,
     + -0.12983531E-02, 0.10273508E-01,-0.24550126E-03,-0.10751795E-02,
     +  0.41289555E-03,-0.74394420E-03,-0.34152795E-01, 0.88979286E-04,
     + -0.17497397E-02,-0.10885214E-02, 0.10505573E-03,-0.11270125E-01,
     +  0.20567942E-02,-0.22139171E-01, 0.53517483E-02, 0.36393213E-02,
     +  0.36681136E-04, 0.20355075E-02,-0.13905944E-02, 0.16327439E-02,
     + -0.57995655E-02,-0.18274874E-02,-0.49492533E-04,-0.38621615E-03,
     +  0.14456632E-03,-0.54090860E-03,-0.20699589E-03,-0.25379716E-02,
     + -0.94035798E-03,-0.15999577E-02,-0.36096370E-02,-0.11566968E-03,
     +  0.64579430E-02,-0.16176143E-02, 0.43175966E-03,-0.29688425E-03,
     +  0.39080651E-02,-0.50112545E-02,-0.37990639E-02, 0.42864932E-02,
     +  0.14291207E-02,-0.12665847E-02,-0.17342571E-02, 0.35690835E-02,
     + -0.41537043E-02, 0.85519644E-03,
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
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      x_sl_q3ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x21        *x51
     5  +coeff(  5)                *x52
     6  +coeff(  6)    *x24            
     7  +coeff(  7)    *x21    *x42    
     8  +coeff(  8)*x11                
      x_sl_q3ex   =x_sl_q3ex   
     9  +coeff(  9)    *x22            
     1  +coeff( 10)            *x42    
     2  +coeff( 11)    *x21*x31*x41    
     3  +coeff( 12)    *x21    *x41    
     4  +coeff( 13)    *x21    *x41*x51
     5  +coeff( 14)    *x23*x31        
     6  +coeff( 15)    *x21*x31        
     7  +coeff( 16)        *x32        
     8  +coeff( 17)        *x31*x41    
      x_sl_q3ex   =x_sl_q3ex   
     9  +coeff( 18)            *x41*x51
     1  +coeff( 19)    *x21*x32        
     2  +coeff( 20)    *x22        *x51
     3  +coeff( 21)    *x21        *x52
     4  +coeff( 22)                *x53
     5  +coeff( 23)    *x21    *x43    
     6  +coeff( 24)    *x23*x31*x42    
     7  +coeff( 25)            *x41    
     8  +coeff( 26)    *x21*x31*x42    
      x_sl_q3ex   =x_sl_q3ex   
     9  +coeff( 27)    *x21    *x42*x51
     1  +coeff( 28)    *x23*x32*x41    
     2  +coeff( 29)    *x23    *x44    
     3  +coeff( 30)        *x31    *x51
     4  +coeff( 31)    *x23            
     5  +coeff( 32)    *x22    *x41    
     6  +coeff( 33)    *x21*x31    *x51
     7  +coeff( 34)        *x31*x41*x51
     8  +coeff( 35)    *x23    *x41    
      x_sl_q3ex   =x_sl_q3ex   
     9  +coeff( 36)    *x21*x32*x41    
     1  +coeff( 37)    *x21*x32*x42    
     2  +coeff( 38)    *x21    *x44    
     3  +coeff( 39)    *x24        *x51
     4  +coeff( 40)    *x21    *x43*x51
     5  +coeff( 41)    *x23*x33        
     6  +coeff( 42)    *x22*x31*x43    
     7  +coeff( 43)    *x24        *x52
     8  +coeff( 44)    *x23*x31*x43    
      x_sl_q3ex   =x_sl_q3ex   
     9  +coeff( 45)    *x22*x31*x41*x53
     1  +coeff( 46)    *x23*x34*x42    
     2  +coeff( 47)        *x32    *x51
     3  +coeff( 48)    *x23        *x51
     4  +coeff( 49)    *x22*x31    *x51
     5  +coeff( 50)    *x22        *x52
     6  +coeff( 51)    *x21*x31*x43    
     7  +coeff( 52)    *x23    *x41*x51
     8  +coeff( 53)    *x23        *x52
      x_sl_q3ex   =x_sl_q3ex   
     9  +coeff( 54)    *x22    *x41*x52
     1  +coeff( 55)    *x24*x32        
     2  +coeff( 56)    *x23    *x43    
     3  +coeff( 57)    *x22    *x44    
     4  +coeff( 58)    *x23*x32*x42    
     5  +coeff( 59)    *x23*x31*x42*x51
     6  +coeff( 60)    *x22*x31*x43*x54
     7  +coeff( 61)        *x33        
     8  +coeff( 62)    *x21*x33        
      x_sl_q3ex   =x_sl_q3ex   
     9  +coeff( 63)    *x22*x31*x41    
     1  +coeff( 64)    *x23*x32        
     2  +coeff( 65)    *x21*x33*x41    
     3  +coeff( 66)    *x23    *x42    
     4  +coeff( 67)*x11            *x51
     5  +coeff( 68)    *x22*x31*x41*x51
     6  +coeff( 69)    *x21    *x42*x52
     7  +coeff( 70)    *x22        *x53
     8  +coeff( 71)    *x23*x31*x41*x51
      x_sl_q3ex   =x_sl_q3ex   
     9  +coeff( 72)    *x23    *x42*x51
     1  +coeff( 73)    *x21*x31*x43*x51
     2  +coeff( 74)    *x23*x34        
     3  +coeff( 75)    *x23*x33*x41    
     4  +coeff( 76)*x11*x21    *x42    
     5  +coeff( 77)    *x23    *x43*x51
     6  +coeff( 78)    *x21*x31*x44*x51
     7  +coeff( 79)    *x23*x32    *x52
     8  +coeff( 80)*x11*x22*x32        
      x_sl_q3ex   =x_sl_q3ex   
     9  +coeff( 81)    *x24*x31*x43    
     1  +coeff( 82)    *x23*x31*x43*x51
     2  +coeff( 83)    *x21*x32*x44*x51
     3  +coeff( 84)    *x23*x31*x42*x52
     4  +coeff( 85)    *x21*x32*x41*x54
     5  +coeff( 86)*x11*x23*x31*x41    
     6  +coeff( 87)*x11*x24*x31*x41    
     7  +coeff( 88)*x11*x23*x31*x43    
     8  +coeff( 89)*x11*x24*x32*x42    
      x_sl_q3ex   =x_sl_q3ex   
     9  +coeff( 90)*x11*x22*x34    *x52
c
      return
      end
      function t_sl_q3ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff( 18)
      data ncoeff/ 17/
      data avdat/ -0.7832501E-03/
      data xmin/
     1 -0.19937E-02,-0.52993E-01,-0.19935E-01,-0.29861E-01,-0.45102E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19971E-02, 0.54222E-01, 0.19987E-01, 0.29928E-01, 0.49800E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.55237496E-02,-0.23087278E-01, 0.10616286E+00, 0.57496675E-02,
     + -0.96267648E-02,-0.44833543E-03, 0.19093701E-02,-0.16622456E-02,
     + -0.15359959E-02,-0.14360204E-03,-0.30103576E-03,-0.45724749E-03,
     +  0.37623435E-03, 0.67301758E-03,-0.53880765E-03, 0.59589336E-03,
     + -0.12024960E-02,
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
      x53 = x52*x5
c
c                  function
c
      t_sl_q3ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x21        *x51
     5  +coeff(  5)                *x52
     6  +coeff(  6)*x11                
     7  +coeff(  7)    *x22            
     8  +coeff(  8)            *x42    
      t_sl_q3ex   =t_sl_q3ex   
     9  +coeff(  9)    *x21    *x42    
     1  +coeff( 10)            *x41    
     2  +coeff( 11)    *x21    *x41    
     3  +coeff( 12)        *x31*x41    
     4  +coeff( 13)            *x41*x51
     5  +coeff( 14)    *x21    *x41*x51
     6  +coeff( 15)    *x21        *x52
     7  +coeff( 16)                *x53
     8  +coeff( 17)    *x24            
c
      return
      end
      function y_sl_q3ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff( 96)
      data ncoeff/ 95/
      data avdat/ -0.4302062E-02/
      data xmin/
     1 -0.19937E-02,-0.52993E-01,-0.19935E-01,-0.29861E-01,-0.45102E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19971E-02, 0.54222E-01, 0.19987E-01, 0.29928E-01, 0.49800E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.19699188E-02,-0.42706713E-01,-0.28502403E-01, 0.29911505E-03,
     + -0.11928904E-01,-0.52597621E-02,-0.51544696E-01, 0.56293737E-02,
     +  0.51978607E-01, 0.45375037E-02, 0.70325811E-02,-0.63244016E-02,
     + -0.11938747E-02,-0.61375299E-02,-0.57643456E-02,-0.39030273E-01,
     + -0.19280355E-03,-0.73262956E-02,-0.33872142E-02, 0.65992456E-02,
     +  0.81320237E-02, 0.24605433E-01, 0.86584445E-02,-0.13493192E-01,
     +  0.43313857E-03,-0.31616536E-02,-0.11826595E-02, 0.37071032E-02,
     +  0.57429308E-03, 0.10053574E-01,-0.87645007E-02,-0.35325647E-02,
     + -0.15897165E-02, 0.33905921E-02,-0.17859176E-02,-0.53957789E-02,
     +  0.18248064E-02, 0.10539616E-01, 0.36550676E-02, 0.20938700E-01,
     +  0.20934099E-02,-0.29231552E-02,-0.16536021E-02,-0.32002197E-02,
     + -0.65372442E-02,-0.23413390E-01,-0.45247607E-01,-0.50838469E-02,
     +  0.51954843E-03, 0.44640963E-03, 0.12182831E-02, 0.51643606E-03,
     + -0.14186039E-03, 0.17124731E-02, 0.45683904E-03,-0.14360036E-02,
     +  0.10003545E-03,-0.13930721E-01, 0.83468594E-02, 0.29199902E-03,
     + -0.17696471E-02,-0.50529307E-02, 0.45046099E-02, 0.16773112E-01,
     + -0.14610225E-03,-0.16817920E-02,-0.13491914E-02,-0.14117517E-01,
     + -0.22381062E-02, 0.33136162E-02, 0.15638877E-01,-0.43727368E-01,
     + -0.63949614E-02, 0.32040544E-02,-0.78492440E-01, 0.87083774E-02,
     + -0.72724201E-01, 0.11650392E-01,-0.42598750E-01, 0.40392812E-01,
     + -0.45544628E-01,-0.41808642E-03,-0.16834614E-02, 0.85541059E-03,
     +  0.22246267E-03, 0.36006582E-02,-0.31682616E-03,-0.77830418E-03,
     +  0.84423536E-03, 0.10500363E-02,-0.94862940E-03,-0.30482188E-02,
     +  0.14916657E-02,-0.35911257E-03,-0.22785892E-02,
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
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      y_sl_q3ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x21            
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x21*x31        
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)        *x31    *x51
      y_sl_q3ex   =y_sl_q3ex   
     9  +coeff(  9)            *x41*x51
     1  +coeff( 10)    *x22            
     2  +coeff( 11)    *x21        *x51
     3  +coeff( 12)                *x52
     4  +coeff( 13)        *x31*x42    
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)    *x22*x31        
     7  +coeff( 16)    *x22    *x41    
     8  +coeff( 17)    *x21*x31    *x51
      y_sl_q3ex   =y_sl_q3ex   
     9  +coeff( 18)    *x21    *x41*x51
     1  +coeff( 19)    *x23            
     2  +coeff( 20)    *x22        *x51
     3  +coeff( 21)    *x21    *x43    
     4  +coeff( 22)    *x22    *x42    
     5  +coeff( 23)    *x23    *x41    
     6  +coeff( 24)    *x22    *x41*x51
     7  +coeff( 25)        *x34*x41    
     8  +coeff( 26)        *x31*x44    
      y_sl_q3ex   =y_sl_q3ex   
     9  +coeff( 27)            *x45    
     1  +coeff( 28)    *x21*x33*x41    
     2  +coeff( 29)        *x32*x42*x51
     3  +coeff( 30)    *x22    *x42*x51
     4  +coeff( 31)    *x22*x33*x41    
     5  +coeff( 32)    *x23*x32*x41    
     6  +coeff( 33)    *x24*x31*x41    
     7  +coeff( 34)    *x23*x31*x41*x52
     8  +coeff( 35)    *x21*x32        
      y_sl_q3ex   =y_sl_q3ex   
     9  +coeff( 36)    *x21*x31*x41    
     1  +coeff( 37)            *x42*x51
     2  +coeff( 38)    *x21*x31*x42    
     3  +coeff( 39)    *x22*x32        
     4  +coeff( 40)    *x22*x31*x41    
     5  +coeff( 41)    *x23*x31        
     6  +coeff( 42)    *x22*x31    *x51
     7  +coeff( 43)    *x24            
     8  +coeff( 44)    *x22    *x43    
      y_sl_q3ex   =y_sl_q3ex   
     9  +coeff( 45)    *x23    *x44    
     1  +coeff( 46)    *x23*x33*x43    
     2  +coeff( 47)    *x24    *x45*x51
     3  +coeff( 48)            *x43    
     4  +coeff( 49)        *x31*x41*x51
     5  +coeff( 50)*x11        *x41    
     6  +coeff( 51)    *x21        *x52
     7  +coeff( 52)                *x53
     8  +coeff( 53)    *x21*x33        
      y_sl_q3ex   =y_sl_q3ex   
     9  +coeff( 54)    *x21    *x42*x51
     1  +coeff( 55)*x11*x21    *x41    
     2  +coeff( 56)    *x23        *x51
     3  +coeff( 57)    *x21    *x44    
     4  +coeff( 58)    *x22*x31*x42    
     5  +coeff( 59)    *x22*x31*x41*x51
     6  +coeff( 60)    *x21    *x45    
     7  +coeff( 61)    *x22*x34        
     8  +coeff( 62)    *x21    *x44*x51
      y_sl_q3ex   =y_sl_q3ex   
     9  +coeff( 63)    *x22*x32*x41*x51
     1  +coeff( 64)    *x22    *x43*x51
     2  +coeff( 65)    *x24    *x42    
     3  +coeff( 66)    *x22    *x42*x52
     4  +coeff( 67)*x11*x22*x31    *x51
     5  +coeff( 68)    *x22    *x45    
     6  +coeff( 69)        *x34*x43*x51
     7  +coeff( 70)    *x23*x31*x44    
     8  +coeff( 71)    *x23    *x45    
      y_sl_q3ex   =y_sl_q3ex   
     9  +coeff( 72)    *x24    *x44    
     1  +coeff( 73)*x11*x23    *x44    
     2  +coeff( 74)*x11*x21    *x44*x52
     3  +coeff( 75)    *x24*x31*x45    
     4  +coeff( 76)    *x22*x35*x45    
     5  +coeff( 77)    *x24*x34*x44    
     6  +coeff( 78)*x11*x23*x35*x42    
     7  +coeff( 79)*x11*x24*x32*x45*x51
     8  +coeff( 80)*x12*x24*x32*x44    
      y_sl_q3ex   =y_sl_q3ex   
     9  +coeff( 81)*x12*x24    *x44*x52
     1  +coeff( 82)            *x42    
     2  +coeff( 83)        *x32*x41    
     3  +coeff( 84)        *x31    *x52
     4  +coeff( 85)*x11*x21            
     5  +coeff( 86)    *x21*x32*x41    
     6  +coeff( 87)*x11*x21*x31        
     7  +coeff( 88)            *x41*x53
     8  +coeff( 89)    *x21*x34        
      y_sl_q3ex   =y_sl_q3ex   
     9  +coeff( 90)    *x22*x33        
     1  +coeff( 91)        *x33    *x52
     2  +coeff( 92)    *x23*x31*x41    
     3  +coeff( 93)    *x22*x32    *x51
     4  +coeff( 94)*x11        *x42*x51
     5  +coeff( 95)    *x21*x31*x41*x52
c
      return
      end
      function p_sl_q3ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.6047979E-03/
      data xmin/
     1 -0.19937E-02,-0.52993E-01,-0.19935E-01,-0.29861E-01,-0.45102E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19971E-02, 0.54222E-01, 0.19987E-01, 0.29928E-01, 0.49800E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.15947319E-03,-0.10618750E-03, 0.14535603E-01,-0.43100161E-02,
     +  0.65222313E-02,-0.28353043E-04,-0.23248792E-02, 0.16072053E-02,
     +  0.19189900E-01, 0.30503495E-03,-0.26180034E-02,-0.44529364E-02,
     + -0.21508783E-01, 0.18659276E-02, 0.35504312E-02, 0.15155208E-01,
     +  0.22035756E-02,-0.21900027E-02, 0.34199439E-02,-0.15935164E-01,
     + -0.38399818E-03, 0.15696757E-02, 0.13074977E-02, 0.15834337E-02,
     +  0.36361525E-03, 0.17357923E-02, 0.11596950E-02,-0.12471813E-02,
     + -0.19565246E-02,-0.26077097E-02,-0.11691768E-01,-0.49521765E-02,
     +  0.18841537E-02, 0.13928404E-01, 0.46041835E-03, 0.69085720E-04,
     +  0.48193726E-03, 0.99782064E-03,-0.12403418E-02,-0.11473270E-02,
     + -0.57153604E-02,-0.40360712E-03, 0.91529330E-02, 0.11234723E-02,
     +  0.16436737E-01,-0.12367778E-02,-0.10129201E-02, 0.18422180E-03,
     + -0.12764711E-01,-0.18188817E-03, 0.43419972E-04,-0.16985650E-03,
     +  0.40102410E-03,-0.52616140E-03,-0.34349150E-03, 0.50002499E-03,
     + -0.88480487E-03,-0.82085817E-03, 0.38006695E-02,-0.93698257E-03,
     +  0.14718180E-01, 0.96578483E-03,-0.17975413E-03, 0.23671603E-02,
     + -0.64094062E-03,-0.17932574E-02, 0.37625392E-02, 0.82838727E-03,
     + -0.26793503E-02, 0.47381306E-02, 0.18013333E-02,-0.47020160E-03,
     +  0.24864709E-02,-0.10604928E-02,-0.59173456E-02, 0.52990941E-02,
     +  0.75972285E-02,-0.30916932E-03,-0.44313069E-02,-0.18612635E-02,
     +  0.50638635E-02, 0.40930882E-02,-0.27359447E-02,-0.26656289E-02,
     + -0.28200990E-02,-0.19781019E-02,-0.60149911E-02,-0.11762178E-02,
     + -0.29673318E-02,-0.60923742E-02,
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
      p_sl_q3ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)*x11                
     7  +coeff(  7)    *x22            
     8  +coeff(  8)    *x21*x31        
      p_sl_q3ex   =p_sl_q3ex   
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)            *x42    
     2  +coeff( 11)    *x21        *x51
     3  +coeff( 12)        *x31    *x51
     4  +coeff( 13)            *x41*x51
     5  +coeff( 14)                *x52
     6  +coeff( 15)    *x22*x31        
     7  +coeff( 16)    *x22    *x41    
     8  +coeff( 17)            *x43    
      p_sl_q3ex   =p_sl_q3ex   
     9  +coeff( 18)    *x22        *x51
     1  +coeff( 19)            *x41*x52
     2  +coeff( 20)    *x22    *x42    
     3  +coeff( 21)    *x23    *x42    
     4  +coeff( 22)    *x24*x31*x41    
     5  +coeff( 23)    *x23            
     6  +coeff( 24)    *x21*x31*x41    
     7  +coeff( 25)        *x32*x41    
     8  +coeff( 26)    *x21    *x42    
      p_sl_q3ex   =p_sl_q3ex   
     9  +coeff( 27)    *x21    *x41*x51
     1  +coeff( 28)            *x42*x51
     2  +coeff( 29)    *x22*x32        
     3  +coeff( 30)    *x23    *x41    
     4  +coeff( 31)    *x22*x31*x41    
     5  +coeff( 32)    *x21    *x43    
     6  +coeff( 33)    *x22    *x41*x51
     7  +coeff( 34)    *x23    *x44    
     8  +coeff( 35)    *x21*x32        
      p_sl_q3ex   =p_sl_q3ex   
     9  +coeff( 36)        *x33        
     1  +coeff( 37)        *x31*x42    
     2  +coeff( 38)        *x31    *x52
     3  +coeff( 39)    *x23*x31        
     4  +coeff( 40)    *x21*x32*x41    
     5  +coeff( 41)    *x21*x31*x42    
     6  +coeff( 42)    *x21*x31    *x52
     7  +coeff( 43)    *x22    *x43    
     8  +coeff( 44)*x11*x22    *x42    
      p_sl_q3ex   =p_sl_q3ex   
     9  +coeff( 45)    *x23*x31*x43    
     1  +coeff( 46)    *x24    *x42*x51
     2  +coeff( 47)    *x22    *x43*x52
     3  +coeff( 48)    *x22*x32*x41*x53
     4  +coeff( 49)*x11*x24    *x44    
     5  +coeff( 50)*x11*x21            
     6  +coeff( 51)*x11    *x31        
     7  +coeff( 52)*x11        *x41    
     8  +coeff( 53)    *x21*x31    *x51
      p_sl_q3ex   =p_sl_q3ex   
     9  +coeff( 54)        *x31*x41*x51
     1  +coeff( 55)                *x53
     2  +coeff( 56)    *x24            
     3  +coeff( 57)            *x41*x53
     4  +coeff( 58)    *x24*x31        
     5  +coeff( 59)    *x22*x32*x41    
     6  +coeff( 60)    *x21*x33*x41    
     7  +coeff( 61)    *x22*x31*x42    
     8  +coeff( 62)        *x31*x44    
      p_sl_q3ex   =p_sl_q3ex   
     9  +coeff( 63)    *x21*x33    *x51
     1  +coeff( 64)    *x21*x31*x41*x52
     2  +coeff( 65)*x11*x23    *x41    
     3  +coeff( 66)*x11*x22*x31*x41    
     4  +coeff( 67)    *x24    *x42    
     5  +coeff( 68)*x11*x22*x31    *x51
     6  +coeff( 69)    *x23    *x41*x52
     7  +coeff( 70)    *x22    *x42*x52
     8  +coeff( 71)    *x21    *x43*x52
      p_sl_q3ex   =p_sl_q3ex   
     9  +coeff( 72)        *x33    *x53
     1  +coeff( 73)*x11*x23    *x42    
     2  +coeff( 74)*x11*x21*x32*x42    
     3  +coeff( 75)    *x22*x31*x44    
     4  +coeff( 76)    *x24    *x44    
     5  +coeff( 77)    *x23*x31*x44    
     6  +coeff( 78)*x12*x23        *x51
     7  +coeff( 79)    *x24*x32*x41*x51
     8  +coeff( 80)*x11*x22*x31*x42*x51
      p_sl_q3ex   =p_sl_q3ex   
     9  +coeff( 81)    *x22*x33*x41*x52
     1  +coeff( 82)*x11*x24    *x43    
     2  +coeff( 83)*x11*x22*x32*x43    
     3  +coeff( 84)*x11*x24*x32    *x51
     4  +coeff( 85)*x11*x24*x31*x41*x51
     5  +coeff( 86)    *x22*x34*x42*x51
     6  +coeff( 87)    *x23*x33*x41*x52
     7  +coeff( 88)*x11*x24*x34        
     8  +coeff( 89)*x11*x23*x33*x42    
      p_sl_q3ex   =p_sl_q3ex   
     9  +coeff( 90)*x11*x23*x32*x43*x51
c
      return
      end
      function l_sl_q3ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff( 16)
      data ncoeff/ 15/
      data avdat/  0.2134420E+02/
      data xmin/
     1 -0.19937E-02,-0.52993E-01,-0.19935E-01,-0.29861E-01,-0.45102E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19971E-02, 0.54222E-01, 0.19987E-01, 0.29928E-01, 0.49800E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.11713376E-03, 0.31109396E+00, 0.30797882E-01,-0.44379979E-02,
     +  0.32031018E-01,-0.56516781E-01, 0.98580765E-02, 0.12771446E-01,
     +  0.73767561E-02,-0.43181904E-01,-0.25861546E-01,-0.23767415E-02,
     +  0.53913337E-02,-0.93278977E-04, 0.49752607E-02,
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
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_sl_q3ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x23    *x42    
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)    *x21        *x51
      l_sl_q3ex   =l_sl_q3ex   
     9  +coeff(  9)                *x52
     1  +coeff( 10)    *x23*x31*x41    
     2  +coeff( 11)    *x21    *x44    
     3  +coeff( 12)            *x41    
     4  +coeff( 13)    *x21*x31        
     5  +coeff( 14)        *x31*x41    
     6  +coeff( 15)            *x42    
c
      return
      end
      function x_sl_fp     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.9311447E-02/
      data xmin/
     1 -0.19999E-02,-0.53425E-01,-0.19919E-01,-0.29893E-01,-0.45450E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19995E-02, 0.54621E-01, 0.19981E-01, 0.29903E-01, 0.49857E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.19617096E-01, 0.44036130E-02, 0.65025371E-04,-0.21782123E-04,
     +  0.62164742E+00, 0.52960210E-01,-0.48269253E-01, 0.21597881E-01,
     + -0.12059635E-01,-0.26911739E-01,-0.45886110E-02,-0.15536460E-01,
     + -0.25379837E-02,-0.39360810E-04, 0.30247956E-02, 0.30453780E-03,
     +  0.43430477E-02, 0.50395988E-02,-0.92357900E-02, 0.19862309E-01,
     +  0.38241358E-02, 0.81116130E-03, 0.93677320E-03,-0.33897688E-02,
     + -0.48030210E-02,-0.25249296E-02,-0.71686081E-03, 0.27738180E-01,
     + -0.32674498E-02, 0.77593041E-03,-0.14215168E-02, 0.47284085E-02,
     + -0.10363095E-02, 0.16360430E-02,-0.25346458E-02, 0.10016531E-01,
     + -0.30742448E-02,-0.29692294E-02,-0.58747181E-02,-0.61151246E-03,
     + -0.56112988E-03,-0.65093650E-03, 0.22323008E-03,-0.11762679E-02,
     +  0.58663450E-03, 0.55470417E-04,-0.16575880E-01, 0.65153912E-02,
     + -0.32037366E-02, 0.56453398E-02,-0.16267559E-01,-0.19025248E-02,
     + -0.13098817E-01,-0.11150754E-01, 0.36985957E-03, 0.75623044E-02,
     +  0.20171304E-02,-0.11219239E-01, 0.13756095E-02, 0.56625973E-03,
     +  0.98306022E-03, 0.39390445E-04, 0.36345053E-03, 0.15312186E-02,
     + -0.13602054E-02,-0.41474295E-02,-0.13500242E-01,-0.14733316E-03,
     + -0.17281321E-02, 0.21598036E-02, 0.62162299E-02,-0.79320197E-03,
     +  0.47593503E-02,-0.53998702E-02,-0.77792434E-02,-0.27723117E-01,
     + -0.14706133E-01, 0.25274415E-03, 0.71646553E-03,-0.24134875E-02,
     +  0.19306449E-02, 0.32812444E-03,-0.21886484E-03,-0.53353999E-02,
     +  0.15342095E-02, 0.20947419E-02,-0.10432374E-03,-0.22440458E-03,
     +  0.35842819E-03,-0.27862829E-03,
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
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      x_sl_fp     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)                *x52
     8  +coeff(  8)    *x22            
      x_sl_fp     =x_sl_fp     
     9  +coeff(  9)            *x42    
     1  +coeff( 10)    *x21    *x42    
     2  +coeff( 11)*x11                
     3  +coeff( 12)    *x21*x31*x41    
     4  +coeff( 13)    *x21        *x52
     5  +coeff( 14)        *x32        
     6  +coeff( 15)            *x41*x51
     7  +coeff( 16)    *x21*x31    *x51
     8  +coeff( 17)    *x21    *x41*x51
      x_sl_fp     =x_sl_fp     
     9  +coeff( 18)                *x53
     1  +coeff( 19)    *x24            
     2  +coeff( 20)    *x21    *x43    
     3  +coeff( 21)    *x23*x32*x41    
     4  +coeff( 22)    *x21*x31        
     5  +coeff( 23)    *x21    *x41    
     6  +coeff( 24)        *x31*x41    
     7  +coeff( 25)    *x21*x32        
     8  +coeff( 26)    *x22*x31*x41    
      x_sl_fp     =x_sl_fp     
     9  +coeff( 27)    *x22    *x42    
     1  +coeff( 28)    *x21*x31*x42    
     2  +coeff( 29)    *x21    *x42*x51
     3  +coeff( 30)        *x31    *x51
     4  +coeff( 31)    *x22    *x41    
     5  +coeff( 32)    *x22        *x51
     6  +coeff( 33)        *x31*x41*x51
     7  +coeff( 34)    *x21*x33        
     8  +coeff( 35)    *x23    *x41    
      x_sl_fp     =x_sl_fp     
     9  +coeff( 36)    *x21*x32*x41    
     1  +coeff( 37)    *x23        *x51
     2  +coeff( 38)    *x22        *x52
     3  +coeff( 39)    *x24        *x51
     4  +coeff( 40)    *x22*x31        
     5  +coeff( 41)        *x32    *x51
     6  +coeff( 42)            *x41*x52
     7  +coeff( 43)    *x22*x31    *x51
     8  +coeff( 44)    *x21*x31*x41*x51
      x_sl_fp     =x_sl_fp     
     9  +coeff( 45)    *x21    *x41*x52
     1  +coeff( 46)*x11*x21            
     2  +coeff( 47)    *x21    *x44    
     3  +coeff( 48)    *x21    *x43*x51
     4  +coeff( 49)    *x23        *x52
     5  +coeff( 50)    *x24    *x42    
     6  +coeff( 51)    *x21*x33*x42    
     7  +coeff( 52)    *x24    *x41*x51
     8  +coeff( 53)    *x21*x34*x42    
      x_sl_fp     =x_sl_fp     
     9  +coeff( 54)    *x21*x33*x43    
     1  +coeff( 55)    *x23*x31*x42*x51
     2  +coeff( 56)    *x24    *x44    
     3  +coeff( 57)    *x22*x32*x44    
     4  +coeff( 58)    *x24*x32*x44    
     5  +coeff( 59)    *x23            
     6  +coeff( 60)    *x22*x32        
     7  +coeff( 61)            *x42*x52
     8  +coeff( 62)*x11    *x31        
      x_sl_fp     =x_sl_fp     
     9  +coeff( 63)    *x22*x33        
     1  +coeff( 64)    *x21*x34        
     2  +coeff( 65)    *x23*x31*x41    
     3  +coeff( 66)    *x23    *x42    
     4  +coeff( 67)    *x21*x31*x43    
     5  +coeff( 68)*x11            *x51
     6  +coeff( 69)    *x22*x31*x41*x51
     7  +coeff( 70)    *x21*x32*x41*x51
     8  +coeff( 71)    *x21*x31*x42*x51
      x_sl_fp     =x_sl_fp     
     9  +coeff( 72)    *x22    *x41*x52
     1  +coeff( 73)    *x24*x31*x41    
     2  +coeff( 74)    *x23*x31*x42    
     3  +coeff( 75)    *x23    *x43    
     4  +coeff( 76)    *x21*x32*x43    
     5  +coeff( 77)    *x21*x31*x44    
     6  +coeff( 78)        *x33*x41*x52
     7  +coeff( 79)    *x22    *x42*x52
     8  +coeff( 80)    *x21    *x42*x53
      x_sl_fp     =x_sl_fp     
     9  +coeff( 81)        *x31*x41*x54
     1  +coeff( 82)*x11*x23            
     2  +coeff( 83)        *x33    *x54
     3  +coeff( 84)    *x23*x34*x41    
     4  +coeff( 85)*x11*x22    *x42*x51
     5  +coeff( 86)*x11*x23    *x42*x52
     6  +coeff( 87)        *x32*x41    
     7  +coeff( 88)        *x31    *x52
     8  +coeff( 89)        *x32    *x52
      x_sl_fp     =x_sl_fp     
     9  +coeff( 90)        *x31    *x53
c
      return
      end
      function t_sl_fp     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff( 19)
      data ncoeff/ 18/
      data avdat/  0.1176425E-02/
      data xmin/
     1 -0.19999E-02,-0.53425E-01,-0.19919E-01,-0.29893E-01,-0.45450E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19995E-02, 0.54621E-01, 0.19981E-01, 0.29903E-01, 0.49857E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.35653585E-02,-0.23732562E-01,-0.49484246E-04,-0.15055630E-03,
     +  0.10668138E+00, 0.57430984E-02,-0.97457115E-02, 0.15948142E-02,
     + -0.45146610E-03,-0.36692666E-03,-0.17641935E-02,-0.80999732E-03,
     +  0.35067947E-03, 0.55016839E-03, 0.79383462E-03, 0.10262458E-02,
     + -0.10804459E-02, 0.13411328E-02,
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
      x53 = x52*x5
c
c                  function
c
      t_sl_fp     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)                *x52
     8  +coeff(  8)    *x22            
      t_sl_fp     =t_sl_fp     
     9  +coeff(  9)*x11                
     1  +coeff( 10)        *x31*x41    
     2  +coeff( 11)            *x42    
     3  +coeff( 12)    *x21        *x52
     4  +coeff( 13)            *x41*x51
     5  +coeff( 14)    *x23            
     6  +coeff( 15)    *x21*x31*x41    
     7  +coeff( 16)                *x53
     8  +coeff( 17)    *x24            
      t_sl_fp     =t_sl_fp     
     9  +coeff( 18)    *x22    *x42    
c
      return
      end
      function y_sl_fp     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff( 96)
      data ncoeff/ 95/
      data avdat/ -0.2163797E-02/
      data xmin/
     1 -0.19999E-02,-0.53425E-01,-0.19919E-01,-0.29893E-01,-0.45450E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19995E-02, 0.54621E-01, 0.19981E-01, 0.29903E-01, 0.49857E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.24499488E-02,-0.67417527E-03,-0.40862828E-01,-0.26038001E-03,
     +  0.68441406E-02,-0.56901185E-05, 0.44037271E-03, 0.40989142E-03,
     + -0.15460377E-02, 0.34127741E-02,-0.71162325E-02,-0.92874002E-02,
     + -0.21619799E-02,-0.11541613E-02,-0.14155227E-02,-0.18525755E-03,
     + -0.20321783E-03, 0.40916173E-03,-0.22998959E-03, 0.63007199E-02,
     +  0.92577818E-03,-0.45252172E-02, 0.38582862E-02, 0.10344232E-01,
     + -0.12953043E-01,-0.42131068E-02, 0.11542011E-02,-0.58996207E-02,
     +  0.68230421E-03, 0.32767807E-02,-0.27292208E-02,-0.19335413E-01,
     + -0.59627718E-02,-0.23369775E-02, 0.70541918E-01, 0.36697529E-01,
     +  0.31723760E-02, 0.31858278E-02, 0.11804031E-02,-0.30742155E-02,
     + -0.61828416E-03,-0.57602623E-02,-0.78013720E-03, 0.10381354E-02,
     + -0.83626084E-01,-0.72286092E-03, 0.38861725E-01, 0.31332497E-01,
     +  0.38315977E-02,-0.41860629E-01, 0.76596034E-02,-0.23939523E-04,
     +  0.52880403E-03, 0.42386801E-03,-0.36794879E-03, 0.16157159E-03,
     + -0.17591766E-02,-0.45334469E-03,-0.10852811E-02,-0.41427364E-03,
     +  0.26591262E-02,-0.83585699E-04, 0.10851456E-03, 0.24806249E-02,
     +  0.18011359E-03,-0.48202105E-01,-0.22418636E-02, 0.95112100E-02,
     +  0.17814107E-02, 0.44361725E-02,-0.18164782E-01, 0.10423620E+00,
     +  0.11755866E-02, 0.91939410E-02,-0.21131053E-02,-0.75559825E-01,
     + -0.28065720E+00,-0.25453019E+00, 0.87960428E-02,-0.16740121E-01,
     +  0.17864033E-03, 0.26048612E-03,-0.38720077E-03,-0.13232004E-02,
     + -0.25380275E-03, 0.44061302E-03,-0.66626829E-03, 0.19217774E-03,
     + -0.68491756E-03, 0.43930244E-03, 0.39889413E-03, 0.11754503E-01,
     + -0.39850478E-03, 0.13254932E-02,-0.56909741E-03,
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
      x35 = x34*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x45 = x44*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      y_sl_fp     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x21            
     5  +coeff(  5)                *x51
     6  +coeff(  6)        *x32        
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)            *x42    
      y_sl_fp     =y_sl_fp     
     9  +coeff(  9)    *x21*x31        
     1  +coeff( 10)    *x21    *x41    
     2  +coeff( 11)        *x31    *x51
     3  +coeff( 12)            *x41*x51
     4  +coeff( 13)    *x22            
     5  +coeff( 14)    *x21        *x51
     6  +coeff( 15)                *x52
     7  +coeff( 16)        *x32*x41    
     8  +coeff( 17)        *x31*x42    
      y_sl_fp     =y_sl_fp     
     9  +coeff( 18)            *x43    
     1  +coeff( 19)            *x42*x51
     2  +coeff( 20)    *x22    *x41    
     3  +coeff( 21)    *x21*x31    *x51
     4  +coeff( 22)    *x21    *x41*x51
     5  +coeff( 23)        *x31    *x52
     6  +coeff( 24)            *x41*x52
     7  +coeff( 25)    *x22    *x42    
     8  +coeff( 26)            *x41*x53
      y_sl_fp     =y_sl_fp     
     9  +coeff( 27)    *x22*x33        
     1  +coeff( 28)    *x22    *x43    
     2  +coeff( 29)    *x24*x31        
     3  +coeff( 30)    *x24    *x41    
     4  +coeff( 31)    *x22*x34        
     5  +coeff( 32)    *x22*x33*x41    
     6  +coeff( 33)    *x22    *x44    
     7  +coeff( 34)    *x23    *x43    
     8  +coeff( 35)    *x22*x32*x43    
      y_sl_fp     =y_sl_fp     
     9  +coeff( 36)    *x22    *x45    
     1  +coeff( 37)    *x22*x35*x41    
     2  +coeff( 38)    *x22*x31        
     3  +coeff( 39)    *x21        *x52
     4  +coeff( 40)    *x22*x31*x41    
     5  +coeff( 41)    *x22*x31    *x51
     6  +coeff( 42)    *x22    *x41*x51
     7  +coeff( 43)    *x23        *x51
     8  +coeff( 44)    *x23    *x41*x51
      y_sl_fp     =y_sl_fp     
     9  +coeff( 45)    *x22*x31*x43    
     1  +coeff( 46)    *x23*x31*x42    
     2  +coeff( 47)    *x22*x33*x42    
     3  +coeff( 48)    *x22*x31*x44    
     4  +coeff( 49)    *x24    *x43*x51
     5  +coeff( 50)    *x22*x32*x45    
     6  +coeff( 51)*x11*x23    *x44*x51
     7  +coeff( 52)*x11*x21            
     8  +coeff( 53)    *x23            
      y_sl_fp     =y_sl_fp     
     9  +coeff( 54)    *x22        *x51
     1  +coeff( 55)                *x53
     2  +coeff( 56)        *x34        
     3  +coeff( 57)            *x43*x51
     4  +coeff( 58)    *x22*x32        
     5  +coeff( 59)        *x31    *x53
     6  +coeff( 60)    *x24            
     7  +coeff( 61)    *x22*x32*x41    
     8  +coeff( 62)    *x21    *x43*x51
      y_sl_fp     =y_sl_fp     
     9  +coeff( 63)    *x23*x32        
     1  +coeff( 64)    *x23    *x42    
     2  +coeff( 65)    *x21*x31*x41*x52
     3  +coeff( 66)    *x22*x32*x42    
     4  +coeff( 67)    *x23*x31*x41*x51
     5  +coeff( 68)    *x22*x34*x41    
     6  +coeff( 69)    *x23*x33*x41    
     7  +coeff( 70)    *x23    *x43*x51
     8  +coeff( 71)    *x22*x33*x43    
      y_sl_fp     =y_sl_fp     
     9  +coeff( 72)    *x22*x31*x45    
     1  +coeff( 73)*x11*x22*x32*x42    
     2  +coeff( 74)    *x22*x34*x43    
     3  +coeff( 75)    *x23    *x42*x53
     4  +coeff( 76)    *x22*x35*x43    
     5  +coeff( 77)    *x22*x34*x44    
     6  +coeff( 78)    *x22*x33*x45    
     7  +coeff( 79)*x11*x23*x32*x44    
     8  +coeff( 80)*x11*x24*x31*x43*x54
      y_sl_fp     =y_sl_fp     
     9  +coeff( 81)    *x21*x32        
     1  +coeff( 82)    *x21*x33        
     2  +coeff( 83)        *x32*x41*x51
     3  +coeff( 84)        *x31*x42*x51
     4  +coeff( 85)    *x21*x32    *x51
     5  +coeff( 86)            *x42*x52
     6  +coeff( 87)    *x23*x31        
     7  +coeff( 88)*x11*x21    *x41    
     8  +coeff( 89)    *x23    *x41    
      y_sl_fp     =y_sl_fp     
     9  +coeff( 90)    *x22        *x52
     1  +coeff( 91)                *x54
     2  +coeff( 92)    *x22*x31*x42    
     3  +coeff( 93)*x11*x21*x32        
     4  +coeff( 94)    *x23*x31*x41    
     5  +coeff( 95)*x11*x21    *x42    
c
      return
      end
      function p_sl_fp     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.1420645E-03/
      data xmin/
     1 -0.19999E-02,-0.53425E-01,-0.19919E-01,-0.29893E-01,-0.45450E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19995E-02, 0.54621E-01, 0.19981E-01, 0.29903E-01, 0.49857E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.42746065E-03, 0.14436061E-01,-0.43622400E-02, 0.66311406E-02,
     + -0.23491955E-02, 0.13496008E-02, 0.18643621E-01, 0.11576367E-03,
     +  0.17371688E-03,-0.28206159E-02,-0.46929633E-02,-0.21671342E-01,
     +  0.18418469E-02, 0.99592458E-03, 0.35492829E-02, 0.14778393E-01,
     +  0.89530123E-03, 0.27481727E-02,-0.26521687E-02, 0.37731847E-02,
     + -0.15916042E-02,-0.47577792E-02,-0.16883092E-02,-0.24417685E-02,
     + -0.30501508E-02,-0.11123642E-02, 0.81659819E-04,-0.88634370E-02,
     +  0.58208029E-02, 0.10319564E-01, 0.33169211E-03, 0.99885568E-03,
     +  0.11889620E-02,-0.15799167E-03,-0.20122298E-02, 0.11109140E-02,
     + -0.38781853E-02,-0.72795311E-02,-0.13680085E-01,-0.22378962E-01,
     +  0.22260232E-01, 0.21150259E-03,-0.25702259E-03, 0.27511395E-02,
     + -0.72898914E-03,-0.46085063E-03,-0.52222018E-02,-0.57821674E-03,
     +  0.13005678E-01, 0.64525469E-02,-0.53984267E-02,-0.12245680E-01,
     + -0.38429755E-02,-0.15807353E-01,-0.24803905E-02, 0.80515370E-02,
     + -0.27655743E-02, 0.45456677E-02,-0.66122161E-02, 0.20159250E-01,
     + -0.19676287E-01,-0.67104811E-04, 0.14589251E-05, 0.16994971E-02,
     + -0.97598287E-03,-0.18499278E-03, 0.12672255E-02,-0.62888017E-03,
     +  0.15849840E-03, 0.83943475E-02,-0.12415784E-02, 0.42482279E-02,
     +  0.44039008E-02, 0.12097076E-02, 0.32940554E-03,-0.16444612E-02,
     + -0.29503493E-01,-0.23337651E-01, 0.11381377E-02, 0.39657145E-02,
     +  0.30353481E-02,-0.77381381E-03, 0.40768128E-03, 0.13455774E-02,
     +  0.64685321E-02, 0.96348412E-02, 0.39127711E-02,-0.74378881E-02,
     +  0.19872885E-02, 0.46573032E-03,
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
      p_sl_fp     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21*x31        
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)        *x31*x41    
      p_sl_fp     =p_sl_fp     
     9  +coeff(  9)            *x42    
     1  +coeff( 10)    *x21        *x51
     2  +coeff( 11)        *x31    *x51
     3  +coeff( 12)            *x41*x51
     4  +coeff( 13)                *x52
     5  +coeff( 14)    *x23            
     6  +coeff( 15)    *x22*x31        
     7  +coeff( 16)    *x22    *x41    
     8  +coeff( 17)    *x21    *x42    
      p_sl_fp     =p_sl_fp     
     9  +coeff( 18)            *x43    
     1  +coeff( 19)    *x22        *x51
     2  +coeff( 20)            *x41*x52
     3  +coeff( 21)    *x23    *x41    
     4  +coeff( 22)    *x22    *x42    
     5  +coeff( 23)    *x21*x33*x41    
     6  +coeff( 24)        *x33*x42    
     7  +coeff( 25)    *x22    *x42*x51
     8  +coeff( 26)    *x24*x31*x41    
      p_sl_fp     =p_sl_fp     
     9  +coeff( 27)    *x24    *x42    
     1  +coeff( 28)    *x23    *x43    
     2  +coeff( 29)    *x24    *x43    
     3  +coeff( 30)    *x24    *x43*x51
     4  +coeff( 31)    *x21*x31    *x51
     5  +coeff( 32)    *x21    *x41*x51
     6  +coeff( 33)        *x31    *x52
     7  +coeff( 34)    *x22*x32        
     8  +coeff( 35)    *x21*x31*x42    
      p_sl_fp     =p_sl_fp     
     9  +coeff( 36)    *x22*x31    *x51
     1  +coeff( 37)    *x22*x31*x41*x51
     2  +coeff( 38)    *x23*x31*x42    
     3  +coeff( 39)    *x22    *x44    
     4  +coeff( 40)    *x22*x33*x43    
     5  +coeff( 41)*x12*x24*x31*x43*x51
     6  +coeff( 42)        *x33        
     7  +coeff( 43)*x11        *x41    
     8  +coeff( 44)    *x21*x31*x41    
      p_sl_fp     =p_sl_fp     
     9  +coeff( 45)            *x42*x51
     1  +coeff( 46)                *x53
     2  +coeff( 47)    *x22*x31*x41    
     3  +coeff( 48)            *x43*x51
     4  +coeff( 49)    *x22*x31*x42    
     5  +coeff( 50)    *x22    *x43    
     6  +coeff( 51)    *x23*x32*x41    
     7  +coeff( 52)    *x22*x33*x41    
     8  +coeff( 53)    *x23*x33*x42    
      p_sl_fp     =p_sl_fp     
     9  +coeff( 54)    *x22*x34*x42    
     1  +coeff( 55)    *x24    *x44    
     2  +coeff( 56)*x11*x23    *x42*x51
     3  +coeff( 57)*x11*x21    *x42*x53
     4  +coeff( 58)*x11*x23    *x44    
     5  +coeff( 59)*x12*x24    *x44    
     6  +coeff( 60)*x12*x24*x34*x44    
     7  +coeff( 61)*x11*x23*x34*x44*x53
     8  +coeff( 62)*x11*x21            
      p_sl_fp     =p_sl_fp     
     9  +coeff( 63)    *x21*x32        
     1  +coeff( 64)        *x31*x42    
     2  +coeff( 65)    *x23*x31        
     3  +coeff( 66)*x11*x21        *x51
     4  +coeff( 67)    *x22    *x41*x51
     5  +coeff( 68)            *x41*x53
     6  +coeff( 69)    *x22*x33        
     7  +coeff( 70)    *x22*x32*x41    
     8  +coeff( 71)*x11*x21    *x42    
      p_sl_fp     =p_sl_fp     
     9  +coeff( 72)    *x23    *x42    
     1  +coeff( 73)    *x21*x32*x42    
     2  +coeff( 74)    *x21    *x43*x51
     3  +coeff( 75)*x11*x22*x32        
     4  +coeff( 76)    *x22*x34        
     5  +coeff( 77)    *x22*x32*x42    
     6  +coeff( 78)    *x22*x31*x43    
     7  +coeff( 79)    *x21*x32*x43    
     8  +coeff( 80)    *x21*x31*x44    
      p_sl_fp     =p_sl_fp     
     9  +coeff( 81)    *x22*x31*x42*x51
     1  +coeff( 82)        *x33*x42*x51
     2  +coeff( 83)*x11*x21*x31    *x52
     3  +coeff( 84)    *x23*x34        
     4  +coeff( 85)    *x23*x33*x41    
     5  +coeff( 86)    *x22*x33*x42    
     6  +coeff( 87)    *x21*x33*x43    
     7  +coeff( 88)    *x22*x31*x44    
     8  +coeff( 89)        *x33*x44    
      p_sl_fp     =p_sl_fp     
     9  +coeff( 90)*x12*x22        *x51
c
      return
      end
      function l_sl_fp     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !, xmean(10)
      dimension coeff( 18)
      data ncoeff/ 17/
      data avdat/ -0.1559083E-01/
      data xmin/
     1 -0.19999E-02,-0.53425E-01,-0.19919E-01,-0.29893E-01,-0.45450E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19995E-02, 0.54621E-01, 0.19981E-01, 0.29903E-01, 0.49857E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.98081063E-02,-0.31619823E+00,-0.33432938E-01, 0.44283322E-02,
     + -0.33274420E-01,-0.24098521E-01, 0.52888338E-01,-0.49244263E-02,
     +  0.33506274E-01,-0.55952333E-02, 0.22456427E-02, 0.38813152E-02,
     + -0.46604360E-02, 0.61948183E-02, 0.56618303E-02,-0.25399506E-01,
     + -0.24844302E-01,
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
      l_sl_fp     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)                *x52
     7  +coeff(  7)    *x21    *x42    
     8  +coeff(  8)    *x21    *x41    
      l_sl_fp     =l_sl_fp     
     9  +coeff(  9)    *x21*x31*x41    
     1  +coeff( 10)    *x21*x33        
     2  +coeff( 11)            *x41    
     3  +coeff( 12)        *x31*x41    
     4  +coeff( 13)    *x21        *x51
     5  +coeff( 14)    *x21*x32        
     6  +coeff( 15)                *x53
     7  +coeff( 16)    *x21    *x43    
     8  +coeff( 17)    *x21*x33*x42    
c
      return
      end

C Forward transfer functions for right hrs with septum based on Rtest_dir.dat
c                     -JJL 8/06/03


c typical call: answer = function(x,5)
c INPUTS: x = 5 or more element array 
c              x(1)=x0  (meters)
c              x(2)=theta0 (really tan(theta0))
c              x(3)=y0   (meters)
c              x(4)=phi0 (really tan(phi0))
c              x(5)=delta (fractional value NOT percent)
c         M=5
c
c OUTPUT: units are the same as inputs
c 
c NOMENCLATURE: function name = prefix + _sr_ +suffix
c           prefixes:     x means xfinal
c                         t means thetafinal
c                         y means yfinal
c                         p means phifinal
c                         l means path length
c     
c           suffixes:     fp means target to focus
c                         q1ex means target to Q1 exit
c                         dent means target to Dipole entrance
c                         dext means target to dipole exit
c                         q3en means target to Q3 entrance
c                         q3ex means target to Q3 exit
c                         ep3  means septum entrance
c                         ep4  means one quarter the way into the septum
c                         ep5  means halfway through the septum
c                         ep6  means three quarters the way through the septum
c                         ep7  means septum exit
c
c          _sr_ is for right hrs with septum
c
c APERTURES:
c    Coordinate systems are different in the spectrometer with septum model compared to the spectrometer
c     without septum. So some apertures even in the body of the spectrometer appear to be different. Numbers
c    given here supercede any numbers given in regard to transfer functions for the spectrometers alone.
c
c     ep3: -0.1486233 < x < -0.08869672
c          -0.110 < y < 0.110
c     ep4: -0.1792231 < x < -0.1089169
c          -0.110 < y < 0.110
c     ep5: -0.2209211 < x < -0.1353789
c          -0.110 < y < 0.110
c     ep6: -0.2763536 < x < -0.1697464
c          -0.110 < y < 0.110
c     ep7: -0.3485396 < x < -0.2156404
c          -0.110 < y < 0.110
c     q1ex is a circle of radius 0.1492 m
c     dent is a trapazoid:
c                                   -5.22008 < x < -4.98099
c             -(-0.192436784*x -0.192436784) < y < -0.192436784*x -0.192436784  
c     dext is also a trapazoid: 
c                                   -0.46188 < x < 0.46188
c                   -(-0.01610808*x + 0.125) < y < -0.01610808*x + 0.125
c     q3en is a circle of radius 0.3 m
c     q3ex is a circle of radius 0.3 m
c


      function x_sr_ep3     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff(  5)
      data ncoeff/  4/
      data avdat/ -0.1164706E+00/
      data xmin/
     1 -0.19975E-02,-0.54936E-01,-0.19979E-01,-0.29994E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54959E-01, 0.19967E-01, 0.29922E-01, 0.49960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.25726720E-02, 0.32485950E-01, 0.75389570E-06, 0.19876383E-01,
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
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
c          set up monomials   functions
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
c
c                  function
c
      x_sr_ep3     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)        *x33        
     4  +coeff(  4)        *x31        
c
      return
      end
      function t_sr_ep3     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff(  7)
      data ncoeff/  6/
      data avdat/ -0.1069279E+00/
      data xmin/
     1 -0.19975E-02,-0.54936E-01,-0.19979E-01,-0.29994E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54959E-01, 0.19967E-01, 0.29922E-01, 0.49960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.11965680E-02,-0.43500770E-06, 0.13630892E-04, 0.29982750E-01,
     +  0.44200543E-04, 0.38676461E-03,
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
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x41 = x4
      x51 = x5
c
c                  function
c
      t_sr_ep3     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
c
      return
      end
      function y_sr_ep3     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff(  7)
      data ncoeff/  6/
      data avdat/  0.2605673E-03/
      data xmin/
     1 -0.19975E-02,-0.54936E-01,-0.19979E-01,-0.29994E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54959E-01, 0.19967E-01, 0.29922E-01, 0.49960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.27264560E-03,-0.59839100E-01,-0.35577011E-06, 0.18868391E-06,
     + -0.19964510E-02,-0.60292681E-04,
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
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x41 = x4
c
c                  function
c
      y_sr_ep3     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x23            
c
      return
      end
      function p_sr_ep3     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff(  5)
      data ncoeff/  4/
      data avdat/  0.2267454E-03/
      data xmin/
     1 -0.19975E-02,-0.54936E-01,-0.19979E-01,-0.29994E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54959E-01, 0.19967E-01, 0.29922E-01, 0.49960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.23819940E-03,-0.54912430E-01, 0.19500340E-03, 0.15295280E-03,
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
      x23 = x22*x2
      x31 = x3
      x41 = x4
c
c                  function
c
      p_sr_ep3     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)    *x23*x31        
c
      return
      end
      function l_sr_ep3     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff(  4)
      data ncoeff/  3/
      data avdat/  0.1089599E+01/
      data xmin/
     1 -0.19975E-02,-0.54936E-01,-0.19979E-01,-0.29994E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54959E-01, 0.19967E-01, 0.29922E-01, 0.49960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.70162540E-03, 0.16447261E-02, 0.52682054E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x41 = x4
      x42 = x41*x4
c
c                  function
c
      l_sr_ep3     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x22            
     3  +coeff(  3)            *x42    
c
      return
      end
      function x_sr_ep4     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff(  4)
      data ncoeff/  3/
      data avdat/ -0.1415698E+00/
      data xmin/
     1 -0.19975E-02,-0.54936E-01,-0.19979E-01,-0.29994E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54959E-01, 0.19967E-01, 0.29732E-01, 0.49960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.29580480E-02, 0.38780430E-01, 0.19813382E-01,
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
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
c          set up monomials   functions
      x31 = x3
      x41 = x4
c
c                  function
c
      x_sr_ep4     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)        *x31        
c
      return
      end
      function t_sr_ep4     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff(  9)
      data ncoeff/  8/
      data avdat/ -0.1268557E+00/
      data xmin/
     1 -0.19975E-02,-0.54936E-01,-0.19979E-01,-0.29994E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54959E-01, 0.19967E-01, 0.29732E-01, 0.49960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.15345194E-02,-0.13169590E-04, 0.19664601E-03, 0.30169720E-01,
     +  0.10404641E-02,-0.29002650E-03, 0.11662253E-04,-0.10712320E-03,
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
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x51 = x5
c
c                  function
c
      t_sr_ep4     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)        *x32        
     8  +coeff(  8)        *x31*x41    
c
      return
      end
      function y_sr_ep4     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff(  9)
      data ncoeff/  8/
      data avdat/  0.2217465E-03/
      data xmin/
     1 -0.19975E-02,-0.54936E-01,-0.19979E-01,-0.29994E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54959E-01, 0.19967E-01, 0.29732E-01, 0.49960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.23665800E-03,-0.71826520E-01, 0.56318140E-07, 0.26447780E-06,
     + -0.19940764E-02, 0.20651500E-03, 0.95339082E-05, 0.98902620E-04,
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
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x41 = x4
c
c                  function
c
      y_sr_ep4     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x23*x31        
     8  +coeff(  8)    *x21*x31        
c
      return
      end
      function p_sr_ep4     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff(  7)
      data ncoeff/  6/
      data avdat/  0.1496922E-03/
      data xmin/
     1 -0.19975E-02,-0.54936E-01,-0.19979E-01,-0.29994E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54959E-01, 0.19967E-01, 0.29732E-01, 0.49960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.16031574E-03,-0.54008760E-01, 0.62500410E-03, 0.13130090E-02,
     +  0.28057102E-03, 0.36220284E-03,
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
      x23 = x22*x2
      x31 = x3
      x41 = x4
      x42 = x41*x4
c
c                  function
c
      p_sr_ep4     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21*x31        
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)    *x21    *x42    
     6  +coeff(  6)    *x23*x31*x41    
c
      return
      end
      function l_sr_ep4     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff(  7)
      data ncoeff/  6/
      data avdat/  0.1308974E+01/
      data xmin/
     1 -0.19975E-02,-0.54936E-01,-0.19979E-01,-0.29994E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54959E-01, 0.19967E-01, 0.29732E-01, 0.49960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.64883050E-03, 0.10586522E-02, 0.19794914E-02,-0.39006231E-03,
     +  0.29099650E-03, 0.68285514E-03,
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
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
c
c                  function
c
      l_sr_ep4     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x22            
     4  +coeff(  4)        *x31*x41    
     5  +coeff(  5)    *x22*x31        
     6  +coeff(  6)        *x33        
c
      return
      end
      function x_sr_ep5     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff(  4)
      data ncoeff/  3/
      data avdat/ -0.1722790E+00/
      data xmin/
     1 -0.19975E-02,-0.54936E-01,-0.19979E-01,-0.29994E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54959E-01, 0.19967E-01, 0.29732E-01, 0.49960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.34938312E-02, 0.44696263E-01, 0.19683990E-01,
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
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
c          set up monomials   functions
      x31 = x3
      x41 = x4
c
c                  function
c
      x_sr_ep5     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)        *x31        
c
      return
      end
      function t_sr_ep5     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff(  9)
      data ncoeff/  8/
      data avdat/ -0.1636614E+00/
      data xmin/
     1 -0.19975E-02,-0.54936E-01,-0.19979E-01,-0.29994E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54959E-01, 0.19967E-01, 0.29732E-01, 0.49960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.21210840E-02,-0.23799480E-04,-0.46609140E-04, 0.29467193E-01,
     +  0.28855800E-02,-0.87882754E-05,-0.15406074E-02, 0.43505911E-04,
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
      x51 = x5
c
c                  function
c
      t_sr_ep5     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)*x11                
     7  +coeff(  7)    *x22            
     8  +coeff(  8)    *x21*x31        
c
      return
      end
      function y_sr_ep5     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff(  9)
      data ncoeff/  8/
      data avdat/  0.2530935E-03/
      data xmin/
     1 -0.19975E-02,-0.54936E-01,-0.19979E-01,-0.29994E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54959E-01, 0.19967E-01, 0.29732E-01, 0.49960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.26953584E-03,-0.83175830E-01,-0.45140510E-06, 0.15610883E-05,
     + -0.19908980E-02, 0.34872000E-03, 0.10074813E-04, 0.15939974E-03,
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
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x41 = x4
c
c                  function
c
      y_sr_ep5     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x23*x31        
     8  +coeff(  8)    *x21*x31        
c
      return
      end
      function p_sr_ep5     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 17)
      data ncoeff/ 16/
      data avdat/  0.1507681E-03/
      data xmin/
     1 -0.19975E-02,-0.54936E-01,-0.19979E-01,-0.29994E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54959E-01, 0.19967E-01, 0.29732E-01, 0.49960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.16207400E-03,-0.53917821E-01,-0.20272022E-05, 0.49577991E-03,
     +  0.49906870E-03,-0.47538971E-04, 0.41928260E-03, 0.15201260E-02,
     +  0.12505110E-02, 0.22944050E-02, 0.46725091E-02, 0.43516140E-02,
     +  0.21824940E-02, 0.47294432E-02, 0.29093330E-02, 0.13238752E-02,
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
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x21 = x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
c
c                  function
c
      p_sr_ep5     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x21*x31        
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)    *x21*x32        
     8  +coeff(  8)    *x21*x31*x41    
      p_sr_ep5     =p_sr_ep5     
     9  +coeff(  9)    *x21    *x42    
     1  +coeff( 10)    *x21*x32*x41    
     2  +coeff( 11)    *x21*x31*x42    
     3  +coeff( 12)    *x21    *x43    
     4  +coeff( 13)    *x21*x32*x42    
     5  +coeff( 14)    *x21*x31*x43    
     6  +coeff( 15)    *x21    *x44    
     7  +coeff( 16)    *x21*x31*x44    
c
      return
      end
      function l_sr_ep5     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff(  6)
      data ncoeff/  5/
      data avdat/  0.1519782E+01/
      data xmin/
     1 -0.19975E-02,-0.54936E-01,-0.19979E-01,-0.29994E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54959E-01, 0.19967E-01, 0.29732E-01, 0.49960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.54557324E-03, 0.33695092E-02, 0.22884630E-02, 0.88453580E-03,
     +  0.19322330E-02,
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
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
c
c                  function
c
      l_sr_ep5     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x22            
     4  +coeff(  4)    *x22*x31        
     5  +coeff(  5)        *x33        
c
      return
      end
      function x_sr_ep6     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff(  4)
      data ncoeff/  3/
      data avdat/ -0.2110171E+00/
      data xmin/
     1 -0.19975E-02,-0.54936E-01,-0.19979E-01,-0.29994E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54959E-01, 0.19967E-01, 0.29220E-01, 0.49960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.35949250E-02, 0.49929410E-01, 0.19489790E-01,
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
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
c          set up monomials   functions
      x31 = x3
      x41 = x4
c
c                  function
c
      x_sr_ep6     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)        *x31        
c
      return
      end
      function t_sr_ep6     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 14)
      data ncoeff/ 13/
      data avdat/ -0.2007603E+00/
      data xmin/
     1 -0.19975E-02,-0.54936E-01,-0.19979E-01,-0.29994E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54959E-01, 0.19967E-01, 0.29220E-01, 0.49960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.25733170E-02,-0.81933630E-04,-0.45774950E-03, 0.29046530E-01,
     +  0.46237830E-02,-0.12769383E-02,-0.47392424E-03, 0.66262400E-04,
     + -0.13738660E-02,-0.46794560E-02,-0.10728694E-01,-0.87227150E-02,
     + -0.29809151E-02,
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
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
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
c
c                  function
c
      t_sr_ep6     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)            *x42    
     8  +coeff(  8)        *x33*x41    
      t_sr_ep6     =t_sr_ep6     
     9  +coeff(  9)        *x31*x41    
     1  +coeff( 10)    *x22    *x41    
     2  +coeff( 11)    *x22    *x42    
     3  +coeff( 12)    *x24*x31*x41    
     4  +coeff( 13)    *x24*x33        
c
      return
      end
      function y_sr_ep6     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 20)
      data ncoeff/ 19/
      data avdat/  0.2892061E-03/
      data xmin/
     1 -0.19975E-02,-0.54936E-01,-0.19979E-01,-0.29994E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54959E-01, 0.19967E-01, 0.29220E-01, 0.49960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.30865200E-03,-0.94699442E-01, 0.46971020E-06, 0.48645442E-05,
     + -0.19834480E-02, 0.48355000E-03, 0.10524840E-02, 0.44425243E-04,
     +  0.17614143E-03,-0.55013670E-03, 0.12433550E-03, 0.15002170E-03,
     +  0.67469052E-03, 0.86995562E-04,-0.10277700E-04, 0.13386090E-02,
     +  0.85970730E-03, 0.21762274E-03,-0.39471942E-03,
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
c
c                  function
c
      y_sr_ep6     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x21    *x42    
     8  +coeff(  8)    *x23*x31        
      y_sr_ep6     =y_sr_ep6     
     9  +coeff(  9)    *x23*x31*x41    
     1  +coeff( 10)    *x21*x33*x42    
     2  +coeff( 11)    *x23    *x43    
     3  +coeff( 12)    *x21*x31        
     4  +coeff( 13)    *x21*x31*x41    
     5  +coeff( 14)    *x21*x33        
     6  +coeff( 15)    *x23    *x41    
     7  +coeff( 16)    *x21*x31*x42    
     8  +coeff( 17)    *x21    *x43    
      y_sr_ep6     =y_sr_ep6     
     9  +coeff( 18)    *x21*x34        
     1  +coeff( 19)    *x21*x32*x42    
c
      return
      end
      function p_sr_ep6     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 26)
      data ncoeff/ 25/
      data avdat/  0.1641316E-03/
      data xmin/
     1 -0.19975E-02,-0.54936E-01,-0.19979E-01,-0.29994E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54959E-01, 0.19967E-01, 0.29220E-01, 0.49960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.17517313E-03,-0.54103010E-01,-0.34445790E-05, 0.13715540E-04,
     +  0.12104370E-03,-0.75596430E-03, 0.25957204E-04,-0.42141910E-04,
     + -0.71577363E-04, 0.56062132E-03, 0.22496990E-02,-0.41677160E-04,
     +  0.15101630E-02,-0.11270600E-04, 0.39476351E-02, 0.83422030E-02,
     + -0.60501050E-04, 0.84579354E-02, 0.50489143E-02, 0.93695390E-02,
     +  0.65457530E-02, 0.30863800E-03, 0.44516450E-04, 0.26603310E-02,
     +  0.12761140E-02,
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
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
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
      p_sr_ep6     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)    *x21*x31        
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21        *x51
      p_sr_ep6     =p_sr_ep6     
     9  +coeff(  9)    *x23            
     1  +coeff( 10)    *x21*x32        
     2  +coeff( 11)    *x21*x31*x41    
     3  +coeff( 12)        *x32*x41    
     4  +coeff( 13)    *x21    *x42    
     5  +coeff( 14)    *x21        *x52
     6  +coeff( 15)    *x21*x32*x41    
     7  +coeff( 16)    *x21*x31*x42    
     8  +coeff( 17)        *x32*x42    
      p_sr_ep6     =p_sr_ep6     
     9  +coeff( 18)    *x21    *x43    
     1  +coeff( 19)    *x21*x32*x42    
     2  +coeff( 20)    *x21*x31*x43    
     3  +coeff( 21)    *x21    *x44    
     4  +coeff( 22)    *x23*x33        
     5  +coeff( 23)    *x21*x32*x43    
     6  +coeff( 24)    *x21*x31*x44    
     7  +coeff( 25)    *x21*x33*x43    
c
      return
      end
      function l_sr_ep6     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff(  6)
      data ncoeff/  5/
      data avdat/  0.1732274E+01/
      data xmin/
     1 -0.19975E-02,-0.54936E-01,-0.19979E-01,-0.29994E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54959E-01, 0.19967E-01, 0.29220E-01, 0.49960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.52108831E-03, 0.54939980E-02, 0.25783370E-02, 0.13462450E-02,
     +  0.28910690E-02,
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
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
c
c                  function
c
      l_sr_ep6     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x22            
     4  +coeff(  4)    *x22*x31        
     5  +coeff(  5)        *x33        
c
      return
      end
      function x_sr_ep7     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff(  9)
      data ncoeff/  8/
      data avdat/ -0.2604133E+00/
      data xmin/
     1 -0.19975E-02,-0.54936E-01,-0.19979E-01,-0.29994E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54959E-01, 0.19967E-01, 0.29220E-01, 0.49960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.45087570E-02, 0.54825901E-01, 0.78956620E-03, 0.24934150E-03,
     + -0.90968930E-03, 0.18677140E-01, 0.23844621E-02,-0.38922110E-02,
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
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
c
c                  function
c
      x_sr_ep7     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x22*x31        
     4  +coeff(  4)        *x33        
     5  +coeff(  5)    *x22*x33        
     6  +coeff(  6)        *x31        
     7  +coeff(  7)                *x51
     8  +coeff(  8)    *x22    *x42    
c
      return
      end
      function t_sr_ep7     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 12)
      data ncoeff/ 11/
      data avdat/ -0.2220073E+00/
      data xmin/
     1 -0.19975E-02,-0.54936E-01,-0.19979E-01,-0.29994E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54959E-01, 0.19967E-01, 0.29220E-01, 0.49960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.31429620E-02,-0.42625630E-04, 0.31703170E-05, 0.29546981E-01,
     +  0.56824300E-02,-0.20011013E-02, 0.41428440E-04,-0.65158060E-02,
     + -0.10814744E-01,-0.14070883E-01,-0.25723173E-02,
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
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
c
c                  function
c
      t_sr_ep7     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x22    *x41    
      t_sr_ep7     =t_sr_ep7     
     9  +coeff(  9)    *x22*x31*x41    
     1  +coeff( 10)    *x22    *x42    
     2  +coeff( 11)    *x22*x31        
c
      return
      end
      function y_sr_ep7     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 25)
      data ncoeff/ 24/
      data avdat/  0.4561184E-03/
      data xmin/
     1 -0.19975E-02,-0.54936E-01,-0.19979E-01,-0.29994E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54959E-01, 0.19967E-01, 0.29220E-01, 0.49960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.47642710E-03,-0.10728790E+00, 0.30013690E-05,-0.19756610E-02,
     + -0.16721343E-03, 0.53434702E-03, 0.21754880E-02, 0.39046080E-04,
     + -0.35328532E-05, 0.89380494E-03, 0.27535611E-02,-0.22835130E-04,
     + -0.10193834E-02,-0.19453621E-03,-0.14468390E-02, 0.15013440E-03,
     +  0.19373423E-03, 0.30433392E-03,-0.28372242E-04, 0.75604191E-04,
     + -0.54864543E-04, 0.23967630E-02, 0.30857560E-02, 0.82873490E-03,
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
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
c
c                  function
c
      y_sr_ep7     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)    *x21*x32        
     7  +coeff(  7)    *x21    *x42    
     8  +coeff(  8)    *x23*x31        
      y_sr_ep7     =y_sr_ep7     
     9  +coeff(  9)    *x23    *x41    
     1  +coeff( 10)    *x21*x32*x41    
     2  +coeff( 11)    *x21    *x43    
     3  +coeff( 12)    *x23*x31*x41    
     4  +coeff( 13)    *x21*x33*x41    
     5  +coeff( 14)    *x23    *x42    
     6  +coeff( 15)    *x21*x32*x42    
     7  +coeff( 16)    *x23*x31*x42    
     8  +coeff( 17)    *x21*x33*x42    
      y_sr_ep7     =y_sr_ep7     
     9  +coeff( 18)    *x23*x34*x41    
     1  +coeff( 19)    *x23*x34*x42    
     2  +coeff( 20)    *x21*x31        
     3  +coeff( 21)    *x21        *x51
     4  +coeff( 22)    *x21*x31*x41    
     5  +coeff( 23)    *x21*x31*x42    
     6  +coeff( 24)    *x21    *x44    
c
      return
      end
      function p_sr_ep7     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 30)
      data ncoeff/ 29/
      data avdat/  0.2245670E-03/
      data xmin/
     1 -0.19975E-02,-0.54936E-01,-0.19979E-01,-0.29994E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54959E-01, 0.19967E-01, 0.29220E-01, 0.49960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.23666140E-03,-0.53710170E-01, 0.21232370E-04, 0.37418340E-05,
     +  0.23260770E-03,-0.46780082E-03, 0.41252322E-04,-0.10180024E-03,
     +  0.16137411E-04,-0.50375281E-04, 0.81289804E-03, 0.32430970E-02,
     + -0.67656663E-04, 0.31188153E-02,-0.41030380E-04,-0.26458974E-04,
     + -0.13430021E-05,-0.31441260E-03, 0.80059384E-04, 0.39861993E-02,
     +  0.10717120E-01,-0.77098630E-04, 0.86225224E-02,-0.53018800E-04,
     +  0.37536950E-02, 0.87895030E-02, 0.58236964E-02, 0.19951321E-03,
     +  0.58506301E-03,
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
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_sr_ep7     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21*x31        
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21        *x51
      p_sr_ep7     =p_sr_ep7     
     9  +coeff(  9)            *x41*x51
     1  +coeff( 10)    *x23            
     2  +coeff( 11)    *x21*x32        
     3  +coeff( 12)    *x21*x31*x41    
     4  +coeff( 13)        *x32*x41    
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)        *x31*x41*x51
     7  +coeff( 16)    *x21        *x52
     8  +coeff( 17)                *x53
      p_sr_ep7     =p_sr_ep7     
     9  +coeff( 18)    *x23*x31        
     1  +coeff( 19)*x11    *x32        
     2  +coeff( 20)    *x21*x32*x41    
     3  +coeff( 21)    *x21*x31*x42    
     4  +coeff( 22)        *x32*x42    
     5  +coeff( 23)    *x21    *x43    
     6  +coeff( 24)    *x21    *x41*x52
     7  +coeff( 25)    *x21*x32*x42    
     8  +coeff( 26)    *x21*x31*x43    
      p_sr_ep7     =p_sr_ep7     
     9  +coeff( 27)    *x21    *x44    
     1  +coeff( 28)    *x23    *x41*x51
     2  +coeff( 29)    *x23*x33        
c
      return
      end
      function l_sr_ep7     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff(  5)
      data ncoeff/  4/
      data avdat/  0.1962772E+01/
      data xmin/
     1 -0.19975E-02,-0.54936E-01,-0.19979E-01,-0.29994E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54959E-01, 0.19967E-01, 0.29220E-01, 0.49960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.78466650E-03, 0.10955550E-01, 0.20659870E-02, 0.52067660E-02,
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
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
c
c                  function
c
      l_sr_ep7     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x22*x31        
     4  +coeff(  4)        *x33        
c
      return
      end
      function x_sr_q1ex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 10)
      data ncoeff/  9/
      data avdat/ -0.2574444E-02/
      data xmin/
     1 -0.19972E-02,-0.54501E-01,-0.19979E-01,-0.29899E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54848E-01, 0.19967E-01, 0.29220E-01, 0.49885E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.30085991E-02, 0.14510660E+00,-0.35522150E-02, 0.32514860E-02,
     + -0.60296720E-02,-0.19442860E-02,-0.73312770E-02,-0.77480940E-02,
     +  0.10292650E-02,
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
c
c                  function
c
      x_sr_q1ex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)    *x21        *x51
     5  +coeff(  5)    *x23    *x42    
     6  +coeff(  6)    *x21*x31        
     7  +coeff(  7)    *x21*x31*x41    
     8  +coeff(  8)    *x21    *x42    
      x_sr_q1ex    =x_sr_q1ex    
     9  +coeff(  9)*x11                
c
      return
      end
      function t_sr_q1ex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 14)
      data ncoeff/ 13/
      data avdat/  0.1792072E-03/
      data xmin/
     1 -0.19972E-02,-0.54501E-01,-0.19979E-01,-0.29899E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54848E-01, 0.19967E-01, 0.29220E-01, 0.49885E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.21492030E-03,-0.91769583E-02,-0.82579744E-03,-0.32899670E-03,
     +  0.29217420E-02,-0.31644100E-02,-0.27577872E-03,-0.22071290E-02,
     + -0.34522500E-03,-0.37294061E-03,-0.20241080E-02,-0.15643793E-02,
     + -0.16467650E-02,
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
c
c                  function
c
      t_sr_q1ex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)    *x21    *x42    
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)    *x21*x31*x41    
      t_sr_q1ex    =t_sr_q1ex    
     9  +coeff(  9)    *x23            
     1  +coeff( 10)    *x21*x32        
     2  +coeff( 11)    *x21    *x43    
     3  +coeff( 12)    *x23*x31*x42    
     4  +coeff( 13)    *x21*x33*x42    
c
      return
      end
      function y_sr_q1ex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 43)
      data ncoeff/ 42/
      data avdat/ -0.8513198E-02/
      data xmin/
     1 -0.19972E-02,-0.54501E-01,-0.19979E-01,-0.29899E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54848E-01, 0.19967E-01, 0.29220E-01, 0.49885E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.44633880E-02,-0.36919920E-04, 0.30132240E-01, 0.15707650E+00,
     +  0.17367031E-01,-0.68896161E-02, 0.41438190E-03,-0.19374421E-02,
     + -0.64761172E-02, 0.11590830E-03,-0.88661123E-03,-0.96695290E-03,
     + -0.90681420E-02, 0.11994220E-02, 0.11783424E-01, 0.66695640E-04,
     + -0.38983260E-01, 0.13444874E-02,-0.12198751E+00,-0.13163261E+00,
     + -0.27367090E-01,-0.10684760E-02,-0.40516260E-02,-0.47314334E-01,
     + -0.95891490E-01, 0.11629020E-01,-0.62598690E-03,-0.51901913E-02,
     + -0.22618500E-01, 0.70611713E-03,-0.25704770E-01,-0.93377614E-02,
     +  0.20769860E-01,-0.51547610E-02, 0.12818680E+00,-0.58523672E-02,
     + -0.17220070E-01,-0.56133340E-02, 0.10203640E+00, 0.55143260E-01,
     + -0.60996674E-01, 0.17670660E-01,
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
      y_sr_q1ex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)            *x42    
     8  +coeff(  8)            *x41*x51
      y_sr_q1ex    =y_sr_q1ex    
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)        *x32*x41    
     2  +coeff( 11)        *x31*x42    
     3  +coeff( 12)            *x43    
     4  +coeff( 13)    *x22*x31*x41    
     5  +coeff( 14)        *x33*x41    
     6  +coeff( 15)    *x22    *x42    
     7  +coeff( 16)    *x24*x31        
     8  +coeff( 17)    *x22    *x43    
      y_sr_q1ex    =y_sr_q1ex    
     9  +coeff( 18)    *x24*x32        
     1  +coeff( 19)    *x22*x31*x43    
     2  +coeff( 20)    *x22    *x44    
     3  +coeff( 21)    *x22*x31*x44    
     4  +coeff( 22)                *x52
     5  +coeff( 23)    *x22*x32        
     6  +coeff( 24)    *x22*x31*x42    
     7  +coeff( 25)    *x22*x32*x42    
     8  +coeff( 26)    *x24*x32*x41    
      y_sr_q1ex    =y_sr_q1ex    
     9  +coeff( 27)        *x31    *x51
     1  +coeff( 28)    *x22*x31        
     2  +coeff( 29)    *x22*x32*x41    
     3  +coeff( 30)    *x22    *x42*x51
     4  +coeff( 31)    *x22*x33*x41    
     5  +coeff( 32)    *x22    *x43*x51
     6  +coeff( 33)    *x24*x31*x42    
     7  +coeff( 34)*x11*x23    *x41*x51
     8  +coeff( 35)    *x22*x32*x44    
      y_sr_q1ex    =y_sr_q1ex    
     9  +coeff( 36)*x11*x24    *x41*x51
     1  +coeff( 37)    *x22*x34*x43    
     2  +coeff( 38)*x12*x23    *x41*x51
     3  +coeff( 39)    *x24*x33*x43    
     4  +coeff( 40)    *x24*x32*x44    
     5  +coeff( 41)    *x22*x34*x44    
     6  +coeff( 42)*x11*x23*x33*x43    
c
      return
      end
      function p_sr_q1ex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 18)
      data ncoeff/ 17/
      data avdat/ -0.3824610E-02/
      data xmin/
     1 -0.19972E-02,-0.54501E-01,-0.19979E-01,-0.29899E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54848E-01, 0.19967E-01, 0.29220E-01, 0.49885E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.18801010E-02, 0.92157530E-02, 0.68227440E-01, 0.92090950E-02,
     + -0.27860074E-02, 0.14326070E-03,-0.82084452E-02,-0.84861350E-02,
     +  0.87969694E-02,-0.22103330E-02,-0.37900132E-02,-0.14807210E-01,
     + -0.17869620E-01,-0.74947632E-01,-0.46234980E-01, 0.17540090E-03,
     + -0.44643480E-02,
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
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
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
c
c                  function
c
      p_sr_q1ex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x22    *x41    
     8  +coeff(  8)    *x22    *x42    
      p_sr_q1ex    =p_sr_q1ex    
     9  +coeff(  9)    *x24*x31*x41    
     1  +coeff( 10)            *x41*x51
     2  +coeff( 11)    *x22*x31        
     3  +coeff( 12)    *x22*x31*x41    
     4  +coeff( 13)    *x22    *x44    
     5  +coeff( 14)    *x24*x31*x43    
     6  +coeff( 15)    *x24    *x44    
     7  +coeff( 16)        *x32        
     8  +coeff( 17)    *x22*x32        
c
      return
      end
      function l_sr_q1ex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 10)
      data ncoeff/  9/
      data avdat/  0.3943433E+01/
      data xmin/
     1 -0.19972E-02,-0.54501E-01,-0.19979E-01,-0.29899E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54848E-01, 0.19967E-01, 0.29220E-01, 0.49885E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.16952070E-02, 0.21620572E-02, 0.48303930E-02, 0.38695900E-02,
     +  0.33554554E-03, 0.33365320E-02, 0.74823520E-03,-0.54391071E-03,
     + -0.26938480E-02,
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
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
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
c
c                  function
c
      l_sr_q1ex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)        *x31*x41    
     6  +coeff(  6)            *x42    
     7  +coeff(  7)            *x41*x51
     8  +coeff(  8)    *x22    *x41    
      l_sr_q1ex    =l_sr_q1ex    
     9  +coeff(  9)    *x24    *x44    
c
      return
      end
      function x_sr_dent    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 14)
      data ncoeff/ 13/
      data avdat/ -0.5098287E+01/
      data xmin/
     1 -0.19972E-02,-0.54501E-01,-0.19979E-01,-0.29899E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19884E-02, 0.54848E-01, 0.19967E-01, 0.29220E-01, 0.49885E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.21937360E-02,-0.12602543E+00, 0.34738490E-03,-0.69291090E-02,
     +  0.21979310E-01, 0.36465651E-02, 0.11515430E-02, 0.16899714E-01,
     +  0.14750410E-02, 0.32325820E-02, 0.28014670E-01, 0.25127730E-01,
     +  0.95651233E-02,
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
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
c
c                  function
c
      x_sr_dent    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)    *x21        *x51
     5  +coeff(  5)    *x21    *x42    
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)    *x21*x31*x41    
      x_sr_dent    =x_sr_dent    
     9  +coeff(  9)*x11                
     1  +coeff( 10)    *x21*x32        
     2  +coeff( 11)    *x21*x31*x42    
     3  +coeff( 12)    *x21    *x43    
     4  +coeff( 13)    *x21*x34*x41    
c
      return
      end
      function t_sr_dent    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 13)
      data ncoeff/ 12/
      data avdat/  0.1298947E+01/
      data xmin/
     1 -0.19972E-02,-0.54501E-01,-0.19979E-01,-0.29899E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19884E-02, 0.54848E-01, 0.19967E-01, 0.29220E-01, 0.49885E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.53297064E-03, 0.55919900E-01,-0.55810520E-02, 0.35877470E-02,
     +  0.30242344E-02,-0.73934951E-02, 0.73339883E-03,-0.89417460E-03,
     + -0.16815784E-02,-0.28431040E-02,-0.16110253E-03,-0.93725612E-02,
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
c
c                  function
c
      t_sr_dent    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x23    *x42    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)    *x21*x31*x41    
     7  +coeff(  7)                *x51
     8  +coeff(  8)*x11                
      t_sr_dent    =t_sr_dent    
     9  +coeff(  9)    *x21*x31        
     1  +coeff( 10)    *x21    *x41    
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x21    *x42    
c
      return
      end
      function y_sr_dent    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 81)
      data ncoeff/ 80/
      data avdat/ -0.5781254E-02/
      data xmin/
     1 -0.19972E-02,-0.54501E-01,-0.19979E-01,-0.29899E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19884E-02, 0.54848E-01, 0.19967E-01, 0.29220E-01, 0.49885E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.33770443E-02,-0.53993810E-03, 0.66151130E-02, 0.11161850E+00,
     +  0.17902344E-01,-0.70137400E-02, 0.28497680E-03, 0.56922701E-02,
     +  0.29516974E-02, 0.25242790E-02, 0.16551030E-01, 0.12690121E-02,
     + -0.13261680E-01,-0.22553620E-02,-0.33018660E-02,-0.36436773E-02,
     + -0.14762711E-02,-0.22325640E-01,-0.73844194E-02,-0.74120662E-02,
     + -0.21701222E-01,-0.49779582E-02, 0.14263550E-02, 0.22941510E-03,
     + -0.82602490E-01, 0.26137900E-01,-0.95567882E-01, 0.19965292E-01,
     + -0.91685093E-02, 0.11017541E-01,-0.15123290E+00,-0.26528320E-02,
     +  0.10629723E-02, 0.13696890E-02, 0.38959140E-02,-0.74797272E-02,
     + -0.14541460E-02,-0.54806380E-02,-0.16219930E-01,-0.92123580E-01,
     + -0.17839441E-01,-0.89974940E-01,-0.77013190E-01,-0.13390061E-01,
     +  0.59387722E-03,-0.31293820E-03,-0.30689220E-03,-0.12907792E-02,
     + -0.36974630E-02, 0.24035380E-03,-0.18573520E-02,-0.60852400E-03,
     +  0.30289310E-02, 0.47103760E-03, 0.24000970E-01, 0.39845321E-01,
     + -0.20750374E-02, 0.46576750E-02,-0.40387650E-02,-0.19914310E-01,
     + -0.16899921E-01, 0.66932400E-01, 0.25745402E-03,-0.12524850E-01,
     +  0.57548941E-02, 0.47138541E-01, 0.22995094E-02,-0.64773983E-02,
     + -0.94764650E-01, 0.94450123E-01,-0.11453174E-01, 0.10085950E-01,
     +  0.74121281E-01,-0.28021184E-02, 0.11840360E+00,-0.10894041E-01,
     + -0.18348911E-01, 0.11603950E-01,-0.36307130E-01, 0.11875580E+00,
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
      y_sr_dent    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)        *x32        
     8  +coeff(  8)    *x21    *x41    
      y_sr_dent    =y_sr_dent    
     9  +coeff(  9)            *x42    
     1  +coeff( 10)        *x31    *x51
     2  +coeff( 11)            *x41*x51
     3  +coeff( 12)                *x52
     4  +coeff( 13)    *x22    *x41    
     5  +coeff( 14)        *x32*x41    
     6  +coeff( 15)        *x31*x42    
     7  +coeff( 16)            *x43    
     8  +coeff( 17)            *x42*x51
      y_sr_dent    =y_sr_dent    
     9  +coeff( 18)    *x22*x31*x41    
     1  +coeff( 19)        *x33*x41    
     2  +coeff( 20)    *x22    *x42    
     3  +coeff( 21)        *x32*x42    
     4  +coeff( 22)            *x44    
     5  +coeff( 23)    *x22    *x41*x51
     6  +coeff( 24)    *x24*x31        
     7  +coeff( 25)    *x22    *x43    
     8  +coeff( 26)    *x24    *x42    
      y_sr_dent    =y_sr_dent    
     9  +coeff( 27)    *x22    *x44    
     1  +coeff( 28)        *x32*x44    
     2  +coeff( 29)    *x24*x31*x42    
     3  +coeff( 30)    *x22*x32*x43    
     4  +coeff( 31)    *x22*x31*x44    
     5  +coeff( 32)    *x22    *x44*x51
     6  +coeff( 33)    *x24*x34        
     7  +coeff( 34)    *x21*x31        
     8  +coeff( 35)        *x31*x41    
      y_sr_dent    =y_sr_dent    
     9  +coeff( 36)    *x22*x31        
     1  +coeff( 37)        *x31*x41*x51
     2  +coeff( 38)    *x22*x32        
     3  +coeff( 39)        *x31*x43    
     4  +coeff( 40)    *x22*x31*x43    
     5  +coeff( 41)    *x22    *x43*x51
     6  +coeff( 42)    *x24*x32*x42    
     7  +coeff( 43)    *x24    *x44    
     8  +coeff( 44)*x12*x21    *x43*x51
      y_sr_dent    =y_sr_dent    
     9  +coeff( 45)    *x21        *x51
     1  +coeff( 46)*x11*x21            
     2  +coeff( 47)        *x33        
     3  +coeff( 48)    *x22        *x51
     4  +coeff( 49)            *x41*x52
     5  +coeff( 50)*x11*x21        *x51
     6  +coeff( 51)        *x31*x42*x51
     7  +coeff( 52)*x11*x22*x31        
     8  +coeff( 53)            *x41*x54
      y_sr_dent    =y_sr_dent    
     9  +coeff( 54)*x11*x22*x31*x41    
     1  +coeff( 55)        *x34*x42    
     2  +coeff( 56)        *x33*x43    
     3  +coeff( 57)    *x24*x31    *x51
     4  +coeff( 58)*x11*x22    *x41*x51
     5  +coeff( 59)*x11*x21    *x42*x51
     6  +coeff( 60)    *x24*x32*x41    
     7  +coeff( 61)    *x22*x33*x42    
     8  +coeff( 62)    *x24    *x43    
      y_sr_dent    =y_sr_dent    
     9  +coeff( 63)        *x34*x43    
     1  +coeff( 64)*x11*x23    *x41*x51
     2  +coeff( 65)*x11*x21*x32*x41*x51
     3  +coeff( 66)    *x22    *x43*x52
     4  +coeff( 67)    *x23*x31    *x53
     5  +coeff( 68)    *x22    *x41*x54
     6  +coeff( 69)    *x22*x34*x42    
     7  +coeff( 70)    *x22*x32*x44    
     8  +coeff( 71)*x11*x24    *x41*x51
      y_sr_dent    =y_sr_dent    
     9  +coeff( 72)    *x22*x34*x43    
     1  +coeff( 73)    *x24*x31*x44    
     2  +coeff( 74)*x11*x21*x32*x44    
     3  +coeff( 75)    *x22*x33*x44    
     4  +coeff( 76)*x12*x21*x31*x42*x51
     5  +coeff( 77)    *x22*x32*x43*x52
     6  +coeff( 78)*x11*x23    *x41*x53
     7  +coeff( 79)    *x22    *x43*x54
     8  +coeff( 80)    *x24*x34*x42    
c
      return
      end
      function p_sr_dent    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 24)
      data ncoeff/ 23/
      data avdat/  0.1392198E-02/
      data xmin/
     1 -0.19972E-02,-0.54501E-01,-0.19979E-01,-0.29899E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19884E-02, 0.54848E-01, 0.19967E-01, 0.29220E-01, 0.49885E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.71128120E-03,-0.67318151E-02,-0.23744640E-01,-0.18127660E-02,
     +  0.85982480E-03,-0.35031030E-02, 0.44305240E-02, 0.73141430E-04,
     +  0.20338030E-03,-0.11086163E-04,-0.46981751E-03, 0.70281820E-03,
     +  0.62512060E-03, 0.32635760E-02, 0.99758581E-04,-0.16703880E-02,
     +  0.33586603E-03,-0.69298333E-03, 0.18663860E-02, 0.30733920E-02,
     +  0.18149350E-02, 0.85733542E-02, 0.84199090E-02,
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
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
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
c
c                  function
c
      p_sr_dent    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)            *x41*x51
     8  +coeff(  8)    *x24    *x42    
      p_sr_dent    =p_sr_dent    
     9  +coeff(  9)    *x21            
     1  +coeff( 10)        *x31*x41    
     2  +coeff( 11)    *x21        *x51
     3  +coeff( 12)        *x31    *x51
     4  +coeff( 13)                *x52
     5  +coeff( 14)    *x22*x31*x41    
     6  +coeff( 15)    *x22*x31        
     7  +coeff( 16)    *x22    *x41    
     8  +coeff( 17)    *x21    *x42    
      p_sr_dent    =p_sr_dent    
     9  +coeff( 18)    *x22        *x51
     1  +coeff( 19)    *x23    *x41    
     2  +coeff( 20)    *x22    *x42    
     3  +coeff( 21)    *x22    *x44    
     4  +coeff( 22)    *x24*x31*x42    
     5  +coeff( 23)    *x24    *x43    
c
      return
      end
      function l_sr_dent    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff(  4)
      data ncoeff/  3/
      data avdat/  0.1075976E+02/
      data xmin/
     1 -0.19972E-02,-0.54501E-01,-0.19979E-01,-0.29899E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19884E-02, 0.54848E-01, 0.19967E-01, 0.29220E-01, 0.49885E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.71754311E-02,-0.23651260E+00, 0.16966061E-01,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
c
c                  function
c
      l_sr_dent    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x22            
c
      return
      end
      function x_sr_dext     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 35)
      data ncoeff/ 34/
      data avdat/ -0.2336794E-02/
      data xmin/
     1 -0.19972E-02,-0.54501E-01,-0.19979E-01,-0.29899E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19884E-02, 0.54848E-01, 0.19967E-01, 0.29220E-01, 0.49885E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.40092320E-02, 0.34526693E+00, 0.13303080E+00, 0.66025620E-02,
     +  0.45855360E-01,-0.36773380E-01, 0.75839660E-02, 0.10891420E-01,
     + -0.10282120E-02,-0.29183880E-01,-0.54791810E-02,-0.33827380E-02,
     + -0.76843863E-02,-0.82020033E-02,-0.81015750E-01,-0.93387514E-01,
     + -0.84693040E-03,-0.39410950E-02,-0.48579390E-02,-0.11765431E-02,
     + -0.10713454E-02,-0.19210143E-02,-0.36073650E-01,-0.51151350E-02,
     + -0.11977761E+00,-0.95100753E-01,-0.19154380E-02, 0.20599694E-01,
     + -0.22464390E-01,-0.54260270E-01,-0.51626130E-02,-0.34355000E-01,
     + -0.62016160E-01,-0.36275330E-01,
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
      x_sr_dext     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)    *x21    *x42    
     7  +coeff(  7)    *x23*x31*x41    
     8  +coeff(  8)    *x22            
      x_sr_dext     =x_sr_dext     
     9  +coeff(  9)    *x21*x31        
     1  +coeff( 10)    *x21*x31*x41    
     2  +coeff( 11)*x11                
     3  +coeff( 12)                *x52
     4  +coeff( 13)    *x21*x32        
     5  +coeff( 14)    *x21    *x41*x51
     6  +coeff( 15)    *x21*x31*x42    
     7  +coeff( 16)    *x21    *x43    
     8  +coeff( 17)    *x21*x34*x41    
      x_sr_dext     =x_sr_dext     
     9  +coeff( 18)            *x42    
     1  +coeff( 19)    *x21*x33        
     2  +coeff( 20)        *x31*x41    
     3  +coeff( 21)            *x41*x51
     4  +coeff( 22)    *x21*x31    *x51
     5  +coeff( 23)    *x21*x32*x41    
     6  +coeff( 24)    *x21    *x42*x51
     7  +coeff( 25)    *x21*x31*x43    
     8  +coeff( 26)    *x21    *x44    
      x_sr_dext     =x_sr_dext     
     9  +coeff( 27)    *x23        *x52
     1  +coeff( 28)    *x23    *x43    
     2  +coeff( 29)    *x21*x32*x43    
     3  +coeff( 30)    *x21*x31*x44    
     4  +coeff( 31)    *x21*x31*x41*x53
     5  +coeff( 32)    *x23*x33*x41    
     6  +coeff( 33)    *x23*x32*x42    
     7  +coeff( 34)    *x21*x34*x42    
c
      return
      end
      function t_sr_dext     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 12)
      data ncoeff/ 11/
      data avdat/  0.5336839E+00/
      data xmin/
     1 -0.19972E-02,-0.54501E-01,-0.19979E-01,-0.29899E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19884E-02, 0.54848E-01, 0.19967E-01, 0.29220E-01, 0.49885E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.53904360E-03,-0.78800424E-01, 0.22719563E-01,-0.33301864E-02,
     +  0.40684253E-02,-0.38496410E-02, 0.14043630E-01, 0.88294290E-03,
     +  0.22366500E-02,-0.12483622E-02, 0.87312562E-02,
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
      t_sr_dext     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)    *x21    *x42    
     8  +coeff(  8)*x11                
      t_sr_dext     =t_sr_dext     
     9  +coeff(  9)    *x21*x31        
     1  +coeff( 10)                *x52
     2  +coeff( 11)    *x21*x31*x41    
c
      return
      end
      function y_sr_dext     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 96)
      data ncoeff/ 95/
      data avdat/  0.4796520E-03/
      data xmin/
     1 -0.19972E-02,-0.54501E-01,-0.19979E-01,-0.29899E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19884E-02, 0.54848E-01, 0.19967E-01, 0.29220E-01, 0.49885E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.55530230E-03,-0.36494333E-01, 0.13115750E-01, 0.17326800E-01,
     + -0.70642520E-02,-0.52925031E-02, 0.77735880E-03,-0.51475290E-01,
     +  0.20992590E-02, 0.16386650E-03,-0.68741482E-02, 0.88135283E-02,
     +  0.54991020E-01, 0.59559860E-02, 0.37386552E-02,-0.94941280E-02,
     + -0.84876391E-03,-0.41887111E-01,-0.30594242E-02, 0.59717064E-02,
     + -0.63011883E-02,-0.81857964E-02,-0.16247543E-02,-0.40236022E-02,
     + -0.53426390E-02,-0.47146340E-02, 0.24162363E-02,-0.14003442E-01,
     +  0.62645480E-02, 0.13198540E-01,-0.79383890E-02, 0.10903120E-01,
     +  0.15805703E-02,-0.45070740E-02,-0.17060970E-02, 0.71103712E-02,
     +  0.69973524E-02,-0.14088303E-01, 0.74969190E-03,-0.65403774E-01,
     + -0.15634470E-02,-0.87189972E-02, 0.25956850E-02, 0.11222950E-01,
     +  0.37422490E-02, 0.24871290E-01, 0.32116010E-01,-0.12568610E+00,
     + -0.39474400E-01, 0.46223890E-01, 0.17731593E-02, 0.76051882E-03,
     +  0.50001140E-02,-0.73558674E-02,-0.10861502E-02,-0.10591790E-02,
     + -0.10359772E-02, 0.18380692E-02,-0.36019750E-02, 0.31306530E-03,
     + -0.41466820E-02,-0.32537630E-02, 0.23359990E-02,-0.78247264E-02,
     +  0.25247340E-02,-0.33821470E-02, 0.42612560E-02,-0.10517240E+00,
     + -0.36064420E-01, 0.62812264E-02,-0.16217300E-01,-0.22357760E-02,
     +  0.65760140E-01,-0.80917462E-01, 0.18220524E-02,-0.57132500E-02,
     + -0.10121642E+00,-0.80271914E-01, 0.20550280E-01, 0.11642792E+00,
     +  0.28414253E-02,-0.58733550E-02, 0.52988223E-01, 0.15439870E+00,
     + -0.50604690E-01,-0.69869750E-01, 0.18300880E-01,-0.10741911E+00,
     + -0.12114670E+00,-0.28346033E-03, 0.20806463E-02,-0.75945150E-02,
     +  0.16358693E-02,-0.32721790E-02, 0.82275290E-03,
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
      y_sr_dext     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21*x31        
     7  +coeff(  7)        *x32        
     8  +coeff(  8)    *x21    *x41    
      y_sr_dext     =y_sr_dext     
     9  +coeff(  9)        *x31*x41    
     1  +coeff( 10)            *x42    
     2  +coeff( 11)    *x21        *x51
     3  +coeff( 12)        *x31    *x51
     4  +coeff( 13)            *x41*x51
     5  +coeff( 14)                *x52
     6  +coeff( 15)    *x23            
     7  +coeff( 16)    *x22*x31        
     8  +coeff( 17)        *x33        
      y_sr_dext     =y_sr_dext     
     9  +coeff( 18)    *x22    *x41    
     1  +coeff( 19)        *x32*x41    
     2  +coeff( 20)    *x21    *x42    
     3  +coeff( 21)            *x43    
     4  +coeff( 22)    *x22        *x51
     5  +coeff( 23)    *x21*x31    *x51
     6  +coeff( 24)    *x21    *x41*x51
     7  +coeff( 25)            *x42*x51
     8  +coeff( 26)            *x41*x52
      y_sr_dext     =y_sr_dext     
     9  +coeff( 27)    *x23    *x41    
     1  +coeff( 28)    *x22*x31*x41    
     2  +coeff( 29)    *x21*x32*x41    
     3  +coeff( 30)    *x22    *x42    
     4  +coeff( 31)        *x32*x42    
     5  +coeff( 32)    *x21    *x43    
     6  +coeff( 33)            *x44    
     7  +coeff( 34)    *x22    *x41*x51
     8  +coeff( 35)            *x43*x51
      y_sr_dext     =y_sr_dext     
     9  +coeff( 36)    *x23*x31*x41    
     1  +coeff( 37)    *x23    *x42    
     2  +coeff( 38)    *x22*x31*x42    
     3  +coeff( 39)        *x33*x42    
     4  +coeff( 40)    *x22    *x43    
     5  +coeff( 41)        *x32*x43    
     6  +coeff( 42)    *x22*x31*x41*x51
     7  +coeff( 43)        *x33*x41*x51
     8  +coeff( 44)    *x22    *x42*x51
      y_sr_dext     =y_sr_dext     
     9  +coeff( 45)    *x24*x32        
     1  +coeff( 46)    *x23*x31*x42    
     2  +coeff( 47)    *x23    *x43    
     3  +coeff( 48)    *x22    *x44    
     4  +coeff( 49)    *x22    *x44*x51
     5  +coeff( 50)    *x24*x31*x43    
     6  +coeff( 51)    *x21*x32        
     7  +coeff( 52)*x11        *x41    
     8  +coeff( 53)    *x21*x31*x41    
      y_sr_dext     =y_sr_dext     
     9  +coeff( 54)        *x31*x42    
     1  +coeff( 55)        *x31*x41*x51
     2  +coeff( 56)        *x31    *x52
     3  +coeff( 57)                *x53
     4  +coeff( 58)    *x24            
     5  +coeff( 59)        *x33*x41    
     6  +coeff( 60)*x11        *x42    
     7  +coeff( 61)        *x31*x43    
     8  +coeff( 62)    *x22*x31    *x51
      y_sr_dext     =y_sr_dext     
     9  +coeff( 63)        *x31*x44    
     1  +coeff( 64)        *x31*x43*x51
     2  +coeff( 65)    *x21*x31    *x53
     3  +coeff( 66)    *x24*x31*x41    
     4  +coeff( 67)    *x23*x32*x41    
     5  +coeff( 68)    *x22*x31*x43    
     6  +coeff( 69)    *x21*x32*x43    
     7  +coeff( 70)    *x23    *x42*x51
     8  +coeff( 71)    *x22    *x43*x51
      y_sr_dext     =y_sr_dext     
     9  +coeff( 72)        *x31*x44*x51
     1  +coeff( 73)    *x24    *x43    
     2  +coeff( 74)    *x22*x31*x44    
     3  +coeff( 75)*x12*x22        *x51
     4  +coeff( 76)*x12*x21    *x41*x51
     5  +coeff( 77)    *x24*x32*x42    
     6  +coeff( 78)    *x22*x34*x42    
     7  +coeff( 79)    *x21*x34*x43    
     8  +coeff( 80)    *x22*x32*x44    
      y_sr_dext     =y_sr_dext     
     9  +coeff( 81)*x12*x21*x31*x41*x51
     1  +coeff( 82)    *x24*x34*x41    
     2  +coeff( 83)    *x24*x31*x44    
     3  +coeff( 84)    *x24*x34*x42    
     4  +coeff( 85)    *x22*x34*x44    
     5  +coeff( 86)    *x24*x34*x43    
     6  +coeff( 87)*x12*x23*x32*x41*x51
     7  +coeff( 88)    *x24*x34*x42*x51
     8  +coeff( 89)    *x24*x33*x43*x51
      y_sr_dext     =y_sr_dext     
     9  +coeff( 90)*x11*x21            
     1  +coeff( 91)    *x23*x31        
     2  +coeff( 92)    *x22*x32        
     3  +coeff( 93)    *x21*x31*x41*x51
     4  +coeff( 94)        *x31*x42*x51
     5  +coeff( 95)    *x21*x31    *x52
c
      return
      end
      function p_sr_dext     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 30)
      data ncoeff/ 29/
      data avdat/  0.9551189E-03/
      data xmin/
     1 -0.19972E-02,-0.54501E-01,-0.19979E-01,-0.29899E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19884E-02, 0.54848E-01, 0.19967E-01, 0.29220E-01, 0.49885E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.37555190E-03,-0.92669430E-02,-0.13063563E-01, 0.12240763E-02,
     + -0.37229070E-03,-0.73650590E-03,-0.56217580E-04,-0.97131934E-02,
     + -0.28184714E-03,-0.12421134E-02, 0.11366330E-02, 0.99695594E-02,
     +  0.12707200E-02, 0.91840770E-03,-0.68805470E-02, 0.18410441E-02,
     + -0.17953274E-02,-0.60153024E-03, 0.15872140E-02,-0.87826041E-03,
     + -0.11972230E-02, 0.24790011E-03,-0.22137491E-02, 0.72861790E-02,
     +  0.69327550E-02,-0.29661363E-03,-0.86757872E-03,-0.77018390E-03,
     + -0.41799770E-03,
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
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      p_sr_dext     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21*x31        
     7  +coeff(  7)        *x32        
     8  +coeff(  8)    *x21    *x41    
      p_sr_dext     =p_sr_dext     
     9  +coeff(  9)            *x42    
     1  +coeff( 10)    *x21        *x51
     2  +coeff( 11)        *x31    *x51
     3  +coeff( 12)            *x41*x51
     4  +coeff( 13)                *x52
     5  +coeff( 14)    *x23            
     6  +coeff( 15)    *x22    *x41    
     7  +coeff( 16)    *x21    *x42    
     8  +coeff( 17)    *x22        *x51
      p_sr_dext     =p_sr_dext     
     9  +coeff( 18)    *x22*x31        
     1  +coeff( 19)    *x21*x31*x41    
     2  +coeff( 20)    *x21    *x41*x51
     3  +coeff( 21)            *x42*x51
     4  +coeff( 22)    *x23    *x41    
     5  +coeff( 23)    *x22    *x41*x51
     6  +coeff( 24)    *x23    *x43    
     7  +coeff( 25)    *x23*x31*x44    
     8  +coeff( 26)        *x31*x41    
      p_sr_dext     =p_sr_dext     
     9  +coeff( 27)            *x43    
     1  +coeff( 28)        *x31*x41*x51
     2  +coeff( 29)            *x41*x52
c
      return
      end
      function l_sr_dext     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 17)
      data ncoeff/ 16/
      data avdat/  0.1734885E+02/
      data xmin/
     1 -0.19972E-02,-0.54501E-01,-0.19979E-01,-0.29899E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19884E-02, 0.54848E-01, 0.19967E-01, 0.29220E-01, 0.49885E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.21031250E-02, 0.49520020E+00, 0.98460420E-01, 0.31669530E-01,
     + -0.16785961E-02, 0.44648122E-01,-0.10148034E+00, 0.14825254E-01,
     + -0.74184370E-02,-0.47976700E-02,-0.80369800E-01,-0.13958903E-01,
     + -0.88777740E-02,-0.11788070E+00,-0.10050160E+00,-0.44408261E-01,
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
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
c
c                  function
c
      l_sr_dext     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)    *x21    *x42    
     8  +coeff(  8)    *x23*x31*x41    
      l_sr_dext     =l_sr_dext     
     9  +coeff(  9)*x11                
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)    *x21*x31*x41    
     3  +coeff( 12)    *x21*x32        
     4  +coeff( 13)    *x21    *x41*x51
     5  +coeff( 14)    *x21*x31*x42    
     6  +coeff( 15)    *x21    *x43    
     7  +coeff( 16)    *x21*x34*x41    
c
      return
      end
      function x_sr_q3en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 39)
      data ncoeff/ 38/
      data avdat/ -0.5595049E-02/
      data xmin/
     1 -0.19972E-02,-0.54501E-01,-0.19979E-01,-0.29899E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19884E-02, 0.54848E-01, 0.19967E-01, 0.29220E-01, 0.49885E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.80596350E-02, 0.21889270E+00,-0.37562483E-03, 0.13866850E+00,
     +  0.18050160E-01, 0.22287773E-02, 0.35813421E-01,-0.23674300E-01,
     +  0.24772910E-02,-0.11296550E-03,-0.58004893E-02,-0.18413914E-01,
     + -0.38600354E-02,-0.42140200E-02,-0.44866753E-02,-0.53588070E-02,
     + -0.67028291E-01,-0.58318010E-01, 0.10215210E-01,-0.16491230E-02,
     + -0.92757411E-01,-0.62358530E-01,-0.26809473E-03,-0.10997240E-02,
     +  0.22667541E-02, 0.22727470E-02,-0.18355700E-02,-0.39635080E-02,
     +  0.63137500E-02,-0.34202670E-01,-0.17277203E-02,-0.45775990E-02,
     + -0.10158820E-01,-0.50502900E-01, 0.10604600E-01,-0.26793092E-01,
     + -0.22188110E-01, 0.15742640E-01,
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
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      x_sr_q3en   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)    *x21    *x42    
      x_sr_q3en   =x_sr_q3en   
     9  +coeff(  9)    *x23*x31*x41    
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)                *x52
     3  +coeff( 12)    *x21*x31*x41    
     4  +coeff( 13)*x11                
     5  +coeff( 14)            *x42    
     6  +coeff( 15)    *x21*x32        
     7  +coeff( 16)    *x21    *x41*x51
     8  +coeff( 17)    *x21*x31*x42    
      x_sr_q3en   =x_sr_q3en   
     9  +coeff( 18)    *x21    *x43    
     1  +coeff( 19)    *x21*x34*x41    
     2  +coeff( 20)        *x31*x41    
     3  +coeff( 21)    *x21*x31*x43    
     4  +coeff( 22)    *x21    *x44    
     5  +coeff( 23)        *x31        
     6  +coeff( 24)            *x41*x51
     7  +coeff( 25)    *x22    *x41    
     8  +coeff( 26)    *x22        *x51
      x_sr_q3en   =x_sr_q3en   
     9  +coeff( 27)    *x21*x31    *x51
     1  +coeff( 28)    *x21*x33        
     2  +coeff( 29)    *x23    *x41    
     3  +coeff( 30)    *x21*x32*x41    
     4  +coeff( 31)    *x21    *x42*x51
     5  +coeff( 32)    *x24    *x41    
     6  +coeff( 33)    *x21*x33*x41    
     7  +coeff( 34)    *x21*x32*x42    
     8  +coeff( 35)    *x23*x31*x42    
      x_sr_q3en   =x_sr_q3en   
     9  +coeff( 36)    *x21*x31*x44    
     1  +coeff( 37)    *x21*x34*x43    
     2  +coeff( 38)    *x23*x31*x44*x51
c
      return
      end
      function t_sr_q3en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 12)
      data ncoeff/ 11/
      data avdat/ -0.6737196E-03/
      data xmin/
     1 -0.19972E-02,-0.54501E-01,-0.19979E-01,-0.29899E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19884E-02, 0.54848E-01, 0.19967E-01, 0.29220E-01, 0.49885E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.20929780E-03,-0.79241600E-01, 0.22705630E-01,-0.45621283E-02,
     +  0.41371220E-02,-0.51315520E-02, 0.14413760E-01, 0.88699190E-03,
     +  0.21576392E-02,-0.14410692E-02, 0.94350170E-02,
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
      t_sr_q3en   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)    *x21    *x42    
     8  +coeff(  8)*x11                
      t_sr_q3en   =t_sr_q3en   
     9  +coeff(  9)    *x21*x31        
     1  +coeff( 10)                *x52
     2  +coeff( 11)    *x21*x31*x41    
c
      return
      end
      function y_sr_q3en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 96)
      data ncoeff/ 95/
      data avdat/  0.1449694E-02/
      data xmin/
     1 -0.19972E-02,-0.54501E-01,-0.19979E-01,-0.29899E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19884E-02, 0.54848E-01, 0.19967E-01, 0.29220E-01, 0.49885E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.93785460E-04, 0.89909640E-03,-0.46059620E-01, 0.21850820E-02,
     +  0.19415610E-01,-0.80955080E-02,-0.33588272E-02, 0.98945530E-03,
     + -0.61476431E-01, 0.52051700E-02, 0.50074132E-02,-0.89644873E-02,
     +  0.10138043E-01, 0.70655220E-01, 0.75132040E-02, 0.28026732E-02,
     + -0.12654421E-01,-0.23342594E-02,-0.52386503E-01,-0.96558680E-02,
     +  0.94336430E-04,-0.16370630E-01,-0.98218042E-02,-0.13281410E-02,
     + -0.55601401E-02,-0.53595510E-02,-0.62436354E-02, 0.61651510E-02,
     + -0.30573340E-01, 0.78057893E-02,-0.88579320E-02,-0.19765650E-01,
     +  0.12258974E-01,-0.13500580E-01,-0.11289110E-01,-0.81402110E-02,
     +  0.16547180E-01, 0.14844173E-01,-0.73073930E-02, 0.14039811E-01,
     + -0.25105122E-01,-0.16491271E-02, 0.16199050E-01,-0.15934740E-02,
     +  0.28335010E-01, 0.30499570E-01,-0.61588510E-02,-0.83888390E-01,
     +  0.91638310E-03, 0.96083470E-03,-0.22308000E-03,-0.19313580E-01,
     + -0.88515100E-03,-0.14937960E-02, 0.33657120E-02,-0.44936700E-02,
     +  0.79337780E-03,-0.27286590E-01, 0.68183990E-02, 0.83414213E-02,
     + -0.26132920E-02,-0.10129192E-01, 0.97907613E-03,-0.18343850E-01,
     + -0.60211290E-01, 0.14366170E-02, 0.41569970E-02,-0.12893943E-01,
     +  0.14385112E-04,-0.20055022E-01, 0.25140710E-02,-0.34826763E-02,
     + -0.64094790E-01,-0.54766014E-02,-0.15511980E-01, 0.86635250E-01,
     +  0.22868970E-02,-0.12977532E-01, 0.12654810E+00,-0.82652300E-02,
     + -0.10482062E-01,-0.47948720E-03,-0.22959813E-02,-0.64027630E-02,
     +  0.65778603E-03,-0.30432364E-02, 0.32783010E-02, 0.31584701E-02,
     +  0.58339770E-02, 0.41856760E-02, 0.12777301E-01, 0.12465610E-01,
     + -0.11827900E-02,-0.69136070E-02, 0.20334450E-02,
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
      y_sr_q3en   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)        *x32        
      y_sr_q3en   =y_sr_q3en   
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)        *x31*x41    
     2  +coeff( 11)            *x42    
     3  +coeff( 12)    *x21        *x51
     4  +coeff( 13)        *x31    *x51
     5  +coeff( 14)            *x41*x51
     6  +coeff( 15)                *x52
     7  +coeff( 16)    *x23            
     8  +coeff( 17)    *x22*x31        
      y_sr_q3en   =y_sr_q3en   
     9  +coeff( 18)        *x33        
     1  +coeff( 19)    *x22    *x41    
     2  +coeff( 20)        *x32*x41    
     3  +coeff( 21)    *x21    *x42    
     4  +coeff( 22)            *x43    
     5  +coeff( 23)    *x22        *x51
     6  +coeff( 24)    *x21*x31    *x51
     7  +coeff( 25)    *x21    *x41*x51
     8  +coeff( 26)            *x42*x51
      y_sr_q3en   =y_sr_q3en   
     9  +coeff( 27)            *x41*x52
     1  +coeff( 28)    *x23    *x41    
     2  +coeff( 29)    *x22*x31*x41    
     3  +coeff( 30)    *x21*x32*x41    
     4  +coeff( 31)    *x22    *x42    
     5  +coeff( 32)        *x32*x42    
     6  +coeff( 33)    *x21    *x43    
     7  +coeff( 34)            *x44    
     8  +coeff( 35)    *x22    *x41*x51
      y_sr_q3en   =y_sr_q3en   
     9  +coeff( 36)            *x43*x51
     1  +coeff( 37)    *x23*x31*x41    
     2  +coeff( 38)    *x23    *x42    
     3  +coeff( 39)    *x22*x31*x42    
     4  +coeff( 40)        *x33*x42    
     5  +coeff( 41)    *x22    *x43    
     6  +coeff( 42)    *x22*x31*x41*x51
     7  +coeff( 43)    *x22    *x42*x51
     8  +coeff( 44)    *x24*x32        
      y_sr_q3en   =y_sr_q3en   
     9  +coeff( 45)    *x23*x31*x42    
     1  +coeff( 46)    *x23    *x43    
     2  +coeff( 47)    *x24*x31    *x51
     3  +coeff( 48)    *x22    *x44*x51
     4  +coeff( 49)    *x21*x32        
     5  +coeff( 50)*x11        *x41    
     6  +coeff( 51)    *x21*x31*x41    
     7  +coeff( 52)        *x31*x42    
     8  +coeff( 53)        *x31*x41*x51
      y_sr_q3en   =y_sr_q3en   
     9  +coeff( 54)                *x53
     1  +coeff( 55)    *x24            
     2  +coeff( 56)        *x33*x41    
     3  +coeff( 57)*x11        *x42    
     4  +coeff( 58)        *x31*x43    
     5  +coeff( 59)    *x21    *x42*x51
     6  +coeff( 60)        *x32*x43    
     7  +coeff( 61)    *x22*x32    *x51
     8  +coeff( 62)        *x31*x43*x51
      y_sr_q3en   =y_sr_q3en   
     9  +coeff( 63)        *x31    *x54
     1  +coeff( 64)    *x21*x32*x43    
     2  +coeff( 65)    *x22    *x44    
     3  +coeff( 66)        *x32*x44    
     4  +coeff( 67)    *x21*x33*x41*x51
     5  +coeff( 68)        *x31*x44*x51
     6  +coeff( 69)    *x22*x32*x43    
     7  +coeff( 70)    *x22*x31*x44    
     8  +coeff( 71)*x12*x22        *x51
      y_sr_q3en   =y_sr_q3en   
     9  +coeff( 72)*x12*x21    *x41*x51
     1  +coeff( 73)    *x22*x31*x43*x51
     2  +coeff( 74)*x12*x21    *x43    
     3  +coeff( 75)    *x24*x31*x43    
     4  +coeff( 76)    *x22*x33*x43    
     5  +coeff( 77)*x12        *x44    
     6  +coeff( 78)    *x24    *x44    
     7  +coeff( 79)    *x22*x32*x44    
     8  +coeff( 80)*x12*x21    *x42*x51
      y_sr_q3en   =y_sr_q3en   
     9  +coeff( 81)    *x22*x34*x42*x51
     1  +coeff( 82)*x11*x21            
     2  +coeff( 83)        *x31    *x52
     3  +coeff( 84)    *x22*x32        
     4  +coeff( 85)*x11    *x31*x41    
     5  +coeff( 86)        *x32*x41*x51
     6  +coeff( 87)    *x23*x32        
     7  +coeff( 88)    *x22*x33        
     8  +coeff( 89)    *x24    *x41    
      y_sr_q3en   =y_sr_q3en   
     9  +coeff( 90)        *x34*x41    
     1  +coeff( 91)    *x21*x31*x43    
     2  +coeff( 92)    *x21    *x44    
     3  +coeff( 93)*x11*x21    *x41*x51
     4  +coeff( 94)        *x32*x42*x51
     5  +coeff( 95)    *x21*x31    *x53
c
      return
      end
      function p_sr_q3en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 30)
      data ncoeff/ 29/
      data avdat/  0.1029573E-02/
      data xmin/
     1 -0.19972E-02,-0.54501E-01,-0.19979E-01,-0.29899E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19884E-02, 0.54848E-01, 0.19967E-01, 0.29220E-01, 0.49885E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.30900631E-03,-0.11264130E-01,-0.12536560E-01, 0.18889630E-02,
     + -0.85437064E-03,-0.11294924E-02,-0.12183130E-03,-0.11806994E-01,
     + -0.67604380E-03,-0.15118750E-02, 0.12700321E-02, 0.12946540E-01,
     +  0.16520984E-02, 0.90901240E-03,-0.10064742E-01, 0.55086380E-03,
     + -0.17493340E-02,-0.32788760E-02, 0.41467220E-03,-0.62088161E-03,
     + -0.11612680E-02,-0.13988020E-02,-0.11917350E-02, 0.69990080E-03,
     + -0.35822254E-02, 0.98596200E-03, 0.16280020E-01, 0.13609220E-01,
     +  0.59933972E-03,
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
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
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
      p_sr_q3en   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21*x31        
     7  +coeff(  7)        *x32        
     8  +coeff(  8)    *x21    *x41    
      p_sr_q3en   =p_sr_q3en   
     9  +coeff(  9)            *x42    
     1  +coeff( 10)    *x21        *x51
     2  +coeff( 11)        *x31    *x51
     3  +coeff( 12)            *x41*x51
     4  +coeff( 13)                *x52
     5  +coeff( 14)    *x23            
     6  +coeff( 15)    *x22    *x41    
     7  +coeff( 16)    *x21    *x42    
     8  +coeff( 17)    *x22        *x51
      p_sr_q3en   =p_sr_q3en   
     9  +coeff( 18)    *x22    *x41*x51
     1  +coeff( 19)    *x21*x33*x41    
     2  +coeff( 20)        *x31*x41    
     3  +coeff( 21)    *x22*x31        
     4  +coeff( 22)            *x43    
     5  +coeff( 23)    *x21    *x41*x51
     6  +coeff( 24)    *x23    *x41    
     7  +coeff( 25)    *x22    *x42*x51
     8  +coeff( 26)    *x23    *x43    
      p_sr_q3en   =p_sr_q3en   
     9  +coeff( 27)    *x23*x31*x43    
     1  +coeff( 28)    *x23    *x44    
     2  +coeff( 29)    *x21*x32        
c
      return
      end
      function l_sr_q3en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 20)
      data ncoeff/ 19/
      data avdat/  0.1836018E+02/
      data xmin/
     1 -0.19972E-02,-0.54501E-01,-0.19979E-01,-0.29899E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19884E-02, 0.54848E-01, 0.19967E-01, 0.29220E-01, 0.49885E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.21838440E-02, 0.32147800E+00, 0.32073870E-01, 0.28518160E-01,
     + -0.32873900E-03, 0.20327840E-01,-0.63310110E-01, 0.83724950E-02,
     + -0.44358642E-02,-0.29988550E-02,-0.50610363E-01, 0.22466250E-02,
     +  0.34943650E-02,-0.85392724E-02,-0.52119302E-02,-0.77154280E-01,
     + -0.65913330E-01,-0.85361780E-02,-0.25100184E-01,
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
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
c
c                  function
c
      l_sr_q3en   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)    *x21    *x42    
     8  +coeff(  8)    *x23*x31*x41    
      l_sr_q3en   =l_sr_q3en   
     9  +coeff(  9)*x11                
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)    *x21*x31*x41    
     3  +coeff( 12)            *x41    
     4  +coeff( 13)            *x42    
     5  +coeff( 14)    *x21*x32        
     6  +coeff( 15)    *x21    *x41*x51
     7  +coeff( 16)    *x21*x31*x42    
     8  +coeff( 17)    *x21    *x43    
      l_sr_q3en   =l_sr_q3en   
     9  +coeff( 18)    *x23*x32*x41    
     1  +coeff( 19)    *x21*x34*x41    
c
      return
      end
      function x_sr_q3ex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 29)
      data ncoeff/ 28/
      data avdat/ -0.1112230E-02/
      data xmin/
     1 -0.19972E-02,-0.54501E-01,-0.19979E-01,-0.29899E-01,-0.46635E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19884E-02, 0.54444E-01, 0.19967E-01, 0.29220E-01, 0.49077E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.95429340E-02, 0.72408670E-01, 0.38153870E-03, 0.31645820E+00,
     +  0.35488083E-01,-0.20998971E-01,-0.36397820E-02, 0.14655410E-02,
     + -0.24533020E-01,-0.32348360E-02, 0.15286350E-01,-0.11285392E-02,
     + -0.57281200E-02,-0.14591872E-01,-0.38719470E-02,-0.45049450E-02,
     + -0.17596630E-02,-0.12771530E-02, 0.15431200E-02,-0.10544710E-02,
     +  0.24967773E-02,-0.12099382E-01,-0.32203060E-01,-0.25489762E-01,
     + -0.31488882E-02,-0.29090870E-01,-0.20850520E-01,-0.11660343E-01,
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
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_sr_q3ex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)                *x52
     7  +coeff(  7)    *x24            
     8  +coeff(  8)    *x21    *x41    
      x_sr_q3ex    =x_sr_q3ex    
     9  +coeff(  9)    *x21    *x42    
     1  +coeff( 10)*x11                
     2  +coeff( 11)    *x22            
     3  +coeff( 12)    *x21*x31        
     4  +coeff( 13)            *x42    
     5  +coeff( 14)    *x21*x31*x41    
     6  +coeff( 15)    *x21*x32        
     7  +coeff( 16)    *x21    *x41*x51
     8  +coeff( 17)        *x31*x41    
      x_sr_q3ex    =x_sr_q3ex    
     9  +coeff( 18)            *x41*x51
     1  +coeff( 19)    *x22        *x51
     2  +coeff( 20)    *x21        *x52
     3  +coeff( 21)                *x53
     4  +coeff( 22)    *x21*x32*x41    
     5  +coeff( 23)    *x21*x31*x42    
     6  +coeff( 24)    *x21    *x43    
     7  +coeff( 25)    *x21    *x42*x51
     8  +coeff( 26)    *x21*x31*x43    
      x_sr_q3ex    =x_sr_q3ex    
     9  +coeff( 27)    *x21    *x44    
     1  +coeff( 28)    *x21*x34*x42    
c
      return
      end
      function t_sr_q3ex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff(  7)
      data ncoeff/  6/
      data avdat/ -0.4451506E-03/
      data xmin/
     1 -0.19972E-02,-0.54501E-01,-0.19979E-01,-0.29899E-01,-0.46635E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19884E-02, 0.54444E-01, 0.19967E-01, 0.29220E-01, 0.49077E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.30094722E-02,-0.24234930E-01, 0.39379490E-03, 0.10775150E+00,
     +  0.56202071E-02,-0.97521350E-02,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x21 = x2
      x41 = x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      t_sr_q3ex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)                *x52
c
      return
      end
      function y_sr_q3ex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 96)
      data ncoeff/ 95/
      data avdat/  0.3068465E-02/
      data xmin/
     1 -0.19972E-02,-0.54501E-01,-0.19979E-01,-0.29899E-01,-0.46635E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19884E-02, 0.54444E-01, 0.19967E-01, 0.29220E-01, 0.49077E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.69063420E-03, 0.10333930E-02,-0.41686080E-01,-0.25507312E-01,
     +  0.12113420E-01,-0.56943000E-02,-0.43119490E-02,-0.18809823E-04,
     + -0.53683303E-01, 0.64017571E-03, 0.99747490E-04,-0.64158720E-02,
     +  0.54643102E-02, 0.52159383E-01, 0.63373770E-02, 0.18448153E-02,
     + -0.63189550E-02,-0.65734220E-03,-0.36228653E-01, 0.16297480E-02,
     + -0.63473230E-02,-0.68743624E-02,-0.37967740E-03,-0.78249853E-02,
     + -0.34793743E-02, 0.52977560E-02,-0.91808550E-02, 0.18108040E-01,
     + -0.15401490E-02, 0.15110281E-01,-0.18306214E-01,-0.65617280E-03,
     +  0.73894471E-02, 0.85552763E-02,-0.12795930E-01,-0.30718792E-01,
     + -0.10095360E-01,-0.55171260E-02, 0.30812384E-01, 0.31207420E-01,
     +  0.40406980E-02,-0.49808580E-01, 0.42360020E-02, 0.26277634E-02,
     + -0.30022021E-02,-0.16778721E-02,-0.61100730E-03,-0.90541522E-02,
     + -0.38728753E-02,-0.52623100E-01,-0.63004190E-01,-0.60393850E-01,
     +  0.45022390E-02, 0.65243830E-02, 0.55177211E-02, 0.28840710E-02,
     + -0.30182240E-01,-0.36746822E-02,-0.13925153E-01,-0.34913800E-02,
     + -0.33789260E-01, 0.43893270E-01,-0.35301250E-01, 0.11684960E-01,
     + -0.50335124E-01,-0.24416801E-03, 0.85943300E-03,-0.46679770E-03,
     + -0.53734844E-02,-0.42245432E-03,-0.10290161E-02, 0.25331211E-03,
     + -0.19901320E-02, 0.46845851E-03,-0.30972480E-02,-0.11492600E-02,
     +  0.82644550E-03, 0.24916520E-02, 0.37417600E-02, 0.51533612E-02,
     + -0.22607810E-02, 0.31132580E-02,-0.11754130E-01, 0.18186680E-02,
     +  0.12123640E-01,-0.19645111E-02,-0.66623762E-02,-0.16556943E-02,
     + -0.88286632E-03, 0.11754783E-01,-0.92584690E-02,-0.83127210E-02,
     + -0.19074860E-02,-0.47894220E-02, 0.70462810E-02,
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
      y_sr_q3ex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)        *x32        
      y_sr_q3ex    =y_sr_q3ex    
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)        *x31*x41    
     2  +coeff( 11)            *x42    
     3  +coeff( 12)    *x21        *x51
     4  +coeff( 13)        *x31    *x51
     5  +coeff( 14)            *x41*x51
     6  +coeff( 15)                *x52
     7  +coeff( 16)    *x23            
     8  +coeff( 17)    *x22*x31        
      y_sr_q3ex    =y_sr_q3ex    
     9  +coeff( 18)        *x33        
     1  +coeff( 19)    *x22    *x41    
     2  +coeff( 20)    *x21    *x42    
     3  +coeff( 21)            *x43    
     4  +coeff( 22)    *x22        *x51
     5  +coeff( 23)    *x21*x31    *x51
     6  +coeff( 24)    *x21    *x41*x51
     7  +coeff( 25)            *x42*x51
     8  +coeff( 26)    *x23    *x41    
      y_sr_q3ex    =y_sr_q3ex    
     9  +coeff( 27)    *x22*x31*x41    
     1  +coeff( 28)    *x21*x32*x41    
     2  +coeff( 29)    *x22    *x42    
     3  +coeff( 30)    *x21    *x43    
     4  +coeff( 31)    *x22    *x41*x51
     5  +coeff( 32)*x11*x22    *x41    
     6  +coeff( 33)    *x23*x31*x41    
     7  +coeff( 34)    *x23    *x42    
     8  +coeff( 35)    *x22*x31*x42    
      y_sr_q3ex    =y_sr_q3ex    
     9  +coeff( 36)    *x22    *x43    
     1  +coeff( 37)    *x22*x31*x41*x51
     2  +coeff( 38)    *x22    *x42*x51
     3  +coeff( 39)    *x23*x31*x42    
     4  +coeff( 40)    *x23    *x43    
     5  +coeff( 41)    *x23*x34        
     6  +coeff( 42)    *x24    *x44*x51
     7  +coeff( 43)    *x21*x31*x41    
     8  +coeff( 44)    *x24            
      y_sr_q3ex    =y_sr_q3ex    
     9  +coeff( 45)    *x22*x32        
     1  +coeff( 46)    *x22*x31    *x51
     2  +coeff( 47)    *x21    *x42*x51
     3  +coeff( 48)    *x22*x32*x41    
     4  +coeff( 49)    *x22*x32    *x51
     5  +coeff( 50)    *x22*x31*x43    
     6  +coeff( 51)    *x21*x32*x43    
     7  +coeff( 52)    *x22    *x44    
     8  +coeff( 53)    *x21*x33*x41*x51
      y_sr_q3ex    =y_sr_q3ex    
     9  +coeff( 54)    *x23    *x42*x51
     1  +coeff( 55)    *x21    *x42*x53
     2  +coeff( 56)    *x22*x32*x43    
     3  +coeff( 57)    *x22*x31*x44    
     4  +coeff( 58)*x12*x21    *x41*x51
     5  +coeff( 59)    *x22*x31*x43*x51
     6  +coeff( 60)*x12*x21        *x52
     7  +coeff( 61)    *x24*x32*x42    
     8  +coeff( 62)    *x21*x34*x43    
      y_sr_q3ex    =y_sr_q3ex    
     9  +coeff( 63)    *x22*x31*x43*x52
     1  +coeff( 64)*x12*x24    *x41*x51
     2  +coeff( 65)    *x22*x32*x42*x54
     3  +coeff( 66)*x11*x21            
     4  +coeff( 67)*x11        *x41    
     5  +coeff( 68)        *x32*x41    
     6  +coeff( 69)        *x31*x42    
     7  +coeff( 70)            *x41*x52
     8  +coeff( 71)                *x53
      y_sr_q3ex    =y_sr_q3ex    
     9  +coeff( 72)*x11        *x42    
     1  +coeff( 73)        *x31*x43    
     2  +coeff( 74)    *x21*x31    *x52
     3  +coeff( 75)    *x21    *x41*x52
     4  +coeff( 76)            *x41*x53
     5  +coeff( 77)*x12*x21            
     6  +coeff( 78)    *x24    *x41    
     7  +coeff( 79)        *x33*x42    
     8  +coeff( 80)    *x21    *x44    
      y_sr_q3ex    =y_sr_q3ex    
     9  +coeff( 81)        *x31*x43*x51
     1  +coeff( 82)    *x21    *x41*x53
     2  +coeff( 83)    *x21*x34*x41    
     3  +coeff( 84)*x11*x22    *x42    
     4  +coeff( 85)    *x24    *x42    
     5  +coeff( 86)*x11*x21*x31*x42    
     6  +coeff( 87)    *x21*x33*x42    
     7  +coeff( 88)        *x34*x42    
     8  +coeff( 89)*x12*x21        *x51
      y_sr_q3ex    =y_sr_q3ex    
     9  +coeff( 90)    *x24    *x41*x51
     1  +coeff( 91)    *x22*x32*x41*x51
     2  +coeff( 92)    *x22*x31*x42*x51
     3  +coeff( 93)*x11*x21    *x41*x52
     4  +coeff( 94)    *x22*x31    *x53
     5  +coeff( 95)    *x21    *x41*x54
c
      return
      end
      function p_sr_q3ex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 30)
      data ncoeff/ 29/
      data avdat/ -0.2974516E-03/
      data xmin/
     1 -0.19972E-02,-0.54501E-01,-0.19979E-01,-0.29899E-01,-0.46635E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19884E-02, 0.54444E-01, 0.19967E-01, 0.29220E-01, 0.49077E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.35036190E-03,-0.33103200E-03, 0.14955300E-01,-0.46445052E-02,
     + -0.60064680E-02, 0.24376020E-02, 0.41709500E-03, 0.51488482E-03,
     +  0.17519583E-01,-0.75080040E-04, 0.27539850E-02,-0.42016454E-02,
     + -0.21942640E-01,-0.19007534E-02,-0.83464011E-03, 0.17053920E-01,
     +  0.38959292E-03, 0.17704420E-02, 0.22879564E-02, 0.23134970E-02,
     + -0.26260891E-02, 0.12194603E-01, 0.10900420E-01,-0.41460423E-03,
     + -0.80567104E-02,-0.12913210E-01, 0.11728240E-01, 0.15601334E-01,
     +  0.43717500E-02,
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
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
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
c
c                  function
c
      p_sr_q3ex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)        *x32        
      p_sr_q3ex    =p_sr_q3ex    
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)            *x42    
     2  +coeff( 11)    *x21        *x51
     3  +coeff( 12)        *x31    *x51
     4  +coeff( 13)            *x41*x51
     5  +coeff( 14)                *x52
     6  +coeff( 15)    *x23            
     7  +coeff( 16)    *x22    *x41    
     8  +coeff( 17)    *x21    *x42    
      p_sr_q3ex    =p_sr_q3ex    
     9  +coeff( 18)            *x43    
     1  +coeff( 19)    *x22        *x51
     2  +coeff( 20)            *x41*x52
     3  +coeff( 21)    *x23    *x41    
     4  +coeff( 22)    *x22*x31*x41    
     5  +coeff( 23)    *x22    *x42    
     6  +coeff( 24)    *x24*x31        
     7  +coeff( 25)    *x23*x31*x41    
     8  +coeff( 26)    *x23    *x42    
      p_sr_q3ex    =p_sr_q3ex    
     9  +coeff( 27)    *x22    *x43*x51
     1  +coeff( 28)    *x24    *x44    
     2  +coeff( 29)    *x22*x31        
c
      return
      end
      function l_sr_q3ex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 21)
      data ncoeff/ 20/
      data avdat/  0.2133930E+02/
      data xmin/
     1 -0.19972E-02,-0.54501E-01,-0.19979E-01,-0.29899E-01,-0.46635E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19884E-02, 0.54444E-01, 0.19967E-01, 0.29220E-01, 0.49077E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.63321064E-03, 0.32105153E+00, 0.30394660E-01, 0.32575763E-01,
     + -0.65469330E-03,-0.62573430E-01, 0.76433440E-02,-0.44968780E-02,
     + -0.30751670E-02, 0.12338270E-01, 0.78210070E-02,-0.49449924E-01,
     +  0.22396490E-02, 0.35996413E-02,-0.84933360E-02,-0.52896010E-02,
     + -0.75872080E-01,-0.65512600E-01,-0.65620220E-02,-0.26024540E-01,
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
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_sr_q3ex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)    *x21    *x42    
     7  +coeff(  7)    *x23*x31*x41    
     8  +coeff(  8)*x11                
      l_sr_q3ex    =l_sr_q3ex    
     9  +coeff(  9)    *x21*x31        
     1  +coeff( 10)    *x21        *x51
     2  +coeff( 11)                *x52
     3  +coeff( 12)    *x21*x31*x41    
     4  +coeff( 13)            *x41    
     5  +coeff( 14)            *x42    
     6  +coeff( 15)    *x21*x32        
     7  +coeff( 16)    *x21    *x41*x51
     8  +coeff( 17)    *x21*x31*x42    
      l_sr_q3ex    =l_sr_q3ex    
     9  +coeff( 18)    *x21    *x43    
     1  +coeff( 19)    *x23*x32*x41    
     2  +coeff( 20)    *x21*x34*x41    
c
      return
      end
      function x_sr_fp     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.9414225E-02/
      data xmin/
     1 -0.19972E-02,-0.54829E-01,-0.99828E-02,-0.29963E-01,-0.45933E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19993E-02, 0.54184E-01, 0.99907E-02, 0.27775E-01, 0.49990E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.17082154E-01, 0.31706710E-02, 0.68658840E-04, 0.10155600E-02,
     +  0.62519884E+00, 0.54543150E-01,-0.48334920E-01, 0.22326970E-01,
     + -0.16351813E-01,-0.10864420E-01,-0.45405332E-02,-0.38314824E-02,
     + -0.23479554E-02,-0.26231980E-02, 0.39855030E-02,-0.29205900E-02,
     +  0.52753314E-02,-0.17898970E-02,-0.14897742E-01,-0.26909130E-01,
     + -0.15818921E-02, 0.45236372E-02,-0.98071983E-02,-0.84522934E-02,
     + -0.31411331E-02,-0.30925171E-01,-0.34642620E-03, 0.23394774E-02,
     + -0.54579292E-03,-0.58766650E-03,-0.66701460E-03,-0.93222810E-03,
     + -0.46441360E-02,-0.14731180E-02, 0.33061301E-02,-0.24017423E-01,
     + -0.56848980E-02,-0.11956120E-01,-0.22474360E-05, 0.22076530E-01,
     +  0.52531860E-03, 0.77435730E-02,-0.23693722E-03,-0.35973920E-02,
     +  0.88762930E-03,-0.34134523E-02,-0.62853324E-03, 0.26523340E-02,
     + -0.80291640E-02,-0.19219690E-03,-0.92119520E-02,-0.42994124E-02,
     + -0.23060990E-03, 0.13056790E-01, 0.59540863E-02,-0.48790540E-02,
     +  0.18265140E-02, 0.10298480E-01,-0.74360840E-02, 0.18213400E-02,
     +  0.50777030E-03, 0.54054064E-02,-0.17876170E-02,-0.20231870E-03,
     +  0.71311442E-04, 0.65126013E-03, 0.80317780E-03, 0.11720004E-02,
     + -0.73045010E-03,-0.92823060E-04,-0.13380510E-02, 0.90851910E-03,
     + -0.55457110E-02, 0.28780780E-02,-0.57601680E-03,-0.15075670E-03,
     +  0.15163222E-02, 0.19927604E-02, 0.59074494E-02,-0.98433904E-03,
     + -0.18264450E-02, 0.87062020E-03,-0.10765620E-01,-0.31017400E-03,
     +  0.38332590E-02, 0.30737340E-02, 0.13160043E-02, 0.50572930E-02,
     + -0.93627130E-02, 0.65799980E-03,
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
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      x_sr_fp     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)                *x52
     8  +coeff(  8)    *x22            
      x_sr_fp     =x_sr_fp     
     9  +coeff(  9)    *x21    *x42    
     1  +coeff( 10)            *x42    
     2  +coeff( 11)*x11                
     3  +coeff( 12)    *x21*x31*x41    
     4  +coeff( 13)    *x21    *x41*x51
     5  +coeff( 14)    *x21        *x52
     6  +coeff( 15)    *x21    *x41    
     7  +coeff( 16)            *x41*x51
     8  +coeff( 17)                *x53
      x_sr_fp     =x_sr_fp     
     9  +coeff( 18)    *x22    *x42    
     1  +coeff( 19)    *x21*x31*x42    
     2  +coeff( 20)    *x21    *x43    
     3  +coeff( 21)        *x31*x41    
     4  +coeff( 22)    *x22        *x51
     5  +coeff( 23)    *x24            
     6  +coeff( 24)    *x21    *x42*x51
     7  +coeff( 25)    *x23    *x42    
     8  +coeff( 26)    *x21    *x44    
      x_sr_fp     =x_sr_fp     
     9  +coeff( 27)        *x31    *x51
     1  +coeff( 28)    *x23            
     2  +coeff( 29)    *x21*x32        
     3  +coeff( 30)        *x31*x41*x51
     4  +coeff( 31)    *x23    *x41    
     5  +coeff( 32)    *x22*x31*x41    
     6  +coeff( 33)    *x21*x32*x41    
     7  +coeff( 34)    *x21    *x41*x52
     8  +coeff( 35)    *x22    *x43    
      x_sr_fp     =x_sr_fp     
     9  +coeff( 36)    *x21*x31*x43    
     1  +coeff( 37)    *x24        *x51
     2  +coeff( 38)    *x21    *x43*x51
     3  +coeff( 39)    *x21*x31    *x53
     4  +coeff( 40)    *x23    *x43    
     5  +coeff( 41)*x11*x23            
     6  +coeff( 42)    *x23*x31*x42*x51
     7  +coeff( 43)            *x43    
     8  +coeff( 44)    *x23        *x51
      x_sr_fp     =x_sr_fp     
     9  +coeff( 45)    *x22    *x41*x51
     1  +coeff( 46)    *x22        *x52
     2  +coeff( 47)    *x21*x33*x41    
     3  +coeff( 48)    *x22*x31*x42    
     4  +coeff( 49)    *x21*x32*x42    
     5  +coeff( 50)*x11            *x51
     6  +coeff( 51)    *x21*x31*x42*x51
     7  +coeff( 52)    *x23        *x52
     8  +coeff( 53)        *x32*x41*x52
      x_sr_fp     =x_sr_fp     
     9  +coeff( 54)    *x23*x31*x42    
     1  +coeff( 55)    *x22    *x44    
     2  +coeff( 56)    *x21*x31*x44    
     3  +coeff( 57)    *x23*x31*x41*x51
     4  +coeff( 58)    *x23    *x42*x51
     5  +coeff( 59)    *x21*x31*x43*x51
     6  +coeff( 60)    *x24    *x43    
     7  +coeff( 61)    *x23*x32*x41*x51
     8  +coeff( 62)    *x21    *x43*x54
      x_sr_fp     =x_sr_fp     
     9  +coeff( 63)    *x22*x33*x43*x51
     1  +coeff( 64)    *x21*x31        
     2  +coeff( 65)        *x32        
     3  +coeff( 66)            *x41*x52
     4  +coeff( 67)        *x31*x41*x52
     5  +coeff( 68)            *x42*x52
     6  +coeff( 69)                *x54
     7  +coeff( 70)*x11        *x41    
     8  +coeff( 71)    *x23*x31*x41    
      x_sr_fp     =x_sr_fp     
     9  +coeff( 72)    *x22*x32*x41    
     1  +coeff( 73)    *x21*x32*x41*x51
     2  +coeff( 74)    *x22    *x42*x51
     3  +coeff( 75)            *x44*x51
     4  +coeff( 76)*x11*x22            
     5  +coeff( 77)    *x24*x31*x41    
     6  +coeff( 78)    *x23*x32*x41    
     7  +coeff( 79)    *x24    *x42    
     8  +coeff( 80)    *x22*x31*x43    
      x_sr_fp     =x_sr_fp     
     9  +coeff( 81)    *x21*x32*x42*x51
     1  +coeff( 82)    *x22        *x54
     2  +coeff( 83)    *x23    *x44    
     3  +coeff( 84)*x11    *x31*x41*x51
     4  +coeff( 85)    *x21*x34*x41*x51
     5  +coeff( 86)    *x21    *x44*x52
     6  +coeff( 87)    *x23        *x54
     7  +coeff( 88)    *x24*x31*x43    
     8  +coeff( 89)    *x23*x31*x44    
      x_sr_fp     =x_sr_fp     
     9  +coeff( 90)*x11*x21*x31*x41*x51
c
      return
      end
      function t_sr_fp     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 18)
      data ncoeff/ 17/
      data avdat/  0.1201644E-02/
      data xmin/
     1 -0.19972E-02,-0.54829E-01,-0.99828E-02,-0.29963E-01,-0.45933E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19993E-02, 0.54184E-01, 0.99907E-02, 0.27775E-01, 0.49990E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.35918470E-02,-0.23732670E-01, 0.45985943E-04, 0.25194630E-03,
     +  0.10735630E+00, 0.59513971E-02,-0.98453790E-02,-0.45717820E-03,
     +  0.76593214E-03,-0.16449372E-02,-0.91194454E-03, 0.56492350E-03,
     + -0.40243414E-03,-0.11544590E-02,-0.71767780E-03, 0.10187060E-02,
     +  0.16581024E-02,
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
      x53 = x52*x5
c
c                  function
c
      t_sr_fp     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)                *x52
     8  +coeff(  8)*x11                
      t_sr_fp     =t_sr_fp     
     9  +coeff(  9)    *x22            
     1  +coeff( 10)            *x42    
     2  +coeff( 11)    *x21        *x52
     3  +coeff( 12)    *x21    *x41    
     4  +coeff( 13)            *x41*x51
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)    *x21    *x41*x51
     7  +coeff( 16)                *x53
     8  +coeff( 17)    *x22    *x42    
c
      return
      end
      function y_sr_fp     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 96)
      data ncoeff/ 95/
      data avdat/  0.2276455E-02/
      data xmin/
     1 -0.19972E-02,-0.54829E-01,-0.99828E-02,-0.29963E-01,-0.45933E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19993E-02, 0.54184E-01, 0.99907E-02, 0.27775E-01, 0.49990E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.95784681E-03,-0.33802680E-03,-0.40157270E-01,-0.64089990E-02,
     +  0.97572483E-04,-0.11466610E-02, 0.31979414E-02,-0.35303083E-02,
     + -0.89742941E-02, 0.22768660E-05, 0.18785974E-02, 0.17067680E-02,
     +  0.10742174E-02, 0.26171433E-03, 0.15914350E-02, 0.13919570E-05,
     + -0.15519261E-03, 0.29308714E-05,-0.90009394E-04, 0.10099550E-01,
     + -0.29143970E-02, 0.98826773E-02,-0.54039760E-03,-0.12078762E-03,
     +  0.69012492E-03, 0.68239540E-04, 0.88765122E-03, 0.40843570E-02,
     +  0.86962210E-02,-0.25360864E-05,-0.13152883E-02,-0.59783593E-02,
     + -0.38525070E-02, 0.79826562E-03,-0.13632790E-02, 0.58312254E-03,
     + -0.84705870E-02,-0.31221930E-01,-0.37916231E-05, 0.44491280E-03,
     + -0.14268523E-02, 0.24850584E-01, 0.72541020E-01, 0.30094090E-01,
     +  0.28210854E-01,-0.88411560E-03, 0.73762121E-03, 0.20029700E-02,
     +  0.54223270E-03, 0.16750190E-02, 0.25745140E-03, 0.32980470E-03,
     + -0.68720530E-03, 0.71179580E-04,-0.44603740E-03, 0.20383340E-02,
     +  0.50597481E-01,-0.18246073E-01, 0.10454940E-01,-0.20445514E-01,
     + -0.86400471E-01,-0.27747321E-03, 0.18702753E-03,-0.72165880E-03,
     +  0.11036790E-03,-0.46096794E-03,-0.40666990E-03,-0.42840890E-03,
     + -0.76572000E-03,-0.11639300E-02,-0.18298634E-02, 0.46527520E-02,
     +  0.29046750E-02, 0.92531510E-02,-0.30948991E-04,-0.84364524E-03,
     +  0.11562842E-01,-0.13706500E-01,-0.36157890E-02,-0.18688340E-02,
     + -0.21617960E-01, 0.60362811E-02,-0.54608490E-02,-0.45875290E-02,
     +  0.12458612E-01, 0.52659963E-02,-0.88400654E-02, 0.20802980E-01,
     + -0.31342470E-02, 0.10645600E-01,-0.13245920E-01,-0.11726690E-01,
     +  0.64544421E-02,-0.16593330E-02, 0.67299650E-02,
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
      x45 = x44*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      y_sr_fp     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)        *x32        
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)        *x31    *x51
      y_sr_fp     =y_sr_fp     
     9  +coeff(  9)            *x41*x51
     1  +coeff( 10)*x11                
     2  +coeff( 11)    *x22            
     3  +coeff( 12)    *x21        *x51
     4  +coeff( 13)        *x32*x41    
     5  +coeff( 14)        *x31*x42    
     6  +coeff( 15)            *x43    
     7  +coeff( 16)    *x21*x32        
     8  +coeff( 17)    *x21    *x42    
      y_sr_fp     =y_sr_fp     
     9  +coeff( 18)            *x42*x51
     1  +coeff( 19)*x11        *x41    
     2  +coeff( 20)    *x22    *x41    
     3  +coeff( 21)    *x21    *x41*x51
     4  +coeff( 22)            *x41*x52
     5  +coeff( 23)    *x23            
     6  +coeff( 24)        *x33*x41    
     7  +coeff( 25)            *x44    
     8  +coeff( 26)    *x21*x33        
      y_sr_fp     =y_sr_fp     
     9  +coeff( 27)    *x21    *x43    
     1  +coeff( 28)    *x22*x31*x41    
     2  +coeff( 29)    *x22    *x42    
     3  +coeff( 30)        *x32    *x52
     4  +coeff( 31)    *x23    *x41    
     5  +coeff( 32)    *x22    *x41*x51
     6  +coeff( 33)            *x41*x53
     7  +coeff( 34)    *x24            
     8  +coeff( 35)            *x45    
      y_sr_fp     =y_sr_fp     
     9  +coeff( 36)    *x22*x33        
     1  +coeff( 37)    *x22*x31*x42    
     2  +coeff( 38)    *x22    *x43    
     3  +coeff( 39)        *x33    *x52
     4  +coeff( 40)        *x31    *x54
     5  +coeff( 41)    *x21        *x54
     6  +coeff( 42)    *x22    *x44    
     7  +coeff( 43)    *x22    *x45    
     8  +coeff( 44)    *x22*x32*x44    
      y_sr_fp     =y_sr_fp     
     9  +coeff( 45)    *x22*x31*x45    
     1  +coeff( 46)    *x21*x31        
     2  +coeff( 47)                *x52
     3  +coeff( 48)    *x22*x31        
     4  +coeff( 49)    *x21*x31    *x51
     5  +coeff( 50)        *x31    *x52
     6  +coeff( 51)    *x21*x31*x42    
     7  +coeff( 52)*x11*x21    *x41    
     8  +coeff( 53)        *x31*x44    
      y_sr_fp     =y_sr_fp     
     9  +coeff( 54)    *x21*x32*x42    
     1  +coeff( 55)            *x43*x52
     2  +coeff( 56)    *x22    *x42*x52
     3  +coeff( 57)    *x22*x31*x44    
     4  +coeff( 58)    *x23    *x45    
     5  +coeff( 59)    *x22*x34*x43    
     6  +coeff( 60)    *x22    *x44*x54
     7  +coeff( 61)*x12*x24*x33*x44*x52
     8  +coeff( 62)    *x21            
      y_sr_fp     =y_sr_fp     
     9  +coeff( 63)*x11*x21            
     1  +coeff( 64)    *x22        *x51
     2  +coeff( 65)                *x53
     3  +coeff( 66)    *x23*x31        
     4  +coeff( 67)    *x22*x31    *x51
     5  +coeff( 68)    *x21    *x41*x52
     6  +coeff( 69)        *x31    *x53
     7  +coeff( 70)        *x34*x41    
     8  +coeff( 71)    *x21    *x41*x53
      y_sr_fp     =y_sr_fp     
     9  +coeff( 72)    *x22*x31*x43    
     1  +coeff( 73)    *x23*x31*x42    
     2  +coeff( 74)    *x22    *x44*x51
     3  +coeff( 75)*x11*x21*x31*x42*x51
     4  +coeff( 76)*x12*x22    *x41    
     5  +coeff( 77)    *x23*x32*x43    
     6  +coeff( 78)    *x23*x31*x44    
     7  +coeff( 79)    *x21    *x45*x52
     8  +coeff( 80)            *x45*x53
      y_sr_fp     =y_sr_fp     
     9  +coeff( 81)    *x24    *x44    
     1  +coeff( 82)    *x23    *x44*x51
     2  +coeff( 83)    *x21    *x44*x53
     3  +coeff( 84)*x11*x22*x31*x42*x51
     4  +coeff( 85)    *x23    *x43*x52
     5  +coeff( 86)*x11*x23*x31*x41*x51
     6  +coeff( 87)    *x22*x34*x42*x51
     7  +coeff( 88)    *x24    *x45    
     8  +coeff( 89)*x11*x22    *x44*x51
      y_sr_fp     =y_sr_fp     
     9  +coeff( 90)    *x22*x32*x42*x53
     1  +coeff( 91)    *x22*x34*x44    
     2  +coeff( 92)    *x23*x34*x43    
     3  +coeff( 93)*x11*x24*x31*x41*x52
     4  +coeff( 94)*x12*x24        *x52
     5  +coeff( 95)*x11*x24*x32*x44    
c
      return
      end
      function p_sr_fp     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.1205423E-02/
      data xmin/
     1 -0.19972E-02,-0.54829E-01,-0.99828E-02,-0.29963E-01,-0.45933E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19993E-02, 0.54184E-01, 0.99907E-02, 0.27775E-01, 0.49990E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.12374633E-02, 0.71601140E-02,-0.52271530E-02,-0.58831010E-02,
     + -0.10946390E-04, 0.19251704E-02, 0.51379553E-03,-0.81407630E-04,
     +  0.19652660E-01, 0.81200472E-04, 0.19306174E-03, 0.31013851E-02,
     + -0.23663500E-02,-0.21109400E-01,-0.19770840E-02, 0.19247272E-03,
     + -0.11227640E-02, 0.17909130E-02,-0.31851960E-03, 0.10988142E-01,
     +  0.77473954E-03, 0.10193000E-02, 0.21010304E-02, 0.70989710E-03,
     +  0.33705640E-02, 0.10172662E-03,-0.38080913E-02,-0.32292670E-02,
     + -0.52930023E-02,-0.16251640E-02,-0.18599034E-02,-0.80861700E-03,
     + -0.48463040E-04,-0.21423120E-02, 0.24644981E-02, 0.16324280E-01,
     +  0.23135700E-01, 0.52649430E-02,-0.31246590E-02,-0.45259950E-02,
     +  0.33335670E-01,-0.94169620E-03, 0.16139302E-02, 0.76720920E-04,
     +  0.70846593E-03, 0.39422000E-02,-0.19208950E-02,-0.14255141E-02,
     + -0.52252173E-03,-0.42531434E-02, 0.39810100E-02,-0.28948270E-02,
     +  0.10876842E-01, 0.42360710E-02, 0.18071403E-01,-0.26216961E-01,
     +  0.22318041E-03, 0.18900963E-02,-0.31122500E-03,-0.62389241E-03,
     + -0.10317372E-03, 0.53762984E-02, 0.22289822E-02,-0.39512193E-02,
     + -0.17491530E-02, 0.32354280E-03,-0.16680350E-02,-0.13115560E-02,
     +  0.20009460E-02, 0.71540991E-02,-0.91117020E-03,-0.22412512E-02,
     +  0.13932010E-02,-0.24071434E-01,-0.14067370E-01,-0.10235590E-01,
     +  0.13131853E-01, 0.19364714E-01, 0.32903372E-02, 0.73544770E-02,
     + -0.43624532E-02,-0.23195890E-02, 0.31371161E-01,-0.15571710E-02,
     +  0.18590260E-02, 0.85943294E-02, 0.19535442E-01,-0.17395073E-01,
     + -0.95191262E-02,-0.23855680E-01,
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
      p_sr_fp     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)        *x32        
      p_sr_fp     =p_sr_fp     
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)        *x31*x41    
     2  +coeff( 11)            *x42    
     3  +coeff( 12)    *x21        *x51
     4  +coeff( 13)        *x31    *x51
     5  +coeff( 14)            *x41*x51
     6  +coeff( 15)                *x52
     7  +coeff( 16)*x11*x21            
     8  +coeff( 17)    *x23            
      p_sr_fp     =p_sr_fp     
     9  +coeff( 18)    *x22*x31        
     1  +coeff( 19)*x11        *x41    
     2  +coeff( 20)    *x22    *x41    
     3  +coeff( 21)    *x21    *x42    
     4  +coeff( 22)            *x43    
     5  +coeff( 23)    *x22        *x51
     6  +coeff( 24)            *x42*x51
     7  +coeff( 25)            *x41*x52
     8  +coeff( 26)*x11*x21    *x41    
      p_sr_fp     =p_sr_fp     
     9  +coeff( 27)    *x23    *x41    
     1  +coeff( 28)    *x22    *x42    
     2  +coeff( 29)    *x21    *x43    
     3  +coeff( 30)            *x44    
     4  +coeff( 31)    *x22    *x41*x51
     5  +coeff( 32)    *x21    *x42*x51
     6  +coeff( 33)            *x43*x51
     7  +coeff( 34)    *x23*x31*x41    
     8  +coeff( 35)    *x23    *x42    
      p_sr_fp     =p_sr_fp     
     9  +coeff( 36)    *x22*x31*x42    
     1  +coeff( 37)    *x22    *x43    
     2  +coeff( 38)    *x22    *x42*x51
     3  +coeff( 39)    *x24*x31*x41    
     4  +coeff( 40)    *x23    *x43    
     5  +coeff( 41)    *x22    *x44    
     6  +coeff( 42)    *x21            
     7  +coeff( 43)    *x21    *x41*x51
     8  +coeff( 44)        *x31*x41*x51
      p_sr_fp     =p_sr_fp     
     9  +coeff( 45)        *x31    *x52
     1  +coeff( 46)    *x22*x31*x41    
     2  +coeff( 47)    *x21*x32*x41    
     3  +coeff( 48)    *x21*x31*x42    
     4  +coeff( 49)        *x31*x43    
     5  +coeff( 50)    *x21    *x44    
     6  +coeff( 51)    *x22*x31*x41*x51
     7  +coeff( 52)    *x23*x31*x42    
     8  +coeff( 53)    *x22*x31*x43    
      p_sr_fp     =p_sr_fp     
     9  +coeff( 54)    *x21*x32*x43    
     1  +coeff( 55)    *x22    *x43*x51
     2  +coeff( 56)    *x24    *x43    
     3  +coeff( 57)    *x21*x31*x41    
     4  +coeff( 58)    *x22*x32        
     5  +coeff( 59)        *x33*x41    
     6  +coeff( 60)    *x23        *x51
     7  +coeff( 61)    *x21    *x41*x52
     8  +coeff( 62)    *x24    *x41    
      p_sr_fp     =p_sr_fp     
     9  +coeff( 63)    *x22*x32*x41    
     1  +coeff( 64)    *x21*x31*x43    
     2  +coeff( 65)    *x21    *x43*x51
     3  +coeff( 66)        *x32    *x53
     4  +coeff( 67)    *x24*x32        
     5  +coeff( 68)    *x24    *x42    
     6  +coeff( 69)*x11*x21*x31*x41*x51
     7  +coeff( 70)    *x22*x31*x42*x51
     8  +coeff( 71)    *x22*x31*x41*x52
      p_sr_fp     =p_sr_fp     
     9  +coeff( 72)            *x43*x53
     1  +coeff( 73)    *x23*x33*x41    
     2  +coeff( 74)    *x24*x31*x42    
     3  +coeff( 75)    *x23    *x44    
     4  +coeff( 76)    *x22*x31*x44    
     5  +coeff( 77)    *x24*x31*x43    
     6  +coeff( 78)    *x24    *x44    
     7  +coeff( 79)*x11*x22*x31*x41*x52
     8  +coeff( 80)    *x23    *x43*x52
      p_sr_fp     =p_sr_fp     
     9  +coeff( 81)    *x21*x32*x42*x53
     1  +coeff( 82)*x12*x24    *x41    
     2  +coeff( 83)    *x24*x31*x44    
     3  +coeff( 84)*x12*x22*x31    *x52
     4  +coeff( 85)*x11*x23*x31*x41*x52
     5  +coeff( 86)    *x21*x34*x44*x51
     6  +coeff( 87)    *x22*x31*x43*x54
     7  +coeff( 88)*x11*x23*x31*x43*x52
     8  +coeff( 89)    *x24*x33*x41*x53
      p_sr_fp     =p_sr_fp     
     9  +coeff( 90)    *x23*x32*x44*x53
c
      return
      end
      function l_sr_fp     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 19)
      data ncoeff/ 18/
      data avdat/ -0.1519356E-01/
      data xmin/
     1 -0.19972E-02,-0.54829E-01,-0.99828E-02,-0.29963E-01,-0.45933E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19993E-02, 0.54184E-01, 0.99907E-02, 0.27775E-01, 0.49990E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.15668330E-01,-0.31905782E+00,-0.33070933E-01, 0.44130400E-02,
     + -0.35053590E-01,-0.33357230E-02,-0.24291200E-01, 0.31800480E-01,
     + -0.33933222E-02,-0.90786663E-04, 0.21768214E-01, 0.60267850E-01,
     + -0.65980880E-02,-0.43015424E-02, 0.54756603E-02, 0.38182243E-01,
     +  0.36959763E-01, 0.17033780E-01,
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
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      l_sr_fp     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)                *x52
     8  +coeff(  8)    *x21    *x42    
      l_sr_fp     =l_sr_fp     
     9  +coeff(  9)            *x41    
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)    *x21*x31*x41    
     3  +coeff( 12)    *x21    *x43    
     4  +coeff( 13)            *x42    
     5  +coeff( 14)    *x21        *x51
     6  +coeff( 15)                *x53
     7  +coeff( 16)    *x21*x31*x42    
     8  +coeff( 17)    *x21    *x44    
      l_sr_fp     =l_sr_fp     
     9  +coeff( 18)    *x22    *x44    
c
      return
      end

