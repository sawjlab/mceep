      function E_TXFIT       (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(10)
      dimension coeff(  3)
      data ncoeff/  2/
      data avdat/ -0.3687619E-02/
      data xmin/
     1 -0.57780E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.50125E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.31657170E-02, 0.89058414E-01,
     +      0.      /
      data ientry/0/
c
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
c
c                  function
c
      E_TXFIT       =avdat
     1  +coeff(  1)    
     2  +coeff(  2)*x11
c
      return
      end
      function e_delta       (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(10)
      dimension coeff( 76)
      data ncoeff/ 75/
      data avdat/ -0.1882823E-02/
      data xmin/
     1 -0.69010E+00,-0.28144E-01,-0.49814E-01,-0.45070E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.54634E+00, 0.24553E-01, 0.49814E-01, 0.45070E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.39884862E-02, 0.50622913E-01,-0.12377950E-02, 0.40097270E-02,
     +  0.39700102E-02,-0.11728850E-02, 0.10129562E-02, 0.12215863E-02,
     +  0.10103060E-03, 0.38528920E-03, 0.25860510E-03, 0.27198652E-03,
     + -0.12176270E-03,-0.56424102E-03, 0.86635900E-03,-0.97910120E-03,
     +  0.71599050E-03,-0.11752060E-02, 0.59314200E-03,-0.86932464E-04,
     + -0.35982830E-03,-0.83722552E-04,-0.86343750E-03,-0.23961050E-04,
     +  0.53753033E-05,-0.17723300E-03,-0.76363590E-03,-0.75737080E-03,
     + -0.91090620E-03,-0.43663560E-04,-0.71769364E-03, 0.34288471E-03,
     +  0.20210870E-04, 0.24764690E-03, 0.16725961E-03,-0.68604580E-03,
     +  0.36006620E-03,-0.79142430E-03, 0.29020872E-03,-0.37886980E-04,
     +  0.93915630E-04, 0.13965122E-03, 0.74772513E-04,-0.23570310E-03,
     + -0.60953380E-03, 0.23991301E-05,-0.47437360E-03, 0.62906343E-04,
     + -0.72657584E-03, 0.38419023E-03,-0.65336270E-03,-0.23461560E-03,
     + -0.11416872E-02,-0.30595992E-03, 0.22321100E-03,-0.33125851E-03,
     +  0.17155780E-03,-0.19372800E-03,-0.23453880E-03, 0.28205500E-04,
     + -0.50900413E-04, 0.30483310E-03, 0.58151203E-04, 0.33781474E-04,
     + -0.62264920E-04, 0.15605430E-03, 0.13921063E-03,-0.17766840E-03,
     + -0.20115260E-03, 0.14353590E-04, 0.19674450E-03, 0.28641452E-03,
     + -0.25983544E-03, 0.20277923E-03,-0.36867454E-03,
     +      0.      /
      data ientry/0/
c
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
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
c
c                  function
c
      e_delta       =avdat
     1  +coeff(  1)                
     2  +coeff(  2)*x11            
     3  +coeff(  3)    *x21        
     4  +coeff(  4)*x12            
     5  +coeff(  5)*x11*x21        
     6  +coeff(  6)    *x22        
     7  +coeff(  7)        *x32    
     8  +coeff(  8)        *x31*x41
      e_delta       =e_delta       
     9  +coeff(  9)            *x42
     1  +coeff( 10)*x13            
     2  +coeff( 11)*x12*x21        
     3  +coeff( 12)*x11*x22        
     4  +coeff( 13)*x11    *x32    
     5  +coeff( 14)    *x21*x32    
     6  +coeff( 15)*x11    *x31*x41
     7  +coeff( 16)    *x21*x31*x41
     8  +coeff( 17)*x11        *x42
      e_delta       =e_delta       
     9  +coeff( 18)    *x21    *x42
     1  +coeff( 19)    *x24        
     2  +coeff( 20)*x12    *x31*x41
     3  +coeff( 21)*x11*x21*x31*x41
     4  +coeff( 22)*x11*x21    *x42
     5  +coeff( 23)    *x22    *x42
     6  +coeff( 24)    *x23        
     7  +coeff( 25)*x13*x21        
     8  +coeff( 26)*x11*x23        
      e_delta       =e_delta       
     9  +coeff( 27)    *x22*x32    
     1  +coeff( 28)    *x22*x31*x41
     2  +coeff( 29)*x11*x22*x31*x41
     3  +coeff( 30)*x12*x21    *x42
     4  +coeff( 31)*x11*x22    *x42
     5  +coeff( 32)    *x23    *x42
     6  +coeff( 33)*x14*x22        
     7  +coeff( 34)*x11*x21*x32    
     8  +coeff( 35)*x11*x24        
      e_delta       =e_delta       
     9  +coeff( 36)    *x23*x32    
     1  +coeff( 37)*x12*x21*x31*x41
     2  +coeff( 38)    *x21*x33*x41
     3  +coeff( 39)    *x24    *x42
     4  +coeff( 40)    *x23*x34    
     5  +coeff( 41)*x14            
     6  +coeff( 42)*x12*x22        
     7  +coeff( 43)        *x32*x42
     8  +coeff( 44)*x12*x23        
      e_delta       =e_delta       
     9  +coeff( 45)    *x21*x34    
     1  +coeff( 46)*x13    *x31*x41
     2  +coeff( 47)    *x21*x32*x42
     3  +coeff( 48)*x11        *x44
     4  +coeff( 49)*x11*x23*x31*x41
     5  +coeff( 50)    *x24*x31*x41
     6  +coeff( 51)*x11*x21*x33*x41
     7  +coeff( 52)*x11*x23    *x42
     8  +coeff( 53)*x11*x21*x32*x42
      e_delta       =e_delta       
     9  +coeff( 54)*x12    *x31*x43
     1  +coeff( 55)*x13*x24        
     2  +coeff( 56)*x14*x21*x32    
     3  +coeff( 57)*x13*x22*x32    
     4  +coeff( 58)*x13    *x31*x43
     5  +coeff( 59)*x14    *x33*x41
     6  +coeff( 60)*x12        *x42
     7  +coeff( 61)*x13    *x32    
     8  +coeff( 62)*x11*x22*x32    
      e_delta       =e_delta       
     9  +coeff( 63)    *x23*x31*x41
     1  +coeff( 64)*x13        *x42
     2  +coeff( 65)*x11    *x32*x42
     3  +coeff( 66)*x11*x23*x32    
     4  +coeff( 67)*x11*x21*x34    
     5  +coeff( 68)    *x22*x34    
     6  +coeff( 69)*x14    *x31*x41
     7  +coeff( 70)*x13*x21*x31*x41
     8  +coeff( 71)*x12*x22*x31*x41
      e_delta       =e_delta       
     9  +coeff( 72)    *x22*x32*x42
     1  +coeff( 73)*x12*x23    *x42
     2  +coeff( 74)*x11*x24    *x42
     3  +coeff( 75)*x11*x22*x32*x42
c
      return
      end
      function e_theta       (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(10)
      dimension coeff( 76)
      data ncoeff/ 75/
      data avdat/  0.2427104E-02/
      data xmin/
     1 -0.69010E+00,-0.28144E-01,-0.49814E-01,-0.45070E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.54634E+00, 0.24553E-01, 0.49814E-01, 0.45070E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.28862934E-02, 0.26194870E-03,-0.60644420E-01, 0.80566080E-02,
     + -0.43439270E-02, 0.42614890E-03, 0.38092610E-03, 0.74502330E-03,
     + -0.60916460E-03,-0.14864461E-02, 0.19036090E-02,-0.70173683E-03,
     + -0.67744503E-03, 0.47386530E-03,-0.14959230E-02,-0.47369170E-04,
     +  0.48181480E-04, 0.27729824E-03,-0.86152820E-03, 0.82670190E-03,
     + -0.38999461E-03, 0.75465090E-03,-0.12472630E-02,-0.17691211E-03,
     + -0.83196560E-03, 0.13777012E-03,-0.59001100E-03,-0.97239244E-04,
     + -0.66026080E-03, 0.75838404E-04, 0.93582460E-03,-0.14742370E-02,
     + -0.58090510E-03,-0.11790260E-03, 0.12113151E-03,-0.18503771E-03,
     +  0.14573951E-03, 0.11403491E-03, 0.52083670E-03, 0.35762650E-03,
     +  0.96164990E-04,-0.14210820E-03, 0.16858090E-04,-0.34694172E-03,
     +  0.19379613E-03,-0.62047560E-03, 0.68812492E-04, 0.29304810E-03,
     + -0.30095193E-02,-0.39137870E-04,-0.25663571E-02,-0.24104882E-02,
     + -0.24849441E-02, 0.90858864E-03,-0.16672342E-02,-0.32029822E-02,
     +  0.10667280E-02,-0.35314880E-02,-0.13071924E-03,-0.16768360E-03,
     + -0.13131740E-03, 0.10484060E-03, 0.20680910E-04,-0.17647470E-03,
     + -0.44763020E-03,-0.65231014E-03, 0.55418262E-03,-0.22222290E-02,
     +  0.10022510E-02,-0.56763580E-03, 0.56457793E-03,-0.54451090E-03,
     +  0.96533420E-03,-0.12792642E-02, 0.75331090E-03,
     +      0.      /
      data ientry/0/
c
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
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
c
c                  function
c
      e_theta       =avdat
     1  +coeff(  1)                
     2  +coeff(  2)*x11            
     3  +coeff(  3)    *x21        
     4  +coeff(  4)*x11*x21        
     5  +coeff(  5)    *x22        
     6  +coeff(  6)        *x32    
     7  +coeff(  7)        *x31*x41
     8  +coeff(  8)            *x42
      e_theta       =e_theta       
     9  +coeff(  9)*x13            
     1  +coeff( 10)*x12*x21        
     2  +coeff( 11)*x11*x22        
     3  +coeff( 12)    *x23        
     4  +coeff( 13)    *x21*x32    
     5  +coeff( 14)*x11    *x31*x41
     6  +coeff( 15)    *x21*x31*x41
     7  +coeff( 16)*x11        *x42
     8  +coeff( 17)*x14            
      e_theta       =e_theta       
     9  +coeff( 18)*x13*x21        
     1  +coeff( 19)*x12*x22        
     2  +coeff( 20)*x11*x23        
     3  +coeff( 21)*x12    *x32    
     4  +coeff( 22)*x11*x21*x32    
     5  +coeff( 23)    *x22*x32    
     6  +coeff( 24)*x12    *x31*x41
     7  +coeff( 25)*x11*x21*x31*x41
     8  +coeff( 26)*x12        *x42
      e_theta       =e_theta       
     9  +coeff( 27)*x11*x21    *x42
     1  +coeff( 28)*x13*x22        
     2  +coeff( 29)*x12*x21*x32    
     3  +coeff( 30)*x13    *x31*x41
     4  +coeff( 31)*x12*x21*x31*x41
     5  +coeff( 32)*x11*x22*x31*x41
     6  +coeff( 33)*x11*x22    *x42
     7  +coeff( 34)*x11        *x44
     8  +coeff( 35)*x14*x22        
      e_theta       =e_theta       
     9  +coeff( 36)*x13*x23        
     1  +coeff( 37)*x12*x24        
     2  +coeff( 38)*x14    *x32    
     3  +coeff( 39)*x12*x22*x32    
     4  +coeff( 40)*x12    *x34    
     5  +coeff( 41)*x14        *x42
     6  +coeff( 42)*x12*x22    *x42
     7  +coeff( 43)    *x23*x34    
     8  +coeff( 44)*x13        *x44
      e_theta       =e_theta       
     9  +coeff( 45)*x12*x21    *x44
     1  +coeff( 46)*x14*x22*x32    
     2  +coeff( 47)*x14*x22    *x42
     3  +coeff( 48)    *x24    *x44
     4  +coeff( 49)*x12            
     5  +coeff( 50)    *x21    *x42
     6  +coeff( 51)    *x23*x32    
     7  +coeff( 52)    *x21*x34    
     8  +coeff( 53)    *x21*x33*x41
      e_theta       =e_theta       
     9  +coeff( 54)*x11*x21*x34    
     1  +coeff( 55)*x11*x23*x31*x41
     2  +coeff( 56)*x11*x21*x32*x42
     3  +coeff( 57)    *x21*x34*x42
     4  +coeff( 58)*x13*x21*x33*x41
     5  +coeff( 59)*x11    *x32    
     6  +coeff( 60)    *x24        
     7  +coeff( 61)    *x22*x31*x41
     8  +coeff( 62)    *x22    *x42
      e_theta       =e_theta       
     9  +coeff( 63)        *x32*x42
     1  +coeff( 64)*x11*x24        
     2  +coeff( 65)    *x23*x31*x41
     3  +coeff( 66)*x13        *x42
     4  +coeff( 67)*x12*x21    *x42
     5  +coeff( 68)    *x21*x32*x42
     6  +coeff( 69)*x11*x23*x32    
     7  +coeff( 70)    *x22*x34    
     8  +coeff( 71)    *x24*x31*x41
      e_theta       =e_theta       
     9  +coeff( 72)*x11*x23    *x42
     1  +coeff( 73)    *x22*x32*x42
     2  +coeff( 74)*x11*x21*x31*x43
     3  +coeff( 75)*x13*x22*x32    
c
      return
      end
      function e_phi         (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(10)
      dimension coeff( 76)
      data ncoeff/ 75/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.69010E+00,-0.28144E-01,-0.49814E-01,-0.45070E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.54634E+00, 0.24553E-01, 0.49814E-01, 0.45070E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.17487643E-01,-0.29655560E-01, 0.98397510E-02, 0.52609070E-02,
     + -0.52526821E-02, 0.72800270E-02,-0.79638370E-03, 0.86906173E-03,
     + -0.33054493E-02, 0.48561770E-02,-0.96981041E-03,-0.36170310E-03,
     + -0.10299570E-02,-0.10739100E-02, 0.11717530E-03, 0.16024122E-02,
     + -0.12350693E-02, 0.37981430E-03,-0.73734030E-03,-0.21286440E-02,
     + -0.27746080E-02, 0.73827170E-03,-0.10559620E-02,-0.50834060E-03,
     + -0.35608332E-03,-0.59472400E-03,-0.80524844E-03, 0.50698581E-03,
     + -0.87376902E-04,-0.16571274E-02, 0.29971050E-02, 0.11621890E-02,
     + -0.33671193E-03, 0.30749131E-03,-0.14976301E-02,-0.12771780E-02,
     + -0.41282970E-03,-0.74392280E-03, 0.12889981E-02,-0.15564570E-03,
     +  0.11562740E-02, 0.20422541E-04,-0.17378691E-03,-0.27197630E-02,
     + -0.18333783E-02,-0.12037600E-02,-0.78850891E-03,-0.17937870E-02,
     +  0.11084770E-02,-0.24426380E-02, 0.15780490E-02,-0.40722161E-03,
     + -0.27267434E-03,-0.54815900E-03,-0.79462080E-03,-0.20813270E-03,
     +  0.82824970E-04,-0.46812290E-02, 0.24091694E-02,-0.14005260E-02,
     + -0.15204910E-02, 0.43287270E-03,-0.82265920E-03, 0.18666980E-02,
     + -0.20100160E-02, 0.10665253E-02, 0.42962073E-02, 0.37922610E-03,
     + -0.94091950E-03, 0.68945984E-03,-0.11669840E-02, 0.36130550E-02,
     +  0.37542170E-03, 0.67739800E-03,-0.36846863E-03,
     +      0.      /
      data ientry/0/
c
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
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
c
c                  function
c
      e_phi         =avdat
     1  +coeff(  1)        *x31    
     2  +coeff(  2)            *x41
     3  +coeff(  3)*x11    *x31    
     4  +coeff(  4)    *x21*x31    
     5  +coeff(  5)*x11        *x41
     6  +coeff(  6)    *x21    *x41
     7  +coeff(  7)*x12    *x31    
     8  +coeff(  8)*x11*x21*x31    
      e_phi         =e_phi         
     9  +coeff(  9)    *x22*x31    
     1  +coeff( 10)*x12        *x41
     2  +coeff( 11)*x11*x21    *x41
     3  +coeff( 12)    *x22    *x41
     4  +coeff( 13)        *x32*x41
     5  +coeff( 14)        *x31*x42
     6  +coeff( 15)            *x43
     7  +coeff( 16)*x11*x22*x31    
     8  +coeff( 17)    *x23*x31    
      e_phi         =e_phi         
     9  +coeff( 18)*x11    *x33    
     1  +coeff( 19)    *x21*x33    
     2  +coeff( 20)*x13        *x41
     3  +coeff( 21)*x11*x22    *x41
     4  +coeff( 22)    *x23    *x41
     5  +coeff( 23)*x11    *x32*x41
     6  +coeff( 24)*x11    *x31*x42
     7  +coeff( 25)    *x21*x31*x42
     8  +coeff( 26)*x11        *x43
      e_phi         =e_phi         
     9  +coeff( 27)*x14    *x31    
     1  +coeff( 28)*x13*x21*x31    
     2  +coeff( 29)*x12    *x33    
     3  +coeff( 30)    *x22*x33    
     4  +coeff( 31)*x14        *x41
     5  +coeff( 32)*x12*x22    *x41
     6  +coeff( 33)*x11*x23    *x41
     7  +coeff( 34)    *x24    *x41
     8  +coeff( 35)    *x22*x32*x41
      e_phi         =e_phi         
     9  +coeff( 36)*x12    *x31*x42
     1  +coeff( 37)*x11*x21*x31*x42
     2  +coeff( 38)        *x33*x42
     3  +coeff( 39)*x12        *x43
     4  +coeff( 40)*x11*x21    *x43
     5  +coeff( 41)    *x22    *x43
     6  +coeff( 42)        *x31*x44
     7  +coeff( 43)*x12*x21*x33    
     8  +coeff( 44)*x14*x21    *x41
      e_phi         =e_phi         
     9  +coeff( 45)*x12*x23    *x41
     1  +coeff( 46)*x11*x24    *x41
     2  +coeff( 47)*x12*x21*x32*x41
     3  +coeff( 48)*x11*x22*x32*x41
     4  +coeff( 49)*x12*x21*x31*x42
     5  +coeff( 50)*x11*x22*x31*x42
     6  +coeff( 51)    *x23*x31*x42
     7  +coeff( 52)*x13        *x43
     8  +coeff( 53)*x12*x21    *x43
      e_phi         =e_phi         
     9  +coeff( 54)    *x23    *x43
     1  +coeff( 55)*x12*x24*x31    
     2  +coeff( 56)*x14    *x33    
     3  +coeff( 57)*x14*x22    *x41
     4  +coeff( 58)*x13*x23    *x41
     5  +coeff( 59)*x12*x24    *x41
     6  +coeff( 60)*x13*x21*x32*x41
     7  +coeff( 61)*x12*x22*x31*x42
     8  +coeff( 62)*x11*x23*x31*x42
      e_phi         =e_phi         
     9  +coeff( 63)    *x24*x31*x42
     1  +coeff( 64)*x12    *x33*x42
     2  +coeff( 65)*x13*x21    *x43
     3  +coeff( 66)*x14*x21*x33    
     4  +coeff( 67)*x14*x23    *x41
     5  +coeff( 68)*x14*x21*x32*x41
     6  +coeff( 69)        *x33    
     7  +coeff( 70)*x13    *x31    
     8  +coeff( 71)*x12*x21*x31    
      e_phi         =e_phi         
     9  +coeff( 72)*x12*x21    *x41
     1  +coeff( 73)    *x21*x32*x41
     2  +coeff( 74)*x11*x23*x31    
     3  +coeff( 75)    *x24*x31    
c
      return
      end
      function e_y00         (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(10)
      dimension coeff( 76)
      data ncoeff/ 75/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.69010E+00,-0.28144E-01,-0.49814E-01,-0.45070E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.54634E+00, 0.24553E-01, 0.49814E-01, 0.45070E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.52208714E-01, 0.39929950E-01,-0.22094262E-01,-0.29987074E-01,
     + -0.31714051E-01,-0.16258660E-01, 0.13512740E-02, 0.27650652E-03,
     +  0.98976800E-02, 0.30010233E-02,-0.11765970E-01,-0.14765980E-01,
     +  0.17623210E-01, 0.39708290E-02, 0.19894510E-02, 0.12103980E-02,
     + -0.87466993E-03,-0.46694830E-02, 0.33423493E-02,-0.14643111E-02,
     +  0.11619210E-02, 0.83493990E-03, 0.21909354E-02, 0.58576660E-02,
     +  0.26590302E-02,-0.24663520E-03, 0.25359560E-02,-0.68681060E-03,
     +  0.14518880E-02, 0.38831902E-03,-0.12668310E-02, 0.14313760E-02,
     +  0.52571980E-02, 0.42008124E-02,-0.38940230E-02, 0.19341380E-02,
     + -0.10813130E-02,-0.21085320E-02, 0.33043194E-02, 0.48266010E-02,
     +  0.53931500E-03, 0.20872962E-02, 0.95260280E-03, 0.21073660E-02,
     + -0.19217660E-02, 0.14807960E-02,-0.12154962E-02,-0.14662470E-02,
     + -0.22175804E-03, 0.15139044E-02, 0.25094180E-02,-0.24052790E-02,
     +  0.62094214E-02,-0.19165310E-02, 0.19272270E-03, 0.64600270E-02,
     +  0.54539710E-03,-0.17769890E-03,-0.38263690E-02, 0.32160030E-02,
     + -0.23764302E-02, 0.10331730E-02, 0.27449650E-02,-0.26372913E-02,
     +  0.24685973E-02, 0.13625590E-02,-0.94646681E-03, 0.14996681E-02,
     +  0.87773670E-03, 0.84332610E-03, 0.11091290E-02, 0.13481281E-02,
     +  0.96075114E-03, 0.35147203E-02, 0.36337173E-02,
     +      0.      /
      data ientry/0/
c
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
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
c
c                  function
c
      e_y00         =avdat
     1  +coeff(  1)        *x31    
     2  +coeff(  2)            *x41
     3  +coeff(  3)*x11    *x31    
     4  +coeff(  4)    *x21*x31    
     5  +coeff(  5)*x11        *x41
     6  +coeff(  6)    *x21    *x41
     7  +coeff(  7)*x12    *x31    
     8  +coeff(  8)*x11*x21*x31    
      e_y00         =e_y00         
     9  +coeff(  9)    *x22*x31    
     1  +coeff( 10)        *x33    
     2  +coeff( 11)*x12        *x41
     3  +coeff( 12)*x11*x21    *x41
     4  +coeff( 13)    *x22    *x41
     5  +coeff( 14)        *x32*x41
     6  +coeff( 15)        *x31*x42
     7  +coeff( 16)            *x43
     8  +coeff( 17)*x13    *x31    
      e_y00         =e_y00         
     9  +coeff( 18)*x11*x22*x31    
     1  +coeff( 19)    *x23*x31    
     2  +coeff( 20)*x11    *x33    
     3  +coeff( 21)*x13        *x41
     4  +coeff( 22)*x12*x21    *x41
     5  +coeff( 23)*x11*x22    *x41
     6  +coeff( 24)    *x23    *x41
     7  +coeff( 25)*x11    *x32*x41
     8  +coeff( 26)    *x21*x32*x41
      e_y00         =e_y00         
     9  +coeff( 27)*x11    *x31*x42
     1  +coeff( 28)    *x21*x31*x42
     2  +coeff( 29)*x11        *x43
     3  +coeff( 30)    *x21    *x43
     4  +coeff( 31)*x14    *x31    
     5  +coeff( 32)*x13*x21*x31    
     6  +coeff( 33)    *x22*x33    
     7  +coeff( 34)*x14        *x41
     8  +coeff( 35)*x13*x21    *x41
      e_y00         =e_y00         
     9  +coeff( 36)*x11*x23    *x41
     1  +coeff( 37)    *x24    *x41
     2  +coeff( 38)*x12    *x32*x41
     3  +coeff( 39)*x11*x21*x32*x41
     4  +coeff( 40)    *x22*x32*x41
     5  +coeff( 41)        *x34*x41
     6  +coeff( 42)*x11*x21*x31*x42
     7  +coeff( 43)    *x22*x31*x42
     8  +coeff( 44)*x12        *x43
      e_y00         =e_y00         
     9  +coeff( 45)*x11*x21    *x43
     1  +coeff( 46)    *x22    *x43
     2  +coeff( 47)*x14*x21*x31    
     3  +coeff( 48)*x12*x21*x33    
     4  +coeff( 49)*x11*x22*x33    
     5  +coeff( 50)    *x23*x33    
     6  +coeff( 51)*x14*x21    *x41
     7  +coeff( 52)*x13*x22    *x41
     8  +coeff( 53)*x11*x22*x32*x41
      e_y00         =e_y00         
     9  +coeff( 54)*x13    *x31*x42
     1  +coeff( 55)*x12*x21*x31*x42
     2  +coeff( 56)*x11*x22*x31*x42
     3  +coeff( 57)*x13        *x43
     4  +coeff( 58)    *x24*x33    
     5  +coeff( 59)*x14*x22    *x41
     6  +coeff( 60)*x13*x23    *x41
     7  +coeff( 61)*x14    *x31*x42
     8  +coeff( 62)*x14*x21*x33    
      e_y00         =e_y00         
     9  +coeff( 63)*x14*x21    *x43
     1  +coeff( 64)    *x23*x31*x44
     2  +coeff( 65)*x12*x21*x31    
     3  +coeff( 66)    *x21*x33    
     4  +coeff( 67)*x11*x23*x31    
     5  +coeff( 68)    *x24*x31    
     6  +coeff( 69)*x12*x22    *x41
     7  +coeff( 70)*x13*x22*x31    
     8  +coeff( 71)*x11*x24    *x41
      e_y00         =e_y00         
     9  +coeff( 72)*x13    *x32*x41
     1  +coeff( 73)*x14*x22*x31    
     2  +coeff( 74)*x13*x21*x31*x42
     3  +coeff( 75)    *x22*x33*x42
c
      return
      end

      function H_TXFIT       (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(10)
      dimension coeff(  3)
      data ncoeff/  2/
      data avdat/ -0.3678100E-02/
      data xmin/
     1 -0.57695E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.50054E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.31561650E-02, 0.88844083E-01,
     +      0.      /
      data ientry/0/
c
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
c
c                  function
c
      H_TXFIT       =avdat
     1  +coeff(  1)    
     2  +coeff(  2)*x11
c
      return
      end
      function h_delta       (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10) !,xmean(10)
      dimension coeff( 76)
      data ncoeff/ 75/
      data avdat/ -0.1228340E-02/
      data xmin/
     1 -0.69147E+00,-0.28492E-01,-0.46587E-01,-0.46772E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.54800E+00, 0.24574E-01, 0.46587E-01, 0.46772E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.46258443E-02, 0.50807181E-01,-0.12846910E-02, 0.40369410E-02,
     +  0.40056360E-02,-0.12421030E-02, 0.86782910E-03, 0.12251740E-02,
     +  0.14077090E-03, 0.37337530E-03, 0.22750960E-03, 0.17283750E-03,
     + -0.13376380E-03,-0.46905182E-03, 0.84913504E-03,-0.98328851E-03,
     +  0.78776110E-03,-0.12475063E-02, 0.60428091E-03,-0.80789950E-04,
     + -0.23781370E-03,-0.11106410E-03,-0.97891400E-03,-0.62463461E-03,
     + -0.14782891E-03, 0.36526480E-03,-0.48124711E-04, 0.98618984E-05,
     + -0.17794610E-03,-0.93285690E-03, 0.37554050E-03, 0.22622190E-04,
     +  0.21213780E-03,-0.54916052E-03,-0.74548860E-03, 0.30468354E-03,
     + -0.66437502E-03, 0.33808051E-03,-0.70316570E-03,-0.84880800E-03,
     +  0.35417923E-03, 0.94932242E-04, 0.20207500E-04, 0.81606860E-04,
     + -0.22095160E-03, 0.37831874E-03,-0.47565830E-03,-0.53897692E-03,
     +  0.11964103E-03, 0.24858150E-03,-0.56795001E-03,-0.21454030E-03,
     + -0.11575910E-02, 0.26509733E-03,-0.31120210E-03,-0.43097380E-04,
     +  0.78313040E-04,-0.27944281E-03,-0.14344362E-04,-0.36188070E-03,
     + -0.55068294E-03,-0.17464880E-03, 0.12205010E-03, 0.33266130E-04,
     +  0.30094470E-04, 0.15961151E-03, 0.90409020E-04,-0.15084030E-04,
     +  0.76407400E-04,-0.74134620E-04, 0.16962680E-03,-0.21226540E-03,
     + -0.24475650E-03,-0.10953430E-03,-0.19985710E-03,
     +      0.      /
      data ientry/0/
c
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
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
c
c                  function
c
      h_delta       =avdat
     1  +coeff(  1)                
     2  +coeff(  2)*x11            
     3  +coeff(  3)    *x21        
     4  +coeff(  4)*x12            
     5  +coeff(  5)*x11*x21        
     6  +coeff(  6)    *x22        
     7  +coeff(  7)        *x32    
     8  +coeff(  8)        *x31*x41
      h_delta       =h_delta       
     9  +coeff(  9)            *x42
     1  +coeff( 10)*x13            
     2  +coeff( 11)*x12*x21        
     3  +coeff( 12)*x11*x22        
     4  +coeff( 13)*x11    *x32    
     5  +coeff( 14)    *x21*x32    
     6  +coeff( 15)*x11    *x31*x41
     7  +coeff( 16)    *x21*x31*x41
     8  +coeff( 17)*x11        *x42
      h_delta       =h_delta       
     9  +coeff( 18)    *x21    *x42
     1  +coeff( 19)    *x24        
     2  +coeff( 20)*x12    *x31*x41
     3  +coeff( 21)*x11*x21*x31*x41
     4  +coeff( 22)*x11*x21    *x42
     5  +coeff( 23)    *x22    *x42
     6  +coeff( 24)*x11*x22    *x42
     7  +coeff( 25)    *x24*x32    
     8  +coeff( 26)    *x24*x31*x41
      h_delta       =h_delta       
     9  +coeff( 27)    *x23        
     1  +coeff( 28)*x13*x21        
     2  +coeff( 29)*x11*x23        
     3  +coeff( 30)*x11*x22*x31*x41
     4  +coeff( 31)    *x23    *x42
     5  +coeff( 32)*x14*x22        
     6  +coeff( 33)*x11*x21*x32    
     7  +coeff( 34)    *x22*x32    
     8  +coeff( 35)    *x22*x31*x41
      h_delta       =h_delta       
     9  +coeff( 36)*x11*x24        
     1  +coeff( 37)    *x23*x32    
     2  +coeff( 38)*x12*x21*x31*x41
     3  +coeff( 39)    *x21*x33*x41
     4  +coeff( 40)*x11*x23*x31*x41
     5  +coeff( 41)    *x24    *x42
     6  +coeff( 42)*x14            
     7  +coeff( 43)        *x34    
     8  +coeff( 44)        *x32*x42
      h_delta       =h_delta       
     9  +coeff( 45)*x12*x23        
     1  +coeff( 46)*x11*x22*x32    
     2  +coeff( 47)    *x21*x34    
     3  +coeff( 48)    *x21*x32*x42
     4  +coeff( 49)*x11*x21*x34    
     5  +coeff( 50)*x12*x22*x31*x41
     6  +coeff( 51)*x11*x21*x33*x41
     7  +coeff( 52)*x11*x23    *x42
     8  +coeff( 53)*x11*x21*x32*x42
      h_delta       =h_delta       
     9  +coeff( 54)    *x22*x32*x42
     1  +coeff( 55)*x12    *x31*x43
     2  +coeff( 56)*x14*x23        
     3  +coeff( 57)*x13*x24        
     4  +coeff( 58)*x14*x21*x32    
     5  +coeff( 59)*x14*x21    *x42
     6  +coeff( 60)*x12*x23    *x42
     7  +coeff( 61)*x11*x22*x32*x42
     8  +coeff( 62)*x13    *x31*x43
      h_delta       =h_delta       
     9  +coeff( 63)*x12*x22        
     1  +coeff( 64)*x12        *x42
     2  +coeff( 65)*x14*x21        
     3  +coeff( 66)*x13*x22        
     4  +coeff( 67)    *x23*x31*x41
     5  +coeff( 68)    *x21*x31*x43
     6  +coeff( 69)*x11        *x44
     7  +coeff( 70)    *x21    *x44
     8  +coeff( 71)*x11*x23*x32    
      h_delta       =h_delta       
     9  +coeff( 72)    *x22*x34    
     1  +coeff( 73)*x14    *x31*x41
     2  +coeff( 74)*x12    *x33*x41
     3  +coeff( 75)*x11*x21*x31*x43
c
      return
      end
      function h_theta       (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10) !,xmean(10)
      dimension coeff( 76)
      data ncoeff/ 75/
      data avdat/  0.1193335E-02/
      data xmin/
     1 -0.69147E+00,-0.28492E-01,-0.46587E-01,-0.46772E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.54800E+00, 0.24574E-01, 0.46587E-01, 0.46772E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.45066480E-02, 0.22988630E-03,-0.60951720E-01, 0.81514883E-02,
     + -0.45520350E-02, 0.79862970E-03,-0.63367830E-03,-0.15306910E-02,
     +  0.17250140E-02,-0.72975320E-03,-0.18164201E-03,-0.71681081E-03,
     +  0.56638660E-03,-0.14038120E-02,-0.14527091E-03, 0.11265811E-03,
     +  0.17903013E-03,-0.51817060E-03, 0.79000890E-03, 0.49360880E-03,
     +  0.13234202E-03,-0.47477840E-03,-0.10149332E-02,-0.46561200E-03,
     +  0.18961150E-03, 0.29080983E-03,-0.77220692E-03, 0.29509340E-03,
     + -0.18731923E-02,-0.18991880E-03,-0.32257750E-04,-0.24887950E-03,
     + -0.37355780E-03,-0.90501731E-03,-0.13129210E-03,-0.78961884E-05,
     +  0.23979100E-04, 0.10775400E-03,-0.12143940E-03,-0.29547800E-03,
     + -0.11386001E-02,-0.36465103E-03, 0.12305751E-02, 0.11242863E-02,
     + -0.50777551E-03,-0.58665940E-03, 0.37900810E-03, 0.28553280E-05,
     +  0.21332412E-04,-0.30898381E-02, 0.26767744E-03, 0.40629721E-03,
     + -0.16054612E-03,-0.13282700E-03,-0.23336160E-02,-0.21909090E-02,
     +  0.24197254E-03,-0.22763340E-02,-0.28624750E-02,-0.63965143E-03,
     + -0.67279880E-02, 0.15513112E-03, 0.71081270E-03,-0.15202853E-02,
     + -0.32722990E-03,-0.60412710E-03,-0.40950200E-03,-0.17266310E-02,
     + -0.31656763E-03, 0.12039412E-02, 0.83241240E-03,-0.97353450E-03,
     +  0.98424253E-03, 0.91998590E-03,-0.54006883E-03,
     +      0.      /
      data ientry/0/
c
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
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
c
c                  function
c
      h_theta       =avdat
     1  +coeff(  1)                
     2  +coeff(  2)*x11            
     3  +coeff(  3)    *x21        
     4  +coeff(  4)*x11*x21        
     5  +coeff(  5)    *x22        
     6  +coeff(  6)            *x42
     7  +coeff(  7)*x13            
     8  +coeff(  8)*x12*x21        
      h_theta       =h_theta       
     9  +coeff(  9)*x11*x22        
     1  +coeff( 10)    *x23        
     2  +coeff( 11)*x11    *x32    
     3  +coeff( 12)    *x21*x32    
     4  +coeff( 13)*x11    *x31*x41
     5  +coeff( 14)    *x21*x31*x41
     6  +coeff( 15)*x11        *x42
     7  +coeff( 16)*x14            
     8  +coeff( 17)*x13*x21        
      h_theta       =h_theta       
     9  +coeff( 18)*x12*x22        
     1  +coeff( 19)*x11*x23        
     2  +coeff( 20)*x11*x21*x32    
     3  +coeff( 21)        *x34    
     4  +coeff( 22)*x12    *x31*x41
     5  +coeff( 23)*x11*x21*x31*x41
     6  +coeff( 24)    *x22*x31*x41
     7  +coeff( 25)        *x33*x41
     8  +coeff( 26)*x12        *x42
      h_theta       =h_theta       
     9  +coeff( 27)*x11*x21    *x42
     1  +coeff( 28)    *x22    *x42
     2  +coeff( 29)*x11*x22*x31*x41
     3  +coeff( 30)*x14*x22        
     4  +coeff( 31)*x13*x23        
     5  +coeff( 32)*x12*x24        
     6  +coeff( 33)*x12*x22*x32    
     7  +coeff( 34)    *x24*x32    
     8  +coeff( 35)*x12    *x34    
      h_theta       =h_theta       
     9  +coeff( 36)*x14        *x42
     1  +coeff( 37)*x12*x22    *x42
     2  +coeff( 38)*x12    *x31*x43
     3  +coeff( 39)*x12        *x44
     4  +coeff( 40)*x13*x24        
     5  +coeff( 41)*x12*x21*x34    
     6  +coeff( 42)    *x23*x34    
     7  +coeff( 43)*x14*x21*x31*x41
     8  +coeff( 44)*x14*x21    *x42
      h_theta       =h_theta       
     9  +coeff( 45)*x13*x22    *x42
     1  +coeff( 46)*x13        *x44
     2  +coeff( 47)*x14*x24        
     3  +coeff( 48)*x14    *x34    
     4  +coeff( 49)*x14*x22    *x42
     5  +coeff( 50)*x12            
     6  +coeff( 51)        *x32    
     7  +coeff( 52)        *x31*x41
     8  +coeff( 53)    *x21    *x42
      h_theta       =h_theta       
     9  +coeff( 54)    *x22*x32    
     1  +coeff( 55)    *x23*x32    
     2  +coeff( 56)    *x21*x33*x41
     3  +coeff( 57)    *x21    *x44
     4  +coeff( 58)*x11*x23*x31*x41
     5  +coeff( 59)*x11*x21*x32*x42
     6  +coeff( 60)    *x21*x34*x42
     7  +coeff( 61)*x11*x21*x33*x43
     8  +coeff( 62)*x13*x22        
      h_theta       =h_theta       
     9  +coeff( 63)*x11*x22*x32    
     1  +coeff( 64)    *x21*x34    
     2  +coeff( 65)    *x23*x31*x41
     3  +coeff( 66)*x13        *x42
     4  +coeff( 67)*x11*x22    *x42
     5  +coeff( 68)    *x21*x32*x42
     6  +coeff( 69)    *x21*x31*x43
     7  +coeff( 70)*x11*x23*x32    
     8  +coeff( 71)*x11*x21*x34    
      h_theta       =h_theta       
     9  +coeff( 72)    *x22*x34    
     1  +coeff( 73)*x12*x22*x31*x41
     2  +coeff( 74)    *x24*x31*x41
     3  +coeff( 75)*x11*x23    *x42
c
      return
      end
      function h_phi         (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10) !,xmean(10)
      dimension coeff( 76)
      data ncoeff/ 75/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.69147E+00,-0.28492E-01,-0.46587E-01,-0.46772E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.54800E+00, 0.24574E-01, 0.46587E-01, 0.46772E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.15343380E-01,-0.30512182E-01, 0.89406110E-02, 0.49360180E-02,
     + -0.47948230E-02, 0.77780922E-02,-0.58445532E-03, 0.55769220E-03,
     + -0.29297540E-02,-0.78356250E-03, 0.47266980E-02, 0.83657410E-03,
     + -0.69641374E-03,-0.12901810E-02,-0.11329964E-02, 0.10245060E-03,
     +  0.11816144E-02,-0.82826990E-03, 0.37434850E-03,-0.40858294E-03,
     + -0.24401640E-02,-0.31085824E-02, 0.42762720E-03,-0.93458121E-03,
     +  0.64999470E-04,-0.82698190E-03, 0.35769710E-03,-0.78134510E-03,
     + -0.10198723E-02, 0.10038820E-02,-0.70115481E-03,-0.14061520E-02,
     +  0.33316253E-02, 0.18622510E-02,-0.22353683E-02, 0.35464980E-03,
     + -0.15361232E-02,-0.84373634E-03, 0.29098891E-03, 0.12482391E-02,
     + -0.13949680E-02, 0.13464730E-02,-0.28540092E-03,-0.24373050E-03,
     +  0.37102690E-03, 0.63984390E-03,-0.23051390E-03, 0.38942770E-03,
     + -0.83995633E-03, 0.14929691E-02,-0.22627792E-02,-0.17435314E-03,
     + -0.53170220E-04,-0.40543410E-03, 0.12805800E-03,-0.45948251E-03,
     +  0.15532350E-02, 0.15211881E-02, 0.89804470E-03,-0.16312670E-02,
     + -0.29222660E-03, 0.53282680E-03, 0.64355940E-04,-0.13356822E-02,
     +  0.85930370E-03,-0.12118790E-02, 0.23000722E-02,-0.31193881E-03,
     +  0.82092842E-03,-0.56001740E-03,-0.34947851E-02, 0.71968240E-03,
     + -0.16322181E-02,-0.59686170E-03,-0.64741950E-03,
     +      0.      /
      data ientry/0/
c
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
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
c
c                  function
c
      h_phi         =avdat
     1  +coeff(  1)        *x31    
     2  +coeff(  2)            *x41
     3  +coeff(  3)*x11    *x31    
     4  +coeff(  4)    *x21*x31    
     5  +coeff(  5)*x11        *x41
     6  +coeff(  6)    *x21    *x41
     7  +coeff(  7)*x12    *x31    
     8  +coeff(  8)*x11*x21*x31    
      h_phi         =h_phi         
     9  +coeff(  9)    *x22*x31    
     1  +coeff( 10)        *x33    
     2  +coeff( 11)*x12        *x41
     3  +coeff( 12)*x11*x21    *x41
     4  +coeff( 13)    *x22    *x41
     5  +coeff( 14)        *x32*x41
     6  +coeff( 15)        *x31*x42
     7  +coeff( 16)            *x43
     8  +coeff( 17)*x11*x22*x31    
      h_phi         =h_phi         
     9  +coeff( 18)    *x23*x31    
     1  +coeff( 19)*x11    *x33    
     2  +coeff( 20)    *x21*x33    
     3  +coeff( 21)*x13        *x41
     4  +coeff( 22)*x11*x22    *x41
     5  +coeff( 23)    *x23    *x41
     6  +coeff( 24)*x11    *x32*x41
     7  +coeff( 25)    *x21*x32*x41
     8  +coeff( 26)*x11    *x31*x42
      h_phi         =h_phi         
     9  +coeff( 27)    *x21*x31*x42
     1  +coeff( 28)*x11        *x43
     2  +coeff( 29)*x14    *x31    
     3  +coeff( 30)*x13*x21*x31    
     4  +coeff( 31)*x12*x22*x31    
     5  +coeff( 32)    *x22*x33    
     6  +coeff( 33)*x14        *x41
     7  +coeff( 34)*x12*x22    *x41
     8  +coeff( 35)*x11*x23    *x41
      h_phi         =h_phi         
     9  +coeff( 36)    *x24    *x41
     1  +coeff( 37)*x11*x21*x32*x41
     2  +coeff( 38)    *x22*x32*x41
     3  +coeff( 39)*x11*x21*x31*x42
     4  +coeff( 40)*x12        *x43
     5  +coeff( 41)*x11*x21    *x43
     6  +coeff( 42)    *x22    *x43
     7  +coeff( 43)*x12*x23*x31    
     8  +coeff( 44)*x13    *x33    
      h_phi         =h_phi         
     9  +coeff( 45)*x12*x21*x33    
     1  +coeff( 46)*x11*x22*x33    
     2  +coeff( 47)*x14*x21    *x41
     3  +coeff( 48)*x12*x23    *x41
     4  +coeff( 49)*x11*x24    *x41
     5  +coeff( 50)*x12*x21*x32*x41
     6  +coeff( 51)*x11*x22*x32*x41
     7  +coeff( 52)*x13        *x43
     8  +coeff( 53)    *x23    *x43
      h_phi         =h_phi         
     9  +coeff( 54)*x13*x23*x31    
     1  +coeff( 55)*x14*x22    *x41
     2  +coeff( 56)*x13*x23    *x41
     3  +coeff( 57)*x12*x24    *x41
     4  +coeff( 58)*x12*x22*x32*x41
     5  +coeff( 59)*x14    *x31*x42
     6  +coeff( 60)*x12*x22*x31*x42
     7  +coeff( 61)*x13*x21    *x43
     8  +coeff( 62)*x14*x23*x31    
      h_phi         =h_phi         
     9  +coeff( 63)*x14*x23    *x41
     1  +coeff( 64)*x14*x21*x32*x41
     2  +coeff( 65)*x13    *x31    
     3  +coeff( 66)*x12*x21*x31    
     4  +coeff( 67)*x12*x21    *x41
     5  +coeff( 68)    *x21    *x43
     6  +coeff( 69)*x11*x23*x31    
     7  +coeff( 70)    *x24*x31    
     8  +coeff( 71)*x13*x21    *x41
      h_phi         =h_phi         
     9  +coeff( 72)*x12    *x32*x41
     1  +coeff( 73)*x12    *x31*x42
     2  +coeff( 74)    *x22*x31*x42
     3  +coeff( 75)    *x23*x33    
c
      return
      end
      function h_y00         (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10) !,xmean(10)
      dimension coeff( 76)
      data ncoeff/ 75/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.69147E+00,-0.28492E-01,-0.46587E-01,-0.46772E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.54800E+00, 0.24574E-01, 0.46587E-01, 0.46772E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.51782100E-01, 0.38850843E-01,-0.20134504E-01,-0.28835970E-01,
     + -0.34928161E-01,-0.17336292E-01, 0.70985534E-03, 0.90171560E-02,
     +  0.22820564E-02,-0.11796310E-01,-0.16142560E-01, 0.19431663E-01,
     +  0.38986850E-02, 0.19339800E-02, 0.14537390E-02,-0.75537500E-03,
     + -0.34106330E-02, 0.38783010E-02,-0.10458381E-02, 0.14375991E-02,
     +  0.10832980E-02, 0.75974650E-03, 0.14030810E-02, 0.65724470E-02,
     +  0.24776251E-02,-0.20288620E-03, 0.28215070E-02,-0.77266973E-04,
     +  0.17721510E-02,-0.12169730E-02, 0.17234080E-02, 0.37366122E-03,
     +  0.48526981E-02, 0.45457342E-02,-0.39764782E-02, 0.18787270E-02,
     + -0.12888820E-02,-0.20477810E-02, 0.19553801E-02, 0.33741070E-02,
     +  0.32260360E-02, 0.16294690E-02, 0.26284314E-02,-0.21260850E-02,
     +  0.19361570E-02,-0.12734160E-02, 0.31915740E-02,-0.19167120E-02,
     +  0.60803150E-02,-0.24044171E-02, 0.15472534E-02, 0.58077932E-02,
     + -0.28211620E-02, 0.13518774E-03, 0.22070760E-02, 0.15971570E-02,
     + -0.49968380E-03,-0.28876160E-02, 0.26809500E-02,-0.67314380E-03,
     +  0.55363350E-02, 0.93509024E-03, 0.20875960E-02,-0.15306700E-02,
     +  0.15597782E-02,-0.78398460E-03, 0.10049162E-02, 0.11932812E-02,
     + -0.13280360E-02,-0.16059570E-02, 0.17910380E-02, 0.34426774E-02,
     +  0.25292013E-02,-0.35794410E-02,-0.35851093E-03,
     +      0.      /
      data ientry/0/
c
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
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
c
c                  function
c
      h_y00         =avdat
     1  +coeff(  1)        *x31    
     2  +coeff(  2)            *x41
     3  +coeff(  3)*x11    *x31    
     4  +coeff(  4)    *x21*x31    
     5  +coeff(  5)*x11        *x41
     6  +coeff(  6)    *x21    *x41
     7  +coeff(  7)*x11*x21*x31    
     8  +coeff(  8)    *x22*x31    
      h_y00         =h_y00         
     9  +coeff(  9)        *x33    
     1  +coeff( 10)*x12        *x41
     2  +coeff( 11)*x11*x21    *x41
     3  +coeff( 12)    *x22    *x41
     4  +coeff( 13)        *x32*x41
     5  +coeff( 14)        *x31*x42
     6  +coeff( 15)            *x43
     7  +coeff( 16)*x13    *x31    
     8  +coeff( 17)*x11*x22*x31    
      h_y00         =h_y00         
     9  +coeff( 18)    *x23*x31    
     1  +coeff( 19)*x11    *x33    
     2  +coeff( 20)    *x21*x33    
     3  +coeff( 21)*x13        *x41
     4  +coeff( 22)*x12*x21    *x41
     5  +coeff( 23)*x11*x22    *x41
     6  +coeff( 24)    *x23    *x41
     7  +coeff( 25)*x11    *x32*x41
     8  +coeff( 26)    *x21*x32*x41
      h_y00         =h_y00         
     9  +coeff( 27)*x11    *x31*x42
     1  +coeff( 28)    *x21*x31*x42
     2  +coeff( 29)*x11        *x43
     3  +coeff( 30)*x14    *x31    
     4  +coeff( 31)*x13*x21*x31    
     5  +coeff( 32)*x12    *x33    
     6  +coeff( 33)    *x22*x33    
     7  +coeff( 34)*x14        *x41
     8  +coeff( 35)*x13*x21    *x41
      h_y00         =h_y00         
     9  +coeff( 36)*x11*x23    *x41
     1  +coeff( 37)    *x24    *x41
     2  +coeff( 38)*x12    *x32*x41
     3  +coeff( 39)*x11*x21*x32*x41
     4  +coeff( 40)    *x22*x32*x41
     5  +coeff( 41)*x11*x21*x31*x42
     6  +coeff( 42)    *x22*x31*x42
     7  +coeff( 43)*x12        *x43
     8  +coeff( 44)*x11*x21    *x43
      h_y00         =h_y00         
     9  +coeff( 45)    *x22    *x43
     1  +coeff( 46)*x14*x21*x31    
     2  +coeff( 47)*x14*x21    *x41
     3  +coeff( 48)*x13*x22    *x41
     4  +coeff( 49)*x11*x22*x32*x41
     5  +coeff( 50)*x13    *x31*x42
     6  +coeff( 51)*x12*x21*x31*x42
     7  +coeff( 52)*x11*x22*x31*x42
     8  +coeff( 53)    *x23*x31*x42
      h_y00         =h_y00         
     9  +coeff( 54)*x13        *x43
     1  +coeff( 55)*x12*x21    *x43
     2  +coeff( 56)*x12*x24*x31    
     3  +coeff( 57)    *x24*x33    
     4  +coeff( 58)*x14*x22    *x41
     5  +coeff( 59)*x13*x23    *x41
     6  +coeff( 60)*x14*x21*x31*x42
     7  +coeff( 61)*x13    *x32*x43
     8  +coeff( 62)*x12    *x31    
      h_y00         =h_y00         
     9  +coeff( 63)*x12*x21*x31    
     1  +coeff( 64)*x11*x23*x31    
     2  +coeff( 65)    *x24*x31    
     3  +coeff( 66)*x12    *x31*x42
     4  +coeff( 67)        *x33*x42
     5  +coeff( 68)*x13*x22*x31    
     6  +coeff( 69)*x11*x24*x31    
     7  +coeff( 70)*x11*x22*x33    
     8  +coeff( 71)*x11*x24    *x41
      h_y00         =h_y00         
     9  +coeff( 72)*x11*x23*x32*x41
     1  +coeff( 73)    *x22*x34*x41
     2  +coeff( 74)*x12*x23*x32*x41
     3  +coeff( 75)*x11*x21*x33    
c
      return
      end



