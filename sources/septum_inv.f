      function sl_txfit    (x,m)
c for left septum at 6 degrees -JJL 3/4/04      
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)
      dimension coeff(  4)
      data ncoeff/  3/
      data avdat/ -0.3438139E-02/
      data xmin/
     1 -0.56268E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49365E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.24954679E-02, 0.90844929E-01,-0.95373677E-03,
     +      0.      /
      data ientry/0/
c
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
c
c                  function
c
      sl_txfit    =avdat
     1  +coeff(  1)    
     2  +coeff(  2)*x11
     3  +coeff(  3)*x12
c
      return
      end
      function sl_delta       (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !,xmean(10)
      dimension coeff( 96)
      data ncoeff/ 95/
      data avdat/  0.1776592E-02/
      data xmin/
     1 -0.68722E+00,-0.26339E-01,-0.41533E-01,-0.29172E-01,-0.19999E-02,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.59420E+00, 0.22418E-01, 0.40099E-01, 0.34626E-01, 0.19995E-02,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.53079491E-02, 0.47928188E-01, 0.38221909E-02, 0.40922058E-02,
     + -0.15057774E-02, 0.37049922E-04, 0.34697444E-03, 0.12938809E-02,
     + -0.93647925E-03, 0.44119582E-03,-0.13279649E-02, 0.12072878E-02,
     + -0.55917411E-03, 0.51839690E-03,-0.11098744E-03, 0.18215622E-03,
     + -0.30075980E-03,-0.32053499E-02, 0.40827626E-02,-0.17350391E-02,
     +  0.15426383E-03,-0.16527968E-03, 0.77360884E-04,-0.20082975E-04,
     + -0.33698916E-02, 0.52599953E-02, 0.17474666E-04,-0.75703912E-03,
     +  0.31349834E-03, 0.52567560E-03, 0.12582098E-03,-0.24859927E-03,
     + -0.91767055E-04,-0.23644810E-03, 0.19775158E-03,-0.23120102E-02,
     +  0.72816233E-02, 0.35133347E-03,-0.75614075E-02, 0.29406196E-02,
     +  0.14205406E-03, 0.35109915E-03,-0.49876858E-03, 0.12070362E-03,
     + -0.36610446E-02, 0.10136695E-01,-0.23933626E-02,-0.93397116E-02,
     +  0.33135726E-02,-0.40063034E-04, 0.63596628E-04,-0.10119873E-04,
     +  0.17804791E-03, 0.26759924E-02,-0.54161006E-04, 0.33093198E-04,
     +  0.38303380E-03,-0.35245167E-03,-0.30358115E-03, 0.10749113E-03,
     +  0.79779849E-04,-0.69648144E-04,-0.39128415E-03,-0.33877420E-03,
     +  0.20972753E-03,-0.19119270E-03,-0.52758271E-03,-0.31828207E-04,
     + -0.45046647E-03, 0.54002483E-03,-0.17833937E-02, 0.56911116E-02,
     + -0.72967508E-02,-0.77820456E-04, 0.50186524E-02, 0.20417884E-03,
     + -0.17333570E-02, 0.19769359E-02,-0.28662682E-02, 0.13707198E-02,
     + -0.23479096E-02,-0.69373148E-03, 0.23005425E-03,-0.12629827E-02,
     + -0.38494781E-03,-0.82792109E-03, 0.11513984E-02, 0.38801234E-04,
     +  0.11039606E-03, 0.19072769E-03, 0.14920680E-03, 0.13620964E-04,
     +  0.11175629E-03,-0.15253609E-03, 0.19315346E-03,
     +      0.      /
      data ientry/0/
c
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
      sl_delta       =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)*x11                
     3  +coeff(  3)*x12                
     4  +coeff(  4)*x11*x21            
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21            
     7  +coeff(  7)                *x51
     8  +coeff(  8)        *x32        
      sl_delta       =sl_delta       
     9  +coeff(  9)        *x31*x41    
     1  +coeff( 10)*x11*x22            
     2  +coeff( 11)    *x21*x32        
     3  +coeff( 12)    *x21*x31*x41    
     4  +coeff( 13)    *x21    *x42    
     5  +coeff( 14)    *x24            
     6  +coeff( 15)        *x31        
     7  +coeff( 16)*x13                
     8  +coeff( 17)    *x23            
      sl_delta       =sl_delta       
     9  +coeff( 18)    *x22*x32        
     1  +coeff( 19)    *x22*x31*x41    
     2  +coeff( 20)    *x22    *x42    
     3  +coeff( 21)            *x42    
     4  +coeff( 22)*x11        *x42    
     5  +coeff( 23)    *x21        *x51
     6  +coeff( 24)*x11*x21*x32        
     7  +coeff( 25)*x11*x21*x33        
     8  +coeff( 26)*x11*x21*x32*x41    
      sl_delta       =sl_delta       
     9  +coeff( 27)    *x25*x31        
     1  +coeff( 28)*x12*x22*x33        
     2  +coeff( 29)*x12*x21            
     3  +coeff( 30)*x11    *x31*x41    
     4  +coeff( 31)        *x31*x42    
     5  +coeff( 32)            *x43    
     6  +coeff( 33)*x11*x23            
     7  +coeff( 34)*x12    *x31*x41    
     8  +coeff( 35)*x12        *x42    
      sl_delta       =sl_delta       
     9  +coeff( 36)    *x21*x33        
     1  +coeff( 37)    *x21*x32*x41    
     2  +coeff( 38)*x11    *x31*x42    
     3  +coeff( 39)    *x21*x31*x42    
     4  +coeff( 40)    *x21    *x43    
     5  +coeff( 41)    *x25            
     6  +coeff( 42)    *x24*x31        
     7  +coeff( 43)    *x23*x32        
     8  +coeff( 44)    *x23    *x42    
      sl_delta       =sl_delta       
     9  +coeff( 45)    *x22*x33        
     1  +coeff( 46)    *x22*x32*x41    
     2  +coeff( 47)*x11*x21*x31*x42    
     3  +coeff( 48)    *x22*x31*x42    
     4  +coeff( 49)    *x22    *x43    
     5  +coeff( 50)*x15        *x41    
     6  +coeff( 51)*x13*x22    *x41    
     7  +coeff( 52)*x11*x24    *x41    
     8  +coeff( 53)    *x25    *x41    
      sl_delta       =sl_delta       
     9  +coeff( 54)*x11*x22*x32*x41    
     1  +coeff( 55)*x12*x25    *x41    
     2  +coeff( 56)*x11    *x31        
     3  +coeff( 57)    *x21*x31        
     4  +coeff( 58)    *x21    *x41    
     5  +coeff( 59)    *x22    *x41    
     6  +coeff( 60)*x12*x22            
     7  +coeff( 61)*x13    *x31        
     8  +coeff( 62)*x11*x22*x31        
      sl_delta       =sl_delta       
     9  +coeff( 63)*x11        *x43    
     1  +coeff( 64)*x12*x23            
     2  +coeff( 65)*x13*x21    *x41    
     3  +coeff( 66)    *x24    *x41    
     4  +coeff( 67)*x11*x22*x32        
     5  +coeff( 68)*x12    *x33        
     6  +coeff( 69)*x12    *x32*x41    
     7  +coeff( 70)*x12    *x31*x42    
     8  +coeff( 71)    *x21*x34        
      sl_delta       =sl_delta       
     9  +coeff( 72)    *x21*x33*x41    
     1  +coeff( 73)    *x21*x32*x42    
     2  +coeff( 74)*x11    *x31*x43    
     3  +coeff( 75)    *x21*x31*x43    
     4  +coeff( 76)*x11        *x44    
     5  +coeff( 77)    *x21    *x44    
     6  +coeff( 78)    *x24*x32        
     7  +coeff( 79)    *x24*x31*x41    
     8  +coeff( 80)    *x24    *x42    
      sl_delta       =sl_delta       
     9  +coeff( 81)*x11*x22*x33        
     1  +coeff( 82)    *x23    *x43    
     2  +coeff( 83)*x13*x24            
     3  +coeff( 84)*x11*x21*x34        
     4  +coeff( 85)*x15*x21*x31        
     5  +coeff( 86)*x15*x21*x32        
     6  +coeff( 87)*x13*x21*x33*x41    
     7  +coeff( 88)            *x41    
     8  +coeff( 89)*x11*x21*x31        
      sl_delta       =sl_delta       
     9  +coeff( 90)    *x22*x31        
     1  +coeff( 91)*x11*x21    *x42    
     2  +coeff( 92)    *x22        *x51
     3  +coeff( 93)*x11*x24            
     4  +coeff( 94)        *x34        
     5  +coeff( 95)        *x33*x41    
c
      return
      end
      function sl_theta       (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !,xmean(10)
      dimension coeff( 96)
      data ncoeff/ 95/
      data avdat/ -0.5502503E-03/
      data xmin/
     1 -0.68722E+00,-0.26339E-01,-0.41533E-01,-0.29172E-01,-0.19999E-02,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.59420E+00, 0.22418E-01, 0.40099E-01, 0.34626E-01, 0.19995E-02,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.49176104E-02,-0.40898475E-03,-0.54596618E-01, 0.83065737E-04,
     +  0.63232658E-02,-0.37534977E-02, 0.18369852E-02,-0.17104602E-02,
     + -0.77408226E-02, 0.87517742E-02,-0.26371765E-02, 0.78646018E-03,
     +  0.15154028E-02,-0.74137142E-03,-0.15861994E-02,-0.39320311E-03,
     + -0.11011544E-02, 0.24735715E-03, 0.19004903E-03,-0.75044105E-03,
     + -0.56155775E-04,-0.42984351E-02,-0.12600703E-01, 0.68168058E-02,
     +  0.40326048E-01,-0.37213364E-02,-0.46371415E-01, 0.17791064E-01,
     + -0.15890766E-02, 0.80764276E-03,-0.37630512E-02, 0.21840583E-02,
     +  0.10816535E-02, 0.12652122E-02, 0.73290564E-03,-0.81623206E-03,
     +  0.11294925E-02,-0.38537596E-03, 0.29288279E-02, 0.10253469E-02,
     +  0.31413487E-02,-0.15446076E-01,-0.18077360E-01,-0.30911153E-02,
     + -0.67000706E-02, 0.31722236E-01, 0.54699421E-01, 0.45931195E-02,
     + -0.19735735E-01,-0.49200688E-01, 0.13822318E-01,-0.18477779E-03,
     + -0.95377141E-03,-0.14021863E-01, 0.12080886E-01,-0.19385670E-02,
     +  0.12525861E-03, 0.12964310E-03,-0.40775473E-03, 0.37630831E-03,
     +  0.27376888E-02, 0.27862698E-03, 0.98649680E-03,-0.94347861E-05,
     +  0.49210689E-02, 0.14202649E-01,-0.10125232E-01,-0.60127261E-02,
     + -0.60502375E-02,-0.21249659E-02,-0.90626295E-03, 0.11096551E-02,
     +  0.64180400E-02, 0.16499422E-01,-0.94275475E-02, 0.19401832E-01,
     + -0.12823422E-01, 0.23006769E-02,-0.45886138E-02,-0.85710064E-02,
     + -0.38781386E-01, 0.10598791E-02, 0.18248342E-01, 0.20040583E-01,
     + -0.26108805E-03,-0.12355099E-02,-0.40607061E-03, 0.16773220E-02,
     +  0.52158313E-03,-0.28594909E-03, 0.11428231E-02, 0.67412219E-03,
     + -0.43570814E-02,-0.45119223E-03,-0.69878827E-03,
     +      0.      /
      data ientry/0/
c
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
c
c                  function
c
      sl_theta       =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)*x11                
     3  +coeff(  3)    *x21            
     4  +coeff(  4)*x12                
     5  +coeff(  5)*x11*x21            
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)    *x21    *x41    
      sl_theta       =sl_theta       
     9  +coeff(  9)    *x21*x32        
     1  +coeff( 10)    *x21*x31*x41    
     2  +coeff( 11)    *x21    *x42    
     3  +coeff( 12)                *x51
     4  +coeff( 13)*x11*x22            
     5  +coeff( 14)    *x22*x31        
     6  +coeff( 15)        *x31*x41    
     7  +coeff( 16)*x12*x21            
     8  +coeff( 17)    *x23            
      sl_theta       =sl_theta       
     9  +coeff( 18)*x11    *x32        
     1  +coeff( 19)    *x22    *x41    
     2  +coeff( 20)            *x43    
     3  +coeff( 21)    *x24            
     4  +coeff( 22)    *x22*x32        
     5  +coeff( 23)    *x21*x33        
     6  +coeff( 24)    *x22*x31*x41    
     7  +coeff( 25)    *x21*x32*x41    
     8  +coeff( 26)    *x22    *x42    
      sl_theta       =sl_theta       
     9  +coeff( 27)    *x21*x31*x42    
     1  +coeff( 28)    *x21    *x43    
     2  +coeff( 29)*x13*x21*x31        
     3  +coeff( 30)*x12*x21*x32        
     4  +coeff( 31)*x13*x21*x32        
     5  +coeff( 32)*x13*x21*x31*x41    
     6  +coeff( 33)        *x32        
     7  +coeff( 34)            *x42    
     8  +coeff( 35)*x11*x21*x31        
      sl_theta       =sl_theta       
     9  +coeff( 36)*x11*x21    *x41    
     1  +coeff( 37)*x11*x23            
     2  +coeff( 38)    *x23*x31        
     3  +coeff( 39)*x11*x21*x32        
     4  +coeff( 40)    *x23    *x41    
     5  +coeff( 41)    *x24*x31        
     6  +coeff( 42)*x11*x21*x33        
     7  +coeff( 43)    *x22*x33        
     8  +coeff( 44)    *x24    *x41    
      sl_theta       =sl_theta       
     9  +coeff( 45)*x11*x22*x31*x41    
     1  +coeff( 46)*x11*x21*x32*x41    
     2  +coeff( 47)    *x22*x32*x41    
     3  +coeff( 48)*x11*x22    *x42    
     4  +coeff( 49)*x11*x21*x31*x42    
     5  +coeff( 50)    *x22*x31*x42    
     6  +coeff( 51)    *x22    *x43    
     7  +coeff( 52)*x11            *x51
     8  +coeff( 53)*x11*x24*x31        
      sl_theta       =sl_theta       
     9  +coeff( 54)*x11*x22*x33        
     1  +coeff( 55)*x11*x22*x32*x41    
     2  +coeff( 56)*x11*x22    *x43    
     3  +coeff( 57)        *x31        
     4  +coeff( 58)*x13*x21            
     5  +coeff( 59)*x12*x22            
     6  +coeff( 60)*x11    *x33        
     7  +coeff( 61)*x11*x21    *x42    
     8  +coeff( 62)*x12    *x33        
      sl_theta       =sl_theta       
     9  +coeff( 63)*x13*x21    *x41    
     1  +coeff( 64)    *x23    *x42    
     2  +coeff( 65)*x11*x21    *x43    
     3  +coeff( 66)    *x21*x31*x43    
     4  +coeff( 67)    *x21    *x44    
     5  +coeff( 68)*x12*x21*x33        
     6  +coeff( 69)*x11*x21*x34        
     7  +coeff( 70)    *x22*x34        
     8  +coeff( 71)*x13*x22    *x41    
      sl_theta       =sl_theta       
     9  +coeff( 72)*x11*x24    *x41    
     1  +coeff( 73)*x12*x21*x32*x41    
     2  +coeff( 74)    *x23*x31*x42    
     3  +coeff( 75)    *x23    *x43    
     4  +coeff( 76)    *x21*x32*x43    
     5  +coeff( 77)    *x21*x31*x44    
     6  +coeff( 78)    *x23*x34        
     7  +coeff( 79)    *x22*x34*x41    
     8  +coeff( 80)*x11*x21*x33*x42    
      sl_theta       =sl_theta       
     9  +coeff( 81)    *x23*x33*x42    
     1  +coeff( 82)*x13*x21*x31*x43    
     2  +coeff( 83)    *x23*x32*x43    
     3  +coeff( 84)    *x23*x34*x43    
     4  +coeff( 85)*x13                
     5  +coeff( 86)        *x32*x41    
     6  +coeff( 87)*x11        *x42    
     7  +coeff( 88)        *x31*x42    
     8  +coeff( 89)*x12*x21*x31        
      sl_theta       =sl_theta       
     9  +coeff( 90)*x13        *x41    
     1  +coeff( 91)*x11*x22    *x41    
     2  +coeff( 92)*x12    *x31*x41    
     3  +coeff( 93)*x11*x21*x31*x41    
     4  +coeff( 94)        *x33*x41    
     5  +coeff( 95)*x12        *x42    
c
      return
      end
      function sl_phi         (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !,xmean(10)
      dimension coeff( 96)
      data ncoeff/ 95/
      data avdat/  0.6279154E-03/
      data xmin/
     1 -0.68722E+00,-0.26339E-01,-0.41533E-01,-0.29172E-01,-0.19999E-02,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.59420E+00, 0.22418E-01, 0.40099E-01, 0.34626E-01, 0.19995E-02,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.37751935E-03,-0.30565031E-01,-0.35010077E-03, 0.39327005E-02,
     +  0.15987980E-02, 0.10980141E-01,-0.14066441E-01, 0.13766055E-02,
     + -0.15502351E-02,-0.56276709E-03,-0.17377811E-02,-0.55683427E-02,
     +  0.52884179E-02,-0.86573064E-02, 0.13793706E-01,-0.77893566E-02,
     +  0.65621189E-02,-0.17165671E-03,-0.35786277E-02,-0.62514157E-02,
     +  0.47467276E-02,-0.21038316E-02, 0.28192871E-02, 0.20892262E-02,
     + -0.31898103E-02, 0.27282181E-03,-0.16864186E-01, 0.15330528E-01,
     + -0.35314392E-02,-0.14758459E-02, 0.67324645E-03, 0.10026610E-01,
     + -0.63437311E-03,-0.55015524E-03,-0.72801998E-03, 0.16077973E-02,
     +  0.46486760E-03,-0.12743771E-02, 0.80754375E-02,-0.13045261E-01,
     + -0.83057899E-02, 0.17502483E-01, 0.56922067E-01,-0.15276895E-01,
     + -0.41227747E-01, 0.78564165E-02,-0.12322154E-02,-0.30832568E-02,
     +  0.27202875E-02,-0.11564890E-01,-0.16050929E-01,-0.45113885E-02,
     + -0.99088484E-02,-0.10100759E-02, 0.32651246E-01,-0.83525106E-01,
     +  0.89292591E-02,-0.43711169E-02,-0.27860743E-02, 0.55570458E-03,
     +  0.61424443E-03,-0.22896597E-01, 0.23786962E-01,-0.26926745E-01,
     + -0.10492810E-02,-0.98492287E-03, 0.74468151E-03, 0.44902074E-02,
     +  0.41470606E-01,-0.44151329E-01, 0.16771080E-01, 0.17476747E-01,
     + -0.39847560E-01, 0.30400550E-01, 0.56495825E-02,-0.58599338E-02,
     + -0.93076052E-02, 0.17613957E-01, 0.16023602E-01, 0.37684079E-01,
     + -0.11861368E-01,-0.34648690E-01, 0.11485470E-01,-0.63648927E-02,
     +  0.33467922E-02, 0.16541802E-02, 0.13400966E-02, 0.52117426E-02,
     +  0.57216105E-02,-0.15878180E-01, 0.75080700E-01,-0.33949379E-01,
     +  0.66277166E-02, 0.10188895E-01, 0.19169299E-03,
     +      0.      /
      data ientry/0/
c
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
      sl_phi         =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)*x11                
     5  +coeff(  5)            *x42    
     6  +coeff(  6)*x11    *x31        
     7  +coeff(  7)*x11        *x41    
     8  +coeff(  8)    *x21    *x41    
      sl_phi         =sl_phi         
     9  +coeff(  9)    *x22            
     1  +coeff( 10)            *x43    
     2  +coeff( 11)    *x21    *x42    
     3  +coeff( 12)    *x22*x31        
     4  +coeff( 13)    *x22    *x41    
     5  +coeff( 14)    *x22*x32        
     6  +coeff( 15)    *x22*x31*x41    
     7  +coeff( 16)    *x22    *x42    
     8  +coeff( 17)    *x21*x31        
      sl_phi         =sl_phi         
     9  +coeff( 18)*x11*x21            
     1  +coeff( 19)        *x33        
     2  +coeff( 20)        *x31*x42    
     3  +coeff( 21)*x12        *x41    
     4  +coeff( 22)*x11*x21    *x41    
     5  +coeff( 23)*x11*x21*x31*x41    
     6  +coeff( 24)*x11*x21    *x42    
     7  +coeff( 25)    *x23*x31        
     8  +coeff( 26)*x12*x21    *x41    
      sl_phi         =sl_phi         
     9  +coeff( 27)        *x34*x41    
     1  +coeff( 28)    *x23*x31*x41    
     2  +coeff( 29)    *x23    *x42    
     3  +coeff( 30)    *x23*x34        
     4  +coeff( 31)*x12                
     5  +coeff( 32)        *x32*x41    
     6  +coeff( 33)*x11    *x32        
     7  +coeff( 34)*x11    *x31*x41    
     8  +coeff( 35)*x12    *x31        
      sl_phi         =sl_phi         
     9  +coeff( 36)*x11*x21*x31        
     1  +coeff( 37)*x11*x22            
     2  +coeff( 38)        *x34        
     3  +coeff( 39)        *x33*x41    
     4  +coeff( 40)            *x44    
     5  +coeff( 41)*x11    *x33        
     6  +coeff( 42)*x11    *x32*x41    
     7  +coeff( 43)    *x21*x32*x41    
     8  +coeff( 44)*x11    *x31*x42    
      sl_phi         =sl_phi         
     9  +coeff( 45)    *x21*x31*x42    
     1  +coeff( 46)*x11        *x43    
     2  +coeff( 47)*x12*x21*x31        
     3  +coeff( 48)*x11*x22    *x41    
     4  +coeff( 49)    *x21*x33*x41    
     5  +coeff( 50)*x11*x21*x33        
     6  +coeff( 51)    *x22*x33        
     7  +coeff( 52)*x11*x22*x32        
     8  +coeff( 53)    *x23*x32        
      sl_phi         =sl_phi         
     9  +coeff( 54)        *x31*x45    
     1  +coeff( 55)    *x21*x35        
     2  +coeff( 56)    *x21*x34*x41    
     3  +coeff( 57)    *x23    *x43    
     4  +coeff( 58)*x11*x23*x32        
     5  +coeff( 59)*x11*x21    *x45    
     6  +coeff( 60)    *x21            
     7  +coeff( 61)    *x21*x32        
     8  +coeff( 62)        *x32*x42    
      sl_phi         =sl_phi         
     9  +coeff( 63)        *x31*x43    
     1  +coeff( 64)    *x21*x33        
     2  +coeff( 65)*x12    *x31*x41    
     3  +coeff( 66)*x13    *x31        
     4  +coeff( 67)*x11*x23            
     5  +coeff( 68)        *x35        
     6  +coeff( 69)        *x32*x43    
     7  +coeff( 70)        *x31*x44    
     8  +coeff( 71)            *x45    
      sl_phi         =sl_phi         
     9  +coeff( 72)*x11    *x33*x41    
     1  +coeff( 73)*x11    *x32*x42    
     2  +coeff( 74)*x11    *x31*x43    
     3  +coeff( 75)    *x21*x31*x43    
     4  +coeff( 76)    *x21    *x44    
     5  +coeff( 77)*x12    *x33        
     6  +coeff( 78)*x12    *x32*x41    
     7  +coeff( 79)*x11*x21*x32*x41    
     8  +coeff( 80)    *x22*x32*x41    
      sl_phi         =sl_phi         
     9  +coeff( 81)*x12    *x31*x42    
     1  +coeff( 82)    *x22*x31*x42    
     2  +coeff( 83)    *x22    *x43    
     3  +coeff( 84)*x12*x21*x31*x41    
     4  +coeff( 85)*x11*x22    *x42    
     5  +coeff( 86)    *x24*x31        
     6  +coeff( 87)*x14        *x41    
     7  +coeff( 88)*x11    *x35        
     8  +coeff( 89)*x11    *x34*x41    
      sl_phi         =sl_phi         
     9  +coeff( 90)*x11    *x33*x42    
     1  +coeff( 91)    *x21*x33*x42    
     2  +coeff( 92)    *x21*x31*x44    
     3  +coeff( 93)*x11        *x45    
     4  +coeff( 94)    *x21    *x45    
     5  +coeff( 95)    *x21        *x51
c
      return
      end
      function sl_y00         (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)  !,xmean(10)
      dimension coeff( 96)
      data ncoeff/ 95/
      data avdat/  0.7871948E-03/
      data xmin/
     1 -0.68722E+00,-0.26339E-01,-0.41533E-01,-0.29172E-01,-0.19999E-02,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.59420E+00, 0.22418E-01, 0.40099E-01, 0.34626E-01, 0.19995E-02,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.36197451E-02,-0.43429405E-03, 0.42599868E-01,-0.63999929E-02,
     + -0.47417460E-02, 0.12637476E-01,-0.10150895E-01,-0.32078527E-01,
     + -0.28597817E-01, 0.15495180E-01, 0.36906474E-02,-0.13898483E-02,
     +  0.20839891E-02,-0.45813750E-02, 0.40904456E-02,-0.71756123E-02,
     + -0.17879005E-02, 0.13316627E-01, 0.20602406E-02, 0.38622266E-02,
     +  0.14896723E-01,-0.11941673E-01,-0.49436660E-02,-0.44956557E-02,
     + -0.23560312E-03, 0.10461499E-01, 0.36944546E-01, 0.11163168E-02,
     +  0.30154779E-02, 0.15122155E-01,-0.29616108E-01, 0.13872006E-01,
     +  0.70655363E-03,-0.14777046E-03, 0.37637260E-01,-0.10215282E+00,
     +  0.88448875E-01,-0.15717980E-01, 0.71397685E-02,-0.11427483E-02,
     + -0.56440579E-02, 0.22897437E-03, 0.12107489E+00,-0.10549717E+00,
     + -0.17964365E-02, 0.47236201E-02, 0.42094253E-02, 0.20846762E-01,
     + -0.86874276E-01, 0.15804084E-01,-0.19980770E-01, 0.10874480E-01,
     +  0.74946419E-02,-0.30106142E-01, 0.45896679E-01,-0.97526900E-01,
     + -0.63062795E-02, 0.93173265E-01,-0.28991856E-01, 0.10609726E-01,
     +  0.26608977E-03,-0.15865391E-01, 0.14847040E-01,-0.16737455E-02,
     +  0.95137465E-03, 0.24826126E-02, 0.57576050E-03,-0.58405962E-01,
     + -0.91399718E-02,-0.38390402E-02, 0.10725462E-01,-0.27279467E-02,
     +  0.41519529E-02, 0.42921158E-04,-0.28685294E-01, 0.13136438E+00,
     + -0.10071784E+00, 0.16538059E-01, 0.16167384E-01, 0.36863352E-02,
     + -0.51915082E-02,-0.26949069E-02,-0.11381987E-01, 0.45152740E-02,
     +  0.21126922E-01, 0.21188729E-01, 0.25454070E-01,-0.19120662E-01,
     + -0.11028523E-01,-0.30926269E-03, 0.10498004E-02, 0.10069782E+00,
     + -0.20471734E+00, 0.16228864E+00,-0.52150249E-01,
     +      0.      /
      data ientry/0/
c
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
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x35 = x34*x3
      x36 = x35*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x45 = x44*x4
      x51 = x5
c
c                  function
c
      sl_y00         =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)*x11                
     5  +coeff(  5)        *x32        
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)*x11    *x31        
      sl_y00         =sl_y00         
     9  +coeff(  9)    *x21*x31        
     1  +coeff( 10)*x11        *x41    
     2  +coeff( 11)    *x21    *x41    
     3  +coeff( 12)*x12                
     4  +coeff( 13)    *x22            
     5  +coeff( 14)            *x43    
     6  +coeff( 15)*x11    *x32        
     7  +coeff( 16)    *x21*x32        
     8  +coeff( 17)*x11    *x31*x41    
      sl_y00         =sl_y00         
     9  +coeff( 18)    *x21*x31*x41    
     1  +coeff( 19)*x11        *x42    
     2  +coeff( 20)*x12    *x31        
     3  +coeff( 21)    *x22*x31        
     4  +coeff( 22)*x12        *x41    
     5  +coeff( 23)*x11*x21    *x41    
     6  +coeff( 24)    *x22    *x41    
     7  +coeff( 25)    *x23            
     8  +coeff( 26)        *x34        
      sl_y00         =sl_y00         
     9  +coeff( 27)            *x44    
     1  +coeff( 28)*x11    *x33        
     2  +coeff( 29)    *x21    *x43    
     3  +coeff( 30)    *x22*x32        
     4  +coeff( 31)    *x22*x31*x41    
     5  +coeff( 32)    *x22    *x42    
     6  +coeff( 33)    *x23*x31        
     7  +coeff( 34)        *x31*x44    
     8  +coeff( 35)    *x22*x33        
      sl_y00         =sl_y00         
     9  +coeff( 36)    *x22*x32*x41    
     1  +coeff( 37)    *x22*x31*x42    
     2  +coeff( 38)    *x22    *x43    
     3  +coeff( 39)*x13*x22    *x41    
     4  +coeff( 40)*x11*x21            
     5  +coeff( 41)    *x21    *x42    
     6  +coeff( 42)*x11*x21*x31        
     7  +coeff( 43)        *x32*x42    
     8  +coeff( 44)        *x31*x43    
      sl_y00         =sl_y00         
     9  +coeff( 45)    *x21*x31*x42    
     1  +coeff( 46)*x11*x21*x32        
     2  +coeff( 47)*x11*x21    *x42    
     3  +coeff( 48)    *x21*x34        
     4  +coeff( 49)    *x21*x33*x41    
     5  +coeff( 50)    *x23*x32        
     6  +coeff( 51)    *x23*x31*x41    
     7  +coeff( 52)    *x23    *x42    
     8  +coeff( 53)*x11*x21*x31*x43    
      sl_y00         =sl_y00         
     9  +coeff( 54)*x11*x21    *x44    
     1  +coeff( 55)    *x23*x33        
     2  +coeff( 56)    *x23*x32*x41    
     3  +coeff( 57)*x11*x22*x31*x42    
     4  +coeff( 58)    *x23*x31*x42    
     5  +coeff( 59)    *x23    *x43    
     6  +coeff( 60)*x12*x22    *x42    
     7  +coeff( 61)    *x21*x32*x44    
     8  +coeff( 62)*x11*x23*x31*x42    
      sl_y00         =sl_y00         
     9  +coeff( 63)*x13*x22*x32        
     1  +coeff( 64)*x11*x24    *x42    
     2  +coeff( 65)    *x21            
     3  +coeff( 66)        *x33        
     4  +coeff( 67)        *x32*x41    
     5  +coeff( 68)        *x33*x41    
     6  +coeff( 69)*x11    *x32*x41    
     7  +coeff( 70)    *x21*x32*x41    
     8  +coeff( 71)*x11    *x31*x42    
      sl_y00         =sl_y00         
     9  +coeff( 72)*x11*x22*x31        
     1  +coeff( 73)    *x23    *x41    
     2  +coeff( 74)*x11*x23            
     3  +coeff( 75)*x11    *x33*x41    
     4  +coeff( 76)    *x21*x32*x42    
     5  +coeff( 77)    *x21*x31*x43    
     6  +coeff( 78)*x11        *x44    
     7  +coeff( 79)    *x21    *x44    
     8  +coeff( 80)*x11*x21*x33        
      sl_y00         =sl_y00         
     9  +coeff( 81)*x12*x21*x32        
     1  +coeff( 82)*x13    *x31*x41    
     2  +coeff( 83)*x11*x22    *x42    
     3  +coeff( 84)*x12    *x34        
     4  +coeff( 85)*x11*x21*x32*x42    
     5  +coeff( 86)    *x22    *x44    
     6  +coeff( 87)*x11*x22*x33        
     7  +coeff( 88)*x11*x22*x32*x41    
     8  +coeff( 89)*x11*x23    *x42    
      sl_y00         =sl_y00         
     9  +coeff( 90)            *x41*x51
     1  +coeff( 91)*x11    *x36        
     2  +coeff( 92)*x11    *x34*x42    
     3  +coeff( 93)*x11    *x33*x43    
     4  +coeff( 94)*x11    *x32*x44    
     5  +coeff( 95)*x11    *x31*x45    
c
      return
      end
      function sr_txfit       (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(  3)
      data ncoeff/  2/
      data avdat/ -0.3435877E-02/
      data xmin/
     1 -0.56270E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49370E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.29719390E-02, 0.90831652E-01,
     +      0.      /
      data ientry/0/
c
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
      sr_txfit       =avdat
     1  +coeff(  1)    
     2  +coeff(  2)*x11
c
      return
      end
      function sr_delta       (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 96)
      data ncoeff/ 95/
      data avdat/  0.1808467E-02/
      data xmin/
     1 -0.68849E+00,-0.25651E-01,-0.36952E-01,-0.26097E-01,-0.19972E-02,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.59465E+00, 0.21663E-01, 0.39547E-01, 0.33592E-01, 0.19993E-02,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.53665390E-02, 0.47898560E-01, 0.38332540E-02, 0.39184770E-02,
     +  0.16522570E-03,-0.14362230E-02, 0.12021110E-02, 0.34039390E-03,
     + -0.22930830E-02, 0.17100874E-04,-0.13411042E-03,-0.98473904E-03,
     +  0.40411470E-03, 0.32437520E-03,-0.34378190E-03, 0.32441872E-02,
     + -0.28367460E-02, 0.49930830E-03, 0.51967351E-03,-0.33154550E-03,
     + -0.40316310E-02, 0.71732460E-02,-0.28395662E-02, 0.31874051E-02,
     + -0.90817473E-02,-0.42950452E-03, 0.78157361E-04,-0.44147763E-03,
     +  0.84380940E-02, 0.36468161E-03,-0.41029570E-02, 0.11084733E-02,
     + -0.37947740E-03, 0.45436170E-02, 0.50299121E-02,-0.11544660E-01,
     +  0.11021663E-01, 0.70630444E-03,-0.58924040E-03,-0.31045850E-02,
     + -0.87476450E-05,-0.26181020E-03, 0.58590091E-03, 0.55629920E-03,
     +  0.73758740E-04, 0.42091422E-04,-0.77415620E-03, 0.90697820E-03,
     + -0.79739620E-03, 0.27977384E-03,-0.48338250E-02,-0.23607292E-02,
     +  0.94952620E-02,-0.88643620E-02, 0.20953812E-02,-0.47534490E-02,
     +  0.78862660E-03, 0.34981850E-02,-0.42419770E-02,-0.30842300E-02,
     +  0.96907170E-04,-0.70463700E-04,-0.21982310E-03,-0.68493123E-03,
     + -0.26917354E-04, 0.51106780E-05,-0.68188900E-04,-0.13638330E-03,
     +  0.15731320E-03, 0.17348043E-03, 0.58417540E-04, 0.39849710E-03,
     +  0.58847960E-03, 0.12313640E-02, 0.24592880E-03, 0.64438220E-04,
     +  0.12909910E-02,-0.29145930E-02, 0.13893423E-02, 0.64788213E-02,
     +  0.57010930E-02,-0.14029750E-02, 0.23427040E-03, 0.52368180E-03,
     +  0.34174450E-02,-0.12546134E-03, 0.79281740E-04, 0.33567514E-03,
     + -0.17260830E-02,-0.24587420E-02,-0.27607030E-02,-0.14122102E-02,
     + -0.10200104E-02, 0.17632550E-03, 0.26846840E-02,
     +      0.      /
      data ientry/0/
c
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
      sr_delta       =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)*x11                
     3  +coeff(  3)*x12                
     4  +coeff(  4)*x11*x21            
     5  +coeff(  5)    *x21            
     6  +coeff(  6)    *x22            
     7  +coeff(  7)        *x32        
     8  +coeff(  8)                *x51
      sr_delta       =sr_delta       
     9  +coeff(  9)    *x21*x32        
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)    *x21    *x41    
     3  +coeff( 12)        *x31*x41    
     4  +coeff( 13)*x13                
     5  +coeff( 14)    *x22*x31        
     6  +coeff( 15)    *x22    *x41    
     7  +coeff( 16)    *x21*x31*x41    
     8  +coeff( 17)    *x21    *x42    
      sr_delta       =sr_delta       
     9  +coeff( 18)    *x24            
     1  +coeff( 19)*x11*x22            
     2  +coeff( 20)*x11*x21*x31        
     3  +coeff( 21)    *x22*x32        
     4  +coeff( 22)    *x22*x31*x41    
     5  +coeff( 23)    *x22    *x42    
     6  +coeff( 24)    *x21*x33        
     7  +coeff( 25)    *x21*x32*x41    
     8  +coeff( 26)    *x23            
      sr_delta       =sr_delta       
     9  +coeff( 27)    *x21        *x51
     1  +coeff( 28)*x11*x21*x32        
     2  +coeff( 29)    *x21*x31*x42    
     3  +coeff( 30)*x11        *x43    
     4  +coeff( 31)    *x21    *x43    
     5  +coeff( 32)            *x44    
     6  +coeff( 33)    *x24*x31        
     7  +coeff( 34)*x11*x21*x33        
     8  +coeff( 35)    *x22*x33        
      sr_delta       =sr_delta       
     9  +coeff( 36)    *x22*x32*x41    
     1  +coeff( 37)    *x22*x31*x42    
     2  +coeff( 38)*x12        *x43    
     3  +coeff( 39)*x11*x21    *x43    
     4  +coeff( 40)    *x22    *x43    
     5  +coeff( 41)        *x31        
     6  +coeff( 42)*x11*x21    *x41    
     7  +coeff( 43)*x11    *x31*x41    
     8  +coeff( 44)        *x32*x41    
      sr_delta       =sr_delta       
     9  +coeff( 45)    *x23*x31        
     1  +coeff( 46)    *x23    *x41    
     2  +coeff( 47)        *x31*x43    
     3  +coeff( 48)*x12*x22*x31        
     4  +coeff( 49)*x11*x22*x32        
     5  +coeff( 50)*x11*x22    *x42    
     6  +coeff( 51)*x11*x21*x32*x41    
     7  +coeff( 52)    *x21*x34        
     8  +coeff( 53)    *x21*x33*x41    
      sr_delta       =sr_delta       
     9  +coeff( 54)    *x21*x32*x42    
     1  +coeff( 55)    *x24*x32        
     2  +coeff( 56)    *x24*x31*x41    
     3  +coeff( 57)*x12*x22    *x42    
     4  +coeff( 58)*x11*x22*x33        
     5  +coeff( 59)*x11*x21*x34        
     6  +coeff( 60)    *x22*x34        
     7  +coeff( 61)    *x22    *x42*x51
     8  +coeff( 62)*x11    *x34*x41    
      sr_delta       =sr_delta       
     9  +coeff( 63)        *x33        
     1  +coeff( 64)        *x31*x42    
     2  +coeff( 65)        *x31    *x51
     3  +coeff( 66)            *x41*x51
     4  +coeff( 67)*x12*x22            
     5  +coeff( 68)*x11*x21*x31*x41    
     6  +coeff( 69)*x11*x21    *x42    
     7  +coeff( 70)    *x25            
     8  +coeff( 71)        *x32    *x51
      sr_delta       =sr_delta       
     9  +coeff( 72)*x13*x21*x31        
     1  +coeff( 73)    *x24    *x41    
     2  +coeff( 74)    *x23    *x42    
     3  +coeff( 75)*x12*x24            
     4  +coeff( 76)            *x43*x51
     5  +coeff( 77)    *x24    *x42    
     6  +coeff( 78)*x11*x22*x32*x41    
     7  +coeff( 79)*x11*x22    *x43    
     8  +coeff( 80)*x11*x21*x33*x41    
      sr_delta       =sr_delta       
     9  +coeff( 81)    *x22*x33*x41    
     1  +coeff( 82)*x12*x24*x31        
     2  +coeff( 83)*x13*x23    *x41    
     3  +coeff( 84)*x11*x25    *x41    
     4  +coeff( 85)    *x21*x31*x44    
     5  +coeff( 86)*x15    *x32        
     6  +coeff( 87)*x12*x23*x31*x41    
     7  +coeff( 88)*x15        *x42    
     8  +coeff( 89)    *x24*x33        
      sr_delta       =sr_delta       
     9  +coeff( 90)*x11*x22*x34        
     1  +coeff( 91)*x12*x21*x33*x41    
     2  +coeff( 92)*x11*x21*x32*x43    
     3  +coeff( 93)*x13*x23*x32        
     4  +coeff( 94)*x15*x24            
     5  +coeff( 95)*x11*x21*x32*x44    
c
      return
      end
      function sr_theta       (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 96)
      data ncoeff/ 95/
      data avdat/ -0.5667348E-03/
      data xmin/
     1 -0.68849E+00,-0.25651E-01,-0.36952E-01,-0.26097E-01,-0.19972E-02,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.59465E+00, 0.21663E-01, 0.39547E-01, 0.33592E-01, 0.19993E-02,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.61547900E-02,-0.39358533E-03,-0.52651510E-01, 0.59574010E-02,
     + -0.33925680E-02, 0.19587000E-03,-0.13377490E-02, 0.19985250E-03,
     + -0.33492830E-02,-0.78905220E-02, 0.21223410E-01, 0.22851410E-02,
     + -0.14980230E-01, 0.79964930E-03, 0.23068140E-02, 0.20305060E-02,
     +  0.39470260E-02, 0.37385553E-01,-0.72060080E-02, 0.26649783E-02,
     +  0.14525401E-01, 0.13030903E-02,-0.33325200E-04,-0.34071632E-04,
     + -0.86914720E-02,-0.23428730E-01, 0.19037820E-01,-0.42094550E-03,
     +  0.10966430E-02,-0.51090940E-01, 0.43324630E-01,-0.20007130E-01,
     + -0.41739661E-02, 0.11548320E-02, 0.13253452E-01, 0.19610760E-01,
     +  0.59815230E-03, 0.28623780E-02,-0.45975264E-01, 0.54024410E-01,
     +  0.70728120E-03,-0.33018920E-03, 0.17372490E-02,-0.32482771E-02,
     + -0.55810830E-03,-0.90955280E-03,-0.14055081E-02, 0.31541343E-02,
     +  0.15402324E-02, 0.98779020E-03, 0.12407740E-02, 0.26594672E-02,
     +  0.20366710E-02,-0.16074473E-01,-0.10687710E-01, 0.40633961E-01,
     + -0.37451882E-04, 0.29387070E-02, 0.78064780E-02,-0.83552412E-02,
     + -0.35327740E-01,-0.15296893E-01, 0.17603063E-01,-0.95378150E-03,
     + -0.40191850E-02,-0.12935152E-03, 0.47210620E-03,-0.68563350E-02,
     + -0.23011960E-01, 0.68733841E-02,-0.13683640E-01,-0.16756820E-01,
     +  0.39540831E-02, 0.14574293E-01, 0.26477000E-04,-0.95003870E-04,
     + -0.49984320E-04, 0.65630650E-03,-0.10058440E-02,-0.24478950E-02,
     +  0.14664750E-02,-0.38381563E-02, 0.39633042E-02,-0.34995350E-03,
     + -0.11048311E-02,-0.15342640E-02,-0.44914680E-02,-0.31546650E-02,
     + -0.16346523E-02,-0.16055973E-02, 0.51038570E-03,-0.14972840E-03,
     +  0.16456860E-03,-0.11907791E-02, 0.62788980E-02,
     +      0.      /
      data ientry/0/
c
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
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
c
c                  function
c
      sr_theta       =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)*x11                
     3  +coeff(  3)    *x21            
     4  +coeff(  4)*x11*x21            
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21*x31        
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)*x14                
      sr_theta       =sr_theta       
     9  +coeff(  9)*x12                
     1  +coeff( 10)    *x21*x32        
     2  +coeff( 11)    *x21*x31*x41    
     3  +coeff( 12)*x11        *x42    
     4  +coeff( 13)    *x21    *x42    
     5  +coeff( 14)                *x51
     6  +coeff( 15)    *x24    *x41    
     7  +coeff( 16)*x11*x21*x31        
     8  +coeff( 17)    *x22*x31        
      sr_theta       =sr_theta       
     9  +coeff( 18)    *x22*x31*x41    
     1  +coeff( 19)    *x22    *x42    
     2  +coeff( 20)*x11*x22*x32        
     3  +coeff( 21)    *x24*x32        
     4  +coeff( 22)*x11*x22            
     5  +coeff( 23)    *x24            
     6  +coeff( 24)    *x23*x31        
     7  +coeff( 25)*x11*x21*x32        
     8  +coeff( 26)    *x22*x32        
      sr_theta       =sr_theta       
     9  +coeff( 27)    *x21*x33        
     1  +coeff( 28)*x11*x22    *x41    
     2  +coeff( 29)    *x23    *x41    
     3  +coeff( 30)    *x21*x32*x41    
     4  +coeff( 31)    *x21*x31*x42    
     5  +coeff( 32)    *x21    *x43    
     6  +coeff( 33)    *x24*x31        
     7  +coeff( 34)    *x23*x32        
     8  +coeff( 35)*x11*x21*x33        
      sr_theta       =sr_theta       
     9  +coeff( 36)    *x22*x33        
     1  +coeff( 37)*x14        *x41    
     2  +coeff( 38)*x11*x23    *x41    
     3  +coeff( 39)    *x22*x32*x41    
     4  +coeff( 40)    *x22*x31*x42    
     5  +coeff( 41)*x11*x21    *x43    
     6  +coeff( 42)        *x31        
     7  +coeff( 43)        *x32        
     8  +coeff( 44)        *x31*x41    
      sr_theta       =sr_theta       
     9  +coeff( 45)*x12*x21            
     1  +coeff( 46)    *x23            
     2  +coeff( 47)    *x22    *x41    
     3  +coeff( 48)        *x32*x41    
     4  +coeff( 49)*x11*x23            
     5  +coeff( 50)        *x33*x41    
     6  +coeff( 51)*x12        *x42    
     7  +coeff( 52)*x11        *x43    
     8  +coeff( 53)            *x44    
      sr_theta       =sr_theta       
     9  +coeff( 54)    *x21*x34        
     1  +coeff( 55)    *x23*x31*x41    
     2  +coeff( 56)    *x21*x33*x41    
     3  +coeff( 57)*x12*x21    *x42    
     4  +coeff( 58)*x11*x22    *x42    
     5  +coeff( 59)    *x23    *x42    
     6  +coeff( 60)*x11*x21*x31*x42    
     7  +coeff( 61)    *x21*x32*x42    
     8  +coeff( 62)    *x22    *x43    
      sr_theta       =sr_theta       
     9  +coeff( 63)    *x21*x31*x43    
     1  +coeff( 64)*x11        *x44    
     2  +coeff( 65)    *x21    *x44    
     3  +coeff( 66)*x11            *x51
     4  +coeff( 67)*x13*x22*x31        
     5  +coeff( 68)    *x23*x33        
     6  +coeff( 69)    *x24*x31*x41    
     7  +coeff( 70)*x11*x22    *x43    
     8  +coeff( 71)*x11*x22*x34        
      sr_theta       =sr_theta       
     9  +coeff( 72)*x11*x23*x32*x41    
     1  +coeff( 73)*x11*x24*x33        
     2  +coeff( 74)    *x24*x33*x41    
     3  +coeff( 75)            *x41    
     4  +coeff( 76)*x11    *x31        
     5  +coeff( 77)*x11        *x41    
     6  +coeff( 78)*x12    *x31        
     7  +coeff( 79)        *x33        
     8  +coeff( 80)*x11*x21    *x41    
      sr_theta       =sr_theta       
     9  +coeff( 81)*x11    *x31*x41    
     1  +coeff( 82)        *x31*x42    
     2  +coeff( 83)            *x43    
     3  +coeff( 84)*x12*x22            
     4  +coeff( 85)*x12    *x32        
     5  +coeff( 86)*x11    *x32*x41    
     6  +coeff( 87)*x11    *x31*x42    
     7  +coeff( 88)        *x31*x43    
     8  +coeff( 89)*x11*x23*x31        
      sr_theta       =sr_theta       
     9  +coeff( 90)*x12*x21*x32        
     1  +coeff( 91)*x12    *x33        
     2  +coeff( 92)        *x31    *x51
     3  +coeff( 93)            *x41*x51
     4  +coeff( 94)*x11*x24*x31        
     5  +coeff( 95)*x11*x22*x33        
c
      return
      end
      function sr_phi         (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 96)
      data ncoeff/ 95/
      data avdat/ -0.3156503E-03/
      data xmin/
     1 -0.68849E+00,-0.25651E-01,-0.36952E-01,-0.26097E-01,-0.19972E-02,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.59465E+00, 0.21663E-01, 0.39547E-01, 0.33592E-01, 0.19993E-02,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.42083943E-03,-0.27565330E-01,-0.24188941E-03,-0.63861790E-02,
     + -0.13437173E-03, 0.67416992E-02, 0.64842533E-02, 0.21034742E-02,
     + -0.75988462E-02, 0.41621520E-02, 0.10618580E-02, 0.24110370E-02,
     + -0.31698440E-01, 0.17392400E-01,-0.19868614E-01,-0.49877970E-02,
     + -0.83347074E-02, 0.68010380E-02, 0.54586282E-02, 0.10590250E-01,
     + -0.14843762E-02, 0.30485913E-02,-0.28913610E-03,-0.14169551E-01,
     +  0.13822440E-01, 0.10149641E-01, 0.79218740E-02, 0.48269554E-02,
     + -0.48301790E-01, 0.21751490E-01, 0.13186170E-01, 0.26844200E-01,
     + -0.12830421E-01,-0.24553100E-01, 0.58456314E-02,-0.42595653E-02,
     + -0.97609350E-02,-0.99692560E-02, 0.46616090E-02,-0.18503414E-02,
     + -0.54209921E-02, 0.14268880E-01,-0.55196800E-03, 0.33676290E-01,
     + -0.23273860E-01,-0.24764174E-02, 0.36334150E-02,-0.28691360E-02,
     +  0.36695500E-01,-0.15906900E-01, 0.35646174E-01, 0.46986630E-02,
     +  0.47021950E-02,-0.18369400E-01, 0.12865351E-01,-0.16341224E-01,
     + -0.20601840E-01, 0.55838744E-02,-0.49735420E-02,-0.40324233E-03,
     +  0.15765890E-02,-0.18141113E-01, 0.40572322E-01, 0.39657563E-02,
     + -0.37959903E-01, 0.14852480E-02,-0.15482810E-01, 0.68649742E-02,
     +  0.77106880E-03,-0.86312050E-03,-0.21172972E-01, 0.17218481E-01,
     + -0.16306610E-01,-0.11570890E-02,-0.40673732E-03, 0.62084900E-03,
     + -0.19433820E-01, 0.13067813E-01, 0.14417494E-01,-0.16758373E-01,
     +  0.11407280E-01,-0.13901030E-02, 0.21381342E-03,-0.81629981E-03,
     +  0.33844443E-02,-0.96629103E-02, 0.83339840E-02,-0.39241020E-01,
     +  0.86547313E-02, 0.83325831E-02,-0.14810110E-01,-0.18851241E-02,
     +  0.27666250E-01, 0.44510043E-02, 0.98047280E-02,
     +      0.      /
      data ientry/0/
c
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
c
c                  function
c
      sr_phi         =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x21            
     6  +coeff(  6)            *x42    
     7  +coeff(  7)*x11    *x31        
     8  +coeff(  8)    *x21*x31        
      sr_phi         =sr_phi         
     9  +coeff(  9)*x11        *x41    
     1  +coeff( 10)    *x21    *x41    
     2  +coeff( 11)*x12                
     3  +coeff( 12)    *x22            
     4  +coeff( 13)        *x31*x42    
     5  +coeff( 14)            *x43    
     6  +coeff( 15)*x11    *x31*x41    
     7  +coeff( 16)*x12    *x31        
     8  +coeff( 17)    *x22*x31        
      sr_phi         =sr_phi         
     9  +coeff( 18)    *x22    *x41    
     1  +coeff( 19)*x12    *x32        
     2  +coeff( 20)*x11*x22    *x44    
     3  +coeff( 21)        *x31*x41    
     4  +coeff( 22)*x11    *x32        
     5  +coeff( 23)    *x21*x32        
     6  +coeff( 24)    *x21*x31*x41    
     7  +coeff( 25)*x11        *x42    
     8  +coeff( 26)    *x21    *x42    
      sr_phi         =sr_phi         
     9  +coeff( 27)*x12        *x41    
     1  +coeff( 28)        *x34        
     2  +coeff( 29)        *x31*x43    
     3  +coeff( 30)            *x44    
     4  +coeff( 31)*x11*x21*x32        
     5  +coeff( 32)    *x22*x32        
     6  +coeff( 33)*x11*x21*x31*x41    
     7  +coeff( 34)    *x22*x31*x41    
     8  +coeff( 35)    *x22    *x42    
      sr_phi         =sr_phi         
     9  +coeff( 36)*x12*x21*x31        
     1  +coeff( 37)*x11*x22    *x41    
     2  +coeff( 38)    *x21*x33*x41    
     3  +coeff( 39)*x12*x21*x32        
     4  +coeff( 40)*x11*x24            
     5  +coeff( 41)        *x33        
     6  +coeff( 42)        *x32*x41    
     7  +coeff( 43)*x11*x21*x31        
     8  +coeff( 44)        *x32*x42    
      sr_phi         =sr_phi         
     9  +coeff( 45)    *x21*x31*x42    
     1  +coeff( 46)*x12    *x31*x41    
     2  +coeff( 47)*x11*x22*x31        
     3  +coeff( 48)        *x33*x42    
     4  +coeff( 49)    *x21*x32*x42    
     5  +coeff( 50)    *x22*x33        
     6  +coeff( 51)    *x22*x32*x41    
     7  +coeff( 52)*x11*x21*x31*x42    
     8  +coeff( 53)*x12        *x43    
      sr_phi         =sr_phi         
     9  +coeff( 54)*x11*x21    *x43    
     1  +coeff( 55)    *x23*x32        
     2  +coeff( 56)*x12*x21    *x42    
     3  +coeff( 57)*x11*x22    *x42    
     4  +coeff( 58)        *x35*x41    
     5  +coeff( 59)        *x34*x42    
     6  +coeff( 60)    *x21        *x51
     7  +coeff( 61)    *x22*x31*x43    
     8  +coeff( 62)    *x23*x33        
      sr_phi         =sr_phi         
     9  +coeff( 63)    *x23*x32*x41    
     1  +coeff( 64)*x11*x22*x31*x42    
     2  +coeff( 65)    *x23*x31*x42    
     3  +coeff( 66)*x11*x22    *x43    
     4  +coeff( 67)    *x24*x32        
     5  +coeff( 68)*x12*x23*x31        
     6  +coeff( 69)    *x21*x31    *x51
     7  +coeff( 70)    *x21    *x41*x51
     8  +coeff( 71)*x12*x22*x33        
      sr_phi         =sr_phi         
     9  +coeff( 72)*x12*x23    *x42    
     1  +coeff( 73)*x11*x23    *x44    
     2  +coeff( 74)        *x32        
     3  +coeff( 75)*x11*x21            
     4  +coeff( 76)*x12*x21            
     5  +coeff( 77)        *x33*x41    
     6  +coeff( 78)*x11    *x32*x41    
     7  +coeff( 79)    *x21*x32*x41    
     8  +coeff( 80)*x11    *x31*x42    
      sr_phi         =sr_phi         
     9  +coeff( 81)*x11        *x43    
     1  +coeff( 82)    *x23    *x41    
     2  +coeff( 83)*x12*x22            
     3  +coeff( 84)*x11*x23            
     4  +coeff( 85)        *x34*x41    
     5  +coeff( 86)        *x31*x44    
     6  +coeff( 87)            *x45    
     7  +coeff( 88)    *x21*x31*x43    
     8  +coeff( 89)*x11        *x44    
      sr_phi         =sr_phi         
     9  +coeff( 90)    *x21    *x44    
     1  +coeff( 91)*x11*x21*x33        
     2  +coeff( 92)*x12    *x32*x41    
     3  +coeff( 93)*x11*x21*x32*x41    
     4  +coeff( 94)*x12*x21*x31*x41    
     5  +coeff( 95)*x11*x22*x31*x41    
c
      return
      end
      function sr_y00         (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 96)
      data ncoeff/ 95/
      data avdat/ -0.2002243E-03/
      data xmin/
     1 -0.68849E+00,-0.25651E-01,-0.36952E-01,-0.26097E-01,-0.19972E-02,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.59465E+00, 0.21663E-01, 0.39547E-01, 0.33592E-01, 0.19993E-02,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.47637610E-02,-0.24139210E-02, 0.39455343E-01, 0.81486190E-02,
     +  0.46124344E-03, 0.51265214E-02,-0.65279980E-03,-0.15615960E-01,
     + -0.25264080E-01,-0.28040240E-01, 0.61612473E-02, 0.37110780E-02,
     + -0.15135674E-02, 0.63533070E-03,-0.35656530E-02, 0.43507350E-01,
     + -0.38047190E-01,-0.55765304E-02,-0.66188923E-02, 0.28967514E-01,
     +  0.27494732E-01,-0.16880463E-01,-0.13259400E-01, 0.51378412E-02,
     +  0.19876644E-01,-0.10137220E-01,-0.19970030E-02,-0.76853511E-02,
     + -0.55652583E-03, 0.33629500E-02,-0.16106910E-01, 0.43093220E-01,
     + -0.43327460E-01,-0.12709903E-01, 0.54017852E-01,-0.16466110E-01,
     + -0.14687750E-01,-0.26682700E-01,-0.44854990E-01, 0.15986840E-01,
     +  0.21064690E-01, 0.58373660E-01,-0.43576300E-02, 0.27699801E-02,
     + -0.90331844E-02, 0.13978273E-01, 0.16647680E-02, 0.15988571E-01,
     + -0.35444583E-01, 0.62353573E-02, 0.13907281E-01,-0.23805502E-01,
     +  0.16860783E-01, 0.42553580E-01, 0.48230011E-01,-0.28190633E-01,
     + -0.57699020E-01,-0.63779100E-01, 0.23718410E-01, 0.12285810E-01,
     +  0.24239584E-01, 0.78746480E-02,-0.19525444E-02, 0.12534132E-01,
     + -0.62728771E-02,-0.92795751E-02, 0.24393620E-01,-0.53190120E-01,
     +  0.74746330E-02, 0.32705810E-01, 0.29072051E-01,-0.19288050E-01,
     + -0.38201780E-01, 0.52746040E-02,-0.95423780E-02, 0.50795054E-03,
     + -0.61507790E-01,-0.12442970E-02, 0.17985640E-01,-0.65956730E-02,
     +  0.69120302E-02, 0.74499710E-02, 0.88582493E-01, 0.89102692E-03,
     + -0.98707510E-02, 0.35903714E-04,-0.69055613E-03, 0.28284694E-01,
     +  0.10543773E-01,-0.32304020E-01,-0.46752393E-01, 0.32931180E-01,
     + -0.11314000E-01, 0.14464933E-01, 0.26202853E-01,
     +      0.      /
      data ientry/0/
c
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
      x36 = x35*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
c
c                  function
c
      sr_y00         =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x21            
     6  +coeff(  6)        *x32        
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)            *x42    
      sr_y00         =sr_y00         
     9  +coeff(  9)*x11    *x31        
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)*x11        *x41    
     3  +coeff( 12)    *x21    *x41    
     4  +coeff( 13)*x12                
     5  +coeff( 14)*x11*x21            
     6  +coeff( 15)    *x22            
     7  +coeff( 16)        *x31*x42    
     8  +coeff( 17)            *x43    
      sr_y00         =sr_y00         
     9  +coeff( 18)*x11    *x32        
     1  +coeff( 19)    *x21*x32        
     2  +coeff( 20)*x11    *x31*x41    
     3  +coeff( 21)    *x21*x31*x41    
     4  +coeff( 22)*x11        *x42    
     5  +coeff( 23)    *x21    *x42    
     6  +coeff( 24)*x12    *x31        
     7  +coeff( 25)    *x22*x31        
     8  +coeff( 26)*x12        *x41    
      sr_y00         =sr_y00         
     9  +coeff( 27)*x11*x21    *x41    
     1  +coeff( 28)    *x22    *x41    
     2  +coeff( 29)*x12*x21            
     3  +coeff( 30)*x11*x22            
     4  +coeff( 31)        *x34        
     5  +coeff( 32)        *x33*x41    
     6  +coeff( 33)        *x32*x42    
     7  +coeff( 34)            *x44    
     8  +coeff( 35)    *x21*x31*x42    
      sr_y00         =sr_y00         
     9  +coeff( 36)*x11        *x43    
     1  +coeff( 37)*x12    *x32        
     2  +coeff( 38)*x11*x21*x32        
     3  +coeff( 39)    *x22*x32        
     4  +coeff( 40)*x12    *x31*x41    
     5  +coeff( 41)*x11*x21*x31*x41    
     6  +coeff( 42)    *x22*x31*x41    
     7  +coeff( 43)    *x22    *x42    
     8  +coeff( 44)*x12*x21*x31        
      sr_y00         =sr_y00         
     9  +coeff( 45)*x11*x22*x31        
     1  +coeff( 46)    *x23*x31        
     2  +coeff( 47)*x12*x21    *x41    
     3  +coeff( 48)*x11*x22    *x41    
     4  +coeff( 49)        *x34*x41    
     5  +coeff( 50)    *x21*x34        
     6  +coeff( 51)    *x21*x33*x41    
     7  +coeff( 52)    *x21*x32*x42    
     8  +coeff( 53)*x12    *x33        
      sr_y00         =sr_y00         
     9  +coeff( 54)*x11*x21*x33        
     1  +coeff( 55)    *x22*x33        
     2  +coeff( 56)*x12    *x32*x41    
     3  +coeff( 57)*x11*x21*x32*x41    
     4  +coeff( 58)    *x22*x32*x41    
     5  +coeff( 59)*x12    *x31*x42    
     6  +coeff( 60)*x11*x21*x31*x42    
     7  +coeff( 61)    *x22*x31*x42    
     8  +coeff( 62)*x11*x21    *x43    
      sr_y00         =sr_y00         
     9  +coeff( 63)*x12*x21*x32        
     1  +coeff( 64)*x12*x21    *x42    
     2  +coeff( 65)*x12*x22    *x41    
     3  +coeff( 66)    *x21*x35        
     4  +coeff( 67)    *x22*x31*x43    
     5  +coeff( 68)    *x23*x32*x41    
     6  +coeff( 69)*x11*x22*x31*x42    
     7  +coeff( 70)    *x23*x31*x42    
     8  +coeff( 71)    *x24*x32        
      sr_y00         =sr_y00         
     9  +coeff( 72)*x12*x22*x31*x41    
     1  +coeff( 73)    *x24*x31*x41    
     2  +coeff( 74)*x11*x23    *x42    
     3  +coeff( 75)    *x25*x31        
     4  +coeff( 76)    *x21        *x51
     5  +coeff( 77)    *x23*x34        
     6  +coeff( 78)    *x21*x31    *x51
     7  +coeff( 79)    *x23*x35        
     8  +coeff( 80)*x13*x22*x34        
      sr_y00         =sr_y00         
     9  +coeff( 81)    *x21    *x43*x51
     1  +coeff( 82)    *x21    *x44*x51
     2  +coeff( 83)    *x25*x36        
     3  +coeff( 84)        *x33        
     4  +coeff( 85)        *x32*x41    
     5  +coeff( 86)*x11*x21*x31        
     6  +coeff( 87)    *x23            
     7  +coeff( 88)        *x31*x43    
     8  +coeff( 89)    *x21*x33        
      sr_y00         =sr_y00         
     9  +coeff( 90)*x11    *x32*x41    
     1  +coeff( 91)    *x21*x32*x41    
     2  +coeff( 92)*x11    *x31*x42    
     3  +coeff( 93)    *x21    *x43    
     4  +coeff( 94)        *x35        
     5  +coeff( 95)        *x33*x42    
c
      return
      end
  
