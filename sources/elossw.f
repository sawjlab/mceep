      subroutine elossw(einc,impeout,eangle,imppout,pangle,model)

c subroutine: elossw
c     author: R.Di Salvo (adapted by L.Todor)
c       date: August 1999
c    purpose: calculate the energy lost by the incident
c             electron and the outgoing electron and proton in
c             the target windows, walls, endcap and spectro's windows. 
c  
c 
c     Modified 4/12/02 to add support for energy loss in
c     polarized He3 target (model 4)
c     Kathy McCormick
c 
c     Modified to add support for the cigar tube target (Model 5).
c     Hassan Ibrahim, December 7, 2004
c

      implicit none

      integer model 
C------------------------------------------------------------------------------
* Ingoing and outgoing particles momenta
C------------------------------------------------------------------------------
      double precision einc,impeout,imppout
C------------------------------------------------------------------------------
* Spectrometers angles
C------------------------------------------------------------------------------
      double precision eangle,pangle
C------------------------------------------------------------------------------
* Thickness of traversed materials
C------------------------------------------------------------------------------
      double precision thAlin,thAlout,thwall,thscatter,thAir,thCap
      double precision thickness,thickendcap
C------------------------------------------------------------------------------
* Density of traversed materials
C------------------------------------------------------------------------------
      double precision rhoAl,rhoGlass,rhoCap,rhoAir
      double precision el1,el2,el3,el4,el5
      double precision prl1,prl2,prl3,prl4,prl5    
      double precision loss_ener
C
      INCLUDE 'eloss.cmn'
C
      
C------------------------------------------------------------------------------
*  thickness and density of the entering and exiting windows 
*  of the target in Al (cm).  If target model = 4 (polarized
*  He3 target, the windows and walls are glass (average thickness
*  assumed to be 130 um)
C------------------------------------------------------------------------------
      rhoAl =  2.7D0    ! g/cm^3 
      rhoglass = 2.76D0   ! g/cm^3

      if (model.eq.2) then
        thAlin  =  0.007112D0  ! cm
        thAlout =  0.011400D0  ! cm
        thAlin  =  thAlin * rhoAl
        thAlout =  thAlout * rhoAl
      elseif (model.eq.4) then
        thAlin = 0.013D0       ! cm
        thAlout = 0.013D0      ! cm
        thAlin  =  thAlin * rhoGlass
        thAlout =  thAlout * rhoGlass
      elseif (model.eq.5) then
        thAlin  =  0.011400D0  ! cm
        thAlout =  0.014000D0  ! cm
        thAlin  =  thAlin * rhoAl
        thAlout =  thAlout * rhoAl  
      else 
        thAlin  =  0.033000D0  ! cm
        thAlout =  0.033000D0  ! cm
        thAlin  =  thAlin * rhoAl
        thAlout =  thAlout * rhoAl
      end if

C------------------------------------------------------------------------------
*  thickness of the walls of the target in Al (cm) + 4 layers superisolation
*  If target model = 4, wall thickness is glass cell wall
C------------------------------------------------------------------------------
      if (model.eq.2) then
        thwall  =  0.02032D0   ! cm
        thwall = thwall * rhoAl
      elseif (model.eq.4) then
        thwall  =  0.14D0      ! cm
        thwall = thwall * rhoGlass
      elseif (model.eq.5) then
        thwall =  0.01270D0   ! cm
        thwall = thwall * rhoAl
      else 
        thwall  =  0.03300D0   ! cm
        thwall = thwall * rhoAl 
      end if

C------------------------------------------------------------------------------
* thickness of the window of scattering chamber in Al (cm)
C------------------------------------------------------------------------------

      thscatter = 0.04064D0  ! cm
      thscatter = thscatter * rhoAl

C------------------------------------------------------------------------------
*  thickness and density of Kapton (cm) 
C------------------------------------------------------------------------------

      thCap = 0.03556D0   ! cm
      rhoCap = 1.4D0    ! g/cm^3

      thCap  = thCap * rhoCap

C------------------------------------------------------------------------------
*  thickness of the air (cm)
C------------------------------------------------------------------------------

      thAir  = 65.1D0    ! cm  
      rhoAir = 0.0012D0 ! g/cm^3

      thAir  = thAir * rhoAir

**********************************************
*    Incoming electron
**********************************************


* Incoming electron energy loss 
      
      if (model.eq.1) then
        elwalk(1)= 0.D0
      else 
        elwalk(1) = loss_ener('e',einc,thAlin,'alu')
      end if

**********************************************
*    Outgoing electron
**********************************************

      if(eangle.eq.0.)  eangle=0.1        ! to calculate thickness
      if(eangle.eq.90.) eangle=89.9       ! to calculate thickendcap

C------------------------------------------------------------------------------
* Outgoing electron energy loss in the walls of the target
C------------------------------------------------------------------------------

      if (model.eq.2 .or. model.eq.4 .or. model.eq.5) then
        thickness = thwall/abs(sin(eangle))
      else
        thickness = thwall
      endif
      el1 = loss_ener('e',impeout,thickness,'alu')

C------------------------------------------------------------------------------
* Outgoing electron energy loss in the endcap of the target
C------------------------------------------------------------------------------

      if (model.eq.2 .or. model.eq.4.or.model.eq.5) then
        thickendcap = thAlout/abs(cos(eangle)) 
      else
        thickendcap = thAlout
      endif
      el2 = loss_ener('e',impeout,thickendcap,'alu')
C------------------------------------------------------------------------------
* Outgoing electron energy loss in the scattering chamber window in Al
C------------------------------------------------------------------------------

      el3 = loss_ener('e',impeout,thscatter,'alu')

C------------------------------------------------------------------------------
* Outgoing electron energy loss in air
C------------------------------------------------------------------------------

      el4 = loss_ener('e',impeout,thAir,'air')

C------------------------------------------------------------------------------
* Outgoing electron energy loss in Capton
C------------------------------------------------------------------------------

      el5 = loss_ener('e',impeout,thCap,'luc')

C------------------------------------------------------------------------------
* Total outgoing electron energy loss
      if (model.eq.1) then
cpeu        elwalk(2)= el4+el5

cpeu - version 3.9
        elwalk(2)= el3+el4+el5   ! need to include scattering chamb. window
        elwalk(3)=elwalk(2) 
cpeu - version 3.9

      else 
*** Through walls :  el1=walls  el3=sc.ch.window el4=air el5=kapton
         elwalk(2) = el1 + el3 + el4 + el5  
*** Through endcap : el2=endcap el3=sc.ch.window el4=air el5=kapton
         elwalk(3) = el2 + el3 + el4 + el5
      end if  
**********************************************
*    Outgoing proton
**********************************************

      if(pangle.eq.0.)  pangle=0.1        ! to calculate thickness
      if(pangle.eq.90.) pangle=89.9       ! to calculate thickendcap

C------------------------------------------------------------------------------
* Outgoing proton energy loss in the walls of the target
C------------------------------------------------------------------------------

      if (model.eq.2 .or. model.eq.4 .or. model.eq.5) then
        thickness = thwall/abs(sin(pangle))
      else
        thickness = thwall
      endif
      prl1 = loss_ener('p',imppout,thickness,'alu')

C------------------------------------------------------------------------------
* Outgoing proton energy loss in the endcap of the target
C------------------------------------------------------------------------------

      if (model.eq.2 .or. model.eq.4 .or. model.eq.5) then
        thickendcap  = thAlout/abs(cos(pangle)) 
      else
        thickendcap  = thAlout 
      endif
      prl2 = loss_ener('p',imppout,thickendcap,'alu')

C------------------------------------------------------------------------------
* Outgoing proton energy loss in the scattering chamber window in Al
C------------------------------------------------------------------------------

      prl3 = loss_ener('p',imppout,thscatter,'alu')

C------------------------------------------------------------------------------
* Outgoing proton energy loss in air
C------------------------------------------------------------------------------

      prl4 = loss_ener('p',imppout,thAir,'air')

C------------------------------------------------------------------------------
* Outgoing proton energy loss in Capton
C------------------------------------------------------------------------------

      prl5 = loss_ener('p',imppout,thCap,'luc')
C------------------------------------------------------------------------------
* Total outgoing proton energy loss
*** Through walls :  prl1=walls  prl3=sc.ch.window prl4=air prl5=capton
C------------------------------------------------------------------------------
      if (model.eq.1) then
cpeu        elwalk(2)= prl4+prl5

cpeu - version 3.9
        elwalk(4)= prl3+prl4+prl5   ! need to include scattering chamb. window
        elwalk(5)=elwalk(4) 
cpeu - version 3.9

      else 
        elwalk(4) = prl1 + prl3 + prl4 + prl5

C------------------------------------------------------------------------------
*** Through endcap : prl2=endcap prl3=sc.ch.window prl4=air prl5=capton
C------------------------------------------------------------------------------
        elwalk(5) = prl2 + prl3 + prl4 + prl5
      end if 
      return
      end


**********************************************************
     
      function loss_ener(part,mom,epaiseur,mat)  

*********************************************************

c======= part= p ou e (proton ou electron)
c======= mom = particle momentum
c======= epaiseur = thickness of the traversed material in g/cm^2
c======= matt = material : 
c======= alu = aluminium
c======= air = air
c======= luc = lucite used for Capton 

      implicit none

      double precision el_al,el_air,el_lu
      double precision pr_al,pr_air,pr_lu
 
      double precision el_mom,pr_mom
      dimension el_mom(30),el_al(30),el_air(30),el_lu(30)
      dimension pr_mom(96),pr_al(96),pr_air(96),pr_lu(96)

      double precision mom
      double precision epaiseur,loss_ener
      character *1 part
      character *3 mat
      integer ie 
      double precision ener

C------------------------------------------------------------------------------
* Electron momentum
C------------------------------------------------------------------------------

      data el_mom/
     @0060.,0080.,
     @0100.,0200.,0300.,0400.,0500.,0600.,0800.,
     @1000.,1200.,1400.,1600.,1800.,
     @2000.,2200.,2400.,2600.,2800.,
     @3000.,3200.,3400.,3600.,3800.,
     @4000.,4200.,4400.,4600.,4800.,
     @5000./

C------------------------------------------------------------------------------
* Electron energy loss in Aluminium (Pierre+Titanic)
C------------------------------------------------------------------------------

      data el_al/
     @1.808,1.831,
     @1.849,1.902,1.933,1.954,1.971,1.984,2.005,
     @2.022,2.034,2.046,2.055,2.064,
     @2.072,2.079,2.085,2.091,2.096,
     @2.101,2.106,2.110,2.114,2.118, 
     @2.122,2.126,2.129,2.132,2.135,
     @2.138/
     
C------------------------------------------------------------------------------
* Electron energy loss in Air (Pierre)
* Starting from 1 GeV the eloss does not change
C------------------------------------------------------------------------------

      data el_air/
     @2.355,2.400,
     @2.433,2.520,2.564,2.593,2.614,2.630,2.655,
     @2.674,2.674,2.674,2.674,2.674,     
     @2.674,2.674,2.674,2.674,2.674,     
     @2.674,2.674,2.674,2.674,2.674,     
     @2.674,2.674,2.674,2.674,2.674,     
     @2.674/
     
C------------------------------------------------------------------------------
* Electron energy loss in Lucite (for Capton) (Pierre)
* Starting from 1 GeV the eloss does not change
C------------------------------------------------------------------------------

      data el_lu/
     @2.089,2.113,
     @2.132,2.189,2.223,2.247,2.265,2.281,2.304,
     @2.323,2.323,2.323,2.323,2.323,      
     @2.323,2.323,2.323,2.323,2.323,     
     @2.323,2.323,2.323,2.323,2.323,     
     @2.323,2.323,2.323,2.323,2.323,     
     @2.323/    

C------------------------------------------------------------------------------
* Proton momentum
C------------------------------------------------------------------------------

      data pr_mom/
     @0310.36,0325.93,0340.86,
     @0355.24,0369.12,0382.57,0395.62,0408.32,0420.69,0432.76,
     @0444.57,0456.12,0467.45,0478.55,0489.46,0500.18,0510.72,
     @0521.10,0531.32,0541.39,0551.33,0561.13,0570.81,0580.38,
     @0589.82,0599.17,0608.41,0617.55,0626.59,0635.55,0644.43,
     @0661.93,0679.13,0696.04,0712.70,0729.11,0745.30,0761.27,
     @0777.05,0792.63,0808.04,0838.37,0868.09,0897.28,0925.98,
     @0954.24,0982.09,1009.57,1036.70,1063.52,1090.05,1116.31,
     @1142.31,1168.07,1193.62,1218.96,1244.10,1269.06,1293.85,
     @1318.47,1342.94,1367.27,1391.46,1415.52,1439.45,1463.26,
     @1486.97,1510.56,1534.06,1557.45,1580.75,1603.97,1627.10,
     @1650.14,1673.11,1696.00,1809.44,1921.38,2032.08,2141.73,
     @2250.48,2358.45,2465.75,2572.46,2678.66,2784.39,2994.68,
     @3203.66,3411.56,3618.56,3824.82,4030.45,4235.55,4440.17,
     @4644.40,4848.27/

C------------------------------------------------------------------------------
* Proton energy loss in Aluminium (Pierre)
C------------------------------------------------------------------------------

      data pr_al/
     @009.621,008.937,008.358,
     @007.861,007.430,007.052,006.718,006.421,006.154,005.913,
     @005.695,005.497,005.315,005.148,004.995,004.852,004.720,
     @004.598,004.483,004.377,004.276,004.182,004.094,004.011,
     @003.932,003.858,003.787,003.721,003.657,003.597,003.540,
     @003.433,003.336,003.246,003.164,003.089,003.019,002.954,
     @002.894,002.838,002.786,002.691,002.607,002.532,002.466,
     @002.406,002.352,002.303,002.259,002.210,002.173,002.139,
     @002.108,002.080,002.053,002.028,002.005,001.894,001.964,
     @001.945,001.928,001.911,001.996,001.882,001.868,001.855,
     @001.843,001.832,001.821,001.811,001.801,001.792,001.983,
     @001.775,001.767,001.760,001.728,001.704,001.684,001.669,
     @001.656,001.647,001.639,001.633,001.628,001.625,001.620,
     @001.619,001.619,001.620,001.622,001.626,001.629,001.633,
     @001.638,001.642/

C------------------------------------------------------------------------------
* Proton energy loss in Air (Pierre)
C------------------------------------------------------------------------------

      data pr_air/
     @010.840,010.050,009.391,
     @008.824,008.333,007.903,007.523,007.185,006.883,006.610,
     @006.363,006.138,005.933,005.744,005.570,005.410,005.261,
     @005.122,004.993,004.873,004.760,004.654,004.554,004.461,
     @004.372,004.288,004.209,004.134,004.063,003.995,003.931,
     @003.811,003.701,003.601,003.509,003.424,003.346,003.273,
     @003.206,003.143,003.084,002.977,002.883,002.800,002.725,
     @002.658,002.598,002.543,002.494,002.448,002.407,002.368,
     @002.333,002.301,002.271,002.243,002.217,002.193,002.171,
     @002.150,002.130,002.112,002.095,002.079,002.063,002.049,
     @002.036,002.023,002.011,002.000,001.989,001.979,001.970,
     @001.961,001.952,001.944,001.910,001.883,001.862,001.846,
     @001.834,001.824,001.817,001.811,001.807,001.804,001.802,
     @001.802,001.805,001.809,001.814,001.819,001.826,001.833,
     @001.840,001.847/

C------------------------------------------------------------------------------
* Proton energy loss in Lucite for Capton (tables)
C------------------------------------------------------------------------------

      data pr_lu/
     @012.260,011.370,010.620,
     @009.972,009.413,008.924,008.492,008.108,007.764,007.455,
     @007.174,006.919,006.686,006.472,006.275,006.093,005.924,
     @005.767,005.621,005.484,005.356,005.236,005.124,005.017,
     @004.917,004.823,004.733,004.648,004.567,004.491,004.418,
     @004.282,004.158,004.045,003.941,003.845,003.757,003.674,
     @003.598,003.527,003.461,003.340,003.234,003.140,003.056,
     @002.980,002.912,002.850,002.794,002.742,002.696,002.652,
     @002.613,002.576,002.542,002.511,002.481,002.454,002.429,
     @002.405,002.382,002.361,002.341,002.322,002.304,002.287,
     @002.271,002.256,002.242,002.228,002.216,002.203,002.192,
     @002.181,002.171,002.161,002.118,002.085,002.058,002.037,
     @002.020,002.006,001.994,001.985,001.978,001.972,001.963,
     @001.959,001.956,001.956,001.957,001.958,001.961,001.964,
     @001.968,001.972/     

      
      if(part.eq.'P'.or.part.eq.'p')then
           ie=1
           if(pr_mom(1).ge.mom)then
                   loss_ener=0.
                   return
                 end if
 11        ie=ie+1               
           if(pr_mom(ie).le.mom.and.ie.eq.96)then
                if(mat.eq.'air'.or.mat.eq.'AIR')then
                   loss_ener=epaiseur*pr_air(ie)
                   return
                 end if
                if(mat.eq.'alu'.or.mat.eq.'ALU')then
                   loss_ener=epaiseur*pr_al(ie)
                   return
                endif
                if(mat.eq.'luc'.or.mat.eq.'LUC')then
                   loss_ener=epaiseur*pr_lu(ie)
                   return
                 else 
                   loss_ener=0.
                   return
                 end if
           end if
           if(pr_mom(ie).lt.mom)then
                   goto 11
           end if

                if(mat.eq.'air'.or.mat.eq.'AIR')then
                   ener=(pr_air(ie-1)-pr_air(ie))
                   ener=ener/(pr_mom(ie-1)-pr_mom(ie))
                   ener=(pr_air(ie-1)+(mom-pr_mom(ie-1))*ener)*epaiseur                     
                   loss_ener=ener
                   return
                end if
                if(mat.eq.'alu'.or.mat.eq.'ALU')then
                   ener=(pr_al(ie-1)-pr_al(ie))
                   ener=ener/(pr_mom(ie-1)-pr_mom(ie))
                   ener=(pr_al(ie-1)+(mom-pr_mom(ie-1))*ener)*epaiseur
                   loss_ener=ener
                   return
                endif
                if(mat.eq.'luc'.or.mat.eq.'LUC')then
                   ener=(pr_lu(ie-1)-pr_lu(ie))
                   ener=ener/(pr_mom(ie-1)-pr_mom(ie))
                   ener=(pr_lu(ie-1)+(mom-pr_mom(ie-1))*ener)*epaiseur
                   loss_ener=ener
                   return
                else 
                   loss_ener=0.
                   return
                end if


       else                  ! cas de lelectron
           ie=1
             if(el_mom(1).ge.mom)then
                  loss_ener=0.
                  return
           end if
 12        ie=ie+1               
           if(el_mom(ie).le.mom.and.ie.eq.30)then
                if(mat.eq.'air'.or.mat.eq.'AIR')then
                   loss_ener=epaiseur*el_air(ie)
                   return
                 end if
                if(mat.eq.'alu'.or.mat.eq.'ALU')then
                   loss_ener=epaiseur*el_al(ie)
                   return
                 end if
                if(mat.eq.'luc'.or.mat.eq.'LUC')then
                   loss_ener=epaiseur*el_lu(ie)                
                 else 
                   loss_ener=0.
                   return
            end if
          end if
         if(el_mom(ie).lt.mom)then
            goto 12
         end if


                if(mat.eq.'air'.or.mat.eq.'AIR')then
                   ener=(el_air(ie-1)-el_air(ie))
                   ener=ener/(el_mom(ie-1)-el_mom(ie))
                   ener=(el_air(ie-1)+(mom-el_mom(ie-1))*ener)*epaiseur  
                   loss_ener=ener
                   return
                 end if
                if(mat.eq.'alu'.or.mat.eq.'ALU')then
                   ener=(el_al(ie-1)-el_al(ie))
                   ener=ener/(el_mom(ie-1)-el_mom(ie))
                   ener=(el_al(ie-1)+(mom-el_mom(ie-1))*ener)*epaiseur
                   loss_ener=ener
                   return
                 end if
                if(mat.eq.'luc'.or.mat.eq.'LUC')then
                   ener=(el_lu(ie-1)-el_lu(ie))
                   ener=ener/(el_mom(ie-1)-el_mom(ie))
                   ener=(el_lu(ie-1)+(mom-el_mom(ie-1))*ener)*epaiseur
                   loss_ener=ener
                   return
                 else 
                    loss_ener=0.
                    goto 999 ! just to provide a path to that statement
                 end if
       end if     
 999   return
      end

