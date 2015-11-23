      subroutine compute_transp_to_cosy(transp_vec,p_cent,c_vec)
c
      implicit none
c
      double precision c_vec(5),transp_vec(6,2),p_cent
      double precision pmom,tg_th,tg_ph,one_over_sqrt
      double precision tkin,central_ke
c
      INCLUDE 'masses.cmn'
*
* compute the input vector for cosy from the particles parameters at 
* target:
*   tar_vec      change to  c_vec
*   | x   (cm)           | x (m)       w.r.t. the central ray (p_o)
*   | theta  (mr)        | a = p_x / p_o
*   | y  (cm)  ----->    | y (m)
*   | phi    (mr)        | b = p_y / p_o
*   | z     (cm)         | t time of flight (=0 for all particles)
*   | p    percent       | relative tkin kinetic energy
      pmom=p_cent*(1.d0+transp_vec(6,2)/100.)
      tg_th = dtan(transp_vec(2,2)/1000.) 
      tg_ph = dtan(transp_vec(4,2)/1000.)
      one_over_sqrt = 1.d0 
     s              / dsqrt(tg_th * tg_th + tg_ph * tg_ph + 1.d0)

* cosy first input coordinate: x at the target:
      c_vec(1) = transp_vec(1,2)/100.

* cosy second input coordinate a = p_x / p_o.
* cosine director along x axis:
*        p_x               tg_th
*  cx = ----- = -----------------------------
*         p      sqrt(tg_th^2 + tg_ph^2 + 1) 
      c_vec(2) = tg_th * one_over_sqrt * pmom / p_cent

* cosy third input coordinate: y at the target:
      c_vec(3) = dble(transp_vec(3,2)/100.)

* cosy fourth input coordinate b = p_y / p_o.
* cosine director along x axis:
*        p_y               tg_ph
*  cy = ----- = -----------------------------
*         p      sqrt(tg_th^2 + tg_ph^2 + 1) 
      c_vec(4) = tg_ph * one_over_sqrt * pmom / p_cent
        
* in this vector for cosy input coordinates, we add a sixth component
* which contains the kinetic energy of the particle:
      tkin = (dsqrt(pmom*pmom + eject_mass*eject_mass) - eject_mass)
      central_ke= (dsqrt(p_cent*p_cent + eject_mass*eject_mass) 
     >              - eject_mass)
      c_vec(5) = (tkin - central_ke) / central_ke  
c
      return
      end
