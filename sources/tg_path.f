       subroutine tg_path(xv,yv,zv,ang,lout,opt,model)
 
c subroutine: tg_path 
c     author: L.Todor
c       date: July 1999
c    purpose: calculate the path length through the target
c             on event basis considering the geometrical 
c             description of the target and the particle angle; 
c             2D at this moment
c       used: for energy loss & bremsstrahlung calculations
c     
c     xv,yv,zv - reaction vertex coordinates
c     ang - particle outgoing angle in xz plane; beam direction z ang=0
c     lout - distance till exiting the target
c     opt - flag set to 1 if particle exits through endcaps (model 2, 4 or 5)
c           0 otherwise  
c     model - target model
c     Modified by R. Suleiman to include the Hall A He tuna can target.
c 
c     Modified to add support for the cigar tube target (Model 5).
c     Hassan Ibrahim, December 7, 2004
c

       implicit none

       integer opt,model 
       double precision lout,d,ang,xv,yv,zv,piover2,arg
c
       integer side_out
       double precision xr,yr,zr,length,tube_diam,cap_radius,cap_height,
     #     theta,phi,dist_out
       double precision arg1,arg2,cs_phi,sn_phi,sn_theta,cs_theta,
     #     tube_radius,z_cent,z_back,z_front,dist_side,dist_base,
     #     dist_front,z_base,dist_cap,arg3,arg4
       logical min_dist
C-----------------------------------------------------------------------
C
	INCLUDE 'eloss.cmn'

       if (model.eq.2 .or. model.eq.4) then       ! LH2/LD2 Cryotarget 
                                                  ! Pol He3 target
         opt=0
         piover2=2.D0*atan(1.D0)
         if (ang.lt.piover2.and.ang.gt.0.D0) then
               lout=(radius-xv)/sin(ang)
               d=(zout-zv)/cos(ang)
               if (d.lt.lout) then 
                   lout=d
                   opt=1
               end if
         else if (ang.gt.(-piover2).and.ang.lt.0.D0) then 
               lout=(radius+xv)/sin(-ang)
               d=(zout-zv)/cos(ang)
               if (d.lt.lout) then 
                   lout=d
                   opt=1
               end if 
         else if (ang.lt.(-piover2)) then
               lout=(radius+xv)/sin(-ang)  
               d=(zv-zin)/abs(cos(ang))
               if (d.lt.lout) then 
                   lout=d
                   opt=1
               end if
         else if (ang.gt.(piover2)) then
               lout=(radius-xv)/sin(ang)
               d=(zv-zin)/abs(cos(ang))
               if (d.lt.lout) then 
                   lout=d
                   opt=1
               end if                 
         end if 

       else if (model.eq.3) then     ! He tuna can

         opt=0
cpeu                              ! protect against rare events (v3.8)
         arg = radius**2-zv**2*sin(ang)**2
cxxx         if(arg .lt. 0.d0 .and. arg .gt. -1.d-5) arg = 0.d0
         if(arg .lt. 0.d0) arg = 0.d0  ! a bit risky, but what the heck
cpeu
         lout = sqrt( arg )
     >          -zv*cos(ang)

      else if (model .eq. 5) then   ! LH2/LD2 Cryotarget cigar tube

c     Changing the variable names

         xr = xv
         yr = yv
         zr = zv
         length = zout-zin
         tube_diam = radius*2.0
         cap_radius =  tube_diam/2.0
         cap_height =  cap_radius

         theta = abs(ang)

         if (ang .ge. 0.d0) then
            phi = 0.d0
         else
            phi = 4.D0*atan(1.D0)
         end if

c     Cigar Tube Target.
c     
c     
c             |<--------------- length ---------------->|
c     
c    _ _       __________________________________ _
c     |       |                                  |  *
c     |       |                                  |    *
c     |       |                                  |     *
c     |       |                                  |      *
c tube_diam   |                    O             |      * <-- end cap
c     |       |                 (0,0,0)          |      *
c     |       |                                  |     *
c     |       |                                  |    *
c    _|_      |__________________________________|_ *
c     
c     
c     
c             |<-----z_front------>|             |<---->| cap_height
c                                  |
c                                  |<--z_base--->|
c                                  |
c                                  |<------z_back------>|
c     
c     

c     xr,yr,zr      reaction point x , y and z coordinates
c     length        length of the target along the beam
c     tube_diam     diameter of the tube
c     cap_radius    radius of the end cap
c     cap_height    height of the end cap
c     theta,phi     particle trajectory (spherical angles wrt beam)
c     dist_out      distance the particle traverses as it exits

c     side_out = 1  particle exits from the side wall
c     side_out = 2  particle exits from the end cap
c     side_out = 3  particle exits from the entrance window

         min_dist = .FALSE.

         cs_phi   = cos(phi)
         sn_phi   = sin(phi)
         cs_theta = cos(theta)
         sn_theta = sin(theta)

         tube_radius = tube_diam/2.d0
         z_back = length/2.d0
         z_front = - z_back
         z_base = z_back - cap_height
         z_cent = z_back - cap_radius
         arg1 = xr*cs_phi+yr*sn_phi
         arg2 = tube_radius**2 - (xr*sn_phi-yr*cs_phi)**2

c   

         if (sn_theta .NE. 0.d0) then
            dist_side = (-arg1+ sqrt(arg2))/sn_theta
            if (dist_side .LT. 0.d0) dist_side = 0.d0
         else
            dist_side = 1.d+10
         endif

c   

         if (cs_theta .NE. 0.d0) then
            dist_front = (z_front-zr)/cs_theta
            dist_base = (z_base-zr)/cs_theta
            if (dist_front .LT. 0.d0) dist_front = 0.d0
            if (dist_base .LT. 0.d0) dist_base = 0.d0
         else
            dist_front = 1.d+10
            dist_base = 1.d+10
         endif

c   

         if (zr .LT. z_base) then
            if (cs_theta .EQ. 1.d0) then
               side_out = 2
            else if (cs_theta .GT. 0.d0) then
               if (dist_side .LE. dist_base) then
                  side_out = 1
               else
                  side_out = 2
               endif
            else if (cs_theta .EQ. 0.d0) then
               side_out = 1
            else
               min_dist = .TRUE.
            endif
c   
         else if (zr .EQ. z_base) then
            if (cs_theta .GT. 0.d0) then
               side_out = 2
            else if (cs_theta .EQ. 0.d0) then
               side_out = 1
            else
               min_dist = .TRUE.
            endif

c   

         else
            if (cs_theta .GE. 0.d0) then
               side_out = 2
            else
               if (dist_side .LT. dist_base) then
                  side_out = 2
               else if (dist_side .EQ. dist_base) then
                  side_out = 1
               else
                  min_dist = .TRUE.
               endif
            endif
         endif

c

         if (side_out .EQ. 1) then
            dist_out = dist_side
         else if (side_out .EQ. 2) then
            arg3 = sn_theta * arg1 + cs_theta * (zr - z_cent)
            arg4 = xr*xr+yr*yr+zr*zr-2*zr*z_cent+z_cent*z_cent
     >                 -cap_radius*cap_radius
            dist_cap = - arg3 + sqrt(arg3*arg3-arg4)
            if  (dist_cap .LT. 0.0D0)  dist_cap = 0.0D0
            dist_out = dist_cap
         else if (min_dist) then
            dist_out = min(dist_side,dist_front)
            if (dist_side .LT. dist_front) then
               side_out = 1
            else
               side_out = 3
            endif
         endif

         lout = dist_out
         if (side_out .eq. 1) opt = 0
         if (side_out .eq. 2 .or. side_out .eq. 3 ) opt = 1

       endif
       return 
       end
