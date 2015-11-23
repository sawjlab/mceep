C
C ----------------------------------------------------------------------
C     SUBROUTINE       GET_VECT
C     AUTHOR:          M.D. Nozar and P.E. Ulmer
C     DATE:            March 31, 1992
C     PURPOSE:
C                 This routine determines the wire chamber coordinates
C                 based on the current Transport vector and stores them
C                 in an array TVEC_CHMBR.  For each wire chamber specified
C                 in the file $mceep_dat/wc_fil (the filename associated
C                 with the
C                 Fortran name wc_fil is specified by the user in the input
C                 file), positional coordinates are calculated from the
C                 intersection point of the ray with that chamber.
C                 The effects of multiple scattering, resolution smearing
C                 and offsets are included based on the parameters
C                 in wc_fil.  The resulting positional coordinates
C                 at each wire chamber are used by the routine TRACK
C                 which reconstructs the Transport vector at the position
C                 corresponding to the point before the TRK element in the
C                 user input file.
C
C                 MULTIPLE SCATTERING is incorporated by adding a gaussian
C                 distributed random number to the angular coordinates of
C                 the ray.  This affects the intersection point of the ray
C                 with subsequent chambers.  The same gaussian width is
C                 assumed for all chambers of a given arm.  Multiple
C                 scattering affects the kinematic bins and not the
C                 cross sections.
C
C                 RESOLUTION SMEARING is included by adding a gaussian
C                 distributed random number to the positional coordinates
C                 defined by the ray's intersection with a given chamber.
C                 This affects only the coordinates for the given chamber
C                 and not the ray itself (in terms of the intersection points
C                 with subsequent chambers).  The same gaussian width is
C                 assumed for all chambers of a given arm.  Resolution
C                 affects the kinematic bins and not the cross sections.
C
C                 OFFSETS(*) are included by addition of a constant to a
C                 specified coordinate for a given chamber
C                 in each arm.  (Only one coordinate can be offset
C                 for each arm for a given run of the program.)
C                 * As of MCEEP Version 3.1, offsets affect the kinematic
C                 bins.  This is necessary for consistency since other
C                 spectrometer elements (like HRS, HRI) are active
C                 only in the kinematics call (ICALL=2).  It doesn't
C                 make sense to apply offsets for the first call since
C                 we need HRS to transport to the focal plane.
C                 Actually, the test on ICALL is now redundant, since
C                 TRK is only called for ICALL=2.
C
C                 The wire chamber parameters are read (in mceep.f) from
C                 the file $mceep_dat/wc_fil in the following order:
C
C                  X_DELTA:      Dispersion in cm/% for a given spectrometer.
C                                (Note that for a realistic traceback
C                                 - as in hrs_inv - X_DELTA is ignored.
C                                 For Transport matrix it is needed though.)
C                  WC_OFF:       A constant (in Transport units) to be
C                                added to a ray coordinate.
C                  WC_OFF_COORD: The positional coordinate to be 
C                                affected by the offset. 
C                  WC_OFF_CHMB:  The chamber at which the offset is to be
C                                applied.
C                  WC_RES:       A gaussian distributed random number
C                                (in centimeters) to be added to the
C                                positional coordinates of the ray.
C                                This parameter is FWHM.
C                  N_WC:         Number of wire chambers specified.
C                                A maximum of five chambers is allowed.
C
C                  For I_WC=1,N_WC: 
C                      WC_C_ANGLE:  
C                                The angle (in deg) between the wire chamber
C                                and the Transport XY plane.  Currently,
C                                it is assumed that the rotation axis is
C                                Transport Y.
C                      WC_W_ANGLE:
C                                The angle (in deg) between the two sets of
C                                wires (assumed to be orthogonal) relative
C                                to the wire chamber frame.  For instance,
C                                this angle is 45 deg for "UV" chambers.
C                      WC_LOC:
C                                The location (in cm) of a W.C. with respect
C                                to the previous element in the input file.
C                      WC_RAD:
C                                x/x0 for wire chamber where:
C                                   x0 = rad. length of "chamber material"
C                                   x  = actual thickness of wire chamber
C                                (The routine sig_mscat determines the
C                                 sigma of a gaussian mult. scatt. distr.
C                                 (wc_msct) based on particle type and
C                                 wc_rad.  A gaussian distributed random
C                                 number based on wc_msct is added to
C                                 the angular coordinates of the ray.)
C
C ----------------------------------------------------------------------
C
       SUBROUTINE GET_VECT(TVEC_TRK,ICALL,IARM)
C
       IMPLICIT          NONE
C
       INCLUDE           'wc.cmn'
C
       INTEGER           I,ICALL,IARM,I_WC
C
       DOUBLE PRECISION  TVEC_TRK(6),TVEC_TRKD(6)
       DOUBLE PRECISION  TVEC_SAV(6),TVEC_UV(6),TVEC_WC(6)
       DOUBLE PRECISION  DELTA_LOC,GASDEVV
       DOUBLE PRECISION  RAY_TMP1(3),RAY_TMP2(3)
C
C ----------------------------------------------------------------------
C      Initialize
C ----------------------------------------------------------------------
C
       RAY_TMP1(1) = 0.D0
       RAY_TMP1(2) = 0.D0
       RAY_TMP1(3) = 1.D0
C
       RAY_TMP2(1) = 0.D0
       RAY_TMP2(2) = 0.D0
       RAY_TMP2(3) = 1.D0
C
C ----------------------------------------------------------------------
C      Loop over the wire chambers for a given arm.
C ----------------------------------------------------------------------
C
       DO I_WC = 1,N_WC(IARM)
C
C ----------------------------------------------------------------------
C        Calculate the distance from this chamber to the previous one.
C        For the first chamber, this is defined to be the distance
C        from the location of the previous element in the input file.
C ----------------------------------------------------------------------
C
         IF (I_WC .EQ. 1) THEN             
            DELTA_LOC = WC_LOC(IARM,1)
         ELSE
            DELTA_LOC = WC_LOC(IARM,I_WC) - WC_LOC(IARM,I_WC-1)
         ENDIF                              
C                                           
C ----------------------------------------------------------------------
C        Drift to the location of the next wire chamber.
C        Normally, for the first chamber, the drift distance should
C        be zero (provided this chamber is located at the position
C        obtained after application of the previous element in the
C        input file).
C ----------------------------------------------------------------------
C
         CALL DRIFT(TVEC_TRK,TVEC_TRKD,DELTA_LOC)
C
C ----------------------------------------------------------------------
C        Rotate (and traceback, if necessary) from the Transport system
C        to the chamber coordinate system and then into the
C        system defined by the wires (for WC_W_ANGLE=45 deg, this gives
C        the UV coordinates, for example).  It is assumed that the
C        chamber consists of two orthogonal wire planes which have zero
C        separation.
C ----------------------------------------------------------------------
C
         CALL ROT_TRPT(TVEC_TRKD,TVEC_SAV,'Y',WC_C_CANGLE(IARM,I_WC),
     #              WC_C_SANGLE(IARM,I_WC),
     #              RAY_TMP1,RAY_TMP2,.FALSE.)
         CALL ROT_TRPT(TVEC_SAV,TVEC_UV,'Z',WC_W_CANGLE(IARM,I_WC),
     #              WC_W_SANGLE(IARM,I_WC),
     #              RAY_TMP1,RAY_TMP2,.FALSE.)
C
C ----------------------------------------------------------------------
C        Add offset to selected coordinate.
C        (As of MCEEP Version 3.1: affects kinematics only (ICALL=2).)
C ----------------------------------------------------------------------
C
         IF((ICALL .EQ. 2) .AND. (I_WC .EQ. WC_OFF_CHMB(IARM)))  
     #     TVEC_UV(WC_OFF_COORD(IARM)) = 
     #               TVEC_UV(WC_OFF_COORD(IARM)) + WC_OFF(IARM)
C
C ----------------------------------------------------------------------
C        Smear positional coordinates with gaussian resolution function.
C        (Affects kinematic bins only (ICALL=2).)
C ----------------------------------------------------------------------
C
         IF(ICALL .EQ. 2) THEN
            TVEC_UV(1) = TVEC_UV(1) + GASDEVV(WC_RES(IARM))
            TVEC_UV(3) = TVEC_UV(3) + GASDEVV(WC_RES(IARM))
         ENDIF
C
C ----------------------------------------------------------------------
C        Rotate back into the chamber coordinate system.
C ----------------------------------------------------------------------
C
         CALL ROT_TRPT(TVEC_UV,TVEC_WC,'Z',WC_W_CANGLE(IARM,I_WC),
     #              -WC_W_SANGLE(IARM,I_WC),
     #              RAY_TMP1,RAY_TMP2,.FALSE.)
C
C ----------------------------------------------------------------------
C        Save the resulting "chamber" vector.
C ----------------------------------------------------------------------
C
         DO I=1,6
           TVEC_CHMBR(I_WC,I) = TVEC_WC(I)
         ENDDO
C
C ----------------------------------------------------------------------
C        Is this the last wire chamber for this arm?  If YES then return.
C ----------------------------------------------------------------------
C
         IF(I_WC .EQ. N_WC(IARM)) GOTO 999
C
C ----------------------------------------------------------------------
C        Add multiple scattering to the angular coordinates.  This does
C        not affect the previously saved "chamber" vector but does
C        affect subsequent "chamber" vectors.  This is because only the
C        ray's orientation is affected which has no impact on the current
C        positional coordinates.
C        (Affects kinematic bins only (ICALL=2).)
C ----------------------------------------------------------------------
C
         IF(ICALL .EQ. 2) THEN
           TVEC_SAV(2) = TVEC_SAV(2) + GASDEVV(WC_MSCT(IARM,I_WC))
           TVEC_SAV(4) = TVEC_SAV(4) + GASDEVV(WC_MSCT(IARM,I_WC))
         ENDIF
C
C ----------------------------------------------------------------------
C        Rotate back (and traceback, if necessary) into the Transport
C        coordinate system and loop back for the next chamber.
C ----------------------------------------------------------------------
C
         CALL ROT_TRPT(TVEC_SAV,TVEC_TRK,'Y',WC_C_CANGLE(IARM,I_WC),
     #              -WC_C_SANGLE(IARM,I_WC),
     #              RAY_TMP1,RAY_TMP2,.FALSE.)
C      
       ENDDO
C
999    RETURN
       END
