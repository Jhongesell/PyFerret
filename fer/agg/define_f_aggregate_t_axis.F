        SUBROUTINE DEFINE_F_AGGREGATE_T_AXIS(nsets, member_sets,
     .                  memb_filenames, axcalf, findices, tline,
     .                  tmatch, t1indices, tdims, status)

*  This software was developed by the Thermal Modeling and Analysis
*  Project(TMAP) of the National Oceanographic and Atmospheric
*  Administration's (NOAA) Pacific Marine Environmental Lab(PMEL),
*  hereafter referred to as NOAA/PMEL/TMAP.
*
*  Access and use of this software shall impose the following
*  obligations and understandings on the user. The user is granted the
*  right, without any fee or cost, to use, copy, modify, alter, enhance
*  and distribute this software, and any derivative works thereof, and
*  its supporting documentation for any purpose whatsoever, provided
*  that this entire notice appears in all copies of the software,
*  derivative works and supporting documentation.  Further, the user
*  agrees to credit NOAA/PMEL/TMAP in any publications that result from
*  the use of this software or in any product that includes this
*  software. The names TMAP, NOAA and/or PMEL, however, may not be used
*  in any advertising or publicity to endorse or promote any products
*  or commercial entity unless specific written permission is obtained
*  from NOAA/PMEL/TMAP. The user also understands that NOAA/PMEL/TMAP
*  is not obligated to provide the user with any support, consulting,
*  training or assistance of any kind with regard to the use, operation
*  and performance of this software nor to provide the user with any
*  updates, revisions, new versions or "bug fixes".
*
*  THIS SOFTWARE IS PROVIDED BY NOAA/PMEL/TMAP "AS IS" AND ANY EXPRESS
*  OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
*  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
*  ARE DISCLAIMED. IN NO EVENT SHALL NOAA/PMEL/TMAP BE LIABLE FOR ANY SPECIAL,
*  INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER
*  RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF
*  CONTRACT, NEGLIGENCE OR OTHER TORTUOUS ACTION, ARISING OUT OF OR IN
*  CONNECTION WITH THE ACCESS, USE OR PERFORMANCE OF THIS SOFTWARE.
*
* Assign the forecasted-time (T) axis from the time axes of the member
* datasets.  If a T-axis to use is not specified, create it from merging
* all the (hopefully overlapping) time axis coordinates in the given
* datasets.  Returns the forecasted-time axis, locations of datasets
* on that axis, and number of timesteps on that axis for each set.
*
* arguments:
*       nsets          - number of forecast datasets (input)
*       member_sets    - dataset numbers for each forecast dataset (input)
*       memb_filenames - array of dataset filenames (input)
*       axcalf         - time-of-forecast (F) axis used (input)
*       findices       - index of each dataset on the time-of-forecast (F) axis (input)
*       tline          - forecasted-times (T) axis to use/used (input/output)
*       tmatch         - method of matching times to tline (input, can be empty);
*                        eg, '@ASN', '@BIN', '@XACT', '@XACT:0.001'
*       t1indices      - index of each first-forecasted-time on the forecasted-time (T) axis (output)
*       tdims          - number of forecasted times for each forecast dataset (output)
*       status         - return status (output)

        IMPLICIT NONE

        include 'tmap_dims.parm'
        include 'ferret.parm'
        include 'errmsg.parm'
        include 'xdset_info.cmn_text'
        external xdset_info_data
        include 'xdset_info.cd_equiv'
        include 'xtm_grid.cmn_text'
        include 'xunits.cmn_text'
        include 'xvariables.cmn'
        include 'netcdf.inc'

* calling argument declarations:
        INTEGER nsets, member_sets(nsets), axcalf, findices(nsets),
     .          tline, t1indices(nsets), tdims(nsets), status
        CHARACTER*512 memb_filenames(nsets)
        CHARACTER*(*) tmatch

* local variable declarations:
        REAL*8  GET_LINE_COORD, TM_WORLD_AX, SECS_FROM_BC,
     .          tlobnd, tcoord, thibnd, delta, origin,
     .          olobnd, ocoord, ohibnd, setorigin,
     .          depsilon, epsil, bnd_epsil
        INTEGER*8 wksize
        INTEGER TM_GET_CALENDAR_ID, TM_UNITS_CAL, TM_LENSTR1,
     .          TM_FIND_LIKE_DYN_LINE,
     .          flen, slen, axcalt, ref_tax, calid, units,
     .          setidx, pre_exist, axlen, axidx, wkspc,
     .          nsteps, ntsteps, igrid, iaxis, idx, previdx
        CHARACTER TM_SECS_TO_DATE*20,
     .          calname*32, buff*256, units_str*64, origin_str*20
        LOGICAL ITSA_TRUEMONTH_AXIS, TM_DFPEQ_TOL,
     .          true_month, irreg_bnds, bad_bnds

        PARAMETER( depsilon = 2.0**(-42) )

        wkspc = unspecified_int4
        axcalt = tline

        IF ( axcalt .LT. 0 ) THEN
           epsil = depsilon
*          Be a little more relaxes about matching boundaries
           bnd_epsil = 16 * depsilon

*          Get the time axis for the first member set
           ref_tax = ds_time_axis(member_sets(1))
           IF ( ref_tax .LE. 0 ) THEN
              flen = TM_LENSTR1(memb_filenames(1))
              CALL ERRMSG(ferr_aggregate_error, status,
     .                    'no time axis for forecast dataset '//
     .                    memb_filenames(1)(1:flen), *5000)
           ENDIF
           calname = line_cal_name(ref_tax)
           calid = TM_GET_CALENDAR_ID(calname)
           IF ( calid .LE. 0 ) THEN
              slen = TM_LENSTR1(calname)
              flen = TM_LENSTR1(memb_filenames(1))
              CALL ERRMSG(ferr_aggregate_error, status,
     .                    'unknown calendar "'//calname(1:slen)//
     .                    '" in forecast dataset '//
     .                    memb_filenames(1)(1:flen), *5000)
           ENDIF
           origin_str = line_t0(ref_tax)
           origin = SECS_FROM_BC(origin_str, calid, status)
           IF ( status .NE. ferr_ok ) THEN
              slen = TM_LENSTR1(origin_str)
              flen = TM_LENSTR1(memb_filenames(1))
              CALL ERRMSG(ferr_aggregate_error, status,
     .                    'invalid time origin "'//origin_str(1:slen)//
     .                    '" in forecast dataset '//
     .                    memb_filenames(1)(1:flen), *5000)
           ENDIF
*          For now, change units of true-month to units of days
           true_month = ITSA_TRUEMONTH_AXIS(ref_tax)
           IF ( true_month ) THEN
              units = pun_day
              units_str = un_name(units)
              true_month = .FALSE.
           ELSE
              units = line_unit_code(ref_tax)
              units_str = line_units(ref_tax)
           ENDIF
*          Change units of year to the appropriate number of days
           units = TM_UNITS_CAL(units, calid, true_month)
           IF  ( units .GE. 0 ) THEN
              slen = TM_LENSTR1(units_str)
              flen = TM_LENSTR1(memb_filenames(1))
              CALL ERRMSG(ferr_aggregate_error, status,
     .                    'unknown time units "'//units_str(1:slen)//
     .                    '" in forecast dataset '//
     .                    memb_filenames(1), *5000)
           ENDIF

*          Assume the time axis of each forecast is the same length and
*          allocate space for each time step in each forecast; normally
*          these axes will contain many overlaps, so should be plenty
           nsteps = line_dim(ref_tax)

           CALL TM_ALLO_TMP_LINE(axcalt, status)
           CALL TM_NEW_LINE_NAME('TF_CAL_T', buff)
           line_name(axcalt)         = buff
           line_name_orig(axcalt)    = buff
           line_dattype(axcalt)      = nf_double
           line_modulo(axcalt)       = .FALSE.
           line_modulo_len(axcalt)   = 0.0
           line_shift_origin(axcalt) = .FALSE.
           line_dim_only(axcalt)     = .FALSE.
           line_cal_name(axcalt)     = calname
           line_unit_code(axcalt)    = units
           line_units(axcalt)        = un_name(units)
           line_tunit(axcalt)        = un_convert(units)
           line_direction(axcalt)    = 'TI'
           line_t0(axcalt)           = origin_str
*          start with irregular axis, then check later if regular
           line_regular(axcalt)      = .FALSE.
           CALL GET_LINE_DYNMEM(nsets * nsteps, axcalt, status)
           IF ( status .NE. ferr_ok )
     .        CALL ERRMSG(ferr_aggregate_error, status,
     .                'out of memory for forecasted-time axis', *5000)

*          Get additional memory for the cell upper bounds until axis is complete
           wksize = nsets * nsteps
           CALL CREATE_DYN_WRKSPC(wksize, wkspc, status)
           IF ( status .NE. ferr_ok ) THEN
*             No memory allocated, so do not attempt to free it
              wkspc = unspecified_int4
              CALL ERRMSG(ferr_aggregate_error, status,
     .                'out of memory for forecasted-time upper bounds',
     .                *5000)
           ENDIF

*          Add all unique time coordinates from all forecast sets
           ntsteps = 0
           irreg_bnds = .FALSE.
           bad_bnds = .FALSE.
           DO setidx = 1,nsets
*             Get the time axis information for this forecast dataset
              iaxis = ds_time_axis(member_sets(setidx))
              IF ( iaxis .LE. 0 ) THEN
                 flen = TM_LENSTR1(memb_filenames(setidx))
                 CALL ERRMSG(ferr_aggregate_error, status,
     .                   'no valid time axis for forecast dataset '//
     .                   memb_filenames(setidx)(1:flen), *5000)
              ENDIF
              calname = line_cal_name(iaxis)
              calid = TM_GET_CALENDAR_ID(calname)
              IF ( calid .LE. 0 ) THEN
                 slen = TM_LENSTR1(calname)
                 flen = TM_LENSTR1(memb_filenames(setidx))
                 CALL ERRMSG(ferr_aggregate_error, status,
     .                   'unknown calendar "'//calname(1:slen)//
     .                   '" for forecast dataset '//
     .                   memb_filenames(setidx)(1:flen), *5000)
              ENDIF
              origin_str = line_t0(iaxis)
              setorigin = SECS_FROM_BC(origin_str, calid, status)
              IF ( status .NE. ferr_ok ) THEN
                 slen = TM_LENSTR1(origin_str)
                 flen = TM_LENSTR1(memb_filenames(setidx))
                 CALL ERRMSG(ferr_aggregate_error, status,
     .                   'invalid time origin "'//origin_str(1:slen)//
     .                   '" for forecast dataset '//
     .                   memb_filenames(setidx)(1:flen), *5000)
              ENDIF

              axidx = 1
              DO idx = 1,line_dim(iaxis)
*                Get this forecasted time (in units of SECS_FROM_BC)
                 tlobnd = TM_WORLD_AX(idx, iaxis, BOX_LO_LIM)
                 tcoord = TM_WORLD_AX(idx, iaxis, BOX_MIDDLE)
                 thibnd = TM_WORLD_AX(idx, iaxis, BOX_HI_LIM)
                 IF ( ITSA_TRUEMONTH_AXIS(iaxis) ) THEN
                    tlobnd = tlobnd * un_convert(pun_day)
                    tcoord = tcoord * un_convert(pun_day)
                    thibnd = thibnd * un_convert(pun_day)
                 ELSE
                    tlobnd = tlobnd * line_tunit(iaxis)
                    tcoord = tcoord * line_tunit(iaxis)
                    thibnd = thibnd * line_tunit(iaxis)
                 ENDIF
                 tlobnd = tlobnd + setorigin
                 tcoord = tcoord + setorigin
                 thibnd = thibnd + setorigin

*                If the coordinate is not in the middle of the cell, cannot be regular
                 IF ( .NOT. TM_DFPEQ_TOL(tcoord-tlobnd,
     .                thibnd-tcoord, bnd_epsil) ) irreg_bnds = .TRUE.

*                Match or add this forecasted time (in units of SECS_FROM_BC)
                 DO WHILE ( axidx .LE. ntsteps )
                    ocoord = GET_LINE_COORD(linemem(axcalt)%ptr,axidx)
*                   If matches an existing value, no changes needed
                    IF ( TM_DFPEQ_TOL(ocoord, tcoord, epsil) ) THEN
*                      Lower limits in lineedg, upper limits in workmem
                       olobnd = GET_LINE_COORD(lineedg(axcalt)%ptr, axidx)
                       ohibnd = GET_LINE_COORD(workmem(wkspc)%ptr, axidx)
*                      If the bounds do not match, cannot use the bounds for the given axes
                       IF ( .NOT. (
     .                     TM_DFPEQ_TOL(olobnd, tlobnd, bnd_epsil) .AND.
     .                     TM_DFPEQ_TOL(ohibnd, thibnd, bnd_epsil) ) )
     .                     bad_bnds = .TRUE.
                       GOTO 340
                    ENDIF
*                   if new value is smaller, insert it here
                    IF ( ocoord .GT. tcoord ) GOTO 320
                    axidx = axidx + 1
                 ENDDO
*                insert the value, moving any value after it
  320            IF ( ntsteps .EQ. nsets*nsteps )
     .              CALL ERRMSG(ferr_aggregate_error, status,
     .                          'unexpected large number of unique '//
     .                          'forecasted times', *5000)
                 ntsteps = ntsteps + 1
*                insert/append the new value
                 previdx = axidx
                 DO WHILE ( previdx .LE. ntsteps )
*                   Lower limits in lineedg, upper limits in workmem, to delay
*                   checking upper n == lower n+1 until all are in place
                    CALL PUT_LINE_COORD(lineedg(axcalt)%ptr, previdx, tlobnd)
                    CALL PUT_LINE_COORD(linemem(axcalt)%ptr, previdx, tcoord)
                    CALL PUT_LINE_COORD(workmem(wkspc)%ptr, previdx, thibnd)
                    tlobnd = olobnd
                    tcoord = ocoord
                    thibnd = ohibnd
                    previdx = previdx + 1
                    IF ( previdx .LT. ntsteps ) THEN
                       olobnd = GET_LINE_COORD(lineedg(axcalt)%ptr, previdx)
                       ocoord = GET_LINE_COORD(linemem(axcalt)%ptr, previdx)
                       ohibnd = GET_LINE_COORD(workmem(wkspc)%ptr, previdx)
                    ENDIF
                 ENDDO
  340            CONTINUE
              ENDDO
           ENDDO

*          Check for problems with the bounds
           IF ( .NOT. bad_bnds ) THEN
              axidx = 1
*             upper limit of previous cell
  350         ohibnd = GET_LINE_COORD(workmem(wkspc)%ptr, axidx)
              axidx = axidx + 1
              IF ( axidx .GT. ntsteps ) GOTO 360
*             lower limit of current cell
              tlobnd = GET_LINE_COORD(lineedg(axcalt)%ptr, axidx)

              IF ( TM_DFPEQ_TOL(tlobnd, ohibnd, bnd_epsil) ) THEN
*                Assign the average in case there is any difference
                 CALL PUT_LINE_COORD(lineedg(axcalt)%ptr, axidx, 0.5*(ohibnd+tlobnd))
              ELSEIF ( ohibnd .GT. tlobnd ) THEN
*                Cells overlap - bounds cannot be used
                 bad_bnds = .TRUE.
                 GOTO 360
              ELSE
*                Gap between cells - add a midpoint coordinate
                 IF ( ntsteps .EQ. nsets*nsteps )
     .              CALL ERRMSG(ferr_aggregate_error, status,
     .                          'unexpected large number of unique '//
     .                          'forecasted times', *5000)
                 ntsteps = ntsteps + 1
*                Cell to add to be specified by tlobnd, tcoord, thibnd
*                Original cell to follow it to be specified by olobnd, ocoord, ohibnd
                 ocoord = GET_LINE_COORD(linemem(axcalt)%ptr, axidx)
                 tcoord = 0.5 * (ohibnd + tlobnd)
                 olobnd = tlobnd
                 tlobnd = ohibnd
                 ohibnd = GET_LINE_COORD(workmem(wkspc)%ptr, axidx)
                 thibnd = olobnd
                 previdx = axidx
                 DO WHILE ( previdx .LE. ntsteps )
*                   Lower limits in lineedg, upper limits in workmem
                    CALL PUT_LINE_COORD(lineedg(axcalt)%ptr, previdx, tlobnd)
                    CALL PUT_LINE_COORD(linemem(axcalt)%ptr, previdx, tcoord)
                    CALL PUT_LINE_COORD(workmem(wkspc)%ptr, previdx, thibnd)
                    tlobnd = olobnd
                    tcoord = ocoord
                    thibnd = ohibnd
                    previdx = previdx + 1
                    IF ( previdx .LT. ntsteps ) THEN
                       olobnd = GET_LINE_COORD(lineedg(axcalt)%ptr, previdx)
                       ocoord = GET_LINE_COORD(linemem(axcalt)%ptr, previdx)
                       ohibnd = GET_LINE_COORD(workmem(wkspc)%ptr, previdx)
                    ENDIF
                 ENDDO
*                Do not need to check the cell just added - continue
*                with the high bound of the cell just repositioned
                 axidx = axidx + 1
              ENDIF
              GOTO 350
           ENDIF

  360      IF ( .NOT. bad_bnds ) THEN
*             Bounds now mesh together - only need the last upper bound from workmem
              ohibnd = GET_LINE_COORD(workmem(wkspc)%ptr, ntsteps)
              CALL PUT_LINE_COORD(lineedg(axcalt)%ptr, ntsteps+1, ohibnd)
           ENDIF
           CALL FREE_WS_DYNMEM(wkspc)
           wkspc = unspecified_int4

           line_dim(axcalt) = ntsteps
*          convert the time in UNITS_FROM_BC to the appropriate time origin
           DO axidx = 1,ntsteps
              ocoord = GET_LINE_COORD(linemem(axcalt)%ptr, axidx)
              ocoord = ocoord - origin
              IF ( ITSA_TRUEMONTH_AXIS(axcalt) ) THEN
                 ocoord = ocoord / un_convert(pun_day)
              ELSE
                 ocoord = ocoord / line_tunit(axcalt)
              ENDIF
              CALL PUT_LINE_COORD(linemem(axcalt)%ptr, axidx, ocoord)
           ENDDO
           IF ( .NOT. bad_bnds ) THEN
              DO axidx = 1,ntsteps+1
                 ocoord = GET_LINE_COORD(lineedg(axcalt)%ptr, axidx)
                 ocoord = ocoord - origin
                 IF ( ITSA_TRUEMONTH_AXIS(axcalt) ) THEN
                    ocoord = ocoord / un_convert(pun_day)
                 ELSE
                    ocoord = ocoord / line_tunit(axcalt)
                 ENDIF
                 CALL PUT_LINE_COORD(lineedg(axcalt)%ptr, axidx, ocoord)
              ENDDO
*             If any coordinate was not in the center of its cell, then
*             the axes cannot be a regular regardless of coordinate spacing
              IF ( irreg_bnds ) GOTO 380
           ELSE
*             Cannot use existing cell bounds; use the "midpoints"
              CALL TM_IRREGAX_MID(linemem(axcalt)%ptr,
     .                            lineedg(axcalt)%ptr, ntsteps)
           ENDIF

*          Check if this is a regular axis
           ocoord = GET_LINE_COORD(linemem(axcalt)%ptr, 1)
           IF ( ntsteps .GT. 1 ) THEN
              tcoord = GET_LINE_COORD(linemem(axcalt)%ptr, ntsteps)
              delta = (tcoord - ocoord) / (ntsteps - 1.0)
           ELSE
              delta = 1.0
           ENDIF
           DO axidx = 2,ntsteps-1
              tcoord = GET_LINE_COORD(linemem(axcalt)%ptr, axidx)
              IF ( .NOT. TM_DFPEQ_TOL(tcoord,
     .             ocoord + (axidx-1.0) * delta, epsil) ) GOTO 380
           ENDDO
           IF ( .NOT. bad_bnds ) THEN
*             Check if the provided cell bounds are at the "midpoints"
              DO axidx = 1,ntsteps+1
                 tlobnd = GET_LINE_COORD(lineedg(axcalt)%ptr, axidx)
                 IF ( .NOT. TM_DFPEQ_TOL(tlobnd,
     .             ocoord + (axidx-1.5) * delta, bnd_epsil) ) GOTO 380
              ENDDO
           ENDIF

*          If we got here, the axis is regular (including bounds)
           CALL FREE_LINE_DYNMEM(axcalt)
           line_regular(axcalt) = .TRUE.
           line_start(axcalt)   = ocoord
           line_delta(axcalt)   = delta

*          If the axis is an exact duplicate of a pre-existing axis, use that axis instead
*          (shared axes are needed in order to E-aggregate multiple FMRCs)
  380      pre_exist = TM_FIND_LIKE_DYN_LINE(axcalt)
           IF ( pre_exist .NE. unspecified_int4 ) THEN
              IF ( .NOT. line_regular(axcalt) )
     .           CALL FREE_LINE_DYNMEM(axcalt)
              CALL TM_DEALLO_DYN_LINE(axcalt)
              axcalt = pre_exist
*             Also set tline so it will not be deleted if an error occurs later
              tline = pre_exist
           ELSE
*             New axis; move the tmp_line to a dyn_line
              CALL TM_RE_ALLO_TMP_LINE(axcalt)
           ENDIF

        ENDIF

*       Map the forecasted times to the created or provided T axis

        axlen = line_dim(axcalt)
        IF ( axlen .LT. 1 )
     .     CALL ERRMSG(ferr_aggregate_error, status,
     .          'invalid forecasted-time (T) axis', *5000)
        IF ( line_direction(axcalt) .NE. 'TI' )
     .     CALL ERRMSG(ferr_aggregate_error, status,
     .          'forecasted-time (T) axis not a time axis', *5000)
        calname = line_cal_name(axcalt)
        calid = TM_GET_CALENDAR_ID(calname)
        IF ( calid .LE. 0 ) THEN
           slen = TM_LENSTR1(calname)
           CALL ERRMSG(ferr_aggregate_error, status,
     .                 'unknown calendar "'//calname(1:slen)//
     .                 '" for provided T axis', *5000)
        ENDIF
        origin_str = line_t0(axcalt)
        origin = SECS_FROM_BC(origin_str, calid, status)
        IF ( status .NE. ferr_ok ) THEN
           slen = TM_LENSTR1(origin_str)
           CALL ERRMSG(ferr_aggregate_error, status,
     .                 'invalid time origin "'//origin_str(1:slen)//
     .                 '" for provided T axis', *5000)
        ENDIF

        IF ( tline .LT. 0 ) THEN
            epsil = depsilon
        ELSEIF ( tmatch .EQ. '@ASN' ) THEN
            GOTO 900
        ELSEIF ( tmatch .EQ. '@BIN' ) THEN
            epsil = -1.0
        ELSEIF ( tmatch .EQ. '@XACT') THEN
            epsil = depsilon
        ELSEIF ( tmatch(1:6) .EQ. '@XACT:' ) THEN
           READ(tmatch(7:), *, ERR=5010) epsil
           IF ( epsil .LT. 0.0 ) GOTO 5010
        ELSEIF ( tmatch .EQ. ' ' ) THEN
           epsil = depsilon
*          Default is to try @XACT; if not successful, try @BIN; if not successful, use @ASN
        ELSE
           GOTO 5010
        ENDIF

*       Check time coordinates from all forecast sets
        DO setidx = 1,nsets
*          get the time axis information for this forecast dataset
           iaxis = ds_time_axis(member_sets(setidx))
           IF ( iaxis .LE. 0 ) THEN
*             @ASN only valid option (assumes a single time step)
              IF ( (tline .GE. 0) .AND. (tmatch .EQ. ' ') ) GOTO 900
              flen = TM_LENSTR1(memb_filenames(setidx))
              CALL ERRMSG(ferr_aggregate_error, status,
     .                    'No valid time axis for forecast dataset '//
     .                    memb_filenames(setidx)(1:flen), *5000)
           ENDIF
           calname = line_cal_name(iaxis)
           calid = TM_GET_CALENDAR_ID(calname)
           IF ( calid .LE. 0 ) THEN
*             @ASN only valid option
              IF ( (tline .GE. 0) .AND. (tmatch .EQ. ' ') ) GOTO 900
              slen = TM_LENSTR1(calname)
              flen = TM_LENSTR1(memb_filenames(setidx))
              CALL ERRMSG(ferr_aggregate_error, status,
     .                    'unknown calendar "'//calname(1:slen)//
     .                    '" for forecast dataset '//
     .                    memb_filenames(setidx)(1:flen), *5000)
           ENDIF
           origin_str = line_t0(iaxis)
           setorigin = SECS_FROM_BC(origin_str, calid, status)
           IF ( status .NE. ferr_ok ) THEN
*             @ASN only valid option
              IF ( (tline .GE. 0) .AND. (tmatch .EQ. ' ') ) GOTO 900
              slen = TM_LENSTR1(origin_str)
              flen = TM_LENSTR1(memb_filenames(setidx))
              CALL ERRMSG(ferr_aggregate_error, status,
     .                    'invalid time origin "'//origin_str(1:slen)//
     .                    '" for forecast dataset '//
     .                    memb_filenames(setidx)(1:flen), *5000)
           ENDIF

           IF ( epsil .GE. 0.0 ) THEN
*             Try matching with @XACT, or with generated axcalt

              axidx = 0
              DO idx = 1,line_dim(iaxis)
*                get this forecasted time (in SECS_FROM_BC)
                 tcoord = TM_WORLD_AX(idx, iaxis, BOX_MIDDLE)
                 IF ( ITSA_TRUEMONTH_AXIS(iaxis) ) THEN
                    tcoord = tcoord * un_convert(pun_day)
                 ELSE
                    tcoord = tcoord * line_tunit(iaxis)
                 ENDIF
                 tcoord = tcoord + setorigin
*                match to a time on the created/provided T axis
  400            axidx = axidx + 1
                 IF ( axidx .GT. axlen ) THEN
                    IF ( (tline .GE. 0) .AND. (tmatch .EQ. ' ') ) THEN
*                      @XACT failed; try @BIN - any previous set assignments made with @XACT
*                      will be the same with @BIN, so no need to exit the setidx DO loop
                       epsil = -1.0
                       GOTO 500
                    ENDIF
                    origin_str = TM_SECS_TO_DATE(tcoord, calid)
                    slen = TM_LENSTR1(origin_str)
                    flen = TM_LENSTR1(memb_filenames(setidx))
                    CALL ERRMSG(ferr_aggregate_error, status,
     .                      'incompatible forecasted-time (T) axis: '//
     .                      'unmatched time "'//origin_str(1:slen)//
     .                      '" in forecast dataset '//
     .                       memb_filenames(setidx)(1:flen), *5000)
                 ENDIF
                 ocoord = TM_WORLD_AX(axidx, axcalt, BOX_MIDDLE)
                 IF ( ITSA_TRUEMONTH_AXIS(axcalt) ) THEN
                    ocoord = ocoord * un_convert(pun_day)
                 ELSE
                    ocoord = ocoord * line_tunit(axcalt)
                 ENDIF
                 ocoord = ocoord + origin
                 IF ( idx .EQ. 1 ) THEN
                    IF ( TM_DFPEQ_TOL(ocoord, tcoord, epsil) ) THEN
*                      First time coordinate matched
                       t1indices(setidx) = axidx
                       tdims(setidx) = line_dim(iaxis)
                    ELSEIF ( tcoord .LT. ocoord ) THEN
*                      Already past first time coordinate - stop the search
                       axidx = axlen
                       GOTO 400
                    ELSE
*                      Keep searching for the first time coordinate
                       GOTO 400
                    ENDIF
                 ELSE
*                   At this time, time coordinates must be consecutive on the created/provided time axis
                    IF ( .NOT. TM_DFPEQ_TOL(ocoord, tcoord, epsil) ) THEN
                       IF ( (tline .GE. 0) .AND. (tmatch .EQ. ' ') ) THEN
*                         @XACT failed; try @BIN - any previous set assignments made with @XACT
*                         will be the same with @BIN, so no need to exit the setidx DO loop
                          epsil = -1.0
                          GOTO 500
                       ENDIF
                       origin_str = TM_SECS_TO_DATE(tcoord, calid)
                       slen = TM_LENSTR1(origin_str)
                       flen = TM_LENSTR1(memb_filenames(setidx))
                       CALL ERRMSG(ferr_aggregate_error, status,
     .                    'incompatible forecasted-time (T) axis: '//
     .                    'non-consecutive times around "'//
     .                    origin_str(1:slen)//'" in forecast dataset '//
     .                    memb_filenames(setidx)(1:flen), *5000)
                    ENDIF
                 ENDIF
              ENDDO

*             Success with @XACT on this dataset - go on to the next dataset
           ENDIF

  500      IF ( epsil .LT. 0.0 ) THEN
*             Try matching with @BIN (either @BIN specified or default matching and @XACT failed)

              axidx = 0
              DO idx = 1,line_dim(iaxis)
*                get this forecasted time (in SECS_FROM_BC)
                 tcoord = TM_WORLD_AX(idx, iaxis, BOX_MIDDLE)
                 IF ( ITSA_TRUEMONTH_AXIS(iaxis) ) THEN
                    tcoord = tcoord * un_convert(pun_day)
                 ELSE
                    tcoord = tcoord * line_tunit(iaxis)
                 ENDIF
                 tcoord = tcoord + setorigin
*                match to a time on the created/provided T axis
  600            axidx = axidx + 1
                 IF ( axidx .GT. axlen ) THEN
                    IF ( tmatch .EQ. ' ' ) THEN
*                      @BIN failed; use @ASN
                       GOTO 900
                    ENDIF
                    origin_str = TM_SECS_TO_DATE(tcoord, calid)
                    slen = TM_LENSTR1(origin_str)
                    flen = TM_LENSTR1(memb_filenames(setidx))
                    CALL ERRMSG(ferr_aggregate_error, status,
     .                      'incompatible forecasted-time (T) axis: '//
     .                      'unmatched time "'//origin_str(1:slen)//
     .                      '" in forecast dataset '//
     .                       memb_filenames(setidx)(1:flen), *5000)
                 ENDIF
                 olobnd = TM_WORLD_AX(axidx, axcalt, BOX_LO_LIM)
                 ohibnd = TM_WORLD_AX(axidx, axcalt, BOX_HI_LIM)
                 IF ( ITSA_TRUEMONTH_AXIS(axcalt) ) THEN
                    olobnd = olobnd * un_convert(pun_day)
                    ohibnd = ohibnd * un_convert(pun_day)
                 ELSE
                    olobnd = olobnd * line_tunit(axcalt)
                    ohibnd = ohibnd * line_tunit(axcalt)
                 ENDIF
                 olobnd = olobnd + origin
                 ohibnd = ohibnd + origin
                 IF ( idx .EQ. 1 ) THEN
                    IF ((tcoord.GE.olobnd).AND.(tcoord.LT.ohibnd)) THEN
*                      First time coordinate matched
                       t1indices(setidx) = axidx
                       tdims(setidx) = line_dim(iaxis)
                    ELSEIF ( tcoord .LT. olobnd ) THEN
*                      Already past first time coordinate - stop the search
                       axidx = axlen
                       GOTO 600
                    ELSE
*                      Keep searching for the first time coordinate
                       GOTO 600
                    ENDIF
                 ELSE
*                   At this time, time coordinates must be consecutive on the created/provided time axis
                    IF ((tcoord.LT.olobnd).OR.(tcoord.GE.ohibnd)) THEN
                       IF ( tmatch .EQ. ' ' ) THEN
*                         @BIN failed; use @ASN
                          GOTO 900
                       ENDIF
                       origin_str = TM_SECS_TO_DATE(tcoord, calid)
                       slen = TM_LENSTR1(origin_str)
                       flen = TM_LENSTR1(memb_filenames(setidx))
                       CALL ERRMSG(ferr_aggregate_error, status,
     .                    'incompatible forecasted-time (T) axis: '//
     .                    'non-consecutive times around "'//
     .                    origin_str(1:slen)//'" in forecast dataset '//
     .                    memb_filenames(setidx)(1:flen), *5000)
                    ENDIF
                 ENDIF
              ENDDO

*             Success with @BIN on this dataset - go on to the next dataset
           ENDIF

        ENDDO

*       Success - all dataset time axes assigned with @XACT or @BIN
        IF ( (tline .GE. 0) .AND. (tmatch .EQ. ' ') ) THEN
           IF ( epsil .GE. 0.0 ) THEN
              CALL WARN('Dataset time axes matched '//
     .                  'to given time axis using @XACT')
           ELSE
              CALL WARN('Dataset time axes matched '//
     .                  'to given time axis using @BIN')
           ENDIF
        ENDIF
        GOTO 1000


*       Blindly associate the datasets in single-step increments onto the given time axis.
  900   DO setidx = 1,nsets
           iaxis = ds_time_axis(member_sets(setidx))
           IF ( iaxis .LE. 0 ) THEN
*             No time axis - assume a single time step
*             (ideally this should work, but may fail elsewhere)
              nsteps = 1
           ELSE
              nsteps = line_dim(iaxis)
           ENDIF
           IF ( (setidx + nsteps) .GT. axlen ) THEN
               flen = TM_LENSTR1(memb_filenames(setidx))
               CALL ERRMSG(ferr_aggregate_error, status,
     .                    'incompatible forecasted-time (T) axis: '//
     .                    'not enough coordinates to associate '//
     .                    'with the time axis of forecast dataset '//
     .                     memb_filenames(setidx)(1:flen), *5000)
           ENDIF
           t1indices(setidx) = setidx
           tdims(setidx) = nsteps
        ENDDO
        IF ( tmatch .EQ. ' ' )
     .     CALL WARN('Dataset time axes matched '//
     .               'to given time axis using @ASN')


* Successful return - match sure tline is assigned as axcalt
 1000   tline = axcalt
        status = ferr_ok
        RETURN

 5000   IF ( (tline .LT. 0) .AND. (axcalt .GE. 0) ) THEN
           CALL FREE_LINE_DYNMEM(axcalt)
           CALL TM_DEALLO_DYN_LINE(axcalt)
        ENDIF
        IF ( wkspc .NE. unspecified_int4 )
     .     CALL FREE_WS_DYNMEM(wkspc)
        RETURN

 5010   idx = TM_LENSTR1(tmatch)
        CALL ERRMSG(ferr_aggregate_error, status,
     .       'Invalid time axis assocation "'
     .       //tmatch(1:idx)//'"', *5000)

        END

