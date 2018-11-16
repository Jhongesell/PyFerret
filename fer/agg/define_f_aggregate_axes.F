        SUBROUTINE DEFINE_F_AGGREGATE_AXES(nsets, member_sets,
     .     torigins, fline, tline, findices, t1indices, tdims, status)

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
* Assign the time-of-forecast (F) axis from the time origins.  If an
* F-axis to use is not specified, create it from these time origins
* using the calendar and units of the time axis of the first member.
*
* Assign the forecasted-time (T) axis from the time axes or member
* sets.  If a T-axis to use is not specified, create it from merging
* all the (hopefully overlapping) times given.
*
* An error is raised if any time does not match a time on the appropriate
* time-of-forecast (F) or forecasted-time (T) axis.
* At this time, an error is also raised if the forecasted times for a forecast
* dataset are not consecutive times on the forecasted-time (T) axis.
*
* arguments:
*       nsets       - number of forecast datasets (input)
*       member_sets - dataset numbers for each forecast dataset (input)
*       torigins    - time axis origin (in SECS_FROM_BC) for each forecast dataset (input)
*       fline       - time-of-forecast (F) axis to use/used (input/output)
*       tline       - forecasted-times (T) axis to use/used (input/output)
*       findices    - index of each time-axis-origin on the time-of-forecast (F) axis (output)
*       t1indices   - index of each first-forecasted-time on the forecasted-time (T) axis (output)
*       tdims       - number of forecasted times for each forecast dataset (output)
*       status      - return status (output)

        IMPLICIT NONE

        include 'tmap_dims.parm'
        include 'ferret.parm'
        include 'xdset_info.cmn_text'
        external xdset_info_data
        include 'xdset_info.cd_equiv'
        include 'xtm_grid.cmn_text'
        include 'xunits.cmn_text'
        include 'errmsg.parm'
        include 'netcdf.inc'

* calling argument declarations:
        INTEGER nsets, member_sets(nsets), fline, tline,
     .          findices(nsets), t1indices(nsets), tdims(nsets), status
        REAL*8  torigins(nsets)

* local variable declarations:
        REAL*8  GET_LINE_COORD, TM_WORLD_AX, SECS_FROM_BC,
     .          udst, tcoord, delta, origin, ocoord, setorigin
        INTEGER TM_GET_CALENDAR_ID, TM_UNITS_CAL, TM_LENSTR1,
     .          TM_FIND_LIKE_DYN_LINE,
     .          axcalf, axcalt, ref_grid, ref_tax, calid, units,
     .          setidx, pre_exist, axlen, axidx,
     .          nsteps, ntsteps, varidx, igrid, iaxis, idx, previdx
        CHARACTER TM_SECS_TO_DATE*20,
     .          calname*32, buff*256, units_str*64, origin_str*20
        LOGICAL ITSA_TRUEMONTH_AXIS, TM_DFPEQ, true_month

        axcalf = fline
        axcalt = tline

*       Find the indices of the time origins on the forecast (F) axis
*       (creating it if not given)
        IF ( axcalf .LE. 0 ) THEN

*          get the time axis information for the first member set
           ref_tax = ds_time_axis(member_sets(1))
           IF ( ref_tax .LE. 0 )
     .        CALL ERRMSG(ferr_aggregate_error, status,
     .                    'no time axis for forecast dataset '
     .                    //ds_name(member_sets(1)), *5000)
           calname = line_cal_name(ref_tax)
           calid = TM_GET_CALENDAR_ID(calname)
           IF ( calid .LE. 0 )
     .        CALL ERRMSG(ferr_aggregate_error, status,
     .                    'unknown calendar "'
     .                    //calname(1:TM_LENSTR1(calname))
     .                    //'" in forecast dataset '
     .                    //ds_name(member_sets(1)), *5000)
*          for now, change units of true-month to units of days
           true_month = ITSA_TRUEMONTH_AXIS(ref_tax)
           IF ( true_month ) THEN
              units = pun_day
              units_str = un_name(units)
              true_month = .FALSE.
           ELSE
              units = line_unit_code(ref_tax)
              units_str = line_units(ref_tax)
           ENDIF
*          change units of year to the appropriate number of days
           units = TM_UNITS_CAL(units, calid, true_month)
           IF  ( units .GE. 0 )
     .        CALL ERRMSG(ferr_aggregate_error, status,
     .                    'unknown time units "'
     .                    //units_str(1:TM_LENSTR1(units_str))
     .                    //'" in a forecast dataset', *5000)

*          Create the forecast (F) axis from the given time origins using
*          the calendar and units of the time axis for the first member set
           CALL TM_ALLO_TMP_LINE(axcalf, status)
           CALL TM_NEW_LINE_NAME('TF_CAL_F', buff)
           line_name(axcalf)         = buff
           line_name_orig(axcalf)    = buff
           line_dattype(axcalf)      = nf_double
           line_modulo(axcalf)       = .FALSE.
           line_modulo_len(axcalf)   = 0.0
           line_shift_origin(axcalf) = .FALSE.
           line_dim_only(axcalf)     = .FALSE.
           line_cal_name(axcalf)     = calname
           line_unit_code(axcalf)    = units
           line_units(axcalf)        = un_name(units)
           line_tunit(axcalf)        = un_convert(units)
           line_direction(axcalf)    = 'FI'
           line_t0(axcalf) = TM_SECS_TO_DATE(torigins(1), calid)
*          start with irregular axis, then check later if regular
           line_regular(axcalf)      = .FALSE.
           CALL GET_LINE_DYNMEM(nsets, axcalf, status)
           IF ( status .NE. ferr_ok )
     .        CALL ERRMSG(ferr_aggregate_error, status,
     .                    'out of memory for forecast axis', *5000)
           udst = un_convert(units)
           DO setidx = 1,nsets
              tcoord = (torigins(setidx) - torigins(1)) / udst
              CALL PUT_LINE_COORD(linemem(axcalf)%ptr, setidx, tcoord)
           ENDDO
           line_dim(axcalf)          = nsets
*          assign the cell boundaries
           CALL TM_IRREGAX_MID(linemem(axcalf)%ptr,
     .                         lineedg(axcalf)%ptr, nsets)

*          check if this is a regular axis
           IF ( nsets .GT. 1 ) THEN
              delta = tcoord / (nsets - 1.0)
           ELSE
              delta = 1.0
           ENDIF
           DO axidx = 2, nsets-1
              tcoord = GET_LINE_COORD(linemem(axcalf)%ptr, axidx)
              IF ( .NOT. TM_DFPEQ(tcoord, (axidx - 1.0) * delta) )
     .           GOTO 100
           ENDDO
*          if we got here, the axis is regular
           CALL FREE_LINE_DYNMEM(axcalf)
           line_regular(axcalf) = .TRUE.
           line_start(axcalf)   = 0.0
           line_delta(axcalf)   = delta

*          If the axis is an exact duplicate of a pre-existing axis, use that axis
*          instead (shared axes are needed in order to E-aggregate multiple FMRCs)
  100      pre_exist = TM_FIND_LIKE_DYN_LINE(axcalf)
           IF ( pre_exist .NE. unspecified_int4 ) THEN
              IF ( .NOT. line_regular(axcalf) )
     .           CALL FREE_LINE_DYNMEM(axcalf)
              axcalf = pre_exist
*             also set fline so it will not be deleted if an error occurs later
              fline = pre_exist
           ELSE
*             move the tmp_line to a dyn_line
              CALL TM_RE_ALLO_TMP_LINE(axcalf)
           ENDIF

*          forecast axis is exactly the time origins
           DO setidx = 1,nsets
              findices(setidx) = setidx
           ENDDO

        ELSE

*          Map the forecast time origins to the provided F axis
           axlen = line_dim(axcalf)
           IF ( axlen .LT. nsets )
     .        CALL ERRMSG(ferr_aggregate_error, status,
     .             'not enough time steps in provided F axis', *5000)
           IF ( line_direction(axcalf) .NE. 'FI' )
     .        CALL ERRMSG(ferr_aggregate_error, status,
     .             'provided F axis not a forecast time axis', *5000)
           calname = line_cal_name(axcalf)
           calid = TM_GET_CALENDAR_ID(calname)
           IF ( calid .LE. 0 )
     .        CALL ERRMSG(ferr_aggregate_error, status,
     .                    'unknown calendar "'
     .                    //calname(1:TM_LENSTR1(calname))
     .                    //'" for provided F axis', *5000)
           origin_str = line_t0(axcalf)
           origin = SECS_FROM_BC(origin_str, calid, status)
           IF ( status .NE. ferr_ok )
     .        CALL ERRMSG(ferr_aggregate_error, status,
     .                    'invalid time origin "'
     .                    //origin_str(1:TM_LENSTR1(origin_str))
     .                    //'" for provided F axis', *5000)
           setidx = 1
           DO axidx = 1,axlen
              tcoord = TM_WORLD_AX(axidx, axcalf, BOX_MIDDLE)
              IF ( ITSA_TRUEMONTH_AXIS(axcalf) ) THEN
                 tcoord = tcoord * un_convert(pun_day)
              ELSE
                 tcoord = tcoord * line_tunit(axcalf)
              ENDIF
              tcoord = tcoord + origin
              IF ( TM_DFPEQ(tcoord, torigins(setidx)) ) THEN
*                found this t0; go to the next t0, or quit if all found
                 findices(setidx) = axidx
                 setidx = setidx + 1
                 IF ( setidx .GT. nsets ) GOTO 200
              ELSEIF ( tcoord .GT. torigins(setidx) ) THEN
*                already past this t0, so not found
                 GOTO 200
              ENDIF
           ENDDO
  200      IF ( setidx .LE. nsets ) THEN
              origin_str = TM_SECS_TO_DATE(torigins(setidx), calid)
              CALL ERRMSG(ferr_aggregate_error, status,
     .                    'incompatible F axis: unmatched time origin "'
     .                    //origin_str(1:TM_LENSTR1(origin_str))
     .                    //'" of forecast dataset '
     .                    //ds_name(member_sets(setidx)), *5000)
           ENDIF

        ENDIF


*       Find the indices of the forecasted times for each forecast
*       on the forecasted-time (T) axis (creating it if not given)
        IF ( tline .LE. 0 ) THEN

*          get the time axis for the first member set
           ref_tax = ds_time_axis(member_sets(1))
           IF ( ref_tax .LE. 0 )
     .         CALL ERRMSG(ferr_aggregate_error, status,
     .                    'no time axis for forecast dataset '
     .                    //ds_name(member_sets(1)), *5000)
           calname = line_cal_name(ref_tax)
           calid = TM_GET_CALENDAR_ID(calname)
           IF ( calid .LE. 0 )
     .        CALL ERRMSG(ferr_aggregate_error, status,
     .                    'unknown calendar "'
     .                    //calname(1:TM_LENSTR1(calname))
     .                    //'" in forecast dataset '
     .                    //ds_name(member_sets(1)), *5000)
           origin_str = line_t0(ref_tax)
           origin = SECS_FROM_BC(origin_str, calid, status)
           IF ( status .NE. ferr_ok )
     .        CALL ERRMSG(ferr_aggregate_error, status,
     .                    'invalid time origin "'
     .                    //origin_str(1:TM_LENSTR1(origin_str))
     .                    //'" in forecast dataset '
     .                    //ds_name(member_sets(1)), *5000)
*          for now, change units of true-month to units of days
           true_month = ITSA_TRUEMONTH_AXIS(ref_tax)
           IF ( true_month ) THEN
              units = pun_day
              units_str = un_name(units)
              true_month = .FALSE.
           ELSE
              units = line_unit_code(ref_tax)
              units_str = line_units(ref_tax)
           ENDIF
*          change units of year to the appropriate number of days
           units = TM_UNITS_CAL(units, calid, true_month)
           IF  ( units .GE. 0 )
     .        CALL ERRMSG(ferr_aggregate_error, status,
     .                    'unknown time units "'
     .                    //units_str(1:TM_LENSTR1(units_str))
     .                    //'" in forecast dataset '
     .                    //ds_name(member_sets(1)), *5000)

*          assume the time axis of each forecast is the same length and
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
     .                    'out of memory for forecast axis', *5000)

*          add all unique time coordinates from all forecast sets
           ntsteps = 0
           DO setidx = 1,nsets
*             get the time axis information for this forecast dataset
              iaxis = ds_time_axis(member_sets(setidx))
              IF ( iaxis .LE. 0 )
     .           CALL ERRMSG(ferr_aggregate_error, status,
     .                       'No valid time axis for forecast dataset '
     .                       //ds_name(member_sets(setidx)), *5000)
              calname = line_cal_name(iaxis)
              calid = TM_GET_CALENDAR_ID(calname)
              IF ( calid .LE. 0 )
     .           CALL ERRMSG(ferr_aggregate_error, status,
     .                       'unknown calendar "'
     .                       //calname(1:TM_LENSTR1(calname))
     .                       //'" for forecast dataset '
     .                       //ds_name(member_sets(setidx)), *5000)
              origin_str = line_t0(iaxis)
              setorigin = SECS_FROM_BC(origin_str, calid, status)
              IF ( status .NE. ferr_ok )
     .            CALL ERRMSG(ferr_aggregate_error, status,
     .                        'invalid time origin "'//
     .                        origin_str(1:TM_LENSTR1(origin_str))
     .                       //'" for forecast dataset '
     .                       //ds_name(member_sets(setidx)), *5000)

              axidx = 1
              DO idx = 1,line_dim(iaxis)
*                get this forecasted time (in units of SECS_FROM_BC)
                 tcoord = TM_WORLD_AX(idx, iaxis, BOX_MIDDLE)
                 IF ( ITSA_TRUEMONTH_AXIS(iaxis) ) THEN
                    tcoord = tcoord * un_convert(pun_day)
                 ELSE
                    tcoord = tcoord * line_tunit(iaxis)
                 ENDIF
                 tcoord = tcoord + setorigin

*                match or add this forecasted time (in units of SECS_FROM_BC)
                 DO WHILE ( axidx .LE. ntsteps )
                    ocoord = GET_LINE_COORD(linemem(axcalt)%ptr, axidx)
*                   if matches an existing value, no changes needed
                    IF ( TM_DFPEQ(ocoord, tcoord) ) GOTO 340
*                   if new value is smaller, insert it here
                    IF ( ocoord .GT. tcoord ) GOTO 320
                    axidx = axidx + 1
                 ENDDO
*                insert the value, moving any value after it
  320            IF ( ntsteps .EQ. nsets*nsteps )
     .              CALL ERRMSG(ferr_aggregate_error, status,
     .                          'unexpected large number of unique '
     .                          //'forecasted times', *5000)
                 ntsteps = ntsteps + 1
*                insert/append the new value
                 previdx = axidx
                 DO WHILE ( previdx .LE. ntsteps )
                    CALL PUT_LINE_COORD(linemem(axcalt)%ptr,
     .                                  previdx, tcoord)
                    tcoord = ocoord
                    previdx = previdx + 1
                    IF ( previdx .LT. ntsteps ) ocoord =
     .                 GET_LINE_COORD(linemem(axcalt)%ptr, previdx)
                 ENDDO
  340            CONTINUE
              ENDDO
           ENDDO

  360      line_dim(axcalt) = ntsteps
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

           CALL TM_IRREGAX_MID(linemem(axcalt)%ptr,
     .                         lineedg(axcalt)%ptr, ntsteps)

*          check if this is a regular axis
           ocoord = GET_LINE_COORD(linemem(axcalt)%ptr, 1)
           IF ( ntsteps .GT. 1 ) THEN
              tcoord = GET_LINE_COORD(linemem(axcalt)%ptr, ntsteps)
              delta = (tcoord - ocoord) / (ntsteps - 1.0)
           ELSE
              delta = 1.0
           ENDIF
           DO axidx = 2,ntsteps-1
              tcoord = GET_LINE_COORD(linemem(axcalt)%ptr, axidx)
              IF ( .NOT. TM_DFPEQ(tcoord, ocoord + (axidx-1.0)*delta) )
     .           GOTO 380
           ENDDO

*          if we got here, the axis is regular
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
*             also set tline so it will not be deleted if an error occurs later
              tline = pre_exist
           ELSE
*             move the tmp_line to a dyn_line
              CALL TM_RE_ALLO_TMP_LINE(axcalt)
           ENDIF

        ENDIF

*       Map the forecasted times to the created/provided T axis
*       (could have been done while creating the T axis, but easier this way)
        axlen = line_dim(axcalt)
        IF ( axlen .LT. 1 )
     .     CALL ERRMSG(ferr_aggregate_error, status,
     .          'invalid provided T axis', *5000)
        IF ( line_direction(axcalt) .NE. 'TI' )
     .     CALL ERRMSG(ferr_aggregate_error, status,
     .          'provided T axis not a time axis', *5000)
        calname = line_cal_name(axcalt)
        calid = TM_GET_CALENDAR_ID(calname)
        IF ( calid .LE. 0 )
     .     CALL ERRMSG(ferr_aggregate_error, status,
     .                 'unknown calendar "'
     .                 //calname(1:TM_LENSTR1(calname))
     .                 //'" for provided T axis', *5000)
        origin_str = line_t0(axcalt)
        origin = SECS_FROM_BC(origin_str, calid, status)
        IF ( status .NE. ferr_ok )
     .     CALL ERRMSG(ferr_aggregate_error, status,
     .                 'invalid time origin "'
     .                 //origin_str(1:TM_LENSTR1(origin_str))
     .                 //'" for provided T axis', *5000)

*       check time coordinates from all forecast sets
        DO setidx = 1,nsets
*          get the time axis information for this forecast dataset
           iaxis = ds_time_axis(member_sets(setidx))
           IF ( iaxis .LE. 0 )
     .        CALL ERRMSG(ferr_aggregate_error, status,
     .                    'No valid time axis for forecast dataset '
     .                    //ds_name(member_sets(setidx)), *5000)
           calname = line_cal_name(iaxis)
           calid = TM_GET_CALENDAR_ID(calname)
           IF ( calid .LE. 0 )
     .        CALL ERRMSG(ferr_aggregate_error, status,
     .                    'unknown calendar "'
     .                    //calname(1:TM_LENSTR1(calname))
     .                    //'" for forecast dataset '
     .                    //ds_name(member_sets(setidx)), *5000)
           origin_str = line_t0(iaxis)
           setorigin = SECS_FROM_BC(origin_str, calid, status)
           IF ( status .NE. ferr_ok )
     .            CALL ERRMSG(ferr_aggregate_error, status,
     .                        'invalid time origin "'//
     .                        origin_str(1:TM_LENSTR1(origin_str))
     .                       //'" for forecast dataset '
     .                       //ds_name(member_sets(setidx)), *5000)

           axidx = 0
           DO idx = 1,line_dim(iaxis)
*             get this forecasted time (in SECS_FROM_BC)
              tcoord = TM_WORLD_AX(idx, iaxis, BOX_MIDDLE)
              IF ( ITSA_TRUEMONTH_AXIS(iaxis) ) THEN
                 tcoord = tcoord * un_convert(pun_day)
              ELSE
                 tcoord = tcoord * line_tunit(iaxis)
              ENDIF
              tcoord = tcoord + setorigin

*             match to a time on the created/provided T axis
  400         axidx = axidx + 1
              IF ( axidx .GT. axlen )
     .           CALL ERRMSG(ferr_aggregate_error, status,
     .                       'incompatible provided T axis: '
     .                       //'unmatched time of forecast dataset '
     .                       //ds_name(member_sets(setidx)), *5000)
              ocoord = TM_WORLD_AX(axidx, axcalt, BOX_MIDDLE)
              IF ( ITSA_TRUEMONTH_AXIS(axcalt) ) THEN
                 ocoord = ocoord * un_convert(pun_day)
              ELSE
                 ocoord = ocoord * line_tunit(axcalt)
              ENDIF
              ocoord = ocoord + origin
              IF ( idx .EQ. 1 ) THEN
*                if not a match, keep searching for the first time coordinate
                 IF ( .NOT. TM_DFPEQ(ocoord, tcoord) ) GOTO 400
                 t1indices(setidx) = axidx
                 tdims(setidx) = line_dim(iaxis)
              ELSE
*                at this time, time coordinates must be consecutive on the created/provided time axis
                 IF ( .NOT. TM_DFPEQ(ocoord, tcoord) )
     .              CALL ERRMSG(ferr_aggregate_error, status,
     .                   'incompatible forecasted time axis: '
     .                   //'non-consecutive times for forecast dataset '
     .                   //ds_name(member_sets(setidx)), *5000)
              ENDIF
           ENDDO
        ENDDO

*       Return the T and F axes used
        tline = axcalt
        fline = axcalf
        RETURN

 5000   CONTINUE
        IF ( (tline .LE. 0) .AND. (axcalt .GT. 0) ) THEN
           CALL FREE_LINE_DYNMEM(axcalt)
           CALL TM_DEALLO_DYN_LINE(axcalt)
        ENDIF
        IF ( (fline .LE. 0) .AND. (axcalf .GT. 0) ) THEN
           CALL FREE_LINE_DYNMEM(axcalf)
           CALL TM_DEALLO_DYN_LINE(axcalf)
        ENDIF
        RETURN

        END

