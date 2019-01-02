      SUBROUTINE DEFINE_F_AGGREGATE_F_AXIS(nsets, member_sets,
     .      memb_filenames, fline, fmatch, findices, status)
*
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
* NOAA/PMEL, Seattle, WA - TMAP

* Assign the time-of-forecast (F) axis from the time origins (T0) of the
* time axis of each dataset.  If an F-axis to use is not specified, create
* it from these time origins using the calendar and units of the time axis
* of the first member.  If appropriate, sort the given datasets using the
* time origins.  Returns the sorted datasets, time-of-forecast axis, and
* locations of the datasets on that axis.

* calling arguments
*   nsets          - number of datasets given to be aggregated (input)
*   member_sets    - array of open dataset numbers (input-output)
*   memb_filenames - array of dataset filenames (input-output)
*   fline          - time-of-forecast (F) axis to use/used (input/output)
*   fmatch         - method of matching times to fline (input, can be empty);
*                    eg, '@ASN', '@BIN', '@XACT', '@XACT:0.001'
*   findices       - index of each dataset on the time-of-forecast (F) axis (output)
*   status         - return status (output)

        IMPLICIT NONE

        include 'tmap_dims.parm'
        include 'ferret.parm'
        include 'errmsg.parm'
        include 'xdset_info.cmn_text'
        external xdset_info_data
        include 'xtm_grid.cmn_text'
        include 'xunits.cmn_text'
        include 'netcdf.inc'

* calling argument declarations:
        INTEGER nsets, member_sets(nsets), fline,
     .          findices(nsets), status
        CHARACTER*512 memb_filenames(nsets)
        CHARACTER*(*) fmatch

* local variable declarations:

        REAL*8  GET_LINE_COORD, TM_WORLD_AX, SECS_FROM_BC,
     .          depsilon, epsil, udst, delta, forigin, fcoord,
     .          torigins(nsets), rindices(nsets), lobnd, hibnd
        INTEGER TM_GET_CALENDAR_ID, TM_UNITS_CAL, TM_LENSTR1,
     .          TM_FIND_LIKE_DYN_LINE,
     .          setidx, dset, taxis, flen, slen, calid, ref_taxis,
     .          ref_calid, axcalf, itmp_memb(nsets), axidx, idx,
     .          units, pre_exist, axlen
        CHARACTER*512 itmp_names(nsets)
        CHARACTER TM_SECS_TO_DATE*20,
     .            calname*32, ref_calname*32, buff*256, units_str*64,
     .            torigin_str*20
        LOGICAL ITSA_TRUEMONTH_AXIS,
     .          true_month

        PARAMETER( depsilon = 2.0**(-42) )

*       Do not modify fline until successful - for error handling
        axcalf = fline
        pre_exist = unspecified_int4

*       If blindly associating to a given axis, do not waste time attempting to sort
        IF ( (axcalf .GE. 0) .AND. (fmatch .EQ. '@ASN') ) GOTO 200

*       Extract the T0 time from each dataset for sorting
        DO setidx = 1, nsets
           dset  = member_sets(setidx)
           taxis = ds_time_axis(dset)
           IF ( taxis .LE. 0 ) THEN
              flen = TM_LENSTR1(memb_filenames(setidx))
              CALL ERRMSG(ferr_aggregate_error, status,
     .                    'no time axis in forecast dataset '//
     .                    memb_filenames(setidx)(1:flen), *5000)
           ENDIF
           calname = line_cal_name(taxis)
           calid = TM_GET_CALENDAR_ID(calname)
           IF ( calid .LE. 0 ) THEN
              slen = TM_LENSTR1(calname)
              flen = TM_LENSTR1(memb_filenames(setidx))
              CALL ERRMSG(ferr_aggregate_error, status,
     .                    'unknown calendar "'//calname(1:slen)//
     .                    '" in forecast dataset '//
     .                    memb_filenames(setidx)(1:flen), *5000)
           ENDIF
           IF ( setidx .EQ. 1 ) THEN
*             Use the T axis of the first given dataset as the template when creating the F axis
              ref_taxis = taxis
              ref_calname = calname
              ref_calid = calid
           ENDIF
           torigin_str = line_t0(taxis)
           torigins(setidx) = SECS_FROM_BC(torigin_str, calid, status)
           IF ( status .NE. ferr_ok ) THEN
              slen = TM_LENSTR1(torigin_str)
              flen = TM_LENSTR1(memb_filenames(setidx))
              CALL ERRMSG(ferr_aggregate_error, status,
     .                    'invalid time origin "'//torigin_str(1:slen)//
     .                    '" in forecast dataset '//
     .                    memb_filenames(setidx)(1:flen), *5000)
           ENDIF
        ENDDO

*       Put indices a REAL*8 array for sorting; assumed HEAP2 "REAL" is actually REAL*8
        DO setidx = 1, nsets
           rindices(setidx) = setidx
        ENDDO
        CALL HEAP2(torigins, rindices, unspecified_val8, nsets)

*       Make sure the time origins are unique within a small allowed error (depsilon) 
*       in units of the F axis to be used.  Convert to seconds for SECS_FROM_BC comparisons.
        epsil = depsilon
        IF ( axcalf .LT. 0 ) THEN
*          For now, change units of true-month to units of days
           true_month = ITSA_TRUEMONTH_AXIS(ref_taxis)
           IF ( true_month ) THEN
              units = pun_day
              units_str = un_name(units)
              true_month = .FALSE.
           ELSE
              units = line_unit_code(ref_taxis)
              units_str = line_units(ref_taxis)
           ENDIF
*          Change units of year to the appropriate number of days
           units = TM_UNITS_CAL(units, ref_calid, true_month)
           IF  ( units .GE. 0 ) THEN
              slen = TM_LENSTR1(units_str)
              CALL ERRMSG(ferr_aggregate_error, status,
     .                 'unknown time units "'//units_str(1:slen)//
     .                 '" in the first given forecast dataset', *5000)
           ENDIF
*          units will be used below for defining the F axis
           epsil = epsil * un_convert(units)
        ELSE
           IF ( ITSA_TRUEMONTH_AXIS(axcalf) ) THEN
              epsil = epsil * un_convert(pun_day)
           ELSE
              epsil = epsil * line_tunit(axcalf)
           ENDIF
        ENDIF

        DO setidx = 2, nsets
           IF ( (torigins(setidx)-torigins(setidx-1)) .LE. epsil ) THEN
*             @ASN only valid option
              IF ( (axcalf .GE. 0 ) .AND. (fmatch .EQ. ' ') ) GOTO 200
              dset  = member_sets(NINT(rindices(setidx)))
              taxis = ds_time_axis(dset)
              torigin_str = line_t0(taxis)
              slen = TM_LENSTR1(torigin_str)
              CALL ERRMSG(ferr_aggregate_error, status,
     .                'duplicate time origin "'//torigin_str(1:slen)//
     .                '" in two forecast datasets', *5000)
           ENDIF
        ENDDO

*       Just save the sorted values in the itmp arrays until sure the sorted values are desired
        DO setidx = 1, nsets
           itmp_memb(setidx) = member_sets(NINT(rindices(setidx)))
        ENDDO
        DO setidx = 1, nsets
           itmp_names(setidx) = memb_filenames(NINT(rindices(setidx)))
        ENDDO

        IF ( axcalf .LT. 0 ) THEN
*          Create the F axis from the given time origins.
*          The value of units was computed above for sorting.

*          Create the forecast (F) axis from the given time origins using the
*          calendar and units of the time axis for the first given dataset but
*          starting at the earliest time origin
           CALL TM_ALLO_TMP_LINE(axcalf, status)
           CALL TM_NEW_LINE_NAME('TF_CAL_F', buff)
           line_name(axcalf)         = buff
           line_name_orig(axcalf)    = buff
           line_dattype(axcalf)      = nf_double
           line_modulo(axcalf)       = .FALSE.
           line_modulo_len(axcalf)   = 0.0
           line_shift_origin(axcalf) = .FALSE.
           line_dim_only(axcalf)     = .FALSE.
           line_cal_name(axcalf)     = ref_calname
           line_unit_code(axcalf)    = units
           line_units(axcalf)        = un_name(units)
           line_tunit(axcalf)        = un_convert(units)
           line_direction(axcalf)    = 'FI'
           line_t0(axcalf) = TM_SECS_TO_DATE(torigins(1), ref_calid)
*          start with irregular axis, then check later if regular
           line_regular(axcalf)      = .FALSE.
           CALL GET_LINE_DYNMEM(nsets, axcalf, status)
           IF ( status .NE. ferr_ok )
     .        CALL ERRMSG(ferr_aggregate_error, status,
     .                    'out of memory for forecast axis', *5000)
           udst = un_convert(units)
           DO setidx = 1,nsets
              fcoord = (torigins(setidx) - torigins(1)) / udst
              CALL PUT_LINE_COORD(linemem(axcalf)%ptr, setidx, fcoord)
           ENDDO
           line_dim(axcalf)          = nsets

*          Assign the cell boundaries
           CALL TM_IRREGAX_MID(linemem(axcalf)%ptr,
     .                         lineedg(axcalf)%ptr, nsets)

*          Check if this is a regular axis
           IF ( nsets .GT. 1 ) THEN
              delta = fcoord / (nsets - 1.0)
           ELSE
              delta = 1.0
           ENDIF
*          Now comparing coordinates in the units of the F axis,
*          so use depsilon for the small allowed error.
           DO axidx = 2, nsets-1
              fcoord = GET_LINE_COORD(linemem(axcalf)%ptr, axidx)
              IF ( ABS(fcoord - (axidx-1.0)*delta) .GT. depsilon )
     .           GOTO 100
           ENDDO

*          If we got here, the axis is regular
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
           ELSE
*             Move the tmp_line to a dyn_line
              CALL TM_RE_ALLO_TMP_LINE(axcalf)
           ENDIF

*          Finish reordering the datasets
           DO setidx = 1, nsets
              member_sets(setidx) = itmp_memb(setidx)
           ENDDO
           DO setidx = 1, nsets
              memb_filenames(setidx) = itmp_names(setidx)
           ENDDO

*          Datasets are contiguous from the start on this generated F axis
           DO setidx = 1,nsets
              findices(setidx) = setidx
           ENDDO

*          Successfully created F axis - done
           GOTO 1000

        ENDIF

*       Map the datasets to the provided F axis

  200   axlen = line_dim(axcalf)
        IF ( axlen .LT. nsets )
     .     CALL ERRMSG(ferr_aggregate_error, status,
     .                 'number of dataset exceeds '//
     .                 'the size of the provided F axis', *5000)

*       The allowed difference for an exact match (epsil) is in units of the given time axis
        IF ( fmatch .EQ. '@ASN' ) THEN
           GOTO 900
        ELSEIF ( fmatch .EQ. '@BIN' ) THEN
           epsil = -1.0
        ELSEIF ( fmatch .EQ. '@XACT' ) THEN
           epsil = depsilon
        ELSEIF ( fmatch(1:6) .EQ. '@XACT:' ) THEN
           READ(fmatch(7:), *, ERR=5010) epsil
           IF ( epsil .LT. 0.0 ) GOTO 5010
        ELSEIF ( fmatch .EQ. ' ' ) THEN
           epsil = depsilon
*          Default is to try @XACT; if not successful, try @BIN; if not successful, use @ASN
        ELSE
           GOTO 5010
        ENDIF

        calname = line_cal_name(axcalf)
        calid = TM_GET_CALENDAR_ID(calname)
        IF ( calid .LE. 0 ) THEN
*          @ASN only valid option
           IF ( fmatch .EQ. ' ' ) GOTO 900
           slen = TM_LENSTR1(calname)
           CALL ERRMSG(ferr_aggregate_error, status,
     .                 'unknown calendar "'//calname(1:slen)//
     .                 '" for provided F axis', *5000)
        ENDIF
        torigin_str = line_t0(axcalf)
        forigin = SECS_FROM_BC(torigin_str, calid, status)
        IF ( status .NE. ferr_ok ) THEN
*          @ASN only valid option
           IF ( fmatch .EQ. ' ' ) GOTO 900
           slen = TM_LENSTR1(torigin_str)
           CALL ERRMSG(ferr_aggregate_error, status,
     .                 'invalid time origin "'//torigin_str(1:slen)//
     .                 '" for provided F axis', *5000)
        ENDIF

        IF ( epsil .GE. 0.0 ) THEN
*          Match time origins with a tolerance of epsil

*          Convert epsil from units of given forecast axis to seconds (for SECS_FROM_BC comparisons)
           IF ( ITSA_TRUEMONTH_AXIS(axcalf) ) THEN
              epsil = epsil * un_convert(pun_day)
           ELSE
              epsil = epsil * line_tunit(axcalf)
           ENDIF

           setidx = 1
           DO axidx = 1,axlen
              fcoord = TM_WORLD_AX(axidx, axcalf, BOX_MIDDLE)
              IF ( ITSA_TRUEMONTH_AXIS(axcalf) ) THEN
                 fcoord = fcoord * un_convert(pun_day)
              ELSE
                 fcoord = fcoord * line_tunit(axcalf)
              ENDIF
              fcoord = fcoord + forigin
              IF ( ABS(fcoord - torigins(setidx)) .LE. epsil ) THEN
*                found this t0; go to the next t0, or quit if all found
                 findices(setidx) = axidx
                 setidx = setidx + 1
                 IF ( setidx .GT. nsets ) THEN
*                   @XACT succeeded - finish reordering the datasets and return
                    DO idx = 1, nsets
                       member_sets(idx) = itmp_memb(idx)
                    ENDDO
                    DO idx = 1, nsets
                       memb_filenames(idx) = itmp_names(idx)
                    ENDDO
                    IF ( fmatch .EQ. ' ' )
     .                 CALL WARN('Forecast axis coordinates matched'//
     .                           ' to time axis origins using @XACT')
                    GOTO 1000
                 ENDIF
              ELSEIF ( fcoord .GT. torigins(setidx) ) THEN
*                already past this t0, so not found
                 GOTO 300
              ENDIF
           ENDDO
  300      IF ( fmatch .NE. ' ' ) THEN
*             @XACT specified and failed - error return
              torigin_str = TM_SECS_TO_DATE(torigins(setidx), calid)
              slen = TM_LENSTR1(torigin_str)
              flen = TM_LENSTR1(itmp_names(setidx))
              CALL ERRMSG(ferr_aggregate_error, status,
     .                'incompatible F axis: unmatched time origin "'//
     .                torigin_str(1:slen)//'" of forecast dataset '//
     .                itmp_names(setidx)(1:flen), *5000)
           ENDIF
*          @XACT failed but not specified; try @BIN
        ENDIF


*       Match time origins by cell boundaries
        setidx = 1
        DO axidx = 1,axlen
           lobnd = TM_WORLD_AX(axidx, axcalf, BOX_LO_LIM)
           hibnd = TM_WORLD_AX(axidx, axcalf, BOX_HI_LIM)
           IF ( ITSA_TRUEMONTH_AXIS(axcalf) ) THEN
              lobnd = lobnd * un_convert(pun_day)
              hibnd = hibnd * un_convert(pun_day)
           ELSE
              lobnd = lobnd * line_tunit(axcalf)
              hibnd = hibnd * line_tunit(axcalf)
           ENDIF
           lobnd = lobnd + forigin
           hibnd = hibnd + forigin
           IF ( (lobnd .LE. torigins(setidx)) .AND.
     .          (hibnd .GT. torigins(setidx)) ) THEN
*             found this t0; go to the next t0, or quit if all found
              findices(setidx) = axidx
              setidx = setidx + 1
              IF ( setidx .GT. nsets ) THEN
*                @BIN succeeded - finish reordering the datasets and return
                 DO idx = 1, nsets
                    member_sets(idx) = itmp_memb(idx)
                 ENDDO
                 DO idx = 1, nsets
                    memb_filenames(idx) = itmp_names(idx)
                 ENDDO
                 IF ( fmatch .EQ. ' ' )
     .              CALL WARN('Forecast axis coordinates matched'//
     .                        ' to time axis origins using @BIN')
                 GOTO 1000
              ENDIF
           ELSEIF ( lobnd .GT. torigins(setidx) ) THEN
*             already past this t0, so not found
              GOTO 400
           ENDIF
        ENDDO
  400   IF ( fmatch .NE. ' ' ) THEN
*          @BIN specified and failed
           torigin_str = TM_SECS_TO_DATE(torigins(setidx), calid)
           slen = TM_LENSTR1(torigin_str)
           flen = TM_LENSTR1(itmp_names(setidx))
           CALL ERRMSG(ferr_aggregate_error, status,
     .             'incompatible F axis: unmatched time origin "'//
     .             torigin_str(1:slen)//'" of forecast dataset '//
     .             itmp_names(setidx)(1:flen), *5000)
        ENDIF
*       @BIN failed but not specified; use @ASN


*       Blindly associate the datasets, as given, contiguously
*       onto the given forecast axis.  Do NOT reorder.
  900   DO setidx = 1,nsets
           findices(setidx) = setidx
        ENDDO
        IF ( fmatch .EQ. ' ' )
     .     CALL WARN('Forecast axis coordinates matched '//
     .               'to time axis origins using @ASN')


* Successful return
 1000   fline = axcalf
        status = ferr_ok
        RETURN

* Error returns
 5000   IF ( (fline .LT. 0) .AND. (axcalf .GE. 0) .AND.
     .       (axcalf .NE. pre_exist) ) THEN
           CALL FREE_LINE_DYNMEM(axcalf)
           CALL TM_DEALLO_DYN_LINE(axcalf)
        ENDIF
        RETURN

 5010   slen = TM_LENSTR1(fmatch)
        CALL ERRMSG(ferr_aggregate_error, status,
     .       'Invalid forecast axis assocation "'
     .       //fmatch(1:slen)//'"', *5000)

        END

