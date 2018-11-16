        SUBROUTINE ASSIGN_T_AGG_AXIS(given_taxis, dset_num, nagfiles,
     .             memb_filename, memb_tstep_start, memb_tstep_end,
     .             memb_nstep_start, memb_nsteps, memb_order, status)

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
*  ARE DISCLAIMED. IN NO EVENT SHALL NOAA/PMEL/TMAP BE LIABLE FOR ANY
*  SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER
*  RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF
*  CONTRACT, NEGLIGENCE OR OTHER TORTUOUS ACTION, ARISING OUT OF OR IN
*  CONNECTION WITH THE ACCESS, USE OR PERFORMANCE OF THIS SOFTWARE.
*

*  Use the given time axis, after validating, as the aggregation axis.
*  Fill in the file time metadata and assign the ordering of the member
*  files according to this time axis.

        IMPLICIT NONE
        include 'tmap_errors.parm'
#include "tmap_dset.parm"
        include 'tmap_dims.parm'
        include 'ferret.parm'
        include 'xdset_info.cmn_text'
        external xdset_info_data
        include 'xdset_info.cd_equiv'
        include 'xprog_state.cmn'
        include 'xtm_grid.cmn_text'
        include 'xunits.cmn_text'
        include 'errmsg.parm'
        include 'xvariables.cmn'
        include 'xcontext.cmn'
        include 'xfr_grid.cmn'
        include 'interp_stack.parm'
        include 'calendar.decl'
        include 'calendar.cmn'
        include 'netcdf.inc'


* arguments:

*       given_taxis      - time axis number to use (input)
*       dset_num         - dataset number 1st file -- to become the T agg (input)
*       nagfiles         - number of files in the aggregation (input)
*       memb_filename    - array of filenames (input)
*       memb_tstep_start - array of first timestep in each file (output)
*       memb_tstep_end   - array of last timestep in each file (output)
*       memb_nstep_start - array of time axis indices of the first step (output)
*       memb_nsteps      - array of number of timesteps in each file (output)
*       memb_order       - sorted ordering index of files (output)
*       status           - return status (output)
*

* Argument definitions
        INTEGER given_taxis, dset_num, nagfiles, memb_nstep_start(nagfiles),
     .          memb_nsteps(nagfiles), memb_order(nagfiles), status
        REAL*8  memb_tstep_start(nagfiles), memb_tstep_end(nagfiles)
        CHARACTER*512 memb_filename(nagfiles)


* local variable declarations:
        LOGICAL TM_DFPEQ_TOL, need_T0_adjust
        INTEGER TM_LENSTR1, TM_GET_CALENDAR_ID,
     .          taxis, agg_nsteps, nsteps, dtype, vlen,
     .          i, ifile, istep, ivar, iline,
     .          cdfid, vid, bounds_id, gg,
     .          agg_unit_code, file_unit_code, agg_cal_id
        REAL*8  TM_WORLD_AX, T0_ADJUST, GET_LINE_COORD,
     .          depsilon, aggstep, tstep, memb_sort(nagfiles)
        CHARACTER LEFINT*20, vname*64, units_str*128, 
     .            calendar_str*128, agg_t0_str*20, file_t0_str*20

        PARAMETER( depsilon = 2.0**(-42) )

*       Get the time axis name from the first file 
*       (assumes rest of the files use the same name for the axis)
        taxis          = ds_time_axis(dset_num)
        IF (taxis .EQ. mnormal) GOTO 5020  ! no time axis?
        vname          = line_name_orig(taxis)
        vlen           = TM_LENSTR1(vname)

        agg_nsteps    = line_dim(given_taxis)
        dtype         = unspecified_int4   ! irrelevant from older netCD
        agg_unit_code = line_unit_code(given_taxis)
        agg_T0_str    = line_T0(given_taxis)
        agg_cal_id    = TM_GET_CALENDAR_ID(line_cal_name(given_taxis))

*       Since putting this on a given time axis, treat the first file same as others
        DO ifile = 1, nagfiles
* ... first read the length of the time axis from the file to be aggregated
           CALL CD_QUIK_NCOPEN(dset_num, memb_filename(ifile), 
     .                    vname(:vlen), cdfid, vid, bounds_id, 
     .                    units_str, calendar_str, nsteps, status)
           IF (status .NE. merr_ok) GOTO 5010
* ... check to see if calendar, units, and T0 are consistent with the agg axis
           CALL CD_GET_T0(memb_filename(ifile), units_str, calendar_str,
     .                    agg_cal_id, file_unit_code, file_t0_str, status)
           IF (status .NE. merr_ok) GOTO 5010
           need_T0_adjust = (file_unit_code .NE. agg_unit_code)
     .                 .OR. (file_T0_str .NE. agg_T0_str)
           IF (need_t0_adjust) THEN
              CALL T0_ADJUST_INIT(agg_cal_id, file_unit_code, 
     .                 file_t0_str, agg_unit_code, agg_t0_str, status)
              IF (status .NE. ferr_ok) GOTO 5020
           ENDIF

* grab a temporary dynamic axis slot for it
           CALL TM_ALLO_TMP_LINE(iline, status)
           CALL GET_LINE_DYNMEM(nsteps, iline, status)
           IF (status .NE. ferr_ok) GOTO 5010

* ... read the array of time steps from the file
           CALL CD_RD_R8_ARR(cdfid, vid, 1, nsteps, dtype,
     .              vname(:vlen), linemem(iline)%ptr, status)
           IF ( status .NE. merr_ok ) GOTO 5010

* ... adjust the time step values to the agg-chosen "<units> since <T0>"
           IF ( need_T0_adjust ) THEN
              DO istep = 1, nsteps
                 tstep = GET_LINE_COORD(linemem(iline)%ptr, istep)
                 CALL PUT_LINE_COORD(linemem(iline)%ptr, istep, 
     .                               T0_ADJUST(tstep))
              ENDDO
           ENDIF

* Verify all the time steps for this file are found and are contiguous on the given axis
           status = ferr_out_of_range
           tstep = GET_LINE_COORD(linemem(iline)%ptr, 1)
           IF ( line_regular(given_taxis) ) THEN
              DO i = 1, agg_nsteps
                 aggstep = line_start(given_taxis) + (i-1)*line_delta(given_taxis)
                 IF ( ABS(tstep - aggstep) .LT. depsilon ) THEN
                    memb_nstep_start(ifile) = i
                    GOTO 100
                 ENDIF
              ENDDO
              GOTO 5020
 100          DO istep = 2, nsteps
                 tstep = GET_LINE_COORD(linemem(iline)%ptr, istep)
                 i = i + 1
                 IF ( i .GT. agg_nsteps ) GOTO 5020
                 aggstep = line_start(given_taxis) + (i-1)*line_delta(given_taxis)
                 IF ( ABS(tstep - aggstep) .GE. depsilon ) GOTO 5020
              ENDDO
           ELSE
              DO i = 1, agg_nsteps
                 aggstep = GET_LINE_COORD(linemem(given_taxis)%ptr, i)
                 IF ( ABS(tstep - aggstep) .LT. depsilon ) THEN
                    memb_nstep_start(ifile) = i
                    GOTO 200
                 ENDIF
              ENDDO
              GOTO 5020
 200          DO istep = 2, nsteps
                 tstep = GET_LINE_COORD(linemem(iline)%ptr, istep)
                 i = i + 1
                 IF ( i .GT. agg_nsteps ) GOTO 5020
                 aggstep = GET_LINE_COORD(linemem(given_taxis)%ptr, i)
                 IF ( ABS(tstep - aggstep) .GE. depsilon ) GOTO 5020
              ENDDO
           ENDIF
           status = ferr_ok

* Assign the time steps for this file
           memb_tstep_start(ifile) = 
     .        GET_LINE_COORD(linemem(iline)%ptr, 1)
           memb_tstep_end(ifile) = 
     .        GET_LINE_COORD(linemem(iline)%ptr, nsteps)
           memb_nsteps(ifile) = nsteps
        ENDDO

* Get the ordering of the files according to the start time
        DO i = 1, nagfiles
           memb_sort(i) = memb_tstep_start(i)
        ENDDO
        CALL SORT_TSERIES_DSETS(nagfiles, memb_order, memb_sort)
* Check that the files do not overlap in time
        DO i = 2, nagfiles
           IF ( memb_tstep_start(memb_order(i)) .LE.
     .          memb_tstep_end(memb_order(i-1)) ) GOTO 5030
        ENDDO

* correct all time-dependent variables in this aggregate dataset
        ds_time_axis(dset_num) = given_taxis
        DO ivar = 1, maxvars
           IF ( ds_var_setnum(ivar) .EQ. dset_num ) THEN
              gg = ds_grid_number(ivar)
              IF ( grid_line(t_dim,gg) .EQ. taxis .OR.
     .             grid_line(t_dim,gg) .EQ. given_taxis ) THEN
                 ds_grid_end(t_dim,ivar) = agg_nsteps
              ENDIF
              IF ( grid_line(t_dim,gg) .EQ. taxis ) THEN
* ... was this grid modified in this loop an instant ago?
*     note: Should use the dset var LIST!   Super-inefficient nested loop
                 DO i = 1, ivar-1
                    IF ( ds_var_setnum(i) .EQ. dset_num ) THEN
                       IF (ds_grid_number(i) .EQ. gg) GOTO 500
!already done
                    ENDIF
                 ENDDO
                 grid_line(t_dim,gg) = given_taxis
* ... "use" the given time axis; "un-use" the one from the file
                 CALL TM_USE_LINE(given_taxis)
                 CALL TM_DEALLO_DYN_LINE(taxis)
              ENDIF
           ENDIF
 500       CONTINUE
        ENDDO
        IF (.NOT.line_regular(taxis)) CALL FREE_LINE_DYNMEM(taxis)

 5000   RETURN
 5010   CALL ERRMSG(ferr_TMAP_error, status, ' ', *5000)
 5020   CALL ERRMSG(ferr_grid_definition, status,
     .              'invalid time axis for the aggregation', *5000)
 5030   CALL ERRMSG(ferr_grid_definition, status,
     .   'Time axes overlap for files '//memb_filename(memb_order(i-1))
     .   //' and '//memb_filename(memb_order(i)), *5000)

        END
