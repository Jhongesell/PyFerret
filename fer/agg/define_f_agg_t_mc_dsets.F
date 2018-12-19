       SUBROUTINE DEFINE_F_AGG_T_MC_DSETS(nsets, memb_sets,
     .            memb_dsetnames, axcalt, t1indices, tdims, status)
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
*  ARE DISCLAIMED. IN NO EVENT SHALL NOAA/PMEL/TMAP BE LIABLE FOR ANY
*  SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER
*  RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF
*  CONTRACT, NEGLIGENCE OR OTHER TORTUOUS ACTION, ARISING OUT OF OR IN
*  CONNECTION WITH THE ACCESS, USE OR PERFORMANCE OF THIS SOFTWARE.
*
*
* Create forecast datasets corresponding to the given forecast dataset
* except defined for the full forecasted time axis axcalt.  Do this by
* create "multifile" netCDF ("MC") T-aggregated datasets for axcalt but
* with only the single file to the appropriate section of the time axis.

* arguments:
*       nsets         - number of forecast datasets (input)
*       memb_sets     - forecast dataset numbers (input)
*       memb_dsetname - forecast dataset filenames (input)
*       axcalt        - forecasted time (T) axis to use for all forecast datasets (input)
*       t1indices     - index on axcalt of the first forecasted time for each forecast dataset (input)
*       tdims         - number of forecasted times for each forecast dataset (input)
*       status        - return status (output)

        IMPLICIT NONE

* Argument definitions
        INTEGER nsets, memb_sets(nsets), axcalt,
     .          t1indices(nsets), tdims(nsets), status
        CHARACTER*512 memb_dsetnames(nsets)

#include "tmap_dset.parm"
        include 'tmap_dims.parm'
        include 'ferret.parm'
        include 'errmsg.parm'
        include 'xtm_grid.cmn_text'
        external xgt_grid_data
        include 'xdset_info.cmn_text'
        external xdset_info_data
        include 'xstep_files.cmn_text'
        external xstep_files_data
        include 'xdset_info.cd_equiv'
        include 'descript.def'

* Internal declarations
        INTEGER ii, dset, taxis, istep, ivar, gg, idx
        REAL*8  GET_LINE_COORD

        DO ii = 1, nsets
           dset = memb_sets(ii)
           taxis = ds_time_axis(dset)

*          find the (already open) dataset in the stepfile records
           istep = 1
           DO WHILE ( sf_setnum(istep) .NE. dset )
              istep = istep + 1
              IF ( istep .GT. maxstepfiles )
     .           CALL ERRMSG(ferr_aggregate_error, status,
     .              'Unexpected failure to find dataset in stepfiles',
     .              *5000)
           ENDDO

*          close the actual dataset - will become the stepfile dataset
           CALL CD_NCCLOSE(dset, status)
           IF ( status .NE. ferr_ok )
     .        CALL ERRMSG(ferr_TMAP_error, status, ' ', *5000)

*          assign the one file (the original dataset) for this stepfile
           sf_name(istep)      = memb_dsetnames(ii)
           sf_ndxs(istep)      = t1indices(ii)
           sf_ndxe(istep)      = t1indices(ii) + tdims(ii) - 1
           sf_regtimes(istep)  = line_regular(axcalt)
           IF ( line_regular(axcalt) ) THEN
              sf_delta(istep)  = line_delta(axcalt)
              sf_1step(istep)  = line_start(axcalt) +
     .              (sf_ndxs(istep) - 1) * line_delta(axcalt)
              sf_lastep(istep) = line_start(axcalt) +
     .              (sf_ndxe(istep) - 1) * line_delta(axcalt)
           ELSE
              sf_delta(istep)  = 0.0
              sf_1step(istep)  =
     .              GET_LINE_COORD(linemem(axcalt)%ptr, sf_ndxs(istep))
              sf_lastep(istep) =
     .              GET_LINE_COORD(linemem(axcalt)%ptr, sf_ndxe(istep))
           ENDIF
*          dataset is not open
           sf_valid(istep)     = .FALSE.
           sf_regvars(istep)   = .TRUE.
           sf_lunit(istep)     = file_not_open
           sf_index(istep)     = ' '
           sf_reclen(istep)    = unspecified_int4

*          associate the dataset number (which was closed) with the stepfile on all of axcalt
           mc_init_stepfile(dset) = istep
           ds_open(dset)       = istep
           ds_type(dset)       = '  MC'
           ds_name(dset)       = memb_dsetnames(ii)
           ds_des_name(dset)   = ds_name(dset)
           ds_title(dset)      =
     .           'T-axis extended '//memb_dsetnames(ii)
           ds_mod_title(dset)  = ds_title(dset)
           ds_time_axis(dset)  = axcalt
           ds_t0time(dset)     = line_t0(axcalt)
           ds_time_unit(dset)  = line_tunit(axcalt)
           ds_cal_name(dset)   = line_cal_name(axcalt)
           ds_regsteps(dset)   = line_regular(axcalt)
           IF ( line_regular(axcalt) ) THEN
             ds_delstep(dset)  = line_delta(axcalt)
             ds_1step(dset)    = line_start(axcalt)
             ds_lastep(dset)   = line_start(axcalt) +
     .              (line_dim(axcalt) - 1) * line_delta(axcalt)
           ELSE
             ds_delstep(dset)  = 0.0
             ds_1step(dset)    = GET_LINE_COORD(linemem(axcalt)%ptr, 1)
             ds_lastep(dset)   =
     .           GET_LINE_COORD(linemem(axcalt)%ptr, line_dim(axcalt))
           ENDIF

* correct all time-dependent variables in this dataset
           DO ivar = 1, maxvars
              IF ( ds_var_setnum(ivar) .EQ. dset ) THEN
                 gg = ds_grid_number(ivar)
                 IF ( grid_line(t_dim,gg) .EQ. taxis .OR.
     .                grid_line(t_dim,gg) .EQ. axcalt ) THEN
                    ds_grid_end(t_dim,ivar) = line_dim(axcalt)
                 ENDIF
                 IF ( grid_line(t_dim,gg) .EQ. taxis ) THEN
* ... was this grid modified in this loop an instant ago?
                    DO idx = 1, ivar-1
                       IF ( ds_var_setnum(idx) .EQ. dset ) THEN
                          IF ( ds_grid_number(idx) .EQ. gg ) GOTO 500
                       ENDIF
                    ENDDO
                    grid_line(t_dim,gg) = axcalt
* ... "use" the given time axis; "un-use" the one from the file
                    CALL TM_USE_LINE(axcalt)
                    CALL TM_DEALLO_DYN_LINE(taxis)
                 ENDIF
              ENDIF
  500         CONTINUE
           ENDDO

        ENDDO

        status = ferr_ok
 5000   RETURN
        END

