        SUBROUTINE ASSIGN_F_AGGREGATE_GRID(agg_dset, fvar, agg_fvar,
     .             axcalf, axcalt, findices, memb_sets, nsets, status)

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
* arguments:
*       agg_dset  - aggregation dataset number (input)
*       fvar      - number of the variable to assign to the aggregate grid (input)
*       agg_fvar  - aggregrate variable number (input)
*       axcalf    - time-of-forecast (F) axis to use for the aggregate grid (input)
*       axcalt    - forecasted-times (T) axis to use for the aggregate grid (input)
*       findices  - index of each time-axis-origin on the time-of-forecast (F) axis (input)
*       memb_sets - numbers of the (MC) forecast datasets using axcalt as their T axis (input)
*       nsets     - number of forecast datasets (input)
*       status    - return status (output)

        IMPLICIT NONE

#include "tmap_dset.parm"
        include 'tmap_dims.parm'
        include 'ferret.parm'
        include 'errmsg.parm'
        include 'xdset_info.cmn_text'
        external xdset_info_data
        include 'xdset_info.cd_equiv'
        include 'xtm_grid.cmn_text'
        include 'xvariables.cmn'
        include 'netcdf.inc'

* calling argument declarations:
        INTEGER agg_dset, fvar, agg_fvar, axcalf, axcalt,
     .          findices(nsets), mc_dsets(nsets), memb_sets(nsets),
     .          nsets, status

* local variable declarations:
        INTEGER STR_SAME, TM_LENSTR1,
     .          var_dset, ref_grid, igrid, agg_grid, idim,
     .          nv, ndim, slen, fvtype, ii, iset, numidx, idx,
     .          digrid, dnv
        CHARACTER vname*128

* initialize
        var_dset = ds_var_setnum(fvar)
        vname    = ds_var_code(fvar)
        ref_grid = ds_grid_number(fvar)
        igrid    = unspecified_int4

* Define aggregate grid
        agg_grid = unspecified_int4
        IF ( ref_grid .NE. unspecified_int4 ) THEN
           CALL TM_ALLO_TMP_GRID(agg_grid, status)
           IF ( status .NE. ferr_ok )
     .        CALL ERRMSG(ferr_aggregate_error, status,
     .                    'Too many grids defined', *5000)
           DO idim = 1, nferdims
              grid_line(idim, agg_grid) = grid_line(idim, ref_grid)
           ENDDO
           grid_line(t_dim, agg_grid) = axcalt
           grid_line(f_dim, agg_grid) = axcalf
* Make up a name for aggregate grid - duplicates are resolved in AGG_CLEAN_GRIDS
           grid_name(agg_grid) = 'GFN'
        ENDIF

* Find an empty ds_var slot and fill with aggregate variable in aggregate dataset.
        nv = 1
        DO WHILE ( STR_SAME(ds_var_code(nv), char_init) .NE. 0 )
           nv = nv + 1
           IF ( nv .GT. maxvars )
     .        CALL ERRMSG(ferr_aggregate_error, status,
     .                    'Too many variables defined', *5000)
        ENDDO
        ds_var_setnum(nv)   = agg_dset
        slen = TM_LENSTR1(vname)
        CALL string_array_modify(ds_var_code_head, nv, vname, slen)
        IF ( ds_var_title(nv) .NE. char_init128 ) THEN
           ds_var_title(nv) = ds_var_title(fvar)
        ELSE
           ds_var_title(nv) = ds_var_code(nv)
        ENDIF
        ds_var_titl_mod(nv) = ds_var_titl_mod(fvar)
        ds_var_units(nv)    = ds_var_units(fvar)
        ds_missing_flag(nv) = ds_missing_flag(fvar)
        ds_bad_flag(nv)     = ds_bad_flag(fvar)
        ds_var_type(nv)     = ds_var_type(fvar)
        cd_scale_factor(nv) = 1.
        cd_add_offset(nv)   = 0.
        cd_scaled(nv)       = .FALSE.

        IF ( ref_grid .NE. unspecified_int4 ) THEN
           ndim = 1
           DO idim = 1, nferdims
              ds_grid_start(idim, nv) = ds_grid_start(idim, fvar)
              ds_grid_end(idim, nv)   = ds_grid_end(idim, fvar)
              ds_ordering(idim, nv)   = ds_ordering(idim, fvar)
              IF ( ds_ordering(ndim, fvar) .GT. 0 ) ndim = ndim + 1
           ENDDO

           ds_grid_start(f_dim, nv) = 1
           ds_grid_end(f_dim, nv) = line_dim(axcalf)
           ds_ordering(ndim, nv) = f_dim  ! last dimension of the grid

*          time axis information should already have been reset to axcalt
*          in DEFINE_F_AGG_T_MC_DSETS

           ds_grid_number(nv) = agg_grid
        ENDIF

* Type should be what is in the original file
        IF ( ds_var_type(nv) .EQ. ptype_string ) THEN
           fvtype = NCCHAR
        ELSE
           slen = TM_LENSTR1(vname)
           CALL CD_GET_VAR_TYPE(var_dset,0,vname(:slen),fvtype,status)
        ENDIF

* Hook the aggregated variable into the linked list structures
        CALL CD_ADD_VAR(agg_dset, vname, agg_fvar, ds_var_title(nv),
     .                  ds_var_units(nv), fvtype, .FALSE.,
     .                  ds_bad_flag(nv), status)
        IF ( status .NE. ferr_ok )
     .     CALL ERRMSG(ferr_aggregate_error, status,
     .                 'Unable to add an aggregate variable', *5000)

* Add aggregate-grid info to the lists in aggregate datset and aggregate-grid list.
* Mark user-variables with a grid number of -999.  Note if there is a LET/D variable
* that duplicates the name of a file var it is the file var that gets registered,
* even though the LET/D var has precedence during evaluation.  This allows for the
* possibility that the LET/D var may be cancelled.  But it also allows the possibility
* that the grid of the fvar is unacceptable -- a chance of a crash for a pathological
* dataset.
        DO ii = 1, nsets
           iset = memb_sets(ii)
           IF ( iset .NE. pdset_dummy ) THEN
              DO nv = 1, maxvars
                 IF ( ds_var_setnum(nv) .NE. iset ) CYCLE
                 IF ( STR_SAME(ds_var_code(nv), vname) .NE. 0 ) CYCLE
                 igrid = ds_grid_number(nv)
                 fvtype = cat_file_var
                 CALL CD_ADD_AGG_VAR_INFO(agg_dset, vname, ii, fvtype,
     .                                iset, igrid, axcalf, nv, status)
                 IF ( status .NE. ferr_ok )
     .              CALL ERRMSG(ferr_aggregate_error, status,
     .                          'Unable to add file var info', *5000)
                 GOTO 100
              ENDDO
*             not a regular file variable - search user-defined variables
* ... add this member data set as a LET/D as part of this aggregated variable
* ... optimize search using deleted_list_get_undel
!             DO nv = 1, max_uvar      ! brute-force equivalent loop
              CALL deleted_list_get_undel(uvar_num_items_head,
     .                                    deleted_list_result(1),
     .                                    max_uvar, numidx)
              DO idx = 1, numidx
                 nv = deleted_list_result(idx)
                 IF ( uvar_num_items(nv) .EQ. uvar_deleted ) CYCLE
                 IF ( uvar_dset(nv) .NE. iset) CYCLE
                 IF ( STR_SAME(uvar_name_code(nv), vname) .NE. 0 ) CYCLE
                 CALL GET_SAVED_UVAR_GRID(nv, iset, igrid)
                 fvtype = cat_user_var
                 CALL CD_ADD_AGG_VAR_INFO(agg_dset, vname, ii, fvtype,
     .                                iset, igrid, axcalf, nv, status)
                 IF ( status .NE. ferr_ok )
     .              CALL ERRMSG(ferr_aggregate_error, status,
     .                          'Unable to add user var info', *5000)
* ... bump the use count for this grid to protect it from being deleted if
*     the user changes the LET/D definition -- grid is needed during the
*     gather process (IS_STRIP & IS_AGGREGATE_GATHER).
*  ==> be sure to un-use (deallo) it when the aggregation gets canceled
*  ==> or if an error occurs prior to the complete initialization of this agg
                 CALL TM_USE_DYN_GRID(igrid)
                 GOTO 100
              ENDDO
           ELSE
* ... add the dummy data set as a part of this variable
*     Note that the first member of the agg by common sense cannot be a dummy
              fvtype = cat_file_var
              CALL CD_ADD_AGG_VAR_INFO(agg_dset, vname, ii, fvtype,
     .                       pdset_dummy, digrid, axcalf, dnv, status)
                 IF ( status .NE. ferr_ok )
     .              CALL ERRMSG(ferr_aggregate_error, status,
     .                          'Unable to add user global var info', *5000)
           ENDIF
* ... grab information from the first dataset to reuse with dummy datasets
  100      IF ( ii .EQ. 1 ) THEN
*    ... OK_AS_AGG_VAR guarantees grid is always valid -- even for LET/D vars
              IF (igrid .EQ. unspecified_int4)
     .           STOP 'DEFINE_F_AGGREGATE_GRID grid err'
              dnv = nv
              digrid = igrid
           ENDIF
        ENDDO

 5000   RETURN
        END
