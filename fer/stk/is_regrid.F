	SUBROUTINE IS_REGRID(  *, status )

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
*
* regrid thae data from the source grid to the destination grid

* programmer - steve hankin
* NOAA/PMEL, Seattle, WA - Tropical Modeling and Analysis Program
* written for VAX computer under VMS operating system
*
* revision 0.0 - 7/25/88
* revision 0.1 -12/30/88 - eliminated P1 and Q1 from regrid call
*			 - renamed DO_REGRID --> DO_LIN_REGRID
* revision 0.2 -  1/9/89 - added xy average regridding
* V200:  8/2/89 - 4D symmetrical version - usually 1 axis per regrid
*	10/16/89 - set up adjustable array limits with NON_ARRAY_SUBSC
*	 11/3/89 - if destination limits are unknown use full axis length
*	 12/1/89 - allow XY average regridding - code split to RCVD_REGRID
* V230:  8/22/92 - added @ASN and @AVE regridding
* V300:  3/16/93 - allow regridding where the source range falls entirely
*                  between points on the destination axis (all bad_data result)
* V312: 5/94 - array "memory" as a calling argument
* V420: 9/95	 - Support for implicit grids: deallocate after it is used
*      11/95	 - added modulo (climatological) regridding
*      12/95	 - added modngd, modmin, etc. (climatological) regridding
*	1/96	 - bug check against attempt to regrid a compressed axis
*	1/29/96	 - regridding @SUM and @NGD are not weighted
*		   (done in ave_regrid_1_axis instead of window_regrid_1_axis)
* V491: 10/97 *sh* - add G=@XACT
* V500:  5/99 *sh* - check if implicit regridding was already done by netCDF
*			read using strides (check cx_delta of src cx)
*		     Allow for the possibility of multiple axes
* V530:  3/01 *sh* - added g=ax@nrst
* V530:	 3/01 *acm* return a status flag from rcvd_regrid, tests for
*                   valid calendar axes.
* V533:	 6/01 *sh* - some regridding ops not allowed for string arrays
*		   - mr_delta changes the handling of stride regridding
* v541:  2/02 *kob*- on hankin request add an initialization line for
*                    cx_delta in loop 500
* V541:  4/02 *sh* - bug fix: delta not removed from result despite 
*                  the fact that result is on a new axis where the delta 
*                  is implicit rather than explicit in the axis definition 
*                  see err540_strides_rev.jnl 
* V552:  4/04 *acm* - bug fix.  For @AVE and related transforms @SUM, etc.
*                     send src_lo_ss and src_hi_ss for the range on the src
*                     axis.  Previous use of mr_lo_ss and mr_hi_ss can result
*                     in incorrect results when a larger then needed chunk
*                     of source data is already in4/15/03 memory
* V554:  3/04 *acm* Fix the check for done_by_strides.
* V603 5/07  *acm* - added MODNBD modulo regridding statistic
* V671:  6/11 *acm*  New 2D linear xy regridding, fixing ticket 1862
*                    The XY linear regridding matches past Ferret behavior, 
*                    not using bilinear regridding for now
*       *acm* 3/12 Add E and F dimensions (use nferdims in tmap_dims.parm)
* V68+ 4/12  *acm*   Fix ticket 1931: the work space for XY linear regridding
*                    was computed incorrectly.
* V696 *acm* 11/15 Added @BIN and @NBIN regridding transforms
* V698 *acm*  3/16 Ticket 2278: add an attribute with time range when
*                  saving a modulo-regridded variable to netCDF
* V702 3/17 *sh* for trac enhancement #2369 -- dynamic memory management 
* V720: 6/17 *acm* Ticket 2346: new @IIN regridding transform
* V741  5/18 *acm* issue 1868: gridding transforms that will grid to same axis

        IMPLICIT NONE
	include 'tmap_dims.parm'
	include	'ferret.parm'
	include	'errmsg.parm'
	include	'interp_stack.parm'
	include	'xprog_state.cmn'
	include	'xcontext.cmn'
	include	'xvariables.cmn'
	include	'xalgebra.cmn'   ! with alg_regrid_code

* calling argument declarations:
	INTEGER	status
* normal exit is for error conditions
* exit #1 is to signal completion of the activity at this stack level

* internal variable declarations:
	LOGICAL	TM_ITS_CDF, 
     .		xy_regrid, valid, valid_y, has_impl_save, done_by_strides,
     .		itsa_string
	INTEGER hi_isp, lo_isp, axis, cx_valid, idim,
     .		rgrd_trans, mrlist(5), wk_size, wk_size2,
     .		dst_cx, src_cx, dst_mr, src_mr, ngd_mr, ave_mr,
     .		axis_list(nferdims), naxis, i, nbd_mr, nxsrc
	INTEGER*8 i8_size
	CHARACTER dtype*4

*        ... interpretation stack setup at the time this is called ...
*	stack level		stack contents
*	isp+2	<--	level of next activity
*	hi_isp	<--	level at which component was obtained
*	lo_isp	<--	level at which regrid was requested (current isp)
* --- end of introductory code ---

* initialize
	lo_isp   = isp                     
	hi_isp   = lo_isp + 1
	src_mr   = is_mr( hi_isp )
	dst_cx   = is_cx( lo_isp )
	src_cx   = is_cx( hi_isp )
	axis	 = is_axis( lo_isp )
	naxis	 = 1
	axis_list(1) = axis
	itsa_string = cx_type(src_cx) .EQ. ptype_string

* what is the regridding transformation for this axis ?
	rgrd_trans = cx_regrid_trans( axis, dst_cx )
	xy_regrid = ( (rgrd_trans .EQ. prgrd_xy_ave) .OR. 
     .                (rgrd_trans .EQ. prgrd_xy_lin) )
	IF ( xy_regrid ) THEN
	   IF ( mode_diagnostic ) WRITE (ttout_lun,*) 'regrid XY'
	   IF (axis.NE.x_dim) CALL ERRMSG
     .		( ferr_internal, status, 'rgrd_xy_regrid', *5000 )	! diagnostic
	   naxis = 2
	   axis_list(2) = y_dim
	ENDIF

* bug check:  has someone attempted to regrid a compressed axis?
* (possible bug in check ... could we have unspec limits for another reason?)
	IF ( cx_category(src_cx) .EQ. cat_user_var
     . .AND. cx_lo_ss(src_cx,axis) .EQ. unspecified_int4 ) THEN
	   IF (  uvar_given(axis,cx_variable(src_cx))
     .	    .NE. uvlim_irrelevant )			 CALL ERRMSG
     .		( ferr_invalid_command, status,
     .		 'attempt to regrid on a compressed axis', *5000 )
	ENDIF

* (5/99) for netCDF variables strides may have been used in READ_TM and
* the regridding already performed there. The src cx then contains the
* limits appropriate to netCDF ... not to Ferret. So, make the source cx
* into a clone of the final result.

* (3/04) do not check for has_impl_grid: when we have done a DEFINE GRID/like=var
* for later comparisons, this grid was stored as a static grid, but we have still 
* read the data using strides, and therefore do not need to do a regridding operation.

	IF ( cx_category(dst_cx) .EQ. cat_file_var ) THEN
	   done_by_strides = 
     .		TM_ITS_CDF(cx_data_set(dst_cx), dtype )
c     .    .AND. cx_has_impl_grid(dst_cx)
     .	  .AND. cx_delta(axis,src_cx) .NE. unspecified_val8
	ELSE
	   done_by_strides = .FALSE.
	ENDIF
	IF (done_by_strides) THEN
	   DO 400 idim = axis+1, nferdims
	      IF (cx_delta(idim,src_cx) .NE. unspecified_val8) THEN
	         naxis = naxis + 1
	         axis_list(naxis) = idim
	      ENDIF
 400	   CONTINUE
	ENDIF

* Make the destination context like the source data context on all except
* the regrid axis.  That one keeps the user-requested size.
* (Other axes may have been modified by information embedded within a
* user-variable definition.  An alternative (untested 8/89) approach is to
* call PASS_AMBIG_LI** for uvars)
	CALL TRANSFER_CONTEXT( src_cx, cx_buff )
	DO 500 i = 1, naxis 
	   idim = axis_list(i)
	   CALL TRANSFER_AXIS( idim, dst_cx, cx_buff )
	   cx_regrid_trans( idim, cx_buff ) = cx_regrid_trans( idim, dst_cx )
	   cx_delta(idim,cx_buff) = unspecified_val8
 500	CONTINUE
	cx_grid( cx_buff ) = cx_grid( dst_cx )
	cx_unstand_grid( cx_buff ) = .TRUE.
	has_impl_save = cx_has_impl_grid(dst_cx)
	CALL TRANSFER_CONTEXT( cx_buff, dst_cx )
	cx_has_impl_grid(dst_cx) = has_impl_save	! preserve this flag

* create a context with the limits that can actually be regridded
* ... if limits were unspecified in request allow full axis size
	IF (done_by_strides) THEN
	   cx_valid = dst_cx
	   valid = .TRUE.
	ELSE
	   cx_valid = cx_buff
	   CALL RCVD_REGRID( src_cx, dst_cx, axis, wk_size,
     .		 cx_lo_ss(cx_valid,axis), cx_hi_ss(cx_valid,axis), status )
	   cx_by_ss( axis, cx_valid ) = .TRUE.
           valid = cx_lo_ss(cx_valid,axis) .LE. cx_hi_ss(cx_valid,axis)  ! 3/93
           IF ( valid ) THEN   ! 3/93
	      CALL CONFINE_AXIS( axis, cx_valid, status )
	      IF ( status .NE. ferr_ok ) GOTO 5000
           ENDIF
	   IF ( xy_regrid ) THEN
	      CALL RCVD_REGRID( src_cx, dst_cx, y_dim, wk_size2,
     .		cx_lo_ss(cx_valid,y_dim),cx_hi_ss(cx_valid,y_dim), status)
	      cx_by_ss( y_dim, cx_valid ) = .TRUE.
              valid_y = cx_lo_ss(cx_valid,y_dim)
     .             .LE. cx_hi_ss(cx_valid,y_dim)
              IF ( valid_y ) THEN
	         CALL CONFINE_AXIS( y_dim, cx_valid, status )
	         IF ( status .NE. ferr_ok ) GOTO 5000
              ENDIF
              valid = valid .AND. valid_y
	   ENDIF
	ENDIF

* create a space for the new result possibly with the edges flagged as bad
	CALL NEW_MR_W_EDGES(  cx_valid, dst_cx, dst_mr, status )
	IF ( status .NE. ferr_ok ) GOTO 5000
	is_mr( lo_isp ) = dst_mr

* cannot actually regrid if the ranges aren't valid -> all bad result
        IF ( .NOT.valid ) GOTO 1000      ! 3/93

* diagnostic message ?
	IF (mode_diagnostic) THEN
	   IF (done_by_strides) THEN
	      CALL DIAGNOSTIC_OUT('nulrgd ',dst_mr,  point_to_mr)
	   ELSE
	      CALL DIAGNOSTIC_OUT('regrid ',dst_mr,  point_to_mr)
	   ENDIF
	ENDIF

* set up 4D array limits as unsubscripted variables in COMMON
	mrlist(1) = src_mr
	mrlist(2) = dst_mr
	CALL NON_ARRAY_SUBSC( mrlist, 2 )

* reserve working array space
	IF ( rgrd_trans .EQ. prgrd_xy_ave) THEN
	   i8_size = wk_size+wk_size2
	ELSEIF (rgrd_trans .EQ. prgrd_xy_lin) THEN
	   nxsrc = cx_hi_ss(src_cx,x_dim)-cx_lo_ss(src_cx,x_dim) + 1
	   i8_size = nxsrc*wk_size2+2*(wk_size+wk_size2)
	ELSE
	   i8_size = 2*wk_size
	ENDIF
	CALL GET_DYN_WORK_SPACE( i8_size, status )
	IF ( status .NE. ferr_ok ) GOTO 5000

* do the regridding
	IF ( done_by_strides ) THEN
	   CALL COPY_GRID( memry(src_mr)%ptr,  src_mr,
     .			   memry(dst_mr)%ptr,  dst_mr )
	ELSEIF ( rgrd_trans .EQ. prgrd_linear_interp ) THEN
	   IF ( itsa_string ) GOTO 5010
	   CALL DO_LIN_REGRID( cx_valid, axis,
     .			memry(src_mr)%ptr, src_mr,
     .			memry(dst_mr)%ptr, dst_mr,
     .			workmem(plegacy_work_buffer)%ptr,
     .			wk_size+1, 2*wk_size )
	ELSEIF ( rgrd_trans .EQ. prgrd_associate ) THEN
	   CALL DO_ASN_REGRID( cx_valid, axis,
     .			memry(src_mr)%ptr, src_mr,
     .			memry(dst_mr)%ptr, dst_mr )
	ELSEIF ( rgrd_trans .EQ. prgrd_average
     .	 .OR.    rgrd_trans .EQ. prgrd_variance
     .	 .OR.    rgrd_trans .EQ. prgrd_sum
     .	 .OR.    rgrd_trans .EQ. prgrd_good_pt ) THEN
	   IF ( itsa_string ) GOTO 5010
	   CALL DO_AVE_REGRID( cx_valid, axis, rgrd_trans,
     .			memry(src_mr)%ptr, src_mr,
     .			memry(dst_mr)%ptr, dst_mr,
     .			workmem(plegacy_work_buffer)%ptr, src_cx, dst_cx,
     .                  wk_size+1 )
	ELSEIF ( rgrd_trans .EQ. prgrd_xy_ave ) THEN
	   IF ( itsa_string ) GOTO 5010
	   CALL DO_XYAVE_REGRID( cx_valid,
     .			memry(src_mr)%ptr, src_mr,
     .			memry(dst_mr)%ptr, dst_mr,
     .			workmem(plegacy_work_buffer)%ptr, src_cx, dst_cx,
     .                  wk_size+1, wk_size+wk_size2 )
	ELSEIF ( rgrd_trans .EQ. prgrd_xy_lin ) THEN
	   IF ( itsa_string ) GOTO 5010
	   CALL DO_XYLIN_REGRID( cx_valid,
     .			memry(src_mr)%ptr, src_mr,
     .			memry(dst_mr)%ptr, dst_mr,
     .			workmem(plegacy_work_buffer)%ptr, src_cx, dst_cx,
     .                  wk_size+1, 2*wk_size+1, 
     .                  2*wk_size+wk_size2+1, 2*(wk_size+wk_size2)+1, 
     .                  nxsrc*wk_size2+2*(wk_size+wk_size2) )
	ELSEIF ( rgrd_trans.EQ.prgrd_min
     .	   .OR. rgrd_trans.EQ.prgrd_max ) THEN
	   IF ( itsa_string ) GOTO 5010
	   CALL DO_WINDOW_REGRID( cx_valid, axis, rgrd_trans,
     .			memry(src_mr)%ptr, src_mr,
     .			memry(dst_mr)%ptr, dst_mr,
     .			workmem(plegacy_work_buffer)%ptr, wk_size+1)
	ELSEIF ( rgrd_trans .EQ. prgrd_xact_pts ) THEN
	   CALL DO_XACT_REGRID( cx_valid, axis,
     .			memry(src_mr)%ptr, src_mr,
     .			memry(dst_mr)%ptr, dst_mr,
     .			workmem(plegacy_work_buffer)%ptr, wk_size+1,
     .			2*wk_size )
	ELSEIF ( rgrd_trans .EQ. prgrd_nearest ) THEN
	   CALL DO_NEAREST_REGRID( cx_valid, axis,
     .			memry(src_mr)%ptr, src_mr,
     .			memry(dst_mr)%ptr, dst_mr,
     .			workmem(plegacy_work_buffer)%ptr, wk_size+1,
     .			2*wk_size )
	ELSEIF ( rgrd_trans .EQ. prgrd_bin
     .	 .OR.    rgrd_trans .EQ. prgrd_nbin ) THEN
	   CALL DO_BIN_REGRID( cx_valid, axis, rgrd_trans,
     .			memry(src_mr)%ptr, src_mr,
     .			memry(dst_mr)%ptr, dst_mr,
     .			workmem(plegacy_work_buffer)%ptr, src_cx, dst_cx,
     .                  wk_size+1 )

	ELSEIF ( rgrd_trans .EQ. prgrd_iin ) THEN
	   CALL DO_IIN_REGRID( cx_valid, axis,
     .			memry(src_mr)%ptr, src_mr,
     .			memry(dst_mr)%ptr, dst_mr,
     .			workmem(plegacy_work_buffer)%ptr, src_cx, dst_cx,
     .                  wk_size+1, status )
     
	   IF ( status .NE. ferr_ok ) GOTO 5020

* modulo regridding
	ELSE
	   IF ( itsa_string ) GOTO 5010
* ... for all modulo regridding cases allocate a temp variable to contain
* ... the modngd (number of valid points) result -- this is used in all
* ... of the modulo calculations.
* ... For modvar (variance) also allocate a temp space for the average
* ... which is needed in the process of computing the variance
* ... For modnbd (number of invalid points) also allocate a temp space for 
* ... the count of number of bad values 
	   ave_mr = dummy_mr	! component 3 default

	   CALL CREATE_TEMP_MEM_VAR(cx_valid,ngd_mr,status)
	   IF ( status .NE. ferr_ok ) GOTO 5100
	   mr_regrid_trans ( axis, ngd_mr ) =  prgrd_modngd

	   CALL CREATE_TEMP_MEM_VAR(cx_valid,nbd_mr,status)
	   IF ( status .NE. ferr_ok ) GOTO 5100
	   mr_regrid_trans ( axis, nbd_mr ) =  prgrd_modnbd

	   IF ( rgrd_trans .EQ. prgrd_modvar ) THEN
	      CALL NEW_MR_W_EDGES(cx_valid,dst_cx,ave_mr,status)
	      IF ( status .NE. ferr_ok ) GOTO 5100
!	      mr_regrid_trans ( axis, ave_mr ) = prgrd_modulo	! bookkeeping
	   ENDIF

* set up 4D array limits and do the regridding
	   mrlist(3) = ngd_mr
	   mrlist(4) = nbd_mr
	   mrlist(5) = ave_mr
	   CALL NON_ARRAY_SUBSC( mrlist, 5 )
	   CALL MODULO_REGRID( cx_valid, axis, rgrd_trans,
     .			memry(src_mr)%ptr, src_mr,
     .			memry(dst_mr)%ptr, dst_mr,
     .			memry(ngd_mr)%ptr, ngd_mr,
     .			memry(nbd_mr)%ptr, nbd_mr,
     .			memry(ave_mr)%ptr, ave_mr )

* ... For modulo regrid in time, add an attribute listing the source dates

	   IF (axis .EQ. t_dim) CALL DEFINE_CLIMAT_ATT (src_cx)

* ... delete temporary variables
	   CALL DELETE_VARIABLE( ngd_mr )
	   CALL DELETE_VARIABLE( nbd_mr )
	   IF ( rgrd_trans .EQ.  prgrd_modvar ) THEN
	      CALL DELETE_VARIABLE( ave_mr )
	   ENDIF

	ENDIF
        CALL  RELEASE_DYN_WORK_SPACE

* indicate we are done with the component
 1000   CALL MR_NOT_IN_USE( src_mr )
	IF (done_by_strides) CALL DELETE_VARIABLE( src_mr ) ! invalid lims in mr

* intermediate grid ? deallocate grid space and delete intermediate data
	IF ( is_phase( lo_isp ) .EQ. prgrd_from_ntrmed ) THEN
	   CALL DEALLO_GRID( status )
	   IF ( status .NE. ferr_ok ) GOTO 5000
	   IF ( mr_protected(src_mr) .EQ. mr_not_protected )
     .			CALL DELETE_VARIABLE ( src_mr )
	ENDIF

* pop off context used by source grid
        cx_stack_ptr = cx_stack_ptr - 1
	RETURN 1

* error exits
 5000	RETURN
 5010	CALL ERRMSG( ferr_invalid_command, status,
     .		'cant use this regrid transform on strings: @'
     .		//alg_regrid_code(rgrd_trans), *5000)
 5020	CALL ERRMSG( ferr_invalid_command, status,
     .   'Units of source and destination axis must match: @'
     .   //alg_regrid_code(rgrd_trans), *5000)
 5100	CALL DELETE_MRS_IN_PROGRESS
	GOTO 5000
	END	





