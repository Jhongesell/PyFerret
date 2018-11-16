        SUBROUTINE CANCEL_ONE_DATA_SET(dset, status)

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
* close a single dataset
* if this dataset is an aggregation remove its special axes

* programmer - steve hankin
* NOAA/PMEL, Seattle, WA - Tropical Modeling and Analysis Program
* written for VAX computer under VMS operating system
*
* V200:	10/20/89 - extracted from XEQ_CANCEL with enhancements to accept /ALL
*		   and data set names
*	 12/8/89 - if saved_dset is cancelled then indicate that
* Unix/RISC port 8/91 *sh*
*     bug fix: calling var, iset, was used as an internal loop index
* Linux Port 1/97 *kob* - Added a preprocessor include for tmap_dset.parm
*			  because it needed a preprocessor.
* V533: *sh* 6/01 - close set after deletions to avoid "%%" diagnostics
* V685 *acm* 4/13 - If its a member of an aggregation, also close the agg
* V698 *sh* 1/16 - complexity removed to the parent routine, CANCEL_DATA_SET
* V698  3/16 *acm* Removed Philander-Seigel code including diagnostic 
*                  variables, pre-defined axes and grids, cat_calc_var.
*            *sh*  deallo special T and F axes even if user has re-defined them
*      11/18 *kms* special T and F axes no longer created; remove (now unused) status return

* calling arguments
*  dset      -  (input) dataset number to be closed
*  status    - (output) return status

        IMPLICIT NONE
	include 'tmap_dims.parm'
	include 'ferret.parm'
	include 'xprog_state.cmn'
	include 'xdset_info.cmn_text'
	external xdset_info_data

* calling argument declarations:
	INTEGER	dset, status

* since ensembles create so much dataset proliferation, need a diagnostic
* message to see what is being canceled
	IF (mode_diagnostic) CALL SPLIT_LIST(pttmode_ops, ttout_lun,
     .		' canceling dset '//ds_des_name(dset), 0 )

* delete any memory variables from this data set
* ( necessary because their pointers will lead to incorrect data set 
*   information and possibly misleading output )
	CALL ERASE_DSET( dset )

* really close and cancel the data set
        CALL TM_CLOSE_SET( dset, status )

        RETURN

	END
