	SUBROUTINE DO_XYAVE_REGRID( cx_lims, src, msrc, dst, mdst,
     .				    work, src_cx, dst_cx, iwork2, itot, 
     .				    aweight )

*
*
*  This software was developed by the Thermal Modeling and Analysis
*  Project(TMAP) of the National Oceanographic and Atmospheric
*  Administration''s (NOAA) Pacific Marine Environmental Lab(PMEL),
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
* call the xy-area regridding routine with appropriate work arrays

* programmer - steve hankin
* NOAA/PMEL, Seattle, WA - Tropical Modeling and Analysis Program
* written for VAX computer under VMS operating system
*
* V200: 12/1/89 - 4D symmetrical version - regrid only 1 axis
* V552:  4/04 *acm* - bug fix.  For @AVE and related transforms @SUM, etc.
*                     send src_lo_ss and src_hi_ss for the range on the src
*                     axis.  
*        2/12 *kms* Add E and F dimensions (use nferdims in tmap_dims.parm)

        IMPLICIT NONE
        include 'tmap_dims.parm'
	include	'ferret.parm'
	include	'xcontext.cmn'

* calling argument declarations:
	INTEGER	cx_lims, axis, msrc, mdst, iwork2, itot, src_cx, dst_cx,
     .		aweight
	REAL	src(*), dst(*), work(itot)

	CALL XYAVE_REGRID( cx_lims,
     .			   cx_lo_ss(cx_lims,x_dim), cx_hi_ss(cx_lims,x_dim),
     .			   cx_lo_ss(cx_lims,y_dim), cx_hi_ss(cx_lims,y_dim),
     .			   src, msrc, dst, mdst, src_cx, dst_cx, 
     .                     work(1), work(iwork2), aweight )

	RETURN
	END
