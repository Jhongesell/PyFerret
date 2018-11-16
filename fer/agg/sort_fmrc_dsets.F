      SUBROUTINE SORT_FMRC_DSETS(nsets, member_sets, member_filenames,
     .                           torigins, status)
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

* Sort the given datasets using the time origin (T0) of the time axis
* of each dataset.  Returns the sorted datasets as well as the time
* origins.  An error is returned if duplicate time origins are found.

* calling arguments
*   nsets            - number of datasets given to be aggregated (input)
*   member_sets      - array of open dataset numbers (input-output)
*   member_filenames - array of dataset filenames (input-output)
*   torigins         - array of sorted time origins (output)
*   status           - return status (output)

        IMPLICIT NONE
        include 'tmap_dims.parm'
        include 'ferret.parm'
        include 'errmsg.parm'
        include 'xdset_info.cmn_text'
        external xdset_info_data
        include 'xtm_grid.cmn_text'

* calling argument declarations:
        INTEGER nsets, member_sets(nsets), status
        CHARACTER*512 member_filenames(nsets)
        REAL*8  torigins(nsets)

* local variable declarations:
        INTEGER TM_LENSTR1, TM_GET_CALENDAR_ID,
     .          iset, dset, taxis, calid, itmp_memb(nsets)
        CHARACTER*512 itmp_names(nsets)
        CHARACTER calname*32, torigin_str*20
        REAL*8  SECS_FROM_BC, rindices(nsets)
        LOGICAL TM_DFPEQ

* extract the T0 time from each dataset
        DO iset = 1, nsets
           dset  = member_sets(iset)
           taxis = ds_time_axis(dset)
           IF ( taxis .LE. 0 )
     .        CALL ERRMSG(ferr_aggregate_error, status,
     .                    'no time axis in a forecast dataset', *5000)
           calname = line_cal_name(taxis)
           calid = TM_GET_CALENDAR_ID(calname)
           IF ( calid .LE. 0 )
     .        CALL ERRMSG(ferr_aggregate_error, status,
     .                    'unknown calendar "'
     .                    //calname(1:TM_LENSTR1(calname))
     .                    //'" in a forecast dataset', *5000)
           torigin_str = line_t0(taxis)
           torigins(iset) = SECS_FROM_BC(torigin_str, calid, status)
           IF ( status .NE. ferr_ok )
     .        CALL ERRMSG(ferr_aggregate_error, status,
     .                    'invalid time origin "'
     .                    //torigin_str(1:TM_LENSTR1(torigin_str))
     .                    //'" in a forecast dataset', *5000)
        ENDDO

* sort the dataset list according to the t0 times
* ... put indices into REAL array
        DO iset = 1, nsets
           rindices(iset) = iset
        ENDDO
*       (this assumes REAL of HEAP2 == REAL*8 of torigins, rindices)
        CALL HEAP2(torigins, rindices, unspecified_val8, nsets)
* ... reorder the dataset numbers to match the sort
        DO iset = 1, nsets
           itmp_memb(iset) = member_sets(NINT(rindices(iset)))
        ENDDO
* ... copy the sorted list back into the input array
        DO iset = 1, nsets
           member_sets(iset) = itmp_memb(iset)
        ENDDO
* ... reorder the dataset filenams to match the sort
        DO iset = 1, nsets
           itmp_names(iset) = member_filenames(NINT(rindices(iset)))
        ENDDO
* ... copy the sorted list back into the input array
        DO iset = 1, nsets
           member_filenames(iset) = itmp_names(iset)
        ENDDO

* check for unique time origins
        DO iset = 2, nsets
           IF ( TM_DFPEQ(torigins(iset-1), torigins(iset)) ) THEN
              dset  = member_sets(iset)
              taxis = ds_time_axis(dset)
              torigin_str = line_t0(taxis)
              CALL ERRMSG(ferr_aggregate_error, status,
     .                    'duplicate time origin "'
     .                    //torigin_str(1:TM_LENSTR1(torigin_str))
     .                    //'" in a forecast dataset', *5000)
           ENDIF
        ENDDO

 5000   RETURN
        END
