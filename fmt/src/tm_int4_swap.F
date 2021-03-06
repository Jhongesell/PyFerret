	SUBROUTINE TM_INT4_SWAP(ival4, author_cpu)

*
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
*
* byte-swap a 4-byte quantity to achieve big/little endian conversion

* programmer - steve hankin
* NOAA/PMEL, Seattle, WA - Tropical Modeling and Analysis Program
* written for VAX/VMS, DEC Ultrix, and SUN
*
* revision 0.0 - 2/13/92
* replaced "elif" syntax with
*	else
*	   if
*  for SGI port	 - kob 4/8/92
* V6.74 3/12 *acm* cleanup ifdefs and unnecessary include files

* calling argument declarations:
	INTEGER ival4, author_cpu

* internal variable declarations:
        LOGICAL swap

      include 'cpu.parm'
#  ifdef sun
      BYTE bb(4),bdum
#  else
      INTEGER*1 bb(4),bdum
#  endif

        INTEGER*4 dum
	EQUIVALENCE (bb(1),dum)

***********   CONVERT between machine types   **************
*** use the pre-processor to select the target word type ***
#ifdef sun
* SUN is the platform in use
      swap = author_cpu .EQ. cptype_vax
     .  .OR. author_cpu .EQ. cptype_dec

#else
#	if unix
* DECstation is the platform in use
      		swap = author_cpu .EQ. cptype_sun

#	else
* VAX is the platform in use
      		swap = author_cpu .EQ. cptype_sun
#	endif
#endif

      IF (swap ) THEN
         dum = ival4
         bdum = bb(1)
         bb(1) = bb(4)
         bb(4) = bdum
         bdum = bb(2)
         bb(2) = bb(3)
         bb(3) = bdum
         ival4 = dum
      ENDIF

      RETURN
      END
