C************************************************************************
C                              AMBER                                   **
C                                                                      **
C                  Copyright (c) 1986, 1991, 1995                      **
C             Regents of the University of California                  **
C                       All Rights Reserved.                           ** 
C                                                                      **
C  This software provided pursuant to a license agreement containing   **
C  restrictions on its disclosure, duplication, and use. This software **
C  contains confidential and proprietary information, and may not be   **
C  extracted or distributed, in whole or in part, for any purpose      **
C  whatsoever, without the express written permission of the authors.  **
C  This notice, and the associated author list, must be attached to    **
C  all copies, or extracts, of this software. Any additional           **
C  restrictions set forth in the license agreement also apply to this  **
C  software.                                                           **
C************************************************************************
C
      subroutine ngfil
c     Author: George Seibel
c     gets unix command line input for program NUCGEN.
*     implicit none
c
c     OUTPUT: (to common)
c
      common /files/ ngin, ngout, ngdat, pdbout, owrite
      character*80   ngin, ngout, ngdat, pdbout
      character owrite
c
c        NUCGEN files as described in amber v 3.0 Rev A documentation
c
c     INTERNAL:
c
      character*80 arg
c        ... temp for each of the whitespace delimited command line words
      integer iarg, narg
c        ... arg pointer, final number of arguments
c
c     --- initialize file names ---
c
      ngin   = 'ngin'
      ngout  = 'ngout'
      ngdat  = 'ngdat'
      pdbout = 'pdbout'
c
c     --- default output file status: 'N'ew
c
      owrite = 'N'
#ifndef NOGETARG
c
c     --- get com line arguments ---
c
# ifdef HITACHI_GETARG
      iarg = 1
# else
      iarg = 0
# endif
      indx = iargc()
   10 continue
           iarg = iarg + 1
           call getarg(iarg,arg)
           if (arg .eq. '-O') then
                owrite = 'U'
           elseif (arg .eq. '-i') then
                iarg = iarg + 1
                call getarg(iarg,ngin)
           elseif (arg .eq. '-o') then
                iarg = iarg + 1
                call getarg(iarg,ngout)
           elseif (arg .eq. '-d') then
                iarg = iarg + 1
                call getarg(iarg,ngdat)
           elseif (arg .eq. '-p') then
                iarg = iarg + 1
                call getarg(iarg,pdbout)
           else
                if (arg .eq. ' ') go to 20
                write(6,'(/,5x,a,a)') 'unknown flag: ',arg
                write(6,9000)
                call mexit(6, 1)
           endif
      if (iarg .lt. indx) go to 10
c
   20 continue
      narg = iarg - 1
      if (narg .lt. 2) then
           write(6,9000)
           call mexit(6, 1)
      endif
c
#endif
      return
 9000 format(/5x,
     .   'usage: nucgen [-O] -i ngin -o ngout -d ngdat -p pdbout')
      end
