$  assign [stewart.TESTDATA] testd
$  assign [stewart.TESTDATA.test1] test1
$!
$!                   Validation tests for MOPAC 6.00
$!                   -------------------------------
$!
$!   TEST A B means "run a MOPAC job with the data A set found in directory B"
$!
$!   First go into the directory which will contain the test-results.
$!
$SD TEST1
$ TEST:== TESTD:@TEST
$DEL/NOCON TEST1:*.*;*
$!
$!  Run these tests in the order given here.  
$!
$!  Test of parameters and single SCF.  If this test fails completely,
$!  MOPAC is corrupt.  If any one test fails, the relevant parameters
$!  are faulty.
$!
$!  To check the results fast, search the output file for "HEAT", this 
$!  should occur three times per element.  The results should agree 
$!  to within one unit in the 3rd decimal place.
$!
$@TESTD:TEST ELEMENTS TESTD:
$!
$!  Test of electronic section.  If this test fails, there is a fault
$!  in COMPFG or in a subroutine called by COMPFG.
$!  
$!  To test 1SCF, GEOMETRY, FORCE, and OLDGEO carry out the check given 
$!  in the comments lines of each data-set.  If no comment, then 
$!  "eyeball" the results.  If nothing looks outrageous, the run was 
$!  successfull.
$!
$@TESTD:TEST 1SCF TESTD:
$!
$!  Test of geometric options. If any test fails, the fault lies in
$!  the relevant geometric subroutine.
$!
$@TESTD:TEST GEOMETRY TESTD:
$!
$!  Test of FORCE calculation.  If this test fails, the fault lies in
$!  subroutine FORCE, or a subroutine called by FORCE.
$!
$@TESTD:TEST FORCE    TESTD:
$!
$!  Test of keywords.  Read the whole output.  Some sections which are
$!  highly repetative, can be examined very superficially.
$!
$!  If any test fails, check the relevant keywords.
$!
$@TESTD:TEST KEYS     TESTD:
$!
$!   Test of use of OLDGEO keyword
$!
$@TESTD:TEST OLDGEO   TESTD:
$!
$!   The following tests are the only tests normal users of MOPAC are
$!   asked to do.  The first test is given in the Manual, the second
$!   runs several very different jobs and essentionally is a 
$!   demonstration of MOPAC
$@TESTD:TEST MNRSD1    TESTD:
$@TESTD:TEST TESTDATA  TESTD:
$!
$!
$!  The tests have been made as non-specific as possible so that MOPAC
$!  should pass on any computer.  On the other hand, the main paths through
$!  MOPAC are examined.  The precision specified is good enough for
$!  all normal research work.
$!
$!**********************************************************************
