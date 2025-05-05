      MODULE test_old
      CONTAINS
      subroutine get_common(OUTT)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /GEOPACK1/ AA(10),SPS,CPS,BB(22)
      COMMON /GEOPACK2/ G(105),H(105),REC(105)
      REAL*8 OUTT(105)

      OUTT = G
      end subroutine get_common
      end module test_old
        
