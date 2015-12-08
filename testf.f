      subroutine TEST(FUNC,BOOLARG,DBLARG,INTARG)

      double precision DBLARG,FUNC
      integer INTARG
      logical BOOLARG
      external FUNC

      BOOLARG = .TRUE.

c remember: argument to function **MUST** be forced to double (8 bytes) here,
c because FORTRAN converts -2.5 to real (4 bytes)
      DBLARG = FUNC(-2.5D0)

      INTARG = 1000

      end

