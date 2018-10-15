c-ml============ some subroutines from cernlib============================================================

      FUNCTION DFRSIN(X)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      

      DIMENSION A(0:16),B(0:15),C1(0:25),C2(0:28)

      PARAMETER (Z1 = 1, R8 = Z1/8, R32 = Z1/32)

      DATA C0 /1.25331 41373 15500 3D0/

      DATA NA,NB,NC1,NC2 /16,15,25,28/

      DATA A( 0) / 0.76435 13866 41860 002D0/
      DATA A( 1) /-0.43135 54754 76601 793D0/
      DATA A( 2) / 0.43288 19997 97266 531D0/
      DATA A( 3) /-0.26973 31033 83871 110D0/
      DATA A( 4) / 0.08416 04532 08769 354D0/
      DATA A( 5) /-0.01546 52448 44613 820D0/
      DATA A( 6) / 0.00187 85542 34398 220D0/
      DATA A( 7) /-0.00016 26497 76188 875D0/
      DATA A( 8) / 0.00001 05739 76563 833D0/
      DATA A( 9) /-0.00000 05360 93398 892D0/
      DATA A(10) / 0.00000 00218 16584 549D0/
      DATA A(11) /-0.00000 00007 29016 212D0/
      DATA A(12) / 0.00000 00000 20373 325D0/
      DATA A(13) /-0.00000 00000 00483 440D0/
      DATA A(14) / 0.00000 00000 00009 865D0/
      DATA A(15) /-0.00000 00000 00000 175D0/
      DATA A(16) / 0.00000 00000 00000 003D0/

      DATA B( 0) / 0.63041 40431 45705 392D0/
      DATA B( 1) /-0.42344 51140 57053 335D0/
      DATA B( 2) / 0.37617 17264 33436 566D0/
      DATA B( 3) /-0.16249 48915 45095 674D0/
      DATA B( 4) / 0.03822 25577 86330 087D0/
      DATA B( 5) /-0.00564 56347 71321 909D0/
      DATA B( 6) / 0.00057 45495 19768 974D0/
      DATA B( 7) /-0.00004 28707 15321 020D0/
      DATA B( 8) / 0.00000 24512 07499 233D0/
      DATA B( 9) /-0.00000 01109 88418 409D0/
      DATA B(10) / 0.00000 00040 82497 317D0/
      DATA B(11) /-0.00000 00001 24498 302D0/
      DATA B(12) / 0.00000 00000 03200 484D0/
      DATA B(13) /-0.00000 00000 00070 324D0/
      DATA B(14) / 0.00000 00000 00001 336D0/
      DATA B(15) /-0.00000 00000 00000 022D0/

      DATA C1( 0) / 0.99056 04793 73497 549D0/
      DATA C1( 1) /-0.01218 35098 31478 997D0/
      DATA C1( 2) /-0.00248 27428 23113 060D0/
      DATA C1( 3) / 0.00026 60949 52647 247D0/
      DATA C1( 4) /-0.00000 10790 68987 406D0/
      DATA C1( 5) /-0.00000 48836 81753 933D0/
      DATA C1( 6) / 0.00000 09990 55266 368D0/
      DATA C1( 7) /-0.00000 00750 92717 372D0/
      DATA C1( 8) /-0.00000 00190 79487 573D0/
      DATA C1( 9) / 0.00000 00090 90797 293D0/
      DATA C1(10) /-0.00000 00019 66236 033D0/
      DATA C1(11) / 0.00000 00001 64772 911D0/
      DATA C1(12) / 0.00000 00000 63079 714D0/
      DATA C1(13) /-0.00000 00000 36432 219D0/
      DATA C1(14) / 0.00000 00000 10536 930D0/
      DATA C1(15) /-0.00000 00000 01716 438D0/
      DATA C1(16) /-0.00000 00000 00107 124D0/
      DATA C1(17) / 0.00000 00000 00204 099D0/
      DATA C1(18) /-0.00000 00000 00090 064D0/
      DATA C1(19) / 0.00000 00000 00025 506D0/
      DATA C1(20) /-0.00000 00000 00004 036D0/
      DATA C1(21) /-0.00000 00000 00000 570D0/
      DATA C1(22) / 0.00000 00000 00000 762D0/
      DATA C1(23) /-0.00000 00000 00000 363D0/
      DATA C1(24) / 0.00000 00000 00000 118D0/
      DATA C1(25) /-0.00000 00000 00000 025D0/

      DATA C2( 0) / 0.04655 77987 37516 4561D0/
      DATA C2( 1) / 0.04499 21302 01239 4140D0/
      DATA C2( 2) /-0.00175 42871 39651 4532D0/
      DATA C2( 3) /-0.00014 65340 02581 0678D0/
      DATA C2( 4) / 0.00003 91330 40863 0159D0/
      DATA C2( 5) /-0.00000 34932 28659 7731D0/
      DATA C2( 6) /-0.00000 03153 53003 2345D0/
      DATA C2( 7) / 0.00000 01876 58200 8529D0/
      DATA C2( 8) /-0.00000 00377 55280 4930D0/
      DATA C2( 9) / 0.00000 00026 65516 5010D0/
      DATA C2(10) / 0.00000 00010 88144 8122D0/
      DATA C2(11) /-0.00000 00005 35500 7671D0/
      DATA C2(12) / 0.00000 00001 31576 5447D0/
      DATA C2(13) /-0.00000 00000 15286 0881D0/
      DATA C2(14) /-0.00000 00000 03394 7646D0/
      DATA C2(15) / 0.00000 00000 02702 0267D0/
      DATA C2(16) /-0.00000 00000 00946 3142D0/
      DATA C2(17) / 0.00000 00000 00207 1565D0/
      DATA C2(18) /-0.00000 00000 00012 6931D0/
      DATA C2(19) /-0.00000 00000 00013 9756D0/
      DATA C2(20) / 0.00000 00000 00008 5929D0/
      DATA C2(21) /-0.00000 00000 00003 1070D0/
      DATA C2(22) / 0.00000 00000 00000 7515D0/
      DATA C2(23) /-0.00000 00000 00000 0648D0/
      DATA C2(24) /-0.00000 00000 00000 0522D0/
      DATA C2(25) / 0.00000 00000 00000 0386D0/
      DATA C2(26) /-0.00000 00000 00000 0165D0/
      DATA C2(27) / 0.00000 00000 00000 0050D0/
      DATA C2(28) /-0.00000 00000 00000 0009D0/

       V=ABS(X)
       IF(V .LT. 8) THEN
       Y=R8*V
       H=2*Y**2-1
       ALFA=H+H
       B1=0
       B2=0
       DO 4 I = NB,0,-1
       B0=B(I)+ALFA*B1-B2
       B2=B1
    4  B1=B0
       H=SQRT(V)*Y*(B0-B2)
      ELSE
       R=1/V
       H=10*R-1
       ALFA=H+H
       B1=0
       B2=0
       DO 5 I = NC1,0,-1
       B0=C1(I)+ALFA*B1-B2
       B2=B1
    5  B1=B0
       S=B0-H*B2
       B1=0
       B2=0
       DO 6 I = NC2,0,-1
       B0=C2(I)+ALFA*B1-B2
       B2=B1
    6  B1=B0
       H=C0-SQRT(R)*(S*COS(V)+(B0-H*B2)*SIN(V))
      END IF
      IF(X .LT. 0) H=-H
      GO TO 9

      ENTRY DFRCOS(X)

       V=ABS(X)
       IF(V .LT. 8) THEN
       H=R32*V**2-1
       ALFA=H+H
       B1=0
       B2=0
       DO 1 I = NA,0,-1
       B0=A(I)+ALFA*B1-B2
       B2=B1
    1  B1=B0
       H=SQRT(V)*(B0-H*B2)
      ELSE
       R=1/V
       H=10*R-1
       ALFA=H+H
       B1=0
       B2=0
       DO 2 I = NC1,0,-1
       B0=C1(I)+ALFA*B1-B2
       B2=B1
    2  B1=B0
       S=B0-H*B2
       B1=0
       B2=0
       DO 3 I = NC2,0,-1
       B0=C2(I)+ALFA*B1-B2
       B2=B1
    3  B1=B0
       H=C0-SQRT(R)*((B0-H*B2)*COS(V)-S*SIN(V))
      END IF
      IF(X .LT. 0) H=-H
    9 DFRSIN=H
      RETURN
      END

c--#include "gen/pilot.h"
      FUNCTION CGAMMA(Z)
 
      
      COMPLEX*8 CGAMMA
      COMPLEX*8 Z,U,V,F,H,S       
      CHARACTER NAME*(*)
      CHARACTER*80 ERRTXT      
      PARAMETER (NAME = 'CGAMMA')

      DIMENSION C(0:15)

      PARAMETER (Z1 = 1, HF = Z1/2)

c--#if defined(CERNLIB_QF2C)
c--#include "gen/gcmpfun.inc"
c--#endif

      DATA PI /3.14159 26535 89793 24D0/
      DATA C1 /2.50662 82746 31000 50D0/

      DATA C( 0) / 41.62443 69164 39068D0/
      DATA C( 1) /-51.22424 10223 74774D0/
      DATA C( 2) / 11.33875 58134 88977D0/
      DATA C( 3) / -0.74773 26877 72388D0/
      DATA C( 4) /  0.00878 28774 93061D0/
      DATA C( 5) / -0.00000 18990 30264D0/
      DATA C( 6) /  0.00000 00019 46335D0/
      DATA C( 7) / -0.00000 00001 99345D0/
      DATA C( 8) /  0.00000 00000 08433D0/
      DATA C( 9) /  0.00000 00000 01486D0/
      DATA C(10) / -0.00000 00000 00806D0/
      DATA C(11) /  0.00000 00000 00293D0/
      DATA C(12) / -0.00000 00000 00102D0/
      DATA C(13) /  0.00000 00000 00037D0/
      DATA C(14) / -0.00000 00000 00014D0/
      DATA C(15) /  0.00000 00000 00006D0/

c----#if !defined(CERNLIB_QF2C)
c----#include "gen/gcmpfun.inc"
c----#endif

      COMPLEX*8 GREAL,GIMAG,XARG,YARG       

      COMPLEX*8  ZARG,GCONJG,GCMPLX                                                  
           GREAL( ZARG)=REAL( ZARG)                                                  
	   GIMAG( ZARG)=AIMAG(ZARG)                                                  
	   GCONJG(ZARG)=CONJG(ZARG)                                                  
c--           GCMPLX(XARG,YARG)= CMPLX(XARG,YARG)     

      U=Z
      X=U
      IF(GIMAG(U) .EQ. 0 .AND. -ABS(X) .EQ. INT(X)) THEN
       F=0
       H=0
       WRITE(ERRTXT,101) X
c-       CALL MTLPRT(NAME,'C305.1',ERRTXT)
      ELSE
       IF(X .GE. 1) THEN
        F=1
        V=U
       ELSEIF(X .GE. 0) THEN
        F=1/U
        V=1+U
       ELSE
        F=1
        V=1-U
       END IF
       H=1
       S=C(0)
       DO 1 K = 1,15
       H=((V-K)/(V+(K-1)))*H
    1  S=S+C(K)*H
       H=V+(4+HF)
       H=C1*EXP((V-HF)*LOG(H)-H)*S
       IF(X .LT. 0) H=PI/(SIN(PI*U)*H)
      ENDIF

c----#if !defined(CERNLIB_DOUBLE)
      CGAMMA=F*H
c----#endif

      RETURN
  101 FORMAT('ARGUMENT EQUALS NON-POSITIVE INTEGER = ',1P,E15.1)
      END

