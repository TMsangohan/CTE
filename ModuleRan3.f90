! ---------------------------------------------------------------------------------------------------------------
! AUTHORS : MIKE BLASKIEWISC, RODERIK BRUCE, MICHAELA SCHAUMANN, TOM MERTENS
! ---------------------------------------------------------------------------------------------------------------
! VERSION 1.0 : 
!   AUTHOR    : TOM MERTENS  
!   DATE      : 11/06/2015
!   COPYRIGHT : CERN
!   REFERENCE : FUNCTION TAKEN FROM NUMERICAL RECIPES IN FORTRAN
!
!   DESCRIPTION : 
!       RANDOM GENERATOR OF REAL NUMBER IN THE INTERVAL -1,1
! ---------------------------------------------------------------------------------------------------------------

MODULE ran
	USE commondata

INTERFACE RAN3
    MODULE PROCEDURE RAN3
END INTERFACE

CONTAINS
FUNCTION RAN3(IDUM)
    USE COMMONDATA
    IMPLICIT NONE

!----------------------------------------------------------------------------------------------------------------
! IDUM : SEED VARIABLE USED AS INPUT TO GENERATE RANDOM NUMBER 
!----------------------------------------------------------------------------------------------------------------
	INTEGER(I4B), INTENT(INOUT)     :: IDUM

	INTEGER(I4B),   PARAMETER       :: MBIG=1000000000,MSEED=161803398,MZ=0
	INTEGER(I4B)                    :: I,II,MJ,MK,K,inext,inextp,IFF
    INTEGER(I4B),   DIMENSION(55)   :: MA

    REAL(SP),       PARAMETER       :: FAC=1.E-9
    REAL(SP)                        :: RAN3
	
	

	!common/inputran3/MA,inext,inextp

    DATA IFF/0/

    IF (IDUM<=0 .OR. IFF==0) THEN
     	IFF=1
     	MJ=MSEED-IABS(IDUM)
     	MJ=MOD(MJ,MBIG)
     	MA(55)=MJ
     	MK=1
     	DO I=1,54
     	  	II=MOD(21*I,55)
     	  	MA(II)=MK
     	  	MK=MJ-MK
     	  	IF(MK<=MZ)MK=MK+MBIG
     	  	MJ=MA(II)
    	ENDDO
     	DO K=1,4
     	  	DO  I=1,55
     	  	 	MA(I)=MA(I)-MA(1+MOD(I+30,55))
     	  	 	IF(MA(I)<MZ)MA(I)=MA(I)+MBIG
     	 	ENDDO
    	ENDDO
     	INEXT=0
     	INEXTP=31
     	IDUM=1
    ENDIF
  
    INEXT=INEXT+1
  
    IF(INEXT==56) INEXT=1
    
    INEXTP=INEXTP+1
    
    IF(INEXTP==56)INEXTP=1
      
    MJ=MA(INEXT)-MA(INEXTP)
    
    IF(MJ<MZ)MJ=MJ+MBIG

    MA(INEXT)=MJ
    RAN3=REAL(MJ*FAC,SP)
   
END FUNCTION RAN3

END MODULE ran