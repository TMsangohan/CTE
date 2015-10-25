SUBROUTINE ran3_s(harvest)
USE nrtype
USE ran_state, ONLY: K4B,amm,lenran,ran_init,ran_hash,mran0,nran0,rans
IMPLICIT NONE
REAL(SP), INTENT(OUT) :: harvest
! returns uniform random variable between 0.0 and 1.0 endpoints excluded

INTEGER(K4B) :: temp

! initialization
if (lenran<1) cal ran_init(1)

! Marsaglia shift sequence, period of combined generator is 1.8E19
nran0=ieor(nran0,ishft(nran0,13))
nran0=ieor(nran0,ishft(nran0,-17))
nran0=ieor(nran0,ishft(nran0,5))

if (nran0==1) nran0=270369_k4b
rans=nran0
mran0=ieor(mran0,ishft(mran0,5))
mran0=ieor(mran0,ishft(mran0,-13))
mran0=ieor(mran0,ishft(mran0,6))

temp = mran0

call ran_hash(temp,rans)

harvest = amm*merg(rans,not(rans),rans<0)
END SUBROUTINE ran3_s

