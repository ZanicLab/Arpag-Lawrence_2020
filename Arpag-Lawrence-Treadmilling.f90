PROGRAM treadmilling

	USE ArpagLawrenceTreadmillingGlobals 


	IMPLICIT NONE


    !!  Simulation variables
        INTEGER ::iargc,numarg, i, j, k, seedin, seedout
        CHARACTER(80) ::charinput, buf, command
        REAL(8) ::prob, probtheta, probcenterx, probcentery, dt=1.0d-5, runtime=0.0d0
        INTEGER(KIND=8) ::maxsteps, step=1
        REAL(8), PARAMETER ::PI=3.141592653589793238462643d0
        REAL(8), DIMENSION(:,:), ALLOCATABLE ::lengthprob

    !!  Initialization variables
        INTEGER(KIND=8) ::noofMTs
        REAL(8) ::dimersize=8.0d-9/13.0d0
        REAL(8) ::L0=0.0d0, L0std=0.0d0, minusinitial=0.0d0, plusinitial=0.0d0

    !!  MT Dynamics variables
        REAL(8) ::vspdummy, vsmdummy, fcatpdummy, fcatmdummy, frespdummy, fresmdummy
        REAL(8) ::growthspeedp=0.0d0, shrinkagespeedp=0.0d0, fcatp=0.0d0, fresp=0.0d0, frenucp=0.0d0
        REAL(8) ::growthspeedm=0.0d0, shrinkagespeedm=0.0d0, fcatm=0.0d0, fresm=0.0d0, frenucm=0.0d0
        REAL(8) ::growthspeedstdevp=0.0d0, shrinkagespeedstdevp=0.0d0, fcatstdevp=0.0d0, fresstdevp=0.0d0, frenucstdevp=0.0d0
        REAL(8) ::growthspeedstdevm=0.0d0, shrinkagespeedstdevm=0.0d0, fcatstdevm=0.0d0, fresstdevm=0.0d0, frenucstdevm=0.0d0
        REAL(8) ::totallength, avepluspos, aveminuspos
        INTEGER(8) ::existingMTcount
        LOGICAL ::mtseed=.FALSE., shrink=.TRUE.
        LOGICAL ::terminateatloss=.FALSE.

    !!  Structures
        TYPE (microtubule), DIMENSION(:), ALLOCATABLE::mt
        !TYPE (motor), DIMENSION(:), ALLOCATABLE::motors

!  Read variables from command line
write(*,*) "read comman line parameters"
	numarg=iargc()
DO i=1,numarg
    CALL getarg(i,charinput)
    IF (charinput(1:6).EQ.'-seed=') THEN
        READ(charinput(7:LEN_TRIM(charinput)),*) seedin
    ELSEIF (charinput(1:9).EQ.'-runtime=') THEN
        READ(charinput(10:LEN_TRIM(charinput)),*) runtime
        maxsteps=INT(60.0d0*runtime/dt)
    ELSEIF (charinput(1:9).EQ.'-noofMTs=') THEN
        READ(charinput(10:LEN_TRIM(charinput)),*) noofMTs
    ELSEIF (charinput(1:4).EQ.'-L0=') THEN
        READ(charinput(5:LEN_TRIM(charinput)),*) L0
        L0=1.0d-6*L0
    ELSEIF (charinput(1:7).EQ.'-L0std=') THEN
        READ(charinput(8:LEN_TRIM(charinput)),*) L0std
        L0std=1.0d-6*L0std
    ELSEIF (charinput(1:7).EQ.'-mtseed') THEN
        mtseed=.TRUE.
    ELSEIF (charinput(1:14).EQ.'-growthspeedp=') THEN
        READ(charinput(15:LEN_TRIM(charinput)),*) growthspeedp
        growthspeedp=1.0d-9*growthspeedp
    ELSEIF (charinput(1:19).EQ.'-growthspeedstdevp=') THEN
        READ(charinput(20:LEN_TRIM(charinput)),*) growthspeedstdevp
        growthspeedstdevp=1.0d-9*growthspeedstdevp
    ELSEIF (charinput(1:17).EQ.'-shrinkagespeedp=') THEN
        READ(charinput(18:LEN_TRIM(charinput)),*) shrinkagespeedp
        shrinkagespeedp=1.0d-9*shrinkagespeedp
    ELSEIF (charinput(1:22).EQ.'-shrinkagespeedstdevp=') THEN
        READ(charinput(23:LEN_TRIM(charinput)),*) shrinkagespeedstdevp
        shrinkagespeedstdevp=1.0d-9*shrinkagespeedstdevp
    ELSEIF (charinput(1:7).EQ.'-fcatp=') THEN
        READ(charinput(8:LEN_TRIM(charinput)),*) fcatp
    ELSEIF (charinput(1:12).EQ.'-fcatstdevp=') THEN
        READ(charinput(13:LEN_TRIM(charinput)),*) fcatstdevp
    ELSEIF (charinput(1:7).EQ.'-fresp=') THEN
        READ(charinput(8:LEN_TRIM(charinput)),*) fresp
    ELSEIF (charinput(1:12).EQ.'-fresstdevp=') THEN
        READ(charinput(13:LEN_TRIM(charinput)),*) fresstdevp
    ELSEIF (charinput(1:9).EQ.'-frenucp=') THEN
        READ(charinput(10:LEN_TRIM(charinput)),*) frenucp
    ELSEIF (charinput(1:14).EQ.'-frenucstdevp=') THEN
        READ(charinput(15:LEN_TRIM(charinput)),*) frenucstdevp
    ELSEIF (charinput(1:14).EQ.'-growthspeedm=') THEN
        READ(charinput(15:LEN_TRIM(charinput)),*) growthspeedm
        growthspeedm=1.0d-9*growthspeedm
    ELSEIF (charinput(1:19).EQ.'-growthspeedstdevm=') THEN
        READ(charinput(20:LEN_TRIM(charinput)),*) growthspeedstdevm
        growthspeedstdevm=1.0d-9*growthspeedstdevm
    ELSEIF (charinput(1:17).EQ.'-shrinkagespeedm=') THEN
        READ(charinput(18:LEN_TRIM(charinput)),*) shrinkagespeedm
        shrinkagespeedm=1.0d-9*shrinkagespeedm
    ELSEIF (charinput(1:22).EQ.'-shrinkagespeedstdevm=') THEN
        READ(charinput(23:LEN_TRIM(charinput)),*) shrinkagespeedstdevm
        shrinkagespeedstdevm=1.0d-9*shrinkagespeedstdevm
    ELSEIF (charinput(1:7).EQ.'-fcatm=') THEN
        READ(charinput(8:LEN_TRIM(charinput)),*) fcatm
    ELSEIF (charinput(1:12).EQ.'-fcatstdevm=') THEN
        READ(charinput(13:LEN_TRIM(charinput)),*) fcatstdevm
    ELSEIF (charinput(1:7).EQ.'-fresm=') THEN
        READ(charinput(8:LEN_TRIM(charinput)),*) fresm
    ELSEIF (charinput(1:12).EQ.'-fresstdevm=') THEN
        READ(charinput(13:LEN_TRIM(charinput)),*) fresstdevm
    ELSEIF (charinput(1:9).EQ.'-frenucm=') THEN
        READ(charinput(10:LEN_TRIM(charinput)),*) frenucm
    ELSEIF (charinput(1:14).EQ.'-frenucstdevm=') THEN
        READ(charinput(15:LEN_TRIM(charinput)),*) frenucstdevm
    ELSEIF (charinput(1:16).EQ.'-terminateatloss') THEN
        terminateatloss=.TRUE.

    ELSE
        WRITE(*,*) "Unknown command line option", charinput
        STOP
    ENDIF
ENDDO


CALL system ('date')

WRITE(*,*) "delta t = ", dt
WRITE(*,*) "Number of MTs = ", noofMTs

CALL system ('rm -rf mtdata')
CALL system ('mkdir mtdata')

write(*,*) "done reading comman line parameters"

!call srand(seedin)
WRITE(*,*) "Seed in", seedin
CALL init_random_seed(seedin,seedout)


OPEN (UNIT=8, FILE="mtdata/totallength.dat", STATUS="REPLACE")

OPEN (UNIT=771, FILE="mtdata/plusgrowth.dat", STATUS="REPLACE")
OPEN (UNIT=772, FILE="mtdata/plusshrinkage.dat", STATUS="REPLACE")

OPEN (UNIT=781, FILE="mtdata/minusgrowth.dat", STATUS="REPLACE")
OPEN (UNIT=782, FILE="mtdata/minusshrinkage.dat", STATUS="REPLACE")

WRITE(*,*) "runtime=", runtime
WRITE(*,*) "maxsteps=", maxsteps


ALLOCATE(mt(noofMTs))

1111  FORMAT(I3.3,'-MTlength.dat') !! time, length, minus tip, plus tip
2222  FORMAT(I3.3,'-MT1D.dat') !! time, length, minus tip, plus tip
3333  FORMAT(I3.3)

! Initialize MT
DO i=1,noofMTs
    !! set length
    prob=-1.0d0
    DO WHILE (prob<0.0d0)
        CALL normal(L0,L0std,seedin, prob)
    ENDDO
    mt(i)%mtlength=prob
    mt(i)%mtlength=dimersize*INT(mt(i)%mtlength/dimersize)

    !! set end positions in 2D
    CALL RANDOM_NUMBER(probtheta)
    CALL RANDOM_NUMBER(probcenterx)
    CALL RANDOM_NUMBER(probcentery)
    mt(i)%theta=2.0d0*PI*probtheta
    mt(i)%minusx=(probcenterx*100.0d-6)-0.0d-6-(DCOS(mt(i)%theta)*mt(i)%mtlength/2.0d0)
    mt(i)%minusy=(probcentery*100.0d-6)-0.0d-6-(DSIN(mt(i)%theta)*mt(i)%mtlength/2.0d0)
    mt(i)%plusx=(probcenterx*100.0d-6)-0.0d-6+(DCOS(mt(i)%theta)*mt(i)%mtlength/2.0d0)
    mt(i)%plusy=(probcentery*100.0d-6)-0.0d-6+(DSIN(mt(i)%theta)*mt(i)%mtlength/2.0d0)

    !! set end positions in 1D
    mt(i)%minusend=0.0d0!-1.0*mt(i)%mtlength/2.0d0
    mt(i)%plusend=mt(i)%minusend+mt(i)%mtlength
    mt(i)%plusinitial=mt(i)%plusend
    mt(i)%minusinitial=mt(i)%minusend

    mt(i)%dimercount=INT(mt(i)%mtlength/dimersize)


!!!! USE 10% ERROR FOR DYNAMIC PARAMETERS
growthspeedstdevp=0.1d0*growthspeedp
shrinkagespeedstdevp=0.1d0*shrinkagespeedp
fcatstdevp=0.1d0*fcatp
fresstdevp=0.1d0*fresp
growthspeedstdevm=0.1d0*growthspeedm
shrinkagespeedstdevm=0.1d0*shrinkagespeedm
fcatstdevm=0.1d0*fcatm
fresstdevm=0.1d0*fresm

    !! set initial DI parameters for plus end
    prob=-1.0d0
    DO WHILE (prob<0.0d0)
        CALL normal(growthspeedp,growthspeedstdevp,seedin, prob)
    ENDDO
    mt(i)%vgp=prob

    prob=-1.0d0
    DO WHILE (prob<0.0d0)
        CALL normal(shrinkagespeedp,shrinkagespeedstdevp,seedin, prob)
    ENDDO
    mt(i)%vsp=prob

    prob=-1.0d0
    DO WHILE (prob<0.0d0)
        CALL normal(fcatp,fcatstdevp,seedin, prob)
    ENDDO
    mt(i)%fcatp=prob

    prob=-1.0d0
    DO WHILE (prob<0.0d0)
        CALL normal(fresp,fresstdevp,seedin, prob)
    ENDDO
    mt(i)%fresp=prob

!! set initial DI parameters for minus end
    prob=-1.0d0
    DO WHILE (prob<0.0d0)
        CALL normal(growthspeedm,growthspeedstdevm,seedin, prob)
    ENDDO
    mt(i)%vgm=prob

    prob=-1.0d0
    DO WHILE (prob<0.0d0)
        CALL normal(shrinkagespeedm,shrinkagespeedstdevm,seedin, prob)
    ENDDO
    mt(i)%vsm=prob

    prob=-1.0d0
    DO WHILE (prob<0.0d0)
        CALL normal(fcatm,fcatstdevm,seedin, prob)
    ENDDO
    mt(i)%fcatm=prob

    prob=-1.0d0
    DO WHILE (prob<0.0d0)
        CALL normal(fresm,fresstdevm,seedin, prob)
    ENDDO
    mt(i)%fresm=prob

    !! initial phase of dynamics
    mt(i)%growingp=.TRUE.
    mt(i)%shrinkingp=.FALSE.

    mt(i)%growingm=.TRUE.
    mt(i)%shrinkingm=.FALSE.


    !! write initial data
    WRITE(buf,1111) i
    buf="mtdata/"//buf(1:16)//char(0)
    OPEN (UNIT=16, FILE=buf, STATUS="REPLACE")
    WRITE(16,*)  0.0d0, mt(i)%minusx, mt(i)%minusy, mt(i)%plusx, mt(i)%plusy, mt(i)%mtlength, mt(i)%theta
    CLOSE(16)

    WRITE(buf,2222) i
    buf="mtdata/"//buf(1:12)//char(0)
    OPEN (UNIT=16, FILE=buf, STATUS="REPLACE")
    WRITE(16,*)  0.0d0, mt(i)%minusend, mt(i)%plusend, mt(i)%mtlength, ((mt(i)%minusend+mt(i)%plusend)/2.0d0)
    CLOSE(16)

ENDDO

! Begin the main loop.
DO WHILE (step.le.maxsteps)

totallength=0.0d0
existingMTcount=0
avepluspos=0.0d0
aveminuspos=0.0d0
DO i=1,noofMTs
IF ((mt(i)%dimercount.GT.1)) THEN
avepluspos=avepluspos+mt(i)%plusend
aveminuspos=aveminuspos+mt(i)%minusend


!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! dynamic instability !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!
    !!! plus end !!!
    !!!!!!!!!!!!!!!!
    ! Growth
    IF (mt(i)%growingp) THEN
    mt(i)%plusgrowthtime=mt(i)%plusgrowthtime+dt
        CALL RANDOM_NUMBER(prob)

        IF ( (mt(i)%vgp.GT.0.0d0) .AND. (prob.LE.(1.0d0-DEXP(-1.0d0*dt*mt(i)%vgp/dimersize))) ) THEN
            mt(i)%plusx=mt(i)%plusx+(DCOS(mt(i)%theta)*dimersize)
            mt(i)%plusy=mt(i)%plusy+(DSIN(mt(i)%theta)*dimersize)
            mt(i)%dimercount=mt(i)%dimercount+1
            mt(i)%mtlength=mt(i)%mtlength+dimersize
            mt(i)%plusend=mt(i)%plusend+dimersize
            mt(i)%plusgrowthlength=mt(i)%plusgrowthlength+dimersize
        ENDIF
    ENDIF

    ! Shrinkage
    IF ((mtseed).AND.((mt(i)%plusend.LE.(mt(i)%plusinitial+1.0d-9)))) mt(i)%shrinkingp=.FALSE.
    IF (mt(i)%dimercount.LE.1) mt(i)%shrinkingp=.FALSE.
    IF(mt(i)%shrinkingp) THEN
        mt(i)%plusshrinkagetime=mt(i)%plusshrinkagetime+dt
        CALL RANDOM_NUMBER(prob)
        vspdummy=mt(i)%vsp
        IF (prob.LE.(1.0d0-DEXP(-1.0d0*dt*vspdummy/dimersize))) THEN
            mt(i)%plusx=mt(i)%plusx-(DCOS(mt(i)%theta)*dimersize)
            mt(i)%plusy=mt(i)%plusy-(DSIN(mt(i)%theta)*dimersize)
            mt(i)%dimercount=mt(i)%dimercount-1
            mt(i)%mtlength=mt(i)%mtlength-dimersize
            mt(i)%plusend=mt(i)%plusend-dimersize
            mt(i)%plusshrinkagelength=mt(i)%plusshrinkagelength+dimersize
        ENDIF
    ENDIF

    ! Catastrophe
    IF (mt(i)%growingp) THEN
        CALL RANDOM_NUMBER(prob)

        fcatpdummy= mt(i)%fcatp
        IF (prob.LE.(fcatpdummy*dt)) THEN
            mt(i)%growingp=.FALSE.
            mt(i)%shrinkingp=.TRUE.

            WRITE(771,*) mt(i)%plusgrowthtime, mt(i)%plusgrowthlength*1.0d9, 1
            mt(i)%plusgrowthtime=0.0d0
            mt(i)%plusgrowthlength=0.0d0

            prob=-1.0d0
            DO WHILE (prob<0.0d0)
                CALL normal(shrinkagespeedp,shrinkagespeedstdevp,seedin, prob)
            ENDDO
            mt(i)%vsp=prob
            prob=-1.0d0
            DO WHILE (prob<0.0d0)
                CALL normal(fresp,fresstdevp,seedin, prob)
            ENDDO
            mt(i)%fresp=prob
        ENDIF
    ENDIF


    ! Rescue
    !IF (mt(i)%shrinkingp) THEN
    IF (.NOT.mt(i)%growingp) THEN
        CALL RANDOM_NUMBER(prob)
        frespdummy= mt(i)%fresp
        IF (prob.LE.(frespdummy*dt)) THEN
            mt(i)%growingp=.TRUE.
            mt(i)%shrinkingp=.FALSE.
            WRITE(772,*) mt(i)%plusshrinkagetime, mt(i)%plusshrinkagelength*1.0d9, 1
            mt(i)%plusshrinkagetime=0.0d0
            mt(i)%plusshrinkagelength=0.0d0
            prob=-1.0d0
            DO WHILE (prob<0.0d0)
                CALL normal(growthspeedp,growthspeedstdevp,seedin, prob)
            ENDDO
            mt(i)%vgp=prob
            prob=-1.0d0
            DO WHILE (prob<0.0d0)
                CALL normal(fcatp,fcatstdevp,seedin, prob)
            ENDDO
            mt(i)%fcatp=prob
        ENDIF
    ENDIF

    !regrow from seeds
    IF ((mtseed).AND.(frenucp.NE.0.0d0).AND.(.NOT.(mt(i)%growingp)).AND.(.NOT.(mt(i)%shrinkingp))) THEN
        CALL RANDOM_NUMBER(prob)
        IF (prob.LE.(dt/frenucp)) THEN
            mt(i)%growingp=.TRUE.
            mt(i)%shrinkingp=.FALSE.
        ENDIF
    ENDIF

    !!!!!!!!!!!!!!!!!
    !!! minus end !!!
    !!!!!!!!!!!!!!!!!
    ! Growth
    IF (mt(i)%growingm) THEN
        mt(i)%minusgrowthtime=mt(i)%minusgrowthtime+dt
        CALL RANDOM_NUMBER(prob)
        IF ( (mt(i)%vgm.GT.0.0d0) .AND. (prob.LE.(1.0d0-DEXP(-1.0d0*dt*mt(i)%vgm/dimersize))) ) THEN
            mt(i)%minusx=mt(i)%minusx-(DCOS(mt(i)%theta)*dimersize)
            mt(i)%minusy=mt(i)%minusy-(DSIN(mt(i)%theta)*dimersize)
            mt(i)%dimercount=mt(i)%dimercount+1
            mt(i)%mtlength=mt(i)%mtlength+dimersize
            mt(i)%minusend=mt(i)%minusend-dimersize
            mt(i)%minusgrowthlength=mt(i)%minusgrowthlength+dimersize
        ENDIF
    ENDIF

    ! Shrinkage
    IF ((mtseed).AND.((mt(i)%minusend.LE.(mt(i)%minusinitial+1.0d-9)))) mt(i)%shrinkingm=.FALSE.
    IF (mt(i)%dimercount.LE.1) mt(i)%shrinkingm=.FALSE.
    IF (mt(i)%shrinkingm) THEN
        mt(i)%minusshrinkagetime=mt(i)%minusshrinkagetime+dt
        CALL RANDOM_NUMBER(prob)
        vsmdummy=mt(i)%vsm
        IF (prob.LE.(1.0d0-DEXP(-1.0d0*dt*vsmdummy/dimersize))) THEN
            mt(i)%minusx=mt(i)%minusx+(DCOS(mt(i)%theta)*dimersize)
            mt(i)%minusy=mt(i)%minusy+(DSIN(mt(i)%theta)*dimersize)
            mt(i)%dimercount=mt(i)%dimercount-1
            mt(i)%mtlength=mt(i)%mtlength-dimersize
            mt(i)%minusend=mt(i)%minusend+dimersize
            mt(i)%minusshrinkagelength=mt(i)%minusshrinkagelength+dimersize
        ENDIF
    ENDIF

    ! Catastrophe
        IF (mt(i)%growingm) THEN
        CALL RANDOM_NUMBER(prob)
        fcatmdummy= mt(i)%fcatm
        IF (prob.LE.(fcatmdummy*dt/1.0d0)) THEN
            mt(i)%growingm=.FALSE.
            mt(i)%shrinkingm=.TRUE.
            WRITE(781,*) mt(i)%minusgrowthtime, mt(i)%minusgrowthlength*1.0d9, 1
            mt(i)%minusgrowthtime=0.0d0
            mt(i)%minusgrowthlength=0.0d0
            prob=-1.0d0
            DO WHILE (prob<0.0d0)
                CALL normal(shrinkagespeedm,shrinkagespeedstdevm,seedin, prob)
            ENDDO
            mt(i)%vsm=prob
            prob=-1.0d0
            DO WHILE (prob<0.0d0)
                CALL normal(fresm,fresstdevm,seedin, prob)
            ENDDO
            mt(i)%fresm=prob
        ENDIF
    ENDIF


    ! Rescue
    IF (.NOT.mt(i)%growingm) THEN
        CALL RANDOM_NUMBER(prob)
        fresmdummy= mt(i)%fresm
        IF (prob.LE.(mt(i)%fresm*dt/1.0d0)) THEN
            mt(i)%growingm=.TRUE.
            mt(i)%shrinkingm=.FALSE.
            WRITE(782,*) mt(i)%minusshrinkagetime, mt(i)%minusshrinkagelength*1.0d9, 1
            mt(i)%minusshrinkagetime=0.0d0
            mt(i)%minusshrinkagelength=0.0d0
            prob=-1.0d0
            DO WHILE (prob<0.0d0)
                CALL normal(growthspeedm,growthspeedstdevm,seedin, prob)
            ENDDO
            mt(i)%vgm=prob
            prob=-1.0d0
            DO WHILE (prob<0.0d0)
                CALL normal(fcatm,fcatstdevm,seedin, prob)
            ENDDO
            mt(i)%fcatm=prob
        ENDIF
    ENDIF

    !regrow from seeds
    IF ((mtseed).AND.(frenucm.NE.0.0d0).AND.(.NOT.(mt(i)%growingm)).AND.(.NOT.(mt(i)%shrinkingm))) THEN
        CALL RANDOM_NUMBER(prob)
        IF (prob.LE.(dt/frenucm)) THEN
            mt(i)%growingm=.TRUE.
            mt(i)%shrinkingm=.FALSE.
        ENDIF
    ENDIF

    IF (mt(i)%dimercount.LE.1) THEN
        mt(i)%growingp=.FALSE.
        mt(i)%shrinkingp=.FALSE.
        mt(i)%growingm=.FALSE.
        mt(i)%shrinkingm=.FALSE.
    ENDIF

    existingMTcount=existingMTcount+1

IF (MOD(step,NINT(2.0d0*NINT(1.0d0/dt))).EQ.0.0d0) THEN
    WRITE(buf,1111) i
    buf="mtdata/"//buf(1:16)//char(0)
    OPEN (UNIT=16, FILE=buf, STATUS="OLD", position="append")
    WRITE(16,*) step*dt, mt(i)%minusx, mt(i)%minusy, mt(i)%plusx, mt(i)%plusy, mt(i)%mtlength, mt(i)%theta
    CLOSE(16)
    WRITE(buf,2222) i
    buf="mtdata/"//buf(1:12)//char(0)
    OPEN (UNIT=16, FILE=buf, STATUS="OLD", position="append")
    WRITE(16,*) step*dt, mt(i)%minusend, mt(i)%plusend, mt(i)%mtlength, ((mt(i)%minusend+mt(i)%plusend)/2.0d0)
    CLOSE(16)
    totallength=totallength+mt(i)%mtlength
ENDIF
ENDIF


IF ( (terminateatloss) .AND. (mt(i)%mtlength.LE.8.0d-9) ) THEN
    EXIT
ENDIF

ENDDO ! noofMTs

IF (MOD(step,NINT(2.0d0*NINT(1.0d0/dt))).EQ.0.0d0) THEN
    IF (existingMTcount.GT.0) THEN
        write(8,*) step*dt, totallength, existingMTcount, totallength/existingMTcount    ! totallength.dat
    ELSE
        write(8,*) step*dt, totallength, existingMTcount, 0.0d0    ! totallength.dat
    ENDIF
ENDIF

IF (MOD(step,(maxsteps/100)).EQ.0) write (*,*) 100*step/maxsteps,"% completed", step*dt, noofMTs
step=step+1



ENDDO ! main (time) loop

DO i=1,noofMTs

if (mt(i)%growingp) then
    WRITE(771,*) mt(i)%plusgrowthtime, mt(i)%plusgrowthlength*1.0d9, 0
endif

if (.NOT.mt(i)%growingp) then
    WRITE(772,*) mt(i)%plusshrinkagetime, mt(i)%plusshrinkagelength*1.0d9, 0 
endif

if (mt(i)%growingm) then
    WRITE(781,*) mt(i)%minusgrowthtime, mt(i)%minusgrowthlength*1.0d9, 0 
endif

if (.NOT.mt(i)%growingm) then
    WRITE(782,*) mt(i)%minusshrinkagetime, mt(i)%minusshrinkagelength*1.0d9, 0 
endif

ENDDO



CLOSE(8)
CLOSE(771)
CLOSE(772)
CLOSE(781)
CLOSE(782)

DEALLOCATE(mt)

CALL system ('date')
STOP

END PROGRAM

SUBROUTINE normal ( mean, stdev, seed, prob )

REAL(8),INTENT(IN) ::mean, stdev
REAL(8),INTENT(OUT) ::prob
REAL(8) ::rone, rtwo,normalab,runiform,xf
REAL(8) , parameter :: rpi = 3.141592653589793D+00
INTEGER ::seed


CALL RANDOM_NUMBER(rone)
CALL RANDOM_NUMBER(rtwo)
x = sqrt ( - 2.0d0 * log ( rone ) ) * cos ( 2.0d0 * rpi * rtwo )

prob = mean + stdev * x

END SUBROUTINE

SUBROUTINE init_random_seed(seedin,seedout)
INTEGER :: i, n, clock
INTEGER, DIMENSION(:), ALLOCATABLE :: seed
INTEGER, INTENT(IN) ::seedin
INTEGER, INTENT(OUT) ::seedout

CALL RANDOM_SEED(size = n)
ALLOCATE(seed(n))

OPEN (UNIT=1, FILE="mtdata/seed.dat", STATUS="REPLACE")

IF (seedin.EQ.0) THEN
CALL SYSTEM_CLOCK(COUNT=clock)
clock=clock+GETPID()
ELSE
clock=seedin
ENDIF
seedout=clock
write(*,*) "clock=",clock
write(*,*) "use clock to repeat the seed"
write(1,*) "clock=",clock
write(1,*) "use clock to repeat the seed"
CLOSE(1)
seed = clock + 37 * (/ (i - 1, i = 1, n) /)
CALL RANDOM_SEED(PUT = seed)

DEALLOCATE(seed)
END SUBROUTINE


