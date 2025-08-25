program prepare_input

implicit none

real :: A(1:69)

integer :: par_f, omni_f, out_f, iostat, nrecs

character(90) :: stormname
character(200) :: in_format, out_format

! TS05 parameters
open(unit=par_f,file='TS05.par')
read(par_f,*) A
close(par_f)

! Read omni file
call get_command_argument(1,stormname)

in_file='./data/omni_data_'//trim(adjustl(stormname))//'_filled_gaps_with_tilt.lst'
out_file='./data/'//trim(adjustl(stormname))//'_TS05_parameters.lst'

in_format = '(2I4,2I3,3F8.2,3F8.1,F7.2,F9.0,F6.2,2I6,2I3,F8.4,F7.2)'
out_format = '(2I4,2I3,3F8.2,3F8.1,F7.2,F9.0,F7.1,2(3X,I2),F8.4,7F7.2)'

nrecs=0
iostat = 0
open(unit=omni_f,file=in_file)
do while ( iostat .eq. 0 )
  read(omni_f,*,iostat=iostat)
  nrecs = nrecs + 1
end do
close(omni_f)

allocate(IYEA(nrecs))
allocate(IDA(nrecs))
allocate(IHOU(nrecs))
allocate(MI(nrecs))
allocate(BXX(nrecs))
allocate(BYY(nrecs))
allocate(BZZ(nrecs))
allocate(VX(nrecs))
allocate(VY(nrecs))
allocate(VZ(nrecs))
allocate(DEN(nrecs))
allocate(T(nrecs))
allocate(PA(nrecs))
allocate(AeIND(nrecs))
allocate(SYMH(nrecs))
allocate(IMFLAG(nrecs))
allocate(SWFLAG(nrecs))
allocate(TILT(nrecs))
allocate(PYDN(nrecs))
allocate(V_SW(nrecs))

open(unit=omni_f,file=in_file)
read(omni_f,*) IYEA,IDA,IHOU,MI,BXX,BYY,BZZ,VX,VY,VZ,DEN,T,PA,AeIND,SYMH,IMFLAG,SWFLAG,TILT,Pdyn
close(omni_f)

V_SW = sqrt(VX**2 + VY**2 + VZ**2)
Pdyn=1.937D-6*DEN*V_SW**2

do IND=1,nrecs

  CALL RECALC_08 (IYEAR(IND),IDAY(IND),IHOUR(IND),MI(IND),0,-400.0,0.0,0.0)

  W1=0.
  W2=0.
  W3=0.
  W4=0.
  W5=0.
  W6=0.

  KEY1=1
  KEY2=1
  KEY3=1
  KEY4=1
  KEY5=1
  KEY6=1

  DO 42 KK=IND,nrecs,-1

  Vnorm=   V_SW(KK)/400.
  Dennorm= DEN(KK)*1.16/5. !  ASSUME HeH=0.04, HENCE 1.16
  Bsnorm= -BZGSM(KK)/5.

                   IF (Bsnorm.LE.0.) THEN
                      Bs1=0.D0
                      Bs2=0.D0
                      Bs3=0.D0
                      Bs4=0.D0
                      Bs5=0.D0
                      Bs6=0.D0
                   ELSE
                      Bs1=Bsnorm**A(53)
                      Bs2=Bsnorm**A(56)
                      Bs3=Bsnorm**A(59)
                      Bs4=Bsnorm**A(62)
                      Bs5=Bsnorm**A(65)
                      Bs6=Bsnorm**A(68)
                   ENDIF

            FAC1=Dennorm**A(51) *Vnorm**A(52) *Bs1
            FAC2=Dennorm**A(54) *Vnorm**A(55) *Bs2
            FAC3=Dennorm**A(57) *Vnorm**A(58) *Bs3
            FAC4=Dennorm**A(60) *Vnorm**A(61) *Bs4
            FAC5=Dennorm**A(63) *Vnorm**A(64) *Bs5
            FAC6=Dennorm**A(66) *Vnorm**A(67) *Bs6

            TAUMT=(IND-KK)*5.

            ARG1=-TAUMT*DT1
            ARG2=-TAUMT*DT2
            ARG3=-TAUMT*DT3
            ARG4=-TAUMT*DT4
            ARG5=-TAUMT*DT5
            ARG6=-TAUMT*DT6

                IF (ARG1.GT.-10..AND.KEY1.EQ.1) THEN
                  W1=W1+FAC1*EXP(ARG1)
                ELSE
                  KEY1=0
                ENDIF

                IF (ARG2.GT.-10..AND.KEY2.EQ.1) THEN
                  W2=W2+FAC2*EXP(ARG2)
                ELSE
                  KEY2=0
                ENDIF

                IF (ARG3.GT.-10..AND.KEY3.EQ.1) THEN
                  W3=W3+FAC3*EXP(ARG3)
                ELSE
                  KEY3=0
                ENDIF

                IF (ARG4.GT.-10..AND.KEY4.EQ.1) THEN
                  W4=W4+FAC4*EXP(ARG4)
                ELSE
                  KEY4=0
                ENDIF

                IF (ARG5.GT.-10..AND.KEY5.EQ.1) THEN
                  W5=W5+FAC5*EXP(ARG5)
                ELSE
                  KEY5=0
                ENDIF

                IF (ARG6.GT.-10..AND.KEY6.EQ.1) THEN
                  W6=W6+FAC6*EXP(ARG6)
                ELSE
                  KEY6=0
                ENDIF


           IF (KEY1.EQ.0.AND.KEY2.EQ.0.AND.KEY3.EQ.0.AND.KEY4.EQ.0.AND.
     +      KEY5.EQ.0.AND.KEY6.EQ.0) GOTO 43


 42         CONTINUE


 43         W1=W1*DT1*5.
            W2=W2*DT2*5.
            W3=W3*DT3*5.
            W4=W4*DT4*5.
            W5=W5*DT5*5.
            W6=W6*DT6*5.




  write(out_f,trim(out_format)) IYEAR,IDAY,IHOUR,MI,BXGSM,TILT,Pdyn,W1,W2,W3,W4,W5,W6
end do



end program prepare_input
