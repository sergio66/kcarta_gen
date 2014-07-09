c Copyright 1997 
c University of Maryland Baltimore County 
c All Rights Reserved

c this is the same as CALXSC.F except that it outputs which of the MXXMOL 
c submatrices the relevant cross-section is stored in (iWhich)
c this is for the KCOMP code, as it considers each gas separately

c simple modification to allow for computation of d/dq,d/dT if
c kJacobian == kSpline == 1 ... this comes in as idQT == 1

c also, to speed up code, this has been split into three blocks that are 
c uncannily similar : if kJacobian = -1
c                     if kJacobian = 0,1 and iGasID .LE. kMaxDQ
c                     if kJacobian = 0,1 and iGasID .GT. kMaxDQ

c NOTE THIS ONLY DOES INTERPOLATION IN TEMPERATURE!!!!!!!!!!!!!!!!!!!
C SO IT DOES NOT CARE WHICH SETS OF LAYERS WE ARE USING

       SUBROUTINE CALXSC(caFName,NFREQ,raFreq,iFstep,iNlay,
     $    raTemp,raXamnt,
     $    NMOL, iaMolid, raaAbs,iWhich, 
     $    idQT, raaDQ, raaDT, iGasID, iDoDQ)

C      Calculate cross-section (XSEC) for a 25 cm-2 region for use
C      with the compressed database etc.

C    INPUT PARAMETERS:
C    type      name      purpose
C    --------  ------    --------------------------
C    CHAR*80   caFName     XSEC binary data file name
C    INTEGER   NFREQ     number of frequency points
C    REAL arr  raFreq    frequency points
C    REAL      iFstep     frequency point spacing
C    INTEGER   iNlay      number of layers
C    REAL arr  raTemp   temperature
C    REAL arr  raXamnt   amount (molecules/cm^2)
C    INTEGER   NMOL      number of XSEC gases desired
C    INT arr   iaMolid   GENLN2 Molecular ID numbers for desired gases

c    INTEGER   idQT    kJacobian=0,1,kSPline=1 => idQT=1 .. calculate d/dq,d/dT
c    INTEGER   iGasID  GENLN2 molcular gas ID of current gas
c    INTEGER   iDODQ   whether or not to compute d/dq for this gas

C    OUTPUT PARAMETERS:
C    type      name    purpose
C    --------  ------  --------------------------
C    REAL arr  XSEC    XSEC output data
C    REAL arr  raaDQ   d/dq (XSEC)
C    REAL arr  raaDT   d/dT (XSEC)

      include '../INCLUDE/kcarta.param'

C      Arguements
       CHARACTER*80 caFName
       INTEGER NFREQ, iNlay, NMOL, iaMolid(MXXMOL), iWhich
       REAL raFreq(kMaxPts), iFstep, raTemp(kProfLayer),
     $    raXamnt(kProfLayer), raaAbs(kMaxPts,kProfLayer),
     $    raaDQ(kMaxPtsJac,kProfLayerJac),
     $    raaDT(kMaxPtsJac,kProfLayerJac)
       INTEGER idQT,iGasID, iDODQ

C      Local Variables
       INTEGER I,J,K,L,iIOUN, IERR, NXSEC, IDXSEC(MXXSEC),NPTS(MXXSEC),
     $    NTEMP(MXXSEC), USEX(MXXSEC), IGAS, INPTS,
     $    INTEMP, ISX, ISXRAW, IEX, IEXRAW, iXsec, ILO, IHI
       REAL FMIN(MXXSEC), FMAX(MXXSEC), XTEMP(MXXTMP), XFMIN,
     $    XFMAX, RAWPTS(MXXTMP,MXXPTS),PTS(MXXTMP,kMaxPts),XJUNK,
     $    DFRAW, FRAWLO, FRAWHI, Q
       REAL dT,XJUNK1

       iIOUN=kCompUnit
       iWhich=0

       OPEN(UNIT=iIOUN,FILE=caFName,STATUS='OLD',FORM='UNFORMATTED',
     $    IOSTAT=IERR)
       IF (IERR .NE. 0) THEN
          WRITE(kStdWarn,1010) IERR, caFName
 1010     FORMAT('Error ',I4,' opening file:',/,A79)
       ENDIF
       kCompUnitOpen=1 

C      Read the datafile header (containing description of data)
       READ(iIOUN) NXSEC
       READ(iIOUN) (IDXSEC(iXsec),iXsec=1,NXSEC)
       READ(iIOUN) (NPTS(iXsec)  ,iXsec=1,NXSEC)
       READ(iIOUN) (NTEMP(iXsec) ,iXsec=1,NXSEC)
       READ(iIOUN) (FMIN(iXsec)  ,iXsec=1,NXSEC)
       READ(iIOUN) (FMAX(iXsec)  ,iXsec=1,NXSEC)

C Decide what XSECs to use for this frequency range
c note that the subroutine called with nmol=1 ie do xsecs individually!!!
c and so iaMolid(1) would be equal to the current GasID
       DO iXsec=1,NXSEC
         USEX(iXsec)=0
         DO J=1,NMOL
           IF (IDXSEC(iXsec) .EQ. iaMolid(J)) THEN
             IF (FMAX(iXsec) .GE. raFreq(1) .AND.
     $          FMIN(iXsec) .LE. raFreq(NFREQ)) THEN
               USEX(iXsec)=J
               iWhich=J
               END IF
             ENDIF
           ENDDO
         ENDDO
c so eventually only ONE of USEX(iI) would be set to > 0

C Loop over the XSEC's in the database. Ignore if it's not
C one of the ones desired.

       DO iXsec=1,NXSEC
          IF (USEX(iXsec) .GT. 0) THEN
C            Read in the mini-header for this particular xsec
             READ(iIOUN) IGAS, INPTS, INTEMP, XFMIN, XFMAX
             IF (IGAS .NE. IDXSEC(iXsec)) THEN
                WRITE(kStdWarn,1050) IGAS, IDXSEC(iXsec)
 1050           FORMAT('Error; reading xsec data for gas ',I2,
     $          ', was expecting gas ',I2)
             ENDIF

C            Note: temperatures in decreasing order (max to min)
             READ(iIOUN) (XTEMP(J),J=1,INTEMP)

C            Read in the data for this xsec
             DO J=1,INTEMP
                READ(iIOUN) (RAWPTS(J,K),K=1,INPTS)
             ENDDO

C            Determine point spacing of raw data
             DFRAW=( XFMAX - XFMIN )/FLOAT( INPTS - 1 )

C            Determine the start points
             XJUNK=raFreq(1) - XFMIN
             IF (XJUNK .GE. 0.0) THEN
               ISX=1
               XJUNK=XJUNK/DFRAW
               ISXRAW=1 + INT(XJUNK)
             ELSE
               XJUNK=-XJUNK/iFstep
               ISX=1 + ( INT(XJUNK)+1 )
               ISXRAW=1
               ENDIF

C            Determine the end points
             XJUNK=XFMAX - raFreq(NFREQ)
             IF (XJUNK .GE. 0.0) THEN
               IEX=NFREQ
               XJUNK=XJUNK/DFRAW
               IEXRAW=INPTS - INT(XJUNK)
             ELSE
               XJUNK=-XJUNK/iFstep
               IEX=NFREQ - ( INT(XJUNK)+1 )
               IEXRAW=INPTS
               ENDIF

             I=ISXRAW
             FRAWLO=XFMIN + FLOAT(I-1)*DFRAW
             FRAWHI=XFMIN + FLOAT(I)*DFRAW

             DO K=ISX,IEX
 10            IF (FRAWHI .LE. raFreq(K) .AND. I .LT. IEXRAW) THEN
                 I=I+1
                 FRAWLO=FRAWHI
                 FRAWHI=XFMIN + (I)*DFRAW
                 GOTO 10
                 ENDIF

               XJUNK=(raFreq(K) - FRAWLO)/DFRAW
               DO J=1,INTEMP
                 PTS(J,K)=XJUNK*(RAWPTS(J,I+1) - RAWPTS(J,I)) +
     $                RAWPTS(J,I)
                 ENDDO
               ENDDO

c ********** this is for kJacobian = -1 ==> idQT = -1  (NO JACOBIANS DONE!)
            IF (idQT .EQ. -1) THEN

C Adjust for temperature and amount
              IF (INTEMP .GT. 1) THEN

c Interpolate xsecs at different temperatures
                DO L=1,iNlay
                  IHI=1
                  ILO=2
C Find which temperatures to use for this layer
 20               IF (raTemp(L) .LT. XTEMP(ILO)
     $                    .AND. ILO .LT. INTEMP) THEN
                    IHI=ILO
                    ILO=ILO + 1
                    GOTO 20
                    ENDIF

                    XJUNK=(raTemp(L) - XTEMP(ILO))/
     $                (XTEMP(IHI) - XTEMP(ILO) )
                    DO K=ISX,IEX
                      raaAbs(K,L)=raXamnt(L)*(
     $                   XJUNK*(PTS(IHI,K) - PTS(ILO,K)) + PTS(ILO,K) )

c check to see if the calculated absorption coeff is > 0 (it could end up
c being negative if the temperature lies OUTSIDE the temp bounds, making
c the interpolation yield a negative value .. if so, use partition fcn 
                      IF (raaAbs(K,L) .LT. 0.0) THEN
C Use partition function to approx temp effects
                        CALL CALQ( IGAS, raTemp(L), Q ,
     $                            idQT, dT)
                        XJUNK1=Q*raXamnt(L)
                        raaAbs(K,L)=XJUNK1*PTS(1,K)
                        END IF

                     ENDDO
                  ENDDO

             ELSE
                DO L=1,iNlay
C Use partition function to approx temp effects
                   CALL CALQ( IGAS, raTemp(L), Q ,
     $                        idQT, dT)

                   XJUNK=Q*raXamnt(L)
                   DO K=ISX,IEX
                      raaAbs(K,L)=XJUNK*PTS(1,K)
                   ENDDO
                ENDDO
             ENDIF
             
c ************** this is for kJacobian = 0,1, iDODQ .GT. 0
c so do d/dT,d/dq
             ELSE IF ((idQT .EQ. 1) .AND. (iDoDQ .GT. 0)) THEN

C Adjust for temperature and amount
              IF (INTEMP .GT. 1) THEN

c Interpolate xsecs at different temperatures
                DO L=1,iNlay
                   IHI=1
                   ILO=2
C Find which temperatures to use for this layer
 21                IF (raTemp(L) .LT. XTEMP(ILO)
     $             .AND. ILO .LT. INTEMP) THEN
                      IHI=ILO
                      ILO=ILO + 1
                      GOTO 21
                   ENDIF

                   XJUNK=(raTemp(L) - XTEMP(ILO))/
     $                (XTEMP(IHI) - XTEMP(ILO) )
                   DO K=ISX,IEX
                      raaAbs(K,L)=raXamnt(L)*(
     $                   XJUNK*(PTS(IHI,K) - PTS(ILO,K)) + PTS(ILO,K) )
                      raaDQ(K,L)=XJUNK*(PTS(IHI,K)-PTS(ILO,K))
     $                              + PTS(ILO,K)
                      raaDT(K,L)=1.0/(XTEMP(IHI) - XTEMP(ILO))*
     $                     raXamnt(L)*(PTS(IHI,K)-PTS(ILO,K))

c check to see if the calculated absorption coeff is > 0 (it could end up
c being negative if the temperature lies OUTSIDE the temp bounds, making
c the interpolation yield a negative value .. if so, use partition fcn 
                      IF (raaAbs(K,L) .LT. 0.0) THEN
C Use partition function to approx temp effects
                        CALL CALQ( IGAS, raTemp(L), Q ,
     $                            idQT, dT)
                        XJUNK1=Q*raXamnt(L)
                        raaAbs(K,L)=XJUNK1*PTS(1,K)
                        raaDQ(K,L)=Q*PTS(1,K)
                        raaDT(K,L)=dT*raXamnt(L)*PTS(1,K)
                        END IF

                   ENDDO
                ENDDO

             ELSE
                DO L=1,iNlay
C Use partition function to approx temp effects
                   CALL CALQ( IGAS, raTemp(L), Q ,
     $                        idQT, dT)

                   XJUNK=Q*raXamnt(L)
                   DO K=ISX,IEX
                      raaAbs(K,L)=XJUNK*PTS(1,K)
                      raaDQ(K,L)=Q*PTS(1,K)
                      raaDT(K,L)=dT*raXamnt(L)*PTS(1,K)
                   ENDDO
                ENDDO
             ENDIF


c ************** this is for kJacobian = 0,1, iDoDQ < 0
c no need to do d/dq
             ELSE IF ((idQT .EQ. 1) .AND. (iDoDQ .LT. 0)) THEN

C Adjust for temperature and amount
              IF (INTEMP .GT. 1) THEN

c Interpolate xsecs at different temperatures
                DO L=1,iNlay
                   IHI=1
                   ILO=2
C Find which temperatures to use for this layer
 22                IF (raTemp(L) .LT. XTEMP(ILO)
     $             .AND. ILO .LT. INTEMP) THEN
                      IHI=ILO
                      ILO=ILO + 1
                      GOTO 22
                   ENDIF

                   XJUNK=(raTemp(L) - XTEMP(ILO))/
     $                (XTEMP(IHI) - XTEMP(ILO) )
                   DO K=ISX,IEX
                      raaAbs(K,L)=raXamnt(L)*(
     $                   XJUNK*(PTS(IHI,K) - PTS(ILO,K)) + PTS(ILO,K) )
                      raaDT(K,L)=1.0/(XTEMP(IHI) - XTEMP(ILO))*
     $                    raXamnt(L)*(PTS(IHI,K)-PTS(ILO,K))

c check to see if the calculated absorption coeff is > 0 (it coild end up
c being negative if the temperature lies OUTSIDE the temp bounds, making
c the interpolation yield a negative value .. if so, use partition fcn 
                      IF (raaAbs(K,L) .LT. 0.0) THEN
C Use partition function to approx temp effects
                        CALL CALQ( IGAS, raTemp(L), Q ,
     $                            idQT, dT)
                        XJUNK1=Q*raXamnt(L)
                        raaAbs(K,L)=XJUNK1*PTS(1,K)
                        raaDT(K,L)=dT*raXamnt(L)*PTS(1,K)
                        END IF
                   ENDDO
                ENDDO

             ELSE
                DO L=1,iNlay
C Use partition function to approx temp effects
                   CALL CALQ( IGAS, raTemp(L), Q ,
     $                        idQT, dT)

                   XJUNK=Q*raXamnt(L)
                   DO K=ISX,IEX
                      raaAbs(K,L)=XJUNK*PTS(1,K)
                      raaDT(K,L)=dT*raXamnt(L)*PTS(1,K)
                   ENDDO
                ENDDO
             ENDIF
           END IF
c ************** end kJacobian options

          ELSE
C Skip over the data
             READ(iIOUN)
             READ(iIOUN)
             DO J=1,NTEMP(iXsec)
                READ(iIOUN)
             ENDDO
          ENDIF

       ENDDO

       CLOSE(iIOUN)
       kCompUnitOpen=-1 

       RETURN
       END
