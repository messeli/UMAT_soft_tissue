      !| MODIFIED FROM THE NOLAN ET AL ABAQ ARTICLE INTO BLEMKER 2005 .
      !| MSA
      !| Nt:In comments b 2000 Holzepfel syntax is used ::: I_4=A.C A

      SUBROUTINE mykstress_calc(DFGRD,XA,G1,G2,XK,XLOFL,XLSTR,
     1                    P1,P2,P3,P4,SIGMAX, 
     2                    XL1,XL2,XL3,XL4,
     3                    XJ,XI1,XI2,XI3,XI4,XI5,      
     4                    XI1B,XI2B,XI3B,XI4B,XI5B,
     5                    CAUSTR,TIME,PROPS) 

      !| call kstress_calc(DFGRD1,XA,G1,G2,XK,XLOFL,XLSTR,
      !|1                    P1,P2,P3,P4,SIGMAX, !|Muscle extra mat para
      !|2                    XL1,XL2,XL3,XL4) !|Tendon extra mater param
      !|3                    XI1,XI2,XI3,XI4,XI5,      !|Invariants
      !|4                    XI1B,XI2B,XI3B,XI4B,XI5B) !|Isochoric invar
      !|5                    CAUSTR,TIME,PROPS) !|Cauchy stress SIGMA

C         INCLUDE 'ABA_PARAM.INC'
        !|:It has this content below in GitHub link:
        !|...KratosLegacyApplications/ConstitutiveLawsApplication/
        !|...custom_external_libraries/umat/ABA_PARAM.INC
C       DOUBLE PRECISION STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
C      1 RPL,DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,
C      2 TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,NSHR,
C      3 PROPS,COORDS,DROT,PNEWDT,CELENT,
C      4 DFGRD0,DFGRD1,NOEL,NPT,KSLAY,KSPT,KSTEP,KINC
C C
C       INTEGER NTENS,NDI,NSTATV,NPROPS,MATERIALNUMBER


        INTENT(IN) :: DFGRD, XA
        INTENT(IN) :: G1,G2,XK,XLOFL,XLSTR
        INTENT(IN) :: P1,P2,P3,P4,SIGMAX
        INTENT(IN) :: XL1,XL2,XL3,XL4
        INTENT(OUT):: XJ,XI1,XI2,XI3,XI4,XI5
        INTENT(OUT):: XI1B,XI2B,XI3B,XI4B,XI5B
        INTENT(OUT):: CAUSTR
        DIMENSION :: TIME(2) !|NEEDED becoz not in ABA_PARAM.INC         
        DIMENSION :: XA(3,1),aCUR(3,1)
        DIMENSION :: DFGRD(3,3),DFGRDB(3,3),DFGRDI(3,3),
     1  A0A(3,3),aoa(3,3),DIAG_A(3,3),
     2  BMAT(3,3),BMAT2(3,3),BMATB(3,3),BMATB2(3,3),
     3  CMAT(3,3),CMAT2(3,3),CMATB(3,3),CMATB2(3,3),CMATI(3,3),
     4  C_A0A(3,3),CB_A0A(3,3),C2_A0A(3,3),CB2A0A(3,3),
     5  DI1BCB(3,3),DI2BCB(3,3),DI3BCB(3,3),DI4BCB(3,3),DI5BCB(3,3),
     7  DW1CB(3,3),S1B(3,3),DEVS1B(3,3),S1(3,3),
     6  DW2CB(3,3),S2B(3,3),DEVS2B(3,3),S2(3,3),
     6  CAUSTR(3,3),
     7  SVOL(3,3),SIG(3,3),SIG1(3,3),SIG2(3,3),SIG3(3,3),SIGVOL(3,3),
     7  EYE(3,3)
       

        EYE = 0.D0
        DO 30 k = 1,3
          EYE(k,k) = 1.D0 
   30   CONTINUE     
        
        DFGRDB = XJ**(-1.D0/3.D0)*DFGRD !|Isochoric defGrad, B for bar
        CALL INVER(DFGRD,DFGRDI) !|Inverse deformation gradient


        BMAT  = MATMUL(DFGRD ,TRANSPOSE(DFGRD) ) !|Left Cauch Green DefT
        BMAT2 = MATMUL(BMAT,BMAT)
        BMATB = MATMUL(DFGRDB,TRANSPOSE(DFGRDB)) !|Left Cauch Greed DefT
        BMATB2= MATMUL(BMATB,BMATB)
C         TRB   = BMAT(1,1)+BMAT(2,2)+BMAT(3,3) 
        CALL TRACE(BMAT,TRB)   
        TRB2  = BMAT2(1,1)+BMAT2(2,2)+BMAT2(3,3)
        TRBB  = BMATB(1,1)+BMATB(2,2)+BMATB(3,3)
        TRBB2 = BMATB2(1,1)+BMATB2(2,2)+BMATB2(3,3)
        !
        CMAT  = MATMUL(TRANSPOSE(DFGRD) ,DFGRD) !|Right Cauch Green DefT
        CMAT2 = MATMUL(CMAT,CMAT)
        CMATB = MATMUL(TRANSPOSE(DFGRDB),DFGRDB) 
        CMATB2= MATMUL(CMATB,CMATB)
        TRC   = CMAT(1,1)+CMAT(2,2)+CMAT(3,3) 
        TRC2  = CMAT2(1,1)+CMAT2(2,2)+CMAT2(3,3)
        TRCB  = CMATB(1,1)+CMATB(2,2)+CMATB(3,3)
        TRCB2 = CMATB2(1,1)+CMATB2(2,2)+CMATB2(3,3)
        CALL INVER(CMAT,CMATI)

        !|AAA> Fibre vector - related calculations
        A0A = MATMUL(XA,TRANSPOSE(XA)) !|Structural Tensor in ref config !TRANSPOSE DOEES NOT WORK with rank 1
C         DO i=1,3 !|IF YOU DEFINE XA WITH RANK 1 : XA(3)
C           DO j=1,3
C             A0A(i,j) = XA(i) * XA(j)
C           END DO
C         END DO            
        
        C_A0A  = MATMUL(CMAT,A0A)
        CB_A0A = MATMUL(CMATB,A0A)
        C2_A0A = MATMUL(CMAT2,A0A)
        CB2A0A = MATMUL(CMATB2,A0A)

        TRCA0A = C_A0A(1,1)+C_A0A(2,2)+C_A0A(3,3)
        TCBA0A = CB_A0A(1,1)+CB_A0A(2,2)+CB_A0A(3,3)
        TC2A0A = C2_A0A(1,1)+C2_A0A(2,2)+C2_A0A(3,3)
        TCB2AA = CB2A0A(1,1)+CB2A0A(2,2)+CB2A0A(3,3)

        aCUR  = MATMUL(DFGRD,XA) !|Fibre vector in current config 
        !|:aCUR is not normalised to 1 !; so it s not aCUR = aCUR/|aCUR|
        aoa = MATMUL(aCUR,TRANSPOSE(aCUR)) !|Structural tensor, cur conf !TRANSPOSE DOEES NOT WORK with rank 1
C         DO i=1,3 !|IF YOU DEFINE aCUR WITH RANK 1 : aCUR(3)
C           DO j=1,3
C             aoa(i,j) = aCUR(i) * aCUR(j)
C           END DO
C         END DO            

        TRaoa = aoa(1,1)+aoa(2,2)+aoa(3,3)
        !|AAA. Fibre vector - related calculations 
  
  
        !|III Invariant calculations
        XJ = +DFGRD(1, 1)*DFGRD(2, 2)*DFGRD(3, 3)
     1       +DFGRD(1, 2)*DFGRD(2, 3)*DFGRD(3, 1)
     2       +DFGRD(1, 3)*DFGRD(2, 1)*DFGRD(3, 2)
     3       -DFGRD(1, 2)*DFGRD(2, 1)*DFGRD(3, 3)
     4       -DFGRD(1, 3)*DFGRD(2, 2)*DFGRD(3, 1)
     5       -DFGRD(1, 1)*DFGRD(2, 3)*DFGRD(3, 2)  !|Jacobian J=det(F)

        XJ_13 = XJ**(-1.D0/3.D0) 
        XJ_23 = XJ**(-2.D0/3.D0)
        XJ_43 = XJ**(-4.D0/3.D0) 

        XI1  = TRC
C         XI1  = TRB
C         XI1B = TRCB
C         XI1B = TRBB
        XI1B = XJ_13*XI1 

        XI2  = 0.5D0*(TRC**2.0-TRC2) 
C         XI2  = 0.5D0*(TRB**2.0-TRB2) 
C         XI2B = 0.5D0*(TRCB**2.0-TRCB2) 
C         XI2B = 0.5D0*(TRBB**2.0-TRBB2) 
        XI2B = XJ_43*XI2

        XI3 = +CMAT(1, 1)*CMAT(2, 2)*CMAT(3, 3)
     1        +CMAT(1, 2)*CMAT(2, 3)*CMAT(3, 1)
     2        +CMAT(1, 3)*CMAT(2, 1)*CMAT(3, 2)
     3        -CMAT(1, 3)*CMAT(3, 1)*CMAT(2, 2)
     4        -CMAT(1, 2)*CMAT(2, 1)*CMAT(3, 3)
     5        -CMAT(1, 1)*CMAT(3, 2)*CMAT(2, 3)
C         XI3 = +BMAT(1, 1)*BMAT(2, 2)*BMAT(3, 3)
C      1        +BMAT(1, 2)*BMAT(2, 3)*BMAT(3, 1)
C      2        +BMAT(1, 3)*BMAT(2, 1)*BMAT(3, 2)
C      3        -BMAT(1, 3)*BMAT(2, 2)*BMAT(3, 1)
C      4        -BMAT(1, 2)*BMAT(2, 1)*BMAT(3, 3)
C      5        -BMAT(1, 1)*BMAT(3, 2)*BMAT(2, 3)
C         XI3 = XJ**2
C         XI3B= +CMATB(1,1)*CMATB(2,2)*CMATB(3,3)
C      1        +CMATB(1,2)*CMATB(2,3)*CMATB(3,1)
C      2        +CMATB(1,3)*CMATB(2,1)*CMATB(3,2)
C      3        -CMATB(1,3)*CMATB(2,2)*CMATB(3,1)
C      4        -CMATB(1,2)*CMATB(2,1)*CMATB(3,3)
C      5        -CMATB(1,1)*CMATB(3,2)*CMATB(2,3) !| = 1.D0
C         XI3B= +BMATB(1,1)*BMATB(2,2)*BMATB(3,3)
C      1        +BMATB(1,2)*BMATB(2,3)*BMATB(3,1)
C      2        +BMATB(1,3)*BMATB(2,1)*BMATB(3,2)
C      3        -BMATB(1,3)*BMATB(2,2)*BMATB(3,1)
C      4        -BMATB(1,2)*BMATB(2,1)*BMATB(3,3)
C      5        -BMATB(1,1)*BMATB(2,3)*BMATB(3,2) !| = 1.D0
        XI3B= 1.D0 !|since isochoric : J=1. So J^2=0

        XI4  = TRCA0A !|I_4 = tr(C A0A) = A.C A  (:= A_i C_ij A_j) 
C         XI4  = TRaoa  !|I_4 = A.F^T F A = a o a 
C         XI4B = TCBA0A !|I_4_bar = tr(C_bar A0A) = A.C_bar A 
        XI4B = XJ_23*I4 

        XI5  = TC2A0A !|I_5 = tr(C^2 A0A) = A.C^2 A  (:= A_i C^2_ij A_j) 
C         XI5B = TCB2AA !|I_5_bar = tr(C_bar^2 A0A) = A.C_bar^2 A 
        XI5B = XJ_43*I5

        XLAM = SQRT(XI4)

        LAMBDB = XJ_13*XLAM
        !|III END Invariant calculations


        !|SSS> STRESS 
        DI1BCB = EYE
        DI2BCB = XI1B*EYE - CMATB
        DI3BCB = 0.0D0
        DI4BCB = A0A
        DI5BCB = MATMUL(XA,TRANSPOSE(MATMUL(CMATB,XA))) 
     +          +MATMUL(MATMUL(TRANSPOSE(CMATB),XA),TRANSPOSE(XA))
C    +          +MATMUL(MATMUL(CMATB,XA),TRANSPOSE(XA)) !|CB symmetric
        !|:??
        !|: Ni.Cjm.Nm  is_elements_of  [N] . [[C] . [N]]^T 
        !|: Nm.Cmi.Nj  is_elements_of  [C]^T . [N] . [N]^T

        !|S_111> SIG1
        DW1I4B = -G1*XI5B/(XI4B**2)
        DW1I5B = +G1/XI4B
        DW1CB = DW1I4B*DI4BCB + DW1I5B*DI5BCB
        S1B = 2.D0*DW1CB
        CALL DEVCUR(CMAT,S1B,DEVS1B)
        S1 = XJ_23*DEVS1B
        SIG1 = XJ**(-1)*DFGRD*S1*TRANSPOSE(DFGRD)
        !|S_111. SIG1
        
        !|S_222> SIG2
        C_W2 = (XI1B*XI4B-XI5B)/(2*SQRT(XI4B))
        CALL ACOSH(C_W2,ACOSH1) 
        C_W2_2 = 2*G2*ACOSH1/SQRT(C_W2**2 - 1)
        DW2I1B = +C_W2_2*(SQRT(XI4B)/2)
        DW2I4B = +C_W2_2*(XI1B*XI4B+XI5B)/(2*XI4B*(3/2))
        DW2I5B = -C_W2_2*(SQRT(XI4B**(3/2)))/2
        DW2CB = DW2I1B*DI1BCB + DW2I4B*DI4BCB + DW2I5B*DI5BCB
        S2B = 2.D0*DW2CB
        CALL DEVCUR(CMAT,S2B,DEVS2B)
        S2 = XJ_23*DEVS2B
        SIG2 = XJ**(-1)*DFGRD*S2*DFGRDI
        !|S_222. SIG2

        !|S_333 SIG3
        !| Piece-wise stretch dependent parameters, Page b15
        IF (PROPS(1).LT.0.0D0) THEN  !|MUSCLE    
          !| MUSCLE - Passive part
          IF (XLAM.LE.XLOFL) THEN 
            FFPAS = 0.D0 !|f_fibre_passive
          ELSE IF (XLAM.GT.XLOFL .AND. XLAM.LT.XLSTR) THEN
            FFPAS = P1*( exp( P2*(XLAM/XLOFL-1) ) - 1 )
          ELSE IF (XLAM.GE.XLSTR) THEN 
            FFPAS = P3*XLAM/XLOFL + P4
          END IF
          
          !| MUSCLE - Active part
          IF (XLAM.LE.(0.6D0*XLOFL)) THEN
            FFACT = 9*(XLAM/XLOFL-0.4D0)**2 !|f_fibre_active
          ELSE IF (XLAM.GT.(0.6D0*XLOFL).AND.XLAM.LT.(1.4D0*XLOFL)) THEN
            FFACT = 1-4*(1-XLAM/XLOFL)**2
          ELSE IF (XLAM.GE.(1.4D0*XLOFL)) THEN 
            FFACT =  9*(XLAM/XLOFL-1.6D0)**2
          END IF

          !| Activation level: TONUS is %1 of total activation  
          ALPHA = 0.01 !|Preset alpha, to be tried first.
          ! ALPHA = TIME(1)*0.01 !|Gradual alpha
          !|:See initial stress in Abaq doc, whether it is useful, 
          !|...instead of tonus definition with alpha.

          FFTOT = FFPAS + ALPHA*FFACT  !|f_fibre_total

          !| Total muscle fibre stress value, IN FIBRE DIRECTION
          SIGFIB = SIGMAX*FFTOT*(XLAM/XLOFL)
       
        ELSE IF (PROPS(1).LT.0.0D0) THEN !|TENDON 
          !| TENDON - direct total response
          IF (XLAM.LE.1.D0) THEN 
            SIGTEN = 0.0D0
          ELSE IF (XLAM.GT.1.0D0.AND.XLAM.LT.XLSTR) THEN 
            SIGTEN = XL1*( EXP(XL2*(XLAM-1)) - 1 )
          ELSE IF (XLAM.GE.XLSTR) THEN
            SIGTEN = XL3*XLAM+XL4
          END IF 
        END IF 

        !| Make the fibre stress value a tensor
        DIAG_A = 0.0D0
        DIAG_A(1,1) = XA(1,1) 
        DIAG_A(2,2) = XA(2,1)
        DIAG_A(3,3) = XA(3,1)
        IF (PROPS(1).GT.0.0D0) THEN !|MUSCLE 
          SIG3 = SIGFIB*DIAG_A
        ELSE IF (PROPS(1).LT.0.0D0) THEN !|TENDON
          SIG3 = SIGTEN*DIAG_A
        END IF
        !|333. SIG3
        
        !|S_VVV> SIGVOL
        PRESS = XK*LOG(XJ)/XJ
        SVOL = XJ*PRESS*CMATI
        SIGVOL = XJ**(-1)*DFGRD*SVOL*DFGRDI
        !|S_VVV. 

        !| TOTAL SIGMA 
        SIG = SIG1+SIG2+SIG3+SIGVOL
        CAUSTR = SIG
        !|SSS. END STRESS
        RETURN
      END SUBROUTINE mykstress_calc


      INCLUDE 'myutilities.for'
      

      !| DEPOT

      !|PARAMETER(ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0)

      !| LOCAL ARRAYS
      !| DFGRD     - the deformation gradient 
      !| BMAT      - left Cauchy-Green deformation tensor
      !| XKIRCH    - Kirchhoff stress
      !| XA     - array containing fibre vectors
      !| B2MAT     - square of the B tensor
      !| XA        - fibre vector in the reference configuration
      !| a_curnt   - fibre vector in the current configuration
      !| aoa       - fibre structural tensor
      !| XANISOK   - anisotropic Kirchhoff stress for a given fibre family
      !| XANISOTOT - total anisotropic Kirchhoff stress
      !| XI4       - array storing the fibre invariant      
