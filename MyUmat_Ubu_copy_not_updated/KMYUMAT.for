      !| MODIFIED FROM THE NOLAN ET AL ABAQ ARTICLE INTO BLEMKER 2005 .
      !| MSA
      !| UMAT DOC http://130.149.89.49:2080/v2016/books/sub/default.htm
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
      !|:
      !|...TIME(1)  Value f step time at d beginJ f d current increment
      !|...TIME(2)  Value f total time at d beginJ f d current increment
      !|...COORDS  currCoords of ds point (Finite deform),not initCoords
      !|...DTIME  Time increment
      !|...NOEL  El number
      !|...NPT  Integration point number
      !|...JSTEP  Step number
      !|...SSE  Specific elastic strain energy
        IMPLICIT REAL*8 (A-H,O-Z)
C         INCLUDE 'ABA_PARAM.INC'
        !|:It has this content below in GitHub link:
        !|...KratosLegacyApplications/ConstitutiveLawsApplication/
        !|...custom_external_libraries/umat/ABA_PARAM.INC
      DOUBLE PRECISION STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,
     2 TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,NSHR,
     3 PROPS,COORDS,DROT,PNEWDT,CELENT,
     4 DFGRD0,DFGRD1,NOEL,NPT,KSLAY,KSPT,KSTEP,KINC
C
      INTEGER NTENS,NDI,NSTATV,NPROPS,MATERIALNUMBER

        CHARACTER*80 CMNAME
        DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1  DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2  STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3  PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4  JSTEP(4) !|Added from abaq 2016 docum


        !| Nt:From UMAT doc of Abaq2016:
        !|...user coding to define DDSDDE, STRESS, STATEV, SSE, SPD, SCD
        !|...and, if necessary, RPL, DDSDDT, DRPLDE, DRPLDT, PNEWDT
        
        double precision CAUSTR,XKIRC1, XA, XIF , G1, G2, XK, XLOFL,
     1  XLSTR
        DIMENSION :: CAUSTR(3,3),XKIRC1(3,3), XA(3,1), XIF(2)
 
        NANISO = 1 !|Drs only 1 fibre family; so transversely isotropic.

        !|XXX> Blemker material properties
        !| Fibre vector 
        XA(1,1) = PROPS(1)
        XA(2,1) = PROPS(2)
        XA(3,1) = PROPS(3)
        XA = XA/(XA(1,1)**2+XA(2,1)**2+XA(3,1)**2)**0.5D0
        !|:PROPS(1:3) are the direction conponents. 
        !|...Not yet decided how to obtain
        !|...Piecewise with manual entry is an option.

        TYPE = PROPS(4) !|Either >0 for muscle or <0 for tendon,:::+1&-1
        
        IF (TYPE.GT.0.0D0) THEN 
        !| Muscle values 
        G1    = 5.00D2 !|G1_muscle , Pa
        G2    = 5.00D2 !|G2_muscle , Pa
        XK    = 1.00D7 !|K_muscle  , Pa
        XLOFL = 1.40D0 !|Stretch for optimum force 
        XLSTR = 1.40D0 !|When passive muscle fibre force becomes linear
        !| Muscle specific
        P1    = 0.05D0 !|A material constant
        P2    = 6.60D0 !|A material constant
        P3    = P1*P2*exp(P2*(XLSTR/XLOFL)-1)
        P4    = P1*(exp(P2*(XLSTR/XLOFL-1))-1) - P3*XLSTR/XLOFL
        SIGMAX = 3.00D5 !|sig_max , Pa , only muscle
        !| Tendon specific
        XL1    = 2.70D6 !|Pa, only tendon
        XL2    = 46.4D0 !|A material constant
        XL3    = XL1*XL2*exp(XL2*(XLSTR-1))
        XL4    = XL1*( exp(XL2*(XLSTR-1)) - 1 ) - XL3*XLSTR

        ELSE IF (TYPE.LT.0.0D0) THEN 
        !| Tendon values
        G1    = 5.00D4 !|G1_tendon , Pa
        G2    = 5.00D4 !|G2_tendon , Pa
        XK    = 1.00D8 !|K_tendon  , Pa
        XLOFL = 1.40D0 !|Stretch for optimum force 
        XLSTR = 1.03D0 !|When passive muscle fibre force becomes linear
        !| Muscle specific
        P1    = 0.05D0 !|A material constant
        P2    = 6.60D0 !|A material constant
        P3    = P1*P2*exp(P2*(XLSTR/XLOFL)-1)
        P4    = P1*(exp(P2*(XLSTR/XLOFL-1))-1) - P3*XLSTR/XLOFL
        SIGMAX= 3.00D5 !|sig_max , Pa , only muscle       
        !| Tendon specific
        XL1    = 2.70D6 !|Pa, only tendon
        XL2    = 46.4D0 !|A material constant
        XL3    = XL1*XL2*exp(XL2*(XLSTR-1))
        XL4    = XL1*( exp(XL2*(XLSTR-1)) - 1 ) - XL3*XLSTR
        END IF
        !|XXX. Blemker material properties
        !|:Nt:Couldve been defined in Props & read automaticly from .inp


        !|YYY> RETURN to the solver what it needs. CALL CALL
        !| CALC STRESS 
        CALL KMYstress_calc(DFGRD1,XA,G1,G2,XK,XLOFL,XLSTR,P1,P2,P3,P4,
     1                      SIGMAX,   !|Muscle extra material param
     2                      XL1,XL2,XL3,XL4, !|Tendon extra mater param
     3                      XI1,XI2,XI3,XI4,XI5,      !|Invariants
     4                      XI1B,XI2B,XI3B,XI4B,XI5B,  !|Isochoric invar
     5                      CAUSTR,TIME,PROPS)    !|Cauchy stress SIGMA
        STATEV(1) = XI1
        STATEV(2) = XI2
        STATEV(3) = XI3
        STATEV(4) = XI4
        STATEV(5) = XI5
        STATEV(6) = XI1B
        STATEV(7) = XI2B
        STATEV(8) = XI3B !|Should always equal to "1.D0"
        STATEV(9) = XI4B
        STATEV(10)= XI5B

        !| Voigt vector of stress 
        CALL KMYMATRIX2VECTOR(CAUSTR, STRESS, nshr)         
        !|STRESS here is in 6x1 vector form

        !| Tangent matrix calculation (:Material gradient/Jacobian)
        DO I=1,NTENS
          DO J=1,NTENS
            DDSDDE(I,J) = 0.D0
          END DO
        END DO          
  
        CALL KMYCTM(XKIRC1,DFGRD1,NTENS,PROPS,XA,NPT, !,ITER
     1              NANISO,DDSDDE,NPROPS,XJ1,NSHR,
     2              G1,G2,XK,XLOFL,XLSTR,P1,P2,P3,P4,
     3              SIGMAX,   
     4              XL1,XL2,XL3,XL4)
        !|YYY. 
        RETURN
      ENDSUBROUTINE UMAT

      
      !| Called Subroutines
      INCLUDE 'KMYSTRESS_CALC.for' !|Subroutine to calc Cauchy stress
      INCLUDE 'KMYPOLARDECOMP.for' !|Polar decomposition 
      INCLUDE 'KMYCTM.for' !|Subroutine to calculate the tangent matrix
      INCLUDE 'KMYUTILITIES.for' !|Helper subroutines      
      

      !| DEPOT

        !| Convert degrees to radians - IS NEEDED ?
C         XPI = 4*ATAN(1.0D0) !|One 4th the angle whose tangent is 1.0D0
C         THETAR = THETAD*XPI/180.D0
  
C         iter=0  !|?? Not used as far as I see in kstress_calc.for
