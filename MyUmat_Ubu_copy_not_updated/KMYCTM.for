C     Subroutine to calculate the consistent tangent matrix using the 
C     perturbation method outlined by Miehe and Sun
C
C     15-Aug-2018
C     David Nolan      
C     Adapted from the above by me.
        
      !|Nt:Var names dt start with ijklmn are implied integers.
      !|...others are implied real's (single precision).
      !|::Source : https://www.ibm.com/

      SUBROUTINE KMYCTM(XKIRC1,DFGRD1,NTENS,PROPS,XAMAT,NPT, 
     1                NANISO,DDSDDE,NPROPS,XJ1,NSHR,
     2                G1,G2,XK,XLOFL,XLSTR,P1,P2,P3,P4,
     3                SIGMAX,   
     4                XL1,XL2,XL3,XL4) !,ITER
      	IMPLICIT REAL*8 (A-H,O-Z)
       ! INCLUDE 'ABA_PARAM.INC'  
              DOUBLE PRECISION STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,
     2 TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,NSHR,
     3 PROPS,COORDS,DROT,PNEWDT,CELENT,
     4 DFGRD0,DFGRD1,NOEL,NPT,KSLAY,KSPT,KSTEP,KINC
C
      INTEGER NTENS,NDI,NSTATV,NPROPS,MATERIALNUMBER

        double precision XKIRC1,DF_IJ,DFGRDP,
     1  XKIRCP,CTM_IJ,CTMCOL,ilist,jlist,
     2  XAMAT,XI4,G1,G2,XK,XLOFL,XLSTR!DFGRD1,DDSDDE,PROPS
		
		DIMENSION :: XKIRC1(3,3),DF_IJ(3,3),DFGRDP(3,3),
     1  XKIRCP(3,3),CTM_IJ(3,3),CTMCOL(NTENS),ilist(6),jlist(6),
     2  XAMAT(3,1),XI4(2),DFGRD1(3,3),DDSDDE(NTENS,NTENS),PROPS(NPROPS)


        eps = 1.0e-08 !|Perturbation parameter
        XKIRC1 = XJ1*CAUSTR !|kirchhoffStr(tau) = J * cauchyStr(sigma)
        !| Nt:Material props are already in XAMAT, and type is not used

        !| The Kirchhoff stress is perturbed six times over the
        !| ...i, j'th components of the deformation gradient. This
        !| ...array lists each of the componets to be perturbed over.
        ilist(1) = 1; jlist(1) = 1; !|11 
        ilist(2) = 2; jlist(2) = 2; !|22
        ilist(3) = 3; jlist(3) = 3; !|33
        ilist(4) = 1; jlist(4) = 2; !|12
        ilist(5) = 1; jlist(5) = 3; !|13
        ilist(6) = 2; jlist(6) = 3; !|23 |See doc Abaq Conventions. 
        
        !| CREATE THE PERTURBED DEFORMATION GRADIENT
        Perturbation: do iter = 1,NTENS          
        !|:Loop over each of the six components of the deformation 
        !|...gradient which are to be perturbed, each time calculating 
        !|...a new column of the tangent matrix.      
          !| Pick out the (i,j)th component f d defGrad to be perturbed
          ii = ilist(iter) ; jj = jlist(iter)     
    
          CALL kdelF(ii, jj, DFGRD1, eps, DF_IJ) !|Delta(F)_i_j : DF_IJ
          DFGRDP = DFGRD1 + DF_IJ !|Perturbed defGrad

          CALL kstress_calc(DFGRDP,G1,G2,XK,XLOFL,XLSTR,P1,P2,P3,P4,
     1                      SIGMAX,   !|Muscle extra material param
     2                      XL1,XL2,XL3,XL4, !|Tendon extra mater param
     3                      XI1,XI2,XI3,XI4,XI5,      !|Invariants
     4                      XI1B,XI2B,XI3B,XI4B,XI5B  !|Isochoric invar
     5                      CAUSTP,XAMAT)    !|PertD Cauchy stress SIGMA
          XKIRCHP = XJ1*CAUSTP !|kirch(tau) = J * cauchy(sigma) 

          CTM_IJ = (XKIRCP - XKIRC1)/XJ1/eps !|cauchy(sig)=kirch(tau)/J
          !|:SEE 1996 Miehe eq2.10>so this is its Eulerian form
          !|:SEE 2014 Nolan eqA.4 > ~same as 1996 Miehe. 
     
          !|FFF FORM THE CONSISTENT TANGENT MATRIX 
          CALL kmatrix2vector(CTM_IJ, CTMCOL, NSHR) !|3x3 Tens > 6x1 vec
          !| Insert this vector into the iter'th column of DDSDDE          
          DO insert = 1,NTENS
            DDSDDE(insert,iter) = CTMCOL(insert)
          END DO
          !|FFF END
        END DO Perturbation
      END SUBROUTINE KMYCTM
      


 
      !| DEPOT - Local arrays
C XKIRC1         - True Kirchhoff stress (i.e. not based on perturbed def. grad.)
C DFGRD1          - True deformation gradient for this increment      
C DF_IJ             - Increment of the perturbed deformation gradient
C DFGRDP      - Perturbed def. grad.
C XKIRCP     - Perturbed Kirchhoff stress
C CTM_IJ             - (:,:,I,J) Components of material jacobian
C CTMCOL          - Above in vector form
C ILIST, JLIST    - Set of index components to be perturbed upon
C DUMSTRSS        - DUMMY STRESS TENSOR
C XAMAT           - array containing reference configuration fibre vectors
C XI4             - array storing the fibre invariant (will not be used in this subroutine)
C DDSDDE          - Voigt notation matrix to store the tangent matrix 
C PROPS           - Properties read in from Abaqus .inp file
