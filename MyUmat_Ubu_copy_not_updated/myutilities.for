      SUBROUTINE mykmatrix2vector(XMAT, VEC, NSHR)
      !|:Matrix to Voigt vector
      !|...11 22 33 12 13 23 <SEE doc Abaq Convensions
C         INCLUDE 'ABA_PARAM.INC'
        INTENT(in) :: XMAT, NSHR
        INTENT(out):: VEC
        DIMENSION xmat(3,3), vec(6)
  
        !| Normal components
        DO i=1,3
          vec(i) = xmat(i,i)
        END DO

        !| Shear components
        vec(4) = xmat(1,2)
        IF (NSHR==3) THEN
          vec(5) = xmat(1,3)
          vec(6) = xmat(2,3)
        END IF
      END SUBROUTINE mykmatrix2vector      


      SUBROUTINE mykdelF(m, n, DFGRD, eps, DF) !|Delta(F) : defG perturb
        !| SEE 1996 Miehe eq2.17 for delta(F) (defGrad pertubation calc)
        !| SEE 2014 Nolan eqA.3  > ~same as above w/o g-mericTensor.  
C         INCLUDE 'ABA_PARAM.INC'
        INTENT (in) :: DFGRD, eps, m, n !|Curr DefG, incr, idx to pert
        INTENT (out):: DF !|Perturbed defGrad increment is output. 
        DIMENSION dyad1(3,3),dyad2(3,3),DFGRD(3,3),DF(3,3)  ! ,DFp1(3,3)
          
        !| Zero the dyad matrices
        DO i = 1,3
          DO j = 1,3
            dyad1(i,j) = zero
            dyad2(i,j) = zero
          END DO
        END DO
        !| Place the 1's in the correct location        
        dyad1(m,n) = 1.0;
        dyad2(n,m) = 1.0;

        DF = 0.5D0 *( MATMUL(dyad1,DFGRD) + MATMUL(dyad2,DFGRD) ) * eps
      END SUBROUTINE mykdelF      


      SUBROUTINE myKPRINTER(TENS, m, n) !|Print out a matrix of any size.
C         INCLUDE 'ABA_PARAM.INC'
        INTENT(in):: TENS, m, n      
        DIMENSION tens(m,n)
        
        write(6,*)
        DO i = 1,m
          DO j = 1,n
            ! write(6,'(e19.9)',advance='no'),tens(i,j) 
            !|:Dss orig, but comma not good for some compilers.
            write(6,'(e19.9)',advance='no') tens(i,j)
          END DO
          write(6,*)
        END DO
        write(6,*)
        RETURN
      END SUBROUTINE myKPRINTER      


      SUBROUTINE ACOSH(X,Y)   
      !| INVERSE HYPERBOLIC FUNCTIONS,SEE EXPLANATION by Arnold Desitter
      !| ...computer-programming-forum.com/49-fortran/7ebc82829d81d040.htm
C         INCLUDE 'ABA_PARAM.INC'
C         DOUBLE PRECISION :: X,Y
        ! DOUBLE PRECISION ATANH, ASINH
        INTENT(IN) :: X
        INTENT(OUT):: Y
        ! ATANH(X)=LOG((X+1.)/(-X+1.))/2
        ! ASINH(X)=LOG(X+SQRT(X**2+1.))
        Y = LOG(X+SQRT((X-1.D0)*(X+1.D0))) !|sO says dss problem fr x~=0
        !|:LOG is natural logarithm
      END SUBROUTINE ACOSH


      SUBROUTINE INVER(A,AI)
C         INCLUDE 'ABA_PARAM.INC'
        INTENT(IN)  :: A
        INTENT(OUT) :: AI
        DIMENSION :: A(3,3),AI(3,3)
        !| TAKE THE INVERSE OF A 3x3 MATRIX
        ! Calculate the inverse determinant of the matrix
        dum1 = + A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2) 
        dum2 = - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1) 
        dum3 = + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1)
        detinv = (dum1+dum2+dum3)**(-1)

        ! Calculate the inverse of the matrix
        AI(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
        AI(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
        AI(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
        AI(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
        AI(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
        AI(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
        AI(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
        AI(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
        AI(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))        
      END SUBROUTINE INVER


      SUBROUTINE TRACE(A,TR_A)
C         INCLUDE 'ABA_PARAM.INC'
        DIMENSION :: A(3,3)
        OUTVAL = A(1,1) + A(2,2) + A(3,3)
      END SUBROUTINE TRACE


      SUBROUTINE DEVCUR(CMAT,AB,DEV_AB) !|B for Bar.
        !| DEVCUR : Deviatoric operator in current config, not reference
        !| AB is the isochoric one, B:bar
        !| A is the full one
C         INCLUDE 'ABA_PARAM.INC'
        DIMENSION :: CMAT(3,3),CMATI(3,3),AB(3,3),DEV_AB(3,3)
        CALL TRACE(MATMUL(TRANSPOSE(AB),CMAT),TR_)
        CALL INVER(CMAT,CMATI)
        DEV_AB = AB -1/3*TR_*CMATI !|Deviatoric part 
      END SUBROUTINE






      !| DEPOT
C       subroutine kdotprod(A, B, dotp, n) !|MATMUL() available fortran 90
C       !|:Dot product of two vectors
C         INCLUDE 'ABA_PARAM.INC'
C         INTENT(in) :: A, B, n
C         INTENT(out):: dotp      
C         DIMENSION A(n), B(n)

C         dotp = 0.0
        
C         DO i = 1,n
C           dotp = dotp + A(i)*B(i)
C         END DO
C       end subroutine kdotprod    

C       SUBROUTINE KMTMS (M, N, L, A, KA, B, KB, C, KC) !|MATMUL available
C       !|:Multiply two real matrices 
C         INCLUDE 'ABA_PARAM.INC'  
C         !| Because IMPLICITE NONE is not used here, 
C         !| ...I,J,K,L,M,N (and what starts with them are integers;
C         !| ...and others are real (single precision)
C         INTENT(in) :: M, N, L, A, KA, B, KB, KC
C         !|:So here KA is an implied integer since it starts with a "K".
C         INTENT(out):: C      
C         DIMENSION A(KA,N), B(KB,L), C(KC,L)
C         !|:A(M,N),B(N,L),C(M,L)
C         !|::KA,KB,KC are they even needed?? 
C         DOUBLE PRECISION W

C         DO 30 J = 1,L
C           DO 20 I = 1,M
C             W = 0.D0 !|Dummy var t store C(I,J)=A(I,:).B(:,J) multiplic.
C             DO 10 K = 1,N
C               W = W + A(I,K) * B(K,J)
C    10       CONTINUE  !|C(I,J) el is done
C             C(I,J) = W
C    20     CONTINUE
C    30   CONTINUE
C         RETURN !|You can RETURN if you like, it is assigned to var RETURN.
C         !|:SEE the website ibm.com/support/knowledgecenter/SSGMCP_5.4.0/
C         !|...applications/developing/rexx/subrou.html
C       END SUBROUTINE KMTMS      
