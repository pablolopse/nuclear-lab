C=======================================================================
C     Fits a single Gaussian peak on a linear background in a
C     channel-vs-counts spectrum, using the "log-parabola" trick:
C     ln(Gaussian) is a parabola, so the peak parameters (centroid,
C     sigma, height, area) come out of a *linear* least-squares fit
C     instead of a nonlinear solver.
C
C     Steps:
C       1. Read the spectrum (channel, counts columns).
C       2. Fit a linear background through two bands flanking the peak.
C       3. Subtract that background inside the peak window, take
C          log(net counts), and fit a parabola to it.
C       4. Recover Gaussian parameters analytically from the parabola
C          coefficients, propagating the background fit's covariance
C          into the log-counts uncertainty.
C
C     Reads: ../data/Mg4-NaI-background-subtracted.txt (two-column
C     spectrum, produced by analysis.ipynb).
C=======================================================================
      PROGRAM MAIN
      IMPLICIT NONE

      INTEGER NCHAN, MA, MBKG
      PARAMETER (NCHAN=512, MA=3, MBKG=2)

      REAL X(NCHAN), Y(NCHAN), SIG(NCHAN)
      REAL XBKG(NCHAN), YBKG(NCHAN), SBKG(NCHAN)
      REAL XPEAK(NCHAN), YPEAK(NCHAN), SPEAK(NCHAN)
      REAL BKG(MBKG), PEAK(MA), COVAR(MA,MA)
      REAL CHISQ, VARF, NET, CENTER, SIGMA, HEIGHT, AREA, PI, X0, T
      INTEGER LISTA(MA)
      INTEGER i, j, k, nbkg, npeak, NPTS, IOS
      INTEGER CHANMIN, CHANMAX, MARGIN

      PI = 4.0*atan(1.0)

C     === CHOOSE THE PEAK: only these two numbers (in channels) ===
      CHANMIN = 14
      CHANMAX = 22
C     (K-40, 1460.8 keV, in background-HPGe.txt)

C     background margin on each side (automatic, no need to touch)
      MARGIN = CHANMAX - CHANMIN

C     reference center (avoids losing precision in single-precision REAL)
      X0 = 0.5*(CHANMIN+CHANMAX)

C     LISTA lists which of the MA parabola coefficients are free; all
C     of them here. Required by the LFIT routine below (Numerical
C     Recipes naming, kept as published).
      DO i=1,MA
        LISTA(i)=i
      END DO

C     --- Read the spectrum (channel, counts) to end of file ---
      open(10, File="../data/Mg4-NaI-background-subtracted.txt")
      NPTS=0
      DO i=1,NCHAN
        read(10,*,iostat=IOS) X(i), Y(i)
        if (IOS.ne.0) goto 5
        if (Y(i).gt.1.0) then
          SIG(i)=sqrt(Y(i))
        else
          SIG(i)=1.0
        end if
        NPTS=NPTS+1
      END DO
 5    close(10)

C     --- Linear background using the bands right next to the peak ---
      nbkg=0
      DO i=1,NPTS
        if ((X(i).ge.CHANMIN-MARGIN .and. X(i).lt.CHANMIN) .or.
     &      (X(i).gt.CHANMAX .and. X(i).le.CHANMAX+MARGIN)) then
          nbkg=nbkg+1
          XBKG(nbkg)=X(i)-X0
          YBKG(nbkg)=Y(i)
          SBKG(nbkg)=SIG(i)
        end if
      END DO
      CALL LFIT(XBKG,YBKG,SBKG,nbkg,BKG,MBKG,LISTA,COVAR,MBKG,MA,
     &          CHISQ)

C     --- log(counts - background) in the peak window, parabola fit ---
      npeak=0
      DO i=1,NPTS
        if (X(i).ge.CHANMIN .and. X(i).le.CHANMAX) then
          T=X(i)-X0
          NET=Y(i)-BKG(1)-BKG(2)*T
          if (NET.gt.0.0) then
            VARF=0.0
            DO j=1,MBKG
              DO k=1,MBKG
                VARF=VARF+COVAR(j,k)*(T**(j-1))*(T**(k-1))
              END DO
            END DO
            npeak=npeak+1
            XPEAK(npeak)=T
            YPEAK(npeak)=log(NET)
            SPEAK(npeak)=sqrt(SIG(i)**2+VARF)/NET
          end if
        end if
      END DO
      CALL LFIT(XPEAK,YPEAK,SPEAK,npeak,PEAK,MA,LISTA,COVAR,MA,MA,
     &          CHISQ)

      if (PEAK(3).ge.0.0) then
        print*, "WARNING: curvature >= 0; window is not over a peak"
      end if

      CENTER = X0 - PEAK(2)/(2.0*PEAK(3))
      SIGMA  = sqrt(-1.0/(2.0*PEAK(3)))
      HEIGHT = exp(PEAK(1)-PEAK(2)**2/(4.0*PEAK(3)))
      AREA   = sqrt(2.0*PI)*HEIGHT*SIGMA

      print*, "Gaussian peak fit (log-parabola method)"
      print*, "Background points:", nbkg
      print*, "Peak points      :", npeak
      print*, "Center channel   :", CENTER
      print*, "Sigma            :", SIGMA
      print*, "FWHM (channels)  :", 2.3548*SIGMA
      print*, "Height           :", HEIGHT
      print*, "Area             :", AREA

      END PROGRAM


      SUBROUTINE FUNCS(x,AFUNC,MA)
C     Basis functions for the parabola fit: 1, x, x^2.
      INTEGER i
      DIMENSION AFUNC(MA)

      DO i=1,MA
        AFUNC(i) = x**(i-1)
      END DO

      END SUBROUTINE


C     LFIT, COVSRT and GAUSSJ below are adapted from "Numerical Recipes
C     in Fortran 77" (Press, Teukolsky, Vetterling & Flannery), not
C     original code. Variable names (e.g. LISTA) follow the book.

c
      SUBROUTINE LFIT(xf,YF,SIG,NF,Af,MA,LISTA,COVAR,MFIT,NCVM,CHISQ)
c
      DIMENSION xF(NF),YF(NF),SIG(NF),Af(MA),LISTA(MFIT),
     *COVAR(NCVM,NCVM),BETA(100),AFUNC(100)
      KK=MFIT+1
      DO 12 J=1,MA
        IHIT=0
        DO 11 K=1,MFIT
          IF (LISTA(K).EQ.J) IHIT=IHIT+1
11      CONTINUE
        IF (IHIT.EQ.0) THEN
          LISTA(KK)=J
          KK=KK+1
        ELSE IF (IHIT.GT.1) THEN
          PAUSE 'Improper set in LISTA'
        ENDIF
12    CONTINUE
      IF (KK.NE.(MA+1)) PAUSE 'Improper set in LISTA'
      DO 14 J=1,MFIT
        DO 13 K=1,MFIT
          COVAR(J,K)=0.
13      CONTINUE
        BETA(J)=0.
14    CONTINUE
      DO 18 I=1,nf
        CALL FUNCS(xF(I),AFUNC,MA)
        YM=YF(I)
        IF(MFIT.LT.MA) THEN
          DO 15 J=MFIT+1,MA
            YM=YM-Af(LISTA(J))*AFUNC(LISTA(J))
15        CONTINUE
        ENDIF
        SIG2I=1./SIG(I)**2
        DO 17 J=1,MFIT
          WT=AFUNC(LISTA(J))*SIG2I
          DO 16 K=1,J
            COVAR(J,K)=COVAR(J,K)+WT*AFUNC(LISTA(K))
16        CONTINUE
          BETA(J)=BETA(J)+YM*WT
17      CONTINUE
18    CONTINUE
      IF (MFIT.GT.1) THEN
        DO 21 J=2,MFIT
          DO 19 K=1,J-1
            COVAR(K,J)=COVAR(J,K)
19        CONTINUE
21      CONTINUE
      ENDIF
      CALL GAUSSJ(COVAR,MFIT,NCVM,BETA,1,1)
      DO 22 J=1,MFIT
        Af(LISTA(J))=BETA(J)
22    CONTINUE
      CHISQ=0.
      DO 24 I=1,NF
        CALL FUNCS(xf(I),AFUNC,MA)
        SUM=0.
        DO 23 J=1,MA
          SUM=SUM+Af(J)*AFUNC(J)
23      CONTINUE
        CHISQ=CHISQ+((Yf(I)-SUM)/SIG(I))**2
24    CONTINUE
      CALL COVSRT(COVAR,NCVM,MA,LISTA,MFIT)
      RETURN
      END
      SUBROUTINE COVSRT(COVAR,NCVM,MA,LISTA,MFIT)
      DIMENSION COVAR(NCVM,NCVM),LISTA(MFIT)
      DO 12 J=1,MA-1
        DO 11 I=J+1,MA
          COVAR(I,J)=0.
11      CONTINUE
12    CONTINUE
      DO 14 I=1,MFIT-1
        DO 13 J=I+1,MFIT
          IF(LISTA(J).GT.LISTA(I)) THEN
            COVAR(LISTA(J),LISTA(I))=COVAR(I,J)
          ELSE
            COVAR(LISTA(I),LISTA(J))=COVAR(I,J)
          ENDIF
13      CONTINUE
14    CONTINUE
      SWAP=COVAR(1,1)
      DO 15 J=1,MA
        COVAR(1,J)=COVAR(J,J)
        COVAR(J,J)=0.
15    CONTINUE
      COVAR(LISTA(1),LISTA(1))=SWAP
      DO 16 J=2,MFIT
        COVAR(LISTA(J),LISTA(J))=COVAR(1,J)
16    CONTINUE
      DO 18 J=2,MA
        DO 17 I=1,J-1
          COVAR(I,J)=COVAR(J,I)
17      CONTINUE
18    CONTINUE
      RETURN
      END
      SUBROUTINE GAUSSJ(A,N,NP,B,M,MP)
      PARAMETER (NMAX=100)
      DIMENSION A(NP,NP),B(NP,MP),IPIV(NMAX),INDXR(NMAX),INDXC(NMAX)
      DO 11 J=1,N
        IPIV(J)=0
11    CONTINUE
      DO 22 I=1,N
        BIG=0.
        DO 13 J=1,N
          IF(IPIV(J).NE.1)THEN
            DO 12 K=1,N
              IF (IPIV(K).EQ.0) THEN
                IF (ABS(A(J,K)).GE.BIG)THEN
                  BIG=ABS(A(J,K))
                  IROW=J
                  ICOL=K
                ENDIF
              ELSE IF (IPIV(K).GT.1) THEN
                PAUSE '1,Singular matrix'
              ENDIF
12          CONTINUE
          ENDIF
13      CONTINUE
        IPIV(ICOL)=IPIV(ICOL)+1
        IF (IROW.NE.ICOL) THEN
          DO 14 L=1,N
            DUM=A(IROW,L)
            A(IROW,L)=A(ICOL,L)
            A(ICOL,L)=DUM
14        CONTINUE
          DO 15 L=1,M
            DUM=B(IROW,L)
            B(IROW,L)=B(ICOL,L)
            B(ICOL,L)=DUM
15        CONTINUE
        ENDIF
        INDXR(I)=IROW
        INDXC(I)=ICOL
        IF (A(ICOL,ICOL).EQ.0.) PAUSE '2,Singular matrix.'
        PIVINV=1./A(ICOL,ICOL)
        A(ICOL,ICOL)=1.
        DO 16 L=1,N
          A(ICOL,L)=A(ICOL,L)*PIVINV
16      CONTINUE
        DO 17 L=1,M
          B(ICOL,L)=B(ICOL,L)*PIVINV
17      CONTINUE
        DO 21 LL=1,N
          IF(LL.NE.ICOL)THEN
            DUM=A(LL,ICOL)
            A(LL,ICOL)=0.
            DO 18 L=1,N
              A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
18          CONTINUE
            DO 19 L=1,M
              B(LL,L)=B(LL,L)-B(ICOL,L)*DUM
19          CONTINUE
          ENDIF
21      CONTINUE
22    CONTINUE
      DO 24 L=N,1,-1
        IF(INDXR(L).NE.INDXC(L))THEN
          DO 23 K=1,N
            DUM=A(K,INDXR(L))
            A(K,INDXR(L))=A(K,INDXC(L))
            A(K,INDXC(L))=DUM
23        CONTINUE
        ENDIF
24    CONTINUE
      RETURN
      END
