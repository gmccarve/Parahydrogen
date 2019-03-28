        implicit real*8 (a-h, o-z)

        INTEGER count, maxj, pts, n, nn, time, TI
        PARAMETER (maxj=1)
        PARAMETER (n=INT((4.0/3.0)*maxj**3+4*maxj**2+(11.0/3.0)*maxj-1))
        PARAMETER (nn=n, pts=1)
        dimension a(nn,nn), evec(nn,nn), eval(nn), fv1(nn), fv2(nn)
        DOUBLE PRECISION coef(nn,nn), coeft(nn,nn), hamil0(nn,nn)
        DOUBLE PRECISION hamil1(nn,nn), fld
        DOUBLE PRECISION hamil2(nn,nn), hamil3(nn,nn)
        DOUBLE PRECISION fieldmat(n+1, pts)

        TI = time()

        a = 0.d0
        hamil1 = 0.d0
        hamil0 = 0.d0
        hamil2 = 0.d0
        coef   = 0.d0
        coeft  = 0.d0
        fieldmat = 0.d0

        call transpose(maxj, n, nn, coef, coeft, hamil0)

        do count = 1, pts

        fld = real(count)
c       fld = 1

        call fullfield(hamil1, fld, maxj, n)

        print *, hamil1

        a = hamil0 + hamil1

        call mult(a, coef, coeft, n, hamil2)

        l = 1

c rss --> ordered
c rs --> non ordered

c        call rss(nn, n, hamil2, eval, l, evec, fv1, fv2, ierr)
        call rss(nn, n, a, eval, l, evec, fv1, fv2, ierr)
c        call rss(nn, n, hamil2, eval, l, evec, fv1, fv2, ierr)

c       write (6, *) 'ierr = ', ierr

c       do i=1, n
c               write(20, *) eval(i)
c               write (6, *) i, eval(i)
c               write (6, *) ''
c               z=0.0d0
c               do j=1, n
c                       write (6, *) j, evec(j, i)
c                       z=z+evec(j, i)**2
c               end do
c               write (6, *) ''
c               write (6, *) 'norm = ', sqrt(z)
c               write (6, *)
c       end do

        fieldmat(1,count) = fld
        fieldmat(2:n+1, count) = eval

        enddo

        do i = 1, n+1
        write(23,*) fieldmat(i, 1:pts)
        enddo

c       print *, time()-TI

        end

ccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine mult(a, coef, coeft, n, mult2)

        DOUBLE PRECISION a(n,n), coef(n,n), coeft(n,n)
        DOUBLE PRECISION mult1(n,n), mult2(n,n)
        INTEGER i, j, k

        mult1 = 0.d0
        mult2 = 0.d0

        do i = 1, n
        do j = 1, n
        do k = 1, n

        mult1(i,j) = mult1(i,j) + coef(i,k)*a(k,j)

        enddo
        enddo
        enddo

        do i = 1, n
        do j = 1, n
        do k = 1, n

        mult2(i,j) = mult2(i,j) + mult1(i,k)*coeft(k,j)

        enddo
        enddo
        enddo

        end


cccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine transpose(maxj, n, nn, coef, coeft, hamil0)

        DOUBLE PRECISION coef(n,n), coeft(n,n), eval(nn)
        DOUBLE PRECISION evec(nn,nn), fv1(nn), fv2(nn), hamil0(n,n)
        INTEGER i, j


        hamil0 = 0.d0

        call nofield(hamil0, maxj, n)

        call rs(nn, n, hamil0, eval, 1, evec, fv1, fv2, ierr)

        coef = evec

        do i = 1, n
        do j = 1, n

        coeft(j,i) = coef(i,j)

        enddo
        enddo

        end

ccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine fullfield(mat, field, jmax, sze)

        INTEGER nmin, nmax, n, i, j, k, m, aa, bb, gg, xx
        DOUBLE PRECISION lilw, bigW, lild, weight, sum1, norm
        INTEGER  kmax, mmax, kmin, mmin
        INTEGER jm1, jm2, jk1, jk2, jmn, jkn, nmk
        DOUBLE PRECISION prod1, prod2, prod3, prod4
        DOUBLE PRECISION prod5, prod6, prod7, prod8
        INTEGER F1, F2, F3, F4, F5, F6, F7, F8
        INTEGER amax, bmax, gmax, jmax
        DOUBLE PRECISION wa, wb, wg, alpha, beta, gama, costheta
        COMPLEX II, BigD, JKM1, JKM2, EE, JKM

        INTEGER jj, kk, mm, nn, nnmin, nnmax
        DOUBLE PRECISION lilww, bigww, lildd, weightt,sum2,normm
        INTEGER jm12, jm22, jk12, jk22, jmn2, jkn2, nmk2
        DOUBLE PRECISION prod12, prod22, prod32, prod42
        DOUBLE PRECISION prod52, prod62, prod72, prod82
        INTEGER F12, F22, F32, F42, F52, F62, F72, F82, test
        COMPLEX BigDD
        INTEGER count1, count2, fle, sze
        DOUBLE PRECISION mat(sze,sze), field, dipole, stark, thresh
        DOUBLE PRECISION coef(sze,sze), coeft(sze,sze)
        DOUBLE PRECISION mult1(sze,sze), mult2(sze,sze), rot(sze,sze)
        DOUBLE PRECISION mat3(sze,sze)
        DOUBLE PRECISION mult3(sze,sze)


        fle   = test + 24
        mat   = 0.0d0
        mult2 = 0.0d0
        mult1 = 0.0d0
        mult3 = 0.0d0

        field  = (field-1)/100.d0
c       field  = 1.d0
        dipole = 1.855d0
        thresh = 0.000001d0

c        write(6,*) field

        count1 = 1
        count2 = 1

        pi = 4.d0 * ATAN(1.d0)
        II = (0,1)

        amax = 10
        bmax = 4
        gmax = 10

        wa = 2*pi/amax
        wg = 2*pi/gmax

        do 20 jj = 0, jmax

        Normm = sqrt((2*jj+1)/(8*pi**2))

        do 20 mm = -jj, jj
c       mm = 0
        do 20 kk = -jj, jj

        nnmin = max(0, (kk-mm))
        nnmax = min((jj-mm), (jj+kk))

        do 20 j = 0, jmax

        Norm = sqrt((2*j+1)/(8*pi**2))

        do 20 m = -j, j
c       m = 0
        do 20 k = -j, j

        sum1 = 0.d0

        nmin = max(0, (k-m))
        nmax = min((j-m), (j+k))

        do 30 aa = 1, amax
        do 30 bb = 1, bmax
        do 30 gg = 1, gmax

        EE = 0.d0

        lilD = 0.d0
        lildd = 0.d0

        alpha = (2*pi*(aa-1))/amax
        gama = (2*pi*(gg-1))/gmax

        if (bb==1) then
                beta = acos(-0.861136311594053)
                wb = 0.347854845137454
        elseif (bb==2) then
                beta = acos(-0.339981043584856)
                wb = 0.652145154862546
        elseif (bb==4) then
                beta = acos(0.861136311594053)
                wb = 0.347854845137454
        elseif (bb==3) then
                beta = acos(0.339981043584856)
                wb = 0.652145154862546
        endif

        do 10 n = nmin, nmax

        prod1 = 1.d0
        prod2 = 1.d0
        prod3 = 1.d0
        prod4 = 1.d0
        prod5 = 1.d0
        prod6 = 1.d0
        prod7 = 1.d0
        prod8 = 1.d0


        jm1 = j+m
        jm2 = j-m
        jk1 = j+k
        jk2 = j-k
        jmn = j-m-n
        jkn = j+k-n
        nmk = n+m-k

        do F1 = 1, jm1
                prod1 = prod1 * real(F1)
        enddo

        do F2 = 1, jm2
                prod2 = prod2 * real(F2)
        enddo

        do F3 = 1,jk1
                prod3 = prod3 * real(F3)
        enddo

        do F4 = 1,jk2
                prod4 = prod4 * real(F4)
        enddo

        do F5 = 1, jmn
                prod5 = prod5 * real(F5)
        enddo

        do F6 = 1, jkn
                prod6 = prod6 * real(F6)
        enddo

        do F7 = 1, nmk
                prod7 = prod7 * real(F7)
        enddo

        do F8 = 1, n
                prod8 = prod8 * real(F8)
        enddo

        lilw = sqrt(prod1 * prod2 * prod3 * prod4) / (prod5 * prod6 *
     + prod7 * prod8)

        bigW = lilw * ((cos(beta/2))**(2*j+k-m-2*n)) *
     + ((-sin(beta/2))**(m-k+2*n))

        lild = lild + ((-1)**n) * bigW

10      continue

        do 40 nn = nnmin, nnmax

        prod12 = 1.d0
        prod22 = 1.d0
        prod32 = 1.d0
        prod42 = 1.d0
        prod52 = 1.d0
        prod62 = 1.d0
        prod72 = 1.d0
        prod82 = 1.d0

        jm12 = jj+mm
        jm22 = jj-mm
        jk12 = jj+kk
        jk22 = jj-kk
        jmn2 = jj-mm-nn
        jkn2 = jj+kk-nn
        nmk2 = nn+mm-kk

        do F12 = 1, jm12
                prod12 = prod12 * real(F12)
        enddo

        do F22 = 1, jm22
                prod22 = prod22 * real(F22)
        enddo

        do F32 = 1,jk12
                prod32 = prod32 * real(F32)
        enddo

        do F42 = 1,jk22
                prod42 = prod42 * real(F42)
        enddo

        do F52 = 1, jmn2
                prod52 = prod52 * real(F52)
        enddo

        do F62 = 1, jkn2
                prod62 = prod62 * real(F62)
        enddo

        do F72 = 1, nmk2
                prod72 = prod72 * real(F72)
        enddo
        do F82 = 1, nn
                prod82 = prod82 * real(F82)
        enddo

        lilww = sqrt(prod12 * prod22 * prod32 * prod42) /
     + (prod52 * prod62 * prod72 * prod82)

        bigWW = lilww * ((cos(beta/2))**(2*jj+kk-mm-2*nn)) *
     + ((-sin(beta/2))**(mm-kk+2*nn))

        lildd = lildd + ((-1)**nn) * bigWW

40      continue


        BigD = exp(-II*(m*alpha + k*gama))*lild

        BigDD = exp(-II*(mm*alpha + kk*gama))*lildd

        JKM1 = CONJG(BigDD)*Normm
        JKM2 = BigD*Norm
        stark = sin(beta) * cos(gama)
        weight = wa*wb*wg

        E = weight*stark*JKM1*JKM2*field*dipole

        sum1 = sum1 + E

30      continue

        if (ABS(sum1) .LT. thresh) then
                sum1 = 0.d0
        endif

        mat(count1, count2) = sum1

        count1 = count1 + 1

        if (count1 == sze+1 ) then
                count1 = 1
                count2 = count2 + 1
        endif

20      continue

        end

ccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccc


        subroutine nofield(rot, jmax, sze)

        INTEGER j, k, m, kk, jmax, sze
        DOUBLE PRECISION rot(sze,sze)

c NIST Water
        A = 27.8761
        B = 14.5074
        C = 9.2877
c Bulthuis Water
c       A = 27.33
c       B = 14.575
c       C = 9.499
c Indole
c       A = 0.18919
c       B = 0.02503
c       C = 0.02210
c D2O
c       A = 15.4204136
c       B = 7.2710285
c       C = 4.8467977

        IULC = 0

        do i=1, sze
        do j=1, sze
                rot(i, j)=0.0d0
                if (i.eq.j) rot(i, j)=dble(i)
        end do
        end do

        do 20 J =   0 , jmax
        do 20 M =  -j , j
        do 30 K =  -j , j
        do 30 KK = -j , j

        val = 0.d0

        if (K == KK) then

        val = (0.5) * (B+C) * (J*(J+1) - K**2) + A*K**2

        elseif (KK == K + 2) then

        val = (0.25) * (B-C) * sqrt(real(J * (J+1) - K * (K+1))) *
     + sqrt(real(J * (J+1) - (K+1) * (K+2)))

        elseif (KK == K - 2) then

        val = (0.25) * (B-C) * sqrt(real(J * (J+1) - K * (K-1))) *
     + sqrt(real(J * (J+1) - (K-1) * (K-2)))

        endif

        rot(IULC+KK+J+1, IULC+K+J+1) = val

        val = 0.d0

30      continue

        IULC = IULC + 2*j+1

20      continue

        end
ccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    Organized eigenvalues/vectors
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine rss(nm,n,a,w,matz,z,fv1,fv2,ierr)
c
      integer n,nm,ierr,matz
      double precision a(nm,n),w(n),z(nm,n),fv1(n),fv2(n)
c
c     this subroutine calls the recommended sequence of
c     subroutines from the eigensystem subroutine package (eispack)
c     to find the eigenvalues and eigenvectors (if desired)
c     of a real symmetric matrix.
c
c     on input
c
c        nm  must be set to the row dimension of the two-dimensional
c        array parameters as declared in the calling program
c        dimension statement.
c
c        n  is the order of the matrix  a.
c
c        a  contains the real symmetric matrix.
c
c        matz  is an integer variable set equal to zero if
c        only eigenvalues are desired.  otherwise it is set to
c        any non-zero integer for both eigenvalues and eigenvectors.
c
c     on output
c
c        w  contains the eigenvalues in ascending order.
c
c        z  contains the eigenvectors if matz is not zero.
c
c        ierr  is an integer output variable set equal to an error
c           completion code described in the documentation for tqlrat
c           and tql2.  the normal completion code is zero.
c
c        fv1  and  fv2  are temporary storage arrays.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      if (n .le. nm) go to 10
      write (6, *) 'IERR!!'
      ierr = 10 * n
      go to 50
c
   10 if (matz .ne. 0) go to 20
c     .......... find eigenvalues only ..........
      call  tred1(nm,n,a,w,fv1,fv2)
*  tqlrat encounters catastrophic underflow on the Vax
*     call  tqlrat(n,w,fv2,ierr)
      call  tql1(n,w,fv1,ierr)
      go to 50
c     .......... find both eigenvalues and eigenvectors ..........
   20 call  tred2(nm,n,a,w,fv1,z)
      call  tql23(nm,n,w,fv1,z,ierr)
   50 return
      end

      subroutine tql23(nm,n,d,e,z,ierr)
c
      integer i,j,k,l,m,n,ii,l1,l2,nm,mml,ierr
      double precision d(n),e(n),z(nm,n)
      double precision c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2,pythag
c
c     this subroutine is a translation of the algol procedure tql23,
c     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
c     wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
c
c     this subroutine finds the eigenvalues and eigenvectors
c     of a symmetric tridiagonal matrix by the ql method.
c     the eigenvectors of a full symmetric matrix can also
c     be found if  tred2  has been used to reduce this
c     full matrix to tridiagonal form.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        d contains the diagonal elements of the input matrix.
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary.
c
c        z contains the transformation matrix produced in the
c          reduction by  tred2, if performed.  if the eigenvectors
c          of the tridiagonal matrix are desired, z must contain
c          the identity matrix.
c
c      on output
c
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct but
c          unordered for indices 1,2,...,ierr-1.
c
c        e has been destroyed.
c
c        z contains orthonormal eigenvectors of the symmetric
c          tridiagonal (or full) matrix.  if an error exit is made,
c          z contains the eigenvectors associated with the stored
c          eigenvalues.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c     calls pythag for  sqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 0
      if (n .eq. 1) go to 1001
c
      do 100 i = 2, n
  100 e(i-1) = e(i)
c
      f = 0.0
      tst1 = 0.0
      e(n) = 0.0
c
      do 240 l = 1, n
         j = 0
         h = abs(d(l)) + abs(e(l))
         if (tst1 .lt. h) tst1 = h
c     .......... look for small sub-diagonal element ..........
         do 110 m = l, n
            tst2 = tst1 + abs(e(m))
            if (tst2 .eq. tst1) go to 120
c     .......... e(n) is always zero, so there is no exit
c                through the bottom of the loop ..........
  110    continue
c
  120    if (m .eq. l) go to 220
  130    if (j .eq. 30) go to 1000
         j = j + 1
c     .......... form shift ..........
         l1 = l + 1
         l2 = l1 + 1
         g = d(l)
         p = (d(l1) - g) / (2.0 * e(l))
         r = pythag(p,1.0d0)
         d(l) = e(l) / (p + sign(r,p))
         d(l1) = e(l) * (p + sign(r,p))
         dl1 = d(l1)
         h = g - d(l)
         if (l2 .gt. n) go to 145
c
         do 140 i = l2, n
  140    d(i) = d(i) - h
c
  145    f = f + h
c     .......... ql transformation ..........
         p = d(m)
         c = 1.0
         c2 = c
         el1 = e(l1)
         s = 0.0
         mml = m - l
c     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c * e(i)
            h = c * p
            r = pythag(p,e(i))
            e(i+1) = s * r
            s = e(i) / r
            c = p / r
            p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
c     .......... form vector ..........
            do 180 k = 1, n
               h = z(k,i+1)
               z(k,i+1) = s * z(k,i) + c * h
               z(k,i) = c * z(k,i) - s * h
  180       continue
c
  200    continue
c
         p = -s * s2 * c3 * el1 * e(l) / dl1
         e(l) = s * p
         d(l) = c * p
         tst2 = tst1 + abs(e(l))
         if (tst2 .gt. tst1) go to 130
  220    d(l) = d(l) + f
  240 continue
c     .......... order eigenvalues and eigenvectors ..........
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
c
         do 260 j = ii, n
            if (d(j) .ge. p) go to 260
            k = j
            p = d(j)
  260    continue
c
         if (k .eq. i) go to 300
         d(k) = d(i)
         d(i) = p
c
         do 280 j = 1, n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
  280    continue
c
  300 continue
c
      go to 1001
c     .......... set error -- no convergence to an
c                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c --- everything from this line down is the actual RS code
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine rs(nm,n,a,w,matz,z,fv1,fv2,ierr)
c
      integer n,nm,ierr,matz
      double precision a(nm,n),w(n),z(nm,n),fv1(n),fv2(n)
c
c     this subroutine calls the recommended sequence of
c     subroutines from the eigensystem subroutine package (eispack)
c     to find the eigenvalues and eigenvectors (if desired)
c     of a real symmetric matrix.
c
c     on input
c
c        nm  must be set to the row dimension of the two-dimensional
c        array parameters as declared in the calling program
c        dimension statement.
c
c        n  is the order of the matrix  a.
c
c        a  contains the real symmetric matrix.
c
c        matz  is an integer variable set equal to zero if
c        only eigenvalues are desired.  otherwise it is set to
c        any non-zero integer for both eigenvalues and eigenvectors.
c
c     on output
c
c        w  contains the eigenvalues in ascending order.
c
c        z  contains the eigenvectors if matz is not zero.
c
c        ierr  is an integer output variable set equal to an error
c           completion code described in the documentation for tqlrat
c           and tql2.  the normal completion code is zero.
c
c        fv1  and  fv2  are temporary storage arrays.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      if (n .le. nm) go to 10
      write (6, *) 'IERR!!'
      ierr = 10 * n
      go to 50
c
   10 if (matz .ne. 0) go to 20
c     .......... find eigenvalues only ..........
      call  tred1(nm,n,a,w,fv1,fv2)
*  tqlrat encounters catastrophic underflow on the Vax
*     call  tqlrat(n,w,fv2,ierr)
      call  tql1(n,w,fv1,ierr)
      go to 50
c     .......... find both eigenvalues and eigenvectors ..........
   20 call  tred2(nm,n,a,w,fv1,z)
      call  tql2(nm,n,w,fv1,z,ierr)
   50 return
      end

      subroutine tql2(nm,n,d,e,z,ierr)
c
      integer i,j,k,l,m,n,ii,l1,l2,nm,mml,ierr
      double precision d(n),e(n),z(nm,n)
      double precision c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2,pythag
c
c     this subroutine is a translation of the algol procedure tql2,
c     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
c     wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
c
c     this subroutine finds the eigenvalues and eigenvectors
c     of a symmetric tridiagonal matrix by the ql method.
c     the eigenvectors of a full symmetric matrix can also
c     be found if  tred2  has been used to reduce this
c     full matrix to tridiagonal form.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        d contains the diagonal elements of the input matrix.
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary.
c
c        z contains the transformation matrix produced in the
c          reduction by  tred2, if performed.  if the eigenvectors
c          of the tridiagonal matrix are desired, z must contain
c          the identity matrix.
c
c      on output
c
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct but
c          unordered for indices 1,2,...,ierr-1.
c
c        e has been destroyed.
c
c        z contains orthonormal eigenvectors of the symmetric
c          tridiagonal (or full) matrix.  if an error exit is made,
c          z contains the eigenvectors associated with the stored
c          eigenvalues.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c     calls pythag for  sqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 0
      if (n .eq. 1) go to 1001
c
      do 100 i = 2, n
  100 e(i-1) = e(i)
c
      f = 0.0
      tst1 = 0.0
      e(n) = 0.0
c
      do 240 l = 1, n
         j = 0
         h = abs(d(l)) + abs(e(l))
         if (tst1 .lt. h) tst1 = h
c     .......... look for small sub-diagonal element ..........
         do 110 m = l, n
            tst2 = tst1 + abs(e(m))
            if (tst2 .eq. tst1) go to 120
c     .......... e(n) is always zero, so there is no exit
c                through the bottom of the loop ..........
  110    continue
c
  120    if (m .eq. l) go to 220
  130    if (j .eq. 30) go to 1000
         j = j + 1
c     .......... form shift ..........
         l1 = l + 1
         l2 = l1 + 1
         g = d(l)
         p = (d(l1) - g) / (2.0 * e(l))
         r = pythag(p,1.0d0)
         d(l) = e(l) / (p + sign(r,p))
         d(l1) = e(l) * (p + sign(r,p))
         dl1 = d(l1)
         h = g - d(l)
         if (l2 .gt. n) go to 145
c
         do 140 i = l2, n
  140    d(i) = d(i) - h
c
  145    f = f + h
c     .......... ql transformation ..........
         p = d(m)
         c = 1.0
         c2 = c
         el1 = e(l1)
         s = 0.0
         mml = m - l
c     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c * e(i)
            h = c * p
            r = pythag(p,e(i))
            e(i+1) = s * r
            s = e(i) / r
            c = p / r
            p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
c     .......... form vector ..........
            do 180 k = 1, n
               h = z(k,i+1)
               z(k,i+1) = s * z(k,i) + c * h
               z(k,i) = c * z(k,i) - s * h
  180       continue
c
  200    continue
c
         p = -s * s2 * c3 * el1 * e(l) / dl1
         e(l) = s * p
         d(l) = c * p
         tst2 = tst1 + abs(e(l))
         if (tst2 .gt. tst1) go to 130
  220    d(l) = d(l) + f
  240 continue

      go to 1001
c     .......... set error -- no convergence to an
c                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end

      subroutine tred2(nm,n,a,d,e,z)
c
      integer i,j,k,l,n,ii,nm,jp1
      double precision a(nm,n),d(n),e(n),z(nm,n)
      double precision f,g,h,hh,scale
c
c     this subroutine is a translation of the algol procedure tred2,
c     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine reduces a real symmetric matrix to a
c     symmetric tridiagonal matrix using and accumulating
c     orthogonal similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        a contains the real symmetric input matrix.  only the
c          lower triangle of the matrix need be supplied.
c
c     on output
c
c        d contains the diagonal elements of the tridiagonal matrix.
c
c        e contains the subdiagonal elements of the tridiagonal
c          matrix in its last n-1 positions.  e(1) is set to zero.
c
c        z contains the orthogonal transformation matrix
c          produced in the reduction.
c
c        a and z may coincide.  if distinct, a is unaltered.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      do 100 i = 1, n
c
         do 80 j = i, n
   80    z(j,i) = a(j,i)
c
         d(i) = a(n,i)
  100 continue
c
      if (n .eq. 1) go to 510
c     .......... for i=n step -1 until 2 do -- ..........
      do 300 ii = 2, n
         i = n + 2 - ii
         l = i - 1
         h = 0.0
         scale = 0.0
         if (l .lt. 2) go to 130
c     .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, l
  120    scale = scale + abs(d(k))
c
         if (scale .ne. 0.0) go to 140
  130    e(i) = d(l)
c
         do 135 j = 1, l
            d(j) = z(l,j)
            z(i,j) = 0.0
            z(j,i) = 0.0
  135    continue
c
         go to 290
c
  140    do 150 k = 1, l
            d(k) = d(k) / scale
            h = h + d(k) * d(k)
  150    continue
c
         f = d(l)
         g = -sign(sqrt(h),f)
         e(i) = scale * g
         h = h - f * g
         d(l) = f - g
c     .......... form a*u ..........
         do 170 j = 1, l
  170    e(j) = 0.0
c
         do 240 j = 1, l
            f = d(j)
            z(j,i) = f
            g = e(j) + z(j,j) * f
            jp1 = j + 1
            if (l .lt. jp1) go to 220
c
            do 200 k = jp1, l
               g = g + z(k,j) * d(k)
               e(k) = e(k) + z(k,j) * f
  200       continue
c
  220       e(j) = g
  240    continue
c     .......... form p ..........
         f = 0.0
c
         do 245 j = 1, l
            e(j) = e(j) / h
            f = f + e(j) * d(j)
  245    continue
c
         hh = f / (h + h)
c     .......... form q ..........
         do 250 j = 1, l
  250    e(j) = e(j) - hh * d(j)
c     .......... form reduced a ..........
         do 280 j = 1, l
            f = d(j)
            g = e(j)
c
            do 260 k = j, l
  260       z(k,j) = z(k,j) - f * e(k) - g * d(k)
c
            d(j) = z(l,j)
            z(i,j) = 0.0
  280    continue
c
  290    d(i) = h
  300 continue
c     .......... accumulation of transformation matrices ..........
      do 500 i = 2, n
         l = i - 1
         z(n,l) = z(l,l)
         z(l,l) = 1.0
         h = d(i)
         if (h .eq. 0.0) go to 380
c
         do 330 k = 1, l
  330    d(k) = z(k,i) / h
c
         do 360 j = 1, l
            g = 0.0
c
            do 340 k = 1, l
  340       g = g + z(k,i) * z(k,j)
c
            do 360 k = 1, l
               z(k,j) = z(k,j) - g * d(k)
  360    continue
c
  380    do 400 k = 1, l
  400    z(k,i) = 0.0
c
  500 continue
c
  510 do 520 i = 1, n
         d(i) = z(n,i)
         z(n,i) = 0.0
  520 continue
c
      z(n,n) = 1.0
      e(1) = 0.0
      return
      end

      subroutine tql1(n,d,e,ierr)
c
      integer i,j,l,m,n,ii,l1,l2,mml,ierr
      double precision d(n),e(n)
      double precision c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2,pythag
c
c     this subroutine is a translation of the algol procedure tql1,
c     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
c     wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
c
c     this subroutine finds the eigenvalues of a symmetric
c     tridiagonal matrix by the ql method.
c
c     on input
c
c        n is the order of the matrix.
c
c        d contains the diagonal elements of the input matrix.
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary.
c
c      on output
c
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct and
c          ordered for indices 1,2,...ierr-1, but may not be
c          the smallest eigenvalues.
c
c        e has been destroyed.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c     calls pythag for  sqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 0
      if (n .eq. 1) go to 1001
c
      do 100 i = 2, n
  100 e(i-1) = e(i)
c
      f = 0.0
      tst1 = 0.0
      e(n) = 0.0
c
      do 290 l = 1, n
         j = 0
         h = abs(d(l)) + abs(e(l))
         if (tst1 .lt. h) tst1 = h
c     .......... look for small sub-diagonal element ..........
         do 110 m = l, n
            tst2 = tst1 + abs(e(m))
            if (tst2 .eq. tst1) go to 120
c     .......... e(n) is always zero, so there is no exit
c                through the bottom of the loop ..........
  110    continue
c
  120    if (m .eq. l) go to 210
  130    if (j .eq. 30) go to 1000
         j = j + 1
c     .......... form shift ..........
         l1 = l + 1
         l2 = l1 + 1
         g = d(l)
         p = (d(l1) - g) / (2.0 * e(l))
         r = pythag(p,1.0d0)
         d(l) = e(l) / (p + sign(r,p))
         d(l1) = e(l) * (p + sign(r,p))
         dl1 = d(l1)
         h = g - d(l)
         if (l2 .gt. n) go to 145
c
         do 140 i = l2, n
  140    d(i) = d(i) - h
c
  145    f = f + h
c     .......... ql transformation ..........
         p = d(m)
         c = 1.0
         c2 = c
         el1 = e(l1)
         s = 0.0
         mml = m - l
c     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c * e(i)
            h = c * p
            r = pythag(p,e(i))
            e(i+1) = s * r
            s = e(i) / r
            c = p / r
            p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
  200    continue
c
         p = -s * s2 * c3 * el1 * e(l) / dl1
         e(l) = s * p
         d(l) = c * p
         tst2 = tst1 + abs(e(l))
         if (tst2 .gt. tst1) go to 130
  210    p = d(l) + f
c     .......... order eigenvalues ..........
         if (l .eq. 1) go to 250
c     .......... for i=l step -1 until 2 do -- ..........
         do 230 ii = 2, l
            i = l + 2 - ii
            if (p .ge. d(i-1)) go to 270
            d(i) = d(i-1)
  230    continue
c
  250    i = 1
  270    d(i) = p
  290 continue
c
      go to 1001
c     .......... set error -- no convergence to an
c                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end

      double precision function pythag(a,b)
      double precision a,b
c
c     finds sqrt(a**2+b**2) without overflow or destructive underflow
c
      double precision p,r,s,t,u
      p = dmax1(abs(a),abs(b))
      if (p .eq. 0.0) go to 20
      r = (dmin1(abs(a),abs(b))/p)**2
   10 continue
         t = 4.0 + r
         if (t .eq. 4.0) go to 20
         s = r/t
         u = 1.0 + 2.0*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 pythag = p
      return
      end

      subroutine tred1(nm,n,a,d,e,e2)
c
      integer i,j,k,l,n,ii,nm,jp1
      double precision a(nm,n),d(n),e(n),e2(n)
      double precision f,g,h,scale
c
c     this subroutine is a translation of the algol procedure tred1,
c     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine reduces a real symmetric matrix
c     to a symmetric tridiagonal matrix using
c     orthogonal similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        a contains the real symmetric input matrix.  only the
c          lower triangle of the matrix need be supplied.
c
c     on output
c
c        a contains information about the orthogonal trans-
c          formations used in the reduction in its strict lower
c          triangle.  the full upper triangle of a is unaltered.
c
c        d contains the diagonal elements of the tridiagonal matrix.
c
c        e contains the subdiagonal elements of the tridiagonal
c          matrix in its last n-1 positions.  e(1) is set to zero.
c
c        e2 contains the squares of the corresponding elements of e.
c          e2 may coincide with e if the squares are not needed.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      do 100 i = 1, n
         d(i) = a(n,i)
         a(n,i) = a(i,i)
  100 continue
c     .......... for i=n step -1 until 1 do -- ..........
      do 300 ii = 1, n
         i = n + 1 - ii
         l = i - 1
         h = 0.0
         scale = 0.0
         if (l .lt. 1) go to 130
c     .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, l
  120    scale = scale + abs(d(k))
c
         if (scale .ne. 0.0) go to 140
c
         do 125 j = 1, l
            d(j) = a(l,j)
            a(l,j) = a(i,j)
            a(i,j) = 0.0
  125    continue
c
  130    e(i) = 0.0
         e2(i) = 0.0
         go to 300
c
  140    do 150 k = 1, l
            d(k) = d(k) / scale
            h = h + d(k) * d(k)
  150    continue
c
         e2(i) = scale * scale * h
         f = d(l)
         g = -sign(sqrt(h),f)
         e(i) = scale * g
         h = h - f * g
         d(l) = f - g
         if (l .eq. 1) go to 285
c     .......... form a*u ..........
         do 170 j = 1, l
  170    e(j) = 0.0
c
         do 240 j = 1, l
            f = d(j)
            g = e(j) + a(j,j) * f
            jp1 = j + 1
            if (l .lt. jp1) go to 220
c
            do 200 k = jp1, l
               g = g + a(k,j) * d(k)
               e(k) = e(k) + a(k,j) * f
  200       continue
c
  220       e(j) = g
  240    continue
c     .......... form p ..........
         f = 0.0
c
         do 245 j = 1, l
            e(j) = e(j) / h
            f = f + e(j) * d(j)
  245    continue
c
         h = f / (h + h)
c     .......... form q ..........
         do 250 j = 1, l
  250    e(j) = e(j) - h * d(j)
c     .......... form reduced a ..........
         do 280 j = 1, l
            f = d(j)
            g = e(j)
c
            do 260 k = j, l
  260       a(k,j) = a(k,j) - f * e(k) - g * d(k)
c
  280    continue
c
  285    do 290 j = 1, l
            f = d(j)
            d(j) = a(l,j)
            a(l,j) = a(i,j)
            a(i,j) = f * scale
  290    continue
c
  300 continue
c
      return
      end

      double precision function epslon (x)
      double precision x
c
c     estimate unit roundoff in quantities of size x.
c
      double precision a,b,c,eps
c
c     this program should function properly on all systems
c     satisfying the following two assumptions,
c        1.  the base used in representing floating point
c            numbers is not a power of three.
c        2.  the quantity  a  in statement 10 is represented to
c            the accuracy used in floating point variables
c            that are stored in memory.
c     the statement number 10 and the go to 10 are intended to
c     force optimizing compilers to generate code satisfying
c     assumption 2.
c     under these assumptions, it should be true that,
c            a  is not exactly equal to four-thirds,
c            b  has a zero for its last bit or digit,
c            c  is not exactly equal to one,
c            eps  measures the separation of 1.0 from
c                 the next larger floating point number.
c     the developers of eispack would appreciate being informed
c     about any systems where these assumptions do not hold.
c
c     this version dated 4/6/83.
c
      a = 4.0/3.0
   10 b = a - 1.0
      c = b + b + b
      eps = abs(c-1.0)
      if (eps .eq. 0.0) go to 10
      epslon = eps*abs(x)
      return
      end
      subroutine tqlrat(n,d,e2,ierr)
c
      integer i,j,l,m,n,ii,l1,mml,ierr
      double precision d(n),e2(n)
      double precision b,c,f,g,h,p,r,s,t,epslon,pythag
c
c     this subroutine is a translation of the algol procedure tqlrat,
c     algorithm 464, comm. acm 16, 689(1973) by reinsch.
c
c     this subroutine finds the eigenvalues of a symmetric
c     tridiagonal matrix by the rational ql method.
c
c     on input
c
c        n is the order of the matrix.
c
c        d contains the diagonal elements of the input matrix.
c
c        e2 contains the squares of the subdiagonal elements of the
c          input matrix in its last n-1 positions.  e2(1) is arbitrary.
c
c      on output
c
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct and
c          ordered for indices 1,2,...ierr-1, but may not be
c          the smallest eigenvalues.
c
c        e2 has been destroyed.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c     calls pythag for  sqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 0
      if (n .eq. 1) go to 1001
c
      do 100 i = 2, n
  100 e2(i-1) = e2(i)
c
      f = 0.0
      t = 0.0
      e2(n) = 0.0
c
      do 290 l = 1, n
         j = 0
         h = abs(d(l)) + sqrt(e2(l))
         if (t .gt. h) go to 105
         t = h
         b = epslon(t)
         c = b * b
c     .......... look for small squared sub-diagonal element ..........
  105    do 110 m = l, n
            if (e2(m) .le. c) go to 120
c     .......... e2(n) is always zero, so there is no exit
c                through the bottom of the loop ..........
  110    continue
c
  120    if (m .eq. l) go to 210
  130    if (j .eq. 30) go to 1000
         j = j + 1
c     .......... form shift ..........
         l1 = l + 1
         s = sqrt(e2(l))
         g = d(l)
         p = (d(l1) - g) / (2.0 * s)
         r = pythag(p,1.0d0)
         d(l) = s / (p + sign(r,p))
         h = g - d(l)
c
         do 140 i = l1, n
  140    d(i) = d(i) - h
c
         f = f + h
c     .......... rational ql transformation ..........
         g = d(m)
         if (g .eq. 0.0) g = b
         h = g
         s = 0.0
         mml = m - l
c     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            i = m - ii
            p = g * h
            r = p + e2(i)
            e2(i+1) = s * r
            s = e2(i) / r
            d(i+1) = h + s * (h + d(i))
            g = d(i) - e2(i) / g
            if (g .eq. 0.0) g = b
            h = g * p / r
  200    continue
c
         e2(l) = s * g
         d(l) = h
c     .......... guard against underflow in convergence test ..........
         if (h .eq. 0.0) go to 210
         if (abs(e2(l)) .le. abs(c/h)) go to 210
         e2(l) = h * e2(l)
         if (e2(l) .ne. 0.0) go to 130
  210    p = d(l) + f
c     .......... order eigenvalues ..........
         if (l .eq. 1) go to 250
c     .......... for i=l step -1 until 2 do -- ..........
         do 230 ii = 2, l
            i = l + 2 - ii
            if (p .ge. d(i-1)) go to 270
            d(i) = d(i-1)
  230    continue
c
  250    i = 1
  270    d(i) = p
  290 continue
c
      go to 1001
c     .......... set error -- no convergence to an
c                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end
      subroutine rebak(nm,n,b,dl,m,z)
c
      integer i,j,k,m,n,i1,ii,nm
      double precision b(nm,n),dl(n),z(nm,m)
      double precision x
c
c     this subroutine is a translation of the algol procedure rebaka,
c     num. math. 11, 99-110(1968) by martin and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 303-314(1971).
c
c     this subroutine forms the eigenvectors of a generalized
c     symmetric eigensystem by back transforming those of the
c     derived symmetric matrix determined by  reduc.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix system.
c
c        b contains information about the similarity transformation
c          (cholesky decomposition) used in the reduction by  reduc
c          in its strict lower triangle.
c
c        dl contains further information about the transformation.
c
c        m is the number of eigenvectors to be back transformed.
c
c        z contains the eigenvectors to be back transformed
c          in its first m columns.
c
c     on output
c
c        z contains the transformed eigenvectors
c          in its first m columns.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      if (m .eq. 0) go to 200
c
      do 100 j = 1, m
c     .......... for i=n step -1 until 1 do -- ..........
         do 100 ii = 1, n
            i = n + 1 - ii
            i1 = i + 1
            x = z(i,j)
            if (i .eq. n) go to 80
c
            do 60 k = i1, n
   60       x = x - b(k,i) * z(k,j)
c
   80       z(i,j) = x / dl(i)
  100 continue
c
  200 return
      end
      subroutine reduc(nm,n,a,b,dl,ierr)
c
      integer i,j,k,n,i1,j1,nm,nn,ierr
      double precision a(nm,n),b(nm,n),dl(n)
      double precision x,y
c
c     this subroutine is a translation of the algol procedure reduc1,
c     num. math. 11, 99-110(1968) by martin and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 303-314(1971).
c
c     this subroutine reduces the generalized symmetric eigenproblem
c     ax=(lambda)bx, where b is positive definite, to the standard
c     symmetric eigenproblem using the cholesky factorization of b.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrices a and b.  if the cholesky
c          factor l of b is already available, n should be prefixed
c          with a minus sign.
c
c        a and b contain the real symmetric input matrices.  only the
c          full upper triangles of the matrices need be supplied.  if
c          n is negative, the strict lower triangle of b contains,
c          instead, the strict lower triangle of its cholesky factor l.
c
c        dl contains, if n is negative, the diagonal elements of l.
c
c     on output
c
c        a contains in its full lower triangle the full lower triangle
c          of the symmetric matrix derived from the reduction to the
c          standard form.  the strict upper triangle of a is unaltered.
c
c        b contains in its strict lower triangle the strict lower
c          triangle of its cholesky factor l.  the full upper
c          triangle of b is unaltered.
c
c        dl contains the diagonal elements of l.
c
c        ierr is set to
c          zero       for normal return,
c          7*n+1      if b is not positive definite.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 0
      nn = iabs(n)
      if (n .lt. 0) go to 100
c     .......... form l in the arrays b and dl ..........
      do 80 i = 1, n
         i1 = i - 1
c
         do 80 j = i, n
            x = b(i,j)
            if (i .eq. 1) go to 40
c
            do 20 k = 1, i1
   20       x = x - b(i,k) * b(j,k)
c
   40       if (j .ne. i) go to 60
            if (x .le. 0.0) go to 1000
            y = sqrt(x)
            dl(i) = y
            go to 80
   60       b(j,i) = x / y
   80 continue
c     .......... form the transpose of the upper triangle of inv(l)*a
c                in the lower triangle of the array a ..........
  100 do 200 i = 1, nn
         i1 = i - 1
         y = dl(i)
c
         do 200 j = i, nn
            x = a(i,j)
            if (i .eq. 1) go to 180
c
            do 160 k = 1, i1
  160       x = x - b(i,k) * a(j,k)
c
  180       a(j,i) = x / y
  200 continue
c     .......... pre-multiply by inv(l) and overwrite ..........
      do 300 j = 1, nn
         j1 = j - 1
c
         do 300 i = j, nn
            x = a(i,j)
            if (i .eq. j) go to 240
            i1 = i - 1
c
            do 220 k = j, i1
  220       x = x - a(k,j) * b(i,k)
c
  240       if (j .eq. 1) go to 280
c
            do 260 k = 1, j1
  260       x = x - a(j,k) * b(i,k)
c
  280       a(i,j) = x / dl(i)
  300 continue
c
      go to 1001
c     .......... set error -- b is not positive definite ..........
 1000 ierr = 7 * n + 1
 1001 return
      end
      subroutine rsg(nm,n,a,b,w,matz,z,fv1,fv2,ierr)
c
      integer n,nm,ierr,matz
      double precision a(nm,n),b(nm,n),w(n),z(nm,n),fv1(n),fv2(n)
c
c     this subroutine calls the recommended sequence of
c     subroutines from the eigensystem subroutine package (eispack)
c     to find the eigenvalues and eigenvectors (if desired)
c     for the real symmetric generalized eigenproblem  ax = (lambda)bx.
c
c     on input
c
c        nm  must be set to the row dimension of the two-dimensional
c        array parameters as declared in the calling program
c        dimension statement.
c
c        n  is the order of the matrices  a  and  b.
c
c        a  contains a real symmetric matrix.
c
c        b  contains a positive definite real symmetric matrix.
c
c        matz  is an integer variable set equal to zero if
c        only eigenvalues are desired.  otherwise it is set to
c        any non-zero integer for both eigenvalues and eigenvectors.
c
c     on output
c
c        w  contains the eigenvalues in ascending order.
c
c        z  contains the eigenvectors if matz is not zero.
c
c        ierr  is an integer output variable set equal to an error
c           completion code described in the documentation for tqlrat
c           and tql2.  the normal completion code is zero.
c
c        fv1  and  fv2  are temporary storage arrays.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
c
   10 call  reduc(nm,n,a,b,fv2,ierr)
      if (ierr .ne. 0) go to 50
      if (matz .ne. 0) go to 20
c     .......... find eigenvalues only ..........
      call  tred1(nm,n,a,w,fv1,fv2)
      call  tqlrat(n,w,fv2,ierr)
      go to 50
c     .......... find both eigenvalues and eigenvectors ..........
   20 call  tred2(nm,n,a,w,fv1,z)
      call  tql2(nm,n,w,fv1,z,ierr)
      if (ierr .ne. 0) go to 50
      call  rebak(nm,n,b,fv2,n,z)
   50 return
      end
