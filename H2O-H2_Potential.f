        INTEGER CNT1, CNT2, CNT3, amax, bmax, gmax, J, K, M, JJ, KK, MM
        INTEGER nmin, nmax, nnmin, nnmax, aa, bb, gg, ROT1, ROT2, jmax
        INTEGER x, y

        PARAMETER (jmax=5)
        PARAMETER (y=INT((4.0/3.0)*jmax**3+4*jmax**2+(11.0/3.0)*jmax+1))

        DOUBLE PRECISION Hamil(y,y), Rot(3,3), Jac(3), NewJac(3)
        DOUBLE PRECISION comHH(3), NewcomHH(3)
        DOUBLE PRECISION COM(3), PlnPt(3), H2Lst(12,3), H2(3)
        DOUBLE PRECISION thresh, pHpH, norm, normm, E1
        DOUBLE PRECISION wa, wb, wg, alpha, beta, gama
        DOUBLE PRECISION LilD, LilDD, E, EH2, BetaVal(20,2)
        DOUBLE PRECISION Dist, Rydberg(16,8), NewAlpha, NewBeta

        DOUBLE COMPLEX II, BigD, BigDD, JKM1, JKM2, sum1

        open(UNIT=0, FILE='Barreto', STATUS='OLD')
        read(0,*)
        do i = 1, 16
        read(0,*) Rydberg(i, 1:8)
        enddo
        close(0)

        pi = 4.d0 * ATAN(1.d0)
        II = (0.d0,1.0d0)

        pHpH = 3.79
        thresh = 1.0E-8

        amax = 20
        bmax = 20
        gmax = 20

        wa = 2 * pi / amax
        wg = 2 * pi / amax

        Jac(1) = 1.d0
        Jac(2) = 0.d0
        Jac(3) = 0.d0
        NewJac = 0.d0

        Hamil = 0.d0
        Rot = 0.d0

        comHH(1) = 0.d0
        comHH(2) = -0.4532164d0
        comHH(3) = 0.d0
        NewcomHH = 0.d0

        COM = 0.d0

        PlnPt(1) = 0.d0
        PlnPt(2) = 1.d0
        PlnPt(3) = 0.d0


c        BetaVal(1,1)  = -0.93246951
c        BetaVal(1,2)  = 0.17132449
c        BetaVal(2,1)  = -0.66120939
c        BetaVal(2,2)  = 0.36076157
c        BetaVal(3,1)  = -0.23861919
c        BetaVal(3,2)  = 0.46791393
c        BetaVal(4,1)  = 0.23861919
c        BetaVal(4,2)  = 0.46791393
c        BetaVal(5,1)  = 0.66120939
c        BetaVal(5,2)  = 0.36076157
c        BetaVal(6,1)  = 0.93246951
c        BetaVal(6,2)  = 0.17132449


        BetaVal(1,1)  = -0.9931285991850949
        BetaVal(1,2)  = 0.017614007139153273
        BetaVal(2,1)  = -0.9639719272779138
        BetaVal(2,2)  = 0.04060142980038622
        BetaVal(3,1)  = -0.9122344282513258
        BetaVal(3,2)  = 0.06267204833410944
        BetaVal(4,1)  = -0.8391169718222188
        BetaVal(4,2)  = 0.08327674157670467
        BetaVal(5,1)  = -0.7463319064601508
        BetaVal(5,2)  = 0.10193011981724026
        BetaVal(6,1)  = -0.636053680726515
        BetaVal(6,2)  = 0.11819453196151825
        BetaVal(7,1)  = -0.5108670019508271
        BetaVal(7,2)  = 0.13168863844917653
        BetaVal(8,1)  = -0.37370608871541955
        BetaVal(8,2)  = 0.14209610931838187
        BetaVal(9,1)  = -0.2277858511416451
        BetaVal(9,2)  = 0.14917298647260366
        BetaVal(10,1) = -0.07652652113349734
        BetaVal(10,2) = 0.15275338713072578
        BetaVal(11,1) = 0.07652652113349734
        BetaVal(11,2) = 0.15275338713072578
        BetaVal(12,1) = 0.2277858511416451
        BetaVal(12,2) = 0.14917298647260366
        BetaVal(13,1) = 0.37370608871541955
        BetaVal(13,2) = 0.14209610931838187
        BetaVal(14,1) = 0.5108670019508271
        BetaVal(14,2) = 0.13168863844917653
        BetaVal(15,1) = 0.636053680726515
        BetaVal(15,2) = 0.11819453196151825
        BetaVal(16,1) = 0.7463319064601508
        BetaVal(16,2) = 0.10193011981724026
        BetaVal(17,1) = 0.8391169718222188
        BetaVal(17,2) = 0.08327674157670467
        BetaVal(18,1) = 0.9122344282513258
        BetaVal(18,2) = 0.06267204833410944
        BetaVal(19,1) = 0.9639719272779138
        BetaVal(19,2) = 0.04060142980038622
        BetaVal(20,1) = 0.9931285991850949
        BetaVal(20,2) = 0.017614007139153273

        H2Lst(1,1) = -pHpH
        H2Lst(1,2) = 0.d0
        H2Lst(1,3) = 0.d0

        H2Lst(2,1) = -pHpH/2.d0
        H2Lst(2,2) = -sqrt(3.d0) * pHpH / 2.d0
        H2Lst(2,3) = 0.d0

        H2Lst(3,1) = pHpH/2.d0
        H2Lst(3,2) = -sqrt(3.d0) * pHpH / 2.d0
        H2Lst(3,3) = 0.d0

        H2Lst(4,1) = pHpH
        H2Lst(4,2) = 0.d0
        H2Lst(4,3) = 0.d0

        H2Lst(5,1) = pHpH/2.d0
        H2Lst(5,2) = sqrt(3.d0) * pHpH / 2.d0
        H2Lst(5,3) = 0.d0

        H2Lst(6,1) = -pHpH/2.d0
        H2Lst(6,2) = sqrt(3.d0) * pHpH / 2.d0
        H2Lst(6,3) = 0.d0

        H2Lst(7,1) = 0.d0
        H2Lst(7,2) = -(sqrt(pHpH**2-(pHpH/2)**2)-sqrt(3.d0)*pHpH/6.d0)
        H2Lst(7,3) = sqrt(6.d0)*pHpH/3.d0

        H2Lst(8,1) = -pHpH/2.d0
        H2Lst(8,2) = sqrt(3.d0) * pHpH / 6.d0
        H2Lst(8,3) = sqrt(6.d0) * pHpH / 3.d0

        H2Lst(9,1) = pHpH/2.d0
        H2Lst(9,2) = sqrt(3.d0) * pHpH / 6.d0
        H2Lst(9,3) = sqrt(6.d0) * pHpH / 3.d0

        H2Lst(10,1) = 0.d0
        H2Lst(10,2) = -(sqrt(pHpH**2-(pHpH/2)**2)-sqrt(3.d0)*pHpH/6.d0)
        H2Lst(10,3) = -sqrt(6.d0) * pHpH/3.d0

        H2Lst(11,1) = -pHpH/ 2.d0
        H2Lst(11,2) = sqrt(3.d0) * pHpH / 6.d0
        H2Lst(11,3) = -sqrt(6.d0) * pHpH / 3.d0

        H2Lst(12,1) = pHpH / 2.d0
        H2Lst(12,2) = sqrt(3.d0) * pHpH / 6.d0
        H2Lst(12,3) = -sqrt(6.d0) * pHpH / 3.d0


        CNT1 = 1
        CNT2 = 1

        do 10 JJ = 0, jmax
        do 10 MM = -JJ, JJ
        do 10 KK = -JJ, JJ

        do 50 J = 0, jmax
        do 50 M = -J, J
        do 50 K = -J, J

        Norm  = sqrt((2*J+1) / (8*pi**2))
        Normm = sqrt((2*JJ+1) / (8*pi**2))

        sum1 = 0.d0

        nnmin = max(0, (KK-MM))
        nnmax = min((JJ-MM), (JJ+KK))

        nmin = max(0, (K-M))
        nmax = min((J-M), (J+K))

        do 20 aa = 1, amax
        do 20 bb = 1, bmax
        do 20 gg = 1, gmax

        NewJac = 0.d0
        NewcomHH = 0.d0

        E = 0.d0
        EH2 = 0.d0

        LilD  = 0.d0
        LilDD = 0.d0
        BigD  = 0.d0
        BigDD = 0.d0

        alpha = (2*pi*(aa-1))/amax
        gama  = (2*pi*(gg-1))/gmax
        beta = acos(BetaVal(bb,1))
        wb = BetaVal(bb,2)

        Rot(1,1) = cos(alpha)*cos(beta)*cos(gama)-sin(alpha)*sin(gama)
        Rot(1,2) = sin(alpha)*cos(beta)*cos(gama)+cos(alpha)*sin(gama)
        Rot(1,3) = -sin(beta) * cos(gama)
        Rot(2,1) = -cos(alpha)*cos(beta)*sin(gama)-sin(alpha)*cos(gama)
        Rot(2,2) = -sin(alpha)*cos(beta)*sin(gama)+cos(alpha)*cos(gama)
        Rot(2,3) = sin(beta) * sin(gama)
        Rot(3,1) = cos(alpha) * sin(beta)
        Rot(3,2) = sin(alpha) * sin(beta)
        Rot(3,3) = cos(beta)

        do 30 ROT1 = 1,3
        do 30 ROT2 = 1,3

        NewJac(ROT1) = NewJac(ROT1) + Jac(ROT2) * Rot(ROT1,ROT2)
        NewcomHH(ROT1)  = NewcomHH(ROT1) + comHH(ROT2) * ROT(ROT1,ROT2)

30      continue

        call LittleD(J, K, M, beta, lilD)
        call LittleD(JJ, KK, MM, beta, LilDD)

        BigD = exp(-(M*alpha + K*gama) * II) * LilD
        BigDD = exp(-(MM*alpha + KK*gama) * II ) * LilDD

c        print *, BigD

        JKM1 = CONJG(BigDD) * Normm
        JKM2 = BigD * Norm

        weight = wa * wb * wg

        do 40 CNT3 = 2, 12

        E1 = 0.d0

        H2(1) = H2Lst(CNT3, 1)
        H2(2) = H2Lst(CNT3, 2)
        H2(3) = H2Lst(CNT3, 3)

        call Angle(NewcomHH, COM, PlnPt, NewAlpha)
        call Angle(NewJac, COM, H2, NewBeta)
        call Distance(H2, COM, Dist)
        call Potential(NewAlpha, NewBeta, Dist, Rydberg, E1)

        EH2 = EH2 + E1

40      continue

        E = weight * JKM1 * JKM2 * EH2

        sum1 = sum1 + E

20      continue

        if (abs(sum1) .LT. thresh) then
        sum1 = 0.d0
        endif

        Hamil(CNT1, CNT2) = sum1

        CNT1 = CNT1 + 1

50      continue

        CNT1 = 1
        CNT2 = CNT2 + 1

10      continue

        do x = 1, y

        write(20, *) Hamil(1:y, x)

        enddo

        end

*********************************************************************
*
*********************************************************************


        subroutine littleD(j, k, m, beta, lild)

        integer j, k, m, n, nmin, nmax
        INTEGER jm1, jm2, jk1, jk2, jmn, jkn, nmk
        DOUBLE PRECISION lild, lilw, bigw, beta
        DOUBLE PRECISION prod1, prod2, prod3, prod4
        DOUBLE PRECISION prod5, prod6, prod7, prod8

        nmin = max(0, (k-m))
        nmax = min((j-m), (j+k))

        lild = 0.d0

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

        call fac(jm1, prod1)
        call fac(jm2, prod2)
        call fac(jk1, prod3)
        call fac(jk2, prod4)
        call fac(jmn, prod5)
        call fac(jkn, prod6)
        call fac(nmk, prod7)
        call fac(  n, prod8)

        lilw = sqrt(prod1 * prod2 * prod3 * prod4) / (prod5 * prod6 *
     + prod7 * prod8)

        bigW = lilw * ((cos(beta/2))**(2*j+k-m-2*n)) *
     + ((-sin(beta/2))**(m-k+2*n))

        lild = lild + ((-1)**n) * bigW

10      continue

        end

*********************************************************************
*
*********************************************************************

        subroutine Potential(a, B2, R, mat, Etot)

        DOUBLE PRECISION R, Etot, mat(16,8)
        DOUBLE PRECISION va000, va022, vb000, vb022, vc000, vc022
        DOUBLE PRECISION VHG, VHP, VHL, VXG, VXP, VZL, VZG, VZP, VXL
        DOUBLE PRECISION VTaG, VTaP, VTbL, VTbP, VTbG, VL
        DOUBLE PRECISION Va, Vb, Vc, B2, pi, a1, a2, a3, a
        INTEGER I, J, K

        pi = 4 * ATAN(1.d0)

        a1 = 0.25 + 0.5*cos(a) + 0.25*cos(2.0*a)
        a2 = 0.5  - 0.5* cos(2.0*a)
        a3 = 0.25 - 0.5*cos(a) + 0.25*cos(2.0*a)


        call Vind(mat, 1, VHG, R)
        call Vind(mat, 2, VHP, R)
        call Vind(mat, 3, VHL, R)
        call Vind(mat, 4, VXG, R)
        call Vind(mat, 5, VXP, R)
        call Vind(mat, 9, VXL, R)
        call Vind(mat, 10, VTaG, R)
        call Vind(mat, 11, VTaP, R)
        call Vind(mat, 12, VTbL, R)
        call Vind(mat, 13, VTbP, R)
        call Vind(mat, 14, VTbG, R)
        call Vind(mat, 15, VL, R)

        va000 = (1.0/9.0) * (2.0*VHL + VL + 2.0*VTaG + 2.0*VTbL+2.0*VXL)
        vb000 = (1.0/9.0) * (2.0*VHP + VL + 2.0*VTaP + 2.0*VTbP+2.0*VXP)
        vc000 = (1.0/9.0) * (2.0*VHG + VL + 2.0*VTaG + 2.0*VTbG+2.0*VXG)

        va022 = (-2*sqrt(5.0)/45.0) * (VHL - VL - 2.0*VTaG + VTbL + VXL)
        vb022 = (-2*sqrt(5.0)/45.0) * (VHP - VL - 2.0*VTaP + VTbP + VXP)
        vc022 = (-2*sqrt(5.0)/45.0) * (VHG - VL - 2.0*VTaG + VTbG + VXG)


        va000 = va000 * a1
        vb000 = vb000 * a2
        vc000 = vc000 * a3

        va022 = va022 * a1
        vb022 = vb022 * a2
        vc022 = vc022 * a3

        Va = va000 + va022 * (sqrt(5.0)/4.0) * (3.0*cos(B2) + 1)

        Vb = vb000 + vb022 * (sqrt(5.0)/4.0) * (3.0*cos(B2) + 1)

        Vc = vc000 + vc022 * (sqrt(5.0)/4.0) * (3.0*cos(B2) + 1)

        Etot = Va + Vb + Vc
        Etot = Etot * 349.75

        end

*********************************************************************
*
*********************************************************************

        subroutine Vind(mat, K, V1, R)

        DOUBLE PRECISION sum1, R, V1, mat(16,8)
        INTEGER I, J, K

        sum1 = 0.d0

        do i = 1, 5

        sum1 = sum1 + mat(K,i) * (R - mat(K,8)) ** i

        enddo

        V1 = -mat(K,7) * (1 + sum1) * exp(-mat(K,1) * (R - mat(K,8))) +
     + mat(K,6)


10      continue

        end


*********************************************************************
*
*********************************************************************

        subroutine fac(n, nn)

        integer n, i
        DOUBLE PRECISION nn

        nn = 1.d0

        do i = 1, n
        nn = nn * real(i)
        enddo

        end

*********************************************************************
*
*********************************************************************

        subroutine angle(P1, P2, P3, angle1)

        DOUBLE PRECISION P1(3), P2(3), P3(3), angle1, d1, d2, d3, Val

        call distance(P1, P2, d1)
        call distance(P1, P3, d2)
        call distance(P2, P3, d3)

        Val = (d1**2 + d3**2 - d2**2)/(2*d1*d3)

        if (Val .LT. -1.0) then
                Val = -1.0
        end if

        if (Val .GT. 1.0) then
                Val = 1.0
        end if

        angle1 = acos(Val)

        end

*********************************************************************
*
*********************************************************************

        subroutine distance(P1, P2, dist)

        DOUBLE PRECISION P1(3), P2(3), dist, A, B, C

        A = (P1(1) - P2(1))**2
        B = (P1(2) - P2(2))**2
        C = (P1(3) - P2(3))**2

        dist = sqrt(A + B + C)

        end
