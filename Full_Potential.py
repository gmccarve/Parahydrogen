import numpy as np
import numpy.linalg as LA
import math
import cmath

#Transpose   : numpy.matrix.transpose(a)
#Multiply    : numpy.matmul(a,b,out)
#Facotrial   : math.factorial(x)
#Diagonalize : numpy.linalg.eig(a)

def fac(x):
    return math.factorial(x)
def cos(x):
    return math.cos(x)
def sin(x):
    return math.sin(x)
def sqrt(x):
    return math.sqrt(x)
def exp(x):
    try:
        return math.exp(x)
    except TypeError:
        return cmath.exp(x)
def acos(x):
    return math.acos(x)

def BarretoData():
    
    Param = np.array([[1.82110830,      -0.68911643,     0.68958598,    -0.0838483647,       0.0126925781,      -3.01419762E-3,         0.18829563,     3.49094606],
                      [1.99840472,	-4.14679440,	 0.69591413,	-0.71312315,        -0.0753054887, 	-1.83877249E-4,		0.00615866,	4.05762898],	
                      [1.82316255,	-2.60212946,	 0.08061900,	 0.65030665,	    -0.0185807957,  	 2.27992315E-3,		0.01575833,	4.06533481],	
                      [1.57801257,	-1.12109621,	 0.80923636,	-0.21623355,	     0.0298999320, 	-1.56642595E-3,		0.35351933,	3.27925573],	
                      [2.35343205,	-0.25199551,	 0.25808430,	-0.0241827175,	    -0.0234774706, 	 2.35364459E-4,		0.03199163,	3.83107720],
                      [2.07730185,	-4.01676002,	 0.64116182,	-0.79382607,	    -0.172110569,	 2.13561279E-3,		0.01982560,	3.66315169],
                      [2.47546373,	-0.94799971,	-0.05054013,	-0.0181169571, 	     0.00947442526, 	 3.32408991E-3,		0.11474494,	3.59636723],
                      [2.20681611,	-0.22328958,	 0.41897355,	 0.12197924,	    -1.83310670E-9,	-2.33746769E-3,		0.35324305,	3.20970135],
                      [1.62963793,	-1.24069967,	 0.92213412,	-0.26527515,	     0.0394768933, 	-1.43091488E-3,		0.46673446,	3.06560853],
                      [1.77462974,	-1.44105532,	 1.01737286,	-0.28604112,	     0.0415765453,	-3.87181069E-4,		0.31942666,	3.22447005],
                      [2.31330927,	-0.36203143,	 0.47312858,	 0.16232889,	    -9.25865897E-9,	-2.37631164E-3,		0.41366228,	3.13952742],	
                      [1.99428335,	-8.70304888,	 0.04377305,	 1.62873895,	     0.395524611,	 4.41058620E-3,		0.01430489,	3.85496247],
                      [1.64451047,	-1.02382763,	 0.62537402,	-0.15048903,	     0.0218892795, 	-7.04592527E-4,		0.41209066,	3.21488250],
                      [1.60351286,	-1.06192799,	 0.77429747,	-0.20795977,	     0.0326880350,	-2.81735049E-3,		0.60250604,	3.11203092],
                      [2.26214546,	-6.36633810,	 4.30604438,	-0.74791939,	    -0.351109419,	 4.97531282E-3,		0.00487686,	3.99151721],
                      [2.10998036,	-2.25756329,	 0.41264115,	-0.27152359,	    -0.130527787,	 2.88966024E-3,		0.06034381,	3.62719190]])
    return Param

def Potential(alpha, Beta2, R):
    
    #Vmat = np.zeros((16))
    v000 = np.zeros((3))
    v022 = np.zeros((3))
    V    = np.zeros((3))

    w1 = 0.25 + 0.5*cos(alpha) + 0.25*cos(2*alpha)
    w2 = 0.5  - 0.5*cos(2*alpha)
    w3 = 0.25 - 0.5*cos(alpha) + 0.25*cos(2*alpha)

    #for K in range(0,16):
    #    Vmat[K] = Rydberg(K, R)

    v000[0] = ((1.0/9.0) * (2.0*Vmat[2] + Vmat[14] + 2.0*Vmat[9]  + 2.0*Vmat[11] + 2.0*Vmat[8])) * w1
    v000[1] = ((1.0/9.0) * (2.0*Vmat[1] + Vmat[14] + 2.0*Vmat[10] + 2.0*Vmat[12] + 2.0*Vmat[4])) * w2
    v000[2] = ((1.0/9.0) * (2.0*Vmat[0] + Vmat[14] + 2.0*Vmat[9]  + 2.0*Vmat[13] + 2.0*Vmat[3])) * w3

    v022[0] = ((-2.0*sqrt(5.0)/45.0) * (Vmat[2] - Vmat[14] - 2.0*Vmat[9]  + Vmat[11] + Vmat[8])) * w1
    v022[1] = ((-2.0*sqrt(5.0)/45.0) * (Vmat[1] - Vmat[14] - 2.0*Vmat[10] + Vmat[12] + Vmat[4])) * w2
    v022[2] = ((-2.0*sqrt(5.0)/45.0) * (Vmat[0] - Vmat[14] - 2.0*Vmat[9]  + Vmat[13] + Vmat[3])) * w3

    V[0] = v000[0] + (v022[0] * (sqrt(5.0)/4.0) * (3.0*cos(Beta2) + 1))
    V[1] = v000[1] + (v022[1] * (sqrt(5.0)/4.0) * (3.0*cos(Beta2) + 1))
    V[2] = v000[2] + (v022[2] * (sqrt(5.0)/4.0) * (3.0*cos(Beta2) + 1))

    Etot = V[0] + V[1] + V[2]
    Etot = Etot * 349.75
    
    return Etot

def Rydberg(K, R):
    
    Val = 0.0
    for x in range(1,6):
        Val = Val + (Ryd[K][x-1] * (R - Ryd[K][7]) ** x)
    V_r = -Ryd[K][6] * (1 + Val) * math.exp(-Ryd[K][0] * (R - Ryd[K][7])) + Ryd[K][5]
    return V_r

def NoField():

    A = 27.8761
    B = 14.5074
    C = 9.2877

    IULC = 0

    RotMat = np.zeros((sze,sze))

    for J in range(0,jmax+1):
        for M in range(-J, J+1):
            for K in range(-J, J+1):
                for KK in range(-J, J+1):

                    val = 0.0

                    if K == KK:
                        val = ((0.50) * (B+C) * (J*(J+1)-K**2) + A*K**2)
                    elif K + 2 == KK:
                        val = ((0.25) * (B-C) * sqrt((J*(J+1) - K * (K+1))) * sqrt(J*(J+1) - (K+1) * (K+2)))
                    elif K - 2 == KK:
                        val = ((0.25) * (B-C) * sqrt((J*(J+1) - K * (K-1))) * sqrt(J*(J+1) - (K-1) * (K-2)))

                    try:
                        RotMat[IULC+KK+J, IULC+K+J] = val
                    except IndexError:
                        pass

            IULC = IULC + 2*J+1
                    
    return RotMat

def FullField(Field):

    Hamil  = np.zeros((sze,sze))
    Rot    = np.zeros((3,3))
    thresh = 1E-10
    pHpH = 3.79

    Jac    = np.array((1.0,0.0,0.0))
    NewJac = np.zeros((3))
    comHH  = np.array((0.0,-0.4532164, 0.0))
    newcomHH = np.zeros((3))
    COM = np.zeros((3))
    PlanePoint = np.array((0.0,1.0,0.0))

    H2List = np.zeros((12,3))

    H2List[0]  = np.array((-pHpH,   0,               0))
    H2List[1]  = np.array((-pHpH/2, -sqrt(3)*pHpH/2, 0))
    H2List[2]  = np.array((pHpH/2,  -sqrt(3)*pHpH/2, 0))
    H2List[3]  = np.array((pHpH,    0,               0))
    H2List[4]  = np.array((pHpH/2,  sqrt(3)*pHpH/2,  0))
    H2List[5]  = np.array((-pHpH/2, sqrt(3)*pHpH/2,  0))
    H2List[6]  = np.array((0,       -(sqrt(pHpH**2 - (pHpH/2)**2) - sqrt(3)*pHpH/6),  sqrt(6)*pHpH/3))
    H2List[7]  = np.array((-pHpH/2, sqrt(3)*pHpH/6,  sqrt(6)*pHpH/3))
    H2List[8]  = np.array((pHpH/2,  sqrt(3)*pHpH/6,  sqrt(6)*pHpH/3))  
    H2List[9] = np.array((0,       -(sqrt(pHpH**2 - (pHpH/2)**2) - sqrt(3)*pHpH/6), -sqrt(6)*pHpH/3))
    H2List[10] = np.array((-pHpH/2, sqrt(3)*pHpH/6,  -sqrt(6)*pHpH/3))
    H2List[11] = np.array((pHpH/2,  sqrt(3)*pHpH/6,  -sqrt(6)*pHpH/3))

    count1 = 0
    count2 = 0

    amax = 6
    bmax = 6
    gmax = 6

    wa = 2*pi/amax
    wg = 2*pi/gmax

    BetaV, BetaW = (np.polynomial.legendre.leggauss(bmax))

    for JJ in range(0, jmax+1):
        for KK in range(-JJ, JJ+1):
            for MM in range(-JJ, JJ+1):
                #print (JJ, KK, MM)
                Normm = sqrt((2*JJ+1)/(8*pi**2))

                for J in range(0, jmax+1):
                    for K in range(-J, J+1):
                        for M in range(-J, J+1):

                            Norm = sqrt((2*J+1)/(8*pi**2))

                            sum1 = 0.0

                            nnmin = max(0, (KK-MM))
                            nnmax = min((JJ-MM), (JJ+KK))

                            nmin = max(0, (K-M))
                            nmax = min((J-M), (J+K))

                            for aa in range(1,amax+1):
                                for bb in range(1,bmax+1):
                                    for gg in range(1,gmax+1):

                                        NewJac = np.zeros((3))
                                        newcomHH = np.zeros((3))

                                        E = 0.0
                                        E_H2 = 0.0

                                        LilD = 0.0
                                        LilDD = 0.0

                                        BigD = 0.0
                                        BigDD = 0.0

                                        alpha = (2*pi*(aa-1))/amax
                                        gamma = (2*pi*(gg-1))/gmax
                                        beta = acos(BetaV[bb-1])
                                        wb = BetaW[bb-1]

                                        Rot[0,0] =  cos(alpha) * cos(beta) * cos(gamma) - sin(alpha) * sin(gamma)
                                        Rot[0,1] =  sin(alpha) * cos(beta) * cos(gamma) + cos(alpha) * sin(gamma)
                                        Rot[0,2] = -sin(beta)  * cos(gamma)
                                        Rot[1,0] = -cos(alpha) * cos(beta) * sin(gamma) - sin(alpha) * cos(gamma)
                                        Rot[1,1] = -sin(alpha) * cos(beta) * sin(gamma) + cos(alpha) * cos(gamma)
                                        Rot[1,2] =  sin(beta)  * sin(gamma)
                                        Rot[2,0] =  cos(alpha) * sin(beta)
                                        Rot[2,1] =  sin(alpha) * sin(beta)
                                        Rot[2,2] =  cos(beta)

                                        for rot1 in range(0,3):
                                            for rot2 in range(0,3):

                                                NewJac[rot1] = NewJac[rot1] + Jac[rot2] * Rot[rot1,rot2]
                                                newcomHH[rot1] = newcomHH[rot1] + comHH[rot2] * Rot[rot1,rot2]

                                        LilD  = LittleD(J, K, M, beta)
                                        LilDD = LittleD(JJ, KK, MM, beta)

                                        BigD = exp(-(M*alpha + K*gamma)*II) * LilD
                                        BigDD = exp(-(MM*alpha + KK*gamma)*II) * LilDD

                                        JKM1 = BigDD.conjugate() * Normm
                                        JKM2 = BigD * Norm
                                        weight = wa * wb * wg
                                        
                                        for item in H2List:
                    
                                            NewAlpha = Angle(newcomHH, COM, PlanePoint)
                                            NewBeta = Angle(NewJac, COM, item)
                                            Dist = Distance(item, COM)
                                            E_H2 = E_H2 +  Potential(NewAlpha, NewBeta, Dist)

                                        #E = weight * JKM1 * JKM2 * 1 * Field * sin(beta) * cos(gamma) #***Stark Effect***#
                                        E = weight * JKM1 * JKM2 * E_H2

                                        sum1 = sum1 + E
                            

                            if abs(sum1) < thresh:
                                sum1 = 0

                            Hamil[count1, count2] = round(sum1.real,6)

                            count1 += 1

                count1 = 0
                count2 += 1
                                   
    return Hamil

def LittleD(J, K, M, beta):

    lilD = 0

    nmin = max(0, (K-M))
    nmax = min((J-M), (J+K))

    for N in range(nmin, nmax+1):

        lilw = sqrt(fac(J+M) * fac(J-M) * fac(J+K) * fac(J-K)) / (fac(J-M-N) * fac(J+K-N) * fac(N+M-K) * fac(N))
        BigW = lilw * ((cos(beta/2.0))**(2*J+K-M-2*N)) * ((-sin(beta/2.0))**(M-K+2*N))

        lilD = lilD + ((-1)**N) * BigW

    return lilD

def Angle(a, b, c):

    d1 = Distance(a,b)
    d2 = Distance(a,c)
    d3 = Distance(b,c)

    Val = (d1**2 + d3**2 - d2**2)/(2*d1*d3)

    if Val > 1.0:
        Val = 1.0
    elif Val < -1.0:
        Val = -1.0

    angle1 = math.acos(Val)

    return angle1

def Distance(a, b):

    A = (a[0] - b[0])**2
    B = (a[1] - b[1])**2
    C = (a[2] - b[2])**2

    return sqrt(A + B + C)

def Stark():

    file = open("C:/Users/gavin/OneDrive/Scripts/Water/val.txt", "w")

    for i in range(0,100):
        i = i / 100
        E = LA.eigvals((FullField(i) + NoField()))
        #Mat = np.sort(E,axis=None)
        file.write(str(i) + '  ')
        for item in E:
            file.write(str(item) + '  ')
        file.write('\n')
        print (i)

    file.close()

    return


#np.set_printoptions(precision=5, linewidth=99999999)
np.set_printoptions(formatter={'float': lambda x: "{0:0.8f}".format(x)},linewidth=99999999)

II = complex(0,1)
pi = math.pi
Ryd = BarretoData()
jmax = 1
sze = int((4.0/3.0) * jmax**3 + 4*jmax**2 + (11.0/3.0) * jmax + 1.0000001)

Vmat = np.zeros((16))
for K in range(0,16):
    Vmat[K] = Rydberg(K, 3.79)

Stark()


#NoField = NoField()

#print (NoField)

#print (LA.eig(NoField))
#FullField = FullField()

#EigVal, EigVec = LA.eig(FullField)

#print (FullField)

#np.save("C:/Users/gavin/OneDrive/Scripts/Water/Matrices/NF_6.npy", NoField)
#np.save("C:/Users/gavin/OneDrive/Scripts/Water/NoField2.npy", NoField)

#FullField = np.load("C:/Users/gavin/OneDrive/Scripts/Water/FullField2.npy")
#NoField = np.load("C:/Users/gavin/OneDrive/Scripts/Water/NoField2.npy")

#FullField = np.load("C:/Users/gavin/OneDrive/Scripts/Water/J=2/FullField.npy")
#NoField = np.load("C:/Users/gavin/OneDrive/Scripts/Water/J=2/NoField.npy")

"""file = open("C:/Users/gavin/OneDrive/Scripts/Water/Values2.txt", "w")

for i in range(0,3):
    #file.write(str(i) + '\n')
    i = (i)/2
    FullMat = i * FullField + NoField

    EigVal, EigVec = LA.eig(FullMat)
    EigVal = np.sort(EigVal)

#    EigVec[abs(EigVec) < 1E-11] = 0.0
    EigVec = abs(EigVec)
    
    file.write(str(i) + '\n')
    for item in EigVal:
        file.write('\t' + str(item) + '\n')
    file.write('\n')

    #file.write(str(i) + '  ' + str(EigVal[0]) + '  ' + str(EigVal[1]) + '  ' + str(EigVal[2]) + '  ' + str(EigVal[3]) + '  ' + str(EigVal[4]) + '  ' +
    #                           str(EigVal[5]) + '  ' + str(EigVal[6]) + '  ' + str(EigVal[7]) + '  ' + str(EigVal[8]) + '  ' + str(EigVal[9]))
    #file.write('\n')


file.close()"""









   






