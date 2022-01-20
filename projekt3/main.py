import math as m

#parametry dla GRS80
a = 6378137
e2 = 0.00669437999013

#dane wierzcholkow A, B, C, D w radianach
fiA = m.radians(52.00)
lambdaA = m.radians(20.45)

fiB = m.radians(51.75)
lambdaB = m.radians(20.75)

fiC = m.radians(52.00)
lambdaC = m.radians(21.25)

fiD = m.radians(51.75)
lambdaD = m.radians(21.25)

def vincent(lambda1, lambda2, fi1, fi2):

    b = a * m.sqrt(1 - e2)
    f = 1 - b / a

    L1 = lambda2 - lambda1

    U1 = m.atan((1 - f) * m.tan(fi1))
    U2 = m.atan((1 - f) * m.tan(fi2))

    L2 = L1

    while True:
        sinS = m.sqrt((m.cos(U2) * m.sin(L2)) ** 2 + (m.cos(U1) * m.sin(U2) - m.sin(U1) * m.cos(U2) * m.cos(L2)) ** 2)
        cosS = m.sin(U1) * m.sin(U2) + m.cos(U1) * m.cos(U2) * m.cos(L2)
        Sigma = m.atan(sinS/cosS)
        sinA = (m.cos(U1) * m.cos(U2) * m.sin(L2))/sinS
        cosA2 = 1 - (sinA) ** 2
        cos2Sm = cosS - ((2 * m.sin(U1) * m.sin(U2))/cosA2)
        C = f/16 * cosA2 * (4 + f * (4 - 3 * cosA2))

        L = L1 + (1 - C) * f * sinA * (Sigma + C * sinS * (cos2Sm + C * (-1 + 2 * ((cos2Sm) ** 2))))

        if (L - L2) < m.radians(0.000001/3600):
            break
        else:
            L2 = L

        u2 = (((a ** 2) - (b ** 2))/(b ** 2)) * cosA2
        A = 1 + (u2/16384) * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)))
        B = (u2/1024) * (256 + u2 * (-128 + u2 * (74 - 47 * u2)))
        deltaSigma = B * sinS * (cos2Sm + (B * (cosS * (-1 + 2 * (cos2Sm ** 2))))/4 - (B * cos2Sm * (-3 + 4 * (sinS ** 2)) * (-3 + 4 * (cos2Sm ** 2)))/6)
        s_12 = b * A * (Sigma - deltaSigma)
        AZ_12 = m.atan((m.cos(U2) * m.sin(L))/((m.cos(U1) * m.sin(U2)) - (m.sin(U1) * m.cos(U2) * m.cos(L))))
        AZ_21 = m.atan((m.cos(U1) * m.sin(L))/(((-m.sin(U1)) * m.cos(U2)) + (m.cos(U1) * m.sin(U2) * m.cos(L)))) + m.pi

        #poprawka azymutow
        y1 = m.cos(U2) * m.sin(L2)
        x1 = m.cos(U1) * m.sin(U2) - m.sin(U1) * m.cos(U2) * m.cos(L2)
        if (y1 > 0 and x1 > 0):
            AZ_12 = m.atan(y1/x1)
        elif (y1 > 0 and x1 < 0):
            AZ_12 = m.atan(y1/x1) + m.pi
        elif (y1 < 0 and x1 < 0):
            AZ_12 = m.atan(y1/x1) + m.pi
        elif (y1 < 0 and x1 > 0):
            AZ_12 = m.atan(y1/x1) + 2 * m.pi

        y2 = m.cos(U1) * m.sin(L2)
        x2 = (-m.sin(U1)) * m.cos(U2) + m.cos(U1) * m.sin(U2) * m.cos(L2)
        if (y2 > 0 and x2 > 0):
            AZ_21 = m.atan(y2/x2) + m.pi
        elif (y2 > 0 and x2 < 0):
            AZ_21 = m.atan(y2/x2) + 2 * m.pi
        elif (y2 < 0 and x2 < 0):
            AZ_21 = m.atan(y2/x2) + 2 * m.pi
        elif (y2 < 0 and x2 > 0):
            AZ_21 = m.atan(y2/x2) + 3 * m.pi

        return s_12, AZ_12, AZ_21

def kivioji(lambda1, fi1, AZ_1, s_12):
    s = s_AD/2
    n = round(s/1000)
    ds = s/n
    for i in range(n):
        M = (a * (1 - e2))/(m.sqrt((1 - e2 * (m.sin(fi1)) ** 2) ** 3))
        N = a/(m.sqrt(1 - e2 * (m.sin(fi1)) ** 2))
        dfi = (m.cos(AZ_1) * ds)/M
        dAZ = (m.sin(AZ_1) * m.tan(fi1) * ds)/N
        pfi = fi1 + (dfi/2)
        Am = AZ_1 + (dAZ/2)
        Nm = a/(m.sqrt(1 - e2 * (m.sin(pfi) ** 2)))
        Mm = (a * (1 - e2))/(m.sqrt((1 - e2 * (m.sin(pfi) ** 2)) ** 3))
        dfm = m.cos(Am) * ds/Mm
        dlm = m.sin(Am) * ds/(Nm * m.cos(pfi))
        dAm = m.sin(Am) * m.tan(pfi) * ds/Nm

        fi1 = fi1 + dfm
        lambda1 = lambda1 + dlm
        AZ_1 = AZ_1 + dAm

    fiP = fi1
    lambdaP = lambda1
    AZ_P = AZ_1

    return  fiP, lambdaP, AZ_P

def pole(lambda1, lambda2, fi1, fi2):
    b = a * m.sqrt(1 - e2)
    e = m.sqrt(e2)
    p1 = b ** 2 * (lambda2 - lambda1) / 2 * \
            ((m.sin(fi1) / (1 - e2 * (m.sin(fi1)) ** 2) + 1 / (2 * e) * m.log((1 + e * m.sin(fi1)) / (1 - e * m.sin(fi1))))
             - (m.sin(fi2) / (1 - e2 * (m.sin(fi2)) ** 2) + 1 / (2 * e) * m.log((1 + e * m.sin(fi2)) / (1 - e * m.sin(fi2)))))
    p = p1/1000000

    return p

def zamiana(s_dz):
    stopnie = int(s_dz)
    minuty = int((s_dz - stopnie) * 60)
    sekundy = str("{:.5f}".format((s_dz - minuty/60 - stopnie) * 3600))
    stopnie, minuty = str(stopnie), str(minuty)

    return stopnie + '° ' + minuty + "' " + sekundy + "'' "

fi_r = (fiA + fiD)/2
fi_s = m.degrees(fi_r)
lamda_r = (lambdaA + lambdaD)/2
lamda_s = m.degrees(lamda_r)

s_AD = vincent(lambdaA, lambdaD, fiA, fiD)[0]
AZ_AD = vincent(lambdaA, lambdaD, fiA, fiD)[1]
AZ_DA = vincent(lambdaA, lambdaD, fiA, fiD)[2]

fiP = kivioji(lambdaA, fiA, AZ_AD, s_AD)[0]
lambdaP = kivioji(lambdaA, fiA, AZ_AD, s_AD)[1]
az_pom = m.degrees(AZ_AD)

fi_srod = m.degrees(fiP)
lam_srod = m.degrees(lambdaP)

s_p = vincent(lamda_r, lambdaP, fi_r, fiP)[0]
s_P = "{:.5}".format(s_p)
az_p = vincent(lamda_r, lambdaP, fi_r, fiP)[1]
az_p2 = vincent(lamda_r, lambdaP, fi_r, fiP)[2]

az_1 = m.degrees(az_p)
az_2 = m.degrees(az_p2)

pole_ABCD = pole(lambdaA, lambdaD, fiA, fiD)

print("Punkt średniej szerokości: (" + zamiana(fi_s) + ", " + zamiana(lamda_s) + ")")
print("Współrzędne punktu środkowego AD: (" + zamiana(fi_srod) + ", " + zamiana(lam_srod) + ")")
print("Różnica odległości pomiędzy tymi punktami wynosi " + str(s_P) + " metrów.")
print("Azymuty w tych punktach: " + zamiana(az_1) + ", " + zamiana(az_2))
print("Pole powierzchni czworokąta wynosi " + str(pole_ABCD) + " kilometrów kwadratowych.")