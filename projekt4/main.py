import math as m
from shapely.geometry.polygon import Polygon

#parametry dla GRS80
a = 6378137
e2 = 0.00669437999013

# dane wierzcholkow A, B, C, D w radianach
fiA = m.radians(52.00)
lambdaA = m.radians(20.45)
fiB = m.radians(51.75)
lambdaB = m.radians(20.75)
fiC = m.radians(52.00)
lambdaC = m.radians(21.25)
fiD = m.radians(51.75)
lambdaD = m.radians(21.25)
fiS = m.radians(51.87581)
lambdaS = m.radians(20.850700)
fiSS = m.radians(51.87500)
lambdaSS = m.radians(20.85000)

b = a * m.sqrt(1 - e2)
e = m.sqrt(e2)
ee2 = ((m.pow(a, 2) - m.pow(b, 2)) / m.pow(b, 2))

Wfi = [fiA, fiB, fiC, fiD, fiS, fiSS]
Wlambda = [lambdaA, lambdaB, lambdaC, lambdaD, lambdaS, lambdaSS]
etykieta = ["Wsp punktu A:", "Wsp punktu B:", "Wsp punktu C:", "Wsp punktu D:", "Wsp punktu środkowego:", "Wsp średniej szerokości:"]

def to_GK(lambda1, fi1):
    L0 = 19 * (m.pi/180)
    N = a/(m.sqrt(1 - e2 * (m.sin(fi1))**2))
    A0 = 1 - (e2/4) - ((3*(e2)**2)/64) - ((5*(e2)**3)/256)
    A2 = (3/8) * (e2 + ((e2)**2)/4 + ((15*(e2)**3)/128))
    A4 = (15/256) * ((e2)**2 + (3*(e2)**3)/4)
    A6 = ((35 * (e2)**3)/3072)
    t = m.tan(fi1)
    L = lambda1 - L0
    n2 = ee2 * (m.cos(fi1) ** 2)
    sigma = a * (A0 * fi1 - A2 * m.sin(2 * fi1) + A4 * m.sin(4 * fi1) - A6 * m.sin(6 * fi1))

    x_GK = sigma + (L**2/2) * N * m.sin(fi1) * m.cos(fi1) * (1 + (L**2/12) * (m.cos(fi1)**2) * (5 - (t**2) + 9*n2 + 4*(n2**2)) + ((L**4)/360) * (m.cos(fi1)**4) * (61 - 58*(t**2) + (t**4) + 270*n2 - 330*n2*(t**2)))
    y_GK = L * N * m.cos(fi1) * (1 + ((L**2)/6) * (m.cos(fi1)**2) * (1 - (t**2) + n2) + ((L**4)/120) * (m.cos(fi1)**4) * (5 - 18*(t**2) + (t**4) + 14*n2 - 58*n2*(t**2)))

    return x_GK, y_GK

def to_92(x_GK, y_GK):
    m0 = 0.9993
    x_92 = m0 * x_GK - 5300000
    y_92 = m0 * y_GK + 500000
    return x_92, y_92

def to_2000(lambda1, fi1):
    if m.degrees(lambda1) < 16.5:
        strefa = 15
    elif 16.5 <= m.degrees(lambda1) < 19.5:
        strefa = 18
    elif 19.5 <= m.degrees(lambda1) < 22.5:
        strefa = 21
    elif m.degrees(lambda1) >= 22.5:
        strefa = 24

    L0 = strefa * (m.pi/180)
    N = a / (m.sqrt(1 - e2 * (m.sin(fi1)) ** 2))
    A0 = 1 - (e2 / 4) - ((3 * (e2) ** 2) / 64) - ((5 * (e2) ** 3) / 256)
    A2 = (3 / 8) * (e2 + ((e2) ** 2) / 4 + ((15 * (e2) ** 3) / 128))
    A4 = (15 / 256) * ((e2) ** 2 + (3 * (e2) ** 3) / 4)
    A6 = ((35 * (e2) ** 3) / 3072)
    t = m.tan(fi1)
    L = lambda1 - L0
    n2 = ee2 * (m.cos(fi1) ** 2)
    sigma = a * (A0 * fi1 - A2 * m.sin(2 * fi1) + A4 * m.sin(4 * fi1) - A6 * m.sin(6 * fi1))

    x_GK = sigma + (L ** 2 / 2) * N * m.sin(fi1) * m.cos(fi1) * (1 + (L ** 2 / 12) * (m.cos(fi1) ** 2) * (5 - (t ** 2) + 9 * n2 + 4 * (n2 ** 2)) + ((L ** 4) / 360) * (m.cos(fi1) ** 4) * (61 - 58 * (t ** 2) + (t ** 4) + 270 * n2 - 330 * n2 * (t ** 2)))
    y_GK = L * N * m.cos(fi1) * (1 + ((L ** 2) / 6) * (m.cos(fi1) ** 2) * (1 - (t ** 2) + n2) + ((L ** 4) / 120) * (m.cos(fi1) ** 4) * (5 - 18 * (t ** 2) + (t ** 4) + 14 * n2 - 58 * n2 * (t ** 2)))

    m0 = 0.999923
    nr_strefy = strefa/3
    x_2000 = m0 * x_GK
    y_2000 = m0 * y_GK + 1000000 * nr_strefy + 500000

    return x_2000, y_2000

wsp_GK = []
wsp_92 = []
wsp_2000 = []

print("Układ GK:")
for i in range (0, len(etykieta)):
    print(etykieta[i], f"{round(to_GK(Wlambda[i], Wfi[i])[0], 3):.3f}",
          f"{round(to_GK(Wlambda[i], Wfi[i])[1], 3):.3f}")
    wsp_GK.extend([to_GK(Wlambda[i], Wfi[i])[0], to_GK(Wlambda[i], Wfi[i])[1]])

print("Układ 2000:")
for i in range (0, len(etykieta)):
    print(etykieta[i], f"{round(to_2000(Wlambda[i], Wfi[i])[0], 3):.3f}",
          f"{round(to_2000(Wlambda[i], Wfi[i])[1], 3):.3f}")
    wsp_2000.extend([to_2000(Wlambda[i], Wfi[i])[0],to_2000(Wlambda[i], Wfi[i])[1]])

print("Układ 1992:")
for i in range (0, len(etykieta)):
    print(etykieta[i],f"{round(to_92(to_GK(Wlambda[i], Wfi[i])[0], to_GK(Wlambda[i], Wfi[i])[1])[0], 3):.3f}",
          round(to_92(to_GK(Wlambda[i], Wfi[i])[0], to_GK(Wlambda[i], Wfi[i])[1])[1], 3))
    wsp_92.extend([to_92(to_GK(Wlambda[i], Wfi[i])[0], to_GK(Wlambda[i], Wfi[i])[1])[0],
                  to_92(to_GK(Wlambda[i], Wfi[i])[0], to_GK(Wlambda[i], Wfi[i])[1])[1]])

def pole(lambda1, lambda2, fi1, fi2):
    p = (b**2 * (lambda2 - lambda1) / 2) * (((m.sin(fi1) / (1 - e2 * (m.sin(fi1)**2))) + (1 / (2 * e)) * m.log(
        (1 + e * m.sin(fi1)) / (1 - e * m.sin(fi1))))
                                        - ((m.sin(fi2) / (1 - e2 * (m.sin(fi2)**2))) + (1 / (2 * e)) * m.log(
        (1 + e * m.sin(fi2)) / (1 - e * m.sin(fi2)))))

    return p

p_GK = Polygon([(wsp_GK[0], wsp_GK[1]), (wsp_GK[2], wsp_GK[3]), (wsp_GK[6], wsp_GK[7]), (wsp_GK[4], wsp_GK[5])])
p_92 = Polygon([(wsp_92[0], wsp_92[1]), (wsp_92[2], wsp_92[3]), (wsp_92[6], wsp_92[7]), (wsp_92[4], wsp_92[5])])
p_2000 = Polygon([(wsp_2000[0], wsp_2000[1]), (wsp_2000[2], wsp_2000[3]), (wsp_2000[6], wsp_2000[7]), (wsp_2000[4], wsp_2000[5])])

print("POLA:")
print("Pole elipsoidalne:", round(pole(lambdaA, lambdaD, fiA, fiD), 12)/1000000, "km2")
print("Pole GK: ", f"{round(p_GK.area,6)/1000000:.12f}", "km2")
print("Pole 2000: ", f"{round(p_2000.area,6)/1000000:.12f}", "km2")
print("Pole 1992: ", f"{round(p_92.area,6)/1000000:.12f}", "km2")

def u92_to_geo(x_92, y_92):
    m0 = 0.9993
    x_GK = ((x_92 + 5300000)/m0)
    y_GK = ((y_92 - 500000)/m0)
    L0 = m.radians(19)

    epsilon = 9999
    A0 = (1 - (e2 / 4) - ((3 * (e2 ** 2)) / 64) - ((5 * m.pow(e2, 3)) / 256))
    A2 = 0.375 * (e2 + (e2 ** 2 / 4) + ((15 * m.pow(e2, 3)) / 128))
    A4 = (15 / 256) * (e2 ** 2 + ((3 * m.pow(e2, 3)) / 4))
    A6 = ((35 * m.pow(e2, 3)) / 3072)
    fi0 = x_GK/(a * A0)

    while epsilon > ((0.00001 / 3600) * (m.pi / 180)):
        fi1 = (fi0 + (x_GK - (a * ((A0 * fi0) - (A2 * m.sin(2 * fi0)) + (A4 * m.sin(4 * fi0)) - (
                    A6 * m.sin(6 * fi0))))) / (a * A0))
        epsilon = abs(fi1 - fi0)
        fi0 = fi1

    n2 = (ee2 * m.pow(m.cos(fi1), 2))
    N = a / (m.sqrt(1 - e2 * m.pow(m.sin(fi1), 2)))
    M = (a * (1 - e2)) / (m.sqrt(pow(1 - e2 * pow(m.sin(fi1), 2), 3)))
    t = m.tan(fi1)

    fi2 = fi1 - (((y_GK ** 2) * t) / (2 * M * N)) * (1 - ((y_GK ** 2) / (12 * (N ** 2))) *
                                                     (5 + 3 * (t ** 2) + n2 - 9 * n2 * (t ** 2) - 4 * (n2 ** 2)) +
                                                     ((y_GK ** 4) / (360 * (N ** 4))) * (61 + 90 * (t ** 2) + 45 * (t ** 4)))

    lambda2 = L0 + (y_GK / (N * m.cos(fi1))) * (1 - ((y_GK ** 2) / (6 * (N ** 2))) * (1 + 2 * (t ** 2) + n2) +
                                                ((y_GK ** 4) / (120 * (N ** 4))) *
                                                (5 + 28 * (t ** 2) + 24 * (t ** 4) + 6 * n2 + 8 * n2 * (t ** 2)))

    fi3 = str(m.degrees(fi2))
    lambda3 = str(m.degrees(lambda2))

    return fi3, lambda3, x_GK, y_GK

def u2000_to_geo(x_2000, y_2000):
    m0 = 0.999923

    if y_2000 < 6000000:
        nr_strefy = 5
    elif 6000000 <= y_2000 < 7000000:
        nr_strefy = 6
    elif 7000000 <= y_2000 < 8000000:
        nr_strefy = 7
    elif y_2000 >= 8000000:
        nr_strefy = 8

    x_GK = x_2000/m0
    y_GK = (y_2000 - (nr_strefy * 1000000) - 500000) / m0

    l0 = m.radians(nr_strefy * 3)

    epsilon = 9999
    A0 = (1 - (e2 / 4) - ((3 * (e2 ** 2)) / 64) - ((5 * m.pow(e2, 3)) / 256))
    A2 = 0.375 * (e2 + (e2 ** 2 / 4) + ((15 * m.pow(e2, 3)) / 128))
    A4 = (15 / 256) * (e2 ** 2 + ((3 * m.pow(e2, 3)) / 4))
    A6 = ((35 * m.pow(e2, 3)) / 3072)
    fi0 = x_GK / (a * A0)

    while epsilon > ((0.00001 / 3600) * (m.pi / 180)):
        fi1 = (fi0 + (x_GK - (a * ((A0 * fi0) - (A2 * m.sin(2 * fi0)) + (A4 * m.sin(4 * fi0)) - (
                    A6 * m.sin(6 * fi0))))) / (a * A0))
        epsilon = abs(fi1 - fi0)
        fi0 = fi1

    N = a / (m.sqrt(1 - e2 * m.pow(m.sin(fi1), 2)))
    M = (a * (1 - e2)) / (m.sqrt(pow(1 - e2 * pow(m.sin(fi1), 2), 3)))
    t = m.tan(fi1)
    n2 = (ee2 * m.pow(m.cos(fi1), 2))

    fi2 = fi1 - (((y_GK ** 2) * t) / (2 * M * N)) * (1 - ((y_GK ** 2) / (12 * (N ** 2))) * (
                5 + 3 * (t ** 2) + n2 - 9 * n2 * (t ** 2) - 4 * (n2 ** 2)) + (
                                                                   (y_GK ** 4) / (360 * (N ** 4))) * (
                                                                   61 + 90 * (t ** 2) + 45 * (t ** 4)))

    lambda2 = l0 + (y_GK / (N * m.cos(fi1))) * (
                1 - ((y_GK ** 2) / (6 * (N ** 2))) * (1 + 2 * (t ** 2) + n2) + ((y_GK ** 4) / (120 * (N ** 4))) * (
                    5 + 28 * (t ** 2) + 24 * (t ** 4) + 6 * n2 + 8 * n2 * (t ** 2)))

    fi3 = m.degrees(fi2)
    lambda3 = m.degrees(lambda2)
    return fi3, lambda3, x_GK, y_GK

def skala_znieksztalcenia_GK(x_GK, y_GK):
    L0 = m.radians(19)
    epsilon = 9999
    A0 = (1 - (e2 / 4) - ((3 * (e2 ** 2)) / 64) - ((5 * m.pow(e2, 3)) / 256))
    A2 = 0.375 * (e2 + (e2**2 / 4) + ((15 * m.pow(e2, 3)) / 128))
    A4 = (15/256) * (e2**2 + ((3 * m.pow(e2, 3)) / 4))
    A6 = ((35 * m.pow(e2, 3)) / 3072)
    fi0 = x_GK / (a * A0)

    while epsilon > ((0.00001 / 3600) * (m.pi / 180)):
        fi1 = (fi0 + (x_GK - (a * ((A0 * fi0) - (A2 * m.sin(2 * fi0)) + (A4 * m.sin(4 * fi0)) - (A6 * m.sin(6 * fi0))))) / (a * A0))
        epsilon = abs(fi1 - fi0)
        fi0 = fi1

    N = a / (m.sqrt(1 - e2 * m.pow(m.sin(fi1), 2)))
    M = (a * (1 - e2)) / (m.sqrt(pow(1 - e2 * pow(m.sin(fi1), 2), 3)))
    t = m.tan(fi1)
    n2 = (ee2 * m.pow(m.cos(fi1), 2))

    fi2 = fi1 - (((y_GK**2) * t) / (2 * M * N)) * (1 - ((y_GK**2) / (12*(N**2))) * (5 + 3*(t**2) + n2 - 9*n2*(t**2) - 4*(n2**2)) + ((y_GK**4) / (360*(N**4))) * (61 + 90*(t**2) + 45*(t**4)))
    lambda2 = L0 + (y_GK / (N * m.cos(fi1))) * (1 - ((y_GK**2) / (6*(N**2))) * (1 + 2*(t**2) + n2) + ((y_GK**4) / (120*(N**4))) * (5 + 28*(t**2) + 24*(t**4) + 6*n2 + 8*n2*(t**2)))

    N1 = a / (m.sqrt(1 - e2 * m.pow(m.sin(fi2), 2)))
    M1 = (a * (1 - e2)) / (m.sqrt(pow(1 - e2 * pow(m.sin(fi2), 2), 3)))
    R = m.sqrt(N1 * M1)

    m_uklad = (1 + (y_GK ** 2 / (2 * (R ** 2))) + (y_GK ** 4 / (24 * (R ** 4))))
    kappa = 1 - m_uklad
    kappa_km = kappa * 1000

    skala_pole = (m_uklad**2)
    kappa_pole = 1 - skala_pole
    kappa_pole_ha = kappa_pole * 10000

    return m_uklad, kappa_km, skala_pole, kappa_pole_ha

print("Skala zniekształcenia GK:")
print("Punkt A: " + str(skala_znieksztalcenia_GK(to_GK(lambdaA, fiA)[0], to_GK(lambdaA, fiA)[1])))
print("Punkt B: " + str(skala_znieksztalcenia_GK(to_GK(lambdaB, fiB)[0], to_GK(lambdaB, fiB)[1])))
print("Punkt C: " + str(skala_znieksztalcenia_GK(to_GK(lambdaC, fiC)[0], to_GK(lambdaC, fiC)[1])))
print("Punkt D: " + str(skala_znieksztalcenia_GK(to_GK(lambdaD, fiD)[0], to_GK(lambdaD, fiD)[1])))
print("Punkt środkowy: " + str(skala_znieksztalcenia_GK(to_GK(lambdaS, fiS)[0], to_GK(lambdaS, fiS)[1])))
print("Punkt średniej szerokości: " + str(skala_znieksztalcenia_GK(to_GK(lambdaSS, fiSS)[0], to_GK(lambdaSS, fiSS)[1])))

def skala_znieksztalcenia(uklad, x, y):
    if uklad == 1992:
        m0 = 0.9993
        fi1, lambda1, x_GK, y_GK = u92_to_geo(x, y)

    elif uklad == 2000:
        m0 = 0.999923
        fi1, lambda1, x_GK, y_GK = u2000_to_geo(x, y)

    N = a / (m.sqrt(1 - e2 * m.pow(m.sin(float(fi1)), 2)))
    M = (a * (1 - e2)) / (m.sqrt(pow(1 - e2 * pow(m.sin(float(fi1)), 2), 3)))
    R = m.sqrt(M * N)
    m_m = (1 + (y_GK**2 / (2*(R**2))) + (y_GK**4 / (24*(R**4))))
    m_uklad = m_m * m0
    kappa = 1 - m_uklad
    kappa_km = kappa * 1000

    m_m2 = (m_m ** 2)
    m0_2 = (m0 ** 2)
    skala_m2 = m_m2 * m0_2
    kappa_pole = 1 - skala_m2
    kappa_pole_ha = kappa_pole * 10000

    return m_uklad, kappa_km, skala_m2, kappa_pole_ha

print("Skala zniekształcenia 1992:")
print("Punkt A: " + str(skala_znieksztalcenia(1992, (to_92(to_GK(lambdaA, fiA)[0], to_GK(lambdaA, fiA)[1])[0], to_92(to_GK(lambdaA, fiA)[0], to_GK(lambdaA, fiA)[1])[1])[0],
                                              (to_92(to_GK(lambdaA, fiA)[0], to_GK(lambdaA, fiA)[1])[0], to_92(to_GK(lambdaA, fiA)[0], to_GK(lambdaA, fiA)[1])[1])[1])))
print("Punkt B: " + str(skala_znieksztalcenia(1992, (to_92(to_GK(lambdaB, fiB)[0], to_GK(lambdaB, fiB)[1])[0], to_92(to_GK(lambdaB, fiB)[0], to_GK(lambdaB, fiB)[1])[1])[0],
                                              (to_92(to_GK(lambdaB, fiB)[0], to_GK(lambdaB, fiB)[1])[0], to_92(to_GK(lambdaB, fiB)[0], to_GK(lambdaB, fiB)[1])[1])[1])))
print("Punkt C: " + str(skala_znieksztalcenia(1992, (to_92(to_GK(lambdaC, fiC)[0], to_GK(lambdaC, fiC)[1])[0], to_92(to_GK(lambdaC, fiC)[0], to_GK(lambdaC, fiC)[1])[1])[0],
                                              (to_92(to_GK(lambdaC, fiC)[0], to_GK(lambdaC, fiC)[1])[0], to_92(to_GK(lambdaC, fiC)[0], to_GK(lambdaC, fiC)[1])[1])[1])))
print("Punkt D: " + str(skala_znieksztalcenia(1992, (to_92(to_GK(lambdaD, fiD)[0], to_GK(lambdaD, fiD)[1])[0], to_92(to_GK(lambdaD, fiD)[0], to_GK(lambdaD, fiD)[1])[1])[0],
                                              (to_92(to_GK(lambdaD, fiD)[0], to_GK(lambdaD, fiD)[1])[0], to_92(to_GK(lambdaD, fiD)[0], to_GK(lambdaD, fiD)[1])[1])[1])))
print("Punkt środkowy: " + str(skala_znieksztalcenia(1992, (to_92(to_GK(lambdaS, fiS)[0], to_GK(lambdaS, fiS)[1])[0], to_92(to_GK(lambdaS, fiS)[0], to_GK(lambdaS, fiS)[1])[1])[0],
                                              (to_92(to_GK(lambdaS, fiS)[0], to_GK(lambdaS, fiS)[1])[0], to_92(to_GK(lambdaS, fiS)[0], to_GK(lambdaS, fiS)[1])[1])[1])))
print("Punkt Średniej szerokości: " + str(skala_znieksztalcenia(1992, (to_92(to_GK(lambdaSS, fiSS)[0], to_GK(lambdaSS, fiSS)[1])[0], to_92(to_GK(lambdaSS, fiA)[0], to_GK(lambdaSS, fiSS)[1])[1])[0],
                                              (to_92(to_GK(lambdaSS, fiSS)[0], to_GK(lambdaSS, fiSS)[1])[0], to_92(to_GK(lambdaSS, fiSS)[0], to_GK(lambdaSS, fiSS)[1])[1])[1])))

print("Skala zniekształcenia 2000:")
print("Punkt A: " + str(skala_znieksztalcenia(2000, (to_2000(lambdaA, fiA)[0], to_2000(lambdaA, fiA)[1])[0], (to_2000(lambdaA, fiA)[0], to_2000(lambdaA, fiA)[1])[1])))
print("Punkt B: " + str(skala_znieksztalcenia(2000, (to_2000(lambdaB, fiB)[0], to_2000(lambdaB, fiB)[1])[0], (to_2000(lambdaB, fiB)[0], to_2000(lambdaB, fiB)[1])[1])))
print("Punkt C: " + str(skala_znieksztalcenia(2000, (to_2000(lambdaC, fiC)[0], to_2000(lambdaC, fiC)[1])[0], (to_2000(lambdaC, fiC)[0], to_2000(lambdaC, fiC)[1])[1])))
print("Punkt D: " + str(skala_znieksztalcenia(2000, (to_2000(lambdaD, fiD)[0], to_2000(lambdaD, fiD)[1])[0], (to_2000(lambdaD, fiD)[0], to_2000(lambdaD, fiD)[1])[1])))
print("Punkt środkowy: " + str(skala_znieksztalcenia(2000, (to_2000(lambdaS, fiS)[0], to_2000(lambdaS, fiS)[1])[0], (to_2000(lambdaS, fiS)[0], to_2000(lambdaS, fiS)[1])[1])))
print("Punkt średniej szerokości: " + str(skala_znieksztalcenia(2000, (to_2000(lambdaSS, fiSS)[0], to_2000(lambdaSS, fiSS)[1])[0], (to_2000(lambdaSS, fiSS)[0], to_2000(lambdaSS, fiSS)[1])[1])))