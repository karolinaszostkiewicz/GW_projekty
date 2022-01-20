import math as m
import numpy as np

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

pkA = [fiA, lambdaA, 100]
pkB = [fiB, lambdaB, 100]
pkC = [fiC, lambdaC, 100]
pkD = [fiD, lambdaD, 100]
pkS = [fiS, lambdaS, 100]
pkSS = [fiSS, lambdaSS, 100]

#parametry GRS80
e2g = 0.00669438002290
ag = 6378137

# parametry Krasowskiego
ak = 6378245
e2k = 0.0066934215520398155
x0 = -33.4297
y0 = 146.5746
z0 = 76.2865
alfa = -0.35867 / 3600
b = -0.05283 / 3600
kappa = 0.8407728e-6
gamma = 0.84354 / 3600

def zamiana(s_dz):
    stopnie = int(s_dz)
    minuty = int((s_dz - stopnie) * 60)
    sekundy = str("{:.5f}".format((s_dz - minuty/60 - stopnie) * 3600))
    stopnie, minuty = str(stopnie), str(minuty)

    return stopnie + '° ' + minuty + "' " + sekundy + "'' "

def xyz(pkt, ag, e2g):
    fi = pkt[0]
    lam = pkt[1]
    h = pkt[2]
    N = ag / (1 - (e2g) * m.sin(fi) ** 2) ** 0.5
    x = (N + h) * m.cos(fi) * m.cos(lam)
    y = (N + h) * m.cos(fi) * m.sin(lam)
    z = (N * (1 - e2g) + h) * m.sin(fi)
    return round(x, 3), round(y, 3), round(z, 3)

def hirv(pkt, a, e2):
    x = pkt[0]
    y = pkt[1]
    z = pkt[2]
    r = m.sqrt(x ** 2 + y ** 2)
    fi0 = m.atan((z/r) * ((1-e2) ** -1))
    N = a / m.sqrt(1 - e2 * m.sin(fi0) ** 2)
    h = (r / m.cos(fi0)) - N
    fi1 = m.atan((z / r) * (1 - e2 * (N / (N + h))) ** -1)
    epsilon = m.radians(0.00005/3600)

    while abs(fi1-fi0) > epsilon:
        fi0 = fi1
        n = a / m.sqrt(1 - e2 * m.sin(fi0) ** 2)
        h = (r / m.cos(fi0)) - n
        fi1 = m.atan((z / r) * (1 - e2 * (n / (n + h))) ** -1)

    lam1 = np.arctan(y/x)
    N = a / m.sqrt(1 - e2 * m.sin(fi1) ** 2)
    h = (r / m.cos(fi1)) - N

    fi1 = np.degrees(fi1)
    lam1= np.degrees(lam1)

    return zamiana(fi1), zamiana(lam1), h

def transform(pkt, x0, y0, z0, kappa, alfa, beta, gamma):
    alfa = m.radians(alfa)
    beta = m.radians(beta)
    gamma = m.radians(gamma)
    pkt_p = np.array([[pkt[0]], [pkt[1]], [pkt[2]]])
    matrix = np.array([[kappa, gamma, -beta],
                       [-gamma, kappa, alfa],
                       [beta, -alfa, kappa]])
    s = np.array([[x0], [y0], [z0]])
    pkt_t = pkt_p + matrix.dot(pkt_p) + s
    xt = pkt_t[0][0]
    yt = pkt_t[1][0]
    zt = pkt_t[2][0]

    return round(xt, 3), round(yt, 3), round(zt, 3)

print("XYZ W GRS80:")
print("Punkt A: xyz w Grs80: " + str(xyz(pkA, ag, e2g)))
print("Punkt B: xyz w Grs80: " + str(xyz(pkB, ag, e2g)))
print("Punkt C: xyz w Grs80: " + str(xyz(pkC, ag, e2g)))
print("Punkt D: xyz w Grs80: " + str(xyz(pkD, ag, e2g)))
print("Punkt środkowy: xyz w Grs80: " + str(xyz(pkS, ag, e2g)))
print("Punkt średniej szerokości: xyz w Grs80: " + str(xyz(pkSS, ag, e2g)))

print("XYZ W ELIPSOIDZIE KRASOWSKIEGO:")
print("Punkt A: xyz w elipsoidzie Krasowskiego: " + str(transform(xyz(pkA, ag, e2g), x0, y0, z0, kappa, alfa, b, gamma)))
print("Punkt B: xyz w elipsoidzie Krasowskiego: " + str(transform(xyz(pkB, ag, e2g), x0, y0, z0, kappa, alfa, b, gamma)))
print("Punkt C: xyz w elipsoidzie Krasowskiego: " + str(transform(xyz(pkC, ag, e2g), x0, y0, z0, kappa, alfa, b, gamma)))
print("Punkt D: xyz w elipsoidzie Krasowskiego: " + str(transform(xyz(pkD, ag, e2g), x0, y0, z0, kappa, alfa, b, gamma)))
print("Punkt środkowy: xyz w elipsoidzie Krasowskiego: " + str(transform(xyz(pkS, ag, e2g), x0, y0, z0, kappa, alfa, b, gamma)))
print("Punkt średniej szerokośći: xyz w elipsoidzie Krasowskiego: " + str(transform(xyz(pkSS, ag, e2g), x0, y0, z0, kappa, alfa, b, gamma)))

print("WSPÓŁRZĘDNE W ELIPSOIDZIE KRASOWSKIEGO:")
print("Punkt A: współrzędne w elipsoidzie Krasowskiego: " + str(hirv(transform(xyz(pkA, ag, e2g), x0, y0, z0, kappa, alfa, b, gamma), ak, e2k)))
print("Punkt B: współrzędne w elipsoidzie Krasowskiego: " + str(hirv(transform(xyz(pkB, ag, e2g), x0, y0, z0, kappa, alfa, b, gamma), ak, e2k)))
print("Punkt C: współrzędne w elipsoidzie Krasowskiego: " + str(hirv(transform(xyz(pkC, ag, e2g), x0, y0, z0, kappa, alfa, b, gamma), ak, e2k)))
print("Punkt D: współrzędne w elipsoidzie Krasowskiego: " + str(hirv(transform(xyz(pkD, ag, e2g), x0, y0, z0, kappa, alfa, b, gamma), ak, e2k)))
print("Punkt środkowy: współrzędne w elipsoidzie Krasowskiego: " + str(hirv(transform(xyz(pkS, ag, e2g), x0, y0, z0, kappa, alfa, b, gamma), ak, e2k)))
print("Punkt średniej szerokości: współrzędne w elipsoidzie Krasowskiego: " + str(hirv(transform(xyz(pkSS, ag, e2g), x0, y0, z0, kappa, alfa, b, gamma), ak, e2k)))