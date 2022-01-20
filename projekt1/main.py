import math as m
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


#współrzędne lotniska Chopina
wsp_fi, wsp_lambda, wsp_h = 52.1786, 20.9559, 213

e2 = 0.00669438002290

#grs80 a parameter
A = 6378137

def x_y_z(fi, lam, h, A, e2):
    fi = np.deg2rad(fi)
    lam = np.deg2rad(lam)
    N = A / (1 - (e2) * m.sin(fi) ** 2) ** 0.5

    x = (N + h) * m.cos(fi) * m.cos(lam)
    y = (N + h) * m.cos(fi) * m.sin(lam)
    z = (N * (1 - e2) + h * m.sin(fi))
    return x, y, z


def geo2neu(F1, L1, H1, F2, L2, H2):
    pkt_1 = x_y_z(F1, L1, H1, A, e2)
    pkt_2 = x_y_z(F2, L2, H2, A, e2)

    F1 = np.deg2rad(F1)
    L1 = np.deg2rad(L1)
    R = np.array(
                 [[-m.sin(F1) * m.cos(L1), -m.sin(L1), m.cos(F1) * m.cos(L1)],
                  [-m.sin(F1) * m.sin(L1), m.cos(L1), m.cos(F1) * m.sin(L1)],
                  [m.cos(F1), 0, m.sin(F1)]
                 ])

    R = R.transpose()
    x = np.array([[pkt_2[0] - pkt_1[0]],
                  [pkt_2[1] - pkt_1[1]],
                  [pkt_2[2] - pkt_1[2]]])
    neu = R @ x
    return neu

# wczytywanie danych
# lot WARSAW, POLAND do PUNTA CANA, DOMINIKANA REPUBLIC
lot = np.loadtxt("lot1.txt")
dane_lotu = pd.DataFrame(lot)

#wyświetlanie trasy lotu w fi lambda i h
x = lot[:, 1]
y = lot[:, 0]

plt.plot(x, y, color='b')
plt.show()


def azymut(up, down, funkcja):
    az = funkcja(up/down)
    if down > 0 and up > 0:
        az = np.rad2deg(az)
    elif down < 0 and up > 0:
        az = np.rad2deg(az + m.pi)
    elif down < 0 and up < 0:
        az = np.rad2deg(az + m.pi)
    else:
        az = np.rad2deg(az + 2 * m.pi)
    if az > 360:
        az -= 360
    elif az < 0:
        az += 360
    return az

def zamiana(s_dz):
    stopnie = int(s_dz)
    minuty = int((s_dz - stopnie) * 60)
    sekundy = str("{:.5f}".format((s_dz - minuty/60 - stopnie) * 3600))
    stopnie, minuty = str(stopnie), str(minuty)

    return stopnie + '° ' + minuty + "' " + sekundy + "'' "

#konwersja na neu
n = []
e = []
u = []

tan_a = []
s = []
cos_z = []


for i in lot:
     i = geo2neu(i[0], i[1], i[2], wsp_fi, wsp_lambda, wsp_h)

     n.append(i[0][0])
     e.append(i[1][0])
     u.append(i[2][0])

print("Azymut A, odległość skośna,  Kąt zenitalny:")
for i in range(1, len(n)):
    os = (n[i]**2 + e[i]**2 + u[i]**2)**0.5

    print(str(zamiana(azymut(e[i], n[i], np.arctan))) + "   " + str(round(os, 3)) + "   " + str(zamiana(azymut(u[i], os, np.arccos))))

fig = plt.figure()
ax = plt.axes(projection="3d")
ax.scatter3D(n, e, u, c=u, cmap='magma')

for x, y, z in zip(n, e, u):
    if z < 0:
        ax.scatter(x, y, z, "blue", s= 200 )
        ax.text(x, y, z, '%s' % ("                  znika za horyzontem"), size=10, zorder=1, color='black')
        print("Samolot zniknął za horyzontem: ", round(x, 3), round(y, 3), round(z, 3))
        break

ax.text(n[0], e[0], u[0], '%s' % ("Punta Cana"), size=10, zorder=1, color='black')
ax.text(n[-1], e[-1], u[-1], '%s' % ("Warszawa"), size=10, zorder=1, color='black')
plt.show()
