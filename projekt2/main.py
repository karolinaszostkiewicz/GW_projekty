import numpy as np
import math as m
import plotly.express as px

# współrzędne Oslo
phi = 59.5446
lambd = 10.4416

# współrzędne Mbandaka
# phi = 0.03
# lambd = 18.16

# współrzędne Sydney
# phi = 33.52
# lambd = 151.12

# rekstascenzja
r = 21.54472

# deklinacja
d = -5.47611

# zamiana jednostki rekstascenzji na radiany
rectanstention = np.deg2rad(r)

# zamiana jednostki deklinacji na radiany
declination = np.deg2rad(d)

# zamiana daty na liczbę dni juliańskich
def GDtoJD(year, month, day, hour):
    if month <= 2:
        year = year - 1
        month = month + 12
    jd = m.floor(365.25*(year + 4716)) + m.floor(30.6001 * (month + 1)) + day + hour / 24 - 1537.5
    return jd

# konwersja liczby dni juliańskich na  Greenwich Mean Sidereal Time
def JDtoGMST(year, month, day, hour):
    T = (GDtoJD(year, month, day, hour) - 2451545) / 36525
    g = 280.46061837 + 360.98564736629 * (GDtoJD(year, month, day, 0) - 2451545) + 0.000387933 * (T ** 2) - (T ** 3) \
        / 38710000
    g = g % 360
    return g

# obliczenie kąta godzinnego
def hour_angle(year, month, day, hour, lambd, alpha):
    g = JDtoGMST(year, month, day, 0)
    UT1 = hour * 1.002737909350795
    S = UT1 * 15 + lambd + g
    t = S - alpha * 15
    return t

# rozwiązanie trójkąta paralaktycznego
def star_azimuth(declination, t, phi):
    phi = np.deg2rad(phi)
    t = np.deg2rad(t)
    numeral = (-(m.cos(declination)) * m.sin(t))
    denom = (m.cos(phi) * m.sin(declination)) - (m.sin(phi) * m.cos(declination) * m.cos(t))
    tgA = np.arctan(numeral / denom)
    tgA = np.rad2deg(tgA)

    if (denom < 0):
        tgA += 180
    elif (numeral < 0):
        tgA += 360
    return tgA

# obliczenie odległości zenitalnej
def zenith_distance(phi, declination, t):
    phi = np.deg2rad(phi)
    t = np.deg2rad(t)
    cos_z = (m.sin(phi) * m.sin(declination)) + (m.cos(phi) * m.cos(declination) * m.cos(t))
    cos_z = np.rad2deg(np.arccos(cos_z))
    return cos_z

# transformacja współrzędnych
def transform(cos_z, tgA):
    cos_z = np.deg2rad(cos_z)
    tgA = np.deg2rad(tgA)
    x = 1 * m.sin(cos_z) * m.cos(tgA)
    y = 1 * m.sin(cos_z) * m.sin(tgA)
    z = 1 * m.cos(cos_z)
    x = np.rad2deg(x)
    y = np.rad2deg(y)
    z = np.rad2deg(z)
    return x, y, z

tgA_array = []
cos_z_array = []

for hour in range(24):
    hA = hour_angle(2022, 5, 1, hour, lambd, rectanstention)
    tgA = star_azimuth(declination, hA, phi)
    cos_z = zenith_distance(phi, declination, hA)

    tgA_array.append(tgA)
    cos_z_array.append(cos_z)

x_array = []
y_array = []
z_array = []

for i in range(24):
    x, y, z = transform(cos_z_array[i], tgA_array[i])

    x_array.append(x)
    y_array.append(y)
    z_array.append(z)

fig = px.scatter_3d(x= x_array, y= y_array, z= z_array, title='Ruch gwiazdy Sadalsuud')
fig.show()



