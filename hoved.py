import numpy as np
import numpy.matlib
from pandas import *


#Bruker studentnr til "Ole Bendik" = xxxx80 -> A = 8 && B = 0

# I-profil (flens, steg, høyde, bredde).
# Rør (radius, tykkelse)

profiler = [
[0.019, 0.0012, 0.6, 0.22],
[0.1, 0.02],
[0.0127, 0.008, 0.36, 0.17],
[0.1, 0.01]
]


def lesinput():

    # Åpner inputfilen
    fid = open("input3.txt", "r")

    # Leser totalt antall punkt
    npunkt = int(fid.readline())       # 'fid.readline()' leser en linje, 'int(...)' gjør at linjen tolkes som et heltall

    # LESER INN XY-KOORDINATER TIL KNUTEPUNKTENE
    # Nodenummer tilsvarer radnummer i "Node-variabel"
    # x-koordinat lagres i kolonne 1, y-koordinat i kolonne 2
    # Grensebetingelse lagres i kolonne 3, 1 = fast innspent og 0 = fri rotasjon
    punkt = np.loadtxt(fid, dtype = int, max_rows = npunkt)     # 'max_rows = npunkt' sorger for at vi bare leser
                                                                # de 'npunkt' neste linjene i tekstfilen

    # Leser antall elementer
    nelem = int(fid.readline())

    # Leser konnektivitet: Sammenheng mellom elementender og knutepunktnummer samt EI for elementene
    # Elementnummer tilsvarer radnummer i "Elem-variabel"
    # Knutepunktnummer for lokal ende 1 lagres i kolonne 1
    # Knutepunktnummer for lokal ende 2 lagres i kolonne 2
    # E-modul for materiale lagres i kolonne 3
    # Tverrsnittstype lagres i kolonne 4, I-profil = 1 og rørprofil = 2
    eleme = np.loadtxt(fid, dtype = np.int64, max_rows = nelem)

    # Leser antall laster som virker på rammen
    nlast = int(fid.readline())

    # Leser antall laster som ikke er ujevnt fordelt
    last = np.loadtxt(fid, dtype = np.float64, max_rows = nlast)
    lastvec = []
    # leser inn lastene og kategoriserer
    for elem in last:
        if elem[0]==1: # punktlast
            lastvec.append([1,elem[1],elem[2],elem[3],elem[4]])
        elif elem[0]==2: # moment
            lastvec.append([2,elem[1],elem[2]])
        elif elem[0]==3: # jevnt fordelt last
            lastvec.append([3,elem[1],elem[2]])

    # leser inn ujevnt fordelte laster
    nujamtlast = int(fid.readline())

    if nujamtlast > 0:
        ujevntfordeltlast = np.loadtxt(fid, dtype = np.float64, max_rows = nujamtlast)

        # Gjør sjekk på om det er kun 1 ujevnlast pga python array indeksering
        if nujamtlast == 1:
            q_start = ujevntfordeltlast[0]
            q_slutt = ujevntfordeltlast[1]
            q_current = q_start
            for i in range(2,ujevntfordeltlast.size):
                lastvec.append([3,ujevntfordeltlast[i],q_current])
                q_old = q_current
                q_current +=((q_slutt-q_start)/(ujevntfordeltlast.size-2))
                lastvec.append([4,ujevntfordeltlast[i],q_current-q_old])

        elif nujamtlast >= 2:
            for ujevn in ujevntfordeltlast:
                q_start = ujevn[0]
                q_slutt = ujevn[1]
                q_current = q_start
                for i in range(2, ujevn.size):
                    lastvec.append([3, ujevn[i], q_current])
                    q_old = q_current
                    q_current += ((q_slutt - q_start) / (ujevn.size - 2))
                    lastvec.append([4, ujevn[i], q_current - q_old])


    # Lukker input-filen
    fid.close()
    return npunkt, punkt, nelem, eleme, nlast, lastvec

def lengder(knutepunkt, element, nelem):
    elementlengder = np.matlib.zeros((nelem, 1))

    # Beregner elementlengder med Pythagoras' laeresetning
    for i in range (0, nelem):
        # OBS! Grunnet indekseringsyntaks i Python-arrays vil ikke denne funksjonen fungere naar vi bare har ett element.
        dx = knutepunkt[int((element[i][0]))][0] - knutepunkt[int((element[i][1]))][0]
        dy = knutepunkt[int((element[i][0]))][1] - knutepunkt[int((element[i][1]))][1]
        elementlengder[i] = np.sqrt((dx**2) + (dy**2))
    return elementlengder

def midtpunktsLaster(nelem, elementlengder, nlast, last, endeM, elem):
    midtLaster =  [0.0]*nelem
    for j in range(nelem):
        hasLoad = False
        for i in range(0, len(last)):
            if last[i][1] == j:
                if last[i][0] == 1: # Punktlast
                    P = last[i][3] * np.cos(last[i][4])
                    a = last[i][2]
                    b = float(elementlengder[int(last[i][1])] - last[i][2])
                    l = float(elementlengder[int(last[i][1])])
                    endeA = float((endeM[j][0] * b) / l)
                    endeB = float((endeM[j][1] * a) / l)

                    midtLaster[int(last[i][1])] = (P*a*b)/l + endeA + endeB
                    hasLoad = True

                elif last[i][0] == 3:  # Jevnlast
                    midtLaster[int(last[i][1])] += float(((last[i][2] * (elementlengder[int(last[i][1])])**2)/8) + \
                                                         (endeM[j][0] * 0.5) + (endeM[j][1] * 0.5))
                    hasLoad = True


                elif last[i][0] == 4:  # Ujevnlast
                    midtLaster[int(last[i][1])] += float(((last[i][2] * (elementlengder[int(last[i][1])])**2)/16))
                    hasLoad = True
        if hasLoad == False:
            midtLaster[j] += endeM[j][0] * 0.5 + endeM[j][1] * 0.5

    return midtLaster



def momentOgLastvektor(npunkt, punkt, nelem, elem, nlast, last, elementlengder):
    mom = np.array([[0 for x in range(2)] for y in range(nelem)]  ) #liste over fim med elementnummer som indeks og ende 1 og 2 som subindeks
    lastVektor = [0.0]*npunkt
    for tempLast in last:

        # Legger til negativt fastinnspenningsmoment i lastvektoren i tillegg

        if tempLast[0] == 1:    # Moment fra punktlast

            mom[int(tempLast[1])][0] -= float(((tempLast[3]*np.cos(tempLast[4])) * (tempLast[2]*(elementlengder[int(tempLast[1])]-tempLast[2])**2))/((elementlengder[int(tempLast[1])])**2))                   # fim for ende a
            lastVektor[int(elem[int(tempLast[1])][0])] += float(((tempLast[3] * np.cos(tempLast[4])) * (tempLast[2] * (elementlengder[int(tempLast[1])] - tempLast[2]) ** 2)) / (elementlengder[int(tempLast[1])]) ** 2)  # fim for ende a

            mom[int(tempLast[1])][1] += float(((tempLast[3]*np.cos(tempLast[4]))*((tempLast[2]**2)*(elementlengder[int(tempLast[1])]-tempLast[2]))) /(elementlengder[int(tempLast[1])])**2)                      # fim for ende b
            lastVektor[int(elem[int(tempLast[1])][1])] -= float(((tempLast[3] * np.cos(tempLast[4])) * ((tempLast[2] ** 2) * (elementlengder[int(tempLast[1])] - tempLast[2]))) / (elementlengder[int(tempLast[1])]) ** 2)  # fim for ende b

        elif tempLast[0] == 2:  # Konsentrerte momenter i knutepunkt
            lastVektor[int(tempLast[1])] += float(tempLast[2])

        elif tempLast[0] == 3:  # Moment fra jevn fordelt last
            mom[int(tempLast[1])][0] -= float((1/12)*tempLast[2]*(elementlengder[int(tempLast[1])]**2)) #fim for ende a
            lastVektor[int(elem[int(tempLast[1])][0])] += float((1 / 12) * tempLast[2] * (elementlengder[int(tempLast[1])]) ** 2)

            mom[int(tempLast[1])][1] += float((1/12)*tempLast[2]*(elementlengder[int(tempLast[1])]**2)) #fim for ende b
            lastVektor[int(elem[int(tempLast[1])][1])] -= float((1 / 12) * tempLast[2] * (elementlengder[int(tempLast[1])]) ** 2)

        elif tempLast[0] == 4:  # Moment fra ujevn fordelt last
            mom[int(tempLast[1])][0] -= float((1/30)*tempLast[2]*(elementlengder[int(tempLast[1])]**2)) #fim for ende a
            lastVektor[int(elem[int(tempLast[1])][0])] += float((1 / 30) * tempLast[2] * (elementlengder[int(tempLast[1])]) ** 2)

            mom[int(tempLast[1])][1] += float((1/20)*tempLast[2]*(elementlengder[int(tempLast[1])]**2)) #fim for ende b
            lastVektor[int(elem[int(tempLast[1])][1])] -= float((1 / 20) * tempLast[2] * (elementlengder[int(tempLast[1])] ** 2))
        else:
            print("Feil i lastvektor")

    return mom, lastVektor

def treghetsMoment(profil):
    if profil == 1 or profil == 3: #returnerer I for I-profil med dimensjoner fra profiler.
        return ((((profiler[int(profil)-1][2]-(2*profiler[int(profil)-1][0]))**3)*profiler[int(profil)-1][1])/12) + 2*((((profiler[int(profil)-1][0]**3)*profiler[int(profil)-1][3])/12)+(((profiler[int(profil)-1][2]-profiler[int(profil)-1][0])/2)**2)*(profiler[int(profil)-1][0]*profiler[int(profil)-1][3]))
    elif profil == 2 or profil == 4: #returnerer I for rør med dimensjoner fra profiler.
        return 2*np.pi*(profiler[int(profil)-1][0]**3)*(profiler[int(profil)-1][1])
    else:
        print("feil profil")
        return 0

def stivhetsMatrise(nelem, elem, elementlengder, npunkt):

    stivhetsmatrise = np.array([[0 for x in range(npunkt)] for y in range(npunkt)], dtype=np.float64)
    counter = 0
    for tempElem in elem:

        p00 = float((4 * tempElem[2] * treghetsMoment(tempElem[3]))/float(elementlengder[counter]))
        p01 = float((2 * tempElem[2] * treghetsMoment(tempElem[3]))/float(elementlengder[counter]))
        p10 = float((2 * tempElem[2] * treghetsMoment(tempElem[3]))/float(elementlengder[counter]))
        p11 = float((4 * tempElem[2] * treghetsMoment(tempElem[3]))/float(elementlengder[counter]))

        stivhetsmatrise[int(tempElem[0])][int(tempElem[0])] += p00
        stivhetsmatrise[int(tempElem[1])][int(tempElem[1])] += p11

        stivhetsmatrise[int(tempElem[1])][int(tempElem[0])] += p01
        stivhetsmatrise[int(tempElem[0])][int(tempElem[1])] += p10
        counter += 1

    return stivhetsmatrise

def randbetingelser(npunkt, punkt, K, b):
    Kn = K
    Bn = b
    for i in range(0, npunkt):
        if punkt[i][2] == 1:
            for j in range(0, npunkt):
                Kn[i][j] = 0    #Setter her rad lik 0 i stivhetsmatrise
                Kn[j][i] = 0    #Setter her kolonne lik 0 i stivhetsmatrise
            Bn[i] = 0           #Nuller ut tilsvarende element i lastvektoren
            Kn[i][i] = 1        #Setter diagonalelementet lik et vilkårlig tall for at matrisen skal bli inverterbar
    return Kn, Bn

def prettyPrint(matrise):
     print(DataFrame(matrise))

def endeMoment(npunkt, punkt, nelem, elem, elementlengder, rot, fim):
    basisM = [[0.0,0.0],[0.0,0.0]]
    basisRot = [0.0,0.0]
    endeM = []
    endeMprint = []

    for i in range(0, nelem):
        basisM[0][0] = float(4 * ((elem[i][2] * treghetsMoment(elem[i][3])) / elementlengder[i]))
        basisM[0][1] = float(2 * ((elem[i][2] * treghetsMoment(elem[i][3])) / elementlengder[i]))
        basisM[1][0] = float(2 * ((elem[i][2] * treghetsMoment(elem[i][3])) / elementlengder[i]))
        basisM[1][1] = float(4 * ((elem[i][2] * treghetsMoment(elem[i][3])) / elementlengder[i]))

        basisRot[0] = rot[int(elem[i][0])]
        basisRot[1] = rot[int(elem[i][1])]

        basisTemp = np.matmul(basisM, basisRot)

        basisTemp[0] += (fim[i][0])
        basisTemp[1] += (fim[i][1])

        basisTemp[0] = round(basisTemp[0], 3)
        basisTemp[1] = round(basisTemp[1], 3)

        endeM.append(basisTemp)

        basisTemp[1] = (-1) * round(basisTemp[1], 3)

        endeMprint.append((basisTemp))

    return endeM, endeMprint

def kritiskBelastedBjelke(endeMoment, midtPunktsLaster, element):

    # Finner den mest kritisk belastede bjelken med hensyn på E-modul og profil. Returnerer forholdet mellom de og flytspenningen.

    flytespenning = 355*(10**6) # Pascal
    currentMax = 0
    currentEle = 0
    distanseFraNoytralakse = 0

    for i, mom in enumerate(endeMoment):
        if element[i][3] == 1 or element[i][3] == 3:
            distanseFraNoytralakse = profiler[int((element[i][3])-1)][2] / 2     #Tversnitt bjelke
        elif element[i][3] == 2 or element[i][3] == 4:
            distanseFraNoytralakse = profiler[(int(element[i][3])-1)][0]         #Radius rør

        for j in mom:
            motstandsmoment = treghetsMoment(element[i][3]) / distanseFraNoytralakse
            if np.abs(j/(motstandsmoment*flytespenning)) > currentMax:
                currentMax = np.abs(j/(motstandsmoment*flytespenning))
                currentEle = i
                print("#Endemoment, ", "ele: ", currentEle,"max: ", currentMax)

    for i, mom in enumerate(midtPunktsLaster):
        if element[i][3] == 1:
            distanseFraNoytralakse = profiler[0][2] / 2
        elif element[i][3] == 2:
            distanseFraNoytralakse = profiler[1][0]

        motstandsmoment = treghetsMoment(element[i][3]) / distanseFraNoytralakse
        if  np.abs(mom / (motstandsmoment * flytespenning)) > currentMax:
            currentMax = np.abs(mom / (motstandsmoment*flytespenning))
            currentEle = i
            print("#Midtpunktslast, ", "ele: ", currentEle, "max: ", currentMax)

    print("Bjelke: ", currentEle, "Max mom: ", currentMax)
    return currentEle, currentMax


def main():
    # -----Leser input-data-----
    npunkt, punkt, nelem, elem, nlast, last = lesinput()

    # -----Regner ut lengder til elementene------
    elementlengder = lengder(punkt, elem, nelem)

    # -----Setter opp lastvektor og fastinnspenningsmoment-----
    fim, b = momentOgLastvektor(npunkt, punkt, nelem, elem, nlast, last, elementlengder)

    # ------Setter opp systemstivhetsmatrisen-----
    K = stivhetsMatrise(nelem, elem, elementlengder, npunkt)

    # ------Innfører randbetingelser------
    Kn, Bn = randbetingelser(npunkt, punkt, K, b)

    # -----Løser ligningssystemet------
    rot = np.linalg.solve(Kn, Bn)

    #------Finner endemoment for hvert element-----
    endemoment, endemomentPrint = endeMoment(npunkt, punkt, nelem, elem, elementlengder, rot, fim)

    #------Finner moment midt på / under punktlast for hvert element-----
    mpLaster = midtpunktsLaster(nelem, elementlengder, nlast, last, endemoment, elem)

    #------Finner elementet med mest kritisk last-----
    kritiskBelastedBjelke(endemoment, mpLaster, elem)

    prettyPrint(endemomentPrint)
    prettyPrint(mpLaster)
    #prettyPrint(Kn)


main()
