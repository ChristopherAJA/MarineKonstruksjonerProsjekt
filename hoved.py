import numpy as np
import numpy.matlib
from pandas import *

#Bruker studentnr til "Ole Bendik" = xxxx80 -> A = 8 && B = 0

profiler = [
[0.0135, 0.0086, 0.4, 0.18],
[0.1, 0.02]
] # I-profil (flens, steg, høyde, bredde).   Rør (radius, tykkelse)


def lesinput():

    # Åpner inputfilen
    fid = open("input_test.txt", "r")

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
    eleme = np.loadtxt(fid, dtype = int, max_rows = nelem)

    # Leser antall laster som virker på rammen
    nlast = int(fid.readline())

    # Leser antall laster som ikke er ujevnt fordelt
    last = np.loadtxt(fid, dtype = float, max_rows = nlast)
    lastvec = []
    # leser inn lastene og kategoriserer
    for elem in last:
        if elem[0]==1:
            lastvec.append([1,elem[1],elem[2],elem[3],elem[4]])
        elif elem[0]==2:
            lastvec.append([2,elem[1],elem[2]])
        elif elem[0]==3:
            lastvec.append([3,elem[1],elem[2]])

    # leser inn ujevnt fordelte laster
    nujamtlast = int(fid.readline())

    if nujamtlast > 0:
        ujevntfordeltlast = np.loadtxt(fid, dtype = float, max_rows = nujamtlast)
        if nujamtlast == 1:
            q_start = ujevntfordeltlast[0]
            q_slutt = ujevntfordeltlast[1]
            q_current = q_start
            q_old = q_start
            for i in range(2,ujevntfordeltlast.size):
                lastvec.append([3,ujevntfordeltlast[i],q_current])
                q_old = q_current
                q_current +=((q_slutt-q_start)/(ujevntfordeltlast.size-2))
                lastvec.append([4,ujevntfordeltlast[i],q_current-q_old])

    # Lukker input-filen
    fid.close()
    return npunkt, punkt, nelem, eleme, nlast, lastvec

def lengder(knutepunkt, element, nelem):
    elementlengder = np.matlib.zeros((nelem, 1))
    # Beregner elementlengder med Pythagoras' laeresetning
    for i in range (0, nelem):
        # OBS! Grunnet indekseringsyntaks i Python-arrays vil ikke denne funksjonen fungere naar vi bare har ett element.
        dx = knutepunkt[element[i, 0], 0] - knutepunkt[element[i, 1], 0]
        dy = knutepunkt[element[i, 0], 1] - knutepunkt[element[i, 1], 1]
        elementlengder[i] = np.sqrt(dx*dx + dy*dy)
    return elementlengder

def midtpunktsLaster(nelem, elementlengder, nlast, last, endeM):
    midtLaster =  [0.0]*nelem
    for i in range(0, len(last)):
        if last[i][0] == 1: # Punktlast
            midtLaster[int(last[i][1])] += float(((last[i][3] * \
                                        np.cos(last[i][4]) * \
                                        last[i][2] * \
                                        (elementlengder[int(last[i][1])] - last[i][2]))/(elementlengder[int(last[i][1])])) + \
                                        (((endeM[i][0] * last[i][2] * \
                                        (elementlengder[int(last[i][1])])) / elementlengder[int(last[i][1])]) - \
                                        (endeM[i][1] * last[i][2]) / elementlengder[int(last[i][1])]))

        elif last[i][0] == 3:  # Jevnlast
            midtLaster[int(last[i][1])] += float(((last[i][2] * (elementlengder[int(last[i][1])])**2)/8) + \
                                                 (endeM[i][0] * 0.5) - (endeM[i][1] * 0.5))

        elif last[i][0] == 4:  # Ujevnlast
            midtLaster[int(last[i][1])] += float(((last[i][2] * (elementlengder[int(last[i][1])])**2)/16) + \
                                                 (endeM[i][0] * 0.5) - (endeM[i][1] * 0.5))
    return midtLaster


def momentOgLastvektor(npunkt, punkt, nelem, elem, nlast, last, elementlengder):
    mom = [0.0]*npunkt #liste over fim med knutepunktnummer som indeks
    lastVektor = [0.0]*npunkt
    for tempLast in last:

        # legger til negativt fastinnspenningsmoment i lastvektoren i tillegg

        if tempLast[0] == 1:    # Moment fra punktlast
            mom[int(elem[int(tempLast[1])][0])] -= float(((tempLast[3]*np.cos(tempLast[4])) * (tempLast[2]*(elementlengder[int(tempLast[1])]-tempLast[2])**2))/((elementlengder[int(tempLast[1])])**2))                   # fim for ende a
            lastVektor[int(elem[int(tempLast[1])][0])] += float(((tempLast[3] * np.cos(tempLast[4])) * (tempLast[2] * (elementlengder[int(tempLast[1])] - tempLast[2]) ** 2)) / (elementlengder[int(tempLast[1])]) ** 2)  # fim for ende a

            mom[int(elem[int(tempLast[1])][1])] += float(((tempLast[3]*np.cos(tempLast[4]))*((tempLast[2]**2)*(elementlengder[int(tempLast[1])]-tempLast[2]))) /(elementlengder[int(tempLast[1])])**2)                      # fim for ende b
            lastVektor[int(elem[int(tempLast[1])][1])] -= float(((tempLast[3] * np.cos(tempLast[4])) * ((tempLast[2] ** 2) * (elementlengder[int(tempLast[1])] - tempLast[2]))) / (elementlengder[int(tempLast[1])]) ** 2)  # fim for ende b

        elif tempLast[0] == 2:  # Konsentrerte momenter i knutepunkt
            lastVektor[int(tempLast[1])] += float(tempLast[2])

        elif tempLast[0] == 3:  # Moment fra jevn fordelt last
            mom[int(elem[int(tempLast[1])][0])] -= float((1/12)*tempLast[2]*(elementlengder[int(tempLast[1])]**2)) #fim for ende a
            lastVektor[int(elem[int(tempLast[1])][0])] += float((1 / 12) * tempLast[2] * (elementlengder[int(tempLast[1])]) ** 2)

            mom[int(elem[int(tempLast[1])][1])] += float((1/12)*tempLast[2]*(elementlengder[int(tempLast[1])]**2)) #fim for ende b
            lastVektor[int(elem[int(tempLast[1])][1])] -= float((1 / 12) * tempLast[2] * (elementlengder[int(tempLast[1])]) ** 2)

        elif tempLast[0] == 4:  # Moment fra ujevn fordelt last
            mom[int(elem[int(tempLast[1])][0])] -= float((1/30)*tempLast[2]*(elementlengder[int(tempLast[1])]**2)) #fim for ende a
            lastVektor[int(elem[int(tempLast[1])][0])] += float((1 / 30) * tempLast[2] * (elementlengder[int(tempLast[1])]) ** 2)

            mom[int(elem[int(tempLast[1])][1])] += float((1/20)*tempLast[2]*(elementlengder[int(tempLast[1])]**2)) #fim for ende b
            lastVektor[int(elem[int(tempLast[1])][1])] -= float((1 / 20) * tempLast[2] * (elementlengder[int(tempLast[1])] ** 2))
        else:
            print("Feil i lastvektor")

    return mom, lastVektor

def treghetsMoment(profil):
    if profil == 1: #returnerer I for I-profil med dimensjoner fra profiler.
        return ((profiler[0][2]**3)*profiler[0][1])/12 + 2*((((profiler[0][0]**3)*profiler[0][3])/12)+(((profiler[0][2]-profiler[0][0])/2)**2)*(profiler[0][0]*profiler[0][3]))
    elif profil == 2: #returnerer I for rør med dimensjoner fra profiler.
        return 2*np.pi*(profiler[1][0]**3)*(profiler[1][1])
    else:
        print("feil profil")
        return 0

def stivhetsMatrise(nelem, elem, elementlengder, npunkt):
    stivhetsmatrise = np.array([[0 for x in range(npunkt)] for y in range(npunkt)], dtype=np.float64)
    counter = 0
    for tempElem in elem:

        p00 = float((4 * tempElem[2] * treghetsMoment(tempElem[3]))/float(elementlengder[counter])) #E-modul må ganges med 10^6 for å få pascal
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
                Kn[i][j] = 0
                Kn[j][i] = 0    #setter her rad og kolonne lik 0 i stivhetsmatrise
            Bn[i] = 0           #Nuller ut tilsvarende element i lastvektoren
            Kn[i][i] = 1337     #Setter diagonalelementet lik et vilkårlig tall for at matrisen skal bli inverterbar
    return Kn, Bn

def prettyPrint(matrise):
     print(DataFrame(matrise))

def endeMoment(npunkt, punkt, nelem, elem, elementlengder, rot, fim):
    basisM = [[0.0,0.0],[0.0,0.0]]
    basisRot = [0.0,0.0]
    endeM = []

    for i in range(0, nelem):
        basisM[0][0] = float(4 * ((elem[i][2] * treghetsMoment(elem[i][3])) / elementlengder[i]))
        basisM[0][1] = float(2 * ((elem[i][2] * treghetsMoment(elem[i][3])) / elementlengder[i]))
        basisM[1][0] = float(2 * ((elem[i][2] * treghetsMoment(elem[i][3])) / elementlengder[i]))
        basisM[1][1] = float(4 * ((elem[i][2] * treghetsMoment(elem[i][3])) / elementlengder[i]))

        basisRot[0] = rot[elem[i][0]]
        basisRot[1] = rot[elem[i][1]]

        basisTemp = np.matmul(basisM, basisRot)

        basisTemp[0] += fim[elem[i][0]]
        basisTemp[1] += fim[elem[i][1]]

        endeM.append(basisTemp)

    return endeM

def kritiskBelastedBjelke(endeMoment, midtPunktsLaster, element):
    flytespenning = 355*(10**6) # Pascal
    currentMax = 0
    currentEle = 0
    distanseFraNoytralakse = 0
    motstandsmoment = 0

    for i, mom in enumerate(endeMoment):
        if element[i][3] == 1:
            distanseFraNoytralakse = profiler[0][2] / 2
        elif element[i][3] == 2:
            distanseFraNoytralakse = profiler[1][0]

        for j in mom:
            motstandsmoment = treghetsMoment(element[i][3]) / distanseFraNoytralakse

            if np.abs(j/(motstandsmoment*flytespenning)) > currentMax:
                currentMax = j/(motstandsmoment*flytespenning)
                currentEle = i

    for i, mom in enumerate(midtPunktsLaster):
        if element[i][3] == 1:
            distanseFraNoytralakse = profiler[0][2] / 2
        elif element[i][3] == 2:
            distanseFraNoytralakse = profiler[1][0]

        motstandsmoment = treghetsMoment(element[i][3]) / distanseFraNoytralakse

        if np.abs(i / (motstandsmoment * flytespenning)) > currentMax:
            currentMax = i / (motstandsmoment*flytespenning)
            currentEle = i

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
    endemoment = endeMoment(npunkt, punkt, nelem, elem, elementlengder, rot, fim)

    #------Finner moment midt på / under punktlast for hvert element-----
    mpLaster = midtpunktsLaster(nelem, elementlengder, nlast, last, endemoment)

    #------Finner elementet med mest kritisk last-----
    kritiskBelastedBjelke(endemoment, mpLaster, elem)

    prettyPrint(rot)

main()
