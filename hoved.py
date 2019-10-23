import numpy as np
import numpy.matlib
from pandas import *

#Bruker studentnr til "Ole Bendik" = xxxx80 -> A = 8 && B = 0

profiler = [
[5.7, 4.1, 100, 55],
[100, 10]
] # flens, steg, høyde, bredde.    radius, tykkelse


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
    # Det anbefales at knutepunktnummerering starter på 0, slik at det samsvarerer med listeindeksering i Python
    # E-modul for materiale lagres i kolonne 3
    # Tverrsnittstype lagres i kolonne 4, I-profil = 1 og rørprofil = 2
    eleme = np.loadtxt(fid, dtype = int, max_rows = nelem)
    #print(eleme)

    # Leser antall laster som virker på rammen
    nlast = int(fid.readline())

    # Leser lastdata
    # Bestem selv verdiene som er nødvendig å lese inn, samt hva verdiene som leses inn skal representere
    last = np.loadtxt(fid, dtype = float, max_rows = nlast)
    lastvec = []     # <-- Forslag til innlesing av lastdata
    for elem in last:
        if elem[0]==1:
            lastvec.append([1,elem[1],elem[2],elem[3],elem[4]])
        elif elem[0]==2:
            lastvec.append([2,elem[1],elem[2]])
        elif elem[0]==3:
            lastvec.append([3,elem[1],elem[2]])

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
    #print(knutepunkt)
    elementlengder = np.matlib.zeros((nelem, 1))
    # Beregner elementlengder med Pythagoras' laeresetning
    for i in range (0, nelem):
        # OBS! Grunnet indekseringsyntaks i Python-arrays vil ikke denne funksjonen fungere naar vi bare har ett element.
        dx = knutepunkt[element[i, 0], 0] - knutepunkt[element[i, 1], 0]
        dy = knutepunkt[element[i, 0], 1] - knutepunkt[element[i, 1], 1]
        elementlengder[i] = np.sqrt(dx*dx + dy*dy)
        #print(elementlengder[i])
    return elementlengder

def momentOgLastvektor(npunkt, punkt, nelem, elem, nlast, last, elementlengder):
    mom = [0.0]*npunkt #liste over fim med knutepunktnummer som indeks
    lastVektor = [0.0]*npunkt
    for tempLast in last:

        if tempLast[0] == 1:    # Moment fra punktlast
            mom[int(elem[int(tempLast[1])][0])] -= float(((tempLast[3]*np.cos(tempLast[4])) *(tempLast[2]*(elementlengder[int(tempLast[1])]-tempLast[2])**2))/(elementlengder[int(tempLast[1])])**2)                                    #fim for ende a
            lastVektor[int(elem[int(tempLast[1])][0])] += float(((tempLast[3] * np.cos(tempLast[4])) * (tempLast[2] * (elementlengder[int(tempLast[1])] - tempLast[2]) ** 2)) / (elementlengder[int(tempLast[1])]) ** 2)  # fim for ende a

            mom[int(elem[int(tempLast[1])][1])] += float(((tempLast[3]*np.cos(tempLast[4]))*((tempLast[2]**2)*(elementlengder[int(tempLast[1])]-tempLast[2]))) /(elementlengder[int(tempLast[1])])**2)                                    #fim for ende b
            lastVektor[int(elem[int(tempLast[1])][1])] -= float(((tempLast[3] * np.cos(tempLast[4])) * ((tempLast[2] ** 2) * (elementlengder[int(tempLast[1])] - tempLast[2]))) / (elementlengder[int(tempLast[1])]) ** 2)  # fim for ende b

        elif tempLast[0] == 2:  # Konsentrerte momenter i knutepunkt
            lastVektor[int(tempLast[1])] += float(tempLast[2])

        elif tempLast[0] == 3:  # Moment fra jevn fordelt last
            mom[int(elem[int(tempLast[1])][0])] -= float((1/12)*tempLast[2]*(elementlengder[int(tempLast[1])]**2)) #fim for ende a
            lastVektor[int(elem[int(tempLast[1])][0])] += float((1 / 12) * tempLast[2] * (elementlengder[int(tempLast[1])] ** 2))

            mom[int(elem[int(tempLast[1])][1])] += float((1/12)*tempLast[2]*(elementlengder[int(tempLast[1])]**2)) #fim for ende b
            lastVektor[int(elem[int(tempLast[1])][1])] -= float((1 / 12) * tempLast[2] * (elementlengder[int(tempLast[1])] ** 2))

        elif tempLast[0] == 4:  # Moment fra ujevn fordelt last
            mom[int(elem[int(tempLast[1])][0])] -= float((1/30)*tempLast[2]*(elementlengder[int(tempLast[1])]**2)) #fim for ende a
            lastVektor[int(elem[int(tempLast[1])][0])] += float((1 / 30) * tempLast[2] * (elementlengder[int(tempLast[1])] ** 2))

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

        stivhetsmatrise[int(tempElem[0])][int(tempElem[1])] += p01
        stivhetsmatrise[int(tempElem[1])][int(tempElem[0])] += p10
        counter += 1

    return stivhetsmatrise

def printMatrix(matrise):
    for rad in matrise:
        for element in rad:
            print(element,end="\t")
        print()



def main():
    # Rammeanalyse

    # -----Leser input-data-----
    npunkt, punkt, nelem, elem, nlast, last = lesinput()
    # -----Regner ut lengder til elementene------
    elementlengder = lengder(punkt, elem, nelem)

    # -----Fastinnspenningsmomentene------
    # Lag funksjon selv


    fim, b = momentOgLastvektor(npunkt, punkt, nelem, elem, nlast, last, elementlengder)

    # -----Setter opp lastvektor-----
    # Lag funksjon selv
    #b = lastvektor(fim, npunkt, punkt, nelem, elem)

    # ------Setter opp systemstivhetsmatrisen-----
    # Lag funksjon selv
    K = stivhetsMatrise(nelem, elem, elementlengder, npunkt)

    # ------Innfører randbetingelser------
    # Lag funksjon selv
    #Kn, Bn = bc(npunkt, punkt, K, b)

    # -----Løser ligningssystemet------
    # Lag funksjon selv
    rot = np.linalg.solve(K,fim)
    # Hint, se side for løsing av lineære systemer i Python

    #------Finner endemoment for hvert element-----
    # Lag funksjon selv
    #endemoment = endeM(npunkt, punkt, nelem, elem, elementlengder, rot, fim)

    #-----Skriver ut hva rotasjonen ble i de forskjellige nodene-----
    #print("Rotasjoner i de ulike punktene:")
    #print(rot)

    #-----Skriver ut hva momentene ble for de forskjellige elementene-----
    #print("Elementvis endemoment:")
    #print(endemoment)
    print(DataFrame(b))

main()
