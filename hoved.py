import numpy as np
import numpy.matlib


def lesinput():

    # Åpner inputfilen
    fid = open("input2.txt", "r")

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
    ujevntfordeltlast = np.loadtxt(fid, dtype = float, max_rows = nujamtlast)
    for elem in ujevntfordeltlast:
        for i in range(2,elem.size()):



    # Lukker input-filen
    fid.close()

    return npunkt, punkt, nelem, eleme, nlast, lastvec



def lengder(knutepunkt, element, nelem):
    #print(knutepunkt)

    elementlengder = np.matlib.zeros((nelem, 1))
    # Beregner elementlengder med Pythagoras' laeresetning
    for i in range (0, nelem-1):
        # OBS! Grunnet indekseringsyntaks i Python-arrays vil ikke denne funksjonen fungere naar vi bare har ett element.
        dx = knutepunkt[element[i, 0], 0] - knutepunkt[element[i, 1], 0]
        dy = knutepunkt[element[i, 0], 1] - knutepunkt[element[i, 1], 1]
        elementlengder[i] = np.sqrt(dx*dx + dy*dy)

    return elementlengder

def moment(npunkt, punkt, nelem, elem, nlast, last, elementlengder):

    mom = []
    mom.append()



def main():
    # Rammeanalyse

    # -----Leser input-data-----
    npunkt, punkt, nelem, elem, nlast, last = lesinput()

    # -----Regner ut lengder til elementene------
    elementlengder = lengder(punkt, elem, nelem)

    # -----Fastinnspenningsmomentene------
    # Lag funksjon selv
    #fim = moment(npunkt, punkt, nelem, elem, nlast, last, elementlengder)

    # -----Setter opp lastvektor-----
    # Lag funksjon selv
    #b = lastvektor(fim, npunkt, punkt, nelem, elem)

    # ------Setter opp systemstivhetsmatrisen-----
    # Lag funksjon selv
    #K = stivhet(nelem, elem, elementlengder, npunkt)

    # ------Innfører randbetingelser------
    # Lag funksjon selv
    #Kn, Bn = bc(npunkt, punkt, K, b)

    # -----Løser ligningssystemet------
    # Lag funksjon selv
    #rot = ...
    # Hint, se side for løsing av lineære systemer i Python

    #------Finner endemoment for hvert element-----
    # Lag funksjon selv
    #endemoment = endeM(npunkt, punkt, nelem, elem, elementlengder, rot, fim)

    #-----Skriver ut hva rotasjonen ble i de forskjellige nodene-----
    print("Rotasjoner i de ulike punktene:")
    print(elem)

    #-----Skriver ut hva momentene ble for de forskjellige elementene-----
    print("Elementvis endemoment:")
    #print(endemoment)

main()
