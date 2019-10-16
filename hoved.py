import numpy as np
import numpy.matlib


def lesinput():

    # Åpner inputfilen
    fid = open("input.txt", "r")

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
    elem = np.loadtxt(fid, dtype = int, max_rows = nelem)

    # Leser antall laster som virker på rammen
    nlast = int(fid.readline())

    # Leser lastdata
    # Bestem selv verdiene som er nødvendig å lese inn, samt hva verdiene som leses inn skal representere
    last = np.loadtxt(fid, dtype = float, max_rows = nlast)     # <-- Forslag til innlesing av lastdata

    # Lukker input-filen
    fid.close()

    return npunkt, punkt, nelem, elem, nlast, last

import numpy as np


def lengder(knutepunkt, element, nelem):

    elementlengder = np.matlib.zeros((nelem, 1))
    # Beregner elementlengder med Pythagoras' laeresetning
    for i in range (0, nelem):
        # OBS! Grunnet indekseringsyntaks i Python-arrays vil ikke denne funksjonen fungere naar vi bare har ett element.
        dx = knutepunkt[element[i, 0], 0] - knutepunkt[element[i, 1], 0]
        dy = knutepunkt[element[i, 0], 1] - knutepunkt[element[i, 1], 1]
        elementlengder[i] = np.sqrt(dx*dx + dy*dy)

    return elementlengder

def moment(npunkt, punkt, nelem, elem, nlast, last, elementlengder):

    mom = np.mathlib.zeros((nelem,1))

def main():
    # Rammeanalyse

    # -----Leser input-data-----
    npunkt, punkt, nelem, elem, nlast, last = lesinput()

    # -----Regner ut lengder til elementene------
    elementlengder = lengder(punkt, elem, nelem)

    # -----Fastinnspenningsmomentene------
    # Lag funksjon selv
    fim = moment(npunkt, punkt, nelem, elem, nlast, last, elementlengder)

    # -----Setter opp lastvektor-----
    # Lag funksjon selv
    b = lastvektor(fim, npunkt, punkt, nelem, elem)

    # ------Setter opp systemstivhetsmatrisen-----
    # Lag funksjon selv
    K = stivhet(nelem, elem, elementlengder, npunkt)

    # ------Innfører randbetingelser------
    # Lag funksjon selv
    Kn, Bn = bc(npunkt, punkt, K, b)

    # -----Løser ligningssystemet------
    # Lag funksjon selv
    rot = ...
    # Hint, se side for løsing av lineære systemer i Python

    #------Finner endemoment for hvert element-----
    # Lag funksjon selv
    endemoment = endeM(npunkt, punkt, nelem, elem, elementlengder, rot, fim)

    #-----Skriver ut hva rotasjonen ble i de forskjellige nodene-----
    print("Rotasjoner i de ulike punktene:")
    print(rot)

    #-----Skriver ut hva momentene ble for de forskjellige elementene-----
    print("Elementvis endemoment:")
    print(endemoment)
