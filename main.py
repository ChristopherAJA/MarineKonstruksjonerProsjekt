import pickle

def konstruksjonsdata():
    knutepunkter = []
    elementer = []
    tall = int(input("Hvor mange knutepunkter?: "))
    for i in range(tall):
        knutepunkter.append(knutepunkt)
        knutepunkter[i].x = input("x-verdi til knutepunkt",i,": ")
        knutepunkter[i].y = input("y-verdi til knutepunkt",i,": ")
        #fast innspent y/n?
        #har moment?
    tall = int(input("Hvor mange elementer?"))
    for i in range(tall):
        elementer.append(element)
        a = int(input("Første knutepunkt festet i element",i,": "))
        elementer[i].a = knutepunkter[a]
        b = int(input("Andre knutepunkt festet i element",i,": "))
        elementer[i].b = knutepunkter[b]
        m = int(input("Materialtype for element",i,": "))
        elementer[i].m = materiale.m #?
        g = int(input("Geometritype for element",i,": "))
        elementer[i].g = geometri.g
        #har punktlast?
