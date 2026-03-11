import random
import pandas as pd
import matplotlib.pyplot as plt


# Generar un nucleòtid aleatori
def generar_cadena():
    return random.choice(['A', 'T', 'C', 'G'])                                      # Bases d'ADN


# Comprovar si un nucleòtid es duplica amb probabilitat q
def es_duplica(q):
    if random.random() > q:        
        return True             # TRUE amb prob -> 1 - q
                                                                              
    return False                # FALSE amb prob -> q


# Fer mutació d'un nucleòtid amb probabilitat p (canal W)
def mutar(nucleotid, p):
    if random.random() > p:                                           
        return nucleotid                    # NO mutació prob -> 1 - p
    
    else:                                                                                
        possibles = ['A', 'T', 'C', 'G']
        possibles.remove(nucleotid)                                
        return random.choice(possibles)     # Mutació uniforme entre es altres 3 amb prob -> p/3


# Fer c cicles de PCR
def cicles_PCR(base_original, c, p, q):
    poblacio = [base_original]
    
    for _ in range(c):                                        
        nova_poblacio = []
        for base in poblacio:
            nova_poblacio.append(base)            # Afegir base original a nova població (sempre sobreviu) --> I NO MUTA!
            
            if es_duplica(q):      
                copia = mutar(base, p)            # Mutar es nucleòtid amb probabilitat p
                nova_poblacio.append(copia)       # Afegir es nucleòtid mutat a nova població
                      
        poblacio = nova_poblacio

    return poblacio


# Comptar quants A, T, C, G hi ha en una llista de nucleòtids
def contador(nucleotids):
    conteig = {'A': 0, 'T': 0, 'C': 0, 'G': 0}

    for idx in nucleotids:
        conteig[idx] += 1 
    return conteig


# Normalitzar es conteig a freqüències
def normalitzar_conteig(conteig):
    total = sum(conteig.values()) 

    if total == 0:  # Evitar divisió per zero que a vegades peta
        return {'A': 0.25, 'T': 0.25, 'C': 0.25, 'G': 0.25}
    
    norm_conteig = {base: valor / total for base, valor in conteig.items()}                                     
    
    return norm_conteig


# Quantificar es valors de freqüència a B bits
def quantificar(freqs, B):
    quantificat = {}
    nivells = 2 ** B

    for base, valor in freqs.items():
        q = round(valor * nivells)

        # cas límit valor = 1
        if q == nivells:
            q = nivells - 1

        quantificat[base] = q / (nivells)  # normalitzar entre [0, 1]

    return quantificat


# Decidir quina base té la freqüència més alta
def decidir_base(freqs):
    max_val = max(freqs.values())
    candidates = [base for base, val in freqs.items() if val == max_val]
    return random.choice(candidates)  # empat → tria aleatòriament


# Mostreig aleatori de nucleòtids
def decidir_base_mostreig_aleatori(poblacio, k, B):
    k = min(k, len(poblacio))
    mostra = random.sample(poblacio, k)
    conteig = contador(mostra)
    freqs = normalitzar_conteig(conteig)
    freqs_q = quantificar(freqs, B)
    return decidir_base(freqs_q)


# Fer c cicles de PCR per cada base d'una cadena i retornar la població resultant per cada base
def poblacio_per_cadena(cadena, c, p, q):
    poblacio = []
    for base in cadena:
        poblacio.append(cicles_PCR(base, c, p, q))
    return poblacio


# Generar codi
def generar_codi(n, length):
    codi = []
    for _ in range(length):               
        paraula = []
        for _ in range(n):         
            paraula.append(generar_cadena())
        codi.append(paraula)
    return codi


# Distància de Hamming entre dues paraules
def distancia_hamming(paraula1, paraula2):
    return sum(el1 != el2 for el1, el2 in zip(paraula1, paraula2))

# Predicció amb distància de Hamming: trobar sa paraula des codi més propera a sa predicció
def predicció_amb_distancia_hamming(prediccio, codi):
    minim_dist = float('inf')
    millor_paraula = None

    for idx in codi:
        dist = distancia_hamming(prediccio, idx)

        if dist < minim_dist:
            minim_dist = dist
            millor_paraula = idx

    return millor_paraula, minim_dist



#############################################################################################################################################################

if __name__ == "__main__":
    c = 15
    p = 0.05
    q = 0.00
    n = 5
    B = 8
    N = 10000
    k = 30
    length_codi = 10

    errors_directes = 0      # Predicció directa incorrecta
    errors_hamming = 0       # Hamming incorrecte
    
    codi = generar_codi(n, length_codi)         # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! DEMANAR, es codi el faig de nou a cada iteració o una sola vegada i a ver quantes endevina amb aquest codi en concret?

    for _ in range(N):
        paraula_triada = random.choice(codi)

        # PCR i predicció directa
        poblacions = poblacio_per_cadena(paraula_triada, c, p, q)
        prediccio = []
        for idx in range(n):
            base = decidir_base_mostreig_aleatori(poblacions[idx], k, B)
            prediccio.append(base)

        # Error predicció directa
        if prediccio != paraula_triada:
            errors_directes += 1

        # Hamming predicció
        paraula_hamming, dist = predicció_amb_distancia_hamming(prediccio, codi)
        if paraula_hamming != paraula_triada:
            errors_hamming += 1

    print(f"Error rate predicció directa: {errors_directes / N:.4f}")
    print(f"Error rate després de Hamming: {errors_hamming / N:.4f}")