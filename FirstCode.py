import random
import pandas as pd
import matplotlib.pyplot as plt


# Generar cadena d'ADN aleatòria de longitud length
def generar_cadena(length):
    cadena = []

    for idx in range(length):                                       # Iterar per tota sa cadena 
        random_nucleotid = random.choice(['A', 'T', 'C', 'G'])      # Escollir un nucleòtid aleatori
        cadena.append(random_nucleotid)                             # Afegir-lo a sa cadena
    return cadena

# Comprovar si un nucleòtid es duplica amb probabilitat q
def es_duplica(q):
    rand = random.random()                                          # Nombre aleatori entre 0 i 1

    if rand > q:                                                    # Si rand > q          
        return True                                                 # SÍ es duplica (prob -> 1 - q)
    
    else:                                                           # NO es duplica (prob -> q)                      
        return False

# Fer mutació d'un nucleòtid amb probabilitat p (canal W)
def mutar(nucleotid, p):
    rand = random.random()                                          # Nombre aleatori entre 0 i 1

    if rand > p:                                                    # Si rand > p          
        return nucleotid                                            # NO mutació (prob -> 1 - p)
    
    else:                                                           # SÍ mutació (prob -> p)                      
        possibles = ['A', 'T', 'C', 'G']
        possibles.remove(nucleotid)                                 # Llevar el nucleòtid actual
        return random.choice(possibles)                             # Mutació uniforme entre es altres 3 (prob -> p/3)


# Realitzar un cicle de PCR duplicant tota sa cadena
############## posar que original tambe muta
def PCR_transition(cadena, p, q):

    cadena_post_PCR = []                                            # Nova cadena després de PCR

    ########for base in cadena:
        #cadena_post_PCR.append(base)                                # Afegir sa cadena original a sa nova cadena

    for base in cadena:
        if es_duplica(q):                                           # Comprovar si es nucleòtid es duplica amb probabilitat q
            base_nova = mutar(base, p)                              # Mutar es nucleòtid amb probabilitat p
            cadena_post_PCR.append(base_nova)                       # Afegir es nucleòtid mutat a sa nova cadena
            cadena_post_PCR.append(base_nova)                       # Afegir es nucleòtid mutat a sa nova cadena (duplicat)

    return cadena_post_PCR                                          # Per tant acabam amb sa cadena duplicada (original) i sa mutada


# Fer c cicles de PCR
def cicles_PCR(cadena, c, p, q):

    for idx in range(c):                                            # Fer c cicles de PCR 
        cadena = PCR_transition(cadena, p, q)                       # Actualitzar sa cadena a sa nova cadena després de PCR
    return cadena


# Comptar quants A, T, C, G hi ha en una llista de nucleòtids
def contador(nucleotids):

    conteig = {'A': 0, 'T': 0, 'C': 0, 'G': 0}                      # Inicialitzar conteig

    for idx in nucleotids:                                          # Iterar per tota sa llista
        conteig[idx] += 1                                           # Incrementar es conteig corresponent
    return conteig

# Normalitzar es conteig a freqüències
def normalitzar_conteig(conteig):

    total = sum(conteig.values())                                   # Num total de nucleòtids
    norm_conteig = conteig.copy()                                   # Copiar es conteig original

    for idx in norm_conteig:                                             # Normalitzar cada valor
        norm_conteig[idx] /= total                                       
    return norm_conteig

# 
#def quantificar(freqs, B):
#    quantificat = {}
#    nivell = 2 ** B
#
#    for idx in freqs:
#        quantificat[idx] = round(freqs[idx] * nivell)
#
#    return quantificat
#

def quantificar(freqs, B):
    quantificat = {}
    nivells = 2 ** B

    for base, valor in freqs.items():
        q = int(valor * nivells)   # truncament

        # cas límit valor = 1
        if q == nivells:
            q = nivells - 1

        quantificat[base] = q

    return quantificat

# Decidir quina base té la freqüència més alta
def decidir_base(frequencies_quantificades):
    return max(frequencies_quantificades, key=frequencies_quantificades.get)        # key=frequencies_quantificades.get retorna es valor associat a cada clau (base) per fer sa comparació


##########################################################################################################################################################
######## 2^B - 1 nromalització per quantificació
if __name__ == "__main__":
    prints = False
    if prints:
        cadena = generar_cadena(1)                                              # cadena inicial de longitud 1
        print("Cadena inicial:", cadena)
        print("Conteig inicial:", contador(cadena))

        resultat = cicles_PCR(cadena, c=20, p=0.05, q=0.01)                    # 5 cicles, p = 0.05, q = 0.01
        print("Longitud final:", len(resultat))                                 # Longitud final de sa cadena

        conteig = contador(resultat)
        print("Conteig final:", conteig)

        frequencies = normalitzar_conteig(conteig)
        print("Freqüències finals:", frequencies)

        frequencies_quantificades = quantificar(frequencies, B=1)
        print("Quantificació:", frequencies_quantificades)
    
    proves = True
    if proves:
        resultats = {}
        p = 0.01                    # probabilitat de mutació
        q = 0.00                    # error de duplicació
        B = 10                       # bits de quantificació
        N = 10                     # nombre d'experiments
        n = 2                     # longitud de sa cadena inicial

        valors_c = range(1, 2)     # de 1 a 30 cicles PCR

        for c in valors_c:
            encerts = 0

            for idx in range(N):
                cadena = generar_cadena(n)
                base_original = decidir_base(contador(cadena))

                resultat = cicles_PCR(cadena, c, p, q)
                frequencies = normalitzar_conteig(contador(resultat))
                frequencies_quantificades = quantificar(frequencies, B)
                base_final = decidir_base(frequencies)

                if base_original == base_final:
                    encerts += 1
                
                print(cadena, base_original, base_final, resultat, frequencies, frequencies_quantificades)
                # canviar cicle pcr 1 per 1, cada nucleotid de sa cadena es iid
                # precisio és precisio de que tota sa cadena estigui be
                # parameter tuning c, n, p , q , b

            precisio = encerts / N
            resultats[c] = precisio
            print(f"Cicles: {c}, Precisió: {precisio:.4f}")


        # Visualitzar resultats
        #df = pd.DataFrame(list(resultats.items()), columns=['Cicles PCR', 'Precisió'])
        #plt.figure(figsize=(10, 6))
        #plt.plot(df['Cicles PCR'], df['Precisió'], marker='o')
        #plt.title('Precisió de la Decisió de la Base vs Nombre de Cicles PCR')
        #plt.xlabel('Nombre de Cicles PCR')
        #plt.ylabel('Precisió')
        #plt.ylim(0, 1)
        #plt.grid()
        #plt.show()

    print(int(1.1))
            
