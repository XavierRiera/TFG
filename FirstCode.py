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

'''
# Realitzar un cicle de PCR duplicant tota sa cadena
def PCR_transition(cadena, p, q, original_muta = False):

    cadena_post_PCR = []                                            # Nova cadena després de PCR

    if original_muta:                                               # Si es nucleòtid original també muta --> demanar: ara està posat que, si muta, s'orignal també muta però a lo mateix. se podria posar que muta però a una base diferent
        for base in cadena:
            if es_duplica(q):
                base_nova = mutar(base, p)                          # Mutar es nucleòtid amb probabilitat p
                cadena_post_PCR.append(base_nova)                   # Afegir es nucleòtid mutat a sa nova cadena
                cadena_post_PCR.append(base_nova)                   # Afegir es nucleòtid mutat a sa nova cadena (duplicat)



    elif not original_muta:                                         # Si es nucleòtid original no muta
        for base in cadena:
            if es_duplica(q):                                       # Comprovar si es nucleòtid es duplica amb probabilitat q
                base_nova = mutar(base, p)                          # Mutar es nucleòtid amb probabilitat p
                cadena_post_PCR.append(base_nova)                   # Afegir es nucleòtid mutat a sa nova cadena
                cadena_post_PCR.append(base)                        # Afegir es nucleòtid original a sa nova cadena


    return cadena_post_PCR                                          # Retornar sa nova cadena després de PCR
'''

# Cicle de PCR per UN SOL nucleòtid 
def PCR_nucleotid(base, p, q, original_muta = False):
    resultat = []

    if es_duplica(q):

        if original_muta:
            base_nova = mutar(base, p)
            resultat.append(base_nova)
            resultat.append(base_nova)

        else:
            base_nova = mutar(base, p)
            resultat.append(base_nova)
            resultat.append(base)

    return resultat

def PCR_cadena(cadena, p, q, original_muta = False):
    resultat = []

    for base in cadena:
        resultat.extend(PCR_nucleotid(base, p, q, original_muta))

    return resultat

# Fer c cicles de PCR
def cicles_PCR(cadena, c, p, q, original_muta = False):

    for idx in range(c):                                            # Fer c cicles de PCR 
        cadena = PCR_cadena(cadena, p, q, original_muta)            # Actualitzar sa cadena a sa nova cadena després de PCR
    return cadena


# Comptar quants A, T, C, G hi ha en una llista de nucleòtids
def contador(nucleotids):

    conteig = {'A': 0, 'T': 0, 'C': 0, 'G': 0}                      # Inicialitzar conteig

    for idx in nucleotids:                                          # Iterar per tota sa llista
        conteig[idx] += 1                                           # Incrementar es conteig des nucleòtid corresponent
    return conteig

# Normalitzar es conteig a freqüències
def normalitzar_conteig(conteig):

    total = sum(conteig.values())                                   # Num total de nucleòtids
    norm_conteig = conteig.copy()                                   # Copiar es conteig original

    for idx in norm_conteig:                                        # Normalitzar cada valor
        norm_conteig[idx] /= total                                       
    return norm_conteig

# Quantificar es valors de freqüència a B bits
def quantificar(freqs, B):
    quantificat = {}
    nivells = 2 ** B

    for base, valor in freqs.items():
        q = int(valor * nivells)   # truncament

        # cas límit valor = 1
        if q == nivells:
            q = nivells - 1

        quantificat[base] = q / (nivells)  # normalitzar entre [0, 1]

    return quantificat

# Decidir quina base té la freqüència més alta
def decidir_base(frequencies_quantificades):
    return max(frequencies_quantificades, key=frequencies_quantificades.get)        # key=frequencies_quantificades.get retorna es valor associat a cada clau (base) per fer sa comparació


##########################################################################################################################################################
'''
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
        p = 0.5                    # probabilitat de mutació
        q = 0.00                    # probabilitat de NO duplicació
        B = 2                       # bits de quantificació
        N = 10000                     # nombre d'experiments
        n = 1                       # longitud de sa cadena inicial
        c = 1                       # nombre de cicles PCR

        valors_c = range(1, 2)     # de 1 a 30 cicles PCR

        
            

        for c in valors_c:

            for idx in range(N):
                cadena = generar_cadena(n)
                print("Cadena inicialASDASSSSSSSSSSSSSSSSSSSSSS:", cadena)
                base_original = decidir_base(contador(cadena))
                #print("Base original:", base_original)

                    

                resultat = cicles_PCR(cadena, c, p, q)
                frequencies = normalitzar_conteig(contador(resultat))
                frequencies_quantificades = quantificar(frequencies, B)
                base_final = decidir_base(frequencies_quantificades)

                encerts = 0
                total = 0

                
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
'''


if __name__ == "__main__":
    p = 0.01            # probabilitat de mutació
    q = 0.0            # probabilitat de NO duplicació
    B = 1              # bits de quantificació
    c = 20              # nombre de cicles PCR
    N = 1000          # nombre d'experiments
    n = 1             # longitud de la cadena original

    cadena_encertada = 0
    cadena_encertada_pre_q = 0
    nucleotids_encertats = 0
    for idx_experiment in range(N):
        if idx_experiment < 2:
            print("EXPERIMENT ", idx_experiment + 1, "\n===================")

        cadena_original = []
        cadena_final_pre_q = []
        cadena_final_post_q = []

        for idx in range(n):
            # Nucleòtid original
            base_original = generar_cadena(1)[0]                                # Generar un nucleòtid aleatori
            cadena_original.append(base_original)                               # Afegir es nucleòtid a sa cadena original  

            # PCR a nucleòtid original
            resultat = cicles_PCR([base_original], c, p, q)                     # Fer c cicles de PCR

            # Decisió de la base després de PCR
            contador_resultat = contador(resultat)                              # conteig    

            # Conteig normalitzat --> entre [0, 1]
            conteig_normalitzat = normalitzar_conteig(contador_resultat)        # normalització

            # Quantificació a B bits
            conteig_quantificat = quantificar(conteig_normalitzat, B)

            # Decisió de la base final --> DESPRES DE QUANTIFICACIÓ
            base_final_post_q = decidir_base(conteig_quantificat)               # IMPORTANT CANVIAR CANVIAR mirar que fer quan hi ha empats #################################################
                             
            # Decisió de la base final --> ABANS DE QUANTIFICACIÓ
            base_final_pre_q = decidir_base(conteig_normalitzat)

            cadena_final_pre_q.append(base_final_pre_q)
            cadena_final_post_q.append(base_final_post_q)


            if idx >= n - 2 and idx_experiment < 2:                             # Printar 2 darrers nucleòtids dels primers 5 experiments
                print("NUCLEÒTID ", idx + 1)
                print("Base original:", base_original)
                print("Conteig post PCR:", contador_resultat)
                print("Conteig normalitzat:", conteig_normalitzat)
                print("Conteig quantificat:", conteig_quantificat)
                print("Base final decidida (post-q):", base_final_post_q)
                print("Base final decidida (pre-q):", base_final_pre_q)
                print("--------------------------------------------------")

        # Encerts per nucleòtid
        for idx in range(n):
            if cadena_final_post_q[idx] == cadena_original[idx]:                # Comprovar nucleòtid a nucleòtid
                nucleotids_encertats += 1

        # Encerts tota sa cadena bé
        if cadena_final_post_q == cadena_original:                              # Comprovar cadena completa
            cadena_encertada += 1
        
        if cadena_final_pre_q == cadena_original:
            cadena_encertada_pre_q += 1

        if idx_experiment < 2:                                                  # Printar cadenes completes de 5 primers experiments
            print("Cadena original:", cadena_original)
            print("Cadena final (pre-q):", cadena_final_pre_q)
            print("Cadena final (post-q):", cadena_final_post_q , "\n")
    
    probabilitat_nucleotid_correcte = nucleotids_encertats / (N * n)            # Quants de nucleòtids s'han encertat en tots es experiments / numero total de nucleòtids en tots es experiments
    print(f"Probabilitat nucleòtid correcte: {probabilitat_nucleotid_correcte:f}")

    probabilitat_tota_cadena_correcta = cadena_encertada / N
    print(f"Probabilitat cadena correcta: {probabilitat_tota_cadena_correcta:f}")

    probabilitat_tota_cadena_correcta_pre_q = cadena_encertada_pre_q / N
    print(f"Probabilitat cadena correcta (pre-q): {probabilitat_tota_cadena_correcta_pre_q:f}")