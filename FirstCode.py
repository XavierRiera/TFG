import random

# Generar cadena d'ADN aleatòria de longitud length
def generar_cadena(length):
    cadena = []

    for idx in range(length):                                       # Iterar per tota sa cadena 
        random_nucleotid = random.choice(['A', 'T', 'C', 'G'])      # Escollir un nucleòtid aleatori
        cadena.append(random_nucleotid)                             # Afegir-lo a sa cadena
    return cadena


# Fer mutació d'un nucleòtid amb probabilitat p (canal W)
def mutar(nucleotid, p):
    rand = random.random()                                          # Nombre aleatori entre 0 i 1

    if rand > p:                                                    # Si rand > p          
        return nucleotid                                            # NO mutació(prob -> 1 - p)
    
    else:                                                           # SÍ mutació (prob -> p)                      
        possibles = ['A', 'T', 'C', 'G']
        possibles.remove(nucleotid)                                 # Llevar el nucleòtid actual
        return random.choice(possibles)                             # Mutació uniforme entre es altres 3 (prob -> p/3)


# Realitzar un cicle de PCR duplicant tota sa cadena
def PCR_transition(cadena, p):

    cadena_post_PCR = []                                            # Nova cadena després de PCR


    for copia in range(2):                                          # Es fa una còpia sencera → per cada cadena surten dues cadenes

        for base in cadena:
            base_nova = mutar(base, p)                              # Mutar es nucleòtid amb probabilitat p
            cadena_post_PCR.append(base_nova)                       # Afegir es nucleòtid mutat a sa nova cadena

    return cadena_post_PCR


# Fer c cicles de PCR
def cicles_PCR(cadena, c, p):

    for idx in range(c):                                            # Fer c cicles de PCR 
        cadena = PCR_transition(cadena, p)                          # Actualitzar sa cadena a sa nova cadena després de PCR
    return cadena


# Comptar quants A, T, C, G hi ha en una llista de nucleòtids
def contador(nucleotids):

    conteig = {'A': 0, 'T': 0, 'C': 0, 'G': 0}                      # Inicialitzar conteig

    for idx in nucleotids:                                          # Iterar per tota sa llista
        conteig[idx] += 1                                           # Incrementar es conteig corresponent
    return conteig


##########################################################################################################################################################

if __name__ == "__main__":
    
    cadena = generar_cadena(1)                                      # cadena inicial de longitud 1
    print("Cadena inicial:", cadena)
    print("Conteig inicial:", contador(cadena))

    resultat = cicles_PCR(cadena, c=5, p=0.05)                      # 5 cicles, p = 0.05
    print("\nLongitud final:", len(resultat))                       # Longitud final de sa cadena
    print("Conteig final:", contador(resultat))
