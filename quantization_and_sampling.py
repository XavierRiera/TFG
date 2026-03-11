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


# Calcular error en funció del nombre de mostres k
def error_rate_vs_k(c, p, q, n, ks, B, N):

    # Inicialitzar diccionaris per comptar errors
    errors_seq = {}             # Error de seqüència (predicció completa incorrecta)
    errors_base = {}            # Error de base (quantes bases estan malament en total)
    for k in ks:
        errors_seq[k] = 0
        errors_base[k] = 0

    # Fer N proves 
    for _ in range(N):

        # Generar cadena de longitud n
        original = []
        for _ in range(n):
            original.append(generar_cadena())

        # Poblacions de cada base després de c cicles de PCR
        poblacions = poblacio_per_cadena(original, c, p, q)

        # Provar cada valor de k
        for k in ks:
            prediccions = []            # Prediccions per aquesta prova i aquest k
            errors = 0                  # Bases encertades

            # Per cada base, decidir amb k samples
            for idx in range(len(original)):
                base_prediccio = decidir_base_mostreig_aleatori(poblacions[idx], k, B)      # Predicció de la base idx amb k mostres i B bits de quantització
                prediccions.append(base_prediccio)

                if base_prediccio != original[idx]:                                         # Si la predicció de la base idx és incorrecta, comptar un error --> error per base
                    errors += 1

            if prediccions != original:                                                     # Si la predicció és incorrecta, comptar un error --> error per seqüència
                errors_seq[k] += 1

            errors_base[k] += errors

    # Calcular error rate per seqüència i per base
    result = {}
    for k in ks:
        seq_err = errors_seq[k] / N             # Error rate seqüència --> quantes seqüències s'han decidit malament sobre es total de proves
        base_err = errors_base[k] / (N * n)     # Error rate base --> quantes bases s'han decidit malament sobre el total de bases
        result[k] = (seq_err, base_err)         # Guardar resultats en diccionari per cada k
    
    return result

#############################################################################################################################################################

if __name__ == "__main__":
    c = 15
    p = 0.05
    q = 0.00
    n = 5
    # Experiment parameters
    N = 10000
    ks = [1, 2, 5, 10, 20, 50, 100]

    # Try several quantization bit-depths (change this list as desired)
    bits_list = [1, 2, 4, 8]

    # Collect results into a DataFrame
    rows = []
    for B in bits_list:
        print(f"Running experiments for B={B} bits...")
        results = error_rate_vs_k(c, p, q, n, ks, B, N)
        for k in ks:
            seq_err, base_err = results[k]
            rows.append({'B': B, 'k': k, 'seq_err': seq_err, 'base_err': base_err})

    df = pd.DataFrame(rows)
    csv_name = "quantization_sampling_results.csv"
    df.to_csv(csv_name, index=False)
    print(f"Saved results CSV: {csv_name}")

    # Plot sequence error vs k with one curve per B
    plt.figure(figsize=(8,5))
    for B in bits_list:
        df_B = df[df['B'] == B].sort_values('k')
        plt.plot(df_B['k'], df_B['seq_err'], marker='o', label=f"B={B}")    

    plt.xlabel("Sample size k")
    plt.ylabel("Sequence error rate")
    plt.title("Sequence Error Rate vs Sample Size for different quantization bits")
    plt.grid(True)
    plt.legend(title='Quantization bits')
    plot_name = "sampling_error_vs_k_by_B.png"
    plt.savefig(plot_name, dpi=150, bbox_inches="tight")
    print(f"Saved plot: {plot_name}")
    plt.show()
