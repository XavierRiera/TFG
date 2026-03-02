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


##########################################################################################################################################################
if __name__ == "__main__":
    c_range = range(1, 16)
    p_range = [0.01, 0.03, 0.05, 0.07, 0.1]
    q_range = [0.0, 0.02, 0.04, 0.06, 0.08, 0.1]
    B_range = range(1, 5)                    
    
    N = 10000                                   
    n = 2   
    
    # Guardar resultats
    results_nucleotid = {}
    results_cadena_post_q = {}
    results_cadena_pre_q = {}
    
    total = len(list(c_range)) * len(p_range) * len(q_range) * len(B_range)
    current = 0
    
    print(f"Testing {total} combinations\n")
    
    for c in c_range:
        for p in p_range:
            for q in q_range:
                for B in B_range:
                    current += 1
                    print(f"[{current}/{total}] Testing: c={c}, p={p}, q={q}, B={B}")
                    
                    cadena_encertada = 0
                    cadena_encertada_pre_q = 0
                    nucleotids_encertats = 0
                    
                    for _ in range(N):
                        cadena_original = []
                        cadena_final_pre_q = []
                        cadena_final_post_q = []

                        for _ in range(n):
                            # Nucleòtid original
                            base_original = generar_cadena()
                            cadena_original.append(base_original)

                            # PCR a nucleòtid original
                            poblacio = cicles_PCR(base_original, c, p, q)

                            # Decisió de la base després de PCR
                            contador_resultat = contador(poblacio)

                            # Conteig normalitzat --> entre [0, 1]
                            conteig_normalitzat = normalitzar_conteig(contador_resultat)

                            # Quantificació a B bits
                            conteig_quantificat = quantificar(conteig_normalitzat, B)

                            # Decisió de la base final
                            base_final_post_q = decidir_base(conteig_quantificat)
                            base_final_pre_q = decidir_base(conteig_normalitzat)  

                            cadena_final_pre_q.append(base_final_pre_q)
                            cadena_final_post_q.append(base_final_post_q)

                        # Encerts per nucleòtid
                        for idx in range(n):
                            if cadena_final_post_q[idx] == cadena_original[idx]:
                                nucleotids_encertats += 1

                        # Encerts tota sa cadena
                        if cadena_final_post_q == cadena_original:
                            cadena_encertada += 1
                        
                        if cadena_final_pre_q == cadena_original:
                            cadena_encertada_pre_q += 1

                    # Calculate error rates (1 - accuracy)
                    error_nucleotid = 1 - (nucleotids_encertats / (N * n))
                    error_cadena_post_q = 1 - (cadena_encertada / N)
                    error_cadena_pre_q = 1 - (cadena_encertada_pre_q / N)
                    
                    # Store results
                    key = (c, p, q, B)
                    results_nucleotid[key] = error_nucleotid
                    results_cadena_post_q[key] = error_cadena_post_q
                    results_cadena_pre_q[key] = error_cadena_pre_q

    
    print("\nDone!\n")

    # --- VISUALITZAR RESULTATS ---

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Error vs cycles (c)
    avg_error_by_c = {}
    for c in c_range:
        errors = [results_cadena_post_q[(c, p, q, B)]
                  for p in p_range for q in q_range for B in B_range]
        avg_error_by_c[c] = sum(errors) / len(errors)

    axes[0].plot(list(avg_error_by_c.keys()), list(avg_error_by_c.values()), marker='o')
    axes[0].set_xlabel('Number of PCR cycles (c)')
    axes[0].set_ylabel('Average Error Rate')
    axes[0].set_title('Error Rate vs PCR Cycles')
    axes[0].set_ylim(0, 1)
    axes[0].grid()

    # Error vs mutation probability (p)
    avg_error_by_p = {}
    for p in p_range:
        errors = [results_cadena_post_q[(c, p, q, B)]
                  for c in c_range for q in q_range for B in B_range]
        avg_error_by_p[p] = sum(errors) / len(errors)

    axes[1].plot(list(avg_error_by_p.keys()), list(avg_error_by_p.values()),
                 marker='s', color='orange')
    axes[1].set_xlabel('Mutation Probability (p)')
    axes[1].set_ylabel('Average Error Rate')
    axes[1].set_title('Error Rate vs Mutation Probability')
    axes[1].set_ylim(0, 1)
    axes[1].grid()

    # Error vs quantization bits (B)
    avg_error_by_B = {}
    for B in B_range:
        errors = [results_cadena_post_q[(c, p, q, B)]
                  for c in c_range for p in p_range for q in q_range]
        avg_error_by_B[B] = sum(errors) / len(errors)

    axes[2].plot(list(avg_error_by_B.keys()), list(avg_error_by_B.values()),
                 marker='^', color='green')
    axes[2].set_xlabel('Quantization Bits (B)')
    axes[2].set_ylabel('Average Error Rate')
    axes[2].set_title('Error Rate vs Quantization Bits')
    axes[2].set_ylim(0, 1)
    axes[2].grid()

    plt.tight_layout()
    plt.savefig('parameter_sweep.png', dpi=150, bbox_inches='tight')
    print("Saved: parameter_sweep.png")

    # --- GUARDAR CSV ---
    df_results = []
    for (c, p, q, B), error in results_cadena_post_q.items():
        df_results.append({
            'cycles': c,
            'mutation_prob': p,
            'non_dup_prob': q,
            'quantization_bits': B,
            'error_rate_post_q': error,
            'error_rate_pre_q': results_cadena_pre_q[(c, p, q, B)],
        })

    df = pd.DataFrame(df_results)
    df.to_csv('pcr_sweep.csv', index=False)
    print("\nSaved: pcr_sweep.csv")
    print("\nTop 10 lowest error combinations:")
    print(df.nsmallest(10, 'error_rate_post_q')[
        ['cycles', 'mutation_prob', 'non_dup_prob', 'quantization_bits', 'error_rate_post_q']
    ].to_string(index=False))

    plt.show()