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
    
    if not es_duplica(q):
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
    
    if total == 0:
        return {'A': 0.25, 'T': 0.25, 'C': 0.25, 'G': 0.25}
    
    norm_conteig = conteig.copy()                                   # Copiar es conteig original

    for idx in norm_conteig:                                        # Normalitzar cada valor
        norm_conteig[idx] /= total                                       
    return norm_conteig

# Quantificar es valors de freqüència a B bits
def quantificar(freqs, B):
    quantificat = {}
    nivells = 2 ** B

    # mirar mes formes de quantificar
    for base, valor in freqs.items():
        # q = int(valor * nivells)   # truncament
        q = round(valor * nivells)   # arrodoniment

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
    p = 0.01                # probabilitat de mutació
    q = 0.0                 # probabilitat de NO duplicació
    B = 2                   # bits de quantificació
    c = 15                  # nombre de cicles PCR
    N = 100                 # nombre d'experiments
    n = 10                  # longitud de la cadena original

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

        # x es B
        # tuning parameters: c, n, p, q, B
        # mostreig aleatori totes secquencies i tenir en compte posició de cada nucleotid
        # figures, començar document, decide 
        # girar probabilitat, mirar probabiitat d'error
'''


if __name__ == "__main__":
    # Parameter tuning ranges
    c_range = range(1, 16)                     # c: 1-15
    p_range = [0.01, 0.03, 0.05, 0.07, 0.1]         # mutation probability
    q_range = [0.0, 0.02, 0.04, 0.06, 0.08, 0.1]    # non-duplication probability
    B_range = range(1, 5)                           # quantization bits: 1-4
    
    N = 100                                         # nombre d'experiments per combinació
    n = 10                                          # longitud de sa cadena original
    
    # Guardar resultats
    results_nucleotid = {}
    results_cadena_post_q = {}
    results_cadena_pre_q = {}
    
    total_combinations = len(list(c_range)) * len(p_range) * len(q_range) * len(B_range)
    current_combination = 0
    
    print(f"Parameter tuning: {total_combinations} combinations to test\n")
    
    for c in c_range:
        for p in p_range:
            for q in q_range:
                for B in B_range:
                    current_combination += 1
                    print(f"[{current_combination}/{total_combinations}] Testing: c={c}, p={p}, q={q}, B={B}")
                    
                    cadena_encertada = 0
                    cadena_encertada_pre_q = 0
                    nucleotids_encertats = 0
                    
                    for idx_experiment in range(N):
                        cadena_original = []
                        cadena_final_pre_q = []
                        cadena_final_post_q = []

                        for idx in range(n):
                            # Nucleòtid original
                            base_original = generar_cadena(1)[0]
                            cadena_original.append(base_original)

                            # PCR a nucleòtid original
                            resultat = cicles_PCR([base_original], c, p, q)

                            # Decisió de la base després de PCR
                            contador_resultat = contador(resultat)

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
    
    print("\n" + "="*60)
    print("PARAMETER TUNING COMPLETED")
    print("="*60 + "\n")
    
    # VISUALITZAR RESULTATS

    # Error vs cicles (c)
    avg_error_by_c = {}
    for c in c_range:
        errors = [results_cadena_post_q[(c, p, q, B)] for p in p_range for q in q_range for B in B_range]
        avg_error_by_c[c] = sum(errors) / len(errors)
    
    plt.figure(figsize=(12, 4))
    
    plt.subplot(1, 3, 1)
    plt.plot(list(avg_error_by_c.keys()), list(avg_error_by_c.values()), marker='o')
    plt.xlabel('Number of PCR cycles (c)')
    plt.ylabel('Average Error Rate')
    plt.title('Error Rate vs PCR Cycles')
    plt.grid()
    
    # Error vs mutation probability (p)
    avg_error_by_p = {}
    for p in p_range:
        errors = [results_cadena_post_q[(c, p, q, B)] for c in c_range for q in q_range for B in B_range]
        avg_error_by_p[p] = sum(errors) / len(errors)
    
    plt.subplot(1, 3, 2)
    plt.plot(list(avg_error_by_p.keys()), list(avg_error_by_p.values()), marker='s', color='orange')
    plt.xlabel('Mutation Probability (p)')
    plt.ylabel('Average Error Rate')
    plt.title('Error Rate vs Mutation Probability')
    plt.grid()
    
    # Error vs quantization bits (B)
    avg_error_by_B = {}
    for B in B_range:
        errors = [results_cadena_post_q[(c, p, q, B)] for c in c_range for p in p_range for q in q_range]
        avg_error_by_B[B] = sum(errors) / len(errors)
    
    plt.subplot(1, 3, 3)
    plt.plot(list(avg_error_by_B.keys()), list(avg_error_by_B.values()), marker='^', color='green')
    plt.xlabel('Quantization Bits (B)')
    plt.ylabel('Average Error Rate')
    plt.title('Error Rate vs Quantization Bits')
    plt.grid()
    
    plt.tight_layout()
    plt.savefig('parameter_sweep_overview.png', dpi=150, bbox_inches='tight')
    print("Saved: parameter_sweep_overview.png\n")
    
    # Heatmap: Error rate for different c and p combinations (averaged over q and B)
    c_values = list(c_range)
    p_values = p_range
    heatmap_data = []
    for c in c_values:
        row = []
        for p in p_values:
            errors = [results_cadena_post_q[(c, p, q, B)] for q in q_range for B in B_range]
            row.append(sum(errors) / len(errors))
        heatmap_data.append(row)
    
    plt.figure(figsize=(10, 8))
    plt.imshow(heatmap_data, cmap='RdYlGn_r', aspect='auto', origin='lower')
    plt.colorbar(label='Error Rate')
    plt.xlabel('Mutation Probability (p)')
    plt.ylabel('Number of PCR cycles (c)')
    plt.title('Error Rate Heatmap: PCR Cycles vs Mutation Probability')
    plt.xticks(range(len(p_values)), p_values)
    plt.yticks(range(len(c_values)), c_values)
    plt.tight_layout()
    plt.savefig('heatmap_c_vs_p.png', dpi=150, bbox_inches='tight')
    print("Saved: heatmap_c_vs_p.png\n")
    
    # Heatmap: Error rate for different q and B combinations
    q_values = q_range
    B_values = list(B_range)
    heatmap_data_qB = []
    for q in q_values:
        row = []
        for B in B_values:
            errors = [results_cadena_post_q[(c, p, q, B)] for c in c_range for p in p_range]
            row.append(sum(errors) / len(errors))
        heatmap_data_qB.append(row)
    
    plt.figure(figsize=(8, 8))
    plt.imshow(heatmap_data_qB, cmap='RdYlGn_r', aspect='auto', origin='lower')
    plt.colorbar(label='Error Rate')
    plt.xlabel('Quantization Bits (B)')
    plt.ylabel('Non-duplication Probability (q)')
    plt.title('Error Rate Heatmap: q vs B')
    plt.xticks(range(len(B_values)), B_values)
    plt.yticks(range(len(q_values)), [f'{q:.2f}' for q in q_values])
    plt.tight_layout()
    plt.savefig('heatmap_q_vs_B.png', dpi=150, bbox_inches='tight')
    print("Saved: heatmap_q_vs_B.png\n")
    
    plt.show()
    
    # Save detailed results to CSV
    df_results = []
    for (c, p, q, B), error in results_cadena_post_q.items():
        df_results.append({
            'cycles': c,
            'mutation_prob': p,
            'non_dup_prob': q,
            'quantization_bits': B,
            'error_rate': error
        })
    
    df = pd.DataFrame(df_results)
    df.to_csv('pcr_parameter_sweep_results.csv', index=False)
    print("Saved: pcr_parameter_sweep_results.csv\n")
    print("Top 10 parameter combinations with lowest error rate:")
    print(df.nsmallest(10, 'error_rate')[['cycles', 'mutation_prob', 'non_dup_prob', 'quantization_bits', 'error_rate']])