# Script pour supprimer la dernière colonne et la troisième colonne de chaque ligne

# Nom de ton fichier
fichier_input = "cs_15N_3.dat"
fichier_output = "cs_15N_3.dat"

# Ouverture du fichier en mode lecture
with open(fichier_input, "r") as f:
    lignes = f.readlines()

# Traitement des lignes pour supprimer la troisième et la dernière colonne
lignes_modifiees = []
for ligne in lignes:
    # On sépare la ligne en colonnes en se basant sur les espaces
    colonnes = ligne.split()
    
    # On garde toutes les colonnes sauf la troisième et la dernière
    if len(colonnes) >= 3:  # S'il y a au moins 3 colonnes
        colonnes_sans_3_et_derniere = colonnes[:2] + colonnes[3:-1]
    else:
        colonnes_sans_3_et_derniere = colonnes  # Si moins de 3 colonnes, on garde la ligne telle quelle
    
    # On reconstruit la ligne sans la troisième et la dernière colonne
    ligne_modifiee = " ".join(colonnes_sans_3_et_derniere) + "\n"
    
    # On ajoute la ligne modifiée à la liste
    lignes_modifiees.append(ligne_modifiee)

# Écriture des lignes modifiées dans un nouveau fichier
with open(fichier_output, "w") as f:
    f.writelines(lignes_modifiees)

print(f"Le fichier a été modifié et sauvegardé sous : cs180deg_15N_2.dat")
