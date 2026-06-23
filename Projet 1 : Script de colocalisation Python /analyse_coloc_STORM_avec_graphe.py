"""
=============================================================
  ANALYSE DE COLOCALISATION STORM
  Auteur   : Stage L3 Génie Bio-Informatique
  Objectif : Déterminer si deux protéines du réticulum
             endoplasmique ou autre (W1 et W2) colocalisent, en
             construisant des réseaux de Delaunay intra-canal
             et en traçant les liaisons inter-canal.
=============================================================

AVANT DE LANCER CE PROGRAMME, INSTALLER LES DÉPENDANCES UNE SEULE FOIS :
  Dans le terminal, taper :
    pip install pandas matplotlib scipy numpy pillow plotly

CE QUE PRODUIT LE PROGRAMME (3 fichiers) :
  1. _1_reference.tif         → nuage de points brut, sans liaisons
  2. _2_colocalisation.tif    → réseaux Delaunay + liaisons jaunes
  3. _graphe_sensibilite.html → graphe interactif à ouvrir dans le navigateur
"""

# ═══════════════════════════════════════════════════════════════
# SECTION 1 — IMPORTS
# Toutes les bibliothèques externes utilisées dans le programme.
# Chacune a un rôle précis :
#   - tkinter       : construire l'interface graphique (fenêtre, boutons, etc.)
#   - pandas        : lire et manipuler les fichiers CSV (tableaux de données)
#   - numpy         : calculs mathématiques sur des tableaux de coordonnées
#   - matplotlib    : dessiner et sauvegarder les images TIF
#   - scipy         : triangulation de Delaunay + recherche de voisins (KD-Tree)
#   - plotly        : générer le graphe HTML interactif
#   - os            : manipuler les chemins de fichiers
#   - threading     : lancer l'analyse en arrière-plan sans bloquer l'interface
#   - datetime      : horodater les fichiers générés
# ═══════════════════════════════════════════════════════════════

import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")   # moteur sans fenêtre → fonctionne partout, même sans écran
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.spatial import Delaunay, cKDTree
import plotly.graph_objects as go
import os, threading, datetime


# ═══════════════════════════════════════════════════════════════
# SECTION 2 — CONSTANTES DE COULEURS ET D'ÉPAISSEUR
# Toutes les couleurs de l'interface et des images sont définies
# ici en un seul endroit. Si tu veux changer une couleur ou
# l'épaisseur des traits, tu n'as qu'à modifier cette section.
#
# Les couleurs sont en format hexadécimal (#RRGGBB).
# LINE_WIDTH contrôle l'épaisseur de TOUS les traits sur les images
# (arêtes vertes W1, arêtes rouges W2, liaisons jaunes).
# Pour des traits plus épais → augmenter (ex: 0.8)
# Pour des traits plus fins  → diminuer (ex: 0.2)
# ═══════════════════════════════════════════════════════════════

# Couleurs de l'interface graphique (fenêtre tkinter)
BG       = "#0f1117"   # fond principal (noir bleuté)
BG2      = "#1a1d27"   # fond secondaire (cartes, zones de texte)
ACCENT   = "#4fc3f7"   # couleur d'accentuation (bleu clair)
ACCENT2  = "#81c784"   # accentuation secondaire (vert clair)
TXT      = "#e8eaf6"   # texte principal (blanc cassé)
TXT2     = "#90a4ae"   # texte secondaire (gris bleuté)
ENTRY_BG = "#252836"   # fond des champs de saisie
BTN_OK   = "#1565c0"   # couleur du bouton "Lancer"
BTN_HOV  = "#1976d2"   # couleur du bouton au survol de la souris

# Couleurs des éléments sur les images TIF générées
C1_COLOR   = "#00e676"   # points W1 (vert vif)
C2_COLOR   = "#ff1744"   # points W2 (rouge vif)
E1_COLOR   = "#00c853"   # arêtes du réseau Delaunay W1 (vert sombre)
E2_COLOR   = "#d50000"   # arêtes du réseau Delaunay W2 (rouge sombre)
LINK_COLOR = "#ffd600"   # liaisons inter-canal W1↔W2 (jaune)

# Épaisseur uniforme de tous les traits sur les images TIF — LIGNE 48
# Changer cette valeur unique modifie simultanément arêtes ET liaisons
LINE_WIDTH = 0.4


# ═══════════════════════════════════════════════════════════════
# SECTION 3 — CHARGEMENT DES DONNÉES
# Deux façons de charger les données :
#   A) charger_csv() : lit deux fichiers CSV séparés (W1 et W2)
#      Format attendu : colonnes "x [nm]" et "y [nm]"
#      C'est le cas de tes fichiers STORM demixés.
#   B) charger_tif() : lit une image TIF RGB
#      Canal vert (index 1) → molécules W1
#      Canal rouge (index 0) → molécules W2
#      Le programme extrait les coordonnées des pixels actifs.
# Les deux fonctions retournent exactement le même format de
# données (df1, df2), donc toute la suite du programme est
# identique quel que soit le mode d'entrée.
# ═══════════════════════════════════════════════════════════════

def charger_csv(path_w1, path_w2):
    """
    Charge deux fichiers CSV STORM (un par canal).
    Vérifie que les colonnes 'x [nm]' et 'y [nm]' existent.
    Retourne deux DataFrames pandas : (df_W1, df_W2).
    """
    df1 = pd.read_csv(path_w1)
    df2 = pd.read_csv(path_w2)
    # Vérification : si les colonnes attendues sont absentes, on arrête
    # avec un message d'erreur clair plutôt qu'un crash incompréhensible
    for df, nom in [(df1, "W1"), (df2, "W2")]:
        if "x [nm]" not in df.columns or "y [nm]" not in df.columns:
            raise ValueError(
                f"Le fichier {nom} n'a pas les colonnes 'x [nm]' / 'y [nm]'.\n"
                f"Colonnes trouvées : {list(df.columns)}"
            )
    return df1, df2


def charger_tif(path_tif):
    """
    Charge une image TIF STORM RGB.
    Extrait les coordonnées des pixels actifs de chaque canal :
      - Canal vert (img[:,:,1]) → W1
      - Canal rouge (img[:,:,0]) → W2
    'seuil' = intensité minimale pour considérer un pixel comme signal
    (filtre le bruit de fond).
    Retourne deux DataFrames dans le même format que charger_csv().
    """
    from PIL import Image as PILImage
    img = np.array(PILImage.open(path_tif))
    if img.ndim == 2:
        # Image en niveaux de gris = pas de séparation possible en deux canaux
        raise ValueError("Image TIF en niveaux de gris. Il faut une image RGB.")
    seuil = 30   # pixels en dessous de 30/255 considérés comme bruit
    ys1, xs1 = np.where(img[:, :, 1] > seuil)   # pixels verts actifs → W1
    ys2, xs2 = np.where(img[:, :, 0] > seuil)   # pixels rouges actifs → W2
    df1 = pd.DataFrame({"x [nm]": xs1.astype(float), "y [nm]": ys1.astype(float)})
    df2 = pd.DataFrame({"x [nm]": xs2.astype(float), "y [nm]": ys2.astype(float)})
    return df1, df2


# ═══════════════════════════════════════════════════════════════
# SECTION 4 — TRIANGULATION DE DELAUNAY FILTRÉE
# La triangulation de Delaunay relie tous les points d'un nuage
# en triangles, de façon à maximiser les angles minimaux
# (on évite les triangles trop "plats" et allongés).
# Pourquoi Delaunay ? Parce que le réticulum endoplasmique forme
# naturellement des réseaux de mailles triangulaires/polygonales.
# Delaunay reproduit visuellement cette structure à partir des
# coordonnées des molécules.
#
# PROBLÈME : Delaunay relie aussi des points très éloignés
# (zones vides de la cellule), ce qui crée de grandes arêtes
# sans sens biologique. On les supprime avec le filtre max_edge_nm.
#
# max_edge_nm = longueur maximale d'une arête conservée (en nm).
# Valeur recommandée : 300-800 nm selon la densité des données.
# ═══════════════════════════════════════════════════════════════

def delaunay_filtered(coords, max_edge_nm):
    """
    Calcule la triangulation de Delaunay sur un ensemble de
    coordonnées 2D (en nm), puis filtre les arêtes trop longues.

    coords      : tableau numpy (N, 2) de coordonnées x,y en nm
    max_edge_nm : longueur maximale des arêtes à conserver (nm)

    Retourne une liste de paires d'indices (i, j) représentant
    les arêtes conservées. Chaque paire signifie : "tracer un
    segment entre coords[i] et coords[j]".
    """
    if len(coords) < 4:
        # Delaunay nécessite au moins 4 points non colinéaires
        return []
    tri = Delaunay(coords)
    edges = {}   # dictionnaire pour éviter les doublons d'arêtes
    for simplex in tri.simplices:
        # Chaque simplex est un triangle défini par 3 indices de points
        for k in range(3):
            a, b = simplex[k], simplex[(k + 1) % 3]
            key = (min(a, b), max(a, b))   # clé unique indépendante de l'ordre
            if key not in edges:
                dist = np.linalg.norm(coords[a] - coords[b])   # longueur de l'arête
                if dist <= max_edge_nm:
                    # On ne garde que les arêtes dans le seuil biologique
                    edges[key] = dist
    return list(edges.keys())


# ═══════════════════════════════════════════════════════════════
# SECTION 5 — CALCUL DU SCORE DE COLOCALISATION
# C'est le cœur mathématique du programme.
# Cette fonction est appelée deux fois :
#   1. pendant l'analyse principale (pour le score final)
#   2. pendant la génération du graphe HTML (pour les courbes)
#
# FORMULE :
#   score_W1 = (nombre de molécules W1 ayant au moins 1 voisin W2
#               dans le rayon) / (total molécules W1) × 100
#   score_W2 = idem dans l'autre sens
#   score_global = √(score_W1 × score_W2)  ← moyenne géométrique
#
# POURQUOI LA MOYENNE GÉOMÉTRIQUE ?
#   Elle pénalise les cas asymétriques. Si score_W1=100% mais
#   score_W2=1%, la moyenne arithmétique donnerait 50.5% (trompeuse),
#   la moyenne géométrique donne √(100×1) = 10% (plus réaliste).
#
# COMMENT ÇA MARCHE TECHNIQUEMENT (KD-Tree) ?
#   Un KD-Tree est une structure de données qui permet de chercher
#   les voisins proches très rapidement (O(log n) au lieu de O(n²)).
#   Sans ça, comparer 5780 points W1 × 53909 points W2 prendrait
#   des heures. Avec KD-Tree, c'est quelques secondes.
#
# LIMITE IMPORTANTE :
#   W2 est très dense (53 909 molécules dans la même zone que
#   5 780 molécules W1). À cette densité, presque n'importe quel
#   point W1 trouvera un voisin W2 dans 250 nm, même par hasard.
#   Le score mesure donc la CO-PRÉSENCE SPATIALE, pas forcément
#   une interaction biologique réelle.
# ═══════════════════════════════════════════════════════════════

def calculer_score(coords1, coords2, rayon_nm):
    """
    Calcule les trois scores de colocalisation pour un rayon donné.
    Retourne : (score_W1, score_W2, score_global, masque_W1, masque_W2)
      - score_W1/W2   : pourcentage de molécules colocalisées (0-100)
      - score_global  : moyenne géométrique des deux scores
      - masque_W1/W2  : tableau booléen (True = molécule colocalisée)
    """
    # Construction des KD-Trees pour recherche rapide de voisins
    tree2 = cKDTree(coords2)   # arbre sur les points W2
    tree1 = cKDTree(coords1)   # arbre sur les points W1

    # Pour chaque point W1, liste des indices W2 dans le rayon
    # query_ball_point retourne une liste de listes d'indices
    m1 = np.array([len(v) > 0 for v in tree2.query_ball_point(coords1, r=rayon_nm)])
    m2 = np.array([len(v) > 0 for v in tree1.query_ball_point(coords2, r=rayon_nm)])

    # Calcul des pourcentages
    sw1 = 100 * m1.sum() / len(coords1)
    sw2 = 100 * m2.sum() / len(coords2)
    sg  = np.sqrt(sw1 * sw2)   # moyenne géométrique
    return sw1, sw2, sg, m1, m2


# ═══════════════════════════════════════════════════════════════
# SECTION 6 — ANALYSE COMPLÈTE
# Fonction principale qui orchestre toutes les étapes :
#   Étape 1 : Delaunay W1 (sur tous les points W1)
#   Étape 2 : Delaunay W2 (sur un sous-échantillon si trop dense)
#   Étape 3 : Calcul des liaisons inter-canal
#   Étape 4 : Calcul du score global
#
# POURQUOI SOUS-ÉCHANTILLONNER W2 POUR DELAUNAY ?
#   W2 a 53 909 molécules. Calculer Delaunay sur autant de points
#   prendrait trop longtemps et produirait des millions d'arêtes
#   illisibles. On tire 15 000 points au hasard (seed=42 pour
#   que ce soit toujours les mêmes points → résultat reproductible).
#
# CORRECTION BUG LIAISONS JAUNES :
#   Dans la version précédente, les liaisons étaient calculées sur
#   TOUS les points W2 (53 909) mais seuls 15 000 étaient affichés.
#   Résultat : des traits jaunes pointaient dans le vide.
#   Correction : les liaisons sont maintenant calculées ET affichées
#   sur le MÊME sous-ensemble de points → cohérence visuelle garantie.
#   Note : le SCORE lui reste calculé sur TOUS les points (plus précis).
# ═══════════════════════════════════════════════════════════════

def analyser_colocalisation(df1, df2, rayon_nm, max_edge_nm, log_fn=None):
    """
    Lance l'analyse complète de colocalisation.
    log_fn : fonction de log optionnelle (pour afficher dans l'interface)
    Retourne un dictionnaire contenant toutes les données nécessaires
    pour générer les images et le graphe.
    """
    def _log(msg):
        if log_fn: log_fn(msg)

    # Conversion des DataFrames en tableaux numpy (plus rapide pour les calculs)
    coords1 = df1[["x [nm]", "y [nm]"]].values   # shape (N1, 2)
    coords2 = df2[["x [nm]", "y [nm]"]].values   # shape (N2, 2)

    # ── Étape 1 : Delaunay W1 ────────────────────────────────
    # On triangule tous les points W1 entre eux.
    # edges1 = liste de paires (i,j) d'indices dans coords1
    _log("Construction du réseau Delaunay W1...")
    edges1 = delaunay_filtered(coords1, max_edge_nm)
    _log(f"  Arêtes W1 : {len(edges1):,}")

    # ── Étape 2 : Delaunay W2 ────────────────────────────────
    # W2 est souvent trop dense pour un Delaunay complet.
    # On sous-échantillonne à MAX_DEL points pour la performance.
    # seed=42 garantit que le même sous-ensemble est choisi à chaque
    # exécution avec les mêmes données → résultats reproductibles.
    MAX_DEL = 15_000
    if len(coords2) > MAX_DEL:
        _log(f"W2 dense ({len(coords2):,} loc.) — sous-échantillon {MAX_DEL:,} pour Delaunay...")
        np.random.seed(42)
        idx_sub = np.random.choice(len(coords2), MAX_DEL, replace=False)
        coords2_del = coords2[idx_sub]   # sous-ensemble pour Delaunay uniquement
    else:
        coords2_del = coords2
    _log("Construction du réseau Delaunay W2...")
    edges2 = delaunay_filtered(coords2_del, max_edge_nm)
    _log(f"  Arêtes W2 : {len(edges2):,}")

    # ── Étape 3 : Liaisons inter-canal ───────────────────────
    # CORRECTION BUG : on définit d'abord les sous-ensembles de points
    # qui seront AFFICHÉS sur l'image, puis on calcule les liaisons
    # UNIQUEMENT entre ces points affichés.
    # Ainsi, chaque trait jaune part d'un point vert visible ET
    # pointe vers un point rouge visible → plus de traits dans le vide.
    _log("Calcul des liaisons inter-canal...")
    MAX_PLOT = 15_000

    np.random.seed(42)
    # Sous-ensemble de W2 affiché sur l'image
    idx2_plot = np.random.choice(len(coords2), min(len(coords2), MAX_PLOT), replace=False)
    coords2_plot = coords2[idx2_plot]

    # Sous-ensemble de W1 affiché sur l'image
    idx1_plot = np.random.choice(len(coords1), min(len(coords1), MAX_PLOT), replace=False)
    coords1_plot = coords1[idx1_plot]

    # KD-Tree construit UNIQUEMENT sur les points W2 affichés
    # → les liaisons ne peuvent pointer que vers des points visibles
    tree2_plot = cKDTree(coords2_plot)
    idx_voisins = tree2_plot.query_ball_point(coords1_plot, r=rayon_nm)

    # Construction de la liste des paires à tracer
    paires = []
    for i, voisins in enumerate(idx_voisins):
        if not voisins:
            continue   # ce point W1 n'a aucun voisin W2 dans le rayon
        # Parmi tous les voisins W2 dans le rayon, on choisit le plus proche
        dists = np.linalg.norm(coords2_plot[voisins] - coords1_plot[i], axis=1)
        j = voisins[np.argmin(dists)]
        # On stocke les coordonnées des deux extrémités du trait
        paires.append([
            coords1_plot[i, 0], coords1_plot[i, 1],   # départ : point W1
            coords2_plot[j, 0], coords2_plot[j, 1]    # arrivée : point W2
        ])
    paires = np.array(paires) if paires else np.empty((0, 4))
    _log(f"  Liaisons tracées : {len(paires):,}")

    # ── Calcul de la distance moyenne de colocalisation ───────
    # Formule : somme des longueurs des traits jaunes / nombre de traits jaunes
    # Chaque trait jaune = une paire [x1, y1, x2, y2]
    # Longueur du trait = sqrt((x2-x1)² + (y2-y1)²)
    if len(paires) > 0:
        longueurs = np.sqrt(
            (paires[:, 2] - paires[:, 0])**2 +
            (paires[:, 3] - paires[:, 1])**2
        )
        dist_moyenne = float(longueurs.mean())
    else:
        dist_moyenne = 0.0
    _log(f"  Distance moyenne de colocalisation : {dist_moyenne:.1f} nm")

    # ── Étape 4 : Score (sur TOUTES les localisations) ───────────
    # On utilise TOUS les points (pas le sous-échantillon) pour
    # que le score soit le plus précis et représentatif possible.
    sw1, sw2, sg, m1, m2 = calculer_score(coords1, coords2, rayon_nm)

    # Verdict selon le score global
    verdict = ("COLOCALISATION SIGNIFICATIVE DETECTEE" if sg >= 20
               else "COLOCALISATION PARTIELLE / FAIBLE" if sg >= 5
               else "PAS DE COLOCALISATION DETECTEE")

    # Message de limite de densité (conditionnel — ratio > 3)
    # Seulement affiché si un canal est plus de 3x plus dense que l autre.
    # Si les deux canaux sont équilibrés : msg_limite = chaine vide.
    n1, n2 = len(coords1), len(coords2)
    ratio = n2 / n1 if n1 > 0 else 1.0
    SEUIL_RATIO = 3.0

    if ratio > SEUIL_RATIO:
        msg_limite = (
            f"LIMITE : W2 est tres dense ({n2:,} loc. vs {n1:,} loc., ratio={ratio:.1f}x).\n"
            f"Le score W2 peut etre surestime par effet de densite.\n"
            f"Validation rigoureuse : test de permutation."
        )
    elif (1.0 / ratio) > SEUIL_RATIO:
        ratio_inv = n1 / n2
        msg_limite = (
            f"LIMITE : W1 est tres dense ({n1:,} loc. vs {n2:,} loc., ratio={ratio_inv:.1f}x).\n"
            f"Le score W1 peut etre surestime par effet de densite.\n"
            f"Validation rigoureuse : test de permutation."
        )
    else:
        msg_limite = ""   # canaux equilibres, aucune limite de densite a signaler

    # On retourne toutes les données dans un dictionnaire
    # pour les passer aux fonctions de génération d images
    return {
        "coords1": coords1,           # tous les points W1
        "coords2": coords2,           # tous les points W2
        "coords1_plot": coords1_plot, # points W1 affichés (sous-ensemble)
        "coords2_plot": coords2_plot, # points W2 affichés (sous-ensemble)
        "coords2_del":  coords2_del,  # points W2 pour Delaunay (sous-ensemble)
        "edges1": edges1,             # arêtes réseau Delaunay W1
        "edges2": edges2,             # arêtes réseau Delaunay W2
        "paires": paires,             # liaisons inter-canal à tracer
        "n_w1": len(coords1), "n_w2": len(coords2),
        "n_w1_coloc": int(m1.sum()), "n_w2_coloc": int(m2.sum()),
        "score_w1": sw1, "score_w2": sw2, "score_global": sg,
        "verdict": verdict,
        "rayon_nm": rayon_nm, "max_edge_nm": max_edge_nm,
        "dist_moyenne": dist_moyenne,
        "msg_limite": msg_limite,
    }


# ═══════════════════════════════════════════════════════════════
# SECTION 7 — GRAPHE HTML INTERACTIF
# Génère un fichier .html ouvrable dans n'importe quel navigateur.
# Le graphe contient deux panneaux côte à côte :
#
# Panneau gauche — Score vs Rayon :
#   Montre comment les scores W1, W2 et global évoluent quand on
#   fait varier le rayon de 10 à 635 nm. La ligne bleue verticale
#   marque le rayon actuellement choisi. La ligne orange horizontale
#   marque le seuil de 20% (colocalisation significative).
#   Utilité : voir si le score est stable ou très sensible au rayon.
#
# Panneau droit — Arêtes vs Seuil Delaunay :
#   Montre combien d'arêtes sont conservées dans les réseaux W1 et W2
#   selon le seuil d'arête max. La ligne bleue verticale marque le
#   seuil actuellement choisi.
#   Utilité : choisir un seuil qui donne des réseaux ni trop denses
#   ni trop clairsemés.
#
# TECHNIQUE :
#   On précalcule les distances au plus proche voisin (tree.query k=1)
#   pour toutes les molécules, ce qui permet de calculer les scores
#   pour tous les rayons en une seule passe → très rapide.
# ═══════════════════════════════════════════════════════════════

def generer_graphe_html(coords1, coords2, rayon_choisi, msg_limite, dossier_sortie, prefixe=""):
    """
    Génère le fichier HTML du graphe interactif — un seul panneau :
    Score de colocalisation (W1, W2, global) en fonction du rayon.

    On précalcule la distance au voisin le plus proche pour chaque
    localisation → permet de calculer les scores pour tous les rayons
    en une seule passe sans refaire un KD-Tree à chaque fois.

    La ligne bleue verticale marque le rayon actuellement choisi.
    La ligne orange horizontale marque le seuil de 20%.
    """

    # ── Calcul des scores pour tous les rayons de 10 à 635 nm ──
    rayons = list(range(10, 650, 15))
    scores_w1_r, scores_w2_r, scores_g_r = [], [], []

    tree2 = cKDTree(coords2)
    tree1 = cKDTree(coords1)
    # Distance de chaque localisation W1 à son voisin W2 le plus proche
    dists_1vers2, _ = tree2.query(coords1, k=1)
    # Distance de chaque localisation W2 à son voisin W1 le plus proche
    dists_2vers1, _ = tree1.query(coords2, k=1)

    for r in rayons:
        # Une localisation est "colocalisée" si son voisin le plus proche est ≤ r
        sw1 = 100 * (dists_1vers2 <= r).sum() / len(coords1)
        sw2 = 100 * (dists_2vers1 <= r).sum() / len(coords2)
        scores_w1_r.append(sw1)
        scores_w2_r.append(sw2)
        scores_g_r.append(np.sqrt(sw1 * sw2))

    # ── Construction du graphe Plotly (un seul panneau) ────────
    fig = go.Figure()

    # Courbe Score W1 (vert)
    fig.add_trace(go.Scatter(
        x=rayons, y=scores_w1_r, name="Score W1",
        mode="lines", line=dict(color="#00e676", width=2.5),
        hovertemplate="Rayon=%{x} nm<br>Score W1=%{y:.1f}%<extra></extra>"
    ))

    # Courbe Score W2 (rouge)
    fig.add_trace(go.Scatter(
        x=rayons, y=scores_w2_r, name="Score W2",
        mode="lines", line=dict(color="#ff1744", width=2.5),
        hovertemplate="Rayon=%{x} nm<br>Score W2=%{y:.1f}%<extra></extra>"
    ))

    # Courbe Score global (jaune pointillé)
    fig.add_trace(go.Scatter(
        x=rayons, y=scores_g_r, name="Score global (moy. géométrique)",
        mode="lines", line=dict(color="#ffd600", width=3, dash="dot"),
        hovertemplate="Rayon=%{x} nm<br>Score global=%{y:.1f}%<extra></extra>"
    ))

    # Ligne verticale bleue = rayon actuellement choisi dans l'interface
    fig.add_vline(
        x=rayon_choisi,
        line=dict(color="#4fc3f7", width=2, dash="dash"),
        annotation_text=f"Rayon actuel : {rayon_choisi} nm",
        annotation_font_color="#4fc3f7",
        annotation_position="top right"
    )

    # Ligne horizontale orange = seuil de significativité (20%)
    fig.add_hline(
        y=20,
        line=dict(color="#f59e0b", width=1.5, dash="longdash"),
        annotation_text="Seuil significatif (20 %)",
        annotation_font_color="#f59e0b",
        annotation_position="bottom right"
    )

    # Mise en forme globale
    fig.update_layout(
        title=dict(
            text="Sensibilité du score de colocalisation en fonction du rayon",
            font=dict(size=16, color="#e8eaf6"), x=0.5
        ),
        xaxis_title="Rayon de colocalisation (nm)",
        yaxis_title="Score (%)",
        yaxis=dict(range=[0, 102]),
        paper_bgcolor="#0f1117",
        plot_bgcolor="#0d1020",
        font=dict(color="#e8eaf6", family="Courier New"),
        legend=dict(bgcolor="#1a1d27", bordercolor="#4fc3f7", borderwidth=1,
                    font=dict(size=11)),
        hovermode="x unified",
        height=520,
        xaxis=dict(gridcolor="#1e2a3a", linecolor="#4fc3f7",
                   tickfont=dict(color="#90a4ae"), title_font=dict(color="#90a4ae")),
        yaxis2=dict(gridcolor="#1e2a3a", linecolor="#4fc3f7",
                    tickfont=dict(color="#90a4ae"), title_font=dict(color="#90a4ae")),
    )

    # Avertissement conditionnel : affiché uniquement si msg_limite n'est pas vide
    # (c'est-à-dire si un canal est plus de 3x plus dense que l'autre)
    if msg_limite:
        fig.add_annotation(
            text=("⚠  " + msg_limite.replace("\n", "<br>")),
            xref="paper", yref="paper", x=0.5, y=-0.15,
            showarrow=False, font=dict(size=10, color="#f59e0b"), align="center"
        )

    # Sauvegarde HTML autonome (librairie Plotly intégrée dans le fichier)
    path_html = os.path.join(dossier_sortie, f"{prefixe}graphe_sensibilite.html")
    fig.write_html(path_html, include_plotlyjs="cdn")
    return path_html


# ═══════════════════════════════════════════════════════════════
# SECTION 8 — GÉNÉRATION DES IMAGES TIF
# Deux images sont produites en format TIF 300 DPI (ou plus si
# tu changes le paramètre dpi) ouvrable dans ImageJ.
#
# Image 1 — Référence :
#   Nuage de points brut. Points W1 en vert, W2 en rouge.
#   Aucune liaison tracée. Sert de vue de base pour comparer
#   avec l'image 2.
#
# Image 2 — Colocalisation Delaunay :
#   - Arêtes vertes sombres : réseau de Delaunay intra W1
#   - Arêtes rouges sombres : réseau de Delaunay intra W2
#   - Traits jaunes : liaisons inter-canal W1↔W2 colocalisées
#   - Points verts/rouges par-dessus les arêtes (zorder=4)
#   - Encadré bas-gauche : score, verdict, et limite
#
# zorder contrôle l'ordre d'empilement des éléments dessinés :
#   zorder=1 → en dessous (arêtes)
#   zorder=3 → au milieu (liaisons jaunes)
#   zorder=4 → au-dessus (points, pour qu'ils soient toujours visibles)
# ═══════════════════════════════════════════════════════════════

def _setup_ax(ax, titre):
    """Configure le style visuel d'un panneau matplotlib (fond sombre, axes colorés)."""
    ax.set_facecolor("#050810")
    ax.set_title(titre, color=TXT, fontsize=10, fontweight="bold", pad=10)
    ax.tick_params(colors=TXT2, labelsize=7)
    for sp in ax.spines.values():
        sp.set_edgecolor(ACCENT); sp.set_linewidth(0.6)
    ax.set_xlabel("X (nm)", color=TXT2, fontsize=8)
    ax.set_ylabel("Y (nm)", color=TXT2, fontsize=8)


def generer_images_tif(res, dossier_sortie, prefixe="", dpi=300):
    """
    Génère les deux fichiers TIF et les sauvegarde dans dossier_sortie.
    res  : dictionnaire retourné par analyser_colocalisation()
    dpi  : résolution en points par pouce (300 = standard, 600 = haute résolution)
           Pour passer à 600 DPI : changer dpi=300 en dpi=600 dans _run()
    """
    # Récupération des données depuis le dictionnaire de résultats
    coords1_plot = res["coords1_plot"]   # points W1 affichés
    coords2_plot = res["coords2_plot"]   # points W2 affichés
    coords1      = res["coords1"]        # tous les points W1 (pour arêtes)
    coords2      = res["coords2"]        # tous les points W2 (non utilisé ici)
    coords2_del  = res["coords2_del"]    # points W2 pour tracer les arêtes Delaunay
    edges1       = res["edges1"]         # arêtes réseau W1
    edges2       = res["edges2"]         # arêtes réseau W2
    paires       = res["paires"]         # liaisons jaunes à tracer

    # Légende commune aux deux images
    leg_base = [
        Line2D([0],[0], marker='o', color='none', markerfacecolor=C1_COLOR,
               markersize=5, label=f"W1  ({res['n_w1']:,} loc.)"),
        Line2D([0],[0], marker='o', color='none', markerfacecolor=C2_COLOR,
               markersize=5, label=f"W2  ({res['n_w2']:,} loc.)"),
    ]

    # ── Image 1 : Référence (nuage brut sans liaisons) ────────
    fig1, ax1 = plt.subplots(figsize=(14, 11), facecolor=BG)
    _setup_ax(ax1, "Image de référence — Distribution spatiale des deux canaux")
    ax1.scatter(coords1_plot[:,0], coords1_plot[:,1], s=1.0, c=C1_COLOR,
                alpha=0.7, linewidths=0, rasterized=True)
    ax1.scatter(coords2_plot[:,0], coords2_plot[:,1], s=0.6, c=C2_COLOR,
                alpha=0.5, linewidths=0, rasterized=True)
    ax1.legend(handles=leg_base, facecolor=BG2, edgecolor=ACCENT,
               labelcolor=TXT, fontsize=9, loc="upper right")
    ax1.set_aspect("equal")   # garantit que 1 nm en X = 1 nm en Y (pas de déformation)
    path1 = os.path.join(dossier_sortie, f"{prefixe}1_reference.tif")
    fig1.savefig(path1, format="tiff", dpi=dpi, bbox_inches="tight", facecolor=BG)
    plt.close(fig1)   # libère la mémoire

    # ── Image 2 : Colocalisation avec réseaux et liaisons ─────
    fig2, ax2 = plt.subplots(figsize=(14, 11), facecolor=BG)
    _setup_ax(ax2, "Réseaux Delaunay + liaisons de colocalisation inter-canal")

    # Arêtes W1 en vert sombre (zorder=1 → en dessous)
    # LINE_WIDTH = épaisseur uniforme définie en Section 2
    for (a, b) in edges1:
        ax2.plot([coords1[a,0], coords1[b,0]], [coords1[a,1], coords1[b,1]],
                 color=E1_COLOR, alpha=0.25, linewidth=LINE_WIDTH, zorder=1)

    # Arêtes W2 en rouge sombre (zorder=1 → en dessous)
    for (a, b) in edges2:
        ax2.plot([coords2_del[a,0], coords2_del[b,0]],
                 [coords2_del[a,1], coords2_del[b,1]],
                 color=E2_COLOR, alpha=0.20, linewidth=LINE_WIDTH, zorder=1)

    # Liaisons jaunes inter-canal (zorder=3 → au-dessus des arêtes)
    # Chaque paire p = [x1, y1, x2, y2] : départ W1 → arrivée W2
    # Les deux points sont dans les sous-ensembles affichés → pas de trait dans le vide
    if len(paires) > 0:
        for p in paires:
            ax2.plot([p[0], p[2]], [p[1], p[3]],
                     color=LINK_COLOR, alpha=0.6, linewidth=LINE_WIDTH, zorder=3)

    # Points par-dessus tout (zorder=4 → toujours visibles)
    ax2.scatter(coords1_plot[:,0], coords1_plot[:,1], s=1.2, c=C1_COLOR,
                alpha=0.8, linewidths=0, rasterized=True, zorder=4)
    ax2.scatter(coords2_plot[:,0], coords2_plot[:,1], s=0.7, c=C2_COLOR,
                alpha=0.6, linewidths=0, rasterized=True, zorder=4)

    # Légende avec les 5 éléments
    leg2 = leg_base + [
        Line2D([0],[0], color=E1_COLOR, linewidth=1.2,
               label=f"Réseau Delaunay W1  ({len(edges1):,} arêtes)"),
        Line2D([0],[0], color=E2_COLOR, linewidth=1.2,
               label=f"Réseau Delaunay W2  ({len(edges2):,} arêtes)"),
        Line2D([0],[0], color=LINK_COLOR, linewidth=1.5,
               label=f"Liaisons coloc. W1<->W2  ({len(paires):,})"),
    ]
    ax2.legend(handles=leg2, facecolor=BG2, edgecolor=ACCENT,
               labelcolor=TXT, fontsize=8, loc="upper right")
    ax2.set_aspect("equal")

    # Encadré de résultats en bas à gauche de l'image
    # Le message de limite (msg_limite) est conditionnel :
    # il ne s'affiche que si un canal est plus de 3x plus dense que l'autre.
    # Si les canaux sont équilibrés, msg_limite est vide → rien n'est affiché.
    base = (
        f"VERDICT : {res['verdict']}\n\n"
        f"Score W1 : {res['score_w1']:.1f}%  ({res['n_w1_coloc']:,}/{res['n_w1']:,} loc.)\n"
        f"Score W2 : {res['score_w2']:.1f}%  ({res['n_w2_coloc']:,}/{res['n_w2']:,} loc.)\n"
        f"Score global (moy. geometrique) : {res['score_global']:.1f}%\n"
        f"Distance moyenne de colocalisation : {res['dist_moyenne']:.1f} nm\n"
        f"  [= somme longueurs traits jaunes / nbre traits jaunes]\n"
        f"Rayon : {res['rayon_nm']} nm  |  Arete Delaunay max : {res['max_edge_nm']} nm"
    )
    # Ajout conditionnel du message de limite
    texte = base + ("\n\n" + res['msg_limite'] if res['msg_limite'] else "")
    ax2.text(0.01, 0.01, texte, transform=ax2.transAxes,
             fontsize=7, color=TXT, verticalalignment="bottom",
             fontfamily="monospace",
             bbox=dict(boxstyle="round,pad=0.6", facecolor=BG2,
                       edgecolor=ACCENT, alpha=0.92))

    path2 = os.path.join(dossier_sortie, f"{prefixe}2_colocalisation_delaunay.tif")
    fig2.savefig(path2, format="tiff", dpi=dpi, bbox_inches="tight", facecolor=BG)
    plt.close(fig2)

    return path1, path2


# ═══════════════════════════════════════════════════════════════
# SECTION 9 — INTERFACE GRAPHIQUE (tkinter)
# Toute la partie visuelle de l'application : fenêtre, onglets,
# champs de saisie, bouton, et zone de log.
#
# La classe App hérite de tk.Tk (la fenêtre principale tkinter).
# Elle contient :
#   _build_ui()  : construit tous les éléments visuels
#   _ligne()     : crée une ligne "label + champ + bouton Parcourir"
#   _param()     : crée une ligne "label + champ numérique + info"
#   _bw1/w2/tif(): ouvrent les boîtes de dialogue de sélection de fichier
#   _log()       : affiche un message horodaté dans la zone de log
#   _lancer()    : valide les entrées et démarre l'analyse
#   _run()       : exécute l'analyse complète dans un thread séparé
#
# POURQUOI UN THREAD SÉPARÉ (_run dans threading.Thread) ?
#   L'analyse prend plusieurs secondes. Si on la lançait dans le
#   thread principal, l'interface se figerait et serait inutilisable
#   pendant tout ce temps. Le thread séparé permet à l'interface
#   de rester réactive (logs qui s'affichent, fenêtre déplaçable)
#   pendant que le calcul tourne en arrière-plan.
# ═══════════════════════════════════════════════════════════════

class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Analyse de Colocalisation STORM — v3")
        self.configure(bg=BG)
        self.resizable(True, True)
        self.minsize(720, 680)

        # Variables tkinter liées aux champs de saisie
        self.path_w1  = tk.StringVar()   # chemin fichier CSV W1
        self.path_w2  = tk.StringVar()   # chemin fichier CSV W2
        self.path_tif = tk.StringVar()   # chemin image TIF
        self.rayon    = tk.StringVar(value="50")   # rayon par défaut : 50 nm
        self.max_edge = tk.StringVar(value="200")   # arête max par défaut : 500 nm
        self.mode     = tk.StringVar(value="csv")   # mode actif : "csv" ou "tif"

        self._build_ui()

    def _build_ui(self):
        """Construit tous les éléments visuels de la fenêtre."""
        # En-tête
        tk.Label(self, text="  Analyse de Colocalisation STORM ",
                 bg=BG, fg=ACCENT, font=("Courier New", 14, "bold")).pack(pady=(16,2))
        tk.Label(self,
                 text="Marquage d'anticorps — Delaunay + Graphe interactif HTML",
                 bg=BG, fg=TXT2, font=("Courier New", 8)).pack(pady=(0,10))

        # Onglets CSV / TIF
        nb = ttk.Notebook(self)
        nb.pack(fill="x", padx=18, pady=4)
        st = ttk.Style(); st.theme_use("clam")
        st.configure("TNotebook",     background=BG2, borderwidth=0)
        st.configure("TNotebook.Tab", background=BG2, foreground=TXT2,
                     padding=[12,6], font=("Courier New",9))
        st.map("TNotebook.Tab",
               background=[("selected", ACCENT)], foreground=[("selected", BG)])

        f_csv = tk.Frame(nb, bg=BG2)
        f_tif = tk.Frame(nb, bg=BG2)
        nb.add(f_csv, text="  Fichiers CSV  ")
        nb.add(f_tif, text="  Image TIF     ")
        # Quand l'utilisateur change d'onglet, on met à jour self.mode
        nb.bind("<<NotebookTabChanged>>",
                lambda e: self.mode.set("csv" if nb.index("current")==0 else "tif"))

        self._ligne(f_csv, "Fichier CSV W1 (ex. Calnexine) :", self.path_w1, self._bw1)
        self._ligne(f_csv, "Fichier CSV W2 (ex. CFTR) :",      self.path_w2, self._bw2)
        self._ligne(f_tif, "Image TIF STORM (RGB : vert=W1, rouge=W2) :", self.path_tif, self._btif)
        tk.Label(f_tif, text="  Canal vert -> W1  |  Canal rouge -> W2",
                 bg=BG2, fg=TXT2, font=("Courier New",8)).pack(anchor="w", padx=14, pady=(2,8))

        tk.Frame(self, bg=ACCENT, height=1).pack(fill="x", padx=18, pady=8)

        # Zone des paramètres numériques
        fp = tk.Frame(self, bg=BG); fp.pack(fill="x", padx=18)
        self._param(fp, "Rayon de colocalisation (nm) :", self.rayon,
                    "distance max W1<->W2 pour tracer une liaison")
        self._param(fp, "Arête Delaunay max (nm) :",      self.max_edge,
                    "longueur max des mailles (300-800 nm recommandé)")

        tk.Frame(self, bg=ACCENT, height=1).pack(fill="x", padx=18, pady=8)

        # Bouton principal
        btn = tk.Button(self, text="  LANCER L'ANALYSE",
                        bg=BTN_OK, fg="white", activebackground=BTN_HOV,
                        font=("Courier New",11,"bold"), relief="flat",
                        cursor="hand2", pady=10, command=self._lancer)
        btn.pack(fill="x", padx=18, pady=2)

        # Effet de survol sur le bouton
        btn.bind("<Enter>", lambda e: btn.config(bg=BTN_HOV))
        btn.bind("<Leave>", lambda e: btn.config(bg=BTN_OK))

        tk.Frame(self, bg=ACCENT, height=1).pack(fill="x", padx=18, pady=8)

        # Zone de log scrollable (affiche les messages de progression)
        self.log = scrolledtext.ScrolledText(
            self, height=12, bg=BG2, fg=TXT, font=("Courier New",9),
            insertbackground=ACCENT, relief="flat", wrap="word", state="disabled")
        self.log.pack(fill="both", expand=True, padx=18, pady=(0,12))

        self._log("Programme prêt. Charge tes fichiers et lance l'analyse.")
        self._log("3 fichiers générés : TIF référence + TIF colocalisation + HTML graphe.")

    def _ligne(self, parent, label, var, cmd):
        """Crée une ligne avec label + champ texte + bouton Parcourir."""
        frm = tk.Frame(parent, bg=BG2); frm.pack(fill="x", padx=12, pady=5)
        tk.Label(frm, text=label, bg=BG2, fg=TXT, font=("Courier New",9)).pack(anchor="w")
        row = tk.Frame(frm, bg=BG2); row.pack(fill="x")
        tk.Entry(row, textvariable=var, bg=ENTRY_BG, fg=ACCENT2,
                 insertbackground=ACCENT2, font=("Courier New",8),
                 relief="flat", bd=3).pack(side="left", fill="x", expand=True, padx=(0,6))
        tk.Button(row, text="Parcourir...", bg=BG, fg=TXT2,
                  activebackground=ACCENT, font=("Courier New",8),
                  relief="flat", cursor="hand2", command=cmd).pack(side="left")

    def _param(self, parent, label, var, info):
        """Crée une ligne avec label + champ numérique + texte d'information."""
        row = tk.Frame(parent, bg=BG); row.pack(fill="x", pady=3)
        tk.Label(row, text=label, bg=BG, fg=TXT, font=("Courier New",9),
                 width=35, anchor="w").pack(side="left")
        tk.Entry(row, textvariable=var, width=8, bg=ENTRY_BG, fg=ACCENT,
                 insertbackground=ACCENT, font=("Courier New",10),
                 relief="flat", bd=4).pack(side="left", padx=8)
        tk.Label(row, text=info, bg=BG, fg=TXT2,
                 font=("Courier New",8)).pack(side="left")

    def _bw1(self):
        """Ouvre la boîte de dialogue pour sélectionner le fichier CSV W1."""
        p = filedialog.askopenfilename(title="CSV W1", filetypes=[("CSV","*.csv"),("Tous","*.*")])
        if p: self.path_w1.set(p)

    def _bw2(self):
        """Ouvre la boîte de dialogue pour sélectionner le fichier CSV W2."""
        p = filedialog.askopenfilename(title="CSV W2", filetypes=[("CSV","*.csv"),("Tous","*.*")])
        if p: self.path_w2.set(p)

    def _btif(self):
        """Ouvre la boîte de dialogue pour sélectionner l'image TIF."""
        p = filedialog.askopenfilename(title="TIF STORM",
            filetypes=[("TIF/TIFF","*.tif *.tiff"),("Tous","*.*")])
        if p: self.path_tif.set(p)

    def _log(self, msg):
        """Affiche un message horodaté dans la zone de log de l'interface."""
        self.log.config(state="normal")
        ts = datetime.datetime.now().strftime("%H:%M:%S")
        self.log.insert("end", f"[{ts}]  {msg}\n")
        self.log.see("end")          # défilement automatique vers le bas
        self.log.config(state="disabled")
        self.update_idletasks()      # force le rafraîchissement de l'interface

    def _lancer(self):
        """
        Vérifie les entrées utilisateur avant de lancer l'analyse.
        Si tout est valide, démarre _run() dans un thread séparé.
        """
        try:
            rayon = float(self.rayon.get())
            max_e = float(self.max_edge.get())
            if rayon <= 0 or max_e <= 0: raise ValueError
        except ValueError:
            messagebox.showerror("Erreur", "Rayon et arête max doivent être des nombres positifs.")
            return
        mode = self.mode.get()
        if mode == "csv" and (not self.path_w1.get() or not self.path_w2.get()):
            messagebox.showerror("Erreur", "Charge les deux fichiers CSV.")
            return
        if mode == "tif" and not self.path_tif.get():
            messagebox.showerror("Erreur", "Charge une image TIF.")
            return
        # Lance l'analyse dans un thread séparé pour ne pas bloquer l'interface
        threading.Thread(target=self._run, args=(mode, rayon, max_e), daemon=True).start()

    def _run(self, mode, rayon, max_edge):
        """
        Exécute l'analyse complète (appelé dans un thread séparé).
        Enchaîne : chargement → analyse → génération TIF → génération HTML.
        """
        self._log("─" * 54)
        self._log(f"Démarrage | Rayon={rayon} nm | Arête max={max_edge} nm")
        try:
            # Chargement selon le mode sélectionné
            if mode == "csv":
                self._log(f"Chargement W1 : {os.path.basename(self.path_w1.get())}")
                self._log(f"Chargement W2 : {os.path.basename(self.path_w2.get())}")
                df1, df2 = charger_csv(self.path_w1.get(), self.path_w2.get())
            else:
                self._log(f"Chargement TIF : {os.path.basename(self.path_tif.get())}")
                df1, df2 = charger_tif(self.path_tif.get())

            self._log(f"W1 : {len(df1):,} localisations  |  W2 : {len(df2):,} localisations")

            # Analyse complète (Delaunay + liaisons + score)
            res = analyser_colocalisation(df1, df2, rayon, max_edge, log_fn=self._log)

            # Affichage des résultats dans le log
            self._log("─" * 54)
            self._log("RÉSULTATS")
            self._log(f"  W1 colocalisées : {res['n_w1_coloc']:,}/{res['n_w1']:,} ({res['score_w1']:.2f}%)")
            self._log(f"  W2 colocalisées : {res['n_w2_coloc']:,}/{res['n_w2']:,} ({res['score_w2']:.2f}%)")
            self._log(f"  Score global    : {res['score_global']:.2f}%")
            self._log(f"  Distance moy. coloc. : {res['dist_moyenne']:.1f} nm")
            self._log(f"  VERDICT : {res['verdict']}")
            self._log("  ATTENTION : score potentiellement surestimé (densité W2 élevée).")

            # Dossier de sortie = même dossier que celui des fichiers chargés
            if mode == "csv":
                dossier = os.path.dirname(os.path.abspath(self.path_w1.get()))
            else:
                dossier = os.path.dirname(os.path.abspath(self.path_tif.get()))
            # Horodatage pour nommer les fichiers de façon unique
            ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            pref = f"COLOC_{ts}_"

            # Génération des deux images TIF
            self._log("Génération des images TIF 300 DPI...")
            p1, p2 = generer_images_tif(res, dossier, prefixe=pref, dpi=300)
            self._log(f"  TIF 1 : {os.path.basename(p1)}")
            self._log(f"  TIF 2 : {os.path.basename(p2)}")

            # Génération du graphe HTML interactif
            self._log("Génération du graphe HTML interactif...")
            ph = generer_graphe_html(
                res["coords1"], res["coords2"],
                rayon, res["msg_limite"], dossier, prefixe=pref
            )
            self._log(f"  HTML  : {os.path.basename(ph)}")
            self._log("Terminé ! Ouvre les TIF dans ImageJ et le HTML dans le navigateur.")
            self._log("─" * 54)

        except Exception as e:
            # En cas d'erreur inattendue, on l'affiche clairement
            self._log(f"ERREUR : {e}")
            messagebox.showerror("Erreur", str(e))


# ═══════════════════════════════════════════════════════════════
# SECTION 10 — POINT D'ENTRÉE
# Ces deux lignes sont le point de départ du programme.
# "if __name__ == '__main__'" signifie : "n'exécute ce bloc que
# si ce fichier est lancé directement" (et non importé par un
# autre script Python).
# app.mainloop() démarre la boucle d'événements tkinter qui
# maintient la fenêtre ouverte et réactive jusqu'à ce que
# l'utilisateur la ferme.
# ═══════════════════════════════════════════════════════════════

if __name__ == "__main__":
    app = App()
    app.mainloop()
