# Livrables de mon stage de licence 3 de Génie Bio-Informatique à la plateforme IMAGEUP de l'université de Poitiers (Mai à Juillet 2026)
Stage portant sur la mise en place d’un logiciel dédié à l’analyse d’images issues de la microscopie haute résolution STORM.

## Projet 1 — Analyse de colocalisation STORM par triangulation de Delaunay

L'objectif ici est de déterminer s'il existe une trame structurelle commune qui va au-delà du simple marquage de deux protéines, c'est-à-dire de quantifier si deuxpopulations de molécules détectées par microscopie STORM (canaux W1 et W2) occupent des positions spatialement proches, ce qui peut indiquer une interaction fonctionnelle.

L'application fonctionne avec deux modes d'entrée : soit deux fichiers CSV de
localisations (un par canal, avec coordonnées x/y en nanomètres), soit une image TIF
STORM bicolore dont les canaux sont extraits automatiquement par seuillage.

Pour chaque canal, une triangulation de Delaunay est calculée afin de reconstruire le
réseau structurel sous-jacent (pertinent notamment pour des organites comme le
réticulum endoplasmique). La colocalisation inter-canal est ensuite quantifiée via une
recherche de voisinage par KD-Tree dans un rayon paramétrable, produisant des scores
de colocalisation (W1, W2, score global), une distance moyenne entre paires
colocalisées, ainsi qu'un graphe interactif de sensibilité montrant l'évolution des
scores selon le rayon choisi.

L'interface graphique (tkinter) a été conçue pour être utilisable sans connaissances
en programmation, avec un guide d'utilisation fourni.

**Stack** : Python, tkinter, scipy (Delaunay, cKDTree), matplotlib, plotly

## Projet 2 — Image Art

Application web locale permettant d'appliquer des effets artistiques (vitrail et mosaïque) à des images de microscopie, développée dans le cadre de mon stage de licence 3 à la plateforme **ImageUP** (Université de Poitiers).

### Contexte

Plusieurs étudiants de M1 GPHY ont développé des macros ImageJ/Fiji produisant des effets visuels de type "vitrail" à partir d'images de microscopie. Ces macros étaient fonctionnelles mais nécessitaient de connaître ImageJ pour être utilisées. L'objectif de ce projet est de rendre ces traitements accessibles via une interface web simple (sans dépendance à Fiji), qui sera hébergée sur le site de la plateforme, afin que n'importe quel membre de la plateforme ou visiteur du site web puisse générer ces rendus en quelques clics.

### Fonctionnalités

- Galerie d'images gérée directement par dépôt de fichiers dans un dossier (aucune interface d'administration nécessaire)
- Sélection d'une image puis d'un effet à appliquer
- Réglage des paramètres de chaque effet via des curseurs interactifs
- Affichage immédiat du résultat dans l'interface
- Téléchargement du résultat au format TIFF

### Effets disponibles

**Vitrail — version Maël Zami** : flou gaussien, quantification des couleurs, détection et épaississement des contours (effet "plomb").

**Vitrail — version Cléo Thury** : segmentation de l'image en zones homogènes (superpixels), couleur moyenne par zone, boost de saturation et choix de teinte.

**Mosaïque** : découpage de l'image en carreaux carrés, chaque carreau prenant la couleur du pixel central avec un boost de brillance ajustable (rendu pixel art).

D'autres versions de l'effet vitrail, portées depuis les macros originales de plusieurs étudiants M1 ont été développées et testées au cours du stage mais ne sont pas activées dans la version retenue par la plateforme.

### Architecture technique

```
Navigateur (HTML/CSS/JavaScript)
        │
        ▼
Serveur Flask (Python)
        │  Pillow + OpenCV pour le traitement d'image
        ▼
Dossier outputs/ — fichiers résultats au format TIFF
```

L'application a initialement été conçue pour piloter Fiji en ligne de commande à partir des macros originales (ImageJ Macro Language), avant d'être entièrement portée en Python (Pillow, OpenCV, NumPy) afin de s'affranchir de la dépendance à Fiji et d'afficher le résultat directement dans l'interface.

### Stack technique

- **Python / Flask** — serveur web et logique de traitement
- **Pillow** — lecture, conversion et quantification des couleurs
- **OpenCV** — flou gaussien, détection de contours, segmentation, morphologie
- **NumPy** — manipulation des tableaux de pixels
- **HTML / CSS / JavaScript** — interface utilisateur
- **Jinja2** — génération dynamique de la galerie d'images

### Installation

```bash
git clone <url-du-depot>
cd vitrail_app
pip install flask pillow opencv-contrib-python numpy
python Image_Art.py
```

L'application est ensuite accessible à l'adresse `http://127.0.0.1:5000`.

### Structure du dépôt

```
vitrail_app/
├── Image_Art.py          # serveur Flask et fonctions de traitement d'image
├── images/                # images disponibles dans la galerie (à gérer manuellement)
├── outputs/               # résultats générés, téléchargeables depuis l'interface
├── macros/                # macros ImageJ originales, conservées à titre d'archive
└── templates/
    └── index.html          # interface utilisateur
```

### Utilisation

1. Lancer `Image_Art.py`
2. Ouvrir `http://127.0.0.1:5000` dans un navigateur
3. Sélectionner une image dans la galerie affichée dans l'interface
4. Choisir un effet (vitrail ou mosaïque)
5. Ajuster les paramètres souhaités
6. Lancer le traitement
7. Télécharger le résultat

Pour ajouter ou retirer des images de la galerie, il suffit d'ajouter ou de supprimer les fichiers correspondants dans le dossier `images/`.

### Auteur

Prudencio AYENAN — stagiaire en Bio-Informatique, à la plateforme ImageUP, Université de Poitiers, 2025-2026

### Remerciements

Les macros ImageJ originales ayant servi de base aux effets vitrail ont été développées par plusieurs étudiants de Master 1 de la plateforme ImageUP : Maël Zami, Cléo Thury, Dorian Esserteau, Kévin Pérès, Samuel Maillé et Rachel Head.
