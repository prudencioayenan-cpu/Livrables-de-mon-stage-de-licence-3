# Livrables de mon stage de licence 3 de Génie Bio-Informatique à la plateforme IMAGEUP de l'université de Poitiers (Mai à Juillet 2026)
Stage portant sur la mise en place d’un logiciel dédié à l’analyse d’images issues de la microscopie haute résolution STORM.

## Projet 1 — Analyse de colocalisation STORM par triangulation de Delaunay

L'objectif ici est de déterminer s'il existe une trame structurelle commune qui va au-delà du simple marquage de deux protéines, c'est-à-dire de quantifier si deux populations de molécules détectées par microscopie STORM (canaux W1 et W2) occupent des positions spatialement proches, ce qui peut indiquer une interaction fonctionnelle, voir une colocalisation entre ces protéines.

L'application fonctionne avec un mode d'entrée : deux fichiers CSV de
localisations (un par canal, avec coordonnées x/y en nanomètres).

Pour chaque canal, une triangulation de Delaunay est calculée afin de reconstruire le
réseau structurel sous-jacent (pertinent notamment pour des organites comme le
réticulum endoplasmique). La colocalisation inter-canal est ensuite quantifiée via une
recherche de voisinage par KD-Tree dans un rayon paramétrable, produisant des scores
de colocalisation (W1, W2, score global), une distance moyenne entre paires
colocalisées, ainsi qu'un graphe interactif de sensibilité montrant l'évolution des
scores selon le rayon choisi.

L'interface graphique a été conçue avec *tkinter* pour être utilisable sans connaissances
en programmation, avec un guide d'utilisation fourni.

**Bibliothèque de code** : Python, tkinter, scipy (Delaunay, cKDTree), matplotlib, plotly

### Installation et Utilisation

## >> 1- Télécharger d'abord Python version 3.14.5 : Python 3.14.5 

## >> 2- Installer les packages requis via le fichier ‘requirements.txt’ 
	Pour installer les packages, se mettre dans le bon dossier dans le terminal windows(‘touche Win +R’ puis cmd ou bien taper directement ‘terminal’ dans la 		barre de recherche du pc), puis une fois dans le terminal taper la commande : pip install *-r requirements.txt* . 

## >> 3- LANCER LE LOGICIEL 
Une fois que Python et les packages sont installés, double cliquer sur le fichier *'analyse_colocalisation_STORM.py'*. Cela ouvre directement le logiciel dans un nouvel onglet.  

## >> 4- COMMENT L'UTILISER ?
Le logiciel présente deux modes d'entrée de fichiers : le mode CSV et le mode TIF 
	- Le mode CSV (celui qui nous intéresse plus ici, puisqu'étant le plus fiable) : vous entrez vos fichiers CSV (le W1 en premier puis le W2 par la suite) ;
	- Le mode TIF : une seule image bicolore (à double marquage) est demandée. 
Cliquez ensuite sur le bouton « Lancer l’analyse ».

## >> 5- ET ALORS QU’EST-CE QUI SE PASSE PAR LA SUITE ? 
Après avoir lancé l’analyse attendez quelques secondes. Dès que vous voyez en bas « Terminé ! Ouvre les TIF dans Image J et le html dans le navigateur », allez dans le dossier dans lequel vous avez pris vos CSV, vous verrez trois fichiers :

	- Le premier : image TIF de référence qui présente le double marquage. C’est l’image de base générée à partir des fichiers CSV.
	- Le deuxième : image TIF qui présente la colocalisation à partir de la triangulation de Delaunay.
	- Le troisième : graphe qui présente l’évolution du score de colocalisation en fonction des paramètres configurés au lancement de l'analyse.


## Projet 2 — Image Art

Application web locale permettant d'appliquer des effets artistiques (vitrail et mosaïque) à des images de microscopie, développée dans le cadre de mon stage de licence 3 à la plateforme **ImageUP** (Université de Poitiers), pour un projet artistique (Arts et Sciences).

### Contexte

Sur demande par ma maitre de stage au cours de l'année scolaire, plusieurs étudiants de master 1 de ma formation (GPHY : Génie Physiologique Biotechnologique et Informatique) ont développé des macros ImageJ/Fiji produisant des effets visuels de type "vitrail" à partir d'images de microscopie. Ces macros étaient fonctionnelles mais nécessitaient de connaître ImageJ pour être utilisées.

L'objectif de ce projet est donc de rendre ces traitements accessibles via une interface web simple (sans dépendance à Fiji), qui sera éventuellement hébergée sur le site web de la plateforme, afin que n'importe quel membre de la plateforme ou visiteur du site web puisse générer ces rendus en quelques clics.

### Fonctionnalités

- Galerie d'images gérée directement par dépôt de fichiers dans un dossier (aucune interface d'administration nécessaire)
- Sélection d'une image puis d'un effet à appliquer
- Réglage des paramètres de chaque effet via des curseurs interactifs
- Affichage immédiat du résultat dans l'interface
- Téléchargement du résultat au format TIFF

### Effets disponibles

**Vitrail 1** : flou gaussien, quantification des couleurs, détection et épaississement des contours (effet "plomb").

**Vitrail 2** : segmentation de l'image en zones homogènes (superpixels), couleur moyenne par zone, boost de saturation et choix de teinte.

**Mosaïque** : découpage de l'image en carreaux carrés, chaque carreau prenant la couleur du pixel central avec un boost de brillance ajustable (rendu pixel art).

D'autres versions de l'effet vitrail, portées depuis les macros originales d'autres étudiants M1 ont été développées et testées au cours du stage mais ne sont pas activées dans la version retenue par la plateforme.

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

### Bibliothèque technique

- **Python / Flask** — serveur web et logique de traitement
- **Pillow** — lecture, conversion et quantification des couleurs
- **OpenCV** — flou gaussien, détection de contours, segmentation, morphologie
- **NumPy** — manipulation des tableaux de pixels
- **HTML / CSS / JavaScript** — interface utilisateur
- **Jinja2** — génération dynamique de la galerie d'images

### Installation et Utilisation

## >> 1- Télécharger d'abord Python version 3.14.5 : Python 3.14.5

## 2- Installer les packages requis via le fichier ‘requirements.txt’
	Pour installer les packages, se mettre dans le bon dossier dans le terminal windows(‘touche Win +R’ puis cmd ou bien taper directement ‘terminal’ dans la 		barre de recherche du pc), puis une fois dans le terminal taper la commande : pip install *-r requirements.txt*.

## >> 3- Lancer l'appli
Une fois que Python et les packages sont installés, double cliquer sur le fichier ‘lancer_ImageArt’, qui est un .bat. Cela ouvre directement le serveur Flask et lance l'application dans le navigateur. Le serveur Flask s’ouvre d’abord via un terminal (puisque c’est lui qui fait tourner l’application, puis l’application elle s’ouvre par la suite quelques secondes d’attente après. 

## >> 4- Dès que vous le navigateur s'ouvre, vous verrez une page avec sécurisée qui vous avertit du fait que la page n'est pas sécurisée etc... Pas de panique c'est juste une mesure de sécurité (certificat auto signé en attendant d’avoir un vrai) https, appliquée sur à page html du logiciel. 
	- Cliquez sur le bouton "Avancé"
	- Puis tout en bas sur "Continuer vers 127.0.0.1:5000 (risqué)"
Pas de panique ce n’est pas du tout riqué… C’est juste une formalité.
Vous verrez s'afficher l'interface de l'application.
Allez, c'est parti !
Bonne analyse vitrail ou mosaique à vous.

Pour ajouter ou retirer des images de la galerie, il suffit d'ajouter ou de supprimer les fichiers correspondants dans le dossier `images/`.

### Auteur

Prudencio AYENAN — stagiaire L3 en Bio-Informatique, à la plateforme ImageUP, Université de Poitiers, Mai-Juillet 2026

