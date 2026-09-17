# Instabilités de Rayleigh-Bénard

Ce projet implémente numériquement les **équations de Navier–Stokes incompressibles en 2D couplées à l'équation de la chaleur** afin de simuler la **convection naturelle de Rayleigh-Bénard**.

Le calcul est réalisé sur un maillage cartésien uniforme et utilise une méthode de **projection** :

1. Calcul de la température ;
2. Prédiction de la vitesse avec prise en compte de la poussée d'Archimède ;
3. Résolution de l'équation de Poisson pour la pression avec une méthode de **Jacobi** ;
4. Correction de la vitesse afin de satisfaire la condition d'incompressibilité.

La direction `x` est **périodique** et les températures sont imposées en haut et en bas du domaine.

---

## Fonctionnalités

Le code permet de simuler les **cellules de convection de Bénard** et de calculer automatiquement le **nombre de Rayleigh**.

Le pas de temps est automatiquement déterminé à partir des contraintes :

* CFL pour la convection ;
* Fourier pour la diffusion thermique et visqueuse.

La température initiale est composée d'un **profil linéaire perturbé par un bruit aléatoire** afin de favoriser l'apparition des instabilités.

Les résultats sont exportés au format `.dat`, compatible notamment avec **Tecplot** et les logiciels de visualisation scientifique.

---

## Structure du projet

```text
.
├── instabilites_rayleigh-benard.f95    # Code Fortran principal
├── input.dat                  # Paramètres de simulation
├── irb3*_tecplot.dat          # Fichiers résultats
└── README.md                  # Documentation
```


## Compilation et exécution

Le code peut être compilé avec `gfortran` :

```bash
gfortran projet_scientifique.f95 -o irb3
```

Puis exécuté avec :

```bash
./irb3
```

La simulation est lancée avec les paramètres définis dans `input.dat`. Il faut ainsi augmenter le nombre de Rayleigh pour observer l'apparition des cellules convectives et inversement. Pour cela il faut directement modifier les paramètres dans le `input.dat`, avec comme rappel :


Les fichiers de résultats sont générés automatiquement toutes les `20` itérations.

Ils contiennent notamment :

* les coordonnées `X` et `Y` ;
* les composantes de vitesse `u` et `v` ;
* la température `Temp` ;
* la pression `p`.

Ces données permettent de visualiser l'évolution de la **convection naturelle et des cellules de Bénard**.
