# Instabilités de Rayleigh-Bénard

Ce projet implémente numériquement les **équations de Navier–Stokes incompressibles en 2D couplées à l'équation de la chaleur** afin de simuler la **convection naturelle de Rayleigh-Bénard**.

Le calcul est réalisé sur un maillage cartésien uniforme et utilise une méthode de **projection** :

1. Calcul de la température ;
2. Prédiction de la vitesse ;
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


Le code doit d'abord être compilé avec un compilateur Fortran, par exemple avec :

```bash
f95 cavite_entrainee.f95
```

Cette commande crée un exécutable nommé a.out.

Il suffit ensuite de l'exécuter :
```bash
./a.out
```

La simulation est lancée avec les paramètres définis dans `input.dat`. Il faut ainsi augmenter le nombre de Rayleigh pour observer l'apparition des cellules convectives et inversement. Pour cela il faut directement modifier les paramètres dans le `input.dat`, avec comme rappel :
```math
Ra = \frac{gL^3\beta\Delta T}{\nu\alpha}
```
Les paramètres sont fixés de tel sorte que le `Ra` équivaut directement au `g`. Il suffit donc juste de modifier seulement le `g` pour changer le `Ra`.

Les fichiers de résultats `*tecplot.dat` sont générés automatiquement toutes les `20` itérations. Ces fichiers peuvent être ouverts avec ParaView pour visualiser les différents champs de la simulation.

Ils contiennent notamment :

* les coordonnées `X` et `Y` ;
* les composantes de vitesse `u` et `v` ;
* la température `Temp` ;
* la pression `p`.

Ces données permettent de visualiser l'évolution de la **convection naturelle et des cellules de Bénard**.
