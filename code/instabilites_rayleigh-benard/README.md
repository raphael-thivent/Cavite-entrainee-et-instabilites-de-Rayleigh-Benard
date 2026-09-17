#  Instabilités de Rayleigh-Bénard

Ce projet implémente numériquement les **équations de Navier–Stokes incompressibles en 2D** pour simuler l'écoulement dans une **cavité entraînée** (*lid-driven cavity*).

Le calcul est réalisé sur un maillage cartésien uniforme et utilise une méthode de **projection** :

1. Prédiction de la vitesse ;
2. Résolution de l'équation de Poisson pour la pression ;
3. Correction de la vitesse afin de satisfaire la condition d'incompressibilité.

Le cas de référence utilisé est notamment celui d'une cavité à **Re = 100**.

---

## Fonctionnalités

Le code permet de choisir entre trois schémas de discrétisation pour le terme convectif :

- **Upwind** — ordre 1
- **Différences centrées** — ordre 2
- **Différences centrées** — ordre 4

La résolution de l'équation de Poisson est effectuée avec une méthode itérative de **Jacobi**.

Le pas de temps est automatiquement déterminé à partir des contraintes :

- CFL pour la convection ;
- Fourier pour la diffusion.

Les résultats sont exportés au format `.dat`, compatible notamment avec **Tecplot**.

---

## Structure du projet

```text
.
├── cavite_entrainee.f95       # Code Fortran principal
├── input.dat        # Paramètres de simulation
├── *tecplot.dat     # Fichiers résultats 
└── README.md        # Documentation
```

## Fonctionnement du code

Les paramètres de la simulation peuvent être facilement modifiés dans le fichier `input.dat`.

Il est notamment possible de modifier :

- le nombre de Reynolds `Re` ;
- le temps de simulation `t_max` ;
- les conditions aux limites `BCE` et `BCL` ;
- les critères de convergence ;
- les coefficients CFL et Fourier ;
- la fréquence d'écriture des résultats ;
- le schéma de discrétisation utilisé.

Le choix du schéma se fait avec `type_discretisation` :

```text
1 → Upwind
2 → Différences centrées d'ordre 2
4 → Différences centrées d'ordre 4
```

## Compilation et exécution

Le code doit d'abord être compilé avec un compilateur Fortran, par exemple avec `gfortran` :

```bash
gfortran cavite_entrainee.f95
```

Cette commande crée un exécutable nommé a.out.

Il suffit ensuite de l'exécuter :
```bash
./a.out
```
La simulation est alors lancée avec les paramètres définis dans input.dat.

À intervalles réguliers, le programme génère automatiquement des fichiers de résultats *tecplot.dat. Ces fichiers peuvent être ouverts avec ParaView pour visualiser les différents champs de la simulation.

Les résultats contiennent notamment :

    - le champ de vitesse ;

    - les composantes de vitesse u et v ;

    - le champ de pression p ;

    - les coordonnées X et Y.

Cela permet notamment de visualiser les champs de vitesse et de pression dans la cavité et d'analyser l'écoulement obtenu.


