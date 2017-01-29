# Notes de cours - Méthodes numeriques

## Système d’équations linéaire : Ax = b
#### Méthodes directes () 
>Décomposition LU, Cholesky ` x = A\b `
    
#### Méthode itérative (forme Qx = b + (Q − A)x)

>Gauss-jacobi : `gjacobi (A,b)`
  
>Gauss-seidel : `gseidel (A,b)`

## Systèmes d’équations non linéaires : Points fixes, Solutions racines
#### Méthode bisection, sur un interval [a,b]

>Pour une fonction f, `bisect (‘f’,a,b)`

#### Méthode Newton : un ou plusieurs variables, avec des valeurs initiales , utilise le jacobien

>Pour une fonction f à 2 variables, avec des valeurs initiales respectives x1 et x2 `newton(’f’,[x1;x2])`

#### Méthode Quasi-Newton : utilise une approximation du jacobien

 * Secant Method : une variable

 * Broyden Method : plusieurs variables, utilise une valeur initiale pour la racine, et une autre pour le Jacobien
 >Pour une fonction f à deux variables, et pour les valeurs initiales x1et x2 des variables `broyden(’f’,[x1;x2])`
    
 >**Note :** Pour ces méthodes, on peut ajouter une backstepping routine, pour éviter les divergences

#### Méthodes exclusives pour Point-fixes
 * Méthode Itération de fonction, pour une valeur initiale x0
 >Pour une fonction g, `fixpoint(’g’,x0)`

 * Complementary Method : utilise le jacobien
Pour résoudre f(x) = 0, pour a < x < b ;  a et b peuvent être Inf

 * Méthode semismooth
 >Pour une fonction f, un intervalle [a,b], et une valeur initiale x0, `ncpsolve(’f’,a,b,x)`
 
 * Méthode minmax
 >Spécifier d’abord l'option 'type' `optset('ncpsolve','type','minmax')`

## Problèmes d’optimisation (recherche du maximum et du minimum)
#### Méthodes sans dérivées :
* Method Golden Search : une variable, optimum local sur un interval [a,b]
>Pour une fonction f, sur un interval [a,b] `golden(’f’,a,b)`

* Méthode Algorithm Nelder-Mead : plusieurs variables, avec des valeurs initiales pour les variables
>Pour une fonction f à deux variables, avec les valeurs initiales x1 et x2, `neldmead(’f’,[x1;x2])`

## Méthode d’intégration et de différentiation
#### Méthode d’intégration

**Calcul de l'aire**
* Méthodes Newton-cotes : calcul de l’aire sous la fonction

**Trapezoid rule :** pour les fonctions discontinues et avec des points d’inflexion
Pour n trapezes, sur un intervalle [a,b], les nodes  et les weights w `[x,w] = qnwtrap(n,a,b)`

**Simpson rule :** pour les fonctions plus lisses `[x,w] = qnwsimp(n,a,b)`

>**Note :** Si w(x)=1, on calcule l’aire sous la fonction

* Méthodes Gaussian quadrature
>Legendre quadrature, pour w(x) = 1 `[x,w] = qnwlege(n,a,b)`

**Calcul de l’espérance**
* Méthodes Gaussian quadrature
Si w(x) = *fonction de densité de probabilité* de x, on calcule l’espérance de la fonction
Pour x suivant une **loi normale (mu, var)**, x les nodes gaussiens et w les weights gaussiens,

`[x,w] = qnwnorm(n,mu, var)`

L’espérance de la fonction est obtenue ensuite par `Somme(w*f(x))`

* Méthodes Intégration Monte-Carlo
Il faut générer pseudoaléatoirement n nodes x d’après la distribution ; les weights `w=1/n` étant identiques

L’espérance de f est obtenue par `Somme(w*f(x))`

* Méthodes Quasi-Monte Carlo
Ici les n nodes x sont déterministes, sur une intervalle [a,b] ;  les weights `w=(b-a)/n` étant identiques

L’espérance de f est obtenue par `Somme(w*f(x))`




