# Méthodes numeriques avec MATLAB - Notes de cours

## Système d’équations linéaire : Ax = b
### Méthodes directes
Décomposition LU, Cholesky ` x = A\b `
    
### Méthode itérative (forme Qx = b + (Q − A)x)

Gauss-jacobi : `gjacobi (A,b)`
  
Gauss-seidel : `gseidel (A,b)`

## Systèmes d’équations non linéaires : Points fixes, Solutions racines
### Méthode bisection, sur un interval [a,b]

Pour une fonction f, `bisect (‘f’,a,b)`

### Méthode Newton : un ou plusieurs variables, avec des valeurs initiales , utilise le jacobien

Pour une fonction f à 2 variables, avec des valeurs initiales respectives x1 et x2 `newton(’f’,[x1;x2])`

### Méthode Quasi-Newton : utilise une approximation du jacobien

 * Secant Method : une variable

 * Broyden Method : plusieurs variables, utilise une valeur initiale pour la racine, et une autre pour le Jacobien
 Pour une fonction f à deux variables, et pour les valeurs initiales x1et x2 des variables `broyden(’f’,[x1;x2])`
    
 >**Note :** Pour ces méthodes, on peut ajouter une backstepping routine, pour éviter les divergences

### Méthodes exclusives pour Point-fixes
 * Méthode Itération de fonction, pour une valeur initiale x0
 Pour une fonction g, 
 
 `fixpoint(’g’,x0)`

 * Complementary Method : utilise le jacobien
 Pour résoudre f(x) = 0, pour *8a < x < b* ;  a et b peuvent être Inf

 * Méthode semismooth
 Pour une fonction f, un intervalle [a,b], et une valeur initiale x0, 
 
 `ncpsolve(’f’,a,b,x)`
 
 * Méthode minmax
 Spécifier d’abord l'option 'type' 
 
 `optset('ncpsolve','type','minmax')`

## Problèmes d’optimisation (recherche du maximum et du minimum)
### Méthodes sans dérivées :
 * Method Golden Search : une variable, optimum local sur un interval [a,b]
 Pour une fonction f, sur un interval [a,b] 
 
 `golden(’f’,a,b)`

 * Méthode Algorithm Nelder-Mead : plusieurs variables, avec des valeurs initiales pour les variables
 Pour une fonction f à deux variables, avec les valeurs initiales x1 et x2, 
 
 `neldmead(’f’,[x1;x2])`

## Méthode d’intégration et de différentiation
### Méthode d’intégration

#### Calcul de l'aire
* Méthodes Newton-cotes : calcul de l’aire sous la fonction
**Trapezoid rule :** pour les fonctions discontinues ayant des points d’inflexion
 
 Pour n trapezes, sur un intervalle [a,b], *n* les nodes et *w* les weights,
 
 `[x,w] = qnwtrap(n,a,b)`
 
**Simpson rule :** pour les fonctions plus lisses
 
 `[x,w] = qnwsimp(n,a,b)`
 
 >Si *w(x)=1*, on calcule l’aire sous la fonction

 * Méthodes Gaussian quadrature
  Legendre quadrature, pour w(x) = 1 

 `[x,w] = qnwlege(n,a,b)`

#### Calcul de l’espérance
* Méthodes Gaussian quadrature
 Pour x suivant une **loi normale (mu, var)**, *n* les nodes gaussiens et *w* les weights gaussiens,
 
 `[x,w] = qnwnorm(n,mu, var)`
 
 >Si w(x) = *fonction de densité de probabilité* de *x*, on calcule l’espérance de la fonction par `Somme(w*f(x))`

* Méthodes Intégration Monte-Carlo
 Il faut générer pseudoaléatoirement *n* nodes *x* d’après la distribution ; les weights *w=1/n* étant identiques

 L’espérance de *f* est obtenue par `Somme(w*f(x))`

* Méthodes Quasi-Monte Carlo
 Ici les *n* nodes *x* sont déterministes, sur une intervalle [a,b] ;  les weights *w=(b-a)/n* étant identiques

 L’espérance de *f* est obtenue par `Somme(w*f(x))`
 
### Méthode de différentiation
Pour une fonction *f*, une approximation *O(h²)* de sa dérivée, autour du point *x0*, est de la forme:
>*f’(x0) = a*f(x0) + b*f(x0 + h) + c*f(x0 + alpha*h) [ + O(h²)]*
 
>où (x0 + h), (x0 + alpha*h) sont deux autres points,
>choisir, alpha = 2, et choisir h quelconque, les paramètres a, b et c, s’obtiennent en résolvant le système suivant :
```
a + b + c = 0
b + cλ = 1/h
b + cλ2 = 0
```

## Initial value problem (IVP)
Prend la forme d’une equation différentielle à résoudre pour une function solution, connaissant les valeurs initiales pour la fonction.

>e.g. à l’ordre 1  *y'(t) = f(t,y(t))*

>ou à l’ordre 2 *y''(t) = f(t,y(t),y'(t))* 

Il faut réécrire parfois le problème sous la forme d’une équation différentielle. Une équation différentielle d’ordre 2 peut se ramener à un système d’équation différentielle d’ordre 1.

Le principe de résolution est de dériver une approximation de taylor de  à l’ordre 1 (Méthode de Euler), ou 2 (Méthode Runge-Kutta 2), ou 4 (Runge-Kutta 4). On choisit un pas h pour subdiviser la période de temps ; plus petit h, mieux c’est en général.

>Méthode de Euler : *y(t+h) = y(t) + y'(t)h*

>Méthode Runge-Kutta 2 : *y(t+h) = y(t) + h[a_1k_1 + a_2k_2]*







