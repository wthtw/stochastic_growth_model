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
 >Spécifier d’abord l'option 'type'
 >`optset('ncpsolve','type','minmax')`

## Problèmes d’optimisation (recherche du maximum et du minimum)
