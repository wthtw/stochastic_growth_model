# Notes de cours - Méthodes numeriques  

## Système d’équations linéaire : Ax = b
* Méthodes directes () 
  > Décomposition LU, Cholesky ` x = A\b `
    
* Méthode itérative (forme Qx = b + (Q − A)x)

  >Gauss-jacobi : `gjacobi (A,b)`
  >Gauss-seidel : `gseidel (A,b)`

## Systèmes d’équations non linéaires : Points fixes, Solutions racines
* Méthode bisection, sur un interval [a,b]

  > Pour une fonction f, `bisect (‘f’,a,b)`

* Méthode Newton : un ou plusieurs variables, avec des valeurs initiales , utilise le jacobien

  > Pour une fonction f à 2 variables, avec des valeurs initiales respectives x1 et x2 `newton(’f’,[x1;x2])`

* Méthode Quasi-Newton : utilise une approximation du jacobien

  * Secant Method : une variable

  * Broyden Method : plusieurs variables, utilise une valeur initiale pour la racine, et une autre pour le Jacobien
    > Pour une fonction f à deux variables, et pour les valeurs initiales x1et x2 des variables
    > `broyden(’f’,[x1;x2])`
    > **Note :** Pour ces méthodes, on peut ajouter une backstepping routine, pour éviter les divergences


