# Notes de cours - Méthodes numeriques  

## Système d’équations linéaire : Ax = b
* Méthodes directes () 
    > Décomposition LU, Cholesky ` x = A\b `
    
* Méthode itérative (forme Qx = b + (Q − A)x)

    >Gauss-jacobi : `gjacobi (A,b)`
    >Gauss-seidel : `gseidel (A,b)`

## Systèmes d’équation non linéaire : Points fixes : x0 tel que f(x0) = x0 / Solutions racines : x0 tel que f(x0) = 0
* Méthode bisection, sur un interval [a,b]

    > Pour une fonction f, `bisect (‘f’,a,b)`

## Méthode Newton : un ou plusieurs variables, avec des valeurs initiales , utilise le jacobien

    > Pour une fonction f à 2 variables, avec des valeurs initiales respectives x1 et x2 `newton(’f’,[x1;x2])`



