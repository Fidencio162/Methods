# Methods
Solución numérica de integrales con el método del trapecio, Simpson (1/3, 3/8)

Para compilar se hace de la siguiente manera: 

• gfortran -c Mathematical.f, lo que da dos archivos: Mathematical.mod y Mathematical.o. Podemos incluir el archivo .o en una librería de programas o bien guardarlo aparte.

• Una vez construido nuestro programa main compilamos y ejecutamos de la siguiente manera: gfortran -o main -J. main.f Mathematical.o.
