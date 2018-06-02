# Optimizaciones posibles

## Factorización LU

Para el cálculo de normales, hay que resolver muchos sistemas de ecuaciones que
utilizan la misma matriz A (formada por las 3 direcciones de luz elegidas al
comienzo del procedimiento), lo único que cambia es el vector indepedendiente.

Aprovechando esto, se puede factorizar A en LU, y resolver con L, para tener
complejidad O(n^2) para cada pixel, en vez de O(n^3) siempre.

## Usar la máscara para reducir los cálculos

Se puede usar la máscara para resolver sólo los sistemas de ecuaciones que
tienen sentido resolver.

No hace falta calcular las normales de los puntos por fuera de la máscara.
Para el caso del cálculo de z, sólo buscar los zs dentro de la máscara Esto
implicaría que la matriz pasaría a ser de 2N x N a 2m x m con m la cantidad de
pixeles dentro de la máscara.

## Tener en cuenta cuando se multiplican matrices simétricas

Se pueden hacer la mitad de los cálculos

## Hacer M matriz banda

Para el mapa de profundidad, armar la matriz M de ecuaciones de tal manera que
quede matriz banda. Las matrices banda son más baratas en espacio y en cantidad
de cálculos.


# Comentarios varios

## Lectura de imágenes PPM

* Imagenes de escala de grises o RGB
* Profundidad de color: 8bit, 16bit, 32bit

## Pasaje de RGB a escala de grises (3 canales a 1 canal)

* Promedio
* Promedio ponderado (buscar referencias)

## Bordes de la matriz M de profundidad

Podemos poner 0 en los casos borde, o tomar otros vectores ortogonales, como
z_{x-1, y} y z_{x, y-1}...

## Procedimiento

### Calibración

```
./calibrate [--mean] IMAGE [IMAGE..]
```

donde `IMAGE` es cada una de las imágenes de la esfera de calibración, junto
con su máscara.  Se espera que una de las imágenes sea la máscara, y que su
nombre de archivo contenga la cadena `mask`.

Por ejemplo, si el directorio contiene 12 archivos (además de la máscara),
devuelve:

```
12
0.403259 0.480808 0.778592
0.0982272 0.163712 0.981606
-0.0654826 0.180077 0.98147
-0.127999 0.431998 0.892745
-0.328606 0.485085 0.810377
-0.110339 0.53593 0.837021
0.239071 0.41439 0.878138
0.0642302 0.417497 0.906406
0.12931 0.339438 0.931698
0.0323953 0.340151 0.939813
0.0985318 0.0492659 0.993914
-0.16119 0.354617 0.921013
```

La primer linea es la cantidad de imagenes analizadas, y las siguientes lineas
contienen los vectores de dirección de luz de cada imagen del conjunto, un
vector linea, y cada componente separada por un espacio.

### Reconstrucción 3D

```
./build mask.ppm img1.ppm img2.ppm img3.ppm < luces
```

donde `luces` es una lista de vectores de dirección de luz correspondientes a
las imágenes elegidas del objeto a reconstruir. Cada vector está separado por
el separador de nueva línea, y cada componente por un espacio. Ej:

```
-0.127999 0.431998 0.892745
-0.328606 0.485085 0.810377
-0.110339 0.53593 0.837021
```


### TODO
  * Al tomar pixeles de los bordes, chequear como cambia la ecuación, probablemente no sea -Nx, etc.
