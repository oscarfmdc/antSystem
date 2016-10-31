#!/bin/bash
# Recopilacion de ejecuciones

# http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/

# m : numero de hormigas
# a : alfa (influencia de la feromona)
# b : beta (influencia de la distancia)
# t: tiempo
# o: valor optimo
# e: Evapoarcion de feromona
# q: Probabilidad de mejor eleccion
# l: Busqueda local -  0: no local search   1: 2-opt   2: 2.5-opt   3: 3-opt
# r : intentos

rm ./results/*
cd ./results

# Problema a280
../src/acotsp -i ../problems/a280.tsp -o 2579 -t 60 -m 280 -b 5

# Problema att48
#../src/acotsp -i ../problems/att48.tsp -o 10628 -t 60 -m 48 -b 5

# Problema berlin52
#../src/acotsp -i ../problems/berlin52.tsp -o 7542 -t 60 -m 52 -b 5

# Problema u1432
../src/acotsp -i ../problems/u1432.tsp -o 152970 -t 60 -m 1022 -b 5

# Problema fnl4461
#../src/acotsp -i ../problems/fnl4461.tsp -o 182566 -t 60 -m 1022 -b 5