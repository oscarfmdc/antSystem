#!/bin/bash
# Resuelve el problema kroAB100 multiobjetivo en 100 segundos
 cd ./moaco/trunk
 ~/bin/moaco_btsp --input ../../btsp/kroAB100.tsp --time 100 --BicriterionAnt > ../../results
