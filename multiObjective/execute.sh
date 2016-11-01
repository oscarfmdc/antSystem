#!/bin/bash
# Resuelve el problema kroAB100 multiobjetivo en 100 segundos
 cd ./moaco/trunk
 
#~/bin/moaco_btsp --input ../../btsp/kroAB100.tsp --time 100 --BicriterionAnt > ../../resultsBicriterion
#~/bin/moaco_btsp --input ../../btsp/kroAB100.tsp --time 100 --PACO > ../../resultsPACO

~/bin/moaco_btsp --input ../../btsp/euclidAB100.tsp --time 200 -m 100 -b 5 --BicriterionAnt > ../../resultsBicriterion
~/bin/moaco_btsp --input ../../btsp/euclidAB100.tsp --time 200 -m 100 -b 5 --PACO > ../../resultsPACO
~/bin/moaco_btsp --input ../../btsp/euclidAB100.tsp --time 200 -m 100 -b 5 --MONACO > ../../resultsMONACO