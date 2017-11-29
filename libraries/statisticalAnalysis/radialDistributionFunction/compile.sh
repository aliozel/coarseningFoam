#!/bin/bash

g++ -o $CFDEM_APP_DIR/radialDistributionFunction radialDistributionFunction.cpp -I../../ann/include/ANN -L../../ann/lib -lANN
