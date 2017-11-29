#!/bin/bash

g++ -o $CFDEM_APP_DIR/twoPointSpatialCorr twoPointSpatialCorr.cpp -I../../ann/include/ANN -L../../ann/lib -lANN
