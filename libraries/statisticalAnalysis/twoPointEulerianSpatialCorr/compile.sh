#!/bin/bash

g++ -o $CFDEM_APP_DIR/twoPointEulerianSpatialCorr twoPointEulerianSpatialCorr.cpp -I../../ann/include/ANN -L../../ann/lib -lANN
