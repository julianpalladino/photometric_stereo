#!/bin/bash

./3dbuild res/mate/ 0 5 10 < experimentacion/calibracion/luces_nuestras_mate_50.txt && ./script/plot_normals mate.n{x,y,z}.txt experimentacion/calibracion/mate_percentil_0.pdf

./3dbuild res/mate/ 0 5 11 < experimentacion/calibracion/luces_nuestras_mate_50.txt && ./script/plot_normals mate.n{x,y,z}.txt experimentacion/calibracion/mate_percentil_25.pdf

./3dbuild res/mate/ 3 8 11 < experimentacion/calibracion/luces_nuestras_mate_50.txt && ./script/plot_normals mate.n{x,y,z}.txt experimentacion/calibracion/mate_percentil_50.pdf

./3dbuild res/mate/ 0 7 11 < experimentacion/calibracion/luces_nuestras_mate_50.txt && ./script/plot_normals mate.n{x,y,z}.txt experimentacion/calibracion/mate_percentil_75.pdf

./3dbuild res/mate/ 3 4 7 < experimentacion/calibracion/luces_nuestras_mate_50.txt && ./script/plot_normals mate.n{x,y,z}.txt experimentacion/calibracion/mate_percentil_100.pdf

./3dbuild res/mate/ 0 4 10 < res/luces.txt && ./script/plot_normals mate.n{x,y,z}.txt experimentacion/calibracion/mate_catedra.pdf
