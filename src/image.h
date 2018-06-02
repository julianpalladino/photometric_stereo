#pragma once
#include "common.h"
#include "matrix.h"

typedef struct {
    int top;
    int right;
    int bottom;
    int left;
} BoundingBox;

// Lee una imagen RGB de tipo PPM, la pasa a escala de grises y devuelve la matriz.
// Para pasar a escala de grises, usa el promedio ponderado BT.709
//
// Flags de entorno:
// * USE_AVG=1: Usa media entre los 3 canales RGB en vez del ponderado
// * SAVE=1: Guarda la imagen resultante en path + "_gray.ppm"
//
Matrix<uchar> read_image_gray(const string& path);

// Calcula el bounding box mínimo para +img+ usando la máscara +mask_img+
//
BoundingBox bounding_box(const Matrix<uchar>& img, Matrix<uchar> mask_img);

// Recorta una imagen en base al bounding box de su máscara
//
Matrix<uchar> trim(const Matrix<uchar>& img, const Matrix<uchar>& mask_img);
