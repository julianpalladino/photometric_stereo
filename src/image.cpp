#include "image.h"
#include "ppmloader.h"

using namespace std;

static void read_image(string filename, uchar** data, int* width, int* height);
static uchar rgb_to_gray(uchar* data, int i, int j, int height, int width, bool use_avg);

// Lee una imagen RGB de tipo PPM, la pasa a escala de grises y devuelve la
// matriz.
//
Matrix<uchar> read_image_gray(const string& path) {
    uchar* data = NULL;
    int height = 0, width = 0;

    // Flags para experimentación
    bool use_avg = get_env_var("AVG") == "1";
    bool save = get_env_var("SAVE") == "1";

    read_image(path, &data, &width, &height);
    Matrix<uchar> image(height, width, 0);

    for (int h = 0; h < height; ++h) {
        for (int w = 0; w < width; ++w) {
            uchar value = rgb_to_gray(data, h, w, height, width, use_avg);
            image[h][w] = value;

            if (save) {
                uchar pos = h*width*3 + w*3;
                data[pos + 0] = data[pos + 1] = data[pos + 2] = value;
            }
        }
    }

    if (save) {
        char comments[100];
        //sprintf(comments, "%s", "Hello world");

        bool ret = SavePPMFile(path.c_str(), data, width, height,
            PPM_LOADER_PIXEL_TYPE_RGB_8B, comments);
        if (!ret) {
            throw runtime_error("No pudo guardar la imagen: " + path);
        }
    }

    delete data;

    return image;
}

// Lee la imagen PPM +filename+ y carga cada canal en el parámetro +data+,
// junto con su tamaño en +width+ y +height+
//
static void read_image(string filename, uchar** data, int* width, int* height) {
    *data = NULL;
    *width = 0;
    *height = 0;
    PPM_LOADER_PIXEL_TYPE pt = PPM_LOADER_PIXEL_TYPE_INVALID;

    bool ret = LoadPPMFile(data, width, height, &pt, filename.c_str());
    if (!ret || width == 0|| height == 0|| pt!=PPM_LOADER_PIXEL_TYPE_RGB_8B){
        throw runtime_error("Fallo al leer la imagen " + filename);
    }
}

BoundingBox bounding_box(const Matrix<uchar>& mask_img) {
    int top    = mask_img.height();
    int right  = 0;
    int bottom = 0;
    int left   = mask_img.width();

    for (int h = 0; h < mask_img.height(); h++) {
        for (int w = 0; w < mask_img.width(); w++) {
            if (mask_img[h][w] != 0) {
                top    = (h < top)    ? h : top;
                right  = (w > right)  ? w : right;
                bottom = (h > bottom) ? h : bottom;
                left   = (w < left)   ? w : left;
            }
        }
    }

    return {top, right, bottom, left};
}

Matrix<uchar> trim(const Matrix<uchar>& img, const Matrix<uchar>& mask_img) {
    BoundingBox bbox = bounding_box(mask_img);
    return img.slice(bbox.top, bbox.bottom, bbox.left, bbox.right);
}

// Devuelve el valor de gris para un pixel determinado
//
// Si +use_avg+ es true, se promedian los 3 valores RGB de manera equitativa.
// Si no, se hace un promeido ponderado basado en cómo el ojo humano percibe
// colores.
//
static uchar rgb_to_gray(uchar* data, int i, int j, int height, int width, bool use_avg) {
    if (i > height || i < 0) {
        throw runtime_error("Posición i inválida");
    }

    if (j > width || j < 0) {
        throw runtime_error("Posición j inválida");
    }

    uchar red   = data[i*width*3 + j*3 + 0];
    uchar green = data[i*width*3 + j*3 + 1];
    uchar blue  = data[i*width*3 + j*3 + 2];

    uchar value;

    if (use_avg) {
      // Promedia los valores de los 3 canales por igual
      value = (red + green + blue) / 3;
    } else {
      // Luma: Promedio ponderado teniendo en cuenta la percepción del ojo
      // humano. Del estándar BT.709 de ITU-R:
      // https://en.wikipedia.org/wiki/Luma_%28video%29
      value = red * 0.2126 + green * 0.7152 + blue * 0.0722;
    }

    return value;
}
