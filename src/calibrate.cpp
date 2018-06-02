#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>
#include "image.h"
#include "matrix.h"

using namespace std;

typedef struct {
    int x;
    int y;
} Point;

const int DEFAULT_WINDOW_SIZE = 5;

void log(const char* fmt, ...) {
    va_list args;

    fprintf(stderr, "[calibrate] ");

    va_start(args, fmt);
    vfprintf(stderr, fmt, args);
    va_end(args);

    fprintf(stderr, "\n");
}

tuple<Point, double> center_and_radius(const Matrix<uchar>& img) {
    int up = img.width();
    int right = img.height();
    int low = 0;
    int left = 0;

    for (int h = 0; h < img.height(); h++) {
        for (int w = 0; w < img.width(); w++) {
            if (img[h][w] != 0) {
                up = (h<up) ? h : up;
                low = (h>low) ? h : low;
                left = (w>left) ? w : left;
                right = (w<right) ? w : right;
            }
        }
    }

    double radius = (left - right) / 2;
    int x = ((right - left) / 2) + left;
    int y = ((low - up) / 2) + up;

    return make_tuple(Point{.x = x, .y = y}, radius);
}

unsigned long sum_window_values(const Matrix<uchar>& img, int row, int col,
        int window_size)
{
    unsigned long sum = 0;

    for (int i = 0; i < window_size; i++) {
        for (int j = 0; j < window_size; j++) {
            sum += img[row + i][col + j];
        }
    }

    return sum;
}

// Calcula el punto de mayor brillo en la imagen, dentro de una ventana de
// tamaño +window_size+.
//
// Busca la ventana tal que la sumatoria de todos sus valores sea máxima en
// toda la imagen, y retorna el punto medio en la ventana.
//
Point most_brightest_point(const Matrix<uchar>& img, int window_size = 5) {
    int max_px = 0, max_py = 0;
    unsigned long max_value = 0;

    for (int i = 0; i < img.height() - window_size; i++) {
        for (int j = 0; j < img.width() - window_size; j++) {
            unsigned long value = sum_window_values(img, i, j, window_size);

            if (value > max_value) {
                max_px = j;
                max_py = i;
                max_value = value;
            }
        }
    }

    return {
        max_px + (window_size / 2),
        max_py + (window_size / 2)
    };
}


int main(int argc, char** argv) {
    if (argc != 2) {
        cerr << "Uso: " << argv[0] << " RUTA\n" \
             << "\n" \
             << "RUTA es la ruta al directorio que contiene las imágenes\n" \
             << "de calibración.\n" \
             << "\n" \
             << "Flags:\n" \
             << "\tAVG:\tUsa promedio sobre RGB (default: 0)\n" \
             << "\tWINDOW:\tTamaño de la ventana (default: 5)\n" \
             << "\tSAVE:\tGuarda imágenes en escala de grises (default: 0)\n" \
             << "\n" \
             << "Ejemplo:\n\tAVG=1 WINDOW=8 " << argv[0] << " res/mate\n" \
             << endl;

        return 1;
    }

    string root = string(argv[1]);
    string basename = get_basename(root);

    // Imprime en pantalla las rutas de las imágenes
    string mask_path = root + "/" + basename + ".mask.ppm";
    log("Máscara: %s", mask_path.c_str());

    vector<string> image_paths;
    image_paths.reserve(12);
    for (int i = 0; i < 12; i++) {
        image_paths.push_back(root + "/" + basename + "." + to_string(i) + ".ppm");
        log("Imagen %d: %s", i, image_paths[i].c_str());
    }

    // Imprime el uso de los flags
    if (get_env_var("AVG") == "1") {
        log("Usa promedio para escala de grises");
    }
    if (get_env_var("SAVE") == "1") {
        log("Guarda archivos temporales para escala de grises");
    }

    int window_size = DEFAULT_WINDOW_SIZE;
    if (!get_env_var("WINDOW").empty()) {
        window_size = stoi(get_env_var("WINDOW"));
    }
    log("Tamaño de la ventana: %d", window_size);

    // Lee imagen de la máscara
    auto mask = read_image_gray(mask_path);

    // Calcula la posición del centro de la esfera en la máscara
    auto center_radius = center_and_radius(mask);
    Point center = get<0>(center_radius);
    double r = get<1>(center_radius);

    log("Radio: %f", r);
    log("Centro: (%d, %d)", center.x, center.y);

    // Imprime la cantidad de imagenes
    cout << image_paths.size() << endl;

    // Para cada imagen, busca el punto de mayor intensidad en una vecindad,
    // y calcula la normal a la superficie en ese punto. La dirección de luz
    // es la normal en dirección opuesta.
    for (int i = 0; i < image_paths.size(); i++) {
        auto path = image_paths[i];
        auto img = read_image_gray(path);

        Point maxP = most_brightest_point(img, window_size);

        log("Máximo %d: (%d, %d)", i, maxP.x, maxP.y);

        // Calcula la normal con el punto más brillante y el centro de la esfera
        double nx, ny, nz;
        nx = -(maxP.x - center.x);
        ny = maxP.y - center.y;
        nz = -sqrt(r*r - nx*nx - ny*ny);

        // Normaliza el vector
        nx /= r; ny /= r; nz /= r;

        cout << -nx << " "
             << -ny << " "
             << -nz << endl;
    }

    return 0;
}
