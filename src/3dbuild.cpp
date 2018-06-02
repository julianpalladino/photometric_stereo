#include <cassert>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <sstream>
#include <string>
#include <vector>
#include "common.h"
#include "matrix.h"
#include "ppmloader.h"
#include "image.h"
#include "algoritmosmatriciales.h"

using namespace std;

void log(const char* fmt, ...) {
    va_list args;

    fprintf(stderr, "[3dbuild] ");

    va_start(args, fmt);
    vfprintf(stderr, fmt, args);
    va_end(args);

    fprintf(stderr, "\n");
}

vector<double> parse_vector(const string line) {
    vector<double> res;

    string token;
    stringstream ss(line);

    while (getline(ss, token, ' ')) {
        if (!token.empty()) {
            res.push_back(stod(token));
        }
    }

    return res;
}

void save_normals_files(const string& basename, const Matrix<normal>& normals) {
    ofstream nx_f, ny_f, nz_f;

    nx_f.open(basename + ".nx.txt");
    ny_f.open(basename + ".ny.txt");
    nz_f.open(basename + ".nz.txt");

    nx_f.precision(6); ny_f.precision(6); nz_f.precision(6);

    for (int i = 0; i < normals.height(); i++) {
        for (int j = 0; j < normals.width(); j++) {
            auto n = normals(i, j);
            nx_f << n.x << ","; ny_f << n.y << ","; nz_f << n.z << ",";
        }
        nx_f << "\n"; ny_f << "\n"; nz_f << "\n";
    }

    nx_f << endl; ny_f << endl; nz_f << endl;
    nx_f.close(); ny_f.close(); nz_f.close();

    log("%s.{nx,ny,nz}.txt guardados", basename.c_str());
}

void save_depth_files(const string& basename, const Matrix<double>& depths) {
    ofstream z_f;

    z_f.open(basename + ".z.txt");
    z_f.precision(6);

    for (int i = 0; i < depths.height(); i++) {
        for (int j = 0; j < depths.width(); j++) {
            double z = depths(i, j);
            z_f << z << ",";
        }
        z_f << "\n";
    }

    z_f << endl;
    z_f.close();

    log("%s.z.txt guardados", basename.c_str());
}

int main(int argc, char** argv) {
    if (argc != 5) {
        cerr << "Uso: " << argv[0] << " RUTA I1 I2 I3 < LUCES\n" \
             << "\n" \
             << "RUTA es la ruta al directorio que contiene las imágenes.\n" \
             << "I1, I2, I3 son los índices de imágenes a usar en el cálculo\n" \
             << "de las normales.\n"
             << "LUCES es la lista de las 12 direcciones de luz obtenidas\n" \
             << "en la calibración.\n" \
             << "\n" \
             << "Flags:\n" \
             << "\tAVG:\tUsa promedio sobre RGB (default: 0)\n" \
             << "\tZ_SIST:\tMétodo para sistema del z: [EG, CL] (default: CL)\n" \
             << "\n" \
             << "Ejemplo:\n\tAVG=1 " << argv[0] << " res/buho 0 4 6 < res/luces.txt\n" \
             << endl;
        return 1;
    }

    // Parsea argumentos
    string root = string(argv[1]);
    string basename = get_basename(root);
    vector<int> idx(3);
    for (int i = 0; i < 3; i++) {
        idx[i] = atoi(argv[i+2]);
    }

    // Carga la máscara
    string mask_path = root + "/" + basename + ".mask.ppm";
    log("Máscara: %s", mask_path.c_str());
    auto mask_img = read_image_gray(mask_path);

    // Carga las 3 imágenes y las recorta usando la máscara
    vector<Matrix<uchar> > imgs;
    for (int i = 0; i < 3; i++) {
        string img_path = root + "/" + basename + "." + to_string(idx[i]) + ".ppm";
        auto img = read_image_gray(img_path);

        if (!get_env_var("TEST_TRIM").empty()) {
            int size = stoi(get_env_var("TEST_TRIM"));
            int half_size = size / 2;
            auto new_img = img.slice(
                    img.height() / 2 - half_size,
                    img.height() / 2 + half_size,
                    img.width() / 2 - half_size,
                    img.width() / 2 + half_size);
            log("[TEST] Tamaño original: %dx%d, imagen recortada: %dx%d",
                    img.width(), img.height(),
                    new_img.width(), new_img.height());
            log("[TEST] x = %d, y = %d", img.height() / 2 - half_size,
                    img.width() / 2 - half_size);
            img = new_img;
            auto new_mask_img = mask_img.slice(
                    mask_img.height() / 2 - half_size,
                    mask_img.height() / 2 + half_size,
                    mask_img.width() / 2 - half_size,
                    mask_img.width() / 2 + half_size);
            mask_img = new_mask_img;
        }

        imgs.push_back(img);
        log("Imagen %d: %s", idx[i], img_path.c_str());
    }

    // Lee stdin
    vector<string> lines;
    for (string line; getline(cin, line);) {
        lines.push_back(line);
    }

    // Parsea los 3 vectores de dirección de luz
    Matrix<double> lights(3, 3, 0);
    for (int i = 0; i < 3; i++) {
        auto s = parse_vector(lines[idx[i] + 1]);
        for (int j = 0; j < 3; j++) {
            lights[i][j] = s[j];
        }
        log("Luz %d: %f %f %f", idx[i], s[0], s[1], s[2]);
    }

    // Cálculo del mapa de normales
    auto normals = obtenerPlanoNormal(lights, imgs);
    log("Mapa normal construido: %dx%d", normals.height(), normals.width());
    save_normals_files(basename, normals);

    // Estimación de la profundidad
    //
    // Construir la matriz M y el vector v
    // Hacer A = M'*M y b = M'*v
    // Resolver sistema Az = b con EG y Cholesky (configurable)
    // Salida esperada: profundidades (z), una línea por cada fila de pixeles de la imagen
    //   y en cada linea, cada z separado por coma.
    auto depths = obtener_profundidades(normals, mask_img);
    save_depth_files(basename, depths);

    return 0;
}
