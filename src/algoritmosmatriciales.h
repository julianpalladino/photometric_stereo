#pragma once
#include <vector>
#include <list>
#include <iostream>
#include <cassert>
#include "matrix.h"
#include "band_matrix.h"
#include "common.h"
#include "dok_matrix.h"

using namespace std;

typedef struct {
    DOKMatrix<double> coeficientes;
    Matrix<double> independiente;
} SistemaDeEcuaciones;

typedef struct {
    double x;
    double y;
    double z;
} normal;

ostream& operator<<(ostream& out, const normal& v) {
    out << "(" << v.x << ", " << v.y << ", " << v.z << ") ";
    return out;
}

bool eq_escalar(double a, double b){
    return (abs(a - b) < eps);
}

template <class T>
ostream& operator<<(ostream& out, const vector<T>& v) {
    for(int i = 0; i < v.size(); i++){
        out << v[i] << " ";
    }
    out << "\n";
    return out;
}

template <class T>
T vector_norm(vector<T>& v) {
    T acum = 0;
    for (int i = 0; i < v.size(); i++) {
        auto vi = v[i];
        acum += vi * vi;
    }
    return sqrt(acum);
}

void print_matrix(Matrix<double>&src){
    for (int i = 0; i < src.height(); ++i)
    {
        for (int j = 0; j < src.width(); ++j)
        {
            printf("%f ", src[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

template <class T, template<class U> class Mat>
void swap_rows(Mat<T>& src, int r1, int r2) {
    T tmp;
    for (int i = 0; i < src.height(); i++) {
        tmp = src(r1, i);
        src(r1, i) = src(r2, i);
        src(r2, i) = tmp;
    }
}

template <class T, template<class U> class Mat>
vector<int> triangular(Mat<T>& src, Mat<T>& ind) {
    const auto& src_r = src;

    auto perm = vector<int>(src_r.height());
    // inicializar permutaciones de 0..n-1
    for (int i = 0; i < src_r.height(); i++){
        perm[i] = i;
    }

    for (int i = 0; i < src_r.height() -1; i++) {
        if (src_r(i, i) == 0) {
            int j = i+1;
            while (j < src_r.height() && src(j, i) == 0) {
                j++;
            }
            // En la columna i son todos 0, no hay nada que hacer.
            if (j == src_r.height()) continue;

            // Tengo una fila que en la columna i no tiene un 0, swapeo la actual con esa.
            swap_rows(src, i, j);

            // Swapeo en el vector de permutaciones
            int tmp = perm[i];
            perm[i] = perm[j];
            perm[j] = tmp;
        }

        for (int j = i+1; j < src_r.height(); j++) {
            if (src_r(i, i) == 0) {
                if (ind(i, 0) == 0) {
                    throw runtime_error("Hay infinitas soluciones");
                }
                throw runtime_error("No hay solucion");
            }

            double coef = src_r(j, i) / src_r(i, i);
            for (int k = i; k < src_r.height(); k++) {
                auto v = src_r(j, k) - src_r(i, k) * coef;
                if (abs(v) > eps) {
                    src(j, k) = v;
                }
            }
            ind(j, 0) = ind(j, 0) - coef * ind(i, 0);
            // Escribo el coeficiente donde iría un 0.
            if (abs(coef) > eps) {
                src(j, i) = coef;
            }
        }
    }

    return perm;
}

template <class T, template<class U> class Mat>
void permutar(Mat<T>& aPermutar, const vector<int>& permutaciones) {
    const Mat<T>& aux = Mat<T>(aPermutar);
    for (int i = 0; i < permutaciones.size(); i++) {
        aPermutar(i, 0) = aux(permutaciones[i], 0);
    }
}

// calcula la suma de los elementos superiores en la matriz, usado en backwards sustitution
template <class T, template<class U> class Mat>
T suma_de_superiores(const Mat<T>& src, const Mat<T>& variables, int i) {
    T result = 0;
    for (int j = i + 1; j < src.height(); j++) {
        result += src(i, j) * variables(j, 0);
    }

    return result;
}

template <class T, template<class U> class Mat>
Mat<T> backward_subtitution(const Mat<T>& src, const Mat<T>& ind, const vector<int>& permutaciones){
    // La matriz ya es triangular superior. Solo falta despejar las variables.
    Mat<T> variables(ind.height(),1,0);
    T suma = T();
    for (int i = variables.height() -1; i >= 0; i--){
        suma = suma_de_superiores(src, variables, i);

        if (src(i, i) == 0) {
            // no tiene solucion
            variables(i, 0) = (ind(i, 0) - suma);
            
            // assert(ind(i, 0) - suma != 0);
        } else {
            variables(i, 0) = (ind(i, 0) - suma) / src(i, i);
        }
    }

    permutar(variables, permutaciones);

    return variables;
}

// calcula la suma de los elementos inferiores en la matriz, usado en forward sustitution
template <class T, template<class U> class Mat>
T suma_de_inferiores(const Mat<T>& src, const Mat<T>& variables, int i) {
    T result = 0;
    for (int j = 0; j < i; j++){
        result += src(i, j) * variables(j, 0);
    }
    return result;
}

template <class T, template<class U> class Mat>
Mat<T> forward_subtitution(const Mat<T>& src, const Mat<T>& ind, const vector<int>& permutaciones) {
    // La matriz ya es triangular inferior. Solo falta despejar las variables.
    Mat<T> variables(ind.height(),1,0);
    T suma = T();
    for (int i = 0; i < variables.height(); i++) {
        suma = suma_de_inferiores(src, variables, i);
        if (src(i, i) == 0) {
            // no tiene solucion
            variables(i, 0) = (ind(i, 0) - suma);
            
      //      assert(ind(i, 0) - suma != 0);
      
        } else {
            variables(i, 0) = (ind(i, 0) - suma) / src(i, i);
        }
    }

    //permutar(variables, permutaciones);

    return variables;
}

template <class T, template<class U> class Mat>
Mat<T> eliminacion_gaussiana(Mat<T>& src, Mat<T>& ind) {
    const vector<int> permutaciones = triangular(src, ind);

    return backward_subtitution(src, ind, permutaciones);
}

//bool es_nula(Matrix<T>& src){

//Recibe una matriz de nx1
template <class T>
double norma_dos(Matrix<T>& src){
  //cout << "calculo la norma de: " <<src << endl;
  double norma = 0;
  for(int i=0; i<src.height(); i++){
    norma += pow(src[i][0], 2);
  }
  return sqrt(norma);
}

vector<int> factorizacion_lu(Matrix<double>& src, Matrix<double>& l, Matrix<double>& u, Matrix<double>& ind){
  auto vectorPerm = triangular(src, ind);

  // en src nos queda L y U mezcladas, las separamos

  // armamos U, son los elementos de la diagonal para arriba, incluida la diagonal.
  for (int i = 0; i < src.height(); ++i){
    for (int j = i; j < src.width(); ++j){
      u[i][j] = src[i][j];
    }
  }

  // armamos L:
  // L tiene los elementos de src desde la diagonal para abajo, sin incluir la diagonal
  for (int i = 0; i < src.height(); ++i){
    for (int j = 0; j < i; ++j){
      l[i][j] = src[i][j];
    }
    l[i][i] = 1;
  }

  return vectorPerm;
}

// Dados L, U y b, calcula x tal que LUx = b
template <class T, template<class U> class Mat>
Mat<T> resolverSistemaLU(Mat<T>& L, Mat<T>& U, vector<int>& P, Mat<T>& b) {
    // Por ahora asumimos que P siempre es la matriz identidad, porque en ninguna
    // imagen la fuente viene desde un lugar perpendicular.  Por lo tanto
    // usamos como vector de permutacion el correspondiente al de la identidad.

    // L.U.x = b
    // Sea y = (U.x)
    // L.y = b
    auto y = forward_subtitution(L, b, P); // L es triang inf -> hacemos forward subt

    // obtuve el y tal que Ly=b
    // Ahora busco el x tal que U.x=y
    auto x = backward_subtitution(U, y, P); // U es triang sup -> hacemos backward subt

    return x;
}

template <class T, template<class U> class Mat>
Mat<T> resolverSistemaLU(Mat<T>& L, Mat<T>& U, Mat<T>& b) {
    // TODO Armar P identidad del tamaño correcto!!
    vector<int> P = {0, 1, 2};
    return resolverSistemaLU(L, U, P, b);
}

template <class T>
Matrix<T> cholesky(const Matrix<T>& B) {
    Matrix<T> res = Matrix<T>(B.height(), B.width());

    for (int j = 0; j < B.height(); j++) {
        T sum = T();
        for (int k = 0; k < j; k++) {
            sum += res(j, k) * res(j, k);
        }

        double resta = B(j, j) - sum;
        if (resta <= 0) {
            throw runtime_error("Matriz no SDP: resta dio " + to_string(resta));
        }
        res(j, j) = sqrt(resta);

        for (int i = j + 1; i < B.height(); i++) {
            //cerr << "(" << i << ", " << j << ")" << endl;
            sum = T();
            for (int k = 0; k < j; k++) {
                sum += res(i, k) * res(j, k);
            }
            res(i, j) = (B(i, j) - sum) / res(j, j);
        }
    }

    return res;
}

template <class T>
BandMatrix<T> cholesky(const BandMatrix<T>& B) {
    BandMatrix<T> res = BandMatrix<T>(B.height(), B.width());

    T sum;
    int p = B.upper_bandwidth();
    int q = B.lower_bandwidth();

    for (int j = 0; j < B.height(); j++) {
        sum = T();
        for (int k = max(0, j - q); k < j; k++) {
            sum += res.at(j, k) * res.at(j, k);
        }

        double resta = B(j, j) - sum;
        if (resta <= 0) {
            throw runtime_error("Matriz no SDP: resta dio " + to_string(resta));
        }
        res(j, j) = sqrt(resta);

        for (int i = j + 1; i <= min(B.height() - 1, j + q); i++) {
            sum = T();
            for (int k = max(0, j - q); k < j; k++) {
                sum += res.at(i, k) * res.at(j, k);
            }
            double r = (B(i, j) - sum) / res.at(j, j);
            if (abs(r) >= eps) {
                res(i, j) = (B(i, j) - sum) / res.at(j, j);
            }
        }
    }

    return res;
}


/*Dada la imagen de la mascara devuelve una lista con los pixeles blancos de la 
  misma y un enum diciendo si el pixel esta en el Borde de la imagen

*/
enum Borde { Not, Right, Bottom, Both};
list<tuple<int, int, Borde>> find_valid_pix(const Matrix<uchar> mask_img, const Matrix<normal>& normales){
    
    int w = mask_img.width();
    int h = mask_img.height();

    list<tuple<int, int, Borde>> res = list<tuple<int, int, Borde>>();
    tuple<int,int, Borde> val;

    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
            // cout << (int)mask_img(i,j) << endl; 
            if (mask_img(i,j) > 0 && normales(i, j).z != 0) {
                bool abajo = i+1 < h && mask_img(i+1,j)>0 && normales(i+1, j).z != 0;
                bool derecha = j+1 < w && mask_img(i,j+1)>0 && normales(i, j+1).z != 0;
                if(abajo && derecha){
                    val = make_tuple(i, j, Not);
                } else if(abajo){
                    val = make_tuple(i, j, Right);                    
                } else if(derecha){
                    val = make_tuple(i, j, Bottom);                    
                }else{
                    val = make_tuple(i, j, Both);                    
                }
                res.push_back(val);
                // cout << "i " << i << " j " << j << endl;
            }
        }
    }
    return res;
}

list<tuple<int,int, Borde>>::iterator find_down(list<tuple<int,int, Borde>>::iterator it){
    auto val = *it;
    int original_i = get<0>(val);
    int original_j = get<1>(val);
   // cout << "orig_i " << original_i << endl;
   // cout << "orig_j " << original_j << endl;
    
    auto new_it  = it;
    int new_i, new_j;
    do{
        new_it++;
        auto new_val = *new_it;
        new_i = get<0>(new_val);
        new_j = get<1>(new_val);
       // cout << "new_i " << new_i << " new_j  " << new_j << endl;
    }while(new_i != original_i+1 || new_j != original_j);
    //cout << "Sali" << endl;    
    return new_it;
}

SistemaDeEcuaciones construir_sistema_de_profundidades_Dok(const Matrix<normal>& normales, list<tuple<int, int, Borde>>& valid_pix) {
    // No calcula el z de los Bordes derecho e inferior de la imagen, porque no
    // tenemos información. En esos caso, z va a ser 0.0 
    int size = valid_pix.size();
    auto M = DOKMatrix<double>(2*size, size);
    auto v = Matrix<double>(2*size, 1, 0);

    tuple<int,int, Borde> val;
    int i,j;
    int count = 0;
    Borde b;
    int count_b = 0;
    for(auto it = valid_pix.begin(); it != valid_pix.end(); it++) {
        val = *it;
        i = get<0>(val);
        j = get<1>(val);
        //cout << count << "/" << size << endl;
        b = get<2>(val);
        //if(b == 3){
        //    count_b++;
        //    cout << count_b <<endl;
        //}
        //cout << b << endl;
        
        auto n = normales(i,j);
        // cout << "N.x: " << n.x << " N.y: " << n.y << endl;
        assert (n.z != 0);

        // primera ecuacion
        M(count, count) = -n.z;
        if(b == Not || b == Bottom){
            M(count, count+1) = n.z;
        }
        v(count, 0) = -n.x;

        // segunda ecuacion
        M(count+size, count) = -n.z;
        if (b == Not || b == Right){
            auto bottom_it = find_down(it);
            int bottom_index  = distance(valid_pix.begin(), bottom_it);
            M(count+ size, bottom_index) = n.z;
        }
        v(count+size, 0) = -n.y;

        count ++;
    }
    return {.coeficientes = M, .independiente = v};
}

/* Dada una matriz de normales <nx, ny, nz> construye un sistema de eq que tiene
  a las profundidades de los píxeles como incógnitas.
  enum Borde { Not, Right, Bottom, Both};
  SistemaDeEcuaciones construir_sistema_de_profundidades(int width, int height, const Matrix<normal>& normales, const Matrix<uchar> mask_img) {
      // No calcula el z de los Bordes derecho e inferior de la imagen, porque no
      // tenemos información. En esos caso, z va a ser 0.0 
      list<tuple<int, int, Borde>> valid_pix = find_valid_pix(mask_img);
      int size =valid_pix.size();
      auto M = BandMatrix<double>(2*size, size);
      auto v = Matrix<double>(2*size, 1, 0);
      
      tuple<int,int, Borde> val;
      int i,j;
      int count = 0;
      Borde b;
      int count_b = 0;
      for(auto it = valid_pix.begin(); it != valid_pix.end(); it++) {
          val = *it;
          i = get<0>(val);
          j = get<1>(val);
          //cout << count << "/" << size << endl;
          b = get<2>(val);
          //if(b == 3){
              //    count_b++;
              //    cout << count_b <<endl;
              //}
              //cout << b << endl;
              
              auto n = normales(i,j);
              
              // primera ecuacion
              M(count, count) = -n.z;
              if(b == Not || b == Bottom){
                  M(count, count+1) = n.z;
                }
                v(count, 0) = -n.x;
                
                // segunda ecuacion
                M(count+size, count) = -n.z;
                if (b == Not || b == Right){
                    auto bottom_it = find_down(it);
                    int bottom_index  = distance(valid_pix.begin(), bottom_it);
                    M(count+ size, bottom_index) = n.z;
                }
                v(count+size, 0) = -n.y;
                
                count ++;
            }
            return {.coeficientes = M, .independiente = v};
        }
        */

Matrix<double> extract_z_from(int width, int height, const BandMatrix<double>& Z_dup,const list<tuple<int, int, Borde>>& valid_pix) {
    double min = Z_dup(0,0);
    for(int i = 0; i < Z_dup.rows(); i++){
        double temp = Z_dup(i,0);
        if(temp < min){min = temp;}
    }
    Matrix<double> Z(height,width, min);
    //Matrix<double> Z(height,width, 0);
    int count = 0;
    for(auto it = valid_pix.begin(); it != valid_pix.end(); it++){
        int i = get<0>(*it);
        int j = get<1>(*it);
        Z(j,i) = Z_dup(count, 0);
        count ++; 
    }
    return Z;
}

//Matrix<double> extract_z_from(int width, int height, const BandMatrix<double>& Z_dup) {
//    size_t n = Z_dup.height() / 2;
//    Matrix<double> Z(height, width);
//
//    for (int y = 0; y < height - 1; y++) {
//        for (int x = 0; x < width - 1; x++) {
//            int row = (y * (width - 1) + x) * 2;
//            int i = x == 0 ? row - 1 : row;
//            // cerr << "i = " << i << ", Z_dup.rows() = " << Z_dup.rows() << endl;
//            double z = Z_dup(i, 0);
//            Z(x, y) = z;
//            // cerr << x << " " << y << " " << z << endl;
//        }
//    }
//
//    return Z;
//}

/* Dada una matriz de normales, calcula las profundidades correspondientes a
 * cada pixel.
 */
Matrix<double> obtener_profundidades(const Matrix<normal>& normales, const Matrix<uchar> mask_img) {
    
    int w = normales.width();
    int h = normales.height();
    list<tuple<int, int, Borde>> valid_pix = find_valid_pix(mask_img, normales);
    
    // Construyo el sistema de ecuaciones Mz = v (13.enunciado.pdf)
    //auto sistema = construir_sistema_de_profundidades(w, h, normales, mask_img);
    cerr << "[3dbuild] Construyendo matriz M de ecuaciones normales" << endl;
    auto sistema = construir_sistema_de_profundidades_Dok(normales, valid_pix);

    DOKMatrix<double>& Md = sistema.coeficientes;
    Matrix<double>& v = sistema.independiente;

    if (Md.has_null_rows()) {
        throw std::runtime_error("M tiene filas nulas, no es de rango completo!");
    }

    // (Guarda M y b en archivos para testeo/comparación con Octave/Matlab)
    if (get_env_var("DEBUG") == "1") {
        Md.save("M.dat");
        v.save("v.dat");
        cerr << "[3dbuild] [DEBUG] M.dat y v.dat escritos" << endl;
        cerr << "[3dbuild] [DEBUG] size(M) == (" << Md.height() << ", " << Md.width() << ")" << endl;
        cerr << "[3dbuild] [DEBUG] size(v) == (" << v.height() << ", " << v.width() << ")" << endl;
    }

    // Multiplico por la transpuesta a M y a v (14.enunciado.pdf),
    // obteniendo Az = b (15.enunciado.pdf)
    
    cerr << "[3dbuild] A = M' * M" << endl;
    auto Ad = Md.mult_transposed();
    cerr << "[3dbuild] b = M' * v" << endl;
    const auto Mdt = Md.transposed();
    auto bd = Mdt.mult_vector(v);
    cerr << "[3dbuild] Mdt size: " << Mdt.height() << "x" << Mdt.width() << endl;
    cerr << "[3dbuild] b size: " << bd.height() << "x" << bd.width() << endl;

    if (get_env_var("DEBUG") == "1") {
        // ..Descomentar sólo si se corre con un valor chico de TEST_TRIM=
        Ad.save("A.dat");
        bd.save("b.dat");
        cerr << "[3dbuild] [DEBUG] A.dat y b.dat escritos" << endl;
        cerr << "[3dbuild] [DEBUG] size(A) == (" << Ad.height() << ", " << Ad.width() << ")" << endl;
        cerr << "[3dbuild] [DEBUG] size(b) == (" << bd.height() << ", " << bd.width() << ")" << endl;
    }

    cerr << "[3dbuild] Ya está listo el sistema Az = b para resolver..." << endl;

    BandMatrix<double> z_dup(bd.height(), bd.width());

    // Convierte las matrices DOK a matrices banda
    auto A = Ad.to_band();
    auto b = bd.to_band();

    cerr << "[3dbuild] A es de " << A.rows() << "x" << A.cols()
         << ", banda " << A.lower_bandwidth() << "," << A.upper_bandwidth() << endl;

    if (A.has_null_rows()) {
        throw std::runtime_error("A tiene filas nulas, no es de rango completo!");
    }

    if (get_env_var("Z_SIST") != "EG"){
        // Como (M_trans * M) es simétrica definida positiva, se puede hacer CL
        auto L = cholesky(A);
        //cerr << "Valor del ultimo lugar de L " <<L(L.rows()-1, L.cols()-1) << endl;
        cerr << "[3dbuild] Listo el Cholesky" << endl;
        auto U = L.transposed();
        //cout << "Valor del ultimo lugar de U " <<U(U.rows()-1, U.cols()-1) << endl;
        //cout << "-------------------------" << endl;
        z_dup = resolverSistemaLU(L, U, b);
        //cout << z_dup << endl;
        cerr << "[3dbuild] LU de Cholesky resuelto" << endl;
    } else {
        cerr << "[3dbuild] Haciendo EG" << endl;
        z_dup = eliminacion_gaussiana(A, b);
        cerr << "[3dbuild] EG resuelto" << endl;
    }

    cerr << "[3dbuild] Mapa de profundidad construido!" << endl;

    return extract_z_from(mask_img.height(), mask_img.width(), z_dup, valid_pix);
}

/*---------------------  Obtencion de normales  ----------------------- */

pair<int, int> obtenerPixel(vector<Matrix<uchar> >& imgs){
  // Busco el primer pixel mas brillante.
  return make_pair(0,0);
}

normal calcularNormalSegunM(Matrix<double>& m){
  assert(m.height() == 3 && m.width() == 1);
  double norma2dem = norma_dos(m);

  normal n = {m[0][0]/norma2dem, m[1][0]/norma2dem, m[2][0]/norma2dem};
  return n;
}



/*
  Dadas 3 imágenes, y las direcciones de luces correspondientes a cada una (almacenadas en la matriz S),
  calcula la normal de cada pixel.

  * Si queremos resolver con LU, Factoriza S en LU
  * Si queremos calcular I0*ro usando un solo pixel:
    * Elegimos el pixel usando obtenerPixel
    * Hacemos el calculo usando obtenerIntensidadRo
    * Lo almacenamos en coef
  Para cada pixel (x, y) de una de las imagenes:
      * Armar vector I de los valores en x, y de las 3 imágenes
      Si queremos resolver con EG:
        * Resolver sistema Lx = I con EG
        * Resolver sistema Um = x con EG
        * Entonces tenemos m
      Si no, queremos resolver con LU:
        * Llamamos a resolverSistemaLU, y que se encargue
        * Entonces tenemos m
      Si queremos calcular n con un solo pixel de referencia para el I0*ro:
        * Hacemos n = m / coef
      Si no, queremos calcular n usando cada pixel como referencia para I0*ro
        * Hacemos n = m / norma(m)
*/
Matrix<normal> obtenerPlanoNormal(Matrix<double>& S, vector<Matrix<uchar> >& imgs){
  bool unPixel = 1;//(get_env_var("1PIXEL") == "1"); // si uso solamente un pixel para sacar I0ro, o si uso todos
  bool factLU = (get_env_var("FACTLU") == "1");

  int dimension = S.height(); //Deberia ser siempre 3, pero para no usar un magic number
  Matrix<double> L(dimension, dimension, 0.0);
  Matrix<double> U(dimension, dimension, 0.0);
  Matrix<double> ind(dimension,  1, 0.0);
  Matrix<double> ind_original(dimension,  1, 0.0);
  Matrix<double> matrixPerm(dimension, dimension, 0.0);
  Matrix<double> m(3,1,0);
  Matrix<double> pix_int(3,1,0); // lo voy a ir usando de termino independiente, con I1, I2, I3.
  int alto = imgs[0].height();
  int ancho = imgs[0].width();
  Matrix<normal> plano_normales(alto,ancho,{1,1,1});
  double norma2dem = 0;
  normal n;
  Matrix<double> S_orig = S;
  vector<int> P;

  if (factLU){
    P = factorizacion_lu(S, L, U, ind);
    S = S_orig; // Hacemos esto porque factorizacion lu destroza la S
  }

  // Recorro el rango de la imagen
  for(int i=0; i<imgs[0].height(); i++){
    for(int j=0; j<imgs[0].width(); j++){
      // Armo vector I de los valores en x, y de las 3 imágenes
      // Para cada pixel tengo que resolver el sistema Lx=(I/coeficienteIRo)
      // Me creo el termino independiente para cada pixel
      for(int k=0; k<3; k++){
        ind[k][0] = imgs[k][i][j];
      }
      // Buscamos m para (i,j)
      if (factLU){
          m = resolverSistemaLU(L, U, P, ind);
      } else {
          m = eliminacion_gaussiana(S, ind);
          S = S_orig;
      }
      // Tenemos m!
      norma2dem = norma_dos(m);
      if (norma2dem == 0) {
        n = {0,0,0};
      }else{
        n = {m[0][0]/norma2dem, m[1][0]/norma2dem, m[2][0]/norma2dem};

      }
      plano_normales[i][j] = n;

    }
  }

  return plano_normales;
}
