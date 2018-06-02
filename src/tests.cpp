#include "assert.h"
#include <iostream>
#include "matrix.h"
#include "band_matrix.h"
#include "dok_matrix.h"
#include "algoritmosmatriciales.h"

// Helpers

/*
void gen_rand_normal() {
    double x, y, z, norm2;
    x = rand(); y = rand(); z = rand();
    norm2 = sqrt(x*x + y*y + z*z);
    return {x: x/norm2, y: y/norm2, z: z/norm2};
}
*/

/// Tests

void test_mult(){
  Matrix<double> id (3,3,0);
  id[0] = {1,0,0};
  id[1] = {0,1,0};
  id[2] = {0,0,1};

  Matrix<double> unos(3,1,0);
  unos[0][0] = 1;
  unos[1][0] = 1;
  unos[2][0] = 1;

  assert(unos.eq_aprox(id.mult(unos)));

  Matrix<double> id_4_dimensiones (4,4,0);
  id_4_dimensiones[0] = {1,0,0,0};
  id_4_dimensiones[1] = {0,1,0,0};
  id_4_dimensiones[2] = {0,0,1,0};
  id_4_dimensiones[3] = {0,0,0,1};

  Matrix<double> unos_4_dimensiones(4,1,0);
  unos_4_dimensiones[0][0] = 1;
  unos_4_dimensiones[1][0] = 1;
  unos_4_dimensiones[2][0] = 1;
  unos_4_dimensiones[3][0] = 1;

  assert(unos_4_dimensiones.eq_aprox(id_4_dimensiones.mult(unos_4_dimensiones)));
}

void test_mult_transposed() {
    Matrix<int> M({ {8,1,6}, {3,5,7}, {4,9,2}, {-1,5,14}, {1,-5,0} });
    auto M_mt = M.mult_transposed();

    Matrix<int> expected_M_mt({ {91,49,63}, {49,157,129}, {63,129,285} });
    assert(M_mt == expected_M_mt);
}


void test_forward_substitution(){
  // -------------- Caso trivial
  Matrix<double> A_1 (3,3,0);
  A_1[0] = {1,0,0};
  A_1[1] = {0,1,0};
  A_1[2] = {0,0,1};

  Matrix<double> ind_1(3,1,0);
  ind_1[0][0] = 1;
  ind_1[1][0] = 2;
  ind_1[2][0] = 3;

  Matrix<double> x_1 (3,1,0);
  x_1[0][0] = 0;
  x_1[1][0] = 0;
  x_1[2][0] = 0;
  vector<int> perm_default_3_dimensiones = {0,1,2};

  x_1 = forward_subtitution(A_1, ind_1, perm_default_3_dimensiones);

  assert(ind_1.eq_aprox(x_1));

  // -------------- Caso trivial, dimension 3
  Matrix<double> A_2 (4,4,0);
  A_2[0] = {1,0,0,0};
  A_2[1] = {0,1,0,0};
  A_2[2] = {0,0,1,0};
  A_2[3] = {0,0,0,1};

  vector<int> perm_default_4_dimensiones = {0,1,2,3};

  Matrix<double> ind_2(4,1,0);
  ind_2[0][0] = 1;
  ind_2[1][0] = 2;
  ind_2[2][0] = 3;
  ind_2[3][0] = 4;

  Matrix<double> x_2_correcto(4,1,0); // lo que deberia dar
  x_2_correcto[0][0] = 1;
  x_2_correcto[1][0] = 2;
  x_2_correcto[2][0] = 3;
  x_2_correcto[3][0] = 4;

  Matrix<double> x_2 (4,1,0);
  x_2[0][0] = 0;
  x_2[1][0] = 0;
  x_2[2][0] = 0;
  x_2[3][0] = 0;

  x_2 = forward_subtitution(A_2, ind_2, perm_default_4_dimensiones);

  assert(x_2.eq_aprox(x_2_correcto));



  // -------------- Caso mas normal
  Matrix<double> A_3 (4,4,0);
  A_3[0] = {1,0,0,0};
  A_3[1] = {2,1,0,0};
  A_3[2] = {0.5,3,1,0};
  A_3[3] = {-1,-0.5,2,1};

  vector<int> perm_default_4_dimensiones_1 = {0,1,2,3};


  Matrix<double> ind_3(4,1,0);
  ind_3[0][0] = 1;
  ind_3[1][0] = 2;
  ind_3[2][0] = 3;
  ind_3[3][0] = 4;

  Matrix<double> x_3_correcto(4,1,0); // lo que deberia dar
  x_3_correcto[0][0] = 1;
  x_3_correcto[1][0] = 0;
  x_3_correcto[2][0] = 2.5;
  x_3_correcto[3][0] = 0;

  Matrix<double> x_3 (4,1,0.0);
  x_3[0][0] = 0.0;
  x_3[1][0] = 0.0;
  x_3[2][0] = 0.0;
  x_3[3][0] = 0.0;


  x_3 = forward_subtitution(A_3, ind_3, perm_default_4_dimensiones_1);

  // verificamos que el test este bien hecho
  assert(ind_3.eq_aprox(A_3.mult(x_3_correcto)));
  assert(x_3.eq_aprox(x_3_correcto));
}

void test_backward_substitution(){
  // -------------- Caso trivial
  Matrix<double> A_1 (3,3,0);
  A_1[0] = {1,0,0};
  A_1[1] = {0,1,0};
  A_1[2] = {0,0,1};

  Matrix<double> ind_1(3,1,0);
  ind_1[0][0] = 1;
  ind_1[1][0] = 2;
  ind_1[2][0] = 3;

  Matrix<double> x_1 (3,1,0);
  x_1[0][0] = 0;
  x_1[1][0] = 0;
  x_1[2][0] = 0;
  vector<int> perm_default_3_dimensiones = {0,1,2};

  x_1 = backward_subtitution(A_1, ind_1, perm_default_3_dimensiones);

  assert(ind_1.eq_aprox(x_1));

  // -------------- Caso trivial, dimension 3
  Matrix<double> A_2 (4,4,0);
  A_2[0] = {1,0,0,0};
  A_2[1] = {0,1,0,0};
  A_2[2] = {0,0,1,0};
  A_2[3] = {0,0,0,1};

  vector<int> perm_default_4_dimensiones = {0,1,2,3};

  Matrix<double> ind_2(4,1,0);
  ind_2[0][0] = 1;
  ind_2[1][0] = 2;
  ind_2[2][0] = 3;
  ind_2[3][0] = 4;

  Matrix<double> x_2_correcto(4,1,0); // lo que deberia dar
  x_2_correcto[0][0] = 1;
  x_2_correcto[1][0] = 2;
  x_2_correcto[2][0] = 3;
  x_2_correcto[3][0] = 4;

  Matrix<double> x_2 (4,1,0);
  x_2[0][0] = 0;
  x_2[1][0] = 0;
  x_2[2][0] = 0;
  x_2[3][0] = 0;

  x_2 = backward_subtitution(A_2, ind_2, perm_default_4_dimensiones);

  assert(x_2.eq_aprox(x_2_correcto));



  // -------------- Caso mas normal
  Matrix<double> A_3 (4,4,0);
  A_3[0] = {6,-2,2,4};
  A_3[1] = {0,-4,2,2};
  A_3[2] = {0,0,2,-5};
  A_3[3] = {0,0,0,-3};

  vector<int> perm_default_4_dimensiones_1 = {0,1,2,3};

  Matrix<double> ind_3(4,1,0);
  ind_3[0][0] = 1;
  ind_3[1][0] = 2;
  ind_3[2][0] = 3;
  ind_3[3][0] = 4;

  Matrix<double> x_3_correcto(4,1,0); // lo que deberia dar
  x_3_correcto[0][0] = 35.0/36.0;
  x_3_correcto[1][0] = -25.0/12.0;
  x_3_correcto[2][0] = -11.0/6.0;
  x_3_correcto[3][0] = -4.0/3.0;

  Matrix<double> x_3 (4,1,0.0);
  x_3[0][0] = 0.0;
  x_3[1][0] = 0.0;
  x_3[2][0] = 0.0;
  x_3[3][0] = 0.0;


  x_3 = backward_subtitution(A_3, ind_3, perm_default_4_dimensiones_1);

  // verificamos que el test este bien hecho
  assert(ind_3.eq_aprox(A_3.mult(x_3_correcto)));
  assert(x_3.eq_aprox(x_3_correcto));
}



// Dado Ax=b, chequea factorizacion PLU y resolucion de Ax=b usando su fact PLU
void test_factorizar_y_resolver_PLU(){
  // ejercicio visto en la teorica
  Matrix<double> A (4,4,0); // matriz
  A[0] = {6,-2,2,4};
  A[1] = {12,-8,6,10};
  A[2] = {3,-13,9,3};
  A[3] = {-6,4,1,-18};

  Matrix<double> ind(4,1,0); // termino independiente
  ind[0][0] = 1;
  ind[1][0] = 2;
  ind[2][0] = 3;
  ind[3][0] = 4;

  Matrix<double> id (4,4,0); // matriz identidad
  id[0] = {1,0,0,0};
  id[1] = {0,1,0,0};
  id[2] = {0,0,1,0};
  id[3] = {0,0,0,1};

  vector<int> P_test = {0, 1, 2, 3};

  Matrix<double> L_test (4,4,0);
  L_test[0] = {1,0,0,0};
  L_test[1] = {2,1,0,0};
  L_test[2] = {0.5,3,1,0};
  L_test[3] = {-1,-0.5,2,1};

  Matrix<double> U_test (4,4,0);
  U_test[0] = {6,-2,2,4};
  U_test[1] = {0,-4,2,2};
  U_test[2] = {0,0,2,-5};
  U_test[3] = {0,0,0,-3};

  vector<double> res_test = {-1/24,5/8,5/4,0.0};

  //chequeamos que el test este bien armado
  assert(A.eq_aprox(L_test.mult(U_test)));

  //_____ Factorizamos A en PLU
  Matrix<double> L (4,4,0);
  Matrix<double> U (4,4,0);


  vector<int> P = factorizacion_lu(A, L, U, ind);

  // chequeamos que P, L y U den aproximadamente lo mismo
  assert(P == P_test);
  assert(L.eq_aprox(L_test));
  assert(U.eq_aprox(U_test));


  // Como P es identidad en este caso, tenemos LUx = b.  Llamamos Ux = y, y resolvemos Ly = b
  //vector<int> perm_default = {1,2,3,4};
  //Matrix<double> y (L.height(), 1, 0);
  //Matrix<double> z = forward_subtitution(L, y, perm_default); // TODO - arreglar


}


// Chequea las funciones de cholesky, transpose y multiplicacionTranspuesta
void test_cholesky_transpose() {
  Matrix<double> A1_test = Matrix<double>(5,5,0);
  A1_test[0] = {5.99246,   0.45174,   0.50921,   0.45632,   0.41578};
  A1_test[1] = {0.45174,   5.02065,   0.40611,   0.52172,   0.11483};
  A1_test[2] = {0.50921,   0.40611,   5.25767,   0.78430,   0.24938};
  A1_test[3] = {0.45632,   0.52172,   0.78430,   5.98219,   0.45188};
  A1_test[4] = {0.41578,   0.11483,   0.24938,   0.45188,   5.2870};


  Matrix<double> L1_test(5,5,0);
  L1_test[0] = {2.44795,   0.00000,   0.00000,   0.00000,   0.00000};
  L1_test[1] = {0.18454,   2.23307,   0.00000,   0.00000,   0.00000};
  L1_test[2] = {0.20802,   0.16467,   2.27756,   0.00000,   0.00000};
  L1_test[3] = {0.18641,   0.21823,   0.31156,   2.40889,   0.00000};
  L1_test[4] = {0.16985,   0.03739,   0.09128,   0.15925,   2.28540};

  Matrix<double> A(3,3,0);
  A[0] = {1.0,   0.00000,   0.00000};
  A[1] = {4.0,   2.0,   0.00000};
  A[2] = {6.0,   2.0,   4.0};

  Matrix<double> B(3,3,0);
  B = A.transpose();

  Matrix<double> C(3,3,0);
  C = A.mult(B);


  Matrix<double> L1 = cholesky(A1_test);
  Matrix<double> A1 = L1.multiplicarPorTraspuesta();
  Matrix<double> A2 = L1.mult(L1.transpose());

  assert(L1.eq_aprox(L1_test));
  assert(A1.eq_aprox(A1_test));
  assert(A1.eq_aprox(A2));
}

void test_obtener_intensidad_fuente_y_ro(){
  double I0ro_correcto = 11.1064;
  Matrix<double> S(3,3,0);
  S[0] = {2,3,4};
  S[1] = {6,14,18};
  S[2] = {4,16,29};

  Matrix<double> S_copia(3,3,0);
  S_copia = S;

  vector<int> P_correcto = {0,1,2};

  Matrix<double> L_correcto(3,3,0);
  L_correcto[0] = {1,0,0};
  L_correcto[1] = {3,1,0};
  L_correcto[2] = {2,2,1};

  Matrix<double> U_correcto(3,3,0);
  U_correcto[0] = {2,3,4};
  U_correcto[1] = {0,5,6};
  U_correcto[2] = {0,0,9};

  Matrix<double> I(3,1,0);
  I[0] = {7};
  I[1] = {3};
  I[2] = {2};

  Matrix<double> I_copia(3,1,0);
  I_copia[0] = {7};
  I_copia[1] = {3};
  I_copia[2] = {2};

  Matrix<double> L(3,3,0);
  Matrix<double> U(3,3,0);

  vector<int> P = factorizacion_lu(S, L, U, I);

  // chequeamos que factorizacion_LU da lo esperado
  assert(S_copia.eq_aprox(L.mult(U)));
  assert(P == P_correcto);
  assert(L.eq_aprox(L_correcto));
  assert(U.eq_aprox(U_correcto));

  //assert(eq_escalar(obtenerIntensidadFuenteRo(L, U, P, I_copia), 11.1064));


}

void test_construir_sistema_de_profundidades_01(){
  auto M = Matrix<normal>(3, 3, {0, 0, 0});
  M[0] = {{ 0,  1,  2}, { 3,  4,  5}, {99, 99, 99}};
  M[1] = {{ 9, 10, 11}, {12, 13, 14}, {99, 99, 99}};
  M[2] = {{99, 99, 99}, {99, 99, 99}, {99, 99, 99}};

  // imagen de 3x3 pero sin bordes -> 2x2 -> N=4, 2*N=8 -> coef es 8x8
  auto coef_expected = BandMatrix<double>(8, 8);
  coef_expected(0, 0) = -2.0;
  coef_expected(0, 1) = 2.0;
  coef_expected(1, 0) = -2.0;
  coef_expected(1, 2) = 2.0;

  coef_expected(2, 1) = -5.0;
  coef_expected(2, 3) = 5.0;
  coef_expected(3, 1) = -5.0;

  coef_expected(4, 4) = -11.0;
  coef_expected(4, 5) = 11.0;
  coef_expected(5, 4) = -11.0;
  coef_expected(5, 6) = 11.0;

  coef_expected(6, 5) = -14.0;
  coef_expected(6, 7) = 14.0;
  coef_expected(7, 5) = -14.0;

  auto indep_expected = Matrix<double>(8, 1, 0);
  indep_expected(0, 0) =   0;
  indep_expected(1, 0) =  -1;
  indep_expected(2, 0) =  -3;
  indep_expected(3, 0) =  -4;
  indep_expected(4, 0) =  -9;
  indep_expected(5, 0) = -10;
  indep_expected(6, 0) = -12;
  indep_expected(7, 0) = -13;

//  auto res = construir_sistema_de_profundidades(3, 3, M);

  /*
  cout << res.coeficientes;
  cout << "\n";
  cout << coef_expected;

  cout << res.independiente;
  cout << "\n";
  cout << indep_expected;
  */

//  assert(res.coeficientes.eq_aprox(coef_expected));
//  assert(res.independiente.eq_aprox(indep_expected));
}


void test_obtener_normales(){
  // va a ser parecido al que calcula intensidad de fuente y ro
  Matrix<double> S(3,3,0);
  S[0] = {2,3,4};
  S[1] = {6,14,18};
  S[2] = {4,16,29};

  double I0ro = 11.1064;

  Matrix<double> S_copia(3,3,0);
  S_copia = S;

  vector<int> P_correcto = {1,2,3};

  Matrix<double> L_correcto(3,3,0);
  L_correcto[0] = {1,0,0};
  L_correcto[1] = {3,1,0};
  L_correcto[2] = {2,2,1};

  Matrix<double> U_correcto(3,3,0);
  U_correcto[0] = {2,3,4};
  U_correcto[1] = {0,5,6};
  U_correcto[2] = {0,0,9};

  Matrix<double> I(3,1,0);
  I[0] = {7};
  I[1] = {3};
  I[2] = {2};

  Matrix<double> I_copia(3,1,0);
  I_copia[0] = {7};
  I_copia[1] = {3};
  I_copia[2] = {2};

  Matrix<double> L(3,3,0);
  Matrix<double> U(3,3,0);

  vector<int> P = factorizacion_lu(S, L, U, I);

  // chequeamos que factorizacion_LU da lo esperado
  assert(S_copia.eq_aprox(L.mult(U)));
  assert(P == P_correcto);
  assert(L.eq_aprox(L_correcto));
  assert(U.eq_aprox(U_correcto));

  vector<Matrix<double>> imgs(3, Matrix<double>(1,1,0)); // las 3 imagenes son de 1 pixel cada una
  imgs[0][0][0] = 7; // imagen 0, en el punto (0,0)
  imgs[1][0][0] = 3; // imagen 1, en el punto (0,0)
  imgs[2][0][0] = 2; // imagen 2, en el punto (0,0)

  //normal n = obtenerNormales(L, U, I0ro, imgs);

  //assert(eq_escalar(n.x, (23/9)));
  //assert(eq_escalar(n.y, (1/3)));
  //assert(eq_escalar(n.z, (2/9)));
}

void test_eliminacion_gaussiana(){
  Matrix<double> x (3,1,0);
  Matrix<double> A (3,1,0);
  A[0] = {2,3,4};
  A[1] = {6,14,18};
  A[2] = {4,16,29};

  Matrix<double> b (3,1,0);
  b[0][0] = {7};
  b[1][0] = {3};
  b[2][0] = {2};

  Matrix<double> x_correcto (3,1,0);
  x_correcto[0][0] = 251.0/30.0;
  x_correcto[1][0] = -34.0/5.0;
  x_correcto[2][0] = 8.0/3.0;

  x = eliminacion_gaussiana(A, b);

  assert(x.eq_aprox(x_correcto));
}

void test_slice() {
    Matrix<int> A(4, 4, 0);

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            A[i][j] = i * 4 + j;
        }
    }
    assert(A.width() == 4);
    assert(A.height() == 4);

    auto B = A.slice(1, 4, 2, 4);
    assert(B.width() == 2);
    assert(B.height() == 3);

    for (int i = 0; i < B.height(); i++) {
        for (int j = 0; j < B.width(); j++) {
            assert(B[i][j] == A[i+1][j+2]);
        }
    }
}


void test_plano_normal_varios_pixel(){
  /* Idea del test:
  - Tengo 3 imagenes de 4 pixeles cada una.
  Img_a = /7  20\
          \4   8/

  Img_b = /3  88\
          \8   32/

  Img_c = /2  123\
          \9   59/

  Y una matriz S, que representa las direcciones de luces (sacada de las primeras 3 direcciones de luces.txt):
  / 0.403259,0.480808,0.778592\
  | 0.0982272,0.163712,0.981606|
  \-0.0654826,0.180077,0.98147/

  El algoritmo debe:

  - Encontrar I0*ro


  - Encontrar m (dimension 3x1) tal que:
                              /img_a[0][0]\                                 /img_a[0][1]\
  S *  m = Intensidad_(0,0) = |img_b[0][0]|     S *  m = Intensidad_(0,0) = |img_b[0][1]|
                              \img_c[0][0]/                                 \img_c[0][1]/

                              /img_a[1][0]\                                 /img_a[1][1]\
  S *  m = Intensidad_(1,0) = |img_b[1][0]|     S *  m = Intensidad_(0,0) = |img_b[1][1]|
                              \img_c[1][0]/                                 \img_c[1][1]/

  - Luego se despeja n = m/(I0*ro)

*/
  // S van a ser mis direcciones de luces, calculadas en la calibracion
  /*Para otro test, puede ser:
  Matrix<double> S(3,3,0);
  S[0] = {0.403259,0.480808,0.778592};
  S[1] = {0.0982272,0.163712,0.981606};
  S[2] = {-0.0654826,0.180077,0.98147};*/

  Matrix<double> S(3,3,0);
  S[0] = {2,3,4};
  S[1] = {6,14,18};
  S[2] = {4,16,29};

  vector<Matrix<uchar>> imgs(3, Matrix<uchar>(2,2,0)); // las 3 imagenes son de 2 pixeles cada una

  Matrix<uchar> img_a = Matrix<uchar>(2,2,0);
  img_a[0] = {7,20};
  img_a[1] = {4,8};

  Matrix<uchar> img_b = Matrix<uchar>(2,2,0);
  img_b[0] = {3,88};
  img_b[1] = {8,32};

  Matrix<uchar> img_c = Matrix<uchar>(2,2,0);
  img_c[0] = {2,123};
  img_c[1] = {9,59};

  imgs = {img_a, img_b, img_c};

  // Las m que tienen que dar
  Matrix<Matrix<double> > m_correctas(2,2,Matrix<double>(3,1,0)); // Matrix de vectores columna

  m_correctas[0][0][0][0] = {251.0/30};
  m_correctas[0][0][1][0] = {-34.0/5};
  m_correctas[0][0][2][0] = {8.0/3};

  m_correctas[0][1][0][0] = {1};
  m_correctas[0][1][1][0] = {2};
  m_correctas[0][1][2][0] = {3};

  m_correctas[1][0][0][0] = {3};
  m_correctas[1][0][1][0] = {-2};
  m_correctas[1][0][2][0] = {1};

  m_correctas[1][1][0][0] = {1};
  m_correctas[1][1][1][0] = {-2};
  m_correctas[1][1][2][0] = {3};

  // Me armo la matriz de normales n_correctas que me debería devolver. Son una por pixel
  Matrix<normal> n_correctas = Matrix<normal>(2,2,{0,0,0});
  n_correctas[0][0] = calcularNormalSegunM(m_correctas[0][0]);
  n_correctas[0][1] = calcularNormalSegunM(m_correctas[0][1]);
  n_correctas[1][0] = calcularNormalSegunM(m_correctas[1][0]);
  n_correctas[1][1] = calcularNormalSegunM(m_correctas[1][1]);

  Matrix<normal> n = obtenerPlanoNormal(S, imgs);

  for (int i = 0; i < n.height(); ++i){
    for (int j = 0; j < n.width(); ++j){
      assert(eq_escalar(n[i][j].x, n_correctas[i][j].x));
      assert(eq_escalar(n[i][j].y, n_correctas[i][j].y));
      assert(eq_escalar(n[i][j].z, n_correctas[i][j].z));
    }
  }
}

void test_matrix() {
    Matrix<double> A(3,3,0);
    A[0] = {0.403259,0.480808,0.778592};
    A[1] = {0.0982272,0.163712,0.981606};
    A[2] = {-0.0654826,0.180077,0.98147};

    Matrix<double> B = A;

    assert(A == B);
}

void test_band_matrix() {
    BandMatrix<int> A(5, 5);

    A(0, 0) = 1;
    A(1, 1) = 1;
    A(2, 2) = 1;

    A(1, 0) = 2;
    A(2, 1) = 4;
    A(3, 2) = 5;
    A(2, 4) = 9;

    assert(A.height() == 5);
    assert(A.width() == 5);

    assert(A(4, 2) != A(2, 4));
    assert(A(2, 4) == 9);
}

void test_band_matrix_symmetric() {
    BandMatrix<int> A(5, 5, true);

    A(0, 0) = 1;
    A(1, 1) = 1;
    A(2, 2) = 1;

    A(1, 0) = 2;
    A(2, 1) = 4;
    A(3, 2) = 5;
    A(2, 4) = 9;

    assert(A.height() == 5);
    assert(A.width() == 5);

    assert(A(4, 2) == A(2, 4));
    assert(A(0, 1) == A(1, 0));
    assert(A(4, 4) == 0);
}

void test_band_matrix_transpose() {
    BandMatrix<int> A(3, 2);
    A(0, 0) = 1;
    A(1, 1) = 2;
    A(0, 1) = 4;

    assert(A.rows() == 3);
    assert(A.cols() == 2);

    auto At = A.transposed();

    assert(At.rows() == A.cols());
    assert(At.cols() == A.rows());

    auto B = At.transposed();

    assert(B.rows() == A.rows());
    assert(B.cols() == A.cols());

    assert(B == A);
    assert(A.eq_aprox(B));
}

void test_band_matrix_mult() {
    BandMatrix<int> A(5, 5);
    for (int i = 0; i < 5; i++) {
        A(i, i) = 1;
    }
    A(1, 0) = 2;
    A(2, 1) = 2;
    A(4, 3) = 2;
    A(1, 2) = 3;
    A(2, 3) = 3;
    A(2, 4) = 4;

    assert(A.lower_bandwidth() == 1);
    assert(A.upper_bandwidth() == 2);

    BandMatrix<int> A_exp(5, 5, true);
    A_exp(0, 0) = 5; A_exp(0, 1) = 2; A_exp(0, 2) = 6;
    A_exp(1, 0) = 2; A_exp(1, 1) = 5; A_exp(1, 2) = 5; A_exp(1, 3) = 6; A_exp(1, 4) = 8;
    A_exp(2, 0) = 6; A_exp(2, 1) = 5; A_exp(2, 2) = 10; A_exp(2, 3) = 3; A_exp(2, 4) = 4;
    A_exp(3, 1) = 6; A_exp(3, 2) = 3; A_exp(3, 3) = 14; A_exp(3, 4) = 14;
    A_exp(4, 1) = 8; A_exp(4, 2) = 4; A_exp(4, 3) = 14; A_exp(4, 4) = 17;

    assert(A.mult_transposed() == A_exp);
}

void test_dok_matrix() {
    DOKMatrix<int> A(5, 5);

    A(0, 0) = 1;
    A(1, 1) = 1;
    A(2, 2) = 1;

    A(1, 0) = 2;
    A(2, 1) = 4;
    A(3, 2) = 5;
    A(2, 4) = 9;

    assert(A.height() == 5);
    assert(A.width() == 5);

    assert(A(4, 2) != A(2, 4));
    assert(A(2, 4) == 9);
}

void test_dok_matrix_transpose() {
    DOKMatrix<int> A(3, 2);
    A(0, 0) = 1;
    A(1, 1) = 2;
    A(0, 1) = 4;

    assert(A.rows() == 3);
    assert(A.cols() == 2);

    auto At = A.transposed();

    assert(At.rows() == A.cols());
    assert(At.cols() == A.rows());

    auto B = At.transposed();

    assert(B.rows() == A.rows());
    assert(B.cols() == A.cols());

    assert(B == A);
    assert(A.eq_aprox(B));
}

void test_dok_matrix_mult() {
    DOKMatrix<int> A(5, 5);
    for (int i = 0; i < 5; i++) {
        A(i, i) = 1;
    }
    A(1, 0) = 2;
    A(2, 1) = 2;
    A(4, 3) = 2;
    A(1, 2) = 3;
    A(2, 3) = 3;
    A(2, 4) = 4;

    DOKMatrix<int> A_exp(5, 5);
    A_exp(0, 0) = 5; A_exp(0, 1) = 2; A_exp(0, 2) = 6;
    A_exp(1, 0) = 2; A_exp(1, 1) = 5; A_exp(1, 2) = 5; A_exp(1, 3) = 6; A_exp(1, 4) = 8;
    A_exp(2, 0) = 6; A_exp(2, 1) = 5; A_exp(2, 2) = 10; A_exp(2, 3) = 3; A_exp(2, 4) = 4;
    A_exp(3, 1) = 6; A_exp(3, 2) = 3; A_exp(3, 3) = 14; A_exp(3, 4) = 14;
    A_exp(4, 1) = 8; A_exp(4, 2) = 4; A_exp(4, 3) = 14; A_exp(4, 4) = 17;



    assert(A.mult_transposed() == A_exp);
}

void test_plano_normal(){
  setenv("FACTLU", "0", true); // setea FACTLU en 0, para probar con EG
  test_plano_normal_varios_pixel();

  setenv("FACTLU", "1", true); // setea FACTLU en 1, para probar con fact LU
  test_plano_normal_varios_pixel();
}

void test_performance_obtenerPlanoNormal(int n){
  Matrix<double> S(3,3,0);
  S[0] = {0.403259,0.480808,0.778592};
  S[1] = {0.0982272,0.163712,0.981606};
  S[2] = {-0.0654826,0.180077,0.98147};
  Matrix<uchar> Unos(n,n,0);
  for(int i=0; i<n;i++){
    for(int j=0; j<n;j++){
      Unos[i][j] = 1;
    }
  }
  Matrix<uchar> Uno2 = Unos;
  Matrix<uchar> Uno3 = Unos;
  vector<Matrix<uchar> > imgs;
  imgs.push_back(Unos);
  imgs.push_back(Uno2);
  imgs.push_back(Uno3);
  clock_t start = clock();
  auto normales = obtenerPlanoNormal(S, imgs);
  clock_t end = clock();
  double segs = (double)(end-start) / CLOCKS_PER_SEC;
  cout << segs << endl;
}

int main() {
  test_cholesky_transpose();
  test_factorizar_y_resolver_PLU();
  test_backward_substitution();
  test_forward_substitution();
  test_mult();
  test_mult_transposed();
  //test_construir_sistema_de_profundidades_01();
  test_obtener_intensidad_fuente_y_ro(); // todo: deprecado, hay que actualizar el test
  test_eliminacion_gaussiana();
  test_slice();
  test_matrix();
  test_band_matrix();
  test_band_matrix_symmetric();
  test_band_matrix_transpose();
  test_band_matrix_mult();
  test_dok_matrix();
  test_dok_matrix_transpose();
  test_dok_matrix_mult();
  test_plano_normal();

  cerr << "Pasaron todos los test!" << endl;

  int n = 100;
  if (!get_env_var("MATRIX_SIZE").empty()) {
    n = stoi(get_env_var("MATRIX_SIZE"));
  }
  test_performance_obtenerPlanoNormal(n);
}
