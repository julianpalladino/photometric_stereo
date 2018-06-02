#pragma once
//#include <cassert>
#include <cmath>
#include <iostream>
#include <cassert>
#include <fstream>
#include <tuple>
#include <vector>
#include "common.h"

using namespace std;

template <class T>
class Matrix {
    vector<vector<T>> M;
    int rows;
    int cols;

  public:
    Matrix(int rows, int cols, T def = T()) :
        rows(rows),
        cols(cols),
        M(rows, vector<T>(cols, def)) { };

    Matrix(const Matrix<T>& rhs) :
        rows(rhs.rows),
        cols(rhs.cols),
        M(rhs.M) { };

    Matrix(const vector<vector<T>>& rhs) :
        rows(rhs.size()),
        cols(rhs[0].size()) {
            M.resize(rows);
            for (int i = 0; i < rows; i++) {
                M[i].resize(cols);
                for (int j = 0; j < cols; j++) {
                    assert(j < (int) rhs[i].size());
                    M[i][j] = rhs[i][j];
                }
            }
        };

    bool operator==(const Matrix<T>& rhs) const;

    Matrix<T> multiplicarPorTraspuesta() const;

    int width() const  { return cols; };
    int height() const { return rows; };

    vector<T>& operator[](int index);
    const vector<T>& operator[](int index) const;

    T& operator()(int i, int j);
    const T operator()(int i, int j) const;
    const T at(int i, int j) const;

    Matrix<T> mult(const Matrix<T>& B) const;
    Matrix<T> transpose() const;
    Matrix<T> mult_transposed() const;
    const tuple<int, int, T> max() const;

    void save(const string& filename) const;

    template <typename U>
    friend ostream& operator<<(ostream& stream, const Matrix<U>& M);

    bool esCuadrada();
    bool esSimetricaYCuadrada();
    bool eq_aprox (const Matrix<T>& B) const;

    Matrix<T> slice(int rows_from, int rows_to, int cols_from, int cols_to) const;
};

template <class T>
Matrix<T> Matrix<T>::transpose() const{
  Matrix<T> res = Matrix<T>(cols, rows, 0);
  for (int i = 0; i < rows; ++i)
  {
    for (int j = 0; j < cols; ++j)
    {
      res[j][i] = M[i][j];
    }
  }
  return res;
}

template <class T>
Matrix<T> Matrix<T>::mult_transposed() const {
    Matrix<T> res(cols, cols, 0);

    for (int k = 0; k < rows; k++) {
        for (int i = 0; i < cols; i++) {
            for (int j = i; j < cols; j++) {
                res(i, j) += this->at(k, i) * this->at(k, j);
            }
        }
    }

    // Por simetría, copio los valores
    for (int i = 0; i < cols; i++) {
        for (int j = 0; j < cols; j++) {
            res(j, i) = res(i, j);
        }
    }

    return res;
}

template<class T>
bool Matrix<T>::operator==(const Matrix<T>& rhs) const {
    return rows == rhs.rows && \
        cols == rhs.cols && \
        M == rhs.M;
}

template <class T>
bool Matrix<T>::eq_aprox(const Matrix<T>& B) const {
  if (height() != B.height() || width() != B.width()) return false;
  for(int i=0; i<rows; i++){
    for(int j=0; j<cols; j++){
      if(abs(M[i][j] - B[i][j]) >= eps){ // ¬(|a-b|<epsilon)
        return false;
      }
    }
  }
  return true;
}

// Busca el máximo valor en la matriz, y devuelve una tupla (fila, col, valor)
// con el valor y su posición.
//
// Complejidad peor caso: O(n^2)
//
template <class T>
const tuple<int, int, T> Matrix<T>::max() const {
    int x, y;
    T v = 0;

    for (int h = 0; h < rows; h++) {
        for (int w = 0; w < cols; w++) {
            if (M[h][w] > v) {
                x = w;
                y = h;
                v = M[h][w];
            }
        }
    }

    return make_tuple(x, y, v);
}

template <class T>
vector<T>& Matrix<T>::operator[](int i) {
    return M[i];
}

template <class T>
const vector<T>& Matrix<T>::operator[](int i) const {
    return M[i];
}

template <class T>
T& Matrix<T>::operator()(int i, int j) {
    //assert(i >= 0 && i < rows);
    //assert(j >= 0 && j < cols);
    return M[i][j];
}

template <class T>
const T Matrix<T>::operator()(int i, int j) const {
    //assert(i >= 0 && i < rows);
    //assert(j >= 0 && j < cols);
    return M[i][j];
}

template <class T>
const T Matrix<T>::at(int i, int j) const {
    //assert(i >= 0 && i < rows);
    //assert(j >= 0 && j < cols);
    return M[i][j];
}

template <class T>
Matrix<T> Matrix<T>::mult(const Matrix<T>& B) const {
    if (cols != B.height()) {
        throw runtime_error("No se puede multiplicar matriz de " \
            "(" + to_string(rows) + ", " + to_string(cols) + ") con " \
            "(" + to_string(B.height()) + ", " + to_string(B.width()));
    }

    Matrix<T> res(rows, B.width(), 0);

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < B.width(); j++) {
          res[i][j] = 0;
            for (int k = 0; k < B.height(); k++) {
              res[i][j] += M[i][k] * B[k][j];
            }
        }
    }

    return res;
}

// Multiplica la matriz consigo misma.
// Considerando que (A * A')' = A'' * A' = A * A'
// Tambien sabemos que la matriz resultante es simetrica
template <class T>
Matrix<T> Matrix<T>::multiplicarPorTraspuesta() const {
  Matrix<T> res(rows, cols, 0);

  for (size_t i = 0; i < rows; i++) {
    for (size_t j = i; j < cols; j++) {
      res[i][j] = 0;
      for (size_t k = 0; k < cols; k++) {
        res[i][j] += M[i][k] * M[j][k];
      }
      res[j][i] = res[i][j];
    }
  }
  return res;
}


template <class T>
ostream& operator<<(ostream& out, const Matrix<T>& M) {
    for (int i = 0; i < M.height(); i++) {
        for (int j = 0; j < M.width(); j++) {
            out << M[i][j] << "\t";
        }
        out << "\n";
    }
    return out;
}

template <class T>
Matrix<T> Matrix<T>::slice(int rows_from, int rows_to, int cols_from, int cols_to) const {
    if (rows_from < 0 || rows_from >= rows || rows_to <= 0 || rows_to > rows || rows_from >= rows_to) {
        throw runtime_error("rows_from o rows_to inválidos");
    }
    if (cols_from < 0 || cols_from >= cols || cols_to <= 0 || cols_to > cols || cols_from >= cols_to) {
        throw runtime_error("cols_from o cols_to inválidos");
    }

    Matrix<T> res(rows_to - rows_from, cols_to - cols_from);

    for (int i = 0; i < res.height(); i++) {
        for (int j = 0; j < res.width(); j++) {
            res[i][j] = M[rows_from + i][cols_from + j];
        }
    }

    return res;
}

template <class T>
void Matrix<T>::save(const string& filename) const {
    ofstream f;
    f.open(filename);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            f << M[i][j] << " ";
        }
        f << "\n";
    }
    f.close();
}
