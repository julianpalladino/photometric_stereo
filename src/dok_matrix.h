#pragma once
#include "common.h"
#include <map>
#include <stdexcept>

using namespace std;



template <class T>
class DOKMatrix {
  public:
    DOKMatrix(int rows, int cols) :
        _rows(rows),
        _cols(cols) {};

    DOKMatrix(const DOKMatrix<T>& rhs) :
        _rows(rhs._rows),
        _cols(rhs._cols),
        _values(rhs._values) {};

    bool operator==(const DOKMatrix<T>& rhs) const;
    bool eq_aprox(const DOKMatrix<T>& rhs) const;

    inline int cols() const { return _cols; };
    inline int rows() const { return _rows; };

    // deprecated
    inline int width() const  { return cols(); };
    inline int height() const { return rows(); };

    T& operator()(int i, int j);
    const T operator()(int i, int j) const;
    const T at(int row, int col) const;

    DOKMatrix<T> transposed() const;
    DOKMatrix<T> mult(const DOKMatrix<T>& M) const;
    DOKMatrix<T> mult_transposed() const;
    DOKMatrix<T> mult_vector(const Matrix<T>& v) const;

    BandMatrix<T> to_band() const;

    void save(const string& filename) const;

    bool has_null_rows() const;

    template <typename U>
    friend ostream& operator<<(ostream& stream, const DOKMatrix<U>& M);

  private:
    int _rows;
    int _cols;
    map<int, map<int, T>> _values;

    bool dict_eq_aprox(const map<int, map<int, T>>& lhs_values,
                       const map<int, map<int, T>>& rhs_values) const;
};

template <class T>
bool DOKMatrix<T>::operator==(const DOKMatrix<T>& rhs) const {
    return _rows == rhs._rows && \
           _cols == rhs._cols && \
           _values == rhs._values;
}

template <class T>
bool DOKMatrix<T>::eq_aprox(const DOKMatrix<T>& rhs) const {
    return _rows == rhs._rows && \
           _cols == rhs._cols && \
           dict_eq_aprox(_values, rhs._values);
}

template <class T>
inline bool DOKMatrix<T>::dict_eq_aprox(
        const map<int, map<int, T>>& lhs_values,
        const map<int, map<int, T>>& rhs_values) const
{
    if (lhs_values.size() != rhs_values.size()) {
        return false;
    }

    for (auto& row_kv: lhs_values) {
        auto rhs_row_kv = rhs_values.find(row_kv.first);
        if (rhs_row_kv == rhs_values.end()) return false;

        for (auto& col_kv: row_kv.second) {
            auto rhs_col_kv = rhs_row_kv->second.find(col_kv.first);
            if (rhs_row_kv == rhs_values.end()) return false;

            if (abs(col_kv.second - rhs_col_kv->second) >= eps) return false;
        }
    }

    return true;
}

template <class T>
T& DOKMatrix<T>::operator()(int i, int j) {
    return _values[i][j];
}

template <class T>
const T DOKMatrix<T>::operator()(int i, int j) const {
    auto row_kv = _values.find(i);
    if (row_kv == _values.end()) return T();
    auto col_kv = row_kv->second.find(j);
    if (col_kv == row_kv->second.end()) return T();
    return col_kv->second;
}

template <class T>
const T DOKMatrix<T>::at(int i, int j) const {
    return (*this)(i, j);
}

template<class T>
DOKMatrix<T> DOKMatrix<T>::transposed() const {
    DOKMatrix<T> res(_cols, _rows);

    for (auto& row_kv: _values) {
        for (auto& col_kv: row_kv.second) {
            auto i = row_kv.first;
            auto j = col_kv.first;
            if (col_kv.second != 0) {
                res(j, i) = col_kv.second;
            }
        }
    }

    return res;
}

template<class T>
DOKMatrix<T> DOKMatrix<T>::mult(const DOKMatrix<T>& M) const {
    DOKMatrix<T> res(_cols, M._rows);

    /*
    for (auto& A_row: _values) {
        for (auto& B_row: M._values) {
            for (auto& B_col: B_row.second) {
                for (auto& A_col: A_row.second) {
                    auto i = A_row.first;
                    auto j = B_col.first;
                    auto k = A_col.first;
                    res(i, j) = this->at(i, k) * M(k, j);
                }
            }
        }
    }
    */

    //map<int, map<int, T>> _values;
    for (auto& A_row: _values) {
        auto i = A_row.first;
        cout << i << endl;
        for (int j = 0; j < M._cols; j++) {
            auto sum = T();
            for (auto& A_col: A_row.second) {
                auto k = A_col.first;
                sum += this->at(i, k) * M(k, j);
            }
            if (sum != T()) {
                res(i, j) = sum;
            }
        }
    }

    return res;
}

template<class T>
DOKMatrix<T> DOKMatrix<T>::mult_transposed() const {
    DOKMatrix<T> res(_cols, _cols);
    const DOKMatrix<T> Mt = transposed();

    vector<pair<int, map<int, T>>> list;
    list.reserve(Mt._values.size());
    for (auto& fila_de_Mt: Mt._values) {
        list.push_back(fila_de_Mt);
    }

    for (auto& fila_Mt: list) {
        auto i = fila_Mt.first;
        for (auto& col_M: list) {
            auto j = col_M.first;
            auto& fila_j = col_M.second;
            auto sum = T();
            for (auto& ik: fila_Mt.second) {
                auto k = ik.first;
                auto jk_it = fila_j.find(k);
                if (jk_it != fila_j.end()) {
                    sum += ik.second * jk_it->second;
                }
            }
            if (abs(sum) > eps) {
                res(i, j) = sum;
            }
        }
    }
    return res;
}

template<class T>
DOKMatrix<T> DOKMatrix<T>::mult_vector(const Matrix<T>& v) const {
    if (v.width() != 1) {
        throw runtime_error("matrix must be a vector (n x 1)");
    }
    if (width() != v.height()) {
        throw runtime_error("can't multiply with vector");
    }

    DOKMatrix<T> res(height(), 1);

    for (auto& row: _values) {
        auto i = row.first;
        auto sum = T();
        for (auto& col: row.second) {
            auto j = col.first;
            sum += this->at(i, j) * v(j, 0);
        }
        if (sum != 0) {
            res(i, 0) = sum;
        }
    }

    return res;
}

// Save as full matrix
template <class T>
void DOKMatrix<T>::save(const string& filename) const {
    ofstream f;
    f.open(filename);
    for (int i = 0; i < _rows; i++) {
        for (int j = 0; j < _cols; j++) {
            f << this->at(i, j) << " ";
        }
        f << "\n";
    }
    f.close();
}

template<class T>
ostream& operator<<(ostream& out, const DOKMatrix<T>& M) {
    for (auto& row: M._values) {
        for (auto& col: row.second) {
            out << row.first << " "<< col.first << ": " << col.second << "\t";
        }
        out << "\n";
    }
    return out;
}

template <class T>
BandMatrix<T> DOKMatrix<T>::to_band() const {
    auto res = BandMatrix<T>(_rows, _cols);
    for (auto& row: _values) {
        auto i = row.first;
        for (auto& col: row.second) {
            auto j = col.first;
            if (col.second != 0) {
                res(i, j) = col.second;
            }
        }
    }
    return res;
}

template <class T>
bool DOKMatrix<T>::has_null_rows() const {
    for (const auto& row: _values) {
        bool zero = true;
        for (const auto& col: row.second) {
            if (abs(col.second) != 0) {
                zero = false;
                break;
            }
        }
        if (zero) {
            return true;
        }
    }
    return false;
}
