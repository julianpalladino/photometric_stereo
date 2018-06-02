#pragma once
#include <algorithm>
#include <cstdlib>
#include <utility>
#include <vector>
#include <fstream>
#include "common.h"

using namespace std;

// Matrix band simétrica
template<class T>
class BandMatrix {
  public:
    BandMatrix(int rows, int cols, bool symmetric = false) :
        rows_(rows),
        cols_(cols),
        symmetric_(symmetric)
        {
            if (symmetric_ && rows_ != cols_) {
                throw runtime_error("symmetric matrix must be square");
            }
        };

    BandMatrix(const BandMatrix<T>& rhs) :
        rows_(rhs.rows_),
        cols_(rhs.cols_),
        symmetric_(rhs.symmetric_),
        lower_bands_(rhs.lower_bands_),
        upper_bands_(rhs.upper_bands_) {}

    bool operator==(const BandMatrix<T>& rhs) const;
    bool eq_aprox(const BandMatrix<T>& rhs) const;

    inline int cols() const { return cols_; };
    inline int rows() const { return rows_; };

    // deprecated
    inline int width() const  { return cols(); };
    inline int height() const { return rows(); };


    T& operator()(int i, int j);
    const T operator()(int i, int j) const;
    const T at(int row, int col) const;

    BandMatrix<T> transposed() const;
    BandMatrix<T> mult_transposed() const;
    BandMatrix<T> mult_vector(const Matrix<T>& v) const;

    int lower_bandwidth() const;
    int upper_bandwidth() const;

    size_t size() const;

    void save(const string& filename) const;

    void save_sparse(
            const string& data_filename,
            const string& rows_filename,
            const string& cols_filename) const;

    bool has_null_rows() const;

    template <typename U>
    friend ostream& operator<<(ostream& stream, const BandMatrix<U>& M);

  private:
    int rows_;
    int cols_;
    bool symmetric_;
    vector<vector<T>> lower_bands_;
    vector<vector<T>> upper_bands_;

    int diagonal(int row, int col) const;
    bool band_eq_aprox(const vector<vector<T>>& lhs_bands, const vector<vector<T>>& rhs_bands) const;
};

template<class T>
bool BandMatrix<T>::operator==(const BandMatrix<T>& rhs) const {
    return rows_ == rhs.rows_ && \
        cols_ == rhs.cols_ && \
        symmetric_ == rhs.symmetric_ && \
        lower_bands_ == rhs.lower_bands_ && \
        upper_bands_ == rhs.upper_bands_;
}

template<class T>
bool BandMatrix<T>::eq_aprox(const BandMatrix<T>& rhs) const {
    return rows_ == rhs.rows_ && \
        cols_ == rhs.cols_ && \
        symmetric_ == rhs.symmetric_ && \
        band_eq_aprox(lower_bands_, rhs.lower_bands_) && \
        band_eq_aprox(upper_bands_, rhs.upper_bands_);
}

template<class T>
T& BandMatrix<T>::operator()(int row, int col) {
    //assert(row < rows_);
    //assert(col < cols_);

    int diag = diagonal(row, col);
    bool use_lower_band = (diag <= 0 || symmetric_);
    vector<vector<T>>& bands = use_lower_band ? lower_bands_ : upper_bands_;
    int i = use_lower_band ? abs(diag) : diag - 1;

    if (i >= bands.size()) {
        bands.resize(i + 1);
    }

    auto j = min(row, col);
    vector<T>& band = bands[i];
    if (j >= band.size()) {
        band.resize(j + 1);
    }

    return band[j];
}

template<class T>
const T BandMatrix<T>::at(int row, int col) const {
    //assert(row < rows_);
    //assert(col < cols_);

    int diag = diagonal(row, col);
    bool use_lower_band = (diag <= 0 || symmetric_);
    const vector<vector<T>>& bands = use_lower_band ? lower_bands_ : upper_bands_;
    int i = use_lower_band ? abs(diag) : diag - 1;

    if (i >= bands.size()) {
        return T();
    }

    auto j = min(row, col);
    const vector<T>& band = bands[i];
    if (j >= band.size()) {
        return T();
    }

    return band[j];
}

template<class T>
const T BandMatrix<T>::operator()(int row, int col) const {
    return at(row, col);
}

template<class T>
BandMatrix<T> BandMatrix<T>::transposed() const {
    if (symmetric_) return (*this);

    BandMatrix<T> res(cols(), rows());

    // lower_bands_ also contain principal diagonal
    res.lower_bands_.resize(upper_bands_.size() + 1);
    res.upper_bands_.resize(max(0, (int)lower_bands_.size() - 1));

    if (lower_bands_.size() > 0) {
        res.lower_bands_[0] = lower_bands_[0];
    }
    for (int i = 1; i < lower_bands_.size(); i++) {
        res.upper_bands_[i-1] = lower_bands_[i];
    }
    for (int i = 0; i < upper_bands_.size(); i++) {
        res.lower_bands_[i+1] = upper_bands_[i];
    }
    return res;
}

template<class T>
static inline pair<int, int> col_range(const BandMatrix<T>& A, int row) {
    return make_pair(max(0, row - A.lower_bandwidth()),
                     min(A.cols() - 1, row + A.upper_bandwidth()));
}

template<class T>
static inline pair<int, int> row_range(const BandMatrix<T>& A, int col) {
    return make_pair(max(0, col - A.upper_bandwidth()),
                     min(A.rows() - 1, col + A.lower_bandwidth()));
}

template<class T>
BandMatrix<T> BandMatrix<T>::mult_transposed() const {
    BandMatrix<T> res(cols(), cols(), true);
    const BandMatrix<T>& A = *this;

    auto p = upper_bandwidth();
    auto q = lower_bandwidth();

    for (int i = 0; i < cols(); i++) {
        for (int j = i; j < min(i + p + q + 1, cols()); j++) {
            T sum = T();
    
            pair<int, int> range1 = row_range(A, i);
            int from1 = get<0>(range1);
            int to1 = get<1>(range1);
            
            pair<int, int> range2 = row_range(A, j);
            int from2 = get<0>(range2);
            int to2 = get<1>(range2);
            
            int k_from = min(from1, from2);
            int k_to = max(to1, to2);
            
            for (int k = k_from; k <= k_to; k++) {
                sum += A(k, i) * A(k, j);
            }

            // Sólo escribo si no es cero, así se mantiene lo más banda posible.
            if (sum != 0) {
                res(i, j) = sum;
            }
        }
    }

    return res;
}

template<class T>
BandMatrix<T> BandMatrix<T>::mult_vector(const Matrix<T>& v) const {
    if (v.width() != 1) {
        throw runtime_error("matrix must be a vector (n x 1)");
    }
    if (width() != v.height()) {
        throw runtime_error("can't multiply with vector");
    }

    BandMatrix<T> res(v.height(), 1);

    auto p = upper_bandwidth();
    auto q = lower_bandwidth();

    for (int i = 0; i < rows(); i++) {
        T sum = T();
        for (int j = max(0, i - q); j <= min(cols() - 1, i + p); j++) {
            sum += (*this)(i, j) * v(j, 0);
        }
        res(i, 0) = sum;
    }

    return res;
}

template<class T>
inline int BandMatrix<T>::lower_bandwidth() const {
    return lower_bands_.size() - 1;
}

template<class T>
inline int BandMatrix<T>::upper_bandwidth() const {
    return symmetric_ ? lower_bands_.size() - 1 : upper_bands_.size();
}

template<class T>
size_t BandMatrix<T>::size() const {
    size_t total = 0;
    for (auto it = lower_bands_.cbegin(); it != lower_bands_.cend(); it++) {
        total += it->size();
    }
    for (auto it = upper_bands_.cbegin(); it != upper_bands_.cend(); it++) {
        total += it->size();
    }
    return total;
}

// Save as full matrix
template <class T>
void BandMatrix<T>::save(const string& filename) const {
    ofstream f;
    f.open(filename);
    for (int i = 0; i < rows(); i++) {
        for (int j = 0; j < cols(); j++) {
            f << (*this)(i, j) << " ";
        }
        f << "\n";
    }
    f.close();
}

// Save sparse matrix in 3 files: data, row indices and column indices
// On Octave/Matlab, use: sparse(r, c, d)
//   where r = row indices, c = column indices, d = data
template<class T>
void BandMatrix<T>::save_sparse(
        const string& data_filename,
        const string& rows_filename,
        const string& cols_filename) const
{
    ofstream data_f, rows_f, cols_f;
    data_f.open(data_filename);
    rows_f.open(rows_filename);
    cols_f.open(cols_filename);
    for (int q = 0; q < lower_bands_.size(); q++) {
        auto& band = lower_bands_[q];
        for (int k = 0; k < band.size(); k++) {
            auto row = q + k;
            auto col = k;
            rows_f << row + 1 << "\n";
            cols_f << col + 1 << "\n";
            data_f << band[k] << "\n";
        }
    }
    // If matrix is symmetric, upper bands are the same as lower bands
    auto& upper = symmetric_ ? lower_bands_ : upper_bands_;
    // Also, lower_bands[0] is the principal diagonal, so skip it here:
    auto p_start = symmetric_ ? 1 : 0;
    for (int p = p_start; p < upper.size(); p++) {
        auto& band = upper[p];
        for (int k = 0; k < band.size(); k++) {
            auto row = k;
            auto col = p + k + 1;
            rows_f << row + 1 << "\n";
            cols_f << col + 1 << "\n";
            data_f << band[k] << "\n";
        }
    }
    data_f.close();
    rows_f.close();
    cols_f.close();
}

template<class T>
ostream& operator<<(ostream& out, const BandMatrix<T>& M) {
    for (size_t i = 0; i < M.rows(); i++) {
        for (size_t j = 0; j < M.cols(); j++) {
            out << M(i, j) << "\t";
        }
        out << "\n";
    }
    return out;
}

// Devuelve índice de la banda. Si es positivo, es una banda superior.
template<class T>
inline int BandMatrix<T>::diagonal(int row, int col) const {
    return col - row;
}

template<class T>
inline bool BandMatrix<T>::band_eq_aprox(const vector<vector<T>>& lhs_bands, const vector<vector<T>>& rhs_bands) const {
    if (lhs_bands.size() != rhs_bands.size()) {
        return false;
    }

    for (size_t i = 0; i < lhs_bands.size(); i++) {
        const vector<T>& lhs_band = lhs_bands[i];
        const vector<T>& rhs_band = rhs_bands[i];

        if (lhs_band.size() != rhs_band.size()) {
            return false;
        }

        for (size_t j = 0; j < lhs_band.size(); j++) {
            if (abs(lhs_band[j] - rhs_band[j]) >= eps) { // ¬(|a-b|<epsilon)
                return false;
            }
        }
    }

    return true;
}

template<class T>
bool BandMatrix<T>::has_null_rows() const {
    auto p = lower_bandwidth();
    auto q = upper_bandwidth();

    for (int i = 0; i < rows_; i++) {
        bool zero = true;
        for (int j = max(0, i - p - 1); j < min(cols_, i + q + 1); j++) {
            if (abs(this->at(i, j)) > 0) {
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
