#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#define epsilon 1e-8

using namespace std;

template <typename T>
void swap(vector<vector<T>> &vec, int i, int j) {
    vector<T> temp = vec[i];
    vec[i] = vec[j];
    vec[j] = temp;
}

template <typename T>
class Matrix {
private:
    int rows;
    int columns;
public:
    vector<vector<T>> matrix;

    Matrix() = default;
    Matrix(int _rows, int _columns) : rows(_rows), columns(_columns) {}


    void setRows(int rows)       { this->rows = rows; } 
    int getRows()const           { return rows; } 

    void setColumns(int columns) { this->columns = columns; }
    int getColumns()const        { return columns; }

    void operator=(const Matrix& rhs_matrix) {
        vector<vector<T>> matrix_data;
        vector<T>         vec;

        for (int i = 0; i < rhs_matrix.getRows(); ++i) {
            for (int j = 0; j < rhs_matrix.getColumns(); ++j) {
                vec.push_back(rhs_matrix.matrix[i][j]);
            }
            matrix_data.push_back(vec);
            vec.clear();
        }
        this->matrix = matrix_data;
    }

    Matrix transposed() {
        Matrix matrix_transposed = Matrix<T>(this->getColumns(), this->getRows());
        vector<vector<T>> matrix_data;
        vector<T>         vec;

        for (int i = 0; i < this->getColumns(); ++i) {
            for (int j = 0; j < this->getRows(); ++j) {
                vec.push_back(this->matrix[j][i]);
            }
            matrix_data.push_back(vec);
            vec.clear();
        }
        matrix_transposed.matrix = matrix_data;

        return matrix_transposed;
    }

    Matrix operator+(const Matrix& rhs_matrix) {
          Matrix return_matrix = Matrix<T>(this->getRows(), this->getColumns());
          vector<vector<T>> matrix_data;
          vector<T>         vec;

          for (int i = 0; i < this->getRows(); ++i) {
              for (int j = 0; j < this->getColumns(); ++j) {
                  vec.push_back(this->matrix[i][j] + rhs_matrix.matrix[i][j]);
              }
              matrix_data.push_back(vec);
              vec.clear();
          }
          return_matrix.matrix = matrix_data;

          return return_matrix;
    }

    Matrix operator-(const Matrix& rhs_matrix) {
        Matrix return_matrix = Matrix<T>(this->getRows(), this->getColumns());
        vector<vector<T>> matrix_data;
        vector<T>         vec;

        for (int i = 0; i < this->getRows(); ++i) {
            for (int j = 0; j < this->getColumns(); ++j) {
                vec.push_back(this->matrix[i][j] - rhs_matrix.matrix[i][j]);
            }
            matrix_data.push_back(vec);
            vec.clear();
        }
        return_matrix.matrix = matrix_data;

        return return_matrix;
    }

    Matrix operator*(const Matrix& rhs_matrix) {
        Matrix return_matrix = Matrix<T>(this->getRows(), rhs_matrix.getColumns());
        vector<vector<T>> matrix_data;
        vector<T>         vec;

        for (int i = 0; i < this->getRows(); ++i) {
            for (int j = 0; j < rhs_matrix.getColumns(); ++j) {
                vec.push_back(0);
                for (int k = 0; k < this->getColumns(); ++k) {
                    vec[j] += this->matrix[i][k] * rhs_matrix.matrix[k][j];
                }
            }
            matrix_data.push_back(vec);
            vec.clear();
        }

        return_matrix.matrix = matrix_data;

        return return_matrix;

    }

    friend istream& operator>>(istream& in, Matrix& this_matrix) {
        vector<T> vec;
        T         x;

        for (int i = 0; i < this_matrix.getRows(); ++i) {
            for (int j = 0; j < this_matrix.getColumns(); ++j) {
                in >> x;
                vec.push_back(x);
            }
            this_matrix.matrix.push_back(vec);
            vec.clear();
        }

        return in;
    }

    friend ostream& operator<<(ostream& out, Matrix& this_matrix) {
        out << fixed;

        for (int i = 0; i < this_matrix.getRows(); ++i) {
            for (int j = 0; j < this_matrix.getColumns(); ++j) {
                if (j != this_matrix.getColumns() - 1) {
                    if (fabs(this_matrix.matrix[i][j]) < epsilon) out << setprecision(4) << 0.0000 << " ";
                    else out << setprecision(4) << this_matrix.matrix[i][j] << " ";
                }
                else {
                    if (fabs(this_matrix.matrix[i][j]) < epsilon) out << setprecision(4) << 0.0000;
                    else out << setprecision(4) << this_matrix.matrix[i][j];
                }
            }
            out << "\n";
        }

        return out;
    }
};

template <typename T>
class SquareMatrix : public Matrix<T> {
private:
    T determinant;
public:
    SquareMatrix(int dim) : Matrix<T> (dim, dim) {}

    void calculateDeterminant() {
        double determinant { 1 };

        toUpperTriangular(this->getRows(), 0, 1, 0, determinant, *this, false);

        for (int i = 0; i < this->getRows(); ++i) {
            determinant *= this->matrix[i][i];
            this->determinant = determinant;
        }
    }
    
    T det() {
        return this->determinant;
    }
};

class IdentityMatrix : public SquareMatrix<double> {
public:
    IdentityMatrix(int dim) : SquareMatrix<double> (dim) {
        vector<double> vec;

        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                if (i == j) vec.push_back(1);
                else vec.push_back(0);
            }
            this->matrix.push_back(vec);
            vec.clear();
        }
    }
};

class EliminationMatrix: public IdentityMatrix {
public:
    EliminationMatrix(Matrix<double> &source_matrix, int dim, int i, int j) : IdentityMatrix(dim) {
        this->matrix[i][j] = -(source_matrix.matrix[i][j] / source_matrix.matrix[j][j]);
    }
};

class PermutationMatrix : public IdentityMatrix {
public:
    PermutationMatrix(int dim, int i, int j) : IdentityMatrix(dim) {
        swap<double>(this->matrix, i, j);
    }
};


template <typename T>
bool upperTriangular(Matrix<T> &source_matrix, int dim) {
    for (int i = 1; i < dim; ++i) {
        for (int j = 0; j < i; ++j) {
            if (fabs(source_matrix.matrix[i][j]) > epsilon) return false;
        }
    }
    return true;
}
template <typename T>
bool lowerTriangular(Matrix<T> &source_matrix, int dim) {
    for (int i = 0; i < dim; ++i) {
        for (int j = i + 1; j < dim; ++j) {
            if (fabs(source_matrix.matrix[i][j]) > epsilon) return false;
        }
    }
    return true;
}

template <typename T>
PermutationMatrix generatePermutation(Matrix<T> &source_matrix, int dim, int row,
                                        int column, T max_val) {
    for (int i = 0; i < dim; ++i) {
        if (source_matrix.matrix[i][column] == max_val) {
            return PermutationMatrix(dim, row, i);
        }
    }
}

template <typename T>
T findPivot(Matrix<T> &a, int dim, int column) {
    T max_val = a.matrix[column][column];
    for (int i = column; i < dim; ++i) {
        if (abs(a.matrix[i][column]) > abs(max_val)) max_val = a.matrix[i][column];
    }
    return max_val;
}

template <typename T>
class ColumnVector : public Matrix<T> {
public:
    ColumnVector (int dim) : Matrix<T> (dim, 1) {}
    ColumnVector operator+(ColumnVector &rhs_column) {
        ColumnVector<T> temp_column = ColumnVector<T>(this->getRows());
        vector<T> vec;
        for (int i = 0; i < this->getRows(); ++i) {
            vec.push_back(this->matrix[i][0] + rhs_column.matrix[i][0]);
            temp_column.matrix.push_back(vec);
            vec.clear();
        }
        return temp_column;
    }
    ColumnVector operator-(ColumnVector &rhs_column) {
        ColumnVector<T> temp_column = ColumnVector<T>(this->getRows());
        vector<T> vec;
        for (int i = 0; i < this->getRows(); ++i) {
            vec.push_back(this->matrix[i][0] - rhs_column.matrix[i][0]);
            temp_column.matrix.push_back(vec);
            vec.clear();
        }
        return temp_column;
    }

    ColumnVector operator*(Matrix<T> &source_matrix) {
        ColumnVector<T> temp_column = ColumnVector<T>(source_matrix.getRows(), 1);
        vector<T> vec;
        for (int i = 0; i < source_matrix.getRows(); ++i) {
            vec.push_back(0);
            for (int j = 0; j < source_matrix.getRows(); ++j) {
                vec[i] += source_matrix[i][j] * this->matrix[j][0];
            }
            temp_column.matrix.push_back(vec);
            vec.clear;
        }
        return temp_column;
    }
    ColumnVector operator=(Matrix<T> &source_matrix) {
        ColumnVector<T> temp_column = ColumnVector<T>(source_matrix.getRows(), 1);
        vector<T> vec;
        for (int i = 0; i < source_matrix.getRows(); ++i) {
            vec.push_back(source_matrix.matrix[i][0]);
            temp_column.matrix.push_back(vec);
            vec.clear();
        }
        return temp_column;
    }

    T norm() {
        T summ_squares = 0;
        for (int i = 0; i < this->getRows(); ++i) {
            summ_squares += this->matrix[i][0] * this->matrix[i][0];
        }
        return sqrt(summ_squares);
    }

};

void toUpperTriangular(int n, int &step, int row, int column, double &determinant, Matrix<double> &a, bool verbose) {
    double max_val;
    Matrix<double> e = SquareMatrix<double> (n);
    Matrix<double> p = SquareMatrix<double> (n);
    while (!upperTriangular(a, n)) {
        if (row == n) {
            column++;
            row = column + 1;
        }
        max_val = findPivot(a, n, column);

        if (max_val == a.matrix[column][column]) {
            if (fabs(max_val) > epsilon) {
                if (fabs(a.matrix[row][column]) > epsilon) {
                    e = EliminationMatrix(a, n, row, column);
                    a = e * a;
                    if (verbose) {
                        cout << "step #" << step << ": elimination\n";
                        cout << a;
                    }
                } else step--;
                row++;
            } else step--;
        } else {
            p = generatePermutation(a, n, row-1, column, max_val);
            a = p * a;
            if (verbose) {
                cout << "step #" << step << ": permutation\n";
                cout << a;
            }
            determinant *= -1;
        }
        step++;
    }
}

void toLowerTriangular(int n, int &step, int row, int column, double &determinant, Matrix<double> &a, bool verbose) {
    double max_val;
    Matrix<double> e = SquareMatrix<double> (n);
    while (!lowerTriangular(a, n)) {
        if (row == (-1)) {
            column--;
            row = column - 1;
        }

        if (fabs(a.matrix[column][column]) > epsilon) {
            if (fabs(a.matrix[row][column]) > epsilon) {
                e = EliminationMatrix(a, n, row, column);
                a = e * a;
                if (verbose) {
                    cout << "step #" << step << ": elimination\n";
                    cout << a;
                }
            } else step--;
            row--;
        } else step--;
        step++;
    }
}
template <typename T>
class AugmentedMatrix : public Matrix<T> {
public:
    AugmentedMatrix(Matrix<T> &source_matrix, int dim) : Matrix<T> (dim, 2*dim) {
        vector<T> vec;
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < 2 * dim; ++j) {
                if (j < dim) {
                    vec.push_back(source_matrix.matrix[i][j]);
                } else {
                    if (i == (j % dim)) vec.push_back(1);
                    else vec.push_back(0);
                }
            }
            this->matrix.push_back(vec);
            vec.clear();
        }
    }
    SquareMatrix<T> inverse() {
        SquareMatrix<T> inverse_matrix = SquareMatrix<T>(this->getRows());
        vector<T> vec;
        for (int i = 0; i < this->getRows(); ++i) {
            for (int j = this->getRows(); j < this->getColumns(); ++j) {
                vec.push_back(this->matrix[i][j]);
            }
            inverse_matrix.matrix.push_back(vec);
            vec.clear();
        }
        return inverse_matrix;
    }
};

template <typename T>
void diagonalNormalization(Matrix<T> &source_matrix) {
    for (int i = 0; i < source_matrix.getRows(); ++i) {
        for (int j = source_matrix.getRows(); j < source_matrix.getColumns(); ++j) {
            source_matrix.matrix[i][j] /= source_matrix.matrix[i][i];
        }
        source_matrix.matrix[i][i] = 1;
    }
}

template <typename T>
bool JacobiApplicable(Matrix<T> &source_matrix) {
    T summ;
    for (int i = 0; i < source_matrix.getRows(); ++i) {
        summ = 0;
        for (int j = 0; j < source_matrix.getColumns(); ++j) {
            if (i != j) summ += fabs(source_matrix.matrix[i][j]);
        }
        if (summ > fabs(source_matrix.matrix[i][i])) return false;
    }
    return true;
}

template <typename T>
Matrix<T> constructMatrix(Matrix<T> &source_matrix) {
    Matrix<T> output = Matrix<T> (source_matrix.getRows(), source_matrix.getColumns());
    vector<T> vec;
    for (int i = 0; i < source_matrix.getRows(); ++i) {
        for (int j = 0; j < source_matrix.getColumns(); ++j) {
            if (i == j) vec.push_back(0);
            else vec.push_back(-(source_matrix.matrix[i][j] / source_matrix.matrix[i][i]));
        }
        output.matrix.push_back(vec);
        vec.clear();
    }
    return output;
}

template <typename T>
ColumnVector<T> constructColumn(ColumnVector<T> &source_vector, Matrix<T> &source_matrix) {
    ColumnVector<T> output = ColumnVector<T> (source_vector.getRows());
    vector<T> vec;
    for (int i = 0; i < source_matrix.getRows(); ++i) {
        vec.push_back(source_vector.matrix[i][0] / source_matrix.matrix[i][i]);
        output.matrix.push_back(vec);
        vec.clear();
    }

    return output;
}

template <typename T>
Matrix<T> constructLowerMatrix(Matrix<T> &source_matrix) {
    Matrix<T> mat = Matrix<T>(source_matrix.getRows(), source_matrix.getColumns());
    vector<T> vec;
    for (int i = 0; i < source_matrix.getRows(); ++i) {
        for (int j = 0; j < source_matrix.getColumns(); ++j) {
            if (i > j) vec.push_back(source_matrix.matrix[i][j]);
            else vec.push_back(0.0000);
        }
        mat.matrix.push_back(vec);
        vec.clear();
    }
    return mat;
}

template <typename T>
Matrix<T> constructUpperMatrix(Matrix<T> &source_matrix) {
    Matrix<T> mat = Matrix<T>(source_matrix.getRows(), source_matrix.getColumns());
    vector<T> vec;
    for (int i = 0; i < source_matrix.getRows(); ++i) {
        for (int j = 0; j < source_matrix.getColumns(); ++j) {
            if (i < j) vec.push_back(source_matrix.matrix[i][j]);
            else vec.push_back(0.0000);
        }
        mat.matrix.push_back(vec);
        vec.clear();
    }
    return mat;
}

int main() {
    return 0;
}
