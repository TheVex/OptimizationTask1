#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

class Matrix {
protected:
    int rows{};
    int columns{};
    vector<vector<double>> matrix;
public:
    Matrix(int n, int m) {
        this->rows = n;
        this->columns = m;
    }

    void matrixInputReader() {
        for (int i = 0; i < rows; i++) {
            matrix.emplace_back();
            for (int j = 0; j < columns; j++) {
                double elem;
                cin >> elem;
                matrix[i].push_back(elem);
            }
        }
    }

    explicit Matrix(const vector<vector<double>>& matrix) {
        this->rows = matrix.size();
        if (this->rows == 0) {
            this->columns = 0;
        }
        else {
            this->columns = matrix[0].size();
        }

        this->matrix = matrix;
    }

    Matrix() = default;

    Matrix& operator=(const Matrix& other) {
        this->rows = other.rows;
        this->columns = other.columns;
        this->matrix = other.matrix;
        return *this;
    }

    Matrix operator+(const Matrix& other) {
        if (this->rows != other.rows || this->columns != other.columns) {
            cout << "Error: the dimensional problem occurred" << endl;
            return {};
        }
        Matrix newMatrix = *this;

        for (int i = 0; i < other.rows; i++) {
            for (int j = 0; j < other.columns; j++) {
                newMatrix.matrix[i][j] += other.matrix[i][j];
            }
        }
        return newMatrix;
    }

    Matrix operator-(const Matrix& other) {
        if (this->rows != other.rows || this->columns != other.columns) {
            cout << "Error: the dimensional problem occurred" << endl;
            return {};
        }
        Matrix newMatrix = *this;
        for (int i = 0; i < other.rows; i++) {
            for (int j = 0; j < other.columns; j++) {
                newMatrix.matrix[i][j] -= other.matrix[i][j];
            }
        }
        return newMatrix;
    }

    Matrix operator*(const Matrix& other) {
        if (this->columns != other.rows) {
            cout << "Error: the dimensional problem occurred" << endl;
            return {};
        }
        Matrix newMatrix(this->rows, other.columns);
        for (int i = 0; i < newMatrix.rows; i++) {
            newMatrix.matrix.emplace_back();
            for (int j = 0; j < newMatrix.columns; j++) {
                newMatrix.matrix[i].push_back(0);
            }
        }
        for (int i = 0; i < newMatrix.rows; i++) {
            for (int j = 0; j < newMatrix.columns; j++) {
                for (int k = 0; k < this->columns; k++) {
                    newMatrix.matrix[i][j] += this->matrix[i][k] * other.matrix[k][j];
                }
            }
        }
        return newMatrix;
    }

    bool operator==(Matrix other) {
        if (this->rows != other.rows || this->columns != other.columns) {
            return false;
        }
        for (int i = 0; i < this->rows; i++) {
            for (int j = 0; j < this->columns; j++) {
                if (matrix[i][j] != other.matrix[i][j]) {
                    return false;
                }
            }
        }
        return true;
    }

    Matrix transpose() {
        Matrix newMatrix(rows, columns);
        for (int i = 0; i < columns; i++) {
            newMatrix.matrix.emplace_back();
            for (int j = 0; j < rows; j++) {
                newMatrix.matrix[i].push_back(matrix[j][i]);
            }
        }
        return newMatrix;
    }

    friend std::ostream& operator<<(ostream& out, const Matrix& matx)
    {
        for (vector<double> row : matx.matrix) {
            for (int i = 0; i < row.size(); i++) {
                cout << row[i];
                if (i != row.size() - 1) {
                    out << " ";
                }
            }
            out << endl;
        }
        return out;
    }

    friend istream& operator>>(istream& in, Matrix& matx)
    {
        for (int i = 0; i < matx.rows; i++) {
            for (int j = 0; j < matx.columns; j++) {
                in >> matx.matrix[i][j];
            }
        }
        return in;
    }

    int getRow() const {
        return this->rows;
    }

    int getColumn() const {
        return this->columns;
    }

    vector<vector<double>> getMatrix() {
        return this->matrix;
    }

    Matrix clone() {
        return Matrix(this->matrix);
    }

    void setMatrixCell(int i, int j, double value) {
        matrix[i][j] = value;
    }

    double getMatrixCell(int i, int j) {
        return matrix[i][j];
    }

};

Matrix solveMaximizationProblem(Matrix C, Matrix A, Matrix b, double epsilon) {
    // This vector will contain iteration table
    vector<vector<double>> it_vector;
    it_vector.emplace_back();

    // Copying values from task to z-row
    for (int i = 0; i < A.getColumn(); i++) {  
        it_vector[0].push_back(-C.getMatrixCell(0, i));
    }
    // Setting solution to 0
    it_vector[0].push_back(0);

    // Completing table with values from constraints and right-hand side vector
    for (int i = 1; i <= A.getRow(); i++) {
        it_vector.emplace_back();
        for (int j = 0; j <= A.getColumn(); j++) {
            if (j == A.getColumn()) {
                it_vector[i].push_back(b.getMatrixCell(i - 1, 0));
            }
            else {
                it_vector[i].push_back(A.getMatrixCell(i - 1, j));
            }
            
        }
    }

    // 0 Iteration table
    Matrix Iteration(it_vector);

    while (true) {
        Matrix NextIteration = Iteration.clone();

        double minZValue = 0;
        int minColumnIndex = 999;

        double minRow = 999999;
        double minRowDelimeter = 999999;
        int minRowIndex = 999;

        // Looking for minimal value for next iteration
        for (int i = 0; i < Iteration.getColumn() - 1; i++) {
            int t = Iteration.getMatrixCell(0, i);
            if (t < 0 && t < minZValue) {
                minZValue = t;
                minColumnIndex = i;
            }
        }
        // If no negative values, then we have done
        if (minZValue >= 0) {
            return Iteration;
        }


        // Defining row with minimal result of function
        for (int i = 1; i < Iteration.getRow(); i++) {
            double v = Iteration.getMatrixCell(i, Iteration.getColumn() - 1) / Iteration.getMatrixCell(i, minColumnIndex);
            if (v < minRow) {
                minRow = v;
                minRowDelimeter = Iteration.getMatrixCell(i, minColumnIndex);
                minRowIndex = i;
            }
        }

        // Dividing new row by its delimeter
        for (int i = 0; i < NextIteration.getColumn(); i++) {
            NextIteration.setMatrixCell(minRowIndex, i, NextIteration.getMatrixCell(minRowIndex, i) / minRowDelimeter);
        }

        // Completing the table
        for (int i = 0; i < NextIteration.getRow(); i++) {
            for (int j = 0; j < NextIteration.getColumn(); j++) {
                if (i == minRowIndex) continue;
                double newValue = Iteration.getMatrixCell(i, j) - Iteration.getMatrixCell(i, minColumnIndex) * NextIteration.getMatrixCell(minRowIndex, j);
                NextIteration.setMatrixCell(i, j, newValue);
            }
        }
        cout << "Debug: " << NextIteration << '\n';
        Iteration = NextIteration.clone();
    }
}


int main() {
    int C_size;
    cin >> C_size;
    Matrix C(1, C_size);
    C.matrixInputReader();

    int A_width, A_length;
    cin >> A_width >> A_length;
    Matrix A(A_width, A_length);
    A.matrixInputReader();

    int b_size;
    cin >> b_size;
    Matrix b(b_size, 1);
    b.matrixInputReader();

    double epsilon;
    cin >> epsilon;

    Matrix ans = solveMaximizationProblem(C, A, b, epsilon);

    cout << "Solution vector is: (";
    for (int i = 0; i < ans.getColumn() - 1; i++) {
        cout << ans.getMatrixCell(0, i);
        if (i < ans.getColumn() - 2) {
            cout << ", ";
        }
    }
    cout << ")\n";
    cout << "Maximum value is: " << ans.getMatrixCell(0, ans.getColumn() - 1);
}
