#include <sys/time.h>
#include "matrix.h"

using std::cin;
using std::cout;
using std::cerr;
using std::string;
using std::ifstream;
using std::endl;

double f(int i, int j, int n) {
    return (i > j) ? n - i : n - j;
    //return i + j == n - 1 ? 1.0 : 0.0;
}

Matrix generate(int n) {
    Matrix A(n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            A[i][j] = f(i, j, n);
    return A;
}

double get_time() {
    struct timeval tv;
    gettimeofday(&tv, 0);
    return tv.tv_usec / 1000000.0 + (double)tv.tv_sec;
}

int main(int argc, char* argv[]) {
    Matrix A, B;
    int n, nthrd;
    double start, end;
    cout << "Enter matrix size: ";
    cin >> n;
    cout << "Enter number of threads: ";
    cin >> nthrd;
    A = generate(n);
    B = A;
    start = get_time();
    A.inverse(nthrd);
    end = get_time();
    A.print(6);
    cout << endl;
    cout << "Spent time: " << end - start << endl;
    cout << "Residual: " << residual(A,B) << endl;
    return 0;
}
