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
    if (argc != 3) {
        cout << "Wrong number of args." << endl;
        return 0;
    }
    int n = atoi(argv[1]);
    int nthrd = atoi(argv[2]);
    if (n < 1 || nthrd < 1 || nthrd > n) {
        cout << "Wrong arg values." << endl;
        return 0;
    }
    cout << "Matrix size: " << n << endl;
    cout << "Number of threads: " << nthrd << endl;
    cout << "Creating matrices." << endl;
    Matrix A, B;
    A = generate(n);
    B = A;
    string s, c;
    double start, end;
    printf("\e[?25l"); /* hide the cursor */
    start = get_time();
    A.inverse(nthrd);
    end = get_time();
    cout << "Spent time: " << end - start << endl;
    cout << "Inversed matrix head: " << endl;
    A.print(6);
    cout << "Residual: " << residual(A, B);
    printf("\e[?25h"); /* show the cursor */
    return 0;
}
