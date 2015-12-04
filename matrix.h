#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <new>
#include <utility>
#include <pthread.h>
#include <vector>

#define eps 1.0e-20
#define it (*this)

using std::vector;
using std::ios;
using std::ifstream;
using std::istream;
using std::ofstream;
using std::ostream;
using std::string;
using std::to_string;
using std::cerr;
using std::cout;
using std::endl;
using std::flush;

class Matrix ;
void synchronize (int total_threads, pthread_cond_t* condvar) ;

class ARGS {
  public:
    Matrix* self;
    Matrix* inv;
    int ncurr;
    int nthrd;
    ARGS() {
        self = NULL;
        inv = NULL;
        ncurr = 0;
        nthrd = 0;
    }
    ~ARGS() {}
    void set(Matrix* s, Matrix* c, int nc, int nt) {
        self = s;
        inv = c;
        ncurr = nc;
        nthrd = nt;
    }
};

class Matrix {
  private:
    int n; // n is usually used as height
    double* matrix;
    void init(int N) {
        n = N;
        try {
            if (matrix != NULL) {
                delete[] matrix;
            }
            matrix = new double[n * n];
        } catch (std::bad_alloc& ba) {
            cerr << "Problems with memory: " << ba.what() << endl;
            exit(1);
        }
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                it[i][j] = 0.0;
            }
        }
    }
    void rotate(int k, int l, Matrix& Q, double* temp, double* qtemp);
  public:
    pthread_mutex_t mutex;
    pthread_cond_t condvar;
    vector<int> processed;
    vector<int> flags;
    friend void QRDecomp(Matrix* self, Matrix* inv,
                         int begin, int end, int nthrd);
    Matrix() {
        matrix = NULL;
        n = 0;
    };
    Matrix(int N) {
        matrix = NULL;
        n = N;
        init(n);
    }
    Matrix(const Matrix& M) {
        matrix = NULL;
        n = M.n;
        init(n);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                it[i][j] = M.matrix[i * n + j];
    }
    ~Matrix() {
        if (matrix != NULL)
            delete[] matrix;
    }
    int size() const { return n; }
    inline double* operator[] (int i) { return (matrix + i * n); }
    friend ifstream& operator >>(ifstream& in, Matrix& a) {
        try {
            in >> a.n;
            a.init(a.n);
            for (int i = 0; i < a.n; ++i) {
                for (int j = 0; j < a.n; ++j) {
                    in >> a[i][j];
                }
            }
        } catch (std::ifstream::failure e) {
            cerr << "Exception with input file." << endl;
            exit(1);
        }
        return in;
    }
    friend istream& operator >>(istream& in, Matrix& a) {
        try {
            in >> a.n;
            a.init(a.n);
            for (int i = 0; i < a.n; ++i) {
                for (int j = 0; j < a.n; ++j) {
                    in >> a[i][j];
                }
            }
        } catch (std::ifstream::failure e) {
            cerr << "Exception with input file." << endl;
            exit(1);
        }
        return in;
    } // rewrite fo rvalue!
    friend ofstream& operator <<(ofstream& out, Matrix& a) {
        try {
            for (int i = 0; i < a.n; ++i) {
                for (int j = 0; j < a.n; ++j) {
                    out.setf(ios::right);
                    out.width(10);
                    out << a[i][j] << '\t';
                }
                out << endl;
            }
        } catch (std::ofstream::failure e) {
            cerr << "Exception with output file." << endl;
            exit(1);
        }
        out << endl;
        return out;
    } // rewrite fo rvalue!
    friend ostream& operator <<(ostream& out, Matrix& a) {
        try {
            for (int i = 0; i < a.n; ++i) {
                for (int j = 0; j < a.n; ++j) {
                    out.setf(ios::right);
                    out.width(10);
                    out << a[i][j] << '\t';
                }
                out << endl;
            }
        } catch (std::ofstream::failure e) {
            cerr << "Exception with output file." << endl;
            exit(1);
        }
        out << endl;
        return out;
    } // rewrite fo rvalue!
    friend ostream& operator <<(ostream& out, Matrix&& a) {
        try {
            for (int i = 0; i < a.n; ++i) {
                for (int j = 0; j < a.n; ++j) {
                    out.setf(ios::right);
                    out.width(10);
                    out << a[i][j] << '\t';
                }
                out << endl;
            }
        } catch (std::ostream::failure e) {
            cerr << "Exception with output file." << endl;
            exit(1);
        }
        out << endl;
        return out;
    }
    friend ofstream& operator <<(ofstream& out, Matrix&& a) {
        try {
            for (int i = 0; i < a.n; ++i) {
                for (int j = 0; j < a.n; ++j) {
                    out.setf(ios::right);
                    out.width(10);
                    out << a[i][j] << '\t';
                }
                out << endl;
            }
        } catch (std::ofstream::failure e) {
            cerr << "Exception with output file." << endl;
            exit(1);
        }
        out << endl;
        return out;
    }
    void unit() {
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                if (i == j)
                    it[i][j] = 1.0;
                else
                    it[i][j] = 0.0;
    }
    Matrix& operator = (Matrix& a) {
        init(a.size());
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                it[i][j] = a[i][j];
        return it;
    }
    Matrix& operator = (Matrix&& a) {
        init(a.size());
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                it[i][j] = a[i][j];
        return it;
    }
    friend Matrix operator * (Matrix& a, Matrix& b) {
        Matrix res;
        res = multiplicate(a, b);
        return res;
    }
    bool isUnit() {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i == j) {
                    if (fabs(it[i][j] - 1) > eps)
                        return 0;
                } else if (fabs(it[i][j]) > eps)
                    return 0;
            }
        }
        return 1;
    }
    friend Matrix multiplicate(Matrix& a, Matrix& b);
    void transpose();
    void inverse(int nthreads);
    void print(int m);
    friend void* threadedInverse(void* args);
    friend void gauss(Matrix& self, Matrix& inv,
                      int begin, int end, int nthrd);
    friend class ARGS;
    friend double residual(Matrix& A, Matrix& B);
    friend double difference(Matrix& A, Matrix& B) {
        double diff = 0.0;
        int n = A.size();
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                diff += fabs(A[i][j] - B[i][j]);
        return diff;
    }
    void request(int k, int l, int ncurr, int nthrd) {
        static int in = 0;
        pthread_mutex_lock (&mutex);
        if (processed[l] >= k) {
            pthread_cond_broadcast (&condvar);
        } else {
            while (processed[l] < k) {
                in++;
                if (in >= nthrd) {
                    in = 1;
                    pthread_cond_broadcast (&condvar);
                }
                pthread_cond_wait (&condvar, &mutex);
            }
        }
        pthread_mutex_unlock (&mutex);
    }
    void report(int l) {
        //static pthread_mutex_t loc_mutex = PTHREAD_MUTEX_INITIALIZER;
        pthread_mutex_lock (&mutex);
        processed[l] += 1;
        pthread_mutex_unlock (&mutex);
    }
    void awake() {
        pthread_cond_broadcast (&condvar);
    }
};
