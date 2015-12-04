#include "matrix.h"

void swap(int& i, int& j) { int tmp = i; i = j; j = tmp; }
void swap(double& i, double& j) { double tmp = i; i = j; j = tmp; }
double max(int i, int j) { return i > j ? i : j; }
double min(int i, int j) { return i < j ? i : j; }

void* threadedInverse(void* args) {
    ARGS* pargs = (ARGS*) args;
    Matrix* self = pargs->self;
    Matrix* inv = pargs->inv;
    int ncurr = pargs->ncurr;
    int nthrd = pargs->nthrd;
    int n = pargs->self->size();
    int begin = ncurr * (n / nthrd);
    int end = ncurr == nthrd - 1 ? n : (ncurr + 1) * (n / nthrd);
    double* temp = new double[2 * n];
    double* qtemp = new double[2 * n];
    string s, c;
    if (ncurr == 0)
        cout << "\tQR decomposition:" << endl;
    static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
    for (int i = 0; i < n - 1; ++i) {
        if (ncurr == 0 )
            cout << '\t' << '[' + s.assign(50 * i / (n - 1), '#') + c.assign(50 - 50 * i / (n - 1), ' ') + ']' +
                 ' ' + to_string(((1 + i) * 100) / (n - 1)) + '%' + '\r' << flush;
        pthread_mutex_lock (&mutex);
        if (self->flags[i] == 0) {
            self->flags[i] = 1;
            pthread_mutex_unlock (&mutex);
            for (int j = i + 1; j < n; ++j) {
                self->request(i, j, ncurr, nthrd);
                self->rotate(i, j, *inv, temp, qtemp);
                self->report(j);
            }
        } else {
            pthread_mutex_unlock (&mutex);
        }
    }
    self->awake();
    synchronize(nthrd, &self->condvar);
    if (ncurr == 0) {
        cout << '\t' << '[' + s.assign(50, '#') + ']' << endl;
        cout << "\tReverse Gauss: " << endl;
    }
    delete[] temp;
    delete[] qtemp;
    synchronize(nthrd, NULL);
    gauss(*self, *inv, begin, end, nthrd);
    return 0;
}

void synchronize (int total_threads, pthread_cond_t* condvar = NULL) {
    static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
    static pthread_cond_t condvar_in = PTHREAD_COND_INITIALIZER;
    static pthread_cond_t condvar_out = PTHREAD_COND_INITIALIZER;
    static int threads_in = 0;
    static int threads_out = 0;
    pthread_mutex_lock (&mutex);
    threads_in++;
    if (threads_in == total_threads -1 && condvar != NULL) 
        pthread_cond_broadcast (condvar); 
    if (threads_in >= total_threads) {
        threads_out = 0;
        pthread_cond_broadcast (&condvar_in);
    } else {
        while (threads_in < total_threads) {
            pthread_cond_wait (&condvar_in, &mutex);
        }
    }
    threads_out ++;
    if (threads_out >= total_threads) {
        threads_in = 0;
        pthread_cond_broadcast (&condvar_out);
    } else {
        while (threads_out < total_threads) {
            pthread_cond_wait (&condvar_out, &mutex);
        }
    }
    pthread_mutex_unlock (&mutex);
}

void gauss(Matrix& self, Matrix& inv, int begin, int end, int nthrd) {
    int n = inv.size();
    string s, c;
    double quot = 0.0;
    for (int j = n - 1; j >= 0; --j) {
        for (int i = min(end - 1, j - 1); i >= begin; --i) {
            quot = self[i][j] / self[j][j];
            for (int k = 0;  k < n; ++k) {
                inv[i][k] -= quot * inv[j][k];
                self[i][k] -= quot * self[j][k];
            }
        }
        synchronize(nthrd, NULL);
        if (begin == 0) {
            cout << '\t' << '[' + s.assign(50 * (n - 1 - j) / n, '#') +
                 c.assign(50 - 50 * (n - 1 - j) / n, ' ') + ']' + ' ' +
                 to_string(((1 + (n - 1 - j)) * 100) / n) + '%' + '\r';
        }
    }
    if (begin == 0)
        cout << '\t' << '[' + s.assign(50, '#') + ']' << endl;
    double denum;
    for (int j = begin; j < end; ++j) {
        denum = self.matrix[j * n + j];
        for (int k = 0; k < n; ++k) {
            self.matrix[j * n + k] /= denum;
            inv.matrix[j * n + k] /= denum;
        }
    }
}

Matrix multiplicate(Matrix& a, Matrix& b) {
    try {
        if (a.size() != b.size())
            throw "Wrong matrix's sizes.";
    } catch (const char* message) {
        cerr << message << endl;
    }
    Matrix res(a.size());
    for (int i = 0; i < a.size(); ++i) {
        for (int j = 0; j < b.size(); ++j) {
            double value = 0;
            for (int k = 0; k < a.size(); ++k)
                value += a[i][k] * b[k][j];
            res[i][j] = value;
        }
    }
    return res;
}

void Matrix::transpose() {
    for (int i = 0; i < n; i++)
        for (int j = 0; j < i; ++j)
            swap(it[j][i], it[i][j]);
    return;
}

void Matrix::rotate(int k, int l, Matrix& q, double* temp, double* qtemp) {
    double denum = sqrt(it[k][k] * it[k][k] + it[l][k] * it[l][k]);
    try {
        if (fabs(denum) < eps)
            throw "Division by zero.";
    } catch (const char* message) {
        return;
    }
    try {
        if (k >= l)
            throw "Wrong values of k and l.";
    } catch (const char* message) {
        cerr << message << endl;
        exit(1);
    }
    double x = it[k][k] / denum;
    double y = -it[l][k] / denum;
    for (int i = 0; i < n; ++i) {
        temp[i] = x * it[k][i] - y * it[l][i];
        temp[n + i] = y * it[k][i] + x * it[l][i];
        qtemp[i] = x * q[k][i] - y * q[l][i];
        qtemp[n + i] = y * q[k][i] + x * q[l][i];
    }
    for (int i = 0; i < n; ++i) {
        it[k][i] = temp[i];
        it[l][i] = temp[n + i];
        q[k][i] = qtemp[i];
        q[l][i] = qtemp[n + i];
    }
}

void Matrix::inverse(int nthreads) {
    Matrix inv(n);
    inv.unit();
    processed.assign(n, 0);
    flags.assign(n - 1, 0);
    pthread_t* threads = new pthread_t[nthreads];
    cout << "Inversing matrix: " << endl;
    pthread_mutex_init(&mutex, 0);
    pthread_cond_init(&condvar, 0);
    ARGS* args = new ARGS[nthreads];
    for (int i = 0; i < nthreads; ++i) {
        args[i].set(this, &inv, i, nthreads);
        if (pthread_create(&threads[i], NULL, threadedInverse, (void*)&args[i])) {
            cerr << "Error creating thread." << endl;
            exit(1);
        }
    }
    for (int i = 0; i < nthreads; ++i)
        if (pthread_join(threads[i], NULL)) {
            cerr << "Can't wait thread." << endl;
            exit(1);
        }
    delete[] threads;
    delete[] args;
    *this = inv;
}

double residual(Matrix& A, Matrix& B) {
    double max = 0.0;
    Matrix R;
    R = A * B;
    Matrix E(A.size());
    E.unit();
    double res = 0.0;
    for (int i = 0; i < R.size(); ++i) {
        res = 0.0;
        for (int j = 0; j < R.size(); ++j) {
            res += fabs(R[i][j] - E[i][j]);
        }
        if (res > max)
            max = res;
    }
    return max;
}

void Matrix::print(int m) {
    for (int i = 0; i < size() && i < m; ++i) {
        for (int j = 0; j < size() && j < m; ++j) {
            cout.setf(ios::right);
            cout.width(10);
            cout << it[i][j] << '\t';
        }
        cout << endl;
    }
}
