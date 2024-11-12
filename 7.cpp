#include <iostream>
#include <iomanip>
#include <algorithm>
#include "..\..\Polinom\Polinom\Polynomial.h" // проект Polinom
#define W 18
using namespace std;

int n = 9; // узлы
double a = 0, b = 1;
double f(double x) { return sin(x); }
double p(double x) { return pow(x, -0.5); }

double integral_pl(Polynomial l) {
    vector<Monomial> m;
    double I = 0;
    m = l.GetMonoms();
    for (auto itr = m.begin(); itr != m.end(); ++itr) {
        double c = itr->GetKoef();
        int n = itr->GetStep();
        I += c / ((double)n + 0.5) * (pow(b, (double)n + 0.5) - pow(a, (double)n + 0.5));  // 0.5, потому что при интегрировании прибавляем ^1, но у нас p=1/sqrt(x)
    }
    return I;
}

int main() {
    setlocale(LC_ALL, "Russian");
    double* A, * x, * fx, h, I, I_calc;
    //cout << "Введите пределы интегрирования и количество узлов: "; 
    //cin >> a >> b >> n;
    h = (b - a) / (n - 1);
    A = new double[n];
    x = new double[n];
    fx = new double[n];

    cout << "Количество узлов: " << n << "\nПроверка для многочлена степени N-1: f = x^" << n - 1 << ".\n";
    for (int i = 0; i < n; ++i) {
        x[i] = a + h * i;
        fx[i] = pow(x[i], n - 1);
    }
    for (int k = 0; k < n; ++k) {
        Polynomial l(1, 0);
        vector<double> v{ 1,0 };
        for (int i = 0; i < k; ++i) {
            v[1] = -1 * x[i];
            Polynomial tmp(v);
            l = l * tmp * (1 / (x[k] - x[i]));
        }
        for (int i = k + 1; i < n; ++i) {
            v[1] = -x[i];
            Polynomial tmp(v);
            l = l * tmp * (1 / (x[k] - x[i]));
        }
        A[k] = integral_pl(l);
    }
    cout << setw(W) << "x_k" << setw(W) << "f_k" << setw(W) << "A_k\n";
    for (int i = 0; i < n; ++i) cout << setw(W) << x[i] << setw(W) << fx[i] << setw(W) << A[i] << "\n";
    cout << "Точное значение интеграла: " << setprecision(12) << (I = integral_pl(Polynomial(1, n - 1))) << "\n";
    I_calc = 0;
    for (int i = 0; i < n; ++i) I_calc += A[i] * pow(x[i], n - 1);
    cout << "Вычисленное значение интеграла: " << setprecision(12) << I_calc << "\n";
    cout << "Погрешность: " << setprecision(12) << fabs(I - I_calc) << "\n";


    cout << "\n\nВычисление интеграла функции sin(x)/sqrt(x):\n";
    for (int i = 0; i < n; ++i) fx[i] = f(x[i]); // x[i] остаются как были
    for (int k = 0; k < n; ++k) {
        Polynomial l(1, 0);
        vector<double> v{ 1,0 };
        for (int i = 0; i < k; ++i) {
            v[1] = -1 * x[i];
            Polynomial tmp(v);
            l = l * tmp * (1 / (x[k] - x[i]));
        }
        for (int i = k + 1; i < n; ++i) {
            v[1] = -x[i];
            Polynomial tmp(v);
            l = l * tmp * (1 / (x[k] - x[i]));
        }
        A[k] = integral_pl(l);
    }
    cout << setw(W) << "x_k" << setw(W) << "f_k" << setw(W) << "A_k\n";
    for (int i = 0; i < n; ++i) cout << setw(W) << x[i] << setw(W) << fx[i] << setw(W) << A[i] << "\n";
    I = 0.6205366034467622036;
    cout << "Точное значение интеграла: " << setprecision(12) << I << "\n";
    I_calc = 0;
    for (int i = 0; i < n; ++i) I_calc += A[i] * f(x[i]);
    cout << "Вычисленное значение интеграла: " << setprecision(12) << I_calc << "\n";
    cout << "Погрешность: " << setprecision(12) << fabs(I - I_calc) << "\n";

    delete[] A;
    delete[] x;
    delete[] fx;
    return 0;
}
