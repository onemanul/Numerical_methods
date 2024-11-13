#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
using namespace std;

double a = 0, b = acos(0.0), e = 1e-12; //pi/2 
double f(double x) { return pow(1 - sin(x) * sin(x) / 2, 0.5); }

double P1(double x) { return x; }
double P2(double x) { return ((3*x*x - 1) / 2); }
double P3(double x) { return ((5*x*x*x - 3*x) / 2); }
double P4(double x) {
    double x2 = x * x, x4 = x2 * x2;
    return ((35*x4 - 30*x2 + 3) / 8);
}
double P5(double x) {
    double x2 = x * x, x4 = x2 * x2;
    return ((63*x4*x - 70*x2*x + 15*x) / 8);
}
double P6(double x) {
    double x2 = x * x, x4 = x2 * x2, x6 = x2*x4;
    return ((231*x6 - 315*x4 + 105*x2 - 5) / 16);
}
double P7(double x) {
    double x2 = x * x, x4 = x2 * x2, x6 = x2 * x4;
    return ((429*x6*x - 693*x4*x + 315*x2*x - 35*x) / 16);
}
double P8(double x) {
    double x2 = x * x, x4 = x2 * x2, x6 = x2 * x4, x8 = x4 * x4;
    return ((6435*x8 - 12012*x6 + 6930*x4 - 1260*x2 + 35) / 128);
}
double Pn_x(double x, int N) {
    switch (N) {
    case 1: return P1(x);
    case 2: return P2(x);
    case 3: return P3(x);
    case 4: return P4(x);
    case 5: return P5(x);
    case 6: return P6(x);
    case 7: return P7(x);
    case 8: return P8(x);
    default: return 1; // PO
    }
}


double x5(double x) { return (6*pow(x,5) + 5*pow(x,4) - 4*pow(x, 3) + 3*x*x - 2*x + 1); } // 6x^5 + 5x^4 - 4x^3 + 3x^2 - 2x + 1, ответ 6
double x7(double x) { return (7*pow(x,7) - 5); } // 7x^7 - 5, ответ -10
double x9(double x) { return (pow(x,9) + pow(x,4)); }// x^9 + x^4, ответ 2/5 
double prov(double x, int N) {
    switch (N) {
    case 3: return x5(x);
    case 4: return x7(x);
    case 5: return x9(x);
    default: return 0;
    }
}

vector<double> sek(vector<double>& v, int N) {
    vector<double> answ;
    double x0, x1, xk;
    for (int i = 0; i < v.size(); i += 2) {
        x0 = v[i];
        x1 = v[i + 1];
        xk = x1 - Pn_x(x1,N) * (x1 - x0) / (Pn_x(x1,N) - Pn_x(x0,N));
        while (abs(xk - x1) > e) {
            x0 = x1;
            x1 = xk;
            xk = x1 - Pn_x(x1,N) * (x1 - x0) / (Pn_x(x1,N) - Pn_x(x0,N));
        }
        answ.push_back(xk);
    }
    return answ;
}
vector<double> intervals(int N) {
    double c;
    vector<double> v{ -1,1 };
    while (v[1] > 0.03125-1) { // точность: 2/2^6, интервал [-1,1]
        for (int i = 1; i < v.size(); i += 2) {
            c = (v[i-1] + v[i]) / 2;
            v.insert(v.begin() + i, c);
        }
    }
    for (int i = 0; i < v.size() - 1; ++i) {
        if (Pn_x(v[i],N) * Pn_x(v[i+1],N) > 0 || !(Pn_x(v[i], N))) { // второе условие - для случая "попадания" в корень, предотвращает двоение интервала
            v.erase(v.begin() + i);
            --i;
        }
        else {
            v.insert(v.begin() + i + 1, v[i+1]);
            ++i;
        }
    }
    v.pop_back();
    return v;
}
vector<double> C_k(vector<double>& t, int N) {
    vector<double> c;
    double x;
    for (int i = 0; i < t.size(); ++i) {
        x = (2 - 2*t[i]*t[i]) / (N*N*Pn_x(t[i],N-1)*Pn_x(t[i], N-1));
        c.push_back(x);
    }
    return c;
}
vector<double> A_k(vector<double>& c, int N) {
    vector<double> A;
    double x;
    for (int i = 0; i < c.size(); ++i) {
        x = (b-a)*c[i] / 2;
        A.push_back(x);
    }
    return A;
}
vector<double> X_k(vector<double>& t) {
    vector<double> X;
    double x;
    for (int i = 0; i < t.size(); ++i) {
        x = (b-a)*t[i]/2 + (b-a)/2;
        X.push_back(x);
    }
    return X;
}

int main() {
    setlocale(LC_ALL, "Russian");
    double sum;
    vector<double> interv, t, C, A, X;

    for (int i = 1; i < 9; ++i) {
        interv = intervals(i);
        t = sek(interv, i);
        C = C_k(t, i);
        cout << "N = " << i << ":\n";
        for (int k = 0; k < C.size(); ++k) {
            cout << "\t Узел: " << fixed << setprecision(6) << setw(10) << t[k];
            cout << "\tКоэффициент: " << fixed << setprecision(12) << C[k] << "\n";
        }
    }
    cout << "\nI = 1.35064388104767550\n\n";
    for (int i = 6; i < 9; ++i) {
        interv = intervals(i);
        t = sek(interv, i);
        C = C_k(t, i);
        A = A_k(C, i);
        X = X_k(t);
        sum = 0;
        for (int k = 0; k < A.size(); ++k) sum += A[k] * f(X[k]);
        cout << "N = " << i << ":\n";
        for (int k = 0; k < C.size(); ++k) {
            cout << "\t Узел: " << fixed << setprecision(6) << setw(10) << t[k];
            cout << "\tКоэффициент для [a,b]: " << fixed << setprecision(12) << A[k] << "\n";
        }
        cout << "Значение интеграла: " << fixed << setprecision(12) << sum << "\n\n";
    }
    cout << "\nБлок проверок.\nОтветы для интервала [-1;1]: степень 5: 6;   степень 7: -10;   степень 9: 2/5\n";
    for (int i = 3; i < 6; ++i) {
        interv = intervals(i);
        t = sek(interv, i);
        C = C_k(t, i);
        sum = 0;
        for (int k = 0; k < C.size(); ++k) sum += C[k] * prov(t[k],i);
        cout << "N = " << i << ":\n\tМногочлен степени " << 2*i-1 << ": " << sum << "\n";
    }

    return 0;
}
