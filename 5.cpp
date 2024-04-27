#include <iostream>
#include <algorithm>
using namespace std;

double funk(double x) { return exp(x) - x; }
double Funk(double a, double b) { return (exp(b)-b*b/2) - (exp(a) - a*a/2); }
double p0(double x) { return 3; }
double P0(double a, double b) { return 3*b - 3*a; }
double p1(double x) { return 3*x+5; }
double P1(double a, double b) { return (3*b*b/2 + 5*b) - (3*a*a/2 + 5*a); }
double p2(double x) { return 2*x*x - 3*x + 7; }
double P2(double a, double b) { return (2*pow(b,3)/3 - 3*b*b/2 + 7*b) - (2*pow(a,3)/3 - 3*a*a/2 + 7*a); }
double p3(double x) { return 5*pow(x,3) - 4*x*x + 3*x - 2; }
double P3(double a, double b) { return (5*pow(b,4)/4 - 4*pow(b,3)/3 + 3*b*b/2 - 2*b) - (5*pow(a,4)/4 - 4*pow(a,3)/3 + 3*a*a/2 - 2*a); }
double p4(double x) { return 6*pow(x,4) + 5*pow(x,3) - 4*x*x + 3*x - 2; }
double P4(double a, double b) { return (6*pow(b,5)/5 + 5*pow(b,4)/4 - 4*pow(b,3)/3 + 3*b*b/2 - 2*b) - (6*pow(a,5)/5 + 5*pow(a,4)/4 - 4*pow(a,3)/3 + 3*a*a/2 - 2*a); }

double f(double x, int c) { 
    switch (c) {
    case 0: return p0(x);
    case 1: return p1(x);
    case 2: return p2(x);
    case 3: return p3(x);
    case 4: return p4(x);
    default: return funk(x);
    }
}
double F(double a, double b, int c) { 
    switch (c) {
    case 0: return P0(a,b);
    case 1: return P1(a,b);
    case 2: return P2(a,b);
    case 3: return P3(a,b);
    case 4: return P4(a,b);
    default: return Funk(a,b);
    }
}

double Lefr_rect(double a, double b, int c) { return (b-a) * f(a, c); }
double Right_rect(double a, double b, int c) { return (b-a) * f(b, c); }
double Middle_rect(double a, double b, int c) { return (b-a) * f((a+b)/2, c); }
double Trapeze(double a, double b, int c) { return (b-a)/2 * (f(a, c)+f(b, c)); }
double Simpson(double a, double b, int c) { return (b-a)/6 * (f(a, c) + 4*f((a+b)/2, c) + f(b, c)); }
double Three_eight(double a, double b, int c) { return (b-a) * (f(a, c)/8 + 3*f((2*a+b)/3,c)/8 + 3*f((a+2*b)/3,c)/8 + f(b, c)/8); }

int main() {
    setlocale(LC_ALL, "Russian");
    double a = 0, b = 1, x, X;
    //cout << "Введите пределы интегрирования: "; 
    //cin >> a >> b;
    for (int c = 0; c < 6; ++c) {
        switch (c) {
        case 0: cout << "\n\nМНОГОЧЛЕН НУЛЕВОЙ СТЕПЕНИ: p = 3"; break;
        case 1: cout << "\n\nМНОГОЧЛЕН ПЕРВОЙ СТЕПЕНИ: p = 3x + 5"; break;
        case 2: cout << "\n\nМНОГОЧЛЕН ВТОРОЙ СТЕПЕНИ: p = 2x^2 - 3x + 7"; break;
        case 3: cout << "\n\nМНОГОЧЛЕН ТРЕТЬЕЙ СТЕПЕНИ: p = 5x^3 - 4x^2 + 3x - 2"; break;
        case 4: cout << "\n\nМНОГОЧЛЕН ЧЕТВЁРТОЙ СТЕПЕНИ: p = 6x^4 + 5x^3 - 4x^2 + 3x - 2"; break;
        default: cout << "\n\nФУНКЦИЯ: f = e^x - x"; break;
        }
        cout << "\nВычисленное значение: " << (X = F(a,b,c)) << "\n";
        cout << "\nКФ левого прямоугольника: " << (x = Lefr_rect(a, b, c)) << "\n\t\t\t\tАбс. факт. погрешность: " << abs(X - x);
        cout << "\nКФ правого прямоугольника: " << (x = Right_rect(a, b, c)) << "\n\t\t\t\tАбс. факт. погрешность: " << abs(X - x);
        cout << "\nКФ среднего прямоугольника : " << (x = Middle_rect(a, b, c)) << "\n\t\t\t\tАбс. факт. погрешность: " << abs(X - x);
        cout << "\nКФ трапеции: " << (x = Trapeze(a, b, c)) << "\n\t\t\t\tАбс. факт. погрешность: " << abs(X - x);
        cout << "\nКФ Симпсона: " << (x = Simpson(a, b, c)) << "\n\t\t\t\tАбс. факт. погрешность: " << abs(X - x);
        cout << "\nКФ 3/8: " << (x = Three_eight(a, b, c)) << "\n\t\t\t\tАбс. факт. погрешность: " << abs(X - x);
    }

return 0;
}
