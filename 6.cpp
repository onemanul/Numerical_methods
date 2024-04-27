#include <iostream>
#include <algorithm>
using namespace std;

double funk(double x) { return exp(x) - x; }
double Funk(double a, double b) { return (exp(b) - b * b / 2) - (exp(a) - a * a / 2); }
double p0(double x) { return 3; }
double P0(double a, double b) { return 3 * b - 3 * a; }
double p1(double x) { return 3 * x + 5; }
double P1(double a, double b) { return (3 * b * b / 2 + 5 * b) - (3 * a * a / 2 + 5 * a); }
double p2(double x) { return 2 * x * x - 3 * x + 7; }
double P2(double a, double b) { return (2 * pow(b, 3) / 3 - 3 * b * b / 2 + 7 * b) - (2 * pow(a, 3) / 3 - 3 * a * a / 2 + 7 * a); }
double p3(double x) { return 5 * pow(x, 3) - 4 * x * x + 3 * x - 2; }
double P3(double a, double b) { return (5 * pow(b, 4) / 4 - 4 * pow(b, 3) / 3 + 3 * b * b / 2 - 2 * b) - (5 * pow(a, 4) / 4 - 4 * pow(a, 3) / 3 + 3 * a * a / 2 - 2 * a); }
double p4(double x) { return 6 * pow(x, 4) + 5 * pow(x, 3) - 4 * x * x + 3 * x - 2; }
double P4(double a, double b) { return (6 * pow(b, 5) / 5 + 5 * pow(b, 4) / 4 - 4 * pow(b, 3) / 3 + 3 * b * b / 2 - 2 * b) - (6 * pow(a, 5) / 5 + 5 * pow(a, 4) / 4 - 4 * pow(a, 3) / 3 + 3 * a * a / 2 - 2 * a); }
double p5(double x) { return 1.27*pow(x,5) + 2.04*x; }
double P5(double a, double b) { return (1.27*pow(b,6)/6 + 1.02*b*b) - (1.27*pow(a,6)/6 + 1.02*a*a); }

double f(double x, int c) {
    switch (c) {
    case 0: return p0(x);
    case 1: return p1(x);
    case 2: return p2(x);
    case 3: return p3(x);
    case 4: return p4(x);
    case 5: return p5(x);
    default: return funk(x);
    }
}
double F(double a, double b, int c) {
    switch (c) {
    case 0: return P0(a, b);
    case 1: return P1(a, b);
    case 2: return P2(a, b);
    case 3: return P3(a, b);
    case 4: return P4(a, b);
    case 5: return P5(a, b);
    default: return Funk(a, b);
    }
}

double Lefr_rect(double a, double b, int m, int c) {
    double S = 0, h = (b-a)/m;
    for (int i = 0; i < m; ++i) {
        S += f(a+i*h, c);
    }
    return h*S;
}
double Right_rect(double a, double b, int m, int c) {
    double S = 0, h = (b-a)/m;
    for (int i = 1; i <= m; ++i) {
        S += f(a+i*h, c);
    }
    return h*S;
}
double Middle_rect(double a, double b, int m, int c) {
    double S = 0, h = (b-a)/m;
    for (int i = 0; i < m; ++i) {
        S += f(a+(i+0.5)*h, c);
    }
    return h*S;
}
double Trapeze(double a, double b, int m, int c) { 
    double S = 0, h = (b-a)/m;
    for (int i = 0; i < m; ++i) {
        S += f(a+i*h, c) + f(a+(i+1)*h, c);
    }
    return S*h/2;
}
double Simpson(double a, double b, int m, int c) { 
    double S = 0, h = (b-a)/m;
    for (int i = 0; i < m; ++i) {
        S += f(a+i*h, c) + 4*f(a+(i+0.5)*h, c) + f(a+(i+1)*h, c);
    }
    return S*h/6;
}

int main() {
    setlocale(LC_ALL, "Russian");
    int m = 100;
    double a = -1, b = 1, x, J, A, w = 1; // w - весовая функция
    //cout << "Введите пределы интегрирования и шаг: "; 
    //cin >> a >> b >> m;
    for (int c = 0; c < 7; ++c) {
        switch (c) {
        case 0: cout << "\n\nМНОГОЧЛЕН НУЛЕВОЙ СТЕПЕНИ: p = 3"; break;
        case 1: cout << "\n\nМНОГОЧЛЕН ПЕРВОЙ СТЕПЕНИ: p = 3x + 5"; break;
        case 2: cout << "\n\nМНОГОЧЛЕН ВТОРОЙ СТЕПЕНИ: p = 2x^2 - 3x + 7"; break;
        case 3: cout << "\n\nМНОГОЧЛЕН ТРЕТЬЕЙ СТЕПЕНИ: p = 5x^3 - 4x^2 + 3x - 2"; break;
        case 4: cout << "\n\nМНОГОЧЛЕН ЧЕТВЁРТОЙ СТЕПЕНИ: p = 6x^4 + 5x^3 - 4x^2 + 3x - 2"; break;
        case 5: cout << "\n\nМНОГОЧЛЕН ПЯТОЙ СТЕПЕНИ: p = 1,27x^5 + 2,04x"; break;
        default: cout << "\n\nФУНКЦИЯ: f = e^x - x"; break;
        }
        J = w*F(a, b, c);
        cout << "\nВычисленное значение: " << J << "\n";
        x = w*Lefr_rect(a, b, m, c);  A = abs(J - x);
        cout << "\nКФ левого прямоугольника: " << x << "\n\t\t\t\tАбс. факт. погрешность: " << A << "\n\t\t\t\tОтн. факт. погрешность: " << A/abs(J);
        x =  w*Right_rect(a, b, m, c);  A = abs(J - x);
        cout << "\nКФ правого прямоугольника: " << x << "\n\t\t\t\tАбс. факт. погрешность: " << A << "\n\t\t\t\tОтн. факт. погрешность: " << A/abs(J);
        x = w*Middle_rect(a, b, m, c); A = abs(J - x);
        cout << "\nКФ среднего прямоугольника : " << x << "\n\t\t\t\tАбс. факт. погрешность: " << A << "\n\t\t\t\tОтн. факт. погрешность: " << A/abs(J);
        x = Trapeze(a, b, m, c); A = abs(J - x);
        cout << "\nКФ трапеции: " << x << "\n\t\t\t\tАбс. факт. погрешность: " << A << "\n\t\t\t\tОтн. факт. погрешность: " << A/abs(J);
        x = Simpson(a, b, m, c); A = abs(J - x);
        cout << "\nКФ Симпсона: " << x << "\n\t\t\t\tАбс. факт. погрешность: " << A << "\n\t\t\t\tОтн. факт. погрешность: " << A/abs(J);
    }

    return 0;
}
