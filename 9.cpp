#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
using namespace std;

double f(double y) { return y*y-2*y; }
double y(double x) { return 2 / (exp(2*x) + 1); }

unsigned int fact(unsigned int n) {
    if (n == 0) return 1;
    return n*fact(n-1);
}
vector<double> Taylor(vector<double>& X, int N, double x0) { // -2
    vector<double> koef_y{ 1, -1, 0, 2, 0, -16 }, Tay;
    double x;
    for (int k = 0; k < X.size(); ++k) {
        x = 0;
        for (int i = 0; i < koef_y.size(); ++i) 
            x += koef_y[i] * pow((X[k]-x0), i) / fact(i);
        Tay.push_back(x);
    }
    return Tay;
}

vector<double> Euler(int N, double h, double y0) { // 0
    vector<double> Eul{ y0 };
    for (int i = 1; i <= N; i++)
        Eul.push_back( Eul[i-1] + h * f(Eul[i-1]) );
    return Eul;
}
vector<double> EulerOne(int N, double h, double y0) { // 0
    vector<double> Eul{ y0 };
    for (int i = 1; i <= N; i++)
        Eul.push_back( Eul[i-1] + h * f(Eul[i-1] + h/2 * f(Eul[i-1])) );
    return Eul;
}
vector<double> EulerTwo(int N, double h, double y0) { // 0
    vector<double> Eul{ y0 };
    for (int i = 1; i <= N; i++)
        Eul.push_back( Eul[i-1] + h/2 * ( f(Eul[i-1]) + f(Eul[i-1] + h*f(Eul[i-1])) ) );
    return Eul;
}

vector<double> Rengue(int N, double h, double y0) { // 0
    vector<double> Ren{ y0 };
    double k1, k2, k3, k4;
    for (int i = 1; i <= N; i++) {
        k1 = h * f( Ren[i-1] );
        k2 = h * f( Ren[i-1] + k1 / 2 );
        k3 = h * f( Ren[i-1] + k2 / 2 );
        k4 = h * f( Ren[i-1] + k3 );
        Ren.push_back( Ren[i-1] + (k1 + 2*k2 + 2*k3 + k4)/6 );
    }
    return Ren;
}
vector<double> Adams(vector<double>& Y, int N, double h) { // -2
    vector<double> Adm;
    double** matrix;
    matrix = new double* [20];
    for (int i = 0; i < 20; i++) matrix[i] = new double[20];

    for (int i = 0; i < 5; i++) matrix[0][i] = Y[i];
    for (int i = 0; i < 5; i++) matrix[1][i] = h * f(matrix[0][i]);
    for (int i = 2; i < 6; i++)
        for (int j = 0; j < 6-i; j++)
            matrix[i][j] = matrix[i-1][j+1] - matrix[i-1][j];

    for (int i = 1; i < 10; i++) {
        matrix[0][4+i] = matrix[0][4 + i - 1] + matrix[1][4 + i - 1] + matrix[2][4 + i - 2] / 2 + 5 * matrix[3][4 + i - 3] / 12 + 3 * matrix[4][4 + i - 4] / 8 + 251 * matrix[5][4 + i - 5] / 720;
        matrix[1][4 + i] = h * f(matrix[0][4+i]);
        matrix[2][4 + i - 1] = matrix[1][4+i] - matrix[1][4 + i - 1];
        matrix[3][4 + i - 2] = matrix[2][4 + i - 1] - matrix[2][4 + i - 2];
        matrix[4][4 + i - 3] = matrix[3][4 + i - 2] - matrix[3][4 + i - 3];
        matrix[5][4 + i - 4] = matrix[4][4 + i - 3] - matrix[4][4 + i - 4];
    }
    for (int i = 3; i <= N; i++) Adm.push_back(matrix[0][i+2]);
    return Adm;
}

void show(vector<double>& X, vector<double>& answ, vector<double>& Y, int beg) {
    for (int i = 0; i < X.size(); ++i) {
        cout << "\t X_" << i+beg << " = " << fixed << setprecision(2) << setw(5) << X[i];
        cout << "\t\tY_" << i+beg << " = " << fixed << setprecision(6) << setw(10) << answ[i];
        cout << "\t\t|yn - y(xn)| = " << fixed << setprecision(8) << setw(12) << fabs(answ[i] - Y[i]) << "\n";
    }
}

int main() {
    setlocale(LC_ALL, "Russian");
    int N = 10;
    double x0 = 0, y0 = 1, h = 0.1;
    vector<double> X{x0}, Y{y0}, Xm2, Ym2, answ, last; // минус 2
    //cout << "Введите шаг и N"; cin >> h >> N;

    for (int i = 1; i <= N; ++i) {
        X.push_back(x0 + i*h);
        Y.push_back(y(X[i]));
    }
    for (int i = -2; i <= N; ++i) {
        Xm2.push_back(x0 + i*h);
        Ym2.push_back(y(Xm2[i+2]));
    }
    last.push_back(Y.back());

    cout << "Задача: y' = y^2 - 2y\nТочное решение для y(0) = 1: y = 2/(e^(2x) + 1)\n";
    cout << "\nТочные значение:\n";
    for (int i = -2; i <= N; ++i) {
        cout << "\t X_" << i << " = " << fixed << setprecision(2) << setw(5) << Xm2[i+2];
        cout << "\t\tY_" << i << " = " << fixed << setprecision(6) << setw(10) << Ym2[i+2] << "\n";
    }

    answ = Taylor(Xm2, N, x0);
    last.push_back(answ.back());
    cout << "\nМетод Тейлора:\n";
    show(Xm2, answ, Ym2, -2);

    answ = Euler(N, h, y0);
    last.push_back(answ.back());
    cout << "\nМетод Эйлера:\n";
    show(X, answ, Y, 0);

    answ = EulerOne(N, h, y0);
    last.push_back(answ.back());
    cout << "\nМетод Эйлера I:\n";
    show(X, answ, Y, 0);

    answ = EulerTwo(N, h, y0);
    last.push_back(answ.back());
    cout << "\nМетод Эйлера II:\n";
    show(X, answ, Y, 0);

    answ = Rengue(N, h, y0);
    last.push_back(answ.back());
    cout << "\nМетод Рунге-Кутта:\n";
    show(X, answ, Y, 0);

    answ = Adams(Ym2, N, h);
    last.push_back(answ.back());
    cout << "\nМетод Адамса:\n";
    for (int i = 3; i <= N; ++i) {
        cout << "\t X_" << i << " = " << fixed << setprecision(2) << setw(5) << X[i];
        cout << "\t\tY_" << i << " = " << fixed << setprecision(6) << setw(10) << answ[i-3];
        cout << "\t\t|yn - y(xn)| = " << fixed << setprecision(8) << setw(12) << fabs(answ[i-3] - Y[i]) << "\n";
    }

    cout << "\nСводная таблица по последним значениям:\n\tМетод \t\t   Y \t\t Абс. погрешность\n";
    cout << "\tТейлор\t\t" << fixed << setprecision(8) << setw(10) << last[1] << "\t" << setprecision(8) << setw(12) << fabs(last[1] - last[0]) << "\n";
    cout << "\tЭйлер\t\t" << fixed << setprecision(8) << setw(10) << last[2] << "\t" << setprecision(8) << setw(12) << fabs(last[2] - last[0]) << "\n";
    cout << "\tЭйлер I\t\t" << fixed << setprecision(8) << setw(10) << last[3] << "\t" << setprecision(8) << setw(12) << fabs(last[3] - last[0]) << "\n";
    cout << "\tЭйлер II\t" << fixed << setprecision(8) << setw(10) << last[4] << "\t" << setprecision(8) << setw(12) << fabs(last[4] - last[0]) << "\n";
    cout << "\tРунге-Кутт\t" << fixed << setprecision(8) << setw(10) << last[5] << "\t" << setprecision(8) << setw(12) << fabs(last[5] - last[0]) << "\n";
    cout << "\tАдамс\t\t" << fixed << setprecision(8) << setw(10) << last[6] << "\t" << setprecision(8) << setw(12) << fabs(last[6] - last[0]) << "\n";

    return 0;
}
