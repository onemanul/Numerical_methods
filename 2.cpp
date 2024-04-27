#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
using namespace std;

double x = 0.65;
double f(double x) { return exp(x) - x; }
bool comp(pair<double, double> a, pair<double, double> b) {
    return abs(a.first - x) < abs(b.first - x);
}
double Newton(vector<pair<double, double>> xf, int n) {
    vector<double> a{ xf[0].second };
    double p = 0, k=2;
    for (int i = 1; i <= n; ++i)
        a.push_back( (xf[i].second - xf[i-1].second) / (xf[i].first - xf[i-1].first) );
    while (k <= n) {
        for (int i = n; i >= k; --i)
            a[i]  = (a[i]-a[i-1])/ (xf[i].first - xf[i-k].first);
        ++k;
    }
    for (int i = 0; i <= n; ++i) {
        k = a[i];
        for (int j = 0; j < i; ++j) k *= (x - xf[j].first);
        p += k;
    }
    return p;
}
double Lagrange(vector<pair<double, double>> xf, int n) {
    double p = 0, a, b;
    for (int i=0; i <= n; ++i) {
        a = 1; b = 1;
        for (int j = 0; j <= n; ++j) {
            if (i != j) {
                a *= x - xf[j].first;
                b *= xf[i].first - xf[j].first;
            }
        }
        p += (a/b)*xf[i].second;
    }
    return p;
}

int main(){
    setlocale(LC_ALL, "Russian");
    int n=7, m=15;
    double a = 0, b = 1, xk, newt, lagr;
    vector<pair<double, double>> xf; 
    /*cout << "Введите: m, концы отрезка [a,b], x - точку интерполирования,\n";
    cout << "n - степень интерполяционного многочлена (n <= m):\n";
    cin >> m >> a >> b >> x >> n;
    while (n > m) {
        cout << "n должно быть не больше m (" << m << ") :";
        cin >> n;
    }*/
    cout << "  k     xk\t f(xk)\n";
    for (int i = 0; i < m + 1; ++i) {
        xk = a + i*(b-a)/m;
        xf.push_back({ xk,f(xk) });
        cout << setw(3) << i << setw(10) << xk << "\t" << setprecision(6) << xf[i].second << "\n";
    }
    sort(xf.begin(), xf.end(), comp);
    cout << "\nОтсортированная таблица\n  k     xk\t f(xk)\n";
    for (int i = 0; i < n + 1; ++i) 
        cout << setw(3) << i << setw(10) << xf[i].first << "\t" << setprecision(6) << xf[i].second << "\n";

    xk = f(x);
    newt = Newton(xf, n);
    lagr = Lagrange(xf, n);
    cout << "\nЗначение в точке: " << setprecision(10) << xk;
    cout << "\nЗначение по Ньютону: " << setprecision(10) << newt;
    cout << "\nАбсолютная погрешность: " << abs(xk-newt);
    cout << "\nЗначение по Лагранжу: " << setprecision(10) << lagr;
    cout << "\nАбсолютная погрешность: " << abs(xk-lagr) << "\n";

return 0;
}
