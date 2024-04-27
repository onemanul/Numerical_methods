#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
using namespace std;

double F = 1.55, e = 1e-12;
double f(double x) { return exp(x) - x; }

bool comp(pair<double, double> a, pair<double, double> b) {
    return abs(a.first - F) < abs(b.first - F);
}
double Newton(vector<pair<double, double>>& xf, int n) {
    vector<double> a{ xf[0].second };
    double p = 0, k = 2;
    for (int i = 1; i <= n; ++i)
        a.push_back((xf[i].second - xf[i - 1].second) / (xf[i].first - xf[i - 1].first));
    while (k <= n) {
        for (int i = n; i >= k; --i)
            a[i] = (a[i] - a[i - 1]) / (xf[i].first - xf[i-k].first);
        ++k;
    }
    for (int i = 0; i <= n; ++i) {
        k = a[i];
        for (int j = 0; j < i; ++j) k *= (F - xf[j].first);
        p += k;
    }
    return p;
}
double first_sp(vector<pair<double, double>> fx, int n) {
    for (int i = 0; i < fx.size(); ++i) swap(fx[i].first, fx[i].second); // fx для f^(-1)
    sort(fx.begin(), fx.end(), comp);
    cout << "\nОтсортированная таблица для первого способа (столбцы поменялись)\n  k     xk\t f(xk)\n";
    for (int i = 0; i < n + 1; ++i)
        cout << setw(3) << i << setw(10) << fx[i].first << "\t" << setprecision(6) << fx[i].second << "\n";
    return Newton(fx, n);
}

vector<double> Newton_Ai(vector<pair<double, double>>& xf, int n) {
    vector<double> a{ xf[0].second };
    double p = 0, k = 2;
    for (int i = 1; i <= n; ++i)
        a.push_back((xf[i].second - xf[i - 1].second) / (xf[i].first - xf[i - 1].first));
    while (k <= n) {
        for (int i = n; i >= k; --i)
            a[i] = (a[i] - a[i - 1]) / (xf[i].first - xf[i - k].first);
        ++k;
    }
    return a;
}
double Pn_x(double x, vector<double>& a, vector<pair<double, double>>& xf) {
    double k, p = 0;
    for (int i = 0; i < a.size(); ++i) {
        k = a[i];
        for (int j = 0; j < i; ++j) k *= (x - xf[j].first);
        p += k;
    }
    return p-F; // Pn(x) - F = 0
}
void intervals(vector<double>& v, vector<double>& a, vector<pair<double, double>>& xf) {
    double c;
    while (v[1] > 0.0625) { // точность: 2/2^5
        for (int i = 1; i < v.size(); i += 2) {
            c = (v[i - 1] + v[i]) / 2;
            v.insert(v.begin() + i, c);
        }
    }
    for (int i = 0; i < v.size() - 1; ++i) {
        if (Pn_x(v[i],a,xf) * Pn_x(v[i+1],a,xf) > 0) {
            v.erase(v.begin() + i);
            --i;
        }
        else {
            v.insert(v.begin() + i + 1, v[i + 1]);
            ++i;
        }
    }
    v.pop_back();
}
vector<double> sek(vector<double>& v, vector<double>& a, vector<pair<double, double>>& xf) {
    vector<double> answ;
    double x0, x1, xk;
    for (int i = 0; i < v.size(); i += 2) {
        x0 = v[i];
        x1 = v[i+1];
        xk = x1 - Pn_x(x1,a,xf)*(x1-x0)/(Pn_x(x1, a, xf) - Pn_x(x0, a, xf));
        while (abs(xk - x1) > e) {
            x0 = x1;
            x1 = xk;
            xk = x1 - Pn_x(x1, a, xf) * (x1 - x0) / (Pn_x(x1, a, xf) - Pn_x(x0, a, xf));
        }
        answ.push_back(xk);
    }
    return answ;
}
vector<double> second_sp(vector<pair<double, double>>& xf, int n) {
    vector<double> a = Newton_Ai(xf, n), interv = {xf[0].first, xf[xf.size()-1].first}; // a и b
    intervals(interv, a, xf);
    return sek(interv, a, xf);
}

int main() {
    setlocale(LC_ALL, "Russian");
    int n = 7, m = 10;
    double a = 0, b = 1, xk, fir; // нет монотонности на (-2;2)
    vector<double> sec;
    vector<pair<double, double>> xf;
    /*cout << "Введите: F - точку для обратной интерполяции,\nn - степень интерполяционного многочлена (n <= m), е - погрешность:\n";
    cin >> F >> n >> e;
    while (n > m) {
        cout << "n должно быть не больше m (" << m << ") :";
        cin >> n;
    }*/
    cout << "  k     xk\t f(xk)\n";
    for (int i = 0; i < m + 1; ++i) {
        xk = a + i * (b - a) / m;
        xf.push_back({ xk,f(xk) });
        cout << setw(3) << i << setw(10) << xk << "\t" << setprecision(6) << xf[i].second << "\n";
    }
    fir = first_sp(xf, n);
    sec = second_sp(xf, n);

    cout << "\nЗначение первым способом: " << setprecision(10) << fir << " ( F = " << F << " )";
    cout << "\nРеальное значение функции в этой точке " << setprecision(10) << (xk=f(fir));
    cout << "\nАбсолютная погрешность: " << abs(xk - F) << "\n";
    for (int i = 0; i < sec.size(); ++i) {
        cout << "\nЗначение вторым способом: " << setprecision(10) << sec[i] << " ( F = " << F << " )";
        cout << "\nРеальное значение функции в этой точке " << setprecision(10) << (xk = f(sec[i]));
        cout << "\nАбсолютная погрешность: " << abs(xk - F) << "\n";
    }

return 0;
}
