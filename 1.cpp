#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

double func(double x) {
    double x2=x*x, x4=x2*x2;
    return (32*x4*x2 - 48*x4 + 18*x2 - 1 - 9/(x+10));
}
double proizv(double x) {
    double x2 = x*x;
    return (32*6*x2*x2*x - 48*4*x2*x + 18*2*x + 9/(x2+20*x+100));
}
double phi(double x, bool var) {
    double x2 = x*x, x4 = x2*x2;
    if (var) return ( (32*x4*x2 - 48*x4 - 1 - 9/(x + 10)) / (-18*x) ); // выразили х из 18х^2. Подходит для 1-го и 3-го корня
    else return ((32 * x4 * x2 + 18 * x2 - 1 - 9 / (x + 10)) / (48 * x2 * x)); // выразили х из -48x^4. Подходит для 2-го корня
}

void intervals(vector<double>& v) {
    double c;
    while (v[1] > 0.0625) { // точность: 2/2^5
        for (int i = 1; i < v.size(); i += 2) {
            c = (v[i - 1] + v[i]) / 2;
            v.insert(v.begin() + i, c);
        }
    }
    for (int i = 0; i < v.size() - 1; ++i) {
        if (func(v[i]) * func(v[i + 1]) > 0) {
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

vector<double> Newton(vector<double>& v) {
    vector<double> answ;
    double x0, xk;
    int k = 1;
    for (int i = 0; i<v.size(); i += 2) {
        x0 = (v[i] + v[i+1]) / 2;
        xk = x0 - func(x0) / proizv(x0);
        while (abs(xk-x0)>0.000001) {
            x0 = xk;
            xk = x0 - func(x0)/proizv(x0);
            ++k;
        }
        answ.push_back(xk);
        cout << "Итераций в Ньютоне для " << i / 2 + 1 << "-го корня: " << k << "\n";
        k = 1;
    }
    return answ;
}

vector<double> Newton_modif(vector<double>& v) {
    vector<double> answ;
    double x0, xk, xk_1;
    int k = 1;
    for (int i = 0; i < v.size(); i += 2) {
        x0 = (v[i] + v[i + 1]) / 2;
        xk_1 = x0;
        xk = xk_1 - func(xk_1) / proizv(x0);
        while (abs(xk - xk_1) > 0.000001) {
            xk_1 = xk;
            xk = xk_1 - func(xk_1) / proizv(x0);
            ++k;
        }
        answ.push_back(xk);
        cout << "Итераций в модифицированном Ньютоне для " << i / 2 + 1 << "-го корня: " << k << "\n";
        k = 1;
    }
    return answ;
}

vector<double> iteracii(vector<double>& v) {
    vector<double> answ;
    double x0, xk;
    bool var;
    int k = 1;
    for (int i = 0; i < v.size(); i += 2) {
        i == 2 ? var = 0 : var = 1; // для нужной функции phi
        x0 = (v[i] + v[i + 1]) / 2;
        xk = phi(x0, var);
        while (0.66*abs(xk-x0)/(1 - 0.66) > 0.001) { // 0.66 - "ручная" оценка, для обеих функций phi
            x0 = xk;
            xk = phi(x0, var);
            ++k;
        }
        answ.push_back(xk);
        cout << "Итераций в методе итераций для " << i/2 + 1 << "-го корня: " << k << "\n";
        k = 1;
    }
    return answ;
}

int main() {
    setlocale(LC_ALL, "Russian");
    double a = 0, b = 2;
    vector<double> interv{a,b}, newt, newt_modif, iter;

    intervals(interv);
    cout << "Первое задание: ";
    for (int i = 0; i < interv.size(); i += 2) cout << "[" << interv[i] << ", " << interv[i+1] << "] ";

    cout << "\n\nВторое задание:\n";
    newt = Newton(interv);
    for (double x : newt) cout << x << " ";

    cout << "\n\nВторое с половиной задание:\n";
    newt_modif = Newton_modif(interv);
    for (double x : newt_modif) cout << x << " ";

    cout << "\n\nТретье задание:\n";
    iter = iteracii(interv);
    for (double x : iter) cout << x << " ";

return 0;
}
