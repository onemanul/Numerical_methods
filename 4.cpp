#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace std;

double h = 0.1, a = 0;
int m = 10; // m >= 2
double f(double x) { return sin(x); }
double pr_1(double x) { return cos(x); }
double pr_2(double x) { return -sin(x); }
double chd_1(double x) { 
    if (x == a) return (-3*f(x) + 4*f(x+h) - f(x+2*h)) / (2*h);
    if (x == a+h*m) return (3*f(x) - 4*f(x-h) + f(x-2*h)) / (2*h);
    return (f(x+h) - f(x-h)) / (2*h); 
}
double chd_2(double x) { 
    if (x == a || x == a + h * m) return INFINITY;
    return (f(x+h) - 2*f(x) + f(x-h)) / (h*h); 
}

int main() {
    setlocale(LC_ALL, "Russian");
    double xk, pr1, pr2, chd1, chd2, abs1, abs2;
    /*cout << "Введите: m - количество значений функции,\na - начальную точку, h - шаг:\n";
    cin >> m >> a >> h;*/
    cout << "\n  xi \t    f(xi) \t  f'(xi)ЧД \t |f'(xi)Т-f'(xi)ЧД| \t отн.погр-ть f'    f''(xi)ЧД \t |f''(xi)Т-f''(xi)ЧД| \t отн.погр-ть f''\n\n";
    for (int i = 0; i < m + 1; ++i) {
        xk = a + i * h;
        chd1 = chd_1(xk);
        chd2 = chd_2(xk);
        pr1 = pr_1(xk);
        pr2 = pr_2(xk);
        abs1 = abs(chd1 - pr1);
        abs2 = abs(chd2 - pr2);
        cout << setw(4) << xk << "\t" << setw(10) << f(xk) << "\t" << setw(10) << chd1 << "\t     ";
        cout << setw(8) << abs1 << "\t\t " << setw(10) << abs1/abs(chd1) << "\t   " << setw(8) << chd2;
        cout << "\t       " << setw(8) << abs2 << "\t       " << setw(10) << abs2/abs(chd2) << "\n";
    }

return 0;
}
