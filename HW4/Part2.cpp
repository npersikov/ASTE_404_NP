#include <iostream>

using namespace std;

// void init(double *a, double *b) {
//   *a = 2; *b = 1;
// }

// int main()
// {
//     // Part 2a
//     // float x = 1.2;
//     // float y = 4.0;
//     // float &r = x;
//     // r = 2*x + y;

//     // cout << x << endl;

//     //Part 2b
//     double a, b;
//     init(&a, &b);
//     double c = a + b;

//     cout << c << endl;

//     return 0;
// }

void init(double &a, double &b) 
{
  a = 2; 
  b = 1;
}

int main()
{
    double a, b;
    init(a, b);
    double c = a + b;

    cout << c << endl;

    return 0;
}
