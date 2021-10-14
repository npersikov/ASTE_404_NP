#include<iostream>

using namespace std;

int main()
{
    // int a[] = {0, 1, 2, 3};
    // int b[] = {2, 3, 4, 5};
    // int *c = new int[4];

    // for(int i = 0; i < 4; i++)
    // {
    //     c[i] = a[i] + b[i];
    //     cout << c[i] << endl;
    // }

    // float a = 0.0, b = 1.0, c = 2.0;
    // float &r = a;
    // float *p = &b;
    // p = &c;
    // *p = 3.0;
    // r = b + *p;

    // cout << c << endl;

    int data[] = {2,1,-2,5,0,3};
    int *p = &data[3];
    p[1] = 4;
    data[2] = 3;

    for(int i = 0; i < 3; i++)
    {
        cout << p[i] << endl;
    }
}