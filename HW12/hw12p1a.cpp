#include <iostream>
#include <thread>

using namespace std;

void doWork(int *counter)
{
    *counter += 1;
}

int main()
{
    for(int s = 0; s < 100000; s++)
    {
        int counter = 0;
        thread t1(doWork, &counter);
        thread t2(doWork, &counter);

        t1.join();
        t2.join();
        if(counter != 2) cout << "*" << endl;
    }
    return 0;
}