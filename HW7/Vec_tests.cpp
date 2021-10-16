#include "gtest/gtest.h"
#include "Vec.h"

TEST(VecTest, VecDot)
{
    double3 a{0,1,2};
    double3 b{0,0,1};
    EXPECT_EQ(dot(a,b), 2);
}

/**
* This class is used to set up for unit testing.
*/
class VecTestClass : public ::testing::Test 
{
    protected:

    VecTestClass() : a{0,1,2}, b{0,0,1} {  }
    ~VecTestClass() override { /* nothing to do */  }

    void SetUp() override 
    {
        // test addition
        double3 c = a+b;
        add_matches = true;  // set default

        for (int i=0;i<3;i++) 
            if (c[i] != a[i] + b[i]) 
                add_matches=false;  // invalidate if component doesn't match
                
        // test subtraction
        double3 d = a-b;
            sub_matches = true;  // set default

        for (int i=0;i<3;i++) 
            if (d[i] != a[i] - b[i]) 
                sub_matches=false;  
    }

    void TearDown() override 
    {
    /* optional code to call after each test prior to destructor */
    }

    double3 a, b; // class data
    bool sub_matches = false;
    bool add_matches = false;
};

/**
* This function tests the addition function to see if it works correctly.
* I am not sure what these arguments are though. One is just a data type and
* one is just a name???
*/
TEST_F(VecTestClass, VecAdd) 
{ // makes a VecTestClass, checks value of add_matches
    EXPECT_EQ(add_matches, true);
}

/**
 *  This tests subtraction
 */
TEST_F(VecTestClass, VecSub)
{
    EXPECT_EQ(sub_matches, true);
}

/**
* @param argc not sure what this one is for tbh
* @param argv not sure wht this is either
*/
int main(int argc, char**argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}