/* 
Author: Trevor Fung
Date: 9/27/19
Purpose: port MATLAB's cmdscale.m to cpp/Arduino. Validating using online compiler
*/

#include <iostream>

using namespace std;

#include <vector>

float* matlab_cmdscale(float** distance_matrix, int n, int m){
    //create identity matrix
    std::vector<vector<float>> eye_n(n);
    for(int i = 0; i < n; ++i){
        eye_n[i] = std::vector<float>(n);
        for(int j = 0; j < n; ++j){
            if(i == j){
                eye_n[i][j] = 1.0;
            }else{
                eye_n[i][j] = 0.0;
            }
        }
    }
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            cout << eye_n[i][j] << " ";
        }
        cout << endl;
    }
}

int main()
{
   // cout << "Hello World" << endl; 
   float** distance_matrix;
   int n = 5;
   int m = 3;
   matlab_cmdscale(distance_matrix, n, m);
   return 0;
}