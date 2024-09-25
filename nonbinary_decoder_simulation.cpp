//
// Created by Siyi Yang on 2/15/24.
//

#include <vector>
#include <time.h>
#include <math.h>
#include "tools.h"
#include "code_construction.h"
#include "quantum_ldpc_nonbinary.h"
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
using namespace std;

vector<vector<int>> readMatrixFromFile(const string& filename) {
    vector<vector<int>> matrix; // Vector to store the matrix
    ifstream file(filename); // Open the file

    if (!file.is_open()) { // Check if the file is opened successfully
        cerr << "Error: Unable to open file " << filename << endl;
        return matrix; // Return an empty matrix
    }

    string line;
    while (getline(file, line)) { // Read each line of the file
        vector<int> row; // Vector to store the elements of the current row
        istringstream iss(line); // Create a string stream from the current line

        int num;
        while (iss >> num) { // Read each integer from the current line
            row.push_back(num); // Add the integer to the current row
        }

        if (!row.empty()) { // Add the row to the matrix if it's not empty
            matrix.push_back(row);
        }
    }

    file.close(); // Close the file
    return matrix; // Return the matrix read from the file
}

//    string filename=argv[1]
//    double p = stod(argv[2]);
//    int num_trials = stoi(argv[3]);
//    int num_iter = stoi(argv[4]);
//    int job_index = stoi(argv[5]);
int main(int argc, char* argv[]) {

    /*
    int l=63;
    int z=7;
    vector<int> a{27,54,0};
    vector<int> b{0,1,6};
    vector<vector<vector<int> > > H=GB(l,z,a,b);
    vector<vector<int> > Hx=H[0];
    vector<vector<int> > Hz=H[1];
    int num_vn=Hx[0].size();

    vector<vector<int> > Hnb=Hx;
    for(int i=0;i<Hz.size();i++) {
        vector<int> t(0);
        for (int j = 0; j < Hz[0].size(); j++) {
            if(Hz[i][j])
                t.push_back(3);
            else
                t.push_back(0);
        }
        Hnb.push_back(t);
    }
     */
    vector<vector<int> > H=readMatrixFromFile(argv[1]);
    int num_vn=H[0].size();
    qCode qsc(H);

    // Read input parameters from command line
    double p = stod(argv[2]);
    int num_trials = stoi(argv[3]);
    int num_iter = stoi(argv[4]);

    // Get job index from command line
    int job_index = stoi(argv[5]);

    vector<vector<double> > lch;

    double p1=p/3.0;
    vector<double> p0({p1,p1,p1});
    lch= get_qlch(p0,num_vn);
    // Set up random number generator with job index as seed
    mt19937 rng(job_index);

    // Set up probability distribution=
    vector<double> probabilities = {1-p, p/3, p/3, p/3};
    discrete_distribution<int> dist(probabilities.begin(), probabilities.end());

    // Initialize variables for FER calculation
    int num_non_conv_errors = 0;
    int num_conv_errors=0;
    int num_success = 0;

    string evname = "ev_" + to_string(job_index) + ".txt";
    ofstream outev(evname);

    int res_type;
    // Loop over trials
    for (int i = 0; i < num_trials; i++) {
        // Generate random vector x
        vector<int> x(num_vn);
        generate(x.begin(), x.end(), [&]() { return dist(rng); });

        // Call LDPC decoder
        bool success = true; // Assume success
        vector<int> result= qsc.decode(x,lch,num_iter);
        res_type=result[num_vn];
        // Update variables for FER calculation
        if (!res_type) {
            num_success++;
            //cout<<i<<": decode success\n";
        }
        else {
            if(res_type==1) {
                num_conv_errors++;
                //cout<<i<<": logical error\n";
            }
            else {
                num_non_conv_errors++;
                //cout<<i<<": non-convergent error\n";
            }

            for(int j=0;j<=num_vn;j++)
                outev <<result[j]<<" ";
            outev<<endl;
        }
    }

    outev.close();
    // Write FER and num_trials to file
    string filename = "job_" + to_string(job_index) + ".txt";
    ofstream outfile(filename);
    outfile << num_conv_errors << " " << num_non_conv_errors << " " << num_trials << endl;
    outfile.close();

    return 0;
}
