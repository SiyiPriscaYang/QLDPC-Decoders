//
// Created by Siyi Yang on 3/30/23.
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

int main(int argc, char *argv[]) {
    if (argc < 4) {
        cerr << "Usage: " << argv[0] << " filename lb rb" << endl;
        return 1;
    }
    string filename = argv[1];
    int lb = stoi(argv[2]);
    int rb = stoi(argv[3]);

    ofstream outfile(filename + ".txt");
    int count = 0; // count of nonzero entries
    int total_conv_errors = 0;
    int total_non_conv_errors = 0;
    int total_trials = 0;
    int num_cols = 0; // number of columns

    bool found_first_file = false;
    int conv_errors, non_conv_errors, trials;
    int index;

    // test with "/Users/prisca/CLionProjects/QLDPC/" at local PC

    for (int i = lb; i <= rb; i++) {
        string filename = "job_" + to_string(i) + ".txt";
        ifstream infile(filename);
        if (infile.is_open()) {
            infile >> conv_errors >> non_conv_errors >> trials;
            if(!found_first_file){
                found_first_file=true;
                break;
            }

            infile.close();
        }
        if(found_first_file)
            break;
    }
    total_trials=trials*(rb-lb+1);

    found_first_file = false;
    for (int i = lb; i <= rb; i++) {
        string infile_name = "ev_" + to_string(i) + ".txt";
        ifstream infile2(infile_name);
        if (infile2.is_open()) {
            string line;
            while (getline(infile2, line)) {
                istringstream iss(line);
                int value;
                int col_count = 0;
                while (iss >> value)
                    col_count++;
                if (!found_first_file) {
                    num_cols = col_count;
                    found_first_file = true;
                    break;
                }
            }
            infile2.close();
        }
        if(found_first_file)
            break;
    }

    num_cols--;

    count=0;
    total_non_conv_errors=0;
    total_conv_errors=0;

    for (int i = lb; i <= rb; i++) {
        string infile_name = "ev_" + to_string(i) + ".txt";
        ifstream infile2(infile_name);
        if (infile2.is_open()) {
            string line;
            while (getline(infile2, line)) {
                istringstream iss(line);
                int value;
                index=0;
                while (iss >> value) {
                    if (value != 0)
                        count++;
                    outfile << value << " ";
                    if(index==num_cols){
                        if(value==1)
                            total_conv_errors++;
                        else
                            total_non_conv_errors++;
                    }
                    index++;
                }
                outfile << endl;
            }
            infile2.close();
        }
    }

    outfile.close();
    int total_errors=total_conv_errors+total_non_conv_errors;
    int bits_errors=count-total_errors;

    double average_fer = ((double) total_errors) / ((double) total_trials);
    cout << "Average FER: " << average_fer << endl;
    cout << "total errors: " << total_errors << endl;
    cout << "total convergent errors: " << total_conv_errors << endl;
    cout << "total non-convergent errors: " << total_non_conv_errors << endl;
    cout << "total trials: " << total_trials << endl;

    double average_ber = (((double) bits_errors) / ((double) total_trials))/((double) num_cols);
    cout << "Average BER: " << average_ber << endl;
    cout << "Total number of erroneous bits: " << bits_errors << endl;
    cout << "Number of columns: " << num_cols << endl;

    return 0;
}
