// Team Members: Sarah Villegas, Ali Fenton, Chanel Aquino
// Date Created: 30 October 2015
// Description: Speech Signal Analysis

//B. Standard deviation is the most different
//C. The zero crossings are the most different
//D. Use the frequency as another measurement
//E. They are not from the same person since the data is too different in percentages.
#include <iostream> // cin, cout
#include <fstream>  // file IO
#include <cmath>    // sqrt(), abs()
#include <iomanip>  // setw()
#include <cstdlib>  // exit()
#include <vector>   

using namespace std;

double mean(vector<double> vct);	// chanel
double averagePower(vector<double> vct);	//chanel
double stdDev(double vct);  //Ali
double avgMagnitude(vector<double> vct); //Ali
double var(vector<double> vct , double& mean);  
int zeroCross(vector<double> vct);


//*******************************
int main()
{
    vector<double> sound1;
    vector<double> sound2;
    
    int 
        index = 0,
        size;
        
   double meanDifference=0, zeroCrossDifference=0, magnitudeDifference=0, powerDifference=0, varDifference=0, stdDifference=0;
   double power1 = 0, power2 = 0, zero1 = 0, zero2 = 0, mag1 =0, mag2 = 0, std1 = 0, std2 = 0;
    double 
        value,
        value1,
        mean1,
        mean2,
        standD1,
        standD2;
    
    // create input/out files
    ifstream fin1, fin2;
    ofstream fout;
    
    // open files
    fin1.open("two_a.txt");
    fin2.open("two_b.txt");
    fout.open("comparison.txt");
    
    // exit program if fail
    if (fin1.fail() || fin2.fail() || fout.fail())
    {
        cout << "ERROR";
        exit(1);
    }    


    while(fin1 >> value)
    {
        sound1.push_back(value);
        
    }
    
    while(fin2 >> value1){
        sound2.push_back(value1);
    }
    
    mean1 = mean(sound1);
    mean2 = mean(sound2);
    standD1 = var(sound1, mean1);
    standD2 = var(sound2, mean2);
    power1 = averagePower(sound1);
    power2 = averagePower(sound2);
    zero1 = zeroCross(sound1);
    zero2 = zeroCross(sound2);
    mag1 = avgMagnitude(sound1);
    mag2 = avgMagnitude(sound2);
    std1 = stdDev(standD1);
    std2 = stdDev(standD2);
    
    powerDifference = (abs(power1 - power2)/ ((power1 + power2)/2)) * 100.0;
    meanDifference = (abs(mean1 - mean2)/ ((mean1 + mean2)/2)) * 100.0;
    zeroCrossDifference = (abs(zero1 - zero2)/ ((zero1 + zero2)/2)) * 100.0;
    magnitudeDifference = (abs(mag1 - mag2)/ ((mag1 + mag2)/2)) * 100.0;
    varDifference = (abs(standD1 - standD2)/ ((standD1 + standD2)/2)) * 100.0;
    stdDifference = (abs(std1 - std2)/ ((std1 + std2)/2)) * 100.0;
    
    fout 
        << "Team Members: Sarah Villegas, Ali Fenton, Chanel Aquino\n" << endl
        << setw(30) << "two_a.txt" << setw(20) << "two_b.txt" << setw(20) << "% difference" << endl
        << "Mean" << setw(27) << mean1 << setw(20) << mean2 << setw(15) << meanDifference << endl
        << "Zero cross" << setw(14) << zeroCross(sound1) << setw(19) << zeroCross(sound2) << setw(23) << zeroCrossDifference <<endl
        << "Average power" << setw(16) << averagePower(sound1) << setw(19) << averagePower(sound2) << setw(17) << powerDifference << endl
        << "Average magnitude"  << setw(12) << avgMagnitude(sound1)  << setw(20) << avgMagnitude(sound2) << setw(17) << magnitudeDifference << endl
        << "Standard Deviation" << setw(12) << stdDev(standD1) << setw(20) << stdDev(standD1) << setw(15) << stdDifference << endl
        << "Variance" << setw(24) << var(sound1, mean1) << setw(20) << var(sound2, mean2) << setw(14) << varDifference << endl;
   
   
       
    // close files
    fin1.close();
    fin2.close();
    fout.close();
    return 0;
}
//*****************Average Magnitude*****************
double avgMagnitude(vector<double> vct){
    double avgMag, sum = 0.0, size = vct.size();
    
    for(int k = 0; k < size; k++){
        vct[k] = abs(vct[k]);
        sum += vct[k];
    }
    
    avgMag = sum / size;
    
    return avgMag;
}

//****************Average Power*****************
double averagePower(vector<double> vct)
{
    double
        avgPower,
        sum = 0.0, 
        size = vct.size();
    
    for(int k = 0; k < size; k++)
    {
        vct[k] *= vct[k];
        sum += vct[k];
    }
    
    avgPower = sum / size;
    
    return avgPower;
}

//*****************Mean (Average)*****************
double mean(vector<double> vct)
{
    double 
        average,
        sum = 0.0,
        size = vct.size();
    
    for (int i = 0; i < size; i++)
    {
        sum += vct[i];
    }
    
    average = sum / size;
    
    return average;
    
}


//*****************Standard Deviation*****************
double stdDev(double vct){
    return sqrt(vct);
}

//*********************************
double var(vector<double> vct, double& mean){
    
    double size = vct.size();
    double results = 0; 
    double squared = 0;  
    
    
    for(int i = 0; i < size; i++ ){
    
        squared = pow(vct[1] - mean ,2);
        results += squared;
   
    }
    results = results / size;

    return results;


}

//**********************Zero Cross**************************
int zeroCross(vector<double> vct){

    int totalZero = 0; 
    int size = vct.size(); 
    for(int i = 1; i < size; i++ ){
    
        if((vct[i - 1]  > 0) && (vct[i ] < 0)){
    
            totalZero++;
        }
        else if((vct[i - 1]  < 0) && (vct[i ] > 0)){
    
            totalZero++;
        }
    }

    return totalZero; 
    
}


