//
//  main.cpp
//  Tourette's_random_walk
//
//  Created by Payton Beeler on 11/14/19.
//  Copyright © 2019 Payton Beeler. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <float.h>
#include <fstream>
#include <sstream>
#include<vector>
#include<iostream>
#include <iomanip>

using namespace std;

#define mode 1 //0 = 0.1 dt; 1 = random dt
#define T 3000  //number of time steps
#define seed (int) 32103253564512457     //random number generator seed

vector<vector<double>>read(string, string);
vector<vector<double>>initialize(int);
vector<vector<double>>walk(vector<vector<double>>, vector<vector<double>>, int);
vector<vector<double>>normalize(vector<vector<double>>, int);
vector<vector<double>>MSD(vector<vector<double>>, int);
vector<vector<double>>ind_MSD(vector<vector<double>>, int);
vector<double>flight_time(vector<vector<double>>, int);
vector<vector<double>>ind_flight_time(vector<vector<double>>, int);
vector<vector<double>>histograms(vector<double>);
vector<vector<double>>ind_histograms(vector<vector<double>>, int);
vector<vector<double>>power_law_fit(vector<vector<double>>, int, string, string);
vector<vector<double>>power_law_fit_msd(vector<vector<double>>, int, string, string);
double psdrand(int );
vector<vector<double>>autocorrelation(vector<vector<double>>, int, string, string);


int main()
{
    string group, visit;
    vector<vector<double>>raw_data;
    vector<vector<double>>walkers;
    vector<vector<double>>normalized;
    vector<vector<double>>msd;
    vector<vector<double>>individual_msd;
    vector<vector<double>>histogram;
    vector<vector<double>>individual_histogram;
    vector<double>flights;
    vector<vector<double>>individual_flights;
    vector<vector<double>>power_law_hist;
    vector<vector<double>>power_law_msd;
    vector<vector<double>>auto_c;
    
    
    int rows;
    int counter;
    
    cout << "Enter group: ";
    cin >> group;
    cout << endl;
    cout << "Enter visit: ";
    cin >> visit;
    
    raw_data=read(group, visit);
    rows=raw_data.size();
    //rows=10;
    walkers=initialize(rows);
    walkers=walk(raw_data, walkers, rows);
    normalized=normalize(walkers, rows);
    auto_c=autocorrelation(walkers,rows,group,visit);
    msd=MSD(walkers, rows);
    individual_msd=ind_MSD(walkers, rows);
    //flights=flight_time(walkers,rows);
    //individual_flights=ind_flight_time(walkers, rows);
    //histogram=histograms(flights);
    //individual_histogram=ind_histograms(individual_flights, rows);
    power_law_hist=power_law_fit(individual_histogram, rows, group, visit);
    power_law_msd=power_law_fit_msd(individual_msd, rows, group, visit);

    
    ofstream outfile ("Positions.csv");
    if (outfile.is_open())
    {
        outfile << "Time, walkers' positions" << endl;
        for (int i=0; i<walkers.size(); i++)
        {
            for (int j=0; j<walkers[i].size(); j++)
            {
                outfile << walkers[i][j] << ", ";
            }
            outfile << endl;
        }
        outfile.close();
    }
    else cout << "Unable to open main file";
    
    ofstream outfile2 ("Normalized_positions.csv");
    if (outfile2.is_open())
    {
        outfile2 << "Time, walkers' positions" << endl;
        for (int i=0; i<normalized.size(); i++)
        {
            for (int j=0; j<normalized[i].size(); j++)
            {
                outfile2 << normalized[i][j] << ", ";
            }
            outfile2 << endl;
        }
        outfile2.close();
    }
    else cout << "Unable to open main file";
    
    ofstream outfile3 ("MSD_all.csv");
    if (outfile3.is_open())
    {
        outfile3 << "Time, MSD, MSD/100" << endl;
        for (int i=0; i<msd.size(); i++)
        {
            outfile3 << msd[i][0] << ", " << msd[i][1]  << ", " << msd[i][1]/100 << endl;
        }
        outfile3.close();
    }
    else cout << "Unable to open main file";
    
    ofstream outfile4 ("Flight_times_all.csv");
    if (outfile4.is_open())
    {
        outfile4 << "Flight Number, dt" << endl;
        for (int i=0; i<flights.size(); i++)
        {
            outfile4 << i+1 << ", " << flights[i] << endl;
        }
        outfile4.close();
    }
    else cout << "Unable to open main file";
    
    
    ofstream outfile5 ("Histogram_all.csv");
    if (outfile5.is_open())
    {
        int counter=1;
        
        outfile5 << "flight duration, number" << endl;
        for (int i=0; i<histogram.size(); i++)
        {
            if (histogram[i][1]==0)
            {
                
            }
            else if (histogram[i][0]==counter)
            {
                outfile5 << histogram[i][0] << ", " << histogram[i][1] << endl;
                counter++;
            } 
        }
        outfile5.close();
    }
    else cout << "Unable to open main file";
    
    ofstream outfile6 ("MSD_individual.csv");
    if (outfile6.is_open())
    {
        outfile6 << "Time, MSD" << endl;
        for (int i=0; i<individual_msd.size(); i++)
        {
            outfile6 << individual_msd[i][0] << ", ";
            
            for (int j=1; j<individual_msd[i].size(); j++)
            {
                outfile6 << individual_msd[i][j]/100 <<", ";
            }
            outfile6 << endl;
        }
        outfile6.close();
    }
    else cout << "Unable to open main file";
    
    
    ofstream outfile7 ("flight_time_individual.csv");
    if (outfile7.is_open())
    {
        for (int i=1; i<rows+1; i++)
        {
            outfile7 << i << ", ";
            
            for (int j=0; j<individual_flights.size(); j++)
            {
                if (individual_flights[j][i]>0)
                {
                    outfile7 << individual_flights[j][i] << ", ";
                }
            }
            outfile7 << endl;
        }
        outfile7.close();
    }
    else cout << "Unable to open main file";
    
    counter = 1;
    
    ofstream outfile8 ("histogram_individual.csv");
    if (outfile8.is_open())
    {
        for (int i=0; i<individual_histogram.size(); i++)
        {
            if (individual_histogram[i][0]==counter)
            {
                outfile8 << individual_histogram[i][0] << ",";
                for (int j=1; j<individual_histogram[i].size(); j++)
                {
                    //if (individual_histogram[i][j]>0)
                    //{
                        outfile8 << individual_histogram[i][j] << ",";
                    //}
                    //else
                    //{
                    //    outfile8 << ",";
                    //}
                }
                outfile8 << endl;
                counter++;
            }
        }
        outfile8.close();
    }
    else cout << "Unable to open main file";
    
    

    
    
    
    
    
}
//------------------------------reads the raw data file and turns it into a vector------------------------------
vector<vector<double>>read(string group, string visit)
{
    vector<vector<double>>v;
    vector<double>temp((T/10)+1);
    
    double data;
    
    
    
    stringstream ss;
    ss << "/Users/paytonbeeler/Desktop/Lab/Data/Levy flight/";
    ss << group;
    ss << "/";
    ss << group;
    ss << "_";
    ss << visit;
    ss << ".txt";
    
    ifstream test;
    test.open(ss.str().c_str());
    
    if (!test.is_open())
    {
        cout <<"FILE ERROR" << endl;
        exit(EXIT_FAILURE);
    }
    
    while (!test.eof())
    {
        for (int i=0; i<(T/10)+1; i++)
        {
            test >> data;
            temp.at(i)=data;
        }
        v.push_back(temp);
    }
    
    
    for (int i=0; i<v.size(); i++)
    {
        for (int j=0; j<5; j++)
        {
            cout << v[i][j] << " ";
        }
        cout << endl;
    }

    
    return v;
}

//------------------------------fills walker vector with timesteps and all zeroes------------------------------
vector<vector<double>>initialize(int row)
{
    vector<vector<double>>v;
    vector<double>temp(row+1);
    
    for (int i=0; i<T; i++)
    {
        temp.at(0)=(0.1*i)+0.1;
        
        for (int j=1; j<row+1; j++)
        {
            temp.at(j)=0;
        }
        
        v.push_back(temp);
    }
    
    return v;
}
//------------------------------------------------walks------------------------------------------------------
vector<vector<double>>walk(vector<vector<double>>old, vector<vector<double>>v, int row)
{
    vector<vector<double>>turns;
    vector<double>temp(row+1);
    int oldrow=0, turnscolumn=1, repeat=1, rand;
    
    cout << endl;
    cout << endl;
    cout << "walking...";
    
    for (int i=0; i<T; i++)
    {
        temp.at(0)=(0.1*i)+0.1;
        
        for (int j=1; j<row+1; j++)
        {
            temp.at(j)=0;
        }
        
        turns.push_back(temp);
    }

    do
    {
        for (int i=0; i<turns.size(); i++)
        {
            for (int j=1; j<(T/10)+1; j++)
            {
                if (turns[i][0]-old[oldrow][j]<0.01 && turns[i][0]-old[oldrow][j]>-1*0.01)
                {
                    repeat=0;
                    
                    if (turns[i][oldrow+1]==0)
                    {
                        turns[i][oldrow+1]=1;
                    }
                    else
                    {
                        while(old[oldrow][j]==old[oldrow][j+repeat])
                        {
                            repeat++;
                        }
                        
                        
                        if (mode==0)
                        {
                            turns[i+repeat][oldrow+1]=1;
                        }
                        
                        if (mode==1)
                        {
                            rand = 0;
                            
                            while (turns[i+rand][oldrow+1]==1)
                            {
                                rand = 10*psdrand(seed);
                            }
                            
                            turns[i+rand][oldrow+1]=1;
                        }
                        
                        
                        
                    }
                }
            }
        }
        
        oldrow++;
        
    } while (oldrow<row);
    
    
    int counter;
    
    for (int j=1; j<row+1; j++)
    {
        counter=-1;
        
        for (int i=0; i<v.size(); i++)
        {
            if (turns[i][j]!=0)
            {
                counter=counter*-1;
            }
            
            if (i==0)
            {
                v[i][j]=v[i][j]+counter;
            }
            
            else
            {
                v[i][j]=v[i-1][j]+counter;
            }
        }
    }
    
    cout << "finished" << endl;
    
    return v;
}
//------------------------------------------------------------normalized displacement------------------------------------------------------
vector<vector<double>>normalize(vector<vector<double>>w, int row)
{
    int max, min;
    
    for (int j=1; j<row+1; j++)
    {
        max=0;
        
        for (int i=0; i<w.size(); i++)
        {
            if (fabs(w[i][j])>max)
            {
                max=fabs(w[i][j]);
            }
        }
        
        for (int i=0; i<w.size(); i++)
        {
            w[i][j]=w[i][j]/max;
        }
    }
    
    return w;
}

//------------------------------------------------------------MSD------------------------------------------------------
vector<vector<double>>MSD(vector<vector<double>>w, int row)
{
    vector<double>temp(2);
    vector<vector<double>>v;
    double sum;
    
    for (int j=0; j<w.size(); j++)
    {
        sum=0;
        
        for (int i=1; i<row+1; i++)
        {
            sum=sum+(w[j][i]*w[j][i]);
        }
        
        temp.at(0)=w[j][0];
        temp.at(1)=sum/row;
        
        v.push_back(temp);
    }

    return v;
}

//----------------------------------------------------individual MSD------------------------------------------------------
vector<vector<double>>ind_MSD(vector<vector<double>>w, int row)
{
    vector<double>temp(row+1);
    vector<vector<double>>v;
    double sum;
    
    for (int j=0; j<w.size(); j++)
    {
        temp.at(0)=w[j][0];
        for (int i=1; i<row+1; i++)
        {
            sum=(w[j][i]*w[j][i]);
            temp.at(i)=sum;
        }
        
        v.push_back(temp);
    }

    return v;
}

//------------------------------------------------------------duration between turns------------------------------------------------------
vector<double>flight_time(vector<vector<double>>w, int row)
{
    vector<double>v;
    vector<double>vcopy;
    vector<vector<double>>copy;
    vector<double>temp(row+1);
    
    
    for (int i=0; i<w.size(); i++)
    {
        for (int j=0; j<w[i].size(); j++)
        {
            temp.at(j)=w[i][j];
        }
        copy.push_back(temp);
    }
    
    for (int i=1; i<row+1; i++)
    {
        for (int j=0; j<w.size(); j++)
        {
            if (j==0)
            {
                copy[j][i]=(w[j][i]-0)/(w[j][0]-0);
            }
            
            else
            {
              //  cout << w[j][i] << ", " << w[j][0] << ", "<< w[j-1][i] << ", "<< w[j-1][0] << endl;
                
                copy[j][i]=(w[j][i]-w[j-1][i])/(w[j][0]-w[j-1][0]);

            }
        }
    }
    
    
    for (int i=0; i<copy.size(); i++)
    {
       // cout << copy[i][0] << ", " << copy[i][1] << endl;
    }
    
    
    for (int i=1; i<row+1; i++)
    {
        for (int j=1; j<copy.size(); j++)
        {
            if (fabs(copy[j][i]-copy[j-1][i])>10)
            {
                v.push_back(copy[j][0]);
            }
        }
    }
    
    for (int i=0; i<v.size(); i++)
    {
        vcopy.push_back(v[i]);
    }
    
    
    for (int i=0; i<v.size(); i++)
    {
        if (i==0 || vcopy[i]<vcopy[i-1])
        {
            v[i]=vcopy[i];
        }
        else
        {
            v[i]=vcopy[i]-vcopy[i-1];
        }
    }
    
    return v;
}
//------------------------------------------------------------ind duration between turns------------------------------------------------------
vector<vector<double>>ind_flight_time(vector<vector<double>>w, int row)
{
    vector<vector<double>>v;
    vector<vector<double>>copy;
    vector<double>temp(row+1);
    vector<double>temp2(row+1);
    
    
    for (int i=0; i<w.size(); i++)
    {
        for (int j=0; j<w[i].size(); j++)
        {
            temp.at(j)=w[i][j];
            
        }
        copy.push_back(temp);
        
        temp2.at(0)=w[i][0];
        for (int j=1; j<w[i].size(); j++)
        {
            temp2.at(j)=0;
        }
        v.push_back(temp2);
    }
    
    for (int i=1; i<row+1; i++)
    {
        for (int j=0; j<w.size(); j++)
        {
            if (j==0)
            {
                copy[j][i]=(w[j][i]-0)/(w[j][0]-0);
            }
            
            else
            {
              //  cout << w[j][i] << ", " << w[j][0] << ", "<< w[j-1][i] << ", "<< w[j-1][0] << endl;
                
                copy[j][i]=(w[j][i]-w[j-1][i])/(w[j][0]-w[j-1][0]);

            }
        }
    }
    
    
    for (int i=0; i<copy.size(); i++)
    {
       // cout << copy[i][0] << ", " << copy[i][3] << endl;
    }
    
    
    for (int i=1; i<copy.size(); i++)
    {
        for (int j=1; j<row+1; j++)
        {
            if (fabs(copy[i][j]-copy[i-1][j])>10)
            {
                v[i][j]=(copy[i][0]);
            }
        }
    }
    
    double last, dt;
    
    for (int i=1; i<row+1; i++)
    {
        last = 0;
        for (int j=0; j<v.size(); j++)
        {
            if (v[j][i]>0)
            {
                dt=v[j][i]-last;
                last=v[j][i];
                v[j][i]=dt;
            }
        }
    }
    
    return v;
}

//------------------------------------------------------------histogram of flight times------------------------------------------------------
vector<vector<double>>histograms(vector<double>dt)
{
    cout << "histograms...";
    vector<vector<double>>v;
    vector<double>temp(2);
    int counter;
    
    for (int i=0; i<T; i++)
    {
        temp.at(0)=(i+1)*0.1;
        temp.at(1)=0;
        v.push_back(temp);
    }
    
    for (int i=0; i<v.size(); i++)
    {
        counter=0;
        
        for (int j=0; j<dt.size(); j++)
        {
            if (dt[j]-v[i][0]<0.001 && dt[j]-v[i][0]>-0.001)
            {
                counter++;
            }
        }
        v[i][1]=counter;
    }
    
    cout << "finsished" << endl;
    
    return v;
}
//-------------------------------random number generator------------------------------------------
double psdrand(int iseed)
{
    int i, j, k, inx;
    double ran_num;
    static const int ndim = 55, m10 = 1000000000, is = 21, ir = 30;
    static const double base = 1.0E9;
    static int jrand, istack[58];
    static bool init = false;
    
    if((!init) || (iseed < 0))
    {
        iseed = abs(iseed);
        istack[ndim] = iseed;
        j = iseed;
        k = 1;
        
        for(i = 1; i <= (ndim - 1); ++i)
        {
            inx = i*is - int((double)(i*is)/(double)(ndim))*ndim;
            istack[inx] = k;
            k = j - k;
            if(k < 0) {k += m10;}
            j = istack[inx];
        }
        
        for(j = 1; j <= 3; ++j)
        {
            for(i = 1; i <= ndim; ++i)
            {
                inx = i + ir - int((double)(i+ir)/(double)(ndim))*ndim;
                istack[i] -= istack[inx+1];
                if(istack[i] < 0) {istack[i] += m10;}
            }
        }
        jrand = 0;
        init = true;
    }
    
    jrand += 1;
    
    if(jrand > ndim)
    {
        for(i = 1; i <= ndim; ++i)
        {
            inx = i + ir - ((int)((double)(i+ir)/(double)(ndim)))*ndim;
            istack[i] -= istack[inx+1];
            if(istack[i] < 0) {istack[i] += m10;}
        }
        jrand = 1;
    }
    
    ran_num = ((double)istack[jrand]) / base;
    
    return ran_num;
}
//------------------------------------------------------------ind histogram of flight times------------------------------------------------------
vector<vector<double>>ind_histograms(vector<vector<double>>flight_times, int rows)
{
    cout << "individual histograms..." << endl;
    
    vector<vector<double>>v;
    vector<double>temp(rows+1);
    int counter;
    double target;
    
    for (int i=0; i<T; i++)
    {
        temp.at(0)=(i+1)*0.1;
        for (int j=1; j<rows+1; j++)
        {
            temp.at(j)=0;
        }
        
        v.push_back(temp);
    }
    

    for (int i=1; i<rows+1; i++)
    {
        cout << i << endl;
        for (int j=0; j<flight_times.size(); j++)
        {
            target=flight_times[j][0];
            counter=0;
            
            for (int k=0; k<flight_times.size(); k++)
            {
                if (abs(flight_times[k][i]-target)<0.01)
                {
                    counter++;
                }
            }
            v[j][i]=counter;
        }
    }
    
    return v;
}
//---------------------------power law fitting of histogram---------------------------------------
vector<vector<double>>power_law_fit(vector<vector<double>>hist, int row, string group, string visit)
{
    vector<vector<double>>out;
    vector<double>x;
    vector<double>y;
    double sum_lnx, sum_lny, sum_lnx_lny, sum_lnx2, sum_lny2, N, a, b, R, total_tics;
    vector<double>temp;
    vector<string>all_pts;
    vector<int>checker;
    vector<string>final_pts;
    string patient;
    int counter=1;
    
    stringstream ss;
    ss << "/Users/paytonbeeler/Desktop/Lab/Data/Levy flight/";
    ss << group;
    ss << "/";
    ss << group;
    ss << "_";
    ss << visit;
    ss << "_names.txt";
    
    ifstream test;
    test.open(ss.str().c_str());
    
    if (!test.is_open())
    {
        cout <<"FILE ERROR" << endl;
        exit(EXIT_FAILURE);
    }
    
    while (!test.eof())
    {
        test >> patient;
        all_pts.push_back(patient);
        checker.push_back(0);
    }
    

    for (int column1=1; column1<row+1; column1++)
    {
        if (checker[column1-1]==0)
        {
            temp.push_back(column1);
            for (int column2=column1; column2<row+1; column2++)
            {
                counter=1;
                if (all_pts[column1-1]==all_pts[column2-1] && checker[column2-1]==0)
                {
                    checker[column1-1]=1;
                    checker[column2-1]=1;

                    for (int i=0; i<hist.size(); i++)
                    {
                        if (hist[i][0]==counter && hist[i][column2]!=0)
                        {
                            //x.push_back(counter);
                            x.push_back(counter+1);
                            y.push_back(hist[i][column2]);
                            counter++;
                        }
                        else if (hist[i][0]==counter && hist[i][column2]==0)
                        {
                            counter++;
                        }
                    }
                    patient = all_pts[column2-1];
                }
            }

            total_tics = 0;
        
            for (int i=0; i<y.size(); i++)
            {
                total_tics = total_tics + y[i];
            }
            //cout << total_tics << endl;
            
            for (int i=0; i<y.size(); i++)
            {
                //y[i] = y[i]/total_tics;
            }
            temp.push_back(total_tics);

        
            sum_lnx=0;sum_lny=0;sum_lnx_lny=0;sum_lnx2=0;sum_lny2=0;N=0;
        
            for (int i=0; i<x.size(); i++)
            {
                sum_lnx=sum_lnx+log(x[i]);
                sum_lny=sum_lny+log(y[i]);
                sum_lnx_lny=sum_lnx_lny+(log(x[i])*log(y[i]));
                sum_lnx2=sum_lnx2+(log(x[i])*log(x[i]));
                sum_lny2=sum_lny2+(log(y[i])*log(y[i]));
                N = N+1;
            }
        
            b = ((sum_lnx_lny)-((1/N)*sum_lnx*sum_lny))/((sum_lnx2)-((1/N)*((sum_lnx)*(sum_lnx))));
            a = exp(((1/N)*sum_lny)-(b*(sum_lnx/N)));
            R = (((sum_lnx_lny)-((1/N)*sum_lnx*sum_lny))*((sum_lnx_lny)-((1/N)*sum_lnx*sum_lny)))/((sum_lnx2-((sum_lnx*sum_lnx)/N))*(sum_lny2-((sum_lny*sum_lny)/N)));
                
            temp.push_back(a);
            temp.push_back(b);
            temp.push_back(R);
            out.push_back(temp);
            final_pts.push_back(patient);
            x.clear();
            y.clear();
            temp.clear();
        }
    }
    
    
    ofstream outfile9 ("histogram_individual_fits.csv");
    if (outfile9.is_open())
    {
        outfile9 << "patient #,total tics,a,b,R^2" << endl;
        for (int i=0; i<out.size(); i++)
        {
            outfile9 << final_pts[i] << "," << out[i][1] << "," << out[i][2] << "," << out[i][3] << "," << out[i][4] << endl;
        }
        outfile9.close();
    }
    else cout << "Unable to open main file";

    
    
    return out;
}
//---------------------------power law fitting of MSD vs t---------------------------------------
vector<vector<double>>power_law_fit_msd(vector<vector<double>>v, int row, string group, string visit)
{
        vector<vector<double>>out;
    vector<double>x;
    vector<double>y;
    double sum_lnx, sum_lny, sum_lnx_lny, sum_lnx2, sum_lny2, N, a, b, R, total_tics;
    vector<double>temp;
    vector<string>all_pts;
    vector<int>checker;
    vector<string>final_pts;
    string patient;
    //int counter=1;
    
    stringstream ss;
    ss << "/Users/paytonbeeler/Desktop/Lab/Data/Levy flight/";
    ss << group;
    ss << "/";
    ss << group;
    ss << "_";
    ss << visit;
    ss << "_names.txt";
    
    ifstream test;
    test.open(ss.str().c_str());
    
    if (!test.is_open())
    {
        cout <<"FILE ERROR" << endl;
        exit(EXIT_FAILURE);
    }
    
    while (!test.eof())
    {
        test >> patient;
        all_pts.push_back(patient);
        checker.push_back(0);
    }
    
    int counter;
    vector<double>total;

    for (int column1=1; column1<row+1; column1++)
    {
        //cout << all_pts[column1-1] << endl;
        
        if (checker[column1-1]==0)
        {
            counter = 1;
            for (int column2=column1+1; column2<row+1; column2++)
            {
                if (all_pts[column1-1]==all_pts[column2-1] && checker[column2-1]==0)
                {
                    checker[column1-1]=1;
                    checker[column2-1]=1;

                    if (counter==1)
                    {
                        counter++;
                        for (int i=0; i<v.size(); i++)
                        {
                            //x.push_back(v[i][0]);
                            total.push_back((v[i][column1]/100)+(v[i][column2]/100));
                        }
                    }
                    else
                    {
                        counter++;
                        checker[column1-1]=1;
                        
                        for (int i=0; i<v.size(); i++)
                        {
                            total[i] = total[i] + (v[i][column2]/100);
                        }
                    }
                    
                    patient = all_pts[column2-1];
                    
                }
                else if (column2==row && counter==1)
                {
                    for (int i=0; i<v.size(); i++)
                    {
                        //x.push_back(v[i][0]);
                        total.push_back((v[i][column1]/100));
                    }
                    patient = all_pts[column1-1];
                }
            }
            
            for (int i=0; i<v.size(); i++)
            {
                if (total[i]!=0)
                {
                    x.push_back(v[i][0]);
                    y.push_back(total[i]/counter);
                }
            }
            
            //cout << patient << " " << counter << endl;
            
            
            sum_lnx=0;sum_lny=0;sum_lnx_lny=0;sum_lnx2=0;sum_lny2=0;N=0;
        
            for (int i=0; i<x.size(); i++)
            {
                sum_lnx=sum_lnx+log(x[i]);
                sum_lny=sum_lny+log(y[i]);
                sum_lnx_lny=sum_lnx_lny+(log(x[i])*log(y[i]));
                sum_lnx2=sum_lnx2+(log(x[i])*log(x[i]));
                sum_lny2=sum_lny2+(log(y[i])*log(y[i]));
                N = N+1;
            }
        
            b = ((sum_lnx_lny)-((1/N)*sum_lnx*sum_lny))/((sum_lnx2)-((1/N)*((sum_lnx)*(sum_lnx))));
            a = exp(((1/N)*sum_lny)-(b*(sum_lnx/N)));
            R = (((sum_lnx_lny)-((1/N)*sum_lnx*sum_lny))*((sum_lnx_lny)-((1/N)*sum_lnx*sum_lny)))/((sum_lnx2-((sum_lnx*sum_lnx)/N))*(sum_lny2-((sum_lny*sum_lny)/N)));
                
            temp.push_back(counter);
            temp.push_back(a);
            temp.push_back(b);
            temp.push_back(R);
            out.push_back(temp);
            final_pts.push_back(patient);
            x.clear();
            y.clear();
            temp.clear();
            total.clear();
        }
    }
    
    
    ofstream outfile9 ("MSD_individual_fits.csv");
    if (outfile9.is_open())
    {
        outfile9 << "patient #,counts,a,b,R^2" << endl;
        for (int i=0; i<out.size(); i++)
        {
            outfile9 << final_pts[i] << "," << out[i][0] << "," << out[i][1] << "," << out[i][2] << "," << out[i][3] << endl;
        }
        outfile9.close();
    }
    else cout << "Unable to open main file";

    
    
    return out;
}

//---------------------------autocorrelation---------------------------------------
vector<vector<double>>autocorrelation(vector<vector<double>>v, int row, string group, string visit)
{
    vector<vector<double>>out;
    vector<string>all_pts;
    vector<int>checker;
    vector<string>final_pts;
    string patient, sessions;
    vector<double>t;
    vector<double>x;
    vector<double>zeros;
    double total, average_x, num_x, variance, autocovariance, lag_d, counter;
    int breaker;
    vector<double>temp;
    vector<vector<double>>autocorrelation_all;
    vector<double>avg_zeros;
    vector<string>pts;
    vector<string>session;
    
    for (int j=0; j<v.size(); j++)
    {
        for (int i=0; i<row+1; i++)
        {
            temp.push_back(0);
        }
        autocorrelation_all.push_back(temp);
        temp.clear();
    }
        
    stringstream ss;
    ss << "/Users/paytonbeeler/Desktop/Lab/Data/Levy flight/";
    ss << group;
    ss << "/";
    ss << group;
    ss << "_";
    ss << visit;
    ss << "_name_session.txt";
    
    ifstream test;
    test.open(ss.str().c_str());
    
    if (!test.is_open())
    {
        cout <<"FILE ERROR" << endl;
        exit(EXIT_FAILURE);
    }
    
    while (!test.eof())
    {
        test >> patient;
        test >> sessions;
        all_pts.push_back(patient);
        session.push_back(sessions);
        checker.push_back(0);
    }
    
    for (int column=1; column<row+1; column++)
    {
        cout << column << endl;
        num_x = 0;
        total=0;
        for (int i=0; i<v.size(); i++)
        {
            t.push_back(v[i][0]);
            x.push_back(v[i][column]/10);
            num_x = num_x+1;
            total=total+(v[i][column]/10);
        }
        average_x=total/num_x;
        
        variance = 0;
        
        for (int i=0; i<x.size(); i++)
        {
            variance = variance + ((x[i]-average_x)*(x[i]-average_x));
        }
        variance = variance/num_x;
        
        for (int lag=1; lag<v.size()-1; lag++)
        {
            autocovariance=0;
            lag_d=lag;
            
            for (int i=0; i<v.size()-lag; i++)
            {
                autocovariance=autocovariance+((x[i]-average_x)*(x[i+lag]-average_x));
            }
            autocovariance=autocovariance/num_x;
            autocorrelation_all[lag][column]=autocovariance/variance;
            autocorrelation_all[lag][0]=lag_d/10;
        }
        
        t.clear();
        x.clear();
        
        breaker=0;
        
        for (int i=1; i<autocorrelation_all.size(); i++)
        {
            if (autocorrelation_all[i][column]<exp(-1) && breaker==0)
            {
                zeros.push_back(autocorrelation_all[i][0]);
                breaker=1;
            }
        }
    }
    
    
//    for (int column1=1; column1<row+1; column1++)
//    {
//        //cout << all_pts[column1-1] << endl;
//
//        if (checker[column1-1]==0)
//        {
//            counter = 1;
//            total=zeros[column1-1];
//
//            for (int column2=column1+1; column2<row+1; column2++)
//            {
//                if (all_pts[column1-1]==all_pts[column2-1] && checker[column2-1]==0)
//                {
//                    checker[column1-1]=1;
//                    checker[column2-1]=1;
//                    counter++;
//                    total=total+zeros[column2-1];
//
//                }
//            }
//            pts.push_back(all_pts[column1-1]);
//            avg_zeros.push_back(total/counter);
//        }
//
//    }

    ofstream outfile9 ("ACF_squared.csv");
    if (outfile9.is_open())
    {
        //outfile9 << "patient #,counts,a,b,R^2" << endl;
        
        for (int i=1; i<autocorrelation_all.size()-1; i++)
        {
            for (int j=0; j<autocorrelation_all[i].size(); j++)
            {
                if (j>0)
                {
                    //outfile9 << autocorrelation_all[i][j]*autocorrelation_all[i][j] << ",";
                    outfile9 << autocorrelation_all[i][j] << ",";
                }
                else
                {
                    outfile9 << autocorrelation_all[i][j] << ",";
                }
            }
            outfile9 << endl;
        }
        outfile9.close();
    }
    else cout << "Unable to open main file";
    
    ofstream outfile ("ACF_e-folding.csv");
    if (outfile.is_open())
    {
        outfile << "Patient #, ACF e-folding" << endl;
        for (int i=0; i<all_pts.size(); i++)
        {
            outfile << all_pts[i] << "," << session[i] << "," << zeros[i] << endl;
        }
        outfile.close();
    }
    else cout << "Unable to open main file";

    
    
    return out;
}
