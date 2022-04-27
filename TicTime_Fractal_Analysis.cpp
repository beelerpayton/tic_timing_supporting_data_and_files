//
//  main.cpp
//  TicTime_Fractal_Analysis
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

#define PI 4.0*atan(1.0)
#define num_qs 150

double psdrand(int );
vector<double>initialize(double, double);
vector<vector<double>>walk(vector<double>, vector<double>, double, int);
vector<vector<double>>Sq(vector<vector<double>>, double);
vector<double>sq(vector<vector<double>>, double, double, double);
vector<double>power_law_fit(vector<vector<double>>);
double qmin (vector<vector<double>>);


int main()
{
    double T, dt;
    int decimal_places=0, breaker=-1, num=1;
    string filename, header, line, output_directory, data_string;
    vector<string>temp_vector;
    vector<string>ID;
    vector<double>raw_data;
    vector<vector<double>>walker;
    vector<vector<double>>trajectories;
    vector<vector<double>>structure_factor;
    vector<double>time;
    vector<double>temp;
    vector<double>power_fit_single(3);
    vector<vector<double> >power_fit_all;
    
    cout << "Max time (in seconds)? ";
    cin >> T;
    cout << "dt (in seconds)? ";
    cin >> dt;
    cout << "Filename? ";
    cin >> filename;
    cout << "Output directory? ";
    cin >> output_directory;
    cout << endl;
    cout << endl;
    
    while(breaker<0)
    {
        if (dt*pow(10,decimal_places)==floor(dt*pow(10,decimal_places)))
        {
            breaker=1;
        }
        else
        {
            decimal_places++;
        }
    }
    
    
    ifstream file;
    file.open(filename);
    
    if (!file.is_open())
    {
        cout <<"FILE ERROR" << endl;
        exit(EXIT_FAILURE);
    }
    
    while (!file.eof())
    {
        getline(file, line);
        
        stringstream temp_filename;
        temp_filename << output_directory;
        temp_filename << "/temp.txt";
        
        ofstream tempout (temp_filename.str().c_str());
        if (tempout.is_open())
        {
            tempout << line;
            tempout.close();
        }
        
        ifstream temp_file;
        temp_file.open(temp_filename.str().c_str());
        temp_vector.clear();
        raw_data.clear();
        structure_factor.clear();
        temp.clear();
        
        while (getline(temp_file, data_string, '\t'))
        {
            temp_vector.push_back(data_string);
        }
    
        ID.push_back(temp_vector[0]);
        cout << num << " " << temp_vector[0] << endl;
        for (int i=1; i<temp_vector.size(); i++)
        {
            raw_data.push_back(atof(temp_vector[i].c_str()));
        }
    
        time=initialize(T, dt);
        walker=walk(raw_data, time, dt, decimal_places);
        
        for (int i=0; i<walker.size(); i++)
        {
            temp.push_back(walker[i][1]);
        }
        trajectories.push_back(temp);
        structure_factor=Sq(walker, dt);
        power_fit_single=power_law_fit(structure_factor);
        power_fit_all.push_back(power_fit_single);
        num++;
    }
    
    stringstream out_filename;
    out_filename << output_directory;
    out_filename << "/Df.csv";
    
    ofstream outfile (out_filename.str().c_str());
    if (outfile.is_open())
    {
        outfile << "Patient ID,Df,R^2" << endl;
        for (int i=0; i<power_fit_all.size(); i++)
        {
            outfile << ID[i] << "," << abs(power_fit_all[i][1]) << "," << power_fit_all[i][2] << endl;
        }
        outfile.close();
    }
    
    stringstream out_filename2;
    out_filename2 << output_directory;
    out_filename2 << "/trajectories.csv";
    
    ofstream outfile2 (out_filename2.str().c_str());
    if (outfile2.is_open())
    {
        outfile2 << "time (s),";
        for (int i=0; i<ID.size(); i++)
        {
            outfile2 << ID[i] << ",";
        }
        outfile2 << endl;
        for (int i=0; i<time.size(); i++)
        {
            outfile2 << time[i] << ",";
            for (int j=0; j<ID.size(); j++)
            {
                outfile2 << trajectories[j][i] << ",";
            }
            outfile2 << endl;
        }
        outfile.close();
    }
}

//================

vector<double>initialize(double T, double dt)
{
    vector<double>v;
    int size = T/dt;
    
    for (int i=0; i<size+1; i++)
    {
        v.push_back(i*dt);
    }
    
    return v;
}

//================

vector<vector<double>>walk(vector<double>raw_data, vector<double>time, double dt, int decimals)
{
    vector<vector<double>>walker;
    vector<double>temp(2);
    vector<int>turns;
    
    for (int i=0; i<time.size(); i++)
    {
        turns.push_back(0);
        temp[0]=0;
        temp[1]=0;
        walker.push_back(temp);
    }
    
    for (int i=0; i<raw_data.size(); i++)
    {
        int repeats=0;
        for (int j=i+1; j<raw_data.size(); j++)
        {
            //cout << raw_data[i] << " " << raw_data[j] << endl;
            if (raw_data[i]==raw_data[j] && raw_data[i]!=0)
            {
                repeats++;
                raw_data[j]=raw_data[j]+(dt*double(repeats));
            }
            else
            {
                break;
            }
        }
    }
    
    for (int i=0; i<time.size(); i++)
    {
        time[i]=time[i]*pow(10,decimals);
        time[i]=round(time[i]);
        time[i]=time[i]*pow(10,-1*decimals);
    }
    
    for (int i=0; i<raw_data.size(); i++)
    {
        raw_data[i]=raw_data[i]*pow(10,decimals);
        raw_data[i]=round(raw_data[i]);
        raw_data[i]=raw_data[i]*pow(10,-1*decimals);
    }
    
    for (int i=0; i<raw_data.size(); i++)
    {
        if (raw_data[i]!=0)
        {
            for (int j=0; j<time.size(); j++)
            {
                if (raw_data[i]==time[j])
                {
                    turns[j]=1;
                }
            }
        }
    }
    
    double velocity=dt;
    
    for (int i=1; i<walker.size(); i++)
    {
        walker[i][0]=time[i];
        walker[i][1]=walker[i-1][1]+velocity;
        if (turns[i]==1)
        {
            velocity=-1*velocity;
        }
    }

    return walker;
}

//================

vector<vector<double>> Sq(vector<vector<double>>walker, double dt)
{
    vector<vector<double>>out;
    vector<vector<double>>unaveraged;
    vector<double>magnitude;
    vector<double>temp;
    vector<double>temp2(2);
    double r_monomer, q, q_max, dq, q_min, Rg;
    int orientations=180;
    
    Rg=qmin(walker);
    q_min=1/Rg;
    
    r_monomer=0.5*sqrt(pow(dt,2)+pow(dt,2));
    q_max=1/r_monomer;
    dq=(q_max-q_min)/num_qs;
    q=q_min;
    
    do
    {
        magnitude.push_back(q);
        q=q+dq;
    } while (q<=q_max);
    
    for (int j=0; j<orientations; j++)
    {
        temp=sq(walker, dq, j, r_monomer);
        unaveraged.push_back(temp);
    }
        
    double sum, avg;

    for (int j=0; j<num_qs; j++)
    {
        sum = 0;
        for (int k=0; k<orientations; k++)
        {
            sum=sum+unaveraged[k][j];
        }
        avg=sum/orientations;
        temp2[0]=magnitude[j];
        temp2[1]=avg;
        out.push_back(temp2);
    }
    
    return out;
}

//================

vector<double> sq(vector<vector<double> >part, double dq, double orient, double r_monomer)
{
    double P, qx, qy, cos_sum=0, sin_sum=0, dot, q, q_min=0.0000001;
    long double s_tot;
    double N=part.size();
    int counter=0, orientations=180;;
    vector<double>back(num_qs);

    q=q_min;
    P=orient*((2*PI)/double(orientations));

    do
    {
        sin_sum=0;
        cos_sum=0;
        qx = q*cos(P);
        qy = q*sin(P);

        for (int i=0; i<part.size(); i++)
        {
            dot=qx*part[i][0]+qy*part[i][1];
            sin_sum=sin_sum+sin(dot);
            cos_sum=cos_sum+cos(dot);
        }

        s_tot=(cos_sum*cos_sum)+(sin_sum*sin_sum);
        s_tot=s_tot/(N*N);
        back[counter]=s_tot;

        counter=counter+1;
        q=q+dq;
    } while (q<=1/r_monomer);

    return back;
}

//=========================

vector<double>power_law_fit(vector<vector<double>>input)
{
    vector<double>out(3);
    double sum_lnx, sum_lny, sum_lnx_lny, sum_lnx2, sum_lny2, N, a, b, R;
    vector<double>x;
    vector<double>y;
    
    for (int j=2; j<input.size(); j++)
    {
        x.push_back(input[j][0]);
        y.push_back(input[j][1]);
    }
    
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
                
    out[0]=a;
    out[1]=b;
    out[2]=R;
    
    return out;
}
//================

double qmin(vector<vector<double>>walker)
{
    double Rg=0, x_com, y_com, xsum=0, ysum=0, dx, dy, r, Rgsum=0;
    int N=0;
    
    for (int i=0; i<walker.size(); i++)
    {
        xsum=xsum+walker[i][0];
        ysum=ysum+walker[i][1];
        N++;
    }
    x_com=xsum/double(N);
    y_com=ysum/double(N);
    for (int i=0; i<walker.size(); i++)
    {
        dx=walker[i][0]-x_com;
        dy=walker[i][1]-y_com;
        r=sqrt((dx*dx)+(dy*dy));
        Rgsum=Rgsum+(r*r);
    }
    Rg=sqrt(Rgsum/double(N));
    return Rg;
}
