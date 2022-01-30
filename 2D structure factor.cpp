//
//  main.cpp
//  2D stucture factor
//
//  Created by Payton Beeler on 6/1/21.
//  Copyright © 2021 Payton Beeler. All rights reserved.
//

#include <ctime>
#include <iostream>
#include <math.h>
#include <iomanip>
#include <vector>
#include <cmath>
#include <fstream>
#include <cstdio>
#include <sstream>

using namespace std;

#define mode 0      //0 = single file; 1 = multiple files (check limits)
#define q_min 0.0000001
#define r_monomer 0.07071067812
#define num_qs 150
#define orientations 180
//#define seed 5078870851
#define PI 4.0*atan(1.0)


vector<vector<double> >build_vector ();
vector<double> sq(vector<vector<double> >, double, double);
double rg (vector<vector<double> >);
vector<vector<double> >build_loop (int);
vector<double>power_law_fit(vector<double>, vector<double>);
double psdrand(int);

int main()
{

    
    double q, s_max=1/r_monomer;
    double dq=s_max/num_qs;
    vector<vector<double> >cluster;
    vector<vector<double> >output;
    vector<vector<double> >final_out;
    vector<vector<double> >power_fit_all;
    vector<double>magnitude;
    vector<double>x;
    vector<double>y;
    vector<double>power_fit_single(3);
    vector<double>averaged;
    vector<double>temp;
    vector<double>temp2(2);
    double Rg;
    int first_file, last_file;
    cout << "First file: ";
    cin >> first_file;
    cout << "Last file: ";
    cin >> last_file;
    cout << "dq: " << dq << endl;
    q=q_min;
    
    
    //srand(seed);
    
    do
    {
        magnitude.push_back(q);
        q=q+dq;
    } while (q<=1/r_monomer);
    

    
    for (int i=first_file; i<=last_file; i++)
    {
        //srand(seed);
        cluster=build_loop(i);
        Rg=rg(cluster);
        cout << i << ", Radius of gyration: " << Rg << endl;
        
        output.clear();
        averaged.clear();
        
        for (int j=0; j<orientations; j++)
        {
            temp=sq(cluster, dq, j);
            output.push_back(temp);
            //cout << j << endl;
        }

        double sum, avg;

        for (int j=0; j<num_qs; j++)
        {
            sum = 0;
            for (int k=0; k<orientations; k++)
            {
                sum=sum+output[k][j];
            }
            avg=sum/orientations;
            averaged.push_back(avg);
        }
        final_out.push_back(averaged);
    }
    
    for (int i=0; i<final_out.size(); i++)
    {
        x.clear();
        y.clear();
        for (int j=2; j<final_out[i].size(); j++)
        {
            x.push_back(magnitude[j]);
            y.push_back(final_out[i][j]);
        }
        power_fit_single=power_law_fit(x,y);
        power_fit_all.push_back(power_fit_single);
    }
    
//    string group, visit, patient;
//    vector<string>all_pts;
//    vector<string>pts;
//    vector<int>checker;
//    vector<double>avg_Df;
//    int rows;
//
//    cout << endl;
//    cout << "Enter group: ";
//    cin >> group;
//    cout << "Enter visit: ";
//    cin >> visit;
//
//    stringstream ss;
//    ss << "/Users/paytonbeeler/Desktop/Lab/Data/Levy flight/";//Amy_Data/";
//    ss << group;
//    ss << "/";
//    ss << group;
//    ss << "_";
//    ss << visit;
//    ss << "_names.txt";
//
//    ifstream test;
//    test.open(ss.str().c_str());
//
//    if (!test.is_open())
//    {
//        cout <<"FILE ERROR" << endl;
//        exit(EXIT_FAILURE);
//    }
//
//    while (!test.eof())
//    {
//        test >> patient;
//        all_pts.push_back(patient);
//        checker.push_back(0);
//    }
//
//    rows=last_file;
//
//    double number, sum;
//
//    for (int column1=0; column1<rows; column1++)
//    {
//        if (checker[column1]==0)
//        {
//            number=0;
//            sum=0;
//            for (int column2=column1; column2<rows; column2++)
//            {
//                if (all_pts[column2]==all_pts[column1] && checker[column2]==0)
//                {
//                    number=number+1;
//                    checker[column2]=1;
//                    sum=sum+0;//power_fit_all[column2][1];
//                }
//            }
//            avg_Df.push_back(sum/number);
//            pts.push_back(all_pts[column1]);
//        }
//    }
    
    
    std::ofstream outfile ("S(q).csv");
    if (outfile.is_open())
    {
        for (int i=0; i<magnitude.size(); i++)
        {
            outfile << magnitude[i] << ",";
        }
        outfile << endl;
        
        for (int i=0; i<final_out.size(); i++)
        {
            for (int j=0; j<final_out[i].size(); j++)
            {
                outfile << final_out[i][j] << ",";
            }
            outfile << endl;
        }
        outfile.close();
    }
        
    std::ofstream outfile2 ("Df.csv");
    if (outfile2.is_open())
    {
        for (int i=0; i<power_fit_all.size(); i++)
        {
            for (int j=0; j<power_fit_all[i].size(); j++)
            {
                outfile2 << abs(power_fit_all[i][j]) << ",";
            }
            outfile2 << endl;
            
        }
        outfile2.close();
    }
    
//    std::ofstream outfile3 ("avg_Df.csv");
//    if (outfile3.is_open())
//    {
//        for (int i=0; i<avg_Df.size(); i++)
//        {
//            outfile3 << pts[i] << "," << avg_Df[i] << endl;
//
//        }
//        outfile3.close();
//    }
}
//-----------------------------------read csv file and turn it into a vector-------------------------------
vector<vector<double> >build_loop (int i)
{
    vector<vector<double> >v;
    int size, row;
    vector<double>temp(2);
    ifstream test;
    
    stringstream ss;
    ss << "/Users/paytonbeeler/Desktop/time series/walker_";
    ss << i;
    ss << ".txt";
    test.open(ss.str().c_str());
    
    if (!test.is_open())
    {
        cout <<"FILE ERROR" << endl;
        exit(EXIT_FAILURE);
    }
    
    double x, y;
    
    while(!test.eof())
    {
        test >> x;
        test >> y;

        temp.at(0)=x;
        temp.at(1)=y;
        
        v.push_back(temp);
    }
    
    size = v.size();
    row = size-1;
    
    if (v.size() > row)
    {
        v.erase(v.begin() + row);
    }
    
    return v;
}
//-------------------------------finds sq at one orientation------------------------------------------
vector<double> sq(vector<vector<double> >part, double dq, double orient)
{
    double P, qx, qy, cos_sum=0, sin_sum=0, dot, q;
    long double s_tot;
    double N=part.size();
    int counter=0;
    vector<double>back(num_qs);
    
    q=q_min;
    //P=2*PI*(double)rand()/RAND_MAX;
    P=orient*((2*PI)/orientations);

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
//------------------------------------------radius of gyration-----------------------------------------
double rg (vector<vector<double> >p)
{
    double xcom, ycom, xsum=0, ysum=0, dx, dy;
    double rgsum=0, RG, distance;
    int mass=p.size();
    
    for (int i=0; i<p.size(); i++)
    {
        xsum=xsum+p[i][0];
        ysum=ysum+p[i][1];
    }
    
    xcom=xsum/mass;
    ycom=ysum/mass;
    
    //cout << xcom << ", " << ycom << endl;
    
    for (int i=0; i<p.size(); i++)
    {
        dx=p[i][0]-xcom;
        dy=p[i][1]-ycom;
        distance=(dx*dx)+(dy*dy);
        rgsum=rgsum+distance;
    }
    RG=rgsum/mass;
    RG=sqrt(RG);
    
    return RG;
}
//---------------------------power law fitting---------------------------------------
vector<double>power_law_fit(vector<double>x, vector<double>y)
{
    vector<double>out(3);
    double sum_lnx, sum_lny, sum_lnx_lny, sum_lnx2, sum_lny2, N, a, b, R;
    
    
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






