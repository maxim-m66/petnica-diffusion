//A planet divided in cells
//Calculating the temperature for every cell
//diffusion between cells

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <cmath>
#include <time.h>

using namespace std;

void srand (unsigned seed);

int main ()

{
    srand(time(NULL));
    int n=40, tt=8015, xx=50;
    double C = 70.0;                    //arbitrary constant converting energy to temp (set by us)
    int i, j, t, x, Nbx, Nwx, Nbsum, Nwsum, Nb[tt], Nw[tt];
    double Tb, Tw, pi, S, d, s, ls, bs, dt, k, Rb, Rw, R3, Txsum, Tsum;
    double Pb[n][n], Pw[n][n], T[n][n], a[n][n], E1[n][n], E2[n][n], Esum[n][n], Ed[n][n], Ef[n][n], Tf[n][n], Tavg[tt], L[tt], Tx[xx];
    dt=10000;
    k=5000000;                         // ration dt/k is important, has to be <1
    Tb=283.15;                          // ideal temp for black daisies
    Tw=313.15;                          // ideal temp for white daisies
    pi=3.141593;
    S=14;                               // standard deviation for temp
    d=175000000000;                     // distance planet - star
    s=5.670367*pow(10, -8);             // Stefan Boltzmann constant
    ls=3.828*pow(10, 26);               // luminosity of sun now

    ifstream lum ("luminositynew.csv");
    if (lum.is_open())
        {
            while (!lum.eof())
            {
                for (i=0;i<tt;i+=1)
                        lum>>L[i];
            }
            lum.close();
        }
    else {cout<<"Could not find file"<<endl; return 0;}

    for (i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            //initial values
            a[i][j]=0.5;
            T[i][j]=pow((L[0]*ls*(1-a[i][j]))/(4*d*d*pi*s),0.25);
            E1[i][j]=T[i][j]*C;
        }
    }
    ofstream temp, numb, numw, ev;
    temp.open ("Final_Daisies.csv");
    numb.open ("Final_Nblack.csv");
    numw.open ("Final_Nwhite.csv");

    for (t=0;t<tt;t++)
    {
        Nbsum=0;
        Nwsum=0;
        Txsum=0;
        for (x=0;x<xx;x++)
        {
            Nbx=0;
            Nwx=0;
            Tsum=0;
            for (i=0;i<n;i++)
            {
                for (j=0;j<n;j++)
                {
                    bs=4;      //number of neighbouring cells

                    if (x%2==0)                                                                         //decides which matrix to use
                    {
                        Esum[i][j]=E1[i-1][j]+E1[i][j-1]+E1[i][j+1]+E1[i+1][j];                         //normal cells
                        if (i==0)    {bs-=1; Esum[i][j]-=E1[i-1][j];}                                   //border cells
                        if (j==0)    {bs-=1; Esum[i][j]-=E1[i][j-1];}
                        if (i==n-1)  {bs-=1; Esum[i][j]-=E1[i+1][j];}
                        if (j==n-1)  {bs-=1; Esum[i][j]-=E1[i][j+1];}
                        Ed[i][j]=(Esum[i][j]-bs*E1[i][j])*(1/bs)*(1-exp(-dt/k));                        //total energy exchange
                        Ef[i][j]=E1[i][j]+(L[t]*ls*(1-a[i][j]))/(4*d*d*pi)-s*pow(T[i][j],4)+Ed[i][j];   //total cell energy
                        E2[i][j]=T[i][j]*C;
                    }
                    else
                    {
                    Esum[i][j]=E2[i-1][j]+E2[i][j-1]+E2[i][j+1]+E2[i+1][j];
                        if (i==0)    {bs-=1; Esum[i][j]-=E2[i-1][j];}
                        if (j==0)    {bs-=1; Esum[i][j]-=E2[i][j-1];}
                        if (i==n-1)  {bs-=1; Esum[i][j]-=E2[i+1][j];}
                        if (j==n-1)  {bs-=1; Esum[i][j]-=E2[i][j+1];}
                        Ed[i][j]=(Esum[i][j]-bs*E2[i][j])*(1/bs)*(1-exp(-dt/k));
                        Ef[i][j]=E2[i][j]+(L[t]*ls*(1-a[i][j]))/(4*d*d*pi)-s*pow(T[i][j],4)+Ed[i][j];
                        E1[i][j]=T[i][j]*C;
                    }
                    T[i][j]=Ef[i][j]/C;

                    Pb[i][j]=exp(-0.5*pow(((T[i][j]-Tb)/S),2));
                    Pw[i][j]=exp(-0.5*pow(((T[i][j]-Tw)/S),2));

                    Rb=((double) rand() / (RAND_MAX));
                    Rw=((double) rand() / (RAND_MAX));

                    if      (T[i][j]<273.15 || T[i][j]>323.15) a[i][j]=0.5;                 //decides if daisy is going to grow or not
                    else if (Pb[i][j]>=Rb && Pw[i][j]<Rw) {a[i][j]=0.25; Nbx+=1;}
                    else if (Pb[i][j]<Rb && Pw[i][j]>=Rw) {a[i][j]=0.75; Nwx+=1;}
                    else if (Pb[i][j]<Rb && Pw[i][j]<Rw) a[i][j]=0.5;
                    else if (Pb[i][j]>Rb && Pw[i][j]>Rw)
                    {
                        R3=((double) rand() / (RAND_MAX));
                        if      (Pb[i][j]/(Pb[i][j]+Pw[i][j])>R3) {a[i][j]=0.25; Nbx+=1;}
                        else if (Pb[i][j]/(Pb[i][j]+Pw[i][j])<R3) {a[i][j]=0.75; Nwx+=1;}
                        else     a[i][j]=0.5;
                    }
                    Tsum+=T[i][j];
                }
            }
            Tx[x]=Tsum/(n*n);                   //avg temp of the planet for every x
            Txsum+=Tx[x];
            Nbsum+=Nbx;                         //number of daisies sum up for all x
            Nwsum+=Nwx;
            temp<<Tx[x]-273.15<<endl;
            numb<<Nbx<<endl;
            numw<<Nwx<<endl;
        }
        Tavg[t]=Txsum/xx;                      //avg temp of the planet for time moment t
        Nb[t]=Nbsum/xx;                        //avg number of daisies for time moment t
        Nw[t]=Nwsum/xx;
    }
    temp.close();
    numb.close();
    numw.close();
    cout<<"Done!"<<endl;

    return 0;
}

//input values - luminosity, albedo, previous energies, current energy of neighbour cells
//output values - new energies, new temperatures

/*algorithm:
    1) importing luminosity data
    2) for every cell of the planet set the starting albedo to 0.5 (uninhabited) and calculate starting temp and energy
    3) for every given time/luminosity value we repeat the cycle 100 times to break down the time intervals
    4) inside every t,x for every cell [i,j] we:
                *use two matrices so we can have both new and old values at any time
                 mod 2 condition tells which matrix to use when
        4.1) add energy from all neighbouring cells
             & through if conditions make sure that the right cells are included (border cells!!)
        4.2) calculate the energy exchanged in diffusion
                the whole energy difference divided by number of neighbour cells times a coeff 1-e^(-dt/k)
        4.4) calculate the whole energy of a cell
                previous energy + energy absorbed - energy radiated + energy exchanged through diffusion
        4.5) find temp after diffusion
        4.6) determine if a daisy grows (see newnodff algorithm)
        4.7) new albedo, sum temp, count the number of black and white daisies
    5) sum up/find average temperatures && number of daisies
    6) store stuff
*/
