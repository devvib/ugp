#include <bits/stdc++.h>
#include "antoine_solver.h"
#include "linear_system_solver.h"
#include "energy_solver.h"
using namespace std;
double R, D, Cond; //Cond == L4
double P_total;
int no_of_comp = 3;
int no_of_stages;
int feed_stage;
double T_Feed;
// double F1, F2,F3;
vector<double> F(no_of_comp+1);
vector<vector<double>> H_param(no_of_comp+1,vector<double>(3));
vector<vector<double>> h_param(no_of_comp+1,vector<double>(3));
vector<double>A(no_of_comp+1),B(no_of_comp+1),C(no_of_comp+1);
  

int main()
{
  cout<<"Enter Reflux_ratio, Distillate, Condensate, Total_pressure, Temperature_of_Feed"<<endl;
  cin>>R>>D>>Cond>>P_total>>T_Feed;
  cout<<"Enter the Flow Rate of Each Component of Feed"<<endl;
  for(int i = 1; i<=no_of_comp; i++)
  {
    cin>>F[i];
  }
  cout<<"Enter the Number of Stages"<<endl;
  cin>>no_of_stages;

  cout<<"Enter the stage at which you're inserting the feed"<<endl;
  cin>>feed_stage;

  cout<<"Enter the 'H' parameters for each component in order of A, B, C"<<endl;
  for(size_t i = 1; i<=no_of_comp; i++)
  {
    double A1, B1, C1;
    cin>>A1>>B1>>C1;
    H_param[i][0] = A1;
    H_param[i][1] = B1;
    H_param[i][2] = C1;
  }

  cout<<"Enter the 'h' parameters for each component in order of A, B, C"<<endl;
  for(size_t i = 1; i<=no_of_comp; i++)
  {
    double A1, B1, C1;
    cin>>A1>>B1>>C1;
    h_param[i][0] = A1;
    h_param[i][1] = B1;
    h_param[i][2] = C1;
  }

  cout<<"Enter the Antoine parameters for each component in order of A, B and C"<<endl;
  for(size_t i = 1; i<=no_of_comp; i++)
  {
    double A1, B1, C1;
    cin>>A1>>B1>>C1;
    A[i] = A1;
    B[i] = B1;
    C[i] = C1;
  }


  vector<double> T(no_of_stages+1),V(no_of_stages+1),
  L(no_of_stages+1);
  vector<vector<double>>Sij,lij,xij;
  double qc, qr;

  // Assumed All the temperatures.
  // for(int i = 1; i<=no_of_stages; i++)
  // {
  //   T[i] = 60;
  // }
  T[1] = 60;
  T[2] = 65;
  T[3] = 70;
  T[4] = 80;

  // Assumed and Calculated V and L.
  V[2] = D+D*R;
  
  V[1] = 0;

  // Assumptions
  for(int i = 3; i<=no_of_stages; i++)
  {
    V[i] = V[2];

  }


  L[1] = D*R;
  L[no_of_stages] = Cond;

  for(int i = 2; i<no_of_stages-1; i++)
  {
    L[i] = V[i+1] + L[i-1] - V[i];
  }
  L[no_of_stages-1] = L[no_of_stages] + V[no_of_stages]; 

  // int ct=0;
  // int stage_no;cin>>stage_no;

 vector<vector<double>>x_dash;
 vector<vector<double>>y_dash;
//  cout<<"L: "<<endl;

   
//       for(auto it:L){
//         cout<<it<<" ";
//      }
//      cout<<endl;
//    cout<<endl<<"V "<<endl;
//       for(auto it:V){
//         cout<<it<<" ";
//      }
    //  cout<<endl;

 

  
  while (true)
  {
   //loop2 needs T,V,L

   
   L[1] = D*R;
  L[no_of_stages] = Cond;

  for(int i = 2; i<no_of_stages-1; i++)
  {
    L[i] = V[i+1] + L[i-1] - V[i];
  }
  L[no_of_stages-1] = L[no_of_stages] + V[no_of_stages]; 

    while (true)
    {
      
  
      // calculation of K (y=K*x);
      vector<vector<double>> K = calculate_K(A,B,C,T, P_total,no_of_stages);
      // cout<<endl<<"k "<<endl;
      // for(auto it:K){
      //   for(auto itt:it)cout<<itt<<" ";
      //   cout<<endl;
      // }


      // calcution of S (S[1][2]=k[1][2]*V[2]/L[2]);
      vector<vector<double>> S = calculate_S(K, V, L,D);
      //  cout<<endl<<"S "<<endl;
      // for(auto it:S){
      //   for(auto itt:it)cout<<itt<<" ";
      //   cout<<endl;
      // }


      // make matrix A,X and B for A*X=C
      vector<double> lf(no_of_comp+1);
      for(int i = 1; i<=no_of_comp; i++)
      {
        lf[i] = F[i];
      }
      vector<vector<double>>l=matrix_solver(S,lf,feed_stage);
      // cout<<"l "<<endl;
      // for(auto it:l){
      //   for(auto itt:it)cout<<itt<<" ";
      //   cout<<endl;
      // }
      // break;
      // we got X all l
      // now we calculate all xij;
      
      vector<vector<double>>x=calculate_x(l);
      x_dash=x;
      // cout<<"x.size "<<x[0].size()<<endl;

      //  Find the temperature using Newton's method
      vector<double> Tnew= temp_solver(T,x,A,B,C,P_total);
      // cout<<"temp: "<<endl;
      // for(auto it:Tnew)cout<<it<<' ';
      // cout<<endl;
      // break;
      //breaking condition for T

      if(accurate(Tnew,T)){Sij=S,lij=l; T=Tnew;xij=x;break;}
      T=Tnew;
      // break;
      // if(ct++>1)break;
    }
    // break;
    // Initializing and calculating vij.
    vector<vector<double> > vij = calculate_vij(Sij,lij);
    // Initializing and calculating Yij.
    vector<vector<double> > Yij = calculate_Yij(vij);
    y_dash=Yij;
    // Initializing and calculating Hij.
    vector<vector<double> > Hij = calculate_Hij(H_param,T);

    // Initializing and calculating hij.
    vector<vector<double> > hij = calculate_hij(h_param,T);

    // Initializing and calculating Hi.
    vector<double> Hi = calculate_Hi(Hij,Yij);

    // Initializing and calculating hi.
    vector<double> hi = calculate_hi(hij,xij);

    // Calculating hfi
    vector<double> hfi = calculatehfi(h_param, T_Feed);


     L[1] = D*R;
  L[no_of_stages] = Cond;

  for(int i = 2; i<no_of_stages-1; i++)
  {
    L[i] = V[i+1] + L[i-1] - V[i];
  }
  L[no_of_stages-1] = L[no_of_stages] + V[no_of_stages]; 
      
    
    // Assuming distillate as 'distillate'
    vector<double> Vnew = calculate_Vnew(Hi, hi, L, V, D, F, hfi,feed_stage);

    qc = V[2]*(Hi[2]-hi[1]);
    qr = (V[no_of_stages]*Hi[no_of_stages])+(L[no_of_stages]*hi[no_of_stages])-((V[no_of_stages]+L[no_of_stages])*hi[no_of_stages-1]);

    if(accurate(Vnew,V)){V = Vnew;break;}
    V = Vnew;
    // break;
  }
  //  L[2] = V[3]-D;
  // L[3] = V[4]+L[4];
  cout<<endl<<"RESULTS: "<<endl;
  cout<<"Temperature of columns : "<<endl;
  for(int i=1;i<=no_of_stages;i++)cout<<T[i]<<" ";
  cout<<endl;
  cout<<"Vapour molar flow rate : "<<endl;
  for(int i=1;i<=no_of_stages;i++)cout<<V[i]<<" ";
  cout<<endl;
  cout<<"Liquid molar flow rate : "<<endl;
  for(int i=1;i<=no_of_stages;i++)cout<<L[i]<<" ";
  cout<<endl;
  cout<<"heat out qc: "<<qc<<endl;
  cout<<"heat given qr: "<<qr<<endl;
  // for(auto it:x_dash[stage_no])cout<<it<<" ";
  // cout<<endl;
  
  for (int j = 1; j <= no_of_stages; j++)
  {
    
  cout<<"mole liq fractions at stage"<<no_of_stages-j+1<<": "<<endl;
  
  for (int i = 0; i <= no_of_comp; i++)
  {
    cout<<x_dash[i][j]<<" ";
  }
  cout<<endl;
  }
  cout<<endl;
  for (int j = 1; j <= no_of_stages; j++)
  {
    
  cout<<"mole vap fractions at stage"<<no_of_stages-j+1<<": "<<endl;
  
  for (int i = 0; i <= no_of_comp; i++)
  {
    cout<<y_dash[i][j]<<" ";
  }
  cout<<endl;
  }
  
  
  
  // for(auto it:x_dash)
  // {
  //   for(auto itt:it)cout<<itt<<" ";
  //   cout<<endl;
    
  // }

 
      
  return 0;
}