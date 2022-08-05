#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <fstream>
#include <algorithm>
#include <math.h>

#include "function.hpp"

using namespace std;
using namespace Eigen;

typedef Vector2d R2;
typedef Matrix<double,2,2> R22;
typedef VectorXi Tri;
typedef Vector4i Edge;

typedef Triplet<double> T;
typedef SparseMatrix<double> SpMat;
typedef VectorXd Vec;
typedef Matrix<double,6,6> R66;
typedef Matrix<double,6,3> R63;
typedef Matrix<double,3,3> R33;

double g(R2 pnt)
{
  if( (pnt[0] == 0) && (pnt[1] >= 1./2) && (pnt[1] <= 1) ){
    return -16*(pnt[1]-1)*(pnt[1]-1./2);
  }
  else{
    return 0.;
  }
}

double psi(const R2& pnt, int i)
{
  assert(i >= 0 && i < 3);
  double result;
  switch(i){
    case 0: result = 1-pnt[0]-pnt[1]; break;
    case 1: result = pnt[0]; break;
    case 2: result = pnt[1]; break;
    default: break;
  }
  return result;
}

R2 Gradpsi(int i)
{
  assert(i >= 0 && i < 3);
  R2 result;
  switch(i){
  case 0:
    result[0] = -1.; result[1] = -1.;
    break;
  case 1:
    result[0] = 1.; result[1] = 0;
    break;
  case 2:
    result[0] = 0; result[1] = 1.;
    break;
  default : break;
  }
  return result;
}

R2 F_T(const R2* T, const R2& u) // envoie un point du triangle de référence sur le triangle T.
{
  R22 B_T;
  B_T(0,0) = (T[1] - T[0])[0];
  B_T(1,0) = (T[1] - T[0])[1];
  B_T(0,1) = (T[2] - T[0])[0];
  B_T(1,1) = (T[2] - T[0])[1];

  return B_T*u+T[0];
}

R2 invF_T(const R2* T, const R2& u) // envoie un point du triangle T sur le triangle de référence.
{
  R22 B_T;
  B_T(0,0) = (T[1] - T[0])[0];
  B_T(1,0) = (T[1] - T[0])[1];
  B_T(0,1) = (T[2] - T[0])[0];
  B_T(1,1) = (T[2] - T[0])[1];

  R22 invB_T = B_T.inverse();

  return invB_T*(u-T[0]);
}
  
double phi(const R2& pnt, int i)
{
  assert(i >= 0 && i < 6);
  double result;
  switch(i){
    
  case 0: result = psi(pnt,i)*(2*psi(pnt,i)-1.); break;
  case 1: result = psi(pnt,i)*(2*psi(pnt,i)-1.); break;
  case 2: result = psi(pnt,i)*(2*psi(pnt,i)-1.); break;
  case 3: result = 4*psi(pnt,0)*psi(pnt,1); break;
  case 4: result = 4*psi(pnt,1)*psi(pnt,2); break;
  case 5: result = 4*psi(pnt,2)*psi(pnt,0); break;
  default: break;

  }
  return result;
}

R2 Gradphi(const R2& pnt, int i)
{
  assert(i >= 0 && i < 6);
  R2 result;
  switch(i){
  case 0:
    result = 4*psi(pnt,0)*Gradpsi(0) - Gradpsi(0);
    break;
  case 1:
    result = 4*psi(pnt,1)*Gradpsi(1) - Gradpsi(1);
    break;
  case 2:
    result = 4*psi(pnt,2)*Gradpsi(2) - Gradpsi(2);
    break;
  case 3:
      result = 4*(Gradpsi(0)*psi(pnt,1)+psi(pnt,0)*Gradpsi(1));
    break;
  case 4:
    result = 4*(Gradpsi(1)*psi(pnt,2)+psi(pnt,1)*Gradpsi(2));
    break;
  case 5:
    result = 4*(Gradpsi(2)*psi(pnt,0)+psi(pnt,2)*Gradpsi(0));
    break;
  default: break;
  }
  return result;
}

R66 MassEle(const R2* T)
{
  R22 B_T;
  B_T(0,0) = (T[1] - T[0])[0];
  B_T(1,0) = (T[1] - T[0])[1];
  B_T(0,1) = (T[2] - T[0])[0];
  B_T(1,1) = (T[2] - T[0])[1];

  R22 invB_T = B_T.inverse();

  double area_T = B_T.determinant()/2;

  R2 alp[7];
  alp[0][0] = 1./3; alp[0][1] = 1./3;
  alp[1][0] = (6-sqrt(15))/21; alp[1][1] = (6-sqrt(15))/21;
  alp[2][0] = (6-sqrt(15))/21; alp[2][1] = (9+2*sqrt(15))/21;
  alp[3][0] = (9+2*sqrt(15))/21; alp[3][1] = (6-sqrt(15))/21;
  alp[4][0] = (6+sqrt(15))/21; alp[4][1] = (6+sqrt(15))/21;
  alp[5][0] = (6+sqrt(15))/21; alp[5][1] = (9-2*sqrt(15))/21;
  alp[6][0] = (9-2*sqrt(15))/21; alp[6][1] = (6+sqrt(15))/21;

  double wei[7];
  wei[0] = 0.225; wei[1] = (155-sqrt(15))/1200; wei[2] = (155-sqrt(15))/1200; wei[3] = (155-sqrt(15))/1200; wei[4] = (155+sqrt(15))/1200; wei[5] = (155+sqrt(15))/1200; wei[6] = (155+sqrt(15))/1200; 

  R66 massele;

  for(int i = 0; i < 6; ++i){
    for(int j = 0; j < 6; ++j){
      massele(i,j) = 0;
    }
  }

  for(int i = 0; i < 6; ++i){
    for(int j = 0; j < 6; ++j){
      for(int k = 0; k < 7; ++k){
	massele(i,j) += area_T*wei[k]*phi(alp[k],i)*phi(alp[k],j);
      }
    }
  }

  return massele;
  
}

SpMat Mass(int Nt, int Nv, int Ne, map<int, Tri> Element, map<int, R2> Nodes)
{
  vector<T> tripletlist;
  for(int t = 0; t < Nt; ++t){
    R2 pnttria[3];
    for(int k = 0; k < 3; ++k){
      int temp = Element[t][k];
      pnttria[k] = Nodes[temp];
    }

    R66 massele = MassEle(pnttria);
    
    for(int i = 0; i < 6; ++i){
      for(int j = 0; j < 6; ++j){
	tripletlist.push_back(T(Element[t][i],Element[t][j],massele(i,j)));
      }
    }
  }

  SpMat result(Nv+Ne,Nv+Ne);
  result.setFromTriplets(tripletlist.begin(), tripletlist.end());
  return result;
}

R66 RigEle(const R2* T)
{
  R22 B_T;
  B_T(0,0) = (T[1] - T[0])[0];
  B_T(1,0) = (T[1] - T[0])[1];
  B_T(0,1) = (T[2] - T[0])[0];
  B_T(1,1) = (T[2] - T[0])[1];

  R22 invB_T = B_T.inverse();

  double area_T = B_T.determinant()/2;

  R2 alp[7];
  alp[0][0] = 1./3; alp[0][1] = 1./3;
  alp[1][0] = (6-sqrt(15))/21; alp[1][1] = (6-sqrt(15))/21;
  alp[2][0] = (6-sqrt(15))/21; alp[2][1] = (9+2*sqrt(15))/21;
  alp[3][0] = (9+2*sqrt(15))/21; alp[3][1] = (6-sqrt(15))/21;
  alp[4][0] = (6+sqrt(15))/21; alp[4][1] = (6+sqrt(15))/21;
  alp[5][0] = (6+sqrt(15))/21; alp[5][1] = (9-2*sqrt(15))/21;
  alp[6][0] = (9-2*sqrt(15))/21; alp[6][1] = (6+sqrt(15))/21;

  double wei[7];
  wei[0] = 0.225; wei[1] = (155-sqrt(15))/1200; wei[2] = (155-sqrt(15))/1200; wei[3] = (155-sqrt(15))/1200; wei[4] = (155+sqrt(15))/1200; wei[5] = (155+sqrt(15))/1200; wei[6] = (155+sqrt(15))/1200; 

  R66 rigele;
  
  for(int i = 0; i < 6; ++i){
    for(int j = 0; j < 6; ++j){
      rigele(i,j) = 0;
    }
  }

  for(int i = 0; i < 6; ++i){
    for(int j = 0; j < 6; ++j){
      for(int k = 0; k < 7; ++k){
	rigele(i,j) += area_T*wei[k]*( ((invB_T.adjoint())*Gradphi(alp[k],i)).adjoint() ) * ( (invB_T.adjoint())*Gradphi(alp[k],j) );
      }
    }
  }

  return rigele;
}

SpMat Rig(int Nt, int Nv, int Ne, map<int, Tri> Element, map<int, R2> Nodes)
{
  vector<T> tripletlist;
  for(int t = 0; t < Nt; ++t){
    R2 pnttria[3];
    for(int k = 0; k < 3; ++k){
      int temp = Element[t][k];
      pnttria[k] = Nodes[temp];
    }

    R66 rigele = RigEle(pnttria);
    
    for(int i = 0; i < 6; ++i){
      for(int j = 0; j < 6; ++j){
	tripletlist.push_back(T(Element[t][i],Element[t][j],rigele(i,j)));
      }
    }
  }

  SpMat result(Nv+Ne,Nv+Ne);
  result.setFromTriplets(tripletlist.begin(), tripletlist.end());
  return result;
}

R63 DivEle(const R2* T, int l)
{
  
  R22 B_T;
  B_T(0,0) = (T[1]-T[0])[0];
  B_T(1,0) = (T[1]-T[0])[1];
  B_T(0,1) = (T[2]-T[0])[0];
  B_T(1,1) = (T[2]-T[0])[1];

  R22 invB_T = B_T.inverse();

  double area_T = B_T.determinant()/2;
  
  R2 alp[7];
  alp[0][0] = 1./3; alp[0][1] = 1./3;
  alp[1][0] = (6-sqrt(15))/21; alp[1][1] = (6-sqrt(15))/21;
  alp[2][0] = (6-sqrt(15))/21; alp[2][1] = (9+2*sqrt(15))/21;
  alp[3][0] = (9+2*sqrt(15))/21; alp[3][1] = (6-sqrt(15))/21;
  alp[4][0] = (6+sqrt(15))/21; alp[4][1] = (6+sqrt(15))/21;
  alp[5][0] = (6+sqrt(15))/21; alp[5][1] = (9-2*sqrt(15))/21;
  alp[6][0] = (9-2*sqrt(15))/21; alp[6][1] = (6+sqrt(15))/21;

  double wei[7];
  wei[0] = 0.225; wei[1] = (155-sqrt(15))/1200; wei[2] = (155-sqrt(15))/1200; wei[3] = (155-sqrt(15))/1200; wei[4] = (155+sqrt(15))/1200; wei[5] = (155+sqrt(15))/1200; wei[6] = (155+sqrt(15))/1200;
  
  R63 divele;
  for(int i = 0; i < 6; ++i){
    for(int j = 0; j < 3; ++j){
      divele(i,j) = 0;
    }
  }
  
  for(int i = 0; i < 6; ++i){
    for(int j = 0; j < 3; ++j){
      for(int k = 0; k < 7; ++k){
	divele(i,j) += area_T*wei[k]*psi(alp[k],j)*( (invB_T.adjoint())*Gradphi(alp[k],i))[l];
      }
    }
  }

  return divele;

}

SpMat Div(int Nt, int Nv, int Ne, map<int, Tri> Element, map<int, R2> Nodes, int l)
{
  vector<T> tripletlist;
  for(int t = 0; t < Nt; ++t){
    R2 pnttria[3];
    for(int k = 0; k < 3; ++k){
      int temp = Element[t][k];
      pnttria[k] = Nodes[temp];
    }
    
    R63 divele = DivEle(pnttria,l);

    for(int i = 0; i < 6; ++i){
      for(int j = 0; j < 3; ++j){
	tripletlist.push_back(T(Element[t][i],Element[t][j],divele(i,j)));
      }
    }
  }

  SpMat result(Nv+Ne,Nv);
  result.setFromTriplets(tripletlist.begin(), tripletlist.end());
  return result;
}

R33 PenEle(const R2* T)
{
  R2 edge0 = T[0]-T[1];
  R2 edge1 = T[0]-T[2];

  double areaT = (1./2)*(edge0[0]*edge1[1]-edge0[1]*edge1[0]);

  R33 penele;
  for(int i = 0; i < 3; ++i){
    for(int j = 0; j < 3; ++j){
      penele(i,j) = 0;
      if(i==j){
	penele(i,j) = areaT/6;
      }
      else{
	penele(i,j) = areaT/12;
      }
    }
  }
  return penele;
}

SpMat Pen(int Nt, int Nv, int Ne, map<int, Tri> Element, map<int, R2> Nodes, int l)
{
  vector<T> tripletlist;
  for(int t = 0; t < Nt; ++t){
    R2 pnttria[3];
    for(int k = 0; k < 3; ++k){
      int temp = Element[t][k];
      pnttria[k] = Nodes[temp];
    }
        
    R33 penele = PenEle(pnttria);

    for(int i = 0; i < 3; ++i){
      for(int j = 0; j < 3; ++j){
	tripletlist.push_back(T(Element[t][i],Element[t][j],penele(i,j)));
      }
    }
  }

  SpMat result(Nv,Nv);
  result.setFromTriplets(tripletlist.begin(), tripletlist.end());
  return result;
}

SpMat AssemSys(mesh& Th, double nu, double alp, double eps, double tgv, Vec& v)
{
  int Nv = Th.Nv();
  int Nt = Th.Nt();
  int Ne = Th.Ne();
  int Neb = Th.Neb();
  int Nlbl = Th.Nlbl();

  map<int, R2> Nodes = Th.Nodes();

  map<int, Tri> Element = Th.Element();

  map<int, Edge> Boundary = Th.Boundary();

  map<int, vector<int> > lbl_nodes = Th.Lbl_Nodes(); // Lbl starts at 1

  map<int, vector<int> > lbl_element = Th.Lbl_Element(); // Lbl starts at 0

  map<int, vector<int> > lbl_boundary = Th.Lbl_Boundary(); // Lbl starts at 1

  map<int, Edge> edge = Th.Edg();
  
  map<int, vector<Edge> > lbl_edge = Th.Lbl_Edge();

  edge = Th.GetEdge();
  
  // for(int i = 0; i < Ne; ++i){
  //   cout << edge[i] << endl << endl;
  // }
  
  map<int, int> lbl_size = Th.Lbl_Size();

  // for(int i = 1; i < 5; ++i){
  //   for(int j = 0; j < lbl_size[i]; ++j){
  //     cout << lbl_edge[i][j] << endl << endl;
  //   }
  // }
  
  Nodes = Th.Nodes();

  Element = Th.Element();
  
  vector<T> tripletlist;
  for(int t = 0; t < Nt; ++t){
    R2 pnttria[3];
    for(int k = 0; k < 3; ++k){
      int temp = Element[t][k];
      pnttria[k] = Nodes[temp];
    }

    R66 massele = MassEle(pnttria);
    R66 rigele = RigEle(pnttria);
    R33 penele = PenEle(pnttria);
    R63 divele1 = DivEle(pnttria,0);
    R63 divele2 = DivEle(pnttria,1);
    
    for(int i = 0; i < 6; ++i){
      for(int j = 0; j < 6; ++j){
        tripletlist.push_back( T(Element[t][i],Element[t][j],nu*rigele(i,j)) );
	tripletlist.push_back( T(Element[t][i] + (Nv+Ne),Element[t][j] + (Nv+Ne),nu*rigele(i,j)) );
	tripletlist.push_back( T(Element[t][i],Element[t][j],alp*massele(i,j)) );
	tripletlist.push_back( T(Element[t][i] + (Nv+Ne),Element[t][j] + (Nv+Ne),alp*massele(i,j)) );
	
      }
    }
    
    for(int i = 0; i < 6; ++i){
      for(int j = 0; j < 3; ++j){
    	tripletlist.push_back( T(2*(Nv+Ne)+Element[t][j],Element[t][i],-divele1(i,j)) );
    	tripletlist.push_back( T(Element[t][i],2*(Nv+Ne)+Element[t][j],-divele1(i,j)) );
    	tripletlist.push_back( T((Nv+Ne)+Element[t][i],2*(Nv+Ne)+Element[t][j],-divele2(i,j)) );
    	tripletlist.push_back( T(2*(Nv+Ne)+Element[t][j],(Nv+Ne)+Element[t][i],-divele2(i,j)) );
      }
    }
    
    for(int i = 0; i < 3; ++i){
      for(int j = 0; j < 3; ++j){
	tripletlist.push_back( T(2*(Nv+Ne)+Element[t][i],2*(Nv+Ne)+Element[t][j],-eps*penele(i,j)) );
      }
    }
  }
  
  SpMat matglob(3*Nv+2*Ne,3*Nv+2*Ne);
  matglob.setFromTriplets(tripletlist.begin(), tripletlist.end());
  
  v.resize(2*(Nv+Ne)+Nv);
  
  for(int i = 1; i < Nlbl+1; ++i){

    // if(i == 3){
    //   break;
    // }
    
    for(int j = 0; j < lbl_size[i]; ++j){
      
      matglob.coeffRef(lbl_edge[i][j][0],lbl_edge[i][j][0]) = tgv;
      matglob.coeffRef(lbl_edge[i][j][1],lbl_edge[i][j][1]) = tgv;
      matglob.coeffRef(lbl_edge[i][j][2],lbl_edge[i][j][2]) = tgv;
      matglob.coeffRef(lbl_edge[i][j][0]+Nv+Ne,lbl_edge[i][j][0]+Nv+Ne) = tgv;
      matglob.coeffRef(lbl_edge[i][j][1]+Nv+Ne,lbl_edge[i][j][1]+Nv+Ne) = tgv;
      matglob.coeffRef(lbl_edge[i][j][2]+Nv+Ne,lbl_edge[i][j][2]+Nv+Ne) = tgv;

      v[lbl_edge[i][j][0]] = 0;
      v[lbl_edge[i][j][1]] = 0;
      v[lbl_edge[i][j][2]] = 0;
      v[lbl_edge[i][j][0]+Nv+Ne] = 0;
      v[lbl_edge[i][j][1]+Nv+Ne] = 0;
      v[lbl_edge[i][j][2]+Nv+Ne] = 0;
      
    }
  }

  for(int j = 0; j < lbl_size[1]; ++j){
    v[lbl_edge[3][j][0]] = tgv;//*g(Nodes[lbl_edge[1][j][0]]);
    v[lbl_edge[3][j][1]] = tgv;//*g(Nodes[lbl_edge[1][j][1]]);
    v[lbl_edge[3][j][2]] = tgv;//*g(Nodes[lbl_edge[1][j][2]]);
  }
  
  return matglob;
}

Vec AssemSnd(const mesh& Th, const Vec& u, double dt)
{
  
  int Nv = Th.Nv();
  int Nt = Th.Nt();
  int Ne = Th.Ne();

  map<int, R2> Nodes = Th.Nodes();
  map<int, Tri> Element = Th.Element();

  Vec result(2*(Nv+Ne)+Nv);

  for(int w = 0; w < 2*(Nv+Ne)+Nv; ++w){result[w] = 0;}
  
  R2 alp[7];
  alp[0][0] = 1./3; alp[0][1] = 1./3;
  alp[1][0] = (6-sqrt(15))/21; alp[1][1] = (6-sqrt(15))/21;
  alp[2][0] = (6-sqrt(15))/21; alp[2][1] = (9+2*sqrt(15))/21;
  alp[3][0] = (9+2*sqrt(15))/21; alp[3][1] = (6-sqrt(15))/21;
  alp[4][0] = (6+sqrt(15))/21; alp[4][1] = (6+sqrt(15))/21;
  alp[5][0] = (6+sqrt(15))/21; alp[5][1] = (9-2*sqrt(15))/21;
  alp[6][0] = (9-2*sqrt(15))/21; alp[6][1] = (6+sqrt(15))/21;
    
  double wei[7];
  wei[0] = 0.225; wei[1] = (155-sqrt(15))/1200; wei[2] = (155-sqrt(15))/1200; wei[3] = (155-sqrt(15))/1200; wei[4] = (155+sqrt(15))/1200; wei[5] = (155+sqrt(15))/1200; wei[6] = (155+sqrt(15))/1200;

  
  for(int t = 0; t < Nt; ++t){

    Tri T = Element[t];
    R2 pnttria[6];
    for(int k = 0; k < 6; ++k){
      pnttria[k] = Nodes[T[k]];
    }
    
    double area_T = ((pnttria[1]-pnttria[0])[0]*(pnttria[2]-pnttria[0])[1] - (pnttria[1]-pnttria[0])[1]*(pnttria[2]-pnttria[0])[0]);
    
    R2 pos_T[7];
    R2 new_pos_T[7];
        
    for(int m = 0; m < 7; ++m){
      
      pos_T[m] = F_T(pnttria,alp[m]);
      new_pos_T[m] = pos_T[m];

      for(int i = 0; i < 6; ++i){

	R2 temp;	
	temp[0] = u[T[i]];
	temp[1] = u[T[i]+Nv+Ne];
	new_pos_T[m] -= dt*temp*phi(alp[m],i);

      }
    }
    
    int new_tri[7];
    
    for(int f = 0; f < 7; ++f){
      new_tri[f] = 0;
    }

    for(int m = 0; m < 7; ++m){

      // int k = 0;
      int current_pnt = T[0];
      int closest_pnt;
      int next_pnt = current_pnt;
      double min_ele = numeric_limits<double>::infinity();
      // double prec;
      map<int,double> dist_to_pnt;//[3*(Th.NbofNeighbors())[current_pnt]];

      dist_to_pnt[current_pnt] = (Nodes[current_pnt]-new_pos_T[m]).norm();//abs((Nodes[current_pnt]-new_pos_T[m]).adjoint()*(Nodes[current_pnt]-new_pos_T[m]));
      dist_to_pnt[next_pnt] = dist_to_pnt[current_pnt];
      
      bool found = false;

      while(found == false){	

	// for(int h = 0; h < 3*(Th.NbofNeighbors())[current_pnt]; ++h){
	//   dist_to_pnt[h] = numeric_limits<double>::infinity();
	// }
      
	for(int n = 0; n < (Th.NbofNeighbors())[current_pnt]; ++n){
	
	  Tri Neighbor = Element[(Th.Neighbors())[current_pnt][n]];
	  R2 pntneighbor[6];
	  pntneighbor[0] = Nodes[Neighbor[0]];
	  pntneighbor[1] = Nodes[Neighbor[1]];
	  pntneighbor[2] = Nodes[Neighbor[2]];
	  pntneighbor[3] = Nodes[Neighbor[3]];
	  pntneighbor[4] = Nodes[Neighbor[4]];
	  pntneighbor[5] = Nodes[Neighbor[5]];
	  
	  R2 s0 = pntneighbor[0]-new_pos_T[m];
	  R2 s1 = pntneighbor[1]-new_pos_T[m];
	  R2 s2 = pntneighbor[2]-new_pos_T[m];
	  
	  double m0 = (s0.adjoint())*s0;
	  double m1 = (s1.adjoint())*s1;
	  double m2 = (s2.adjoint())*s2;
	  
	  double mini[] = {m0,m1,m2};
	  
	  // prec = min_ele;
	  
	  min_ele = *min_element(mini,mini+3);
	  
	  int y = -1;
	  
	  for(int i = 0; i < 3; ++i){
	    if (min_ele == mini[i]){
	      y = i;
	    }
	  }

	  if( mini[y] <= dist_to_pnt[next_pnt] ){
	    next_pnt = Neighbor[y];
	    dist_to_pnt[next_pnt] = mini[y];
	  }
	  
	  // dist_to_pnt[Neighbor[y]] = mini[y];
  
	  R2 u = invF_T(pntneighbor,new_pos_T[m]);
	  	  
	  if ( (0 <= u[0]) && (u[0] <= 1) && (0 <= u[1]) && (u[1] <= 1) ) {
	    found = true;

	    // cout << "We have found our triangle !! it's the triangle nb :" << endl;
	    // cout << (Th.Neighbors())[current_pnt][n] << endl;
	    
	    // cout << (Th.NbofNeighbors())[current_pnt] << endl;
	    // cout << current_pnt << endl;
	    // cout << n << endl;
	    // cout << (Th.Neighbors())[current_pnt][n] << endl << endl;

	    new_tri[m] = (Th.Neighbors())[current_pnt][n];

	    // int trash;
	    // cin >> trash;

	    
	    break;
	  }

	  // if (min_ele <= prec){
	  //   k = current_pnt;
	  //   current_pnt = Neighbor[y];
	  // }

	}
	
	if(found == true){
	  break;
	}

	
	
	// double dist_to_pnt_ele = *min_element(dist_to_pnt, dist_to_pnt+3*(Th.NbofNeighbors())[current_pnt]);

	// int c = -1;
	
	// for(int i = 0; i < (Th.NbofNeighbors())[current_pnt]; ++i){
	//   if (dist_to_pnt_ele == dist_to_pnt[i]){
	//     c = i;
	//   }
	// }
	
	// k = c;
	
	if(next_pnt != current_pnt){
	  current_pnt = next_pnt;
	}
	else{
	  // cout << "le point recherché est en dehors du maillaige !" << endl;
	  // cout << new_pos_T[m] << endl << endl;

	  R2 border_point = Nodes[current_pnt];  

	  R2 proj_pot[3*(Th.NbofNeighbors())[current_pnt]];

	  for(int y = 0; y < (Th.NbofNeighbors())[current_pnt]; ++y){
	    Tri A = Element[(Th.Neighbors())[current_pnt][y]];

	    for(int e = 0; e < 3; ++e){
	      double ps;
	      if(Nodes[A[e]] != border_point){
		ps = abs(((Nodes[A[e]]-border_point).adjoint())*(new_pos_T[m]-border_point));
		proj_pot[3*y+e] = border_point + ps*(Nodes[A[e]] - border_point)/sqrt((Nodes[A[e]] - border_point).adjoint()*(Nodes[A[e]] - border_point));
	      }
	      else{
		proj_pot[3*y+e] = border_point;
	      }
	    }
	    
	  }

	  double min_proj[3*(Th.NbofNeighbors())[current_pnt]];

	  for(int e = 0; e < 3*(Th.NbofNeighbors())[current_pnt]; ++e){
	    min_proj[e] = ((proj_pot[e]-new_pos_T[m]).adjoint())*(proj_pot[e]-new_pos_T[m]);
	  }

	  double minim = *min_element(min_proj,min_proj+3*(Th.NbofNeighbors())[current_pnt]);

	  int z = -1;
	  
	  for(int s = 0; s < 3*(Th.NbofNeighbors())[current_pnt]; ++s){
	    if(min_proj[s] == minim){
	      z = s;
	    }
	  }

	  int q;
	  q = z/3; // ((Th.NbofNeighbors())[current_pnt]);
	  
	  // r = z%(3*(Th.NbofNeighbors())[current_pnt]);
	  
	  new_pos_T[m] = proj_pot[z];
	  new_tri[m] = (Th.Neighbors())[current_pnt][q];

	  // cout << new_pos_T[m] << endl << endl;
	  // cout << new_tri[m] << endl << endl;

	  // cin >> trash;

	  found = true;
	  
	}
	


      }

	    
	
    }

    // for(int i = 0; i < 7; ++i){
    //   cout << new_pos_T[i] << endl << endl;
    //   cout << new_tri[i] << endl;
    // }

    // for(int f = 0; f < 7; ++f){
    //   cout << new_tri[f] << endl;
    // }
    // cout << endl;

    for(int i = 0; i < 6; ++i){
      for(int k = 0; k < 7; ++k){
	for(int j = 0; j < 6; ++j){
	  // R2 pnt[6] = {Nodes[Element[new_tri[k]][0]], Nodes[Element[new_tri[k]][1]], Nodes[Element[new_tri[k]][2]], Nodes[Element[new_tri[k]][3]], Nodes[Element[new_tri[k]][4]], Nodes[Element[new_tri[k]][5]]};
	  
	  result[T[i]] += (1/dt)*area_T*wei[k]*u[Element[new_tri[k]][j]]*phi(alp[k],j)*phi(alp[k],i);
	  result[T[i]+Nv+Ne] += (1/dt)*area_T*wei[k]*u[Element[new_tri[k]][j]+Nv+Ne]*phi(alp[k],j)*phi(alp[k],i);
	}
      }
    }
    
  }
  
  return result;

}

void SndBorderCond(mesh& Th, Vec& v, double tgv)
{
  int Nlbl = Th.Nlbl();
  int Nv = Th.Nv();
  int Ne = Th.Ne();
  map<int, R2> Nodes = Th.Nodes();
  map<int, int> lbl_size = Th.Lbl_Size();
  map<int, vector<Edge> > lbl_edge = Th.Lbl_Edge();
  
  for(int i = 1; i < Nlbl+1; ++i){

    // if (i == 3){
    //   break;
    // }
    
    for(int j = 0; j < lbl_size[i]; ++j){

      v[lbl_edge[i][j][0]] = 0;
      v[lbl_edge[i][j][1]] = 0;
      v[lbl_edge[i][j][2]] = 0;
      v[lbl_edge[i][j][0]+Nv+Ne] = 0;
      v[lbl_edge[i][j][1]+Nv+Ne] = 0;
      v[lbl_edge[i][j][2]+Nv+Ne] = 0;
      
    }
  }

  for(int j = 0; j < lbl_size[1]; ++j){
    v[lbl_edge[3][j][0]] = tgv;//*g(Nodes[lbl_edge[1][j][0]]);
    v[lbl_edge[3][j][1]] = tgv;//*g(Nodes[lbl_edge[1][j][1]]);
    v[lbl_edge[3][j][2]] = tgv;//*g(Nodes[lbl_edge[1][j][2]]);
  }

}

void SavetoPlot(const mesh Th, const Vec& x, const char* filename)
{
  ofstream f(filename, ofstream::out);

  int Nbt = Th.Nt();
  int Nbv = Th.Nv();
  int Nbe = Th.Ne();

  map<int, Tri> Element = Th.Element();

  for(int t = 0; t < Nbt; ++t){
    
    for(int j = 0; j < 6; ++j){
      if(j < 3){
	f << x[Element[t][j]] << " ";
      }
      else{
	f << x[Element[t][(j+1)%3+3]] << " ";
      }
    }
    
    for(int j = 0; j < 6; ++j){
      if(j < 3){
	f << x[Element[t][j]+Nbv+Nbe] << " ";
      }
      else{
	f << x[Element[t][(j+1)%3+3]+Nbv+Nbe] << " ";
      }
    }
        
    for(int j = 0; j < 3; ++j){
      f << x[Element[t][j]+2*(Nbv+Nbe)] << " ";
    }
    
    f << endl;
    
  }
}


int StructMatGlob(const SpMat& M, bool coef)
{
  int result = 0;
  ofstream f("matrice.dat");
  for(int i = 0; i < M.rows(); ++i){
    for(int j = 0; j < M.cols(); ++j){
	if(abs(M.coeff(i,j)) > pow(10,-12)){
	  result += 1;
	  f << i << "\t" << j << "\t" << M.coeff(i,j) << endl;
	  if(coef == false){
	    cout << "x";
	  }
	  else{
	    cout << M.coeff(i,j) << " ";
	  }
	}
	else{
	  if(coef == false){
	    cout << "o";
	  }
	  else{
	    cout << 0 << " ";
	  }
	}
    }
    cout << endl;
    f << endl;
  }
  cout << result << endl;
  return result;
}
