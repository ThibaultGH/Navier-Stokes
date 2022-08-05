#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

#include "function.hpp"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <assert.h>

using namespace std;
using namespace Eigen;

typedef Eigen::Triplet<double> T;
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Matrix<double,6,6> RigMat;
typedef Eigen::Matrix<double,6,3> DivMat;
typedef Eigen::Matrix<double,3,3> PenMat;
typedef Eigen::VectorXd VEigen;

typedef vector<double> vecd;
typedef vector<int> veci;

void Dim(int& Nodes, int& Element, int& Edge, int& Boundary, const char* filename)
{
  ifstream myflux(filename, ios::in);

  if(myflux){
    
    string word;

    myflux >> word;
    Nodes = atoi(word.c_str());
    myflux >> word;
    Element = atoi(word.c_str());
    myflux >> word;
    Boundary = atoi(word.c_str());
  }
  Edge = (3*Element-Boundary)/2 + Boundary;
}

void Load(VEigen* Nodes, veci* Element, veci* Boundary, const char* filename)
{
  ifstream myflux(filename, ios::in);

  if(myflux){

    string word;
    int NumNodes;
    int NumElement;
    int NumBoundary;

    myflux >> word;
    NumNodes = atoi(word.c_str());
    myflux >> word;
    NumElement = atoi(word.c_str());
    myflux >> word;
    NumBoundary = atoi(word.c_str());
    
    
    for(int i = 0; i < NumNodes; ++i){
      Nodes[i].resize(2);
      myflux >> word;
      Nodes[i][0] = atof(word.c_str());
      myflux >> word;
      Nodes[i][1] = atof(word.c_str());
      myflux >> word;
    }
    
    for(int i = 0; i < NumElement; ++i){
      Element[i].resize(3);
      myflux >> word;
      Element[i][0] = atoi(word.c_str())-1;
      myflux >> word;
      Element[i][1] = atoi(word.c_str())-1;
      myflux >> word;
      Element[i][2] = atoi(word.c_str())-1;
      myflux >> word;
    }

    for(int i = 0; i < NumBoundary; ++i){
      Boundary[i].resize(3);
      myflux >> word;
      Boundary[i][0] = atoi(word.c_str())-1;
      myflux >> word;
      Boundary[i][1] = atoi(word.c_str())-1;
      myflux >> word;
      Boundary[i][2] = atoi(word.c_str());
    }
  }
}

void Trian(int Nbv, int Nbt, veci* Element, veci* trian)
{
  for(int i = 0; i < Nbv; ++i){trian[i].push_back(0);}
  for(int t = 0; t < Nbt; ++t){
    for(int i = 0; i < 3; ++i){
      trian[Element[t][i]].push_back(t);
      trian[Element[t][i]][0] += 1;
    }
  }
}

/*
void Edge(int Nbv, int Nbt, int Nbe, VEigen* Nodes, veci* Element, veci* edge)
{
  








}
*/

void Edge(int Nbv, int Nbt, int Nbe, veci* Element, veci* edge, veci* Numedge)
{
  veci first(Nbv,-1);
  veci next;
  
  int v0,v1;
  int I;
  int found;
  int Na = 0;
  
  for(int T = 0; T < Nbt; ++T){
    for(int i = 0; i < 3; ++i){
      
      v0 = min(Element[T][i], Element[T][(i+1)%3]);
      v1 = max(Element[T][i], Element[T][(i+1)%3]);
      
      Numedge[T].resize(3);
      Numedge[T][i] = Na;
      
      I = first[v0];
      found = 0;
      
      while( (found == 0) && (I != -1) ){
	
	if ( (v0 == edge[I][0]) && (v1 == edge[I][1]) ){
	  found = 1;
	  Numedge[T][i] = I;
	  edge[I][2] = 0;
	}
	I = next[I];
      }
      
      if(found == 0){

	edge[Na].resize(3);
	
	edge[Na][0] = v0;
	edge[Na][1] = v1;
	edge[Na][2] = 1;

	next.push_back(first[v0]);
	
	first[v0] = Na;
	
	Na += 1;
      }
    }
  } 
}

veci IdNodesBound(int Nbv, int Nbe, int Nbeb, veci* edge, veci* Boundary)
{
  veci idnodesbound(2*Nbeb);
  veci doublon(Nbv,0);
  //for(int i = 0; i < Nbv; ++i){doublon[i] = 0;}
  int k = 0;
  
  for(int i = 0; i < Nbeb; ++i){
    for(int j = 0; j < 2; ++j){
      if(doublon[Boundary[i][j]] == 0){
	idnodesbound[k] = Boundary[i][j];
	k += 1;
	doublon[Boundary[i][j]] = 1;
      }
    }
  }
  
  for(int i = 0; i < Nbe; ++i){
    if(edge[i][2] == 1){
      idnodesbound[k] = Nbv+i;
      k += 1;
    }
  }
  return idnodesbound;
}

veci IdNodesLbl(veci& idnodesbound, VEigen* MidPnt, int lbl)
{
  veci result;
  int n = idnodesbound.size();
  switch(lbl){
  case 0:
    for(int i = 0; i < n; ++i){
      if(abs(MidPnt[idnodesbound[i]][1]-1) < pow(10,-8)){
	result.push_back(idnodesbound[i]);
      }
    }
  case 1:
    for(int i = 0; i < n; ++i){
      if(abs(3*MidPnt[idnodesbound[i]][1]+MidPnt[idnodesbound[i]][0]-8.5) < pow(10,-8)){
	result.push_back(idnodesbound[i]);
      }
    }
  case 2:
    for(int i = 0; i < n; ++i){
      if(abs(MidPnt[idnodesbound[i]][1]-20) < pow(10,-8)){
	result.push_back(idnodesbound[i]);
      }
    }
  default: break;
  }
  return result;
}
 

void ndofP2(int& Nbv, int& Nbe, VEigen* Nodes, veci* edge, VEigen* MidPnt)
{  
  for(int i = 0; i < Nbv+Nbe; ++i){
    if(i<Nbv){
      MidPnt[i].resize(2);
      MidPnt[i] = Nodes[i];
    }
    else{
      MidPnt[i].resize(2);
      MidPnt[i] = 1./2*(Nodes[edge[i-Nbv][0]]+Nodes[edge[i-Nbv][1]]);
    }
  }
}

void dofP2(int Nbt, int Nbv, veci* Element, veci* Numedge, veci* ddl)
{
  for(int i = 0; i < Nbt; ++i){
    ddl[i].resize(6);
    for(int j = 0; j < 6; ++j){
      if(j<3){
	ddl[i][j] = Element[i][j];
      }
      else{
	ddl[i][j] = Nbv + Numedge[i][j-3];
      }
    }
  }
}

bool IsSymetric(const SpMat& M)
{
  bool result = true;
  int m = M.rows();
  int n = M.cols();
  for(int i = 0; i < m; ++i){
    for(int j = 0; j < n; ++j){
      if(abs(M.coeff(i,j) - M.coeff(j,i)) > pow(10,-8)){
	result = false;
      }
    }
  }
  return result;
}

double g(VEigen pnt)
{
  if( (pnt[0] == 0) && (pnt[1] >= 1./2) && (pnt[1] <= 1) ){
    return -16*(pnt[1]-1)*(pnt[1]-1./2);
  }
  else{
    return 0.;
  }
}

VEigen Gradg(VEigen pnt)
{
  VEigen result(2);
  result(0) = 0;
  result(1) = -32*pnt[1]+24;
  return result;
}

double lambda(const VEigen& pnt, int i)
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

VEigen Gradlambda(int i)
{
  assert(i >= 0 && i < 3);
  VEigen result(2);
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

double divlambda(int i)
{
  assert( (i >= 0) && ( i < 3 ) );
  double result;
  switch(i){
  case 0: result = -2; break;
  case 1: result = 1; break;
  case 2: result = 1; break;
  default: break;
  }
  return result;
}

double phi(const VEigen& pnt, int i)
{
  assert(i >= 0 && i < 6);
  double result;
  switch(i){
  case 0: result = lambda(pnt,i)*(2*lambda(pnt,i)-1.); break;
  case 1: result = lambda(pnt,i)*(2*lambda(pnt,i)-1.); break;
  case 2: result = lambda(pnt,i)*(2*lambda(pnt,i)-1.); break;
  case 3: result = 4*lambda(pnt,0)*lambda(pnt,1); break;
  case 4: result = 4*lambda(pnt,1)*lambda(pnt,2); break;
  case 5: result = 4*lambda(pnt,2)*lambda(pnt,0); break;
  default: break;
  }
  return result;
}

VEigen Gradphi(const VEigen& pnt, int i)
{
  assert(i >= 0 && i < 6);
  VEigen result(2);
  switch(i){
  case 0:
    result = 4*lambda(pnt,0)*Gradlambda(0) - Gradlambda(0);
    break;
  case 1:
    result = 4*lambda(pnt,1)*Gradlambda(1) - Gradlambda(1);
    break;
  case 2:
    result = 4*lambda(pnt,2)*Gradlambda(2) - Gradlambda(2);
    break;
  case 3:
      result = 4*(Gradlambda(0)*lambda(pnt,1)+lambda(pnt,0)*Gradlambda(1));
    break;
  case 4:
    result = 4*(Gradlambda(1)*lambda(pnt,2)+lambda(pnt,1)*Gradlambda(2));
    break;
  case 5:
    result = 4*(Gradlambda(2)*lambda(pnt,0)+lambda(pnt,2)*Gradlambda(0));
    break;
  default: break;
  }
  return result;
}

double dxdyphi(const VEigen& pnt, int i, int j)
{
  assert( (i >= 0) && (i < 6) );
  double result;
  if(j==1){
    switch(i){
    case 0: result = 4*lambda(pnt,0)*Gradlambda(0)[0]-Gradlambda(0)[0]; break;
    case 1: result = 4*lambda(pnt,1)*Gradlambda(1)[0]-Gradlambda(1)[0]; break;
    case 2: result = 4*lambda(pnt,2)*Gradlambda(2)[0]-Gradlambda(2)[0]; break;
    case 3: result = 4*(Gradlambda(0)[0]*lambda(pnt,1)+lambda(pnt,0)*Gradlambda(0)[0]); break;
    case 4: result = 4*(Gradlambda(1)[0]*lambda(pnt,2)+lambda(pnt,1)*Gradlambda(2)[0]); break;
    case 5: result = 4*(Gradlambda(2)[0]*lambda(pnt,0)+lambda(pnt,2)*Gradlambda(0)[0]); break;
    default: break;
    }
  }
  else{
    switch(i){
    case 0: result = 4*lambda(pnt,0)*Gradlambda(0)[1]-Gradlambda(0)[1]; break;
    case 1: result = 4*lambda(pnt,1)*Gradlambda(1)[1]-Gradlambda(1)[1]; break;
    case 2: result = 4*lambda(pnt,2)*Gradlambda(2)[1]-Gradlambda(2)[1]; break;
    case 3: result = 4*(Gradlambda(0)[1]*lambda(pnt,1)+lambda(pnt,0)*Gradlambda(0)[1]); break;
    case 4: result = 4*(Gradlambda(1)[1]*lambda(pnt,2)+lambda(pnt,1)*Gradlambda(2)[1]); break;
    case 5: result = 4*(Gradlambda(2)[1]*lambda(pnt,0)+lambda(pnt,2)*Gradlambda(0)[1]); break;
    default: break;
    }
  }
  return result;
}

RigMat MassEle(const VEigen* Nodes)
{
  Matrix<double,2,2> B_T;
  B_T(0,0) = (Nodes[1]-Nodes[0])[0];
  B_T(1,0) = (Nodes[1]-Nodes[0])[1];
  B_T(0,1) = (Nodes[2]-Nodes[0])[0];
  B_T(1,1) = (Nodes[2]-Nodes[0])[1];

  double areaT = B_T.determinant()/2;

  Vector2d alp0(1./3,1./3);
  Vector2d alp1((6-sqrt(15))/21,(6-sqrt(15))/21);
  Vector2d alp2((6-sqrt(15))/21,(9+2*sqrt(15))/21);
  Vector2d alp3((9+2*sqrt(15))/21,(6-sqrt(15))/21);
  Vector2d alp4((6+sqrt(15))/21,(6+sqrt(15))/21);
  Vector2d alp5((6+sqrt(15))/21,(9-2*sqrt(15))/21);
  Vector2d alp6((9-2*sqrt(15))/21,(6+sqrt(15))/21);

  double weight0 = 0.225;
  double weight1 = (155-sqrt(15))/1200;
  double weight2 = (155+sqrt(15))/1200;
  
  RigMat massele;
  for(int i = 0; i < 6; ++i){
    for(int j = 0; j < 6; ++j){
      massele(i,j) = weight0*areaT*phi(alp0,i)*phi(alp0,j) + weight1*areaT*phi(alp1,i)*phi(alp1,j) + weight1*areaT*phi(alp2,i)*phi(alp2,j) + weight1*areaT*phi(alp3,i)*phi(alp3,j) + weight2*areaT*phi(alp4,i)*phi(alp4,j) + weight2*areaT*phi(alp5,i)*phi(alp5,j) + weight2*areaT*phi(alp6,i)*phi(alp6,j);
    }
  }

  return massele;
}

SpMat Mass(VEigen* MidPnt, veci* ddl, veci* Element, int Nbv, int Nbt, int Nbe)
{
  vector<T> tripletlist;
  for(int t = 0; t < Nbt; ++t){
    VEigen pnttria[6];
    for(int k = 0; k < 6; ++k){
      pnttria[k] = MidPnt[ddl[t][k]];
    }

    RigMat massele = MassEle(pnttria);
    
    for(int i = 0; i < 6; ++i){
      for(int j = 0; j < 6; ++j){
	tripletlist.push_back(T(ddl[t][i],ddl[t][j],massele(i,j)));
      }
    }
  }

  SpMat result(Nbv+Nbe,Nbv+Nbe);
  result.setFromTriplets(tripletlist.begin(), tripletlist.end());
  return result;
}

RigMat RigEle(const VEigen* Nodes)
{
  Matrix<double,2,2> B_T;
  B_T(0,0) = (Nodes[1]-Nodes[0])[0];
  B_T(1,0) = (Nodes[1]-Nodes[0])[1];
  B_T(0,1) = (Nodes[2]-Nodes[0])[0];
  B_T(1,1) = (Nodes[2]-Nodes[0])[1];

  //cout << B_T << endl;

  Matrix<double,2,2> invB_T = B_T.inverse();

  //cout << invB_T << endl;

  double areaT = B_T.determinant()/2;

  VEigen alp[7]; for(int i = 0; i < 7; ++i){alp[i].resize(2);}
  VEigen alp0(2);
  alp0[0] = 1./3; alp0[1] = 1./3;
  alp[0] = alp0;
  VEigen alp1(2);
  alp1[0] = (6-sqrt(15))/21; alp1[1] = (6-sqrt(15))/21;
  alp[1] = alp1;
  VEigen alp2(2);
  alp2[0] = (6-sqrt(15))/21; alp2[1] = (9+2*sqrt(15))/21;
  alp[2] = alp2;
  VEigen alp3(2);
  alp3[0] = (9+2*sqrt(15))/21; alp3[1] = (6-sqrt(15))/21;
  alp[3] = alp3;
  VEigen alp4(2);
  alp4[0] = (6+sqrt(15))/21; alp4[1] = (6+sqrt(15))/21;
  alp[4] = alp4;
  VEigen alp5(2);
  alp5[0] = (6+sqrt(15))/21; alp5[1] = (9-2*sqrt(15))/21;
  alp[5] = alp5;
  VEigen alp6(2);
  alp6[0] = (9-2*sqrt(15))/21; alp6[1] = (6+sqrt(15))/21;
  alp[6] = alp6;
  
  VEigen wei(7);
  wei[0] = 0.225; wei[1] = (155-sqrt(15))/1200; wei[2] = (155-sqrt(15))/1200; wei[3] = (155-sqrt(15))/1200; wei[4] = (155+sqrt(15))/1200; wei[5] = (155+sqrt(15))/1200; wei[6] = (155+sqrt(15))/1200; 

  RigMat rigele;
  for(int i = 0; i < 6; ++i){
    for(int j = 0; j < 6; ++j){
      rigele(i,j) = 0;
    }
  }

  VEigen nbr1(2);
  VEigen nbr2(2);
  
  for(int i = 0; i < 6; ++i){
    for(int j = 0; j < 6; ++j){
      for(int k = 0; k < 7; ++k){
	//cout << areaT*wei[k]*( (invB_T*Gradphi(alp[k],i)).adjoint() ) * ( invB_T*Gradphi(alp[k],j) ) << endl << endl;
	//nbr1 += invB_T*Gradphi(alp[k],i); nbr2 += invB_T*Gradphi(alp[k],j);
	rigele(i,j) += areaT*wei[k]*( ((invB_T.adjoint())*Gradphi(alp[k],i)).adjoint() ) * ( (invB_T.adjoint())*Gradphi(alp[k],j) );
      }
      //cout << i << " " << j << endl << nbr1.adjoint()*nbr2 << endl << endl;
      //cout << i << " " << j << " " << rigele(i,j) << endl << endl;
    }
  }      
      /*
	C'est très étrange, lorsqu'on mets : 
	rigele(i,j) = Gradphi(alp1,i).adjoint()*Gradphi(alp1,j);
	Il n'y a aucun probleme, le membre de droite est un nombre, on peut l'affecter à rigele(i,j), mais lorsqu'on mets :
	rigele(i,j) = Gradphi(alp1,i).adjoint()*Gradphi(alp1,j)+Gradphi(alp2,i).adjoint()*Gradphi(alp2,j);
	On a une erreur à la compilation qui nous dit que le membre de droite n'est pas un nombre et donc il ne peut pas affecter à rigele(i,j)...
	Ce qui est très pénible du coup..
      
      
      rigele(i,j) = areaT*weight0*((invB_T*Gradphi(alp0,i))).adjoint()*(invB_T*Gradphi(alp0,j));
      rigele(i,j) += areaT*weight1*((invB_T*Gradphi(alp1,i))).adjoint()*(invB_T*Gradphi(alp1,j));
      rigele(i,j) += areaT*weight1*((invB_T*Gradphi(alp2,i))).adjoint()*(invB_T*Gradphi(alp2,j));
      rigele(i,j) += areaT*weight1*((invB_T*Gradphi(alp3,i))).adjoint()*(invB_T*Gradphi(alp3,j));
      rigele(i,j) += areaT*weight2*((invB_T*Gradphi(alp4,i))).adjoint()*(invB_T*Gradphi(alp4,j));
      rigele(i,j) += areaT*weight2*((invB_T*Gradphi(alp5,i))).adjoint()*(invB_T*Gradphi(alp5,j));
      rigele(i,j) += areaT*weight2*((invB_T*Gradphi(alp6,i))).adjoint()*(invB_T*Gradphi(alp6,j));      */
  
  
  //cout << rigele << endl << endl;
  return rigele;
}

SpMat Rig(VEigen* MidPnt, veci* ddl, veci* Element, int Nbv, int Nbt, int Nbe)
{
  vector<T> tripletlist;
  for(int t = 0; t < Nbt; ++t){
    VEigen pnttria[6];
    for(int k = 0; k < 6; ++k){
      pnttria[k] = MidPnt[ddl[t][k]];
    }
    
    RigMat rigele = RigEle(pnttria);
    
    for(int i = 0; i < 6; ++i){
      for(int j = 0; j < 6; ++j){
	tripletlist.push_back(T(ddl[t][i],ddl[t][j],rigele(i,j)));
      }
    }
  }

  SpMat result(Nbv+Nbe,Nbv+Nbe);
  result.setFromTriplets(tripletlist.begin(), tripletlist.end());
  
  return result;
}

DivMat DivEle(const VEigen* Nodes, int l)
{
  Matrix<double,2,2> B_T;
  B_T(0,0) = (Nodes[1]-Nodes[0])[0];
  B_T(1,0) = (Nodes[1]-Nodes[0])[1];
  B_T(0,1) = (Nodes[2]-Nodes[0])[0];
  B_T(1,1) = (Nodes[2]-Nodes[0])[1];

  Matrix<double,2,2> invB_T = B_T.inverse();

  double areaT = B_T.determinant()/2;

  Vector2d alp0(1./3,1./3);
  Vector2d alp1((6-sqrt(15))/21,(6-sqrt(15))/21);
  Vector2d alp2((6-sqrt(15))/21,(9+2*sqrt(15))/21);
  Vector2d alp3((9+2*sqrt(15))/21,(6-sqrt(15))/21);
  Vector2d alp4((6+sqrt(15))/21,(6+sqrt(15))/21);
  Vector2d alp5((6+sqrt(15))/21,(9-2*sqrt(15))/21);
  Vector2d alp6((9-2*sqrt(15))/21,(6+sqrt(15))/21);

  double weight0 = 0.225;
  double weight1 = (155-sqrt(15))/1200;
  double weight2 = (155+sqrt(15))/1200;
  
  DivMat divele;
  for(int i = 0; i < 6; ++i){
    for(int j = 0; j < 3; ++j){
	divele(i,j) = weight0*areaT*lambda(alp0,j)*(invB_T*Gradphi(alp0,i))[l] + weight1*areaT*(lambda(alp1,j)*(invB_T*Gradphi(alp1,i))[l] + lambda(alp2,j)*(invB_T*Gradphi(alp2,i))[l] + lambda(alp3,j)*(invB_T*Gradphi(alp3,i))[l]) + weight2*areaT*(lambda(alp4,j)*(invB_T*Gradphi(alp4,i))[l] + lambda(alp5,j)*(invB_T*Gradphi(alp5,i))[l] + lambda(alp6,j)*(invB_T*Gradphi(alp6,i))[l]);      
    }
  }

  return divele;
}

SpMat Div(VEigen* MidPnt, veci* ddl, veci* Element, int Nbv, int Nbt, int Nbe, int l)
{
  vector<T> tripletlist;
  for(int t = 0; t < Nbt; ++t){
    VEigen pnttria[6];
    for(int k = 0; k < 6; ++k){
      pnttria[k] = MidPnt[ddl[t][k]];
    }
    
    DivMat divele = DivEle(pnttria,l);

    for(int i = 0; i < 6; ++i){
      for(int j = 0; j < 3; ++j){
	tripletlist.push_back(T(ddl[t][i],ddl[t][j],divele(i,j)));
      }
    }
  }

  SpMat result(Nbv+Nbe,Nbv);
  result.setFromTriplets(tripletlist.begin(), tripletlist.end());
  return result;
}

PenMat PenEle(const VEigen* Nodes)
{
  VEigen edge0 = Nodes[0]-Nodes[1];
  VEigen edge1 = Nodes[0]-Nodes[2];

  double areaT = (1./2)*(edge0[0]*edge1[1]-edge0[1]*edge1[0]);

  PenMat penele;
  for(int i = 0; i < 3; ++i){
    for(int j = 0; j < 3; ++j){
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

SpMat Pen(VEigen* MidPnt, veci* ddl, veci* Element, int Nbv, int Nbt, int Nbe)
{
  vector<T> tripletlist;
  for(int t = 0; t < Nbt; ++t){
    VEigen pnttria[6];
    for(int k = 0; k < 6; ++k){
      pnttria[k] = MidPnt[ddl[t][k]];
    }
        
    PenMat penele = PenEle(pnttria);

    for(int i = 0; i < 3; ++i){
      for(int j = 0; j < 3; ++j){
	tripletlist.push_back(T(ddl[t][i],ddl[t][j],penele(i,j)));
      }
    }
  }

  SpMat result(Nbv,Nbv);
  result.setFromTriplets(tripletlist.begin(), tripletlist.end());
  return result;
}


SpMat AssemSys(VEigen* MidPnt, veci* ddl, veci* Element, const veci& IdNodesBound, const veci& IdNodesLbl, int Nbv, int Nbt, int Nbe, double eps, double nu, double alp, VEigen& RightHandSide)
{
  vector<T> tripletlist;
  for(int t = 0; t < Nbt; ++t){
    VEigen pnttria[6];
    for(int k = 0; k < 6; ++k){
      pnttria[k] = MidPnt[ddl[t][k]];
    }

    RigMat massele = MassEle(pnttria);
    RigMat rigele = RigEle(pnttria);
    PenMat penele = PenEle(pnttria);
    DivMat divele1 = DivEle(pnttria,0);
    DivMat divele2 = DivEle(pnttria,1);

    // cout << rigele << endl << endl;

    for(int i = 0; i < 6; ++i){
      for(int j = 0; j < 6; ++j){
        tripletlist.push_back( T(ddl[t][i],ddl[t][j],nu*rigele(i,j)) );
	tripletlist.push_back( T(ddl[t][i] + (Nbv+Nbe),ddl[t][j] + (Nbv+Nbe),nu*rigele(i,j)) );
	tripletlist.push_back( T(ddl[t][i],ddl[t][j],alp*massele(i,j)) );
	tripletlist.push_back( T(ddl[t][i] + (Nbv+Nbe),ddl[t][j] + (Nbv+Nbe),alp*massele(i,j)) );
      }
    }
    
    for(int i = 0; i < 6; ++i){
      for(int j = 0; j < 3; ++j){
	tripletlist.push_back( T(2*(Nbv+Nbe)+ddl[t][j],ddl[t][i],-divele1(i,j)) );
	tripletlist.push_back( T(ddl[t][i],2*(Nbv+Nbe)+ddl[t][j],-divele1(i,j)) );
	tripletlist.push_back( T((Nbv+Nbe)+ddl[t][i],2*(Nbv+Nbe)+ddl[t][j],-divele2(i,j)) );
	tripletlist.push_back( T(2*(Nbv+Nbe)+ddl[t][j],(Nbv+Nbe)+ddl[t][i],-divele2(i,j)) );
      }
    }
    
    for(int i = 0; i < 3; ++i){
      for(int j = 0; j < 3; ++j){
	tripletlist.push_back( T(2*(Nbv+Nbe)+ddl[t][i],2*(Nbv+Nbe)+ddl[t][j],-eps*penele(i,j)) );
      }
    }
  }


  SpMat matglob(3*Nbv+2*Nbe,3*Nbv+2*Nbe);
  matglob.setFromTriplets(tripletlist.begin(), tripletlist.end());

  //for(int i = 2*(Nbv+Nbe); i < 2*(Nbv+Nbe)+Nbv; ++i){
  //  matglob.coeffRef(i,i) = -eps;}

  for(int i = 0; i < IdNodesBound.size(); ++i){
    matglob.coeffRef(IdNodesBound[i],IdNodesBound[i]) = pow(10,10);
    matglob.coeffRef(IdNodesBound[i]+Nbv+Nbe,IdNodesBound[i]+Nbv+Nbe) = pow(10,10);
    RightHandSide(IdNodesBound[i]) = 0;
    RightHandSide(IdNodesBound[i]+Nbv+Nbe) = 0;
  }

  for(int i = 0; i < IdNodesLbl.size(); ++i){
    //cout << g(MidPnt[IdNodesLbl[i]]) << endl;
    RightHandSide(IdNodesLbl[i]) = pow(10,10);//*g(MidPnt[IdNodesLbl[i]]);
    //RightHandSide(IdNodesLbl[i]+Nbv+Nbe) = pow(10,30);
  }
  
  return matglob;
}

/*
SpMat AssemSysNS(veci* trian, double dt, const VEigen& u_n, VEigen* MidPnt, veci* ddl, veci* Element, const veci& IdNodesBound, const veci& IdNodesLbl, int Nbv, int Nbt, int Nbe, double eps, double nu, VEigen& RightHandSide)
{
  vector<T> tripletlist;
  for(int t = 0; t < Nbt; ++t){
    VEigen pnttria[6];
    for(int k = 0; k < 6; ++k){
      pnttria[k] = MidPnt[ddl[t][k]];
    }

    RigMat massele = MassEle(pnttria);
    RigMat rigele = RigEle(pnttria);
    DivMat divele1 = DivEle(pnttria,0);
    DivMat divele2 = DivEle(pnttria,1);
    Matrix<double,7,2> righthandsideele;
    righthandsideele = RightHandSideEle(dt,u_n,pnttria,Nbv,Nbe,MidPnt,t,ddl,trian);
    
    for(int i = 0; i < 6; ++i){
      for(int j = 0; j < 6; ++j){
        tripletlist.push_back( T(ddl[t][i],ddl[t][j],nu*rigele(i,j)) );
	tripletlist.push_back( T(ddl[t][i] + (Nbv+Nbe),ddl[t][j] + (Nbv+Nbe),nu*rigele(i,j)) );
	tripletlist.push_back( T(ddl[t][i],ddl[t][j],(1./dt)*massele(i,j)) );
	tripletlist.push_back( T(ddl[t][i] + (Nbv+Nbe),ddl[t][j] + (Nbv+Nbe),(1./dt)*massele(i,j)) );
      }
    }
    
    for(int i = 0; i < 6; ++i){
      for(int j = 0; j < 3; ++j){
	tripletlist.push_back( T(2*(Nbv+Nbe)+ddl[t][j],ddl[t][i],-divele1(i,j)) );
	tripletlist.push_back( T(ddl[t][i],2*(Nbv+Nbe)+ddl[t][j],-divele1(i,j)) );
	tripletlist.push_back( T((Nbv+Nbe)+ddl[t][i],2*(Nbv+Nbe)+ddl[t][j],-divele2(i,j)) );
	tripletlist.push_back( T(2*(Nbv+Nbe)+ddl[t][j],(Nbv+Nbe)+ddl[t][i],-divele2(i,j)) );
      }
    }
    
    
    for(int i = 0; i < 3; ++i){
      for(int j = 0; j < 3; ++j){
	tripletlist.push_back( T(2*(Nbv+Nbe)+ddl[t][i],2*(Nbv+Nbe)+ddl[t][j],-eps*massele(i,j)) );
      }
    }

    for(int i = 0; i < 6; ++i){
      RightHandSide[Element[t][i]] += righthandsideele(i,0);
      RightHandSide[Element[t][i]+Nbv+Nbe] += righthandsideele(i,1);
    }
  
    
  }


  SpMat matglob(3*Nbv+2*Nbe,3*Nbv+2*Nbe);
  matglob.setFromTriplets(tripletlist.begin(), tripletlist.end());

  //for(int i = 2*(Nbv+Nbe); i < 2*(Nbv+Nbe)+Nbv; ++i){
  //  matglob.coeffRef(i,i) = -eps;}

  for(int i = 0; i < IdNodesBound.size(); ++i){
    matglob.coeffRef(IdNodesBound[i],IdNodesBound[i]) = pow(10,7);
    matglob.coeffRef(IdNodesBound[i]+Nbv+Nbe,IdNodesBound[i]+Nbv+Nbe) = pow(10,7);
  }

  for(int i = 0; i < IdNodesBound.size(); ++i){
    RightHandSide(IdNodesBound[i]) = 0;
  }
  
  for(int i = 0; i < IdNodesLbl.size(); ++i){
    //cout << g(MidPnt[IdNodesLbl[i]]) << endl;
    RightHandSide(IdNodesLbl[i]) = pow(10,7);//*g(MidPnt[IdNodesLbl[i]]);
    //RightHandSide(IdNodesLbl[i]+Nbv+Nbe) = pow(10,30);
  }
  
  
  
  return matglob;
}
*/

int StructMatGlob(const SpMat& M, bool coef)
{
  int result = 0;
  ofstream f("matrice.dat");
  for(int i = 0; i < M.rows(); ++i){
    for(int j = 0; j < M.cols(); ++j){
      if(coef == 0){
	if(abs(M.coeff(i,j)) > pow(10,-12)){
	  result += 1;
	  cout << "x";
	}
	else{
	  cout << "o";
	}
      }
      else{
	if(abs(M.coeff(i,j)) > pow(10,-7)){
	  result += 1;
	  f << i << "\t" << j << "\t" << M.coeff(i,j) << endl;
	}
      }
    }
    f << endl;
  }
  cout << result << endl;
  return result;
}

void SavetoPlot(int Nbe, int Nbv, int Nbt, const veci* ddl, const VEigen& x, const char* filename)
{
  ofstream f(filename);

  for(int t = 0; t < Nbt; ++t){
    
    for(int j = 0; j < 6; ++j){
      if(j < 3){
	f << x[ddl[t][j]] << " ";
      }
      else{
	f << x[ddl[t][(j+1)%3+3]] << " ";
      }
    }
    
    for(int j = 0; j < 6; ++j){
      if(j < 3){
	f << x[ddl[t][j]+Nbv+Nbe] << " ";
      }
      else{
	f << x[ddl[t][(j+1)%3+3]+Nbv+Nbe] << " ";
      }
    }
        
    for(int j = 0; j < 3; ++j){
      f << x[ddl[t][j]+2*(Nbv+Nbe)] << " ";
    }



    

    /*
    for(int j = 0; j < 6; ++j){
      f << x[ddl[t][j]] << " ";
    }
    for(int j = 0; j < 6; ++j){
      f << x[ddl[t][j]+Nbv+Nbe] << " ";
    }
    for(int j = 0; j < 3; ++j){
      f << x[ddl[t][j]+2*(Nbv+Nbe)] << " ";
    }
    /*
    for(int j = 0; j < 3; ++j){
      //cout << ddl[t][(j+1)%3+3]+Nbv+Nbe << endl;
      f << x[ddl[t][(j+1)%3+3]] << " ";
    }
    for(int j = 0; j < 3; ++j){
      //cout << ddl[t][(j+1)%3+3]+Nbv+Nbe << endl;
      f << x[ddl[t][(j+1)%3+3]+(Nbv+Nbe)] << " ";
    }
    */
    f << endl;
  }

  /*
  for(int t = 0; t < Nbt; ++t){
    for(int j =0; j < 3; ++j){
      f << x[ddl[t][j]] << " " << x[ddl[t][j]+Nbv+Nbe] << " " << x[ddl[t][j]+2*(Nbv+Nbe)] << " ";
    }
    for(int j = 0; j < 3; ++j){
      f << x[ddl[t][(j+1)%3+3]] << " ";
    }
    for(int j = 0; j < 3; ++j){
      f << x[ddl[t][(j+1)%3+3]+Nbv+Nbe] << " ";
    }
    f << endl;
  }
  */
}

/*
void SavetoPlot(int Nbv, int Nbe, int Nbt, const VEigen* MidPnt, const veci* ddl, const VEigen& x)
{
  ofstream u1("u1cpp.txt");
  ofstream u2("u2cpp.txt");
  ofstream u("ucpp.txt");
  ofstream p("pcpp.txt");

  u1 << Nbv+Nbe << "\n";
  u2 << Nbv+Nbe << "\n";
  u << Nbv+Nbe << "\n";
  p << Nbv << "\n";
  
  for(int t = 0; t < Nbt; ++t){
    for(int i = 0; i < 6; ++i){
      u1 << MidPnt[ddl[t][i]][0] << "\t"  << MidPnt[ddl[t][i]][1] << "\t" << x[ddl[t][i]] << "\n";
      u2 << MidPnt[ddl[t][i]][0] << "\t" << MidPnt[ddl[t][i]][1] << "\t" << x[ddl[t][i]+Nbv+Nbe] << "\n";
      u << MidPnt[ddl[t][i]][0] << "\t" << MidPnt[ddl[t][i]][1] <<"\t" << x[ddl[t][i]] << "\t" << x[ddl[t][i]+Nbv+Nbe] << "\n";
      if(i<3){
	p << MidPnt[ddl[t][i]][0] << "\t" << MidPnt[ddl[t][i]][1] << "\t" << x[ddl[t][i]+2*(Nbv+Nbe)] << "\n";
      }
    }
    u1 << endl << endl;
    u2 << endl << endl;
    u << endl << endl;
    p << endl << endl;
  }
}
*/

/*
int WhereWasI(int Nbv, int Nbe, VEigen* Nodes, VEigen& a_l, int n_t, veci* Element, veci* trian, VEigen& u_n, double& dt)
{
  VEigen transla(2);

  // On calcul le translaté de a_l en fonction de la vitesse au pas de temps précédent
  
  transla[0] = (u_n[Element[n_t][0]]+u_n[Element[n_t][1]]+u_n[Element[n_t][2]])/3;
  transla[1] = (u_n[Element[n_t][0]+Nbv+Nbe]+u_n[Element[n_t][1]+Nbv+Nbe]+u_n[Element[n_t][2]+Nbv+Nbe])/3;
  VEigen a_l_tr(2);
  a_l_tr = a_l - dt*transla;
  ////////////////////////////////////////////////////////////////////////////////////
  
  for(int i = 0; i < 3; ++i){
    
    int k = trian[Element[n_t][i]][0]; // Pour chaque sommet du triangle du quel est issu a_l, on regarde le nombre de voisins de ce sommets k

    for(int j = 1; j < k+1; ++j){

      // Pour chaque voisins, on vas checker voir si le translater de a_l est dans ce voisin ou pas

      int num_neighb = trian[Element[n_t][i]][j]; // numéro du triangle voisin que l'on regarde
      Matrix<double,2,2> B_T;
      B_T(0,0) = (Nodes[Element[num_neighb][1]]-Nodes[Element[num_neighb][0]])[0];
      B_T(1,0) = (Nodes[Element[num_neighb][1]]-Nodes[Element[num_neighb][0]])[1];
      B_T(0,1) = (Nodes[Element[num_neighb][2]]-Nodes[Element[num_neighb][0]])[1];
      B_T(1,1) = (Nodes[Element[num_neighb][2]]-Nodes[Element[num_neighb][0]])[1];
      Matrix<double,2,2> inv_B_T = B_T.inverse();
      VEigen a = inv_B_T*(a_l_tr-Nodes[Element[num_neighb][0]]); // projeté du a_l translaté sur le triangle de référence
      if( (phi(a,0) >= 0) & (phi(a,1) >= 0) ){ // test si on est dans le bon triangle voisin
	return num_neighb;
      }
    }
  }
}


Matrix<double,6,2> RightHandSideEle(double dt, const VEigen& u_n, const VEigen* Nodes, int Nbv, int Nbe, const VEigen* NodesP2, int n_t, const veci* Element, veci* trian)
{
  Matrix<double,2,2> B_T;
  B_T(0,0) = (Nodes[1]-Nodes[0])[0];
  B_T(1,0) = (Nodes[1]-Nodes[0])[1];
  B_T(0,1) = (Nodes[2]-Nodes[0])[0];
  B_T(1,1) = (Nodes[2]-Nodes[0])[1];

  double areaT = B_T.determinant()/2;

  Matrix<double,2,2> inv_B_T = B_T.inverse();

  VEigen alp[7]; for(int i = 0; i < 7; ++i){alp[i].resize(2);}
  VEigen alp0(2);
  alp0[0] = 1./3; alp0[1] = 1./3;
  VEigen alp1(2);
  alp1[0] = (6-sqrt(15))/21; alp1[1] = (6-sqrt(15))/21;
  VEigen alp2(2);
  alp2[0] = (6-sqrt(15))/21; alp2[1] = (9+2*sqrt(15))/21;
  VEigen alp3(2);
  alp3[0] = (9+2*sqrt(15))/21; alp3[1] = (6-sqrt(15))/21;
  VEigen alp4(2);
  alp4[0] = (6+sqrt(15))/21; alp4[1] = (6+sqrt(15))/21;
  VEigen alp5(2);
  alp5[0] = (6+sqrt(15))/21; alp5[1] = (9-2*sqrt(15))/21;
  VEigen alp6(2);
  alp6[0] = (9-2*sqrt(15))/21; alp6[1] = (6+sqrt(15))/21;
  
  alp[0] = alp0;
  alp[1] = alp1;
  alp[2] = alp2;
  alp[3] = alp3;
  alp[4] = alp4;
  alp[5] = alp5;
  alp[6] = alp6;

  veci wherewasi(7);
  for(int i = 0; i < 7; ++i){wherewasi[i] = WhereWasI(Nbv,Nbe,NodesP2,B_T*alp[i]+Nodes[0],n_t,Element,trian,u_n,dt);}

  VEigen w_k(7);
  w_k[0] = 0.225; w_k[1] = (155-sqrt(15))/1200; w_k[2] = (155-sqrt(15))/1200; w_k[3] = (155-sqrt(15))/1200; w_k[4] = (155+sqrt(15))/1200; w_k[5] = (155+sqrt(15))/1200; w_k[6] = (155+sqrt(15))/1200; 
  // double weight0 = 0.225;
  // double weight1 = (155-sqrt(15))/1200;
  // double weight2 = (155+sqrt(15))/1200;
  
  Matrix<double,6,2> result;
  for(int i = 0; i < 6; ++i){
    for(int k = 0; k < 6; ++k){
      for(int j = 0; j < 5; ++j){
	result(i,0) += w_k[k]*u_n[Element[wherewasi[k]][j]]*phi(B_T*alp[k]+Nodes[0],j)*phi(B_T*alp[k]+Nodes[0],i);
	result(i,1) += w_k[k]*u_n[Element[wherewasi[k]][j]]+Nbv+Nbe]*phi(B_T*alp[k]+Nodes[0],j)*phi(B_T*alp[k]+Nodes[0],i);
      }
    }
    result(i,0) *= areaT;
  }
  return result;
}
*/
