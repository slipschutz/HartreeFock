

#include <iostream>
#include <vector>
#include "TMath.h"
#include "TMatrixD.h"
#include "TVectorD.h"

using namespace std;

int nStates=8;
int nParticles=4;
int HFIterations=100;


class Operator {

public:
  Operator();
  virtual double Apply(vector<bool>& ){cout<<"Nope"<<endl;}
  virtual void Print(){cout<<"No Print for you"<<endl;}
  
private:
};

Operator::Operator() 
{}

class Term {
public:
  double Coef;
  vector <Operator*> theOperators;
  int GetNumOps(){return theOperators.size();}
};


vector <Term> thePotential;

class Raise : public Operator{
public:
  Raise();
  Raise(int);
  double Apply(vector<bool>& s){
    if (sigma ==-1){
      hiddenP=p;
    }else if (sigma==1){
      hiddenP=4+p;
    }
    if (s[hiddenP] ==0){
      s[hiddenP]=1;
      return 1;
    }else if (s[hiddenP] ==1){
      return 0;
    }
  }
 void Print(){
    cout<<"R"<<p<<sigma;
  }
  int sigma;
  int p;
private:
  int hiddenP;
};

class Lower : public Operator{
public:
  Lower();
  Lower(int);
  double Apply(vector<bool>& s){

    if (sigma ==-1){
      hiddenP=p;
    }else if (sigma==1){
      hiddenP=4+p;
    }
    if (s[hiddenP] ==1){
      s[hiddenP]=0;
      return 1;
    }else if (s[hiddenP] ==0){
      return 0;
    }
  }
  void Print(){
    cout<<"L"<<p<<sigma;
  }
  int sigma;
  int p;
private:
  int hiddenP;
};

Lower::Lower()
{}

Raise::Raise()
{}

Lower::Lower(int n){
  if(n<=3){
    p=n;
    sigma=-1;
  }else{
    p=n-4;
    sigma=1;
  }
}


Raise::Raise(int n){
  if(n<=3){
    p=n;
    sigma=-1;
  }else{
    p=n-4;
    sigma=1;
  }
}


void BuildPotential(){
  
  for (int i=0;i<nParticles;i++){
    for (int j=0;j<nParticles;j++){
      for (int sigma=-1;sigma<2;sigma=sigma+2){
	Term aTerm;
	
	Lower *low1 = new Lower();
	low1->sigma = -sigma;
	low1->p=i;

	Lower *low2= new Lower();;
	low2->sigma=-sigma;
	low2->p=j;

	Raise *raise1=new Raise();
	raise1->sigma=sigma;
	raise1->p=j;
	Raise *raise2=new Raise();
	raise2->sigma=sigma;
	raise2->p=i;

	aTerm.theOperators.push_back(low1);
	aTerm.theOperators.push_back(low2);
	aTerm.theOperators.push_back(raise1);
	aTerm.theOperators.push_back(raise2);
	thePotential.push_back(aTerm);
      }
    }
  }


}

    // for (int a=0;a<theCurrentTerm.theOperators.size();a++){
    //   theCurrentTerm.theOperators[theCurrentTerm.theOperators.size()-1-a]->Print();
    //   cout<<" ";
    // }


double matrixElement(int alpha,int beta,int gamma,int delta){
  vector <bool> theState(8,0);
  double total=0;
  for (int i=0;i<thePotential.size();i++){
    Term theCurrentTerm =thePotential[i];
    //need to add on the operators assciated with gamma and delta
    //they should be at the beginning of the vector of operators
    //should be gamma alpha POTENTIAL beta delta
    Raise firstRaiseForState(beta);
    Raise secondRaiseForState(delta);
    Lower firstLowerForState(alpha);
    Lower secondLowerForState(gamma);
    theCurrentTerm.theOperators.insert(theCurrentTerm.theOperators.begin(),&firstRaiseForState);
    theCurrentTerm.theOperators.insert(theCurrentTerm.theOperators.begin(),&secondRaiseForState);
    theCurrentTerm.theOperators.push_back(&firstLowerForState);
    theCurrentTerm.theOperators.push_back(&secondLowerForState);
    
    int numOps =theCurrentTerm.theOperators.size();
    bool isAZero=false;
    for (int q=0;q<numOps;q++){
      double ret=theCurrentTerm.theOperators[q]->Apply(theState);
      if (ret == 0){
	q=numOps+1;//end loop
	isAZero=true;
      }
    }
    theState.clear();
    theState.resize(8,0);
    if (isAZero==false){
      total++;
    }
  }

  return total*(-1.0/6.0);
}

double h0(int alpha,int gamma){

  if (alpha != gamma){
    return 0;
  } 

  if (gamma <=3){
    return -1;
  }else {
    return 1;
  }
}

void hartreeFockrun() {
    double interaction;
    // --------------- Setting up the HF-hamiltonian using C = 1 as guess
    TMatrixD H(nStates,nStates);
    TVectorD E(nStates);
    TVectorD ePrev(nStates);
    TMatrixD C(nStates,nStates);
    for (int i=0;i<nStates;i++){
      for (int j=0;j<nStates;j++){
	if (i==j){
	  C[i][j]=1;
	}else{
	  C[i][j]=0;
	}
      }
    }
    TVectorD diff(nStates);

    double threshold =0.00000001;

    // Hartree-Fock loop
    int hfIt = 0;
    while (hfIt < HFIterations) {
        cout << "iteration = " << hfIt << endl;

        for (int alpha = 0; alpha < nStates; alpha++) {
            for (int gamma = 0; gamma < nStates; gamma++) {
                interaction = 0;
                for (int p = 0; p < nParticles; p++) {
                    for (int beta = 0; beta < nStates; beta++) {
                        for (int delta = 0; delta < nStates; delta++) {
			  interaction += C[p][gamma] * C[p][delta] * matrixElement(alpha, gamma, beta, delta);//switched alpha and beta
                        }
                    }
                }
                H[alpha][gamma] =h0(alpha, gamma) + interaction;
		H[gamma][alpha]=H[alpha][gamma];  
            }
        }
	H.Print();
        //Computing the HF one-body energies
        C=H.EigenVectors(E);
        // Transposing the vectors
        C.Transpose(C);
        hfIt++;
        // Convergence test

	diff = E - ePrev;
	cout<<"Diff is "<<TMath::Abs(diff.Max())<<endl;
        if (TMath::Abs(diff.Max()) < threshold){
	  cout<<"FINISHED"<<endl;
	  break;
	}
        ePrev = E;
    }
    //    double E0 = calcEnergy(C);
    //    cout << "Final energy E = " << E0 << " after " << hfIt << " iterations, error < " << threshold << endl;*/
}





int main(){
  BuildPotential();

  hartreeFockrun();

  //  matrixElement(4,0,5,1);
  return 0;

  for (int i=0;i<thePotential.size();i++){

    int num=thePotential[i].theOperators.size();
    for (int j=0;j<num;j++){
      thePotential[i].theOperators[num-j-1]->Print();
      cout<<" ";
    }
    cout<<endl;
  }

  return 0;  


  return 0;
  TMatrixD aMatrix(2,2);
  aMatrix[0][0]=1;
  aMatrix[0][1]=13;
  aMatrix[1][0]=1;
  aMatrix[1][1]=3;

  TMatrixD TransposedMatrix(2,2);
  TransposedMatrix.Transpose(aMatrix);
  //  TMatrixD theVecs = aMatrix.EigenVectors(theValues);

  cout<<"A MATRIX IS "<<endl;
  aMatrix.Print();
  cout<<"TRANSPOSED IS "<<endl;

  TransposedMatrix.Print();

  
  
  return 0;
}
