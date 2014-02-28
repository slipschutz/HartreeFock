

#include <iostream>
#include <vector>

#include "TMatrixD.h"
#include "TVectorD.h"

using namespace std;

int nStates=5;
int HFIterations=100;

double MatrixElements(int alpha,int beta,int gamma,int delta){

  return 1;
}

double h0(int alpha,int gamma){

  return 1;
}

void hartreeFockrun() {
    double interaction;
    // --------------- Setting up the HF-hamiltonian using C = 1 as guess, Armadillo is used for vectors
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
    TVectorD diff;



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
                            interaction += C[p][beta] * C[p][delta] * matrixElement(alpha, beta, gamma, delta);
                        }
                    }
                }
                H[alpha][gamma] =h0(alpha, gamma) + interaction;
		H[gamma][alpha]=H[alpha][gamma];  
            }
        }
        //Computing the HF one-body energies
        C=H.EigenVectors(E);
        // Transposing the vectors
        C = C.Transpose(C);
        hfIt++;
        // Convergence test
        diff = E - ePrev;
        if (abs(diff.max()) < threshold)
            break;
        ePrev = E;
    }
    //    double E0 = calcEnergy(C);
    cout << "Final energy E = " << E0 << " after " << hfIt << " iterations, error < " << threshold << endl;*/
}




int main(){
  hartreeFockrun();
  TMatrixD aMatrix(2,2);
  aMatrix[0][0]=1;
  aMatrix[0][1]=13;
  aMatrix[1][0]=1;
  aMatrix[1][1]=3;

  TVectorD theValues;

  TMatrixD theVecs = aMatrix.EigenVectors(theValues);

  aMatrix.Print();

  theValues.Print();
  
  return 0;
}
