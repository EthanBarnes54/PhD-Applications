/* 

I don't normally put comments onto my code, I aim to use variable naims whihc make it explicitely clear as to what th programme is doing; 
but I feel it is convenient to do so for your cosideration and understanding, as this problem may be alien to you. I as yet, have not finalised the construction of the programme 
which I will run alongside my MATLAb code to provide speedup via parallel computation. As such, the below is a translation, where appropriate of my project code into C. 
I believe this demosntrates my ability to learn new programming syntax and apply universal coding methods across platforms. 

The current phase of my reseacrh proejct is to emulate findings of a previous piece of research which can be found here: https://arxiv.org/abs/2108.13554 

As a surface level overview; I am performing random walks through Hilbert space to investigate the adiabtaic evolution of pure-multiqubit states, to ascertain 
the behaviour and nature of (adiabatic) "quantum" computation given an evolution/ run time.

The programme below is translated from matlab, and describes the process implemented in the above piece of research to model this adaiabtic evolution, 
"A random walk through Hilbert Space".

*/


// Including all relevant packages and Libraries

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <string.h>

// Defining constant(s) which will be used throughout

#define PI 3.141592653589793 

// Defining the "Protoypes" (Functions) whihc are integral to the simulation at hand

void State_Generator(double *Parameters, int Dimensions, double complex *State);
void Partial_Differentiator(double *Parameters, int Dimensions, double complex *Derivative_States);
void Random_Walk_Vector(int Dimensions, double *RWV);
void Metric_Generating_Function(int Dimensions, double complex *Original_State, double *Initial_Parameters, double *MGF);
void Step_Generator(double Step_Size, int Dimensions, double *Initial_Parameters, double complex *First_State, double *Step);
int Random_Walk(int Number_of_Qubits, int Number_of_Steps, double Step_Size, double Precision);

// Below are the definitions of the implemented functions for simulation

//A 'random' number generator, which generates a normally distributed random number belonging to [-1,1]
double Random_Number_Generator(){

    return 2.0 * ((double)randn() / RAND_MAX) - 1.0;
}

// A fucntion to construct a pure quantum state according to state vector representation on the "Bloch Sphere" (N-dimensional Spherical Coordinates)

void State_Generator(double *Parameters, int Dimensions, double complex *State){

    double Array_of_Sines[Dimensions - 1];

    for (int i = 0; i <= Dimensions - 1; i++){

        Array_of_Sines[i] = 1.0;
    }

    for (int i = 0; i < Dimensions - 1; i++){

        State[i] = cexp(I * Parameters[Dimensions - 1 + i]) * cos(Parameters[i]);
        for (int j = 0; j < i; j++) {
            State[i] *= sin(Parameters[j]);
        }
        Array_of_Sines[i] = sin(Parameters[i]);
    }

    State[Dimensions-1] = 1.0;
    for (int i = 0; i < Dimensions; i++){

        State[Dimensions-1] *= Array_of_Sines[i];
    }
}

// A fucntion which performs the derivatives (Analytically) of a pure quantum state, returning every derivative wrt. each parameter, as a D x 2(D-1) array as required by the problem

void Partial_Differentiator(double *Parameters, int Dimensions, double complex *Derivative_States){

    double *Array_of_Sines = (double *)calloc(Dimensions, sizeof(double));
    double complex *Derivative_State = (double complex *)calloc(Dimensions, sizeof(double complex));
    double Product_of_Sines = 1.0;

        for (int i = 0; i < Dimensions; i++) {
            for (int j = 0; j < 2 * (Dimensions - 1); j++) {

             Derivative_States[i * (2 * (Dimensions - 1)) + j] = 0.0 + 0.0 * I;
            }
        }

    for (int i = 0; i < Dimensions - 1; i++) {
        
        for (int k = 0; k < Dimensions - 1; k++) {
            Array_of_Sines[k] = 1.0;
        }
 
        for (int k = 0; k <= i; k++) {
            Array_of_Sines[k] = sin(Parameters[k]);
        }

        double Product_of_Sines = 1.0;
        for (int k = 0; k <= i; k++) {
            Product_of_Sines *= Array_of_Sines[k];
        }
        Derivative_State[i] = -1 * cexp(I * Parameters[i + (Dimensions - 1)]) * Product_of_Sines;

        Array_of_Sines[i] = 1.0;

        for (int j = i + 1; j < Dimensions - 1; j++) {
            Product_of_Sines = 1.0;
            for (int k = 0; k < Dimensions - 1; k++) {
                Product_of_Sines *= Array_of_Sines[k];
            }
            Derivative_State[j] = cexp(I * Parameters[j + (Dimensions - 1)]) * Product_of_Sines * cos(Parameters[i]) * cos(Parameters[j]);

            Array_of_Sines[j] = sin(Parameters[j]);
        }

        Product_of_Sines = 1.0;
        for (int k = 0; k < Dimensions - 1; k++) {
            Product_of_Sines *= Array_of_Sines[k];
        }
        Derivative_State[Dimensions - 1] = Product_of_Sines * cos(Parameters[i]);

        for (int k = 0; k < Dimensions; k++) {
            Derivative_States[k * (2 * (Dimensions - 1)) + i] = Derivative_State[k];
        }
    }
   
    for (int i = Dimensions - 1; i < 2 * (Dimensions - 1); i++) {
        for (int k = i - (Dimensions - 1); k < Dimensions; k++) {
            Array_of_Sines[k] = 1.0;
        }

        double Product_of_Sines = 1.0;

        for (int k = i - (Dimensions - 1); k < Dimensions; k++) {
            Product_of_Sines *= Array_of_Sines[k];
        }
         Derivative_States[(i - (Dimensions - 1)) * (2 * (Dimensions - 1)) + i] = I * cexp(I * Parameters[i]) * Product_of_Sines * cos(Parameters[i - (Dimensions - 1)]);
    
        Array_of_Sines[i - (Dimensions - 1)] = sin(Parameters[i - (Dimensions - 1)]);
    }

    free(Array_of_Sines);
    free(Derivative_State);
}

// A function which constructs the Fubini-Study Metric (Used to calculate the curvature of the projective quantum space) required for the contruction of the random parameter evolutions

void Metric_Generating_Function(int Dimensions, double complex *Original_State, double *Initial_Parameters, double *MGF){
    
    double *FS_Matrix = (double *)calloc(Dimensions * Dimensions, sizeof(double));
    double complex *Partial_Derivative_Matrix = (double complex *)calloc(Dimensions * 2 * (Dimensions - 1), sizeof(double complex));
    
    Partial_Differentiator(Initial_Parameters, Dimensions, Partial_Derivative_Matrix);

    for (int i = 0; i < 2 * (Dimensions - 1); i++) {
        for (int j = 0; j < 2 * (Dimensions - 1); j++){

            double complex Curvature_Elements = 0.0;
            double complex Complex_Projection = 0.0;
            double complex Real_Projection = 0.0;

            for (int k = 0; k < Dimensions; k++) {
                Curvature_Elements += Partial_Derivative_Matrix[k * 2 * (Dimensions - 1) + i] * conj(Partial_Derivative_Matrix[k * 2 * (Dimensions - 1) + j]);
                Complex_Projection += Partial_Derivative_Matrix[k * 2 * (Dimensions - 1) + i] * conj(Original_State[k]);
                Real_Projection += Original_State[k] * conj(Partial_Derivative_Matrix[k * 2 * (Dimensions - 1) + j]);
            }

            FS_Matrix[i * 2 * (Dimensions - 1) + j] = creal(Curvature_Elements - (Complex_Projection * Real_Projection));
        }
    }

    for (int i = 0; i < 2 * (Dimensions - 1); i++){
        for (int j = 0; j < 2 * (Dimensions - 1); j++) {
            MGF[i * 2 * (Dimensions - 1) + j] = 0.5 * (FS_Matrix[i * 2 * (Dimensions - 1) + j] + FS_Matrix[j * 2 * (Dimensions - 1) + i]);
        }
    }

    free(FS_Matrix);
    free(Partial_Derivative_Matrix);
}

// A function that constructs a vector of normally distributed random direction, acting as the "direction" for the random walker

void Random_Walk_Vector(int Dimensions, double *RWV){

    double sum = 0.0;

    for (int i = 0; i < 2 * (Dimensions - 1); i++){

        RWV[i] = Random_Number_Generator();
        sum += RWV[i] * RWV[i];
    }
    for (int i = 0; i < 2 * (Dimensions - 1); i++){

        RWV[i] /= sqrt(sum);
    }
}

// A non-functional prototype which when implemented into matlab generates evolutions to the set of parameter which describe 
// a quantum state using Eigen-decomposition of the Fubini-Study Metric (This is not finished as Eigen-decomposition is native to matlab; and as such is not required for the project).

void Step_Generator(double Step_Size, int Dimensions, double *Initial_Parameters, double complex *First_State, double *Step){

    int Array_Dimensions = 2*(Dimensions - 1);
    double RWV[Array_Dimensions];
    double Metric[Array_Dimensions * Array_Dimensions];
    double Eigen_Vectors[Array_Dimensions * Array_Dimensions];
    double Eigen_Values[Array_Dimensions];  

    Random_Walk_Vector(Dimensions, RWV);
    Metric_Generating_Function(Dimensions, First_State, Initial_Parameters, Metric);

    // Placeholder for actual Metric Generating Function and Eigen decomposition Which will defnitely be performed in MATLAB as it is native

    for (int i = 0; i < Array_Dimensions; i++){

        Step[i] = (Step_Size / sqrt(Eigen_Values[i])) * RWV[i];
    }
}

// The final prototye which actually performs the random walk given a certain number of steps, step size and number of qubits.

int Random_Walk(int Number_of_Qubits, int Number_of_Steps, double Step_Size, double Precision){

    int Dimension_of_Hilbert_Space = (int)pow(2, Number_of_Qubits);
    double complex State_Vectors[Number_of_Steps][Dimension_of_Hilbert_Space];
    double State_Parameters[2 * (Dimension_of_Hilbert_Space - 1)][Number_of_Steps];

    for (int i = 0; i < 2 * (Dimension_of_Hilbert_Space - 1); i++){

        State_Parameters[i][0] = 0.0;
    }
    State_Generator(State_Parameters[0], Dimension_of_Hilbert_Space, State_Vectors[0]);

    for (int t = 1; t <= Number_of_Steps; t++){

        double Step[2 * (Dimension_of_Hilbert_Space - 1)];
        Step_Generator(Step_Size, Dimension_of_Hilbert_Space, State_Parameters[t - 1], State_Vectors[t - 1], Step);

        for (int i = 0; i <= 2 * (Dimension_of_Hilbert_Space - 1); i++){

            State_Parameters[i][t] = State_Parameters[i][t - 1] + Step[i];
        }

        State_Generator(State_Parameters[t], Dimension_of_Hilbert_Space, State_Vectors[t]);
    }

    double State_Overlap = 0.0;

    for (int i = 0; i < (sizeof(State_Vectors[0])/sizeof(State_Vectors[0][0])); i++){
        
        State_Overlap += State_Vectors[Number_of_Steps][i] * State_Vectors[0][i];
    }

    double Final_Distance = acos(cabs(State_Overlap));

    return (fabs(PI / 2 - Final_Distance) <= Precision) ? 1 : 0 ;
}

// A simple fucntion for the plotting of the calcualted results. 

void plot_results(const char *filename) {
    FILE *gnuplot = popen("gnuplot -persistent", "w");
    if (!gnuplot) {
        printf("Error: Gnuplot not found.\n");
        return;
    }

    fprintf(gnuplot, "set title 'Critical Step Size vs Number of Steps'\n");
    fprintf(gnuplot, "set xlabel 'Number of Steps'\n");
    fprintf(gnuplot, "set ylabel 'Critical Step Size'\n");
    fprintf(gnuplot, "set grid\n");
    fprintf(gnuplot, "plot '%s' using 1:2 with linespoints title 'Critical Step Size'\n", filename);

    pclose(gnuplot);
}

// Main Calcualtion function which performs the random walk and collates the data (Step size vs. Number of steps required) for plotting and analysis

int main() {
    srand(time(NULL));

    int Number_of_Qubits = 1;
    int Number_of_Walks = 5;

    double Step_Size_Increment = 0.005; 
    double Precision = ((PI/2) - acos(1 / sqrt(pow(2, Number_of_Qubits))));

    const char *data_file = "critical_step_size.dat";
    FILE *fp = fopen(data_file, "w");

    if(!fp) {
        perror("Failed to open file for writing");
        return EXIT_FAILURE;
    }

    printf("Running Random Walk Simulation...\n");

    for(int j = 0; j < Number_of_Walks; j++){
        for (int Number_of_Steps = 3; Number_of_Steps <= 10; Number_of_Steps++){

            double Step_Size;
            double Critical_Step_Size = -1.0;

            for (Step_Size = 0.0; Step_Size <= PI * 0.5; Step_Size += Step_Size_Increment){

                int Result = Random_Walk(Number_of_Qubits, Number_of_Steps, Step_Size, Precision);

                if(Result){
                    break;
                }
                else{}
            }

            if (Critical_Step_Size >= 0.0) {

                fprintf(fp, "%d %f\n", Number_of_Steps, Critical_Step_Size);
            }
        }
    }
    fclose(fp);

    plot_results(data_file);

    return 0;
}