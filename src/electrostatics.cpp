#include "finiteDifferenceSolve.h"
#include "analyticalSolutions.h"
#include "SolvedElectrostaticSystem.h"
#include "UnsolvedElectrostaticSystem.h"
#include <iostream>
#include <cmath>

//double problem1analytical();
//double problem2analytical();
void solveProblem1();
void solveProblem2();

/* Currently asks the user to select problem 1 or 2 and to input the relevant paramters.
 * The output (plot data for the analytical solution, plot data for the numerical solution, 
 * and plot data for the difference between the two solutions) is saved to the specified files.
 *
 * It will need to be generalised / changed to allow solving of arbitrary electrostatic
 * systems.
 */
int main() {
    int problemNumber;

    std::cout << "\nWhich system would you like to generate plot data for?\n" <<
        "1. The system with the two concentical cylinders.\n" <<
        "2. The system with the cylinder between two plates.\n" <<
        "Enter '1' or '2': ";
    std::cin >> problemNumber;
    std::cout << "\n";

    switch(problemNumber) {
        case 1:
            solveProblem1();
            break;
        case 2:
            solveProblem2();
            break;
    }
    std::cout << "\n";
    return 0;
}


void solveProblem1() {
    std::cout << "Solving system 1 using grid size ixj.\n";

    // Get the dimensions of the problem
    int iMin;
    int iMax;
    int jMin;
    int jMax;
    std::cout <<"Enter minimum value for i: ";
    std::cin >> iMin;
    std::cout <<"Enter maximum value for i: ";
    std::cin >> iMax;
    std::cout <<"Enter minimum value for j: ";
    std::cin >> jMin;
    std::cout << "Enter maximum value for j: ";
    std::cin >> jMax;

    // Get the necessary parameters - inner cylinder is A, outer is B
    double radiusA;
    double potentialA;
    double radiusB;
    double potentialB;
    std::cout << "Enter radius for inner cylinder: ";
    std::cin >> radiusA;
    std::cout << "Enter potential for inner cylinder: ";
    std::cin >> potentialA; 
    std::cout << "Enter radius for outer cylinder: ";
    std::cin >> radiusB;
    std::cout << "Enter potential for outer cylinder: ";
    std::cin >> potentialB;


    // Numerical Solution
    // Setup the unsolved system with the boundary conditions, and then solve it
    electrostatics::UnsolvedElectrostaticSystem unsolvedSystem(iMin, iMax, jMin, jMax);
    unsolvedSystem.setBoundaryCircle(0, 0, radiusA, potentialA);
    unsolvedSystem.setBoundaryRing(0, 0, radiusB, potentialB);

    electrostatics::SolvedElectrostaticSystem solvedSystemNumerical(iMin, iMax, jMin, jMax);
    electrostatics::finiteDifferenceSolve(unsolvedSystem, solvedSystemNumerical);
    solvedSystemNumerical.saveFile("numericalProblem1");
    solvedSystemNumerical.saveFieldGNUPlot("numericalFieldProblem1");


    // Analytical solution
    electrostatics::SolvedElectrostaticSystem systemAnalytical(iMin, iMax, jMin, jMax);
    for(int i=iMin; i<=iMax; i++) {
        for(int j=jMin; j<=jMax; j++) {
            systemAnalytical.setPotentialIJ(i, j,
                    electrostatics::analyticalProblem1(i, j, radiusA, radiusB, potentialA, potentialB));
        }
    }
    systemAnalytical.saveFile("analyticalProblem1");


    // Difference between analytical and numerical solutions
    electrostatics::SolvedElectrostaticSystem solutionComparison(iMin, iMax, jMin, jMax);
    solvedSystemNumerical.compareTo(systemAnalytical, solutionComparison);
    solutionComparison.saveFile("differenceProblem1");
}


void solveProblem2() {
    // Get the dimensions of the problem
    int iMin;
    int iMax;
    int jMin;
    int jMax;
    std::cout <<"Enter minimum value for i: ";
    std::cin >> iMin;
    std::cout <<"Enter maximum value for i: ";
    std::cin >> iMax;
    std::cout <<"Enter minimum value for j: ";
    std::cin >> jMin;
    std::cout << "Enter maximum value for j: ";
    std::cin >> jMax;

    // Get the necessary parameters
    double cylinderRadius;
    double cylinderPotential;
    double leftPotential;
    double rightPotential;
    std::cout << "Enter radius for cylinder: ";
    std::cin >> cylinderRadius;
    cylinderPotential = 0;
    std::cout << "Enter potential for left plate: ";
    std::cin >> leftPotential;
    std::cout << "Enter potential for right plate: ";
    std::cin >> rightPotential;


    // Numerical Solution
    // Setup the unsolved system with the boundary conditions, and then solve it
    electrostatics::UnsolvedElectrostaticSystem unsolvedSystem(iMin, iMax, jMin, jMax);
    unsolvedSystem.setBoundaryCircle(0, 0, cylinderRadius, cylinderPotential);
    unsolvedSystem.setLeftBoundary(leftPotential);
    unsolvedSystem.setRightBoundary(rightPotential);

    electrostatics::SolvedElectrostaticSystem solvedSystemNumerical(iMin, iMax, jMin, jMax);
    electrostatics::finiteDifferenceSolve(unsolvedSystem, solvedSystemNumerical);
    solvedSystemNumerical.saveFile("numericalProblem2");
    solvedSystemNumerical.saveFieldGNUPlot("numericalFieldProblem2");

    // Analytical solution
    electrostatics::SolvedElectrostaticSystem systemAnalytical(iMin, iMax, jMin, jMax);
    // Calculate the potential at each point using the analytical solution
    double uniformField =electrostatics::uniformField(iMin, iMax, leftPotential, rightPotential);
    for(int i=iMin; i<=iMax; i++) {
        for(int j=jMin; j<=jMax; j++) {
            systemAnalytical.setPotentialIJ(i, j, 
                    electrostatics::analyticalProblem2(i, j, cylinderRadius, uniformField));
        }
    }
    systemAnalytical.saveFile("analyticalProblem2");


    // Difference between analytical and numerical solutions
    electrostatics::SolvedElectrostaticSystem solutionComparison(iMin, iMax, jMin, jMax);
    solvedSystemNumerical.compareTo(systemAnalytical, solutionComparison);
    solutionComparison.saveFile("differenceProblem2");
}
