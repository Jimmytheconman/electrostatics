#include <Eigen/Dense>
#include <stdexcept>
#include "ElectrostaticSystem.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <list>
#include <vector>
#include "SolvedElectrostaticSystem.h"

namespace electrostatics {

/* Constructors */

SolvedElectrostaticSystem::SolvedElectrostaticSystem(int iMin, int iMax, int jMin, int jMax) :
    ElectrostaticSystem(iMin, iMax, jMin, jMax) {
        fieldFound = false;
        maxField = 0;
}


/* Methods */

void SolvedElectrostaticSystem::findField() {
    fieldX = doubleGrid::Zero(iMax-iMin+1, jMax-jMin+1);    // X components of field
    fieldY = doubleGrid::Zero(iMax-iMin+1, jMax-jMin+1);    // Y components of field
    field = doubleGrid::Zero(iMax-iMin+1, jMax-jMin+1);     // Magnitude of field
    for(int i=iMin; i<=iMax; i++) {
        for(int j=jMin; j<=jMax; j++) {
            // Components in i direction
            if(i==iMax) fieldX(i-iMin, j-jMin) = -(getPotentialIJ(i, j) - getPotentialIJ(i-1, j));
            else if(i==iMin) fieldX(i-iMin, j-jMin) = -(getPotentialIJ(i+1, j)  - getPotentialIJ(i, j));
            else fieldX(i-iMin, j-jMin) = -((getPotentialIJ(i+1, j) - getPotentialIJ(i-1, j))/2);

            // Components in j direction
            if(j==jMax) fieldY(i-iMin, j-jMin) = -(getPotentialIJ(i, j) - getPotentialIJ(i, j-1));
            else if(j==jMin) fieldY(i-iMin, j-jMin) = -(getPotentialIJ(i, j+1)  - getPotentialIJ(i, j));
            else fieldY(i-iMin, j-jMin) = -((getPotentialIJ(i, j+1) - getPotentialIJ(i, j-1))/2);

            // Full field
            double currentField = sqrt( pow(fieldX(i-iMin, j-jMin), 2) +
                    pow(fieldY(i-iMin, j-jMin), 2) );
            field(i-iMin, j-jMin) = currentField;

            if(currentField > maxField) maxField = currentField;
        }
    }
    fieldFound = true;
}

void SolvedElectrostaticSystem::saveFieldGNUPlot(std::string fileName) {
    if(!fieldFound) findField();
    std::ofstream outputFile;
    outputFile.open(fileName.c_str());
    for(int i=iMin; i<=iMax; i++) {
        for(int j=jMin; j<=jMax; j++) {
            outputFile << i << " " << j << " " << fieldX(i-iMin, j-jMin)/maxField << " " << fieldY(i-iMin, j-jMin)/maxField <<"\n";
        }
        outputFile << "\n";
    }
    outputFile.close();
}

} // namespace electrostatics
