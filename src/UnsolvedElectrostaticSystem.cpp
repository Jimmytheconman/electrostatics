#include <Eigen/Dense>
#include <stdexcept>
#include <cmath>
#include "ElectrostaticSystem.h"
#include "UnsolvedElectrostaticSystem.h"

namespace electrostatics {

/* Constructors */

UnsolvedElectrostaticSystem::UnsolvedElectrostaticSystem(int iMin, int iMax, int jMin, int jMax) :
    ElectrostaticSystem(iMin, iMax, jMin, jMax) {
        boundaryConditionPositions = boolGrid(potentials.rows(), potentials.cols());
        boundaryConditionPositions.fill(false);
}


/* Methods */

bool UnsolvedElectrostaticSystem::isBoundaryConditionIJ(int i, int j) const {
    if(i>iMax || i<iMin || j>jMax || j<jMin) throw std::out_of_range("Error: Trying to get element out of range!");
    return boundaryConditionPositions(i-iMin, j-jMin);
}
bool UnsolvedElectrostaticSystem::isBoundaryConditionK(long k) const {
    if(k>kMax || k<0) throw std::out_of_range("Error: Trying to get element out of range!");
    int* ij = k2ij(k);
    return boundaryConditionPositions(ij[0]-iMin, ij[1]-jMin);
}

void UnsolvedElectrostaticSystem::setBoundaryConditionIJ(int i, int j, bool isBoundaryCondition) {
    if(i>iMax || i<iMin || j>jMax || j<jMin) throw std::out_of_range("Error: Trying to set element out of range!");
    boundaryConditionPositions(i-iMin, j-jMin) = isBoundaryCondition;
}
void UnsolvedElectrostaticSystem::setBoundaryConditionK(long k, bool isBoundaryCondition) {
    if(k>kMax || k<0) throw std::out_of_range("Error: Trying to set element out of range!");
    int* ij = k2ij(k);
    boundaryConditionPositions(ij[0]-iMin, ij[1]-jMin) = isBoundaryCondition;
}

void UnsolvedElectrostaticSystem::setBoundaryPoint(int i, int j, double potential) {
    if(i>iMax || i<iMin || j>jMax || j<jMin) throw std::out_of_range("Error: Trying to set element out of range!");
    setBoundaryConditionIJ(i, j, true);
    setPotentialIJ(i, j, potential);
}

void UnsolvedElectrostaticSystem::setBoundaryRing(int centreI, int centreJ, double radius, double potential) {
    int iOffset, jOffset;
    // i offsets from the centre of the ring
    for(iOffset=ceil(radius); iOffset>=0; iOffset--) {
        // Will only happen for the most extreme i values
        if(iOffset > radius) jOffset = 0;
        else jOffset = round(sqrt( pow(radius,2)-pow(iOffset,2) ));
        // Try / catch is not the proper way to do this, but it is by far the most convenient
        try{
            setBoundaryPoint(centreI+iOffset, centreJ+jOffset, potential);
            setBoundaryPoint(centreI+iOffset, centreJ-jOffset, potential);
            setBoundaryPoint(centreI-iOffset, centreJ+jOffset, potential);
            setBoundaryPoint(centreI-iOffset, centreJ-jOffset, potential);
        } catch (const std::out_of_range& e) {}
    }
    // j offsets from the centre of the ring
    for(jOffset=ceil(radius); jOffset>=0; jOffset--) {
        // Will only happen for the most extreme j values
        if(jOffset > radius) iOffset = 0;
        else iOffset = round(sqrt( std::pow(radius,2)-pow(jOffset,2) ));
        // Try / catch is not the proper way to do this, but it is by far the most convenient
        try{
            setBoundaryPoint(centreI+iOffset, centreJ+jOffset, potential);
            setBoundaryPoint(centreI+iOffset, centreJ-jOffset, potential);
            setBoundaryPoint(centreI-iOffset, centreJ+jOffset, potential);
            setBoundaryPoint(centreI-iOffset, centreJ-jOffset, potential);
        } catch (const std::out_of_range& e) {}
    }
}

void UnsolvedElectrostaticSystem::setBoundaryCircle(int centreI, int centreJ, double radius, double potential) {
    for(int i=centreI-ceil(radius); i<=centreI+ceil(radius); i++) {
        for(int j=centreJ-ceil(radius); j<=centreJ+ceil(radius); j++) {
            // If point is inside the circle and on the grid, set it to potential
            if((sqrt( pow(i-centreI,2)+pow(j-centreJ,2) ) <= radius) && 
                    i>=iMin && i<=iMax && j>=jMin && j<=jMax) {
                setBoundaryPoint(i, j, potential);
            }
        }
    }
}

void UnsolvedElectrostaticSystem::setBoundaryLine(int i1, int j1, int i2, int j2, double potential) {
    if(i1>iMax || i1<iMin || j1>jMax || j1<jMin || i2>iMax || i2<iMin || j2>jMax || j2<jMin)  { 
        throw std::out_of_range("Error: Trying to set elements out of range!");
    }
    // Constant x value
    if((i2-i1) == 0) {
        // Swap so j is increasing
        if(j2<j1) {
            int jTemp = j1;
            j1 = j2;
            j2 = jTemp;
        }
        for(int j=j1; j<=j2; j++) setBoundaryPoint(i1, j, potential);
    }
    // Constant y value
    else if((j2-j1) == 0) {
        // Swap so i is increasing
        if(i2<i1) {
            int iTemp = i1;
            i1 = i2;
            i2 = iTemp;
        }
        for(int i=i1; i<=i2; i++) setBoundaryPoint(i, j1, potential);
    }
    // Any other gradient
    else {
        // Swap points so i is increasing
        if(i2<i1) {
            int temp = i1;
            i1 = i2;
            i2 = temp;
            temp = j1;
            j1 = j2;
            j2 = temp;
        }
        double x = i1;
        double y = j1;
        double gradient = (j2-j1)/(double)(i2-i1);
        while(x<=i2) {
            setBoundaryPoint(round(x), round(y), potential);
            x += (gradient<=1)?(1):(1/gradient);
            y += (gradient<=1)?(gradient):(1);
        }
    }
}

void UnsolvedElectrostaticSystem::setLeftBoundary(double potential) { 
    setBoundaryLine(iMin, jMin, iMin, jMax, potential);
}
void UnsolvedElectrostaticSystem::setRightBoundary(double potential) { 
    setBoundaryLine(iMax, jMin, iMax, jMax, potential);
}
void UnsolvedElectrostaticSystem::setTopBoundary(double potential) { 
    setBoundaryLine(iMin, jMax, iMax, jMax, potential);
}
void UnsolvedElectrostaticSystem::setBottomBoundary(double potential) { 
    setBoundaryLine(iMin, jMin, iMax, jMin, potential);
}

void UnsolvedElectrostaticSystem::setBoundaryRectangle(int left, int right,
        int top, int bottom, double potential) {
    setBoundaryLine(left, top, right, top, potential);
    setBoundaryLine(left, bottom, right, bottom, potential);
    setBoundaryLine(left, bottom, left, top, potential);
    setBoundaryLine(right, bottom, right, top, potential);
}

  void UnsolvedElectrostaticSystem::setBoundaryQuad(int x1, int x2, int x3, int x4, int y1, int y2, int y3, int y4, double potential) {
    setBoundaryLine(x1,y1,x2,y2);
    setBoundaryLine(x2,y2,x3,y3);
    setBoundaryLine(x3,y3,x4,y4);
    setBoundaryLine(x4,y4,x1,y1);
  }


} // namespace electrostatics
