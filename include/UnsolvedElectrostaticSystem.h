/**
 * A class to represent an unsolved electrostatic system.
 *
 * Inherits from ElectrostaticSystem, adding a matrix (from eigen, dense) to
 * specify which points of the grid are boundary conditions and which are 
 * unspecified points to calculate the potential at.
 * 
 * eg:
 *
 * potentials:      and         boundaryCondition:
 * 6 7 8                    1 1 1
 * 3 4 5                    0 1 0
 * 0 1 2                    0 0 1
 *
 * would mean that positions 2, 4, 6, 7 and 8 of the potentials are boundary 
 * conditions, while the other positions are unknown points.
 */

#ifndef UNSOLVEDELECTROSTATICSYSTEM_H
#define UNSOLVEDELECTROSTATICSYSTEM_H

#include <Eigen/Dense>
#include "ElectrostaticSystem.h"

namespace electrostatics {

// Matrix to represent a grid of boolean values    
typedef Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> boolGrid; 

class UnsolvedElectrostaticSystem : public ElectrostaticSystem{
    protected:
        /* Boolean grid to mark the positions of the boundary conditions. */
        boolGrid boundaryConditionPositions; 

    public:
        /* Constructor */
        UnsolvedElectrostaticSystem(int iMin, int iMax, int jMin, int jMax);
        

        /* Methods */

        /* Test if position (i, j) or (k) is a boundary condition. */
        bool isBoundaryConditionIJ(int i, int j) const;
        bool isBoundaryConditionK(long k) const;

        /* Set/unset position (i, j) or (k) as a boundary condition. */
        void setBoundaryConditionIJ(int i, int j, bool isBoundaryCondition);
        void setBoundaryConditionK(long k, bool isBoundaryCondition);

        /* Set a point as a boundary condition with specified potential. */
        void setBoundaryPoint(int i, int j, double potential);

        /* Set ring as boundary condition with specified potential.
         * Note: This creates a ring only.
         * For a filled circle use setBoundaryCircle()
         */
        void setBoundaryRing(int centreI, int centreJ, double radius, double potential);

        /* Set circle as boundary condition with specified potential.
         * Note: This creates a filled circle.
         * For a ring only use setBoundaryRing()
         */
        void setBoundaryCircle(int centreI, int centreJ, double radius, double potential);

        /* Set a line as a boundary condition - from (i1,j1) to (i2,j2). */
        void setBoundaryLine(int i1, int j1, int i2, int j2, double potential);

        /* Set the edges of the system as a boundary conditions. */
        void setLeftBoundary(double potential);
        void setRightBoundary(double potential);
        void setTopBoundary(double potential);
        void setBottomBoundary(double potential);

        /* Set a rectangle as a boundary condition.
         * Top left point: (left, top).
         * Bottom left point: (left, bottom).
         * Top right point: (right, top).
         * Bottom right point: (right, bottom).
         * NOTE: Rectange is not filled!
         */
        void setBoundaryRectangle(int left, int right, int top, int bottom, double potential);
};

/* Set a quadrilateral as a boundary condition. 
 * With vertices defined by (x1,y1),(x2,y2),(x3,y3),(x4,y4).
 * Reccomend sketching quadrilateral on paper to help define correct points for desired shape
 * NOTE: Quadrilateral is not filled! 
 */
 void setBoundaryQuad(int x1, int x2, int x3, int x4, int y1, int y2, int y3, int y4, double potential);

} // namespace electrostatics
#endif
