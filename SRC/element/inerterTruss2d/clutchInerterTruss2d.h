/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

#ifndef clutchInerterTruss2d_h
#define clutchInerterTruss2d_h

#include <Element.h>
#include <Matrix.h>

class Node;
class Channel;
class UniaxialMaterial;

class clutchInerterTruss2d : public Element
{
  public:
    // constructors
    clutchInerterTruss2d(int tag, int dim,
        int Nd1, int Nd2,
        double J1, double alpha1, double D1,
        double J2, double alpha2, double D2, int gFlag);
    clutchInerterTruss2d();

    // destructor
    ~clutchInerterTruss2d();

    const char *getClassType(void) const {return "clutchInerterTruss2d";};

    // public methods to obtain information about dof & connectivity
    int getNumExternalNodes(void) const;
    const ID &getExternalNodes(void);
    Node **getNodePtrs(void);

    int getNumDOF(void);
    void setDomain(Domain *theDomain);

    // public methods to set the state of the element
    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);
    int update(void);

    // public methods to obtain stiffness, mass, damping and residual information
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);
    const Matrix &getDamp(void);
    const Matrix &getMass(void);

    //void zeroLoad(void);
    //int addLoad(ElementalLoad *theLoad, double loadFactor);
    //int addInertiaLoadToUnbalance(const Vector &accel);

    const Vector &getResistingForce(void);
    const Vector &getResistingForceIncInertia(void);

    // public methods for element output
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    int displaySelf(Renderer &, int mode, float fact, const char **displayModes=0, int numModes=0);
    void Print(OPS_Stream &s, int flag =0);

    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &eleInformation);

  protected:

  private:
    // private attributes - a copy for each object of the class
    ID  externalNodes;      // contains the tags of the end nodes

    double J1;              //
    double alpha1;          //
    double D1;              //
    double J2;              //
    double alpha2;          //
    double D2;              //
    int gFlag;              // geometric nonlinearity flag

    int state;              //
    double time;            //
    double theta1_dot;      //
    double theta2_dot;      //
    int state_commit;           //
    double time_commit;         //
    double theta1_dot_commit;   //
    double theta2_dot_commit;   //

    double b;               //
    double c;               //

    double Lo;              // initial length
    double Lc;              // current length

    double so;              // initial sin of angle
    double sc;              // current sin of angle
    double co;              // initial cos of angle
    double cc;              // current cos of angle

    double dc;              // current relative displacement
    double vc;              // current relative velocity
    double ac;              // current relative acceleration

    Node *theNodes[2];      // node pointers
    int numDIM;             // number of dimensions for truss
    int numDOF;	            // number of dof for truss

    Matrix *theMatrix;      // pointer to objects matrix (a class wide Matrix)
    Vector *theVector;      // pointer to objects vector (a class wide Vector)

    // static data - single copy for all objects of the class
    static Matrix clutchInerterTruss2dM4;   // class wide matrix for 4*4
    static Matrix clutchInerterTruss2dM6;   // class wide matrix for 6*6
    static Vector clutchInerterTruss2dV4;   // class wide Vector for size 4
    static Vector clutchInerterTruss2dV6;   // class wide Vector for size 6
};

#endif
