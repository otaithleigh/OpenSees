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

#include "clutchInerterTruss2d.h"
#include <elementAPI.h>
#include <G3Globals.h>

#include <Information.h>
#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <Message.h>
#include <FEM_ObjectBroker.h>
#include <UniaxialMaterial.h>
#include <Renderer.h>
#include <ElementResponse.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

// initialise the class wide variables
Matrix clutchInerterTruss2d::clutchInerterTruss2dM4(4,4);
Matrix clutchInerterTruss2d::clutchInerterTruss2dM6(6,6);
Vector clutchInerterTruss2d::clutchInerterTruss2dV4(4);
Vector clutchInerterTruss2d::clutchInerterTruss2dV6(6);

// Documentation
// Truss Element with Inertance
//
// element clutchInerterTruss2d $tag $iNode $jNode $J1 $alpha1 $D1 <$J2 $alpha2 $D2> <-geomNonlinear $gFlag>
//
// Required Input Parameters:
//   $tag                   integer tag identifying the element
//   $iNode, $jNode         end nodes
//   $J1                    inertia of the flywheel active in compression (units: mass times length squared)
//   $alpha1                ratio of rotation of the flywheel active in compression to relative displacement of end nodes (units: rad per length)
//   $D1                    damping coefficient of the flywheel active in compression (units: mass times length squared per second)
//
// Optional Input:
//   $J2, $alpha2, $D2      corresponding values of J, alpha, and D for the flywheel active in extension
//                              the default is for these values to equal those for the flywheel active in compression
//   $gFlag                 geometric nonlinearity flag
//                              gFlag = 0, use geometric linear formulation (default)
//                              gFlag = 1, use geometric nonlinear formulation
//
// Element Notes:
//   To avoid conflict with mass matrix used to develop forces in the ground motion,
//   inertance is defined with a restoring force. Accordingly, a nonlinear solver is
//   necessary when using this element, even for linear problems.
//
// Developed by Mark D. Denavit, University of Tennessee, Knoxville

void * OPS_clutchInerterTruss2d() {

    Element *theTruss = 0;

    int numRemainingArgs = OPS_GetNumRemainingInputArgs();
    if (numRemainingArgs == 0) { // parallel processing
        theTruss = new clutchInerterTruss2d();
        return theTruss;
    }

    if (numRemainingArgs < 6) {
        opserr << "ERROR - clutchInerterTruss2d not enough args provided, want: element clutchInerterTruss2d tag? iNode? jNode? J1? alpha1? D1?\n";
    }

    // get the id and end nodes
    int iData[3];
    double dData[3];
    char *sData  = new char[40];
    int numData;
    int ndm = OPS_GetNDM();


    numData = 3;
    if (OPS_GetIntInput(&numData, iData) != 0) {
        opserr << "WARNING invalid element data\n";
        return 0;
    }
    int eleTag = iData[0];
    int Nd1 = iData[1];
    int Nd2 = iData[2];

    numData = 3;
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
        opserr << "WARNING error reading element area for element" << eleTag << endln;
        return 0;
    }
    double J1 = dData[0];
    double alpha1 = dData[1];
    double D1 = dData[2];

    // Set Default Values for Optional Input
    double J2 = J1;
    double alpha2 = alpha1;
    double D2 = D1;
    int gFlag = 0;

    // If three more doubles are next then they are J2, alpha2, and D2
    numData = 3;
    if (OPS_GetDoubleInput(&numData, dData) == 0) {
        J2 = dData[0];
        alpha2 = dData[1];
        D2 = dData[2];
    }

    // Loop through remaining arguments to get optional input
    while ( OPS_GetNumRemainingInputArgs() > 0 ) {
        if ( OPS_GetStringCopy(&sData) != 0 ) {
            opserr << "WARNING invalid input";
            return 0;
        }

        if ( strcmp(sData,"-geomNonlinear") == 0 ) {
            numData = 1;
            if (OPS_GetInt(&numData, &gFlag) != 0) {
                opserr << "WARNING: Invalid input, want: -geomNonlinear $gFlag in element clutchInerterTruss2d " << eleTag;
                return 0;
            }
        } else {
            opserr << "WARNING unknown option " << sData << "\n";
        }
    }

    // now create the truss and add it to the Domain
    theTruss = new clutchInerterTruss2d(eleTag, ndm, Nd1, Nd2, J1, alpha1, D1, J2, alpha2, D2, gFlag);

    if (theTruss == 0) {
        opserr << "WARNING ran out of memory creating element with tag " << eleTag << endln;
        return 0;
    }

    return theTruss;
}

// constructor:
//  responsible for allocating the necessary space needed by each object
//  and storing the tags of the clutchInerterTruss2d end nodes.

clutchInerterTruss2d::clutchInerterTruss2d(int tag, int dim,
           int Nd1, int Nd2,
           double iJ1, double ialpha1, double iD1,
           double iJ2, double ialpha2, double iD2, int igFlag):
    Element(tag, 0), externalNodes(2),
    numDIM(dim), numDOF(0),
    theMatrix(0), theVector(0),
    J1(iJ1), alpha1(ialpha1), D1(iD1),
    J2(iJ2), alpha2(ialpha2), D2(iD2), gFlag(igFlag),
    state(1), time(0.0), theta1_dot(0.0), theta2_dot(0.0),
    state_commit(1), time_commit(0.0), theta1_dot_commit(0.0), theta2_dot_commit(0.0),
    Lo(0.0), Lc(0.0), so(0.0), sc(0.0), co(0.0), cc(0.0)
{
    // ensure the connectedExternalNode ID is of correct size & set values
    if (externalNodes.Size() != 2) {
      opserr << "FATAL clutchInerterTruss2d::clutchInerterTruss2d - " <<  tag << "failed to create an ID of size 2\n";
      exit(-1);
    }

    externalNodes(0) = Nd1;
    externalNodes(1) = Nd2;

    // set node pointers to NULL
    for (int i=0; i<2; i++)
        theNodes[i] = 0;
}

// constructor:
//   invoked by a FEM_ObjectBroker - blank object that recvSelf needs
//   to be invoked upon
clutchInerterTruss2d::clutchInerterTruss2d():
    Element(0, 0), externalNodes(2),
    numDIM(0), numDOF(0),
    theMatrix(0), theVector(0),
    J1(0.0), alpha1(0.0), D1(0.0),
    J2(0.0), alpha2(0.0), D2(0.0), gFlag(0),
    state(1), time(0.0), theta1_dot(0.0), theta2_dot(0.0),
    state_commit(1), time_commit(0.0), theta1_dot_commit(0.0), theta2_dot_commit(0.0),
    Lo(0.0), Lc(0.0), so(0.0), sc(0.0), co(0.0), cc(0.0)
{
    // ensure the connectedExternalNode ID is of correct size
    if (externalNodes.Size() != 2) {
        opserr << "FATAL clutchInerterTruss2d::clutchInerterTruss2d - failed to create an ID of size 2\n";
        exit(-1);
    }

    for (int i=0; i<2; i++)
        theNodes[i] = 0;
}

//  destructor
//     delete must be invoked on any objects created by the object
//     and on the matertial object.
clutchInerterTruss2d::~clutchInerterTruss2d() {
    // invoke the destructor on any objects created by the object
    // that the object still holds a pointer to
}


int
clutchInerterTruss2d::getNumExternalNodes(void) const {
    return 2;
}

const ID &
clutchInerterTruss2d::getExternalNodes(void) {
    return externalNodes;
}

Node **
clutchInerterTruss2d::getNodePtrs(void) {
  return theNodes;
}

int
clutchInerterTruss2d::getNumDOF(void) {
    return numDOF;
}


// method: setDomain()
//    to set a link to the enclosing Domain and to set the node pointers.
//    also determines the number of dof associated
//    with the clutchInerterTruss2d element, we set matrix and vector pointers,
//    allocate space for t matrix, determine the length
//    and set the transformation matrix.
void
clutchInerterTruss2d::setDomain(Domain *theDomain)
{
    // check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
        theNodes[0] = 0;
        theNodes[1] = 0;
        Lo = 0;
        Lc = 0;
        return;
    }

    // first set the node pointers
    int Nd1 = externalNodes(0);
    int Nd2 = externalNodes(1);
    theNodes[0] = theDomain->getNode(Nd1);
    theNodes[1] = theDomain->getNode(Nd2);

    // if can't find both - send a warning message
    if ((theNodes[0] == 0) || (theNodes[1] == 0)) {
        if (theNodes[0] == 0)
            opserr <<"clutchInerterTruss2d::setDomain() - clutchInerterTruss2d" << this->getTag() << " node " << Nd1 <<
                "does not exist in the model\n";
        else
            opserr <<"clutchInerterTruss2d::setDomain() - clutchInerterTruss2d" << this->getTag() << " node " << Nd2 <<
                "does not exist in the model\n";

        // fill this in so don't segment fault later
        numDOF = 4;
        theMatrix = &clutchInerterTruss2dM4;
        theVector = &clutchInerterTruss2dV4;

        return;
    }

    // now determine the number of dof and the dimesnion
    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();

    // if differing dof at the ends - print a warning message
    if (dofNd1 != dofNd2) {
        opserr <<"WARNING clutchInerterTruss2d::setDomain(): nodes " << Nd1 << " and " << Nd2 <<
            "have differing dof at ends for clutchInerterTruss2d " << this->getTag() << endln;

        // fill this in so don't segment fault later
        numDOF = 4;
        theMatrix = &clutchInerterTruss2dM4;
        theVector = &clutchInerterTruss2dV4;

        return;
    }

    // call the base class method
    this->DomainComponent::setDomain(theDomain);

    // now set the number of dof for element and set matrix and vector pointer
    if (numDIM == 2 && dofNd1 == 2) {
        numDOF = 4;
        theMatrix = &clutchInerterTruss2dM4;
        theVector = &clutchInerterTruss2dV4;
    } else if (numDIM == 2 && dofNd1 == 3) {
        numDOF = 6;
        theMatrix = &clutchInerterTruss2dM6;
        theVector = &clutchInerterTruss2dV6;
    } else {
        opserr <<"WARNING clutchInerterTruss2d::setDomain cannot handle " << numDIM << " dofs at nodes in " <<
            dofNd1  << " problem\n";

        numDOF = 4;
        theMatrix = &clutchInerterTruss2dM4;
        theVector = &clutchInerterTruss2dV4;
        return;
    }


    // now determine the length, cosines and fill in the transformation
    // NOTE t = -t(every one else uses for residual calc)
    const Vector &end1Crd = theNodes[0]->getCrds();
    const Vector &end2Crd = theNodes[1]->getCrds();

    // Set undeformed and initial length
    double Lx = end2Crd(0)-end1Crd(0);
    double Ly = end2Crd(1)-end1Crd(1);
    Lo = sqrt(Lx*Lx + Ly*Ly);

    if (Lo == 0.0) {
        opserr <<"WARNING clutchInerterTruss2d::setDomain() - clutchInerterTruss2d " << this->getTag() << " has zero length\n";
        return;
    }

    co = Lx/Lo;
    so = Ly/Lo;

    // Set current to initial
    Lc = Lo;
    cc = co;
    sc = so;
    b = 0.0;
    c = 0.0;

    return;
}


int
clutchInerterTruss2d::commitState() {
    int err;

    // call element commitState to do any base class stuff
    if ((err = this->Element::commitState()) != 0) {
        opserr << "clutchInerterTruss2d::commitState () - failed in base class";
        return err;
    }

    // commit the element state variables
    state_commit = state;
    time_commit = time;
    theta1_dot_commit = theta1_dot;
    theta2_dot_commit = theta2_dot;

    return 0;
}

int
clutchInerterTruss2d::revertToLastCommit() {
    state = state_commit;
    time = time_commit;
    theta1_dot = theta1_dot_commit;
    theta2_dot = theta2_dot_commit;

    return 0;
}

int
clutchInerterTruss2d::revertToStart() {
    state = 0;
    time = 0.0;
    time_commit = 0.0;
    theta1_dot = 0.0;
    theta2_dot = 0.0;
    theta1_dot_commit = 0.0;
    theta2_dot_commit = 0.0;
    dc = 0.0;
    vc = 0.0;
    ac = 0.0;
    b = 0.0;
    c = 0.0;

    return 0;
}

int
clutchInerterTruss2d::update(void) {
    // Nodal displacements, velocities, and accelerations
    const Vector &end1Disp  = theNodes[0]->getTrialDisp();
    const Vector &end2Disp  = theNodes[1]->getTrialDisp();
    const Vector &end1Vel   = theNodes[0]->getTrialVel();
    const Vector &end2Vel   = theNodes[1]->getTrialVel();
    const Vector &end1Accel = theNodes[0]->getTrialAccel();
    const Vector &end2Accel = theNodes[1]->getTrialAccel();


    if (gFlag == 0) { // geometric linear
        dc =   (end2Disp(0)-end1Disp(0))*co +   (end2Disp(1)-end1Disp(1))*so;
        vc =     (end2Vel(0)-end1Vel(0))*co +     (end2Vel(1)-end1Vel(1))*so;
        ac = (end2Accel(0)-end1Accel(0))*co + (end2Accel(1)-end1Accel(1))*so;

    } else { // geometric nonlinear
        const Vector &end1Crd   = theNodes[0]->getCrds();
        const Vector &end2Crd   = theNodes[1]->getCrds();
        double Lx = (end2Crd(0)+end2Disp(0))-(end1Crd(0)+end1Disp(0));
        double Ly = (end2Crd(1)+end2Disp(1))-(end1Crd(1)+end1Disp(1));
        Lc = sqrt(Lx*Lx+Ly*Ly);
        cc = Lx/Lc;
        sc = Ly/Lc;
        dc = Lc-Lo;
        vc =     (end2Vel(0)-end1Vel(0))*cc +     (end2Vel(1)-end1Vel(1))*sc;
        ac = (end2Accel(0)-end1Accel(0))*cc + (end2Accel(1)-end1Accel(1))*sc;
    }

    // state determination
    time = (this->getDomain())->getCurrentTime();
    if (alpha1*vc < theta1_dot_commit*exp(-D1*(time-time_commit)/J1)) {
        state = 1;
        theta1_dot = alpha1*vc;
        theta2_dot = theta2_dot_commit*exp(-D2*(time-time_commit)/J2);
        b = J1*alpha1*alpha1;
        c = D1*alpha1*alpha1;
    } else if (alpha2*vc > theta2_dot_commit*exp(-D2*(time-time_commit)/J2)) {
        state = 2;
        theta1_dot = theta1_dot_commit*exp(-D1*(time-time_commit)/J1);
        theta2_dot = alpha2*vc;
        b = J2*alpha2*alpha2;
        c = D2*alpha2*alpha2;
    } else {
        state = 3;
        theta1_dot = theta1_dot_commit*exp(-D1*(time-time_commit)/J1);
        theta2_dot = theta2_dot_commit*exp(-D2*(time-time_commit)/J2);
        b = 0.0;
        c = 0.0;
    }

    if ( isnan(theta1_dot) || isnan(theta2_dot) ) {
        return -1;
    }

    return 0;
}


const Matrix &
clutchInerterTruss2d::getTangentStiff(void) {
    Matrix &K = *theMatrix;
    K.Zero();
    return K;
}


const Matrix &
clutchInerterTruss2d::getInitialStiff(void) {
    Matrix &K = *theMatrix;
    K.Zero();
    return K;
}

const Matrix &
clutchInerterTruss2d::getDamp(void) {
    Matrix &C = *theMatrix;
    C.Zero();

    if (Lo == 0.0) { // - problem in setDomain() no further warnings
        return C;
    }

    // add in defined damping
    if (c != 0.0) {
        // transformation matrix
        static Matrix T(4,4);
        T.Zero();
        double cosq, sinq;
        if (gFlag == 0) {
            cosq = co; // geometric linear
            sinq = so;
        } else {
            cosq = cc; // geometric nonlinear
            sinq = sc;
        }
        T(0,0) =  cosq;
        T(1,0) = -sinq;
        T(0,1) =  sinq;
        T(1,1) =  cosq;
        T(2,2) =  cosq;
        T(3,2) = -sinq;
        T(2,3) =  sinq;
        T(3,3) =  cosq;

        static Matrix cb(4,4);
        cb.Zero();
        cb(0,0) =  c;
        cb(0,2) = -c;
        cb(2,0) = -c;
        cb(2,2) =  c;

        static Matrix ce(4,4);
        ce.Zero();
        ce.addMatrixTripleProduct(0.0, T, cb, 1.0);

        int numDOF2 = numDOF/2;
        for (int i = 0; i < numDIM; i++) {
            for (int j = 0; j < numDIM; j++) {
                C(i,j)                 += ce(i,j);
                C(i+numDOF2,j)         += ce(i+numDIM,j);
                C(i,j+numDOF2)         += ce(i,j+numDIM);
                C(i+numDOF2,j+numDOF2) += ce(i+numDIM,j+numDIM);
            }
        }
    }

    return C;
}


const Matrix &
clutchInerterTruss2d::getMass(void) {
    // zero the matrix
    Matrix &M = *theMatrix;
    M.Zero();

    // check for quick return
    if (Lo == 0.0) { // - problem in setDomain() no further warnings
        return M;
    }

    if (b != 0.0) {
        // transformation matrix
        static Matrix T(4,4);
        T.Zero();
        double cosq, sinq;
        if (gFlag == 0) {
            cosq = co; // geometric linear
            sinq = so;
        } else {
            cosq = cc; // geometric nonlinear
            sinq = sc;
        }
        T(0,0) =  cosq;
        T(1,0) = -sinq;
        T(0,1) =  sinq;
        T(1,1) =  cosq;
        T(2,2) =  cosq;
        T(3,2) = -sinq;
        T(2,3) =  sinq;
        T(3,3) =  cosq;

        static Matrix bb(4,4);
        bb.Zero();
        bb(0,0) =  b;
        bb(0,2) = -b;
        bb(2,0) = -b;
        bb(2,2) =  b;

        static Matrix be(4,4);
        be.addMatrixTripleProduct(0.0, T, bb, 1.0);

        int numDOF2 = numDOF/2;
        for (int i = 0; i < numDIM; i++) {
            for (int j = 0; j < numDIM; j++) {
                M(i,j)                 += be(i,j);
                M(i+numDOF2,j)         += be(i+numDIM,j);
                M(i,j+numDOF2)         += be(i,j+numDIM);
                M(i+numDOF2,j+numDOF2) += be(i+numDIM,j+numDIM);
            }
        }
    }

    return M;
}


const Vector &
clutchInerterTruss2d::getResistingForce() {
    Vector &P = *theVector;
    P.Zero();

    if (Lo == 0.0) { // - problem in setDomain() no further warnings
        return P;
    }

    double force = b*ac + c*vc;

    double cosq, sinq;
    if (gFlag == 0) {
        cosq = co; // geometric linear
        sinq = so;
    } else {
        cosq = cc; // geometric nonlinear
        sinq = sc;
    }

    int numDOF2 = numDOF/2;
    P(0)         = -force*cosq;
    P(1)         = -force*sinq;
    P(numDOF2+0) =  force*cosq;
    P(numDOF2+1) =  force*sinq;

    return P;
}


const Vector &
clutchInerterTruss2d::getResistingForceIncInertia() {
    Vector &P = *theVector;
    P = this->getResistingForce();
    return P;
}

int
clutchInerterTruss2d::sendSelf(int commitTag, Channel &theChannel) {
    opserr << "Error: clutchInerterTruss2d::sendSelf -- not yet implemented\n";
    return -1;
}

int
clutchInerterTruss2d::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker) {
    opserr << "Error: clutchInerterTruss2d::recvSelf -- not yet implemented\n";
    return -1;
}

int
clutchInerterTruss2d::displaySelf(Renderer &theViewer, int displayMode, float fact,
		   const char **displayModes, int numModes) {
    opserr << "Error: clutchInerterTruss2d::displaySelf -- not yet implemented\n";
    return -1;
}

void
clutchInerterTruss2d::Print(OPS_Stream &s, int flag) {
    s << "Element: " << this->getTag() << "\n";
    s << " type: clutchInerterTruss2d\n";
    s << " iNode: " << externalNodes(0) << "\n";
    s << " jNode: " << externalNodes(1) << "\n";
    s << " J1: " << J1 << "\n";
    s << " alpha1: " << alpha1 << "\n";
    s << " D1: " << D1 << "\n";
    s << " J2: " << J2 << "\n";
    s << " alpha2: " << alpha2 << "\n";
    s << " D2: " << D2 << "\n";
}

Response*
clutchInerterTruss2d::setResponse(const char **argv, int argc, OPS_Stream &output) {
    Response *theResponse = 0;

    output.tag("ElementOutput");
    output.attr("eleType","clutchInerterTruss2d");
    output.attr("eleTag",this->getTag());
    output.attr("node1",externalNodes[0]);
    output.attr("node2",externalNodes[1]);

    //
    // we compare argv[0] for known response types for the clutchInerterTruss2d
    //

    if (strcmp(argv[0],"force") == 0 ||
        strcmp(argv[0],"forces") == 0 ||
        strcmp(argv[0],"globalForce") == 0 ||
        strcmp(argv[0],"globalForces") == 0) {

        char outputData[10];
        int numDOFperNode = numDOF/2;
        for (int i=0; i<numDOFperNode; i++) {
            sprintf(outputData,"P1_%d", i+1);
            output.tag("ResponseType", outputData);
        }
        for (int j=0; j<numDOFperNode; j++) {
            sprintf(outputData,"P2_%d", j+1);
            output.tag("ResponseType", outputData);
        }
        theResponse =  new ElementResponse(this, 1, Vector(numDOF));

    } else if (strcmp(argv[0],"axialForce") == 0 ||
               strcmp(argv[0],"basicForce") == 0 ||
               strcmp(argv[0],"localForce") == 0 ||
               strcmp(argv[0],"basicForces") == 0) {

        output.tag("ResponseType", "N");
        theResponse =  new ElementResponse(this, 2, Vector(1));

    } else if (strcmp(argv[0],"defo") == 0 ||
               strcmp(argv[0],"deformation") == 0 ||
               strcmp(argv[0],"deformations") == 0 ||
               strcmp(argv[0],"basicDefo") == 0 ||
               strcmp(argv[0],"basicDeformation") == 0 ||
               strcmp(argv[0],"basicDeformations") == 0) {

        output.tag("ResponseType", "U");
        theResponse = new ElementResponse(this, 3, Vector(1));

    } else if (strcmp(argv[0],"state") == 0)  {

        theResponse = new ElementResponse(this, 4, Vector(1));

    } else if (strcmp(argv[0],"flywheelvelocity") == 0)  {

        theResponse = new ElementResponse(this, 5, Vector(2));

    }

    output.endTag();

    return theResponse;
}

int
clutchInerterTruss2d::getResponse(int responseID, Information &eleInfo) {
    static Vector fVec(1);
    static Vector fVec2(2);
    static Matrix kVec(1,1);

    switch (responseID) {
        case 1:
            return eleInfo.setVector(this->getResistingForce());

        case 2:
            fVec(0) = b*ac + c*vc;
            return eleInfo.setVector(fVec);

        case 3:
            fVec(0) = dc;
            return eleInfo.setVector(fVec);

        case 4:
            fVec(0) = state;
            return eleInfo.setVector(fVec);

       case 5:
            fVec2(0) = theta1_dot;
            fVec2(1) = theta2_dot;
            return eleInfo.setVector(fVec2);

        default:
            return 0;
    }
}
