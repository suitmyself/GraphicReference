/********************************************************************
*   RigidBody.cpp
*   Richard Chaney
*   Version 4.0 - June 2000
*
*   Class for a 3D rigid body with 4th order Runge-Kutta integrator.
*
*   To set up a rigid body do as follows :
*       1.  Call setGravity() to set g (gravity is always along the z axis)
*       2.  Call setBox(), setCylinder() or setSphere() to create a rigid body primitive
*       3.  Use combine() to join primitives together if needed
*       4.  Call diagonalise() to make moment of inertia matrix diagonal. It also initialise various bits.
*       5.  Don't forget to transform your graphical object and collision detection object by the transform :
*               newCoord = diagRotPos * oldCoord + diagLinPos
*       6.  Call initialise() with initial parameters
*
*   Bodies can be integrated with Euler or 4th order Runge-Kutta. Runge-Kutta is recommended as it is
*   much more accurate for a given step size, however the step takes 4 times as long to calculate.
*
*   To update it with Euler method do :
*       1.  Call setStepSize() to set size of step (applies to all rigid bodies)
*       2.  Call RKInit()
*       3.  Use linPos, linVel, rotPos and rotVel to determine forces on body.
*           Add forces and torques with provided functions, or add them to linForce and rotForce directly.
*       4.  Calculate constraint forces if necessary and add to linForce and rotForce
*       5.  Call eulerStep()
*       6.  Draw body at linPos, linVel, rotPos and rotVel.
*
*   To update it with Runge-Kutta method do :
*       1.  Call setStepSize() to set size of step (applies to all rigid bodies)
*       2.  Call RKInit()
*       3.  Use linPos, linVel, rotPos and rotVel to determine forces on body.
*           Add forces and torques with provided functions, or add them to linForce and rotForce directly.
*       4.  Calculate constraint forces if necessary and add to linForce and rotForce
*       5.  Call RKStep1()
*       6.  Repeat step 2 and 3.
*       7.  Call RKStep2()
*       8.  Repeat step 2 and 3.
*       9.  Call RKStep3()
*       10. Repeat step 2 and 3.
*       11. Call RKStep4()
*       12. Draw body at linPos, linVel, rotPos and rotVel.
*
*   To go back to the previous position after a step is made call undo(). This is
*   useful when a collision is detected.
********************************************************************/

#include <math.h>
#include "RigidBody.h"
#include <myassert.h>
#include <float.h>

float RigidBody::g = 0.0;           // gravity acceleration, z accel = -g
float RigidBody::RKh = 0.0;
float RigidBody::RKh2 = 0.0;
float RigidBody::RKh4 = 0.0;
float RigidBody::RKh6 = 0.0;
									// Runke-Kutta step sizes, same for
									// all rigid bodies.

/********************************************************************
*   setGravity
*
*   params :
*       float g               -   gravity (m/s/s)
*
*   Set up the gravity. You must call this before using rigid bodies.
********************************************************************/

void RigidBody::setGravity(float g)
{
	RigidBody::g = g;
}

/********************************************************************
*   setStepSize
*
*   params :
*       float h               -   step size in secs
*
*   Set up the step size. Call this before making the step.
*   Applies to all bodies.
********************************************************************/

void RigidBody::setStepSize(float h)
{
	RKh = h;
	RKh2 = h / 2.0;
	RKh4 = h / 4.0;
	RKh6 = h / 6.0;
}

/***********************************************************************************
*   clear
*
*   Set the body as zero mass. May be useful if combine() is being used.
***********************************************************************************/

void RigidBody::clear(void)
{
	this->mass = 0.0;
    cOfM = Vector3(0.0, 0.0, 0.0);
    Ix = Iy = Iz = 0.0;
	Ixx = Iyy = Izz = 0.0;
    Ixy = Iyz = Izx = 0.0;
}

/***********************************************************************************
*   setBox
*
*   params :
*       float xRadius, yRadius, zRadius         -   half the size of the box
*       float mass                              -   total mass
*
*   Set the body as a rectangular box centred on the origin.
***********************************************************************************/

void RigidBody::setBox(float xRadius, float yRadius, float zRadius, float mass)
{
	this->mass = mass;
    cOfM = Vector3(0.0, 0.0, 0.0);
    Ix = Iy = Iz = 0.0;
	Ixx = mass / 3.0 * (yRadius * yRadius + zRadius * zRadius);
	Iyy = mass / 3.0 * (zRadius * zRadius + xRadius * xRadius);
	Izz = mass / 3.0 * (xRadius * xRadius + yRadius * yRadius);
    Ixy = Iyz = Izx = 0.0;
}

/***********************************************************************************
*   setCylinder
*
*   params :
*       float radius, height                    -   size
*       float mass                              -   total mass
*
*   Set the body as a cylinder centred on the origin. Cylinder goes along X axis from
*   -height/2 to +height/2.
***********************************************************************************/

void RigidBody::setCylinder(float radius, float height, float mass)
{
	this->mass = mass;
    cOfM = Vector3(0.0, 0.0, 0.0);
    Ix = Iy = Iz = 0.0;
	Ixx = mass * radius * radius / 2.0;
	Iyy = Izz = mass / 12.0 * (3 * radius * radius + height * height);
    Ixy = Iyz = Izx = 0.0;
}

/***********************************************************************************
*   setSphere
*
*   params :
*       float radius                            -   size
*       float mass                              -   total mass
*
*   Set the body as a sphere centred on the origin.
***********************************************************************************/

void RigidBody::setSphere(float radius, float mass)
{
	this->mass = mass;
    cOfM = Vector3(0.0, 0.0, 0.0);
    Ix = Iy = Iz = 0.0;
	Ixx = Iyy = Izz = mass * radius * radius * 2.0 / 5.0;
    Ixy = Iyz = Izx = 0.0;
}

/***********************************************************************************
*   combine
*
*   params :
*       RigidBody & other                       -   other body to combine with this one
*       Matrix3x3 & rotPos                      -   orientation of other body
*       Vector3 & linPos                        -   position of other body origin
*
*   Combine this body with another body. Other body has the following transform :
*       Rt = rotPos * Ro + linPos
*   where :
*       Ro = coord on other body
*       Rt = coord on this body
*   =>
*       Xt = Xo * m[0] + Yo * m[1] + Zo * m[2] + pos[0]
*       Yt = Xo * m[3] + Yo * m[4] + Zo * m[5] + pos[1]
*       Zt = Xo * m[6] + Yo * m[7] + Zo * m[8] + pos[2]
***********************************************************************************/

void RigidBody::combine(RigidBody & other, Matrix3x3 & rotPos, Vector3 & linPos)
{
// first combine centre of mass of both bodies into one.
// New centre of mass = (thisCOM * thisMass + otherCOM * otherMass) / totalMass
    Vector3 otherCOfM;                              // other COM in this body coords
    otherCOfM.transform(other.cOfM, rotPos, linPos);
    float totalMass = mass + other.mass;
    cOfM.mult(mass / totalMass);
    cOfM.mac(otherCOfM, other.mass / totalMass);
    mass = totalMass;

// now combine other moment of inertia with this one
    float *m = rotPos.m;
    Ix += other.Ix * m[0] + other.Iy * m[1] + other.Iz * m[2] + other.mass * linPos.m[0];
    Iy += other.Ix * m[3] + other.Iy * m[4] + other.Iz * m[5] + other.mass * linPos.m[1];
    Iz += other.Ix * m[6] + other.Iy * m[7] + other.Iz * m[8] + other.mass * linPos.m[2];
    float oMxx = (other.Iyy + other.Izz - other.Ixx) / 2.0;
    float oMyy = (other.Izz + other.Ixx - other.Iyy) / 2.0;
    float oMzz = (other.Ixx + other.Iyy - other.Izz) / 2.0;
                                                    // sigma of m * x^2, m * y^2 and m * z^2 (other body)
    float tMxx = oMxx * m[0] * m[0] +               // increment to sigma of m * x^2, m * y^2 and m * z^2 (this body)
        oMyy * m[1] * m[1] +
        oMzz * m[2] * m[2] +
        other.mass * linPos.m[0] * linPos.m[0] +
        2.0 * other.Ixy * m[0] * m[1] +
        2.0 * other.Iyz * m[1] * m[2] +
        2.0 * other.Izx * m[2] * m[0] +
        2.0 * other.Ix * m[0] * linPos.m[0] +
        2.0 * other.Iy * m[1] * linPos.m[0] +
        2.0 * other.Iz * m[2] * linPos.m[0];
    float tMyy = oMxx * m[3] * m[3] +
        oMyy * m[4] * m[4] +
        oMzz * m[5] * m[5] +
        other.mass * linPos.m[1] * linPos.m[1] +
        2.0 * other.Ixy * m[3] * m[4] +
        2.0 * other.Iyz * m[4] * m[5] +
        2.0 * other.Izx * m[5] * m[3] +
        2.0 * other.Ix * m[3] * linPos.m[1] +
        2.0 * other.Iy * m[4] * linPos.m[1] +
        2.0 * other.Iz * m[5] * linPos.m[1];
    float tMzz = oMxx * m[6] * m[6] +
        oMyy * m[7] * m[7] +
        oMzz * m[8] * m[8] +
        other.mass * linPos.m[2] * linPos.m[2] +
        2.0 * other.Ixy * m[6] * m[7] +
        2.0 * other.Iyz * m[7] * m[8] +
        2.0 * other.Izx * m[8] * m[6] +
        2.0 * other.Ix * m[6] * linPos.m[2] +
        2.0 * other.Iy * m[7] * linPos.m[2] +
        2.0 * other.Iz * m[8] * linPos.m[2];
    Ixx += tMyy + tMzz;
    Iyy += tMzz + tMxx;
    Izz += tMxx + tMyy;

    Ixy += oMxx * m[0] * m[3] +
        other.Ixy * m[1] * m[3] +
        other.Izx * m[2] * m[3] +
        other.Ix * linPos.m[0] * m[3] +
        other.Ixy * m[0] * m[4] +
        oMyy * m[1] * m[4] +
        other.Iyz * m[2] * m[4] +
        other.Iy * linPos.m[0] * m[4] +
        other.Izx * m[0] * m[5] +
        other.Iyz * m[1] * m[5] +
        oMzz * m[2] * m[5] +
        other.Iz * linPos.m[0] * m[5] +
        other.Ix * linPos.m[1] * m[0] +
        other.Iy * linPos.m[1] * m[1] +
        other.Iz * linPos.m[1] * m[2] +
        other.mass * linPos.m[0] * linPos.m[1];

    Iyz += oMxx * m[3] * m[6] +
        other.Ixy * m[4] * m[6] +
        other.Izx * m[5] * m[6] +
        other.Ix * linPos.m[1] * m[6] +
        other.Ixy * m[3] * m[7] +
        oMyy * m[4] * m[7] +
        other.Iyz * m[5] * m[7] +
        other.Iy * linPos.m[1] * m[7] +
        other.Izx * m[3] * m[8] +
        other.Iyz * m[4] * m[8] +
        oMzz * m[5] * m[8] +
        other.Iz * linPos.m[1] * m[8] +
        other.Ix * linPos.m[2] * m[3] +
        other.Iy * linPos.m[2] * m[4] +
        other.Iz * linPos.m[2] * m[5] +
        other.mass * linPos.m[1] * linPos.m[2];

    Izx += oMxx * m[6] * m[0] +
        other.Ixy * m[7] * m[0] +
        other.Izx * m[8] * m[0] +
        other.Ix * linPos.m[2] * m[0] +
        other.Ixy * m[6] * m[1] +
        oMyy * m[7] * m[1] +
        other.Iyz * m[8] * m[1] +
        other.Iy * linPos.m[2] * m[1] +
        other.Izx * m[6] * m[2] +
        other.Iyz * m[7] * m[2] +
        oMzz * m[8] * m[2] +
        other.Iz * linPos.m[2] * m[2] +
        other.Ix * linPos.m[0] * m[6] +
        other.Iy * linPos.m[0] * m[7] +
        other.Iz * linPos.m[0] * m[8] +
        other.mass * linPos.m[2] * linPos.m[0];
}

/***********************************************************************************
*   diagonalise
*
*   Find the transform to convert body into diagonalised body coords. This makes
*   the moment of inertia matrix diagonal, and the centre of mass at the origin.
*
*   The transform from old body coords to diagonalised body coords is :
*       Rd = diagRotPos * Ro + diagLinPos
***********************************************************************************/

void RigidBody::diagonalise(void)
{
    Matrix3x3 identity;
    identity.identity();
    Vector3 cOfM(-Ix / mass, -Iy / mass, -Iz / mass);   // centre of mass in old coords
    RigidBody old = *this;
    clear();
    combine(old, identity, cOfM);                       // shift centre of mass to origin

    Matrix3x3 I(                                    // moment of inertia (old coords)
        Ixx, -Ixy, -Izx,
        -Ixy, Iyy, -Iyz,
        -Izx, -Iyz, Izz);
    Vector3 d;
    Matrix3x3::eigenvectors(I, diagRotPos, d);      // find diagonalising transform with eigen vectors
    diagRotPos.transpose();
    diagLinPos.mult(diagRotPos, cOfM);
    diagIxx = d.x;                                  // eigen values are Ixx, Iyy, Izz
    diagIyy = d.y;
    diagIzz = d.z;

// now precalculate some useful numbers for speed ups
    oneOverDiagIxx = 1.0 / diagIxx;
    oneOverDiagIyy = 1.0 / diagIyy;
    oneOverDiagIzz = 1.0 / diagIzz;
    oneOverMass = 1.0 / mass;
    oneOverDiagIxx = 1.0 / diagIxx;
    oneOverDiagIyy = 1.0 / diagIyy;
    oneOverDiagIzz = 1.0 / diagIzz;
    minusMG = -mass * g;
    diagIyyMinusIzzOverIxx = (diagIyy - diagIzz) / diagIxx;
    diagIzzMinusIxxOverIyy = (diagIzz - diagIxx) / diagIyy;
    diagIxxMinusIyyOverIzz = (diagIxx - diagIyy) / diagIzz;
    diagIzzMinusIyy = diagIzz - diagIyy;
    diagIxxMinusIzz = diagIxx - diagIzz;
    diagIyyMinusIxx = diagIyy - diagIxx;
}

/***********************************************************************************
*   initialise
*
*   params :
*       Vector3 & oLinPos, & oLinVel                -   initial pos and vel (world coords)
*       Quaternion & oQRotPos                       -   initial rotation (world coords).
*                                                       will automatically be normalised.
*       Vector3 & oRotVel                           -   initial angular vel (old body coords)
*       float linKD, rotKD                          -   damping constants. will get scaled by the mass.
*
*   Set up the initial position of the body. The positions are for old body coords (these are automaticlly
*   transformed into diagonalised body coords).
*
*   The transform from old body coords to world coords is :
*   (1) Rw = oRotPos * Ro + oLinPos
*   The transform from old body coords to diagonalised body coords is :
*       Rd = diagRotPos * Ro + diagLinPos
*   (2) Ro =  trans(diagRotPos) * (Rd - diagLinPos)
*   The transform from diagonalised body coords to world coords is :
*   (3) Rw = rotPos * Rd + linPos
*
*   Now sub (2) in (1) :
*       Rw = oRotPos * trans(diagRotPos) * (Rd - diagLinPos) + oLinPos
*   (4) Rw = oRotPos * trans(diagRotPos) * Rd - oRotPos * trans(diagRotPos) * diagLinPos + oLinPos
*
*   Compare (4) with (3) :
*       rotPos = oRotPos * trans(diagRotPos)
*       linPos = oLinPos - oRotPos * trans(diagRotPos) * diagLinPos
*       rotVel = diagRotPos * oRotVel
*
*   Initial linVel :
*       COM is Rd = 0, so from (2) :
*       COM_Ro = trans(diagRotPos) * (-diagLinPos)
*       velocity of COM (old body coords) = oRotVel x COM_Ro
*       velocity of COM (world coords) =  oRotPos * (oRotVel x COM_Ro) + oLinVel
***********************************************************************************/

void RigidBody::initialise(Vector3 & oLinPos, Vector3 & oLinVel, Quaternion & oQRotPos, Vector3 & oRotVel,
    float linKD, float rotKD)
{
    Matrix3x3 oRotPos;
    oRotPos.fromQuatL2(oQRotPos, oQRotPos.mod2());    // find matrix of initial rotation
    Matrix3x3 M1;
    M1.transpose(diagRotPos);                       // M1 = trans(diagRotPos)
    rotPos.mult(oRotPos, M1);                       // rotPos = oRotPos * trans(diagRotPos)
    qRotPos.fromMatrix(rotPos);

    Vector3 t1;
    t1.mult(rotPos, diagLinPos);                    // t1 = oRotPos * trans(diagRotPos) * diagLinPos
    linPos.sub(oLinPos, t1);                        // linPos = oLinPos - oRotPos * trans(diagRotPos) * diagLinPos
    rotVel.mult(diagRotPos, oRotVel);

    Vector3 COM;                                    // centre of mass in old body coords
    COM.transMult(diagRotPos, diagLinPos);
    COM.neg();                                      // COM = trans(diagRotPos) * (-diagLinPos)
    t1.cross(oRotVel, COM);                         // t1 = oRotVel x COM_Ro
    linVel.mult(oRotPos, t1);
    linVel.add(oLinVel);                            // velocity of COM (world coords)

    this->linKD = -linKD * mass;                    // scale damping by mass
    rotKDx = -rotKD * diagIxx;
    rotKDy = -rotKD * diagIyy;
    rotKDz = -rotKD * diagIzz;
}

/***********************************************************************************
*   RKInit
*
*   Initialise for Runge-Kutta or Euler integration. Call this first. Then find linForce and rotForce
*   from linPos, linVel, rotPos, rotVel.
***********************************************************************************/

void RigidBody::RKInit(void)
{
	linPos0 = linPos;                               // save position before step is made
    qRotPos0 = qRotPos;
    rotPos0 = rotPos;
    linVel0 = linVel;
    rotVel0 = rotVel;

	linForce.x = linKD * linVel.x;                  // built in damping and gravity
    linForce.y = linKD * linVel.y;
    linForce.z = linKD * linVel.z + minusMG;
	rotForce.x = rotKDx * rotVel.x;
	rotForce.y = rotKDy * rotVel.y;
	rotForce.z = rotKDz * rotVel.z;
}

/***********************************************************************************
*   eulerStep
*
*   One and only step for Euler integration.
***********************************************************************************/

void RigidBody::eulerStep(void)
{
    linPos.x += RKh * linVel.x;
    linPos.y += RKh * linVel.y;
    linPos.z += RKh * linVel.z;

    linVel.x += linForce.x * oneOverMass * RKh;
    linVel.y += linForce.y * oneOverMass * RKh;
    linVel.z += linForce.z * oneOverMass * RKh;

    qRotPos.w -= RKh2 * (qRotPos0.x * rotVel.x + qRotPos0.y * rotVel.y + qRotPos0.z * rotVel.z);
    qRotPos.x += RKh2 * (qRotPos0.w * rotVel.x - qRotPos0.z * rotVel.y + qRotPos0.y * rotVel.z);
    qRotPos.y += RKh2 * (qRotPos0.z * rotVel.x + qRotPos0.w * rotVel.y - qRotPos0.x * rotVel.z);
    qRotPos.z += RKh2 * (qRotPos0.x * rotVel.y - qRotPos0.y * rotVel.x + qRotPos0.w * rotVel.z);

    float l2 = qRotPos.x * qRotPos.x + qRotPos.y * qRotPos.y + qRotPos.z * qRotPos.z + qRotPos.w * qRotPos.w;
    if (l2 < 0.9999 || l2 > 1.0001) {               // normalise to unit length (approximately)
        float n = (l2 + 1.0) / (2.0 * l2);
        qRotPos.x *= n;
        qRotPos.y *= n;
        qRotPos.z *= n;
        qRotPos.w *= n;
        l2 *= n * n;
    }
    rotPos.fromQuatL2(qRotPos, l2);                 // make new matrix from quaternion.
                                                    // this will be orthogonal.

    rotVel1.x = rotVel.x + (rotVel.y * rotVel.z * diagIyyMinusIzzOverIxx + rotForce.x * oneOverDiagIxx) * RKh;
    rotVel1.y = rotVel.y + (rotVel.z * rotVel.x * diagIzzMinusIxxOverIyy + rotForce.y * oneOverDiagIyy) * RKh;
    rotVel1.z = rotVel.z + (rotVel.x * rotVel.y * diagIxxMinusIyyOverIzz + rotForce.z * oneOverDiagIzz) * RKh;
    rotVel = rotVel1;
                                                    // calculate angular velocity with Eulers equations
}

/***********************************************************************************
*   RKStep1
*
*   Runge Kutta step 1.
***********************************************************************************/

void RigidBody::RKStep1(void)
{
    linAcc1.x = linForce.x * oneOverMass * RKh2;
    linAcc1.y = linForce.y * oneOverMass * RKh2;
    linAcc1.z = linForce.z * oneOverMass * RKh2;
    rotAcc1.x = (rotVel.y * rotVel.z * diagIyyMinusIzzOverIxx + rotForce.x * oneOverDiagIxx) * RKh2;
    rotAcc1.y = (rotVel.z * rotVel.x * diagIzzMinusIxxOverIyy + rotForce.y * oneOverDiagIyy) * RKh2;
    rotAcc1.z = (rotVel.x * rotVel.y * diagIxxMinusIyyOverIzz + rotForce.z * oneOverDiagIzz) * RKh2;

    linVel1 = linVel;
    rotVel1 = rotVel;

    linPos.x += RKh2 * linVel1.x;
    linPos.y += RKh2 * linVel1.y;
    linPos.z += RKh2 * linVel1.z;
    qRotPos.w -= RKh4 * (qRotPos0.x * rotVel1.x + qRotPos0.y * rotVel1.y + qRotPos0.z * rotVel1.z);
    qRotPos.x += RKh4 * (qRotPos0.w * rotVel1.x - qRotPos0.z * rotVel1.y + qRotPos0.y * rotVel1.z);
    qRotPos.y += RKh4 * (qRotPos0.z * rotVel1.x + qRotPos0.w * rotVel1.y - qRotPos0.x * rotVel1.z);
    qRotPos.z += RKh4 * (qRotPos0.x * rotVel1.y - qRotPos0.y * rotVel1.x + qRotPos0.w * rotVel1.z);
    float l2 = qRotPos.x * qRotPos.x + qRotPos.y * qRotPos.y + qRotPos.z * qRotPos.z + qRotPos.w * qRotPos.w;
    if (l2 < 0.9999 || l2 > 1.0001) {               // normalise to unit length (approximately)
        float n = (l2 + 1.0) / (2.0 * l2);
        qRotPos.x *= n;
        qRotPos.y *= n;
        qRotPos.z *= n;
        qRotPos.w *= n;
        l2 *= n * n;
    }
    rotPos.fromQuatL2(qRotPos, l2);                 // make new matrix from quaternion.
                                                    // this will be orthogonal.
    linVel.x += linAcc1.x;
    linVel.y += linAcc1.y;
    linVel.z += linAcc1.z;
    rotVel.x += rotAcc1.x;
    rotVel.y += rotAcc1.y;
    rotVel.z += rotAcc1.z;

	linForce.x = linKD * linVel.x;                  // built in damping and gravity
    linForce.y = linKD * linVel.y;
    linForce.z = linKD * linVel.z + minusMG;
	rotForce.x = rotKDx * rotVel.x;
	rotForce.y = rotKDy * rotVel.y;
	rotForce.z = rotKDz * rotVel.z;
}

/***********************************************************************************
*   RKStep2
*
*   Runge Kutta step 2.
***********************************************************************************/

void RigidBody::RKStep2(void)
{
    linAcc2.x = linForce.x * oneOverMass * RKh2;
    linAcc2.y = linForce.y * oneOverMass * RKh2;
    linAcc2.z = linForce.z * oneOverMass * RKh2;
    rotAcc2.x = (rotVel.y * rotVel.z * diagIyyMinusIzzOverIxx + rotForce.x * oneOverDiagIxx) * RKh2;
    rotAcc2.y = (rotVel.z * rotVel.x * diagIzzMinusIxxOverIyy + rotForce.y * oneOverDiagIyy) * RKh2;
    rotAcc2.z = (rotVel.x * rotVel.y * diagIxxMinusIyyOverIzz + rotForce.z * oneOverDiagIzz) * RKh2;

    linVel2 = linVel;
    rotVel2 = rotVel;

    linPos.x = linPos0.x + RKh2 * linVel2.x;
    linPos.y = linPos0.y + RKh2 * linVel2.y;
    linPos.z = linPos0.z + RKh2 * linVel2.z;
    qRotPos.w = qRotPos0.w - RKh4 * (qRotPos0.x * rotVel2.x + qRotPos0.y * rotVel2.y + qRotPos0.z * rotVel2.z);
    qRotPos.x = qRotPos0.x + RKh4 * (qRotPos0.w * rotVel2.x - qRotPos0.z * rotVel2.y + qRotPos0.y * rotVel2.z);
    qRotPos.y = qRotPos0.y + RKh4 * (qRotPos0.z * rotVel2.x + qRotPos0.w * rotVel2.y - qRotPos0.x * rotVel2.z);
    qRotPos.z = qRotPos0.z + RKh4 * (qRotPos0.x * rotVel2.y - qRotPos0.y * rotVel2.x + qRotPos0.w * rotVel2.z);
    float l2 = qRotPos.x * qRotPos.x + qRotPos.y * qRotPos.y + qRotPos.z * qRotPos.z + qRotPos.w * qRotPos.w;
    if (l2 < 0.9999 || l2 > 1.0001) {               // normalise to unit length (approximately)
        float n = (l2 + 1.0) / (2.0 * l2);
        qRotPos.x *= n;
        qRotPos.y *= n;
        qRotPos.z *= n;
        qRotPos.w *= n;
        l2 *= n * n;
    }
    rotPos.fromQuatL2(qRotPos, l2);                 // make new matrix from quaternion.
                                                    // this will be orthogonal.
    linVel.x = linVel1.x + linAcc2.x;
    linVel.y = linVel1.y + linAcc2.y;
    linVel.z = linVel1.z + linAcc2.z;
    rotVel.x = rotVel1.x + rotAcc2.x;
    rotVel.y = rotVel1.y + rotAcc2.y;
    rotVel.z = rotVel1.z + rotAcc2.z;

	linForce.x = linKD * linVel.x;                 // built in damping and gravity
    linForce.y = linKD * linVel.y;
    linForce.z = linKD * linVel.z + minusMG;
	rotForce.x = rotKDx * rotVel.x;
	rotForce.y = rotKDy * rotVel.y;
	rotForce.z = rotKDz * rotVel.z;
}

/***********************************************************************************
*   RKStep3
*
*   Runge Kutta step 3.
***********************************************************************************/

void RigidBody::RKStep3(void)
{
    linAcc3.x = linForce.x * oneOverMass * RKh;
    linAcc3.y = linForce.y * oneOverMass * RKh;
    linAcc3.z = linForce.z * oneOverMass * RKh;
    rotAcc3.x = (rotVel.y * rotVel.z * diagIyyMinusIzzOverIxx + rotForce.x * oneOverDiagIxx) * RKh;
    rotAcc3.y = (rotVel.z * rotVel.x * diagIzzMinusIxxOverIyy + rotForce.y * oneOverDiagIyy) * RKh;
    rotAcc3.z = (rotVel.x * rotVel.y * diagIxxMinusIyyOverIzz + rotForce.z * oneOverDiagIzz) * RKh;

    linVel3 = linVel;
    rotVel3 = rotVel;

    linPos.x = linPos0.x + RKh * linVel3.x;
    linPos.y = linPos0.y + RKh * linVel3.y;
    linPos.z = linPos0.z + RKh * linVel3.z;
    qRotPos.w = qRotPos0.w - RKh2 * (qRotPos0.x * rotVel3.x + qRotPos0.y * rotVel3.y + qRotPos0.z * rotVel3.z);
    qRotPos.x = qRotPos0.x + RKh2 * (qRotPos0.w * rotVel3.x - qRotPos0.z * rotVel3.y + qRotPos0.y * rotVel3.z);
    qRotPos.y = qRotPos0.y + RKh2 * (qRotPos0.z * rotVel3.x + qRotPos0.w * rotVel3.y - qRotPos0.x * rotVel3.z);
    qRotPos.z = qRotPos0.z + RKh2 * (qRotPos0.x * rotVel3.y - qRotPos0.y * rotVel3.x + qRotPos0.w * rotVel3.z);
    float l2 = qRotPos.x * qRotPos.x + qRotPos.y * qRotPos.y + qRotPos.z * qRotPos.z + qRotPos.w * qRotPos.w;
    if (l2 < 0.9999 || l2 > 1.0001) {               // normalise to unit length (approximately)
        float n = (l2 + 1.0) / (2.0 * l2);
        qRotPos.x *= n;
        qRotPos.y *= n;
        qRotPos.z *= n;
        qRotPos.w *= n;
        l2 *= n * n;
    }
    rotPos.fromQuatL2(qRotPos, l2);                 // make new matrix from quaternion.
                                                    // this will be orthogonal.
    linVel.x = linVel1.x + linAcc3.x;
    linVel.y = linVel1.y + linAcc3.y;
    linVel.z = linVel1.z + linAcc3.z;
    rotVel.x = rotVel1.x + rotAcc3.x;
    rotVel.y = rotVel1.y + rotAcc3.y;
    rotVel.z = rotVel1.z + rotAcc3.z;

	linForce.x = linKD * linVel.x;                  // built in damping and gravity
    linForce.y = linKD * linVel.y;
    linForce.z = linKD * linVel.z + minusMG;
	rotForce.x = rotKDx * rotVel.x;
	rotForce.y = rotKDy * rotVel.y;
	rotForce.z = rotKDz * rotVel.z;
}

/***********************************************************************************
*   RKStep4
*
*   Runge Kutta step 4.
***********************************************************************************/

void RigidBody::RKStep4(void)
{
    linPos.x = linPos0.x + (linVel1.x + 2.0 * (linVel2.x + linVel3.x) + linVel.x) * RKh6;
    linPos.y = linPos0.y + (linVel1.y + 2.0 * (linVel2.y + linVel3.y) + linVel.y) * RKh6;
    linPos.z = linPos0.z + (linVel1.z + 2.0 * (linVel2.z + linVel3.z) + linVel.z) * RKh6;
    Vector3 rotVel4;
    rotVel4.x = (rotVel1.x + rotVel.x) * 0.5 + rotVel2.x + rotVel3.x;
    rotVel4.y = (rotVel1.y + rotVel.y) * 0.5 + rotVel2.y + rotVel3.y;
    rotVel4.z = (rotVel1.z + rotVel.z) * 0.5 + rotVel2.z + rotVel3.z;
    qRotPos.w = qRotPos0.w - RKh6 * (qRotPos0.x * rotVel4.x + qRotPos0.y * rotVel4.y + qRotPos0.z * rotVel4.z);
    qRotPos.x = qRotPos0.x + RKh6 * (qRotPos0.w * rotVel4.x - qRotPos0.z * rotVel4.y + qRotPos0.y * rotVel4.z);
    qRotPos.y = qRotPos0.y + RKh6 * (qRotPos0.z * rotVel4.x + qRotPos0.w * rotVel4.y - qRotPos0.x * rotVel4.z);
    qRotPos.z = qRotPos0.z + RKh6 * (qRotPos0.x * rotVel4.y - qRotPos0.y * rotVel4.x + qRotPos0.w * rotVel4.z);
    float l2 = qRotPos.x * qRotPos.x + qRotPos.y * qRotPos.y + qRotPos.z * qRotPos.z + qRotPos.w * qRotPos.w;
    if (l2 < 0.9999 || l2 > 1.0001) {               // normalise to unit length (approximately)
        float n = (l2 + 1.0) / (2.0 * l2);
        qRotPos.x *= n;
        qRotPos.y *= n;
        qRotPos.z *= n;
        qRotPos.w *= n;
        l2 *= n * n;
    }
    rotPos.fromQuatL2(qRotPos, l2);                 // make new matrix from quaternion.

  	linVel.x = linVel1.x + (linAcc1.x + 2.0 * linAcc2.x + linAcc3.x + linForce.x * oneOverMass * RKh2) * (1.0 / 3.0);
  	linVel.y = linVel1.y + (linAcc1.y + 2.0 * linAcc2.y + linAcc3.y + linForce.y * oneOverMass * RKh2) * (1.0 / 3.0);
  	linVel.z = linVel1.z + (linAcc1.z + 2.0 * linAcc2.z + linAcc3.z + linForce.z * oneOverMass * RKh2) * (1.0 / 3.0);

  	rotVel4.x = rotVel1.x + (rotAcc1.x + 2.0 * rotAcc2.x + rotAcc3.x + (rotVel.y * rotVel.z * diagIyyMinusIzzOverIxx + rotForce.x * oneOverDiagIxx) * RKh2) * (1.0 / 3.0);
  	rotVel4.y = rotVel1.y + (rotAcc1.y + 2.0 * rotAcc2.y + rotAcc3.y + (rotVel.z * rotVel.x * diagIzzMinusIxxOverIyy + rotForce.y * oneOverDiagIyy) * RKh2) * (1.0 / 3.0);
  	rotVel4.z = rotVel1.z + (rotAcc1.z + 2.0 * rotAcc2.z + rotAcc3.z + (rotVel.x * rotVel.y * diagIxxMinusIyyOverIzz + rotForce.z * oneOverDiagIzz) * RKh2) * (1.0 / 3.0);
    rotVel = rotVel4;
}

/***********************************************************************************
*   addBodyTorque
*
*   params :
*       Vector3 & torque    -   torque, diagonalised body coords
*
*   Add a torque given in diagonalised body coords.
***********************************************************************************/

void RigidBody::addBodyTorque(Vector3 & torque)
{
	rotForce.add(torque);
}

/***********************************************************************************
*   addWorldTorque
*
*   params :
*       Vector3 & torque    -   torque, world coords
*
*   Add a torque given in world coords.
***********************************************************************************/

void RigidBody::addWorldTorque(Vector3 & torque)
{
	Vector3 torqueBody;                     // torque in body coords
	torqueBody.transMult(rotPos, torque);
	rotForce.add(torqueBody);
}

/***********************************************************************************
*   addBodyBodyForce
*
*   params :
*       Vector3 & force     -   direction of force, diag body coords
*       Vector3 & pos       -   position of application, diag body coords
*   returns :
*       void
*
*   Force is given in diag body coords, point is diag body coords
***********************************************************************************/

void RigidBody::addBodyBodyForce(Vector3 & force, Vector3 & pos)
{
	Vector3 torqueBody;             // torque in body coords
	torqueBody.cross(pos, force);
	rotForce.add(torqueBody);
	Vector3 forceWorld;
	forceWorld.mult(rotPos, force);
	linForce.add(forceWorld);      // include the linear force
}

/***********************************************************************************
*   addWorldBodyForce
*
*   params :
*       Vector3 & force     -   direction of force, world coords
*       Vector3 & pos       -   position of application, diag body coords
*   returns :
*       void
*
*   Force is given in world coords, point is body coords
***********************************************************************************/

void RigidBody::addWorldBodyForce(Vector3 & force, Vector3 & pos)
{
	Vector3 forceBody;              // force in body coords
	forceBody.transMult(rotPos, force);
	Vector3 torqueBody;             // torque in body coords
	torqueBody.cross(pos, forceBody);
	rotForce.add(torqueBody);
	linForce.add(force);           // include the linear force
}

/***********************************************************************************
*   addWorldWorldForce
*
*   params :
*       Vector3 & force     -   direction of force, world coords
*       Vector3 & pos       -   position of application, world coords
*
*   Force is given in world coords, point is world coords
***********************************************************************************/

void RigidBody::addWorldWorldForce(Vector3 & force, Vector3 & pos)
{
	Vector3 torque;                 // torque in world coords
    Vector3 posCOM;                 // pos relative to COM
    posCOM.sub(pos, linPos);
	torque.cross(posCOM, force);
	Vector3 torqueBody;                     // torque in body coords
	torqueBody.transMult(rotPos, torque);
	rotForce.add(torqueBody);
	linForce.add(force);           // include the linear force
}

/********************************************************************
*   findWorldPos
*
*   params :
*       Vector3 & bodyPos       -   point in diag body coords
*       Vector3 & worldPos      -   returns position of point in world coords
*
*   Find the world coords of a point on the body.
********************************************************************/

void RigidBody::findWorldPos(Vector3 & bodyPos, Vector3 & worldPos)
{
    worldPos.transform(bodyPos, rotPos, linPos);
}

/********************************************************************
*   findBodyPos
*
*   params :
*       Vector3 & worldPos       -   position of point in world coords
*       Vector3 & bodyPos        -   returns point in diag body coords
*
*   Find the body coords of a point in world coords.
********************************************************************/

void RigidBody::findBodyPos(Vector3 & worldPos, Vector3 & bodyPos)
{
    bodyPos.invTransform(worldPos, rotPos, linPos);
}

/********************************************************************
*   findBodyWorldVel
*
*   params :
*       Vector3 & bodyPos       -   point in diag body coords
*       Vector3 & worldVel      -   returns velocity of point in world coords
*
*   Find the velocity in world coords of a point on the body.
********************************************************************/

void RigidBody::findBodyWorldVel(Vector3 & bodyPos, Vector3 & worldVel)
{
	Vector3 bodyVel;                            // velocity in body coords
    bodyVel.cross(rotVel, bodyPos);
	worldVel.transform(bodyVel, rotPos, linVel);
}

/********************************************************************
*   findWorldPosVel
*
*   params :
*       Vector3 & bodyPos       -   point in diag body coords
*       Vector3 & worldPos      -   returns position of point in world coords
*       Vector3 & worldVel      -   returns velocity of point in world coords
*
*   Find the world coords and world velocity of a point on the body.
********************************************************************/

void RigidBody::findWorldPosVel(Vector3 & bodyPos, Vector3 & worldPos, Vector3 & worldVel)
{
    worldPos.transform(bodyPos, rotPos, linPos);
	Vector3 bodyVel;                            // velocity in body coords
    bodyVel.cross(rotVel, bodyPos);
	worldVel.transform(bodyVel, rotPos, linVel);
}

/********************************************************************
*   findWorldWorldVel
*
*   params :
*       Vector3 & worldPos      -   point in world coords
*       Vector3 & worldVel      -   returns velocity of point in world coords
*
*   Find the velocity in world coords of a point given in world coords.
********************************************************************/

void RigidBody::findWorldWorldVel(Vector3 & worldPos, Vector3 & worldVel)
{
	Vector3 pos;                                // pos relative to COM in world coords
	pos.sub(worldPos, linPos);
    Vector3 worldRotVel;                        // angular velocity in world coords
    worldRotVel.mult(rotPos, rotVel);
    worldVel.cross(worldRotVel, pos);
	worldVel.add(linVel);
}

/********************************************************************
*   energy
*
*   returns :
*       float                       -   total energy
*
*   Find the total energy of the body (for testing constraints).
*   Total energy = 1/2 * mv^2 + 1/2 * Iw^2 + mgz
********************************************************************/

float RigidBody::energy(void)
{
    float Iw2 = diagIxx * rotVel.x * rotVel.x + diagIyy * rotVel.y * rotVel.y + diagIzz * rotVel.z * rotVel.z;
    float v2 = linVel.x * linVel.x + linVel.y * linVel.y + linVel.z * linVel.z;
    return 0.5 * (mass * v2 + Iw2) + mass * g * linPos.z;
}

/********************************************************************
*   momentum
*
*   params :
*       Vector3 & lin, & rot                    -   returns linear and angular momentum in world coords
*
*   Find the momentum in world coords (for testing).
*   Lo = Io.wo = R.Ib.Rt.R.wb = R.Ib.wb
********************************************************************/

void RigidBody::momentum(Vector3 & lin, Vector3 & rot)
{
    lin.mult(linVel, mass);
    Vector3 temp(rotVel.x * diagIxx, rotVel.y * diagIyy, rotVel.z * diagIzz);
                                                    // momentum in body coords (Ib.wb)
    rot.mult(rotPos, temp);                         // Lo = R.Ib.wb
}

/********************************************************************
*   undo
*
*   Go back in time to start of previous step. For collision detection.
********************************************************************/

void RigidBody::undo(void)
{
	linPos = linPos0;
    qRotPos = qRotPos0;
    rotPos = rotPos0;
    linVel = linVel0;
    rotVel = rotVel0;
}

