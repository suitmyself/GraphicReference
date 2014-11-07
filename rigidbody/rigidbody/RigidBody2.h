#ifndef RIGIDBODY_H
#define RIGIDBODY_H

#include "Matrix.h"

class RigidBody
{
// Runge-Kutta integration stuff
    Vector3 linPos0, linVel0, rotVel0;
    Quaternion qRotPos0;
    Matrix3x3 rotPos0;
    Vector3 linVel1, rotVel1, linAcc1, rotAcc1;
    Vector3 linVel2, rotVel2, linAcc2, rotAcc2;
    Vector3 linVel3, rotVel3, linAcc3, rotAcc3;

public:
// properties of body (these are public in case you need access, but don't change them)
    float mass, oneOverMass;                        // mass
    float minusMG;                                  // -mass * g
    Vector3 cOfM;                                   // centre of mass (in old body coords)
    float Ix, Iy, Iz;                               // Ix = sum of mx
    float Ixx, Iyy, Izz;                            // Ixx = sum of m(y^2 + z^2)
    float Ixy, Iyz, Izx;                            // Ixy = sum of mxy
                                                    // these are in old body coords

    Matrix3x3 diagRotPos;                           // transform from old body coords to diagonalised body coords.
    Vector3 diagLinPos;                             // Rd = diagRotPos * Rb + diagLinPos
    float diagIxx, oneOverDiagIxx;                  // transformed moment of inertia (diagonalised body coords)
    float diagIyy, oneOverDiagIyy;                  // also 1 / Ixx * h for Euler integrater
    float diagIzz, oneOverDiagIzz;
    float diagIyyMinusIzzOverIxx;
    float diagIzzMinusIxxOverIyy;
    float diagIxxMinusIyyOverIzz;
                                                    // (Iyy - Izz) / Ixx * h, etc..
    float diagIzzMinusIyy;
    float diagIxxMinusIzz;
    float diagIyyMinusIxx;
    float linKD;                                    // linear damping (-ve and multiplied by mass)
    float rotKDx, rotKDy, rotKDz;                   // angular damping (-ve and multiplied by Ixx, etc..)

	static float g;                                 // gravity, z accel = -g
    static float RKh, RKh2, RKh4, RKh6;             // step sizes

// current state of body
	Vector3 linPos, linVel;                         // linear position and velocity
									                // of centre of mass (Rd = 0). (world coords)
	Matrix3x3 rotPos;                               // rotational position. transform from
									                // diagonalised body coord to world coords is
									                // Rw = rotPos * Rd + linPos
    Quaternion qRotPos;                             // quaternion of above matrix (matrix and quaternion are kept the same)
	Vector3 rotVel;                                 // rotational velocity (diagonalised body coords)

// current forces on body
	Vector3 linForce, rotForce;                     // current forces

// functions to set up rigid body mass distribution
	void setBox(float xRadius, float yRadius, float zRadius, float mass);
    void setCylinder(float radius, float height, float mass);
    void setSphere(float radius, float mass);
    void clear(void);
    void combine(RigidBody & other, Matrix3x3 & rotPos, Vector3 & linPos);

// functions for initialising rigid body
    void diagonalise(void);
    void initialise(Vector3 & oLinPos, Vector3 & oLinVel, Quaternion & oQRotPos, Vector3 & oRotVel,
        float linKD, float rotKD);
	static void setGravity(float g);
	static void setStepSize(float h);

// functions to update from one frame to another
    void RKInit(void);
    void eulerStep(void);
    void RKStep1(void);
    void RKStep2(void);
    void RKStep3(void);
    void RKStep4(void);
    void undo(void);

// functions to add force and torque
	void addBodyTorque(Vector3 & torque);
	void addWorldTorque(Vector3 & torque);
	void addWorldWorldForce(Vector3 & force, Vector3 & pos);
	void addWorldBodyForce(Vector3 & force, Vector3 & pos);
	void addBodyBodyForce(Vector3 & force, Vector3 & pos);

// functions to get current position and velocity of points
	void findWorldPos(Vector3 & bodyPos, Vector3 & worldPos);
	void findBodyPos(Vector3 & worldPos, Vector3 & bodyPos);
	void findBodyWorldVel(Vector3 & bodyPos, Vector3 & worldVel);
	void findWorldWorldVel(Vector3 & worldPos, Vector3 & worldVel);
    void findWorldPosVel(Vector3 & bodyPos, Vector3 & worldPos, Vector3 & worldVel);

// functions for energy and momentum (for testing)
    float energy(void);
    void momentum(Vector3 & lin, Vector3 & rot);
};

#endif


