/***********************************************************************
*   matrix.h
***********************************************************************/

#ifndef MATRIX_H
#define MATRIX_H
#include <math.h>

#define EIGENROTATE(a,i,j,k,l)      g=a.n[i][j];h=a.n[k][l];a.n[i][j]=g-s*(h+g*tau);\
                                    a.n[k][l]=h+s*(g-h*tau);
                                    // macro for Matrix3x3::eigenvectors() function.
#define FABS(X)                     asm {and X,0x7fffffff;}                         // X = fabs(X)

class Matrix3x3;
class Quaternion;

/***********************************************************************************
*   class Vector3
*
*   A general 3D vector and some methods for it. It uses float precision data.
*   This vectors is assumed to be a column vector when multiplied by matricies.
***********************************************************************************/

class Vector3
{
public:
    union {
        struct {
        	float x, y, z;                             // elements here
        };
        float m[3];
    };

	Vector3() { }                                   // blank constructor
    Vector3(float x, float y, float z);   // constructor with 3 values

   void add(Vector3 & A);
   void add(Vector3 & A, Vector3 & B);
   void sub(Vector3 & A);
   void sub(Vector3 & A, Vector3 & B);
   void mult(Matrix3x3 & M, Vector3 & A);
   void mult2(Matrix3x3 & M1, Vector3 & A1, Matrix3x3 & M2, Vector3 & A2);
   void mac(Matrix3x3 & M, Vector3 & A);
   void transMult(Matrix3x3 & M, Vector3 & A);
   void mult(Vector3 & A, float scalar);
   void mult(float scalar);                // multiply vector by scalar
   void mac(Vector3 & A, Vector3 & B, float scalar);
   void mac(Vector3 & A, float scalar);
   void msc(Vector3 & A, Vector3 & B, float scalar);
   void msc(Vector3 & A, float scalar);
   void mac2(Vector3 & A, Vector3 & B, float scalarB, Vector3 & C, float scalarC);
   void msc2(Vector3 & A, Vector3 & B, float scalarB, Vector3 & C, float scalarC);
   void lerp(Vector3 & A, Vector3 & B, float t);
   void herp(Vector3 & p0, Vector3 & p1, Vector3 & p2, Vector3 & p3, float t);
   int normalise(Vector3 & A);
   int normalise(void);                 // scale to length = 1
   void neg(void);
   void neg(Vector3 & A);
   float mod(void);                          // length of vector
   float mod2(void);                         // length ^ 2 of vector
   float dot(Vector3 & A);
   void cross(Vector3 & A, Vector3 & B);
   float dist2(Vector3 & A);
   float dist(Vector3 & A);
   void transform(Vector3 & r, Matrix3x3 & rotPos, Vector3 & linPos);
   void invTransform(Vector3 & r, Matrix3x3 & rotPos, Vector3 & linPos);
   static void transformArray(Vector3 *source, Vector3 *dest, int N, Matrix3x3 & rotPos, Vector3 & linPos);
   static void invTransformArray(Vector3 *source, Vector3 *dest, int N, Matrix3x3 & rotPos, Vector3 & linPos);
   static void transformNormalArray(Vector3 *source, Vector3 *dest, int N, Matrix3x3 & rotPos);
   static void invTransformNormalArray(Vector3 *source, Vector3 *dest, int N, Matrix3x3 & rotPos);

    static Vector3 ZERO;
    static Vector3 I;
    static Vector3 J;
    static Vector3 K;
    static Vector3 MI;
    static Vector3 MJ;
    static Vector3 MK;
};

/***********************************************************************************
*   class Matrix3x3
*
*   A general 3 x 3 matrix and some methods for it. It uses float precision data.
*
*   2 ways to think of the matrix are :
*
*       |   m[0]        m[1]        m[2]    |
*       |   m[3]        m[4]        m[5]    |           1D array
*       |   m[6]        m[7]        m[8]    |
*
*       |   n[0][0]     n[0][1]     n[0][2] |
*       |   n[1][0]     n[1][1]     n[1][2] |           2D array
*       |   n[2][0]     n[2][1]     n[2][2] |
***********************************************************************************/

class Matrix3x3
{
public:
    union {
    	float m[9];                         // elements here
        float n[3][3];
    };

	Matrix3x3() { }                             // blank constructor
    Matrix3x3(float m11, float m12, float m13, float m21, float m22, float m23,
	float m31, float m32, float m33);
													// constructor with 9 values
	void identity(void);
	void neg(void);
	void neg(Matrix3x3 & A);
	void add(Matrix3x3 & A);
	void add(Matrix3x3 & A, Matrix3x3 & B);
	void sub(Matrix3x3 & A);
	void sub(Matrix3x3 & A, Matrix3x3 & B);
	void mult(float s);
	void mult(Matrix3x3 & A, float s);
	void mac(Matrix3x3 & A, Matrix3x3 & B, float s);
	void mac(Matrix3x3 & A, float s);
	void transpose(void);
	void transpose(Matrix3x3 & A);
	void mult(Matrix3x3 & A, Matrix3x3 & B);
	void multAtB(Matrix3x3 & A, Matrix3x3 & B);
	void multABt(Matrix3x3 & A, Matrix3x3 & B);
	float det(void);
	int inverse(Matrix3x3 & M);
	void exp(Matrix3x3 & A);

	void xAxisRotation(float thetaX);
	void yAxisRotation(float thetaY);
	void zAxisRotation(float thetaZ);
														// set matrix as x, y, z axis
														// rotations.
	void arbAxisRotation(Vector3 & axis, float theta);
										  // set matrix as rotation about
										  // the given axis.
   void rpyRotation(float roll, float pitch, float yaw);
														// constructor from roll, pitch, yaw angles
	void star(Vector3 & a);
    void fromQuat(Quaternion & q);                      // create from quaternion
    void fromQuatL2(Quaternion & q, float l2);
    static int eigenvectors(Matrix3x3 & M, Matrix3x3 & V, Vector3 & d);
                                                        // find eigen vectors of symmetric matrix
    static Matrix3x3 ZERO;
    static Matrix3x3 I;
};

/******************************************************************************
*   class Quaternion
*
*   Quaternion for doing rotations.
******************************************************************************/

class Quaternion
{
public:
    union {
        struct {
            float x, y, z, w;
        };
        float m[4];
    };

    Quaternion() {}
    Quaternion(float w, float x, float y, float z);
    void rotation(Vector3 & axis, float angle);  // define as a rotation
    void fromMatrix(Matrix3x3 & M);
    void mult(Quaternion & A, Quaternion & B);
    void normalise(void);
    float mod(void);
    float mod2(void);

    static Quaternion I;
};

#endif


