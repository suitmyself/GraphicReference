/***********************************************************************
*   matrix.cpp
*
*   Version 2.0
*
*   Classes in here :
*       Vector3         -   3D column vector
*       Matrix3x3       -   3 x 3 matrix
*       Quaternion      -   quaternion for rotation stuff
***********************************************************************/

#include <math.h>
#include "Matrix.h"

/******************************************************************************
*   Define some useful static constants
******************************************************************************/

Vector3 Vector3::ZERO = Vector3(0.0, 0.0, 0.0);
Vector3 Vector3::I = Vector3(1.0, 0.0, 0.0);
Vector3 Vector3::J = Vector3(0.0, 1.0, 0.0);
Vector3 Vector3::K = Vector3(0.0, 0.0, 1.0);
Vector3 Vector3::MI = Vector3(-1.0, 0.0, 0.0);
Vector3 Vector3::MJ = Vector3(0.0, -1.0, 0.0);
Vector3 Vector3::MK = Vector3(0.0, 0.0, -1.0);
Matrix3x3 Matrix3x3::ZERO = Matrix3x3(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
Matrix3x3 Matrix3x3::I = Matrix3x3(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
Quaternion Quaternion::I = Quaternion(1.0, 0.0, 0.0, 0.0);

/***********************************************************************************
*   Matrix3x3::Matrix3x3
*
*   params :
*       float m11, m12, m13, m21, m22, m23, m31, m32, m33  -   initial value
*   returns :
*       none
*
*   Constructor with initial value.
***********************************************************************************/

Matrix3x3::Matrix3x3(float m11, float m12, float m13, float m21, float m22, float m23,
	float m31, float m32, float m33)
{
	m[0] = m11;
	m[1] = m12;
	m[2] = m13;
	m[3] = m21;
	m[4] = m22;
	m[5] = m23;
	m[6] = m31;
	m[7] = m32;
	m[8] = m33;
}

/***********************************************************************************
*   Matrix3x3::identity
*
*   params :
*       void
*   returns :
*       void
*
*   Set the matrix to identity matrix
***********************************************************************************/

void Matrix3x3::identity(void)
{
	m[0] = 1.0;
	m[1] = 0.0;
	m[2] = 0.0;
	m[3] = 0.0;
	m[4] = 1.0;
	m[5] = 0.0;
	m[6] = 0.0;
	m[7] = 0.0;
	m[8] = 1.0;
}

/***********************************************************************************
*   Matrix3x3::neg
*
*   this = -this
***********************************************************************************/

void Matrix3x3::neg(void)
{
	m[0] = -m[0];
	m[1] = -m[1];
	m[2] = -m[2];
	m[3] = -m[3];
	m[4] = -m[4];
	m[5] = -m[5];
	m[6] = -m[6];
	m[7] = -m[7];
	m[8] = -m[8];
}

/***********************************************************************************
*   Matrix3x3::neg
*
*   params :
*       Matrix3x3 & A       -   input A
*
*   this = -A
***********************************************************************************/

void Matrix3x3::neg(Matrix3x3 & A)
{
	m[0] = -A.m[0];
	m[1] = -A.m[1];
	m[2] = -A.m[2];
	m[3] = -A.m[3];
	m[4] = -A.m[4];
	m[5] = -A.m[5];
	m[6] = -A.m[6];
	m[7] = -A.m[7];
	m[8] = -A.m[8];
}

/***********************************************************************
*  Matrix3x3::add
*
*   params :
*       Matrix3x3 & A       -   matrix to add
*   returns :
*       void
*
*  Add another matrix to this matrix :
*       this = this + A
***********************************************************************/

void Matrix3x3::add(Matrix3x3 & A)
{
	m[0] += A.m[0];
	m[1] += A.m[1];
	m[2] += A.m[2];
	m[3] += A.m[3];
	m[4] += A.m[4];
	m[5] += A.m[5];
	m[6] += A.m[6];
	m[7] += A.m[7];
	m[8] += A.m[8];
}

/***********************************************************************
*   Matrix3x3::mac
*
*   params :
*       Matrix3x3 & A       -   input A
*       Matrix3x3 & B       -   input B
*       float s             -   scalar
*
*   this = A + B * s
***********************************************************************/

void Matrix3x3::mac(Matrix3x3 & A, Matrix3x3 & B, float s)
{
	m[0] = A.m[0] + B.m[0] * s;
	m[1] = A.m[1] + B.m[1] * s;
	m[2] = A.m[2] + B.m[2] * s;
	m[3] = A.m[3] + B.m[3] * s;
	m[4] = A.m[4] + B.m[4] * s;
	m[5] = A.m[5] + B.m[5] * s;
	m[6] = A.m[6] + B.m[6] * s;
	m[7] = A.m[7] + B.m[7] * s;
	m[8] = A.m[8] + B.m[8] * s;
}

/***********************************************************************
*   Matrix3x3::mac
*
*   params :
*       Matrix3x3 & A       -   input A
*       float s             -   scalar
*
*   this = this + A * s
***********************************************************************/

void Matrix3x3::mac(Matrix3x3 & A, float s)
{
    m[0] += A.m[0] * s;
    m[1] += A.m[1] * s;
    m[2] += A.m[2] * s;
    m[3] += A.m[3] * s;
    m[4] += A.m[4] * s;
    m[5] += A.m[5] * s;
    m[6] += A.m[6] * s;
    m[7] += A.m[7] * s;
    m[8] += A.m[8] * s;
}

/***********************************************************************************
*   Matrix3x3::add
*
*   params :
*       Matrix3x3 & A       -   input A
*       Matrix3x3 & B       -   input B
*   returns :
*       void
*
*   Add two matricies into this matrix :
*       this = A + B
***********************************************************************************/

void Matrix3x3::add(Matrix3x3 & A, Matrix3x3 & B)
{
	m[0] = A.m[0] + B.m[0];
	m[1] = A.m[1] + B.m[1];
	m[2] = A.m[2] + B.m[2];
	m[3] = A.m[3] + B.m[3];
	m[4] = A.m[4] + B.m[4];
	m[5] = A.m[5] + B.m[5];
	m[6] = A.m[6] + B.m[6];
	m[7] = A.m[7] + B.m[7];
	m[8] = A.m[8] + B.m[8];
}

/***********************************************************************
*  Matrix3x3::sub
*
*   params :
*       Matrix3x3 & A       -   matrix to subtract
*   returns :
*       void
*
*  Subtract another matrix from this matrix :
*       this = this - A
***********************************************************************/

void Matrix3x3::sub(Matrix3x3 & A)
{
	m[0] -= A.m[0];
	m[1] -= A.m[1];
	m[2] -= A.m[2];
	m[3] -= A.m[3];
	m[4] -= A.m[4];
	m[5] -= A.m[5];
	m[6] -= A.m[6];
	m[7] -= A.m[7];
	m[8] -= A.m[8];
}

/***********************************************************************************
*   Matrix3x3::sub
*
*   params :
*       Matrix3x3 & A       -   input A
*       Matrix3x3 & B       -   input B
*   returns :
*       void
*
*   Add two matricies into this matrix :
*       this = A - B
***********************************************************************************/

void Matrix3x3::sub(Matrix3x3 & A, Matrix3x3 & B)
{
	m[0] = A.m[0] - B.m[0];
	m[1] = A.m[1] - B.m[1];
	m[2] = A.m[2] - B.m[2];
	m[3] = A.m[3] - B.m[3];
	m[4] = A.m[4] - B.m[4];
	m[5] = A.m[5] - B.m[5];
	m[6] = A.m[6] - B.m[6];
	m[7] = A.m[7] - B.m[7];
	m[8] = A.m[8] - B.m[8];
}

/***********************************************************************
*  Matrix3x3::mult
*
*   params :
*       float s        -       scalar to multiply by
*   returns :
*       void
*
*  Multiply this matrix by a scalar :
*       this = this * s
***********************************************************************/

void Matrix3x3::mult(float s)
{
	m[0] *= s;
	m[1] *= s;
	m[2] *= s;
	m[3] *= s;
	m[4] *= s;
	m[5] *= s;
	m[6] *= s;
	m[7] *= s;
	m[8] *= s;
}

/***********************************************************************
*  Matrix3x3::mult
*
*   params :
*         Matrix3x3 & A     -           input A
*       float s           -       scalar to multiply by
*   returns :
*       void
*
*  Multiply a matrix by a scalar into this matrix :
*       this = A * s
***********************************************************************/

void Matrix3x3::mult(Matrix3x3 & A, float s)
{
	m[0] = A.m[0] * s;
	m[1] = A.m[1] * s;
	m[2] = A.m[2] * s;
	m[3] = A.m[3] * s;
	m[4] = A.m[4] * s;
	m[5] = A.m[5] * s;
	m[6] = A.m[6] * s;
	m[7] = A.m[7] * s;
	m[8] = A.m[8] * s;
}

/***********************************************************************************
*   Matrix3x3::mult
*
*   params :
*       Matrix3x3 & A       -   input A
*       Matrix3x3 & B       -   input B
*
*   Multiply two matricies into this matrix :
*       this = A * B
***********************************************************************************/

void Matrix3x3::mult(Matrix3x3 & A, Matrix3x3 & B)
{
	m[0] = A.m[0] * B.m[0] + A.m[1] * B.m[3] + A.m[2] * B.m[6];
	m[1] = A.m[0] * B.m[1] + A.m[1] * B.m[4] + A.m[2] * B.m[7];
	m[2] = A.m[0] * B.m[2] + A.m[1] * B.m[5] + A.m[2] * B.m[8];
	m[3] = A.m[3] * B.m[0] + A.m[4] * B.m[3] + A.m[5] * B.m[6];
	m[4] = A.m[3] * B.m[1] + A.m[4] * B.m[4] + A.m[5] * B.m[7];
	m[5] = A.m[3] * B.m[2] + A.m[4] * B.m[5] + A.m[5] * B.m[8];
	m[6] = A.m[6] * B.m[0] + A.m[7] * B.m[3] + A.m[8] * B.m[6];
	m[7] = A.m[6] * B.m[1] + A.m[7] * B.m[4] + A.m[8] * B.m[7];
	m[8] = A.m[6] * B.m[2] + A.m[7] * B.m[5] + A.m[8] * B.m[8];
}

/***********************************************************************************
*   Matrix3x3::multAtB
*
*   params :
*       Matrix3x3 & A       -   input A
*       Matrix3x3 & B       -   input B
*
*   Multiply two matricies into this matrix, transposing A :
*       this = transpose(A) * B
***********************************************************************************/

void Matrix3x3::multAtB(Matrix3x3 & A, Matrix3x3 & B)
{
	m[0] = A.m[0] * B.m[0] + A.m[3] * B.m[3] + A.m[6] * B.m[6];
	m[1] = A.m[0] * B.m[1] + A.m[3] * B.m[4] + A.m[6] * B.m[7];
	m[2] = A.m[0] * B.m[2] + A.m[3] * B.m[5] + A.m[6] * B.m[8];
	m[3] = A.m[1] * B.m[0] + A.m[4] * B.m[3] + A.m[7] * B.m[6];
	m[4] = A.m[1] * B.m[1] + A.m[4] * B.m[4] + A.m[7] * B.m[7];
	m[5] = A.m[1] * B.m[2] + A.m[4] * B.m[5] + A.m[7] * B.m[8];
	m[6] = A.m[2] * B.m[0] + A.m[5] * B.m[3] + A.m[8] * B.m[6];
	m[7] = A.m[2] * B.m[1] + A.m[5] * B.m[4] + A.m[8] * B.m[7];
	m[8] = A.m[2] * B.m[2] + A.m[5] * B.m[5] + A.m[8] * B.m[8];
}

/***********************************************************************************
*   Matrix3x3::multABt
*
*   params :
*       Matrix3x3 & A       -   input A
*       Matrix3x3 & B       -   input B
*
*   Multiply two matricies into this matrix, transposing B:
*       this = A * transpose(B)
***********************************************************************************/

void Matrix3x3::multABt(Matrix3x3 & A, Matrix3x3 & B)
{
	m[0] = A.m[0] * B.m[0] + A.m[1] * B.m[1] + A.m[2] * B.m[2];
	m[1] = A.m[0] * B.m[3] + A.m[1] * B.m[4] + A.m[2] * B.m[5];
	m[2] = A.m[0] * B.m[6] + A.m[1] * B.m[7] + A.m[2] * B.m[8];
	m[3] = A.m[3] * B.m[0] + A.m[4] * B.m[1] + A.m[5] * B.m[2];
	m[4] = A.m[3] * B.m[3] + A.m[4] * B.m[4] + A.m[5] * B.m[5];
	m[5] = A.m[3] * B.m[6] + A.m[4] * B.m[7] + A.m[5] * B.m[8];
	m[6] = A.m[6] * B.m[0] + A.m[7] * B.m[1] + A.m[8] * B.m[2];
	m[7] = A.m[6] * B.m[3] + A.m[7] * B.m[4] + A.m[8] * B.m[5];
	m[8] = A.m[6] * B.m[6] + A.m[7] * B.m[7] + A.m[8] * B.m[8];
}

/***********************************************************************
*  Matrix3x3::transpose
*
*   params :
*       void
*   returns :
*       void
*
*  Transpose the matrix in place :
*       this = trans(this)
***********************************************************************/

void Matrix3x3::transpose(void)
{
	float t;
	t = m[1];
	m[1] = m[3];
	m[3] = t;
	t = m[2];
	m[2] = m[6];
	m[6] = t;
	t = m[5];
	m[5] = m[7];
	m[7] = t;
}

/***********************************************************************************
*   Matrix3x3::transpose
*
*   params :
*       Matrix3x3 & A       -   input A
*   returns :
*       void
*
*   Transpose the matrix A into another matrix :
*       this = trans(A)
***********************************************************************************/

void Matrix3x3::transpose(Matrix3x3 & A)
{
	m[0] = A.m[0];
   m[1] = A.m[3];
   m[2] = A.m[6];
   m[3] = A.m[1];
   m[4] = A.m[4];
   m[5] = A.m[7];
   m[6] = A.m[2];
   m[7] = A.m[5];
   m[8] = A.m[8];
}

/***********************************************************************************
*   Matrix3x3::det
*
*   params :
*       void
*   returns :
*       float              -   determinant of matrix
*
*   Find the determinant of this matrix.
***********************************************************************************/

float Matrix3x3::det(void)
{
	return m[0] * (m[4] * m[8] - m[5] * m[7])
		- m[1] * (m[3] * m[8] - m[5] * m[6])
		+ m[2] * (m[3] * m[7] - m[4] * m[6]);
}

/***********************************************************************************
*   Matrix3x3::inverse
*
*   params :
*       Matrix3x3 & M       -   input M
*   returns :
*       int                 -   0 if OK, 20 if singular
*
*   Make this matrix the inverse of the given matrix using determinants :
*       this = inv(M)
***********************************************************************************/

int Matrix3x3::inverse(Matrix3x3 & M)
{
	float detM = M.det();
	if (fabs(detM) < 1e-20)                     // cant invert this
		return 20;

   float oneOverDetM = 1.0 / detM;     // one over det M
	m[0] = oneOverDetM * (M.m[4] * M.m[8] - M.m[5] * M.m[7]);
	m[1] = -oneOverDetM * (M.m[1] * M.m[8] - M.m[2] * M.m[7]);
	m[2] = oneOverDetM * (M.m[1] * M.m[5] - M.m[2] * M.m[4]);
	m[3] = -oneOverDetM * (M.m[3] * M.m[8] - M.m[5] * M.m[6]);
	m[4] = oneOverDetM * (M.m[0] * M.m[8] - M.m[2] * M.m[6]);
	m[5] = -oneOverDetM * (M.m[0] * M.m[5] - M.m[2] * M.m[3]);
	m[6] = oneOverDetM * (M.m[3] * M.m[7] - M.m[4] * M.m[6]);
	m[7] = -oneOverDetM * (M.m[0] * M.m[7] - M.m[1] * M.m[6]);
	m[8] = oneOverDetM * (M.m[0] * M.m[4] - M.m[1] * M.m[3]);

   return 0;
}

/***********************************************************************************
*   Matrix3x3::exp
*
*   params :
*       Matrix3x3 & A       -   input A
*   returns :
*       void
*
*   Make this the exponential of matrix A, using a power series :
*       this = exp(A)
*
*   The constant numTerms is how many terms in the series.
*    This is SLOW so beware.
***********************************************************************************/

void Matrix3x3::exp(Matrix3x3 & A)
{
	const int numTerms = 100;

	Matrix3x3 term;                        // next term in series
   Matrix3x3 temp;

	identity();                                     // clear this matrix to I
	term.identity();

	for (int div = 1; div <= numTerms; div++) {
		temp.mult(term, A);
		term.mult(temp, float(1.0 / div));
													// find next term = term * A / div
		add(term);
	}
}

/***********************************************************************
*   Matrix3x3::xAxisRotation
*
*   params :
*       float thetaX       -       angle (radians)
*   returns :
*       none
*
*   Set up the matrix as a rotation about X axis
*       [ 1    0   0 ]
*   M = [ 0    c  -s ]
*       [ 0    s   c ]
***********************************************************************/

void Matrix3x3::xAxisRotation(float thetaX)
{
	m[4] = m[8] = cos(thetaX);
	m[7] = sin(thetaX);
	m[5] = -m[7];
	m[1] = m[2] = m[3] = m[6] = 0.0;
	m[0] = 1.0;
}

/***********************************************************************
*   Matrix3x3::yAxisRotation
*
*   params :
*       float thetaY       -       angle (radians)
*   returns :
*       none
*
*   Set up the matrix as a rotation about Y axis
*       [  c   0   s ]
*   M = [  0   1   0 ]
*       [ -s   0   c ]
***********************************************************************/

void Matrix3x3::yAxisRotation(float thetaY)
{
	m[0] = m[8] = cos(thetaY);
	m[2] = sin(thetaY);
	m[6] = -m[2];
	m[1] = m[3] = m[5] = m[7] = 0.0;
	m[4] = 1.0;
}

/***********************************************************************
*   Matrix3x3::zAxisRotation
*
*   params :
*       float thetaZ       -       angle (radians)
*   returns :
*       none
*
*   Set up the matrix as a rotation about Z axis.
*       [ c   -s   0 ]
*   M = [ s    c   0 ]
*       [ 0    0   1 ]
***********************************************************************/

void Matrix3x3::zAxisRotation(float thetaZ)
{
	m[0] = m[4] = cos(thetaZ);
	m[3] = sin(thetaZ);
	m[1] = -m[3];
	m[2] = m[5] = m[6] = m[7] = 0.0;
	m[8] = 1.0;
}

/***********************************************************************
*   Matrix3x3::arbAxisRotation
*
*   params :
*       Vector3 & axis      -   axis to rotate around.
*       float theta        -   angle (radians), clockwise when axis goes away
*
*   Set up the matrix as a rotation about a given axis.
***********************************************************************/

void Matrix3x3::arbAxisRotation(Vector3 & axis, float theta)
{
    Quaternion Q;                                   // use a quaternion to make the matrix
    Q.rotation(axis, theta);
    fromQuat(Q);
}

/***********************************************************************************
*   Matrix3x3::rpyRotation
*
*   params :
*       float roll     -       roll (radians)
*       float pitch    -       pitch (radians)
*       float yaw      -       yaw (radians)
*   returns :
*       none
*
*   Set matrix as 3 concatenated rotations using roll, pitch and yaw angles. The order is :
*
*       roll around global X axis
*       pitch around global Y axis
*       yaw around global Z axis
*
*   which is equivalent to :
*
*       yaw around object Z axis
*       pitch around object Y axis
*       roll around object X axis
***********************************************************************************/

void Matrix3x3::rpyRotation(float roll, float pitch, float yaw)
{
	Matrix3x3 rollM, pitchM, yawM;
   rollM.xAxisRotation(roll);
	pitchM.yAxisRotation(pitch);
   yawM.zAxisRotation(yaw);

	Matrix3x3 temp;
    temp.mult(yawM, pitchM);
	mult(temp, rollM);
														// total matrix is yawM * pitchM * rollM
}

/***********************************************************************************
*   Matrix3x3::star
*
*   params :
*       Vector3 & a             -   input vector
*
*  Set matrix as the * or ~ operator on the vector :
*
*          [  0   -z    y ]
*   this = [  z    0   -x ]
*          [ -y    x    0 ]
***********************************************************************************/

void Matrix3x3::star(Vector3 & a)
{
	m[0] = m[4] = m[8] = 0.0;
    m[3] = a.z;
    m[1] = -m[3];
    m[2] = a.y;
    m[6] = -m[2];
    m[7] = a.x;
    m[5] = -m[7];
}

/***********************************************************************************
*   Vector3::Vector3
*
*   params :
*       float x, y, z      -       initial value
*   returns :
*       none
*
*   Constructor with initial value.
***********************************************************************************/

Vector3::Vector3(float x, float y, float z)
{
	this->x = x;
	this->y = y;
	this->z = z;
};

/***********************************************************************************
*   Vector3::add
*
*   params :
*       Vector3 & A     -       input A
*       Vector3 & B     -       input B
*   returns :
*       void
*
*   Adds 2 vectors into this vector :
*       this = A + B
***********************************************************************************/

void Vector3::add(Vector3 & A, Vector3 & B)
{
	x = A.x + B.x;
   y = A.y + B.y;
   z = A.z + B.z;
}

/***********************************************************************************
*   Vector3::add
*
*   params :
*       Vector3 & A     -       input A
*   returns :
*       void
*
*   Adds a vector into this vector :
*       this = this + A
***********************************************************************************/

void Vector3::add(Vector3 & A)
{
	x += A.x;
   y += A.y;
   z += A.z;
}

/***********************************************************************************
*   Vector3::mac
*
*   params :
*       Vector3 & A, & B                 -  input A, B
*       float scalar                     -  mult by this
*
*   this = A + B * scalar
***********************************************************************************/

void Vector3::mac(Vector3 & A, Vector3 & B, float scalar)
{
    x = A.x + B.x * scalar;
    y = A.y + B.y * scalar;
    z = A.z + B.z * scalar;
}

/***********************************************************************************
*   Vector3::mac
*
*   params :
*       Vector3 & A                      -  input A
*       float scalar                     -  mult by this
*
*   this = this + A * scalar
***********************************************************************************/

void Vector3::mac(Vector3 & A, float scalar)
{
    x += A.x * scalar;
    y += A.y * scalar;
    z += A.z * scalar;
}

/***********************************************************************************
*   Vector3::msc
*
*   params :
*       Vector3 & A, & B                 -  input A, B
*       float scalar                     -  mult by this
*
*   this = A - B * scalar
***********************************************************************************/

void Vector3::msc(Vector3 & A, Vector3 & B, float scalar)
{
    x = A.x - B.x * scalar;
    y = A.y - B.y * scalar;
    z = A.z - B.z * scalar;
}

/***********************************************************************************
*   Vector3::msc
*
*   params :
*       Vector3 & A                      -  input A
*       float scalar                     -  mult by this
*
*   this = this - A * scalar
***********************************************************************************/

void Vector3::msc(Vector3 & A, float scalar)
{
    x -= A.x * scalar;
    y -= A.y * scalar;
    z -= A.z * scalar;
}

/***********************************************************************************
*   Vector3::mac2
*
*   params :
*       Vector3 & A                      -  input A
*       Vector3 & B                      -  input B
*       float scalarB                    -  mult B by this
*       Vector3 & C                      -  input C
*       float scalarC                    -  mult C by this
*
*   this = A + B * scalarB + C * scalarC
***********************************************************************************/

void Vector3::mac2(Vector3 & A, Vector3 & B, float scalarB, Vector3 & C, float scalarC)
{
    x = A.x + B.x * scalarB + C.x * scalarC;
    y = A.y + B.y * scalarB + C.y * scalarC;
    z = A.z + B.z * scalarB + C.z * scalarC;
}

/***********************************************************************************
*   Vector3::msc2
*
*   params :
*       Vector3 & A                      -  input A
*       Vector3 & B                      -  input B
*       float scalarB                    -  mult B by this
*       Vector3 & C                      -  input C
*       float scalarC                    -  mult C by this
*
*   this = A - B * scalarB - C * scalarC
***********************************************************************************/

void Vector3::msc2(Vector3 & A, Vector3 & B, float scalarB, Vector3 & C, float scalarC)
{
    x = A.x - B.x * scalarB - C.x * scalarC;
    y = A.y - B.y * scalarB - C.y * scalarC;
    z = A.z - B.z * scalarB - C.z * scalarC;
}

/***********************************************************************************
*   Vector3::lerp
*
*   params :
*       Vector3 & A, & B                -  input A, B
*       float t                         -  fraction between A and B
*
*   Linear interpolate between two vectors :
*       this = A + t * (B - A)
***********************************************************************************/

void Vector3::lerp(Vector3 & A, Vector3 & B, float t)
{
    x = A.x + t * (B.x - A.x);
    y = A.y + t * (B.y - A.y);
    z = A.z + t * (B.z - A.z);
}

/***********************************************************************************
*   Vector3::herp
*
*   params :
*       Vector3 & p0, & p1, & p2, & p3      -   4 control points
*       float t                             -  fraction between p1 and p2
*
*   Hermite interpolate between p1 and p2. p0 and p3 are used for finding gradient at p1 and p2.
*       this = p0 * (2t^2 - t^3 - t)/2
*            + p1 * (3t^3 - 5t^2 + 2)/2
*            + p2 * (4t^2 - 3t^3 + t)/2
*            + p3 * (t^3 - t^2)/2
***********************************************************************************/

void Vector3::herp(Vector3 & p0, Vector3 & p1, Vector3 & p2, Vector3 & p3, float t)
{
    float t2 = t * t;
    float t3 = t2 * t;
    float kp0 = (2.0 * t2 - t3 - t) / 2.0;
    float kp1 = (3.0 * t3 - 5.0 * t2 + 2.0) / 2.0;
    float kp2 = (4.0 * t2 - 3.0 * t3 + t) / 2.0;
    float kp3 = (t3 - t2) / 2.0;
    x = p0.x * kp0 + p1.x * kp1 + p2.x * kp2 + p3.x * kp3;
    y = p0.y * kp0 + p1.y * kp1 + p2.y * kp2 + p3.y * kp3;
    z = p0.z * kp0 + p1.z * kp1 + p2.z * kp2 + p3.z * kp3;
}

/***********************************************************************************
*   Vector3::dist
*
*   params :
*       Vector3 & A                         -   input A
*   returns :
*       float                               -   distance this to A
*
*   Find distance between this and A.
***********************************************************************************/

float Vector3::dist(Vector3 & A)
{
    float dx = A.x - x;
    float dy = A.y - y;
    float dz = A.z - z;
    return sqrt(dx * dx + dy * dy + dz * dz);
}

/***********************************************************************************
*   Vector3::dist2
*
*   params :
*       Vector3 & A                         -   input A
*   returns :
*       float                               -   distance (this to A)^2
*
*   Find distance^2 between this and A (faster).
***********************************************************************************/

float Vector3::dist2(Vector3 & A)
{
    float dx = A.x - x;
    float dy = A.y - y;
    float dz = A.z - z;
    return (dx * dx + dy * dy + dz * dz);
}

/***********************************************************************************
*   Vector3::sub
*
*   params :
*       Vector3 & A     -       input A
*       Vector3 & B     -       input B
*   returns :
*       void
*
*   Subtracts 2 vectors into this vector :
*       this = A - B
***********************************************************************************/

void Vector3::sub(Vector3 & A, Vector3 & B)
{
	x = A.x - B.x;
   y = A.y - B.y;
   z = A.z - B.z;
}

/***********************************************************************************
*   Vector3::sub
*
*   params :
*       Vector3 & A     -       input A
*   returns :
*       void
*
*   Subtracts a vector from this vector :
*       this = this - A
***********************************************************************************/

void Vector3::sub(Vector3 & A)
{
	x -= A.x;
   y -= A.y;
   z -= A.z;
}

/***********************************************************************************
*   Vector3::neg
*
*   params :
*       void
*   returns :
*       void
*
*  Negate this vector :
*       this = -this
***********************************************************************************/

void Vector3::neg(void)
{
   x = -x;
   y = -y;
   z = -z;
}

/***********************************************************************************
*   Vector3::neg
*
*   params :
*       Vector & A      -   input A
*   returns :
*       void
*
*  Negate input into this vector :
*       this = -A
***********************************************************************************/

void Vector3::neg(Vector3 & A)
{
   x = -A.x;
   y = -A.y;
   z = -A.z;
}

/***********************************************************************************
*   Vector3::mult
*
*   params :
*       Matrix3x3 & M   -       matrix M
*       Vector3 & A     -       input A
*   returns :
*       void
*
*   Pre-Multiply the column vector A by the matrix M :
*       this = M * A
***********************************************************************************/

void Vector3::mult(Matrix3x3 & M, Vector3 & A)
{
	x = A.x * M.m[0] + A.y * M.m[1] + A.z * M.m[2];
	y = A.x * M.m[3] + A.y * M.m[4] + A.z * M.m[5];
	z = A.x * M.m[6] + A.y * M.m[7] + A.z * M.m[8];
}

/***********************************************************************************
*   Vector3::mult2
*
*   params :
*       Matrix3x3 & M1  -       matrix M1
*       Vector3 & A1    -       input A1
*       Matrix3x3 & M2  -       matrix M2
*       Vector3 & A2    -       input A2
*
*   this = M1 * A1 + M2 * A2
***********************************************************************************/

void Vector3::mult2(Matrix3x3 & M1, Vector3 & A1, Matrix3x3 & M2, Vector3 & A2)
{
	x = A1.x * M1.m[0] + A1.y * M1.m[1] + A1.z * M1.m[2] + A2.x * M2.m[0] + A2.y * M2.m[1] + A2.z * M2.m[2];
	y = A1.x * M1.m[3] + A1.y * M1.m[4] + A1.z * M1.m[5] + A2.x * M2.m[3] + A2.y * M2.m[4] + A2.z * M2.m[5];
	z = A1.x * M1.m[6] + A1.y * M1.m[7] + A1.z * M1.m[8] + A2.x * M2.m[6] + A2.y * M2.m[7] + A2.z * M2.m[8];
}

/***********************************************************************************
*   Vector3::mac
*
*   params :
*       Matrix3x3 & M   -       matrix M
*       Vector3 & A     -       input A
*
*   this = this + M * A
***********************************************************************************/

void Vector3::mac(Matrix3x3 & M, Vector3 & A)
{
	x += A.x * M.m[0] + A.y * M.m[1] + A.z * M.m[2];
	y += A.x * M.m[3] + A.y * M.m[4] + A.z * M.m[5];
	z += A.x * M.m[6] + A.y * M.m[7] + A.z * M.m[8];
}

/***********************************************************************************
*   Vector3::transMult
*
*   params :
*       Matrix3x3 & M   -       matrix M
*       Vector3 & A     -       input A
*   returns :
*       void
*
*   Pre-Multiply the column vector A by the matrix transpose(M) :
*       this = transpose(M) * A
***********************************************************************************/

void Vector3::transMult(Matrix3x3 & M, Vector3 & A)
{
	x = A.x * M.m[0] + A.y * M.m[3] + A.z * M.m[6];
	y = A.x * M.m[1] + A.y * M.m[4] + A.z * M.m[7];
	z = A.x * M.m[2] + A.y * M.m[5] + A.z * M.m[8];
}

/***********************************************************************************
*   Vector3::mult
*
*   params :
*       Vector3 & A     -       input A
*       float scalar   -       scalar to multiply by
*   returns :
*       void
*
*   Multiply the vector A by the scalar into this :
*       this = A * scalar
***********************************************************************************/

void Vector3::mult(Vector3 & A, float scalar)
{
	x = A.x * scalar;
   y = A.y * scalar;
   z = A.z * scalar;
}

/***********************************************************************************
*   Vector3::mult
*
*   params :
*       float scalar   -       scalar to multiply by
*   returns :
*       void
*
*   Multiply this vector by the scalar :
*       this = this * scalar
***********************************************************************************/

void Vector3::mult(float scalar)
{
	x *= scalar;
   y *= scalar;
   z *= scalar;
}

/***********************************************************************************
*   Vector3::normalise
*
*   params :
*       Vector3 & A     -       input A
*   returns :
*       int                 -           0 if OK, 20 if zero length.
*
*   Normalise the given vector to length 1.0 into this :
*       this = norm(A)
***********************************************************************************/

int Vector3::normalise(Vector3 & A)
{
    float modA2 = A.x * A.x + A.y * A.y + A.z * A.z;
    if (modA2 < 1e-10)                          // zero length vector
		return 20;
    float oneOverMod = 1.0 / sqrt(modA2);
    x = A.x * oneOverMod;
    y = A.y * oneOverMod;
    z = A.z * oneOverMod;
    return 0;
}

/***********************************************************************************
*   Vector3::normalise
*
*   params :
*       void
*   returns :
*       int                 -           0 if OK, 20 if zero length.
*
*   Normalise the vector to length 1.0 :
*       this = norm(this)
***********************************************************************************/

int Vector3::normalise(void)
{
    float modA2 = x * x + y * y + z * z;
    if (modA2 < 1e-10)                          // zero length vector
		return 20;
    float oneOverMod = 1.0 / sqrt(modA2);
    x *= oneOverMod;
    y *= oneOverMod;
    z *= oneOverMod;
    return 0;
}

/***********************************************************************************
*   Vector3::mod
*
*   params :
*       void
*   returns :
*       float          -       length of this
*
*   Find modulus (length) of this vector.
***********************************************************************************/

float Vector3::mod(void)
{
	return sqrt(x * x + y * y + z * z);
}

/***********************************************************************************
*   Vector3::mod2
*
*   params :
*       void
*   returns :
*       float          -       length of this ^ 2
*
*   Find modulus ^ 2 of this vector.
***********************************************************************************/

float Vector3::mod2(void)
{
	return x * x + y * y + z * z;
}

/***********************************************************************************
*   Vector3::dot
*
*   params :
*       Vector3 & A     -       input A
*   returns :
*       float          -       this.A
*
*   Dot product of this vector and A.
***********************************************************************************/

float Vector3::dot(Vector3 & A)
{
	return A.x * x + A.y * y + A.z * z;
}

/***********************************************************************************
*   Vector3::cross
*
*   params :
*       Vector3 & A     -       input A
*       Vector3 & B     -       input B
*   returns :
*       void
*
*   Vector product of A and B into this vector :
*       this = A x B
***********************************************************************************/

void Vector3::cross(Vector3 & A, Vector3 & B)
{
	x = A.y * B.z - A.z * B.y;
	y = A.z * B.x - A.x * B.z;
	z = A.x * B.y - A.y * B.x;
}

/***********************************************************************
*  Matrix3x3::fromQuat
*
*   params :
*       Quaternion & q                  -   quaternion to make into matrix
*
*   Make the matrix a rotation matrix from the quaternion [x y z w].
*   w = cos(theta / 2)
*   [x y z] = axis * sin(theta / 2)
*
*   This assumes the quaternion is normalised. If not the matrix will not be
*   orthogonal.
***********************************************************************/

void Matrix3x3::fromQuat(Quaternion & q)
{
    float x2, y2, z2;
    x2 = q.x * q.x;
    y2 = q.y * q.y;
    z2 = q.z * q.z;
    float xy, wz, xz, wy, yz, wx;
    xy = q.x * q.y;
    wz = q.w * q.z;
    xz = q.x * q.z;
    wy = q.w * q.y;
    yz = q.y * q.z;
    wx = q.w * q.x;

    m[0] = 1.0 - 2.0 * (y2 + z2);
    m[1] = 2.0 * (xy - wz);
    m[2] = 2.0 * (xz + wy);
    m[3] = 2.0 * (xy + wz);
    m[4] = 1.0 - 2.0 * (z2 + x2);
    m[5] = 2.0 * (yz - wx);
    m[6] = 2.0 * (xz - wy);
    m[7] = 2.0 * (yz + wx);
    m[8] = 1.0 - 2.0 * (x2 + y2);
}

/***********************************************************************
*  Matrix3x3::fromQuatL2
*
*   params :
*       Quaternion & q                  -   quaternion to make into matrix
*       float l2                        -   length^2 of quaternion
*
*   Make the matrix a rotation matrix from the quaternion [x y z w].
*   w = cos(theta / 2)
*   [x y z] = axis * sin(theta / 2)
*
*   Will always be orthogonal even if quaternion not normalised because the
*   matrix is automatically scaled by l2.
***********************************************************************/

void Matrix3x3::fromQuatL2(Quaternion & q, float l2)
{
    float x2, y2, z2;
    x2 = q.x * q.x;
    y2 = q.y * q.y;
    z2 = q.z * q.z;
    float xy, wz, xz, wy, yz, wx;
    xy = q.x * q.y;
    wz = q.w * q.z;
    xz = q.x * q.z;
    wy = q.w * q.y;
    yz = q.y * q.z;
    wx = q.w * q.x;

    float scale = 2.0 / l2;                         // scale to make it orthogonal
    m[0] = 1.0 - scale * (y2 + z2);
    m[1] = scale * (xy - wz);
    m[2] = scale * (xz + wy);
    m[3] = scale * (xy + wz);
    m[4] = 1.0 - scale * (z2 + x2);
    m[5] = scale * (yz - wx);
    m[6] = scale * (xz - wy);
    m[7] = scale * (yz + wx);
    m[8] = 1.0 - scale * (x2 + y2);
}

/***********************************************************************
*   Quaternion::Quaternion
*
*   params :
*       float w, x, y, z                -   initial value
*
*   Constructor with initial value.
***********************************************************************/

Quaternion::Quaternion(float w, float x, float y, float z)
{
    this->w = w;
    this->x = x;
    this->y = y;
    this->z = z;
}

/***********************************************************************
*   Quaternion::rotation
*
*   params :
*       Vector3 & axis              -   rotation axis
*       float angle                 -   angle (radians)
*
*   Make this quaternion a rotation. +ve rotation is clockwise around axis when axis points away.
*   Quaternion is :
*   w = cos(theta / 2)
*   [x y z] = axis * sin(theta / 2)
***********************************************************************/

void Quaternion::rotation(Vector3 & axis, float angle)
{
    Vector3 nAxis;
    nAxis.normalise(axis);
    float angle2 = angle / 2.0;
    w = cos(angle2);
    float s = sin(angle2);
    x = nAxis.x * s;
    y = nAxis.y * s;
    z = nAxis.z * s;
}

/***********************************************************************
*   Quaternion::mult
*
*   params :
*       Quaternion & A, & B         -   multiply these
*
*   this = A * B
*   Note : Rotation by q1 followed by q2 gives a total rotation of (q2 * q1).
***********************************************************************/

void Quaternion::mult(Quaternion & A, Quaternion & B)
{
    w = A.w * B.w - A.x * B.x - A.y * B.y - A.z * B.z;
    x = A.w * B.x + A.x * B.w + A.y * B.z - A.z * B.y;
    y = A.w * B.y - A.x * B.z + A.y * B.w + A.z * B.x;
    z = A.w * B.z + A.x * B.y - A.y * B.x + A.z * B.w;
}

/***********************************************************************
*   Quaternion::fromMatrix
*
*   params :
*       Matrix3x3 & M       -   matrix to make into quaternion
*
*   Find quaternion from 3x3 matrix with numerically stable method.
***********************************************************************/

void Quaternion::fromMatrix(Matrix3x3 & M)
{
    float qw2 = 0.25 * (M.m[0] + M.m[4] + M.m[8] + 1.0);
    float qx2 = qw2 - 0.5 * (M.m[4] + M.m[8]);
    float qy2 = qw2 - 0.5 * (M.m[8] + M.m[0]);
    float qz2 = qw2 - 0.5 * (M.m[0] + M.m[4]);      // quaternion components squared
    int i = (qw2 > qx2 ) ?                          // maximum magnitude component
        ((qw2 > qy2) ? ((qw2 > qz2) ? 0 : 3) : ((qy2 > qz2) ? 2 : 3)) :
        ((qx2 > qy2) ? ((qx2 > qz2) ? 1 : 3) : ((qy2 > qz2) ? 2 : 3));
    float tmp;

// compute signed quaternion components using numerically stable method
    switch(i) {
        case 0:
            w = sqrt(qw2);
            tmp = 0.25 / w;
            x = (M.m[7] - M.m[5]) * tmp;
            y = (M.m[2] - M.m[6]) * tmp;
            z = (M.m[3] - M.m[1]) * tmp;
            break;
        case 1:
            x = sqrt(qx2);
            tmp = 0.25 / x;
            w = (M.m[7] - M.m[5]) * tmp;
            y = (M.m[1] + M.m[3]) * tmp;
            z = (M.m[6] + M.m[2]) * tmp;
            break;
        case 2:
            y = sqrt(qy2);
            tmp = 0.25 / y;
            w = (M.m[2] - M.m[6]) * tmp;
            x = (M.m[1] + M.m[3]) * tmp;
            z = (M.m[5] + M.m[7]) * tmp;
            break;
        case 3:
            z = sqrt(qz2);
            tmp = 0.25 / z;
            w = (M.m[3] - M.m[1]) * tmp;
            x = (M.m[2] + M.m[6]) * tmp;
            y = (M.m[5] + M.m[7]) * tmp;
            break;
    }

    if (i && w < 0.0) {                             // make w always +ve
        w = -w;
        x = -x;
        y = -y;
        z = -z;
    }

    tmp = 1.0 / sqrt(w * w + x * x + y * y + z * z);
    w *= tmp;                                      // normalise it to be safe
    x *= tmp;
    y *= tmp;
    z *= tmp;
}

/***********************************************************************
*  eigenvectors
*
*   params :
*       Matrix3x3 & M                   -   matrix to find eigen vectors for (symmetric)
*       Matrix3x3 & V                   -   unit length eigen vectors returned here (in columns)
*       Vector3 & d                     -   eigen values returned here
*   returns :
*       int                             -   0 OK, 20 error
*
*   Find the eigen vectors and values for a matrix using the Jacobi rotation method.
*   Assumes the matrix is symmetric, i.e. the bottom half is ignored.
***********************************************************************/

int Matrix3x3::eigenvectors(Matrix3x3 & M, Matrix3x3 & V, Vector3 & d)
{
    Matrix3x3 M2 = M;                       // make a copy of M (this gets trashed)
    int j, iq, ip, i;
	float tresh, theta, tau, t, sm, s, h, g, c, b[3], z[3];

    V.identity();
	for (ip=0;ip<3;ip++) {
		b[ip]=d.m[ip]=M2.n[ip][ip];
	    z[ip]=0.0;
    }

	for (i=0;i<50;i++) {
        sm=0.0;
		for (ip=0;ip<2;ip++) {
		    for (iq=ip+1;iq<3;iq++)
	            sm += M2.n[ip][iq];
        }
        if (sm == 0.0) {
		    return 0;
        }
		if (i < 3)
		    tresh=0.2*fabs(sm)/(3*3);
		else
	        tresh=0.0;
		for (ip=0;ip<2;ip++) {
		    for (iq=ip+1;iq<3;iq++) {
		        g=100.0*fabs(M2.n[ip][iq]);

                if (i > 4 && (float)(fabs(d.m[ip])+g) == (float)fabs(d.m[ip])
                    && (float)(fabs(d.m[iq])+g) == (float)fabs(d.m[iq]))
                    M2.n[ip][iq]=0.0;

            	else if (fabs(M2.n[ip][iq]) > tresh) {
                	h=d.m[iq]-d.m[ip];
                    if ((float)(fabs(h)+g) == (float)fabs(h))
                    	t=(M2.n[ip][iq])/h;
                    else {
                    	theta=0.5*h/(M2.n[ip][iq]);
                        t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
                    	if (theta < 0.0) t = -t;
                    }
	                c=1.0/sqrt(1+t*t);
                	s=t*c;
	                tau=s/(1.0+c);
                    h=t*M2.n[ip][iq];
                	z[ip] -= h;
                	z[iq] += h;
                	d.m[ip] -= h;
                    d.m[iq] += h;
                	M2.n[ip][iq] = 0.0;
	                for (j=0;j<=ip-1;j++) {
                    	EIGENROTATE(M2,j,ip,j,iq)
                    }
                	for (j=ip+1;j<=iq-1;j++) {
                    	EIGENROTATE(M2,ip,j,j,iq)
                    }
                	for (j=iq+1;j<3;j++) {
                    	EIGENROTATE(M2,ip,j,iq,j)
                    }
                   	for (j=0;j<3;j++) {
                    	EIGENROTATE(V,j,ip,j,iq)
                    }
			    }
            }
        }
    	for (ip=0;ip<3;ip++) {
        	b[ip] += z[ip];
	        d.m[ip] = b[ip];
	        z[ip]=0.0;
        }
    }
    return 20;
}

/***********************************************************************************
*   Quaternion::normalise
*
*   Normalise to length 1.0 :
*       this = norm(this)
***********************************************************************************/

void Quaternion::normalise(void)
{
    float oneOverMod = 1.0 / sqrt(x * x + y * y + z * z + w * w);
    w *= oneOverMod;
    x *= oneOverMod;
    y *= oneOverMod;
    z *= oneOverMod;
}

/***********************************************************************************
*   Quaternion::mod
*
*   returns :
*       float                           -   length of it
***********************************************************************************/

float Quaternion::mod(void)
{
    return sqrt(x * x + y * y + z * z + w * w);
}

/***********************************************************************************
*   Quaternion::mod2
*
*   returns :
*       float                           -   length of it ^ 2
***********************************************************************************/

float Quaternion::mod2(void)
{
    return x * x + y * y + z * z + w * w;
}

/***********************************************************************************
*   Vector3::transform
*
*   params :
*       Vector3 & r                     -   point in coord system 1
*       Matrix3x3 & rotPos              -   rotation from system 1 to 2
*       Vector3 & linPos                -   translation from system 1 to 2
*
*   Transform a point to a new coords system.
*
*   this = rotPos * r + linPos
***********************************************************************************/

void Vector3::transform(Vector3 & r, Matrix3x3 & rotPos, Vector3 & linPos)
{
	x = r.x * rotPos.m[0] + r.y * rotPos.m[1] + r.z * rotPos.m[2] + linPos.x;
	y = r.x * rotPos.m[3] + r.y * rotPos.m[4] + r.z * rotPos.m[5] + linPos.y;
	z = r.x * rotPos.m[6] + r.y * rotPos.m[7] + r.z * rotPos.m[8] + linPos.z;
}

/***********************************************************************************
*   Vector3::transformArray
*
*   params :
*       Vector3 *source                 -   points in coord system 1
*       Vector3 *dest                   -   points in coord system 2
*       int N                           -   number of points
*       Matrix3x3 & rotPos              -   rotation from system 1 to 2
*       Vector3 & linPos                -   translation from system 1 to 2
*
*   Transform a list of points to a new coords system.
*
*   dest[i] = rotPos * source[i] + linPos
***********************************************************************************/

void Vector3::transformArray(Vector3 *source, Vector3 *dest, int N, Matrix3x3 & rotPos, Vector3 & linPos)
{
    float sx, sy, sz;
    for (int i = N; i--;) {
        sx = source->x;
        sy = source->y;
        sz = source->z;
    	dest->x = sx * rotPos.m[0] + sy * rotPos.m[1] + sz * rotPos.m[2] + linPos.x;
	    dest->y = sx * rotPos.m[3] + sy * rotPos.m[4] + sz * rotPos.m[5] + linPos.y;
    	dest->z = sx * rotPos.m[6] + sy * rotPos.m[7] + sz * rotPos.m[8] + linPos.z;
        dest++;
        source++;
    }
}

/***********************************************************************************
*   Vector3::invTransform
*
*   params :
*       Vector3 & r                     -   point in coord system 1
*       Matrix3x3 & rotPos              -   rotation from system 1 to 2
*       Vector3 & linPos                -   translation from system 1 to 2
*
*   Inverse of transform().
*
*   this = trans(rotPos) * (r - linPos)
***********************************************************************************/

void Vector3::invTransform(Vector3 & r, Matrix3x3 & rotPos, Vector3 & linPos)
{
    float sx, sy, sz;
    sx = r.x - linPos.x;
    sy = r.y - linPos.y;
    sz = r.z - linPos.z;
	x = sx * rotPos.m[0] + sy * rotPos.m[3] + sz * rotPos.m[6];
	y = sx * rotPos.m[1] + sy * rotPos.m[4] + sz * rotPos.m[7];
	z = sx * rotPos.m[2] + sy * rotPos.m[5] + sz * rotPos.m[8];
}

/***********************************************************************************
*   Vector3::invTransformArray
*
*   params :
*       Vector3 *source                 -   points in coord system 1
*       Vector3 *dest                   -   points in coord system 2
*       int N                           -   number of points
*       Matrix3x3 & rotPos              -   rotation from system 1 to 2
*       Vector3 & linPos                -   translation from system 1 to 2
*
*   Inverse of transformArray().
*
*   dest[i] = trans(rotPos) * (source[i] - linPos)
***********************************************************************************/

void Vector3::invTransformArray(Vector3 *source, Vector3 *dest, int N, Matrix3x3 & rotPos, Vector3 & linPos)
{
    float sx, sy, sz;
    for (int i = N; i--;) {
        sx = source->x - linPos.x;
        sy = source->y - linPos.y;
        sz = source->z - linPos.z;
    	dest->x = sx * rotPos.m[0] + sy * rotPos.m[3] + sz * rotPos.m[6];
	    dest->y = sx * rotPos.m[1] + sy * rotPos.m[4] + sz * rotPos.m[7];
    	dest->z = sx * rotPos.m[2] + sy * rotPos.m[5] + sz * rotPos.m[8];
        dest++;
        source++;
    }
}

/***********************************************************************************
*   Vector3::transformNormalArray
*
*   params :
*       Vector3 *source                 -   normals in coord system 1
*       Vector3 *dest                   -   normals in coord system 2
*       int N                           -   number of normals
*       Matrix3x3 & rotPos              -   rotation from system 1 to 2
*
*   Transform a list of normals to a new coords system.
*
*   dest[i] = rotPos * source[i]
***********************************************************************************/

void Vector3::transformNormalArray(Vector3 *source, Vector3 *dest, int N, Matrix3x3 & rotPos)
{
    float sx, sy, sz;
    for (int i = N; i--;) {
        sx = source->x;
        sy = source->y;
        sz = source->z;
    	dest->x = sx * rotPos.m[0] + sy * rotPos.m[1] + sz * rotPos.m[2];
	    dest->y = sx * rotPos.m[3] + sy * rotPos.m[4] + sz * rotPos.m[5];
    	dest->z = sx * rotPos.m[6] + sy * rotPos.m[7] + sz * rotPos.m[8];
        dest++;
        source++;
    }
}

/***********************************************************************************
*   Vector3::invTransformNormalArray
*
*   params :
*       Vector3 *source                 -   normals in coord system 1
*       Vector3 *dest                   -   normals in coord system 2
*       int N                           -   number of normals
*       Matrix3x3 & rotPos              -   rotation from system 1 to 2
*
*   Inverse of transformNormalArray()
*
*   dest[i] = trans(rotPos) * source[i]
***********************************************************************************/

void Vector3::invTransformNormalArray(Vector3 *source, Vector3 *dest, int N, Matrix3x3 & rotPos)
{
    float sx, sy, sz;
    for (int i = N; i--;) {
        sx = source->x;
        sy = source->y;
        sz = source->z;
    	dest->x = sx * rotPos.m[0] + sy * rotPos.m[3] + sz * rotPos.m[6];
	    dest->y = sx * rotPos.m[1] + sy * rotPos.m[4] + sz * rotPos.m[7];
    	dest->z = sx * rotPos.m[2] + sy * rotPos.m[5] + sz * rotPos.m[8];
        dest++;
        source++;
    }
}

