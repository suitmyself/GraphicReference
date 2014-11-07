/************************************************************************
* main.cpp
* piggy 2010-08-27
* to test the "RigidBody" class

************************************************************************/

#include "RigidBody.h"
#include <iostream>
using namespace std;

int main()
{
	float h = 0.01; // step size = 0.01 secs
	float g = 9.81; // gravity = 9.81 m/s/s
	float linKD = 0.0; // linear damping
	float rotKD = 0.0; // rotational damping
	RigidBody::setStatics(h, g); // set up step size and g for all
	                                        // rigid bodies
	RigidBody mainBody, tempBody; // declare two bodies. We will combine
	                                                  // the mass distributions into one
	                                                 // as an example.
	mainBody.setBox(4.0, 2.0, 1.0, 1000.0); // set up a box :
																// x radius = 4.0 metres
																// y radius = 2.0 metres
																// z radius = 1.0 metres
																// mass = 1000.0 kg
	tempBody.setCylinder(1.0, 2.0, 500.0); // set up a cylinder (along x axis) :
																// radius = 1.0 metres
																// height = 2.0 metres
																// mass = 500.0 kg
	Vector3 tempLinPos(5.0, 0.0, 0.0); // position and orientation of cylinder

	Matrix3x3 tempRotPos;              // relative to main body.(��һ�����Ҫ)
	tempRotPos.identity();
	//mainBody.combine(tempBody, tempRotPos, tempLinPos);
																								// merge the cylinder into the main body
																								// to create one body
	mainBody.diagonalise(); // diagonalise its moment of inertia and
										// move the centre of mass to origin

	Vector3 linPos(0.0, -10.0, 10.0); // initial position
	Vector3 linVel(0.0, 0.0, 0.0);  // initial velocity
	Quaternion qRotPos;             // initial rotation. Set up the
	//Vector3 rotAxis(1.0, 1.0, 0.0); // rotation as 0.3 radians around the
	//qRotPos.rotation(rotAxis, 0.3); // axis (1.0, 1.0, 0.0).
	//Vector3 rotVel(1.0, 2.0, -0.5); // initial angular velocity
	Vector3 rotVel(0, 0, 0); 
	mainBody.initialise(linPos, linVel, qRotPos, rotVel, linKD, rotKD);
																										// initialise state of body
	float time = 0.0;
	const int numSteps = 5; // take 5 steps between drawing graphics

	int step=1;
	while (step<=10) 
	{
		for (int i = numSteps; i--;) 
		{
			mainBody.RKInit(); // start Runge-Kutta
			mainBody.addWorldBodyForce(Vector3(1000,0,0),Vector3(0,0,0));//forces(); 
			mainBody.RKStep1();
			mainBody.addWorldBodyForce(Vector3(1000,0,0),Vector3(0,0,0));//forces();
			mainBody.RKStep2();
			mainBody.addWorldBodyForce(Vector3(1000,0,0),Vector3(0,0,0));//forces();
			mainBody.RKStep3();
			mainBody.addWorldBodyForce(Vector3(1000,0,0),Vector3(0,0,0));//forces();
			mainBody.RKStep4();
			time += RigidBody::RKh;
		}

		//drawBody(mainBody.linPos, mainBody.rotPos);
		// draw graphic at current position
		// with your own code.


		//�������
		cout<<"step:  "<<step<<endl;
		cout<<"��������:   "<<mainBody.mass<<endl;
		cout<<"��������������ϵ�е���������:   ("<<mainBody.linPos.x<<","<<mainBody.linPos.y<<","<<mainBody.linPos.z<<")"<<endl;
		//cout<<"��������������ϵ�е�ת������:   ("<<mainBody.rotPos.x<<","<<mainBody.rotPos.y<<","<<mainBody.rotPos.z<<")"<<endl;
		cout<<"��������������ϵ�е������ٶ�:   ("<<mainBody.linVel.x<<","<<mainBody.linVel.y<<","<<mainBody.linVel.z<<")"<<endl;
		cout<<"��������������ϵ�еĽ��ٶ�:  ("<<mainBody.rotVel.x<<","<<mainBody.rotVel.y<<","<<mainBody.rotVel.z<<")"<<endl;
		cout<<"��������������ϵ�е�����:  "<<mainBody.energy()<<endl;
		Vector3 linMomentum,rotMomentum;
		mainBody.momentum(linMomentum,rotMomentum);
		cout<<"����ƽ������:  ("<<linMomentum.x<<","<<linMomentum.y<<","<<linMomentum.z<<")"<<endl;
		cout<<"����ת���Ƕ���:  ("<<rotMomentum.x<<","<<rotMomentum.y<<","<<rotMomentum.z<<")"<<endl;

		cout<<"�����ת������:  ("<<mainBody.rotPos.n[0][0]<<","<<mainBody.rotPos.n[0][1]<<","<<mainBody.rotPos.n[0][2]<<") "
			               <<"("<<mainBody.rotPos.n[1][0]<<","<<mainBody.rotPos.n[1][1]<<","<<mainBody.rotPos.n[1][2]<<") " 
						   <<"("<<mainBody.rotPos.n[2][0]<<","<<mainBody.rotPos.n[2][1]<<","<<mainBody.rotPos.n[2][2]<<") "
						   <<endl;
		step++;
	}
	system("pause");

	return 0;
}