
#include "Attitude.h"
#include "libmesh/libmesh.h"

using namespace Eigen;
using namespace libMesh;

Attitude::Attitude(Real attack, Real sideslip, Real pitch, Real yaw, Real roll) :
	_attack(attack),
	_sideslip(sideslip),
	_pitch(pitch),
	_yaw(yaw),
	_roll(roll)
{

}

Eigen::Quaterniond Attitude::bodyFromWind()
{
	Quaterniond q_attack(AngleAxisd(_attack, Vector3d::UnitY()));
	Quaterniond q_sideslip (AngleAxisd(pi-_sideslip,  Vector3d::UnitZ()));
	return q_attack*q_sideslip;
}

Eigen::Quaterniond Attitude::earthFromBody()
{
	Quaterniond q_pitch(AngleAxisd(_pitch, Vector3d::UnitY()));
	Quaterniond q_yaw(AngleAxisd(_yaw, Vector3d::UnitZ()));
	Quaterniond q_roll(AngleAxisd(_roll, Vector3d::UnitX()));
	return q_roll*q_pitch*q_yaw;
}

Eigen::Quaterniond Attitude::earthFromWind()
{
	return earthFromBody()*bodyFromWind();
}
