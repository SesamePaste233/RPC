#pragma once
#include"Optimizer.h"
#include"CommonHeader.h"
#include"BaseType.h"

#define PI 3.1415926
namespace types {

	typedef struct {
		double X, Y, Z;
	}XYZ;

	typedef struct {
		double B, L, H;
	}BLH;

	typedef struct {
		double p, r, y;
	}PRY;//pitch roll yaw
}

namespace tf {
	Eigen::Matrix3d toRotationMatrix(types::PRY PRY);

	types::PRY toPRY(Eigen::Matrix3d R);

	Eigen::Matrix3d toRotationMatrix(types::Quaternion Q);

	types::BLH XYZ2BLH(types::XYZ XYZ);

	types::XYZ BLH2XYZ(types::BLH BLH);

}


