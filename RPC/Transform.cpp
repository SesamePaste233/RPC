#include "Transform.h"

Eigen::Matrix3d tf::toRotationMatrix(types::PRY PRY)
{
	Eigen::Matrix3d R;
	Eigen::Matrix3d Rp;
	Eigen::Matrix3d Rr;
	Eigen::Matrix3d Ry;
	double p = PRY.p;
	double r = PRY.r;
	double y = PRY.y;
	Rp(0, 0) = cos(p);
	Rp(0, 1) = 0;
	Rp(0, 2) = sin(p);
	Rp(1, 0) = 0;
	Rp(1, 1) = 1;
	Rp(1, 2) = 0;
	Rp(2, 0) = -sin(p);
	Rp(2, 1) = 0;
	Rp(2, 2) = cos(p);

	Rr(0, 0) = 1;
	Rr(0, 1) = 0;
	Rr(0, 2) = 0;
	Rr(1, 0) = 0;
	Rr(1, 1) = cos(r);
	Rr(1, 2) = -sin(r);
	Rr(2, 0) = 0;
	Rr(2, 1) = sin(r);
	Rr(2, 2) = cos(r);

	Ry(0, 0) = cos(y);
	Ry(0, 1) = -sin(y);
	Ry(0, 2) = 0;
	Ry(1, 0) = sin(y);
	Ry(1, 1) = cos(y);
	Ry(1, 2) = 0;
	Ry(2, 0) = 0;
	Ry(2, 1) = 0;
	Ry(2, 2) = 1;
	return Rp*Rr*Ry;
}

types::PRY tf::toPRY(Eigen::Matrix3d R)
{
	types::PRY PRY;
	PRY.p = atan2(R(0, 2), R(2, 2));
	PRY.r = atan2(-R(1, 2), sqrt(R(1, 0)* R(1, 0) + R(1, 1)* R(1, 1)));
	PRY.y = atan2(R(1, 0), R(1, 1));
	return PRY;
}

Eigen::Matrix3d tf::toRotationMatrix(types::Quaternion Q) {
	Eigen::Matrix3d R;
	double x = Q.x;
	double y = Q.y;
	double z = Q.z;
	double w = Q.w;
	R(0, 0) = 1 - 2 * y * y - 2 * z * z;
	R(0, 1) = 2 * x * y - 2 * z * w;
	R(0, 2) = 2 * x * z + 2 * y * w;
	R(1, 0) = 2 * x * y + 2 * z * w;
	R(1, 1) = 1 - 2 * x * x - 2 * z * z;
	R(1, 2) = 2 * y * z - 2 * x * w;
	R(2, 0) = 2 * x * z - 2 * y * w;
	R(2, 1) = 2 * y * z + 2 * x * w;
	R(2, 2) = 1 - 2 * x * x - 2 * y * y;
	return R;
}

types::BLH tf::XYZ2BLH(types::XYZ XYZ)
{
	types::BLH BLH;
	double a = 6378137;//地球长半轴
	double b = 6356755;//地球短半轴
	double e2 = 0.00669437999013;//(a * a - b * b) / (b * b);//第一偏心率的平方
	double X = XYZ.X;
	double Y = XYZ.Y;
	double Z = XYZ.Z;

	double L = std::atan2(Y, X);
	//由于B、H计算时两者函数相关，故先取B0，迭代计算
	double newB = atan(Z / sqrt(X * X + Y * Y));
	double B = 0;
	double H = 0;
	do {
		B = newB;
		double N = a / sqrt(1 - e2 * sin(B) * sin(B));//卯酉圈半径
		H = Z / sin(B) - N * (1 - e2);
		newB = atan((Z * N + Z * H) / (sqrt(X * X + Y * Y) * (N - N * e2 + H)));
	} while (fabs(B - newB) > 1e-10);

	BLH.B = newB * 180 / PI;
	BLH.H = H;
	BLH.L = L * 180 / PI;

	return BLH;
}


types::XYZ tf::BLH2XYZ(types::BLH BLH)
{
	types::XYZ XYZ;
	double a = 6378137;
	double e2 = 0.00669437999013;
	double B = BLH.B * PI / 180;
	double L = BLH.H * PI / 180;
	double W = sqrt(1 - e2 * sin(B) * sin(B));
	double N = a / W;

	XYZ.X= (N + BLH.H) * cos(B) * cos(L);
	XYZ.Y= (N + BLH.H) * cos(B) * sin(L);
	XYZ.Z = (N * (1 - e2) + BLH.H) * sin(B);
	return XYZ;
}
