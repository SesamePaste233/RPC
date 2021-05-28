#include "AlgorithmBase.h"
using namespace tf;
std::vector<Eigen::Vector3d> utils::genRectGridWithLeveledHeight(int cols, int rows, int dims, double c_min, double c_max, double r_min, double r_max, double h_min, double h_max)
{
	std::vector<Eigen::Vector3d> pts;

	double d_c = (c_max - c_min) / double(cols);
	double d_r = (r_max - r_min) / double(rows);
	double d_h = (h_max - h_min) / double(dims);

	for (double h = h_min + d_h / 2;h < h_max;h += d_h) {
		for (double r = r_min + d_r / 2;r < r_max;r += d_r) {
			for (double c = c_min + d_c / 2;c < c_max;c += d_c) {
				Eigen::Vector3d pt;
				pt << c, r, h;
				pts.push_back(pt);
			}
		}
	}

	return pts;
}

types::Quaternion utils::interpolate(types::Quaternion p, types::Quaternion q, double t0, double t1, double t)
{
	//该函数用来进行四元数球面插值
	//t0时刻四元数为p；t1时刻四元数为q
	types::Quaternion r;
	//用点乘计算两个四元数夹角的cos值
	double k1 = 0;
	double k2 = 0;
	double cosAngle = p.w * q.w + p.x * q.x + p.y * q.y + p.z * q.z;
	double sinAngle;
	double Angle;
	//如果点乘为负，则翻转一个四元数以取得短的4D弧
	if (cosAngle < 0) {
		q.x = -q.x;
		q.y = -q.y;
		q.z = -q.z;
		cosAngle = -cosAngle;
	}

	if (cosAngle > 0.99999) {
		//重合
		k1 = (t1 - t) / (t1 - t0);
		k2 = (t - t0) / (t1 - t0);
	}

	else {
		sinAngle = sqrt(1 - cosAngle * cosAngle);
		Angle = atan2(sinAngle, cosAngle);
		k1 = sin((t1 - t) / (t1 - t0) * Angle) / sinAngle;
		k2 = sin((t - t0) / (t1 - t0) * Angle) / sinAngle;
	}
	r.w = p.w * k1 + q.w * k2;
	r.x = p.x * k1 + q.x * k2;
	r.y = p.y * k1 + q.y * k2;
	r.z = p.z * k1 + q.z * k2;
	return r;
}

std::vector<types::GPSDataRaw> utils::synchronize(std::vector<types::GPSDataRaw> GPS, std::vector<types::AttDataRaw> Att) {
	//GPS数据间隔1s，Att数据间隔0.25s
	std::vector<types::GPSDataRaw> newGPS;
	types::GPSDataRaw oneGPS;
	for (int i = 0; i < Att.size(); i++) {
		int j = int(i / 4);
		if (i % 4 == 0) {
			//当时间为整数时不用内插
			newGPS.push_back(GPS[j]);
		}
		else {

			oneGPS.header = Att[i].header;
			if (j > 0 && j < 99) {
				std::vector<double> v;
				std::vector<double> w;
				w.push_back(GPS[j - 1].header.timeCode);
				w.push_back(GPS[j].header.timeCode);
				w.push_back(GPS[j + 1].header.timeCode);
				w.push_back(GPS[j + 2].header.timeCode);
				v.push_back(GPS[j - 1].PX);
				v.push_back(GPS[j].PX);
				v.push_back(GPS[j + 1].PX);
				v.push_back(GPS[j + 2].PX);
				oneGPS.PX = lagrangeInterpolation(w, v, Att[i].header.timeCode);

				std::vector<double> v1;
				v1.push_back(GPS[j - 1].PY);
				v1.push_back(GPS[j].PY);
				v1.push_back(GPS[j + 1].PY);
				v1.push_back(GPS[j + 2].PY);
				oneGPS.PY = lagrangeInterpolation(w, v1, Att[i].header.timeCode);

				std::vector<double> v2;
				v2.push_back(GPS[j - 1].PZ);
				v2.push_back(GPS[j].PZ);
				v2.push_back(GPS[j + 1].PZ);
				v2.push_back(GPS[j + 2].PZ);
				oneGPS.PZ = lagrangeInterpolation(w, v2, Att[i].header.timeCode);

				std::vector<double> v3;
				v3.push_back(GPS[j - 1].VX);
				v3.push_back(GPS[j].VX);
				v3.push_back(GPS[j + 1].VX);
				v3.push_back(GPS[j + 2].VX);
				oneGPS.VX = lagrangeInterpolation(w, v3, Att[i].header.timeCode);

				std::vector<double> v4;
				v4.push_back(GPS[j - 1].VY);
				v4.push_back(GPS[j].VY);
				v4.push_back(GPS[j + 1].VY);
				v4.push_back(GPS[j + 2].VY);
				oneGPS.VY = lagrangeInterpolation(w, v4, Att[i].header.timeCode);

				std::vector<double> v5;
				v5.push_back(GPS[j - 1].VZ);
				v5.push_back(GPS[j].VZ);
				v5.push_back(GPS[j + 1].VZ);
				v5.push_back(GPS[j + 2].VZ);
				oneGPS.VZ = lagrangeInterpolation(w, v5, Att[i].header.timeCode);

			}
			newGPS.push_back(oneGPS);
		}

	}
	return newGPS;
}

double utils::lagrangeInterpolation(std::vector<double> v, std::vector<double> w, double x)
{
	double output = 0;
	double L;
	int k = v.size();
	for (int i = 0; i < k; i++) {
		L = 1;
		for (int j = 0; j < k; j++) {

			if (j != i) {
				L *= (x - v[j]) / (v[i] - v[j]);
			}

		}
		output += L * w[i];
	}
	return output;
}

types::PRY utils::interpolate(Eigen::Matrix3d R1, Eigen::Matrix3d R2, double t1, double t2, double T)
{
	types::PRY newPRY;
	types::PRY PRY1 = toPRY(R1);
	types::PRY PRY2 = toPRY(R2);
	double p = PRY1.p + (PRY2.p - PRY1.p) * (T - t1) / (t2 - T);
	double r = PRY1.r + (PRY2.r - PRY1.r) * (T - t1) / (t2 - T);
	double y = PRY1.y + (PRY2.y - PRY1.y) * (T - t1) / (t2 - T);
	newPRY.p = p;
	newPRY.r = r;
	newPRY.y = y;
	return newPRY;
}

Eigen::Vector3d utils::interpolate(std::vector<types::GPSDataRaw> gps_round, double time)
{
	Eigen::Vector3d XYZ;
	std::vector<double> time_seq;
	std::vector<double> dataX_seq;
	std::vector<double> dataY_seq;
	std::vector<double> dataZ_seq;
	if (time< gps_round[0].header.timeCode || time > gps_round[gps_round.size() - 1].header.timeCode) {
		std::cout << "time is not include in your time_seq!" << std::endl;
	}
	for (int i = 0; i < gps_round.size(); i++) {
		time_seq.push_back(gps_round[i].header.timeCode);
		dataX_seq.push_back(gps_round[i].PX);
		dataY_seq.push_back(gps_round[i].PY);
		dataZ_seq.push_back(gps_round[i].PZ);
	}
	XYZ(0, 0) = lagrangeInterpolation(time_seq, dataX_seq, time);
	XYZ(1, 0) = lagrangeInterpolation(time_seq, dataY_seq, time);
	XYZ(2, 0) = lagrangeInterpolation(time_seq, dataZ_seq, time);
	return XYZ;
}

