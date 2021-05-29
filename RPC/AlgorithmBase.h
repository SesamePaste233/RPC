#pragma once
#include "CommonHeader.h"
#include "BaseType.h"
#include "Transform.h"

namespace utils {
	/// <summary>
	/// �����������н�ȡ�ض�ʱ�䷶Χ�����ݵ��㷨
	/// </summary>
	/// <typeparam name="T">��������</typeparam>
	/// <param name="t1">��Сʱ��</param>
	/// <param name="t2">���ʱ��</param>
	/// <param name="data_seq">��������</param>
	/// <returns>ʱ������Ĵ�����������</returns>
	template<class T>
	std::vector<T> splitDataByTime(double t1, double t2, std::vector<T> data_seq) {
		if (t1 >= t2)return data_seq;
		std::sort(data_seq.begin(), data_seq.end(), [](T a, T b) {return a.header.timeCode < b.header.timeCode;});
		auto iter1 = std::find_if(data_seq.begin(), data_seq.end(), [&](T a) {return a.header.timeCode >= t1;});
		auto iter2 = std::find_if(data_seq.begin(), data_seq.end(), [&](T a) {return a.header.timeCode > t2;});
		std::vector<T> output;
		output.assign(iter1, iter2-1);
		return output;
	}

	/// <summary>
	/// ����ʱ������������������ض�ʱ�丽�������ݵĿ��������㷨
	/// </summary>
	/// <typeparam name="T">��������</typeparam>
	/// <param name="t">Ѱ�ҵ�ʱ��</param>
	/// <param name="data_seq">��������</param>
	/// <param name="range_off_set">��Ҫ�ķ�Χ</param>
	/// <returns></returns>
	template<class T>
	auto findFloorAndCeil(double t, std::vector<T> data_seq,int range_off_set = 0) {
		//std::sort(data_seq.begin(), data_seq.end(), [](T a, T b) {return a.header.timeCode < b.header.timeCode;});

		int size = data_seq.size();
		double delta_t = (data_seq[size - 1].header.timeCode - data_seq[0].header.timeCode) / double(size - 1);

		double st = t - data_seq[0].header.timeCode;
		auto begin = data_seq.begin();
		auto end = data_seq.end();

		int radius = 3;
		if (delta_t != 0) {
			if (st / delta_t - radius < 0)return std::vector<T>();
			begin += st / delta_t - radius;
			if (st / delta_t + radius >= data_seq.size()) {
				end = data_seq.begin() + data_seq.size() - 1;
			}
			else end = begin + 2 * radius;
		}
		auto iter = data_seq.begin();
		do {
			iter = std::find_if(begin, end, [&](T a) {return a.header.timeCode > t;});
			if (iter != data_seq.end())break;
			try {
				begin -= radius;
				end += radius;
			}
			catch (...) {
				return std::vector<T>();
			}
		} while (iter == data_seq.end());

		std::vector<T> output;
		if ((iter - 1)->header.timeCode == t) {
			for (auto i = iter - 1 - range_off_set;i != iter + range_off_set;i++) {
				output.push_back(*i);
			}
		}
		for (auto i = iter - 1 - range_off_set;i != iter + 1 + range_off_set;i++) {
			output.push_back(*i);
		}
		return output;
	}


	template<class T>
	auto _find_floor_and_ceil(double v, std::vector<T> data_seq, int range_off_set = 0) {
		//std::sort(data_seq.begin(), data_seq.end(), [](T a, T b) {return a.header.timeCode < b.header.timeCode;});

		int size = data_seq.size();
		double delta_t = (data_seq[size - 1].header.timeCode - data_seq[0].header.timeCode) / double(size - 1);

		double st = t - data_seq[0].header.timeCode;
		auto begin = data_seq.begin();
		auto end = data_seq.end();

		int radius = 3;
		if (delta_t != 0) {
			if (st / delta_t - radius < 0)return std::vector<T>();
			begin += st / delta_t - radius;
			if (st / delta_t + radius >= data_seq.size()) {
				end = data_seq.begin() + data_seq.size() - 1;
			}
			else end = begin + 2 * radius;
		}
		auto iter = data_seq.begin();
		do {
			iter = std::find_if(begin, end, [&](T a) {return a.header.timeCode > t;});
			if (iter != data_seq.end())break;
			try {
				begin -= radius;
				end += radius;
			}
			catch (...) {
				return std::vector<T>();
			}
		} while (iter == data_seq.end());

		std::vector<T> output;
		if ((iter - 1)->header.timeCode == t) {
			for (auto i = iter - 1 - range_off_set;i != iter + range_off_set;i++) {
				output.push_back(*i);
			}
		}
		for (auto i = iter - 1 - range_off_set;i != iter + 1 + range_off_set;i++) {
			output.push_back(*i);
		}
		return output;
	}

	std::vector<Eigen::Vector3d> genRectGridWithLeveledHeight(int cols, int rows, int dims, double c_min, double c_max, double r_min, double r_max, double h_min, double h_max);

	types::Quaternion interpolate(types::Quaternion p, types::Quaternion q, double t0, double t1, double t);

	//t1Ϊ��һ��ʱ�̣�t2Ϊ��һ��ʱ�̣�T�Ǳ�ʱ��
	types::PRY interpolate(Eigen::Matrix3d R1, Eigen::Matrix3d R2, double t1, double t2, double T);

	Eigen::Vector3d interpolate(std::vector<types::GPSDataRaw> gps_round, double time);

	//��ȡ����attʱ�̵�λ�ú��ٶ�
	std::vector<types::GPSDataRaw> synchronize(std::vector<types::GPSDataRaw> GPS, std::vector<types::AttDataRaw> Att);

	//@brief ����ʱ������time_seq ������ٶ�����w��kΪ���г��ȣ�xΪ��ȡʱ��
	//@param time_seq ʱ������
	double lagrangeInterpolation(std::vector<double> time_seq, std::vector<double> w, double x);

}