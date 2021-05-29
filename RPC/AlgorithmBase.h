#pragma once
#include "CommonHeader.h"
#include "BaseType.h"
#include "Transform.h"

namespace utils {
	/// <summary>
	/// 对于数据序列截取特定时间范围的数据的算法
	/// </summary>
	/// <typeparam name="T">数据类型</typeparam>
	/// <param name="t1">最小时间</param>
	/// <param name="t2">最大时间</param>
	/// <param name="data_seq">数据序列</param>
	/// <returns>时间升序的待求数据序列</returns>
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
	/// 对于时间升序的数据序列求特定时间附近的数据的快速搜索算法
	/// </summary>
	/// <typeparam name="T">数据类型</typeparam>
	/// <param name="t">寻找的时间</param>
	/// <param name="data_seq">数据序列</param>
	/// <param name="range_off_set">需要的范围</param>
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

	//t1为上一个时刻，t2为下一个时刻，T是本时刻
	types::PRY interpolate(Eigen::Matrix3d R1, Eigen::Matrix3d R2, double t1, double t2, double T);

	Eigen::Vector3d interpolate(std::vector<types::GPSDataRaw> gps_round, double time);

	//获取所有att时刻的位置和速度
	std::vector<types::GPSDataRaw> synchronize(std::vector<types::GPSDataRaw> GPS, std::vector<types::AttDataRaw> Att);

	//@brief 输入时间序列time_seq 坐标和速度序列w，k为序列长度，x为所取时刻
	//@param time_seq 时间序列
	double lagrangeInterpolation(std::vector<double> time_seq, std::vector<double> w, double x);

}