#pragma once
#include "BaseType.h"
#include "CommonHeader.h"
#include "Transform.h"
#include "AlgorithmBase.h"

typedef class RPM
{
public:
	std::vector<types::AttDataRaw> attitudes;
	std::vector<types::CCDAngleRaw> ccd_angles;
	std::vector<types::GPSDataRaw> gps_data;
	std::vector<types::ImagingTimeRaw> imaging_times;
	std::vector<types::RotationRaw> J2000_WGS84_rotations;

	types::PRY pry_sen2body;

	double WGS84_A, WGS84_B;

	RPM() = default;

	inline bool feed(std::vector<types::AttDataRaw> attitudes,
		std::vector<types::CCDAngleRaw> ccd_angles,
		std::vector<types::GPSDataRaw> gps_data,
		std::vector<types::ImagingTimeRaw> imaging_times,
		std::vector<types::RotationRaw> J2000_WGS84_rotations
		) 
	{
		this->attitudes = attitudes;
		this->ccd_angles = ccd_angles;
		this->gps_data = gps_data;
		this->imaging_times = imaging_times;
		this->J2000_WGS84_rotations = J2000_WGS84_rotations;
		return true;
	};

	Eigen::Vector3d predict(double u, double v, double h);

	Eigen::Vector3d predict(Eigen::Vector3d pt_with_h);

	Eigen::Vector2d forward(Eigen::Vector3d obj_pt);

protected:

	Eigen::Vector2d reproject(double line_id,Eigen::Vector3d obj_pt);

	bool getExtrinsicElems(double line_id, Eigen::Vector3d& translation, Eigen::Matrix3d& rotation);

	inline double getImagingTime(double line_id) {
		double line_below = std::floor(line_id), line_above = std::ceil(line_id);
		double t1 = imaging_times[line_below-1].time;
		double t2 = imaging_times[line_above-1].time;
		return t1 + (t2 - t1) * (line_id - line_below);
	};

}RigorousPhysicalModel;

typedef class RFM {
public:
	Eigen::Vector2d X_RECT_COEFF;
	Eigen::Vector2d Y_RECT_COEFF;
	Eigen::Vector2d Z_RECT_COEFF;
	Eigen::Vector2d LINE_RECT_COEFF;
	Eigen::Vector2d SAMPLE_RECT_COEFF;
	Eigen::VectorXd LINE_NUM_COEFF;
	Eigen::VectorXd SAMP_NUM_COEFF;
	Eigen::VectorXd LINE_DEN_COEFF;
	Eigen::VectorXd SAMP_DEN_COEFF;

	RFM() {
		this->SAMP_NUM_COEFF = Eigen::VectorXd(20);
		this->LINE_NUM_COEFF = Eigen::VectorXd(20);
		this->SAMP_DEN_COEFF = Eigen::VectorXd(20);
		this->LINE_DEN_COEFF = Eigen::VectorXd(20);
		this->SAMP_DEN_COEFF.fill(1), this->LINE_DEN_COEFF.fill(1);
	};

	void solve(std::vector<Eigen::Vector3d> BLH_pts,std::vector<Eigen::Vector2d> img_pts);

	Eigen::MatrixXd forward(std::vector<types::BLH> BLH);

	bool readCoeffFromFile(std::string file_name);

	inline friend std::ostream& operator<<(std::ostream& os,const RFM a) {
		os << "X_RECT_COEFF = " << std::endl << "{" << std::endl << a.X_RECT_COEFF << std::endl << "}" << std::endl;
		os << "Y_RECT_COEFF = " << std::endl << "{" << std::endl << a.Y_RECT_COEFF << std::endl << "}" << std::endl;
		os << "Z_RECT_COEFF = " << std::endl << "{" << std::endl << a.Z_RECT_COEFF << std::endl << "}" << std::endl;
		os << "LINE_RECT_COEFF = " << std::endl << "{" << std::endl << a.LINE_RECT_COEFF << std::endl << "}" << std::endl;
		os << "SAMPLE_RECT_COEFF = " << std::endl << "{" << std::endl << a.SAMPLE_RECT_COEFF << std::endl << "}" << std::endl;
		os << "LINE_NUM_COEFF = " << std::endl << "{" << std::endl << a.LINE_NUM_COEFF << std::endl << "}" << std::endl;
		os << "LINE_DEN_COEFF = " << std::endl << "{" << std::endl << a.LINE_DEN_COEFF << std::endl << "}" << std::endl;
		os << "SAMP_NUM_COEFF = " << std::endl << "{" << std::endl << a.SAMP_NUM_COEFF << std::endl << "}" << std::endl;
		os << "SAMP_DEN_COEFF = " << std::endl << "{" << std::endl << a.SAMP_DEN_COEFF << std::endl << "}" << std::endl;
		return os;
	}

protected:
	static Eigen::VectorXd normalize(Eigen::VectorXd value, double& offset, double& scale);

}RationalFunctionModel;