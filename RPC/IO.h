#pragma once
#include "CommonHeader.h"
#include "BaseType.h"
namespace io {
	bool readFromRawFile(std::string file_name, std::vector<types::AttDataRaw>& att_data_list);
	bool readFromRawFile(std::string file_name, std::vector<types::GPSDataRaw>& gps_data_list);
	bool readFromRawFile(std::string file_name, std::vector<types::CCDAngleRaw>& ccd_angle_list);
	bool readFromRawFile(std::string file_name, std::vector<types::ImagingTimeRaw>& imaging_time_list);
	bool readFromRawFile(std::string file_name, std::vector<types::RotationRaw>& rotations);
}