#include "IO.h"
using namespace std;
using namespace types;
bool io::readFromRawFile(std::string file_name, std::vector<types::AttDataRaw>& att_data_list)
{
    std::ifstream ifs(file_name);
    if(!ifs.is_open())
        return false;
    
    att_data_list.clear();

    int count = 1;

    for (int i = 0;i < count && !ifs.eof();i++) {
        AttDataRaw att_data;
        std::string line;
        std::getline(ifs, line, '{');
        if (i == 0) {
            auto n = line.substr(line.find("groupNumber") + 14);
            count = atof(n.substr(0, n.find(";")).c_str());
            if (!count)return false;
        }
        std::getline(ifs, line, '=');
        ifs >> att_data.header.timeCode;
        std::getline(ifs, line, '"');
        std::getline(ifs, att_data.header.dateTime, '"');
        std::getline(ifs, line, '=');
        ifs >> att_data.eulor1;
        std::getline(ifs, line, '=');
        ifs >> att_data.eulor2;
        std::getline(ifs, line, '=');
        ifs >> att_data.eulor3;
        std::getline(ifs, line, '=');
        ifs >> att_data.roll_velocity;
        std::getline(ifs, line, '=');
        ifs >> att_data.pitch_velocity;
        std::getline(ifs, line, '=');
        ifs >> att_data.yaw_velocity;
        std::getline(ifs, line, '=');
        ifs >> att_data.quaternion.x;
        std::getline(ifs, line, '=');
        ifs >> att_data.quaternion.y;
        std::getline(ifs, line, '=');
        ifs >> att_data.quaternion.z;
        std::getline(ifs, line, '=');
        ifs >> att_data.quaternion.w;

        att_data_list.push_back(att_data);
    }
    return true;
}

bool io::readFromRawFile(std::string file_name, std::vector<types::GPSDataRaw>& gps_data_list)
{
    std::ifstream ifs(file_name);
    if (!ifs.is_open())
        return false;

    gps_data_list.clear();

    int count = 1;

    for (int i = 0;i < count && !ifs.eof();i++) {
        GPSDataRaw gps_data;
        std::string line;
        std::getline(ifs, line, '{');
        if (i == 0) {
            auto n = line.substr(line.find("groupNumber") + 14);
            count = atof(n.substr(0, n.find(";")).c_str());
            if (!count)return false;
        }
        std::getline(ifs, line, '=');
        ifs >> gps_data.header.timeCode;
        std::getline(ifs, line, '"');
        std::getline(ifs, gps_data.header.dateTime, '"');
        std::getline(ifs, line, '=');
        ifs >> gps_data.PX;
        std::getline(ifs, line, '=');
        ifs >> gps_data.PY;
        std::getline(ifs, line, '=');
        ifs >> gps_data.PZ;
        std::getline(ifs, line, '=');
        ifs >> gps_data.VX;
        std::getline(ifs, line, '=');
        ifs >> gps_data.VY;
        std::getline(ifs, line, '=');
        ifs >> gps_data.VZ;

        gps_data_list.push_back(gps_data);
    }
    return true;
}

bool io::readFromRawFile(std::string file_name, std::vector<types::CCDAngleRaw>& ccd_angle_list) {
    std::ifstream ifs(file_name);
    if (!ifs.is_open())
        return false;

    ccd_angle_list.clear();

    int count = 0;
    ifs >> count;
    for (int i = 0;i < count && !ifs.eof();i++) {
        CCDAngleRaw ccd_angle;
        ifs >> ccd_angle.sample_id >> ccd_angle.PsiX >> ccd_angle.PsiY;
        ccd_angle_list.push_back(ccd_angle);
    }
    return true;
}

bool io::readFromRawFile(std::string file_name, std::vector<types::ImagingTimeRaw>& imaging_time_list) {
    std::ifstream ifs(file_name);
    if (!ifs.is_open())
        return false;

    imaging_time_list.clear();
    string desc;
    getline(ifs, desc, '\n');
    for (int i = 0;!ifs.eof();i++) {
        ImagingTimeRaw imaging;
        ifs >> imaging.line_id >> imaging.time >> imaging.deltaTime;
        imaging.header.timeCode = imaging.time;
        imaging_time_list.push_back(imaging);
    }
    if (!imaging_time_list.empty())imaging_time_list.pop_back();
    return true;
}

bool io::readFromRawFile(string file_name, std::vector<types::RotationRaw>& rotations) {
    std::ifstream ifs(file_name);
    if (!ifs.is_open())
        return false;
    rotations.clear();
    for (int i = 0;!ifs.eof();i++) {
        RotationRaw r;
        ifs >> r.header.timeCode
            >> r.R(0, 0) >> r.R(0, 1) >> r.R(0, 2)
            >> r.R(1, 0) >> r.R(1, 1) >> r.R(1, 2)
            >> r.R(2, 0) >> r.R(2, 1) >> r.R(2, 2);
        rotations.push_back(r);
    }
    if (!rotations.empty())rotations.pop_back();
}

bool io::readFromRawFile(std::string file_name, Eigen::MatrixXd& dem)
{
    std::ifstream ifs(file_name);
    if (!ifs.is_open())
        return false;
    double cols = dem.cols();
    double rows = dem.rows();
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            ifs >> dem(i, j);
        }
    }
}
