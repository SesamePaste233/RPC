#pragma once
#include "CommonHeader.h"

namespace types {
    typedef struct {
        double timeCode;
        std::string dateTime;
    }Header;

    typedef struct {
        double x, y, z, w;
    }Quaternion;

	typedef struct {
        Header header;
        double eulor1,eulor2,eulor3;
        double roll_velocity,pitch_velocity,yaw_velocity;
        Quaternion quaternion;
	}AttDataRaw;

    typedef struct {
        Header header;
        double PX, PY, PZ;
        double VX, VY, VZ;
    }GPSDataRaw;

    typedef struct {
        long sample_id;
        double PsiX;
        double PsiY;
    }CCDAngleRaw;

    typedef struct {
        Header header;
        long line_id;
        double time, deltaTime;
    }ImagingTimeRaw;

    typedef struct {
        Header header;
        Eigen::Matrix3d R;
    }RotationRaw;


}