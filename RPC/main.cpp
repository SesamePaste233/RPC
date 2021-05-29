#include "Optimizer.h"
#include "IO.h"
#include "AlgorithmBase.h"
#include "Model.h"

int main(){
	std::vector<types::AttDataRaw> att_datas;
	std::vector<types::GPSDataRaw> gps_datas;
	std::vector<types::CCDAngleRaw> ccd_datas;
	std::vector<types::ImagingTimeRaw> imaging_times;
	std::vector<types::RotationRaw> J2000_WGS84;
	io::readFromRawFile("DX_ZY3_NAD_att.txt", att_datas);
	io::readFromRawFile("DX_ZY3_NAD_gps.txt", gps_datas);
	io::readFromRawFile("NAD.cbr", ccd_datas);
	io::readFromRawFile("DX_ZY3_NAD_imagingTime.txt", imaging_times);
	io::readFromRawFile("jw.txt", J2000_WGS84);

	att_datas = utils::splitDataByTime(131862404.25, 131862408.25, att_datas);
	gps_datas = utils::splitDataByTime(131862402, 131862413, gps_datas);

	RPM rpm;
	rpm.feed(att_datas, ccd_datas, gps_datas, imaging_times, J2000_WGS84);
	rpm.pry_sen2body = { -0.000511776876952,0.001828916699906,0.003770429577750 };
	rpm.WGS84_A = 6378137, rpm.WGS84_B = 6356752.3;

	/*
	auto img_pts = utils::genRectGridWithLeveledHeight(10, 10, 5, 0, 5378, 0, 8192, 45, 65);

	std::vector<Eigen::Vector2d> ctrl_pts_img;
	std::vector<Eigen::Vector3d> ctrl_pts_obj;
	std::vector<types::BLH> BLHarray;

	for (auto pt : img_pts) {
		auto XYZ = rpm.predict(pt);
		types::XYZ xyz = { XYZ[0],XYZ[1],XYZ[2] };
		auto BLH = tf::XYZ2BLH(xyz);
		ctrl_pts_img.push_back(pt.head(2));
		ctrl_pts_obj.push_back(Eigen::Vector3d(BLH.B, BLH.L, BLH.H));
		BLHarray.push_back(BLH);
	}
	*/

	/*
	RFM rfm;
	rfm.solve(ctrl_pts_obj, ctrl_pts_img);

	Eigen::MatrixXd ctrl_img_reprojection;
	ctrl_img_reprojection = rfm.forward(BLHarray);
	double reproj_error_x = 0, reproj_error_y = 0;
	for (int i = 0; i < ctrl_pts_obj.size(); i++) {
		Eigen::Vector2d ctrl_img_reprojection_pt = ctrl_img_reprojection.row(i);
		std::cout << "Id: " << i << "\t Error (x,y): ";
		std::cout << std::setprecision(8) << std::fixed << ctrl_img_reprojection(i, 0) - ctrl_pts_img[i](0) << '\t';
		std::cout << std::setprecision(8) << std::fixed << ctrl_img_reprojection(i, 1) - ctrl_pts_img[i](1) << std::endl;
		reproj_error_x += pow(ctrl_img_reprojection(i, 0) - ctrl_pts_img[i](0), 2);
		reproj_error_y += pow(ctrl_img_reprojection(i, 1) - ctrl_pts_img[i](1), 2);
	}
	reproj_error_x = sqrt(reproj_error_x / double(ctrl_pts_obj.size()));
	reproj_error_y = sqrt(reproj_error_y / double(ctrl_pts_obj.size()));

	std::cout << "------------------------------------------------------------" << std::endl;
	std::cout << "像方重投影误差(像素): " << std::setprecision(8) << std::fixed << reproj_error_x <<'\t'<< reproj_error_y << std::endl;

	std::ofstream ofs("zy3_coeffs.txt");
	ofs << rfm;
	*/
	return 0;
}