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

	
	auto img_pts_ctrl = utils::genRectGridWithLeveledHeight(10, 10, 5, 0, 5378, 0, 8192, 45, 65);

	std::vector<Eigen::Vector2d> ctrl_pts_img;
	std::vector<Eigen::Vector3d> ctrl_pts_obj;
	std::vector<types::BLH> BLHarray;

	for (auto pt : img_pts_ctrl) {
		auto XYZ = rpm.predict(pt);
		types::XYZ xyz = { XYZ[0],XYZ[1],XYZ[2] };
		auto BLH = tf::XYZ2BLH(xyz);
		ctrl_pts_img.push_back(pt.head(2));
		ctrl_pts_obj.push_back(Eigen::Vector3d(BLH.B, BLH.L, BLH.H));
		BLHarray.push_back(BLH);
	}
	

	
	RFM rfm;
	//rfm.readCoeffFromFile("zy3_coeffs.txt");
	
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

	Eigen::MatrixXd DEM(206, 327);

	io::readFromRawFile("dem.txt", DEM);

	double North = 35.96528996;
	double West = 114.60535575;
	double dLat = 0.0007999999999;
	double dLon = 0.0007999999999;
	auto ground_pts_BLH = utils::genRectGridbyDEM(DEM, North, West, dLat, dLon, 20, 20);
	std::vector<Eigen::Vector3d> obj_pts;
	std::vector<Eigen::Vector2d> img_pts;
	std::vector<types::BLH> obj_pts_BLH;


	for (auto pt : ground_pts_BLH) {
		auto XYZ = tf::BLH2XYZ({ pt(0),pt(1),pt(2) });
		auto img_pt = rpm.forward(Eigen::Vector3d(XYZ.X, XYZ.Y, XYZ.Z));
		if (img_pt(0) == 0 || img_pt(1) == 0) {
			continue;
		}
		obj_pts_BLH.push_back({ pt(0),pt(1),pt(2) });
		img_pts.push_back(img_pt);
	}

	Eigen::MatrixXd check_uv = rfm.forward(obj_pts_BLH);
	
	double pix_error_x = 0;
	double pix_error_y = 0;
	for (int i = 0;i < img_pts.size();i++) {
		pix_error_x += (check_uv(i, 0) - img_pts[i](0)) * (check_uv(i, 0) - img_pts[i](0));
		pix_error_y += (check_uv(i, 1) - img_pts[i](1)) * (check_uv(i, 1) - img_pts[i](1));
		std::cout << "检查点像方误差(像素): " << std::setprecision(8) << std::fixed << check_uv(i, 0) - img_pts[i](0) << '\t' << check_uv(i, 1) - img_pts[i](1) << std::endl;
	}
	pix_error_x = sqrt(pix_error_x / img_pts.size());
	pix_error_y = sqrt(pix_error_y / img_pts.size());
	double pix_error_xy = sqrt(pix_error_x * pix_error_x + pix_error_y * pix_error_y);
	std::cout << "------------------------------------------------------------" << std::endl;
	std::cout << "检查点像方误差(像素): " << std::setprecision(8) << std::fixed << pix_error_x << '\t' << pix_error_y << '\t' << pix_error_xy << std::endl;
	/*
	std::vector<Eigen::Vector2d> check_pts_img = img_pts;
	std::vector<Eigen::Vector3d> check_pts_obj;
	std::vector<types::BLH> BLHarray;
	
	for (auto pt : img_pts) {
		auto XYZ = rpm.predict(pt);
		types::XYZ xyz = { XYZ(0),XYZ(1),XYZ(2) };
		auto BLH = tf::XYZ2BLH(xyz);
		check_pts_obj.push_back(Eigen::Vector3d(BLH.B, BLH.L, BLH.H));
		BLHarray.push_back(BLH);
	}
	
	double check_error_B = 0;
	double check_error_L = 0;
	double check_error_H = 0;
	for (int i = 0;i < obj_pts_BLH.size();i++) {
		check_error_B += (obj_pts_BLH[i].B - BLHarray[i].B) * (obj_pts_BLH[i].B - BLHarray[i].B);
		check_error_L += (obj_pts_BLH[i].L - BLHarray[i].L) * (obj_pts_BLH[i].L - BLHarray[i].L);
		check_error_H += (obj_pts_BLH[i].H - BLHarray[i].H) * (obj_pts_BLH[i].H - BLHarray[i].H);
	}
	check_error_B = sqrt(check_error_B / obj_pts_BLH.size());
	check_error_L = sqrt(check_error_L / obj_pts_BLH.size());
	check_error_H = sqrt(check_error_H / obj_pts_BLH.size());
	std::cout << "------------------------------------------------------------" << std::endl;
	std::cout << "检查点物方误差: " << std::setprecision(8) << std::fixed << check_error_B << '\t' << check_error_L << '\t' << check_error_H <<std::endl;
	return 0;
	*/
}