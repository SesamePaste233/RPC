#include "Model.h"

#define _X(a,b) ((a).cwiseProduct(b))

Eigen::Vector3d RPM::predict(double u, double v, double h)
{
    double PsiX = this->ccd_angles[v - 1].PsiX;
    double PsiY = this->ccd_angles[v - 1].PsiY;

    Eigen::Vector3d t;
    Eigen::Matrix3d R;
    this->getExtrinsicElems(u, t, R);

    Eigen::Vector3d position;
    position << tan(PsiY), tan(PsiX), -1;

    //视方向
    position = R * position;

    //求解椭球交点
    double A = this->WGS84_A + h, B = this->WGS84_B + h;

    double a = (position(0) * position(0) + position(1) * position(1)) / (A * A) + position(2) * position(2) / (B * B);
    double b = 2 * ((t(0) * position(0) + t(1) * position(1)) / (A * A) + t(2) * position(2) / (B * B));
    double c = (t(0) * t(0) + t(1) * t(1)) / (A * A) + t(2) * t(2) / (B * B) - 1;
    double t1 = (-b - sqrt(b * b - 4 * a * c)) / (2 * a), t2 = (-b + sqrt(b * b - 4 * a * c)) / (2 * a);
    double m = std::max(t1, t2);

    return t + m * position;
}

Eigen::Vector3d RPM::predict(Eigen::Vector3d pt_with_h)
{
    return this->predict(pt_with_h(0), pt_with_h(1), pt_with_h(2));
}


Eigen::Vector2d RPM::forward(Eigen::Vector3d obj_pt)
{
    int low = 0, high = this->imaging_times[imaging_times.size() - 1].line_id;
    std::pair<double, double> ls(low, this->reproject(low, obj_pt)(0));
    std::pair<double, double> le(high, this->reproject(high, obj_pt)(0));
    double new_l = 0;
    int iter = 0;
    bool found = false;
    for (double l = (low + high) / 2.f;iter<20;iter++,l = new_l) {
        std::pair<double, double> lm(l, this->reproject(l, obj_pt)(0));
        if (lm.second * ls.second < 0) {
            le = lm;
        }
        else if (lm.second * le.second < 0) {
            ls = lm;
        }
        else return Eigen::Vector2d(0, 0);
        new_l = (ls.first + le.first) / 2.f;
        if (abs(l - new_l) < 10) {
            found = true;
            break;
        }
    }

    if (found) {
        double min = 10000;
        int line_id = (ls.first+le.first)/2;
        for (int l = ls.first;l < le.first;l++) {
            double dx = abs(this->reproject(l, obj_pt)(0));
            if (dx < min) {
                min = dx;
                line_id = l;
            }
        }

        int sample_id = 0;
        double dy = atan(this->reproject(line_id-1, obj_pt)(1));
        auto near = utils::_find_floor_and_ceil<types::CCDAngleRaw>(dy, this->ccd_angles, [](types::CCDAngleRaw a){return a.PsiX;});
        if (near.size() == 1) {
            sample_id = near[0].sample_id;
        }
        else if(near.size() == 2) {
            if (abs(near[0].PsiX - dy) < abs(near[1].PsiX - dy)) {
                sample_id = near[0].sample_id;
            }
            else {
                sample_id = near[1].sample_id;
            }
        }
        else {
            return Eigen::Vector2d(0, 0);
        }

        return Eigen::Vector2d(line_id + 1, sample_id + 1);
    }

    return Eigen::Vector2d(0, 0);
}

Eigen::Vector2d RPM::reproject(double line_id, Eigen::Vector3d obj_pt)
{
    Eigen::Vector3d t;Eigen::Matrix3d R;
    getExtrinsicElems(line_id + 1, t, R);
    Eigen::Vector3d Psi = R.transpose() * (obj_pt - t);
    double PsiY = -Psi(0) / Psi(2), PsiX = -Psi(1) / Psi(2);
    return Eigen::Vector2d(PsiY, PsiX);
}

bool RPM::getExtrinsicElems(double line_id, Eigen::Vector3d& translation, Eigen::Matrix3d& rotation)
{
    double time = this->getImagingTime(line_id);
    auto R_s2b = tf::toRotationMatrix(this->pry_sen2body);

    auto floor_ceil_q = utils::findFloorAndCeil(time, this->attitudes);
    auto quaternion = floor_ceil_q[0].quaternion;
    if (floor_ceil_q.size() != 1) {
        quaternion = utils::interpolate(floor_ceil_q[0].quaternion, floor_ceil_q[1].quaternion, floor_ceil_q[0].header.timeCode, floor_ceil_q[1].header.timeCode, time);
    }

    auto floor_ceil_R = utils::findFloorAndCeil(time, this->J2000_WGS84_rotations);
    types::PRY q_j2w;
    if (floor_ceil_R.size() != 1) {
        q_j2w = utils::interpolate(floor_ceil_R[0].R, floor_ceil_R[1].R, floor_ceil_R[0].header.timeCode, floor_ceil_R[1].header.timeCode, time);
    }
    else {
        q_j2w = tf::toPRY(floor_ceil_R[0].R);
    }

    auto R_b2j = tf::toRotationMatrix(quaternion);
    auto R_j2w = tf::toRotationMatrix(q_j2w);

    rotation = R_j2w * R_b2j * R_s2b;
    
    translation = utils::interpolate(gps_data, time);

    return true;
}

void RFM::solve(std::vector<Eigen::Vector3d> BLH_pts, std::vector<Eigen::Vector2d> img_pts)
{
    int n = BLH_pts.size();
    Eigen::VectorXd X(n), Y(n), Z(n), Line(n), Sample(n);
    for (int i = 0;i < n;i++) {
        X(i) = BLH_pts[i](0);
        Y(i) = BLH_pts[i](1);
        Z(i) = BLH_pts[i](2);
        Line(i) = img_pts[i](0);
        Sample(i) = img_pts[i](1);
    }

    double offset, scale;
    X = this->normalize(X, offset, scale);
    this->X_RECT_COEFF << offset, scale;
    Y = this->normalize(Y, offset, scale);
    this->Y_RECT_COEFF << offset, scale;
    Z = this->normalize(Z, offset, scale);
    this->Z_RECT_COEFF << offset, scale;
    Line = this->normalize(Line, offset, scale);
    this->LINE_RECT_COEFF << offset, scale;
    Sample = this->normalize(Sample, offset, scale);
    this->SAMPLE_RECT_COEFF << offset, scale;

    Eigen::VectorXd coeff(78);
    Eigen::MatrixXd A(2, 78);
    Eigen::VectorXd L(2);

    /*std::cout << X << std::endl << std::endl;
    std::cout << Y << std::endl << std::endl;
    std::cout << Z << std::endl << std::endl;
    std::cout << Line << std::endl << std::endl;
    std::cout << Sample << std::endl << std::endl;*/

    A.fill(0);
    coeff.fill(0);

    Eigen::MatrixXd diff(n, 19);
    diff << Z, Y, X, _X(Z, Y), _X(Z, X), _X(Y, X), _X(Z, Z), _X(Y, Y), _X(X, X), _X(_X(Z, Y), X), _X(_X(Z, Z), Y), _X(_X(Z, Z), X), _X(_X(Y, Y), Z), _X(_X(Y, Y), X), _X(_X(Z, X), X), _X(_X(X, Y), X), _X(_X(Z, Z), Z), _X(_X(Y, Y), Y), _X(_X(X, X), X);

    //std::cout<< std::setprecision(10) << diff;
    //Eigen::VectorXd Line_coeff(39);
    //Eigen::VectorXd Samp_coeff(39);
    //Eigen::MatrixXd M(n, 39);
    //Eigen::MatrixXd N(n, 39);
    //for (int i = 0;i < n;i++) {
    //    M.row(i) << 1, diff.row(i), -Line(i) * diff.row(i);
    //    N.row(i) << 1, diff.row(i), -Sample(i) * diff.row(i);
    //}
    //Eigen::VectorXd B(500, 1), D(500, 1);
    //Samp_coeff.fill(0), Line_coeff.fill(0);
    //B.fill(0), D.fill(0);
    //

    //for (int iteration = 0;iteration < 20;iteration++) {
    //    B = diff * Line_coeff.tail(19) + Eigen::VectorXd::Ones(500), D = diff * Samp_coeff.tail(19) + Eigen::VectorXd::Ones(500);
    //    Eigen::MatrixXd Wr = Eigen::VectorXd::Ones(500).cwiseQuotient(B).asDiagonal(), Wc = Eigen::VectorXd::Ones(500).cwiseQuotient(D).asDiagonal();
    //    //std::cout << Wr << std::endl;
    //    Eigen::MatrixXd Nbb = M.transpose() * Wr * Wr * M;
    //    Eigen::VectorXd Wb = M.transpose() * Wr * Wr * Line;
    //    Eigen::MatrixXd Ndd = N.transpose() * Wc * Wc * N;
    //    Eigen::VectorXd Wd = N.transpose() * Wc * Wc * Sample;

    //    Line_coeff = opt::pinv(Nbb) * Wb;
    //    Samp_coeff = opt::pinv(Ndd) * Wd;
    //    std::cout << Line_coeff << std::endl;
    //}

    opt::LeastSquareSolver solver(A, L, coeff, n);

    auto objectFunction = [&](void*)->bool {
        int i = solver.i();
        A.fill(0);
        double r = Line(i), c = Sample(i);
        double B = diff.row(i).dot(coeff.segment(20,19)) + 1;
        double D = diff.row(i).dot(coeff.segment(59,19)) + 1;

        A(0, 0) = 1 / B;
        A.row(0).segment(1, 19) = diff.row(i) / B, A.row(0).segment(20, 19) =  -r * diff.row(i) / B;
        A(1, 39) = 1 / D;
        A.row(1).segment(40, 19) = diff.row(i) / D, A.row(1).segment(59, 19) = -c * diff.row(i) / D;

        //std::cout << A << std::endl;

        L(0) = r / B;

        L(1) = c / D;

        return solver.is_enough_observation();
    };

    auto condition = [&]()->bool {
        std::cout << solver.sigma() << std::endl;
        return solver.sigma() < 1e-4;
    };
    solver.addQuitCondition(condition);
    solver.setObjectFunction(objectFunction);

    solver.dropX = true;

    solver.maxError() = 1e-6;

    solver.maxIteration() = 30;

    auto flag = solver.solve();

    if (flag == opt::QuitFlag::MaxIterReached) {
        std::cout << "Optimization failed." << std::endl;
    }
    else if (flag == opt::QuitFlag::MinErrorReached) {
        std::cout << "Solved with minimun error reached in iteration " << solver.iteration << " ." << std::endl;
    }
    else if (flag == opt::QuitFlag::QuitCondSatisfied) {
        std::cout << "Solved with condition satisfied in iteration " << solver.iteration << " ." << std::endl;
    }

    this->LINE_NUM_COEFF = coeff.head(20);
    this->LINE_DEN_COEFF.segment(1, 19) = coeff.segment(20, 19);
    this->SAMP_NUM_COEFF = coeff.segment(39, 20);
    this->SAMP_DEN_COEFF.segment(1, 19) = coeff.segment(59, 19);
}

Eigen::VectorXd RFM::normalize(Eigen::VectorXd value, double& offset, double& scale)
{
    offset = value.mean();
    scale = std::max(abs(value.maxCoeff() - offset), abs(value.minCoeff() - offset));
    Eigen::VectorXd ones(value.size());
    ones.fill(1);
    return (value - ones * offset) / scale;
}

Eigen::MatrixXd RFM::forward(std::vector<types::BLH> BLH)
{
    int n = BLH.size();
    Eigen::MatrixXd LS(n, 2);
    Eigen::VectorXd X(n);
    Eigen::VectorXd Y(n);
    Eigen::VectorXd Z(n);
    Eigen::VectorXd ones(n);
    ones.fill(1);
    for (int i = 0; i < n; i++) {
        X(i) = BLH[i].B;
        Y(i) = BLH[i].L;
        Z(i) = BLH[i].H;
    }
    //正则化BLH
    Eigen::VectorXd normX;
    double X_offset, X_scale;
    Eigen::VectorXd normY;
    double Y_offset, Y_scale;
    Eigen::VectorXd normZ;
    double Z_offset, Z_scale;
    normX = (X - ones * this->X_RECT_COEFF(0)) / X_RECT_COEFF(1);
    normY = (Y - ones * this->Y_RECT_COEFF(0)) / Y_RECT_COEFF(1);
    normZ = (Z - ones * this->Z_RECT_COEFF(0)) / Z_RECT_COEFF(1);
    Eigen::MatrixXd A(n, 20);
    for (int i = 0; i < n; i++) {
        A(i, 0) = ones(i);
        A(i, 1) = normZ(i);
        A(i, 2) = normY(i);
        A(i, 3) = normX(i);
        A(i, 4) = normZ(i) * normY(i);
        A(i, 5) = normZ(i) * normX(i);
        A(i, 6) = normY(i) * normX(i);
        A(i, 7) = normZ(i) * normZ(i);
        A(i, 8) = normY(i) * normY(i);
        A(i, 9) = normX(i) * normX(i);
        A(i, 10) = normZ(i) * normY(i) * normX(i);
        A(i, 11) = normZ(i) * normZ(i) * normY(i);
        A(i, 12) = normZ(i) * normZ(i) * normX(i);
        A(i, 13) = normY(i) * normY(i) * normZ(i);
        A(i, 14) = normY(i) * normY(i) * normX(i);
        A(i, 15) = normX(i) * normX(i) * normZ(i);
        A(i, 16) = normX(i) * normX(i) * normY(i);
        A(i, 17) = normZ(i) * normZ(i) * normZ(i);
        A(i, 18) = normY(i) * normY(i) * normY(i);
        A(i, 19) = normX(i) * normX(i) * normX(i);
    }

    Eigen::VectorXd line(n);
    Eigen::VectorXd sample(n);
    line = (A * LINE_NUM_COEFF).cwiseQuotient(A * LINE_DEN_COEFF);
    sample = (A * SAMP_NUM_COEFF).cwiseQuotient(A * SAMP_DEN_COEFF);

    //反正则化
    for (int i = 0; i < n; i++) {
        LS(i, 0) = line(i) * LINE_RECT_COEFF(1) + LINE_RECT_COEFF(0);
        LS(i, 1) = sample(i) * SAMPLE_RECT_COEFF(1) + SAMPLE_RECT_COEFF(0);
        //        std::cout << LS(i, 0) << "   " << LS(i, 1) << std::endl;
    }
    return LS;
}

bool RFM::readCoeffFromFile(std::string file_name)
{
    std::ifstream ifs(file_name);
    if (!ifs.is_open())
        return false;

    int count = 9;
    bool read[9] = { false };
    int i = 0;
    for (i = 0;i < count && !ifs.eof();) {
        std::string line;
        std::getline(ifs, line, '{');
        double value = 0;
        if (!read[0] && line.find("X_RECT_COEFF") != line.npos) {
            for (int j = 0;j < 2;j++) {
                ifs >> value;
                this->X_RECT_COEFF(j) = value;
            }
            i++, read[0] = true;
        }
        else if (!read[1] && line.find("Y_RECT_COEFF") != line.npos) {
            for (int j = 0;j < 2;j++) {
                ifs >> value;
                this->Y_RECT_COEFF(j) = value;
            }
            i++, read[1] = true;
        }
        else if (!read[2] && line.find("Z_RECT_COEFF") != line.npos) {
            for (int j = 0;j < 2;j++) {
                ifs >> value;
                this->Z_RECT_COEFF(j) = value;
            }
            i++, read[2] = true;
        }
        else if (!read[3] && line.find("LINE_RECT_COEFF") != line.npos) {
            for (int j = 0;j < 2;j++) {
                ifs >> value;
                this->LINE_RECT_COEFF(j) = value;
            }
            i++, read[3] = true;
        }
        else if (!read[4] && line.find("SAMPLE_RECT_COEFF") != line.npos) {
            for (int j = 0;j < 2;j++) {
                ifs >> value;
                this->SAMPLE_RECT_COEFF(j) = value;
            }
            i++, read[4] = true;
        }
        else if (!read[5] && line.find("LINE_NUM_COEFF") != line.npos) {
            for (int j = 0;j < 20;j++) {
                ifs >> value;
                this->LINE_NUM_COEFF(j) = value;
            }
            i++, read[5] = true;
        }
        else if (!read[6] && line.find("LINE_DEN_COEFF") != line.npos) {
            for (int j = 0;j < 20;j++) {
                ifs >> value;
                this->LINE_DEN_COEFF(j) = value;
            }
            i++, read[6] = true;
        }
        else if (!read[7] && line.find("SAMP_NUM_COEFF") != line.npos) {
            for (int j = 0;j < 20;j++) {
                ifs >> value;
                this->SAMP_NUM_COEFF(j) = value;
            }
            i++, read[7] = true;
        }
        else if (!read[8] && line.find("SAMP_DEN_COEFF") != line.npos) {
            for (int j = 0;j < 20;j++) {
                ifs >> value;
                this->SAMP_DEN_COEFF(j) = value;
            }
            i++, read[8] = true;
        }
    }
    if (i == 9) {
        return true;
    }

    return false;
}
