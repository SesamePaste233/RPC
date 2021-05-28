#include "Model.h"

Eigen::Vector3d RPM::predict(double u, double v, double h)
{
    auto v1 = floor(v),v2 = ceil(v);
    double PsiX1 = this->ccd_angles[v1].PsiX;
    double PsiY1 = this->ccd_angles[v1].PsiY;
    double PsiX2 = this->ccd_angles[v2].PsiX;
    double PsiY2 = this->ccd_angles[v2].PsiY;
    double PsiX = (PsiX2 - PsiX1) * (v - v1) + PsiX1;
    double PsiY = (PsiY2 - PsiY1) * (v - v1) + PsiY1;

    Eigen::Vector3d t;
    Eigen::Matrix3d R;
    this->getExtrinsicElems(u, t, R);

    Eigen::Vector3d position;
    position << tan(PsiY), tan(PsiX), -1;

    position = R * position;

    //Çó½âÍÖÇò½Çµã
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

void RPM::getExtrinsicElems(double line_id, Eigen::Vector3d& translation, Eigen::Matrix3d& rotation)
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
    
    auto floor_ceil_gps = utils::findFloorAndCeil(time, this->gps_data, 2);
    translation = utils::interpolate(floor_ceil_gps, time);
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

    A.fill(0);
    coeff.fill(0);

    opt::LeastSquareSolver solver(A, L, coeff, n);

    auto objectFunction = [&](void*)->bool {
        int i = solver.i();
        A.fill(0);
        double x = X(i), y = Y(i), z = Z(i), r = Line(i), c = Sample(i);
        Eigen::VectorXd diff(19);
        diff << z, y, x, z* y, z* x, y* x, z* z, y* y, x* x, z* y* x, z* z* y, z* z* x, y* y* z, y* y* x, x* x* z, x* x* y, z* z* z, y* y* y, x* x* x;

        double B = diff.dot(coeff.segment(20,19)) + 1;
        double D = diff.dot(coeff.segment(59,19)) + 1;

        A(0, 0) = 1;
        A.row(0).segment(1, 19) = diff / B, A.row(0).segment(20, 19) =  -r * diff / B;
        A(1, 39) = 1;
        A.row(1).segment(40, 19) = diff / D, A.row(1).segment(59, 19) = -r * diff / D;

        //std::cout << A << std::endl;

        L(0) = r / B;

        L(1) = c / D;

        return solver.is_enough_observation();
    };

    solver.setObjectFunction(objectFunction);

    solver.dropX = true;

    solver.maxError() = 1e-6;

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
