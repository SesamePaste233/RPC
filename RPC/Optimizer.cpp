#include "Optimizer.h"
#define MAX(a,b) ((a)>(b)?(a):(b))


Eigen::VectorXd sqrt(Eigen::VectorXd value) {
	int n = value.size();
	Eigen::VectorXd new_value(n);
	for (int i = 0;i < n;i++) {
		new_value(i) = sqrt(value(i));
	}
	return new_value;
}


opt::QuitFlag opt::LeastSquareSolver::solve(void* data)
{
	this->is_solved = false;
	if (is_linear) {
		this->X.fill(0);
	}
	if (is_solved && !is_feasable) {
		return QuitFlag::UnfeasableError;
	}
	try {
		bool meet_condition = false;
		if (method == Method::GaussNewton) {
			Hessian.resize(Hessian_dims, Hessian_dims);
			for (iteration = 0; iteration < max_iteration && !is_solved && !meet_condition;iteration++) {
				Hessian.fill(0);
				Gradiant.fill(0);

				bool flag = true;
				do {//Calculate hessian matrix
					flag = !obj_fcn(data);
					Hessian += Jacobi.transpose() * Jacobi;
					Gradiant += Jacobi.transpose() * L;

					_sigma += (Jacobi * X - L).transpose() * (Jacobi * X - L);
					observation_count++;
				} while (flag);
				observation_number = observation_count;
				observation_count = 0;

				if (is_linear && iteration != 0) {
					is_solved = true;
					continue;
				}

				deltaX = pinv(Hessian) * Gradiant;//Update equation for Gauss-Newton Method

				//std::cout << deltaX <<std::endl << std::endl;

				//std::cout << Hessian << std::endl << Gradiant << std::endl;

				Eigen::VectorXd tempX(X.size());
				tempX.fill(0);
				if (dropX) {
					tempX = X;
					X.fill(0);
				}
				X += deltaX;

				meet_condition = false;
				for (auto cond : quit_conditions) {
					meet_condition = cond();
				}
				if (meet_condition) {
					return QuitFlag::QuitCondSatisfied;
				}

				for (int j = 0;j < X.rows();j++) {
					if (abs(deltaX(j)-tempX(j)) > max_error) {
						is_solved = false;
						_sigma = 0;
						break;
					}
					else {
						is_solved = true;
					}
				}

				if (is_solved) {
					_sigma = 0;
					do {//Calculate hessian matrix
						flag = !obj_fcn(data);
						Hessian += Jacobi.transpose() * Jacobi;
						_sigma += (Jacobi * deltaX - L).transpose() * (Jacobi * deltaX - L);
						observation_count++;
					} while (flag);
					observation_number = observation_count;
					observation_count = 0;

					Eigen::MatrixXd Q = Hessian.inverse();
					float sigma0 = this->sigma();
					sigma_x = sqrt(Q.diagonal()) * sigma0;
				}
			}
		}
		else if (method == Method::LevenbergMarquardt) {
			Hessian.resize(Hessian_dims, Hessian_dims);
			double mu = 1, tau = 1E-3, rho = 1, f_old = 0, v = 2;
			auto g = Gradiant;
			auto h = Hessian;
			for (iteration = 0; iteration < max_iteration && !is_solved && !meet_condition;iteration++) {
				Hessian.fill(0);
				Gradiant.fill(0);

				double f_new = 0;
				bool flag = true;

				do {//Calculate hessian matrix
					flag = !obj_fcn(data);
					Hessian += Jacobi.transpose() * Jacobi;
					Gradiant += Jacobi.transpose() * L;
					f_new += L.array().square().sum();

					_sigma += (Jacobi * X - L).transpose() * (Jacobi * X - L);
					observation_count++;
				} while (flag);
				observation_number = observation_count;
				observation_count = 0;

				if (is_linear && iteration != 0) {
					is_solved = true;
					continue;
				}

				if (iteration == 0) {
					mu = tau * Hessian.diagonal().maxCoeff();
					f_old = f_new;
				}
				else {
					rho = (f_old - f_new) / (deltaX.transpose() * g + 0.5 * deltaX.transpose() * h * deltaX).value();
					if (rho > 0) {
						f_old = f_new;
						mu = mu * MAX(1.0 / 3.0, pow(1 - (2 * rho - 1), 3)), v = 2;
					}
					else {
						X -= deltaX;
						Gradiant = g;
						Hessian = h;
						mu = mu * v, v = v * 2;
					}
				}

				g = Gradiant, h = Hessian;

				deltaX = (Hessian + mu * Eigen::MatrixXd::Identity(Hessian.rows(), Hessian.cols())).inverse() * Gradiant;//Update equation for LM Method
				
				X += deltaX;


				meet_condition = false;
				for (auto cond : quit_conditions) {
					meet_condition = cond();
				}
				if (meet_condition) {
					return QuitFlag::QuitCondSatisfied;
				}

				for (int j = 0;j < X.rows();j++) {
					if (abs(deltaX(j)) > max_error) {
						is_solved = false;
						_sigma = 0;
						break;
					}
					else {
						is_solved = true;
					}
				}
				if (is_solved) {
					_sigma = 0;
					do {//Calculate hessian matrix
						flag = !obj_fcn(data);
						_sigma += (Jacobi * deltaX - L).transpose() * (Jacobi * deltaX - L);
						observation_count++;
					} while (flag);
					observation_number = observation_count;
					observation_count = 0;

					Eigen::MatrixXd Q = Hessian.inverse();
					float sigma0 = this->sigma();
					sigma_x = sqrt(Q.diagonal()) * sigma0;
				}
			}
		}
		else if (method == Method::SparseLM) {

			Eigen::SparseMatrix<double> sparse_Hessian(Hessian_dims,Hessian_dims);

			double mu = 1, tau = 1E-3, rho = 1, f_old = 0, v = 2;
			auto g = Gradiant;
			auto h = sparse_Hessian;
			for (iteration = 0; iteration < max_iteration && !is_solved && !meet_condition;iteration++) {
				Gradiant.fill(0);
				sparse_Hessian.setZero();

				double f_new = 0;
				bool flag = true;
				do {//Calculate hessian matrix
					flag = !obj_fcn(data);

					Eigen::SparseMatrix<double,Eigen::RowMajor> sparse_Jacobi(Jacobi.rows(),Jacobi.cols());
					for (int r = 0;r < Jacobi.rows();r++) {
						for (int c = 0;c < Jacobi.cols();c++) {
							double j = Jacobi(r, c);
							if (j!=0)
								sparse_Jacobi.insert(r, c) = j;
						}
					}

					sparse_Hessian += sparse_Jacobi.transpose() * sparse_Jacobi;

					Gradiant += Jacobi.transpose() * L;

					f_new += L.array().square().sum();

					_sigma += (Jacobi * X - L).transpose() * (Jacobi * X - L);

					observation_count++;
				} while (flag);
				observation_number = observation_count;
				observation_count = 0;

				if (is_linear && iteration != 0) {
					is_solved = true;
					continue;
				}

				if (iteration == 0) {
					mu = tau * sparse_Hessian.diagonal().maxCoeff();
					f_old = f_new;
				}else {
					rho = (f_old - f_new) / (deltaX.transpose() * g + 0.5 * deltaX.transpose() * h * deltaX).value();
					if (rho > 0) {
						f_old = f_new;
						mu = mu * MAX(1.0 / 3.0, pow(1 - (2 * rho - 1), 3)), v = 2;
					}
					else {
						X -= deltaX;
						Gradiant = g;
						sparse_Hessian = h;
						mu = mu * v, v = v * 2;
					}
				}

				g = Gradiant, h = sparse_Hessian;

				Eigen::SparseMatrix<double> eye(Hessian_dims, Hessian_dims);
				for (int i = 0;i < Hessian_dims;i++) {
					eye.insert(i, i) = mu;
				}

				Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;

				solver.compute(sparse_Hessian + eye);
				if (solver.info() != Eigen::Success) {
					std::cout << "ERROR: SparseLM failed to compute Hessian matrix at iteration: " << iteration << "." << std::endl;
					throw;
				}

				deltaX = solver.solve(Gradiant);//Update equation for LM Method

				if (solver.info() != Eigen::Success) {
					std::cout << "ERROR: SparseLM failed to compute delta-X at iteration: " << iteration << "." << std::endl;
					throw;
				}

				X += deltaX;


				meet_condition = false;
				for (auto cond : quit_conditions) {
					meet_condition = cond();
				}
				if (meet_condition) {
					return QuitFlag::QuitCondSatisfied;
				}

				for (int j = 0;j < X.rows();j++) {
					if (abs(deltaX(j)) > max_error) {
						is_solved = false;
						_sigma = 0;
						break;
					}
					else {
						is_solved = true;
					}
				}
			}
		}
		else if (method == Method::SparseQR) {
			if (function_num == 0) {
				std::cout << "ERROR: at SparseQR. Please set function number before calculation." << std::endl;
				throw;
			}
			Eigen::SparseMatrix<double,Eigen::RowMajor> sparse_Jacobi(function_num,Hessian_dims);
			Eigen::VectorXd B(function_num);
			for (iteration = 0; iteration < max_iteration && !is_solved && !meet_condition;iteration++) {
				Gradiant.fill(0);
				sparse_Jacobi.setZero();
				bool flag = true;
				int fcn_n = 0;
				do {//Calculate sparse Jacobian matrix
					flag = !obj_fcn(data);
					for (int r = 0;r < Jacobi.rows();r++) {
						for (int c = 0;c < Jacobi.cols();c++) {
							double j = Jacobi(r, c);
							if(j!=0)
								sparse_Jacobi.insert(fcn_n + r, c) = j;
						}
						B(fcn_n + r) = L(r);
					}
					fcn_n += Jacobi.rows();

					_sigma += (Jacobi * X - L).transpose() * (Jacobi * X - L);
					observation_count++;
				} while (flag);
				observation_number = observation_count;
				observation_count = 0;

				if (is_linear && iteration != 0) {
					is_solved = true;
					continue;
				}

				Eigen::SparseQR<Eigen::SparseMatrix<double>,Eigen::COLAMDOrdering<int>> solver;

				solver.compute(sparse_Jacobi);

				deltaX = solver.solve(B);

				X += deltaX;

				meet_condition = false;
				for (auto cond : quit_conditions) {
					meet_condition = cond();
				}
				if (meet_condition) {
					return QuitFlag::QuitCondSatisfied;
				}

				for (int j = 0;j < X.rows();j++) {
					if (abs(deltaX(j)) > max_error) {
						is_solved = false;
						_sigma = 0;
						break;
					}
					else {
						is_solved = true;
					}
				}
			}
		}
	}
	catch (...) {
		return QuitFlag::UnknownError;
	}
	if (is_solved) {
		return QuitFlag::MinErrorReached;
	}
	else {
		return QuitFlag::MaxIterReached;
	}
}

float opt::LeastSquareSolver::sigma()
{
	float redundancy = observation_number * L.size() - X.size();

	return sqrt(_sigma/redundancy);
}

std::list<int> randomNIntNumber(int n, int min, int max,bool allow_same) {
	std::uniform_int_distribution<int> u(min,max);
	std::default_random_engine generator;

	generator.seed(clock());

	std::list<int> list;
	while (list.size() < n) {
		list.push_back(u(generator));
		if (!allow_same) {
			list.sort();
			list.unique();
		}
	}

	return list;
}

opt::QuitFlag opt::RANSAC::solve(void* data)
{
	lss_solver.isLinear() = is_linear;
	if (is_linear)method = Method::GaussNewton;
	lss_solver.method = method;
	int iterations = 0;
	float max_error = 10000000;

	int k = init_estimate_observation_ratio * total_observations;

	int min_observations = X.size() / L.size() + 1;

	if (inliner_tolerance == 0) {
		min_observations++;
	}

	if (k < min_observations)k = min_observations;

	while (iterations < this->max_iteration) {
		this->possible_inliners = randomNIntNumber(k, 0, total_observations - 1, false);
		inliner_count = possible_inliners.begin();

		auto new_inliner = possible_inliners;

		lss_solver.solve();

		float estimated_error = 0;
		float tolerance = 0;
		if (inliner_tolerance == 0) {
			for (auto i:possible_inliners) {
				estimated_error += evaluate_observation(i);
			}
			estimated_error /= (float)possible_inliners.size();
			tolerance = 6 * estimated_error;
		}
		else {
			tolerance = inliner_tolerance;
		}

		for (int i = 0;i < total_observations;i++) {
			if (find(possible_inliners.begin(), possible_inliners.end(), i) == possible_inliners.end()) {
				float error = evaluate_observation(i);
				if (error < tolerance) {
					new_inliner.push_back(i);
				}
			}
		}
		

		possible_inliners = new_inliner;

		lss_solver.solve();
		float error = 0;
		for (auto i : possible_inliners) {
			error += evaluate_observation(i);
		}
		error /= (float)possible_inliners.size();
		
		if (error < max_error) {
			max_error = error;
			best_model = X;
			best_error = error;
			best_inliner = new_inliner;
		}
		if (possible_inliners.size() == total_observations) {
			break;
		}
		iterations++;
	}

	return opt::QuitFlag::MaxIterReached;
}
