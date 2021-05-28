#pragma once
#include "CommonHeader.h"
#include <random>
#include <time.h>
#include "Eigen/Sparse"
#include "Eigen/SparseCholesky"
#include "Eigen/SparseQR"


namespace opt {
		
	enum class QuitFlag
	{
		MinErrorReached = 0,
		QuitCondSatisfied = 1,
		MaxIterReached = 2,
		UnknownError = -1,
		UnfeasableError = -2
	};

	enum class Method {
		GaussNewton=1,
		LevenbergMarquardt=2,
		SparseLM=3,
		SparseQR=4
	};

	class LeastSquareSolver {//Least Square Optimization Solver. In form of V = Jacobi * X - L
	protected:
		std::function<bool(void*)> obj_fcn;

		std::vector<std::function<bool()>> quit_conditions;

		//Read-only
		Eigen::MatrixXd& Jacobi;
		Eigen::VectorXd& L;
		Eigen::VectorXd& X;

		//Status indicators
		bool is_solved;
		bool is_feasable;

		double max_error;
		int max_iteration;
		unsigned int observation_count;
		unsigned int observation_number;

		float _sigma;

	public:
		Eigen::MatrixXd Hessian;
		int Hessian_dims;

		Eigen::MatrixXd Gradiant;

		Method method;

		int function_num;

		int iteration;

		bool is_linear;

		bool dropX;

		Eigen::VectorXd sigma_x;

		Eigen::VectorXd deltaX;

		LeastSquareSolver(Eigen::MatrixXd& A, Eigen::VectorXd& L, Eigen::VectorXd& X, int observations = 1) :Jacobi(A), L(L), X(X), max_error(3E-5), max_iteration(20), iteration(0), observation_number(observations), _sigma(0), dropX(false){
			deltaX.resizeLike(X);
			Hessian_dims = X.rows();
			Gradiant.resizeLike(X);
			Hessian.fill(0);
			Gradiant.fill(0);
			function_num = 0;
			method = Method::GaussNewton;
			is_solved = false;
			is_feasable = true;
			is_linear = false;
		}

		inline bool& isLinear() {
			return is_linear;
		}

		void setObjectFunction(std::function<bool(void*)> obj_fcn) {
			this->obj_fcn = obj_fcn;
			is_feasable = true;
		}

		void addQuitCondition(std::function<bool()> quit) {
			this->quit_conditions.push_back(quit);
			is_feasable = true;
		}

		void clearQuitCondition() {
			quit_conditions.clear();
			is_feasable = true;
		}

		double& maxError() {
			is_feasable = true;
			return max_error;
		}

		int& maxIteration() {
			is_feasable = true;
			return max_iteration;
		}

		void setObservationNum(unsigned int i) {
			this->observation_number = i;
		}

		int i() {
			return observation_count;
		}

		bool is_enough_observation() {
			if (observation_count >= observation_number - 1) {
				return true;
			}
			else {
				return false;
			}
		}

		QuitFlag solve(void* data = nullptr); /**/

		float sigma();
	};

	class RANSAC {

	protected:
		std::function<bool(void*)> obj_fcn;

		std::function<float(int)> evaluate_observation;

		Eigen::MatrixXd& Jacobi;
		Eigen::VectorXd& L;//error
		Eigen::VectorXd& X;//model

		Eigen::VectorXd X_init;

		std::list<int> possible_inliners;

		int max_iteration;

		float inliner_tolerance;

		float init_estimate_observation_ratio;

		int total_observations;

		opt::LeastSquareSolver lss_solver;

		std::list<int>::iterator inliner_count;

		Eigen::VectorXd best_model;

		float best_error;

		std::list<int> best_inliner;

		bool is_linear;

		Method method;

	public:

		RANSAC(Eigen::MatrixXd& A, Eigen::VectorXd& L, Eigen::VectorXd& X, int observation_count) :
			Jacobi(A), L(L), X(X),
			lss_solver(A, L, X),
			inliner_tolerance(0),
			max_iteration(10),
			X_init(X),
			total_observations(observation_count),
			is_linear(false),
			method(opt::Method::LevenbergMarquardt)
		{
			lss_solver.isLinear() = is_linear;
			lss_solver.method = method;
		};

		inline bool& isLinear() {
			return is_linear;
		}

		inline int i() {
			int i = *inliner_count;
			inliner_count++;
			return i;
		}

		inline void setInitEstiRatio(float i) {
			if (i < 1 && i>0) {
				this->init_estimate_observation_ratio = i;
			}
		}

		bool is_enough_observation() {
			if (inliner_count == possible_inliners.end()) {
				inliner_count = possible_inliners.begin();
				return true;
			}
			return false;
		}

		void setObjectFunction(std::function<bool(void*)> obj_fcn) {
			this->obj_fcn = obj_fcn;
			this->lss_solver.setObjectFunction(obj_fcn);
		}

		void setMaxIteration(int max_iter) {
			this->max_iteration = max_iter;
		}

		void setEvaluationFunction(std::function<float(int)> eval) {
			this->evaluate_observation = eval;
		}

		QuitFlag solve(void* data = nullptr);

		inline Eigen::VectorXd bestModel() {
			return best_model;
		};

		inline float error() {
			return best_error;
		}

		std::list<int> bestInliner() {
			return best_inliner;
		}
	};
}