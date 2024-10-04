#include "compute_g.h"
#include <igl/cotmatrix.h>
#include <igl/grad.h>
#include <Eigen/Sparse>
#include <limits>
#include <cmath>

Eigen::VectorXd compute_g(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    const Eigen::VectorXd& f,
    const Eigen::MatrixXd& GF,
    const Eigen::MatrixXd& B,
    const Eigen::SparseMatrix<double>& G)
{

    std::cout << "V dimensions: " << V.rows() << "x" << V.cols() << std::endl;
    std::cout << "F dimensions: " << F.rows() << "x" << F.cols() << std::endl;
    std::cout << "f dimensions: " << f.rows() << "x" << f.cols() << std::endl;
    std::cout << "GF dimensions: " << GF.rows() << "x" << GF.cols() << std::endl;
    std::cout << "B dimensions: " << B.rows() << "x" << B.cols() << std::endl;
    std::cout << "G dimensions: " << G.rows() << "x" << G.cols() << std::endl;

    // Initialize g with random values
    Eigen::VectorXd g = Eigen::VectorXd::Random(V.rows()) * 0.01;

    // Set g(B) = 0
    for (int i = 0; i < B.rows(); ++i) {
        // Find the index of this point in V
        for (int j = 0; j < V.rows(); ++j) {
            if (B.row(i) == V.row(j)) {
                g(j) = 0;
                break;  // We found the matching vertex, no need to continue the inner loop
            }
        }
    }

    // Compute cotangent Laplacian
    Eigen::SparseMatrix<double> L;
    igl::cotmatrix(V, F, L);
    std::cout << "L dimensions: " << L.rows() << "x" << L.cols() << std::endl;

    // Gradient descent parameters
    double learning_rate = 1e-6;
    int max_iterations = 10000;
    double tolerance = 1e-6;

    double prev_objective = std::numeric_limits<double>::max();

    for (int iter = 0; iter < max_iterations; ++iter) {
        // Compute gradient of g
        Eigen::MatrixXd grad_g = Eigen::Map<const Eigen::MatrixXd>((G*g).eval().data(),F.rows(),3);

        // std::cout << "grad_g dimensions: " << grad_g.rows() << "x" << grad_g.cols() << std::endl;

        // Compute objective function
        Eigen::VectorXd objective = (GF.array() * grad_g.array() - 1).square().rowwise().sum();

        // std::cout << "objective dimensions: " << objective.rows() << "x" << objective.cols() << std::endl;

        double total_objective = objective.sum();

        // Compute gradient of objective function
        Eigen::MatrixXd temp = (GF.array() * (GF.array() * grad_g.array() - 1)).matrix().reshaped();
        Eigen::VectorXd grad_objective = 2 * G.transpose() * Eigen::Map<Eigen::VectorXd>(temp.data(), F.rows() * 3, 1);
        // std::cout << "grad_objective dimensions: " << grad_objective.rows() << "x" << grad_objective.cols() << std::endl;

        // Update g
        g -= learning_rate * grad_objective;

        // Set g(B) = 0 again
        for (int i = 0; i < B.rows(); ++i) {
            // Find the index of this point in V
            for (int j = 0; j < V.rows(); ++j) {
                if (B.row(i) == V.row(j)) {
                    g(j) = 0;
                    break;  // We found the matching vertex, no need to continue the inner loop
                }
            }
        }

        // Check for NaN values
        if (g.array().isNaN().any()) {
            std::cout << "NaN values detected in g at iteration " << iter << std::endl;
            break;
        }

        // Print debug information
        if (iter % 100 == 0) {
            std::cout << "Iteration " << iter << ", Objective: " << total_objective << std::endl;
            std::cout << "g min: " << g.minCoeff() << ", g max: " << g.maxCoeff() << std::endl;
        }

        // Check for convergence
        if (iter > 0 && std::abs(total_objective - prev_objective) < tolerance) {
            std::cout << "Converged at iteration " << iter << std::endl;
            break;
        }

        prev_objective = total_objective;
    }

    return g;
}
