#include "compute_g.h"
#include <igl/cotmatrix.h>
#include <igl/grad.h>

Eigen::VectorXd compute_g(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    const Eigen::VectorXd& f,
    const Eigen::MatrixXd& GF,
    const Eigen::MatrixXd& B)
{
    // Initialize g with random values
    Eigen::VectorXd g = Eigen::VectorXd::Random(V.rows());

    // Set g(B) = 0
    for (int i = 0; i < B.rows(); ++i) {
        int closest_vertex = -1;
        double min_distance = std::numeric_limits<double>::max();
        for (int j = 0; j < V.rows(); ++j) {
            double distance = (V.row(j) - B.row(i)).norm();
            if (distance < min_distance) {
                min_distance = distance;
                closest_vertex = j;
            }
        }
        if (closest_vertex != -1) {
            g(closest_vertex) = 0;
        }
    }

    // Compute cotangent Laplacian
    Eigen::SparseMatrix<double> L;
    igl::cotmatrix(V, F, L);

    // Compute gradient operator
    Eigen::SparseMatrix<double> G;
    igl::grad(V, F, G);

    // Gradient descent parameters
    double learning_rate = 0.01;
    int max_iterations = 1000;
    double tolerance = 1e-6;

    for (int iter = 0; iter < max_iterations; ++iter) {
        // Compute gradient of g
        Eigen::MatrixXd grad_g = G * g;

        // Compute objective function
        Eigen::VectorXd objective = (GF.array() * grad_g.array() - 1).square().rowwise().sum();
        double total_objective = objective.sum();

        // Compute gradient of objective function
        Eigen::VectorXd grad_objective = 2 * G.transpose() * (GF.array() * (GF.array() * grad_g.array() - 1)).matrix();

        // Update g
        g -= learning_rate * grad_objective;

        // Set g(B) = 0 again
        for (int i = 0; i < B.rows(); ++i) {
            int closest_vertex = -1;
            double min_distance = std::numeric_limits<double>::max();
            for (int j = 0; j < V.rows(); ++j) {
                double distance = (V.row(j) - B.row(i)).norm();
                if (distance < min_distance) {
                    min_distance = distance;
                    closest_vertex = j;
                }
            }
            if (closest_vertex != -1) {
                g(closest_vertex) = 0;
            }
        }

        // Check for convergence
        if (iter > 0 && std::abs(total_objective - prev_objective) < tolerance) {
            break;
        }

        prev_objective = total_objective;
    }

    return g;
}
