#ifndef COMPUTE_G_H
#define COMPUTE_G_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

Eigen::VectorXd compute_g(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    const Eigen::VectorXd& f,
    const Eigen::MatrixXd& GF,
    const Eigen::MatrixXd& B,
    const Eigen::SparseMatrix<double>& G);

#endif
