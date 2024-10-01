#include <Eigen/Core>
#include <Eigen/SparseCore>
namespace igl
{
  template <typename Scalar>
  struct HeatGeodesicsData;
}

bool update(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const double t,
  const double x,
  const double y,
  const Eigen::Matrix4f& model,
  const Eigen::Matrix4f& proj,
  const Eigen::Vector4f& viewport,
  const igl::HeatGeodesicsData<double>& data,
  Eigen::VectorXd& D,
  const Eigen::SparseMatrix<double>& G,
  Eigen::MatrixXd& GF,
  Eigen::MatrixXd& BaryCenter,
  Eigen::MatrixXd& J_Delta_F_arrow,
  Eigen::MatrixXd& Cut_Path,
  Eigen::MatrixXd & V_cut,
  Eigen::MatrixXi & F_cut);
