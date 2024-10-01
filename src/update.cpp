#include "update.h"
#include <igl/heat_geodesics.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/grad.h>
#include <igl/barycenter.h>
#include <igl/per_vertex_normals.h>

bool update(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const double t,
  const double x,
  const double y,
  const Eigen::Matrix4f& view,
  const Eigen::Matrix4f& proj,
  const Eigen::Vector4f& viewport,
  const igl::HeatGeodesicsData<double>& data,
  Eigen::VectorXd& D,
  const Eigen::SparseMatrix<double>& G,
  Eigen::MatrixXd& GF,
  Eigen::MatrixXd& BaryCenter,
  Eigen::MatrixXd& J_Delta_F_arrow,
  Eigen::MatrixXd& Cut_Path)
{
  int fid;
  Eigen::Vector3f bc;
  igl::barycenter(V,F,BaryCenter);
  // Cast a ray in the view direction starting from the mouse position
  if(igl::unproject_onto_mesh(Eigen::Vector2f(x,y), view,
    proj, viewport, V, F, fid, bc))
  {
    // if big mesh, just use closest vertex. Otherwise, blend distances to
    // vertices of face using barycentric coordinates.
    if(F.rows()>100000)
    {
      // 3d position of hit
      const Eigen::RowVector3d m3 =
        V.row(F(fid,0))*bc(0) + V.row(F(fid,1))*bc(1) + V.row(F(fid,2))*bc(2);
      int cid = 0;
      Eigen::Vector3d(
          (V.row(F(fid,0))-m3).squaredNorm(),
          (V.row(F(fid,1))-m3).squaredNorm(),
          (V.row(F(fid,2))-m3).squaredNorm()).minCoeff(&cid);
      const int vid = F(fid,cid);
      igl::heat_geodesics_solve(data,(Eigen::VectorXi(1,1)<<vid).finished(),D);
    }else
    {
      D = Eigen::VectorXd::Zero(V.rows());
      for(int cid = 0;cid<3;cid++)
      {
        const int vid = F(fid,cid);
        Eigen::VectorXd Dc;
        igl::heat_geodesics_solve(data,(Eigen::VectorXi(1,1)<<vid).finished(),Dc);
        D += Dc*bc(cid);
      }
    }

    // Find the vertex for the smallest and largest D
    Eigen::MatrixXd extreme_vertices(2, 3);
    int min_index, max_index;
    D.minCoeff(&min_index);
    D.maxCoeff(&max_index);
    extreme_vertices.row(0) = V.row(min_index);
    extreme_vertices.row(1) = V.row(max_index);

    std::cout << "Extreme vertices:" << std::endl;
    std::cout << "Min: " << extreme_vertices.row(0) << std::endl;
    std::cout << "Max: " << extreme_vertices.row(1) << std::endl;

    GF = Eigen::Map<const Eigen::MatrixXd>((G*D).eval().data(),F.rows(),3);
    Eigen::MatrixXd N_vertices;
    igl::per_vertex_normals(V, F, N_vertices);

    Eigen::MatrixXd GF_rotated = Eigen::MatrixXd::Zero(F.rows(), 3);
    for (int i = 0; i < F.rows(); i++) {
        for (int j = 0; j < 3; j++) {
            int vid = F(i, j);
            Eigen::Vector3d normal = N_vertices.row(vid);
            Eigen::Vector3d gradient = GF.row(i);

            // Project gradient onto tangent plane
            Eigen::Vector3d tangent_gradient = gradient - gradient.dot(normal) * normal;

            // Rotate by π/2 in tangent plane
            Eigen::Vector3d rotated = normal.cross(tangent_gradient);

            GF_rotated.row(i) += rotated / 3.0; // Average over vertices of the face
        }
    }

    // Use GF_rotated instead of GF for visualization
    const Eigen::VectorXd GF_mag = GF_rotated.rowwise().norm();

    // Average edge length divided by average gradient (for scaling)
    const double max_size = igl::avg_edge_length(V,F) / GF_mag.mean();
    J_Delta_F_arrow = BaryCenter+max_size*GF_rotated;
    return true;
  }
  return false;
}
