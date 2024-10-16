#include "update.h"
#include "compute_g.h"
#include <igl/heat_geodesics.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/grad.h>
#include <igl/barycenter.h>
#include <igl/per_vertex_normals.h>
#include <igl/cut_mesh.h>
#include <igl/boundary_loop.h>

Eigen::MatrixXi path_to_edges(const Eigen::MatrixXd& Cut_Path) {
    Eigen::MatrixXi edges(Cut_Path.rows() - 1, 2);
    for (int i = 0; i < Cut_Path.rows() - 1; ++i) {
        edges(i, 0) = i;
        edges(i, 1) = i + 1;
    }
    return edges;
}

Eigen::VectorXd computePerVertexGradients(
    const Eigen::MatrixXd& V, const Eigen::MatrixXi& F,
    const Eigen::MatrixXd& perFaceGradients) {

    Eigen::VectorXd perVertexGradients = Eigen::VectorXd::Zero(V.rows() * 3);
    Eigen::VectorXd vertexAreas = Eigen::VectorXd::Zero(V.rows());

    for (int i = 0; i < F.rows(); i++) {
        Eigen::Vector3d faceGradient = perFaceGradients.row(i);
        // Compute face area
        Eigen::Vector3d v1 = V.row(F(i, 1)) - V.row(F(i, 0));
        Eigen::Vector3d v2 = V.row(F(i, 2)) - V.row(F(i, 0));
        double faceArea = 0.5 * (v1.cross(v2)).norm();

        for (int j = 0; j < 3; j++) {
            int vId = F(i, j);
            perVertexGradients.segment<3>(vId * 3) += faceGradient * faceArea;
            vertexAreas(vId) += faceArea;
        }
    }

    // Normalize
    for (int i = 0; i < V.rows(); i++) {
        if (vertexAreas(i) > 0) {
            perVertexGradients.segment<3>(i * 3) /= vertexAreas(i);
        }
    }

    return perVertexGradients;
}

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
  Eigen::MatrixXd& Cut_Path,
  Eigen::MatrixXd & V_cut,
  Eigen::MatrixXi & F_cut,
  Eigen::MatrixXd & B,
  Eigen::VectorXd & g)
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

    // Find the shortest path from min to max vertex
    Eigen::VectorXi prev(V.rows());
    Eigen::VectorXd min_dist(V.rows());
    std::priority_queue<std::pair<double, int>,
                        std::vector<std::pair<double, int>>,
                        std::greater<std::pair<double, int>>> pq;

    // Initialize distances
    min_dist.setConstant(std::numeric_limits<double>::infinity());
    min_dist(min_index) = 0;
    pq.push({0, min_index});

    // Dijkstra's algorithm
    while (!pq.empty()) {
        int u = pq.top().second;
        pq.pop();

        if (u == max_index) break;

        for (int i = 0; i < F.rows(); i++) {
        for (int j = 0; j < 3; j++) {
            if (F(i, j) == u) {
            for (int k = 0; k < 3; k++) {
                int v = F(i, k);
                double weight = (V.row(u) - V.row(v)).norm();
                if (min_dist(u) + weight < min_dist(v)) {
                min_dist(v) = min_dist(u) + weight;
                prev(v) = u;
                pq.push({min_dist(v), v});
                }
            }
            }
        }
        }
    }

    // Reconstruct path
    std::vector<int> path;
    for (int v = max_index; v != min_index; v = prev(v)) {
        path.push_back(v);
    }
    path.push_back(min_index);
    std::reverse(path.begin(), path.end());

    // Convert path to Cut_Path
    Cut_Path.resize(path.size(), 3);
    for (int i = 0; i < path.size(); i++) {
        Cut_Path.row(i) = V.row(path[i]);
    }

    std::cout << "Cut Path:" << std::endl;
    for (int i = 0; i < Cut_Path.rows(); i++) {
        std::cout << Cut_Path.row(i) << std::endl;
    }

    // Convert Cut_Path to edges
    Eigen::MatrixXi cut_edges = path_to_edges(Cut_Path);

    std::cout << "Cut Edges:" << std::endl;
    for (int i = 0; i < cut_edges.rows(); i++) {
        std::cout << "Edge " << i << ": " << cut_edges(i, 0) << " - " << cut_edges(i, 1) << std::endl;
    }

    // Create cuts matrix
    Eigen::MatrixXi cuts = Eigen::MatrixXi::Zero(F.rows(), 3);
    for (int i = 0; i < cut_edges.rows(); ++i) {
        for (int j = 0; j < F.rows(); ++j) {
            for (int k = 0; k < 3; ++k) {
                int v1 = F(j, k);
                int v2 = F(j, (k + 1) % 3);
                if ((v1 == cut_edges(i, 0) && v2 == cut_edges(i, 1)) ||
                    (v1 == cut_edges(i, 1) && v2 == cut_edges(i, 0))) {
                    cuts(j, k) = 1;
                }
            }
        }
    }

    // Cut the mesh
    V_cut = V;
    F_cut = F;
    Eigen::VectorXi I;
    igl::cut_mesh(V_cut, F_cut, cuts, I);

    std::cout << "Mesh cut along the specified edges." << std::endl;
    std::cout << "Original mesh: " << V.rows() << " vertices, " << F.rows() << " faces." << std::endl;
    std::cout << "Cut mesh: " << V_cut.rows() << " vertices, " << F_cut.rows() << " faces." << std::endl;

    // Create D_cut by mapping D values to the cut mesh vertices
    Eigen::VectorXd D_cut(V_cut.rows());
    for (int i = 0; i < V_cut.rows(); ++i) {
        D_cut(i) = D(I(i));
    }

    std::cout << "Original D size: " << D.size() << std::endl;
    std::cout << "D_cut size: " << D_cut.size() << std::endl;

    // TODO: Calculate B: the longest connected path on Cut_Path along which D_cut is strictly monotone
    B = Cut_Path;

    std::cout << "Longest monotone path size: " << B.rows() << std::endl;
    std::cout << "Cut path size: " << Cut_Path.rows() << std::endl;

    // Output the values for D_cut along B
    std::cout << "D_cut values along B:" << std::endl;
    for (int i = 0; i < B.rows(); ++i) {
        // Find the index of the vertex in V_cut that matches B.row(i)
        int vertex_index = -1;
        for (int j = 0; j < V_cut.rows(); ++j) {
            if (V_cut.row(j) == B.row(i)) {
                vertex_index = j;
                break;
            }
        }
        if (vertex_index != -1) {
            std::cout << "Vertex " << i << ": " << D_cut(vertex_index) << std::endl;
        } else {
            std::cout << "Vertex " << i << ": Not found in V_cut" << std::endl;
        }
    }
    GF = Eigen::Map<const Eigen::MatrixXd>((G*D).eval().data(),F.rows(),3);
    Eigen::MatrixXd GF_vertex = Eigen::Map<Eigen::MatrixXd>(computePerVertexGradients(V, F, GF).data(), V.rows(), 3);
    std::cout << "Shape of: " << GF_vertex.rows() << " x " << GF_vertex.cols() << std::endl;
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

    g = compute_g(V, F, D, GF_rotated, B, G);

    std::cout << "Computed g. Size: " << g.size() << std::endl;
    std::cout << "g values:" << std::endl;
    for (int i = 0; i < g.size(); ++i) {
        std::cout << "g[" << i << "] = " << g[i] << std::endl;
    }
    // const Eigen::VectorXd g_mag = g.rowwise().norm();
    // const double max_size_g = igl::avg_edge_length(V,F) / g_mag.mean();
    // g = BaryCenter+max_size_g*g;

    // Average edge length divided by average gradient (for scaling)
    const double max_size_jf = igl::avg_edge_length(V,F) / GF_mag.mean();
    J_Delta_F_arrow = BaryCenter+max_size_jf*GF_rotated;
    // D = D_cut;
    std::cout << "Shape of G: " << G.rows() << " x " << G.cols() << std::endl;
    std::cout << "Shape of GF: " << GF.rows() << " x " << GF.cols() << std::endl;
    std::cout << "Shape of GF_rotated: " << GF_rotated.rows() << " x " << GF_rotated.cols() << std::endl;

    // Write D and g values to a CSV file
    std::ofstream csvFile("D_and_g_values.csv");
    csvFile << "D,g" << std::endl;
    for (int i = 0; i < D.size(); ++i) {
        csvFile << D[i] << "," << g[i] << std::endl;
    }
    csvFile.close();
    std::cout << "D and g values written to D_and_g_values.csv" << std::endl;
    return true;
  }
  return false;
}
