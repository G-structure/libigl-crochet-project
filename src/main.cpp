#include "isolines_colormap.h"
#include "update.h"
#include <igl/read_triangle_mesh.h>
#include <igl/heat_geodesics.h>
#include <igl/avg_edge_length.h>
#include <igl/opengl/glfw/Viewer.h>
#include <iostream>
#include <igl/avg_edge_length.h>
#include <igl/barycenter.h>
#include <igl/grad.h>

int main(int argc, char *argv[])
{
  Eigen::MatrixXi F;
  Eigen::MatrixXd V;
  igl::read_triangle_mesh( argc>1?argv[1]: TUTORIAL_SHARED_PATH "/sphere.obj",V,F);
  double t = std::pow(igl::avg_edge_length(V,F),2);

  // Precomputation
  igl::HeatGeodesicsData<double> data;
  const auto precompute = [&]()
  {
    if(!igl::heat_geodesics_precompute(V,F,t,data))
    {
      std::cerr<<"Error: heat_geodesics_precompute failed."<<std::endl;
      exit(EXIT_FAILURE);
    };
  };
  precompute();

  igl::opengl::glfw::Viewer viewer;
  bool down_on_mesh = false;
  const auto update = [&]()->bool
  {
    const double x = viewer.current_mouse_x;
    const double y = viewer.core().viewport(3) - viewer.current_mouse_y;
    Eigen::VectorXd D;
    Eigen::MatrixXd GF;
    Eigen::SparseMatrix<double> G;
    igl::grad(V,F,G);
    if(::update(
      V,F,t,x,y,
      viewer.core().view,viewer.core().proj,viewer.core().viewport,
      data,
      D,
      G,
      GF))
    {
      ## TODO need to move GF_mag, max_size into the update function
      ## TODO need to fix error where gradent arrows are not being cleared on update
      viewer.data().set_data(D);
      const Eigen::VectorXd GF_mag = GF.rowwise().norm();
      // Average edge length divided by average gradient (for scaling)
      const double max_size = igl::avg_edge_length(V,F) / GF_mag.mean();
      // Draw a black segment in direction of gradient at face barycenters
      Eigen::MatrixXd BC;
      igl::barycenter(V,F,BC);
      const Eigen::RowVector3d black(0,0,0);
      viewer.data().add_edges(BC,BC+max_size*GF, black);
      return true;
    }
    return false;
  };
  viewer.callback_mouse_down =
    [&](igl::opengl::glfw::Viewer& viewer, int, int)->bool
  {
    if(update())
    {
      down_on_mesh = true;
      return true;
    }
    return false;
  };
  viewer.callback_mouse_move =
    [&](igl::opengl::glfw::Viewer& viewer, int, int)->bool
    {
      if(down_on_mesh)
      {
        update();
        return true;
      }
      return false;
    };
  viewer.callback_mouse_up =
    [&down_on_mesh](igl::opengl::glfw::Viewer& viewer, int, int)->bool
  {
    down_on_mesh = false;
    return false;
  };
  std::cout<<R"(Usage:
  [click]  Click on shape to pick new geodesic distance source
  ,/.      Decrease/increase t by factor of 10.0
  D,d      Toggle using intrinsic Delaunay discrete differential operators

)";

  viewer.callback_key_pressed =
    [&](igl::opengl::glfw::Viewer& /*viewer*/, unsigned int key, int mod)->bool
  {
    switch(key)
    {
    default:
      return false;
    case 'D':
    case 'd':
      data.use_intrinsic_delaunay = !data.use_intrinsic_delaunay;
      std::cout<<(data.use_intrinsic_delaunay?"":"not ")<<
        "using intrinsic delaunay..."<<std::endl;
      precompute();
      update();
      break;
    case '.':
    case ',':
      t *= (key=='.'?10.0:0.1);
      precompute();
      update();
      std::cout<<"t: "<<t<<std::endl;
      break;
    }
    return true;
  };

  // Show mesh
  viewer.data().set_mesh(V, F);
  viewer.data().set_data(Eigen::VectorXd::Zero(V.rows()));
  viewer.data().set_colormap(isolines_colormap());
  viewer.data().show_lines = false;
  viewer.launch();

}
