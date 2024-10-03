#include "isolines_colormap.h"
#include "update.h"
#include <igl/read_triangle_mesh.h>
#include <igl/heat_geodesics.h>
#include <igl/avg_edge_length.h>
#include <igl/opengl/glfw/Viewer.h>
#include <iostream>
#include <igl/avg_edge_length.h>
#include <igl/grad.h>
#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>

int main(int argc, char *argv[])
{
  Eigen::MatrixXi F;
  Eigen::MatrixXd V;
  Eigen::VectorXd D;
  Eigen::VectorXd g;
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

  // Attach a menu plugin
  igl::opengl::glfw::imgui::ImGuiPlugin plugin;
  viewer.plugins.push_back(&plugin);
  igl::opengl::glfw::imgui::ImGuiMenu menu;
  plugin.widgets.push_back(&menu);

  bool down_on_mesh = false;
  const auto update = [&]()->bool
  {
    const double x = viewer.current_mouse_x;
    const double y = viewer.core().viewport(3) - viewer.current_mouse_y;
    Eigen::MatrixXd GF;
    Eigen::SparseMatrix<double> G;
    Eigen::MatrixXd BaryCenter;
    Eigen::MatrixXd J_Delta_F_arrow;
    Eigen::MatrixXd Cut_Path;
    Eigen::MatrixXd V_cut;
    Eigen::MatrixXi F_cut;
    Eigen::MatrixXd B;

    igl::grad(V,F,G);
    if(::update(
      V,F,t,x,y,
      viewer.core().view,viewer.core().proj,viewer.core().viewport,
      data,
      D,
      G,
      GF,
      BaryCenter,
      J_Delta_F_arrow,
      Cut_Path,
      V_cut,
      F_cut,
      B,
      g))
    {
      viewer.data().clear();
      viewer.data().set_mesh(V, F);

      viewer.data().set_data(D);
      viewer.data().clear_edges();
      // Draw a black segment in direction of gradient at face barycenters
      const Eigen::RowVector3d black(0,0,0);
      const Eigen::RowVector3d red(1,0,0);
      const Eigen::RowVector3d green(0,1,0);
      const Eigen::RowVector3d blue(0,0,1);
      const Eigen::RowVector3d yellow(1,1,0);
      const Eigen::RowVector3d magenta(1,0,1);
      const Eigen::RowVector3d cyan(0,1,1);
      viewer.data().add_edges(BaryCenter,J_Delta_F_arrow, black);
      // viewer.data().add_edges(BaryCenter,g, green);

      viewer.data().add_edges(Cut_Path.topRows(Cut_Path.rows()-1), Cut_Path.bottomRows(Cut_Path.rows()-1), Eigen::RowVector3d(1,0,0));
      // viewer.data().add_edges(B.topRows(B.rows()-1), B.bottomRows(B.rows()-1), Eigen::RowVector3d(1,1,0));

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

  // Draw additional windows
  menu.callback_draw_custom_window = [&]()
  {
    // Define next window position + size
    ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 10), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(200, 160), ImGuiCond_FirstUseEver);
    ImGui::Begin(
        "New Window", nullptr,
        ImGuiWindowFlags_NoSavedSettings
    );

    // Expose the same variable directly ...
    ImGui::PushItemWidth(-80);
    // ImGui::DragScalar("double", ImGuiDataType_Double, &doubleVariable, 0.1, 0, 0, "%.4f");
    ImGui::PopItemWidth();

    // static std::string str = "bunny";
    // ImGui::InputText("Name", str);

    ImGui::End();
  };

  // Show mesh
  viewer.data().set_mesh(V, F);
  viewer.data().set_data(Eigen::VectorXd::Zero(V.rows()));
  viewer.data().set_colormap(isolines_colormap());
  viewer.data().show_lines = false;
  viewer.launch();

}
