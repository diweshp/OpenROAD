/////////////////////////////////////////////////////////////////////////////
//
// BSD 3-Clause License
//
// Copyright (c) 2019, The Regents of the University of California
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//   list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its
//   contributors may be used to endorse or promote products derived from
//   this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
///////////////////////////////////////////////////////////////////////////////

#include <string>
#include <utility>
#include <vector>
#include <boost/functional/hash.hpp>
#include <boost/icl/interval.hpp>
#include <boost/icl/interval_set.hpp>
#include <boost/multi_array.hpp>
#include <set>
#include <unordered_map>

#include "grt/GRoute.h"
#include "grt/RoutePt.h"
#include "odb/db.h"
#include "MorseDataType.h"
#include "grt/GRoute.h"
#include "odb/geom.h"
#include "stt/SteinerTreeBuilder.h"

namespace odb {
class dbDatabase;
class dbTechLayerDir;
class dbTechLayer;
}  // namespace odb

namespace stt {
class SteinerTreeBuilder;
}


using boost::multi_array;
using boost::icl::interval;
using boost::icl::interval_set;
namespace grt {

   struct MorseParent3D
{
  int16_t layer;
  int x, y;
};

struct MorseCostParams
{
  const float logistic_coef;
  const float cost_height;
  const int slope;

  MorseCostParams(const float logistic_coef,
             const float cost_height,
             const int slope)
      : logistic_coef(logistic_coef), cost_height(cost_height), slope(slope)
  {
  }
};

class MorseRoute
{
 public:
  MorseRoute(odb::dbDatabase* db);
  int x_corner() const { return x_corner_; }
  int y_corner() const { return y_corner_; }
  int tile_size() const { return tile_size_; }
  void setNumAdjustments(int nAdjustments);
  void setMaxNetDegree(int deg);
  void initEdges();
  void clear();
  void clearNets();
  void setGridsAndLayers(int x, int y, int nLayers);
  void addVCapacity(short verticalCapacity, int layer);
  void addHCapacity(short horizontalCapacity, int layer);
  void setLowerLeft(int x, int y);
  void setTileSize(int size);
  void addLayerDirection(int layer_idx, const odb::dbTechLayerDir& direction);
  void setVerbose(bool v);
  void setCriticalNetsPercentage(float u);
  float getCriticalNetsPercentage() { return critical_nets_percentage_; };
  //void setMakeWireParasiticsBuilder(AbstractMakeWireParasitics* builder);
  void setOverflowIterations(int iterations);
  void setCongestionReportIterStep(int congestion_report_iter_step);
  void setCongestionReportFile(const char* congestion_file_name);
  void setGridMax(int x_max, int y_max);
  void getCongestionNets(std::set<odb::dbNet*>& congestion_nets);
  void computeCongestionInformation();
  void setRegularX(bool regular_x) { regular_x_ = regular_x; }
  void setRegularY(bool regular_y) { regular_y_ = regular_y; }
  void incrementEdge3DUsage(int x1, int y1, int x2, int y2, int layer);
  void setLastColVCapacity(short cap, int layer)
  {
    last_col_v_capacity_3D_[layer] = cap;
  }
  void setLastRowHCapacity(short cap, int layer)
  {
    last_row_h_capacity_3D_[layer] = cap;
  }
  const std::vector<int16_t>& getLastColumnVerticalCapacities()
  {
    return last_col_v_capacity_3D_;
  }
  const std::vector<int16_t>& getLastRowHorizontalCapacities()
  {
    return last_row_h_capacity_3D_;
  }
  const std::vector<short>& getVerticalCapacities() { return v_capacity_3D_; }
  const std::vector<short>& getHorizontalCapacities() { return h_capacity_3D_; }
  int getAvailableResources(int x1, int y1, int x2, int y2, int layer);
  int getEdgeCapacity(int x1, int y1, int x2, int y2, int layer);
  const multi_array<Edge3D, 3>& getHorizontalEdges3D() { return h_edges_3D_; }
  const multi_array<Edge3D, 3>& getVerticalEdges3D() { return v_edges_3D_; }
  void updateEdge2DAnd3DUsage(int x1,
                              int y1,
                              int x2,
                              int y2,
                              int layer,
                              int used);
  MorseNet* addNet(odb::dbNet* db_net,
                bool is_clock,
                int driver_idx,
                int cost,
                int min_layer,
                int max_layer,
                float slack,
                std::vector<int>* edge_cost_per_layer);

  
private:
  typedef std::tuple<int, int, int> Tile;
odb::dbDatabase* _db;
int max_degree_;
std::vector<MorseNet*> nets_;
  std::vector<int> cap_per_layer_;
  std::vector<int> usage_per_layer_;
  std::vector<int> overflow_per_layer_;
  std::vector<int> max_h_overflow_;
  std::vector<int> max_v_overflow_;
  int overflow_iterations_;
  int congestion_report_iter_step_;
  std::string congestion_file_name_;
  std::vector<odb::dbTechLayerDir> layer_directions_;
  int x_range_;
  int y_range_;
  int num_adjust_;
  int v_capacity_;
  int h_capacity_;
  int x_grid_;
  int y_grid_;
  int x_grid_max_;
  int y_grid_max_;
  int x_corner_;
  int y_corner_;
  int tile_size_;
  int enlarge_;
  int costheight_;
  int ahth_;
  int num_layers_;
  int total_overflow_;  // total # overflow
  bool has_2D_overflow_;
  int grid_hv_;
  bool verbose_;
  float critical_nets_percentage_;
  int via_cost_;
  int mazeedge_threshold_;
  float v_capacity_lb_;
  float h_capacity_lb_;
  bool regular_x_;
  bool regular_y_;
  std::vector<short> v_capacity_3D_;
  std::vector<short> h_capacity_3D_;
  std::vector<short> last_col_v_capacity_3D_;
  std::vector<short> last_row_h_capacity_3D_;
  std::vector<double> cost_hvh_;       // Horizontal first Z
  std::vector<double> cost_vhv_;       // Vertical first Z
  std::vector<double> cost_h_;         // Horizontal segment cost
  std::vector<double> cost_v_;         // Vertical segment cost
  std::vector<double> cost_lr_;        // Left and right boundary cost
  std::vector<double> cost_tb_;        // Top and bottom boundary cost
  std::vector<double> cost_hvh_test_;  // Vertical first Z
  std::vector<double> cost_v_test_;    // Vertical segment cost
  std::vector<double> cost_tb_test_;   // Top and bottom boundary cost
  std::vector<double> h_cost_table_;
  std::vector<double> v_cost_table_;
  std::vector<int> xcor_;
  std::vector<int> ycor_;
  std::vector<int> dcor_;
  // Maze 3D variables
  multi_array<MorseDirection, 3> directions_3D_;
  multi_array<int, 3> corr_edge_3D_;
  multi_array<MorseParent3D, 3> pr_3D_;
  std::vector<bool> pop_heap2_3D_;
  std::vector<int*> src_heap_3D_;
  std::vector<int*> dest_heap_3D_;
  multi_array<int, 3> d1_3D_;
  multi_array<int, 3> d2_3D_;

  multi_array<MorseEdge, 2> v_edges_;       // The way it is indexed is (Y, X)
  multi_array<MorseEdge, 2> h_edges_;       // The way it is indexed is (Y, X)
  multi_array<MorseEdge3D, 3> h_edges_3D_;  // The way it is indexed is (Layer, Y, X)
  multi_array<MorseEdge3D, 3> v_edges_3D_;  // The way it is indexed is (Layer, Y, X)
  multi_array<int, 2> corr_edge_;
  multi_array<short, 2> parent_x1_;
  multi_array<short, 2> parent_y1_;
  multi_array<short, 2> parent_x3_;
  multi_array<short, 2> parent_y3_;
  multi_array<bool, 2> hv_;
  multi_array<bool, 2> hyper_v_;
  multi_array<bool, 2> hyper_h_;
  multi_array<bool, 2> in_region_;

  std::unordered_map<odb::dbNet*, int> db_net_id_map_;  // db net -> net id
  std::vector<MorseOrderNetEdge> net_eo_;
  std::vector<std::vector<int>>
      gxs_;  // the copy of xs for nets, used for second FLUTE
  std::vector<std::vector<int>>
      gys_;  // the copy of xs for nets, used for second FLUTE
  std::vector<std::vector<int>>
      gs_;  // the copy of vertical sequence for nets, used for second FLUTE
  std::vector<std::vector<MorseSegment>> seglist_;  // indexed by netID, segID
  std::vector<MorseOrderNetPin> tree_order_pv_;
  std::vector<MorseOrderTree> tree_order_cong_;
  std::set<std::pair<int, int>> h_used_ggrid_;
  std::set<std::pair<int, int>> v_used_ggrid_;
  std::vector<int> net_ids_;

  std::vector<MorseStTree> sttrees_;  // the Steiner trees
  std::vector<MorseStTree> sttrees_bk_;

  utl::Logger* logger_;
  stt::SteinerTreeBuilder* stt_builder_;
  //AbstractMakeWireParasitics* parasitics_builder_;

  //std::unique_ptr<DebugSetting> debug_;

  std::unordered_map<Tile, interval_set<int>, boost::hash<Tile>>
      vertical_blocked_intervals_;
  std::unordered_map<Tile, interval_set<int>, boost::hash<Tile>>
      horizontal_blocked_intervals_;
  
};

}  // namespace grt
