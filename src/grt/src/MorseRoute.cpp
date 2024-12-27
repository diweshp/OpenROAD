#include "MorseRoute.h"

#include <algorithm>
#include <cmath>
#include <unordered_set>

namespace grt {
MorseRoute::MorseRoute(odb::dbDatabase* db)
    : db_(db),
      max_degree_(0),
      overflow_iterations_(0),
      congestion_report_iter_step_(0),
      x_range_(0),
      y_range_(0),
      num_adjust_(0),
      v_capacity_(0),
      h_capacity_(0),
      x_grid_(0),
      y_grid_(0),
      x_grid_max_(0),
      y_grid_max_(0),
      x_corner_(0),
      y_corner_(0),
      tile_size_(0),
      enlarge_(0),
      costheight_(0),
      ahth_(0),
      num_layers_(0),
      total_overflow_(0),
      has_2D_overflow_(false),
      grid_hv_(0),
      verbose_(false),
      critical_nets_percentage_(10),
      via_cost_(0),
      mazeedge_threshold_(0),
      v_capacity_lb_(0),
      h_capacity_lb_(0),
      regular_x_(false),
      regular_y_(false)
{
}

void MorseRoute::clear()
{
  clearNets();

  num_adjust_ = 0;
  v_capacity_ = 0;
  h_capacity_ = 0;
  total_overflow_ = 0;
  has_2D_overflow_ = false;

  h_edges_.resize(boost::extents[0][0]);
  v_edges_.resize(boost::extents[0][0]);
  seglist_.clear();

  gxs_.clear();
  gys_.clear();
  gs_.clear();

  tree_order_pv_.clear();
  tree_order_cong_.clear();

  h_edges_3D_.resize(boost::extents[0][0][0]);
  v_edges_3D_.resize(boost::extents[0][0][0]);

  parent_x1_.resize(boost::extents[0][0]);
  parent_y1_.resize(boost::extents[0][0]);
  parent_x3_.resize(boost::extents[0][0]);
  parent_y3_.resize(boost::extents[0][0]);

  net_eo_.clear();

  xcor_.clear();
  ycor_.clear();
  dcor_.clear();

  hv_.resize(boost::extents[0][0]);
  hyper_v_.resize(boost::extents[0][0]);
  hyper_h_.resize(boost::extents[0][0]);
  corr_edge_.resize(boost::extents[0][0]);

  in_region_.resize(boost::extents[0][0]);

  v_capacity_3D_.clear();
  h_capacity_3D_.clear();

  cost_hvh_.clear();
  cost_vhv_.clear();
  cost_h_.clear();
  cost_v_.clear();
  cost_lr_.clear();
  cost_tb_.clear();
  cost_hvh_test_.clear();
  cost_v_test_.clear();
  cost_tb_test_.clear();

  vertical_blocked_intervals_.clear();
  horizontal_blocked_intervals_.clear();
}

void MorseRoute::clearNets()
{
  if (!sttrees_.empty()) {
    sttrees_.clear();
  }

  for (MorseNet* net : nets_) {
    delete net;
  }
  nets_.clear();
  net_ids_.clear();
  seglist_.clear();
  db_net_id_map_.clear();
}

void MorseRoute::clearNetRoute(odb::dbNet* db_net)
{
  if (db_net_id_map_.find(db_net) != db_net_id_map_.end()) {
    const int net_id = db_net_id_map_[db_net];
    clearNetRoute(net_id);
  }
}

void MorseRoute::getNetId(odb::dbNet* db_net, int& net_id, bool& exists)
{
  auto itr = db_net_id_map_.find(db_net);
  exists = itr != db_net_id_map_.end();
  net_id = exists ? itr->second : 0;
}

void MorseRoute::clearNetRoute(const int netID)
{
  // clear used resources for the net route
  //releaseNetResources(netID);

  // clear stree
  sttrees_[netID].nodes.clear();
  sttrees_[netID].edges.clear();
}
void MorseRoute::setGridsAndLayers(int x, int y, int nLayers)
{
  x_grid_ = x;
  y_grid_ = y;
  num_layers_ = nLayers;
  layer_directions_.resize(num_layers_);
  if (std::max(x_grid_, y_grid_) >= 1000) {
    x_range_ = std::max(x_grid_, y_grid_);
    y_range_ = std::max(x_grid_, y_grid_);
  } else {
    x_range_ = 1000;
    y_range_ = 1000;
  }

  v_capacity_3D_.resize(num_layers_);
  h_capacity_3D_.resize(num_layers_);
  last_col_v_capacity_3D_.resize(num_layers_);
  last_row_h_capacity_3D_.resize(num_layers_);

  for (int i = 0; i < num_layers_; i++) {
    v_capacity_3D_[i] = 0;
    h_capacity_3D_[i] = 0;
    last_col_v_capacity_3D_[i] = 0;
    last_row_h_capacity_3D_[i] = 0;
  }

  hv_.resize(boost::extents[y_range_][x_range_]);
  hyper_v_.resize(boost::extents[y_range_][x_range_]);
  hyper_h_.resize(boost::extents[y_range_][x_range_]);
  corr_edge_.resize(boost::extents[y_range_][x_range_]);

  in_region_.resize(boost::extents[y_range_][x_range_]);

  cost_hvh_.resize(x_range_);  // Horizontal first Z
  cost_vhv_.resize(y_range_);  // Vertical first Z
  cost_h_.resize(y_range_);    // Horizontal segment cost
  cost_v_.resize(x_range_);    // Vertical segment cost
  cost_lr_.resize(y_range_);   // Left and right boundary cost
  cost_tb_.resize(x_range_);   // Top and bottom boundary cost

  cost_hvh_test_.resize(y_range_);  // Vertical first Z
  cost_v_test_.resize(x_range_);    // Vertical segment cost
  cost_tb_test_.resize(x_range_);   // Top and bottom boundary cost

  // maze3D variables
  directions_3D_.resize(boost::extents[num_layers_][y_grid_][x_grid_]);
  corr_edge_3D_.resize(boost::extents[num_layers_][y_grid_][x_grid_]);
  pr_3D_.resize(boost::extents[num_layers_][y_grid_][x_grid_]);

  int64 total_size = static_cast<int64>(num_layers_) * y_range_ * x_range_;
  pop_heap2_3D_.resize(total_size, false);

  // allocate memory for priority queue
  total_size = static_cast<int64>(y_grid_) * x_grid_ * num_layers_;
  src_heap_3D_.resize(total_size);
  dest_heap_3D_.resize(total_size);

  d1_3D_.resize(boost::extents[num_layers_][y_range_][x_range_]);
  d2_3D_.resize(boost::extents[num_layers_][y_range_][x_range_]);
}

void MorseRoute::addVCapacity(short verticalCapacity, int layer)
{
  v_capacity_3D_[layer - 1] = verticalCapacity;
  v_capacity_ += v_capacity_3D_[layer - 1];
}

void MorseRoute::addHCapacity(short horizontalCapacity, int layer)
{
  h_capacity_3D_[layer - 1] = horizontalCapacity;
  h_capacity_ += h_capacity_3D_[layer - 1];
}

void MorseRoute::setLowerLeft(int x, int y)
{
  x_corner_ = x;
  y_corner_ = y;
}

void MorseRoute::setTileSize(int size)
{
  tile_size_ = size;
}

void MorseRoute::addLayerDirection(int layer_idx,
                                   const odb::dbTechLayerDir& direction)
{
  layer_directions_[layer_idx] = direction;
}

void MorseRoute::initEdges()
{
  const float LB = 0.9;
  v_capacity_lb_ = LB * v_capacity_;
  h_capacity_lb_ = LB * h_capacity_;

  // allocate memory and initialize for edges

  h_edges_.resize(boost::extents[y_grid_][x_grid_ - 1]);
  v_edges_.resize(boost::extents[y_grid_ - 1][x_grid_]);

  v_edges_3D_.resize(boost::extents[num_layers_][y_grid_][x_grid_]);
  h_edges_3D_.resize(boost::extents[num_layers_][y_grid_][x_grid_]);

  for (int i = 0; i < y_grid_; i++) {
    for (int j = 0; j < x_grid_; j++) {
      // 2D edge initialization
      if (j < x_grid_ - 1) {
        h_edges_[i][j].cap = h_capacity_;
        h_edges_[i][j].usage = 0;
        h_edges_[i][j].est_usage = 0;
        h_edges_[i][j].red = 0;
        h_edges_[i][j].last_usage = 0;
      }

      // 3D edge initialization
      for (int k = 0; k < num_layers_; k++) {
        h_edges_3D_[k][i][j].cap = h_capacity_3D_[k];
        h_edges_3D_[k][i][j].usage = 0;
        h_edges_3D_[k][i][j].red = 0;
      }
    }
  }
  for (int i = 0; i < y_grid_; i++) {
    for (int j = 0; j < x_grid_; j++) {
      // 2D edge initialization
      if (i < y_grid_ - 1) {
        v_edges_[i][j].cap = v_capacity_;
        v_edges_[i][j].usage = 0;
        v_edges_[i][j].est_usage = 0;
        v_edges_[i][j].red = 0;
        v_edges_[i][j].last_usage = 0;
      }

      // 3D edge initialization
      for (int k = 0; k < num_layers_; k++) {
        v_edges_3D_[k][i][j].cap = v_capacity_3D_[k];
        v_edges_3D_[k][i][j].usage = 0;
        v_edges_3D_[k][i][j].red = 0;
      }
    }
  }
}

void MorseRoute::setNumAdjustments(int nAdjustments)
{
  num_adjust_ = nAdjustments;
}

void MorseRoute::setMaxNetDegree(int deg)
{
  max_degree_ = deg;
}

MorseNet* MorseRoute::addNet(odb::dbNet* db_net,
                             bool is_clock,
                             int driver_idx,
                             int cost,
                             int min_layer,
                             int max_layer,
                             float slack,
                             std::vector<int>* edge_cost_per_layer)
{
  int netID;
  bool exists;
  MorseNet* net;
  getNetId(db_net, netID, exists);
  if (exists) {
    net = nets_[netID];
    clearNetRoute(netID);
    seglist_[netID].clear();
  } else {
    net = new MorseNet;
    nets_.push_back(net);
    netID = nets_.size() - 1;
    db_net_id_map_[db_net] = netID;
    // at most (2*num_pins-2) nodes -> (2*num_pins-3) segs_ for a net
    seglist_.emplace_back();
    sttrees_.emplace_back();
    gxs_.emplace_back();
    gys_.emplace_back();
    gs_.emplace_back();
  }
  net->reset(db_net,
             is_clock,
             driver_idx,
             cost,
             min_layer,
             max_layer,
             slack,
             edge_cost_per_layer);
  net_ids_.push_back(netID);

  return net;
}
void MorseRoute::setVerbose(bool v)
{
  verbose_ = v;
}

void MorseRoute::setCriticalNetsPercentage(float u)
{
  critical_nets_percentage_ = u;
}

/*void FastRouteCore::setMakeWireParasiticsBuilder(
    AbstractMakeWireParasitics* builder)
{
  parasitics_builder_ = builder;
}*/

void MorseRoute::setOverflowIterations(int iterations)
{
  overflow_iterations_ = iterations;
}

void MorseRoute::setCongestionReportIterStep(int congestion_report_iter_step)
{
  congestion_report_iter_step_ = congestion_report_iter_step;
}

void MorseRoute::setCongestionReportFile(const char* congestion_file_name)
{
  congestion_file_name_ = congestion_file_name;
}

void MorseNet::addPin(int x, int y, int layer)
{
  pin_x_.push_back(x);
  pin_y_.push_back(y);
  pin_l_.push_back(layer);
}

}  // namespace grt