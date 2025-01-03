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

void MorseRoute::setGridMax(int x_max, int y_max)
{
  x_grid_max_ = x_max;
  y_grid_max_ = y_max;
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
void MorseNet::reset(odb::dbNet* db_net,
                  bool is_clock,
                  int driver_idx,
                  int edge_cost,
                  int min_layer,
                  int max_layer,
                  float slack,
                  std::vector<int>* edge_cost_per_layer)
{
  db_net_ = db_net;
  is_critical_ = false;
  is_clock_ = is_clock;
  driver_idx_ = driver_idx;
  edge_cost_ = edge_cost;
  min_layer_ = min_layer;
  max_layer_ = max_layer;
  slack_ = slack;
  edge_cost_per_layer_.reset(edge_cost_per_layer);
  pin_x_.clear();
  pin_y_.clear();
  pin_l_.clear();
}
NetRouteMap MorseRoute::getRoutes()
{
  NetRouteMap routes;
  for (const int& netID : net_ids_) {
    odb::dbNet* db_net = nets_[netID]->getDbNet();
    GRoute& route = routes[db_net];
    std::unordered_set<GSegment, GSegmentHash> net_segs;

    const auto& treeedges = sttrees_[netID].edges;
    const int num_edges = sttrees_[netID].num_edges();

    for (int edgeID = 0; edgeID < num_edges; edgeID++) {
      const MorseTreeEdge* treeedge = &(treeedges[edgeID]);
      if (treeedge->len > 0 || treeedge->route.routelen > 0) {
        int routeLen = treeedge->route.routelen;
        const std::vector<short>& gridsX = treeedge->route.gridsX;
        const std::vector<short>& gridsY = treeedge->route.gridsY;
        const std::vector<short>& gridsL = treeedge->route.gridsL;
        int lastX = tile_size_ * (gridsX[0] + 0.5) + x_corner_;
        int lastY = tile_size_ * (gridsY[0] + 0.5) + y_corner_;
        int lastL = gridsL[0];
        for (int i = 1; i <= routeLen; i++) {
          const int xreal = tile_size_ * (gridsX[i] + 0.5) + x_corner_;
          const int yreal = tile_size_ * (gridsY[i] + 0.5) + y_corner_;

          GSegment segment
              = GSegment(lastX, lastY, lastL + 1, xreal, yreal, gridsL[i] + 1);

          lastX = xreal;
          lastY = yreal;
          lastL = gridsL[i];
          if (net_segs.find(segment) == net_segs.end()) {
            if (segment.init_layer != segment.final_layer) {
              GSegment invet_via = GSegment(segment.final_x,
                                            segment.final_y,
                                            segment.final_layer,
                                            segment.init_x,
                                            segment.init_y,
                                            segment.init_layer);
              if (net_segs.find(invet_via) != net_segs.end()) {
                continue;
              }
            }

            net_segs.insert(segment);
            route.push_back(segment);
          }
        }
      }
    }
  }

  return routes;
}

int MorseRoute::getOverflow2D(int* maxOverflow)
{
  // check 2D edges for invalid usage values
  //check2DEdgesUsage();

  // get overflow
  int H_overflow = 0;
  int V_overflow = 0;
  int max_H_overflow = 0;
  int max_V_overflow = 0;
  int hCap = 0;
  int vCap = 0;
  int numedges = 0;

  int total_usage = 0;

  for (const auto& [i, j] : h_used_ggrid_) {
    total_usage += h_edges_[i][j].est_usage;
    const int overflow = h_edges_[i][j].est_usage - h_edges_[i][j].cap;
    hCap += h_edges_[i][j].cap;
    if (overflow > 0) {
      H_overflow += overflow;
      max_H_overflow = std::max(max_H_overflow, overflow);
      numedges++;
    }
  }

  for (const auto& [i, j] : v_used_ggrid_) {
    total_usage += v_edges_[i][j].est_usage;
    const int overflow = v_edges_[i][j].est_usage - v_edges_[i][j].cap;
    vCap += v_edges_[i][j].cap;
    if (overflow > 0) {
      V_overflow += overflow;
      max_V_overflow = std::max(max_V_overflow, overflow);
      numedges++;
    }
  }

  const int max_overflow = std::max(max_H_overflow, max_V_overflow);
  total_overflow_ = H_overflow + V_overflow;
  *maxOverflow = max_overflow;

  if (total_usage > 800000) {
    ahth_ = 30;
  } else {
    ahth_ = 20;
  }

  if (logger_->debugCheck(GRT, "checkRoute3D", 1)) {
    logger_->report("Overflow report.");
    logger_->report("Total hCap               : {}", hCap);
    logger_->report("Total vCap               : {}", vCap);
    logger_->report("Total usage              : {}", total_usage);
    logger_->report("Max H overflow           : {}", max_H_overflow);
    logger_->report("Max V overflow           : {}", max_V_overflow);
    logger_->report("Max overflow             : {}", max_overflow);
    logger_->report("Number of overflow edges : {}", numedges);
    logger_->report("H   overflow             : {}", H_overflow);
    logger_->report("V   overflow             : {}", V_overflow);
    logger_->report("Final overflow           : {}\n", total_overflow_);
  }

  return total_overflow_;
}

void MorseRoute::InitEstUsage()
{
  for (int i = 0; i < y_grid_; i++) {
    for (int j = 0; j < x_grid_ - 1; j++) {
      h_edges_[i][j].est_usage = 0;
    }
  }

  for (int i = 0; i < y_grid_ - 1; i++) {
    for (int j = 0; j < x_grid_; j++) {
      v_edges_[i][j].est_usage = 0;
    }
  }
}

NetRouteMap MorseRoute::run()
{
  if (netCount() == 0) {
    return getRoutes();
  }

  v_used_ggrid_.clear();
  h_used_ggrid_.clear();

  int tUsage;
  int cost_step;
  int maxOverflow = 0;
  int minoflrnd = 0;
  int bwcnt = 0;

  // Init grid variables when debug mode is actived
  if (debug_->isOn()) {
    //fastrouteRender()->setGridVariables(tile_size_, x_corner_, y_corner_);
  }

  // TODO: check this size
  int max_degree2 = 2 * max_degree_;
  xcor_.resize(max_degree2);
  ycor_.resize(max_degree2);
  dcor_.resize(max_degree2);
  net_eo_.reserve(max_degree2);

  int THRESH_M = 20;
  const int ENLARGE = 15;  // 5
  const int ESTEP1 = 10;   // 10
  const int ESTEP2 = 5;    // 5
  const int ESTEP3 = 5;    // 5
  int CSTEP1 = 2;          // 5
  const int CSTEP2 = 2;    // 3
  const int CSTEP3 = 5;    // 15
  const int COSHEIGHT = 4;
  int L = 0;
  int VIA = 2;
  const int Ripvalue = -1;
  const bool goingLV = true;
  const bool noADJ = false;
  const int thStep1 = 10;
  const int thStep2 = 4;
  const int LVIter = 3;
  const int mazeRound = 500;
  int bmfl = BIG_INT;
  int minofl = BIG_INT;
  float logistic_coef = 0;
  int slope;
  int max_adj;

  // call FLUTE to generate RSMT and break the nets into segments (2-pin nets)

  via_cost_ = 0;
  gen_brk_RSMT(false, false, false, false, noADJ);
  getOverflow2D(&maxOverflow);
  


  //  past_cong = getOverflow2Dmaze( &maxOverflow);

  InitEstUsage();

  return routes;
}
int MorseRoute::getAvailableResources(int x1,
                                         int y1,
                                         int x2,
                                         int y2,
                                         int layer)
{
  const int k = layer - 1;
  int available_cap = 0;
  if (y1 == y2) {  // horizontal edge
    available_cap = h_edges_3D_[k][y1][x1].cap - h_edges_3D_[k][y1][x1].usage;
  } else if (x1 == x2) {  // vertical edge
    available_cap = v_edges_3D_[k][y1][x1].cap - v_edges_3D_[k][y1][x1].usage;
  } else {
    logger_->error(
        GRT,
        213,
        "Cannot get available resources: edge is not vertical or horizontal.");
  }
  return available_cap;
}

int MorseRoute::getEdgeCapacity(int x1, int y1, int x2, int y2, int layer)
{
  const int k = layer - 1;

  if (y1 == y2) {  // horizontal edge
    return h_edges_3D_[k][y1][x1].cap;
  } else if (x1 == x2) {  // vertical edge
    return v_edges_3D_[k][y1][x1].cap;
  } else {
    logger_->error(
        GRT,
        214,
        "Cannot get edge capacity: edge is not vertical or horizontal.");
    return 0;
  }
}

int MorseRoute::getEdgeCapacity(MorseNet* net,
                                   int x1,
                                   int y1,
                                   MorseEdgeDirection direction)
{
  int cap = 0;

  // get 2D edge capacity respecting layer restrictions
  for (int l = net->getMinLayer(); l <= net->getMaxLayer(); l++) {
    if (direction == MorseEdgeDirection::Horizontal) {
      cap += h_edges_3D_[l][y1][x1].cap;
    } else {
      cap += v_edges_3D_[l][y1][x1].cap;
    }
  }

  return cap;
}

void MorseRoute::fluteNormal(const int netID,
                                const std::vector<int>& x,
                                const std::vector<int>& y,
                                const int acc,
                                const float coeffV,
                                Tree& t)
{
  const int d = x.size();

  if (d == 2) {
    t.deg = 2;
    t.length = abs(x[0] - x[1]) + abs(y[0] - y[1]);
    t.branch.resize(2);
    t.branch[0].x = x[0];
    t.branch[0].y = y[0];
    t.branch[0].n = 1;
    t.branch[1].x = x[1];
    t.branch[1].y = y[1];
    t.branch[1].n = 1;
  } else if (d == 3) {
    t.deg = 3;
    int x_max, x_min, x_mid;
    if (x[0] < x[1]) {
      if (x[0] < x[2]) {
        x_min = x[0];
        x_mid = std::min(x[1], x[2]);
        x_max = std::max(x[1], x[2]);
      } else {
        x_min = x[2];
        x_mid = x[0];
        x_max = x[1];
      }
    } else {
      if (x[0] < x[2]) {
        x_min = x[1];
        x_mid = x[0];
        x_max = x[2];
      } else {
        x_min = std::min(x[1], x[2]);
        x_mid = std::max(x[1], x[2]);
        x_max = x[0];
      }
    }
    int y_max, y_min, y_mid;
    if (y[0] < y[1]) {
      if (y[0] < y[2]) {
        y_min = y[0];
        y_mid = std::min(y[1], y[2]);
        y_max = std::max(y[1], y[2]);
      } else {
        y_min = y[2];
        y_mid = y[0];
        y_max = y[1];
      }
    } else {
      if (y[0] < y[2]) {
        y_min = y[1];
        y_mid = y[0];
        y_max = y[2];
      } else {
        y_min = std::min(y[1], y[2]);
        y_mid = std::max(y[1], y[2]);
        y_max = y[0];
      }
    }

    t.length = abs(x_max - x_min) + abs(y_max - y_min);
    t.branch.resize(4);
    t.branch[0].x = x[0];
    t.branch[0].y = y[0];
    t.branch[0].n = 3;
    t.branch[1].x = x[1];
    t.branch[1].y = y[1];
    t.branch[1].n = 3;
    t.branch[2].x = x[2];
    t.branch[2].y = y[2];
    t.branch[2].n = 3;
    t.branch[3].x = x_mid;
    t.branch[3].y = y_mid;
    t.branch[3].n = 3;
  } else {
    std::vector<int> xs(d);
    std::vector<int> ys(d);

    std::vector<int> tmp_xs(d);
    std::vector<int> tmp_ys(d);
    std::vector<int> s(d);
    pnt* pt = new pnt[d];
    std::vector<pnt*> ptp(d);

    for (int i = 0; i < d; i++) {
      pt[i].x = x[i];
      pt[i].y = y[i];
      ptp[i] = &pt[i];
    }

    if (d < 1000) {
      for (int i = 0; i < d - 1; i++) {
        int minval = ptp[i]->x;
        int minidx = i;
        for (int j = i + 1; j < d; j++) {
          if (minval > ptp[j]->x) {
            minval = ptp[j]->x;
            minidx = j;
          }
        }
        std::swap(ptp[i], ptp[minidx]);
      }
    } else {
      std::stable_sort(ptp.begin(), ptp.end(), orderx);
    }

    for (int i = 0; i < d; i++) {
      xs[i] = ptp[i]->x;
      ptp[i]->o = i;
    }

    // sort y to find s[]
    if (d < 1000) {
      for (int i = 0; i < d - 1; i++) {
        int minval = ptp[i]->y;
        int minidx = i;
        for (int j = i + 1; j < d; j++) {
          if (minval > ptp[j]->y) {
            minval = ptp[j]->y;
            minidx = j;
          }
        }
        ys[i] = ptp[minidx]->y;
        s[i] = ptp[minidx]->o;
        ptp[minidx] = ptp[i];
      }
      ys[d - 1] = ptp[d - 1]->y;
      s[d - 1] = ptp[d - 1]->o;
    } else {
      std::stable_sort(ptp.begin(), ptp.end(), ordery);
      for (int i = 0; i < d; i++) {
        ys[i] = ptp[i]->y;
        s[i] = ptp[i]->o;
      }
    }

    gxs_[netID].resize(d);
    gys_[netID].resize(d);
    gs_[netID].resize(d);

    for (int i = 0; i < d; i++) {
      gxs_[netID][i] = xs[i];
      gys_[netID][i] = ys[i];
      gs_[netID][i] = s[i];

      tmp_xs[i] = xs[i] * 100;
      tmp_ys[i] = ys[i] * ((int) (100 * coeffV));
    }

    t = stt_builder_->makeSteinerTree(tmp_xs, tmp_ys, s, acc);

    for (auto& branch : t.branch) {
      branch.x /= 100;
      branch.y /= ((int) (100 * coeffV));
    }

    delete[] pt;
  }
}

void MorseRoute::fluteCongest(const int netID,
                                 const std::vector<int>& x,
                                 const std::vector<int>& y,
                                 const int acc,
                                 const float coeffV,
                                 Tree& t)
{
  const float coeffH = 1;
  const int d = x.size();

  if (d == 2) {
    t.deg = 2;
    t.length = abs(x[0] - x[1]) + abs(y[0] - y[1]);
    t.branch.resize(2);
    t.branch[0].x = x[0];
    t.branch[0].y = y[0];
    t.branch[0].n = 1;
    t.branch[1].x = x[1];
    t.branch[1].y = y[1];
    t.branch[1].n = 1;
  } else if (d == 3) {
    t.deg = 3;
    int x_max, x_min, x_mid;
    if (x[0] < x[1]) {
      if (x[0] < x[2]) {
        x_min = x[0];
        x_mid = std::min(x[1], x[2]);
        x_max = std::max(x[1], x[2]);
      } else {
        x_min = x[2];
        x_mid = x[0];
        x_max = x[1];
      }
    } else {
      if (x[0] < x[2]) {
        x_min = x[1];
        x_mid = x[0];
        x_max = x[2];
      } else {
        x_min = std::min(x[1], x[2]);
        x_mid = std::max(x[1], x[2]);
        x_max = x[0];
      }
    }
    int y_max, y_min, y_mid;
    if (y[0] < y[1]) {
      if (y[0] < y[2]) {
        y_min = y[0];
        y_mid = std::min(y[1], y[2]);
        y_max = std::max(y[1], y[2]);
      } else {
        y_min = y[2];
        y_mid = y[0];
        y_max = y[1];
      }
    } else {
      if (y[0] < y[2]) {
        y_min = y[1];
        y_mid = y[0];
        y_max = y[2];
      } else {
        y_min = std::min(y[1], y[2]);
        y_mid = std::max(y[1], y[2]);
        y_max = y[0];
      }
    }

    t.length = abs(x_max - x_min) + abs(y_max - y_min);
    t.branch.resize(4);
    t.branch[0].x = x[0];
    t.branch[0].y = y[0];
    t.branch[0].n = 3;
    t.branch[1].x = x[1];
    t.branch[1].y = y[1];
    t.branch[1].n = 3;
    t.branch[2].x = x[2];
    t.branch[2].y = y[2];
    t.branch[2].n = 3;
    t.branch[3].x = x_mid;
    t.branch[3].y = y_mid;
    t.branch[3].n = 3;
  } else {
    std::vector<int> xs(d);
    std::vector<int> ys(d);
    std::vector<int> nxs(d);
    std::vector<int> nys(d);
    std::vector<int> x_seg(d - 1);
    std::vector<int> y_seg(d - 1);
    std::vector<int> s(d);

    for (int i = 0; i < d; i++) {
      xs[i] = gxs_[netID][i];
      ys[i] = gys_[netID][i];
      s[i] = gs_[netID][i];
    }

    // get the new coordinates considering congestion
    for (int i = 0; i < d - 1; i++) {
      x_seg[i] = (xs[i + 1] - xs[i]) * 100;
      y_seg[i] = (ys[i + 1] - ys[i]) * 100;
    }

    const int height = ys[d - 1] - ys[0] + 1;  // # vertical grids the net span
    const int width = xs[d - 1] - xs[0] + 1;  // # horizontal grids the net span

    for (int i = 0; i < d - 1; i++) {
      int usageH = 0;
      for (int k = ys[0]; k <= ys[d - 1]; k++)  // all grids in the column
      {
        for (int j = xs[i]; j < xs[i + 1]; j++)
          usageH += h_edges_[k][j].est_usage_red();
      }
      if (x_seg[i] != 0 && usageH != 0) {
        x_seg[i]
            *= coeffH * usageH / ((xs[i + 1] - xs[i]) * height * h_capacity_);
        x_seg[i] = std::max(1, x_seg[i]);  // the segment len is at least 1 if
                                           // original segment len > 0
      }
      int usageV = 0;
      for (int j = ys[i]; j < ys[i + 1]; j++) {
        for (int k = xs[0]; k <= xs[d - 1]; k++)  // all grids in the row
          usageV += v_edges_[j][k].est_usage_red();
      }
      if (y_seg[i] != 0 && usageV != 0) {
        y_seg[i]
            *= coeffV * usageV / ((ys[i + 1] - ys[i]) * width * v_capacity_);
        y_seg[i] = std::max(1, y_seg[i]);  // the segment len is at least 1 if
                                           // original segment len > 0
      }
    }

    nxs[0] = xs[0];
    nys[0] = ys[0];
    for (int i = 0; i < d - 1; i++) {
      nxs[i + 1] = nxs[i] + x_seg[i];
      nys[i + 1] = nys[i] + y_seg[i];
    }

    t = stt_builder_->makeSteinerTree(nxs, nys, s, acc);

    // map the new coordinates back to original coordinates
    for (auto& branch : t.branch) {
      branch.x = mapxy(branch.x, xs, nxs, d);
      branch.y = mapxy(branch.y, ys, nys, d);
    }
  }
}

bool MorseRoute::netCongestion(const int netID)
{
  for (const MorseSegment& seg : seglist_[netID]) {
    const int ymin = std::min(seg.y1, seg.y2);
    const int ymax = std::max(seg.y1, seg.y2);

    // remove L routing
    if (seg.xFirst) {
      for (int i = seg.x1; i < seg.x2; i++) {
        const int cap = getEdgeCapacity(
            nets_[netID], i, seg.y1, MorseEdgeDirection::Horizontal);
        if (h_edges_[seg.y1][i].est_usage >= cap) {
          return true;
        }
      }
      for (int i = ymin; i < ymax; i++) {
        const int cap
            = getEdgeCapacity(nets_[netID], seg.x2, i, MorseEdgeDirection::Vertical);
        if (v_edges_[i][seg.x2].est_usage >= cap) {
          return true;
        }
      }
    } else {
      for (int i = ymin; i < ymax; i++) {
        const int cap
            = getEdgeCapacity(nets_[netID], seg.x1, i, MorseEdgeDirection::Vertical);
        if (v_edges_[i][seg.x1].est_usage >= cap) {
          return true;
        }
      }
      for (int i = seg.x1; i < seg.x2; i++) {
        const int cap = getEdgeCapacity(
            nets_[netID], i, seg.y2, MorseEdgeDirection::Horizontal);
        if (h_edges_[seg.y2][i].est_usage >= cap) {
          return true;
        }
      }
    }
  }
  return false;
}

void MorseRoute::gen_brk_RSMT(const bool congestionDriven,
                                 const bool reRoute,
                                 const bool genTree,
                                 const bool newType,
                                 const bool noADJ)
{
  Tree rsmt;
  int numShift = 0;

  int wl = 0;
  int wl1 = 0;
  int totalNumSeg = 0;

  const int flute_accuracy = 2;

  for (const int& netID : net_ids_) {
    MorseNet* net = nets_[netID];

    int d = net->getNumPins();


    // check net alpha because FastRoute has a special implementation of flute
    // TODO: move this flute implementation to SteinerTreeBuilder
    const float net_alpha = stt_builder_->getAlpha(net->getDbNet());
    if (net_alpha > 0.0) {
      rsmt = stt_builder_->makeSteinerTree(
          net->getDbNet(), net->getPinX(), net->getPinY(), net->getDriverIdx());
    } else {
      float coeffV = 1.36;

      if (congestionDriven) {
        // call congestion driven flute to generate RSMT
        bool cong;
        //coeffV = noADJ ? 1.2 : coeffADJ(netID);
        cong = netCongestion(netID);
        if (cong) {
          fluteCongest(netID,
                       net->getPinX(),
                       net->getPinY(),
                       flute_accuracy,
                       coeffV,
                       rsmt);
        } else {
          fluteNormal(netID,
                      net->getPinX(),
                      net->getPinY(),
                      flute_accuracy,
                      coeffV,
                      rsmt);
        }
        /*if (d > 3) {
          numShift += edgeShiftNew(rsmt, netID);
        }*/
      } else {
        // call FLUTE to generate RSMT for each net
        //if (noADJ || HTreeSuite(netID)) {
          coeffV = 1.2;
        //}
        fluteNormal(netID,
                    net->getPinX(),
                    net->getPinY(),
                    flute_accuracy,
                    coeffV,
                    rsmt);
      }
    }
    /*if (debug_->isOn() && debug_->steinerTree_
        && net->getDbNet() == debug_->net_) {
      steinerTreeVisualization(rsmt, net);
    }*/

    /*if (genTree) {
      copyStTree(netID, rsmt);
    }

    if (net->getNumPins() != rsmt.deg) {
      d = rsmt.deg;
    }

    if (congestionDriven) {
      for (int j = 0; j < sttrees_[netID].num_edges(); j++) {
        wl1 += sttrees_[netID].edges[j].len;
      }
    }*/

    for (int j = 0; j < rsmt.branchCount(); j++) {
      const int x1 = rsmt.branch[j].x;
      const int y1 = rsmt.branch[j].y;
      const int n = rsmt.branch[j].n;
      const int x2 = rsmt.branch[n].x;
      const int y2 = rsmt.branch[n].y;

      wl += abs(x1 - x2) + abs(y1 - y2);

      if (x1 != x2 || y1 != y2) {  // the branch is not degraded (a point)
        // the position of this segment in seglist
        seglist_[netID].push_back(Segment());
        auto& seg = seglist_[netID].back();
        if (x1 < x2) {
          seg.x1 = x1;
          seg.x2 = x2;
          seg.y1 = y1;
          seg.y2 = y2;
        } else {
          seg.x1 = x2;
          seg.x2 = x1;
          seg.y1 = y2;
          seg.y2 = y1;
        }

        seg.netID = netID;
      }
    }  // loop j

    totalNumSeg += seglist_[netID].size();

    /*if (reRoute) {
      // update the est_usage due to the segments in this net
      newrouteL(
          netID,
          RouteType::NoRoute,
          true);  // route the net with no previous route for each tree edge
    }*/
  }  // loop i

  debugPrint(logger_,
             GRT,
             "rsmt",
             1,
             "Wirelength: {}, Wirelength1: {}\nNumber of segments: {}\nNumber "
             "of shifts: {}",
             wl,
             wl1,
             totalNumSeg,
             numShift);
}

void MorseRoute::updateDbCongestion(int min_routing_layer,
                                       int max_routing_layer)
{
  if (h_edges_3D_.num_elements() == 0) {  // no information
    return;
  }
  auto block = db_->getChip()->getBlock();
  auto db_gcell = block->getGCellGrid();
  if (db_gcell)
    db_gcell->resetGrid();
  else
    db_gcell = odb::dbGCellGrid::create(block);

  db_gcell->addGridPatternX(x_corner_, x_grid_, tile_size_);
  db_gcell->addGridPatternY(y_corner_, y_grid_, tile_size_);
  auto db_tech = db_->getTech();
  for (int k = min_routing_layer - 1; k <= max_routing_layer - 1; k++) {
    auto layer = db_tech->findRoutingLayer(k + 1);
    if (layer == nullptr) {
      continue;
    }

    const uint8_t capH = h_capacity_3D_[k];
    const uint8_t capV = v_capacity_3D_[k];
    const uint8_t last_row_capH = last_row_h_capacity_3D_[k];
    const uint8_t last_col_capV = last_col_v_capacity_3D_[k];
    bool is_horizontal
        = layer_directions_[k] == odb::dbTechLayerDir::HORIZONTAL;
    for (int y = 0; y < y_grid_; y++) {
      for (int x = 0; x < x_grid_; x++) {
        if (is_horizontal) {
          if (!regular_y_ && y == y_grid_ - 1) {
            db_gcell->setCapacity(layer, x, y, last_row_capH);
          } else {
            db_gcell->setCapacity(layer, x, y, capH);
          }
        } else {
          if (!regular_x_ && x == x_grid_ - 1) {
            db_gcell->setCapacity(layer, x, y, last_col_capV);
          } else {
            db_gcell->setCapacity(layer, x, y, capV);
          }
        }
        if (x == x_grid_ - 1 && y == y_grid_ - 1 && x_grid_ > 1
            && y_grid_ > 1) {
          uint8_t blockageH = h_edges_3D_[k][y][x - 1].red;
          uint8_t blockageV = v_edges_3D_[k][y - 1][x].red;
          uint8_t usageH = h_edges_3D_[k][y - 1][x - 1].usage + blockageH;
          uint8_t usageV = v_edges_3D_[k][y - 1][x - 1].usage + blockageV;
          db_gcell->setUsage(layer, x, y, usageH + usageV);
        } else {
          uint8_t blockageH = h_edges_3D_[k][y][x].red;
          uint8_t blockageV = v_edges_3D_[k][y][x].red;
          uint8_t usageH = h_edges_3D_[k][y][x].usage + blockageH;
          uint8_t usageV = v_edges_3D_[k][y][x].usage + blockageV;
          db_gcell->setUsage(layer, x, y, usageH + usageV);
        }
      }
    }
  }
}

}  // namespace grt