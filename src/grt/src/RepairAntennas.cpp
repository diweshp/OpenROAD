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

#include "RepairAntennas.h"

#include <omp.h>

#include <algorithm>
#include <limits>
#include <map>
#include <string>
#include <unordered_set>
#include <vector>

#include "Net.h"
#include "Pin.h"
#include "grt/GlobalRouter.h"
#include "utl/Logger.h"

#include <boost/pending/disjoint_sets.hpp>

namespace grt {

using utl::GRT;

RepairAntennas::RepairAntennas(GlobalRouter* grouter,
                               ant::AntennaChecker* arc,
                               dpl::Opendp* opendp,
                               odb::dbDatabase* db,
                               utl::Logger* logger)
    : grouter_(grouter),
      arc_(arc),
      opendp_(opendp),
      db_(db),
      logger_(logger),
      unique_diode_index_(1),
      illegal_diode_placement_count_(0)
{
  block_ = db_->getChip()->getBlock();
  while (block_->findInst(
      fmt::format("ANTENNA_{}", unique_diode_index_).c_str())) {
    unique_diode_index_++;
  }
}

bool RepairAntennas::checkAntennaViolations(
    NetRouteMap& routing,
    const std::vector<odb::dbNet*>& nets_to_repair,
    int max_routing_layer,
    odb::dbMTerm* diode_mterm,
    float ratio_margin,
    const int num_threads)
{
  for (odb::dbNet* db_net : nets_to_repair) {
    antenna_violations_[db_net];
  }

  bool destroy_wires = !grouter_->haveDetailedRoutes();

  makeNetWires(routing, nets_to_repair, max_routing_layer);
  arc_->initAntennaRules();
  omp_set_num_threads(num_threads);
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < nets_to_repair.size(); i++) {
    odb::dbNet* db_net = nets_to_repair[i];
    checkNetViolations(db_net, diode_mterm, ratio_margin);
  }

  if (destroy_wires) {
    destroyNetWires(nets_to_repair);
  }

  // remove nets with zero violations
  for (auto it = antenna_violations_.begin();
       it != antenna_violations_.end();) {
    if (it->second.empty()) {
      it = antenna_violations_.erase(it);
    } else {
      ++it;
    }
  }

  logger_->info(
      GRT, 12, "Found {} antenna violations.", antenna_violations_.size());
  return !antenna_violations_.empty();
}

void RepairAntennas::checkNetViolations(odb::dbNet* db_net,
                                        odb::dbMTerm* diode_mterm,
                                        float ratio_margin)
{
  if (!db_net->isSpecial() && db_net->getWire()) {
    std::vector<ant::Violation> net_violations
        = arc_->getAntennaViolations(db_net, diode_mterm, ratio_margin);
    if (!net_violations.empty()) {
      antenna_violations_[db_net] = std::move(net_violations);
      debugPrint(logger_,
                 GRT,
                 "repair_antennas",
                 1,
                 "antenna violations {}",
                 db_net->getConstName());
    }
  }
}

void RepairAntennas::makeNetWires(
    NetRouteMap& routing,
    const std::vector<odb::dbNet*>& nets_to_repair,
    int max_routing_layer)
{
  std::map<int, odb::dbTechVia*> default_vias
      = grouter_->getDefaultVias(max_routing_layer);

  for (odb::dbNet* db_net : nets_to_repair) {
    if (!db_net->isSpecial() && !db_net->isConnectedByAbutment()
        && !grouter_->getNet(db_net)->isLocal()
        && !grouter_->isDetailedRouted(db_net)) {
      makeNetWire(db_net, routing[db_net], default_vias);
    }
  }
}

odb::dbWire* RepairAntennas::makeNetWire(
    odb::dbNet* db_net,
    GRoute& route,
    std::map<int, odb::dbTechVia*>& default_vias)
{
  odb::dbWire* wire = odb::dbWire::create(db_net);
  if (wire) {
    Net* net = grouter_->getNet(db_net);
    odb::dbTech* tech = db_->getTech();
    odb::dbWireEncoder wire_encoder;
    wire_encoder.begin(wire);
    RoutePtPinsMap route_pt_pins = findRoutePtPins(net);
    std::unordered_set<GSegment, GSegmentHash> wire_segments;
    int prev_conn_layer = -1;
    for (GSegment& seg : route) {
      int l1 = seg.init_layer;
      int l2 = seg.final_layer;
      auto [bottom_layer, top_layer] = std::minmax(l1, l2);

      odb::dbTechLayer* bottom_tech_layer
          = tech->findRoutingLayer(bottom_layer);
      odb::dbTechLayer* top_tech_layer = tech->findRoutingLayer(top_layer);

      if (std::abs(seg.init_layer - seg.final_layer) > 1) {
        debugPrint(logger_,
                   GRT,
                   "check_antennas",
                   1,
                   "invalid seg: ({}, {})um to ({}, {})um",
                   block_->dbuToMicrons(seg.init_x),
                   block_->dbuToMicrons(seg.init_y),
                   block_->dbuToMicrons(seg.final_x),
                   block_->dbuToMicrons(seg.final_y));

        logger_->error(GRT,
                       68,
                       "Global route segment for net {} not "
                       "valid. The layers {} and {} "
                       "are not adjacent.",
                       net->getName(),
                       bottom_tech_layer->getName(),
                       top_tech_layer->getName());
      }
      if (wire_segments.find(seg) == wire_segments.end()) {
        int x1 = seg.init_x;
        int y1 = seg.init_y;

        if (seg.isVia()) {
          if (bottom_layer >= grouter_->getMinRoutingLayer()) {
            if (bottom_layer == prev_conn_layer) {
              wire_encoder.newPath(bottom_tech_layer, odb::dbWireType::ROUTED);
              prev_conn_layer = std::max(l1, l2);
            } else if (top_layer == prev_conn_layer) {
              wire_encoder.newPath(top_tech_layer, odb::dbWireType::ROUTED);
              prev_conn_layer = std::min(l1, l2);
            } else {
              // if a via is the first object added to the wire_encoder, or the
              // via starts a new path and is not connected to previous wires
              // create a new path using the bottom layer and do not update the
              // prev_conn_layer. this way, this process is repeated until the
              // first wire is added and properly update the prev_conn_layer
              wire_encoder.newPath(bottom_tech_layer, odb::dbWireType::ROUTED);
            }
            wire_encoder.addPoint(x1, y1);
            wire_encoder.addTechVia(default_vias[bottom_layer]);
            addWireTerms(net,
                         route,
                         x1,
                         y1,
                         bottom_layer,
                         bottom_tech_layer,
                         route_pt_pins,
                         wire_encoder,
                         default_vias,
                         false);
            wire_segments.insert(seg);
          }
        } else {
          // Add wire
          int x2 = seg.final_x;
          int y2 = seg.final_y;
          if (x1 != x2 || y1 != y2) {
            odb::dbTechLayer* tech_layer = tech->findRoutingLayer(l1);
            addWireTerms(net,
                         route,
                         x1,
                         y1,
                         l1,
                         tech_layer,
                         route_pt_pins,
                         wire_encoder,
                         default_vias,
                         true);
            wire_encoder.newPath(tech_layer, odb::dbWireType::ROUTED);
            wire_encoder.addPoint(x1, y1);
            wire_encoder.addPoint(x2, y2);
            addWireTerms(net,
                         route,
                         x2,
                         y2,
                         l1,
                         tech_layer,
                         route_pt_pins,
                         wire_encoder,
                         default_vias,
                         true);
            wire_segments.insert(seg);
            prev_conn_layer = l1;
          }
        }
      }
    }
    wire_encoder.end();

    return wire;
  } else {
    logger_->error(
        GRT, 221, "Cannot create wire for net {}.", db_net->getConstName());
    // suppress gcc warning
    return nullptr;
  }
}

RoutePtPinsMap RepairAntennas::findRoutePtPins(Net* net)
{
  RoutePtPinsMap route_pt_pins;
  for (Pin& pin : net->getPins()) {
    int conn_layer = pin.getConnectionLayer();
    odb::Point grid_pt = pin.getOnGridPosition();
    RoutePt route_pt(grid_pt.x(), grid_pt.y(), conn_layer);
    route_pt_pins[route_pt].pins.push_back(&pin);
    route_pt_pins[route_pt].connected = false;
  }
  return route_pt_pins;
}

void RepairAntennas::addWireTerms(Net* net,
                                  GRoute& route,
                                  int grid_x,
                                  int grid_y,
                                  int layer,
                                  odb::dbTechLayer* tech_layer,
                                  RoutePtPinsMap& route_pt_pins,
                                  odb::dbWireEncoder& wire_encoder,
                                  std::map<int, odb::dbTechVia*>& default_vias,
                                  bool connect_to_segment)
{
  std::vector<int> layers;
  layers.push_back(layer);
  if (layer == grouter_->getMinRoutingLayer()) {
    layer--;
    layers.push_back(layer);
  }

  for (int l : layers) {
    auto itr = route_pt_pins.find(RoutePt(grid_x, grid_y, l));
    if (itr != route_pt_pins.end() && !itr->second.connected) {
      for (const Pin* pin : itr->second.pins) {
        itr->second.connected = true;
        int conn_layer = pin->getConnectionLayer();
        std::vector<odb::Rect> pin_boxes = pin->getBoxes().at(conn_layer);
        odb::Point grid_pt = pin->getOnGridPosition();
        odb::Point pin_pt = grid_pt;
        // create the local connection with the pin center only when the global
        // segment doesn't overlap the pin
        if (!pinOverlapsGSegment(grid_pt, conn_layer, pin_boxes, route)) {
          int min_dist = std::numeric_limits<int>::max();
          for (const odb::Rect& pin_box : pin_boxes) {
            odb::Point pos = grouter_->getRectMiddle(pin_box);
            int dist = odb::Point::manhattanDistance(pos, pin_pt);
            if (dist < min_dist) {
              min_dist = dist;
              pin_pt = pos;
            }
          }
        }

        if (conn_layer >= grouter_->getMinRoutingLayer()) {
          wire_encoder.newPath(tech_layer, odb::dbWireType::ROUTED);
          wire_encoder.addPoint(grid_pt.x(), grid_pt.y());
          wire_encoder.addPoint(pin_pt.x(), grid_pt.y());
          wire_encoder.addPoint(pin_pt.x(), pin_pt.y());
        } else {
          odb::dbTech* tech = db_->getTech();
          odb::dbTechLayer* min_layer
              = tech->findRoutingLayer(grouter_->getMinRoutingLayer());

          if (connect_to_segment && tech_layer != min_layer) {
            // create vias to connect the guide segment to the min routing
            // layer. the min routing layer will be used to connect to the pin.
            wire_encoder.newPath(tech_layer, odb::dbWireType::ROUTED);
            wire_encoder.addPoint(grid_pt.x(), grid_pt.y());
            for (int l = min_layer->getRoutingLevel();
                 l < tech_layer->getRoutingLevel();
                 l++) {
              wire_encoder.addTechVia(default_vias[l]);
            }
          }

          if (min_layer->getDirection() == odb::dbTechLayerDir::VERTICAL) {
            makeWire(wire_encoder,
                     min_layer,
                     grid_pt,
                     odb::Point(grid_pt.x(), pin_pt.y()));
            wire_encoder.addTechVia(
                default_vias[grouter_->getMinRoutingLayer()]);
            makeWire(wire_encoder,
                     min_layer,
                     odb::Point(grid_pt.x(), pin_pt.y()),
                     pin_pt);
          } else {
            makeWire(wire_encoder,
                     min_layer,
                     grid_pt,
                     odb::Point(pin_pt.x(), grid_pt.y()));
            wire_encoder.addTechVia(
                default_vias[grouter_->getMinRoutingLayer()]);
            makeWire(wire_encoder,
                     min_layer,
                     odb::Point(pin_pt.x(), grid_pt.y()),
                     pin_pt);
          }

          // create vias to reach the pin
          for (int i = min_layer->getRoutingLevel() - 1; i >= conn_layer; i--) {
            wire_encoder.addTechVia(default_vias[i]);
          }
        }
      }
    }
  }
}

void RepairAntennas::makeWire(odb::dbWireEncoder& wire_encoder,
                              odb::dbTechLayer* layer,
                              const odb::Point& start,
                              const odb::Point& end)
{
  wire_encoder.newPath(layer, odb::dbWireType::ROUTED);
  wire_encoder.addPoint(start.x(), start.y());
  wire_encoder.addPoint(end.x(), end.y());
}

bool RepairAntennas::pinOverlapsGSegment(
    const odb::Point& pin_position,
    const int pin_layer,
    const std::vector<odb::Rect>& pin_boxes,
    const GRoute& route)
{
  // check if pin position on grid overlaps with the pin shape
  for (const odb::Rect& box : pin_boxes) {
    if (box.overlaps(pin_position)) {
      return true;
    }
  }

  // check if pin position on grid overlaps with at least one GSegment
  for (const odb::Rect& box : pin_boxes) {
    for (const GSegment& seg : route) {
      if (seg.init_layer == seg.final_layer &&  // ignore vias
          seg.init_layer == pin_layer) {
        int x0 = std::min(seg.init_x, seg.final_x);
        int y0 = std::min(seg.init_y, seg.final_y);
        int x1 = std::max(seg.init_x, seg.final_x);
        int y1 = std::max(seg.init_y, seg.final_y);
        odb::Rect seg_rect(x0, y0, x1, y1);

        if (box.intersects(seg_rect)) {
          return true;
        }
      }
    }
  }

  return false;
}

void RepairAntennas::destroyNetWires(
    const std::vector<odb::dbNet*>& nets_to_repair)
{
  for (odb::dbNet* db_net : nets_to_repair) {
    odb::dbWire* wire = db_net->getWire();
    if (wire) {
      odb::dbWire::destroy(wire);
    }
  }
}

void RepairAntennas::repairAntennas(odb::dbMTerm* diode_mterm)
{
  int site_width = -1;
  r_tree fixed_insts;
  odb::dbTech* tech = db_->getTech();

  illegal_diode_placement_count_ = 0;
  diode_insts_.clear();

  auto rows = block_->getRows();
  for (odb::dbRow* db_row : rows) {
    odb::dbSite* site = db_row->getSite();
    if (site->getClass() == odb::dbSiteClass::PAD) {
      continue;
    }
    int width = site->getWidth();
    if (site_width == -1) {
      site_width = width;
    }

    if (site_width != width) {
      logger_->warn(GRT, 27, "Design has rows with different site widths.");
    }
  }

  setInstsPlacementStatus(odb::dbPlacementStatus::FIRM);
  getFixedInstances(fixed_insts);

  bool repair_failures = false;
  for (auto const& net_violations : antenna_violations_) {
    odb::dbNet* db_net = net_violations.first;
    auto violations = net_violations.second;

    bool inserted_diodes = false;
    for (ant::Violation& violation : violations) {
      debugPrint(logger_,
                 GRT,
                 "repair_antennas",
                 2,
                 "antenna {} insert {} diodes",
                 db_net->getConstName(),
                 violation.diode_count_per_gate * violation.gates.size());
      if (violation.diode_count_per_gate > 0) {
        for (odb::dbITerm* gate : violation.gates) {
          for (int j = 0; j < violation.diode_count_per_gate; j++) {
            odb::dbTechLayer* violation_layer
                = tech->findRoutingLayer(violation.routing_level);
            insertDiode(db_net,
                        diode_mterm,
                        gate,
                        site_width,
                        fixed_insts,
                        violation_layer);
            inserted_diodes = true;
          }
        }
      } else
        repair_failures = true;
    }
    if (inserted_diodes)
      grouter_->addDirtyNet(db_net);
  }
  if (repair_failures)
    logger_->warn(GRT, 243, "Unable to repair antennas on net with diodes.");
}

void RepairAntennas::legalizePlacedCells()
{
  opendp_->detailedPlacement(0, 0, "");
  // After legalize placement, diodes and violated insts don't need to be FIRM
  setInstsPlacementStatus(odb::dbPlacementStatus::PLACED);
}

void RepairAntennas::insertDiode(odb::dbNet* net,
                                 odb::dbMTerm* diode_mterm,
                                 odb::dbITerm* gate,
                                 int site_width,
                                 r_tree& fixed_insts,
                                 odb::dbTechLayer* violation_layer)
{
  odb::dbMaster* diode_master = diode_mterm->getMaster();
  std::string diode_inst_name
      = "ANTENNA_" + std::to_string(unique_diode_index_++);
  odb::dbInst* diode_inst
      = odb::dbInst::create(block_, diode_master, diode_inst_name.c_str());

  bool place_vertically
      = violation_layer->getDirection() == odb::dbTechLayerDir::VERTICAL;
  bool legally_placed = setDiodeLoc(
      diode_inst, gate, site_width, place_vertically, fixed_insts);

  odb::Rect inst_rect = diode_inst->getBBox()->getBox();

  legally_placed = legally_placed && diodeInRow(inst_rect);

  if (!legally_placed)
    illegal_diode_placement_count_++;

  // allow detailed placement to move diodes with geometry out of the core area,
  // or near macro pins (can be placed out of row), or illegal placed diodes
  const odb::Rect& core_area = block_->getCoreArea();
  odb::dbInst* sink_inst = gate->getInst();
  if (core_area.contains(inst_rect) && !sink_inst->getMaster()->isBlock()
      && legally_placed) {
    diode_inst->setPlacementStatus(odb::dbPlacementStatus::FIRM);
  } else {
    diode_inst->setPlacementStatus(odb::dbPlacementStatus::PLACED);
  }

  odb::dbITerm* diode_iterm
      = diode_inst->findITerm(diode_mterm->getConstName());
  diode_iterm->connect(net);
  diode_insts_.push_back(diode_inst);

  // Add diode to the R-tree of fixed instances
  int fixed_inst_id = fixed_insts.size();
  box b(point(inst_rect.xMin(), inst_rect.yMin()),
        point(inst_rect.xMax(), inst_rect.yMax()));
  value v(b, fixed_inst_id);
  fixed_insts.insert(v);
  fixed_inst_id++;
}

void RepairAntennas::getFixedInstances(r_tree& fixed_insts)
{
  int fixed_inst_id = 0;
  for (odb::dbInst* inst : block_->getInsts()) {
    odb::dbPlacementStatus status = inst->getPlacementStatus();
    if (status == odb::dbPlacementStatus::FIRM
        || status == odb::dbPlacementStatus::LOCKED) {
      odb::dbBox* instBox = inst->getBBox();
      box b(point(instBox->xMin(), instBox->yMin()),
            point(instBox->xMax(), instBox->yMax()));
      value v(b, fixed_inst_id);
      fixed_insts.insert(v);
      fixed_inst_id++;
    }
  }
}

void RepairAntennas::setInstsPlacementStatus(
    odb::dbPlacementStatus placement_status)
{
  for (auto const& violation : antenna_violations_) {
    for (int i = 0; i < violation.second.size(); i++) {
      for (odb::dbITerm* gate : violation.second[i].gates) {
        if (!gate->getMTerm()->getMaster()->isBlock()) {
          gate->getInst()->setPlacementStatus(placement_status);
        }
      }
    }
  }

  for (odb::dbInst* diode_inst : diode_insts_) {
    diode_inst->setPlacementStatus(placement_status);
  }
}

bool RepairAntennas::setDiodeLoc(odb::dbInst* diode_inst,
                                 odb::dbITerm* gate,
                                 int site_width,
                                 const bool place_vertically,
                                 r_tree& fixed_insts)
{
  const int max_legalize_itr = 50;
  bool place_at_left = true;
  bool place_at_top = false;
  int left_offset = 0, right_offset = 0;
  int top_offset = 0, bottom_offset = 0;
  int horizontal_offset = 0, vertical_offset = 0;
  bool legally_placed = false;

  int inst_loc_x, inst_loc_y, inst_width, inst_height;
  odb::dbOrientType inst_orient;
  getInstancePlacementData(
      gate, inst_loc_x, inst_loc_y, inst_width, inst_height, inst_orient);

  odb::dbBox* diode_bbox = diode_inst->getBBox();
  int diode_width = diode_bbox->xMax() - diode_bbox->xMin();
  int diode_height = diode_bbox->yMax() - diode_bbox->yMin();
  odb::dbInst* sink_inst = gate->getInst();

  // Use R-tree to check if diode will not overlap or cause 1-site spacing with
  // other fixed cells
  int legalize_itr = 0;
  while (!legally_placed && legalize_itr < max_legalize_itr) {
    if (place_vertically) {
      computeVerticalOffset(inst_height,
                            top_offset,
                            bottom_offset,
                            place_at_top,
                            vertical_offset);
    } else {
      computeHorizontalOffset(diode_width,
                              inst_width,
                              site_width,
                              left_offset,
                              right_offset,
                              place_at_left,
                              horizontal_offset);
    }
    diode_inst->setOrient(inst_orient);
    if (sink_inst->isBlock() || sink_inst->isPad() || place_vertically) {
      int x_center = inst_loc_x + horizontal_offset + diode_width / 2;
      int y_center = inst_loc_y + vertical_offset + diode_height / 2;
      odb::Point diode_center(x_center, y_center);
      odb::dbOrientType orient = getRowOrient(diode_center);
      diode_inst->setOrient(orient);
    }
    diode_inst->setLocation(inst_loc_x + horizontal_offset,
                            inst_loc_y + vertical_offset);

    legally_placed = checkDiodeLoc(diode_inst, site_width, fixed_insts);
    legalize_itr++;
  }

  return legally_placed;
}

void RepairAntennas::getInstancePlacementData(odb::dbITerm* gate,
                                              int& inst_loc_x,
                                              int& inst_loc_y,
                                              int& inst_width,
                                              int& inst_height,
                                              odb::dbOrientType& inst_orient)
{
  odb::dbInst* sink_inst = gate->getInst();
  odb::Rect sink_bbox = getInstRect(sink_inst, gate);
  inst_loc_x = sink_bbox.xMin();
  inst_loc_y = sink_bbox.yMin();
  inst_width = sink_bbox.xMax() - sink_bbox.xMin();
  inst_height = sink_bbox.yMax() - sink_bbox.yMin();
  inst_orient = sink_inst->getOrient();
}

bool RepairAntennas::checkDiodeLoc(odb::dbInst* diode_inst,
                                   const int site_width,
                                   r_tree& fixed_insts)
{
  const odb::Rect& core_area = block_->getCoreArea();
  const int left_pad = opendp_->padLeft(diode_inst);
  const int right_pad = opendp_->padRight(diode_inst);
  odb::dbBox* instBox = diode_inst->getBBox();
  box box(point(instBox->xMin() - ((left_pad + right_pad) * site_width) + 1,
                instBox->yMin() + 1),
          point(instBox->xMax() + ((left_pad + right_pad) * site_width) - 1,
                instBox->yMax() - 1));

  std::vector<value> overlap_insts;
  fixed_insts.query(bgi::intersects(box), std::back_inserter(overlap_insts));

  return overlap_insts.empty() && core_area.contains(instBox->getBox());
}

void RepairAntennas::computeHorizontalOffset(const int diode_width,
                                             const int inst_width,
                                             const int site_width,
                                             int& left_offset,
                                             int& right_offset,
                                             bool& place_at_left,
                                             int& offset)
{
  if (place_at_left) {
    offset = -(diode_width + left_offset * site_width);
    left_offset++;
    place_at_left = false;
  } else {
    offset = inst_width + right_offset * site_width;
    right_offset++;
    place_at_left = true;
  }
}

void RepairAntennas::computeVerticalOffset(const int inst_height,
                                           int& top_offset,
                                           int& bottom_offset,
                                           bool& place_at_top,
                                           int& offset)
{
  if (place_at_top) {
    offset = top_offset * inst_height;
    top_offset++;
    place_at_top = false;
  } else {
    offset = -(bottom_offset * inst_height);
    bottom_offset++;
    place_at_top = true;
  }
}

odb::Rect RepairAntennas::getInstRect(odb::dbInst* inst, odb::dbITerm* iterm)
{
  const odb::dbTransform transform = inst->getTransform();

  odb::Rect inst_rect;

  if (inst->getMaster()->isBlock()) {
    inst_rect.mergeInit();
    odb::dbMTerm* mterm = iterm->getMTerm();
    if (mterm != nullptr) {
      for (odb::dbMPin* mterm_pin : mterm->getMPins()) {
        for (odb::dbBox* mterm_box : mterm_pin->getGeometry()) {
          odb::Rect rect = mterm_box->getBox();
          transform.apply(rect);

          inst_rect = rect;
        }
      }
    }
  } else {
    inst_rect = inst->getBBox()->getBox();
  }

  return inst_rect;
}

bool RepairAntennas::diodeInRow(odb::Rect diode_rect)
{
  int diode_height = diode_rect.dy();
  for (odb::dbRow* row : block_->getRows()) {
    odb::Rect row_rect = row->getBBox();
    int row_height = row_rect.dy();

    if (row_rect.contains(diode_rect) && diode_height == row_height) {
      return true;
    }
  }

  return false;
}

odb::dbOrientType RepairAntennas::getRowOrient(const odb::Point& point)
{
  odb::dbOrientType orient;
  for (odb::dbRow* row : block_->getRows()) {
    odb::Rect row_rect = row->getBBox();

    if (row_rect.overlaps(point)) {
      orient = row->getOrient();
    }
  }

  return orient;
}

odb::dbMTerm* RepairAntennas::findDiodeMTerm()
{
  for (odb::dbLib* lib : db_->getLibs()) {
    for (auto master : lib->getMasters()) {
      if (master->getType() == odb::dbMasterType::CORE_ANTENNACELL) {
        for (odb::dbMTerm* mterm : master->getMTerms()) {
          if (diffArea(mterm) > 0.0)
            return mterm;
        }
      }
    }
  }
  return nullptr;
}

// copied from AntennaChecker
double RepairAntennas::diffArea(odb::dbMTerm* mterm)
{
  double max_diff_area = 0.0;
  std::vector<std::pair<double, odb::dbTechLayer*>> diff_areas;
  mterm->getDiffArea(diff_areas);
  for (auto [diff_area, layer] : diff_areas) {
    max_diff_area = std::max(max_diff_area, diff_area);
  }
  return max_diff_area;
}

//////////////////////////////////

void addSegments(GRoute& route,const int& init_x,const int& init_y,const int& final_x,const int& final_y,const int& layer_level){
  // create vias
  route.push_back(GSegment(init_x, init_y, layer_level, init_x, init_y, layer_level + 1));
  route.push_back(GSegment(init_x, init_y, layer_level + 1, init_x, init_y, layer_level + 2));
  route.push_back(GSegment(final_x, final_y, layer_level, final_x, final_y, layer_level + 1));
  route.push_back(GSegment(final_x, final_y, layer_level + 1, final_x, final_y, layer_level + 2));
  // divide segment (new segment is created with size = required_size)
  route.push_back(GSegment(init_x, init_y, layer_level + 2, final_x, final_y, layer_level + 2));
}

int getSegmentPos(std::vector<GSegment*>& segments, int& req_size, int bridge_size, bool isHorizontal, bool inStart)
{
  if (isHorizontal) {
    sort(segments.begin(),
         segments.end(),
         [](const GSegment* seg1, const GSegment* seg2) {
           return seg1->init_x < seg2->init_x;
         });
  } else {
    sort(segments.begin(),
         segments.end(),
         [](const GSegment* seg1, const GSegment* seg2) {
           return seg1->init_y < seg2->init_y;
	 });
  }
  int size_acum = 0;
  int cand = -1;
  if (inStart) {
    for (int pos = 0; pos < segments.size(); pos++) {
      size_acum += segments[pos]->length();
      if (size_acum > req_size && size_acum >= req_size + bridge_size) {
        req_size -= (size_acum - segments[pos]->length());
        cand = pos;
	break;
      }
    }
  } else {
    for (int pos = segments.size()-1; pos >= 0; pos--) {
      size_acum += segments[pos]->length();
      if (size_acum > req_size && size_acum >= req_size + bridge_size) {
        req_size -= (size_acum - segments[pos]->length());
        cand = pos;
	break;
      }
    }
  }
  return cand;
}

int divideSegment(std::vector<GSegment*>& segments, GRoute& route, odb::dbTechLayer* violation_layer, const int& tile_size, const double& ratio, const int& init_c, const int& final_c, const double& init_area, const double& final_area)
{
  int length = 0;
  bool is_horizontal = violation_layer->getDirection() == odb::dbTechLayerDir::HORIZONTAL;
  for (const auto& it : segments) {
    length += it->length();
  }
  const int n_tiles = length/tile_size;
  int req_tiles = int(n_tiles/ratio) * 0.8;
  if (init_c != 0 && final_c != 0) {
    req_tiles = std::max(2, int(req_tiles * 0.15));//(init_area / (init_area+final_area));
  }
  int req_size = req_tiles * tile_size;
  const int bridge_size = 2 * tile_size;
  const int layer_level = violation_layer->getRoutingLevel();
  const int bridges_number = std::ceil(double(n_tiles - req_tiles) / double(req_tiles + 2));
  const int m_size = std::ceil(double(n_tiles - (2 * bridges_number)) / (bridges_number + 1)); //tiles
  std::cerr << "tiles of segment: " << n_tiles <<  " tiles required: " << req_tiles << " size required: " << req_tiles * tile_size << "\n";
  std::cerr << "jumpers need: " << bridges_number << " size of subseg: " << m_size << "\n";

  int jumper_count = 0;
  // place bridge in segment begin
  if (final_c == 0 || (init_c != 0 && final_c != 0 && init_area > 0.0)) {
    jumper_count++;
    int pos = getSegmentPos(segments, req_size, bridge_size, is_horizontal, true);
    if (pos == -1) std::cerr << "Not found\n";
    auto & seg = segments[pos];
    if (is_horizontal) {
      // get jumper position
      const int bridge_init_x = seg->init_x + req_size;
      const int bridge_final_x = bridge_init_x + bridge_size;
      // create vias
      std::cerr << "Last size: " << route.size() << " ";
      addSegments(route, bridge_init_x, seg->init_y, bridge_final_x, seg->final_y, layer_level);
      route.push_back(GSegment(seg->init_x, seg->init_y, layer_level, bridge_init_x, seg->init_y, layer_level));
      seg->init_x = bridge_final_x;
      std::cerr << "New size: " << route.size() << std::endl;
    } else {
      // Get jumper position
      const int bridge_init_y = seg->init_y + req_size;
      const int bridge_final_y = bridge_init_y + bridge_size;
      // create vias
      std::cerr << "Last size: " << route.size() << " ";
      addSegments(route, seg->init_x, bridge_init_y, seg->final_x, bridge_final_y, layer_level);
      route.push_back(GSegment(seg->init_x, seg->init_y, layer_level, seg->init_x, bridge_init_y, layer_level));
      seg->init_y = bridge_final_y;
      std::cerr << "New size: " << route.size() << std::endl;
    }
  }
  // if need place other in segment end
  if (init_c == 0 || (final_c != 0 && init_c != 0 && final_area > 0.0)) {
    jumper_count++;
    int pos = getSegmentPos(segments, req_size, bridge_size, is_horizontal, false);
    if (pos == -1) std::cerr << "Not found\n";
    auto & seg = segments[pos];
    if (is_horizontal) {
      // Get jumper position
      const int bridge_init_x = seg->final_x - req_size - bridge_size;
      const int bridge_final_x = bridge_init_x + bridge_size;
      // create vias
      std::cerr << "Last size: " << route.size() << " ";
      addSegments(route, bridge_init_x, seg->init_y, bridge_final_x, seg->final_y, layer_level);
      route.push_back(GSegment(seg->init_x, seg->init_y, layer_level, bridge_init_x, seg->init_y, layer_level));
      seg->init_x = bridge_final_x;
      std::cerr << "New size: " << route.size() << std::endl;
    } else {
      // Get jumper position
      const int bridge_init_y = seg->final_y - req_size - bridge_size;
      const int bridge_final_y = bridge_init_y + bridge_size;
      // create vias
      std::cerr << "Last size: " << route.size() << " ";
      addSegments(route, seg->init_x, bridge_init_y, seg->final_x, bridge_final_y, layer_level);
      route.push_back(GSegment(seg->init_x, seg->init_y, layer_level, seg->init_x, bridge_init_y, layer_level));
      seg->init_y = bridge_final_y;
      std::cerr << "New size: " << route.size() << std::endl;
    }
  }
  return jumper_count;
}

int manhattanDistance(int x1, int y1, int x2, int y2)
{
  return std::abs(x1-x2) + std::abs(y1-y2);
}

double gateArea(odb::dbMTerm* mterm)
{       
  double max_gate_area = 0;
  if (mterm->hasDefaultAntennaModel()) {
    odb::dbTechAntennaPinModel* pin_model = mterm->getDefaultAntennaModel();
    std::vector<std::pair<double, odb::dbTechLayer*>> gate_areas;
    pin_model->getGateArea(gate_areas);
  
    for (const auto& [gate_area, layer] : gate_areas) {
      max_gate_area = std::max(max_gate_area, gate_area);
    }
  }
  return max_gate_area;
}

void RepairAntennas::getPinNumberNearEndPoint(std::vector<GSegment*>& segments, const std::vector<odb::dbITerm*>& gates, int& init_c, int& final_c, double& init_area, double& final_area)
{
  int seg_init_x, seg_init_y, seg_final_x, seg_final_y;
  seg_init_x = seg_init_y = std::numeric_limits<int>::max();;
  seg_final_x = seg_final_y = 0;
  for (const auto& seg : segments) {
    seg_init_x = std::min(seg_init_x, seg->init_x);
    seg_init_y = std::min(seg_init_y, seg->init_y);
    seg_final_x = std::max(seg_final_x, seg->final_x);
    seg_final_y = std::max(seg_final_y, seg->final_y);
  }
  for (const auto& iterm : gates) {
    odb::dbInst* sink_inst = iterm->getInst();
    odb::Rect sink_bbox = getInstRect(sink_inst, iterm);
    const int x_min = sink_bbox.xMin();
    const int y_min = sink_bbox.yMin();
    const int x_max = sink_bbox.xMax();
    const int y_max = sink_bbox.yMax();
    const int dist1 = std::min(std::min(manhattanDistance(x_min, y_min, seg_init_x, seg_init_y),
		                        manhattanDistance(x_min, y_max, seg_init_x, seg_init_y)),
			       std::min(manhattanDistance(x_max, y_min, seg_init_x, seg_init_y),
				        manhattanDistance(x_max, y_max, seg_init_x, seg_init_y)));
    const int dist2 = std::min(std::min(manhattanDistance(x_min, y_min, seg_final_x, seg_final_y),
		                        manhattanDistance(x_min, y_max, seg_final_x, seg_final_y)),
			       std::min(manhattanDistance(x_max, y_min, seg_final_x, seg_final_y),
				        manhattanDistance(x_max, y_max, seg_final_x, seg_final_y)));
    //std::cerr << x_min << " " << y_min << " " << x_max << " " << y_max << "\n";
    if (dist1 < dist2) {
      init_c++;
      init_area += gateArea(iterm->getMTerm());
    }
    else {
      final_c++;
      final_area += gateArea(iterm->getMTerm());
    }
  }
  std::cerr << "Init area " << init_area << " final area " << final_area << "\n";
  std::cerr << "Init count " << init_c << " final count " << final_c << "\n";
}

struct SegInfo
{
  int id = -1;
  GSegment* seg;
  odb::Rect rect;
  std::vector<std::pair<odb::dbTechLayer*,int>> low_adj;
  std::unordered_set<std::string> gates;
  SegInfo(int id_, GSegment* seg_, odb::Rect rect_)
  {
    id = id_;
    seg = seg_;
    rect = rect_;
  }
};

void getSegmentsWithOverlap(SegInfo& seg_info, const std::vector<SegInfo>& low_segs, odb::dbTechLayer* layer)
{
  int index = 0;
  for (const SegInfo& seg_it : low_segs) {
    if (seg_info.rect.overlaps(seg_it.rect)) {
      seg_info.low_adj.push_back({layer, index});
    }
    index++;
  }
}

SegmentByViolation RepairAntennas::getSegmentsWithViolation(odb::dbNet* db_net, GRoute& route, int& max_layer, std::map<int, int>& layer_with_violation)
{
  odb::dbTech* tech = db_->getTech();
  int min_layer = 1;
  // get segments information
  std::unordered_map<odb::dbTechLayer*, std::vector<SegInfo>> segment_by_layer;
  int seg_count = 0;
  for (GSegment& seg: route) {
    if (std::max(seg.final_layer, seg.init_layer) <= max_layer) {
      odb::dbTechLayer* tech_layer = tech->findRoutingLayer(std::min(seg.init_layer, seg.final_layer));
      if (seg.isVia()) {
        segment_by_layer[tech_layer->getUpperLayer()].push_back(SegInfo(seg_count, &seg, grouter_->globalRoutingToBox(seg)));
        // add one segment in upper and lower layer (to connect stacked vias)
	segment_by_layer[tech->findRoutingLayer(seg.init_layer)].push_back(SegInfo(seg_count++, nullptr, grouter_->globalRoutingToBox(seg)));
        segment_by_layer[tech->findRoutingLayer(seg.final_layer)].push_back(SegInfo(seg_count++, nullptr, grouter_->globalRoutingToBox(seg)));
      } else {
        segment_by_layer[tech_layer].push_back(SegInfo(seg_count, &seg, grouter_->globalRoutingToBox(seg)));
      }
      seg_count++;
    }
  }

  // Set adjacents segments (neighbor)
  for (auto& layer_it : segment_by_layer) {
     odb::dbTechLayer* tech_layer = layer_it.first;
     // find segments in same layer that touch segment
     for (SegInfo& seg_info : layer_it.second) {
       getSegmentsWithOverlap(seg_info, segment_by_layer[tech_layer], tech_layer);
     }
     // find segments in lower layer
     odb::dbTechLayer* lower_layer = tech_layer->getLowerLayer();
     if (lower_layer && segment_by_layer.find(lower_layer) != segment_by_layer.end()) {
       for (SegInfo& seg_info : layer_it.second) {
         getSegmentsWithOverlap(seg_info, segment_by_layer[lower_layer], lower_layer);
       }
     }
     //std::cerr << "layer " << layer_it.first->getConstName() << " " << layer_it.second.size() << " segments\n";
  }

  // iterate all instance pins
  std::unordered_map<std::string, std::unordered_set<int>> seg_connected;
  for (odb::dbITerm* iterm : db_net->getITerms()) {
    odb::dbMTerm* mterm = iterm->getMTerm();
    std::string pin_name = fmt::format("  {}/{} ({})",
                                       iterm->getInst()->getConstName(),
                                       mterm->getConstName(),
                                       mterm->getMaster()->getConstName());
    odb::dbInst* inst = iterm->getInst();
    const odb::dbTransform transform = inst->getTransform();
    for (odb::dbMPin* mterm : mterm->getMPins()) {
      for (odb::dbBox* box : mterm->getGeometry()) {
        odb::dbTechLayer* tech_layer = box->getTechLayer();
        if (tech_layer->getType() != odb::dbTechLayerType::ROUTING) {
          continue; 
        } 
        // get lower and upper layer (avoid vias)
        odb::dbTechLayer* upper_layer = tech_layer->getUpperLayer();
        odb::dbTechLayer* lower_layer = tech_layer->getLowerLayer(); 
    
        odb::Rect pin_rect = box->getBox();
        transform.apply(pin_rect);

	for (const SegInfo& it: segment_by_layer[tech_layer]) {
          if (it.rect.overlaps(pin_rect)) {
            seg_connected[pin_name].insert(it.id);
	   }
        }
        // optional
        if (upper_layer) {
	  for (const SegInfo& it: segment_by_layer[upper_layer]) {
            if (it.rect.overlaps(pin_rect)) {
              seg_connected[pin_name].insert(it.id);
	    }
	  }
	}

	if (lower_layer) {
          for (const SegInfo& it: segment_by_layer[lower_layer]) {
            if (it.rect.overlaps(pin_rect)) {
              seg_connected[pin_name].insert(it.id);
	    }
	  }
	}
      }
    }
  }

  auto violations = antenna_violations_[db_net];
  SegmentByViolation segment_with_violations(violations.size());

  // Run DSU
  std::vector<int> dsu_parent(seg_count);
  std::vector<int> dsu_size(seg_count);
  for (int i = 0; i < seg_count; i++) {
    dsu_size[i] = 1;
    dsu_parent[i] = i;
  }

  boost::disjoint_sets<int*, int*> dsu(&dsu_size[0], &dsu_parent[0]);

  odb::dbTechLayer* lower_layer;
  odb::dbTechLayer* layer_iter = tech->findRoutingLayer(min_layer);
  for (; layer_iter; layer_iter = layer_iter->getUpperLayer() ) {
    // iterate each node of this layer to union set
    for (auto& seg_it : segment_by_layer[layer_iter]) {
      int id_u = seg_it.id;
      // get neighbors and union
      for (const auto& lower_it : seg_it.low_adj) {
        int id_v = segment_by_layer[lower_it.first][lower_it.second].id;
        // if they are on different sets then union
        if (dsu.find_set(id_u) != dsu.find_set(id_v)) {
          dsu.union_set(id_u, id_v);
        }
      }
    }
    int iter = layer_iter->getRoutingLevel();
    // verify if layer has violation and save segments with same set
    if (layer_with_violation.find(iter) != layer_with_violation.end()) {
      for (auto& seg_it : segment_by_layer[layer_iter]) {
	int id_u = seg_it.id;
	//std::cerr << "Seg in layer " << layer_iter->getConstName() << "\n";
        for (const auto& iterm : violations[layer_with_violation[iter]].gates) {
          odb::dbMTerm* mterm = iterm->getMTerm();
          std::string pin_name = fmt::format("  {}/{} ({})",
                                       iterm->getInst()->getConstName(),
                                       mterm->getConstName(),
                                       mterm->getMaster()->getConstName());
	  //std::cerr << pin_name << " has " << seg_connected[pin_name].size() << " segments connected\n";
          bool is_conected = false;
	  for (const int& nbr_id : seg_connected[pin_name]) {
            if (dsu.find_set(id_u) == dsu.find_set(nbr_id)) {
              //std::cerr << "Seg found by layer " << layer_iter->getConstName() << "\n";
	      is_conected = true;
              break;
	    }
	  }
	  if (is_conected) {
            if (seg_it.seg) {
              segment_with_violations[layer_with_violation[iter]].push_back(seg_it.seg);
              std::cerr << "Seg with violation found in layer " << layer_iter->getConstName() << "\n";
	    }
            break;
	  }
	}
      }
    } 
  }
  return segment_with_violations;
}

void RepairAntennas::jumperInsertion(NetRouteMap& routing, const int tile_size)
{
  odb::dbTech* tech = db_->getTech();
  int total_jumper = 0;
  for (auto const& net_violations : antenna_violations_) {
    odb::dbNet* db_net = net_violations.first;
    const auto violations = net_violations.second;
    std::map<int, int> routing_layer_with_violations;
    std::cerr << "Net " << db_net->getConstName() << "\n";

    int max_layer = 1;
    // Iterate found layers with violations
    int violation_id = 0;
    for (const ant::Violation& violation : violations) {
      odb::dbTechLayer* violation_layer
                = tech->findRoutingLayer(violation.routing_level);
      odb::dbTechLayer* upper_layer = tech->findRoutingLayer(violation.routing_level + 2);
      if (upper_layer && violation_layer->getType() == odb::dbTechLayerType::ROUTING) { // avoid other layers 
        routing_layer_with_violations[violation.routing_level] = violation_id;
        max_layer = std::max(max_layer, violation.routing_level);	
        std::cerr << "Layer " << violation_layer->getConstName() << " ratio: " << violation.ratio << " n_pins: " << violation.gates.size() << "\n";
      }
      violation_id++;
    }

    // Iterate route of net (get segments)
    GRoute& route = routing[db_net];

    SegmentByViolation segment_with_violations;
    if (routing_layer_with_violations.size()) {
      segment_with_violations = getSegmentsWithViolation(db_net, route, max_layer, routing_layer_with_violations);
    }
    int init_c, final_c;
    double init_area, final_area;
    // add jumpers in segments for each layer
    for (auto it : routing_layer_with_violations) {
      std::cerr << "Violation in layer " << tech->findRoutingLayer(it.first)->getConstName() << ", segments found: " << segment_with_violations[it.second].size() << "\n";
      init_c = 0;
      final_c = 0;
      init_area = 0.0;
      final_area = 0.0;
      getPinNumberNearEndPoint(segment_with_violations[it.second], violations[it.second].gates, init_c, final_c, init_area, final_area);
      total_jumper += divideSegment(segment_with_violations[it.second], route, tech->findRoutingLayer(it.first), tile_size, violations[it.second].ratio, init_c, final_c, init_area, final_area);
    }
  }
  std::cerr << "Total jumper inserted " << total_jumper << "\n";
}

}  // namespace grt
