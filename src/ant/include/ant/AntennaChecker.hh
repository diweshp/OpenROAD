// BSD 3-Clause License
//
// Copyright (c) 2020, MICL, DD-Lab, University of Michigan
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

#pragma once

#include <map>
#include <unordered_set>
#include <queue>

//#include <boost/polygon/polygon.hpp>

#include "odb/db.h"
#include "odb/dbWireGraph.h"
#include "utl/Logger.h"

namespace grt {
class GlobalRouter;
}

namespace ant {

using std::vector;

using odb::dbInst;
using odb::dbITerm;
using odb::dbMTerm;
using odb::dbNet;
using odb::dbSWire;
using odb::dbTechLayer;
using odb::dbWire;
using odb::dbWireGraph;

using utl::Logger;

struct PARinfo;
struct ARinfo;
struct AntennaModel;

///////////////////////////////////////
//namespace gtl = boost::polygon;
struct PinType;
struct GraphNode;

struct InfoType{
  double PAR;
  double PSR;
  double diff_PAR;
  double diff_PSR;
  double area;
  double side_area;
  double iterm_gate_area;
  double iterm_diff_area;

  double CAR;
  double CSR;
  double diff_CAR;
  double diff_CSR;

  std::vector<dbITerm*> iterms;

  InfoType& operator +=(const InfoType& a)
  {
    PAR += a.PAR;
    PSR += a.PSR;
    diff_PAR += a.diff_PAR;
    diff_PSR += a.diff_PSR;
    area += a.area;
    side_area += a.side_area;
    return *this;
  }
  InfoType () {
    PAR = 0.0;
    PSR = 0.0;
    diff_PAR = 0.0;
    diff_PSR = 0.0;

    area = 0.0;
    side_area = 0.0;
    iterm_gate_area = 0.0;
    iterm_diff_area = 0.0;

    CAR = 0.0;
    CSR = 0.0;
    diff_CAR = 0.0;
    diff_CSR = 0.0; 
  }
};

typedef std::unordered_map<dbTechLayer*, InfoType> LayerInfoVector;
typedef std::vector<GraphNode*> GraphNodeVector;
///////////////////////////////////////

class GlobalRouteSource
{
 public:
  virtual ~GlobalRouteSource() = default;

  virtual bool haveRoutes() = 0;
  virtual void makeNetWires() = 0;
  virtual void destroyNetWires() = 0;
};

struct Violation
{
  int routing_level;
  vector<odb::dbITerm*> gates;
  int diode_count_per_gate;
};

class AntennaChecker
{
 public:
  AntennaChecker();
  ~AntennaChecker();

  void init(odb::dbDatabase* db,
            GlobalRouteSource* global_route_source,
            utl::Logger* logger);

  // net nullptr -> check all nets
  int checkAntennas(dbNet* net = nullptr, bool verbose = false);
  int antennaViolationCount() const;

  void findMaxWireLength();

  vector<Violation> getAntennaViolations2(dbNet* net,
                                         odb::dbMTerm* diode_mterm,
                                         float ratio_margin);
  vector<Violation> getAntennaViolations(dbNet* net,
                                         odb::dbMTerm* diode_mterm,
                                         float ratio_margin);
  void initAntennaRules();
  void setReportFileName(const char* file_name);

 private:
  bool haveRoutedNets();
  double dbuToMicrons(int value);

  dbWireGraph::Node* findSegmentRoot(dbWireGraph::Node* node_info,
                                     int wire_level);
  dbWireGraph::Node* findSegmentStart(dbWireGraph::Node* node);
  bool ifSegmentRoot(dbWireGraph::Node* node, int wire_level);

  void findWireBelowIterms(dbWireGraph::Node* node,
                           double& iterm_gate_area,
                           double& iterm_diff_area,
                           int wire_level,
                           std::set<dbITerm*>& iv,
                           std::set<dbWireGraph::Node*>& nv);
  std::pair<double, double> calculateWireArea(
      dbWireGraph::Node* node,
      int wire_level,
      std::set<dbWireGraph::Node*>& nv,
      std::set<dbWireGraph::Node*>& level_nodes);

  double getViaArea(dbWireGraph::Edge* edge);
  dbTechLayer* getViaLayer(dbWireGraph::Edge* edge);
  std::string getViaName(dbWireGraph::Edge* edge);
  double calculateViaArea(dbWireGraph::Node* node, int wire_level);
  dbWireGraph::Edge* findVia(dbWireGraph::Node* node, int wire_level);

  void findCarPath(dbWireGraph::Node* node,
                   int wire_level,
                   dbWireGraph::Node* goal,
                   vector<dbWireGraph::Node*>& current_path,
                   vector<dbWireGraph::Node*>& path_found);

  void calculateParInfo(PARinfo& par_info);
  double getPwlFactor(odb::dbTechLayerAntennaRule::pwl_pair pwl_info,
                      double ref_val,
                      double def);

  vector<PARinfo> buildWireParTable(
      const vector<dbWireGraph::Node*>& wire_roots);
  vector<ARinfo> buildWireCarTable(
      const vector<PARinfo>& PARtable,
      const vector<PARinfo>& VIA_PARtable,
      const vector<dbWireGraph::Node*>& gate_iterms);
  vector<PARinfo> buildViaParTable(
      const vector<dbWireGraph::Node*>& wire_roots);
  vector<ARinfo> buildViaCarTable(
      const vector<PARinfo>& PARtable,
      const vector<PARinfo>& VIA_PARtable,
      const vector<dbWireGraph::Node*>& gate_iterms);

  vector<dbWireGraph::Node*> findWireRoots(dbWire* wire);
  void findWireRoots(dbWire* wire,
                     // Return values.
                     vector<dbWireGraph::Node*>& wire_roots,
                     vector<dbWireGraph::Node*>& gate_iterms);

  std::pair<bool, bool> checkWirePar(const ARinfo& AntennaRatio,
                                     bool report,
                                     bool verbose,
                                     std::ofstream& report_file);
  std::pair<bool, bool> checkWireCar(const ARinfo& AntennaRatio,
                                     bool par_checked,
                                     bool report,
                                     bool verbose,
                                     std::ofstream& report_file);
  bool checkViaPar(const ARinfo& AntennaRatio,
                   bool report,
                   bool verbose,
                   std::ofstream& report_file);
  bool checkViaCar(const ARinfo& AntennaRatio,
                   bool report,
                   bool verbose,
                   std::ofstream& report_file);

  void checkNet(dbNet* net,
                bool report_if_no_violation,
                bool verbose,
                std::ofstream& report_file,
                // Return values.
                int& net_violation_count,
                int& pin_violation_count);
  void checkGate(dbWireGraph::Node* gate,
                 vector<ARinfo>& CARtable,
                 vector<ARinfo>& VIA_CARtable,
                 bool report,
                 bool verbose,
                 std::ofstream& report_file,
                 // Return values.
                 bool& violation,
                 std::unordered_set<dbWireGraph::Node*>& violated_gates);
  bool checkViolation(const PARinfo& par_info, dbTechLayer* layer);
  bool antennaRatioDiffDependent(dbTechLayer* layer);

  void findWireRootIterms(dbWireGraph::Node* node,
                          int wire_level,
                          vector<dbITerm*>& gates);
  double diffArea(dbMTerm* mterm);
  double gateArea(dbMTerm* mterm);

  vector<std::pair<double, vector<dbITerm*>>> parMaxWireLength(dbNet* net,
                                                               int layer);
  vector<std::pair<double, vector<dbITerm*>>> getViolatedWireLength(
      dbNet* net,
      int routing_level);

  odb::dbDatabase* db_{nullptr};
  odb::dbBlock* block_{nullptr};
  int dbu_per_micron_{0};
  GlobalRouteSource* global_route_source_{nullptr};
  utl::Logger* logger_{nullptr};
  std::map<odb::dbTechLayer*, AntennaModel> layer_info_;
  int net_violation_count_{0};
  float ratio_margin_{0};
  std::string report_file_name_;

  static constexpr int max_diode_count_per_gate = 10;

  /////////////////////////////////////////////////////////////////////////////////
  std::unordered_map<odb::dbTechLayer*, GraphNodeVector> node_by_layer_map_;
  std::unordered_map<std::string, LayerInfoVector> info_;
  std::vector<Violation> antenna_violations_;
  int node_count_;
  dbTechLayer *min_layer_;
  // dsu variables
  std::vector<int> dsu_parent_, dsu_size_;

  // DSU functions
  void init_dsu();
  int find_set(int v);
  void union_set(int u, int v);
  bool dsu_same(int u, int v);

  bool isValidGate(dbMTerm* mterm);
  void buildLayerMaps(dbNet* net);
  void checkNet(dbNet* net, bool verbose, bool report, std::ofstream& report_file, dbMTerm* diode_mterm, float ratio_margin, int& net_violation_count, int& pin_violation_count);
  void saveGates(dbNet* db_net);
  void calculateAreas();
  void calculatePAR();
  void calculateCAR();
  int checkInfo(dbNet* db_net, bool verbose, bool report, std::ofstream& report_file, dbMTerm* diode_mterm, float ratio_margin);
  void calculateViaPar(dbTechLayer* tech_layer, InfoType& info);
  void calculateWirePar(dbTechLayer* tech_layer, InfoType& info);
  std::pair<bool, bool> checkPAR(dbTechLayer* tech_layer, const InfoType& info, bool verbose, bool report, std::ofstream& report_file);
  std::pair<bool, bool> checkPSR(dbTechLayer* tech_layer, const InfoType& info, bool verbose, bool report, std::ofstream& report_file);
  bool checkCAR(dbTechLayer* tech_layer, const InfoType& info, bool verbose, bool report, std::ofstream& report_file);
  bool checkCSR(dbTechLayer* tech_layer, const InfoType& info, bool verbose, bool report, std::ofstream& report_file); 
  /////////////////////////////////////////////////////////////////////////////////
};

}  // namespace ant
