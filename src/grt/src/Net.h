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

#pragma once

#include <string>
#include <utility>
#include <vector>

#include "Pin.h"
#include "grt/GRoute.h"
#include "odb/db.h"

namespace grt {

class Net
{
 public:
  Net(odb::dbNet* net, bool has_wires);
  odb::dbNet* getDbNet() const { return net_; }
  const std::string getName() const;
  const char* getConstName() const;
  odb::dbSigType getSignalType() const;
  void addPin(Pin& pin);
  void deleteSegment(int seg_id, GRoute& routes);
  std::vector<Pin>& getPins() { return pins_; }
  int getNumPins() const { return pins_.size(); }
  float getSlack() const { return slack_; }
  void setSlack(float slack) { slack_ = slack; }
  void setHasWires(bool in) { has_wires_ = in; }
  void setSegmentParent(std::vector<uint16_t> segment_parent)
  {
    segment_parent_ = std::move(segment_parent);
  }
  std::vector<uint16_t> getSegmentParent() const { return segment_parent_; }
  std::vector<std::vector<uint16_t>> getSegmentGraph();
  bool isLocal();
  void destroyPins();
  bool hasWires() const { return has_wires_; }
  bool hasStackedVias(odb::dbTechLayer* max_routing_layer);
  odb::Rect computeBBox();

 private:
  int getNumBTermsAboveMaxLayer(odb::dbTechLayer* max_routing_layer);

  odb::dbNet* net_;
  std::vector<Pin> pins_;
  float slack_;
  bool has_wires_;
  std::vector<uint16_t> segment_parent_;
};

}  // namespace grt
