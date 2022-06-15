///////////////////////////////////////////////////////////////////////////////
// BSD 3-Clause License
//
// Copyright (c) 2021, Andrew Kennings
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

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
#include <map>
#include <vector>

#include "rectangle.h"

namespace dpo {
class Node;
}

namespace dpo {

////////////////////////////////////////////////////////////////////////////////
// Forward declarations.
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Classes.
////////////////////////////////////////////////////////////////////////////////
class RoutingParams
{
 public:
  class EdgeAdjust
  {
   public:
    int irow_ = -1;
    int icol_ = -1;
    int ilayer_ = -1;
    int jrow_ = -1;
    int jcol_ = -1;
    int jlayer_ = -1;
    double rcap_ = 0.0;

    EdgeAdjust() = default;

    EdgeAdjust(int irow,
               int icol,
               int ilayer,
               int jrow,
               int jcol,
               int jlayer,
               double rcap)
    {
      init(irow, icol, ilayer, jrow, jcol, jlayer, rcap);
    }
    EdgeAdjust(const EdgeAdjust& other)
    {
      init(other.irow_,
           other.icol_,
           other.ilayer_,
           other.jrow_,
           other.jcol_,
           other.jlayer_,
           other.rcap_);
    }
    EdgeAdjust& operator=(const EdgeAdjust& other)
    {
      if (this != &other) {
        init(other.irow_,
             other.icol_,
             other.ilayer_,
             other.jrow_,
             other.jcol_,
             other.jlayer_,
             other.rcap_);
      }
      return *this;
    }
    virtual ~EdgeAdjust() {}

    void init(int irow,
              int icol,
              int ilayer,
              int jrow,
              int jcol,
              int jlayer,
              double rcap)
    {
      irow_ = irow;
      icol_ = icol;
      ilayer_ = ilayer;
      jrow_ = jrow;
      jcol_ = jcol;
      jlayer_ = jlayer;
      rcap_ = rcap;
    }
  };

 public:
  RoutingParams()
      : grid_x_(0),
        grid_y_(0),
        num_layers_(0),
        default_layer_(1),
        origin_x_(0),
        origin_y_(0),
        tile_size_x_(0),
        tile_size_y_(0),
        blockage_porosity_(0.0),
        num_ni_terminals_(0),
        num_route_blockages_(0),
        num_edge_adjusts_(0),
        Xlowerbound_(0.0),
        Xupperbound_(0.0),
        Ylowerbound_(0.0),
        Yupperbound_(0.0),
        XpitchGcd_(0.0),
        YpitchGcd_(0.0),
        hasObs_(0),
        numRules_(0)
  {
    v_capacity_.erase(v_capacity_.begin(), v_capacity_.end());
    h_capacity_.erase(h_capacity_.begin(), h_capacity_.end());
    wire_width_.erase(wire_width_.begin(), wire_width_.end());
    wire_spacing_.erase(wire_spacing_.begin(), wire_spacing_.end());
    via_spacing_.erase(via_spacing_.begin(), via_spacing_.end());

    edge_adjusts_.erase(edge_adjusts_.begin(), edge_adjusts_.end());
  }

  void postProcess();

  // Get spacing between two objects.
  double get_spacing(int layer,
                     double xmin1,
                     double xmax1,
                     double ymin1,
                     double ymax1,
                     double xmin2,
                     double xmax2,
                     double ymin2,
                     double ymax2);
  double get_spacing(int layer, double width, double parallel);
  double get_maximum_spacing(int layer);

 public:
  int grid_x_;
  int grid_y_;
  int num_layers_;

  int default_layer_;

  double origin_x_;
  double origin_y_;

  std::vector<double> v_capacity_;
  std::vector<double> h_capacity_;
  std::vector<double> wire_width_;
  std::vector<double> wire_spacing_;
  std::vector<double> via_spacing_;
  std::vector<int> layer_dir_;

  double tile_size_x_;
  double tile_size_y_;

  double blockage_porosity_;

  int num_ni_terminals_;
  int num_route_blockages_;

  // Stuff for edge adjustements (ICCAD12).
  int num_edge_adjusts_;
  std::vector<EdgeAdjust> edge_adjusts_;

  // Map for routing blockages...  We have the name of the node and a vector of
  // layers with which it interferes...
  std::map<Node*, std::vector<unsigned>*> blockage_;

  // Other blockages which are simply specified by rectangles on a layer...
  std::vector<std::vector<Rectangle>> layerBlockages_;

  // Added to get information from LEF/DEF...
  double Xlowerbound_;
  double Xupperbound_;
  double Ylowerbound_;
  double Yupperbound_;
  double XpitchGcd_;
  double YpitchGcd_;
  int hasObs_;
  std::vector<std::vector<std::vector<unsigned>>> obs_;

  // Stuff for routing rules...  These vectors should all be the same length...
  int numRules_;
  std::vector<std::vector<double>> ruleWidths_;
  std::vector<std::vector<double>> ruleSpacings_;

  // Stuff for spacing tables...  Only one spacing table per layer...
  std::vector<std::vector<double>> spacingTableWidth_;
  std::vector<std::vector<double>> spacingTableLength_;
  std::vector<std::vector<std::vector<double>>> spacingTable_;
};

}  // namespace dpo
