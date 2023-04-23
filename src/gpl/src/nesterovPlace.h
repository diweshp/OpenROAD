///////////////////////////////////////////////////////////////////////////////
// BSD 3-Clause License
//
// Copyright (c) 2018-2020, The Regents of the University of California
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
// ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
///////////////////////////////////////////////////////////////////////////////

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "point.h"
#include "nesterovBase.h"

namespace utl {
class Logger;
}

namespace odb {
class dbInst;
}

namespace gpl {

class PlacerBase;
class PlacerBaseCommon;
class Instance;
class RouteBase;
class TimingBase;
class Graphics;

class NesterovPlace
{
 public:
  NesterovPlace();
  NesterovPlace(const NesterovPlaceVars& npVars,
                std::shared_ptr<PlacerBaseCommon> pbc,
                std::shared_ptr<NesterovBaseCommon> nbc,
                std::vector<std::shared_ptr<PlacerBase>>& pbVec,
                std::vector<std::shared_ptr<NesterovBase>>& nbVec,
                std::shared_ptr<RouteBase> rb,
                std::shared_ptr<TimingBase> tb,
                utl::Logger* log);
  ~NesterovPlace();

  // return iteration count
  int doNesterovPlace(int start_iter = 0);


  void updateWireLengthCoef(float overflow);


  void updateNextIter(const int iter);

  void updateDb();

  float getWireLengthCoefX() const { return wireLengthCoefX_; }
  float getWireLengthCoefY() const { return wireLengthCoefY_; }
  float getDensityPenalty() const { return densityPenalty_; }

  void setTargetOverflow(float overflow) { npVars_.targetOverflow = overflow; }
  void setMaxIters(int limit) { npVars_.maxNesterovIter = limit; }

 private:
  std::shared_ptr<PlacerBaseCommon> pbc_;
  std::shared_ptr<NesterovBaseCommon> nbc_;
  std::vector<std::shared_ptr<PlacerBase>> pbVec_;
  std::vector<std::shared_ptr<NesterovBase>> nbVec_;
  utl::Logger* log_;
  std::shared_ptr<RouteBase> rb_;
  std::shared_ptr<TimingBase> tb_;
  NesterovPlaceVars npVars_;
  std::unique_ptr<Graphics> graphics_;

  // // SLP is Step Length Prediction.
  // //
  // // y_st, y_dst, y_wdst, w_pdst
  // std::vector<FloatPoint> curSLPCoordi_;
  // std::vector<FloatPoint> curSLPWireLengthGrads_;
  // std::vector<FloatPoint> curSLPDensityGrads_;
  // std::vector<FloatPoint> curSLPSumGrads_;

  // // y0_st, y0_dst, y0_wdst, y0_pdst
  // std::vector<FloatPoint> nextSLPCoordi_;
  // std::vector<FloatPoint> nextSLPWireLengthGrads_;
  // std::vector<FloatPoint> nextSLPDensityGrads_;
  // std::vector<FloatPoint> nextSLPSumGrads_;

  // // z_st, z_dst, z_wdst, z_pdst
  // std::vector<FloatPoint> prevSLPCoordi_;
  // std::vector<FloatPoint> prevSLPWireLengthGrads_;
  // std::vector<FloatPoint> prevSLPDensityGrads_;
  // std::vector<FloatPoint> prevSLPSumGrads_;

  // // x_st and x0_st
  // std::vector<FloatPoint> curCoordi_;
  // std::vector<FloatPoint> nextCoordi_;

  // // save initial coordinates -- needed for RD
  // std::vector<FloatPoint> initCoordi_;

  // densityPenalty stor
  std::vector<float> densityPenaltyStor_;

  float wireLengthGradSum_;
  float densityGradSum_;


  // opt_phi_cof
  float densityPenalty_;

  // base_wcof
  float baseWireLengthCoef_;

  // wlen_cof
  float wireLengthCoefX_;
  float wireLengthCoefY_;

  // phi is described in ePlace paper.
  float sumOverflow_;
  float sumOverflowUnscaled_;

  // half-parameter-wire-length
  int64_t prevHpwl_;

  float isDiverged_;
  float isRoutabilityNeed_;

  std::string divergeMsg_;
  int divergeCode_;

  int recursionCntWlCoef_;
  int recursionCntInitSLPCoef_;

  void cutFillerCoordinates();

  void init();
  void reset();
};
}  // namespace gpl
