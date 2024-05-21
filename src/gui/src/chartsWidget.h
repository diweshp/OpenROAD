///////////////////////////////////////////////////////////////////////////////
// BSD 3-Clause License
//
// Copyright (c) 2021, The Regents of the University of California
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

#include <QDockWidget>
#include <QLabel>

#ifdef ENABLE_CHARTS
#include <QComboBox>
#include <QPushButton>
#include <QString>
#include <QtCharts>
#include <memory>

#include "gui/gui.h"
#include "staGuiInterface.h"

namespace sta {
class Pin;
class dbSta;
class Clock;
}  // namespace sta
#endif

namespace gui {
#ifdef ENABLE_CHARTS

using ITermBTermPinsLists = std::pair<StaPins, StaPins>;

enum StartEndPathType
{
  RegisterToRegister,
  RegisterToIO,
  IOToRegister,
  IOToIO,
};

struct SlackHistogramData
{
  StaPins constrained_pins;
  std::set<sta::Clock*> clocks;
};

struct Buckets
{
  std::deque<std::vector<const sta::Pin*>> positive;
  std::deque<std::vector<const sta::Pin*>> negative;
};

class HistogramView : public QChartView
{
  Q_OBJECT

 public:
  HistogramView(QChart* chart, QWidget* parent);

  virtual void mousePressEvent(QMouseEvent* event) override;

 signals:
  void barIndex(int bar_index);
};

#endif

class ChartsWidget : public QDockWidget
{
  Q_OBJECT

 public:
  ChartsWidget(QWidget* parent = nullptr);
#ifdef ENABLE_CHARTS
  void setSTA(sta::dbSta* sta);
  void setLogger(utl::Logger* logger)
  {
    logger_ = logger;
  }

 signals:
  void endPointsToReport(const std::set<const sta::Pin*>& report_pins);

 private slots:
  void changeMode();
  void changeStartEndFilter();
  void showToolTip(bool is_hovering, int bar_index);
  void emitEndPointsInBucket(int bar_index);

 private:
  enum Mode
  {
    SELECT,
    SLACK_HISTOGRAM
  };

  static std::string toString(enum StartEndPathType);

  void setSlackHistogram();
  void setModeMenu();
  void setStartEndFiltersMenu();
  void setBucketInterval();
  void setBucketInterval(float bucket_interval)
  {
    bucket_interval_ = bucket_interval;
  }
  void setNegativeCountOffset(int neg_count_offset)
  {
    neg_count_offset_ = neg_count_offset;
  }
  void setDecimalPrecision(int precision_count)
  {
    precision_count_ = precision_count;
  }

  SlackHistogramData fetchSlackHistogramData();
  void fetchConstrainedPins(StaPins& end_points, bool set_mix_max);
  TimingPathList fetchPathsBasedOnStartEnd(const StaPins& end_points,
                                           const StartEndPathType path_type);
  StaPins getEndPointsFromPaths(const TimingPathList& paths);
  ITermBTermPinsLists separatePinsIntoBTermsAndITerms(const StaPins& pins);

  void populateBuckets(const StaPins& end_points);
  void populateBarSets(QBarSet& neg_set, QBarSet& pos_set);

  int computeSnapBucketInterval(float exact_interval);
  float computeSnapBucketDecimalInterval(float minimum_interval);
  int computeNumberofBuckets(int bucket_interval,
                             float max_slack,
                             float min_slack);

  void setXAxisConfig(int all_bars_count, const std::set<sta::Clock*>& clocks);
  void setYAxisConfig();
  QString createXAxisTitle(const std::set<sta::Clock*>& clocks);
  int computeMaxYSnap();
  int computeNumberOfDigits(int value);
  int computeFirstDigit(int value, int digits);
  int computeYInterval();
  void clearBarSets();

  void clearChart();

  utl::Logger* logger_;
  sta::dbSta* sta_;
  std::unique_ptr<STAGuiInterface> stagui_;

  QComboBox* mode_menu_;
  QComboBox* filters_menu_;
  QChart* chart_;
  HistogramView* display_;
  QValueAxis* axis_x_;
  QValueAxis* axis_y_;

  std::unique_ptr<Buckets> buckets_;
  int prev_filter_index_;

  float max_slack_;
  float min_slack_;

  const int default_number_of_buckets_ = 15;
  int largest_slack_count_ = 0;  // Used to configure the y axis.
  int precision_count_ = 0;      // Used to configure the x labels.

  float bucket_interval_ = 0;
  int neg_count_offset_ = 0;
#endif
  QLabel* label_;
};

}  // namespace gui
