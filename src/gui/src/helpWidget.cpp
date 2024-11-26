/////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2024, The Regents of the University of California
// All rights reserved.
//
// BSD 3-Clause License
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

#include "helpWidget.h"

#include <dirent.h>
#include <sys/stat.h>

#include <QFile>
#include <QHBoxLayout>
#include <QListWidgetItem>
#include <QTextStream>
#include <QVBoxLayout>

namespace gui {

HelpWidget::HelpWidget(QWidget* parent)
    : QDockWidget("Help Browser", parent),
      category_selector_(new QComboBox(this)),
      help_list_(new QListWidget(this)),
      viewer_(new QTextBrowser(this))
{
  setObjectName("help_viewer");

  has_help_ = false;

  viewer_->setOpenExternalLinks(true);

  QHBoxLayout* layout = new QHBoxLayout;

  QVBoxLayout* select_layout = new QVBoxLayout;
  select_layout->addWidget(category_selector_);
  select_layout->addWidget(help_list_, /* stretch*/ 1);

  layout->addLayout(select_layout);
  layout->addWidget(viewer_, /* stretch*/ 1);

  QWidget* container = new QWidget(this);
  container->setLayout(layout);
  setWidget(container);

  // add connect
  connect(category_selector_,
          qOverload<int>(&QComboBox::currentIndexChanged),
          this,
          &HelpWidget::changeCategory);
  connect(help_list_,
          &QListWidget::currentItemChanged,
          this,
          &HelpWidget::showHelpInformation);
}

void HelpWidget::init(const std::string& path)
{
  category_selector_->clear();

  // validate path
  struct stat dir_info;
  if (stat(path.c_str(), &dir_info) == 0) {
    // check for html
    const std::string html_path = path + "/html";
    if (stat(html_path.c_str(), &dir_info) != 0) {
      return;
    }

    // add html1, html2, html3
    auto add_help = [this, &html_path](const std::string& category,
                                       const std::string dir) {
      struct stat dir_info;
      const std::string doc_path = html_path + "/" + dir;
      if (stat(doc_path.c_str(), &dir_info) == 0) {
        category_selector_->addItem(QString::fromStdString(category),
                                    QString::fromStdString(doc_path));
        has_help_ = true;
      }
    };

    add_help("application", "html1");
    add_help("commands", "html2");
    add_help("messages", "html3");
  }

  if (hasHelp()) {
    category_selector_->setCurrentIndex(1);
  }
}

void HelpWidget::changeCategory()
{
  help_list_->clear();

  const std::string path
      = category_selector_->currentData().toString().toStdString();

  DIR* dir = opendir(path.c_str());

  if (dir == nullptr) {
    return;
  }

  struct dirent* ent;
  while ((ent = readdir(dir)) != nullptr) {
    const std::string help_path = path + "/" + ent->d_name;

    if (help_path.find(".html") == std::string::npos) {
      continue;
    }

    const std::string name = std::string(ent->d_name, strlen(ent->d_name) - 5);

    auto* qitem = new QListWidgetItem(QString::fromStdString(name));
    qitem->setData(Qt::UserRole, QString::fromStdString(help_path));
    help_list_->addItem(qitem);
  }

  closedir(dir);

  help_list_->sortItems();
}

void HelpWidget::showHelpInformation(QListWidgetItem* item)
{
  if (item == nullptr) {
    return;
  }
  const QString path = item->data(Qt::UserRole).toString();
  QFile html_file(path);
  if (html_file.open(QIODevice::ReadOnly)) {
    QTextStream in(&html_file);
    viewer_->setHtml(in.readAll());
  } else {
    viewer_->setHtml("Unable to open: " + path);
  }
}

void HelpWidget::selectHelp(const std::string& item)
{
  const int start_index = category_selector_->currentIndex();

  for (int i = 0; i < category_selector_->count(); i++) {
    category_selector_->setCurrentIndex(i);

    const auto help_items
        = help_list_->findItems(QString::fromStdString(item), Qt::MatchExactly);
    if (!help_items.empty()) {
      raise();
      show();

      help_list_->setCurrentItem(help_items[0]);
      return;
    }
  }

  category_selector_->setCurrentIndex(start_index);
}

}  // namespace gui
