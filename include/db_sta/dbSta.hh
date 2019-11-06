// StaDb, OpenSTA on OpenDB
// Copyright (c) 2019, Parallax Software, Inc.
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#ifndef STA_DB_OPENSTADB_H
#define STA_DB_OPENSTADB_H

#include "opendb/db.h"
#include "Sta.hh"

namespace sta {

class dbNetwork;

using odb::dbDatabase;

class dbSta : public Sta
{
public:
  dbSta(dbDatabase *db);
  dbDatabase *db() { return db_; }
  virtual void makeComponents();
  dbNetwork *dbNetwork();
  void readDbAfter();
  virtual LibertyLibrary *readLiberty(const char *filename,
				      Corner *corner,
				      const MinMaxAll *min_max,
				      bool infer_latches);

protected:
  virtual void makeNetwork();
  virtual void makeSdcNetwork();

  dbDatabase *db_;
};

} // namespace
#endif
