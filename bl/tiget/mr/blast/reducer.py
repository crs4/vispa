# BEGIN_COPYRIGHT
# 
# Copyright (C) 2013-2014 CRS4.
# 
# This file is part of vispa.
# 
# vispa is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
# 
# vispa is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
# 
# You should have received a copy of the GNU General Public License along with
# vispa.  If not, see <http://www.gnu.org/licenses/>.
# 
# END_COPYRIGHT
import pydoop.pipes as pp


class Reducer(pp.Reducer):
  
  def __init__(self, ctx):
    super(reducer, self).__init__(ctx)

  def reduce(self, ctx):
    pass
