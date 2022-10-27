#
# Copyright (c) 2018-2020 Helge Eichhorn and the AstroBase.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#

export iau2006, iau2006a

abstract type IAUModel end
abstract type IAU2006Model <: IAUModel end


########
# 2006 #
########

struct IAU2006A <: IAU2006Model end

"""
    iau2006a

The singleton instance of type `IAU2006a`, representing the IAU 2006A family of models.
"""
const iau2006a = IAU2006A()

struct IAU2006B <: IAU2006Model end

"""
    iau2006a

The singleton instance of type `IAU2006b`, representing the IAU 2006B family of models.
"""
const iau2006b = IAU2006B()