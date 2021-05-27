# Copyright 2019 TerraPower, LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Deals with elements of the periodic table.
"""

import os
from typing import List, Dict

from armi import context
from armi.utils.units import HEAVY_METAL_CUTOFF_Z
from armi.nucDirectory import chart

# Re-export Element since most of the code still expects it to be defined here
from armi.nucDirectory.chart import Element


# These remain for backwards compatibility of the old global Elements and NuclideBases.
# Since they are essentially references back to the global container, they are managed
# by it; do not set them here
byZ: Dict[int, Element] = {}
byName: Dict[str, Element] = {}
bySymbol: Dict[str, Element] = {}
