# Copyright 2021 TerraPower, LLC
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
Implementation of a chart of the nuclides class.
"""
from collections.abc import Iterable
import pathlib

from armi import context

class Element:
    r"""
    Represents an element, defined by its atomic number.

    Attributes
    ----------
    z : int
        atomic number, number of protons

    symbol : str
        element symbol

    name : str
        element name

    isotopes : list of Nuclide objects
        The represented isotopes of the element
    """

    def __init__(self, z, symbol, name):
        r"""
        Creates an instance of an Element.

        Parameters
        ----------
        z : int
            atomic number, number of protons

        symbol : str
            element symbol

        name: str
            element name

        """
        self.z = z
        self.symbol = symbol
        self.name = name
        self.standardWeight = None
        self.isotopes = []

        if other is not None and other == self:
            raise Exception(
                "Element with atomic weight {} already exists".format(self)
            )

    def __repr__(self):
        return "<Element {} {}>".format(self.symbol, self.z)

    def __eq__(self, other):
        return (
            self.z == other.z
            and self.symbol == other.symbol
            and self.name == other.name
        )

    def __hash__(self):
        return hash(self.name)

    def __iter__(self):
        for nuc in self.isotopes:
            yield nuc

    def append(self, nuclide):
        self.isotopes.append(nuclide)

    def isNaturallyOccurring(self):
        r"""
        Calculates the total natural abundance and if this value is zero returns False.
        If any isotopes are naturally occurring the total abundance will be >0 so it will return True
        """
        totalAbundance = 0.0
        for nuc in self.isotopes:
            totalAbundance += nuc.abundance
        return totalAbundance > 0.0

    def getNaturalIsotopics(self):
        """
        Return the nuclide bases of any naturally-occurring isotopes of this element.

        Notes
        -----
        Some elements have no naturally-occurring isotopes (Tc, Pu, etc.). To
        allow this method to be used in loops it will simply return an
        empty list in these situations.
        """
        return [nuc for nuc in self.isotopes if nuc.abundance > 0.0 and nuc.a > 0]

    def isHeavyMetal(self):
        return self.z > HEAVY_METAL_CUTOFF_Z

class ChartOfTheNuclides:
    """
    A representation of the chart of the nuclides.

    This is a container class that stores information about elements and the isotopes
    that belong to them. It is sometimes useful to have convenient access to properties
    of the contained elements (symbol, atomic number, natural abundance vector of
    isotopes, atomic mass, etc.) and of the specific nuclides (mass number, atomic mass,
    individual natural abundance, etc.). However, these should not be represented
    separately because they need to be mutually-consistent (e.g., the natural
    abundance of the C-14 nuclide should be the same as the Carbon element's idea of the
    isotope's abundance, and the atomic mass of Carbon should be the abundance-averaged
    atomic masses of its isotopes, etc.)

    These also should not be maintained as global state for all of the regular reasons
    why global state is problematic, and because it is entirely reasonable for one to
    desire multiple charts of the nuclides that make different approximations (e.g.,
    replacing certain nuclides with a lumped fission product, or neglecting certain
    nuclides altogether).
    """

    def __init__(self, nuclides: Iterable[Nuclide]):
        self.nuclides = []
        self.elements = []

        with open(pathlib.Path(context.RES) / "elements.dat", "r") as f:
            for line in f:
                # read z, symbol, and name
                lineData = line.split()
                z = int(lineData[0])
                sym = lineData[1].upper()
                name = lineData[2].lower()
                self.elements.append(Element(z, sym, name))


