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
from typing import List, Iterator
import pathlib

import numpy

from armi import context
from armi.nucDirectory import nuclideBases
from armi.utils.units import HEAVY_METAL_CUTOFF_Z


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

    isotopes : list of Nuclides
        Nuclides that make up the isotopes of this element. This includes a Nuclide that
        represents the natural element itself, with an abundance of 0.0.
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
        self.isotopes = set()

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
        self.isotopes.add(nuclide)

    def isNaturallyOccurring(self):
        r"""
        Return true if any isotopes have a greater-than-zero natural abundance.
        """
        return any(nuc.abundance > 0 for nuc in self.isotopes)

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

    def updateWeight(self):
        """
        Re-compute the element's standard weight based on its isotopes.
        """
        totalAbundance = sum(nuc.abundance for nuc in self.isotopes)
        if totalAbundance > 0.0:
            self.standardWeight = (
                sum(nuc.weight * nuc.abundance for nuc in self.isotopes)
                / totalAbundance
            )


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

    class Elements:
        def __init__(self):
            self.byZ = dict()
            self.byName = dict()
            self.bySymbol = dict()
            self._elements = set()

            with open(pathlib.Path(context.RES) / "elements.dat", "r") as f:
                for line in f:
                    # read z, symbol, and name
                    lineData = line.split()
                    z = int(lineData[0])
                    sym = lineData[1].upper()
                    name = lineData[2].lower()
                    self.add(Element(z, sym, name))

        def __iter__(self) -> Iterator[Element]:
            return iter(self._elements)

        def add(self, element: Element):
            self._elements.add(element)
            self.byZ[element.z] = element
            self.byName[element.name] = element
            self.bySymbol[element.symbol] = element

        def update(self, elements: List[Element]):
            for element in elements:
                self.add(element)

        def setIsotopics(self, nuclides):
            for nuc in nuclides:
                self.byZ[nuc.z].append(nuc)

        def clear(self):
            self.byZ.clear()
            self.byName.clear()
            self.bySymbol.clear()
            self._elements = set()

    class Nuclides:
        """
        Rudimentary container for Nuclide instances

        This is nested within the ``ChartOfTheNuclides`` class because of the
        interrelated nature of Nuclides (their references to their Elements) and
        Elements (with their references to isotope Nuclides).
        """

        def __init__(self):
            self._nuclides = set()
            self.byName = dict()
            self.byDbName = dict()
            self.byLabel = dict()

        def __iter__(self):
            return iter(self._nuclides)

        def add(self, nuclide):
            """
            Insert a new Nuclide instance.

            Raises ValueError if the nuclide already exists.
            """
            if nuclide in self._nuclides:
                raise ValueError("Nuclide `{}` already in directory".format(nuclide))

            self._nuclides.add(nuclide)
            self.byName[nuclide.name] = nuclide
            self.byDbName[nuclide.getDatabaseName] = nuclide
            self.byLabel[nuclide.label] = nuclide

        def update(self, other):
            """
            Fold entries from another set of Nuclides into this one.

            Nuclide collisions between self and other will produce exceptions.
            """
            for nuclide in other:
                self.add(other)

        def where(self, predicate):
            r"""Get all :py:class:`Nuclides <Nuclide>` matching a condition.

            Returns an iterator of :py:class:`Nuclides <Nuclide>` matching the specified condition.

            Attributes
            ----------

            predicate: lambda
                A lambda, or function, accepting a :py:class:`Nuclide` as a parameter

            Examples
            --------

            >>> from armi.nucDirectory import nuclideBases
            >>> nucDir = nuclideBases.NuclideDirectory()
            >>> ...
            >>> [nn.name for nn in nucDir.where(lambda nb: 'Z' in nb.name)]
            ['ZN64', 'ZN66', 'ZN67', 'ZN68', 'ZN70', 'ZR90', 'ZR91', 'ZR92', 'ZR94', 'ZR96', 'ZR93', 'ZR95', 'ZR']

            >>> # in order to get length, convert to list
            >>> isomers90 = list(nucDir.where(lambda nb: nb.a == 95))
            >>> len(isomers90)
            3
            >>> for iso in isomers: print(iso)
            <NuclideBase MO95: Z:42, A:95, S:0, label:MO2N, mc2id:MO95 5>
            <NuclideBase NB95: Z:41, A:95, S:0, label:NB2N, mc2id:NB95 5>
            <NuclideBase ZR95: Z:40, A:95, S:0, label:ZR2N, mc2id:ZR95 5>

            """
            for nuc in self._nuclides:
                if predicate(nuc):
                    yield (nuc)

        def single(self, predicate):
            r"""Get a single :py:class:`Nuclide` meeting the specified condition.

            Similar to :py:func:`where`, this function uses a lambda input to filter
            the :py:attr:`Nuclide instances <instances>`. If there is not 1 and only
            1 match for the specified condition, an exception is raised.

            Examples
            --------

            >>> from armi.nucDirectory import nuclideBases
            >>> nucDir = nuclideBases.NuclideDirectory()
            >>> ...
            >>> nucDir.single(lambda nb: nb.name == 'C')
            <NaturalNuclideBase C: Z:6, w:12.0107358968, label:C, mc2id:C    5>

            >>> nucDir.single(lambda nb: nb.z == 95 and nb.a == 242 and nb.state == 1)
            <NuclideBase AM242M: Z:95, A:242, S:1, label:AM4C, mc2id:AM242M>

            """
            matches = [nuc for nuc in self._nuclides if predicate(nuc)]
            if len(matches) != 1:
                raise ValueError(
                    "Expected one nuclide matching the predicate, but got `{}`".format(
                        matches
                    )
                )
            return matches[0]

    def __init__(self, nuclides: List[nuclideBases.Nuclide]):
        self.nuclides = self.Nuclides()
        self.elements = self.Elements()

        with open(pathlib.Path(context.RES) / "elements.dat", "r") as f:
            for line in f:
                # read z, symbol, and name
                lineData = line.split()
                z = int(lineData[0])
                sym = lineData[1].upper()
                name = lineData[2].lower()
                self.elements.add(Element(z, sym, name))

        for nuclide in nuclides:
            self.nuclides.add(nuclide)

        # make sure that the elements have data consistent with their
        # isotopes/abundances
        self._connectElements()

    def _connectElements(self) -> None:
        for nuc in self.nuclides:
            self.elements.byZ[nuc.z].append(nuc)

        for element in self.elements:
            element.updateWeight()


def factory():
    # pylint:disable=import-location
    # import here, since we need to associate element's collections with ours, while
    # elements is importing us for the actual definition of the Element class. Kind of
    # bonkers, but elements still exists mostly to support old code that used to use the
    # global collections.
    #
    # TODO consider moving the Element class defn back to elements. I think i only put
    # it here because i was still working around elements being in charge of its
    # collections.
    from armi.nucDirectory import elements

    global chart
    chart = ChartOfTheNuclides([])

    elements.byZ = chart.elements.byZ
    elements.byName = chart.elements.byName
    elements.bySymbol = chart.elements.bySymbol
