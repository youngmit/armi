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

r"""
This module serves as a backwards-compatibility layer to the old nuclideBases
implementation.

The below docs corresponded to the "old" nuclideBases approach. These should be migrated
to new locations, as appropriate.

The nuclideBases module classes for providing *base* nuclide information, such as Z, A, state and energy release.

For details on how to setup the RIPL-3 data files to process extra nuclide data see:
:doc:`/user/user_install`

The nuclide class structure is outlined in :ref:`nuclide-bases-class-diagram`.

.. _nuclide-bases-class-diagram:

.. pyreverse:: armi.nucDirectory.nuclideBases
    :align: center
    :width: 90%

    Class inheritance diagram for :py:class:`Nuclide`.

See Also
--------
armi.nucDirectory.nuclideBases._addNuclideToIndices : builds this object
    
Examples
--------
>>> nuclideBases.byName['U235']
<NuclideBase U235: Z:92, A:235, S:0, label:U_6J, mc2id:U-2355>

>>> nuclideBases.byLabel['U_6J']
<NuclideBase U235: Z:92, A:235, S:0, label:U_6J, mc2id:U-2355>

Retrieve U-235 by the MC**2-v2 ID.

>>> nuclideBases.byMccId['U-2355']
<NuclideBase U235: Z:92, A:235, S:0, label:U_6J, mc2id:U-2355>
U235_7

Can get the same nuclide by the MC**2-v3 ID.

>>> nuclideBases.byMccId['U235_7']
<NuclideBase U235: Z:92, A:235, S:0, label:U_6J, mc2id:U-2355>

Retrieve U-235 by the MCNP ID.

>>> nuclideBases.byMcnpId['92235']
<NuclideBase U235: Z:92, A:235, S:0, label:U235, mc2id:U-2355>
U235_7

Retrieve U-235 by the AAAZZZS ID.

>>> nuclideBases.byAAAZZZSId['2350920']
<NuclideBase U235: Z:92, A:235, S:0, label:U235, mc2id:U-2355>
U235_7

.. exec::
    from tabulate import tabulate
    from armi.nucDirectory import nuclideBases

    attributes = ['name',
                  'a',
                  'z',
                  'state',
                  'weight',
                  'label',
                  'type']

    def getAttributes(nuc):
        return [
            '``{}``'.format(nuc.name),
            '``{}``'.format(nuc.a),
            '``{}``'.format(nuc.z),
            '``{}``'.format(nuc.state),
            '``{}``'.format(nuc.weight),
            '``{}``'.format(nuc.label),
            ':py:class:`~armi.nucDirectory.nuclideBases.{}`'.format(nuc.__class__.__name__),
        ]

    sortedNucs = sorted(nuclideBases.instances, key=lambda nb: (nb.z, nb.a))

    return create_table(tabulate(tabular_data=[getAttributes(nuc) for nuc in sortedNucs],
                                 headers=attributes,
                                 tablefmt='rst'),
                        caption='List of nuclides')

"""

import os
import pathlib
from typing import Optional, Union, List, Dict
import zlib

import yaml

import armi
from armi.nucDirectory import elements
from armi.nucDirectory.nuclides import Nuclide
from armi.nucDirectory.nuclides import LumpNuclideBase
from armi.nucDirectory.nuclides import DummyNuclideBase
from armi.localization import errors
from armi import runLog
from armi.utils.units import HEAVY_METAL_CUTOFF_Z

# used to prevent multiple applications of burn chains, which would snowball
# unphysically. This is a bit of a crutch for the global state that is the nuclide
# directory.
_burnChainImposed = False

instances = []

# Dictionary of Nuclides by the Nuclide.name for fast indexing
byName: Dict[str, Union[Nuclide, elements.Element]] = {}

byDBName = {}

byLabel = {}

byMccId = {}

byMcnpId = {}

byAAAZZZSId = {}

_riplEnvironVariable = "ARMI_RIPL_PATH"
RIPL_PATH = None


def getIsotopics(nucName: str) -> List[Nuclide]:
    """Expand elemental nuc name to isotopic nuc bases."""
    from armi.nucDirectory import chart
    obj = chart.byName[nucName]
    if isinstance(obj, (LumpNuclideBase, DummyNuclideBase)):
        # skip lumped fission products or dumps
        return []
    elif isinstance(obj, elements.Element):
        isotopics = obj.getNaturalIsotopics()
    else:
        isotopics = [obj]
    return isotopics


def nucNameFromDBName(dbName):
    """
    Return the nuc name of the given param name if the param name has a corresponding nuc name.

    If there is no nuc with that param name return None.
    """
    try:
        return byDBName[dbName].name
    except KeyError:
        return None


def isMonoIsotopicElement(name):
    """Return true if this is the only naturally occurring isotope of its element"""
    base = byName[name]
    return (
        base.abundance > 0
        and len([e for e in base.element.isotopes if e.abundance > 0]) == 1
    )


def where(predicate):
    directory.where(predicate)


def single(predicate):
    directory.single(predicate)


def changeLabel(nuclideBase, newLabel):
    nuclideBase.label = newLabel
    byLabel[newLabel] = nuclideBase


def __readRiplNuclides():
    """
    Initialize all nuclides with experimentally-measured masses.

    This includes roughly 4000 nuclides and should represent anything we ever
    want to model. This builds the large set of NuclideBases available.

    RIPL is the Reference Input Parameter Library (RIPL-3), which can be found at
    https://www-nds.iaea.org/RIPL-3/.
    """
    from armi.nuclearDataIO import ripl

    elements.clearNuclideBases()
    for z, a, symbol, mass, _err in ripl.readFRDMMassFile(
        os.path.join(armi.context.RES, "ripl-mass-frdm95.dat")
    ):
        if z == 0 and a == 1:
            # skip the neutron
            continue
        element = elements.bySymbol[symbol.upper()]
        NuclideBase(element, a, mass, 0, 0, None)


def __readRiplAbundance():
    """
    Read natural abundances of any natural nuclides.

    This adjusts already-existing NuclideBases and Elements with the new information.
    """
    from armi.nuclearDataIO import ripl

    with open(os.path.join(armi.context.RES, "ripl-abundance.dat")) as ripl_abundance:
        for _z, a, sym, percent, _err in ripl.readAbundanceFile(ripl_abundance):
            nb = byName[sym + "{}".format(a)]
            nb.abundance = percent / 100.0


def __readMc2Nuclides():
    """
    Read nuclides as defined in the MC2 library.

    Notes
    -----
    This assigns MC2 labels and often adds metastable versions of nuclides
    that have already been added from RIPL.
    """
    with open(os.path.join(armi.context.RES, "mc2Nuclides.yaml"), "r") as mc2Nucs:
        mc2Nuclides = yaml.load(mc2Nucs, Loader=yaml.FullLoader)
    # now add the mc2 specific nuclideBases, and correct the mc2Ids when a > 0 and state = 0
    for name, data in mc2Nuclides.items():
        z = data["z"]
        a = data["a"]
        state = data["state"]
        iid = data["id"]

        element = elements.byZ[z] if z > 0 else None
        if z == 0:
            weight = data["weight"]
            if "LFP" in name:
                LumpNuclideBase(name, z, iid, weight)
            # Allows names like REGXX to be in the ISOTXS file (macroscopic/region xs considered as 'lumped' parameters)
            elif "LREGN" in name:
                LumpNuclideBase(name, z, iid, weight)
            else:
                DummyNuclideBase(name, iid, weight)
        # I am considering removal of the concept of "natural" nuclides. These really
        # are just Elements, and the only reason that i can think of maintaining them
        # like this is for some form of unified concept of  "what can be used as a key
        # in a numberDensity array". Might as well allow Elements to function like this.
        # If we want to get more formal about what it means for something to function as
        # a key in a number-density dictionary, let's call a spade a spade and invent
        # `NumberDensityKey` (essentially I think what INuclide was originally intended
        # to be? The name probably didn't do it any favors...).
        # Natrually, some lattice physics codes/nuclear
        # data libraries may only have some of the Elements represented directly, but
        # really its on them to sort that out, not ARMI.
        # elif a == 0:
        #     NaturalNuclideBase(name, element, iid)
        else:
            # state == 0 nuclide *should* already exist
            needToAdd = True
            if state == 0:
                clide = [
                    nn
                    for nn in element.isotopes
                    if nn.z == z and nn.a == a and nn.state == state
                ]
                if len(clide) > 1:
                    raise ValueError(
                        "More than 1 nuclide meets specific criteria: {}".format(clide)
                    )
                needToAdd = len(clide) == 0
                if not needToAdd and iid:
                    clide[0].mc2id = iid
                    byMccId[iid] = clide[0]
            # state != 0, nuclide should not exist, create it
            if needToAdd:
                NuclideBase(
                    element, a, element.standardWeight or float(a), 0.0, state, iid
                )

    # special case AM242. Treat the metastable as the main one and specify ground state as AM242G.
    # This is a typical approach in many reactor codes including MCNP since you almost always
    # are interested in AM242M.
    am242g = byName["AM242"]
    am242g.name = "AM242G"
    am242 = byName["AM242M"]
    am242.name = "AM242"
    am242.weight = am242g.weight  # use RIPL mass for metastable too
    byName[am242.name] = am242
    byDBName[am242.getDatabaseName()] = am242
    byName["AM242G"] = am242g
    byDBName[byName["AM242G"].getDatabaseName()] = am242g


def getDepletableNuclides(activeNuclides, obj):
    """Get nuclides in this object that are in the burn chain."""
    return sorted(set(activeNuclides) & set(obj.getNuclides()))


def imposeBurnChain(burnChainStream):
    """
    Apply transmutation and decay information to each nuclide.

    Notes
    -----
    You cannot impose a burn chain twice. Doing so would require that you clean out the
    transmutations and decays from all the module-level nuclide bases, which generally
    requires that you rebuild them. But rebuilding those is not an option because some
    of them get set as class-level attributes and would be orphaned. If a need to change
    burn chains mid-run re-arises, then a better nuclideBase-level burnchain cleanup
    should be implemented so the objects don't have to change identity.

    Notes
    -----
    We believe the transmutation information would probably be better stored on a
    less fundamental place (e.g. not on the NuclideBase).

    See Also
    --------
    armi.nucDirectory.transmutations : describes file format
    """
    global _burnChainImposed  # pylint: disable=global-statement
    if _burnChainImposed:
        # the only time this should happen is if in a unit test that has already
        # processed conftest.py and is now building a Case that also imposes this.
        runLog.warning("Burn chain already imposed. Skipping reimposition.")
        return
    _burnChainImposed = True
    burnData = yaml.load(burnChainStream, Loader=yaml.FullLoader)
    for nucName, burnInfo in burnData.items():
        nuclide = byName[nucName]
        # think of this protected stuff as "module level protection" rather than class.
        nuclide._processBurnData(burnInfo)  # pylint: disable=protected-access


def factory():
    """
    Reads data files to instantiate the :py:class:`Nuclides <Nuclide>`.

    Reads NIST, MC**2 and burn chain data files to instantiate the :py:class:`Nuclides <Nuclide>`.
    Also clears and fills in the
    :py:data:`~armi.nucDirectory.nuclideBases.instances`,
    :py:data:`byName`, :py:attr:`byLabel`, and
    :py:data:`byMccId` module attributes. This method is automatically run upon
    loading the module, hence it is not usually necessary to re-run it unless there is a
    change to the data files, which should not happen during run time, or a *bad*
    :py:class`Nuclide` is created.

    Notes
    -----
    This may cannot be run more than once. NuclideBase instances are used throughout the ARMI
    ecosystem and are even class attributes in some cases. Re-instantiating them would orphan
    any existing ones and break everything.

    Nuclide labels from MC2-2, MC2-3, and MCNP are currently handled directly.
    Moving forward, we plan to implement a more generic labeling system so that
    plugins can provide code-specific nuclide labels in a more extensible fashion.
    """
    # this intentionally clears and reinstantiates all nuclideBases
    global instances  # pylint: disable=global-statement
    if len(instances) == 0:
        # make sure the elements actually exist...
        elements.factory()
        del instances[:]  # there is no .clear() for a list
        byName.clear()
        byDBName.clear()
        byLabel.clear()
        byMccId.clear()
        byMcnpId.clear()
        byAAAZZZSId.clear()
        __readRiplNuclides()
        __readRiplAbundance()
        # load the mc2Nuclide.json file. This will be used to supply nuclide IDs
        __readMc2Nuclides()
        elements.deriveNaturalWeights()
        __readRiplDecayData()
        # reload the thermal scattering library with the new nuclideBases too
        # pylint: disable=import-outside-toplevel; cyclic import
        from . import thermalScattering

        thermalScattering.factory()


def __readRiplDecayData():
    """
    Read in the RIPL-3 decay data files and update nuclide bases.

    Notes
    -----
    This makes an assumption that the RIPL-3 data files have a
    `z???.dat` naming convention and assumes that there are 118
    total data files in the package.

    The processing is skipped if the ``ARMI_RIPL_PATH`` environment
    variable has not been set.

    Raises
    ------
    ValueError
        If the ``ARMI_RIPL_PATH`` is defined, but set incorrectly.
    """
    from armi.nuclearDataIO import ripl

    global RIPL_PATH

    riplPath = os.environ.get(_riplEnvironVariable, None)
    if riplPath is None:
        return None

    path = pathlib.Path(riplPath)
    if not path.exists() or not path.is_dir():
        raise ValueError(f"`{_riplEnvironVariable}`: {path} is invalid.")

    # Check for all (.dat) data files within the directory. These
    # are ordered from z000.dat to z117.dat. If all files do not
    # exist then an exception is thrown for the missing data files.
    numRIPLDataFiles = 118
    dataFileNames = ["z{:>03d}.dat".format(i) for i in range(0, numRIPLDataFiles)]
    missingFileNames = []
    for df in dataFileNames:
        expectedDataFilePath = os.path.abspath(os.path.join(path, df))
        if not os.path.exists(expectedDataFilePath):
            missingFileNames.append(df)
    if missingFileNames:
        raise ValueError(
            f"There are {len(missingFileNames)} missing RIPL data files in `{_riplEnvironVariable}`: {path}.\n"
            f"The following data files were expected: {missingFileNames}"
        )

    ripl.makeDecayConstantTable(directory=path)
    RIPL_PATH = path


def _renormalizeElementRelationship():
    for nuc in instances:
        if nuc.element is not None:
            nuc.element = elements.byZ[nuc.z]
            nuc.element.append(nuc)


elements.nuclideRenormalization = _renormalizeElementRelationship


def _addNuclideToIndices(nuc):
    instances.append(nuc)
    byName[nuc.name] = nuc
    byDBName[nuc.getDatabaseName()] = nuc
    byLabel[nuc.label] = nuc
    if nuc.mc2id:
        byMccId[nuc.mc2id] = nuc
    mc3 = nuc.getMcc3Id()
    if mc3:
        byMccId[mc3] = nuc

    mcnp = nuc.getMcnpId()
    if mcnp is not None:
        byMcnpId[mcnp] = nuc
    try:
        byAAAZZZSId[nuc.getAAAZZZSId()] = nuc
    except AttributeError:
        pass


def initReachableActiveNuclidesThroughBurnChain(numberDensityDict, activeNuclides):
    """
    March through the depletion chain and find all nuclides that can be reached by depleting nuclides passed in.

    This limits depletion to the smallest set of nuclides that matters.

    Parameters
    ----------
    numberDensityDict : dict
        Starting number densities.

    activeNuclides : OrderedSet
        Active nuclides defined on the reactor blueprints object. See: armi.reactor.blueprints.py
    """
    missingActiveNuclides = set()
    memo = set()
    difference = set(numberDensityDict).difference(memo)
    while any(difference):
        nuclide = difference.pop()
        memo.add(nuclide)
        # Skip the nuclide if it is not `active` in the burn-chain
        if not nuclide in activeNuclides:
            continue
        nuclideObj = byName[nuclide]
        for interaction in nuclideObj.trans + nuclideObj.decays:
            try:
                # Interaction nuclides can only be added to the number density
                # dictionary if they are a part of the user-defined active nuclides
                productNuclide = interaction.getPreferredProduct(activeNuclides)
                if productNuclide not in numberDensityDict:
                    numberDensityDict[productNuclide] = 0.0
            except KeyError:
                # Keep track of the first production nuclide
                missingActiveNuclides.add(interaction.productNuclides)

        difference = set(numberDensityDict).difference(memo)

    if missingActiveNuclides:
        _failOnMissingActiveNuclides(missingActiveNuclides)


def _failOnMissingActiveNuclides(missingActiveNuclides):
    """Raise ValueError with notification of which nuclides to include in the burn-chain."""
    msg = "Missing active nuclides in loading file. Add the following nuclides:"
    for i, nucList in enumerate(missingActiveNuclides, 1):
        msg += "\n {} - ".format(i)  # Index of
        for j, nuc in enumerate(nucList, 1):
            delimiter = " or " if j < len(nucList) else ""
            msg += "{}{}".format(nuc, delimiter)
    raise ValueError(msg)
