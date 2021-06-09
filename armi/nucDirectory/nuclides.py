from typing import Optional

from armi import runLog
from armi.nucDirectory import transmutations
from armi.nucDirectory import elements


# lookup table from https://t2.lanl.gov/nis/data/endf/endfvii-n.html
BASE_ENDFB7_MAT_NUM = {
    "PM": 139,
    "RA": 223,
    "AC": 225,
    "TH": 227,
    "PA": 229,
    "NP": 230,
    "PU": 235,
    "AM": 235,
    "CM": 240,
    "BK": 240,
    "CF": 240,
    "TC": 99,
}


class Auto:
    pass


class Nuclide:
    r"""
    Represents an individual nuclide

    This is meant to be a lookup class, with just one per unique nuclide
    in the problem. That means there won't be copies of these in every block
    that has these nuclideBases.

    Attributes
    ----------
    z : int
        Number of protons.

    a : int
        Number of nucleons.

    state : int
        Indicates excitement, 1 is more excited than 0.

    abundance : float
        Isotopic fraction of a naturally occurring nuclide. The sum of all nuclide
        abundances for a naturally occurring element should be 1.0. This is atom
        fraction, not mass fraction.

    name : str
        ARMI's unique name for the given nuclide.

    label : str
        ARMI's unique 4 character label for the nuclide.
        These are not human readable, but do not lose any information.
        The label is effectively the
        :attr:`Element.symbol `armi.nucDirectory.elements.Element.symbol`
        padded to two characters, plus the mass number (A) in base-26 (0-9, A-Z).
        Additional support for meta-states is provided by adding 100 * the state
        to the mass number (A).

    mc2id : str
        The unique id used by MC**2 version 2.

    nuSF : float
        Neutrons released per spontaneous fission.
        This should probably be moved at some point.
    """

    fissile = ["U235", "PU239", "PU241", "AM242M", "CM244", "U233"]
    TRANSMUTATION = "transmutation"
    DECAY = "decay"
    SPONTANEOUS_FISSION = "nuSF"

    def __init__(
        self,
        element: Optional[elements.Element],
        z: Optional[int],
        a,
        state,
        weight,
        abundance,
        name: Union[Auto, str],
        label,
        mc2id,
    ):
        r"""
        Create an instance of an Nuclide.

        .. warning::
            Do not call this constructor directly; use the factory instead.

        """
        if not (element is None) ^ (z is None):
            raise ValueError("Either `element` or `z` must be provided, but not both")
        self.z = z or element.z
        self.a = a
        self.state = state
        self.decays = []
        self.trans = []
        self.weight = weight
        self.abundance = abundance
        self.name = name if name is not Auto else self._createName(element, a, state)
        self.label = (
            label if label is not Auto else self._createLabel(element, a, state)
        )
        self.element = element
        self.mc2id = mc2id
        self.nuSF = 0.0

    def __hash__(self):
        return hash((self.a, self.z, self.state))

    def __reduce__(self):
        return fromName, (self.name,)

    def __repr__(self):
        return "<Nuclide {}: Z:{}, A:{}, S:{}, label:{}, mc2id:{}>".format(
            self.name, self.z, self.a, self.state, self.label, self.mc2id
        )

    def __lt__(self, other):
        return (self.z, self.a, self.state) < (other.z, other.a, other.state)

    def _processBurnData(self, burnInfo):
        """
        Process YAML burn transmutation, decay, and spontaneous fission data for this nuclide.

        This clears out any existing transmutation/decay information before processing.

        Parameters
        ----------
        burnInfo: list
            List of dictionaries containing burn information for the current nuclide
        """
        self.decays = []
        self.trans = []
        for nuclideBurnCategory in burnInfo:
            # Check that the burn category has only one defined burn type
            if len(nuclideBurnCategory) > 1:
                raise ValueError(
                    "Improperly defined ``burn-chain`` of {}. {} should be a single burn type.".format(
                        self, nuclideBurnCategory.keys()
                    )
                )
            nuclideBurnType = list(nuclideBurnCategory.keys())[0]
            if nuclideBurnType == self.TRANSMUTATION:
                self.trans.append(
                    transmutations.Transmutation(
                        self, nuclideBurnCategory[nuclideBurnType]
                    )
                )
            elif nuclideBurnType == self.DECAY:
                self.decays.append(
                    transmutations.DecayMode(self, nuclideBurnCategory[nuclideBurnType])
                )
            elif nuclideBurnType == self.SPONTANEOUS_FISSION:
                self.nuSF = nuclideBurnCategory[nuclideBurnType]
            else:
                raise Exception(
                    "Undefined Burn Data {} for {}. Expected {}, {}, or {}."
                    "".format(
                        nuclideBurnType,
                        self,
                        self.TRANSMUTATION,
                        self.DECAY,
                        self.SPONTANEOUS_FISSION,
                    )
                )

    def getDecay(self, decayType):
        r"""Get a :py:class:`~armi.nucDirectory.transmutations.DecayMode`.

        Retrieve the first :py:class:`~armi.nucDirectory.transmutations.DecayMode`
        matching the specified decType.

        Parameters
        ----------
        decType: str
            Name of decay mode e.g. 'sf', 'alpha'

        Returns
        -------
        decay : :py:class:`DecayModes <armi.nucDirectory.transmutations.DecayMode>`

        """
        for d in self.decays:
            if d.type == decayType:
                return d

    def isFissile(self):
        r"""Determine if the nuclide is fissile.

        Determines if the nuclide is fissle.

        Returns
        -------
        answer: bool
            True if the :py:class:`Nuclide` is fissile, otherwise False.
        """
        return self.name in self.fissile

    def getDatabaseName(self):
        """Get the name of the nuclide used in the database (i.e. "nPu239")"""
        return "n{}".format(self.name.capitalize())

    def isHeavyMetal(self):
        return self.z > HEAVY_METAL_CUTOFF_Z

    @staticmethod
    def _createName(element, a, state):
        # state is either 0 or 1, so some nuclides will get an M at the end
        metaChar = ["", "M"]
        return "{}{}{}".format(element.symbol, a, metaChar[state])

    @staticmethod
    def _createLabel(element, a, state):
        """
        Make label for nuclide base.

        The logic causes labels for things with A<10 to be zero padded like H03 or tritium
        instead of H3. This avoids the metastable tritium collision which would look
        like elemental HE. It also allows things like MO100 to be held within 4 characters,
        which is a constraint of the ISOTXS format if we append 2 characters for XS type.
        """
        # len(e.symbol) is 1 or 2 => a % (either 1000 or 100)
        #                         => gives exact a, or last two digits.
        # the division by 10 removes the last digit.
        firstTwoDigits = (a % (10 ** (4 - len(element.symbol)))) // 10
        # the last digit is either 0-9 if state=0, or A-J if state=1
        lastDigit = "0123456789" "ABCDEFGHIJ"[(a % 10) + state * 10]

        return "{}{}{}".format(element.symbol, firstTwoDigits, lastDigit)

    def isNatural(self):
        """
        A Nuclide is "natural" if it has a zero mass number.

        These are useful as a sort of Element-in-Nuclide's-clothing to define number
        densities of a whole element's worth of natrual isotopics. In future it may make
        sense to be more explicit about this and simply support Elements as keys in
        number densities.
        """
        return self.a == 0

    def getNaturalIsotopics(self):
        """Gets the natural isotopics root :py:class:`~elements.Element`.

        Gets the naturally occurring nuclides for this nuclide.

        Returns
        -------
        nuclides: list
            List of :py:class:`Nuclides <Nuclide>`

        See Also
        --------
        :meth:`Nuclide.getNaturalIsotopics`
        """
        return self.element.getNaturalIsotopics()

    def getMcc3Id(self):
        """Gets the MC**2-v3 nuclide ID.

        Returns
        -------
        name: str
            The MC**2 ID: ``AM42M7``, ``B10__7``, etc.

        See Also
        --------
        :meth:`Nuclide.getMcc3Id`
        """
        base = ""
        if self.state > 0:
            base = "{}{}M".format(self.element.symbol, self.a % 100)
        else:
            base = "{}{}".format(self.element.symbol, self.a)
        return "{:_<5}7".format(base)

    def getMcnpId(self) -> Optional[str]:
        """
        Gets the MCNP label for this nuclide

        Returns
        -------
        id : str
            The MCNP ID e.g. ``92235``, ``94239``, ``6000``

        """
        z, a = self.z, self.a

        if z == 95 and a == 242:
            # Am242 has special rules
            if self.state != 1:
                # MCNP uses base state for the common metastable state AM242M , so AM242M is just 95242
                # AM242 base state is called 95642 (+400) in mcnp.
                # see https://mcnp.lanl.gov/pdf_files/la-ur-08-1999.pdf
                # New ACE-Formatted Neutron and Proton Libraries Based on ENDF/B-VII.0
                a += 300 + 100 * max(self.state, 1)
        elif self.state > 0:
            # in general mcnp adds 300 + 100*m to the Z number for metastables. see above source
            a += 300 + 100 * self.state

        return "{z:d}{a:03d}".format(z=z, a=a)

    def getAAAZZZSId(self):
        """
        Gets the AAAZZZS label for this nuclide

        Returns
        -------
        id : str
            The MCNP ID e.g. ``2350920``, ``2390940``, ``1200600``

        """

        aaa = "{}".format(self.a)
        zzz = "{0:03}".format(self.z)
        s = "1" if self.state > 0 else "0"

        return "{}{}{}".format(aaa, zzz, s)

    def getSerpentId(self):
        """
        Returns the SERPENT style ID for this nuclide.

        Returns
        -------
        id: str
            The ID of this nuclide based on it's elemental name, weight,
            and state, eg ``U-235``, ``Te-129m``,
        """
        symbol = self.element.symbol.capitalize()
        return "{}-{}{}".format(symbol, self.a, "m" if self.state else "")

    def getEndfMatNum(self):
        """
        Gets the ENDF MAT number

        MAT numbers are defined as described in section 0.4.1 of the NJOY manual.
        Basically, it's Z * 100 + I where I is an isotope number. I=25 is defined
        as the lightest known stable isotope of element Z, so for Uranium,
        Z=92 and I=25 refers to U234. The values of I go up by 3 for each
        mass number, so U235 is 9228. This leaves room for three isomeric
        states of each nuclide.

        Returns
        -------
        id : str
            The MAT number e.g. ``9237`` for U238

        """
        z, a = self.z, self.a
        if self.element.symbol in BASE_ENDFB7_MAT_NUM:
            # no stable isotopes (or other special case). Use lookup table
            smallestStableA = BASE_ENDFB7_MAT_NUM[self.element.symbol]
        else:
            naturalIsotopes = self.getNaturalIsotopics()
            if naturalIsotopes:
                smallestStableA = min(
                    ni.a for ni in naturalIsotopes
                )  # no guarantee they were sorted
            else:
                raise KeyError("Nuclide {0} is unknown in the MAT number lookup")

        isotopeNum = (a - smallestStableA) * 3 + self.state + 25
        mat = z * 100 + isotopeNum
        return "{0}".format(mat)


class NaturalNuclideBase(Nuclide):
    def __init__(self, name, element, mc2id):
        self.element = element
        Nuclide.__init__(
            self,
            element,
            element.z,
            0,
            0,
            sum([nn.weight * nn.abundance for nn in element.getNaturalIsotopics()]),
            0.0,  # keep abundance 0.0 to not interfere with the isotopes
            name,
            name,
            mc2id,
        )
        self.element.append(self)

    def __repr__(self):
        return "<NaturalNuclideBase {}: Z:{}, w:{}, label:{}, mc2id:{}>" "".format(
            self.name, self.z, self.weight, self.label, self.mc2id
        )

    def getNaturalIsotopics(self):
        r"""Gets the natural isotopics root :py:class:`~elements.Element`.

        Gets the naturally occurring nuclides for this nuclide.

        Returns
        -------
        nuclides: list
            List of :py:class:`Nuclides <Nuclide>`.

        See Also
        --------
        :meth:`Nuclide.getNaturalIsotopics`
        """
        return self.element.getNaturalIsotopics()

    def getMcc3Id(self):
        r"""Gets the MC**2-v3 nuclide ID.

        Returns
        -------
        id: str
            The MC**2 ID: ``FE___7``, ``C____7``, etc.

        See Also
        --------
        :meth:`Nuclide.getMcc3Id`
        """
        return "{:_<5}7".format(self.element.symbol)

    def getAAAZZZSId(self):
        """Gets the AAAZZZS ID for a few elements.

        Notes
        -----
        the natural nuclides 'C' and 'V' do not have isotopic nuclide data for MC2 so sometimes they tag along in the
        list of active nuclides. This method is designed to fail in the same as if there was not getAAAZZZSId method
        defined.
        """
        if self.element.symbol == "C":
            return "120060"
        elif self.element.symbol == "V":
            return "510230"
        else:
            return None

    def getSerpentId(self):
        """Gets the SERPENT ID for this natural nuclide.

        Returns
        -------
        id: str
            SERPENT ID: ``C-nat``, `Fe-nat``
        """
        return "{}-nat".format(self.element.symbol.capitalize())

    def getEndfMatNum(self):
        """Get the ENDF mat number for this element."""
        if self.z != 6:
            runLog.warning(
                "The only elemental in ENDF/B VII.1 is carbon. "
                "ENDF mat num was requested for the elemental {} and will not be helpful "
                "for working with ENDF/B VII.1. Try to expandElementalsToIsotopics".format(
                    self
                )
            )
        return "{0}".format(self.z * 100)


class LumpNuclideBase(Nuclide):
    """
    Lump nuclides are used for lumped fission products.

    See Also
    --------
    armi.physics.neutronics.fissionProduct model:
        Describes what nuclides LumpNuclideBase is expend to.
    """

    def __init__(self, name, z, mc2id, weight):
        Nuclide.__init__(self, z, 0, 0, weight, 0.0, name, name[1:], mc2id)

    def __repr__(self):
        return "<LumpNuclideBase {}: Z:{}, w:{}, label:{}, mc2id:{}>" "".format(
            self.name, self.z, self.weight, self.label, self.mc2id
        )

    def getNaturalIsotopics(self):
        r"""Gets the natural isotopics, an empty iterator.

        Gets the naturally occurring nuclides for this nuclide.

        Returns
        -------
        empty: iterator
            An empty generator

        See Also
        --------
        :meth:`Nuclide.getNaturalIsotopics`
        """
        return
        yield

    def getMcc3Id(self):
        r"""Gets the MC**2-v3 nuclide ID.

        Returns
        -------
        name: str
            The MC**2 ID: ``LFP38``, etc.

        See Also
        --------
        :meth:`Nuclide.getMcc3Id`
        """
        return self.mc2id


class DummyNuclideBase(Nuclide):
    """
    Dummy nuclides are used nuclides which transmute into isotopes that are not defined in blueprints.

    Notes
    -----
    If DMP number density is not very small, cross section may be artifically depressed.
    """

    def __init__(self, name, mc2id, weight):
        Nuclide.__init__(
            self,
            None,
            0,
            0,
            0,
            weight,
            0.0,
            name,
            "DMP" + name[4],
            mc2id,  # z  # a  # state
        )

    def __repr__(self):
        return "<DummyNuclideBase {}: Z:{}, w:{}, label:{}, mc2id:{}>" "".format(
            self.name, self.z, self.weight, self.label, self.mc2id
        )

    def getNaturalIsotopics(self):
        r"""Gets the natural isotopics, an empty iterator.

        Gets the naturally occurring nuclides for this nuclide.

        Returns
        -------
        empty: iterator
            An empty generator

        See Also
        --------
        :meth:`Nuclide.getNaturalIsotopics`
        """
        return
        yield

    def getMcc3Id(self):
        r"""Gets the MC**2-v3 nuclide ID.

        Returns
        -------
        name: str
            The MC**2 ID: ``DUMMY`` for all.

        See Also
        --------
        :meth:`Nuclide.getMcc3Id`
        """
        return "DUMMY"

    def getMcnpId(self):
        return None


def fromName(name):
    r"""Get a nuclide from its name."""
    matches = [nn for nn in instances if nn.name == name]
    if len(matches) != 1:
        raise errors.nuclides_TooManyOrTooFew_number_MatchesForNuclide_name(
            len(matches), name
        )
    return matches[0]


