/*
 * ElementImpl.cpp
 *
 * Copyright (c) 2011,2012 Mathias Wilhelm
 * Copyright (c) 2011,2012 Marc Kirchner
 *
 * This file is part of the Mass Spectrometry Toolkit (MSTK).
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "MSTK/aas/ElementImpl.hpp"
#include "MSTK/common/Error.hpp"

namespace mstk {
namespace aas {
namespace elements {

const Size nEntries = 107;

// TODO are there other important heavy isotopes?
const String symbols[] = { "H", "He", "Li", "Be", "B", "C", "N", "O", "F",
                                "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl",
                                "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn",
                                "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As",
                                "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb",
                                "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
                                "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La",
                                "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
                                "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta",
                                "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl",
                                "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac",
                                "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk",
                                "Cf", "Es", "Fm", "Md", "No", "Lr", "2H",
                                "13C", "15N", "18O" };

const Size nIsotopes[] = { 2, 2, 2, 1, 2, 2, 2, 3, 1, 3, 1, 3, 1, 3, 1, 4,
                                2, 3, 3, 6, 1, 5, 2, 4, 1, 4, 1, 5, 2, 5, 2, 5,
                                1, 6, 2, 6, 2, 4, 1, 5, 1, 7, 1, 7, 1, 6, 2, 8,
                                2, 10, 2, 8, 1, 9, 1, 7, 2, 4, 1, 7, 1, 7, 2,
                                7, 1, 7, 1, 6, 1, 7, 2, 6, 2, 5, 2, 7, 2, 6, 1,
                                7, 2, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1,
                                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };

const Size atomicNumbers[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
                                    14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,
                                    25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35,
                                    36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46,
                                    47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57,
                                    58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68,
                                    69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
                                    80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90,
                                    91, 92, 93, 94, 95, 96, 97, 98, 99, 100,
                                    101, 102, 103, 1, 12, 14, 16 };

const Double masses[] = { 1.0078246, 2.0141021, 3.01603, 4.00260,
                               6.015121, 7.016003, 9.012182, 10.012937,
                               11.009305, 12.0000000, 13.0033554, 14.0030732,
                               15.0001088, 15.9949141, 16.9991322, 17.9991616,
                               18.9984032, 19.992435, 20.993843, 21.991383,
                               22.989767, 23.985042, 24.985837, 25.982593,
                               26.981539, 27.976927, 28.976495, 29.973770,
                               30.973762, 31.972070, 32.971456, 33.967866,
                               35.967080, 34.9688531, 36.9659034, 35.967545,
                               37.962732, 39.962384, 38.963707, 39.963999,
                               40.961825, 39.962591, 41.958618, 42.958766,
                               43.955480, 45.953689, 47.952533, 44.955910,
                               45.952629, 46.951764, 47.947947, 48.947871,
                               49.944792, 49.947161, 50.943962, 49.946046,
                               51.940509, 52.940651, 53.938882, 54.938047,
                               53.939612, 55.934939, 56.935396, 57.933277,
                               58.933198, 57.935346, 59.930788, 60.931058,
                               61.928346, 63.927968, 62.939598, 64.927793,
                               63.929145, 65.926034, 66.927129, 67.924846,
                               69.925325, 68.925580, 70.924700, 69.924250,
                               71.922079, 72.923463, 73.921177, 75.921401,
                               74.921594, 73.922475, 75.919212, 76.919912,
                               77.9190, 79.916520, 81.916698, 78.918336,
                               80.916289, 77.914, 79.916380, 81.913482,
                               82.914135, 83.911507, 85.910616, 84.911794,
                               86.909187, 83.913430, 85.909267, 86.908884,
                               87.905619, 88.905849, 89.904703, 90.905644,
                               91.905039, 93.906314, 95.908275, 92.906377,
                               91.906808, 93.905085, 94.905840, 95.904678,
                               96.906020, 97.905406, 99.907477, 98.0,
                               95.907599, 97.905287, 98.905939, 99.904219,
                               100.905582, 101.904348, 103.905424, 102.905500,
                               101.905634, 103.904029, 104.905079, 105.903478,
                               107.903895, 109.905167, 106.905092, 108.904757,
                               105.906461, 107.904176, 109.903005, 110.904182,
                               111.902758, 112.904400, 113.903357, 115.904754,
                               112.904061, 114.903880, 111.904826, 113.902784,
                               114.903348, 115.901747, 116.902956, 117.901609,
                               118.903310, 119.902200, 121.903440, 123.905274,
                               120.903821, 122.904216, 119.904048, 121.903054,
                               122.904271, 123.902823, 124.904433, 125.903314,
                               127.904463, 129.906229, 126.904473, 123.905894,
                               125.904281, 127.903531, 128.904780, 129.903509,
                               130.905072, 131.904144, 133.905395, 135.907214,
                               132.905429, 129.906282, 131.905042, 133.904486,
                               134.905665, 135.904553, 136.905812, 137.905232,
                               137.90711, 138.906347, 135.907140, 137.905985,
                               139.905433, 141.909241, 140.907647, 141.907719,
                               142.909810, 143.910083, 144.912570, 145.913113,
                               147.916889, 149.920887, 145.0, 143.911998,
                               146.914895, 147.914820, 148.917181, 149.917273,
                               151.919729, 153.922206, 150.919847, 152.921225,
                               151.919786, 153.920861, 154.922618, 155.922118,
                               156.923956, 157.924099, 159.927049, 158.925342,
                               155.925277, 157.924403, 159.925193, 160.926930,
                               161.926795, 162.928728, 163.929171, 164.930319,
                               161.928775, 163.929198, 165.930290, 166.932046,
                               167.932368, 169.935461, 168.934212, 167.933894,
                               169.934759, 170.936323, 171.936378, 172.938208,
                               173.938859, 175.942564, 174.940770, 175.942679,
                               173.940044, 175.941406, 176.943217, 177.943696,
                               178.945812, 179.946545, 179.947462, 180.947992,
                               179.946701, 181.948202, 182.950220, 183.950928,
                               185.954357, 184.952951, 186.955744, 183.952488,
                               185.953830, 186.955741, 187.955860, 188.958137,
                               189.958436, 191.961467, 190.960584, 192.962917,
                               189.959917, 191.961019, 193.962655, 194.964766,
                               195.964926, 197.967869, 196.966543, 195.965807,
                               197.966743, 198.968254, 199.968300, 200.970277,
                               201.970617, 203.973467, 202.972320, 204.974401,
                               203.973020, 205.974440, 206.975872, 207.976627,
                               208.980374, 209.0, 210.0, 222.0, 223.0, 226.025,
                               227.028, 232.038054, 231.0359, 234.040946,
                               235.043924, 238.050784, 237.048, 244.0, 243.0,
                               247.0, 247.0, 251.0, 252.0, 257.0, 258.0, 259.0,
                               260.0, 2.0141021, 13.0033554, 15.0001088,
                               17.9991616 };

const Double frequencies[] = { 0.99985, 0.00015, 0.00000138, 0.99999862,
                                    0.075, 0.925, 1.0, 0.199, 0.801, 0.988930,
                                    0.011070, 0.996337, 0.003663, 0.997590,
                                    0.000374, 0.002036, 1.0, 0.9048, 0.0027,
                                    0.0925, 1.0, 0.7899, 0.1000, 0.1101, 1.0,
                                    0.9223, 0.0467, 0.0310, 1.0, 0.9502,
                                    0.0075, 0.0421, 0.0002, 0.755290, 0.244710,
                                    0.00337, 0.00063, 0.99600, 0.932581,
                                    0.000117, 0.067302, 0.96941, 0.00647,
                                    0.00135, 0.02086, 0.00004, 0.00187, 1.0,
                                    0.080, 0.073, 0.738, 0.055, 0.054, 0.00250,
                                    0.99750, 0.04345, 0.83790, 0.09500,
                                    0.02365, 1.0, 0.0590, 0.9172, 0.0210,
                                    0.0028, 1.0, 0.6827, 0.2610, 0.0113,
                                    0.0359, 0.0091, 0.6917, 0.3083, 0.486,
                                    0.279, 0.041, 0.188, 0.006, 0.60108,
                                    0.39892, 0.205, 0.274, 0.078, 0.365, 0.078,
                                    1.0, 0.009, 0.091, 0.076, 0.236, 0.499,
                                    0.089, 0.5069, 0.4931, 0.0035, 0.0225,
                                    0.116, 0.115, 0.570, 0.173, 0.7217, 0.2783,
                                    0.0056, 0.0986, 0.0700, 0.8258, 1.0,
                                    0.5145, 0.1122, 0.1715, 0.1738, 0.0280,
                                    1.0, 0.1484, 0.0925, 0.1592, 0.1668,
                                    0.0955, 0.2413, 0.0963, 1.0, 0.0554,
                                    0.0186, 0.127, 0.126, 0.171, 0.316, 0.186,
                                    1.0, 0.0102, 0.1114, 0.2233, 0.2733,
                                    0.2646, 0.1172, 0.51839, 0.48161, 0.0125,
                                    0.0089, 0.1249, 0.1280, 0.2413, 0.1222,
                                    0.2873, 0.0749, 0.043, 0.957, 0.0097,
                                    0.0065, 0.0036, 0.1453, 0.0768, 0.2422,
                                    0.0858, 0.3259, 0.0463, 0.0579, 0.574,
                                    0.426, 0.00095, 0.0259, 0.00905, 0.0479,
                                    0.0712, 0.1893, 0.3170, 0.3387, 1.0,
                                    0.0010, 0.0009, 0.0191, 0.264, 0.041,
                                    0.212, 0.269, 0.104, 0.089, 1.0, 0.00106,
                                    0.00101, 0.0242, 0.06593, 0.0785, 0.1123,
                                    0.7170, 0.00090, 0.99910, 0.0019, 0.0025,
                                    0.8843, 0.1113, 1.0, 0.2713, 0.1218,
                                    0.2380, 0.0830, 0.1719, 0.0576, 0.0564,
                                    1.0, 0.031, 0.150, 0.113, 0.138, 0.074,
                                    0.267, 0.227, 0.478, 0.522, 0.0020, 0.0218,
                                    0.1480, 0.2047, 0.1565, 0.2484, 0.2186,
                                    1.0, 0.0006, 0.0010, 0.0234, 0.189, 0.255,
                                    0.249, 0.282, 1.0, 0.0014, 0.0161, 0.336,
                                    0.2295, 0.268, 0.149, 1.0, 0.0013, 0.0305,
                                    0.143, 0.219, 0.1612, 0.318, 0.127, 0.9741,
                                    0.0259, 0.00162, 0.05206, 0.18606, 0.27297,
                                    0.13629, 0.35100, 0.00012, 0.99988, 0.0012,
                                    0.263, 0.1428, 0.307, 0.286, 0.3740,
                                    0.6260, 0.0002, 0.0158, 0.016, 0.133,
                                    0.161, 0.264, 0.410, 0.373, 0.627, 0.0001,
                                    0.0079, 0.329, 0.338, 0.253, 0.072, 1.0,
                                    0.0015, 0.100, 0.169, 0.231, 0.132, 0.298,
                                    0.0685, 0.29524, 0.70476, 0.014, 0.241,
                                    0.221, 0.524, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                    1.0, 1.0, 1.0, 0.000055, 0.00720, 0.992745,
                                    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };

ElementImpl::ElementImplKeyType ElementImpl::freeId = nEntries + 1;
std::map<ElementImpl::ElementImplSymbolType, ElementImpl::ElementImplKeyType> ElementImpl::elementMapping;

ElementImpl::ElementImplKeyType ElementImpl::getNextId()
{
    return freeId++;
}

const std::map<ElementImpl::ElementImplSymbolType,
        ElementImpl::ElementImplKeyType>& ElementImpl::getDefaultMapping()
{
    if (elementMapping.empty()) {
        for (Size i = 0; i < nEntries; ++i) {
            elementMapping.insert(
                std::pair<ElementImpl::ElementImplSymbolType,
                        ElementImpl::ElementImplKeyType>(symbols[i], i + 1));
        }
    }
    return elementMapping;
}

ElementImpl::ElementImplKeyType ElementImpl::getDefaultKeyForElementSymbol(
    const ElementImplSymbolType& symbol)
{
    const std::map<ElementImplSymbolType, ElementImplKeyType>& map =
            getDefaultMapping();
    std::map<ElementImplSymbolType, ElementImplKeyType>::const_iterator it =
            map.find(symbol);
    if (it != map.end()) {
        return it->second;
    }
    throw mstk::LogicError("ElementImpl::getDefaultKeyForElementSymbol(): Cannot find given symbol in default element mapping.");
}

ElementImpl::ElementImpl(const ElementImpl::ElementImplKeyType& id) :
        id_(id), isotopes_()
{
    ElementImpl::ElementImplKeyType idc = id - 1;
    if (idc >= 0 && idc < nEntries) {
        symbol_ = symbols[idc];
        atomicNumber_ = atomicNumbers[idc];
        //        symbol_ = symbols[id];
        //        atomicNumber_ = atomicNumbers[id];
        Size k = 0;
        for (Size i = 0; i < idc; ++i) {
            k += nIsotopes[i];
        }
        for (Size j = 0; j < nIsotopes[idc]; ++j, ++k) {
            isotopes_.push_back(Isotope(masses[k], frequencies[k]));
        }
    } else {
    	throw mstk::LogicError("Element::Element(): Unkown element id.");
    }
}

ElementImpl::ElementImpl(const ElementImplKeyType& id, const String& symbol,
    const Size& atomicNumber) :
        id_(id), symbol_(symbol), atomicNumber_(atomicNumber), isotopes_()
{
}

ElementImpl::ElementImpl(const ElementImplKeyType& id, const String& symbol,
    const Size& atomicNumber, const std::vector<Isotope>& isotopes) :
        id_(id), symbol_(symbol), atomicNumber_(atomicNumber), isotopes_(
            isotopes)
{
}

const String& ElementImpl::getSymbol() const
{
    return symbol_;
}

const Size& ElementImpl::getAtomicNumber() const
{
    return atomicNumber_;
}

const std::vector<Isotope>& ElementImpl::getIsotopes() const
{
    return isotopes_;
}

void ElementImpl::addIsotope(const Isotope& i)
{
    isotopes_.push_back(i);
}

void ElementImpl::addIsotope(const Double& mass, const Double& frequency)
{
    addIsotope(Isotope(mass, frequency));
}

void ElementImpl::clearIsotopes()
{
    isotopes_.clear();
}

void ElementImpl::setIsotopes(const std::vector<Isotope>& isotopes)
{
    isotopes_ = isotopes;
}

Size ElementImpl::getNumberOfStandardElements()
{
    return nEntries;
}

bool ElementImpl::operator==(const ElementImpl& e) const
{
    return e.id_ == id_ && e.symbol_ == symbol_
            && e.atomicNumber_ == atomicNumber_ && e.isotopes_ == isotopes_;
}

bool ElementImpl::operator!=(const ElementImpl& e) const
{
    return !(operator ==(e));
}

ElementImpl& ElementImpl::operator=(const ElementImpl& rhs)
{
    if (this != &rhs) {
        id_ = rhs.id_;
        symbol_ = rhs.symbol_;
        atomicNumber_ = rhs.atomicNumber_;
        isotopes_ = rhs.isotopes_;
    }
    return *this;
}

std::ostream& operator<<(std::ostream& os, const ElementImpl& o)
{
    os << o.getId() << ":" << o.getSymbol() << "(" << o.getAtomicNumber()
            << ") - " << o.getIsotopes();
    return os;
}

} // namespace elements
} // namespace aas
} // namespace mstk
