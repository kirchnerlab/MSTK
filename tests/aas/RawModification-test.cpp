/*
 * RawModification-test.cpp
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

#include <MSTK/aas/RawModification.hpp>
#include <MSTK/aas/RawAminoAcid.hpp>
#include <MSTK/common/Error.hpp>

#include "unittest.hxx"

#include <iostream>
#include <vector>

using namespace mstk;
using namespace aas;
using namespace aas::modifications;
using namespace aas::aminoAcids;
using namespace aas::stoichiometries;

/** Short description.
 * Long description.
 */
struct RawModificationTestSuite : vigra::test_suite
{
    /** Constructor.
     * The RawModificationTestSuite constructor adds all RawModification tests to
     * the test suite. If you write an additional test, add the test
     * case here.
     */
    RawModificationTestSuite() :
            vigra::test_suite("RawModification")
    {
        add(testCase(&RawModificationTestSuite::testRawModification));
        add(testCase(&RawModificationTestSuite::testRawModificationFW));
        add(testCase(&RawModificationTestSuite::testRawModificationRef));
        add(testCase(&RawModificationTestSuite::testAddRawModification));
        add(testCase(&RawModificationTestSuite::testAddRawModificationRef));
        add(
            testCase(&RawModificationTestSuite::testOverrideUninitializedRawModification));
        add(
            testCase(&RawModificationTestSuite::testOverrideInitializedRawModification));
        add(testCase(&RawModificationTestSuite::testCreateModification));
    }

    void testRawModification()
    {
        RawModificationImpl::RawModificationImplKeyType k = "Deamidated";
        RawModificationImpl deamidated(k);

        shouldEqual(deamidated.getId(), k);
        shouldEqual(deamidated.getName(), "Deamidated");
        shouldEqual(deamidated.getFullName(), "Deamidation");
        shouldEqual(deamidated.isVerified(), false);
        std::vector<String> altNames = deamidated.getAltNames();
        size_t altNameSize = 2;
        shouldEqual(altNames.size(), altNameSize);
        shouldEqual(altNames[0], "phenyllactyl from N-term Phe");
        shouldEqual(altNames[1], "Citrullination");

        aas::elements::Element H(1);
        aas::elements::Element N(7);
        aas::elements::Element O(8);
        Stoichiometry s = deamidated.getStoichiometry(), expected_s;
        expected_s.set(H, -1);
        expected_s.set(N, -1);
        expected_s.set(O, 1);
        shouldEqual(s, expected_s);

        const std::vector<Specificity>& specificities =
                deamidated.getSpecificities();
        size_t numberOfSpecificities = 4;
        shouldEqual(specificities.size(), numberOfSpecificities);
        Specificity spec0(RawAminoAcid('Q'), Specificity::ANYWHERE,
            Specificity::ARTEFACT);
        Specificity spec1(RawAminoAcid('R'), Specificity::ANYWHERE,
            Specificity::POST_TRANSLATIONAL);
        Specificity spec2(RawAminoAcid('N'), Specificity::ANYWHERE,
            Specificity::ARTEFACT);
        Specificity spec3(RawAminoAcid('F'), Specificity::PROTEIN_N_TERM,
            Specificity::POST_TRANSLATIONAL);
        shouldEqual(specificities[0], spec0);
        shouldEqual(specificities[1], spec1);
        shouldEqual(specificities[2], spec2);
        shouldEqual(specificities[3], spec3);

        RawModificationImpl test("id", "name", "fullname", false);
        String name = "test name";
        String fullName = "test full name";
        Bool verified = false;
        Stoichiometry st;
        st.set(H, 10);
        st.set(N, 5);
        st.set(O, 2);
        std::vector<Specificity> specs;
        specs.push_back(
            Specificity(RawAminoAcid('A'), Specificity::ANYWHERE,
                Specificity::CHEMICAL_DERIVATIVE));
        std::vector<String> alts;
        alts.push_back("Alt name");

        test.setName(name);
        shouldEqual(test.getName(), name);
        test.setFullName(fullName);
        shouldEqual(test.getName(), name);
        test.setVerified(verified);
        shouldEqual(test.isVerified(), verified);
        test.setStoichiometry(st);
        shouldEqual(test.getStoichiometry(), st);
        test.setSpecificities(specs);
        test.addSpecificity(
            Specificity(RawAminoAcid('C'), Specificity::ANY_C_TERM,
                Specificity::POST_TRANSLATIONAL));
        specs.push_back(
            Specificity(RawAminoAcid('C'), Specificity::ANY_C_TERM,
                Specificity::POST_TRANSLATIONAL));
        shouldEqual(test.getSpecificities(), specs);
        test.setAltNames(alts);
        shouldEqual(test.getAltNames().size(), alts.size());
        shouldEqual(test.getAltNames()[0], alts[0]);

        RawModificationImpl heavy("Propionyl:13C(3)");
        Stoichiometry esheavy;
        esheavy.set(
            elements::Element(
                elements::ElementImpl::getDefaultKeyForElementSymbol("13C")),
            3);
        esheavy.set(
            elements::Element(
                elements::ElementImpl::getDefaultKeyForElementSymbol("O")), 1);
        esheavy.set(
            elements::Element(
                elements::ElementImpl::getDefaultKeyForElementSymbol("H")), 4);
        shouldEqual(heavy.getStoichiometry(), esheavy);

        RawModificationImpl tmp1 = test;
        tmp1 = deamidated;
        shouldEqual(tmp1, deamidated);
        shouldEqual(tmp1 == test, false);
        shouldEqual(tmp1 != test, true);
        shouldEqual(tmp1 != deamidated, false);
    }

    void testRawModificationFW()
    {
        RawModification m1("ESP");
        RawModification m2("Oxidation");
        RawModification m3("TMT");

        shouldEqual(m1 < m1, false);
        shouldEqual(m1 < m2, true);
        shouldEqual(m1 < m3, true);

        shouldEqual(m1 <= m1, true);
        shouldEqual(m1 <= m2, true);
        shouldEqual(m1 <= m3, true);

        shouldEqual(m1 > m1, false);
        shouldEqual(m1 > m2, false);
        shouldEqual(m1 > m3, false);

        shouldEqual(m1 >= m1, true);
        shouldEqual(m1 >= m2, false);
        shouldEqual(m1 >= m3, false);
    }

    void testRawModificationRef()
    {
        RawModificationImpl::RawModificationImplKeyType k1 = "Carbamyl", k2 =
                "Carboxymethyl";

        // create two different raw modifications without using the factory
        RawModificationImpl e_2(k1);
        RawModificationImpl e_3(k2);

        // create the same raw modifications again, but using the flyweight implementation
        RawModification er_2(k1);
        RawModification er_3(k2);

        // testing equality
        shouldEqual(er_2, e_2);
        shouldEqual(er_3, e_3);

        // testing whether the id extractor works correctly
        shouldEqual(er_2.get_key(), er_2.get().getId());
        shouldEqual(er_3.get_key(), er_3.get().getId());

        // retrieve the amino acids again
        RawModification er_t2(k1);
        RawModification er_t3(k2);

        // test whether the old references and new references point to the same
        // object
        shouldEqual(&er_2.get(), &er_t2.get());
        shouldEqual(&er_3.get(), &er_t3.get());
    }

    void testAddRawModification()
    {
        // setting up test data
        aas::elements::Element H(1);
        aas::elements::Element C(6);
        aas::elements::Element O(8);

        // create raw modification
        RawModificationImpl::RawModificationImplKeyType k = "customID";
        String name = "Name";
        String name2 = "Othername";
        String fullName = "fullName";
        Bool verified = false;
        std::vector<String> altNames;
        Stoichiometry stoichiometry;
        std::vector<Specificity> specificities;

        // adding not existing element
        shouldEqual(
            addRawModification(k, name, fullName, altNames, stoichiometry, specificities, verified),
            true);
        RawModificationImpl e(k, name, fullName, verified);
        e.setAltNames(altNames);
        e.setStoichiometry(stoichiometry);
        e.setSpecificities(specificities);

        RawModification e_r(k);
        shouldEqual(e_r, e);

        // a second try to add a raw modification with the same id but different properties should fail
        shouldEqual(
            addRawModification(k, name2, fullName, altNames, stoichiometry, specificities, verified),
            false);
        // and the raw modification added in the first place should stay the same
        shouldEqual(RawModification(k), e);
        RawModificationImpl e2(k, name2, fullName, verified);
        e2.setAltNames(altNames);
        e2.setStoichiometry(stoichiometry);
        e2.setSpecificities(specificities);
        shouldEqual(RawModification(k) == e2, false);
    }

    void testAddRawModificationRef()
    {
        // first test should result in an out_of_range exceptions since the
        // default raw modification table does no contain a raw modification with id "unkown"
        bool thrown = false;
        try {
            RawModification test("unkown");
        } catch (mstk::LogicError& e) {
            thrown = true;
        };

        shouldEqual(thrown, true);

        // create an arbitrary raw modification
        RawModificationImpl::RawModificationImplKeyType k1 = "unkown";
        String name = "Name";
        String fullName = "fullName";
        Bool verified = false;
        RawModificationImpl t(k1, name, fullName, verified);
        // add raw modification to the flyweight table
        RawModification tr(t);

        // test equality of the original amino acid and the const ref
        shouldEqual(tr, t);
        shouldEqual(t.getId(), tr.get_key());

        // retrieve the raw modification directly from the flyweight table
        RawModification tr_t(k1);

        // test equality of two const refs of the same raw modification
        shouldEqual(tr, tr_t);
        // test reference pointer to ensure it is the same object
        shouldEqual(&tr.get(), &tr_t.get());
    }

    void testOverrideUninitializedRawModification()
    {
        // create different raw modification "Deamidated"
        RawModificationImpl::RawModificationImplKeyType k1 = "Phospho";
        String name = "Name";
        String fullName = "fullName";
        Bool verified = false;
        RawModificationImpl t(k1, name, fullName, verified);
        // store it in flyweight table
        RawModification tr(t);
        // retrieve element id k1
        RawModification tr_t(k1);

        // getSymbol returns "Name" instead of "Phospho"
        shouldEqual(tr_t.get().getName(), name);
        shouldEqual(tr_t.get().getName() != RawModificationImpl(k1).getName(),
            true);

    }

    void testOverrideInitializedRawModification()
    {
        RawModificationImpl::RawModificationImplKeyType k1 = "Phospho";
        RawModification tr_1(k1);

        RawModificationImpl t(k1, "other name", "full name", false);
        // flyweight checks the id of the raw modification t and recognizes it, since it
        // was retrieved(initialized) earlier
        // as a consequence, the flyweight factory returns the known object
        // instead of overriding it with the given one!
        RawModification tr_2(t);

        shouldEqual(tr_2.get().getName() != "other name", true);
    }

    void testCreateModification()
    {
        // create all modifications
        RawModificationImpl::RawModificationImplKeyType names[] = {
                "Acetyl", "Amidated", "Biotin", "Carbamidomethyl", "Carbamyl",
                "Carboxymethyl", "Deamidated", "ICAT-G", "ICAT-G:2H(8)",
                "Met->Hse", "Met->Hsl", "ICAT-D:2H(8)", "ICAT-D", "NIPCAM",
                "PEO-Iodoacetyl-LC-Biotin", "Phospho", "Dehydrated",
                "Propionamide", "Pyridylacetyl", "Pyro-carbamidomethyl",
                "Glu->pyro-Glu", "Gln->pyro-Glu", "SMA", "Cation:Na",
                "Pyridylethyl", "Methyl", "Oxidation", "Dimethyl", "Trimethyl",
                "Methylthio", "Sulfo", "Hex", "Lipoyl", "HexNAc", "Farnesyl",
                "Myristoyl", "PyridoxalPhosphate", "Palmitoyl",
                "GeranylGeranyl", "Phosphopantetheine", "FAD", "Tripalmitate",
                "Guanidinyl", "HNE", "Glucuronyl", "Glutathione",
                "Acetyl:2H(3)", "Propionyl", "Propionyl:13C(3)", "GIST-Quat",
                "GIST-Quat:2H(3)", "GIST-Quat:2H(6)", "GIST-Quat:2H(9)",
                "Succinyl", "Succinyl:2H(4)", "Succinyl:13C(4)",
                "probiotinhydrazide", "Pro->pyro-Glu", "His->Asn", "His->Asp",
                "Trp->Hydroxykynurenin", "Delta:H(4)C(3)", "Delta:H(4)C(2)",
                "Cys->Dha", "Arg->GluSA", "Trioxidation", "Iminobiotin", "ESP",
                "ESP:2H(10)", "NHS-LC-Biotin", "EDT-maleimide-PEO-biotin",
                "IMID", "IMID:2H(4)", "Lysbiotinhydrazide",
                "Propionamide:2H(3)", "Nitro", "ICAT-C", "Delta:H(2)C(2)",
                "Trp->Kynurenin", "Lys->Allysine", "ICAT-C:13C(9)",
                "FormylMet", "Nethylmaleimide", "OxLysBiotinRed", "IBTP",
                "OxLysBiotin", "OxProBiotinRed", "OxProBiotin", "OxArgBiotin",
                "OxArgBiotinRed", "EDT-iodoacetyl-PEO-biotin", "GlyGly",
                "Formyl", "ICAT-H", "ICAT-H:13C(6)", "Cation:K", "Thioacyl",
                "Fluoro", "Fluorescein", "Iodo", "Diiodo", "Triiodo",
                "Myristoleyl", "Pro->Pyrrolidinone", "Myristoyl+Delta:H(-4)",
                "Benzoyl", "Hex(5)HexNAc(2)", "Dansyl", "a-type-ion",
                "Amidine", "HexNAc(1)dHex(1)", "HexNAc(2)", "Hex(3)",
                "HexNAc(1)dHex(2)", "Hex(1)HexNAc(1)dHex(1)",
                "HexNAc(2)dHex(1)", "Hex(1)HexNAc(2)",
                "Hex(1)HexNAc(1)NeuAc(1)", "HexNAc(2)dHex(2)",
                "Hex(1)HexNAc(2)Pent(1)", "Hex(1)HexNAc(2)dHex(1)",
                "Hex(2)HexNAc(2)", "Hex(3)HexNAc(1)Pent(1)",
                "Hex(1)HexNAc(2)dHex(1)Pent(1)", "Hex(1)HexNAc(2)dHex(2)",
                "Hex(2)HexNAc(2)Pent(1)", "Hex(2)HexNAc(2)dHex(1)",
                "Hex(3)HexNAc(2)", "Hex(1)HexNAc(1)NeuAc(2)",
                "Hex(3)HexNAc(2)P(1)", "Delta:S(-1)Se(1)", "NBS:13C(6)",
                "Methyl:2H(3)13C(1)", "Dimethyl:2H(6)13C(2)", "NBS",
                "Delta:H(1)O(-1)18O(1)", "QAT", "BHT",
                "Delta:H(4)C(2)O(-1)S(1)", "DAET", "Pro->Pyrrolidone",
                "Label:13C(9)", "Label:13C(9)+Phospho", "Label:13C(6)", "HPG",
                "2HPG", "QAT:2H(3)", "Label:18O(2)", "AccQTag",
                "Dimethyl:2H(4)", "EQAT", "EQAT:2H(5)", "Ethanedithiol",
                "NEIAA:2H(5)", "Delta:H(6)C(6)O(1)", "Delta:H(4)C(3)O(1)",
                "Delta:H(2)C(3)", "Delta:H(4)C(6)", "Delta:H(8)C(6)O(2)",
                "ADP-Ribosyl", "NEIAA", "iTRAQ4plex", "Crotonaldehyde",
                "Bromo", "Amino", "Argbiotinhydrazide", "Label:18O(1)",
                "Label:13C(6)15N(2)", "Thiophospho", "SPITC", "IGBP",
                "Cytopiloyne", "Cytopiloyne+water", "Label:13C(6)15N(4)",
                "Label:13C(9)15N(1)", "Label:2H(3)", "Label:13C(5)15N(1)",
                "PET", "CAF", "Xlink:SSD", "Nitrosyl", "AEBS", "Ethanolyl",
                "Label:13C(6)15N(2)+Dimethyl", "HMVK", "Ethyl", "CoenzymeA",
                "Methyl+Deamidated", "Delta:H(5)C(2)", "Methyl:2H(2)",
                "SulfanilicAcid", "SulfanilicAcid:13C(6)", "Biotin-PEO-Amine",
                "Trp->Oxolactone", "Biotin-HPDP", "Delta:Hg(1)", "IodoU-AMP",
                "CAMthiopropanoyl", "IED-Biotin", "dHex", "Methyl:2H(3)",
                "Carboxy", "Bromobimane", "Menadione", "DeStreak",
                "dHex(1)Hex(3)HexNAc(4)", "dHex(1)Hex(4)HexNAc(4)",
                "dHex(1)Hex(5)HexNAc(4)", "Hex(3)HexNAc(4)", "Hex(4)HexNAc(4)",
                "Hex(5)HexNAc(4)", "Cysteinyl", "Lys-loss", "Nmethylmaleimide",
                "CyDye-Cy3", "DimethylpyrroleAdduct", "Delta:H(2)C(5)",
                "Delta:H(2)C(3)O(1)", "Nethylmaleimide+water",
                "Methyl+Acetyl:2H(3)", "Xlink:B10621", "DTBP", "FP-Biotin",
                "Thiophos-S-S-biotin", "Can-FP-biotin", "HNE+Delta:H(2)",
                "Thrbiotinhydrazide", "Methylamine", "Diisopropylphosphate",
                "Isopropylphospho", "ICPL:13C(6)", "CarbamidomethylDTT",
                "ICPL", "Deamidated:18O(1)", "Arg->Orn", "Cation:Cu[I]",
                "Dehydro", "Diphthamide", "Hydroxyfarnesyl", "Diacylglycerol",
                "Carboxyethyl", "Hypusine", "Retinylidene",
                "Lys->AminoadipicAcid", "Cys->PyruvicAcid", "Ammonia-loss",
                "Phycocyanobilin", "Phycoerythrobilin", "Phytochromobilin",
                "Heme", "Molybdopterin", "Quinone", "Glucosylgalactosyl",
                "GPIanchor", "PhosphoribosyldephosphoCoA", "GlycerylPE",
                "Triiodothyronine", "Thyroxine", "Tyr->Dha", "Didehydro",
                "Cys->Oxoalanine", "Ser->LacticAcid", "GluGlu",
                "Phosphoadenosine", "Glu", "Hydroxycinnamyl", "Glycosyl",
                "FMNH", "Archaeol", "Phenylisocyanate",
                "Phenylisocyanate:2H(5)", "Phosphoguanosine", "Hydroxymethyl",
                "MolybdopterinGD+Delta:S(-1)Se(1)", "Dipyrrolylmethanemethyl",
                "PhosphoUridine", "Glycerophospho", "Carboxy->Thiocarboxy",
                "Sulfide", "PyruvicAcidIminyl", "Delta:Se(1)",
                "MolybdopterinGD", "Dioxidation", "Octanoyl", "PhosphoHexNAc",
                "PhosphoHex", "Palmitoleyl", "Cholesterol",
                "Didehydroretinylidene", "CHDH", "Methylpyrroline",
                "Hydroxyheme", "MicrocinC7", "Cyano", "Diironsubcluster",
                "Amidino", "FMN", "FMNC", "CuSMo", "Hydroxytrimethyl", "Deoxy",
                "Microcin", "Decanoyl", "GluGluGlu", "GluGluGluGlu", "HexN",
                "Xlink:DMP-s", "Xlink:DMP", "NDA", "SPITC:13C(6)",
                "TMAB:2H(9)", "TMAB", "FTC", "AEC-MAEC", "BADGE",
                "Label:2H(4)", "Hep", "CyDye-Cy5", "DHP", "BHTOH",
                "IGBP:13C(2)", "Nmethylmaleimide+water", "PyMIC",
                "LG-lactam-K", "BisANS", "Piperidine", "Diethyl",
                "LG-Hlactam-K", "Dimethyl:2H(4)13C(2)", "C8-QAT", "Hex(2)",
                "LG-lactam-R", "Withaferin", "Biotin:Thermo-88317",
                "CLIP_TRAQ_2", "LG-Hlactam-R", "Maleimide-PEO2-Biotin",
                "Sulfo-NHS-LC-LC-Biotin", "FNEM", "PropylNAGthiazoline",
                "Dethiomethyl", "iTRAQ4plex114", "iTRAQ4plex115", "Dibromo",
                "LeuArgGlyGly", "CLIP_TRAQ_3", "CLIP_TRAQ_4",
                "Biotin:Cayman-10141", "Biotin:Cayman-10013", "Ala->Ser",
                "Ala->Thr", "Ala->Asp", "Ala->Pro", "Ala->Gly", "Ala->Glu",
                "Ala->Val", "Cys->Phe", "Cys->Ser", "Cys->Trp", "Cys->Tyr",
                "Cys->Arg", "Cys->Gly", "Asp->Ala", "Asp->His", "Asp->Asn",
                "Asp->Gly", "Asp->Tyr", "Asp->Glu", "Asp->Val", "Glu->Ala",
                "Glu->Gln", "Glu->Asp", "Glu->Lys", "Glu->Gly", "Glu->Val",
                "Phe->Ser", "Phe->Cys", "Phe->Xle", "Phe->Tyr", "Phe->Val",
                "Gly->Ala", "Gly->Ser", "Gly->Trp", "Gly->Glu", "Gly->Val",
                "Gly->Asp", "Gly->Cys", "Gly->Arg", "dNIC", "His->Pro",
                "His->Tyr", "His->Gln", "NIC", "His->Arg", "His->Xle",
                "Xle->Ala", "Xle->Thr", "Xle->Asn", "Xle->Lys", "Lys->Thr",
                "Lys->Asn", "Lys->Glu", "Lys->Gln", "Lys->Met", "Lys->Arg",
                "Lys->Xle", "Xle->Ser", "Xle->Phe", "Xle->Trp", "Xle->Pro",
                "Xle->Val", "Xle->His", "Xle->Gln", "Xle->Met", "Xle->Arg",
                "Met->Thr", "Met->Arg", "Met->Lys", "Met->Xle", "Met->Val",
                "Asn->Ser", "Asn->Thr", "Asn->Lys", "Asn->Tyr", "Asn->His",
                "Asn->Asp", "Asn->Xle", "Pro->Ser", "Pro->Ala", "Pro->His",
                "Pro->Gln", "Pro->Thr", "Pro->Arg", "Pro->Xle", "Gln->Pro",
                "Gln->Lys", "Gln->Glu", "Gln->His", "Gln->Arg", "Gln->Xle",
                "Arg->Ser", "Arg->Trp", "Arg->Thr", "Arg->Pro", "Arg->Lys",
                "Arg->His", "Arg->Gln", "Arg->Met", "Arg->Cys", "Arg->Xle",
                "Arg->Gly", "Ser->Phe", "Ser->Ala", "Ser->Trp", "Ser->Thr",
                "Ser->Asn", "Ser->Pro", "Ser->Tyr", "Ser->Cys", "Ser->Arg",
                "Ser->Xle", "Ser->Gly", "Thr->Ser", "Thr->Ala", "Thr->Asn",
                "Thr->Lys", "Thr->Pro", "Thr->Met", "Thr->Xle", "Thr->Arg",
                "Val->Phe", "Val->Ala", "Val->Glu", "Val->Met", "Val->Asp",
                "Val->Xle", "Val->Gly", "Trp->Ser", "Trp->Cys", "Trp->Arg",
                "Trp->Gly", "Trp->Xle", "Tyr->Phe", "Tyr->Ser", "Tyr->Asn",
                "Tyr->His", "Tyr->Asp", "Tyr->Cys", "BDMAPP", "NA-LNO2",
                "NA-OA-NO2", "ICPL:2H(4)", "CarboxymethylDTT", "iTRAQ8plex",
                "Label:13C(6)15N(1)", "Label:2H(9)13C(6)15N(2)",
                "HNE-Delta:H(2)O", "4-ONE", "O-Dimethylphosphate",
                "O-Methylphosphate", "Diethylphosphate", "Ethylphosphate",
                "O-pinacolylmethylphosphonate", "Methylphosphonate",
                "O-Isopropylmethylphosphonate", "iTRAQ8plex:13C(6)15N(2)",
                "DTT_ST", "Ethanolamine", "TMT6plex", "DTT_C", "TMT2plex",
                "TMT", "ExacTagThiol", "ExacTagAmine", "NO_SMX_SEMD",
                "4-ONE+Delta:H(-2)O(-1)", "NO_SMX_SMCT", "NO_SMX_SIMD",
                "Malonyl", "3sulfo", "trifluoro", "TNBS", "Biotin-phenacyl",
                "DTT_C:2H(6)", "lapachenole", "Label:13C(5)", "maleimide",
                "IDEnT", "DTT_ST:2H(6)", "Met-loss", "Met-loss+Acetyl",
                "Menadione-HQ", "Carboxymethyl:13C(2)", "NEM:2H(5)",
                "Gly-loss+Amide", "TMPP-Ac", "Label:13C(6)+GlyGly", "Arg->Npo",
                "Label:2H(4)+Acetyl", "Pentylamine", "Biotin:Thermo-21345",
                "Dihydroxyimidazolidine", "DFDNB", "Cy3b-maleimide",
                "Hex1HexNAc1", "AEC-MAEC:2H(4)", "BMOE", "Biotin:Thermo-21360",
                "Label:13C(6)+Acetyl", "Label:13C(6)15N(2)+Acetyl", "EQIGG",
                "cGMP", "cGMP+RMP-loss", "mTRAQ", "Arg2PG",
                "Label:2H(4)+GlyGly", "Label:13C(8)15N(2)",
                "Label:13C(1)2H(3)", "ZGB", "MG-H1", "G-H1",
                "Label:13C(6)15N(2)+GlyGly", "ICPL:13C(6)2H(4)",
                "DyLight-maleimide", "mTRAQ:13C(3)15N(1)",
                "Methyl-PEO12-Maleimide", "MDCC", "QQQTGG", "QEQTGG",
                "HydroxymethylOP", "Biotin:Thermo-21325",
                "Label:13C(1)2H(3)+Oxidation", "Bodipy", "Biotin-PEG-PRA",
                "Met->Aha", "Label:15N(4)", "pyrophospho", "Met->Hpg",
                "4AcAllylGal", "DimethylArsino", "Lys->CamCys", "Phe->CamCys",
                "Leu->MetOx", "Lys->MetOx", "Galactosyl", "SMCC-maleimide",
                "Bacillosamine", "MTSL", "HNE-BAHAH", "Ethoxyformyl",
                "Methylmalonylation", "AROD", "Cys->methylaminoAla",
                "Cys->ethylaminoAla", "Label:13C(4)15N(2)+GlyGly",
                "ethylamino", "MercaptoEthanol", "Atto495Maleimide",
                "AMTzHexNAc2", "EthylAmide", "VFQQQTGG", "VIEVYQEQTGG",
                "Chlorination", "dichlorination", "DNPS", "SulfoGMBS",
                "DimethylamineGMBS", "Label:15N(2)2H(9)", "LG-anhydrolactam",
                "LG-pyrrole", "LG-anhyropyrrole", "3-deoxyglucosone",
                "Cation:Li", "Cation:Ca[II]", "Cation:Fe[II]", "Cation:Ni[II]",
                "Cation:Zn[II]", "Cation:Ag", "Cation:Mg[II]", "2-succinyl",
                "Propargylamine", "Phosphopropargyl", "SUMO2135", "SUMO3549",
                "Chlorpyrifos", "BITC", "Carbofuran", "PEITC", "thioacylPA",
                "maleimide3", "maleimide5", "Puromycin", "glucosone",
                "Label:13C(6)+Dimethyl", "cysTMT", "cysTMT6plex",
                "ISD_z+2_ion", "Ammonium", "BHAc", "Biotin:Sigma-B1267",
                "Label:15N(1)", "Label:15N(2)", "Label:15N(3)", "sulfo+amino",
                "AHA-Alkyne", "AHA-Alkyne-KDDDD", "EGCG1", "EGCG2",
                "Label:13C(6)15N(4)+Methyl", "Label:13C(6)15N(4)+Dimethyl",
                "Label:13C(6)15N(4)+Methyl:2H(3)13C(1)",
                "Label:13C(6)15N(4)+Dimethyl:2H(6)13C(2)",
                "SecCarbamidomethyl", "Thiazolidine", "DEDGFLYMVYASQETFG",
                "Biotin:Invitrogen-M1602", "Xlink:DSS", "DMPO", "glycidamide",
                "Ahx2+Hsl", "ICDID", "ICDID:2H(6)", "Xlink:EGS", "Xlink:DST",
                "Xlink:DTSSP", "Xlink:SMCC", "2-nitrobenzyl", "Xlink:DMP-de",
                "Xlink:EGScleaved", "SecNEM", "SecNEM:2H(5)", "Thiadiazole",
                "Biotin:Thermo-88310", "TAMRA-FP", "Biotin:Thermo-21901+H2O",
                "Deoxyhypusine", "Acetyldeoxyhypusine", "Acetylhypusine",
                "Ala->Cys", "Ala->Phe", "Ala->His", "Ala->Xle", "Ala->Lys",
                "Ala->Met", "Ala->Asn", "Ala->Gln", "Ala->Arg", "Ala->Trp",
                "Ala->Tyr", "Cys->Ala", "Cys->Asp", "Cys->Glu", "Cys->His",
                "Cys->Xle", "Cys->Lys", "Cys->Met", "Cys->Asn", "Cys->Pro",
                "Cys->Gln", "Cys->Thr", "Cys->Val", "Asp->Cys", "Asp->Phe",
                "Asp->Xle", "Asp->Lys", "Asp->Met", "Asp->Pro", "Asp->Gln",
                "Asp->Arg", "Asp->Ser", "Asp->Thr", "Asp->Trp", "Glu->Cys",
                "Glu->Phe", "Glu->His", "Glu->Xle", "Glu->Met", "Glu->Asn",
                "Glu->Pro", "Glu->Arg", "Glu->Ser", "Glu->Thr", "Glu->Trp",
                "Glu->Tyr", "Phe->Ala", "Phe->Asp", "Phe->Glu", "Phe->Gly",
                "Phe->His", "Phe->Lys", "Phe->Met", "Phe->Asn", "Phe->Pro",
                "Phe->Gln", "Phe->Arg", "Phe->Thr", "Phe->Trp", "Gly->Phe",
                "Gly->His", "Gly->Xle", "Gly->Lys", "Gly->Met", "Gly->Asn",
                "Gly->Pro", "Gly->Gln", "Gly->Thr", "Gly->Tyr", "His->Ala",
                "His->Cys", "His->Glu", "His->Phe", "His->Gly", "His->Lys",
                "His->Met", "His->Ser", "His->Thr", "His->Val", "His->Trp",
                "Xle->Cys", "Xle->Asp", "Xle->Glu", "Xle->Gly", "Xle->Tyr",
                "Lys->Ala", "Lys->Cys", "Lys->Asp", "Lys->Phe", "Lys->Gly",
                "Lys->His", "Lys->Pro", "Lys->Ser", "Lys->Val", "Lys->Trp",
                "Lys->Tyr", "Met->Ala", "Met->Cys", "Met->Asp", "Met->Glu",
                "Met->Phe", "Met->Gly", "Met->His", "Met->Asn", "Met->Pro",
                "Met->Gln", "Met->Ser", "Met->Trp", "Met->Tyr", "Asn->Ala",
                "Asn->Cys", "Asn->Glu", "Asn->Phe", "Asn->Gly", "Asn->Met",
                "Asn->Pro", "Asn->Gln", "Asn->Arg", "Asn->Val", "Asn->Trp",
                "Pro->Cys", "Pro->Asp", "Pro->Glu", "Pro->Phe", "Pro->Gly",
                "Pro->Lys", "Pro->Met", "Pro->Asn", "Pro->Val", "Pro->Trp",
                "Pro->Tyr", "Gln->Ala", "Gln->Cys", "Gln->Asp", "Gln->Phe",
                "Gln->Gly", "Gln->Met", "Gln->Asn", "Gln->Ser", "Gln->Thr",
                "Gln->Val", "Gln->Trp", "Gln->Tyr", "Arg->Ala", "Arg->Asp",
                "Arg->Glu", "Arg->Asn", "Arg->Val", "Arg->Tyr", "Arg->Phe",
                "Ser->Asp", "Ser->Glu", "Ser->His", "Ser->Lys", "Ser->Met",
                "Ser->Gln", "Ser->Val", "Thr->Cys", "Thr->Asp", "Thr->Glu",
                "Thr->Phe", "Thr->Gly", "Thr->His", "Thr->Gln", "Thr->Val",
                "Thr->Trp", "Thr->Tyr", "Val->Cys", "Val->His", "Val->Lys",
                "Val->Asn", "Val->Pro", "Val->Gln", "Val->Arg", "Val->Ser",
                "Val->Thr", "Val->Trp", "Val->Tyr", "Trp->Ala", "Trp->Asp",
                "Trp->Glu", "Trp->Phe", "Trp->His", "Trp->Lys", "Trp->Met",
                "Trp->Asn", "Trp->Pro", "Trp->Gln", "Trp->Thr", "Trp->Val",
                "Trp->Tyr", "Tyr->Ala", "Tyr->Glu", "Tyr->Gly", "Tyr->Lys",
                "Tyr->Met", "Tyr->Pro", "Tyr->Gln", "Tyr->Arg", "Tyr->Thr",
                "Tyr->Val", "Tyr->Trp", "Tyr->Xle", "AHA-SS", "AHA-SS_CAM",
                "Biotin:Thermo-33033", "Biotin:Thermo-33033-H",
                "2-monomethylsuccinyl", "Saligenin", "Cresylphosphate",
                "CresylSaligeninPhosphate", "Ub-Br2", "Ub-VME", "Ub-amide",
                "Ub-fluorescein", "2-dimethylsuccinyl", "Gly", "pupylation",
                "Label:13C(4)", "HCysteinyl", "Label:13C(4)+Oxidation",
                "UgiJoullie", "HCysThiolactone", "UgiJoullieProGly",
                "Bipyridine", "Furan", "Difuran", "BMP-piperidinol",
                "UgiJoullieProGlyProGly", "Arg-loss", "Arg", "IMEHex(2)NeuAc",
                "Butyryl", "Dicarbamidomethyl", "Dimethyl:2H(6)", "GGQ",
                "QTGG", "Label:13C(3)15N(1)", "Label:13C(3)",
                "Label:13C(4)15N(1)", "Label:2H(10)", "Label:2H(4)13C(1)",
                "Lys", "mTRAQ:13C(6)15N(2)", "NeuAc", "NeuGc", "Propyl",
                "Propyl:2H(6)", "Propiophenone" };
        Size n = 927;
        std::cout << std::endl;
        for (Size i = 0; i < n; ++i) {
            RawModification m(names[i]);
            shouldEqual(m.get_key(), names[i]);
        }
    }

};

/** The main function that runs the tests for class RawModification.
 * Under normal circumstances you need not edit this.
 */
int main()
{
    RawModificationTestSuite test;
    int success = test.run();
    std::cout << test.report() << std::endl;
    return success;
}

