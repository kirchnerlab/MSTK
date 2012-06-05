/*!
\page Complex1 Complex example 1
\dontinclude examples/aas/Complex1.cpp

See examples/aas/Complex1.cpp

\section main Create an amino acid sequence using custom amino acids with custom modifications and a custom element configuration

\subsection cce create/add custom element

\skip ElementImpl
\until addElement

\subsection csc create/add stoichiometry configuration

\skip StoichiometryConfigImplKeyType
\until addStoichiometryConfig

\subsection ccrm create/add custom raw modification

\skip RawModificationImplKeyType
\until addRawModification

\subsection cm create a modification using the custom raw modification and custom stoichiometry configuration

\skipline customMod

\subsection craa create/add raw amino acid

\skip RawAminoAcidImplKeyType
\until addRawAminoAcid

\subsection ccm create custom modification

\skip acetylKey
\until POST_TRANSLATIONAL

\subsection caas create amino acid sequence and apply custom modification

\skip AminoAcidSequence
\until applyModificationAtPosition

*/
