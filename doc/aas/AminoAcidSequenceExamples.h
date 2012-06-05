/*!
\page AminoAcidSequenceExamples Examples for AminoAcidSequence
\dontinclude examples/aas/AminoAcidSequence.cpp

See examples/aas/AminoAcidSequence.cpp


\section caas Create an amino acid sequence

\skipline aas

\section asm Add a standard modification

\skip mod1
\until applyModificationAtPosition

\skip mod2
\until applyFixedModifications

\section acm Add a custom modification

\skip RawModificationImplKeyType
\until applyModificationAtPosition

\section rm Remove a (custom) modification
\subsection bmk by modification key

\skipline aas.remove

\subsection bmr by modification ref

\skipline aas.remove(mod1

\section asc Apply stoichiometry config

\skip StoichiometryConfigImplKeyType
\until applyModificationStoichiometryConfig

\section rs Retrieve the stoichiometry

\skipline getStoichiometry

*/
