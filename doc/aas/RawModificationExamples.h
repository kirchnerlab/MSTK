/*!
\page RawModificationExamples Examples for RawModification
\dontinclude examples/aas/RawModifications.cpp

See examples/aas/RawModifications.cpp

\section rsm Retrieve a standard modification

\skipline RawModification mod

\section ccm Create a custom modification

\skip RawModificationImpl
\until setStoichiometry

\section acs Add a custom specificity

\skip Specificity
\until addSpecificity

\section acm Add a custom modification
\subsection acmm manually

\skipline customModRef

\subsection acmbcf by convenience functions

\skip addRawModification
\until }

\skip std::vector
\until }

\section rcm Retrieve a custom modification

\skipline RawModification retrieved

\section rs Retrieve the stoichiometry

\skipline Stoichiometry

\section osm Override a standard modification

\subsection osmfs from scratch

\skip mod1
\until } 

\subsection osmusm using standard modification

\skip mod2
\until }

*/
