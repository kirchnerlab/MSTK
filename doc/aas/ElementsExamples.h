/*!
\page ElementsExamples Examples for Elements
\dontinclude examples/aas/Elements.cpp

See examples/aas/Elements.cpp

\section rse Retrieve a standard element

\skipline Element C

\section cce Create a custom element

\skip ElementImpl::ElementImplKeyType
\until addIsotope

\section ace Add a custom element

\subsection acem manually

\skipline customElementRef

\subsection acebcf by convenience functions

\skip addElement
\until }

\skip addElement
\until }

\section rce Retrieve a custom element

\skipline retrievedCustomElementRef

\section ose Override a standard element

\subsection osefs from scratch

\skip ElementImpl
\until }

\subsection oseuse using standard element

\skip ElementImpl O
\until }

*/
