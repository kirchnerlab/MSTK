/*!
\page StoichiometryConfigExamples Examples for StoichiometryConfigs
\dontinclude examples/aas/StoichiometryConfigs.cpp

See examples/aas/StoichiometryConfigs.cpp

\section rdsc Retrieve default stoichiometry configuration

\skip defaultStoichiometry
\until DEFAULT_ELEMENT_CONFIG

\section ccsc Create a custom stoichiometry configuration

\skip StoichiometryConfigImpl
\until insertElement(Element

\section acsc Add a custom stoichiometry configuration

\subsection acscm manually

\skipline customConfigRef

\subsection acscbcf by convenience functions

\skip addStoichiometryConfig
\until }

\skip addStoichiometryConfig
\until }

\section odsc Override the default stoichiometry configuration

\subsection odscfs from scratch

not recommended

\subsection odscfs using standard mapping

\skip defaultConfig
\until }

*/
