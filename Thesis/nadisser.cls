\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{nadisser}[2018/07/09 My custom class for my PhD thesis]

\RequirePackage{etoolbox}

%prepend english and ngerman to global options
\ifx\@empty\@classoptionslist
   \preto\@classoptionslist{ngerman,english}
\else
   \preto\@classoptionslist{ngerman,english,}
\fi

\DeclareOption*{
   \PassOptionsToClass{\CurrentOption}{nathesis}
}

\ProcessOptions \relax

\LoadClass{nathesis}[2018/07/09]

%\newcommand\rulethickness{0.4pt}
%\newcommand\ruleshort{\rule{3cm}{\rulethickness}}
%\newcommand\rulemiddle{\rule{5cm}{\rulethickness}}
%\newcommand\rulelong{\rule{8cm}{\rulethickness}}

\newcommand\birthplace[1]{\def\@birthplace{#1}}
%\newcommand\promoter[2]{\def\@promoter{#1} \def\@promoteraff{#2}}
%\newcommand\copromoter[2]{\def\@copromoter{#1} \def\@copromoteraff{#2}}

\renewcommand\maketitle{%
   \begin{titlepage}
      \begin{otherlanguage}{ngerman}
         \vspace*{\fill}
         \parbox[t][0pt]{20em}{}
         \begin{center}
            \large
            {\Huge\textbf\@title}\\[2em]
            ACADEMISCH PROEFSCHRIFT\\ \vspace*{1em}
            ter verkrijging van de graad van doctor\\ \vspace*{0.5em}
            aan de Universiteit van Amsterdam\\ \vspace*{0.5em} 
            op gezag van de Rector Magnificus\\ \vspace*{0.5em}
            prof. dr. ir. K.I.J. Maex\\ \vspace*{0.5em}
            ten overstaan van een door het College voor Promoties ingestelde commissie,\\ \vspace*{0.5em}
	    in het openbaar te verdedigen in de Online\\ \vspace*{0.5em}
            op dinsdag 26 januari 2021, te 12:00 uur\\ \vspace*{0.5em}
	    door \@author \\
	    geboren te \@birthplace
         \end{center}
         \vspace*{\fill}
         \newpage
         \thispagestyle{empty}
         %\vspace*{\fill}
	 \hspace*{-1.8cm}
	 \begin{tabular}{p{3cm} p{5cm} p{8cm}}
	    \textbf{\textit{Promotiecommissie}} & &\\
	    &&\\
            \textit{Promotor:}& dr. T. Bereau \qquad & Universiteit van Amsterdam \\
	    &&\\
          \textit{Copromotor:}  &prof. dr. E.~J. Meijer \qquad &  Universiteit van Amsterdam\\
		&&\\
	    \textit{Overige leden:} & prof. dr. P.~G. Bolhuis   \qquad & Universiteit van Amsterdam\\
			   & prof. dr. K. Kremer       \qquad & Max Planck Institute for Polymer Research\\
			   & prof. dr. D. Bonn		\qquad & Universiteit van Amsterdam\\
		 	   & prof. dr.~ir. M. Dijkstra 	\qquad & Universiteit Utrecht\\
			   & dr. E. Lerner		\qquad & Universiteit van Amsterdam\\
         \end{tabular}
	
	 \hspace*{-1.8cm}
         Faculteit der Natuurwetenschappen, Wiskunde en Informatica
 	
      \end{otherlanguage}
   \end{titlepage}%
}
\endinput
