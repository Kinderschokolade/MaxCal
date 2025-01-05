\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{nathesis}[2018/07/09 My custom class for my PhD thesis]

%\RequirePackage{xstring}
%
%\def\AfterEqualSignImpl#1=#2\relax{#2}
%\def\AfterEqualSign#1{\expandafter\AfterEqualSignImpl#1\relax}
%
%\newcommand\BabelOptDefault{
%   \ClassWarning{nathesis}{No babelopt has been set. Defaulting to babelopt=english}
%   \PassOptionsToPackage{english}{babel}
%}
%
%\DeclareOption{opt}{\renewcommand\content{OPTION}}
%\DeclareOption*{
%   \IfBeginWith\CurrentOption{babelopt=}{
%      \PassOptionsToPackage{\AfterEqualSign\CurrentOption}{babel}
%      \renewcommand\BabelOptDefault{}
%   }{
%      \PassOptionsToClass{\CurrentOption}{scrreprt}
%   }
%}
%\ProcessOptions \relax
%
%\BabelOptDefault
%\let\BabelOptDefault\undefined

\PassOptionsToClass{twoside,DIV=calc,headsepline}{scrreprt}

\DeclareOption*{
   \PassOptionsToClass{\CurrentOption}{scrreprt}
}

\ProcessOptions \relax

\LoadClass{scrreprt}
\RequirePackage{graphicx}
\RequirePackage{geometry}
\RequirePackage{mathtools}
\RequirePackage{babel}
\RequirePackage[hidelinks]{hyperref}
\RequirePackage{cleveref}
\pagestyle{headings}

\endinput
