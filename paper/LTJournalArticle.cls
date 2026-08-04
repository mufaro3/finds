%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LaTeX Templates Journal Article
% LaTeX Class
% Version 2.0 (February 7, 2023)
%
% This class originates from:
% https://www.LaTeXTemplates.com
%
% Author:
% Vel (vel@latextemplates.com)
%
% License:
% CC BY-NC-SA 4.0 (https://creativecommons.org/licenses/by-nc-sa/4.0/)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%	CLASS CONFIGURATION

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{LTJournalArticle}
              [2023/02/07 LaTeX Templates Journal Article Class v2.0]

% Required for conditional logic and easily changing commands
\usepackage{etoolbox}

% Create toggle for a class option
\newtoggle{unnumberedsections}

% Default value for the class option
\settoggle{unnumberedsections}{false}

% Set the class option toggle if the class option was used in the template
\DeclareOption{unnumberedsections}{\settoggle{unnumberedsections}{true}}

% Pass through any extra options specified to the base class
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}

% Process class options
\ProcessOptions\relax

\LoadClass[twocolumn]{article} % Load the base class

%	REQUIRED PACKAGES AND MISC CONFIGURATIONS

\usepackage{graphicx}
\graphicspath{{figures/}}

% Force footnotes to the bottom of the page and to the left margin
\usepackage[bottom, hang]{footmisc}

% Horizontal space between the footnote marker and text
\setlength{\footnotemargin}{6pt}

%	MARGINS

\usepackage[
  top=2.5cm,
  bottom=2.5cm,
  left=2cm,
  right=2cm,
  footskip=1cm,
  headsep=0.75cm,
  columnsep=20pt
]{geometry}

%	FONTS

% Required for inputting international characters
\usepackage[utf8]{inputenc}

% Output font encoding for international characters
\usepackage[T1]{fontenc}

% Use the Palatino font
%\usepackage[sc]{mathpazo}

% Increase line spacing slightly
\linespread{1.05}

% Slightly tweak font spacing for aesthetics
\usepackage{microtype}

%	HEADERS AND FOOTERS

% Required for customizing headers and footers
\usepackage{fancyhdr}

% Enable custom headers and footers
\pagestyle{fancy}

% Top horizontal rule thickness
\renewcommand{\headrulewidth}{0pt}

% Clear default headers/footers
\fancyhf{}

\fancyhead[RO]{\small\textit{\runninghead}} % Right-odd page header
\fancyhead[LE]{\small\textit{\runninghead}} % Left-even page header

\fancyfoot[RO]{\small\textbf{\thepage}} % Right-odd page footer
\fancyfoot[LO]{\footnotesize\footertext} % Left-odd page footer
\fancyfoot[LE]{\small\textbf{\thepage}} % Left-even page footer
\fancyfoot[RE]{\footnotesize\footertext} % Left-even page footer

%	TITLE SECTION CUSTOMIZATION

% Required for customizing the title section
\usepackage{titling}

% Move the title up
\setlength{\droptitle}{-4\baselineskip}

% Article title pre-formatting
\pretitle{\begin{center}\huge\bfseries}

% Article title post-formatting
\posttitle{\end{center}}

% Left margin for the first \thanks line
\setlength{\thanksmarkwidth}{3pt}

% Left margin for the second and onwards \thanks line
\setlength{\thanksmargin}{-3pt}

% Set the headers and footers style for the first page to empty
\patchcmd{\maketitle}{plain}{empty}{}{}

%	ABSTRACT CUSTOMIZATION

\usepackage{abstract} % Allows abstract customization

% Set the "Abstract" text to bold
\renewcommand{\abstractnamefont}{\normalfont\bfseries\vspace{0.5\baselineskip}}

% Set the abstract itself to small italic text
\renewcommand{\abstracttextfont}
             {\vspace{-0.5\baselineskip}\normalfont\small\itshape}

%	BIBLIOGRAPHY

\usepackage[
  backend=biber,
  citestyle=numeric,
  bibstyle=numeric,
  sorting=none
]{biblatex}

%	TABLES

% Required for better horizontal rules in tables
\usepackage{booktabs}

% Required for manipulating table columns
\usepackage{array}

% Define a new right-aligned paragraph column type
\newcolumntype{R}[1]{>{\raggedleft\arraybackslash}p{#1}}

% Define a new left-aligned (no justification) paragraph column type
\newcolumntype{L}[1]{>{\raggedright\arraybackslash}p{#1}}

% Define a new centered paragraph column type
\newcolumntype{C}[1]{>{\centering\arraybackslash}p{#1}}

%	CAPTIONS

% Required for customizing captions
\usepackage{caption}

% Vertical whitespace between figures/tables and the caption (default is 10pt)
\captionsetup{skip=6pt}

% Define caption font style
\captionsetup{labelfont={bf,small}, textfont={it,small}}

%	LISTS

% Required for list customization
\usepackage{enumitem}

% Customize spacing around and inside lists
\setlist{noitemsep}

%	LINKS AND COLOR

% Required for links

\usepackage{xcolor}
\definecolor{crimson}{rgb}{0.75, 0.0, 0.2}
\usepackage{hyperref}

\hypersetup{
	colorlinks=true,
	urlcolor=crimson,
	linkcolor=crimson,
	citecolor=crimson,
}


%	CUSTOM COMMANDS

\newcommand{\runninghead}[1]{\renewcommand{\runninghead}{#1}}
\newcommand{\footertext}[1]{\renewcommand{\footertext}{#1}}
