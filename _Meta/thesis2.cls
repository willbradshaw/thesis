%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thesis document class
% Simplified from Krishna Kumar's thesis2 class
% (github.com/kks32/phd-thesis-template/)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ************************** Class Identification *****************************
\newcommand\fileversion{w}
\newcommand\filedate{2019-04-01}
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{thesis2}

% ******************************************************************************
% **************************** Class Definition ********************************
% ******************************************************************************

% *********************** Define a Print/Online Version ************************
\newif\if@print\@printfalse
\DeclareOption{print}{\@printtrue}

% ****************** Chapter Mode - To print only selected chapters ************
\newif \ifdefineChapter\defineChapterfalse

\newif\ifPHD@chapter\PHD@chapterfalse
\DeclareOption{chapter}{
  \PHD@chaptertrue
  \ClassWarning{thesis2}{You have chosen an option that generates only 
    selected chapters with references, if this was intentional, ignore this
    warning. If not, please remove the option `chapter' from the document class
    and recompile. Good luck with your writing!}
}

\ProcessOptions\relax%

% *************************** Pre-defined Options ******************************

% Font Size
\newcommand\PHD@ptsize{11pt}
\DeclareOption{10pt}{\renewcommand\PHD@ptsize{10pt}}%
\DeclareOption{11pt}{\renewcommand\PHD@ptsize{11pt}}%
\DeclareOption{12pt}{\renewcommand\PHD@ptsize{12pt}}%
\PassOptionsToClass{\PHD@ptsize}{book}%
\message{thesis2: Font size is set as \PHD@ptsize}

% Page Size
\newcommand\PHD@papersize{a4paper} % Set Default as a4paper
\DeclareOption{a4paper}{\renewcommand\PHD@papersize{a4paper}}%
\DeclareOption{a5paper}{\renewcommand\PHD@papersize{a5paper}}%
\DeclareOption{letterpaper}{\renewcommand\PHD@papersize{letterpaper}}%

\PassOptionsToClass{\PHD@papersize}{book}%

% Column layout
\newif\if@oneside\@onesidefalse
\DeclareOption{oneside}{\PassOptionsToClass{\CurrentOption}{book}\@onesidetrue}%
\DeclareOption{twoside}{\PassOptionsToClass{\CurrentOption}{book}}%

% Line Numbers
\newif\ifPHD@lineno\PHD@linenotrue
\DeclareOption{lineno}{\PHD@linenotrue}
\DeclareOption{nolineno}{\PHD@linenofalse}


% Generates Warning for unknown options
\DeclareOption*{
  \ClassWarning{thesis2}{Unknown or non-standard option
'\CurrentOption'. I'll see if I can load it from the book class. If you get a
warning unused global option(s): `\CurrentOption` then the option is not
supported!}
  \PassOptionsToClass{\CurrentOption}{book}
}

% Determine whether to run pdftex or dvips
\ProcessOptions\relax%
\newif\ifsetDVI\setDVIfalse
\ifx\pdfoutput\undefined
  % we are not running PDFLaTeX
  \setDVItrue
  \LoadClass[dvips,fleqn,openright]{book}%
\else % we are running PDFLaTeX
  \ifnum \pdfoutput>0
    %PDF-Output
    \setDVIfalse
    \LoadClass[pdftex,fleqn,openright]{book}%
  \else
    %DVI-output
    \setDVItrue
    \LoadClass[fleqn,openright]{book}%
  \fi
\fi

% ******************************************************************************
% ******************************** Packages ************************************
% ******************************************************************************
% ************************** Layout and Formatting *****************************

\RequirePackage{lscape}   % Supports Landscape Layout
\RequirePackage{setspace} % Define line spacing in paragraph
\RequirePackage{calc}     % calculate vertical spacing

% ************************* Conditional Statements *****************************
\RequirePackage{ifthen}   % Conditional statements
\RequirePackage{ifpdf}    % Check for pdfLaTeX
\RequirePackage{ifxetex}  % XeLaTeX

% *********************** Table of Contents & Appendices ***********************
% add Bibliography, List of figures and tables to contents
\RequirePackage[nottoc]{tocbibind}
% Add appendices
\RequirePackage[title,titletoc]{appendix}

\renewcommand{\contentsname}{Table of contents}
\renewcommand{\listfigurename}{List of figures}
\renewcommand{\listtablename}{List of tables}

% Hides Contents appearing from TOC, but adds it to bookmarks
\let\temptableofcontents\tableofcontents
\renewcommand{\tableofcontents}{
  \cleardoublepage
  \pdfbookmark[0]{\contentsname}{Contents}
  \temptableofcontents
}

% *************************** Graphics and Figures *****************************
\RequirePackage[usenames, dvipsnames]{color}
\ifpdf
  % Convert eps figures to pdf
  \RequirePackage{epstopdf}
  \RequirePackage[pdftex,demo]{graphicx}
  \DeclareGraphicsExtensions{.png, .jpg, .pdf}
  \pdfcompresslevel=9
  \graphicspath{{Figs/Raster/}{Figs/}}
\else
  \RequirePackage[demo]{graphicx}
  \DeclareGraphicsExtensions{.eps, .ps}
  \graphicspath{{Figs/Vector/}{Figs/}}
\fi



% ******************************* Time Stamp ***********************************
% Compute the timestamp based on an idea of
% Tim Piessens \texttt{<Tim.Piessens@esat.kuleuven.ac.be>}.

\RequirePackage{datetime}

\newcount\PHD@hour \newcount\PHD@minute
\PHD@hour=\time
\divide \PHD@hour by 60
\PHD@minute=\time
\count@=\PHD@hour
\multiply \count@ by -60
\advance \PHD@minute by \count@
\newcommand{\PHD@daytime}{%
  \today\space--\space%
  \ifnum\PHD@hour=0 00\else\ifnum\PHD@hour<10 0\fi%
  \number\PHD@hour\fi:\ifnum\PHD@minute<10 0\fi\number\PHD@minute
}

% ******************************* Bibliography *********************************

\RequirePackage[backend=biber, citestyle=numeric-comp, style=nature, intitle=true, sorting=none, natbib=true, maxnames=3, minnames=1, url=false, giveninits=true, sortcites=true, date=year, doi=false,isbn=false]{biblatex}
\renewcommand{\intitlepunct}{\addspace\nopunct}
\renewcommand{\bibname}{References}

% ********************************** Fonts **********************************
\RequirePackage{textcomp}
\RequirePackage{mathptmx} % times roman, including math (where possible)

\ifxetex
  % XeLaTeX
  \usepackage{amsmath}
  \usepackage{fontspec}
  \usepackage[]{unicode-math}
  \setmainfont[
    Extension = .otf,
    UprightFont = *-regular,
    BoldFont = *-bold,
    ItalicFont = *-italic,
    BoldItalicFont = *-bolditalic,
  ]{xits}
  
  \setmathfont[ 
    Extension = .otf,
    BoldFont = *bold,
  ]{xits-math}
  
\else
  % default: pdfLaTeX
  \RequirePackage[utf8]{inputenc}
  \RequirePackage[T1]{fontenc}
  
  % If building with PDFLaTeX, use microtype spacing adjustments
  \RequirePackage[final]{microtype}
  
  \input{glyphtounicode}
  \pdfglyphtounicode{f_f}{FB00}
  \pdfglyphtounicode{f_i}{FB01}
  \pdfglyphtounicode{f_l}{FB02}
  \pdfglyphtounicode{f_f_i}{FB03}
  \pdfglyphtounicode{f_f_l}{FB04}
  \pdfgentounicode=1

  \RequirePackage{amsfonts}
  \RequirePackage{amsmath}
  \RequirePackage{amssymb}
\fi

% Don't break enumeration (etc.) across pages in an ugly manner
\clubpenalty=10000
\widowpenalty=10000

%******************************* Print / Online ********************************
% Defines a print / online version to define page-layout and hyperrefering

% Moved below other usepackage definitions to fix PDFLaTeX footnote warning
% Based on the suggestion by John Plaice

\ifsetDVI
  \special{papersize=\the\paperwidth,\the\paperheight}
  \RequirePackage[unicode=true, hidelinks]{hyperref}  % extra option dvips
\else
  \RequirePackage[unicode=true, hidelinks]{hyperref}
  \pdfpagewidth=\the\paperwidth
  \pdfpageheight=\the\paperheight
\fi

\if@print
  % For Print version
  \hypersetup{
    final=true,
    plainpages=false,
    pdfstartview=FitV,
    pdftoolbar=true,
    pdfmenubar=true,
    bookmarksopen=true,
    bookmarksnumbered=true,
    breaklinks=true,
    linktocpage,
    colorlinks=true,
    linkcolor=black,
    urlcolor=black,
    citecolor=black,
    anchorcolor=black
  }
    \ifsetDVI
      % Odd and Even side Margin for binding and set viewmode for PDF
      \RequirePackage[dvips,paper=\PHD@papersize,hmarginratio=1:1,
        vmarginratio=1:1,scale=0.75,bindingoffset=5mm]{geometry}
    \else
      \RequirePackage[pdftex,paper=\PHD@papersize,hmarginratio=1:1,
        vmarginratio=1:1,scale=0.75,bindingoffset=5mm]{geometry}
    \fi

  \if@twoside
    \hypersetup{pdfpagelayout=TwoPageRight}
  \else
    \hypersetup{pdfpagelayout=OneColumn}
  \fi

\else
  % For PDF Online version
  \hypersetup{
    final=true,
    plainpages=false,
    pdfstartview=FitV,
    pdftoolbar=true,
    pdfmenubar=true,
    bookmarksopen=true,
    bookmarksnumbered=true,
    breaklinks=true,
    linktocpage,
    colorlinks=true,
    linkcolor=blue,
    urlcolor=blue,
    citecolor=blue,
    anchorcolor=green
  }

      % No Margin staggering on Odd and Even side
      \ifsetDVI
        \RequirePackage[paper=\PHD@papersize,hmarginratio=1:1,
          vmarginratio=1:1,scale=0.75]{geometry} % dvips
      \else
        \RequirePackage[pdftex,paper=\PHD@papersize,hmarginratio=1:1,
          vmarginratio=1:1,scale=0.75]{geometry}
      \fi

    \hypersetup{pdfpagelayout=OneColumn}
\fi

% ************************ URL Package and Definition **************************
\RequirePackage{url}
% Redefining urlstyle to use smaller fontsize in References with URLs
\newcommand{\url@leostyle}{%
 \@ifundefined{selectfont}{\renewcommand{\UrlFont}{\sffamily}}
 {\renewcommand{\UrlFont}{\normalsize}}}
\urlstyle{leo}

% option to split urls over multiple lines for latex >> DVIPS >> PDF option
% PDFLaTeX does it automatically.

\ifxetex % XeTeX doesn't require breakurl only dvi requires break url
\else
  \ifsetDVI
    \RequirePackage{breakurl} % to split the url over multiple lines
  \fi
\fi  

% ******************************************************************************
% **************************** Pre-defined Settings ****************************
% ******************************************************************************

% *************************** Setting PDF Meta-Data ****************************
\ifpdf
\AtBeginDocument{
  \hypersetup{
    pdftitle = {\@title},
    pdfauthor = {\@author},
    pdfsubject={\@subject},
    pdfkeywords={\@keywords}
  }
}
\fi

% ******************************** Line Spacing ********************************
% Set spacing as 1.5 line spacing for the PhD Thesis
% In practice, fortunately, nobody knows really what “one-and-a-half spaced
% type” means exactly (in terms of millimetres baseline distance). The following
% LaTeX setting has routinely been considered perfectly acceptable:

\renewcommand\baselinestretch{1.2}

% ************************** TOC and Hide Sections *****************************
\newcommand{\nocontentsline}[3]{}
\newcommand{\tochide}[2]{
	\bgroup\let
	\addcontentsline=\nocontentsline#1{#2}
	\egroup}
% Removes pagenumber appearing from TOC
\addtocontents{toc}{\protect\thispagestyle{empty}}


% ***************************** Header Formatting ******************************
% Custom Header with Chapter Number, Page Number and Section Numbering

\RequirePackage{fancyhdr} % Define custom header

% Set Fancy Header Command is defined to Load FancyHdr after Geometry is defined
\newcommand{\setFancyHdr}{
\pagestyle{fancy}
% Default Style: Sets Page Number at the Top (LE/RO) with Chapter/Section Name
% on LO/RE and an empty footer
\renewcommand{\chaptermark}[1]{\markboth {##1}{}}
\renewcommand{\sectionmark}[1]{\markright{\thesection\ ##1}}
\fancyhf{}
\fancyhead[LO]{\nouppercase \rightmark}
\fancyhead[LE,RO]{\bfseries\thepage}
\fancyhead[RE]{\nouppercase \leftmark}
}

\setlength{\headheight}{14.5pt}
%\renewcommand{\headrulewidth}{0.5pt}
%\renewcommand{\footrulewidth}{0pt}
\fancypagestyle{plain}{
  \fancyhead{}
  \renewcommand{\headrulewidth}{0pt}
}

\setFancyHdr

% **************** Clear Header Style on the Last Empty Odd pages **************
\renewcommand{\cleardoublepage}{\clearpage\if@twoside \ifodd\c@page\else%
	\hbox{}%
	\thispagestyle{empty}  % Empty header styles
	\newpage%
	\if@twocolumn\hbox{}\newpage\fi\fi\fi}


% ******************************************************************************
% **************************** Macro Definitions *******************************
% ******************************************************************************
% These macros are used to declare arguments needed for the
% construction of the title page and other preamble.

% Subtitle (optional)
\newcommand{\@subtitle}{}
\newcommand{\subtitle}[1]{\renewcommand{\@subtitle}{#1}}

% The year and term the degree will be officially conferred
\newcommand{\@degreedate}{\the\day.\space\monthname[\the\month]\space\the\year}
\newcommand{\degreedate}[1]{\renewcommand{\@degreedate}{#1}}
\newcommand{\@degreeyear}{\the\year}
\newcommand{\degreeyear}[1]{\renewcommand{\@degreeyear}{#1}}

% The full (unabbreviated) name of the degree
\newcommand{\@degreetitle}{}
\newcommand{\degreetitle}[1]{\renewcommand{\@degreetitle}{#1}}

% The name of your department(eg. Engineering, Maths, Physics)
\newcommand{\@dept}{}
\newcommand{\dept}[1]{\renewcommand{\@dept}{#1}}

% The name of your college (eg. King's)
\newcommand{\@college}{}
\newcommand{\college}[1]{\renewcommand{\@college}{#1}}

% The name of your University
\newcommand{\@university}{}
\newcommand{\university}[1]{\renewcommand{\@university}{#1}}

% Defining the crest
\newcommand{\@crest}{}
\newcommand{\crest}[1]{\renewcommand{\@crest}{#1}}

% Defining the college crest
\newif\ifPHD@collegeshield\PHD@collegeshieldfalse
\newcommand{\@collegeshield}{}
\newcommand{\collegeshield}[1]{\renewcommand{\@collegeshield}{#1}\PHD@collegeshieldtrue}

% Supervisor
\newif\ifPHD@supervisor\PHD@supervisorfalse
\newcommand{\@supervisor}{}
\newcommand{\supervisor}[1]{\renewcommand{\@supervisor}{#1}\PHD@supervisortrue}

% Supervisor Title (Supervisor - Default, can be changed)
\newcommand{\@supervisorrole}{\textbf{Gutachter:} }
\newcommand{\supervisorrole}[1]{\renewcommand{\@supervisorrole}{#1}}

% Supervisor width to align
\newcommand{\@supervisorlinewidth}{0.5\textwidth}
\newcommand{\supervisorlinewidth}[1]{\renewcommand{\@supervisorlinewidth}{#1}}

% Advisor
\newif\ifPHD@advisor\PHD@advisorfalse
\newcommand{\@advisor}{}
\newcommand{\advisor}[1]{\renewcommand{\@advisor}{#1}\PHD@advisortrue}

% Advisor Title (Advisor - Default, can be changed)
\newcommand{\@advisorrole}{Advisor: }
\newcommand{\advisorrole}[1]{\renewcommand{\@advisorrole}{#1}}

% Advisor width to align
\newcommand{\@advisorlinewidth}{0.25\textwidth}
\newcommand{\advisorlinewidth}[1]{\renewcommand{\@advisorlinewidth}{#1}}

% Submission Text
\newcommand{\submissiontext}{\Large{Inaugural-Disseration\\
zur\\
Erlangung des Doktorgrades\\
der Mathematisch-Naturwissenschaftlichen Fakult\"at\\
der Universit\"at zu K\"oln\\
vorgelegt von}\\}

\newcommand{\pruefungstext}{\textbf{Tag der m\"undlichen Pr\"ufung:} 6. Juni 2019}

% keywords (These keywords will appear in the PDF meta-information
% called `pdfkeywords`.)
\newcommand{\@keywords}{}
\newcommand{\keywords}[1]{\renewcommand{\@keywords}{#1}}

% subjectline (This subject will appear in the PDF meta-information
% called `pdfsubject`.)
\newcommand{\@subject}{}
\newcommand{\subject}[1]{\renewcommand{\@subject}{#1}}


% These macros define an environment for front matter that is always
% single column even in a double-column document.
\newenvironment{alwayssingle}{%
       \@restonecolfalse\if@twocolumn\@restonecoltrue\onecolumn
       \else\newpage\fi}
       {\if@restonecol\twocolumn\else\newpage\fi}

% Set single column even under two column layout
\newcommand{\setsinglecolumn}{
\if@twocolumn
   \onecolumn
\else
\fi
}


% ******************************************************************************
% *************************** Front Matter Layout ******************************
% ******************************************************************************

% ******************************** Title Page **********************************
\renewcommand{\maketitle}{

% To compute the free vertical space in Title page
\computeTitlePageSpacing

\thispagestyle{empty}
\begin{center}

\begin{singlespace}

% University Crest Long if college crest is defined
\ifPHD@collegeshield%
	{\usebox{\PHD@crest}}
	\vspace{.15\PHD@titlepagespacing}
\fi

% Title
{\usebox{\PHD@Title}} % subtitle is defined


% Crest

\ifPHD@collegeshield%
	\vspace{.15\PHD@titlepagespacing}
	{\usebox{\PHD@collegecrest}}
	\vspace{.1\PHD@titlepagespacing}
\else
	\vspace{.25\PHD@titlepagespacing}
	{\usebox{\PHD@crest}}
	\vspace{.15\PHD@titlepagespacing}
\fi

% Submission Text
\vspace*{1em}
{\onehalfspacing
\Large{Inaugural-Disseration\\
zur\\
Erlangung des Doktorgrades\\
der Mathematisch-Naturwissenschaftlichen Fakult\"at\\
der Universit\"at zu K\"oln\\
\vspace*{2em}
vorgelegt von}\vspace{0.75em}\\}
{\LARGE\bfseries{\@author}}\vspace{0.75em}
{\usebox{\PHD@dept}}\\
\vspace*{5em}
{\Large{K\"{o}ln \@degreeyear{}}}
\vspace{.2\PHD@titlepagespacing}

\clearpage
% Supervisor
\ifPHD@supervisor%
{\usebox{\PHD@supervisor}}
\vspace*{0.5em}
\fi

% Advisor
\ifPHD@advisor%
{\usebox{\PHD@advisor}}
\vspace*{0.5em}
\fi

\vspace*{1em}
{\usebox{\PHD@pruefung}}
\vspace*{2em}

\end{singlespace}

\end{center}
}

% ******************************* Declaration **********************************
% The declaration environment puts a large, bold, centered
% "Declaration" label at the top of the page.

% Add date and name under signature space
\newcommand{\declare}[1]{\begin{flushright}\vspace{#1}\@author{}\\\@degreedate{}\end{flushright}\vfill}

% ***************************** Chapter Mode ***********************************
% The chapter mode allows user to only print particular chapters with references
% All other options are disabled by default
% To include only specific chapters without TOC, LOF, Title and Front Matter
% To send it to supervisior for changes

\ifPHD@chapter
    \defineChaptertrue
    % Disable the table of contents, figures, tables, index and nomenclature
    \renewcommand{\maketitle}{}
    \renewcommand{\tableofcontents}{}
    \renewcommand{\listoffigures}{}
    \renewcommand{\listoftables}{}
    \renewcommand{\printnomencl}{}
\else
    \defineChapterfalse
\fi

% ******************** To compute empty space in title page *******************
% Boxes below are used to space different contents on the title page
% Based on https://github.com/cambridge/thesis

\newcommand{\computeTitlePageSpacing}{


% Title Box
\newsavebox{\PHD@Title}
\begin{lrbox}{\PHD@Title}
  \begin{minipage}[c]{0.98\textwidth}
    \centering \Huge \bfseries{\@title}
    \ifthenelse{\equal{\@subtitle}{}}{
      % subtitle is not defined
    }{
      \\
      \centering \Large {\@subtitle}
    } % subtitle is defined
  \end{minipage}
\end{lrbox}

% University Crest Box
\newsavebox{\PHD@crest}
\begin{lrbox}{\PHD@crest}
  {\@crest \par}
\end{lrbox}

\newsavebox{\PHD@collegecrest}
\begin{lrbox}{\PHD@collegecrest}
  \ifPHD@collegeshield%
    {\@collegeshield}
  \else
    % College shield is undefined
  \fi
\end{lrbox}

% Author Box
\newsavebox{\PHD@author}
\begin{lrbox}{\PHD@author}
  \begin{minipage}[c]{\textwidth}
    \centering \Large \bfseries{\@author}
      \vspace{0.5em}   
  \end{minipage}
\end{lrbox}

% Supervisor Box
\newsavebox{\PHD@supervisor}
\begin{lrbox}{\PHD@supervisor}
  \begin{minipage}[c]{\textwidth}
    \ifthenelse{\equal{\@supervisor}{}}{
      % supervisor is not defined
    }{
      \begin{center}
        \makebox[\widthof{supervisorrole} + \@supervisorlinewidth][c]{
          \Large
          \@supervisorrole \parbox[t]{\@supervisorlinewidth}{\@supervisor}
        }
      \end{center}
      \vspace{0.5em}
    } % supervisor is defined
  \end{minipage}
\end{lrbox}

% Advisor Box
\newsavebox{\PHD@advisor}
\begin{lrbox}{\PHD@advisor}
  \begin{minipage}[c]{\textwidth}
    \ifthenelse{\equal{\@advisor}{}}{
      % advisor is not defined
    }{
      \begin{center}
        \makebox[\widthof{advisorrole} + \@advisorlinewidth][c]{
          \Large
          \@advisorrole \parbox[t]{\@advisorlinewidth}{\@advisor}
        }
      \end{center}
      \vspace{0.5em}
    } % advisor is defined
  \end{minipage}
\end{lrbox}

% Department Box
\newsavebox{\PHD@dept}
\begin{lrbox}{\PHD@dept}
  \begin{minipage}[c]{\textwidth}
    \centering {\Large \@dept \par}
	       {\large \@university \par}
  \end{minipage}
\end{lrbox}

% Submission Box
\newsavebox{\PHD@submission}
\begin{lrbox}{\PHD@submission}
  \begin{minipage}[c]{\textwidth}
    \begin{center}
      \large \submissiontext \par
      \large \textit {\@degreetitle} \par
    \end{center}
  \end{minipage}
\end{lrbox}

\newsavebox{\PHD@pruefung}
\begin{lrbox}{\PHD@pruefung}
  \begin{minipage}[c]{\textwidth}
    %\begin{center}
      \Large \pruefungstext \par
    %\end{center}
  \end{minipage}
\end{lrbox}

% College and Date Box
\newsavebox{\PHD@collegedate}
\begin{lrbox}{\PHD@collegedate}
\begin{minipage}[c]{\textwidth}
  \large
  \ifthenelse{\equal{\@college}{}}{
    % college is not defined
    \begin{center}
      \@degreedate
    \end{center}
  }{
    % college is defined
    \begin{minipage}[b]{0.49\textwidth}
      \flushleft\@college
    \end{minipage}
    \begin{minipage}[b]{0.49\textwidth}
      \flushright \@degreedate
    \end{minipage}
  } % college is defined
\end{minipage}
\end{lrbox}

%  Now to compute the free vertical space
\newlength{\PHD@titlepagespacing}
\setlength{\PHD@titlepagespacing}{ \textheight %
			- \totalheightof{\usebox{\PHD@Title}}
			- \totalheightof{\usebox{\PHD@crest}}
			- \totalheightof{\usebox{\PHD@collegecrest}}
			- \totalheightof{\usebox{\PHD@author}}
                   	- \totalheightof{\usebox{\PHD@supervisor}}
                   	- \totalheightof{\usebox{\PHD@advisor}}
			- \totalheightof{\usebox{\PHD@dept}}
			- \totalheightof{\usebox{\PHD@submission}}
			- \totalheightof{\usebox{\PHD@collegedate}}
}
}

