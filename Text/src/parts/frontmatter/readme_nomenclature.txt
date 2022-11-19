Hinweise zur Benutzung der automatischen Nomenklatur.

1. 	Nomenklatur Block im Dokumenten Header einkommentieren:
		%==================================================================
		%Einheiten
		\usepackage{units}
		%%Nomenklauturpaket
		\usepackage[intoc]{nomenclTRi_A4}
		\renewcommand{\nomname}{}
		%Definition der Einheitenspalte
		\newcommand{\nomenunit}[1]{\hskip-0em\parbox{4em}{$\left[ #1 \right]$}}
		
		\makenomenclature
		%==================================================================

2. 	Zeile im Dokument zum Einbinden der Nomenklatur einkommentieren:
		%\input{src/parts/frontmatter/nomenclature.tex}

3. 	Die Variablen für die Nomenklatur können irgendwo im Dokument stehen. Beispiele befinden sich in ./src/parts/frontmatter/nomenclature.tex.

4. 	Die folgende Anweisung gilt für MikTex in Kombination mit TexWorks als Editor:
	- einmal mit pdflatex.exe kompilieren
	- einmal mit makeindex.exe und folgenden Argumenten komplilieren:
		makeindex.exe
		  $basename.nlo
		  -s
		  nomencl.ist
		  -o
		  $basename.nls
	- erneut mit pdflatex.exe kompilieren
	