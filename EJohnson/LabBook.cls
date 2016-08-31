%%
%% This is file `labboode.tex',
%% generated with the docstrip utility.
%%
%% The original source files were:
%%
%% labbook.dtx  (with options: `labboode')
%% 
%% Copyright (C) 2002 Frank Küster.
%% 
%% All rights reserved.
\ProvidesFile{labboode.tex}%
            [2003/05/20 v1.1 LaTeX2e labbook document class]
\documentclass[a4paper]{article}
\usepackage{doc}
\usepackage{ngerman}
\usepackage[T1]{fontenc}
\usepackage[latin1]{inputenc}

\MacroIndent=2em

\MakeShortVerb{\|}

\begin{document}

\title{labbook.cls, eine \LaTeX-Klasse zum F"uhren von experimentellen
  Laborjournalen, Version 1.1}

\author{Frank K"uster}
\date{2003-05-20}
\maketitle

\begin{abstract}
  Diese Klasse soll es erm"oglichen, Laborb"ucher mit chronologischen
  Aufzeichnungen "uber Experimente mit \LaTeX\ zu f"uhren. Die neuen
  Seiten werden einfach jeden Tag hinten angeheftet, aus den
  Gliederungselementen wird ein Experiment-Index generiert, um das
  Auffinden von bestimmten Experimenten zu erleichtern. Durch
  unterschiedliche Seitennummerierung kann der Vorspann (Index,
  Inhalt, Abkürzungsverzeichnis) unabhängig vom Hauptteil verl"angert
  werden.  Die Klasse basiert auf der KOMA-Script-Klasse
  |scrbook.cls|. Es k"onnen daher alle Merkmale dieser Klasse
  verwendet werden -- die Lekt"ure des |scrguide|, der
  KOMA-Script-Dokumentation, wird daher sehr empfohlen.
\end{abstract}

\tableofcontents

\section{Generelles}

\subsection{Rechtliches}

Die F"uhrung von Laborb"uchern ist teilweise gesetzlich geregelt, oft
stellen Institutionen zur Forschungsf"orderung weitere Anforderungen.
In der Regel ist es daher \emph{nicht} erlaubt, ein Laborjournal
lediglich in elektronischer Form zu f"uhren; teilweise sind gebundene
B"ucher vorgeschrieben. Allerdings scheint es auch dort "ublich zu
sein, die Notizen, die man w"ahrend der Experimente auf losen
Bl"attern, in Computerdateien oder einem "`Schmierbuch"' gef"uhrt hat,
mehr oder weniger regelm"a"sig ins gebundene Laborjournal zu
"ubertragen.

Aus diesem Grund scheint es mir vertretbar, mit einer elektronischen
und einer parallelen ausgedruckten Fassung des Laborjournals zu
arbeiten. Bei Bedarf kann die ausgedruckte Fassung geheftet oder
gebunden werden, eine auf CD-ROM gebrannte Version sollte einem Buch
an Beweiskraft nicht nachstehen. Dem stehen die gro"sen Vorteile v.a. einer
PDF-Fassung gegen"uber: Volltextsuche, einfache Einbindung von
Graphen und Bildern, M"oglichkeit der Verlinkung von Datendateien etc.

\subsection{Merkmale}
\label{sec:features}

\begin{description}
\item[Basisklasse] |labbook.cls| basiert auf der KOMA-Script-Klasse
  |scrbook.cls| und bietet daher alle Möglichkeiten dieser Klasse, vom
  Seitenlayout über Kopfzeilengestaltung und Absatzauszeichnung bis
  hin zu Variationen des Layouts von \textit{floats} und
  Randbemerkungen. Die Lektüre des scrguide wird sehr empfohlen!
\item[Gliederung] Die Gliederung erfolgt zun"achst chronologisch,
  statt |\chapter| wird |\labday| verwendet. Innerhalb eines Tages
  wird nach Experimenten untergliedert -- das kann eine vollständige
  Messung von der Planung bis zur Auswertung sein, aber auch ein
  Arbeitsgang eines mehrere Tage dauernden Versuchs. |\experiment|
  steht auf der Ebene von |\section|, zusätzlich steht
  |\subexperiment| auf der |\subsection|-Ebene zur Verfügung. Darunter
  gibt es die "ublichen Gliederungsebenen, auch wenn eine sehr tiefe
  Schachtelung wenig sinnvoll erscheint -- man beachte das
  |\minisec|-Kommando von KOMA-Script.
\item[Index und Inhaltsverzeichnis] Wegen der chronologischen Gliederung ist ein
  traditionelles Inhaltsverzeichnis nicht ausreichend, besonders wenn
  Versuche in einem zusammenh"angenden Experiment an verschiedenen
  Tagen durchgef"uhrt werden. Daher wird zusätzlich ein Index
  ausgegeben. Der Index wird aus den Inhaltsverzeichnis-Einträgen der
  Abschnittsbefehle erstellt. Als Argument für |\experiment| und
  |\subexperiment| kann man auch vorher definierte Abkürzungen
  verwenden. Dadurch soll die konsistente Erstellung des Index
  erleichtert werden. Außerdem können einem Experiment mehrere
  Indexeinträge zugeordnet werden.
\item[hyperref-Integration] Die Klasse arbeitet gut mit dem Paket
  |hyperref.sty| zur Erstellung von pdf-Dateien mit
  Navigationsfunktionen zusammen.
\end{description}

\section{Benutzung}
\label{sec:usage}

\subsection{Initialisierung}
\label{sec:options}

Wenn man hyperref verwenden möchte \emph{muss} man die Klassenoption
|hyperref| angeben:
\begin{verbatim}
\documentclass[hyperref]{labbook}
\end{verbatim}
Danach kann man hyperref und andere Pakete in der jeweils passenden
Reihenfolge laden. Wird hyperref geladen, ohne die Option anzugeben,
dann überschreibt es einige Änderungen, die in dieser Klasse in
\LaTeX-Interna eingefügt werden. Umgekehrt werden diese Änderungen gar
nicht erst durchgeführt, wenn man die Option |hyperref| angibt, aber
das Laden des Paketes vergisst! Bei der Verwendung von hyperref
sollte man außerdem die KOMA-Option |idxtotoc| in Erwägung ziehen: Der
Index wird nämlich nicht nur in das Inhaltsverzeichnis aufgenommen,
sondern auch in die PDF-Bookmarks.

Ansonsten werden derzeit einfach alle Optionen an scrbook
durchgereicht. Die Verwendung von |openany| bietet sich an, damit
Labortage auf jeder Seite beginnen können.

\section{Gliederungsbefehle}
\label{sec:gliederungsbefehle}

subsection{\textbackslash\texttt{labday}}
label{sec:labday}

\DescribeMacro{\labday} Der Befehl |\labday| kann verwendet werden, um
eine nicht nummerierte "Uberschrift zu erzeugen. Ihr Text (bzw. der des
optionalen Arguments) wird ins Inhaltsverzeichnis aufgenommen und als
lebender Kolumnentitel gesetzt.  In der Regel wird man einfach das
Datum, vielleicht mit dem Wochentag, verwenden.  |\labday| kümmert
sich auch um die Indexeinträge der untergeordneten Gliederungspunkte;
man sollte es nur innerhalb von mainmatter verwenden. Intern wird
|\addchap| aufgerufen, es gibt keine Sternform und kein optionales
Argument.

\subsection{\textbackslash\texttt{experiment} und
  \textbackslash\texttt{subexperiment}: Einfache Verwendung}
\label{sec:exper-und-subexp}

\DescribeMacro{\experiment} Innerhalb eines Tages k"onnen nummerierte
"Uberschriften für einzelne Experimente mit
\begin{verbatim}
\experiment[<Kurzform>]{<Langform>}
\end{verbatim}
erzeugt werden. Ihr Text (bzw. der des optionalen Arguments) wird
nicht nur ins Inhaltsverzeichnis und den Seitenkopf geschrieben,
sondern auch in den Index. Dieser erleichtert die Orientierung,
besonders wenn Experimente mehrfach oder "uber mehrere Tage hinweg
durchgef"uhrt werden. Die Indexeintr"age geben nicht nur die Seiten
des Beginns jedes Experiments an, sondern den ganzen Seitenbereich,
und fassen zusammenhängende Bereiche des gleichen Experiments an
unterschiedlichen Tagen zusammen.

Beachten Sie, dass man im optionalen Argument keine Kommas verwenden
darf, weil es auch eine durch Kommas abgetrennte Liste enthalten kann
(s.\,u.). Wenn Sie doch ein Komma benötigen, schließen Sie das
optionale Argument nochmals in geschweifte Klammern ein:
\begin{verbatim}
\experiment[{eins, zwei, drei}]{Die Eins, die zwei und die Drei}
\end{verbatim}

\DescribeMacro{\subexperiment} Das selbe gilt f"ur die darunter
liegenden |\subexperiment|-Eintr"age. Sie sind f"ur Untergliederungen
der Art "`Planung, Durchf"uhrung, Auswertung"' oder "`Herstellung,
Reinigung, Messung"' vorgesehen.

\subsection{Fortgeschrittene Verwendung: Abk"urzungen}
\label{sec:abkurzungen}

\DescribeMacro{\newexperiment} Wenn man für das gleiche Experiment an
unterschiedlichen Stellen leicht unterschiedliche Formulierungen (oder
Schreibungen) wählt, erhält man unterschiedliche Indexeinträge. Um
dies zu vermeiden, kann man sich Abk"urzungen für häufig auftretende
Experimente definieren. Dies geschieht mit dem Makro

|\newexperiment{<abk>}{<Langform>}{<Kurzform>}|

Dabei bezeichnet |<abk>| die Abk"urzung, mit der man auf die
definierte |<Langform>| und |<Kurzform>| zugreifen kann (letztere für
Index, Inhaltsverzeichnis und Kolumnentitel). Der Text der Abk"urzung
ist dabei \emph{ohne} vorangestellten Backslash einzugeben, also als
|\experiment{abk}|! Die Abk"urzungen d"urfen außer der Tilde
(\textasciitilde), dem Komma und Leerzeichen sämtliche Zeichen
enthalten.

\DescribeMacro{\newsubexperiment} Für |\subexperiment| steht ein
analoges Makro, |\newsubexperiment|, zur Verfügung.

Bei jeder Verwendung eines dieser Makros wird "uberpr"uft, ob die
gew"unschte Abk"urzung bereits verwendet wurde, und gegebenenfalls ein
Fehler erzeugt. Gleichlautende Abk"urzungen f"ur |\experiment| und
|\subexperiment| sind grunds"atzlich m"oglich, d"urften aber nur
Verwirrung stiften (beim Benutzer, nicht bei \TeX). Wird das dritte
Argument dieser Makros leer gelassen, so wird die Langform auch f"ur
Index und Inhaltsverzeichnis verwendet.

|\experiment| und |\subexperiment| k"onnen wie |\section| mit oder
ohne optionales Argument verwendet werden, und zwar sowohl mit v"ollig
frei gew"ahlten "Uberschriften (Kurz- und Langform) als auch mit einer
Abk"urzung. Die Kombination einer Abk"urzung im optionalen Argument
und einer variierten, freien Langform im Hauptargument ist m"oglich,
die Verwendung eines frei formulierten optionalen Arguments mit einer
Abk"urzung im Pflichargument f"uhrt jedoch zu einem Fehler, denn die
Abk"urzung legt ja bereits einen Index- und Verzeichniseintrag fest.
Die Verwendung zweier unterschiedlicher Abk"urzungen in optionalem und
Pflichtargument ist m"oglich, wenn sie beide zum selben Indexeintrag
f"uhren.

\subsection{Mehrfache Indexeinträge}

Manchmal führt man zusammenpassende Arbeitschritte für verschiedene
Experimente parallel aus, das verkompliziert die Indexeinträge.
Stellen Sie sich z.\,B. vor, Sie hätten mit einer Screening-Methode
die Substanzen A152 und B96 aus einer kombinatorischen Bibliothek als
vielversprechende Medikamente identifiziert. Als nächstes möchten Sie
deren exakte Struktur, Zusammensetzung oder Sequenz verifizieren und
eine Präparationsmethode für etwas größere Mengen entwickeln.
Wahrscheinlich können Sie Zeit sparen, indem Sie einige dieser
Schritte für beide Substanzen parallel ausführen, aber dann erhalten
Sie einen Index mit Einträgen wie "`A152 und B96, Sequenzierung"' --
und zwei Monate später sollen Sie sich erinnern, ob sie B96 zusammen
mit A152 sequenziert haben oder stattdessen in der folgenden Woche mit
A43 und C12. Daher wäre es praktisch, für das Experiment
"`Sequenzierung von A152 und B96"' zwei Indexeinträge zu erhalten,
nämlich "A152, Sequenzierung"' und "B96, Sequenzierung"'. Und genau
das können Sie einfach erzeugen.

\DescribeMacro{\subexperiment}\DescribeMacro{\experiment}
Denn die Syntax von |\experiment| und |\subexperiment| erlaubt es,
eine durch Kommata abgetrennte Liste im optionalen Argument zu
verwenden. Das erste Element wird für das Inhaltsverzeichnis und die
Kolumnentitel verwendet, und die folgenden Elemente liefern
Indexeinträge. Angenommen, Sie haben die Abkürzungen A152-seq und
B96-seq definiert, dann brauchen Sie nur zu schreiben:
\begin{verbatim}
\experiment[Sequenzierung A152 und B96, A152-seq,
B96-seq]{Sequenzierung der Inhibitor-Kandidaten A152 und B96}
\end{verbatim}
und Sie erhalten was Sie wünschen. Leerzeichen vor und nach den
Kommata werden ignoriert. Sie können überall Abkürzungen oder freien
Text wählen, obwohl für die Indexeinträge wahrscheinlich nur
Abkürzungen wirklich sinnvoll sind.

\section{Implementierung}

Die Implementierung ist auf Englisch dokumentiert und kann durch
Übersetzung von |labbook.dtx| mit latex erhalten werden.

\PrintIndex

\end{document}
\endinput
%%