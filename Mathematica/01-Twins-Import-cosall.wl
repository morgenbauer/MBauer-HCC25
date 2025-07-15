(* ::Package:: *)

ClearAll@ImportPMC;
Options[ImportPMC] = {"Sheet" -> 1, "ColumnMap" -> {}};

ImportPMC::usage = "ImportPMC[fName] expects to import an Excel spreadsheet whose first tab is the data with the first line being headers";
ImportPMC[fPath_String, OptionsPattern[]]:=Module[{dirPath, tab, normalizedTab, sheetNames, headersFromExcel, expectedHeaders, missingHeaders, attributes, validOptSheet, dataDir},
	dataDir = FileNameTake[fPath, {1,-2}];
	(* Verify : file paths and file name *)
	If[!DirectoryQ[dataDir],
		Return[Failure["DirectoryNotFound",
		<|"Message"->"The specified directory \"" <> dataDir <> "\" does not exist"|>]]];
  
	If[!FileExistsQ[fPath],
		Return[Failure["FileNotFound",
		<|"Message"->"The specified file \""<>FileNameTake@fPath<>"\" was not found in the directory: " <> dataDir|>]]];
  
	(* Verify : Optional Sheets argument makes sense *)

	sheetNames = Import[fPath, {"XLSX","Sheets"}];
	(*TODO : Add logic to process one or more sheets names, ordinal values, etc.*)
	(* Verify : Optional Sheets argument makes sense *)

	tab = Import[fPath, {"XLSX","Tabular", OptionValue["Sheet"]},"HeaderLines"->1];
     
	expectedHeaders = {"PR", "Mesor", "Amplitude", "Tau", "PID", "Variable"};
	(*TODO : verify we have minimum set of required columns for each sheet*)

	normalizedTab=RenameColumns[tab, OptionValue["ColumnMap"]]
];

(* TESTING *)
$fName = "cosall.xlsx";
$dataPath = FileNameJoin[{NotebookDirectory[], "..", "GC", $fName}];
$tab=ImportPMC[$dataPath]


Export[FileNameJoin[{NotebookDirectory[],FileBaseName[$fName]<>".wxf"}] , $tab]
