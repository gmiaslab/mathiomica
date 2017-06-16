(* ::Package:: *)
(* Wolfram Language Package *)
(* Created by the Wolfram Workbench Nov 25, 2015 *)
(*The MIT License (MIT)

Copyright (c) 2016 George I. Mias, G. Mias Lab, Department of Biochemistry and Molecular Biology, Michigan State University, East Lansing 48824.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.*)
(* TODO Consider Mouse annotations as example *)
(* TODO Consider Drug enrichment analysis for KEGG *)
(* TODO Consider KEGG analysis for ath and mouse including KEGG and GO analyses and visualization - need UCSC gene tables*)
(* TODO Omics Tutorials, Guide, Overview at the very end *)
(* TODO generate example multi-omics data for KEGGPathwayVisual *)

BeginPackage["MathIOmica`",{"HierarchicalClustering`",
    "DatabaseLink`",
    "JLink`",
    "WebServices`",
    "RLink`"}]

(* ::Section:: *)
(*#####GlobalConstants#####*)

ConstantMathIOmicaDataDirectory::usage ="ConstantMathIOmicaDataDirectory is a global variable pointing to the MathIOmica data directory.";

ConstantMathIOmicaExamplesDirectory::usage ="ConstantMathIOmicaExamplesDirectory is a global variable pointing to the MathIOmica example data directory.";

ConstantMathIOmicaExampleVideosDirectory::usage ="ConstantMathIOmicaExampleVideosDirectory is a global variable pointing to the MathIOmica example videos directory."

ConstantGeneDictionary::usage ="ConstantGeneDictionary is a global gene/protein dictionary variable typically created by GetGeneDictionary."

(* ::Section:: *)
(*#####AnnotationsAndEnumerations#####*)

OBOGODictionary::usage = "OBOGODictionary[] is an Open Biomedical Ontologies (OBO) Gene Ontology (GO) vocabulary dictionary generator."

GOAnalysis::usage = "GOAnalysis[data] calculates input data over-representation analysis for Gene Ontology (GO) categories."

GOAnalysisAssigner::usage = "GOAnalysisAssigner[] creates identifier to gene ontology (GO) accession associations, restricted to required background set, downloading the data if necessary."

GetGeneDictionary::usage = "GetGeneDictionary[] creates an ID/accession dictionary from a UCSC search - typically of gene annotations."

GeneTranslation::usage = "GeneTranslation[inputList,targetIDList,geneDictionary] uses geneDictionary to convert inputList IDs to different annotations as indicated by targetIDList."

KEGGAnalysisAssigner::usage = "KEGGAnalysisAssigner[] creates KEGG: Kyoto Encyclopedia of Genes and Genomes pathway associations, restricted to required background set, downloading the data if necessary."

KEGGDictionary::usage = "KEGGDictionary[] creates a dictionary from KEGG: Kyoto Encyclopedia of Genes and Genomes terms - typically association of pathways and members therein."

KEGGAnalysis::usage = "KEGGAnalysis[data] calculates input data over-representation analysis for KEGG: Kyoto Encyclopedia of Genes and Genomes pathways."

MassMatcher::usage = "MassMatcher[data, accuracy] assigns putative mass identification to input data based on monoisotopic mass (using MathIOmica's mass dictionary), using the accuracy in parts per million."

MassDictionary::usage = "MassDictionary[] loads MathIOmica's current mass dictionary."

OmicsObjectUniqueMassConverter::usage="OmicsObjectUniqueMassConverter[omicsObject, massAccuracy] assigns a unique putative mass identification to each of omicsObject's inner association keys, using the massAccuracy in parts per million."; 

EnrichmentReportExport::usage = "EnrichmentReportExport[results] exports results from enrichment analyses to Excel spreadsheets."

(* ::Section:: *)
(*#####ClassificationAndClustering#####*)
MissingDataCreator::usage = "MissingDataCreator[data, setSamples] fills in Missing tags in a paired dataset for which the first component is not a member of a given sample list."

TimeSeriesClassification::usage = "TimeSeriesClassification[data, setTimes] identifies classes of temporal behavior in time series datasets that run over a specified set of time points."

TimeSeriesClusters::usage = "TimeSeriesClusters[data] performs clustering of time series data using two tiers of hierarchical clustering to identify groups and subgroups in the data.";

TimeSeriesSingleClusters::usage = "TimeSeriesSingleClusters[data] performs clustering of a time series data using a single tier of hierachrical clustering to identify groups in the data.";

MatrixClusters::usage = "MatrixClusters[data] performs hierarchical clustering in both dimensions of the provided data matrix.";


(* ::Section:: *)
(*#####Databases#####*)

UCSCBrowserSQL::usage = "UCSCBrowserSQL[query] performs a MySQL string query on the UCSC Genome Browser database tables."

(* ::Section:: *)
(*#####DataProcessing#####*)

OmicsObjectCreator::usage = "OmicsObjectCreator[outerLabels, innerLabels, measurements, metadata] creates an OmicsObject for use with MathIOmica."

FileSelector::usage = "FileSelector[variable] allows assignment of multiple file names and first lines to variable.";

DataImporter::usage = "DataImporter[associationName] provides a graphical interface to extract data and create OmicsObject variable for associations of information.";

DataImporterDirect::usage = "DataImporterDirect[positionsList, fileList, headerLines] creates an OmicsObject association by importing data from the files at the paths specified in the fileList, using the positionsList to select columns, and importing data by skipping a number of headerLines.";    

DataImporterDirectLabeled::usage = "DataImporterDirectLabeled[sampleRules,fileList,headerLines,headerColumnAssociations] creates an OmicsObject association for variableName, by imporing data from the files at the paths specified in fileList, using sampleRules as rules from labels to column-headers being imported for each file, and the headerColumnAssociations list of associations to associate column headers to column positions for each file.";
    
BootstrapGeneral::usage = "BootstrapGeneral[omicsObject, numberResampled] performs a resampling of the omicsObject data with replacement, and generates a new association structure with numbering corresponding to the numberResampled of new identities.";

QuantileEstimator::usage = "QuantileEstimator[data, timepoints] obtains the quantile estimator following bootstrap for time series.";

FilteringFunction::usage = "FilteringFunction[omicsObject, cutoff] filters OmicsObject data by a chosen comparison (by default greater or equal) to a cutoff.";

FilterMissing::usage = "FilterMissing[omicsObject, percentage] filters out data from omicsObject, retaining data across the datasets with a percentage of data points not missing."; 

QuantileNormalization::usage = "QuantileNormalization[data] performs quantile normalization of data.";

BoxCoxTransform::usage = "BoxCoxTransform[data,lambda] computes the Box-Cox transformation for a given parameter \[Lambda].";

ApplyBoxCoxTransform::usage = "ApplyBoxCoxTransform[data] for a given data set, computes the Box-Cox transformation at the maximum likelihood \[Lambda] parameter."

StandardizeExtended::usage = "StandardizeExtended[inputList, subtract, divide] allows standardization of data that may include Missing values with specified transformations."

Applier::usage = "Applier[function, inputData] applies function to OmicsObject, association or list inputData components."

ApplierList::usage = "ApplierList[function, inputData] applies applies function to list of lists from an association, nested association or components or a matrix inputData.";

Returner::usage = "Returner[originalAssociation, update] returns a modified originalAssociation updated at a specified position by the single association update, e.g. from Applier or ApplierList result.";
    
ConstantAssociator::usage = "ConstantAssociator[inputAssociation, associationAddition] adds multi-key constant to an OmicsObject (or an association of associations) inputAssociation, with each addition specified in a single association associationAddition, of form <|addition1-> Value1,addition2-> Value2,...|>.";

EnlargeInnerAssociation::usage = "EnlargeInnerAssociation[omicsObjectList] combines a list of OmicsObject (associations of associations)  omicsObjectList elements  by enlarging the inner associations - inner association Keys must be different.";

EnlargeOuterAssociation::usage = "EnlargeOuterAssociation[omicsObjectList] combines a list, omicsObjectList, of OmicsObject (or associations of associations) elements to a combined output by enlarging the outer associations - outer association keys must be different.";

LowValueTag::usage = "LowValueTag[omicsObject, valueCutoff] takes an omicsObject and tags values in specified position as Missing[] based on provided valueCutoff.";

MeasurementApplier::usage="MeasurementApplier[function,omicsObject] applies a function to the measurement list of an omicsObject, ignoring missing values.";

CreateTimeSeries::usage = "CreateTimeSeries[omicsObject] creates a time series list across an OmicsObject using outer keys as times.";
    
TimeExtractor::usage = "TimeExtractor[omicsObject] extracts a list of sorted times from an OmicsObject's outer keys.";
    
SeriesApplier::usage = "SeriesApplier[function,data] applies a given function to data, an association of lists, implementing masking for Missing values.";

SeriesCompare::usage = "SeriesCompare[series1, series2] merges the values of two associations of series (lists) by pointwise operation (by default subtraction) on the values of each matching pair of keys.";

SeriesInternalCompare::usage = "SeriesInternalCompare[associationOfLists] compares each value in each list of associationOfLists to an internal reference value in the list, if the reference point itself is not Missing.";

ConstantSeriesClean::usage = "ConstantSeriesClean[associationOfLists] removes constant list series from associationOfLists values.";

JoinNestedAssociations::usage = "JoinNestedAssociations[associationList] merges the nested associationList by joining the inner associations.";

(* ::Section:: *)
(*#####HypothesisTesting#####*)

BenjaminiHochbergFDR::usage = "BenjaminiHochbergFDR[pValues] calculates for a list of pValues the Benjamini Hochberg approach false discovery rates (FDR)."

(* ::Section:: *)
(*#####MassSpec#####*)

MSViewer::usage = "MSViewer[file] is a spectrum viewer for mass spectrometry .mzXML or .mzML data, in the provided file path."

(* ::Section:: *)
(*#####SpectralAnalysis#####*)

LombScargle::usage = "LombScargle[data, setTimes] calculates the Lomb-Scargle power spectrum for time series data that runs over specified setTimes.";
InverseAutocovariance::usage = "InverseAutocovariance[data, setTimes] calculates the inverse Fourier transform of the autocovariance for a time series that runs over specified set times, using a Lomb-Scargle approach.";
Autocorrelation::usage = "Autocorrelation[data, setTimes] calculates  the normalized autocorrelation for a time series over specified set times, using a Lomb-Scargle based inverse autocovariance.";

(* ::Section:: *)
(*#####Visualization#####*)

Heatmapper::usage = "Heatmapper[data] produces a heamap of the input data.";

TimeSeriesDendrogramHeatmap::usage = "TimeSeriesDendrogramHeatmap[data] generates a dendrogram and heatmap plot for time series data clusters."

TimeSeriesDendrogramsHeatmaps::usage = "TimeSeriesDendrogramsHeatmaps[data] generates a dendrogram and heatmap plots for all classified classes of time series data clusters."

TimeSeriesSingleDendrogramHeatmap::usage = "TimeSeriesSingleDendrogramHeatmap[data] generates a dendrogram and heatmap plot for time series data without subgrouping for classified time series data."

TimeSeriesSingleDendrogramsHeatmaps::usage = "TimeSeriesSingleDendrogramsHeatmaps[data] generate a dendrogram and heatmap plot for time series data without subgrouping for all classified time series data."

MatrixDendrogramHeatmap::usage = "MatrixDendrogramHeatmap[data] generates dendrograms and heatmap plots for a matrix."

MatrixDendrogramsHeatmaps::usage = "MatrixDendrogramsHeatmaps[data] generates dendrograms and heatmap plots for a matrix, all classifications."

KEGGPathwayVisual::usage = "KEGGPathwayVisual[pathway] generates a visual representation for a KEGG: Kyoto Encyclopedia of Genes and Genomes pathway."

(* ::Section:: *)
(* Function Options *)

AdditionalFilter::usage="AdditionalFilter is an option for various MathIOmica functions, that provides additional filtering that may be applied to the standard function output structure to be returned.";
(*GOAnalysis, KEGGAnalysis*)

AnalysisType::usage="AnalysisType is an option for various MathIOmica functions that specify different analysis methods that may be used.";
(*KEGGAnalysis, KEGGPathwayVisual*)

AppendString::usage="AppendString is an option for EnrichmentReportExport that will append a string to the exported file name following the class name. If a string is not provided, by default the current date will be appended.";
(*EnrichmentReportExport*)

AssociationLabels::usage="AssociationLabels is an option to add labels for the associations created by MathIOmica functions, such as DataImporter, i.e. first association keys in an OmicsObject."
(*DataImporter, DataImporterDirect, DataImporterDirectLabeled*)

AugmentDictionary::usage="AugmentDictionary is an option for MathIOmica functions, that provides a choice whether or not to augment the current ConstantGeneDictionary global variable or create a new one.";
(*GOAnalysis, KEGGAnalysis, KEGGPathwyVisual*)

AutocorrelationCutoffs::usage="AutocorrelationCutoffs is an option for TimeSeriesClassification that specifies a list of cutoffs Cutoffs, for \"Autocorrelation\" and \"InterpolatedAutocorrelation\" methods, for different lags that will be used to filter out data series for which the lags are not within cutoffs. The list length corresponds to cuttofs at different lags, with the ith lag cutoff provided as the ith index.";
(*TimeSeriesClassification*)

AutocorrelationLogic::usage="AutocorrelationLogic is an option for TimeSeriesClassification, which specifies whether or not to return the autocorrelation logic list for each signal, with the default set to False. If set to True, a logic vector is returned indicating whether or not at a particular lag the autocorrelation for a signal is above or below the AutocorrelationCutoffs.";
(*TimeSeriesClassification*)

AutocorrelationOptions::usage="AutocorrelationOptions is an option for MathIOmica functions that use internally Autocorrelation, and specifies a list of options to be passed to this internal Autocorrelation function.";
(*QuantileEstimator,TimeSeriesClassification*)

Averaging::usage="Averaging is an option for QuantileNormalization that specifies with what value to replace the same-rank elements across samples.";
(*QuantileNormalization*)

BackgroundSet::usage="BackgroundSet is an option for various MathIOmica functions, that involve considering pathways/groups/sets and that provides a list of IDs (e.g. gene accessions) that should be considered as the background for the calculation.";
(*GOAnalysis, GOAnalysisAssigner, KEGGAnalysis, KEGGAnalysisAssigner*)

BlendColors::usage="BlendColors is an option for KEGGPathwayVisual that provides a list of colors to be used in coloring intensities provided and is used by the IntensityFunction as its first argument. The colors must be provided as RGBColor[] specification.";
(*KEGGPathwayVisual*)

ClusterLabeling::usage="ClusterLabeling is an option for MathIOmica's clustering functions, specifies an additional label string to append to each cluster being computed, to prepend to the inbuilt cluster group labels."
(*MatrixClusters, TimeSeriesClusters,TimeSeriesSingleClusters*)

ColorBlending::usage="ColorBlending is an option for MathIOmica's heatmap generating functions, that provides a color scheme for the plot generated.";
(*Heatmapper,MatrixDendrogramHeatmap,TimeSeriesDendrogramHeatmap,TimeSeriesSingleDendrogramHeatmap*)

ColorSelection::usage="ColorSelection is an option for KEGGPathwayVisual that assigns foreground and background colors in the KEGG pathway through an association. The Keys point to labels for multi-omics data, and the values \"bg\" and \"fg\" can point to background and foreground representations respectively for each key.";
(*KEGGPathwayVisual*)

ColumnLabels::usage="ColumnLabels is an option for MatrixClusters that assigns labels to the columns of the input data matrix.";
(*MatrixClusters*)

CompareFunction::usage="CompareFunction is an option used by various MathIOmica comparison functions, that specifies a function used in implementing the comparisons.";
(*SeriesCompare, SeriesInternalCompare*)

ComparisonIndex::usage="ComparisonIndex is an option used by SeriesInternalCompare, that specifies the list position for the list value that will be used as a reference data point.";
(*SeriesInternalCompare*)

ComponentIndex::usage="ComponentIndex is an option for MathIOmica functions, such as Applier, that allows selection of which component of a list to use in an association or OmicsObject input or output values."
(*Applier, ApplyBoxCoxTransform, FilteringFunction, LowValueTag,MeasurementApplier,QuantileNormalization*)

DataTransforms::usage="DataTransforms is an option for MatrixClusters, that specifies transformations to be applied to the input matrix data for clustering, in the form {transformation for horizontal, transformation for vertical}.";
(*DataTransforms*)

DefaultColors::usage="DefaultColors is an option for KEGGPathwayVisual that provides a list of rules for setting the colors to be used as default values for the foreground \"fg\" and background \"bg\" respectively in the generated pathways. The colors must be provided as RGBColor[] specification.";
(*KEGGPathwayVisual*)

DeleteRule::usage="DeleteRule allows the customization of how to select values for the reference data point for which its key should be deleted.  The DeleteRule value takes the structure deleteRuleOptionValue={MatchQ first argument, MatchQ second argument}. The MatchQ function referred to here is implemented by SeriesInternalCompare internally, and uses the deleteRuleOptionValue as: MatchQ[deleteRuleOptionValue[[1]][reference comparison value], deleteRuleOptionValue[[2]]].  The default removes the corresponding key if the value used for reference in the comparison is actually Missing, i.e. the comparison reference point has Head that matches Missing.";
(*SeriesInternalCompare*)

DelimiterList::usage="DelimiterList is an option for DataImporter that allows the user to select the delimiters in importing text files.";
(*DataImporter*)

DendrogramColor::usage="DendrogramColor is an option for MathIOmica's dendrogram/heatmap generating functions, that provides color(s) option for highlighting the dendrogram(s).";
(*MatrixDendrogramHeatmatp,TimeSeriesDendrogramHeatmap,TimeSeriesSingleDendrogramHeatmap*)

DendrogramPlotOptions::usage="DendrogramPlotOptions is an option for MathIOmica functions that use internally the DendrogramPlot function, and provides a list of options to be passed to this internal DendrogramPlot function.";
(*MatrixClusters,TimeSeriesClusters,TimeSeriesSingleClusters*)

ExportMovieOptions::usage="ExportMovieOptions is an option for KEGGPathwayVisual that provides a list of options for the Export function used internally, to export the pathway list when Intensities have been provided for a time series representation of data.";
(*KEGGPathwayVisual*)

FileExtend::usage="FileExtend is an option for KEGGPathwayVisual that provides a string to be appended to the file name if the ResultsFormat is set to \"Movie\".";
(*KEGGPathwayVisual*)

FileURL::usage="FileURL is an option for OBOGODictionary that provides the location of the Open Biomedical Ontologies (OBO) Gene Ontology (GO) file in case this will be downloaded from the web.";
(*OBOGODictionary*)

FilterSignificant::usage="FilterSignificant is an option for MathIOmica functions can be set to True to filter data based on whether the analysis result is statistically significant, or if set to False to return all membership computations.";
(*KEGGAnalysis, GOAnalysis*)

FrameName::usage="FrameName is an option for various MathIOmica plots that specifies a label for the plot's frame.";
(*MatrixDendrogramHeatmap,TimeSeriesDendrogramHeatmap*)

FrequenciesOnly::usage="FrequenciesOnly is an option for LombScargle that specifies whether to return only the computation frequencies as an association."
(*LombScargle*)

FunctionOptions::usage="FunctionOptions is an option for MatrixDendrogramsHeatmaps and TimeSeriesDendrogramsHeatmaps that specifies a list of options to be passed to the internal MatrixDendrogramHeatmap and TimeSeriesDendrogramHeatmap respectively.";
(*MatrixDendrogramsHeatmaps,TimeSeriesDendrogramsHeatmaps,TimeSeriesSingleDendrogramsHeatmaps*)

GeneDictionary::usage="GeneDictionary is an option for MathIOmica functions, that points to an existing variable to use as a gene dictionary in annotations. If set to None the default ConstantGeneDictionary will be used.";
(*GOAnalysis, KEGGAnalysis,KEGGPathwayVisual*)

GetGeneDictionaryOptions::usage="GetGeneDictionaryOptions is an option for MathIOmica functions that use GetGeneDictionary internally, and that specifies a list of options that will be passed to this internal GetGeneDictionary function.";
(*GOAnalysis, KEGGAnalysis,KEGGPathwayVisual*)

GOAnalysisAssignerOptions::usage="GOAnalysisAssignerOptions is an option for GOAnalysis that specifies a list of options that will be passed to the internal GOAnalysisAssigner function.";
(*GOAnalysis*)

GOFileColumns::usage="GOFileColumns is an option for GOAnalysisAssigner that selects which columns to use for IDs and GO:accessions respectively from the downloaded GO annotation file, used when ImportDirectly is set to True to obtain a new GO association file.";
(*GOAnalysisAssigner*)

GOFileName::usage="GOFileName is an option for GOAnalysisAssigner that provides the name for the specific GO file to download from the GOURL if option ImportDirectly is set to True.";
(*GOAnalysisAssigner*)

GOURL::usage="GOURL is an option for GOAnalysisAssigner that provides the provides the location (base URL) where the GO association annotation files are downloaded from.";
(*GOAnalysisAssigner*)

GroupSubSize::usage="GroupSubSize is an option for TimeSeriesDendrogramHeatmap that provides a list to specify the relative size of group and subgroup reference columns in the plot.";
(*TimeSeriesDendrogramHeatmap*)

HorizontalAxisName::usage="HorizontalAxisName is an option for various MathIOmica plot generating functions that provides a label for a plot's horizontal axis.";
(*MatrixDendrogramHeatmap,TimeSeriesDendrogramHeatmap,TimeSeriesSingleDendrogramHeatmap*)

HorizontalLabels::usage="HorizontalLabels is an option for various MathIOmica heatmap generating functions, that provides labels for the horizontal axis, for each column in the heatmap.";
(*MatrixDendrogramHeatmap,TimeSeriesDendrogramHeatmap,TimeSeriesSingleDendrogramHeatmap*)

HorizontalSelection::usage="HorizontalSelection is an option for various MathIOmica functions, such as Applier, that allows for horizontal selection across components for a single level association with multi-list values."
(*Applier, ApplyBoxCoxTransform*)

HypothesisFunction::usage="HypothesisFunction is an option for various MathIOmica functions, that allows the choice of function for implementing multiple hypothesis testing considerations."
(*GOAnalysis, KEGGAnalysis*)

IgnorePattern::usage="IgnorePattern is an option for MeasurementApplier specifying a pattern of values to delete prior to applying the function to the measurement list."
(*MeasurementApplier*)

ImportDirectly::usage="ImportDirectly is an option for various MathIOmica functions to import various data directly from the internet.";
(*GetGeneDictionary, GOAnalysisAssigner,KEGGAnalysisAssigner, KEGGDictionary, OBOGODictionary*)

IndexColor::usage="IndexColor is an option for various MathIOmica dendrogram/heatmap generating plots that provides a color gradient scheme for the index color legends that indicate groupings from the clustering.";
(*MatrixDendrogramHeatmap,TimeSeriesDendrogramHeatmap,TimeSeriesSingleDendrogramHeatmap*)

InputID::usage="InputID is an option for various MathIOmica functions that specifies the kind of identifiers/accessions used as input.";
(*GetGeneDictionary, GeneTranslation, GOAnalysis, KEGGAnalysis,KEGGPathwayVisual*)

Intensities::usage="Intensities is an option for KEGGPathwayVisual that specifies may be used to provide a set of intensities that will be used for coloring components of the pathway. The intensities are provided as an association for each ID as single values, or as a list of values in the case of series data. Intensities must be scaled from -1 to 1, or selected such that the IntensityFunction can convert them to a number between 0 to 1.";
(*KEGGPathwayVisual*)

IntensityFunction::usage="IntensityFunction is an option for KEGGPathwayVisual that specifies a function of two arguments that allows customization of the coloring for the intensities. The IntensityFunction value can be any function which outputs a color, I(#1,#2), (where #1 is the BlendColors option value, and #2 is an intensity vector, that has values typically ranging from [-1,1])."
(*KEGGPathwayVisual*)

InterpolationDeltaT::usage="InterpolationDeltaT is an option for various MathIOmica functions that utilize an interpolation function, and specifies the time step used to grid the time window over which calculations will be performed.";
(*QuantileEstimator,TimeSeriesClassification*)

InterpolationOptions::usage="InterpolationOptions is an option for various MathIOmica functions that utilize an internal Interpolation function, and specifies a list of options to be passed to this internal Interpolation function";
(*QuantileEstimator,TimeSeriesClassification*)

InverseSelection::usage="InverseSelection is an option for ConstantSeriesClean that can be set to True to invert the selection/filtering process and return the constant series instead of the non-constant ones.";
(*ConstantSeriesClean*)

JavaGBs::usage="JavaGBs is an option for MathIOmica functions that utilize Java to assign much memory to use for Java in GBs."
(*GetGeneDictionary, DataImporter*)

KeyModifiers::usage="KeyModifiers is an option for EnlargeInnerAssociation that selects whether keys will be modified.";
(*EnlargeinnerAssociation*)

KEGGAnalysisAssignerOptions::usage="KEGGAnalysisAssignerOptions is an option for MathIOmica functions that use internally KEGGAnalysisAssigner, and specifies a list of options that will be passed to this internal KEGGAnalysisAssigner function."
(*KEGGAnalysis, KEGGAnalysisVisual*)

KEGGDatabase::usage="KEGGDatabase is an option for MathIOmica functions that use the KEGG database, and indicates which KEGG database to use as the target database."
(*KEGGAnalysis,KEGGAnalysisVisual*)

KEGGDictionaryOptions::usage="KEGGDictionaryOptions is an option for MathIOmica functions that utilize KEGGDictionary internally, and specifies a list of options to be passed to the internal KEGGDictionary function that provides the KEGG annotations.";
(*KEGGAnalysis*)

KEGGDictionaryVariable::usage="KEGGDictionaryVariable is an option for MathIOmica functions that use a KEGG dictionary, and provides a KEGG annotation variable. If set to None, KEGGDictionary will be used internally to automatically generate the default KEGG annotation.";
(*KEGGAnalysis*)

KEGGMolecular::usage="KEGGMolecular is an option for various MathIOmica functions that perform molecular analysis using the KEGG database and specifies which database to use for molecular analysis (default is the compound database: cpd).";
(*KEGGAnalysis,KEGGAnalysisVisual*)

KEGGOrganism::usage="KEGGOrganism is an option for various MathIOmica functions that perform analysis using the KEGG database and indicates which organism (org) to use for \"Genomic\" type of analysis (default is human analysis: org=\"hsa\").";
(*KEGGAnalysis,KEGGAnalysisVisual*)

KEGGUCSCSplit::usage="KEGGUCSCSplit is an option for various MathIOmica functions that specifies a two component list, {True/False, label}. If the first component is set to True the initially imported KEGG IDs, identified by the second component label,  are split on + string to fix nomenclature issues, retaining the string following +.";
(*GetGeneDictionary*)

KEGGQuery1::usage="KEGGQuery1 is an option for various MathIOmica functions that make KEGG API calls, and sets string query1 in http://rest.kegg.jp/link/<> query1 <> / <> query2. Typically this will be used as the target database to find related entries by using database cross-references.";
(*KEGGAnalysisAssigner, KEGGDictionary*)

KEGGQuery2::usage="KEGGQuery2 is an option for various MathIOmica functions that make KEGG API calls, and sets string query2 in http://rest.kegg.jp/link/<> query1 <> / <> query2. Typically this will be used as the source database to find related entries by using database cross-references.";
(*KEGGAnalysisAssigner, KEGGDictionary*)

Labels::usage="Labels is an option for various MathIOmica functions that provides a string list for how keys in a created association will be named."; 
(*KEGGAnalysisAssigner*)

LabelFunction::usage="LabelFunction is an option for BootstrapGeneral, indicating which function to use to generate the labels for the simulated data.";
(*BootstrapGeneral*)

LengthFilter::usage="LengthFilter is an option for MathIOmica functions that perform computations of membership in pathways/ontologies/groups/sets, that allows the selection of how many members each category can have, as typically restricted by the LengthFilterFunction."
(*GOAnalysisAssigner, KEGGAnalysisAssigner*)

LengthFilterFunction::usage="LengthFilterFunction is an option for MathIOmica functions that perform computations of membership in pathways/ontologies/groups/sets, that specifies which function to use to filter the number of members a reported category has compared to the number typically provided by LengthFilter.";
(*GOAnalysisAssigner, KEGGAnalysisAssigner*)

LinkageMeasure::usage="LinkageMeasure is an option for MathIOmica's clustering functions that specifies which linkage measure(s) to use in computing fusion coefficients.";
(*MatrixClusters,TimeSeriesClusters,TimeSeriesSingleClusters*)

ListIndex::usage="ListIndex is an option for MathIOmica functions, such as Applier that allows selection of which list to use in the association or OmicsObject input or output values."
(*Applier, ApplyBoxCoxTransform, FilteringFunction, LowValueTag,MeasurementApplier,QuantileNormalization*)

LombScargleCutoff::usage="LombScargleCutoff is an option for TimeSeriesClassification that provides a cutoff value for the \"LombScargle\" Method, for identifying the highest intensity observed in the power spectrum that is greater than this value."
(*TimeSeriesClassification*)

LombScargleOptions::usage="LombScargleOptions is an option for various MathIOmica functions that utilize an internal LombScargle function, and specifies a list of options to be passed to this internal LombScargle function";
(*QuantileEstimator,TimeSeriesClassification*)

MassDictionaryVariable::usage="MassDictionaryVariable is an option for MassMatcher that can provide a mass dictionary variable. If set to None, MathIOmica's inbuilt mass dictionary (MassDictionary) will be loaded and used.";
(*MassMatcher*)

MassMatcherOptions::usage="MassMatcherOptions is an option for MathIOmica functions that internally use MassMatcher, and provides a list of options to be passed to this internal MassMatcher function.";
(*OmicsObjectUniqueMassConverter*)

MathIOmicaDataDirectory::usage="MathIOmicaDataDirectory is an option for various MathIOmica functions option specifying the directory where the default MathIOmica package data is stored.";
(*GetGeneDictionary, GOAnalysisAssigner, KEGGAnalysis, KEGGAnalysisAssigner, KEGGAnalysisVisual,MSViewer,MassDictionary,MassMatcher,OBOGODictionary*)

MemberSet::usage="MemberSet is an option for KEGGAnalysisVisual that selects which members of the pathway are to be considered. The choices are:
All: return the pathway only.
{list of identifiers}: a list of identifiers that will be highlighted. If ORA is set to True the list must be the output from an over representation analysis, and the identifiers will be selected from the last list, second sublist.
Only IDs that are found to match in the pathway are colored.
An internal gene dictionary (see GetGeneDictionary) is used to convert IDs to KEGG IDs."
(*KEGGAnalysisVisual*)

MergeFunction::usage="MergeFunction is an option for EnlargeInnerAssociation to use in the merging of the input omicObjectList  OmicsObject variables' inner associations.";
(*EnlargeInnerAssociation*)

MinimumPoints::usage="MinimumPoints is an option for FilterMissing to select the minimum number of datapoints to keep.";
(*FilterMissing*)

MissingMask::usage="MissingMask is an option for SeriesApplier that specifies the value to be used to mask Missing data.";
(*SeriesApplier*)

MissingReplacement::usage="MissingReplacement is an option for functions that handle missing data (e.g. EnlargeInnerAssociation) that indicates what value missing data should be globally assigned.";
(*EnlargeInnerAssociation, EnlargeOuterAssociation*)

MissingSubstitution::usage="MissingSubstitution is an option for MathIOmica's clustering functions that specifies a substitution rule for Missing (or other pattern) of data.";
(*MatrixClusters*)

MissingValueColor::usage="MissingValueColor is an option for KEGGPathwayVisual, that specifies a color to be used when Intensities are provided to represent values that are tagged as Missing[]. The color must be provided as RGBColor[] specification.";
(*KEGGPathwayVisual*)

MolecularInputID::usage="MolecularInputID is an option for MathIOmica functions that perform molecular analyses, and specifies a string list to indicate the kind of ID to use for the input molecule entries.";
(*KEGGAnalysis,KEGGPathwayVisual*)

MolecularOutputID::usage="MolecularOutputID is an option for MathIOmica functions that perform molecular analyses, and provides a string a string to indicate the kind of ID to convert input molecule entries. The default is typically cpd, consistent with use of the cpd database as the default molecular analysis.";
(*KEGGAnalysis,KEGGPathwayVisual*)

MolecularSpecies::usage="MolecularSpecies is an option for MathIOmica functions that perform molecular analyses, and specifies the kind of molecular input.";
(*KEGGAnalysis,KEGGPathwayVisual,MassDictionary, MassMatcher*)

MovieFilePath::usage="MovieFilePath is an option for KEGGPathwayVisual that indicates the path (including file name) where, if ResultsFormat is set to \"Movie\"  the movie generated will be saved. The default value None will generate a file named after the pathway with extension set by the FileExtend option in the current directory.";
(*KEGGPathwayVisual*)

MultipleList::usage="MultipleList is an option for MathIOmica functions that have IDs as inputs, and that specifies whether the input accessions list constituted a multi-omics list input that is annotated so."
(*GOAnalysis, KEGGAnalysis*)

MultipleListCorrection::usage="MultipleListCorrection is an option for MathIOmica functions that have IDs as inputs that may come from multi-omics inputs, and that specifies whether or not to correct for multi-omics analysis. The choices are None, Automatic, or a custom number.";
(*GOAnalysis, KEGGAnalysis*)

NonUCSC::usage="NonUCSC is an option for MathIOmica functions that use KEGG identifiers in dictionaries, and specifies if UCSC browser was used in determining an internal GeneDictionary used in ID translations, where the KEGG identifiers for genes are number strings (e.g. 4790).The NonUCSC option can be set to True if standard KEGG accessions are used in a user provided GeneDictionary variable, in the form OptionValue[KEGGOrganism] <>:<>numberString, e.g. hsa:4790.";
(*KEGGAnalysis,KEGGPathwayVisual*)

NormalizeIntensities::usage="NormalizeIntensities is an option for LombScarle that specifies whether the intensities list should be normalized or not.";
(*LombScargle*)

OBODictionaryVariable::usage="OBODictionaryVariable is an option for GOAnalysis that can provide a GO annotation variable. If set to None, OBOGODictionary will be used internally to automatically generate the default GO annotation.";
(*GOAnalysis*)

OBOFile::usage="OBOFile is an option for OBOGODictionary that specifies the name of the file that will be used to construct the ontology dictionary, or if not present saved when dowloaded from the web.";
(*OBOGODictionary*)

OBOGODictionaryOptions::usage="OBOGODictionaryOptions is an option for GOAnalysis that specifies a list of options to be passed to the internal OBOGODictionary function that provides the GO annotations.";
(*GOAnalysis*)

OntologyLengthFilter::usage="OntologyLengthFilter is an option for GOAnalysis that can be used to set the value for which terms to consider in the computation, by excluding GO terms that have fewer items compared to the OntologyLengthFilter value. It is used by the internal GOAnalysisAssigner function.";
(*GOAnalysis*)

ORA::usage="ORA is an option for KEGGPathwayVisual, that can be set to True or False depending on whether the input is from an over representation analysis (e.g. output from KEGGAnalysis), or not respectively."; 
(*KEGGPathwayVisual*)

OtherReplacement="OtherReplacement is an option for LowValueTag that provides a replacement rule for any other kind of replacement in the data.";
(*LowValueTag*)

OutputDirectory::usage="OutputDirectory is an option for EnrichmentReportExport that specifies the location of a directory to output the Excel spreadsheets generated. If it is set to None the NotebookDirectory[] will be used as a default output directory."
(*EnrichmentReportExport*)

OutputID::usage="OutputID is an option for various MathIOmica functions that use input IDs, and takes a string value that specifies what kind of IDs/accessions to convert the input IDs/accession numbers in the function's analysis."
(*GOAnalysis, KEGGAnalysis,KEGGPathwayVisual*)

OversamplingRate::usage="OversamplingRate is an option for LombScargle, that specifies the rate at which to oversample the input time series using zero-padding.";
(*LombScarlge*)

PairReturn::usage="PairReturn is an option for InverseAutocorrelation and LombScargle that specifies whether paired data should be returned as two lists or as a list of pairs.";
(*Inverse Autocorrelation, LombScargle*)

PathwayLengthFilter::usage="PathwayLengthFilter is an option for MathIOmica functions that use pathways, and specifies which pathways to consider in the computation, by excluding pathways that have fewer items compared to the PathwayLengthFilter value.";
(*KEGGAnalysis,KEGGAnalysis*)

PrintDendrograms::usage="PrintDendrgorams is an option for MathIOmica's clustering functions that indicates whether or not to print dendrograms for the clustering computed.";
(*MatrixClusters,TimeSeriesClusters,TimeSeriesSingleClusters*)

PSImsOBOFile::usage="PSImsOBOFile is an option for MSViewer that points to an Open Biomedical Ontologies (OBO) file for the mass spectromety information found in the ConstantMathIOmicaDataDirectory.";
(*MSViewer*)

pValueCutoff::usage="pValueCutoff is an option for various MathIOmica functions that utilize p-Values, that provides a cutoff p-value for (adjusted) p-values to assess statistical significance.";
(*GOAnalysis, KEGGAnalysis*)

QuantileValue::usage="QuantileValue is an option for QuantileEstimator that specifies which quantile to extract."; 
(*QuantileEstimator*)

Reference::usage="Reference is an option for FilterMissing that allows selection of a reference outer key for which to remove a dataset if the reference point has Missing values.";
(*FilterMissing*)

ReportFilter::usage="ReportFilter is an option for MathIOmica functions that use pathways/ontologies/groups, and provides a cutoff for membership in ontologies/pathways/groups in selecting which terms/categories to return. It is typically used in conjunction with ReportFilterFunction.";
(*GOAnalysis, KEGGAnalysis*)

ReportFilterFunction::usage="ReportFilterFunction is an option for GOAnalysis that specifies what operator form will be used to compare against ReportFilter option value in selecting which terms/categories to return."; 
(*GOAnalysis*)

ResultsFormat::usage="ResultsFormat is an option for KEGGPathwayVisual that provides provides a choice of output format, the choices are:
\"URL\": returns a URL of the pathway,
\"Figure\": returns figure output(s) for the pathway,
\"Movie\": in the case of series data returns a movie/animation of the series pathway snapshots.";
(*KEGGPathwayVisual*)

ReturnAllSpikes::usage="ReturnAllSpikes is an option for TimeSeriesClassification that specifies whether each signal may maintain unique membership to each spike class, or be allowed to belong to multiple classes. The option is used in \"Autocorrelation\" and \"InterpolatedAutocorrelation\"  Method selections. If set to False, first spike maxima are classified, and only signals found not to belong to spike maxima are then considered for membership in the spike minima class."
(*TimeSeriesClassification*)

ReturnData::usage="ReturnData is an option for TimeSeriesClassification that can take values:
True: the function will return input keys to data associations in the classification.
False: the function will only return the keys of the input data in the classification."
(*TimeSeriesClassification*)

ReturnDendrograms::usage="ReturnDendrograms is an option for MathIOmica's clustering functions that specifies whether or not to to return the dendrograms as output.";
(*MatrixClusters,TimeSeriesClusters,TimeSeriesSingleClusters*)

ReturnDropped::usage="ReturnDropped is an option for ConstantSeriesClean that can be set to True to return the keys of the constant series in addition to the filtered list. The data is returned in an association:
<|Data->non-constant data,Dropped keys-> keys of dropped values |>";
(*ConstantSeriesClean*)

ReturnModels::usage="ReturnModels is an option for TimeSeriesClassification, that indicates whether or not to return the models as well as the classification information for the input data. The data is returned as an association with the  key \"TimeSeriesClasses\" for classification groups and one of the following: (i) \"Models\" for model-based Method, (ii) \"LombScargle\" for periodograms in the \"LombScargle\" Method, (iii) \"Autocorrelations\" for any autocorrelation based Method."
(*TimeSeriesClassification*)

SamplingFunction::usage="SamplingFunction is an option for BootstrapGeneral for selecting a sampling strategy function.";
(*BootstrapGeneral*)

SampleKind::usage="SampleKind is an option for DataImporterDirect and DataImporterDirectLabeled, that assigns a label for the sample kind, e.g. {\"Protein\"}."
(*DataImporterDirect,DataImporterDirectLabeled*)

ScaleShift::usage="ScaleShift is an option for various MathIOmica plot generating functions to reset the blend of the colors used overall. The option is a real positive number, and is used as a multiplier for the internal Blend function's second argument.";
(*Heatmapper, MatrixDendrogramHeatmap,TimeSeriesDendrogramHeatmap,TimeSeriesSingleDendrogramHeatmap*)

SelectionFunction::usage="SelectionFunction is an option for FilteringFunction to select which function will be used in filtering.";
(*FilteringFunction*)

ShowPlots::usage="ShowPlots is an option for FilterMissing to select whether or not to show summary plots.";
(*FilterMissing*)

SignificanceCriterion::usage="SignificanceCriterion is an option for MathIOmica's clustering functions that specifies the method used in determining the number of groups in the clustering.";
(*MatrixClusters,TimeSeriesClusters,TimeSeriesSingleClusters*)

SimilarityVectors::usage="SimilarityVectors is an option for MatrixClusters that can provide a list of similarity vectors to be used for clustering.";
(*MatrixClusters*)

SingleAssociationLabel::usage="SingleAssociationLabel is an option for MathIOmica's clustering functions that specifies a label to use in case a list is provided to name the class of data produced.";
(*MatrixClusters,TimeSeriesClusters,TimeSeriesSingleClusters*)

SingleColorPlace::usage="SingleColorPlace is an option for KEGGPathwayVisual that selects in the case of a single identifier input whether to place the color to the foreground, (\"fg\") or background (\"bg\" set by default)."; 
(*KEGGPathwayVisual*)

Species::usage="Species is an option for various MathIOmica functions, that specifies the species considered in the calculation, by default corresponding to human.";
(*GeneTranslation, GetGeneDictionary, GOAnalysis, GOAnalysisAssigner,KEGGPathwayVisual*)

SpectrumFunction::usage="SpectrumFunction is an option for Autocorrelation that allows selection of the specific function to use to generate a periodogram on which inverse Fourier transform will be performed to obtain the autocorrelation function."
(*Autocorrelation*)

SpikeCutoffs::usage="SpikeCutoffs is an option for TimeSeriesClassification that provides an association with number, n, of data points as keys, and values corresponding to cutoffs, in the form <|n-> {Maximum Spike nth Cutoff,Minimum Spike nth Cutoff}|> used to call spike maxima and minima for a time series with n number of datapoints. The values are provided by the user - the default values are only place-holders and should be replaced by actual values. The association must have corresponding keys for all possible lengths of input datasets. i.e. all possible lengths of series constructed by excluding Missing or other non-numeric values."; 
(*TimeSeriesClassification*)

StandardHighlight::usage="StandardHighlight is an option for KEGGPathwayVisual that provides a list of rules for setting the highlight colors for the IDs represented in the pathway (when no intensities are provided). The list specifies color rules for foregroung, \"fg\", and background, \"bg\", respectively. The colors must be provided as RGBColor[] specification.";
(*KEGGPathwayVisual*)

SubclusteringDistanceFunction::usage="SubclusteringDistanceFunction is an option for TimeSeriesClusters that specifies which DistanceFunction will be used in calculating the similarities between different time series in the second tier of clustering.";
(*TimeSeriesClusters*)

TestFunction::usage="TestFunction is an option for MathIOmica functions that perform statistical tests, and provides a function used to calculate p-values.";
(*GOAnalysis, KEGGAnalysis*)

Ties::usage="Ties is an option for QuantileNormalization that indicates how ties should be handled.";
(*QuantileNormalization*)

UCSCSQLString::usage="UCSCSQLString is an option for MathIOmica functions, such as GetGeneDictionary, that use a MySQL connection to the UCSC browser, providing an association to be used to obtain data from the UCSC Browser tables. The key of the association must match the Species option value used (default: human). The value for the species corresponds to the actual MySQL command used.";
(*GetGeneDictionary*)

UCSCSQLSelectLabels::usage="UCSCSQLSelectLabels is an option for GetGeneDictionary, that provides an association to be used to assign key labels for the data improted from the UCSC Browser tables. The key of the association must match the Species option value used (default: human). The value is a multi component string list corresponding to the matrices in the MathIOmica data file, or the tables used in the MySQL query provided by UCSCSQLString.";
(*GetGeneDictionary*)

UpperFrequencyFactor::usage="UpperFrequencyFactor is an option for LombScargle and InverseAutocovariance for scaling the upper Nyquist cutoff frequency."
(*Autocorrelation, InverseAutocorrelation, LombScargle*)

ValueReplacement::usage="ValueReplacement is an option for LowValueTag that specifies a value for replacement of tagged data points."
(*LowValueTag*)

VerticalLabels::usage="VerticalLabels is an option for various MathIOmica heatmap generating plots, that specifies labels to be used for vertical axis labeling of each row.";
(*MatrixDendrogramHeatmap,TimeSeriesDendrogramHeatmap,TimeSeriesSingleDendrogramHeatmap*)

VocabularyVariable::usage="VocabularyVariable is an option for MSViewer that provides a vocabulary variable to use for mass spectrometry vocabularies in lieu of using the PSImsOBOFile option. The content should have the form of an OBO file imported as lines: Import[filename,\"Lines\"].";
(*MSViewer*)

Begin["`Private`"]
(* ::Section:: *)
(*#####Intro#####*)
Print["MathIOmica (", Hyperlink["http://mathiomica.org"], "),", 
 Style[" by ", Italic], 
 Hyperlink[Style["G. Mias Lab", Italic], 
  "http://georgemias.org"], "."];

(* ::Section:: *)
(*#####GlobalConstants#####*)

(* definition of ConstantMathIOmicaDataDirectory *)
ConstantMathIOmicaDataDirectory = If[ ! DirectoryQ[
       FileNameJoin[
        Flatten[{FileNameSplit[$UserBaseDirectory], "Applications", 
          "MathIOmica", "MathIOmicaData"}]]],
                               CreateDirectory[
                                FileNameJoin[
                                 Flatten[{FileNameSplit[$UserBaseDirectory], "Applications", 
                                   "MathIOmica", "MathIOmicaData"}]]],
                               FileNameJoin[
                                Flatten[{FileNameSplit[$UserBaseDirectory], "Applications", 
                                  "MathIOmica", "MathIOmicaData"}]]
                           ];
(* definition of ConstantMathIOmicaExamplesDirectory *)

ConstantMathIOmicaExamplesDirectory=If[ ! DirectoryQ[
       FileNameJoin[
        Flatten[{FileNameSplit[$UserBaseDirectory], "Applications", 
          "MathIOmica", "MathIOmicaData","ExampleData"}]]],
                                   CreateDirectory[
                                    FileNameJoin[
                                     Flatten[{FileNameSplit[$UserBaseDirectory], "Applications", 
                                       "MathIOmica", "MathIOmicaData","ExampleData"}]]],
                                   FileNameJoin[
                                    Flatten[{FileNameSplit[$UserBaseDirectory], "Applications", 
                                      "MathIOmica", "MathIOmicaData","ExampleData"}]]
                               ];

ConstantMathIOmicaExampleVideosDirectory = If[ ! DirectoryQ[
       FileNameJoin[
        Flatten[{FileNameSplit[$UserBaseDirectory], "Applications", 
          "MathIOmica", "MathIOmicaData","ExampleVideos"}]]],
                               CreateDirectory[
                                FileNameJoin[
                                 Flatten[{FileNameSplit[$UserBaseDirectory], "Applications", 
                                   "MathIOmica", "MathIOmicaData","ExampleVideos"}]]],
                               FileNameJoin[
                                Flatten[{FileNameSplit[$UserBaseDirectory], "Applications", 
                                  "MathIOmica", "MathIOmicaData","ExampleVideos"}]]
                           ];
    
(* ::Section:: *)
(*#####AnnotationsAndEnumerations#####*)

(****Gene Ontology****)

(* ::Function:: *)
(* f:OBOGODictionary *)
(***Options***)
Options[OBOGODictionary] = {FileURL-> "http://purl.obolibrary.org/obo/go/go-basic.obo",
    ImportDirectly -> False,
    MathIOmicaDataDirectory -> 
    If[ ! DirectoryQ[
       FileNameJoin[
        Flatten[{FileNameSplit[$UserBaseDirectory], "Applications", 
          "MathIOmica", "MathIOmicaData"}]]],
        CreateDirectory[
         FileNameJoin[
          Flatten[{FileNameSplit[$UserBaseDirectory], "Applications", 
            "MathIOmica", "MathIOmicaData"}]]],
        FileNameJoin[
         Flatten[{FileNameSplit[$UserBaseDirectory], "Applications", 
           "MathIOmica", "MathIOmicaData"}]]
    ],
    OBOFile -> "goBasicObo.txt"};
(***Function***)
OBOGODictionary[OptionsPattern[]] :=
    Module[ {dir = FileNameSplit[OptionValue[MathIOmicaDataDirectory]], 
      importTrueOrFalse = OptionValue[ImportDirectly],fileURL = OptionValue[FileURL],
      oboFile = OptionValue[OBOFile], fileGOOBO, tempPrint, inputFile, 
      outDictionary},
     (*import the GO OBO file*)
        fileGOOBO = FileNameJoin[Flatten[{dir, oboFile}]];
        (*we check if the OBO file Exist, if not, 
        attempt to download and create it*)
        If[ !FileExistsQ[fileGOOBO],
            tempPrint = 
             PrintTemporary[
              "Did Not Find Annotation Files, Attempting to Download..."];
            importTrueOrFalse = True
        ];
        If[ importTrueOrFalse,
            If[ FileExistsQ[fileGOOBO],
                DeleteFile[fileGOOBO]
            ];(*deleted existing file*)
            URLSave[fileURL, 
             FileNameJoin[Flatten[{dir, oboFile}]]];
            (*cleaned up archive*)
            NotebookDelete[tempPrint];
            If[ FileExistsQ[fileGOOBO],
                Print["Created Annotation Files at ", fileGOOBO],
                Return["Did Not Find Annotation Files, Aborting Process"]
            ]
        ];
        inputFile = Import[fileGOOBO, "Lines"];
        (*find keys "accessions (id):" and "name:" and "namespace" but \
      extract their corresponding values in a list and map them to their \
      corresponding [Term] positions, 
        once the "accessions (id):" and its corresponding "name:" in a \
      list, make an association between them so you can serch this \
      association using the key "accessions (id):" to get the value "name:" \
      and "namespace"*)
        outDictionary = 
         Association[#[[
              1]] -> {#[[2]], #[[3]]} & /@ (StringTrim[inputFile[[#]], 
               "id: " | "name: " | 
                "namespace: "] & /@ ({#[[1]] + 1, #[[1]] + 2, #[[1]] + 
                  3} & /@ Position[inputFile, "[Term]"]))];
        Return[outDictionary]
    ];

(* ::Function:: *)
(* f:GOAnalysisAssigner *)
(*Function will download and \
create gene associations and restrict to required background set*)
(***Options***)
Options[GOAnalysisAssigner] = {BackgroundSet -> All,
   GOFileColumns -> {2, 5},
   GOFileName -> None,
   GOURL-> "http://geneontology.org/gene-associations/",
   ImportDirectly -> False, 
   LengthFilter -> None, 
   LengthFilterFunction -> GreaterEqualThan,
   MathIOmicaDataDirectory -> 
   If[ ! DirectoryQ[
      FileNameJoin[
       Flatten[{FileNameSplit[$UserBaseDirectory], "Applications", 
         "MathIOmica", "MathIOmicaData"}]]],
       CreateDirectory[
        FileNameJoin[
         Flatten[{FileNameSplit[$UserBaseDirectory], "Applications", 
           "MathIOmica", "MathIOmicaData"}]]],
       FileNameJoin[
        Flatten[{FileNameSplit[$UserBaseDirectory], "Applications", 
          "MathIOmica", "MathIOmicaData"}]]
   ],
   Species -> "human"};
(***Function***)
GOAnalysisAssigner[OptionsPattern[]] :=
    Module[ {dir = FileNameSplit[OptionValue[MathIOmicaDataDirectory]], 
      importTrueOrFalse = OptionValue[ImportDirectly],
      backGeneSet = OptionValue[BackgroundSet], 
      species = OptionValue[Species],
      lengthFilter = OptionValue[LengthFilter], 
      lengthFilterFunction = OptionValue[LengthFilterFunction],
      goFileName = OptionValue[GOFileName],
      goColumns = OptionValue[GOFileColumns],
      goURL = OptionValue[GOURL],
      file,
      localFile,
      localZipFile,
      fileGOAssociations, goData, identifierAssoc, geneOntAssoc, 
      returning, tempPrint},
     (*local variables*)(*if the user asked us to import directly,
     import directly from GO website,otherwise,
     get it from a directory they specify*)
        file = 
         If[ MatchQ[goFileName, None],
             "gene_association.goa_" <> species <> ".gz",
             goFileName
         ];
        localFile = 
         FileNameJoin[Flatten[{dir, "gene_association.goa_" <> species}]];
        localZipFile = 
         FileNameJoin[
          Flatten[{dir, "gene_association.goa_" <> species <> ".gz"}]];
        fileGOAssociations = 
         FileNameJoin[Flatten[{dir, #}]] & /@ {species <> "GeneOntAssoc", 
           species <> "IdentifierAssoc"};
        (*we check if the Annotations Exist, if not, 
        attempt to download and create them*)
        If[ !(And @@ (FileExistsQ[#] & /@ fileGOAssociations)),
            tempPrint = 
             PrintTemporary[
              "Did Not Find Annotation Files, Attempting to Download..."];
            importTrueOrFalse = True
        ];
        If[ importTrueOrFalse,
            If[ FileExistsQ[localFile],
                DeleteFile[localFile]
            ];(*deleted existing file*)
            URLSave[goURL <> file, 
             FileNameJoin[localZipFile]];
            ExtractArchive[localZipFile, DirectoryName@localFile];
            DeleteFile[localZipFile];
            (*cleaned up archive*)
            goData = 
             DeleteCases[
              Import[localFile, 
               "TSV"], _?(StringMatchQ[#[[1]], "!" ~~ ___] &)];
            (*remove comments by "!" lines*)
            geneOntAssoc = 
             Association@(#[[1, 2]] -> (Union@ #[[All, 1]]) & /@ 
                Query[All /* (GatherBy[#, Last] &), goColumns]@goData);
            Put[Join[Date[], {geneOntAssoc}], fileGOAssociations[[1]]];
            identifierAssoc = 
             Association@(#[[1, 1]] -> (Union@ #[[All, 2]]) & /@ 
                Query[All /* (GatherBy[#, First] &), goColumns]@goData);
            Put[Join[Date[], {identifierAssoc}], fileGOAssociations[[2]]];
            (*save created annotations*)
            NotebookDelete[tempPrint];
            DeleteFile[localFile];
            If[ And @@ (FileExistsQ[#] & /@ fileGOAssociations),
                Print["Created Annotation Files at ", fileGOAssociations],
                Return["Did Not Find Annotation Files, Aborting Process"]
            ],
            (geneOntAssoc = (Get[
                 fileGOAssociations[[
                  1]]])[[-1]];(*otherwise we get from the user specified \
        directory*)
             identifierAssoc = (Get[
                 fileGOAssociations[[2]]])[[-1]];)
        ];(*and again,
        the reverse one*)
        If[ !MatchQ[backGeneSet, 
          All],(*using provided \
     background list to create annotation projection to limited background \
     space*)
            identifierAssoc = 
             Query[backGeneSet /* DeleteMissing]@identifierAssoc;
            geneOntAssoc = 
             GroupBy[#, Last -> First, Union] &@
                Flatten[#, 1] &@(Tuples[{{#[[1]]}, #[[2]]}] & /@ 
                Normal@identifierAssoc)
        (*DeleteCases[{}]@(Query[Union@Flatten@
            Values[identifierAssoc],DeleteCases[Except[x_/;MemberQ[Keys[
            identifierAssoc],x]]]]@geneOntAssoc)*)];
        If[ !MatchQ[lengthFilter, None],
            geneOntAssoc = 
             Query[Select[(lengthFilterFunction[lengthFilter][Length[#]]) &]]@
              geneOntAssoc;
            identifierAssoc = 
             GroupBy[#, Last -> First, Union] &@
                Flatten[#, 1] &@(Tuples[{{#[[1]]}, #[[2]]}] & /@ 
                Normal@geneOntAssoc)
        ];
        returning = 
         Association[
          species -> 
           AssociationThread[{"IDToGO", "GOToID"}, {identifierAssoc, 
             geneOntAssoc}]];
        Return[returning]
    ];

(* ::Function:: *)
(* f:GOAnalysis *)
(***Options***)
Options[GOAnalysis] = {AdditionalFilter -> None (*Select[MatchQ[#[[3,1,2]],"biological_process"]&]*),
    AugmentDictionary -> True,
    BackgroundSet -> All,
    FilterSignificant -> True,
    GeneDictionary -> None,
    GetGeneDictionaryOptions -> {},
    GOAnalysisAssignerOptions -> {},
    HypothesisFunction ->   (Query["Results"]@
       BenjaminiHochbergFDR [#1, SignificanceLevel -> #2] &),
    InputID -> {"UniProt ID", "Gene Symbol"},
    MultipleListCorrection -> None (*Correct for multiple lists, e.g protein+RNA*),
    MultipleList -> False (*whether input is multiple omics or single - for non-omics-object inputs*),
    OBOGODictionaryOptions -> {},
    OBODictionaryVariable -> None,
    OntologyLengthFilter -> 2,
    OutputID -> "UniProt ID",
    pValueCutoff -> 0.05,
    ReportFilter -> 1,
    ReportFilterFunction -> GreaterEqualThan,
    Species -> "human",
    TestFunction -> (N[
       1 - CDF[HypergeometricDistribution[#1, #2, #3], #4 - 1]] &)};
(*Input can be: clustering object*)
(***Function***)
GOAnalysis[dataIn_, OptionsPattern[]] :=
    Module[ {data = dataIn,
    getGeneDictionaryOptions = OptionValue[GetGeneDictionaryOptions],
    augmentDictionary = OptionValue[AugmentDictionary],
    inputID = OptionValue[InputID],
    outputID = OptionValue[OutputID],
    goAnalysisAssignerOptions = OptionValue[GOAnalysisAssignerOptions],
    background = OptionValue[BackgroundSet],
    species = OptionValue[Species],
    lengthFilter = OptionValue[OntologyLengthFilter],
    rptFilter = OptionValue[ReportFilter],
    rptFilterFn = OptionValue[ReportFilterFunction],
    pValCut = OptionValue[pValueCutoff],
    testFn = OptionValue[TestFunction],
    hypothesisFn = OptionValue[HypothesisFunction],
    filterSig = OptionValue[FilterSignificant],
    OBODictVar = OptionValue[OBODictionaryVariable],
    oboGOOptions = OptionValue[OBOGODictionaryOptions],
    multiCorrect = OptionValue[MultipleListCorrection],
    multiList = OptionValue[MultipleList],
    addFilter = OptionValue[AdditionalFilter],
    geneDict = OptionValue[GeneDictionary],
    multiCorr,
    OBODict,
    goAssignment,
    listToggle = False,
    membersWithAssociations,
    testCats,
    countsAll (*counts of all GO categories*),
    totalGenes,
    totalCategories,
    ontologyResultsHCct,
    ontologyResultsFltrd,
    returning},
  (*Obtain OBO dictionary. 
  If externally defined use user definition for OBODict Var*)
        OBODict = 
         If[ MatchQ[OBODictVar, None],
             OBOGODictionary[Sequence @@ oboGOOptions],
             OBODictVar
         ];
        (*Obtain gene dictionary - 
        if it exists can either augment with new information or species or \
      create new, if not exist then create variable*)
        If[ ValueQ[ConstantGeneDictionary],(*variable Exists*)
            If[ augmentDictionary,(*augment*)
                ConstantGeneDictionary = 
                Join[ConstantGeneDictionary, 
                If[ MatchQ[geneDict, None],
                    GetGeneDictionary[Sequence @@ getGeneDictionaryOptions],
                    geneDict
                ]],
    (*replace*)
                ConstantGeneDictionary = 
                If[ MatchQ[geneDict, None],
                    GetGeneDictionary[Sequence @@ getGeneDictionaryOptions],
                    geneDict
                ]
            ],
            (*create/load UCSC based translation dictionary -
            NB global variable or use specified variable*)
            ConstantGeneDictionary =
              If[ MatchQ[geneDict, None],
                  GetGeneDictionary[Sequence @@ getGeneDictionaryOptions],
                  geneDict
              ]
        ];
        (*get the right GO terms for the background requested and correct \
      species*)
        goAssignment = If[ MatchQ[goAnalysisAssignerOptions, {}],
          (*If no specific options for function use background, 
          species request, length request*)
                           GOAnalysisAssigner[BackgroundSet -> background, 
                            Species -> species , LengthFilter -> lengthFilter],
                           GOAnalysisAssigner[Sequence @@ goAnalysisAssignerOptions]
                       ];
        countsAll = Query[species, "GOToID", All, Length]@goAssignment;
        totalGenes = Query[species, "IDToGO", Length]@goAssignment;
        totalCategories = Query[species, "GOToID", Length]@goAssignment;
        (*If the input is a list*)
        If[ MatchQ[Head[data], List],
            listToggle = True;(*LIST entered*)
            data = Association@{"Input List" -> data}
        ];
        If[ MatchQ[Head[data], Association],
         (*check if association of associations*)
            If[ MatchQ[Head@FirstCase[Values[data], Except[Missing]], 
              Association],
             (*clustering object*)
                multiCorr = 
                 Switch[multiCorrect, None, 1, Automatic, 
                  Max@Values@
                    Flatten@Query[All /* Values /* Normal, "GroupAssociations", 
                       All, All /* Tally /* Length, -1]@data, _, multiCorrect];
                (*generate an association of all members with named association*)
                membersWithAssociations =
                Query[All, "GroupAssociations", All,(*Extract the translations, 
                match the ones that intersect for different modes (e.g. RNA/
                Protein label (last element)must be same)*)(Union@#[[All, 
                         
                         1]] ->(*query the GO database in same step*)(Query[# \
                /* Values /* Flatten /* Union /* 
                            DeleteMissing]@(Query[species, "IDToGO"]@
                            goAssignment) &@Union@ Flatten[#[[All, 2]]]) & /@ 
                     Gather[#, ((MatchQ[#1[[1, -1]], #2[[1, -1]]]) && 
                         IntersectingQ[#1[[2]], #2[[2]]]) &] &@({#, (Flatten@
                         Values@GeneTranslation[#[[{1}]], 
                           outputID, ConstantGeneDictionary, InputID -> inputID, 
                           Species -> species])} & /@ #) &) /* Association /* 
                DeleteCases[{Missing[]} | {}]]@data;
            (*testing Category groupings*)
                testCats = (Query[All, 
                    All, {(*Counts in the list*)
                       Sequence@{AssociationThread[#1 -> 
                             ConstantArray[#2, Length[#1]]],
                           (*Counts in the Family*)(Query[species, "GOToID", #1, 
                              Length]@goAssignment), Counts@#1, 
                           Query[#1]@OBODict} &[Flatten@Values[#], Length[#]], 
                       GroupBy[Last -> First]@
                          Flatten[#, 1] &@(Tuples[{{#[[1]]}, #[[2]]}] & /@ 
                          Normal[#])} & /* 
                     Merge[Identity] /* ({(testFn[#[[1]], multiCorr*#[[2]], 
                            multiCorr*totalGenes, #[[3]]]), {#[[1]], 
                           multiCorr*#[[2]], multiCorr*totalGenes, #[[3]]}, #[[
                           4 ;;]]} & /@ # &)]@membersWithAssociations);
                (*Hypothesis Testing*)
                ontologyResultsHCct = 
                 Query[All, 
                   Returner[#, 
                     Applier[hypothesisFn[#, pValCut, totalCategories] &, #, 
                      ListIndex -> 1], ListIndex -> 1] &]@testCats;
                (*Filters*)
                (*Length filter*)
                ontologyResultsFltrd = 
                 Query[All, All, Select[rptFilterFn[rptFilter][#[[2, 4]]] &]]@
                  ontologyResultsHCct;
                returning = 
                 Query[All, 
                   All, (If[ filterSig,
                             Select[#[[1, -1]] &],
                             All
                         ]) /* 
                    SortBy[#[[1, 1]] &]]@ontologyResultsFltrd;
                If[ !MatchQ[addFilter, None],
                    returning = Query[All, All, addFilter]@returning
                ],
                (*Association of Lists*)
                If[ multiList,
                 (*Multi Omics Associations of Lists input*)
                    multiCorr =
                     Switch[multiCorrect, None, 1, Automatic, 
                      Max@Values@Query[All, All /* Tally /* Length, -1]@data, _, 
                      multiCorrect];
                    membersWithAssociations =
                     Query[All,(*Extract the translations, 
                       match the ones that intersect for different modes (e.g. RNA/
                       Protein label (last element)must be same)*)(Union@#[[All, 
                                 1]] ->(*query the GO database in same \
                    step*)(Query[# /* Values /* Flatten /* Union /* 
                                   DeleteMissing]@(Query[species, "IDToGO"]@
                                   goAssignment) &@Union@ Flatten[#[[All, 2]]]) & /@ 
                             Gather[#, ((MatchQ[#1[[1, -1]], #2[[1, -1]]]) && 
                                 IntersectingQ[#1[[2]], #2[[2]]]) &] &@({#, (Flatten@
                                 Values@GeneTranslation[#[[{1}]], 
                                   outputID, ConstantGeneDictionary, InputID -> inputID, 
                                   Species -> species])} & /@ #) &) /* Association /*
                         DeleteCases[{Missing[]} | {}]]@data,
                    (*Single Omics Associations of Lists Input*)
                    multiCorr = 
                     Switch[multiCorrect, None, 1, _, 
                      multiCorrect](*might have put in mixed IDs correspocting to \
               different omics types still, give option to user*);
                    membersWithAssociations =
                     Query[All,(*Extract the translations, 
                       match the ones that intersect for different modes (e.g. RNA/
                       Protein label (last element)must be same)*)(Union@#[[All, 
                                1]] ->(*query the GO database in same step*)(Query[# \
                    /* Values /* Flatten /* Union /* 
                                   DeleteMissing]@(Query[species, "IDToGO"]@
                                   goAssignment) &@Union@ Flatten[#[[All, 2]]]) & /@ 
                            Gather[#, (IntersectingQ[#1[[2]], #2[[2]]]) &] &@({If[ ListQ[#],
                                                                                   #,
                                                                                   #
                                                                               ] , (Flatten@
                                Values@GeneTranslation[If[ ListQ[#],
                                                           #[[{1}]],
                                                           {#}
                                                       ], 
                                  outputID, ConstantGeneDictionary, InputID -> inputID, 
                                  Species -> species])} & /@ #) &)/*Association/*
                       DeleteCases[{Missing[]}|{}]]@data
                ];
                testCats = (Query[
                    All, {(*Counts in the list*)
                       Sequence@{AssociationThread[#1 -> 
                             ConstantArray[#2, Length[#1]]],
                           (*Counts in the Family*)(Query[species, "GOToID", #1, 
                              Length]@goAssignment), Counts@#1, 
                           Query[#1]@OBODict} &[Flatten@Values[#], Length[#]], 
                       GroupBy[Last -> First]@
                          Flatten[#, 1] &@(Tuples[{{#[[1]]}, #[[2]]}] & /@ 
                          Normal[#])} & /* 
                     Merge[Identity] /* ({testFn[#[[1]], multiCorr*#[[2]], 
                           multiCorr*totalGenes, #[[3]]], {#[[1]], 
                           multiCorr*#[[2]], multiCorr*totalGenes, #[[3]]}, #[[
                           4 ;;]]} & /@ # &)]@membersWithAssociations);
                ontologyResultsHCct = 
                 Query[Returner[#, 
                     Applier[hypothesisFn[#, pValCut, totalCategories] &, #, 
                      ListIndex -> 1], ListIndex -> 1] &]@testCats;
                (*Filters*)
                (*Length filter*)
                ontologyResultsFltrd = 
                 Query[All, Select[rptFilterFn[rptFilter][#[[2, 4]]] &]]@
                  ontologyResultsHCct;
                returning = 
                 Query[All, (If[ filterSig,
                                 Select[#[[1, -1]] &],
                                 All
                             ]) /* 
                    SortBy[#[[1, 1]] &]]@ontologyResultsFltrd;
                If[ !MatchQ[addFilter, None],
                    returning = Query[All, addFilter]@returning
                ];
                If[ listToggle,
                    returning = 
                     Values[returning][[1]]
                (*If a single list was provided, 
                    return the association for Gene Ontologies*)]
            ]
        ];
        Return[returning]
    ];
  

(* ::Function:: *)
(* f:GetGeneDictionary *)
(***Options***)
Options[GetGeneDictionary] = {ImportDirectly -> False,
   JavaGBs -> 8,
   KEGGUCSCSplit -> {True,"KEGG Gene ID"} (*Fix nomenclature for UCSC*)
   (*Update these to match any change in query*),
   MathIOmicaDataDirectory -> 
    If[ ! DirectoryQ[
       FileNameJoin[
        Flatten[{FileNameSplit[$UserBaseDirectory], "Applications", 
          "MathIOmica", "MathIOmicaData"}]]],
        CreateDirectory[
         FileNameJoin[
          Flatten[{FileNameSplit[$UserBaseDirectory], "Applications", 
            "MathIOmica", "MathIOmicaData"}]]],
        FileNameJoin[
         Flatten[{FileNameSplit[$UserBaseDirectory], "Applications", 
           "MathIOmica", "MathIOmicaData"}]]
    ],
   Species -> "human",
   UCSCSQLSelectLabels -> 
    Association@{"human" -> {"UCSC ID", "UniProt ID", "Gene Symbol", 
        "RefSeq ID", "NCBI Protein Accession", "Ensembl ID", 
        "KEGG Gene ID", "HGU133Plus2 Affymetrix ID"}}, 
   UCSCSQLString -> 
    Association@{"human" -> 
       "SELECT hg19.kgXref.kgID, hg19.kgXref.spID, \
hg19.kgXref.geneSymbol, hg19.kgXref.refseq, hg19.kgXref.protAcc, \
hg19.knownToEnsembl.value, hg19.knownToKeggEntrez.keggEntrez, \
hg19.knownToU133Plus2.value FROM hg19.kgXref LEFT JOIN \
hg19.knownToEnsembl ON hg19.kgXref.kgID = hg19.knownToEnsembl.name \
LEFT JOIN hg19.knownToKeggEntrez ON hg19.kgXref.kgID = \
hg19.knownToKeggEntrez.name LEFT JOIN hg19.knownToU133Plus2 ON \
hg19.kgXref.kgID = hg19.knownToU133Plus2.name"} (*query for UCSC SQL \
server*)};
(***Function***)
GetGeneDictionary[OptionsPattern[]] :=
    Module[ {geneUCSCTable = 
       FileNameJoin[
        Flatten[{FileNameSplit[OptionValue[MathIOmicaDataDirectory]], 
          OptionValue[Species] <> "GeneUCSCTable"}]],
      importTrueOrFalse = OptionValue[ImportDirectly],
      gigs = ToString[OptionValue[JavaGBs]],
      species = OptionValue[Species],
      ucscSQLString = OptionValue[UCSCSQLString], 
      ucscSQLSelectLabels = OptionValue[UCSCSQLSelectLabels],
      keggSplit = OptionValue[KEGGUCSCSplit],
      ucscDatabase, termTable, tempPrint,
       returning},(*defined local variables*)(*if the user asked us to \
  import directly,import directly with SQL,otherwise,
     get it from a directory they specify*)
        If[ !FileExistsQ[geneUCSCTable],
            tempPrint = 
             PrintTemporary[
              "Did Not Find Gene Translation Files, Attempting to Download from UCSC..."];
            importTrueOrFalse = True
        ];
        If[ importTrueOrFalse,
            ((*run Java through \
            Mathematica*)(*InstallJava[];*)
            
            ReinstallJava[(*CommandLine\[Rule]"java",*)
            JVMArguments -> "-Xmx" <> gigs <> "g"];
            ucscDatabase = 
             OpenSQLConnection[
              JDBC["MySQL(Connector/J)", "genome-mysql.cse.ucsc.edu"], 
              "Username" -> "genomep", 
              "Password" -> "password"];(*get the database from UCSC*)
            If[ MatchQ[ucscDatabase,$Failed],
                Return["Could not establish connection to UCSC. Please try again or add MathIOmica's dictionary manually at "<> geneUCSCTable]
            ];
            termTable = 
            Transpose @
              SQLExecute[ucscDatabase, 
               ucscSQLString[species]] /. {Null | "" -> 
               Missing[]};(*get all the terms we are going to need*)
           \
     (*import with SQL the combined tables,and export with a time stamp*)
            Put[Join[Date[], {termTable}], geneUCSCTable];
            NotebookDelete[tempPrint];
            (*close SQL connection*)
            CloseSQLConnection[ucscDatabase];
            If[ FileExistsQ[geneUCSCTable],
                Print["Created Annotation Files at ", geneUCSCTable],
                Return["Did Not Find Annotation Files, Aborting Process"]
            ]),
            \
            (termTable = (Get[geneUCSCTable])[[-1]];)
        (*otherwise,
     get it from the directory*)];
       (*finally,return what we need*)
        returning = 
         Association[
          species -> 
           AssociationThread[ucscSQLSelectLabels[species] -> termTable]];
        If[ keggSplit[[1]],
            returning[species, keggSplit[[2]]] = 
             If[ MissingQ[#],
                 #,
                 StringSplit[#, {"+"}][[2]]
             ] & /@ 
              returning[species, keggSplit[[2]]]
        ];
        Return[returning]
    ];  

(* ::Function:: *)
(* f:GeneTranslation *)
(***Options***)
Options[GeneTranslation] = { InputID -> None, Species -> "human"};
(***Function***)
GeneTranslation[InputList_, TargetIDList_, GeneDictionary_, 
   OptionsPattern[]] :=
    Module[ {geneDir = GeneDictionary, targetID = TargetIDList, 
      inputID = OptionValue[InputID], species = OptionValue[Species], 
      inputList = InputList, returning},(*LocalVariables*)(*now generate the transformed output set*)
        If[ MatchQ[inputID, None],
            returning = 
             Query[species, targetID, 
               Merge[Query[species, All, Position[#]]@geneDir & /@ inputList, 
                Identity]]@geneDir,
            returning = (Query[species, 
                  targetID, # /* (AssociationThread[
                      
                      inputList -> ((*If[MatchQ[Head[#],List],
                        DeleteMissing[#],#]&@*)(DeleteMissing[
                            Flatten[Union[#]]] /. {} -> 
                            Missing[]) & /@ #)] &)]@
                 geneDir) &@ (Flatten[#] & /@ 
                Query[species /* Values /* Transpose, inputID, 
                  Position[#] & /@ inputList]@
                 geneDir)
        ];(*query of position for given list of identifiers*)
        Return[returning]
    ];

(* ::Function::*)
(* f:KEGGAnalysisAssigner *)
(***Options***) 
Options[KEGGAnalysisAssigner] = {BackgroundSet -> All,
    ImportDirectly -> False,
    KEGGQuery1 -> "pathway",
    KEGGQuery2 -> "hsa",
    Labels -> {"IDToPath", "PathToID"},
    LengthFilter -> None,
    LengthFilterFunction -> GreaterEqualThan,
    MathIOmicaDataDirectory -> 
    If[ ! DirectoryQ[
      FileNameJoin[
       Flatten[{FileNameSplit[$UserBaseDirectory], "Applications", 
         "MathIOmica", "MathIOmicaData"}]]],
        CreateDirectory[
         FileNameJoin[
          Flatten[{FileNameSplit[$UserBaseDirectory], "Applications", 
            "MathIOmica", "MathIOmicaData"}]]],
        FileNameJoin[
         Flatten[{FileNameSplit[$UserBaseDirectory], "Applications", 
           "MathIOmica", "MathIOmicaData"}]]
    ]
  };
(***Function***)
KEGGAnalysisAssigner[OptionsPattern[]] :=
    Module[ {dir = FileNameSplit[OptionValue[MathIOmicaDataDirectory]], 
      importTrueOrFalse = OptionValue[ImportDirectly],
      backMemberSet = OptionValue[BackgroundSet],
      query1 = OptionValue[KEGGQuery1],
      query2 = OptionValue[KEGGQuery2],
      lengthFilter = OptionValue[LengthFilter], 
      lengthFilterFunction = OptionValue[LengthFilterFunction],
      labels = OptionValue[Labels],
      importedIDs, idToPath, pathToID, fileKEGGAssociations, tempPrint, 
      returning},(*local variables defined*)(*if the user asked us to \
   import directly,import directly from KEGG website,otherwise,
     get it from a directory they specify*)
        fileKEGGAssociations = 
         FileNameJoin[Flatten[{dir, #}]] & /@ {query1 <> "_" <> query2 <> 
            "KEGGMemberToPathAssociation", 
           query1 <> "_" <> query2 <> "KEGGPathToMemberAssociation"};
        If[ !(And @@ (FileExistsQ[#] & /@ fileKEGGAssociations)),
            tempPrint = 
             PrintTemporary[
              "Did Not Find Annotation Files, Attempting to Download..."];
            importTrueOrFalse = True
        ];
        If[ importTrueOrFalse,
            If[ FileExistsQ[
              FileNameJoin[Flatten[{dir, query1 <> "_" <> query2 <> ".tsv"}]]],
                DeleteFile[
                FileNameJoin[Flatten[{dir, query1 <> "_" <> query2 <> ".tsv"}]]]
            ];
            URLSave[
             "http://rest.kegg.jp/link/" <> query1 <> 
              If[ MatchQ[query2, ""],
                  "",
                  "/" <> query2
              ], 
             FileNameJoin[Flatten[{dir, query1 <> "_" <> query2 <> ".tsv"}]]
             ];
            importedIDs = 
             Import[FileNameJoin[
               Flatten[{dir, query1 <> "_" <> query2 <> ".tsv"}]]];
            (*gene to pathway association*)
            idToPath = 
             Association@(#[[1, 1]] -> (Union@ #[[All, 2]]) & /@ 
                Query[All /* (GatherBy[#, First] &)]@importedIDs);
            Put[Join[Date[], {idToPath}], fileKEGGAssociations[[1]]];
            (*pathway to gene association*)
            pathToID = 
             Association@(#[[1, 2]] -> (Union@ #[[All, 1]]) & /@ 
                Query[All /* (GatherBy[#, Last] &)]@importedIDs);
            (*time stamp and save associations*)
            Put[Join[Date[], {pathToID}], fileKEGGAssociations[[2]]];
            NotebookDelete[tempPrint];
            DeleteFile[
             FileNameJoin[Flatten[{dir, query1 <> "_" <> query2 <> ".tsv"}]]];
            If[ And @@ (FileExistsQ[#] & /@ fileKEGGAssociations),
                Print["Created Annotation Files at ", fileKEGGAssociations],
                Return["Did Not Find Annotation Files, Aborting Process"]
            ],
            \
            (idToPath = (Get[
            fileKEGGAssociations[[
            1]]])[[-1]];(*otherwise import the necessary associations \
from web*)
             pathToID = (Get[fileKEGGAssociations[[2]]])[[-1]];)
        ];
        If[ !MatchQ[backMemberSet, 
          All],(*using provided \
      background list to create annotation projection to limited background \
      space*)
            idToPath = Query[backMemberSet /* DeleteMissing]@idToPath;
            pathToID = 
             GroupBy[#, Last -> First, Union] &@
                Flatten[#, 1] &@(Tuples[{{#[[1]]}, #[[2]]}] & /@ 
                Normal@idToPath)
        ];
        If[ !MatchQ[lengthFilter, None],
            pathToID = 
             Query[Select[(lengthFilterFunction[lengthFilter][Length[#]]) &]]@
              pathToID;
            idToPath = 
             GroupBy[#, Last -> First, Union] &@
                Flatten[#, 1] &@(Tuples[{{#[[1]]}, #[[2]]}] & /@ 
                Normal@pathToID)
        ];
        returning = 
         Association[
          query2 -> AssociationThread[labels, {idToPath, pathToID}]];
        Return[returning]
    ]

(* ::Function:: *)
(* f:KEGGDictionary *)
(***Options***)
Options[KEGGDictionary] = {ImportDirectly -> False,
    KEGGQuery1 -> "pathway", 
    KEGGQuery2 -> "hsa",
    MathIOmicaDataDirectory -> If[ ! DirectoryQ[
      FileNameJoin[
       Flatten[{FileNameSplit[$UserBaseDirectory], "Applications", 
         "MathIOmica", "MathIOmicaData"}]]],
                                   CreateDirectory[
                                    FileNameJoin[
                                     Flatten[{FileNameSplit[$UserBaseDirectory], "Applications", 
                                       "MathIOmica", "MathIOmicaData"}]]],
                                   FileNameJoin[
                                    Flatten[{FileNameSplit[$UserBaseDirectory], "Applications", 
                                      "MathIOmica", "MathIOmicaData"}]]
                               ]
          };
(***Function***)
KEGGDictionary[OptionsPattern[]] :=
    Module[ {dir = FileNameSplit[OptionValue[MathIOmicaDataDirectory]], 
      importTrueOrFalse = OptionValue[ImportDirectly],
      query1 = OptionValue[KEGGQuery1],
      query2 = OptionValue[KEGGQuery2],
      fileKEGGDict, tempPrint, importedDictionary, 
      associationKEGG},(*local variables defined*)(*if the user asked us \
   to import directly,import directly from KEGG website,otherwise,
     get it from a directory they specify*)
        fileKEGGDict = 
         FileNameJoin[
          Flatten[{dir, {query1 <> "_" <> query2 <> "_KEGGDictionary"}}]];
        If[ (FileExistsQ[fileKEGGDict]),
            associationKEGG = Get[fileKEGGDict][[-1]],
            tempPrint = 
             PrintTemporary[
              "Did Not Find Annotation Files, Attempting to Download..."];
            importTrueOrFalse = True
        ];
        If[ importTrueOrFalse,
            If[ FileExistsQ[
              FileNameJoin[Flatten[{dir, query1 <> "_" <> query2 <> ".tsv"}]]],
                DeleteFile[
                FileNameJoin[Flatten[{dir, query1 <> "_" <> query2 <> ".tsv"}]]]
            ];
            URLSave[
             "http://rest.kegg.jp/list/" <> query1 <> 
              If[ MatchQ[query2, ""],
                  "",
                  "/" <> query2
              ], 
             FileNameJoin[Flatten[{dir, query1 <> "_" <> query2 <> ".tsv"}]]
             ];
            importedDictionary = 
             Import[FileNameJoin[
               Flatten[{dir, query1 <> "_" <> query2 <> ".tsv"}]]];
            associationKEGG = 
             AssociationThread[#[[1]] -> #[[2]] &@
               Transpose@importedDictionary];
            Put[Join[Date[], {associationKEGG}], fileKEGGDict];
            NotebookDelete[tempPrint];
            DeleteFile[
             FileNameJoin[Flatten[{dir, query1 <> "_" <> query2 <> ".tsv"}]]];
            If[ FileExistsQ[fileKEGGDict],
                Print["Created Annotation Files at ", fileKEGGDict],
                Return["Did Not Find Annotation Files, Aborting Process"]
            ]
        ];
        Return[associationKEGG]
    ]

(* ::Function:: *)
(* f:KEGGAnalysis *)
(***Options***)
Options[KEGGAnalysis] = {AdditionalFilter -> None (*Select[MatchQ[#[[3,1,2]],"biological_process"]&]*),
    AnalysisType -> "Genomic" (*options are "Genome", "Molecular","All"*),
    AugmentDictionary -> True,
    BackgroundSet -> All,
    FilterSignificant -> True,
    GeneDictionary -> None,
    GetGeneDictionaryOptions -> {},
    HypothesisFunction ->   (Query["Results"]@
       BenjaminiHochbergFDR [#1, SignificanceLevel -> #2] &),
    InputID -> {"UniProt ID", "Gene Symbol"},
    KEGGAnalysisAssignerOptions -> {},
    KEGGDatabase -> "pathway",
    KEGGDictionaryOptions -> {},
    KEGGDictionaryVariable -> None,
    KEGGMolecular -> "cpd",
    KEGGOrganism -> "hsa",
    MathIOmicaDataDirectory -> 
    If[ ! DirectoryQ[
       FileNameJoin[
        Flatten[{FileNameSplit[$UserBaseDirectory], "Applications", 
          "MathIOmica", "MathIOmicaData"}]]],
        CreateDirectory[
         FileNameJoin[
          Flatten[{FileNameSplit[$UserBaseDirectory], "Applications", 
            "MathIOmica", "MathIOmicaData"}]]],
        FileNameJoin[
         Flatten[{FileNameSplit[$UserBaseDirectory], "Applications", 
           "MathIOmica", "MathIOmicaData"}]]
    ],
    MolecularInputID -> {"cpd"},
    MolecularOutputID -> "cpd",
    MolecularSpecies -> "compound",
    MultipleListCorrection -> None (*Correct for multiple lists, e.g protein+RNA*),
    MultipleList -> 
    False (*whether input is multiple omics or single - for non-omics-object inputs*),
    NonUCSC -> False,
    OutputID -> "KEGG Gene ID",
    PathwayLengthFilter -> 2,
    pValueCutoff -> 0.05,
    ReportFilter -> 1,
    ReportFilterFunction -> GreaterEqualThan,
    Species -> "human" (*Used in GeneDictionary*),
    TestFunction -> (N[
       1 - CDF[HypergeometricDistribution[#1, #2, #3], #4 - 1]] &)
    };
(*Input can be: clustering object*)
(***Function***)

KEGGAnalysis[dataIn_, OptionsPattern[]] :=
    Module[ {data = dataIn,
    type = OptionValue[AnalysisType],
    getGeneDictionaryOptions = OptionValue[GetGeneDictionaryOptions],
    augmentDictionary = OptionValue[AugmentDictionary],
    inputID = OptionValue[InputID],
    outputID = OptionValue[OutputID],
    molInID = OptionValue[MolecularInputID],
    molOutID = OptionValue[MolecularOutputID],
    keggAnalysisAssignOpt = OptionValue[KEGGAnalysisAssignerOptions],
    background = OptionValue[BackgroundSet],
    keggOrg = OptionValue[KEGGOrganism],
    keggMol = OptionValue[KEGGMolecular],
    keggDb = OptionValue[KEGGDatabase],
    lengthFilter = OptionValue[PathwayLengthFilter],
    rptFilter = OptionValue[ReportFilter],
    rptFilterFn = OptionValue[ReportFilterFunction],
    pValCut = OptionValue[pValueCutoff],
    testFn = OptionValue[TestFunction],
    hypothesisFn = OptionValue[HypothesisFunction],
    filterSig = OptionValue[FilterSignificant],
    keggDictVar = OptionValue[KEGGDictionaryVariable],
    keggDictOpt = OptionValue[KEGGDictionaryOptions],
    multiCorrect = OptionValue[MultipleListCorrection],
    multiList = OptionValue[MultipleList],
    addFilter = OptionValue[AdditionalFilter],
    geneDict = OptionValue[GeneDictionary],
    species = OptionValue[Species],
    molSpecies = OptionValue[MolecularSpecies],
    nonUCSC = OptionValue[NonUCSC],
    dirIn = OptionValue[MathIOmicaDataDirectory],
    dir,
    fileMolDict,
    multiCorr,
    keggDict,
    pathAssign,
    listToggle = False,
    membersWithAssociations,
    testCats,
    countsAll (*counts of all KEGG categories*),
    totalMembers,
    totalCategories,
    keggResultsHCct,
    keggResultsFltrd,
    returning},
        dir = FileNameSplit[dirIn];
        Switch[type,
         "Genomic",
         (*Gene Identifier based analysis*)
         (*Obtain OBO dictionary. 
         If externally defined use user definition for OBODict Var*)
         
         keggDict = 
          If[ MatchQ[keggDictVar, None],
              KEGGDictionary[Sequence @@ keggDictOpt],
              keggDictVar
          ];
         (*Obtain gene dictionary - 
         if it exists can either augment with new information or species or \
      create new, if not exist then create variable*)
         If[ ValueQ[ConstantGeneDictionary],(*variable Exists*)
             If[ augmentDictionary,(*augment*)
                 ConstantGeneDictionary = 
                 Join[ConstantGeneDictionary, 
                 If[ MatchQ[geneDict, None],
                     GetGeneDictionary[Sequence @@ getGeneDictionaryOptions],
                     geneDict
                 ]],
     (*replace*)
                 ConstantGeneDictionary = 
                 If[ MatchQ[geneDict, None],
                     GetGeneDictionary[Sequence @@ getGeneDictionaryOptions],
                     geneDict
                 ]
             ],
             (*create/load UCSC based translation dictionary -
             NB global variable or use specified variable*)
             \
ConstantGeneDictionary = 
              If[ MatchQ[geneDict, None],
                  GetGeneDictionary[Sequence @@ getGeneDictionaryOptions],
                  geneDict
              ]
         ];
         (*get the right KEGG terms for the background requested and \
      correct species*)
         pathAssign = If[ MatchQ[keggAnalysisAssignOpt, {}],
           (*If no specific options for function use background, 
           species request, length request*)
                          KEGGAnalysisAssigner[BackgroundSet -> background, 
                           KEGGQuery1 -> keggDb, KEGGQuery2 -> keggOrg , 
                           LengthFilter -> lengthFilter],
                          KEGGAnalysisAssigner[Sequence @@ keggAnalysisAssignOpt]
                      ],
         "Molecular",
         (*molecular analysis*)
         inputID = molInID;
         outputID = molOutID;
         species = molSpecies;
         nonUCSC = True;
         keggOrg = keggMol;
         multiCorrect = None;
         keggDict = 
          If[ MatchQ[keggDictVar, None],
              KEGGDictionary[
               If[ MatchQ[keggDictOpt, {}],
                   Sequence @@ {KEGGQuery1 -> "pathway", KEGGQuery2 -> ""},
                   Sequence @@ keggDictOpt
               ]],
              keggDictVar
          ];
         (*Obtain gene dictionary - 
         if it exists can either augment with new information or species or \
      create new, if not exist then create variable*)
         fileMolDict = 
          FileNameJoin[Flatten[{dir, {"MathIOmicaMolecularDictionary"}}]];
         If[ ValueQ[ConstantGeneDictionary],(*variable Exists*)
             If[ augmentDictionary,(*augment*)
                 ConstantGeneDictionary = 
                 Join[ConstantGeneDictionary, 
                 If[ FileExistsQ[fileMolDict],
                     Get[fileMolDict],
                     Return["Could not find annotation file at " <> fileMolDict <> 
                       " Please either obtain an annotation file from \
mathiomica.org or provide a GeneDictionary option variable."]
                 ]],
     (*replace*)
                 ConstantGeneDictionary = 
                 If[ MatchQ[geneDict, None],
                     If[ FileExistsQ[fileMolDict],
                         Get[fileMolDict],
                         Return["Could not find annotation file at " <> fileMolDict <> 
                           " Please either obtain an annotation file from \
mathiomica.org or provide a GeneDictionary option variable."]
                     ],
                     geneDict
                 ]
             ],
             (*create/load UCSC based translation dictionary -
             NB global variable or use specified variable*)
             \
ConstantGeneDictionary = 
              If[ MatchQ[geneDict, None],
                  If[ FileExistsQ[fileMolDict],
                      Get[fileMolDict],
                      Return["Could not find annotation file at " <> fileMolDict <> 
                        " Please either obtain an annotation file from \
mathiomica.org or provide a GeneDictionary option variable."]
                  ],
                  geneDict
              ]
         ];
         (*get the right KEGG terms for the background requested and \
      correct species*)
         pathAssign = If[ MatchQ[keggAnalysisAssignOpt, {}],
           (*If no specific options for function use background, 
           species request, length request*)
                          KEGGAnalysisAssigner[BackgroundSet -> background, 
                           KEGGQuery1 -> keggDb, KEGGQuery2 -> keggOrg , 
                           LengthFilter -> lengthFilter],
                          KEGGAnalysisAssigner[Sequence @@ keggAnalysisAssignOpt]
                      ],
         All,
         returning = 
          Association @@ (# ->  
               KEGGAnalysis[data, 
                GetGeneDictionaryOptions -> getGeneDictionaryOptions,
                AugmentDictionary -> augmentDictionary,
                InputID -> inputID,
                OutputID -> outputID,
                MolecularInputID -> molInID,
                MolecularOutputID -> molOutID,
                KEGGAnalysisAssignerOptions -> keggAnalysisAssignOpt,
                BackgroundSet -> background,
                KEGGOrganism -> keggOrg,
                KEGGMolecular -> keggMol,
                KEGGDatabase -> keggDb,
                PathwayLengthFilter -> lengthFilter,
                ReportFilter -> rptFilter,
                ReportFilterFunction -> rptFilterFn,
                pValueCutoff -> pValCut,
                TestFunction -> testFn,
                HypothesisFunction ->   hypothesisFn,
                FilterSignificant -> filterSig,
                KEGGDictionaryVariable -> keggDictVar,
                KEGGDictionaryOptions -> keggDictOpt,
                
                MultipleListCorrection -> 
                 multiCorrect (*Correct for multiple lists, e.g protein+
                RNA*),
                
                MultipleList -> 
                 multiList (*whether input is multiple omics or single - 
                for non-omics-object inputs*),
                AdditionalFilter -> addFilter (*Select[MatchQ[#[[3,1,2]],
                "biological_process"]&]*),
                GeneDictionary -> geneDict,
                Species -> species (*Used in GeneDictionary*),
                MolecularSpecies -> molSpecies,
                NonUCSC -> nonUCSC,
                AnalysisType -> # (*options are "Genome", "Molecular",
                "All"*),
                MathIOmicaDataDirectory -> dirIn] & /@ {"Molecular", 
              "Genomic"});
         Return[returning],
         __,
         (*non match*)
         
         Return["AnalysisType\[Rule]" <> 
           If[ StringQ[type],
               "\"" <> type <> "\"",
               ToString@type
           ] <> 
           " is not a valid  choice."]];
        countsAll = Query[keggOrg, "PathToID", All, Length]@pathAssign;
        totalMembers = Query[keggOrg, "IDToPath", Length]@pathAssign;
        totalCategories = Query[keggOrg, "PathToID", Length]@pathAssign;
        (*If the input is a list*)
        If[ MatchQ[Head[data], List],
            listToggle = True;(*LIST entered*)
            data = Association@{"Input List" -> data}
        ];
        If[ MatchQ[Head[data], Association],
         (*check if association of associations*)
            If[ MatchQ[Head@FirstCase[Values[data], Except[Missing]], 
              Association],
             (*clustering object*)
                multiCorr = 
                 Switch[multiCorrect, None, 1, Automatic, 
                  Max@Values@
                    Flatten@Query[All /* Values /* Normal, "GroupAssociations", 
                       All, All /* Tally /* Length, -1]@data, _, multiCorrect];
                (*generate an association of all members with named association*)
                membersWithAssociations =
                Query[All, "GroupAssociations", All,(*Extract the translations, 
                match the ones that intersect for different modes (e.g. RNA/
                Protein label (last element)must be same)*)(Union@#[[All, 
                          1]] ->(*query the KEGG database in same \
                step*)(Query[# /* Values /* Flatten /* Union /* 
                            DeleteMissing]@(Query[keggOrg, "IDToPath"]@
                            pathAssign) &@Union@ Flatten[#[[All, 2]]]) & /@ 
                     Gather[#, ((MatchQ[#1[[1, -1]], #2[[1, -1]]]) && 
                         IntersectingQ[#1[[2]], #2[[2]]]) &] &@({#, 
                       If[ nonUCSC,
                           #,
                           If[ MissingQ[#],
                               #,
                               keggOrg <> ":" <> #
                           ]
                       ] & /@ (Flatten@
                          Values@GeneTranslation[#[[{1}]], 
                            outputID, ConstantGeneDictionary, InputID -> inputID, 
                            Species -> species])} & /@ #) &) /* Association /*
                 DeleteCases[{Missing[]} | {}]]@data;
            (*testing Category groupings*)
                testCats = (Query[All, 
                    All, {(*Counts in the list*)
                       Sequence@{AssociationThread[#1 -> 
                             ConstantArray[#2, Length[#1]]],
                           (*Counts in the Family*)(Query[keggOrg, "PathToID", #1,
                               Length]@pathAssign), Counts@#1, 
                           Query[#1]@keggDict} &[Flatten@Values[#], Length[#]], 
                       GroupBy[Last -> First]@
                          Flatten[#, 1] &@(Tuples[{{#[[1]]}, #[[2]]}] & /@ 
                          Normal[#])} & /* 
                     Merge[Identity] /* ({(testFn[#[[1]], multiCorr*#[[2]], 
                            multiCorr*totalMembers, #[[3]]]), {#[[1]], 
                           multiCorr*#[[2]], 
                           multiCorr*totalMembers, #[[3]]}, #[[4 ;;]]} & /@ # &)]@
                   membersWithAssociations);
                (*Hypothesis Testing*)
                keggResultsHCct = 
                 Query[All, 
                   Returner[#, 
                     Applier[hypothesisFn[#, pValCut, totalCategories] &, #, 
                      ListIndex -> 1], ListIndex -> 1] &]@testCats;
                (*Filters*)
                (*Length filter*)
                keggResultsFltrd = 
                 Query[All, All, 
                   Select[(*If[MatchQ[Values[#],{}],#,*)
                    rptFilterFn[rptFilter][#[[2, 4]]](*]*)&]]@keggResultsHCct;
                returning = 
                 Query[All, 
                   All, (If[ filterSig,
                             Select[#[[1, -1]] &],
                             All
                         ]) /* 
                    SortBy[#[[1, 1]] &]]@keggResultsFltrd;
                If[ !MatchQ[addFilter, None],
                    returning = Query[All, All, addFilter]@returning
                ],
                
                (*Association of Lists*)
                If[ multiList,
                 (*Multi Omics Associations of Lists input*)
                    multiCorr =
                     Switch[multiCorrect, None, 1, Automatic, 
                      Max@Values@Query[All, All /* Tally /* Length, -1]@data, _, 
                      multiCorrect];
                    membersWithAssociations =
                     Query[All,(*Extract the translations, 
                       match the ones that intersect for different modes (e.g. RNA/
                       Protein label (last element)must be same)*)(Union@#[[All, 
                                  1]] ->(*query the KEGG database in same \
                    step*)(Query[# /* Values /* Flatten /* Union /* 
                                   DeleteMissing]@(Query[keggOrg, "IDToPath"]@
                                   pathAssign) &@Union@ Flatten[#[[All, 2]]]) & /@ 
                             Gather[#, ((MatchQ[#1[[1, -1]], #2[[1, -1]]]) && 
                                 IntersectingQ[#1[[2]], #2[[2]]]) &] &@({#, 
                          
                                If[ nonUCSC,
                                    #,
                                    If[ MissingQ[#],
                                        #,
                                        keggOrg <> ":" <> #
                                    ]
                                ] & /@ (Flatten@
                                  Values@
                                   GeneTranslation[#[[{1}]], 
                                   outputID, ConstantGeneDictionary, InputID -> inputID, 
                                   Species -> species])} & /@ #) &) /* Association /*
                         DeleteCases[{Missing[]} | {}]]@data,
                    (*Single Omics Associations of Lists Input*)
                    multiCorr = 
                     Switch[multiCorrect, None, 1, _, 
                      multiCorrect](*might have put in mixed IDs correspocting to \
               different omics types still, give option to user*);
                    membersWithAssociations =
                     Query[All,(*Extract the translations, 
                       match the ones that intersect for different modes (e.g. RNA/
                       Protein label (last element)must be same)*)(Union@#[[All, 
                                 1]] ->(*query the KEGG database in same \
                    step*)(Query[# /* Values /* Flatten /* Union /* 
                                   DeleteMissing]@(Query[keggOrg, "IDToPath"]@
                                   pathAssign) &@Union@ Flatten[#[[All, 2]]]) & /@ 
                            Gather[#, (IntersectingQ[#1[[2]], #2[[2]]]) &] &@({If[ ListQ[#],
                                                                                   #,
                                                                                   #
                                                                               ] , 
                              
                              If[ nonUCSC,
                                  #,
                                  If[ MissingQ[#],
                                      #,
                                      keggOrg <> ":" <> #
                                  ]
                              ] & /@ (Flatten@
                                 Values@GeneTranslation[If[ ListQ[#],
                                                            #[[{1}]],
                                                            {#}
                                                        ], 
                                   outputID, ConstantGeneDictionary, InputID -> inputID, 
                                   Species -> species])} & /@ #) &)/*Association/*
                       DeleteCases[{Missing[]}|{}]]@data
                ];
                testCats = (Query[
                    All, {(*Counts in the list*)
                       Sequence@{AssociationThread[#1 -> 
                             ConstantArray[#2, Length[#1]]],
                           (*Counts in the Family*)(Query[keggOrg, "PathToID", #1,
                               Length]@pathAssign), Counts@#1, 
                           Query[#1]@keggDict} &[Flatten@Values[#], Length[#]], 
                       GroupBy[Last -> First]@
                          Flatten[#, 1] &@(Tuples[{{#[[1]]}, #[[2]]}] & /@ 
                          Normal[#])} & /* 
                     Merge[Identity] /* ({testFn[#[[1]], multiCorr*#[[2]], 
                           multiCorr*totalMembers, #[[3]]], {#[[1]], 
                           multiCorr*#[[2]], 
                           multiCorr*totalMembers, #[[3]]}, #[[4 ;;]]} & /@ # &)]@
                   membersWithAssociations);
                keggResultsHCct = 
                 Query[Returner[#, 
                     Applier[hypothesisFn[#, pValCut, totalCategories] &, #, 
                      ListIndex -> 1], ListIndex -> 1] &]@testCats;
                (*Filters*)
                (*Length filter*)
                keggResultsFltrd = 
                 Query[All, Select[rptFilterFn[rptFilter][#[[2, 4]]] &]]@
                  keggResultsHCct;
                returning = 
                 Query[All, (If[ filterSig,
                                 Select[#[[1, -1]] &],
                                 All
                             ]) /* 
                    SortBy[#[[1, 1]] &]]@keggResultsFltrd;
                If[ !MatchQ[addFilter, None],
                    returning = Query[All, addFilter]@returning
                ];
                If[ listToggle,
                    returning = 
                     Values[returning][[1]]
                (*If a single list was provided, 
                    return the association for Gene Ontologies*)]
            ]
        ];
        Return[returning]
    ];
  
(* ::Function:: *)
(* f:MassMatcher *)
(***Options***)
Options[MassMatcher] = {MassDictionaryVariable -> None,MolecularSpecies->"cpd"}
(***Function***)
MassMatcher[data_, accuracy_, OptionsPattern[]] :=
    Module[ {x = data,
        ppm = accuracy*(10^-6), 
        dictIn = OptionValue[MassDictionaryVariable],
        molSpecies = OptionValue[MolecularSpecies],
        dict,
        returning},
        dict = If[ MatchQ[dictIn,None],
                   MassDictionary[],
                   dictIn
               ];
        returning = Keys@Query[Select[x (1 - ppm) < # < x (1 + ppm) &]]@dict[molSpecies];
        Return[returning]
    ];
    
(* ::Function:: *)
(* f:MassDictionary:: *)
(*** Options ***)
Options[MassDictionary] = {MathIOmicaDataDirectory -> 
    If[ ! DirectoryQ[
       FileNameJoin[
        Flatten[{FileNameSplit[$UserBaseDirectory], "Applications", 
          "MathIOmica", "MathIOmicaData"}]]],
        CreateDirectory[
         FileNameJoin[
          Flatten[{FileNameSplit[$UserBaseDirectory], "Applications", 
            "MathIOmica", "MathIOmicaData"}]]],
        FileNameJoin[
         Flatten[{FileNameSplit[$UserBaseDirectory], "Applications", 
           "MathIOmica", "MathIOmicaData"}]]
    ]};
MassDictionary[OptionsPattern[]] :=
    Module[ {dir = OptionValue[MathIOmicaDataDirectory],fileMassDict},
        fileMassDict = 
              FileNameJoin[Flatten[{dir, {"MathIOmicaMassDictionary"}}]];
        If[ FileExistsQ[fileMassDict],
            Get[fileMassDict],
            Return["Could not find MathIOmica's mass dictionary at " <> fileMassDict <> 
              " Please either obtain a mass dictionary file from mathiomica.org or provide a custom file at the above location."]
        ]
    ];

(* ::Function:: *)
(* f:OmicsObjectUniqueMassConverter:: *)
(*** Options ***)
Options[OmicsObjectUniqueMassConverter] := {MassMatcherOptions -> {}};
(***Function***)
OmicsObjectUniqueMassConverter[omicsObject_, massAccuracy_, 
   OptionsPattern[]] := 
  Module[{data = omicsObject, ppm = massAccuracy, 
    massMatcherOpts = OptionValue[MassMatcherOptions], keyMapper, 
    returning},
   keyMapper = 
    Association@(( #[[2]] -> 
          If[MatchQ[#[[1]], {}], #[[2]], 
           If[Length[#[[1]]] == 1, 
            Flatten[{#[[1]], #[[2]]}], #[[2]]]]) & /@ \
(({MassMatcher[#[[1]], ppm, 
             Sequence @@ massMatcherOpts], #} &) /@ 
         Query[1, Keys]@data));
   returning = Query[All, KeyMap[keyMapper]]@data;
   Return[returning]];


(* ::Function:: *)
(* f:EnrichmentReportExport *)
(***Options***)
Options[EnrichmentReportExport] = {AppendString-> "",OutputDirectory -> None};
(***Function***)
EnrichmentReportExport[results_, OptionsPattern[]] :=
    Module[ {res = results, dir = OptionValue[OutputDirectory],appendString=OptionValue[AppendString]},
        If[ (! MemberQ[Keys@NotebookInformation[], "FileName"]) && (MatchQ[
            dir, None]),
            Print["Please save current notebook and try again, or select a \
directory by setting the OutputDirectory option."];
            NotebookSave[],
            If[MatchQ[appendString,""],appendString=(DateString[] // StringReplace[{" " -> "_", ":" -> "_"}])];
            Export[FileNameJoin[{If[ MatchQ[dir, None],
                                     NotebookDirectory[],
                                     dir
                                 ], #[[1]] <> "_"<>appendString<>".xlsx"}], "Sheets" -> #[[2]], "Rules"] & /@ 
             Query[Transpose@{Keys[#], Values[#]} &, Normal, 
               Flatten[#, 2] & /@ Transpose@{Keys[#], Values[#]} &]@(res)
        ]
    ];

(* ::Section:: *)
(*#####CLassificationAndClustering#####*)

(* ::Function:: *)
(* f:MissingDataCreator*)
(*Fill In of Missing Data*)
(***Function***)
MissingDataCreator[data_, setSamples_] :=
    Module[ {dataIn = data, pSamples = setSamples, 
      inputAssociation, returning},
     (*local variables created*)
        inputAssociation = 
         AssociationThread[Sequence @@ Transpose[dataIn]];
        returning = (Query[(Key[#] & /@ pSamples) /* Values]@
            inputAssociation) /. _Missing -> 
           Missing[];(*substitution to ensure forward compatibility + 
        other functions*)
        Return[returning]
    ];
  
(* ::Function:: *)
(* f:TimeSeriesClassification*)
(*Identifies classes of temporal behavior in time series datasets.*)
(***Options***)
Options[TimeSeriesClassification] = {Method -> 
    "LombScargle"(*Which method,
   currently other option"Autocorrelation"*), 
   LombScargleOptions -> {PairReturn -> 
      False,NormalizeIntensities->True}  (*Method Specific Options to be Passed to function*), 
   AutocorrelationOptions -> {UpperFrequencyFactor -> 
      1},(*Autocorrelation Options*)
   AutocorrelationLogic -> False,(*Common Options*)
   ReturnAllSpikes -> False, ReturnData -> True, 
   SpikeCutoffs -> <|1->{0.99,-0.99},2->{0.99,-0.99}|>,
   ReturnModels -> False,
   AutocorrelationCutoffs -> {0},
   LombScargleCutoff -> 0,
   InterpolationOptions -> {},
   InterpolationDeltaT -> "Auto"};
(***Function***)
TimeSeriesClassification[data_, setTimes_, OptionsPattern[]] :=
    Module[ {dataInPreClean = N[data], dataIn, 
      associationToggle = AssociationQ[data], 
      inputSetTimes = N[setTimes], 
      autoC = OptionValue[AutocorrelationCutoffs], 
      lsC = OptionValue[LombScargleCutoff], 
      spikeC = OptionValue[SpikeCutoffs], 
      spikesAll = OptionValue[ReturnAllSpikes], 
      uniqueAutocorr = OptionValue[AutocorrelationLogic], 
      dataReturn = OptionValue[ReturnData], 
      method = OptionValue[Method],
      (*InterpolatedAutocorrelation variables*)
      iDeltaT = OptionValue[InterpolationDeltaT],
      interOptions = OptionValue[InterpolationOptions], iTimePts, 
      inputTLength,
      (*autocorrelation variables*)
      autoCorrOpts = OptionValue[AutocorrelationOptions],
      returnModels = OptionValue[ReturnModels],
      autocorrelationList, selectedAutocorrs,(*periodogram variables*)
      lsOpts = OptionValue[LombScargleOptions], lsAll, lsClass, 
      spikeData, spikeMax, spikeMin,(*spike variables*)spikeMaxTest, 
      spikeMinTest, spikesUniq, allGroups, classes, groups, models, 
      counter, 
      returning},(*using local variables*)(*check if the input is an \
  association,if not create one*)
        If[ ! associationToggle,
            dataInPreClean = 
             AssociationThread[Range[Length[#]], #] &@dataInPreClean
        ];
        (*Clean the time series*)
        dataIn = ConstantSeriesClean[dataInPreClean];
        (*If pairs of {times,
        values} given create missing data format based on input full set \
     times*)
        If[ Length[dataIn[[1, 1]]] != 0,
            dataIn = 
             Query[All, MissingDataCreator[#,inputSetTimes] &]@dataIn
        ];
    (*select procedure based on method*)
        Switch[method,(****)"Autocorrelation"(****), Print["Method \[Rule] \"Autocorrelation\" "];
                                                     autocorrelationList = 
                                                      Query[All, (Autocorrelation[#, inputSetTimes, 
                                                           Sequence @@ autoCorrOpts] &)]@dataIn;
                                                     (*calculate the autocorrelation list for each \
                                                 timeseries*)(*select Unique autocorrelations based on provided \
                                                 significance cutoffs*)
                                                     If[ uniqueAutocorr,
                                                         selectedAutocorrs = 
                                                         Query[All, 
                                                         MapThread[
                                                         GreaterEqual[#1, #2] &, {#[[2 ;; 1 + Length[autoC]]], 
                                                         autoC}] &]@autocorrelationList;
                                                         Return[selectedAutocorrs],
                                                         selectedAutocorrs = 
                                                          Query[All, 
                                                            MapThread[
                                                               GreaterEqual[#1, #2] &, {#[[2 ;; 1 + Length[autoC]]], 
                                                                autoC}] & /* (FirstPosition[#, True] &)]@
                                                           autocorrelationList
                                                     ];
                                                     (*if no unique autocorrelations are significant,
                                                     test the remainder data for spike behavior:
                                                     N.B.spike max given priority over spike min unless spikesAll \
                                                 option is toggled*)
                                                     If[ spikesAll,(*check for spike Max*)
                                                         spikeMaxTest = 
                                                          DeleteCases[#, False] &@
                                                           Query[Key[#] & /@ Keys@Select[#, MissingQ] &@
                                                              selectedAutocorrs, ((Max[#] > (spikeC[Length[#]][[1]])) &@
                                                                 Cases[#, _Real] &) /* (# /. {True -> "SpikeMax"} &)]@
                                                            dataIn;
                                                         (*check for spike Min*)
                                                         spikeMinTest = 
                                                          DeleteCases[#, False] &@
                                                           Query[Key[#] & /@ Keys@Select[#, MissingQ] &@
                                                              selectedAutocorrs, ((Min[#] < (spikeC[Length[#]][[2]])) &@
                                                                 Cases[#, _Real] &) /* (# /. {True -> "SpikeMin"} &)]@
                                                            dataIn;
                                                         allGroups = (Query[All, 
                                                             Keys]@(Join @@ {GroupBy[spikeMaxTest, # &], 
                                                               GroupBy[spikeMinTest, # &], 
                                                               GroupBy[DeleteMissing@selectedAutocorrs, # &]})),
                                                         spikesUniq = 
                                                          DeleteMissing@
                                                           Query[Key[#] & /@ Keys@Select[#, MissingQ] &@
                                                              selectedAutocorrs, ({Max[#] > (spikeC[Length[#]][[1]]), 
                                                                   Min[#] < (spikeC[Length[#]][[2]])} &@
                                                                 Cases[#, _Real] &) /* (FirstPosition[#, 
                                                                 True] &) /* (# /. {{1} -> "SpikeMax", {2} -> 
                                                                   "SpikeMin"} &)]@dataIn;
                                                         allGroups = (Query[All, 
                                                             Keys]@(Join @@ {GroupBy[spikesUniq, # &], 
                                                               GroupBy[DeleteMissing@selectedAutocorrs, # &]}))
                                                     ];
                                                     returning = 
                                                      If[ dataReturn,
                                                          If[ returnModels,
                                                              AssociationThread[{"TimeSeriesClasses", 
                                                                "Autocorrelations"}, {KeyMap[
                                                                  If[ MemberQ[{"SpikeMax", "SpikeMin"}, #],
                                                                      #,
                                                                      "Lag" <> ToString[First@#]
                                                                  ] &]@
                                                                 KeySort@Query[All, All /* (Association[#] &), 
                                                                    Association@(# -> {Query[Key[#]]@autocorrelationList, 
                                                                         Query[Key[#]]@dataIn} &)]@allGroups, 
                                                                autocorrelationList}],
                                                              KeyMap[If[ MemberQ[{"SpikeMax", "SpikeMin"}, #],
                                                                         #,
                                                                         "Lag" <> ToString[First@#]
                                                                     ] &]@
                                                               KeySort@
                                                                Query[All, All /* (Association[#] &), 
                                                                  Association@(# -> {Query[Key[#]]@autocorrelationList, 
                                                                       Query[Key[#]]@dataIn} &)]@allGroups
                                                          ],
                                                          KeyMap[If[ MemberQ[{"SpikeMax", "SpikeMin"}, #],
                                                                     #,
                                                                     "Lag" <> ToString[First@#]
                                                                 ] &]@KeySort@allGroups
                                                      ];,(****)
         "LombScargle"(****), Print["Method \[Rule] \"LombScargle\" "];
                              lsAll = ({LombScargle[#, inputSetTimes, 
                                     Sequence @@ lsOpts][[2]], #} & /@ dataIn);
                              lsClass = 
                               KeySort@GroupBy[
                                  Last -> First]@(Query[All /* Select[Max[#[[1]]] > lsC &], 
                                    All /* ({#, 
                                        "f" <> ToString@(Sequence @@ Ordering[#[[1]], -1])} &)]@
                                   lsAll);
                              spikeData = 
                               KeyDrop[lsAll, 
                                Query[All /* Values /* Keys /* (Key[#] & /@ Flatten[#, 1] &)]@
                                 lsClass];
                              spikeMax = 
                               Association@("SpikeMax" -> 
                                  Query[All /* DeleteMissing /* 
                                     Select[(Max[#] > (spikeC[Length[#]][[1]])) &@
                                        Cases[#[[2]], _Real] &]]@spikeData);
                              If[ spikesAll,
                                  spikeMin = 
                                   Association@("SpikeMin" -> 
                                      Query[All /* DeleteMissing /* 
                                         Select[(Min[#] < (spikeC[Length[#]][[2]])) &@
                                            Cases[#[[2]], _Real] &]]@spikeData);
                                  Print[spikeMin];(*return all spikes with multi-class belonging*),
                                  spikeMin = 
                                  Association@("SpikeMin" -> 
                                     Query[All /* DeleteMissing /* 
                                        Select[(Min[#] < (spikeC[Length[#]][[2]])) &@
                                           Cases[#[[2]], _Real] &]]@(KeyDrop[
                                         Key[#] & /@ Flatten[#, 1] &@Keys@Values@spikeMax]@
                                        spikeData))
                              ];
                              allGroups = Join[spikeMax, spikeMin, lsClass];
                              (*return unique spikes,priority for max*)
                              returning = 
                               If[ dataReturn,
                                   If[ returnModels,
                                       AssociationThread[{"TimeSeriesClasses", 
                                         "LombScargle"}, {DeleteCases[<||>]@allGroups, lsAll}],
                                       allGroups
                                   ],
                                   Query[All, Keys]@allGroups
                               ];,
         "TimeSeriesModelAggregate", Print["Method \[Rule] \"TimeSeriesModelAggregate\" "];
                                     classes = {(If[ !MatchQ[Head[#],TimeSeriesModelFit],
                                                     {#[{"AIC", "BestFit", "BestFitParameters", 
                                                     "CandidateSelectionTable"}], #},
                                                     {{"Invalid",Missing[],{},{{Null, {Null, {"Missing"}}}}},#}
                                                 ] &@(Quiet@TimeSeriesModelFit@
                                             DeleteCases[{x_, y_} /; (MissingQ[x] || MissingQ[y])]@
                                              Transpose[{inputSetTimes, #}])), #} & /@ dataIn;
                                     groups = 
                                      GroupBy[First -> Last]@
                                       Query[All, {ToString@Head@#[[1, 1, 2]], {Chop@#[[2]], 
                                            Chop@#[[2]]}} &]@classes;
                                     models = Query[All, 1, 2]@classes;
                                     returning = 
                                     If[ dataReturn,
                                         If[ returnModels,
                                             AssociationThread[{"TimeSeriesClasses", "Models"}, {groups, 
                                               models}],
                                             groups
                                         ],
                                         Query[All,Keys]@groups
                                     ];,
         "TimeSeriesModelDetailed",
         (*Time Series Models with Degree subdivision*)
         Print["Method \[Rule] \"TimeSeriesModelDetailed\" "];
         classes = {(If[ !MatchQ[Head[#],TimeSeriesModelFit],
                         {#[{"AIC", "BestFit", "BestFitParameters", 
                         "CandidateSelectionTable"}], #},
                         {{"Invalid",Missing[],{},{{Null, {Null, {"Missing"}}}}},#}
                     ] &@(Quiet@TimeSeriesModelFit@
         DeleteCases[{x_, y_} /; (MissingQ[x] || MissingQ[y])]@
         Transpose[{inputSetTimes, #}])), #} & /@ dataIn;
         models = Query[All, 1, 2]@classes;
         counter = 1;
         groups = KeySort@(GroupBy[First -> Last]@
             Query[All, {(Normal@#[[1, 1, 4]])[[1, 2, 2, 
                   1]], {#[[1, 1, 3]] /. {{} -> {counter++}}, Chop@#[[2]]}} &]@
              classes);
         returning = 
          If[ dataReturn,
              If[ returnModels,
                  AssociationThread[{"TimeSeriesClasses", "Models"}, {groups, 
                    models}],
                  groups
              ],
              Query[All,Keys]@groups
          ];,
         "InterpolatedAutocorrelation",
         Print["Method \[Rule] \"InterpolatedAutocorrelation\" "];
         iTimePts = 
          Table[i, {i, Min[#], Max[#], 
              If[ MatchQ[iDeltaT, "Auto"],
                  (Max[#] - Min[#])/(Length[#] - 1),
                  iDeltaT
              ]}] &@inputSetTimes;
         inputTLength = Length@iTimePts;
         autocorrelationList = 
          CorrelationFunction[#, {0, 
                inputTLength - 
            1}] &@((Interpolation[#, Sequence @@ interOptions] &@
                 N@DeleteCases[{x_, y_} /; (MissingQ[x] || MissingQ[y])]@
                   Transpose[{inputSetTimes, #}])[iTimePts]) & /@ dataIn;
         (*calculate the autocorrelation list for each \
     timeseries*)(*select Unique autocorrelations based on provided \
     significance cutoffs*)
         If[ uniqueAutocorr,
             selectedAutocorrs = 
             Query[All, 
             MapThread[
             GreaterEqual[#1, #2] &, {#[[2 ;; 1 + Length[autoC]]], 
             autoC}] &]@autocorrelationList;
             Return[selectedAutocorrs],
             selectedAutocorrs = 
              Query[All, 
                MapThread[
                   GreaterEqual[#1, #2] &, {#[[2 ;; 1 + Length[autoC]]], 
                    autoC}] & /* (FirstPosition[#, True] &)]@
               autocorrelationList
         ];
         (*if no unique autocorrelations are significant,
         test the remainder data for spike behavior:
         N.B.spike max given priority over spike min unless spikesAll \
     option is toggled*)
         If[ spikesAll,(*check for spike Max*)
             spikeMaxTest = 
              DeleteCases[#, False] &@
               Query[Key[#] & /@ Keys@Select[#, MissingQ] &@
                  selectedAutocorrs, ((Max[#] > (spikeC[Length[#]][[1]])) &@
                     Cases[#, _Real] &) /* (# /. {True -> "SpikeMax"} &)]@
                dataIn;
             (*check for spike Min*)
             spikeMinTest = 
              DeleteCases[#, False] &@
               Query[Key[#] & /@ Keys@Select[#, MissingQ] &@
                  selectedAutocorrs, ((Min[#] < (spikeC[Length[#]][[2]])) &@
                     Cases[#, _Real] &) /* (# /. {True -> "SpikeMin"} &)]@
                dataIn;
             allGroups = (Query[All, 
                 Keys]@(Join @@ {GroupBy[spikeMaxTest, # &], 
                   GroupBy[spikeMinTest, # &], 
                   GroupBy[DeleteMissing@selectedAutocorrs, # &]})),
             spikesUniq = 
              DeleteMissing@
               Query[Key[#] & /@ Keys@Select[#, MissingQ] &@
                  selectedAutocorrs, ({Max[#] > (spikeC[Length[#]][[1]]), 
                       Min[#] < (spikeC[Length[#]][[2]])} &@
                     Cases[#, _Real] &) /* (FirstPosition[#, 
                     True] &) /* (# /. {{1} -> "SpikeMax", {2} -> 
                       "SpikeMin"} &)]@dataIn;
             allGroups = (Query[All, 
                 Keys]@(Join @@ {GroupBy[spikesUniq, # &], 
                   GroupBy[DeleteMissing@selectedAutocorrs, # &]}))
         ];
         returning = 
          If[ dataReturn,
              If[ returnModels,
                  AssociationThread[{"TimeSeriesClasses", 
                    "Autocorrelations"}, {KeyMap[
                      If[ MemberQ[{"SpikeMax", "SpikeMin"}, #],
                          #,
                          "Lag" <> ToString[First@#]
                      ] &]@
                     KeySort@Query[All, All /* (Association[#] &), 
                        Association@(# -> {Query[Key[#]]@autocorrelationList, 
                             Query[Key[#]]@dataIn} &)]@allGroups, 
                    autocorrelationList}],
                  KeyMap[If[ MemberQ[{"SpikeMax", "SpikeMin"}, #],
                             #,
                             "Lag" <> ToString[First@#]
                         ] &]@
                   KeySort@Query[All, All /* (Association[#] &), 
                      Association@(# -> {Query[Key[#]]@autocorrelationList, 
                           Query[Key[#]]@dataIn} &)]@allGroups
              ],
              KeyMap[If[ MemberQ[{"SpikeMax", "SpikeMin"}, #],
                         #,
                         "Lag" <> ToString[First@#]
                     ] &]@KeySort@allGroups
          ];,
         (****)_(****), 
         Return["Not a valid Method, Please choose Method\[Rule]option by \
selecting from one of these options: {Autocorrelation,LombScargle}, \
Default is Method\[Rule]LombScargle"]];
        Return[DeleteCases[<||>]@returning]
    ];

(* ::Function:: *)
(* f:TimeSeriesClusters *)
(*Hierarchical Clustering of a Time Series.*)
(***Options***)
Options[TimeSeriesClusters] = {
    ClusterLabeling -> "",
    DendrogramPlotOptions -> {},
    DistanceFunction -> EuclideanDistance,
    LinkageMeasure -> "Average",
    SignificanceCriterion -> "Silhouette",
    PrintDendrograms -> False, 
    ReturnDendrograms -> False,
    SingleAssociationLabel -> "1",
    SubclusteringDistanceFunction -> EuclideanDistance};
(***Function***)
TimeSeriesClusters[data_, OptionsPattern[]] :=
    Module[ {dataIn = data, 
      distFun1 = OptionValue[DistanceFunction], 
      distFun2 = OptionValue[SubclusteringDistanceFunction], 
      linkMeas = OptionValue[LinkageMeasure], 
      sigTest = OptionValue[SignificanceCriterion], 
      clusLabels = OptionValue[ClusterLabeling], 
      singleALabel = OptionValue[SingleAssociationLabel], 
      printDendro = OptionValue[PrintDendrograms], 
      dendroOptions = OptionValue[DendrogramPlotOptions], 
      dendroReturn = OptionValue[ReturnDendrograms], keysIn, dataVecs, 
      clusVec, numFirstClusters, agglomerativeCluster, 
      initialSplitClusters, subsplitClusVectors, numSubsplits, 
      intermediateClusters, subsplitClusters, clusteredData, 
      clusterAssocnOuts, auxiliaryLabels, 
      returning},(*defined Local Variables*)(*Cluster based on \
  similarity vectors, display output based on the data vectors: 
     auxilary variable clusVec*)
        If[ AssociationQ[dataIn],(*association Input*)
            If[ MemberQ[Keys@dataIn,"TimeSeriesClasses"],
                dataIn = dataIn["TimeSeriesClasses"]
            ];
            If[ AssociationQ[dataIn[[1]]],
                (*multiple Associations*)
                keysIn = Keys[dataIn];
                dataIn = Values[dataIn],
                keysIn = {singleALabel};
                dataIn = {dataIn}
            (*single association*)
            ],
            If[ AssociationQ[dataIn[[1]]],
                keysIn = Range[Length[dataIn]](*list of Associations*),
                keysIn = {singleALabel};
                dataIn = {AssociationThread[Range[Length[dataIn]], dataIn]}
            (*single list input*)
                ]
        ];
        clusVec = 
         Query[All, 
           KeyValueMap[((#2[[1]] /. 
                 Missing[] -> 0) -> (#2[[2]] -> #1)) &]]@dataIn;
        (*extract all the second values, i.e. the time series*)
        dataVecs = Values[dataIn][[All, All, 2]];
        (*we can use the FindClusters function,with the Silhouette option,
        to determine the number of clusters-
        but this FindClusters is NOT a cluster object*)
        numFirstClusters = 
         Length[FindClusters[#, DistanceFunction -> distFun1, 
             Method -> {"Agglomerate", "Linkage" -> linkMeas, 
               "SignificanceTest" -> sigTest}]] & /@ clusVec;
        (*now run agglomerative hierarchical clustering with all user \
     specified inputs*)
        agglomerativeCluster = 
         Agglomerate[#, DistanceFunction -> distFun1, 
            Linkage -> linkMeas] & /@ clusVec;
        (*and split into pre-determined number of clusters*)
        initialSplitClusters = 
         MapThread[
          If[ #2 == 1,
              {#1},
              ClusterSplit[#1, #2]
          ] &, {agglomerativeCluster,
            numFirstClusters}];
        (*now we take split clusters and recluster based on sign (make \
     sure to use case delete to convert to 0 any times with missing \
     points),and another user specified distance function*)
        (*we will \
     cluster based on caseDeletion of "Missing" but display output based \
     on data vectors*)
        (*Return[numFirstClusters];*)
        subsplitClusVectors = 
         If[ MatchQ[Head[#], Cluster],
             ClusterFlatten[# /. Missing[] -> 0][[All, 1]] -> 
              ClusterFlatten[#],
             {# /. Missing[] -> 0}[[All, 
             1]] -> {#}
         ] & /@ # & /@ initialSplitClusters;
        (*find separate list of subsplit levels to highlight*)
        
        (*Return[initialSplitClusters]*);
        numSubsplits = 
         Length[FindClusters[#, DistanceFunction -> distFun2, 
               Method -> {"Agglomerate", "Linkage" -> linkMeas, 
                 "SignificanceTest" -> sigTest}]] & /@ # & /@ 
          subsplitClusVectors;
        (*do agglomerative subclustering for each split cluster*)
        intermediateClusters = 
         Agglomerate[#, DistanceFunction -> distFun2, 
              Linkage -> linkMeas] & /@ # & /@ subsplitClusVectors;
        (*and again subsplit based on pre-
        determined number of subsplit clusters,unless there's only one*)
        subsplitClusters = 
         If[ MatchQ[Head[#[[1]]], Cluster],
             ClusterSplit,
             {#1} &
         ][#[[
               1]], #[[2]]] & /@ Transpose[#] & /@ 
          Transpose[{intermediateClusters, numSubsplits}];
        (*generate output clustered dataset for heatmap display based on \
     subsplits*)
        clusteredData = 
         Flatten[((If[ MatchQ[Head[#], Cluster],
                       ClusterFlatten[#],
                       {#}
                   ]) & /@ Flatten[#, 1]), 1] & /@ 
          subsplitClusters;
        clusterAssocnOuts = 
         Table[Association@{clusLabels <> "G" <> ToString[i] <> "S" <> 
                ToString[#] -> 
               If[ MatchQ[Head[subsplitClusters[[j]][[i, #]]], Cluster],
                   ClusterFlatten[subsplitClusters[[j]][[i, #]]][[All, 
                    2]],
                   {subsplitClusters[[j]][[i, #]][[2]]}
               ]} & /@ 
           Range[Length[subsplitClusters[[j]][[i]]]], {j, 1, 
           Length[subsplitClusters]}, {i, 1, 
           Length[subsplitClusters[[j]]]}];
        If[ printDendro || dendroReturn,
            auxiliaryLabels = 
             Table[clusLabels <> "G" <> ToString[i], {j, 1, 
               Length[subsplitClusters]}, {i, 1, 
               Length[subsplitClusters[[j]]]}];
            If[ dendroReturn,
                Return,
                Print
            ][
             AssociationThread[keysIn, 
              Association@(#[[1]] ->  
                     AssociationThread[#[[2]] -> #[[
                        3]]] & /@ #) & /@ (Transpose[#] & /@ 
                 Transpose[{auxiliaryLabels, Keys[clusterAssocnOuts], 
                   If[ MatchQ[Head[#], Cluster],
                       DendrogramPlot[#, 
                        Sequence @@ dendroOptions],
                       {#}
                   ] & /@ # & /@ # & /@ 
                    subsplitClusters}])]]
        ];
        returning = 
         AssociationThread[keysIn, 
          AssociationThread[{"Cluster", "InitialSplitCluster", 
              "IntermediateClusters", "SubsplitClusters", "Data", 
              "GroupAssociations"}, #] & /@ 
           Transpose[{agglomerativeCluster, initialSplitClusters, 
             intermediateClusters, subsplitClusters, clusteredData, 
             Association[#] & /@ clusterAssocnOuts}]];
        Return[returning]
    ];

(* ::Function:: *)
(* f:TimeSeriesSingleClusters *)
(*Hierarchical clustering of a time series without sub-clustring*)
(***Options***)
Options[TimeSeriesSingleClusters] = \
{ClusterLabeling -> "",
    DendrogramPlotOptions -> {},
    DistanceFunction -> EuclideanDistance, 
    LinkageMeasure -> "Average",
    PrintDendrograms -> False, 
    ReturnDendrograms -> False,
    SignificanceCriterion -> "Silhouette",
    SingleAssociationLabel -> "1" 
    };
(***Function***)
TimeSeriesSingleClusters[data_, OptionsPattern[]] :=
    Module[ {dataIn = data, 
      distFun1 = OptionValue[DistanceFunction], 
      linkMeas = OptionValue[LinkageMeasure], 
      sigTest = OptionValue[SignificanceCriterion], 
      clusLabels = OptionValue[ClusterLabeling], 
      singleALabel = OptionValue[SingleAssociationLabel], 
      printDendro = OptionValue[PrintDendrograms], 
      dendroOptions = OptionValue[DendrogramPlotOptions],
      dendroReturn = OptionValue[ReturnDendrograms], keysIn, dataVecs, 
      clusVec, numFirstClusters, agglomerativeCluster, 
      initialSplitClusters, clusterAssocnOuts, 
      returning},(*the layer above transforms all the external variables \
   into internal variables and defines all the other variables we are \
   going to need*)(*we can cluster based on the similarity vectors,
     but display the output based on the data vectors-
     to do so we define auxilary clusVec*)
        If[ AssociationQ[dataIn],(*association Input*)
            If[ MemberQ[Keys@dataIn, "TimeSeriesClasses"],
                dataIn = dataIn["TimeSeriesClasses"]
            ];
            If[ AssociationQ[dataIn[[1]]],
                       (*Multiple Associations*)
                keysIn = Keys[dataIn];
                dataIn = Values[dataIn],
                keysIn = {singleALabel};
                dataIn = {dataIn}
            (*"SingleAssociation"*)
                ],
            If[ AssociationQ[dataIn[[1]]],
                keysIn = Range[Length[dataIn]](*list of Associations*),
                keysIn = {singleALabel};
                dataIn = {AssociationThread[Range[Length[dataIn]], dataIn]}
            (*single list input*)
                ]
        ];
        clusVec = 
         Query[All, 
           KeyValueMap[((#2[[1]] /. Missing[] -> 0) -> (#2[[2]] -> #1)) &]]@
          dataIn;
        dataVecs = Values[dataIn][[All, All, 2]];
        (*we can use the FindClusters function,with the Silhouette option,
        to determine the number of clusters-
        but this FindClusters is NOT a cluster object*)
        numFirstClusters = 
         Length[FindClusters[#, DistanceFunction -> distFun1, 
             Method -> {"Agglomerate", "Linkage" -> linkMeas, 
               "SignificanceTest" -> sigTest}]] & /@ clusVec;
        (*now run agglomerative hierarchical clustering with all user \
      specified inputs*)
        agglomerativeCluster = 
         Agglomerate[#, DistanceFunction -> distFun1, 
            Linkage -> linkMeas] & /@ clusVec;
        (*and split into pre-determined number of clusters*)
        initialSplitClusters = 
         MapThread[
          If[ #2 == 1,
              {#1},
              ClusterSplit[#1, #2]
          ] &, {agglomerativeCluster, 
           numFirstClusters}];
        clusterAssocnOuts = 
         Table[Association@{clusLabels <> "G" <> ToString[i] -> 
             If[ MatchQ[Head[initialSplitClusters[[j]][[i]]], Cluster],
                 ClusterFlatten[initialSplitClusters[[j]][[i]]][[All, 
                  2]],
                 {initialSplitClusters[[j]][[i]][[2]]}
             ]}, {j, 1, 
           Length[initialSplitClusters]}, {i, 1, 
           Length[initialSplitClusters[[j]]]}];
        (*returning*)
        If[ printDendro || dendroReturn,
            If[ dendroReturn,
                Return,
                Print
            ][
             AssociationThread[
              keysIn, (AssociationThread[#[[1]] -> #[[2]]] & /@ 
                Transpose[{Keys[clusterAssocnOuts], 
                  If[ MatchQ[Head[#], Cluster],
                      DendrogramPlot[#, 
                       Sequence @@ dendroOptions],
                      {#}
                  ] & /@ # & /@ 
                   initialSplitClusters}])]]
        ];
        returning = 
         AssociationThread[keysIn, 
          AssociationThread[{"Cluster", "InitialSplitCluster", 
              "GroupAssociations"}, #] & /@ 
           Transpose[{agglomerativeCluster, initialSplitClusters, 
             Association[#] & /@ clusterAssocnOuts}]];
        Return[returning]
    ];

(* ::Function:: *)  
(* f:MatrixClusters *)
(*Hierarchical clustering,in both dimensions,of a data matrix*)
(***Options***)
Options[MatrixClusters] = {ClusterLabeling -> "",
    ColumnLabels -> {},
    DataTransforms -> {Identity, Identity},
    DendrogramPlotOptions -> {{}, {}},
    DistanceFunction -> EuclideanDistance,
    LinkageMeasure -> "Average",
    MissingSubstitution -> (Missing[] -> 0),
    PrintDendrograms -> False, 
    ReturnDendrograms -> False,
    SignificanceCriterion -> "Silhouette",
    SimilarityVectors -> {},
    SingleAssociationLabel -> "1"
    };
(***Function***)
MatrixClusters[data_, OptionsPattern[]] :=
    Module[ {dataIn = data, 
        simVecs = OptionValue[SimilarityVectors], 
        distFun = OptionValue[DistanceFunction], 
      linkMeas = OptionValue[LinkageMeasure], 
      sigTest = OptionValue[SignificanceCriterion],
      singleALabel = OptionValue[SingleAssociationLabel], 
      printDendro = OptionValue[PrintDendrograms],
      dendroOptions = OptionValue[DendrogramPlotOptions],
      dendroReturn = OptionValue[ReturnDendrograms],
      columnLabels = OptionValue[ColumnLabels],
      clusLabels = OptionValue[ClusterLabeling],
      dataTransforms = OptionValue[DataTransforms],
      missingSubstitution = OptionValue[MissingSubstitution], keysIn, 
      distFunH, distFunV, linkMeasH, linkMeasV, sigTestsH, sigTestsV, 
      horizClusVec, vertVecs, vertClusVec, numClustersHoriz, 
      agglomerativeClusterHoriz, splitClustersHoriz, numClustersVert, 
      agglomerativeClusterVert, splitClustersVert, clusterAssocnOutsH, 
      clusterAssocnOutsV, returning, simVecsH, 
      simVecsV},(*the layer above transforms all the external variables \
  into internal variables and defines all the other variables we are \
  going to need*)(*we can use the FindClusters function,
     with the Silhouette option,to determine the number of clusters-
     but this FindClusters is NOT a cluster object*)
        If[ Length[distFun] < 2,
            distFunH = distFunV = First@Flatten[{distFun}],
            distFunH = distFun[[1]];
            distFunV = distFun[[2]]
        ];
        If[ Length[linkMeas] < 2,
            linkMeasH = linkMeasV = First@Flatten[{linkMeas}],
            linkMeasH = linkMeas[[1]];
            linkMeasV = linkMeas[[2]]
        ];
        If[ Length[sigTest] < 2,
            sigTestsH = sigTestsV = First@Flatten[{sigTest}],
            sigTestsH = sigTest[[1]];
            sigTestsV = sigTest[[2]]
        ];
        If[ AssociationQ[dataIn],(*association Input*)
            If[ AssociationQ[dataIn[[1]]],
             (*Multiple Associations*)
                keysIn = Keys[dataIn];
                dataIn = Values[dataIn],
                keysIn = {singleALabel};
                dataIn = {dataIn}
            (*"SingleAssociation"*)
                ],
            If[ AssociationQ[dataIn[[1]]],
                keysIn = Range[Length[dataIn]](*list of Associations*),
                keysIn = {singleALabel};
                dataIn = {AssociationThread[Range[Length[dataIn]], dataIn]}
            (*single list input*)
                ]
        ];
        If[ MatchQ[simVecs, {}],
         (*data is used for clustering*)
            horizClusVec = 
             Query[All, 
               KeyValueMap[((dataTransforms[[
                      1]]@(#2 /. missingSubstitution)) -> (#2 -> #1)) &]]@
              dataIn;
            vertVecs = (Query[All, Values /* Transpose]@dataIn) /. 
              missingSubstitution;(*Extract a list of transposed vectors*)
            vertClusVec = (dataTransforms[[2]]@#[[1]]) -> #[[2]] & /@ 
               Transpose[#] & /@ 
             Transpose[{vertVecs, 
               If[ columnLabels == {},
                   Range[Length[#] & /@ vertVecs],
                   columnLabels
               ]}],
            simVecsH = simVecs[[1]];
            simVecsV = simVecs[[2]];
            If[ MatchQ[simVecsH, {}],
                horizClusVec = 
                 Query[All, 
                   KeyValueMap[((dataTransforms[[
                          1]]@(#2 /. missingSubstitution)) -> (#2 -> #1)) &]]@
                  dataIn,
                If[ AssociationQ[simVecsH],(*association similarity vectors Input*)
                    If[ AssociationQ[simVecsH[[1]]],
               (*Multiple Associations*)
                        simVecsH = Values[simVecsH],
                        (*"SingleAssociation"*)
                        simVecsH = {simVecsH}
                    ],
                    If[ !AssociationQ[simVecsH[[1]]],
                        simVecsH = {AssociationThread[Range[Length[simVecsH]], 
                           simVecsH]}
                    (*single list input*)
                        ]
                ];
                (*data is combined with simvecs based on associations using \
           Merge and then Queried to create the vectors for the clustering*)
                horizClusVec = 
                Query[All, 
                  KeyValueMap[((#2[[2]] /. 
                        missingSubstitution) -> (#2[[1]] -> #1)) &]]@(Merge[#,Identity] & /@ Transpose@{dataIn, simVecsH})
            ];
            (*Vertical*)
            If[ MatchQ[simVecsV, {}],
             (*Vertical Empty, use Data*)
                vertVecs = (Query[All, Values /* Transpose]@dataIn) /. 
                  missingSubstitution;(*Extract a list of transposed vectors*)
                vertClusVec = (dataTransforms[[2]]@#[[1]]) -> #[[2]] & /@ 
                 Transpose[#] & /@ 
                Transpose[{vertVecs, 
                 If[ columnLabels == {},
                     Range[Length[#] & /@ vertVecs],
                     columnLabels
                 ]}],
                If[ AssociationQ[simVecsV],(*association similarity vectors Input*)
                    If[ AssociationQ[simVecsV[[1]]],
               (*Multiple Associations*)
                        vertVecs = (Values@Values@simVecsV) /. missingSubstitution;
                        If[ MatchQ[columnLabels, {}],
                            columnLabels = Keys@Values@simVecsV
                        ],
                        (*"SingleAssociation"*)
                        vertVecs = {Values@simVecsV} /. missingSubstitution;
                        If[ MatchQ[columnLabels, {}],
                            columnLabels = {Keys@simVecsV}
                        ]
                    ],
                    If[ AssociationQ[simVecsV[[1]]],
                     (*list of Associations*)
                        vertVecs = (Values@simVecsV) /. missingSubstitution;
                        If[ MatchQ[columnLabels, {}],
                            columnLabels = Keys@simVecsV
                        ],
                        vertVecs = {simVecsV} /. missingSubstitution;
                        If[ MatchQ[columnLabels, {}],
                            columnLabels = Range[Length[#] & /@ vertVecs]
                        ]
                    (*single list input*)
                        ]
                ];
                vertClusVec = #[[1]] -> #[[2]] & /@ Transpose[#] & /@ 
                  Transpose[{vertVecs, 
                    If[ MatchQ[columnLabels, {}],
                        Range[Length[#] & /@ vertVecs],
                        columnLabels
                    ]}]
            (*Use simvecs*)
                ]
        ];
        (*find the number of horizontal clusters*)
        numClustersHoriz = 
         Length@FindClusters[#, DistanceFunction -> distFunH, 
             Method -> {"Agglomerate", "Linkage" -> linkMeasH, 
               "SignificanceTest" -> sigTestsH}] & /@ horizClusVec;
        (*find the number of vertical clusters*)
        numClustersVert = 
         Length@FindClusters[#, DistanceFunction -> distFunV, 
             Method -> {"Agglomerate", "Linkage" -> linkMeasV, 
               "SignificanceTest" -> sigTestsV}] & /@ vertClusVec;
        (*now run agglomerative hierarchical clustering with all user \
     specified inputs*)
        agglomerativeClusterHoriz = 
         Agglomerate[#, DistanceFunction -> distFunH, 
            Linkage -> linkMeasH] & /@ horizClusVec;
        agglomerativeClusterVert = 
         Agglomerate[#, DistanceFunction -> distFunV, 
            Linkage -> linkMeasV] & /@ vertClusVec;
        (*and split into pre-determined number of clusters*)
        splitClustersHoriz = 
         MapThread[
          If[ #2 == 1,
              {#1},
              ClusterSplit[#1, #2]
          ] &, {agglomerativeClusterHoriz, 
           numClustersHoriz}];
        splitClustersVert = 
         MapThread[
          If[ #2 == 1,
              {#1},
              ClusterSplit[#1, #2]
          ] &, {agglomerativeClusterVert, 
           numClustersVert}];
        clusterAssocnOutsH = 
         Table[Association@{clusLabels <> "H:G" <> ToString[i] -> 
             If[ MatchQ[Head[splitClustersHoriz[[j]][[i]]], Cluster],
                 ClusterFlatten[splitClustersHoriz[[j]][[i]]][[All, 
                  2]],
                 {splitClustersHoriz[[j]][[i]][[2]]}
             ]}, {j, 1, 
           Length[splitClustersHoriz]}, {i, 1, 
           Length[splitClustersHoriz[[j]]]}];
        clusterAssocnOutsV = 
         Table[Association@{clusLabels <> "V:G" <> ToString[i] -> 
             If[ MatchQ[Head[splitClustersVert[[j]][[i]]], Cluster],
                 ClusterFlatten[
                  splitClustersVert[[j]][[i]]],
                 {splitClustersVert[[j]][[
                 i]]}
             ]}, {j, 1, Length[splitClustersVert]}, {i, 1, 
           Length[splitClustersVert[[j]]]}];
        (*printDendrograms*)
        If[ printDendro || dendroReturn,
            If[ dendroReturn,
                Return,
                Print
            ][
             AssociationThread[keysIn, 
              Join @@ {#[[1]], #[[
                   2]]} & /@ (Transpose@{(AssociationThread[#[[1]] -> #[[
                        2]]] & /@ 
                    Transpose[{Keys[clusterAssocnOutsH], 
                      If[ MatchQ[Head[#], Cluster],
                          DendrogramPlot[#, 
                           Sequence @@ dendroOptions[[1]]],
                          {#}
                      ] & /@ # & /@ 
                       splitClustersHoriz}]), (AssociationThread[#[[1]] -> #[[
                        2]]] & /@ 
                    Transpose[{Keys[clusterAssocnOutsV], 
                      If[ MatchQ[Head[#], Cluster],
                          DendrogramPlot[#, 
                           Sequence @@ dendroOptions[[2]]],
                          {#}
                      ] & /@ # & /@ 
                       splitClustersVert}])})]]
        ];
        returning = 
         AssociationThread[keysIn, 
          AssociationThread[{"RowCluster", "ColumnCluster", 
              "RowSplitClusters", "ColumnSplitClusters", 
              "GroupAssociationsRows", "GroupAssociationsColumns"}, #] & /@
            Transpose[{agglomerativeClusterHoriz, agglomerativeClusterVert,
              splitClustersHoriz, splitClustersVert, 
             Association[#] & /@ clusterAssocnOutsH, 
             Association[#] & /@ clusterAssocnOutsV}]];
        Return[returning]
    ];
(*#####Databases#####*)
(* ::Function:: *)
(* f:UCSCBrowserSQL*)
(*UCSC SQL query*)
(***Function***)
UCSCBrowserSQL[inSQLString_] :=
    Module[ {inputString = inSQLString, ucscDatabase, 
      outputData},(*the layer above transforms all the external \
  variables into internal variables and defines all the other variables \
  we are going to need*)(*get the database from UCSC*)
        ucscDatabase = 
         OpenSQLConnection[
          JDBC["MySQL(Connector/J)", "genome-mysql.cse.ucsc.edu"], 
          "Username" -> "genomep", "Password" -> "password"];
        (*get the output from the user input string*)
        outputData = SQLExecute[ucscDatabase, inputString];
        (*close SQL connection*)
        CloseSQLConnection[ucscDatabase];
   (*and finally return it*)
        Return[outputData]
    ];
    
(* ::Section:: *)
(*#####DataProcessing#####*)

(****Data Import****)

(* ::Function:: *)
(* f:OmicsObjectCreator:: *)
(***function***)
OmicsObjectCreator[outerLabels_, innerLabels_, measurements_, 
  metadata_] :=
    Module[ {outLabels = outerLabels, inLabels = innerLabels, 
      dataIn = measurements, meta = metadata},
        AssociationThread[
          outLabels -> 
           KeyUnion@(AssociationThread[#] & /@ (#[[1]] -> #[[2]] & /@ \
        (Transpose@{inLabels, 
                   Transpose[#] & /@ 
                    Transpose@{dataIn, meta}})))] /. _Missing -> 
          Missing[]
    ]


(* ::Function:: *)
(* f:FileSelector *)
Options[FileSelector] = {DelimiterList -> {","}};
(*delimiters for possible text import*)
(***Attributes***)
SetAttributes[FileSelector, HoldFirst];
(*invalid format messager*)
(***Messages***)
FileSelector::argx = 
  "Not a valid format. Please re-evaluate Input Cell, and use a \
.csv, .xlsx, or delimited .txt file with headers";
(***Function***)
FileSelector[variableName_, OptionsPattern[]] :=
    Module[ {files, importedFiles, 
      delimiters = OptionValue[DelimiterList]},
     (*import button*)
        Button["File Import: Select your file(s)", 
         files = SystemDialogInput[
           "FileOpen", {NotebookDirectory[], {"All Files" -> {"*"},"Comma Separated Values \
(.csv)" -> {"*.csv"}, "Excel (.xlsx)" -> {"*.xlsx"}, 
             "Text (.txt)" -> {"*.txt"}}}, 
           WindowTitle -> 
            "Import Data Files"];(*open a file dialogue to select \
      appropriate files*)
         If[ files =!= $Canceled && files =!= $Failed,
          (*import first line as header to use*)
             importedFiles = 
              Switch[FileFormat[#], "CSV", Import[#, {"Data", 1, All}], "TSV", 
                 Import[#, {"Data", 1, All}], "XLSX", 
                 Import[#, {"Data", 1, 1}], "Text", 
                 StringSplit[ReadString[#, "\n"], delimiters], _, 
                 Message[FileSelector::argx, ""]] & /@ Flatten[{files}];
             variableName = Transpose[{importedFiles, Flatten[{files}]}],
             importedFiles = ""
         ], Method -> "Queued"]
    ];
    
(* ::Function:: *)
(* f:DataImporterDirectLabeled *)
(***Options***)
Options[DataImporterDirectLabeled] = {AssociationLabels -> None, 
  SampleKind -> {}};
(***Function***)
DataImporterDirectLabeled[sampleRules_, fileList_, headerLines_, 
  headerColumnAssociations_, OptionsPattern[]] :=
    Module[ {varName = sampleRules, fiNames = fileList, 
      hdrLines = headerLines, interA = headerColumnAssociations,
      aLabels = OptionValue[AssociationLabels], 
      sampleKind = Flatten@{OptionValue[SampleKind]}, returning, dataIn, 
      associationIn},
        DistributeDefinitions[varName,interA,fiNames,hdrLines,sampleKind];
        dataIn = 
         Parallelize@MapIndexed[Switch[FileFormat[fiNames[[#2[[1]]]]],
             "CSV", 
             Import[fiNames[[#2[[1]]]], {"Data", All, 
               If[ Length[#] == 1,
                   Flatten[#, 1],
                   #
               ] &@
                MapIndexed[interA[[#2[[1]]]][#1] &, varName, {5}][[All, 
                  2]][[#2[[1]]]]}, HeaderLines -> hdrLines],
             "TSV", 
             Import[fiNames[[#2[[1]]]], {"Data", All, 
               If[ Length[#] == 1,
                   Flatten[#, 1],
                   #
               ] &@
                
                MapIndexed[interA[[#2[[1]]]][#1] &, varName, {5}][[All, 
                  2]][[#2[[1]]]]}, HeaderLines -> hdrLines],
             "XLSX", 
             Import[fiNames[[#2[[1]]]], {"Data", 1, All, 
                If[ Length[#] == 1,
                    Flatten[#, 1],
                    #
                ] &@
                 MapIndexed[interA[[#2[[1]]]][#1] &, varName, {5}][[All, 
                   2]][[#2[[1]]]]}][[(1 + hdrLines) ;;]],
             "Text", 
             Import[fiNames[[#2[[1]]]], {"Data", All, 
               If[ Length[#] == 1,
                   Flatten[#, 1],
                   #
               ] &@
                MapIndexed[interA[[#2[[1]]]][#1] &, varName, {5}][[All, 
                  2]][[#2[[1]]]]}, HeaderLines -> hdrLines],
             _,
             Message[FileSelector::argx, ""];
             Print["Invalid File Format ", 
              FileFormat[fiNames[[#2[[1]]]]]]] &, interA];
        (*dataIn imported correct columns based on input file type*)
        \
      (*create the association using input association names for samples*)
        associationIn = 
        Parallelize[
         AssociationThread[#[[1]] -> #[[2]]] & /@ 
          Transpose[{MapIndexed[interA[[#2[[1]]]][#1] &, varName, {5}][[
             All, 1]], 
            Map[If[ Depth[#] == 
                
                4,
                    {Association@
                    Map[Join[#[[1]], sampleKind] -> #[[2 ;;]] &, #]},
                    Association[#] & /@ 
                     Map[Join[#[[1]], sampleKind] -> #[[2 ;;]] &, 
                      Transpose[#], {2}]
                ] &, dataIn]}]];
       (*create main result by unionizing keys across the input sample \
     points*)
        returning = 
         Query[All, All, 
           If[ MissingQ[#],
               Missing[],
               If[ Count[Flatten[{#[[1]]}], ""] != 0,
                   Missing[],
                   #
               ]
           ] &]@(AssociationThread[
              Keys[#] -> KeyUnion[Values[#]]] &@(Join @@ associationIn));
        If[ MatchQ[aLabels, None],
            Return[returning],
            (*if there are labels then rename the output keys*)
            returning = 
             KeyMap[AssociationThread[Keys[returning], aLabels]]@returning
        ];
        Return[returning]
    ];

(* ::Function:: *)
(* f:DataImporterDirect *)
(*Direct Association Creator from file Indexes and File Names*)
(*take as input direct column specification (multi-file) to extract features*)
(***Options***)
Options[DataImporterDirect] = {AssociationLabels -> None, 
  SampleKind -> {}};
(***Function***)
DataImporterDirect[positionsList_, fileList_, headerLines_, 
  OptionsPattern[]] :=
    Module[ {varName = positionsList, fiNames = fileList, 
      hdrLines = headerLines,
      aLabels = OptionValue[AssociationLabels], 
      sampleKind = Flatten@{OptionValue[SampleKind]}, returning, dataIn, 
      associationIn},
        DistributeDefinitions[varName,fiNames,hdrLines,sampleKind];
        dataIn = 
         Parallelize@MapIndexed[Switch[FileFormat[fiNames[[#2[[1]]]]],
             "CSV", 
             Import[fiNames[[#2[[1]]]], {"Data", All, 
               If[ Length[#] == 1,
                   Flatten[#, 1],
                   #
               ] &@varName[[#2[[1]]]]}, 
              HeaderLines -> hdrLines],
             "TSV", 
             Import[fiNames[[#2[[1]]]], {"Data", All, 
               If[ Length[#] == 1,
                   Flatten[#, 1],
                   #
               ] &@varName[[#2[[1]]]]}, 
              HeaderLines -> hdrLines],
             "XLSX", 
             Import[fiNames[[#2[[1]]]], {"Data", 1, All, 
                If[ Length[#] == 1,
                    Flatten[#, 1],
                    #
                ] &@
                 varName[[#2[[1]]]]}][[(1 + hdrLines) ;;]],
             "Text", 
             Import[fiNames[[#2[[1]]]], {"Data", All, 
               If[ Length[#] == 1,
                   Flatten[#, 1],
                   #
               ] &@varName[[#2[[1]]]]}, 
              HeaderLines -> hdrLines],
             _,
             Message[FileSelector::argx, ""];
             Print["Invalid File Format ", 
              FileFormat[fiNames[[#2[[1]]]]]]] &, 
           Range[(Length[#] &@varName)]];
        (*dataIn imported correct columns based on input file type*)
        \
      (*create the association using indexed association names for files <i>
        f<file number>*)
        associationIn = 
         Parallelize[
          AssociationThread[#[[1]] -> #[[2]]] & /@ 
           Transpose[{MapIndexed[
              Table["Other" <> ToString[i] <> "f" <> ToString[#2[[1]]], {i, 
                 1, Length[#1]}] &, varName], 
             Map[If[ Depth[#] == 4,
                     {Association@
                     Map[Join[#[[1]], sampleKind] -> #[[2 ;;]] &, #]},
                     Association[#] & /@ 
                      Map[Join[#[[1]], sampleKind] -> #[[2 ;;]] &, 
                       Transpose[#], {2}]
                 ] &, dataIn]}]];
        (*create main result by unionizing keys across the input sample \
      points*)
        returning = 
         Query[All, All, 
           If[ MissingQ[#],
               Missing[],
               If[ Count[Flatten[{#[[1]]}], ""] != 0,
                   Missing[],
                   #
               ]
           ] &]@(AssociationThread[
              Keys[#] -> KeyUnion[Values[#]]] &@(Join @@ associationIn));
        If[ MatchQ[aLabels, None],
            Return[returning],
            (*if there are labels then rename the output keys*)
            returning = 
             KeyMap[AssociationThread[Keys[returning], aLabels]]@returning
        ];
        Return[returning]
    ]  

(* ::Function:: *)
(* f:PanelImport *)
(*Panel to import data: the panel lets the user select columns if there is a header line in the files*)
PanelImport::usage = "PanelImport[variableName, input] imports data Using a Panel, used with DataImporter.";
(***Options***)
Options[PanelImport] = {VariableIndex -> ""};
(***Attributes***)
SetAttributes[PanelImport, HoldFirst];
(***Function***)
PanelImport[variableName_, input_, OptionsPattern[]] :=
    DynamicModule[{inputAssociationNames = 
       MapIndexed[#1 -> #1 &, #] &@(input), usedInternal = input, 
      fButtonNotPressed = True, identifier = {1, 1, 1}, 
      variableIndex = ToString@OptionValue[VariableIndex], 
      labels = {"Other" <> ToString@OptionValue[VariableIndex] <> "1"}, 
      counter = 1},
     Panel[TabView[{"Sample Labeling" -> 
         TabView[{"Current Selection of labels" -> 
            Column[{Dynamic@
               Grid@ MapIndexed[{ToString[#1], 
                   InputField[Dynamic[labels[[#2[[1]]]]], String]} &, 
                 labels], 
              Button["Add Custom Field", counter++;
                                         AppendTo[labels, 
                                          "Other" <> variableIndex <> ToString[counter]]], 
              Button["Remove Last Field", 
               labels = Setting[labels[[1 ;; -2]]]],
              DynamicWrapper["", usedInternal = Join[input, labels]]}], 
           "Select from file header" -> 
            Grid@Transpose@{{Dynamic@
                 CheckboxBar[Dynamic[labels], 
                  MapIndexed[#1 -> #1 &, #] &@usedInternal, 
                  Appearance -> "Row"]}}}] , 
        "Measurement Selection" -> Row@{
           Dynamic[
            identifier = 
             ConstantArray["Pending", {Length[Setting[labels]], 3}];
            ""],
            Dynamic@
            TabView[MapIndexed[#1 -> 
                TabView[{"Unique Id" -> 
                   CheckboxBar[Dynamic@identifier[[#2[[1]], 1]], 
                    inputAssociationNames, Appearance -> "Row"],
                  
                  "Measurements" -> 
                   CheckboxBar[Dynamic@identifier[[#2[[1]], 2]], 
                    inputAssociationNames, Appearance -> "Row"],
                  
                  "Extras" -> 
                   CheckboxBar[Dynamic@identifier[[#2[[1]], 3]], 
                    inputAssociationNames, Appearance -> "Row"]}, 
                 ImageSize -> Automatic] &, labels]]}, 
        "Finalize" -> 
         Column@{"Push Button to Finalize Selection", 
           Dynamic@{labels, identifier},
           Button["Set Input Measurement", 
            fButtonNotPressed = Dynamic[False];
            variableName = 
             DeleteCases[Setting[labels -> identifier], "Pending", 4]], 
           Dynamic@If[ Setting@fButtonNotPressed,
                       "Selection Pending",
                       "Selection Made"
                   ]}}]] ];    

(* ::Function:: *)
(* f:DataImporter *)
(*Importer of Data*)
(*graphical interface to extract data and create \
variable for association of information*)
(***Options***)
Options[DataImporter] = {AssociationLabels -> None,
	DelimiterList -> {","},
	JavaGBs -> None};
(***Attributes***)
SetAttributes[DataImporter, HoldFirst];
(***Function***)
DataImporter[associationName_, OptionsPattern[]] :=
    If[ ! MemberQ[Keys@NotebookInformation[], "FileName"],
        Print["Please retry after Notebook has been saved "];
        NotebookSave[],
        Style[DynamicModule[{variableName, fileNames, 
           delimiters = OptionValue[DelimiterList], 
           aLabels = OptionValue[AssociationLabels], 
           javaReinstall = OptionValue[JavaGBs], 
           idButtonNotPressed = True, importButtonNotPressed = True, first, 
           intermediateA, global, radioButton, 
           radioButtonManual, headerLines = 1, sampleKind = {}, 
           globalSelection = {{{"identifier1"}, {"measurement1"}, \
        {"metadata1"}}, {{"identifier2"}, {"measurement2"}, {"metadata2"}}}, 
           individualSelection = {{{{"identifier1f1"}, {"measurement1f1"}, \
        {"metadata1f1"}}, {{"identifier2f1"}, {"measurement2f1"}, \
        {"metadata2f1"}}, {{"identifier3f1"}, {"measurement3f1"}, \
        {"metadata3f1"}}}, {{{"identifier1f2"}, {"measurement1f2"}, \
        {"metadata1f2"}}, {{"identifier2f2"}, {"measurement2f2"}, \
        {"metadata2f2"}}}}},
          If[ !MatchQ[javaReinstall, None],(*JLink` introduced in package*)
              ReinstallJava[
               JVMArguments -> "-Xmx" <> ToString@javaReinstall <> "g"]
          ];
          Panel[
           Grid@{{Style["MathIOmica Data Importer", Bold]}, {Panel@Dynamic[
                
                If[ Setting[importButtonNotPressed],
                    If[ Setting[idButtonNotPressed],
                        Column[{DynamicWrapper[
                           Style["Selection of Files", Italic], first], 
                          FileSelector[first, DelimiterList -> delimiters],
                          If[ ListQ[first],
                              DynamicWrapper["", 
                               idButtonNotPressed = Dynamic[False]]
                          ]}](*file selection completed*),
                        Column[{DynamicWrapper[
                           Style["Select Method for Choosing Data Columns", 
                            Italic], first], variableName = Range[Length[first]];
                                             fileNames = first[[All, 2]];
                                             Panel@RadioButtonBar[Dynamic[radioButton],
                                               (*Global Selection of Columns*)
                                               \
                                             {Grid[{{Style[
                                                   "Selections for All Files Will be Uniform based \
on first selected File", Italic]},
                                                   {"PathName: " <> first[[1, 2]]},
                                                   {PanelImport[global, first[[1, 1]]]},
                                                   {Row@{"Number of Header Lines: ", 
                                                   InputField[Dynamic[headerLines]]}},
                                                   {Row@{"Sample Kind (e.g. \"RNA\"): ", 
                                                   InputField[
                                                   Dynamic[sampleKind]]}}, {Grid[{{Button[
                                                   "Reset File Selection", 
                                                   globalSelection = {{{"identifier1"}, \
                                             {"measurement1"}, {"metadata1"}}, {{"identifier2"}, {"measurement2"}, \
                                             {"metadata2"}}};
                                                   individualSelection = {{{{"identifier1f1"}, \
                                                   {"measurement1f1"}, {"metadata1f1"}}, {{"identifier2f1"}, \
                                                   {"measurement2f1"}, {"metadata2f1"}}, {{"identifier3f1"}, \
                                                   {"measurement3f1"}, {"metadata3f1"}}}, {{{"identifier1f2"}, \
                                                   {"measurement1f2"}, {"metadata1f2"}}, {{"identifier2f2"}, \
                                                   {"measurement2f2"}, {"metadata2f2"}}}};
                                                   (*reset variables for re-selection*)
                                                   ClearAll[first, radioButton, variableName];
                                                   idButtonNotPressed = Dynamic[True], 
                                                   Method -> "Queued"],
                                                   Button["Import and Create Dataset",
                                                   
                                                   intermediateA = 
                                                   Association@
                                                   MapIndexed[#1 -> #2[[1]] &, #[[1]]] & /@ first;
                                                   (*create indexing of variables*)
                                                   variableName = #[[1]] -> #[[2]] & /@ 
                                                   Transpose[{Table[# <> "f" <> ToString[i] & /@ 
                                                   Keys[global], {i, 1, Length[fileNames]}], 
                                                   ConstantArray[Values[global], Length[fileNames]]}];
                                                (*create global variable*)
                                                   associationName = 
                                                   Setting@
                                                   DataImporterDirectLabeled[Setting@variableName, 
                                                   Setting@fileNames, Setting@headerLines, 
                                                   Setting@intermediateA, 
                                                   AssociationLabels -> Setting@aLabels, 
                                                   SampleKind -> Setting@sampleKind];
                                                   importButtonNotPressed = Dynamic[False], 
                                                   Method -> "Queued"
                                                   ]}}, Spacings -> 2]}}, Alignment -> Left] -> 
                                                 "Global Selection",
                                                (*Itemized Selection of Columns per File*)
                                                Grid[{{( 
                                                   MenuView[#, ImageSize -> 800] &@
                                                   MapIndexed[(#2[[1]] -> FileNameTake[#1[[2]]]) -> 
                                                   Grid[{{"PathName: " <> #1[[2]]}, { 
                                                   PanelImport[variableName[[#2[[1]]]], #1[[1]], 
                                                   VariableIndex -> #2[[1]]]}}, Alignment -> Left] &,
                                                    first])},
                                                   (*input for number of lines to skip as headers*)
                                             \
                                                  {Row@{"Number of Header Lines: ", 
                                                   InputField[Dynamic[headerLines]]}},
                                                   {Row@{"Sample Kind (e.g. \"RNA\"): ", 
                                                   InputField[Dynamic[sampleKind]]}},
                                                   (*reset file selection button*){Grid[{{Button[
                                                   "Reset File Selection", 
                                                   globalSelection = {{{"identifier1"}, \
                                             {"measurement1"}, {"metadata1"}}, {{"identifier2"}, {"measurement2"}, \
                                             {"metadata2"}}};
                                                   individualSelection = {{{{"identifier1f1"}, \
                                                   {"measurement1f1"}, {"metadata1f1"}}, {{"identifier2f1"}, \
                                                   {"measurement2f1"}, {"metadata2f1"}}, {{"identifier3f1"}, \
                                                   {"measurement3f1"}, {"metadata3f1"}}}, {{{"identifier1f2"}, \
                                                   {"measurement1f2"}, {"metadata1f2"}}, {{"identifier2f2"}, \
                                                   {"measurement2f2"}, {"metadata2f2"}}}};
                                                   ClearAll[first, radioButton, variableName];
                                                   idButtonNotPressed = Dynamic[True], 
                                                   Method -> "Queued"],
                                                   Button["Import and Create Dataset",            
                                                   intermediateA = 
                                                   Association@
                                                   MapIndexed[#1 -> #2[[1]] &, #[[1]]] & /@ first;
                                                   associationName = 
                                                   Setting@
                                                   DataImporterDirectLabeled[Setting@variableName, 
                                                   Setting@fileNames, Setting@headerLines, 
                                                   Setting@intermediateA, 
                                                   AssociationLabels -> Setting@aLabels, 
                                                   SampleKind -> Setting@sampleKind];
                                                   importButtonNotPressed = Dynamic[False], 
                                                   Method -> "Queued"
                                                   ]}}, Spacings -> 2]}}, Alignment -> Left] -> 
                                                 "Itemized Selection",
                                                (*Expert Manual Selection*)
                                                Panel@Grid[{{RadioButtonBar[
                                                   Dynamic[
                                                   radioButtonManual], {Grid[{{"Import columns \
numbering ", InputField[Dynamic[globalSelection]]},
                                                   {Row@{"Number of Header Lines: ", 
                                                   InputField[
                                                   Dynamic[
                                                   headerLines]]}}, {Row@{"Sample Kind (e.g. \
\"RNA\"): ", InputField[
                                                   Dynamic[sampleKind]]}}, {Grid[{{Button[
                                                   "Reset File Selection", 
                                                   globalSelection = {{{"identifier1"}, \
                                             {"measurement1"}, {"metadata1"}}, {{"identifier2"}, {"measurement2"}, \
                                             {"metadata2"}}};
                                                   individualSelection = {{{{"identifier1f1"}, \
                                                   {"measurement1f1"}, {"metadata1f1"}}, {{"identifier2f1"}, \
                                                   {"measurement2f1"}, {"metadata2f1"}}, {{"identifier3f1"}, \
                                                   {"measurement3f1"}, {"metadata3f1"}}}, {{{"identifier1f2"}, \
                                                   {"measurement1f2"}, {"metadata1f2"}}, {{"identifier2f2"}, \
                                                   {"measurement2f2"}, {"metadata2f2"}}}};
                                                   ClearAll[first, radioButton, variableName];
                                                   idButtonNotPressed = Dynamic[True], 
                                                   Method -> "Queued"],
                                                   Button["Import and Create Dataset",
                                                   
                                                   associationName = 
                                                   Setting@
                                                   DataImporterDirect[
                                                   Setting[
                                                   ConstantArray[globalSelection, 
                                                   Length[fileNames]]], Setting@fileNames, 
                                                   Setting@headerLines, 
                                                   AssociationLabels -> Setting@aLabels, 
                                                   SampleKind -> Setting@sampleKind];
                                                   importButtonNotPressed = Dynamic[False], 
                                                   Method -> "Queued"
                                                   ]}}, Spacings -> 2]}}, Alignment -> Left] -> 
                                                   "Global Selection", 
                                                   Grid[{{"Import columns numbering ", 
                                                   InputField[Dynamic[individualSelection]]},
                                                   {Row@{"Number of Header Lines: ", 
                                                   InputField[
                                                   Dynamic[
                                                   headerLines]]}}, {Row@{"Sample Kind (e.g. \
\"RNA\"): ", InputField[
                                                   Dynamic[sampleKind]]}}, {Grid[{{Button[
                                                   "Reset File Selection", 
                                                   globalSelection = {{{"identifier1"}, \
                                             {"measurement1"}, {"metadata1"}}, {{"identifier2"}, {"measurement2"}, \
                                             {"metadata2"}}};
                                                   individualSelection = {{{{"identifier1f1"}, \
                                                   {"measurement1f1"}, {"metadata1f1"}}, {{"identifier2f1"}, \
                                                   {"measurement2f1"}, {"metadata2f1"}}, {{"identifier3f1"}, \
                                                   {"measurement3f1"}, {"metadata3f1"}}}, {{{"identifier1f2"}, \
                                                   {"measurement1f2"}, {"metadata1f2"}}, {{"identifier2f2"}, \
                                                   {"measurement2f2"}, {"metadata2f2"}}}};
                                                   ClearAll[first, radioButton, variableName, 
                                                   associationName];
                                                   idButtonNotPressed = Dynamic[True], 
                                                   Method -> "Queued"],
                                                   Button["Import and Create Dataset",
                                                   
                                                   associationName = 
                                                   Setting@
                                                   DataImporterDirect[
                                                   Setting@individualSelection, Setting@fileNames, 
                                                   Setting@headerLines, 
                                                   AssociationLabels -> Setting@aLabels, 
                                                   SampleKind -> Setting@sampleKind];
                                                   importButtonNotPressed = Dynamic[False], 
                                                   Method -> "Queued"
                                                   ]}}, Spacings -> 2]}}, Alignment -> Left] -> 
                                                   " Individual Selection"}]}, {Dynamic[
                                                   radioButtonManual]}}, Alignment -> Left] -> 
                                                 "Manual Selection"}
                                               ], Dynamic[radioButton]}
                         ]
                    ],(*confirmation of variable creation*)
                    Column@{"CreatedAssociation: ", 
                      ToString@SymbolName@Unevaluated@associationName }
                ]
                ]}}
           ]
          ], DynamicEvaluationTimeout -> Infinity]
    ];


(****Bootstap Assistive Functions****)

(* ::Function:: *)
(* f:BootstrapGeneral *)
(***Options***)
Options[BootstrapGeneral] = {LabelFunction -> Range, 
  SamplingFunction -> RandomChoice};
(***Function***)
BootstrapGeneral[data_, sampling_, OptionsPattern[]] :=
    Module[ {inputdata = data, number = sampling, 
      samplingFunction = OptionValue[SamplingFunction], 
      labelFunction = OptionValue[LabelFunction], returning}, 
     (*perform a resampling of the data with replacement, 
     and generate a new association structure with numbering \
   corresponding to new identities - 
     this is necessary for constructing the time series, or using non-
     standard identities*)
        returning = (Query[All, 
            All /* Values /* (AssociationThread[
                labelFunction[number] -> samplingFunction[#, number]] &)]@
           inputdata)/.(_Missing-> Missing[]);
        Return[returning]
    ]

(* ::Function:: *)
(* f:QuantileEstimator *)
(*Note - data is from constant time series, inputs as "name, timeseries"*)
(***Options***)
Options[QuantileEstimator] = {AutocorrelationOptions -> {}, 
   InterpolationDeltaT -> "Auto",
   InterpolationOptions -> {},
   LombScargleOptions -> {PairReturn -> False,NormalizeIntensities->True}, 
   Method -> "LombScargle"(*Autocorrelation*),
   QuantileValue -> 0.95};
(***Function***)
QuantileEstimator[data_, timepoints_, OptionsPattern[]] :=
    Module[ {dataIn = data,
      timepts = timepoints,
      q = OptionValue[QuantileValue],
      optsAutoCorr = OptionValue[AutocorrelationOptions],
      optsLS = OptionValue[LombScargleOptions],
      type = OptionValue[Method],
      iDeltaT = OptionValue[InterpolationDeltaT],
      interOptions = OptionValue[InterpolationOptions],
      dataValues,
      iTimePts,
      returning},
        dataValues = If[ AssociationQ[dataIn],
                         Values@dataIn,
                         dataIn
                     ];
        Switch[type,
         "Autocorrelation",
         returning = 
          Quantile[Cases[N[#],_Real],q] & /@ ((Transpose@(Autocorrelation[#, timepts, 
                   Sequence @@ optsAutoCorr] & /@ (dataValues)))[[2 ;;]]),
         "LombScargle",
         returning = 
          Quantile[Cases[N[#],_Real],q] &@
        (Max@(LombScargle[#, timepts, 
                    Sequence @@ optsLS])[[2]] & /@ (dataValues)),
        "Spikes",
        returning = (Query[All, 
        All /* Transpose /* ({Quantile[Cases[N@#[[1]],_Real], q], 
        Quantile[Cases[N@#[[2]],_Real], 1-q]} &), ({Max[#], Min[#]} &@
        Cases[#, _Real] &)]@ 
        GroupBy[(dataValues), Length[DeleteCases[#, Except[_Real]]] &]),
        "InterpolatedAutocorrelation",
        iTimePts = Table[i, {i, Min[#], Max[#], 
              If[ MatchQ[iDeltaT, "Auto"],
                  (Max[#] - Min[#])/(Length[#] - 1),
                  iDeltaT
              ]}] &@timepts;
        inputTLength = Length@iTimePts;
        returning = Quiet@(Quantile[Cases[N[#],_Real],q] & /@ (Transpose@((CorrelationFunction[#, {1, 
               inputTLength - 
            1}] &@((Interpolation[#, Sequence @@ interOptions] &@
                N@DeleteCases[{x_, y_} /; (MissingQ[x] || MissingQ[y])]@
                  Transpose[{timepts, #}])[iTimePts])) & /@ dataValues))),
         __,
         Return[
          "Type\[Rule]" <> ToString[type] <> 
           "is not a valid method. Please make another selection"]
         ];
        Return[returning]
    ];


(****Data processing & filtering****)

(* ::Function:: *)
(* f:FilteringFunction *)
(***Options***)
Options[FilteringFunction] = {ComponentIndex -> Missing[],
	ListIndex -> Missing[], 
	SelectionFunction -> GreaterEqual};
(***Function***)
FilteringFunction[omicsObject_, cutoff_, OptionsPattern[]] :=
    Module[ {data = omicsObject, outputData, 
      filterFunction = OptionValue[SelectionFunction], 
      filterValue = cutoff, list = OptionValue[ListIndex], 
      component = OptionValue[ComponentIndex]},
        outputData = 
         AssociationThread[
            Keys[#] -> 
             Query[All /* Values /* KeyUnion, 
               Select[If[ MissingQ[#],
                          Identity,
                          filterFunction[#[[
                            Sequence @@ DeleteMissing[{list, component}]]], 
                           filterValue]
                      ] &]]@#] &@data;
        Return[outputData/._Missing->Missing[]]
    ];

(* ::Function:: *)
(* f:FilterMissing *)
(***Options***)
Options[FilterMissing] = {MinimumPoints -> 3, Reference -> {}, 
   ShowPlots -> True};
(*Reference for dropping missing in the reference*)
(***Function***)
FilterMissing[omicsObject_, percentage_, OptionsPattern[]] :=
    Module[ {data = omicsObject, pointMinima = OptionValue[MinimumPoints], 
      percent = percentage, outputData, outputData2, counts, 
      histogramToggle = OptionValue[ShowPlots], testSelection, 
      referencing = Flatten[{OptionValue[Reference]}]},
        counts = Merge[Values[data], Count[Flatten[#], _Missing] &];
        If[ histogramToggle,
            Print[#] & /@ {Histogram[counts, 
               PlotLabel -> "Number of Missing Data Points per Component", 
               AxesLabel -> {"Number of Missing Points", "Counts"}, 
               ImageSize -> Medium], {"Missing -> Counts: ", KeySort@Counts[counts]}, 
              PieChart[KeySort@Counts[counts], ChartLabels -> Values[KeySort@Counts[counts]], 
               ChartLegends -> Keys[KeySort@Counts[counts]], 
               PlotLabel -> "Pie Chart of number of missing components", 
               ImageSize -> Medium]}
        ];
        outputData = 
         Query[All, 
           KeyDrop[Keys[
               Query[Select[(# > 
                     Min[(1 - percent)*Length[data], 
                      Length[data] - pointMinima]) &]]@counts]][#] &]@data;
        If[ MatchQ[referencing, {}],
            Return[outputData],
            testSelection = 
             Query[referencing /* Values /* (Flatten[#, 1] &) /* 
                Union /* (Map[Key, #, 1] &), Select[MissingQ[#] &] /* Keys]@
              outputData;
            outputData2 = Query[All, KeyDrop[testSelection][#] &]@outputData
        ];
        Return[outputData2]
    ]

(* ::Function:: *)
(* f:BoxCoxTransform *)
(***Function***)
BoxCoxTransform = 
 Compile[{{x, _Real}, {lamda, _Real}}, 
  Piecewise[{{((x)^lamda - 1)/lamda, Chop[Abs[lamda]] > 0.}}, Log[x]],
   "RuntimeOptions" -> {"EvaluateSymbolically" -> False}
  (*this is to avoid errors within modules where local variables are \
defined*)];

(* ::Function:: *)
(* f:NormalLogLikelihoodBoxCoxBase *)
NormalLogLikelihoodBoxCoxBase = 
  Compile[{{data, _Real, 1}, {\[Lambda], _Real}},
   (*compile as a function of the data as an array, 
   and lambda as the box cox variable, 
   further defining the function for the log likelihood*)
   N[-(N[Length[data]]/
         2)*(Log[(Variance[
          BoxCoxTransform[N[#], \[Lambda]] & /@ 
           data])]) + ((\[Lambda] - 1)*Plus @@ Log[N[data]]) ], 
   "RuntimeOptions" -> {"EvaluateSymbolically" -> 
      False}(*this is to avoid errors within modules where local \
variables are defined*)];

(* ::Function:: *)
(* f:MaximumLikelihoodBoxCoxTransformBase *) 
MaximumLikelihoodBoxCoxTransformBase[data_] :=
    Module[ {inputData = data, lamdaHat, tempPrint, x, outputData},
     (*maximize the log likelihood, 
     extract exponent and map BoxCoxTransform to Real Number data*)
        tempPrint = 
         PrintTemporary[
          "Computing \!\(\*OverscriptBox[\(\[Lambda]\), \(^\)]\) "];
        lamdaHat = 
         x /. (NMaximize[
               NormalLogLikelihoodBoxCoxBase[Cases[N[#], _Real], x], 
               x][[-1]] &@inputData);
        outputData = 
         Map[If[ MatchQ[#1, _Real],
                 BoxCoxTransform,
                 #1 &
             ][N[#], 
            lamdaHat] &, inputData];
        NotebookDelete[tempPrint];
        Print["Calculated Box-Cox parameter \!\(\*OverscriptBox[\(\[Lambda]\
\), \(^\)]\) = ", lamdaHat];
        Return[outputData]
    ];

(* ::Function:: *)
(* f:ApplyBoxCoxTransform *) 
Options[ApplyBoxCoxTransform] = {ComponentIndex -> Missing[],
	HorizontalSelection -> False,
	ListIndex -> Missing[]}; 
ApplyBoxCoxTransform[data_, OptionsPattern[]] :=
    Module[ {inputData = N[data], list = OptionValue[ListIndex], 
      component = OptionValue[ComponentIndex], 
      horizontal = OptionValue[HorizontalSelection], outputData}(*internal variables*),
        outputData = 
         Returner[inputData, 
          Applier[MaximumLikelihoodBoxCoxTransformBase, inputData, 
           ListIndex -> list, ComponentIndex -> component, 
           HorizontalSelection -> horizontal], ListIndex -> list, 
          ComponentIndex -> component, HorizontalSelection -> horizontal];
        Return[outputData]
    ];
  
(* ::Function:: *)
(* f:StandardizeExtended *) 
StandardizeExtended[inputList_, subtract_, divide_] :=
    Module[ {input = inputList, output, sub, div, subtractFn = subtract, 
      divideFn = divide},
        sub = subtractFn[Cases[N[input], _Real]];
        div = divideFn[Cases[N[input], _Real]];
        output = Map[If[ MatchQ[N[#], _Real],
                         (# - sub)/div,
                         #
                     ] &, input];
        Return[output]
    ]

(* ::Function:: *)
(* f:QuantileNormalizationBase *) 
(***Options***)
Options[QuantileNormalizationBase] = {Averaging -> Mean, Ties -> Mean};
(***Function***)
QuantileNormalizationBase[data_, OptionsPattern[]] :=
    Module[ {inputData = N[data], means, replacementSet, outputData, 
      averagingFunction = OptionValue[Averaging], 
      tiesFunction = OptionValue[Ties]},
     (*create a list of the ordering in each row*)
     (*merge across \
   associations by averaging same ranks*)
        means = Association@ 
          MapIndexed[#2[[1]] -> averagingFunction[#] &, 
           Transpose[Sort[#] & /@ inputData]];
        If[ MatchQ[ToString[tiesFunction],"None"],
         (*if ties are to be assigned randomly*)
            outputData = 
             Map[means[#] &, 
              Keys[#] & /@ (AssociationThread[Ordering@Ordering[#] -> #] & /@ 
                 inputData), {2}],
            (*else create averaging of ties within each list*)
            If[ AllTrue[inputData, DuplicateFreeQ],
                outputData = 
                 Map[means[#] &, 
                  Keys[#] & /@ (AssociationThread[Ordering@Ordering[#] -> #] & /@ 
                     inputData), {2}],
                replacementSet = 
                 Merge[
                    Association[#[[2]] -> means[#[[1]]] ] & /@ 
                     Transpose[{Ordering@Ordering[#], #}], tiesFunction] & /@ 
                  inputData;
                (*for each list the ties can be different, 
                so each replacement is threaded*)
                outputData = MapThread[Lookup, {replacementSet, inputData}]
            ]
        ];
        Return[outputData]
    ]

(* ::Function:: *)
(* f:QuantileNormalization *)
(***Options***)
Options[QuantileNormalization] = {Averaging -> Mean,
	ComponentIndex -> Missing[],
	ListIndex -> Missing[],
	Ties -> Mean};
(***Function***)   
QuantileNormalization[data_, OptionsPattern[]] :=
    Module[ {inputData = N[data], outputData, 
      averagingFunction = OptionValue[Averaging], 
      tiesFunction = OptionValue[Ties], list = OptionValue[ListIndex], 
      component = OptionValue[ComponentIndex]},
        outputData = 
         Returner[#, 
            ApplierList[
             QuantileNormalizationBase[#, Averaging -> averagingFunction, 
               Ties -> tiesFunction] &, 
             Applier[# &, #, ListIndex -> list, 
              ComponentIndex -> component]], ListIndex -> list, 
            ComponentIndex -> component] &@inputData;
        Return[outputData]
    ]

(* ::Function:: *)
(* f:Applier *)
(***Options***)
Options[Applier] = {ComponentIndex -> Missing[],
	HorizontalSelection -> False,
	ListIndex -> Missing[]};
(***Function***)  
Applier[function_, inputData_, OptionsPattern[]] :=
    Module[ {functional = function, data = N[inputData], 
      list = OptionValue[ListIndex], 
      component = OptionValue[ComponentIndex], 
      horizontal = OptionValue[HorizontalSelection], outputData, 
      dataTest},
        dataTest = FirstCase[data, Except[_Missing]];
        If[ AssociationQ[data],
            If[ ! (ListQ[dataTest] || 
                AssociationQ[
                 dataTest])(*check if this is an association of numbers*),
                outputData = 
                 AssociationThread[Keys[data], functional[Values[data]]],
                If[ horizontal == True,
                    outputData = 
                     AssociationThread[Keys[data], 
                      functional[
                       Query[Sequence @@ 
                          DeleteMissing[{All /* Values, list, component}]]@data]],
                    outputData = 
                     If[ ListQ[dataTest](*check if this is a list*),
                         Query[Sequence @@ 
                            DeleteMissing[{All /* (Map[functional, #] &), All, list, 
                              component}]]@data,
                         If[ AssociationQ[dataTest](*check if this is an Association*),
                             AssociationThread[
                                Keys[#] -> 
                                 Query[Sequence @@ 
                                    DeleteMissing[{All /* Values /* functional, list, 
                                      component}]]@#] & /@ data,
                             {"Not Association", 
                             Map[# &, data, 1]}
                         ]
                     ]
                ]
            ],
            outputData = If[ ListQ[dataTest],
                             Map[functional, data, {1}],
                             functional@data
                         ];
        ];
        Return[outputData]
    ];
     
(* ::Function:: *)
(* f:ApplierList ::*)
(***Function***)
ApplierList[function_, inputData_] :=
    Module[ {data = N[inputData], outputData, f = function, dataTest, 
      dataTestFirst},
        dataTest = FirstCase[data, Except[_Missing]];
        If[ AssociationQ[data](*Check if this an association*),
            If[ AssociationQ[
              dataTest](*check if first component is an association*),
                outputData = 
                 AssociationThread[#1, #2] & @@ {Keys[#], (MapThread[
                       AssociationThread[#1 -> #2] &, {Keys[#], f@Values[#]} &@
                        Values[#]])} &@data,
                dataTestFirst = FirstCase[dataTest, Except[_Missing]];
                If[ ListQ[dataTestFirst],
                    outputData = 
                     AssociationThread[Keys[#] -> { f @@ Values[#]}] &@data,
                    outputData = AssociationThread[Keys[#] -> f@Values[#]] &@data
                ]
            ],
            If[ AssociationQ[dataTest],
                outputData = AssociationThread[#[[1]] -> #[[2]]] & /@Transpose[{Keys@data,f[Values[data]]}],
                outputData = f[data]
            ]
        ];
        Return[outputData]
    ];
  
(* ::Function:: *)
(* f:Returner *)
(***Options***)
Options[Returner] = {ComponentIndex -> Missing[],
	HorizontalSelection -> False,
	ListIndex -> Missing[]};
(***Function***)   
Returner[originalAssociation_, update_, OptionsPattern[]] :=
    Module[ {data = originalAssociation, modifiedData = update, 
      list = OptionValue[ListIndex], 
      component = OptionValue[ComponentIndex], 
      horizontal = OptionValue[HorizontalSelection], outputData, 
      dataTest, dataTestFirst},
        dataTest = FirstCase[data, Except[_Missing | <||>]];
        If[ AssociationQ[data],
            If[ ! (ListQ[dataTest] || 
                AssociationQ[
                 dataTest])(*check if this is an association of numbers*),
                outputData = modifiedData,
                If[ horizontal == True,
                    outputData = Merge[
                       ReplacePart[{Sequence @@ 
                             DeleteMissing[{list, component}]} -> #[[2]]][#[[
                          1]]] &][{data, modifiedData}],
                    dataTestFirst = FirstCase[dataTest, Except[_Missing]];
                    outputData = If[ ListQ[dataTest],
                                     If[ ListQ[dataTestFirst],
                                         Merge[
                                           
                                           MapThread[
                                             ReplacePart[{Sequence @@ 
                                                   DeleteMissing[{list, 
                                                     component}]} -> #2][#1] &, #] &][{data, 
                                           modifiedData}],
                                         modifiedData
                                     ],
                                     If[ ! (ListQ[dataTestFirst] || AssociationQ[dataTestFirst]),
                                         outputData = modifiedData,
                                         MapThread[Merge[{#1, #2},
                                            
                                            ReplacePart[{Sequence @@ 
                                                  DeleteMissing[{list, component}]} -> #[[2]]][#[[
                                               1]]] &] &, {data, modifiedData}]
                                     ]
                                 ]
                ]
            ],
            outputData = update
        ];
        Return[outputData]
    ];
  
(* ::Function:: *)
(* f:ConstantAssociation *)  
(*Adds multi key constant to association, each addition specified in \
a single associaiton of form <|addition1\[Rule] Value1,addition2 \
\[Rule] Value2,...|>*)
(***Function***)
ConstantAssociator[inputAssociation_, constant_] :=
    Module[ {dataIn = inputAssociation, mod = constant, keysIn, dataOut},
     (*get association In keys values*)
        keysIn = Keys@dataIn[[1]];
        (*extract keys, create an array to add to each, 
        add resulting association to original*)
        dataOut = 
         AssociateTo[dataIn, 
          AssociationThread@(Keys[mod] -> ( 
              AssociationThread[
                 keysIn -> ConstantArray[#, Length[keysIn]]] & /@ 
               Values[mod]))];
        Return[dataOut]
    ]

(* ::Function:: *)
(* f:EnlargeInnerAssociation *)
(*enlarges associations of multiple inputs - care must be taken to \
have different identities for each association *)
(*Inside Keys must \
be different*)
(***Options***)
Options[EnlargeInnerAssociation] = {KeyModifiers -> None,
	MergeFunction -> (Join[Sequence @@ #] &),
	MissingReplacement -> Missing[]};
(***Function***) 
EnlargeInnerAssociation[omicsObjectList_, OptionsPattern[]] :=
    Module[ {dataIn = omicsObjectList, missingSub = OptionValue[MissingReplacement], 
      mergeF = OptionValue[MergeFunction], 
      keyMods = OptionValue[KeyModifiers], dataOut},
     (*If modifiers included key map the modifiers first*)
        If[ !MatchQ[keyMods, None],
            dataIn = 
             MapThread[ 
              Query[All, 
                 KeyMap[Append[#2]]]@(Query[All, 
                   KeyMap[If[ ListQ[#],
                              #,
                              {#}
                          ] &]][#1]) &, {dataIn, 
               keyMods}]
        ];
        (*Create Missing[] association elements if the associations are \
      unequal, and then Merge by joining the associations of each sample*)
        dataOut = 
        Merge[MapThread[
          If[ MatchQ[#2, <||>],
              #1,
              ConstantAssociator[#1, #2 /. _Missing -> 
                 missingSub]
          ] &, {dataIn, 
           Query[All, Select[MatchQ[#, _Missing] &]]@KeyUnion[dataIn]}], 
         mergeF];
        Return[dataOut]
    ];
  
  
(* ::Function:: *)
(* f:EnlargeOuterAssociation *)
(*add new data to same type of association*)
(*outside Keys must be \
different*)
(***Options***)
Options[EnlargeOuterAssociation] = {MissingReplacement-> Missing[]};
(***Function***)
EnlargeOuterAssociation[omicsObjectList_,OptionsPattern[]] :=
    Module[ {dataIn = omicsObjectList, missingSub = OptionValue[MissingReplacement],outerKeys, values, dataOut},
     (*If modifiers included key map the modifiers first*)
        outerKeys = Join[Sequence @@ Keys[#] & /@ dataIn];
        values = KeyUnion@Join[Sequence @@ Values[#] & /@ dataIn];
        dataOut = AssociationThread[outerKeys, values/. _Missing-> missingSub];
        Return[dataOut]
    ];
  
(* ::Function:: *)
(* f:LowValueTag *)
(*function takes in omics association and tags values in list1, index \
1 or otherwise as Missing[] based on cutoff provided*)
(***Options***)
Options[LowValueTag] = {ComponentIndex -> 1, 
	ListIndex -> 1,
	OtherReplacement -> _Missing :> Missing[],
	ValueReplacement -> Missing[]};
(***Function***) 
LowValueTag[omicsObject_, valueCutoff_, OptionsPattern[]] :=
    Module[ {data = omicsObject, cutoff = valueCutoff, 
      listIn = OptionValue[ListIndex], 
      comIn = OptionValue[ComponentIndex], 
      replacer = OptionValue[ValueReplacement], 
      otherReplacer = OptionValue[OtherReplacement], returning},
        returning = 
         Returner[data, 
          Applier[(# /. {x_ /; x <= cutoff :> replacer, otherReplacer}) &, data, ListIndex -> listIn, 
           ComponentIndex -> comIn], ListIndex -> listIn, 
          ComponentIndex -> comIn];
        Return[returning]
    ];

(* ::Function:: *)
(* f:MeasurementApplier *)
(*function takes in omics object and applies a function to its measurement list, ignoring Missing values*)
(***Options***)
Options[MeasurementApplier] := {ComponentIndex -> All,
  IgnorePattern -> _Missing,
  ListIndex -> 1};

(***function***)
MeasurementApplier[function_, omicsObject_, OptionsPattern[]] := 
  Module[{fn = function, data = omicsObject, 
    deletePattern = OptionValue[IgnorePattern], 
    list = OptionValue[ListIndex], 
    component = OptionValue[ComponentIndex], returning},
   returning = 
    Returner[data, 
     Applier[Map[
        If[MissingQ[#], Missing[], 
          If[MatchQ[#, {}], Missing[], {fn[#]}] &@
           DeleteCases[#, deletePattern]] &, #] &, data, 
      ListIndex -> list, ComponentIndex -> component], ListIndex -> 1];
   Return[returning]];


(* ::Function:: *)
(* f:CreateTimeSeries *)
(***Function***)  
CreateTimeSeries[
   dataIn_] :=
    (Query[All /* Values /* Merge[Flatten[#] &], All, 
        1]@(KeySortBy[#, ToExpression[#] &] &@
         dataIn)) /. {"KeyAbsent" -> Missing[], _Missing -> Missing[]};

(* ::Function:: *)
(* f:TimeExtractor *)
(***Function***)
TimeExtractor[dataIn_] :=
    Sort[ToExpression@Keys@dataIn];

(* ::Function:: *)
(* f:SeriesApplier *)
(***Options***)
Options[SeriesApplier] = {MissingMask -> 0};
(***Function***) 
SeriesApplier[function_,dataIn_, OptionsPattern[]] :=
    Module[ {data = dataIn, functionF = function, 
      missingMask = OptionValue[MissingMask], returning},
        returning = Query[All, (Apply[If[ NumericQ[#1],
                                          #2,
                                          #1
                                      ] &, 
             Transpose@{#,  
               functionF[# /. Missing[] -> missingMask]}, {1}]) &]@(data /. 
           "KeyAbsent" -> Missing[]);
        Return[returning]
    ];
 
(* ::Function:: *)
(* f:SeriesCompare *)
(*merge two time series by pointwise binary operation, default is \
subtraction*)
(***Options***)
Options[SeriesCompare] = {CompareFunction -> (If[ MatchQ[Head[#1], Missing] || MatchQ[Head[#2], Missing],
                                                  Missing[],
                                                  (#1 - #2)
                                              ] &)};
(***Function***)
SeriesCompare[series1_, series2_, OptionsPattern[]] :=
    Module[ {seriesMain = series1, seriesControl = series2, returning, 
      binaryF = OptionValue[CompareFunction]},
        returning = 
         Merge[KeyIntersection[{seriesMain, seriesControl}], 
          MapThread[binaryF, #] &];
        Return[returning]
    ];

(* ::Function:: *)
(* f:SeriesInternalCompare *)  
(*compare all points to an internal reference, If refernce is absent \
Delete corresponding series*)
(***Options***)
Options[SeriesInternalCompare] = {CompareFunction -> (If[ MatchQ[Head[#1], Missing] || MatchQ[Head[#2], Missing],
                          Missing[],
                          Table[If[ MatchQ[Head[#1[[i]]], Missing],
                                    Missing[],
                                    #1[[i]] - #2
                                ], {i, Length[#1]}]
                      ](*/._Plus\[Rule] 
     Missing[]*)&),
     ComparisonIndex -> 1,
     DeleteRule -> {Head, Missing}}; 
(***Function***)
SeriesInternalCompare[dataIn_, OptionsPattern[]] :=
    Module[ {data = dataIn, comIn = OptionValue[ComparisonIndex], 
      delRule = OptionValue[DeleteRule], 
      binaryF = OptionValue[CompareFunction], filterMissing, 
      returning},
        filterMissing = 
         DeleteCases[data, 
          x_ /; MatchQ[delRule[[1]][x[[comIn]]], delRule[[2]]]];
        returning = Query[All, binaryF[#, #[[comIn]]] &]@filterMissing;
        Return[returning]
    ]

(* ::Function:: *)
(* f:ConstantSeriesClean *)
(***Options***)
Options[ConstantSeriesClean] = {InverseSelection -> False,ReturnDropped -> False};
(***Function***)   
ConstantSeriesClean[inputData_, OptionsPattern[]] :=
    Module[ {dataIn = inputData, data, constantKeys, 
      returnDropped = OptionValue[ReturnDropped], 
      inverseSelection = OptionValue[InverseSelection], returning},
        data = If[ AssociationQ[dataIn],(*already an association*)
                   dataIn,
    (*create an association*)
                   AssociationThread@( Range[Length[#]] -> #) &@(If[ Length[dataIn[[1]]] == 0,
                                                                     List,
                                                                     Identity
                                                                 ]@dataIn)
               ];
        constantKeys = (Query[
            All /* (# &) /* (Select[# == inverseSelection &] /* 
               Keys), (Tally[DeleteMissing[#]] &) /* ( Length[#] > 1 &)]@
           data);
        returning = 
        If[ returnDropped == False,
            If[ !MatchQ[constantKeys, {}],
                Print["Removed series and returning filtered list. If you would like a list of removed keys run \
the command ConstantSeriesClean[data,ReturnDropped \[Rule] \
True]."]
            ];
            KeyDrop[data, constantKeys],
            <|
            "Data" -> KeyDrop[data, constantKeys], 
            "Dropped Keys" -> constantKeys|>
        ];
        Return[If[ AssociationQ[dataIn],
                   returning,
                   If[ Length[dataIn[[1]]] == 0,
                       First@Values[returning],
                       Values[returning]
                   ]
               ]]
    ];

JoinNestedAssociations[data_] :=
    Module[ {dataIn = data,returning},
        returning = Merge[dataIn,Join@@#&];
        Return[returning]
    ];


(* :: Section:: *)
(*#####HypothesisTesting#####*)

(* ::Function:: *)
(*BenjaminiHochbergFDR*)
(***Options***)
Options[BenjaminiHochbergFDR] = {SignificanceLevel -> 0.05}; 
(***Function***)
BenjaminiHochbergFDR[pValues_, OptionsPattern[]] :=
    Module[ {pVals = N[pValues], 
      sigLevel = N[OptionValue[SignificanceLevel]], nTests, sortedpVals, 
      sortingIDs, weightedpVals, adjustedpVals, qVals, cutoffqValue, 
      pValqValAssociation, pValCutoff, 
      returning},(*count number of hypotheses tested*)
        nTests = Length[pVals];
        (*sort the pValues in order*)
        sortedpVals = Sort[pVals];
  (*generate a sorting ID by ordering*)
        sortingIDs = Ordering[Ordering[pVals]];
        (*adjust p values to weighted p values*)
        weightedpVals = nTests*sortedpVals/Range[nTests];
        (*flatten the weighted p-
        values to smooth out any local maxima and get adjusted p-vals*)
        adjustedpVals = Table[Min[weightedpVals[[i ;;]]], {i, 1, nTests}];
        (*finally,generate the qVals by reordering*)
        qVals = adjustedpVals[[sortingIDs]];
        (*create an association from which to identify correspondence \
      between p-values and q-values*)(*Print[{qVals,pVals}]*);
        pValqValAssociation = AssociationThread[qVals -> pVals];
        (*get the cutoff adjusted q-value*)
        cutoffqValue = 
         SelectFirst[Reverse[adjustedpVals], # <= sigLevel &];
        (*Print[cutoffqValue];*)(*then and identify corresponding cutoff p-
        value*)
        pValCutoff = 
        If[ ToString[Head[cutoffqValue]] == ToString[Missing],
            0;
            cutoffqValue = 0,
            pValqValAssociation[cutoffqValue]
        ];
  (*return it all*)
        returning = 
         Association@{"Results" -> 
            Transpose@{pVals, qVals, # <= cutoffqValue & /@ qVals}, 
           "p-Value Cutoff" -> pValCutoff, "q-Value Cutoff" -> cutoffqValue};
        (*get q-vals and q-
        val cutoff*)
        (*now test the qVals for being above or below \
      significance level-
        return "true" if enriched and "false" if not*)
        Return[returning]
    ];


(* ::Section:: *)
(*#####MassSpec#####*)


(* ::Subsection:: *)
(* MzML *)

(* ::Subsubsection:: *)
(* MetaData *)

(* ::Function:: *)
(* f:MzmlRawMetaData:: *)
(***Function***) 
MzmlRawMetaData[mzMLFile_] :=
    Module[ {inputFile = mzMLFile, emptyList, filePath, allMetaData, testXML},
        emptyList = {};
        filePath = inputFile;
        (*OpenRead is used to open the input file and form a <Stream>*)
        inputFile = OpenRead[filePath];
        (*SetStreamPosition will set the data streaming position to zero \
      and this way it will stream from the begining of the file. 
        While loop and Read are used to keep reading the mzML file.*)
        SetStreamPosition[filePath, 0];
        While[True,
         allMetaData = Read[filePath, String];
         (*StringCases to search for "scan" tag and <..> 
         is to keep keep representing "scan" until it matches and if no \
      matching, use AppendTo to append into the emptylist. 
         If statement and Break to stop once "scan" matchs.*)
         If[ StringCases[allMetaData, WordCharacter ..][[1]] == "spectrum",
             Break[],
             AppendTo[emptyList, allMetaData]
         ]];
        (*These </spectrumList></run></mzML></indexedmzML> 
        are the closing tags for the streamed elements that we need to \
      convert to the mathematica "XML" format.*)
        emptyList = 
         AppendTo[emptyList, "</spectrumList></run></mzML></indexedmzML>"];
        (*Import as "XML" for data conversion to XML format*)
        testXML = ImportString[StringJoin[emptyList], "XML"];
        Close[inputFile];
        Return[testXML]
    ]

(* ::Function:: *)
(* f:MzmlVersion:: *)
(***Function***) 
MzmlVersion[mzMLSchemaVersion_] :=
    Module[ {
      mzMLSchema = mzMLSchemaVersion,
      mzMLSchemaTag, allTogether},
     
     (*First convert to string, 
     then we use StringCases and string pattern matching to look for \
   "mzML_" followed by a number <NumberString> 
     to find the mzXML version format. 
     "mzML_" will always be followed by a number that corresponds to the \
   mzML schema. This will be output in the <MetaData>.*)
        mzMLSchemaTag = Grid[StringCases[ToString[mzMLSchema],
           ___ ~~ "mzML" ~~ x : NumberString ~~ ___ :> {"mzML:", x}]];
        
        (*allTogether is to output a readable mzML schema version *)
        allTogether = 
         Insert[Insert[
           Grid[Transpose[
             Transpose[{{Style["mzML Schema:", Bold, 
                 Underlined]}, {mzMLSchemaTag}}]]], Alignment -> Left, 2], 
          Spacings -> {2, {2, {0.7}, 0.9}}, 2];
        Return[allTogether]
    ]

(* ::Function:: *)
(* f:MzmlCvList:: *)
(***Function***) 
MzmlCvList[cvListXMLElement_] :=
    Module[ {
      cvListElement,
      (*cvList is equal to the only input and that is MzmlRawMetaData \
    output*)
      cvList = cvListXMLElement,
      cvListTag1, cvListTag2, cvListTag3, allTogether},
     
     (*extract Element <cv> 
     which contain the information about file ontology including 3 \
   required attributes, URI, FullName and id.*)
        cvListElement = Cases[cvList, XMLElement["cv", _, _], Infinity];
        
        (*extracting all <id>'s, <version>'s, <fullName> and <URI> 
        that are inside cvListElement*)
        cvListTag1 =
         {Cases[#, x_ /; Length[x] > 1 && x[[1]] == "id" :> x[[2]], All],
            (*since the attribute <version> is optional, 
            so we have an if statment here*)
            
            If[ Cases[#, 
               x_ /; Length[x] > 1 && x[[1]] == "version" :> x[[2]], 
               All] == {},
                {"Unknown"},
                Cases[#, x_ /; Length[x] > 1 && x[[1]] == "version" :> x[[2]], 
                 All]
            ], Cases[#, 
             x_ /; Length[x] > 1 && x[[1]] == "fullName" :> x[[2]], All], 
            Cases[#, x_ /; Length[x] > 1 && x[[1]] == "URI" :> x[[2]], 
             All]} & /@ cvListElement;
        
        (*capitalized ID: and Full Name: 
        labeling for version and URI respectively*)
        cvListTag2 =
         {StringInsert[#[[1]], "ID: ", 1], 
            StringInsert[#[[2]], "Version: ", 1], 
            StringInsert[#[[3]], "Full Name: ", 1], 
            StringInsert[#[[4]], "URI: ", 1]} & /@ cvListTag1;
        cvListTag3 = 
         Insert[Grid[
           Partition[Grid[#, Alignment -> Left] & /@ cvListTag2, 1], 
           Dividers -> {Center, {True, True, True}}], Alignment -> Left, 
          2];
        
        (*allTogether is to output a readable mzML schema version *)
        allTogether =
         (*incase Element <cv> is absent from the file*)
         
         If[ Cases[cvList, XMLElement["cv", _, _], Infinity] == {},
             Insert[Insert[
               Grid[Transpose[
                 Transpose[{{Style["Ontology Information:", Bold, 
                     Underlined]}, {"Not Found"}}]]], Alignment -> Left, 2], 
              Spacings -> {2, {2, {0.7}, 0.9}}, 2],
             Insert[Insert[
               Grid[Transpose[
                 Transpose[{{Style["Ontology Information:", Bold, 
                     Underlined]}, {cvListTag3}}]]], Alignment -> Left, 2], 
              Spacings -> {2, {2, {0.7}, 0.9}}, 2]
         ];
        Return[allTogether]
    ]

(* ::Function:: *)
(* f:MzmlFileDescriptionFileContent:: *)
(***Function***) 
MzmlFileDescriptionFileContent[
  fileDescriptionfileContentXMLElement1_] :=
    Module[ {fileDescriptionElement,
      fileDescription = fileDescriptionfileContentXMLElement1,
      fileContentTag1, fileContentTag2, fileContentTag3, fileContentTag4, allTogether},
     
     (*fileDescriptionElement is equal to the only input and that is the \
   MzmlRawMetaData output*)
        fileDescriptionElement =
         (*extract element <fileDescription> 
         where information pertaining to the entire mzML file*)
         
         Cases[fileDescription,
          XMLElement["fileDescription", _, _], Infinity];
        
        (*use Cases to find the extract subelement "fileContent" located \
      inside element <fileDescription> 
        which summarizes the different types of spectra that can be \
      expected in the file. *)
        fileContentTag1 = Cases[fileDescriptionElement,
          XMLElement["fileContent", _, _], Infinity];
        
        (*cvParam holds additional data or annotation. 
        Only controlled values are allowed here, 
        while userParam holds uncontrolled user parameters*)
        fileContentTag2 =
         (*We use an if statment, 
         since subelements "cvParam" and "userParam" are not requered and \
        there could be many of them*)
         
         If[ Cases[fileContentTag1, 
            XMLElement["cvParam" | "userParam", _, _], 
            Infinity] == {},
             {"Unknown"},
             Cases[fileContentTag1, XMLElement["cvParam" | "userParam", _, _], 
              Infinity]
         ];
        
        (* The attribute "name" is the actual required name for the \
      parameter,from the referred-to controlled vocabulary. 
        The attribute "name" is the actual name for the parameter,
        from the referred-to controlled vocabulary. 
        The attribute "value" is the optional value for the parameter;
        may be absent if not appropriate. And if absent we output "Unknown"*)
        fileContentTag3 = {Cases[#, 
             x_ /; Length[x] > 1 && x[[1]] == "name" :> x[[2]], All], 
            If[ Cases[#, x_ /; Length[x] > 1 && x[[1]] == "value" :> x[[2]], 
               All] == {""},
                {"Unknown"},
                Cases[#, x_ /; Length[x] > 1 && x[[1]] == "value" :> x[[2]], 
                 All]
            ]} & /@ fileContentTag2;
        
        (*capitalized Name: and Value: labeling*)
        fileContentTag4 = 
         Insert[Grid[
           Flatten[{StringInsert[#[[1]], "Name: ", 1], 
               StringInsert[#[[2]], " Value: ", 1]}, 1] & /@ 
            fileContentTag3], Alignment -> Left, 2];
        
        (*allTogether is to output a readable mzML schema version *)
        allTogether =
         (*if statment in case the subelements "cvParam" and "userParam" \
        are absent*)
         
         If[ Cases[fileContentTag1, 
            XMLElement["cvParam" | "userParam", _, _], Infinity] == {},
             Insert[Insert[
               Grid[Transpose[
                 Transpose[{{Style["File Content:", Bold, 
                     Underlined]}, {"Not Found"}}]]], Alignment -> Left, 2], 
              Spacings -> {2, {2, {0.7}, 0.9}}, 2],
             Insert[Insert[
               Grid[Transpose[
                 Transpose[{{Style["File Content:", Bold, 
                     Underlined]}, {fileContentTag4}}]]], Alignment -> Left, 
               2], Spacings -> {2, {2, {0.7}, 0.9}}, 2]
         ];
        Return[allTogether]
    ];
  
(* ::Function:: *)
(* f:MzmlFileDescriptionSourceFileList:: *)
(***Function***) 
MzmlFileDescriptionSourceFileList[
  FileDescriptionSourceFileListXMLElement1_] :=
    Module[ {fileDescriptionElement,
      fileDescription = FileDescriptionSourceFileListXMLElement1,
      sourceFileListTag1, sourceFileListTag2, sourceFileListTag3, 
      sourceFileListTag4, sourceFileListTag5, sourceFileListTag6, 
      sourceFileListTag7, sourceFileListTag8},
     
     (*fileDescriptionElement is equal to the only input and that is the \
   MzmlRawMetaData output*)
        fileDescriptionElement =
         (*extract element <fileDescription> 
         where information pertaining to the entire mzML file*)
         
         Cases[fileDescription,
          XMLElement["fileDescription", _, _], Infinity];
        
        (*
        sourceFileListTag1 is equal to the soureFile subelement and \
      that is where list and descriptions of the source files this mzML \
      document was generated or derived from*)
        sourceFileListTag1 =
         (*we use if statement in case the subelement "sourceFile" is \
        absent*) If[ Cases[fileDescriptionElement,
            XMLElement["sourceFile", _, _], Infinity] == {},
                     {"Unknown"},
                     Cases[fileDescriptionElement,
                      XMLElement["sourceFile", _, _], Infinity]
                 ];
        
        
        (* The attribute "id" is the required identifier for the current \
      file. The attribute "location" is the required URI-
        formatted location where the file was retrieved.*)
        sourceFileListTag2 = {Cases[#, 
             x_ /; Length[x] > 1 && x[[1]] == "id" :> x[[2]], All], 
            Cases[#, x_ /; Length[x] > 1 && x[[1]] == "location" :> x[[2]], 
             All]} & /@ sourceFileListTag1;
        
        (*capitalized ID: and Location: labeling*)
        sourceFileListTag3 = 
         Insert[Grid[
           Flatten[{StringInsert[#[[1]], "ID: ", 1], 
               StringInsert[#[[2]], "Location: ", 1]}, 1] & /@ 
            sourceFileListTag2, Dividers -> Center], Alignment -> Left, 2];
        
        
        (*cvParam holds additional data or annotation. 
        Only controlled values are allowed here, 
        while userParam holds uncontrolled user parameters*)
        sourceFileListTag4 =
         (*We use an if statment, 
         since subelements "cvParam" and "userParam" are not requered and \
        there could be many of them*)
         If[ Cases[sourceFileListTag1, 
            XMLElement["cvParam" | "userParam", _, _], 
            Infinity] == {},
             {"Unknown"},
             Cases[sourceFileListTag1, 
              XMLElement["cvParam" | "userParam", _, _], Infinity]
         ];
        
        
        (* The attribute "name" is the actual required name for the \
      parameter,from the referred-to controlled vocabulary. 
        The attribute "name" is the actual name for the parameter,
        from the referred-to controlled vocabulary. 
        The attribute "value" is the optional value for the parameter;
        may be absent if not appropriate. And if absent we output "Unknown"*)
        sourceFileListTag5 = {Cases[#, 
           x_ /; Length[x] > 1 && x[[1]] == "name" :> x[[2]], All], 
          If[ Cases[#, x_ /; Length[x] > 1 && x[[1]] == "value" :> x[[2]], 
             All] == {""},
              {"Unknown"},
              Cases[#, x_ /; Length[x] > 1 && x[[1]] == "value" :> x[[2]], 
               All]
          ]} & /@ sourceFileListTag4;
      
      
      (*capitalized Name: and Value: labeling*)
        sourceFileListTag6 = 
         Insert[Grid[
           Flatten[{StringInsert[#[[1]], "Name: ", 1], 
               StringInsert[#[[2]], "Value: ", 1]}, 1] & /@ 
            sourceFileListTag5, Dividers -> Center], Alignment -> Left, 2];
        sourceFileListTag7 = 
         Insert[Grid[
           Transpose[
            Transpose[{{sourceFileListTag3}, {sourceFileListTag6}}]], 
           Dividers -> Center], Alignment -> Left, 2];
        
        (*sourceFileListTag8 is to output a readable mzML schema version *)
        sourceFileListTag8 =
        (*if statment in case the subelements "cvParam" and "userParam" \
        are absent*)
        If[ Cases[sourceFileListTag1, 
          XMLElement["cvParam" | "userParam", _, _], Infinity] == {},
            Insert[Grid[
              Transpose[
               Transpose[{{Style["Source File: ", Bold, 
                   Underlined]}, {"Not Found"}}]]], Alignment -> Left, 2],
            Insert[Grid[
              Transpose[
               Transpose[{{Style["Source File: ", Bold, 
                   Underlined]}, {sourceFileListTag7}}]]], Alignment -> Left, 
             2]
        ];
        Return[sourceFileListTag8]
    ];  

(* ::Function:: *)
(* f:MzmlReferenceParamGroup:: *)
(***Function***) 
MzmlReferenceParamGroup[referenceableParamGroupXMLElement_] :=
    Module[ {referenceableParamGroupElement,
      referenceableParamGroup = referenceableParamGroupXMLElement,
      referenceableParamGroupTag1, referenceableParamGroupTag2, 
      referenceableParamGroupTag3, referenceableParamGroupTag4},
     
     (*Using Cases to extract the file element "referenceableParamGroup"*)
        referenceableParamGroupElement = Cases[referenceableParamGroup,
        XMLElement["referenceableParamGroup", _, _], Infinity];
      
      
      (*cvParam holds additional data or annotation. 
      Only controlled values are allowed here, 
      while userParam holds uncontrolled user parameters*)
        referenceableParamGroupTag1 = 
         Cases[referenceableParamGroupElement, 
          XMLElement["cvParam" | "userParam", _, _], Infinity];
        
        (*The attribute "name" is the actual name for the parameter,
        from the referred-to controlled vocabulary. 
        The attribute "value" is the optional value for the parameter;
        may be absent if not appropriate. And if absent we output "Unknown"*)
        referenceableParamGroupTag2 = {Cases[#, 
           x_ /; Length[x] > 1 && x[[1]] == "name" :> x[[2]], All], 
          If[ Cases[#, x_ /; Length[x] > 1 && x[[1]] == "value" :> x[[2]], 
             All] == {""},
              {"Unknown"},
              Cases[#, x_ /; Length[x] > 1 && x[[1]] == "value" :> x[[2]], 
               All]
          ]} & /@ referenceableParamGroupTag1;
      
      
      (*capitalized Name: and Value: labeling*)
        referenceableParamGroupTag3 = 
         Insert[Grid[
           Flatten[{StringInsert[#[[1]], "Name: ", 1], 
               StringInsert[#[[2]], " Value: ", 1]}, 1] & /@ 
            referenceableParamGroupTag2], Alignment -> Left, 2];
        
        
        (*we use if statement in case the subelement \
      "referenceableParamGroup" is absent*)
        referenceableParamGroupTag4 = If[ Cases[referenceableParamGroup,
            XMLElement["referenceableParamGroup", _, _], Infinity] == {},
                                          Insert[Insert[
                                            Grid[Transpose[
                                              Transpose[{{Style["Instrument Params: ", Bold, 
                                                  Underlined]}, {"Not Found"}}]], 
                                             Dividers -> {{False, False}, {False, True, True}}], 
                                            Alignment -> Left, 2], Spacings -> {2, {2, {0.7}, 0.9}}, 2],
                                          Insert[Insert[
                                            Grid[Transpose[
                                              Transpose[{{Style["Instrument Params: ", Bold, 
                                                  Underlined]}, {referenceableParamGroupTag3}}]], 
                                             Dividers -> {{False, False}, {False, True, True}}], 
                                            Alignment -> Left, 2], Spacings -> {2, {2, {0.7}, 0.9}}, 2]
                                      ];
        Return[referenceableParamGroupTag4]
    ]

(* ::Function:: *)
(* f:MzmlSoftware:: *)
(***Function***)
MzmlSoftware[softwareXMLElement_] :=
    Module[ {softwareElement,
      software = softwareXMLElement,
      softwareTag1, softwareTag2, allTogether},
     
     (*Using Cases to extract the file element "software"*)
        softwareElement = Cases[software,
          XMLElement["software", _, _], Infinity];
        
        (*The attribute "id" is an identifier for the software that is \
      unique across all software types used in this file. 
        The attribute "version" is the software version. 
        While the attribute "name" is the actual name for the parameter,
        from the referred-to controlled vocabulary. 
        The attribute "value" is the optional value for the parameter;
        may be absent if not appropriate. 
        And if the attributes "name" and "value" are absent we output \
      "Unknown"*)
        softwareTag1 = 
         Flatten[{Cases[#, x_ /; Length[x] > 1 && x[[1]] == "id" :> x[[2]], 
              All], Cases[#, 
              x_ /; Length[x] > 1 && x[[1]] == "version" :> x[[2]], All], 
             If[ Cases[#, x_ /; Length[x] > 1 && x[[1]] == "name" :> x[[2]], 
                All] == {""},
                 {"Unknown"},
                 Cases[#, x_ /; Length[x] > 1 && x[[1]] == "name" :> x[[2]], 
                  All]
             ], If[ Cases[#, x_ /; Length[x] > 1 && x[[1]] == "value" :> x[[2]], 
                All] == {""},
                    {"Unknown"},
                    Cases[#, x_ /; Length[x] > 1 && x[[1]] == "value" :> x[[2]], 
                     All]
                ]}, 1] & /@ softwareElement;
        
        
        (*capitalized ID, Name and Value labeling*)
        softwareTag2 = 
         Insert[Grid[
           Transpose[{StringInsert[#[[1]], "ID: ", 1], 
               StringInsert[#[[2]], "Version: ", 1], 
               StringInsert[#[[3]], "Name: ", 1], 
               StringInsert[#[[4]], "Value: ", 1]} & /@ softwareTag1]], 
          Alignment -> Left, 2];
        
        
        (*we use if statement in case the subelement "software" is absent*)
        allTogether = If[ Cases[software,
          XMLElement["software", _, _], Infinity] == {},
                          Insert[Insert[
                            Grid[Transpose[
                              Transpose[{{Style["Instrument Software List: ", Bold, 
                                  Underlined]}, {"Not Found"}}]]], Alignment -> Left, 2], 
                           Spacings -> {2, {2, {0.7}, 0.9}}, 2],
                          Insert[
                           Insert[Grid[
                             Transpose[
                              Transpose[{{Style["Instrument Software List: ", Bold, 
                                  Underlined]}, {softwareTag2}}]]], Alignment -> Left, 2], 
                           Spacings -> {2, {2, {0.7}, 0.9}}, 2]
                      ];
        Return[allTogether]
    ]

(* ::Function:: *)
(* f:MzmlInstrumentConfiguration:: *)
(***Function***)
MzmlInstrumentConfiguration[instrumentConfigurationXMLElement_] :=
    Module[ {
      instrumentConfiguration = instrumentConfigurationXMLElement,
      EmptyFileContent, SourceAllTogether, AnalyzerAllTogether, 
      DetectorAllTogether, overAllInstrumentConfigurationTag},
     
     (*Extracting the Element <instrumentConfigurationList> 
     where there will be a list and descriptions of the instrument \
   configurations. 
     The attribute "id" is required and that is identifier for the \
   instrument configuration. Also,
     At least one instrument configuration MUST be specified,
     even if it is only to specify that the instrument is unknown. 
     In that case,
     the "instrument model" term is used to indicate the unknown \
   instrument in the instrumentConfiguration.*)
        EmptyFileContent[Content_] :=
            Insert[Grid[
              Flatten[{StringInsert[#[[1]], "Name: ", 1], 
                  StringInsert[#[[2]], " Value: ", 1]}, 1] & /@ (If[ Cases[Content, XMLElement["cvParam" | "userParam", _, _], 
                   Infinity] != {},
                                                                     Flatten[{
                                                                         
                                                                         If[ Cases[#, 
                                                                            x_ /; Length[x] > 1 && x[[1]] == "name" :> x[[2]], 
                                                                            All] == {""},
                                                                             {"Unknown"},
                                                                             Cases[#, 
                                                                              x_ /; Length[x] > 1 && x[[1]] == "name" :> x[[2]], 
                                                                              All]
                                                                         ], If[ Cases[#, 
                                                                            x_ /; Length[x] > 1 && x[[1]] == "value" :> x[[2]], 
                                                                            All] == {""},
                                                                                {"Unknown"},
                                                                                Cases[#, 
                                                                                 x_ /; Length[x] > 1 && x[[1]] == "value" :> x[[2]], 
                                                                                 All]
                                                                            ]}, 1] & /@ (Cases[Content,
                                                                        XMLElement["instrumentConfigurationList", _, _], 
                                                                        Infinity])
                                                                 ])], Alignment -> Left, 2];
        
        (*Extracting the subelement "source" with in the Element <
        instrumentConfigurationList>. 
        This subelement indicate the mass spectroscopy ionization source \
      component.*)
        SourceAllTogether[Content1_] :=
            Insert[Insert[
              Grid[Transpose[
                Transpose[{{Style["Ionization Source: ", Bold, 
                    Underlined]}, ({Insert[
                     Grid[Flatten[{StringInsert[#[[1]], "Name: ", 1], 
                          StringInsert[#[[2]], " Value: ", 1]}, 
                         
                         1] & /@ (Flatten[{If[ Cases[#, 
                             x_ /; Length[x] > 1 && x[[1]] == "name" :> x[[2]],
                              All] == {""},
                                               {"Unknown"},
                                               Cases[#, 
                                               x_ /; Length[x] > 1 && x[[1]] == "name" :> x[[2]],
                                                All]
                                           ], 
                            If[ Cases[#, 
                             x_ /; Length[x] > 1 && x[[1]] == "value" :> 
                             x[[2]], All] == {""},
                                {"Unknown"},
                                Cases[#, 
                                x_ /; Length[x] > 1 && x[[1]] == "value" :> 
                                x[[2]], All]
                            ]}, 1] & /@ (If[ Cases[(If[ Cases[#,
                             XMLElement["source", _, _], 
                             Infinity] == {},
                                                        {"Unknown"},
                                                        Cases[#,
                                                        XMLElement["source", _, _], Infinity]
                                                    ] & /@ 
                             Content1), 
                             XMLElement["cvParam" | "userParam", _, _], 
                             Infinity] == {},
                                             {"Unknown"},
                                             Cases[(If[ Cases[#,
                                             XMLElement["source", _, _], 
                                             Infinity] == {},
                                                        {"Unknown"},
                                                        Cases[#,
                                                        XMLElement["source", _, _], Infinity]
                                                    ] & /@ 
                                             Content1), 
                                             XMLElement["cvParam" | "userParam", _, _], 
                                             Infinity]
                                         ]))], Alignment -> Left, 2]})}]]], 
              Alignment -> Left, 2], Spacings -> {2, {2, {0.7}, 0.9}}, 2];
        
        
        (*Extracting the subelement "analyzer" with in the Element <
        instrumentConfigurationList>. 
        This subelement indicate mass analyzer (or mass filter) component.*)
        AnalyzerAllTogether[Content2_] :=
            Insert[Insert[
              Grid[Transpose[
                Transpose[{{Style["Mass Analyzer: ", Bold, 
                    Underlined]}, {(Insert[
                     Grid[Flatten[{StringInsert[#[[1]], "Name: ", 1], 
                          StringInsert[#[[2]], " Value: ", 1]}, 
                         1] & /@ (Flatten[{If[ Cases[#, 
                             x_ /; Length[x] > 1 && x[[1]] == "name" :> x[[2]],
                              All] == {""},
                                               {"Unknown"},
                                               Cases[#, 
                                               x_ /; Length[x] > 1 && x[[1]] == "name" :> x[[2]],
                                                All]
                                           ], 
                            If[ Cases[#, 
                             x_ /; Length[x] > 1 && x[[1]] == "value" :> 
                             x[[2]], All] == {""},
                                {"Unknown"},
                                Cases[#, 
                                x_ /; Length[x] > 1 && x[[1]] == "value" :> 
                                x[[2]], All]
                            ]}, 1] & /@ (If[ Cases[(If[ Cases[#,
                             XMLElement["analyzer", _, _], 
                             Infinity] == {},
                                                        {"Unknown"},
                                                        Cases[#,
                                                        XMLElement["analyzer", _, _], Infinity]
                                                    ] & /@ 
                             Content2), 
                             XMLElement["cvParam" | "userParam", _, _], 
                             Infinity] == {},
                                             {"Unknown"},
                                             Cases[(If[ Cases[#,
                                             XMLElement["analyzer", _, _], 
                                             Infinity] == {},
                                                        {"Unknown"},
                                                        Cases[#,
                                                        XMLElement["analyzer", _, _], Infinity]
                                                    ] & /@ 
                                             Content2), 
                                             XMLElement["cvParam" | "userParam", _, _], 
                                             Infinity]
                                         ]))], Alignment -> Left, 2])}}]]], 
              Alignment -> Left, 2], Spacings -> {2, {2, {0.7}, 0.9}}, 2];
      
      
      
      (*Extracting the subelement "detector" with in the Element <
      instrumentConfigurationList>. 
      This subelement indicate the detector component.*)
        DetectorAllTogether[Content3_] :=
            Insert[Insert[
              Grid[Transpose[
                Transpose[{{Style["Detector: ", Bold, 
                    Underlined]}, {(Insert[
                     Grid[Flatten[{StringInsert[#[[1]], "Name: ", 1], 
                          StringInsert[#[[2]], " Value: ", 1]}, 
                         1] & /@ (Flatten[{If[ Cases[#, 
                             x_ /; Length[x] > 1 && x[[1]] == "name" :> x[[2]],
                              All] == {""},
                                               {"Unknown"},
                                               Cases[#, 
                                               x_ /; Length[x] > 1 && x[[1]] == "name" :> x[[2]],
                                                All]
                                           ], 
                            If[ Cases[#, 
                             x_ /; Length[x] > 1 && x[[1]] == "value" :> 
                             x[[2]], All] == {""},
                                {"Unknown"},
                                Cases[#, 
                                x_ /; Length[x] > 1 && x[[1]] == "value" :> 
                                x[[2]], All]
                            ]}, 1] & /@ (If[ Cases[(If[ Cases[#,
                             XMLElement["detector", _, _], 
                             Infinity] == {},
                                                        {"Unknown"},
                                                        Cases[#,
                                                        XMLElement["detector", _, _], Infinity]
                                                    ] & /@ 
                             Content3), 
                             XMLElement["cvParam" | "userParam", _, _], 
                             Infinity] == {},
                                             {"Unknown"},
                                             Cases[(If[ Cases[#,
                                             XMLElement["detector", _, _], 
                                             Infinity] == {},
                                                        {"Unknown"},
                                                        Cases[#,
                                                        XMLElement["detector", _, _], Infinity]
                                                    ] & /@ 
                                             Content3), 
                                             XMLElement["cvParam" | "userParam", _, _], 
                                             Infinity]
                                         ]))], Alignment -> Left, 2])}}]]], 
              Alignment -> Left, 2], Spacings -> {2, {2, {0.7}, 0.9}}, 2];
        
        
        
        (*overAllInstrumentConfigurationTag is to output a readable format*)
        overAllInstrumentConfigurationTag =
        (*we use if statement in case the subelement "componentList" is \
        absent*)
        If[ Cases[instrumentConfiguration,
          XMLElement["componentList", _, _], Infinity] == {},
            Insert[Insert[
              Grid[Transpose[
                Transpose[{{Style["Instrument Configuration: ", Bold, 
                    Underlined]}, {EmptyFileContent[
                    instrumentConfiguration]}}]]], Alignment -> Left, 2], 
             Spacings -> {2, {2, {0.7}, 0.9}}, 2],
            Insert[Insert[
              Grid[Transpose[
                Transpose[{{Style["Instrument Configuration: ", Bold, 
                    Underlined]}, {SourceAllTogether[
                    instrumentConfiguration]}, {DetectorAllTogether[
                    instrumentConfiguration]}, {AnalyzerAllTogether[
                    instrumentConfiguration]}}]]], Alignment -> Left, 2], 
             Spacings -> {2, {2, {0.7}, 0.9}}, 2]
        ];
        Return[overAllInstrumentConfigurationTag]
    ]

(* ::Function:: *)
(* f:MzmlDataProcessing:: *)
(***Function***)
MzmlDataProcessing[dataProcessingXMLElement_] :=
    Module[ {
      dataProcessingElement,
      dataProcessing = dataProcessingXMLElement,
      dataProcessingTag1, dataProcessingTag2, dataProcessingElement1, 
      overAlldataProcessingTag},
     
     (*Using Cases to extract the file subelement "processingMethod" \
   found in Element <dataProcessing>*)
        dataProcessingElement = If[ Cases[#,
              XMLElement["processingMethod", _, _], 
              Infinity] == {},
                                    {"Unknown"},
                                    Cases[#,
                                    XMLElement["processingMethod", _, _], Infinity]
                                ] &@
          dataProcessing;
        
        
        (*cvParam holds additional data or annotation. 
        Only controlled values are allowed here, 
        while userParam holds uncontrolled user parameters*)
        dataProcessingElement1 = 
         If[ Cases[dataProcessingElement, XMLElement["cvParam", _, _], 
            Infinity] == {},
             {"Unknown"},
             Cases[dataProcessingElement, XMLElement["cvParam", _, _], 
              Infinity]
         ];
        
        
        (*The attribute "name" is the actual name for the parameter,
        from the referred-to controlled vocabulary. 
        The attribute "value" is the optional value for the parameter;
        may be absent if not appropriate. And if absent we output "Unknown"*)
        dataProcessingTag1 = 
        Flatten[{Cases[#, 
            x_ /; Length[x] > 1 && x[[1]] == "name" :> x[[2]], All], 
           If[ Cases[#, x_ /; Length[x] > 1 && x[[1]] == "value" :> x[[2]],
               All] == {""},
               {"Unknown"},
               Cases[#, x_ /; Length[x] > 1 && x[[1]] == "value" :> x[[2]], 
                All]
           ]}, 1] & /@ dataProcessingElement1;
      
      (*capitalized Name: and Value: labeling*)
        dataProcessingTag2 = 
         Insert[Grid[{StringInsert[#[[1]], "Name: ", 1], 
              StringInsert[#[[2]], "Value: ", 1]} & /@ dataProcessingTag1], 
          Alignment -> Left, 2];
        
        (*overAlldataProcessingTag is to output a readable "Data \
      Processing: " *)
        overAlldataProcessingTag = 
         Insert[Insert[
           Grid[Transpose[
             Transpose[{{Style["Data Processing: ", Bold, 
                 Underlined]}, {dataProcessingTag2}}]]], Alignment -> Left, 
           2], Spacings -> {2, {2, {0.7}, 0.9}}, 2];
        Return[overAlldataProcessingTag]
    ]

(* ::Function:: *)
(* f:MzmlSpectrumList:: *)
(***Function***)
MzmlSpectrumList[spectrumListXMLElement_] :=
    Module[ {
      spectrumListElement,
      spectrumList = spectrumListXMLElement,
      spectrumListTag1, spectrumListTag2, spectrumListTag3, 
      spectrumListTag4, spectrumListTag5, overAllspectrumListTag},
     
     (*Using Cases to extract the file Element <run>*)
        spectrumListElement = 
         Cases[spectrumList, XMLElement["run", _, _], Infinity];
        
        (* The attribute "id" is a unique identifier for this run. 
        And the attribute "defaultInstrumentConfigurationRef" is the \
      attribute that MUST reference the 'id' of the default instrument \
      configuration.*)
        spectrumListTag1 = 
         Flatten[{Cases[#, x_ /; Length[x] > 1 && x[[1]] == "id" :> x[[2]], 
              All], Cases[#, 
              x_ /; Length[x] > 1 && 
                 x[[1]] == "defaultInstrumentConfigurationRef" :> x[[2]], 
              All]}, 1] & /@ spectrumListElement;
        
        (*capitalized ID, Name and Value labeling*)
        spectrumListTag2 =
         Insert[
          Grid[Transpose[{StringInsert[#[[1]], "ID: ", 1], 
               StringInsert[#[[2]], "Default Instrument Configuration: ", 
                1]} & /@ spectrumListTag1]], Alignment -> Left, 2];
        
        (*The subelement "spectrumList" is where all mass spectra and the \
      acquisitions underlying them are described and attached here. 
        Subsidiary data arrays are also both described and attached here.*)
        spectrumListTag3 = If[ Cases[#,
            XMLElement["spectrumList", _, _], 
            Infinity] == {},
                               {"Unknown"},
                               Cases[#,
                               XMLElement["spectrumList", _, _], Infinity]
                           ] &@
        spectrumListElement;
      
      (*the attribute "count" is the number of spectra defined in this \
    mzML file. 
      And the attribute "defaultDataProcessingRef" is used as a reference \
    of the default data processing for the spectrum list.*)
        spectrumListTag4 = 
         Flatten[{Cases[#, 
              x_ /; Length[x] > 1 && x[[1]] == "count" :> x[[2]], All], 
             Cases[#, 
              x_ /; Length[x] > 1 && x[[1]] == "defaultDataProcessingRef" :>
                x[[2]], All]}, 1] & /@ spectrumListTag3;
        
        (*adding  "Number of Spectra: " labeling to the final spectra count*)
        spectrumListTag5 =
        Insert[
        Grid[Transpose[{StringInsert[#[[1]], "Number of Spectra: ", 1], 
             StringInsert[#[[2]], "Default Data Processing: ", 1]} & /@ 
           spectrumListTag4]], Alignment -> Left, 2];
        overAllspectrumListTag = 
         Insert[Insert[
           Grid[Transpose[
             Transpose[{{Style["Instrument Scans: ", Bold, 
                 Underlined]}, {spectrumListTag2}, {spectrumListTag5}}]]], 
           Alignment -> Left, 2], Spacings -> {2, {2, {0.7}, 0.9}}, 2];
        Return[overAllspectrumListTag]
    ]

(* ::Function:: *)
(* f:MzmlMetaData:: *)
(***Function***)
MzmlMetaData[input_] :=
    Module[ {
      (*overall is equal to the only input to this function*)
      
      overall = input,
      metaDataOverall},
        metaDataOverall = Grid[Partition[
           (*this is the viewer title joined with the other functions*)
          
            Join[{Insert[
              Grid[Transpose[
                Partition[{Style["MathIOmica MSViewer", 35, Bold], \!\(\*
        GraphicsBox[
        TagBox[RasterBox[CompressedData["
        1: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
        jpmo8wy6M5VMP9O86qJ/D6nlkJYf1sdHfxUz1VfC36OuGSGfCumHsqwcao2s
        PNq1CuxQ3ATYd/Dw0N9qwS+r1rfh/1HHVmX/T1css4vOGMkyyNdH766lB8ZG
        lhQyoy1mKT+d/D3kT/Fn5jVV/a2RlaprtWRnAvtDkfeg5ZfVuF4t2U7+n65Y
        ZhedMVLvc1dd5VgLfljK7yEd1pK/h+IP6pi9xb4dVlyz2s4E9ofEx+T4Wcgx
        mMK0beEIceZo1Zl0iuWcaZ510dlD0l48YTrO7zXscWfXqjvV+9nLS79VyL+j
        6bCyj1POCz/qHrcSl2qZ2VecAfY698oHO+UP+lXsYZ3jrHIsUz5LkPtgtd59
        pvnbRWeePydJtNcMfhTMFVHv57DgsmZmWYfViglSHgF0Qtnfc8p0QhUuL2Km
        HCBDngFmTKvvQ/RmV+uyScxUm6nlpzoluuwZ4GsXXk8/Fs8/BXS07zXISk9P
        fbqW/idwhniClg6rlUMj+3so5+6UxgOle6C8djlvdnBsVHGUZrwkPKFIde1y
        XEcdy9TqM3IiOsiJnI2ng79dOD11WFTXUtC+s5YuUJH6+y9ydatcot7HbrrK
        55mpnkvOUbdU90w5p+o6jVOaP6NaLzlmAt3TkM8Ou7gTLkPLVjBzWa+uDyVd
        Fn7ZY/qQTxJ7al6cav4eE6tnev+fTXQMfmnxh2qDqT7YGqLv6cZM/lD6G6/3
        UPbxk+p9HOBXcRkzr+eyFLMkXfBYMfpTZqNJayj7fzryf5Bjr74fRb+9n5nO
        Fa3rJz+yg4XrP9a+18KeFl9PN2/VeNW83jONibMMj1pYlOvzqUa/h6CeEvWS
        6BINkl/vqSK5P0d36TcMcQbN+EhUcSnTrn3OZp1zZygOCB32tOe0qdazky6L
        uWLq+3HVVX7OzPPs1fKe/LJyjUkf1tlndSy9RsaeGkcdvPXy0mfI5O2tz4yJ
        LB16FP7KPLbEXxmvMk67MGoZj2ZY9PTShyM2AR1SkVd3ObtWPaElt9QE/4yL
        a+VTBnKrvEfnXrkM5OOjr/L3q8gVexPUV+wr6mcl97Qi3vfBZzr9hmvVD8w8
        R51kJfVOl3MJosR71Xua7LNTHgPUwKWZLqucKZGa+YQx2VVM2z7GPci6LNnH
        fZipZx5h05Juc7Gylh6KrE5T+DpF5quyngePxVctIj6Dv8p+GR8cWB6vwV+Z
        tzJmZZyqMdoJn2caM//PmKRz9UJPT70b6o3Bp5Ph0YmQAbtulU8rNuLNCj9n
        B/pX5MfHFEKW9WfmfSAHKvtohfrzHh76R5h5T9gMgUeaVYl9LccAybaUazMo
        l+C0xBmkNSZdlvyyBl1Wy8aMjiidxyzrsrhfquGm+6DzhbBp2Of+vhUeClZK
        kbeI+NLJ4u5keavg/j4fb/2EhJgC2A1yf0/Ca2/pumWMynL0X4NPDUxeoJx3
        Y7Ty284EYZ4IeKroTXODA8pHKNdnp2B3r/p9AX76DcykvxIm8xTKF4/Y16nM
        lGfqz0x1lnIt42mrmWKdZaZZXh7mcXa+r4pbpfshuW8p394gM2Mj8519vMuz
        Xd30SxWsK2cq4kmVZ5yXROBfSGB5kThHYBfbCKwSTvuIM5L03WPi80zj6P8B
        j5pz6c4cVZrIuZI7OOlx/v6pfl9sZAnilklCnmQLPI4QlC+eo9w7kjPq/o9q
        n88plZkauJRzDGz8fPUT1ffl4al/mZnPGOtUnxYVlpnu61U43dWt+HoHp9LX
        hjmWK+tUoayX3rBmpjU807w0J8ylCAsuh08A/mQ7gdNBzNSfl+RoL9a519k/
        DptqTCrrE3K2yMej41LPHbDXNN7n7Z63NTwkvYkZ9VZgEXX/peKxUDwvy0vC
        pbr2n+yy0+JzkL7nAmaeYzAgIqRMs8+muO4sJunkAT7p83XO2bvt7HLfHjKk
        gNvZjeBD7Yu5vUMpJ1w6dMKlddhU9smb6H8ik1ZvvxMlxW7ZotjR8GE5iXOS
        MAo5OkBglHqmkH6ric8zja1ThUnISEeNuQFH4c1zyHVT9MZLo8LKqqPCSmti
        olKhP2Jd4eOELwI1xjQTMjwhJi8/JrK0EYQ5HYF+FZsU2XCzu7t+n1afY01c
        OleK878zLiFDBw/O5zY2OXzQoKyfXR0zHvP3TrvGQ+czVvn9CoFNyMwMIWco
        dol6Kdm+VM9nJn3WYqzBSlzKNqYhX1ar7gu5sqGBCTVeuvRNdrbpr/bvl8EH
        DMjkNgOzlfvN4UOG5HNb4HJoEbcfpuDSocyIS6fjw6Wbe+ULCi8e8vXRr1Ps
        +mXR4aW1Co0RfA0SfA1W8TVMnG8RgiIVPhfGRJY0gccKf3e6e1TuPxHMx0Vn
        4Ptg88OP5SjOTFtxbg5gpt69xCPya53T2LSAyaOuF+w75TzbLHojw3/py0xz
        HkOY+TxWxLuhW4Gf8gzlY1FcfHRhrYLxGUH+Fdcrv/eEaZbHsXGJfTh4cJ6C
        y2w+cEAWH9A/g2Pv9uubzocNTX/DxyNjW2hQQqOQNSnMVGOM+6E4CemyFGcg
        uXlMv+CJ7Adae2YexyQbcwjOPfX9DR1a+CPdD8iAy/6ZfODALD5IOYtw77a2
        hQKXJYb1MMpMvZkuq5yBbyh68R7gD2dqbHQ6fLry3E46V2XMga9REm+Pxd8E
        9XMKxmf5+OjvPZZsxesJMQXQ2b008CnPn5Blp4zNc1KnlTGp7HvNXonS+fVs
        cGDFOGbU9WQ8Bgu+0VzkWIkfSYKSBQ0XlHIUovfQZ+g7EqMjk3NiI4vbgvzL
        t7i7V7zs5Kz/0bjHVPtW0d8GD8rtwGV/Ay477+GhQ7Pf9/PKvSYyNLWMmfsy
        Kd/HRuI72TSW/IInrOfK688syEs3XUWnXPyhCtY631OG4Z5sOnCp6LIKLocK
        XDo6l/9Xpyt/OMC3/LKo0NIawT/CIOk2JP9k/Mm4I2wlCjoZ/nbwOCoiMU85
        35chLmtp3ymy+21FTscL3hA+cXbSvCYbZuprpI4HnXPYlPeEv1/FQEu6K3za
        isyqE2uBNfGR8Ii9HC14RjhMltYfdhB8EulWUJqgVIlS6Hmta7a1KzToc9in
        wKBF2TIox6DvQrY4OI54x8+neGVsZC6+n/w/NKeL/IIUY5P9gscVW5P3hsbz
        hEuDfenvUx7t4lp27TCH0iMGTKn1AUX2yToA8AmcGnRZ5SzCmWQ/tOBzD7fC
        vcH+hUtiI7NGCKyRTUGP4YLUGJTxRzxV4414azV/I8MTRgQHlnfqjd9hb3rq
        H2fG8wIynPAJHlHND/FIK1Z7TvmCpGu9QLFhtmmth6LXvpYQkwscApOeYk3U
        eEyUeCbzCLZb5imkDNVjVlxMVq3WdQd4597h7JjzlhqX6v07SNm/ZraY8JG4
        Qa74lU1DLSQz6Uyyz4Fi4UeLral9uFq52h25Gt5eep0z+vg56z+Cvgl7EJjE
        NWndo/q+bAZk/OLskPlSoG/G1eEhSQ3M6BOSc2jjBeZk3fNYGLSEvQyJThl/
        46KLp2rNuAcJ/zrOD+rLjX1J/XMRD6I52oRN+bw8J2xNJp3TyBWwhEllX9Ic
        DJxP/mJNIgU/CY+pzByHWRpkCWdadDTe0exmQ36AotNeqr5ub2/9m8pr9Qq1
        ODvpZgX5pu7UuaS/NHBA+q/Yv2p972i+S0VX+MHdQ787LLi8nJlmzlJ8TfYL
        yjJU7R9U+wk78t5Cg8oHGfvCw/ct+ZgJl45lBr+Nln8Gsn7QoIyv3F3SHgj2
        T16pfB90mjEKjVZIz4y+Z8hJwifl7Ml2AvFPPlPVfNDipxZZ4ufx8LjjexR7
        sllrP4o6N5wfNKubZpsBm87MeH7K2LTYC/9M4+8YmDTISsW2vkpLd42NKoY8
        JExSTUaMWJtkZsLj0bCYwTrrpKnHSZbOaXwv/DW5yr6+WX3tis6Nvo44Wycq
        hBr/GQrNAoUGJG/w0mU8YWeb9W0nO0zBpSX/CAi5N6iLFGcVxddoFq1az9WS
        ox25pm7u+hLl+24w94tKMVlnwmW50Zeq4dsKDy65iRlxCNk4XhB8WePE86OY
        0fdM+MwT6ybzSxMbzHx+vfq9aoyp7Q0tm0Ot/6rtG7MzPTykbIMWNuNjCnDe
        Jok9SHmH0GngR3dhJmwetff2mcbgsWQls9Br1ctLv4ppYzKJmfo70jpS3ol8
        bhKfiBdm/pujEL0nmZn7EtQ6lSHfxcdH38lXGRFSdpfyWotC7QphvgcwOUcQ
        /p6m0KTw4OGX+njmPjhsWOFXdpLf0hRPsByDR28sZe9gz8P3YEmOUo6KAafB
        AaXeLq4VS51c9B8afcmd8Sj7mHENwCVkuJNz+XedcBlSdjszYhL32iYR/g18
        IiZUrRDyaeHXorgt4kMU98yWHtU4VONW5ithT/bPyTy2RFp8lvnbgVHMTVLf
        c6B/xW3it/EZwqYsN0mnlfNCKF/rrMYmk2SlpV6rcdEjcM+wJ2VMJoq1o/o/
        macyHolfWHfykct2TbRYT/I1yBTNzGMrsh+QeNqhO+vc9e+qMaPoQZczy7ic
        LZ5rF++BbKmJicic6eNV/KijU9nPJxLngwz1861YHB+TjXXSwuhAP5+yPGeX
        8puM32v8bgdnCffOnXGJ151dK97w9y1dFheVkxsSVL5G/ds+3vq3xD2QXkA0
        SaFWZpKdNcyITejiMjZJflKusIxL0nvVeg7pwYQrdQxEzWfitczfGOl9sn2L
        7yOfYVpYcNm16nt2d9e/x7T74hM2ZV/6OTPjjKlkJfKWO9mVxl6r7sxU5681
        i474SbwkPCYL7JAdgDWDLQBskx9QpmCJyEdIpI6ZEU8Jp0laemBMVCpwN17s
        T8jGmcyIRxDJS8Jlg7Rvoe8VRocXLfHyKn/mRPLVEPuHDBWz4Z2jI9L9vT2L
        Jw1zKPnYYLMK/bjDdiV8dshlIzm76D/x9q5YEx9dgLX2F2sSo9gUnWq+FN59
        w4zyEVjEPc8QNFXc3wRmxCbsTticlcyk0wKbecx0tsr6qYxDkmUkC4m3FKuU
        MUi8Dpf4rcVzOc9Ezd8Ovz7qZrTWmpnOEbn/C/Ueg34HPU+dS3nacpxPtazE
        dWr1QMb+YsY6RKqRimfmvaq0MJkieBYn1jdCrD2+I0DsMeDcTzz6iN/wEeQr
        ve7HTHE1iq2pcxai4qILa7RwqbxWy4x4g8yYLPbsTLFnp4vnsJ+BXbU8Qe40
        ZElOZHh8iXJmX+PmXnn4WLg0k6EuFZ8rturPBptV2K0dsX0DPsvM8Kng8Qd3
        j4o9kSGlI8Wa+Iv77bAd0G/dwh6dzEz2M509uNep0j2SPku4JHlJeCSSbQ/S
        T5OYuVwEf2OZSQ7KZy5hkHhOfCeeq4n4LOeFhYvvxW8kWLhnuR5IPQ8R3yf3
        uD4nZCYzx6Vm7w3RA5n6x8lzPAiTWBcZk6RPxIh1DRY8wdrj/PIUawUZ7CaR
        TiJ38R5P8Rkib4mPcg5KaExEcaNaF1Ts4teZ0d+Bvdgs9ib2Lul47eI5vAbs
        yv4R6HiES/keMxXduAm9Aiz58DVJuTb4UsmnRPgk+enkXPZmkH/J7OjI1Bhm
        OnsophgpMJAk1jjXxbXyZ/VvJMTkX8uMmJyr0Dzx2GFDMxMua5jRR32jQjco
        tFWhdxTix0EvKXSNQlcLIhuE5CLpPYGCRzIOvSX+ewo+E9FznuJ9fsyEUUOe
        inImvtPZxiqayUzYxHmCc4R6Gcp9AeEDkvszyb7ZswaXTKXD4joV23qp+r6R
        k8XMZSXpr+raRbWOT3ElwiPwBh+ZMzOeX7C/HKVHIifxHhexnoRT4hnhk+Sq
        4YxV7LYNalvN37fiUWbEWLXYh9DlIDfJJ4J/N4vXIFe1MEl2VzYz94dkQm5h
        JoGiQ757IvgEDg2+JfuiH3W6orujwnMhG7Vy2ihPCrKK+oPhmoq8vPWHOu/R
        wm3MhEuQLC/xCN/Qi+z48HeihL6ZOMNlGenHTLoQ8Y/OXlfBYyeJiO8gnNUe
        TNIZkOujcRa1iH2Yy8z1WeqlIs/SVs84o/qDcw6XmCnOzGdeqfVXS3a3n+AD
        5WFAx4cfBPGExVbwf65YZ1CHDPX2LLtRbaeFBJbtUF7DjB3E72qYEX+wsVBX
        AtkBGQl7i/Q6+CmBxwJmLiPJ9yH7mzvI1zNtla1t7neQgceLTyeXyt+C/Mtv
        jopIhB4pxxRk34r8m1QrivspC/SvUNWiV2JuCeYjAH/k2wIu9yr0sRXrfTL0
        nOAxYVOWk9gPwKQzM53LFGeyF3tkmCC87iI+521Bj6U4tiwz5Z4NR+sHc9bp
        skxlWyrU3dNLX6HhT8D5Ks+8IllJMoTqcBNYZz+1MzPlr8HuprjBQit4Plfi
        cYd+6+xS+qK5L0XPoyOKF4m9DDsKNmMNM+KQYu7AI2QV2ZIF0j2R70odQyXf
        fRbwaDMw/UvK6QMhRwH5Qlr581oEfTQooOKuuJhMnBGFKhohEa4Pchxnh16x
        dfcQHomCA8pfYCbf1l4r1vhU0RGFVjOjvPMRvCJsAieO0v6wFTRE+psw6hQX
        kxWs5XNj5jaVLCeSmPl8RHXN3rmCyx5B/hWxnfwWrlU/MiPeEpmpzk9eAy1b
        W51zQbF26s0zzwpez2cmPbbDV+ToWPpfky/FiM3YyBHt4jqxr4E9YHO0IPxN
        eTCF4n0k93E/6joI3H9ygF94vrtrxvaBAzN+Qt6bkUy4NOXdphty+hydK/4+
        XhkKGRgbnT5RujZQhbj2Cum5qrio4tVqH5ePdwVsxP8o9IsV63s66BNm1FVI
        dhI2sUconkR5GaBB0t8GORoUUD5O7Q+HvcmOLiuoxp16wpz2voanGJcd9bda
        +yU+ekQ166zDUm283LNC7oMs50HJdeU4p+ZawWPoR7BfyK9g8NEOFbUSchwi
        MiyhQFwvHiFrsL9pbwOP0Anzmbnug/OFfBkdtYMhQTHpHq7Z1w8elP2jIc92
        QKaRkC/en/Bpwuggm/SvfDzSdri6eDTHxxQsRz7g8eJTOR+fjo1Og02IdccZ
        MlJQpaCq+OiCpTIm3d0V/SBy+P9OEY5OF+1hplxz4JN8DlSrJfcksBGEv209
        3Cuu63QO+ejvZZ39kLL/R+4LSDYm5f/Ic5/OGhuTmduWHb3X0NtKvU8Qx2bm
        PZDVdqVsY8t9kCEnCZOUhwaaYwVvL2NG3FDtQ3hUWEYedEe7jjiEEZuCN7hO
        2GaQmdSfgPw6edI9xInvpXgO4T7Q3S3v2iFD8n5EXcYgQ71YjrFmTMGngQYY
        8+GB0UE2GV/5e6VuZCadmWza5oSY/KsUXfP548VnREjZnYH+Ya3iu+CTghwl
        bFbS/lTsVO7tFXCmMXe89AYz6ifqWhBnZl6zI+NysJOL/mN1bgfqq5kplk52
        JvkE6JylHkdavQ3PKlwyDZ8PE7hErzn1/kANMrOMyyRmPquYfNIUx6UcbvJ/
        WYtL5ARSPA3rHhcakDXBUAdiW9iBTVfXkpeZqRdlrtgLwKIcR89gpnOF+Ee2
        kLevZ/7MoUML/4vvRT67kfINebRGjOYY6o6B0WH2Wa+FBqbCbgZmIOcQj6Bc
        VWAL+iliNFMVuXZZcGD5S8dpf/6CHDsFn/gewqVBvxV780zj7GQIPUNxZkWI
        Nce+8RR7B9iBXks9twYH+pcVauVCJsTkQQ+ieLrce4zsK+wP8nl4sXMYl9Hh
        pcFa+wM5Jsy8v5PWvVvyexEmqT7RGlwiZgZ+UH5KeqBv5jxj7aQp79zVtfgV
        ZsrbzRGfKZQ+K/elJEwazm8fz8Lp9vZFn+B7OkiRx6jjNGHUiE8Hh9yDYUEZ
        8LMA89CTIdcQc0GOArCEWAzsXMrDIX/pbAWfl4cElh+jl4boJaDgMyyoDDke
        wKQhvzU2suiT04CZ/y+C3wJ+ccrPwfpTrZYzM9mdQ1xcy3ap85QVWYE6Ieo7
        KssMeW8mMfM5M8CleibbOYFLhWwUXfZ19f7w9dHvVd27lj/aS6ypfO9yXSph
        0xpcXs+McoPijCM83DJuMqtrVvDj4zkCtgzlCUJeFjAVngXf6Nr9g/3zRw0b
        NuKQWdxf0YnxSM+BgFMnp8L94cG5beJ7gHtgvkxcG3AJ/RWxUchKnGlyLQvu
        n2KMcwU+9x0Nl5Qr6+pW8VVkWBHi+hutxMXZQMBmi+CDPLsBtiBkp2NUWEGU
        7MszxcDKsXYRrPOsxGPh8rT36D5F2CS/TwcutfohGvw/Mdmj2InjUu6BRHqs
        NX4f5Khgz3fEHT1c0nZRHaWxtrmA+/kU3shMcS2K/VHuDl07bJOI8ODssmHD
        Cl+008qTU5GzS8lLESEFDYLfieJ78DvAPOxX6Jk1zCgvG5kplxw6LMUw5Dy5
        jnw5BZ/LO+u35rjE3gwNLDnTeDqV9JNYnzRmytOBbQ/Z6ebsUnwr+Quk/po/
        xESlAY/ApdynW441a+1NNS7PKj3WAi6pj8xARW/31Oq7puD1AXb0M4l6Gch6
        rNzXm37HmjjJbmaUP5BV0IPqbIekvNGvb5rB90J1lCEB+VcyExYLpMcccd3D
        w4KHl7k45d7dUQst96WScstBzi5lL0SEFOMcgL4F/SmWmXyBpCdTnBTvw7lB
        9uUEsfegz8LPOk3QdEHTxPOG3MCYyJRFZv4hqbYk0K9UuVa7M42lU01fMKO/
        GXsKOqkBS5GhmenqPn7ghbtHOfKW4M+Q87XJviRcyvl45PeRc37OZn+sWfyS
        mXru23l766/Ukplx0UXTmOlcMusfzEx+H3lWpFm/fWaUn9bg8k5m1AOhE2Kf
        N9sOSnqrb+9UQ4wCMhM+mYjQLLyHdFeK0Rvys328g4rdXTO3DhqU/aO6J5Wx
        R4HpbHZ2KT8UYcwhp5xcqqmR9wP5laBXQ7+GXlHDVP5YZp4HOFEiygtsZabc
        wMb4mJzFhviKwKVOV3Eu+V1PlN4T65cr1jXG0SHvMZyZamzGReVi74UzU72h
        HMNU95xXz7Q4l+KX1A8RZ4ih7xp62WjJTE9P/UFmPi+SYgxyToEcvyU9vjsz
        zduabwX/ELNawow5Q9AH2/v1Sfymd6/hnGQm7MyIkExgH76YYmbywxYEeKcu
        GWST8TliG4aeVKKHo9wrDph0cq74JMCvDLKM8qihW1HfFPI1UNwMejJwCXkJ
        XEJmApvVzGhrIrcIGIWPFvptvQUaJ94zRhA+Wx0XXbDUy0t/KDY683Rh4jeF
        XhCEHFvYrpsEbRR0s0LPC/rtNF3HLczoO8sL8kudQvzpwKZiY+h0JahzJ0xS
        Dhr59wiX6hxZ9UwL7MuzPd9H9v10zMLAPaCWREtmoq8DszyjWB0nIl2WMAns
        L7CCdw8qhFrnZQohfjWrx0Xx/JKeSdwoMzMMeGOmmnxQcZBfwmTbIamvI+bf
        YYvamPxEHTqsY+kPfj6lS5kp91buvSnLScpfhQ4LOUzyspSZ5+dUSlQlUaWF
        19TPUx7Eaiv3vBah3gbnHHKgtwnaotBmZvSvXcdMuATPUa9yjaDd4vOn+ppm
        ODvpxg4ckPYF+ET9/IBNe/uin2IiM3G+yv2kZH+PXNeUxDrH7+SZZucSLjtm
        YShkHx+TFaTVaxf1TZERiSVMYxYGM6+rkXur0JxD/L3ACr49ptBasU+Bzfnd
        u8XyXt0TeZ9LUjp0WSZifDo3r1oXh5R7+vRO4X37pHXIVMKl3P/Ox7P4muiI
        FNwH1QFS3WOEuD+5tkPO/SLblWKk8iPFS4nI1qValQJmbv/KubFFzJQX+/kp
        2PNE/1UIM22BLeTtQRaiPxB8ZYRN4BJ4pJqudYKw9ldJtEuhj07htb3lPCzx
        Ieg/4JcBmwON2PT2zMe1JTHTuSjX5Wvlrctzk5xZ5xgJ9fo5KzCpgU0556e3
        uHbcg3NwYMVMLZnp461/lnWuJZF1eSdm6rMrzzXEozV5608w47mNPXK5s5PH
        2ou7xfCe3RM46bIClyP9vZIu79sn5Wdg0oTL9E4zAnSuhffERGQCC3JNL9U9
        kr5E5zP1Tuno98VMPl/6m/zVcg8Oyn0friKtvkUdufHiu7adwn2PfIv7FLqH
        GW119MgBPncyIy63MiMmIR+BRzoDkc9xhUIrNQjPP3mqrtHZPp736pHIe19i
        wuaQIdlfMJPtSFjMU1EuM9fhKDf2nMjB08Cllo05WNwL5iK+qIXNqLCyVcy8
        Jyn5ZaH/Uf8GmrlI2ISebA0un2ZG/QrYXO3pFnrjRd2iec+L4xVcJhv4aD80
        813bQclvgK+QoSDouB247G/E5VD73HfDgrPhOwLmqA+G3HNGrn8k34Jsx5jV
        YzLLPY2ob4OcdxsuUUcerngf6Wn4XZwXP1q71y+++OK/e/fui1wLzAF9SCH4
        1qHD3qEQ8hWAS+iz0GGBSZx7yL+EzbBcoUuZMTd5keDfQvE3CPb+UvEdVuPS
        w82P97gorgObOGtDA1LwO9SPSK61Id8BxcAymMm2JF+kbFuetbUkx6nLDhL3
        osNMCEv6LPp7MnOZSbMjKW+D+pHJs2MXWcG3Z8TegY61VsHlTcBlh42pYK9/
        //S/LlEwCpyCrwZ8Cnlp6Hk8MPMnf59M7LkcZjpfCUeJzLw3m6V+uFr9xTpy
        A5mpxw3V8AcyU/8US71TkG8k27P4jqut3ec9evTgvbp5Kn/r+OB+Hh8yoy3w
        MDPKTfhRgEvosVuZ0deznhllIPAILC5gpv6B01Uk90yCr227tdcLcrSLMfC0
        V48k7uY0HNcKPxj84nJ9TakgqtGjmiDwQI5byrG7c6mPiBqXsi6Le/IKCqhY
        qSUzUauPun3WWWZir+Gcgm+aZnsMFGRNXfSzYv9AZq7zdA3ZaY5LY7wEfxux
        aY5Ld5e0vR46H/hGKXed6tzVvW21sKjuD0c94YBHyEWKdVPuO+XaUs8MuSbY
        iXWu03djpt4MwCmwbHWuXc9uHsqjWwdddJH7T3a2jgeY0cYELmFjQoeF/go5
        CT0I9QHzmRGLiK1SXz3kSYyXqJmZemNC9wBOb7f2mr3cgznsk7594r8bZu88
        i5n6olHdLNXZ6AUvqWc1eCLrsFRHAb2N8rXPWh32OHRZ8ss6invzw9w7LWz6
        +VY8yMxjJrQmdFZB7lK/RsjOJVbwDLXwsLeAy/UOtuHPG3EZ14FL4E+NS3vb
        1LeC/BMWMGNcEXyFr1OuhyY9SO6lqqWrUi8B3KssG3HPJBeBRdg0VMcE7Ml1
        hh052YKGiOftxFo5ic/WWru/7QYCk64KuQhyM8hN5fL44N6+H/bu3Qcyc5dY
        U9gGyMcgTEIeIr4K/CGOM0asn1x7NlqsJ16jnGBg9Flrr71vryju4xUGnXqO
        uBbqHUb9cGuZqa69kJlmC5Mf1pAzxEy9fc5aP+wxcNkpV1bcE854b8xRsjQL
        TbE1r2Dm+T9y72tnZqqvAz6XWsEvxM/QK2qzm5P/Hd27h/12Ubco3l3CJWQj
        4bJvn+G/eOuSIQ+QG9TATP3HaW+Bp6QLyf5SdV+fDGbqD4d7lHvw416pl5G6
        5pfqlmjOMfVix2Nf8Tf920a811Z89iYr1gn6628XMbc/jXh0lnDprhD0Wh/e
        87yAn12dfGBvQn+Ff9Xg42ZGuYcch3pmPktBvU7kQ4bOVCHWFHiBfP3amut3
        c/ZGrgHsW9i2FK9GThThE9dWLa6tSPAJvIEOSzWXzqxzTuy51HNdrsXEtZP/
        B/fUITNjIkrrLdU+YM4LM49n0pkl9xQZJtb5ZPmF2PeNg/v6v84uDOLndwvj
        ZrjsbcQl7JJhtkmvODq4TmWmPBvSwaivTw0z2SxlYm/JPgTyq5INSvqqFh7V
        vYxoXirNfTPrt87MY7o9mSmGBKwCv8Dn91asE+htBZsH+16s+8Yclx4GTCpq
        NjeuYSh3Ghp10MZmMHRY2P7TxXqNFesDPMo9VmR7nPzKlGNBPYiAla1WXv9X
        zGjrQoavYEb7h2Qn9fcDH6sE3/D7iYI3NCNRznEhWXlW67DHITMplgn9y2Bn
        KhSIec5auBR+oCbWub8u2ZrOgpadLK/69xv4Wo9uvt+x82hPGXHZQ8Jl/36p
        v/p5xm0WvGsXNImZ8mrJVoHsxHkLnbajPoWZ+qmSjCSdFfYjzuIwcV9yrz9n
        Zl7XSzNSCY+EQ7m2hqibWPMe4n1Y9zFW7mlgGj0lYUe+bDPQ7p0Lz3P73SQr
        /Tg7P5Cfp2Dywm4RBltucP/4LwN8o68W69Uo1qaMmXo6yOtAfi3yMePfCeI9
        1LcPeP7CyvuAHwkxVMhywibqHiDPcXZAZo4UfEsX1wD7nmICFKs7p2SlCpdq
        mUn5srB/HMT+M/Rdx8xmLWxijoTo8ZYk+BfOTL2v3cV3rDhZPtn19zPuKQmX
        3VS49PNKv5+ZcsMpP5xyw6lnLNkp1DOW/Ae4drnmBGdMPDPZkNTrj+7HhZnr
        q2o8Uv0pndNa87xozaneBti80sr9jL4AwOV+ZoxZ7ldk54FBl3h/ply6cf26
        BfMLuoWLc83o+4Su4e2eBpsTcqiSmeSQ3NMhjJl6xBLJOVF4H2Qq5Os1Vt4H
        4tUUu1nDTNicLXjZIPhHuIQuI/cnoJgl5cOeE7LyBGWms7jfgOjI4cO1+uyS
        jzYqIgl6T6LYz9RjxFvs58tPlk/Hg8vQwCTwkfo0yoRcd5p7gLO2Sew/0mWh
        f8GulHNGqJcB6ayklxMeqbZe7i2mhUcZi3QGyiSvO87EAye7Rgr+/mBG/xjh
        ErSPGW3zpx2Gue3vdV7gL5CV3bpFcsqXolgh8i6cHLLfDA9Jgq8Vco/q34mX
        wCH16CGeUs9trFEoM9XBYU2tib9iDgd8xVuYiFkzo18KMZkpgofAJfX7jBLX
        QDos+XuAy3NGVh6nzKQ8A9ynm9ibwQmxeUWWeo5L2ExgJj8Q9RNddbJ8Oh5c
        hgQmwbe4QKL5gqgOhXwH4CnkJc0jwfXKPZuimfnc0456XWaafXqieOwgjXU/
        X1r3k5YxTsO8fmbGOC9s8ZcEvSCee1yhB/r163efw5Cwt9Q5GcDkIFErZ2dX
        8HNYUC58dFSrSrPs6Iyl+I+OmXrmU4wHMitCfO4JK+7nO2bMZ0cvYMRxYG8i
        xwh+IJo9M0rwDjI9TPDKmRl1GDm/x6y/Oju3cGlJZlI800nwwjAbAb3yLGHT
        z6cCOSWkB1KMHfrOmpPlUwcumT/veV7QT/0uDv9aA5fgIeWg4JHyVIBNyE3g
        kuRlDTP6ZElfo/6GlFtJPlbCI82LPpa+elQ8aqw9yUx8NtqKfcxt+/pDPgGD
        kJnAI+Qk4hbAB2L0yPNBLsF2P6/Iuwb0SzoC/bUjT1w1B9TLo/iByPBE6KTU
        lwy8l/tNUizImZl65NNMRmDzpPUjQcgVRCwHfnjos6sEXwmXOFeh5yQx47lB
        NdBqf885JSuPQ2ZSPHOQ2JfkAzLMR4gOL51lqQeGr08FbD3ZRgOfkHPJyypG
        8+TM0VwXUcD9Mpp5aPFcXjx6Ji8tH3VMXA7oFfCJcubvGnhJ+GedcZm4Mysr
        d13duIn3ZeWPeT2hoOWT4tqFHxePnnWotX3xPeMaJsB2k30G0LXg54FciLOx
        GRQ1Y9aldWOb5l3T2HbplpYpq9a1TL5i0aqrdtRfceVG8J1iHCeNx3cOf9Jr
        w6Y9kVdds2fEug135F91za1pex94zk6sefPJ7mHk9vTt6QcZ85TAIggYRf4q
        cvDuE/scMRjIn6vs7IZd6eqU/oahdtVCPaqLa/nh2Mj8GmbEJfm86LzCfqDe
        6SBnZsSrl3hvqhWY5OLa94hrRtwa5zrO3BlirUYK/iWIPYnrory7PszUL+Ss
        y1G3AptyDhB8QLI+C95A/kVGhJR1mqFOPTB8zeUmZNHVo2sb+GWrVvL5s8dw
        /8BoruCIZ9Uu5+6NN3Nd9bV8VE2LJo/sB/r96WTvhb2GmORNalz265uq4HrC
        AZeo0m9tvRN+0nmFfD12VOmr+SVjX3EpW/6Xrmbj3+41G/6uGLf06QEDBsAf
        C7sSZy1sk7ix4xrTEwpa99qEjfgsbnjmF2kpaY9OnLxocnzxtAdcRq3/06ny
        yj/0zVft2HT9jbaC39SH4UKmOo+1MHnn3U+6ZZVNXtY3KPctXUjqJ4nZNbuH
        pLe9azdi4ZGhRYt/HJY/91B55dgHT3YPuzh7ES6fEPv5SfE3YRK5OJA9W5gx
        xw864QKFpsdE5m0JCq36PSBqFLcJLOWDYhu4bVwTD0tt5fEpE3l4RIOi15Ys
        Ztq6vZOES5KdBpvH2SnvTuWsswaXsLXvkK4buixszJkCl1UCl9hfOAegz9mL
        /XoJOwf9PUfBpZwDJOcaQH+zZSZ9FvYiZEi0gs3LNHvTmLBJ9RMbZs+Zy3/8
        /hn++ccb+Ni60TwqKppPa2/gIeO3cbcJe7hv3fV8dE1TJx716tULfgCc99DD
        dtn0Dv9cxqV/TB13SWn6eWhg5hEvv+j3aqur9iybEfXNytlR3zXW1z3oO2bD
        b25td3Bd6228qu2aZ5T9ArsSNmXClGkLRrqlN79ml9zw0XXXX/rFj9/f8vWO
        zYvgFxqYm5s/sLJ+8Vb3lt3cffyuv8Nq1x+8/a5H7SU8HlVGfv3ND+cXVC9c
        MDCu9vW+3ilv2XvFPDVrxszsZdPS7rp2VdvY0eNmzXIqWfq559itv4/QT/jr
        ZPewp18id3CK/fX888+HHfmYQo8yo+5KOXew1Si3B7kgSzIyspcXFtfd6RCQ
        cXiQe+wX3v7pPyTFFPHQ0NE8PKmRD82Yyp3KLufOFau5T+48npY6Afgm/w/5
        wSh2C3IkXDo7F+5BjWt4SJI1uHybmfIFce2Isy6TcAm/MWwQ+OgCxLXItuXF
        Mn/ONMZOocwk/z3ps3JM05NJc6QVbC7TwiX6YPh4lz8YERaPvIxNxcXFfN8z
        a/jr+6/lHh6e/KW7SviutfG8rmEUd510D9cp5N1wAy8tG6nmEWoiUJ+Ec393
        n4vCfiJchkTX8aG587lnZtuPfpHZn4SHR65rn9SyY/W8qF+uWhD98/SpE25o
        a2t9zHvC7r/dpj3E3Sfv5a1zN2OfDh/fMqkguOKyV4tqF75nmzT23V03rfzi
        95/v/PrB2ybBF9j39ScWrFkzL/a/BdVNv+sm3sXdWm/nUQ3XH7r/4RfsjoZH
        0JdffX9+RvWylUNyZ742MKry1dC4vE0HDrza+/5dU2vWzA//cO3C8A+fffTm
        yOUr1ocFVS37Nqdw9N8nu4e9QrP5Rb4juG9o1peKPgBsPiQweTcznmW4H+T8
        w695aVFx5TaP+FHv2gdmf6jzjd1XVl65MC01vWbSuNRPJ9Sk8MS4eu7pPZrH
        pE3gnrWbua5xJ/eqvZ6n5c/+ctysHZsXrr591sTJM2meuDyzzcXJecQe1Jyj
        pjkiJN0aXL7PjHVpuH7YmOslXMK+hM6DsxUxEj9xDbbMFB+5iJn8aucsLlXY
        VPvvezOppwEz6jBkRxj6eQb46Xdr4RI9Bl3dyt5T9ssOWvMZU9sMj7duaeaP
        70zjMydX8VHNy7nXzIcNFDP+Bu7i4irz6FXBI5yftyH+BlzqdAXcrWwV9ymY
        y1ctbf9uUuu4x52cnA3xrnmTRry6ZGH7HXFxCfMb6mtX6kdV/uw69QHuMeNh
        rpuy9895SzfNnjWz5aGn7xjzwYzpLYfD02pv33PLxoIXHl5z81OP3gb/T883
        Hir45NZrEn5bOSflkaiqS/e5KZ93m3o/j5lw88GDr7/fTwuPRDm1l87Xjb72
        5aqGhS/7JZbtOHDgtUvw/JsvP+C0ak74e1fMDnv5uSfvG/DAfbcOHFsV9k10
        7MnLlmEOobxHYDnvFVTxt3Nw/hF7e/tHxHrhHEP9M+Wmr8gp0N/kmzHhA/eE
        2vd9glOeSkhIgl5fPnLUyPor5kZ9vGJWzFcZqQ2/gnfhYXW8qLSIezXdzD2U
        M9Njwm18WMHiv9zyZ78zefplkFnko3WPDBse4ug44kHYqKg5h70aEpBlDS6R
        u0++KuAS+jf0WNiX8N3Bb4ccEPilfMW5AN8c1XTRPL1/Ei5lP6G63oRy9HBO
        Qp/pmJOJua1auETvnNio3G9ozcc31BoeL774Yp6Tk8PLi3K5l5c3T2i6jics
        foqHzH2cN7Uvk3n0GjPOqsJeuwO47Nc7gfvmLeBupat4dNEU/sEztUc+fn7M
        98sXNcG2umzF4ingKWIkiFuOnzO1bVtew1Vves58hAfPe5z7T9v715SZEz96
        ak/JR1fMyzm4eMky3Ee/+267Kee5Jx7CfV389N3z2+7ZVvvw/bdvSrzn3scG
        xbVsfz9IuTbfWY/yMQvuvs7SOk5fuHmEV8ON+92qN+xfNLfplc9fqHz48dva
        Z9Prd26eN1E6/7rNm1i8vHfvPic996fnxb58sHPc733DR/P+UbXcLariyODB
        Q3D//5H29BUFhRU7ggtnvx+Y1fJuc0PN4antDasUW4Ly2fKWzRm9a82K9vUF
        +WU1Hp76gyEhY/j4URE8O69cweStOM+4m6I3eKROoXkEBpszPDQxyt6+4E3q
        gU24DPLLsQaX36hwSfYl/LHIK0BOUrrYe+CXk4TLs65H7CnG5oXMVG9CuWLQ
        Z+XYiSEXiImcrOCA8m1auIyJyO9Y89aWsWY8GKUvNDxW17Xx3Cte5KM3vcoT
        5z+i6LNV9B7Yl7BVYePcDVyGpLRzz8r13CF/Ebf3L+QbLi//6IHbl7zS0twA
        W2TJ0kVToKvibIUPdtz8We3X1je0jPMbf+P3wH7hVft5cuv1PwRkjj30yD23
        AieG3ifPPXFn4qGXH/YS937hkw/dEfLso3fhPs+7etOt8WFzH/0jY+ULPGTO
        o39s3LkvS71+23Y95eM5dsvzXvXbXxqSOeXF5YsmLPvoqeK3vn1l5AP0nrVL
        R86Tcbnu0ukTrNi/HHmv/Xr6fesRkvXd4MR6PiRpPE/IHI38b5KVa7OzCzcF
        5k55P2TEnHe94qte3XF15fMv3FV2+OaNLehbYKh5WzCt4dqysnL4qQ114IF+
        FbcXZE3godFjuV1qG9e138t9ZzzIvafexxOGT0b/invCgjKqhg7News9IAy9
        smwLOvrT+3vnWoNLzOvDObxbhcsZEi7TJFxSb4K+7Byo6TpF2NSyNWVsQr/v
        yPeIDC1bcTRcTmkztx9LCzN5dlY2rx8ZzLNb1/LFj3/Jp931Ea+bt9MgU5kx
        v4xspns8vYt5dN3V3KvqGm6fNpW7uhbyebMaH2WiDkHRXZfNmNYKXkJWNqWm
        pjUtmtsG2ypr9qL1SxRM8XE3vMmbd73NRzcvvf7JB3fWMFE7e/dt230+PPxU
        rODp+fueuWvYa/vug7523r4ntg0vbZzzev32N3njzW/zjPkP//f9j7/vS+v2
        7gdfXRI9Yec9EZNuedG9dtPz6SOapr3zxms9X3p8+Yyfvn97kAmXlUsIl2PG
        jO3eVp8/31pc9u3p901AYPDTdnHV3L9kKnfJnqLIOb2h/0pycup1YflT3w4p
        XvCeV3Ld2yWlo5a0tTVNmdOW+uoNmxYtYKK2u7IkZdvMmdPgTzH0m8vJySvU
        l1Q/Zh8//s+hsY3cIW8O9535EE+59GkeMeVOxQatQTzlT/ThMZChj11+x8yI
        QN8ca+Zx/sJM9uU2gUvUPkBe1gtcIhYTwUzyknBJObH/VFzKtqacC0RxzaFi
        PeA392fSjOnI0NKVlnA5eeI4Mx7UVOUb8NdcX8rrxozkRWte5Ne8/hNv2/0u
        Hzd+Kt7zFjP6/B90cHC8P7h0Bdc3TVdweTUPCK/htkOy+KzpLcCsodY+NDR8
        0YiiYvj2DfURlZVVbevXXrq4X79+adMXrl2X1Hbjr3lXvsSvOfg9b9n66ps3
        bpoPeWnIpdyx7Tr711+6PY7u//GH7un76N7/hGJdGkf6vj662IvnzNrxx1UH
        f+RNN7/D29c8s5TWrWnhra2e0x9+fvTc9QeCR115wxdffnvh3Xt2Dnxq79ps
        eX13bJxbS2fejVs3BlQVxzx3svt3wAAbAy5d7EOPRIXY/ODp5/V3WOUsntm2
        nHvmT//dw8Pz+pjs+pfDy5e+l1TU+F5Bcd12saezJzePgCyl/Nfoouyw+/Pz
        CyjvKbk4P2lrfnHRhz45U990SGj+JSi+iYePWsMzVzzHJ/7nMM8Zt8GgF4G/
        wKRhJoXAZ4BvDmIx1tSSApfwW5E/FvFvxC+RUwl5WSrhEjaus8ClHLv8R+FS
        Q2aqc4HkPD0554ByPQy9VgP9im7QwmVTY62ynwZ0/Lt21AhznE5cxsfvOsxv
        fPcn3rTucT6preELBbePKa89VD9x6Ztpk7bxUaOKeEJhM3d3y+ee7iP4xAlj
        7xZ8g+yBz84gKxWqqa2taW5tba6bNWPyrJr6se/MmTL6S78pd/MqRV9ed+A7
        vmTd+ruZ6Yzt8OGp+fna/sf63rJlSXhN/aTtY7a9wTe9+ROvuvblP57a/5n3
        HXsP+kTMffTpxLl7n6udetm+m3Y+5I7PfPXl5+f/+vO3fSzpInOmNhZVl8ed
        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
        7CkPJMRk8yHCrzN0aJ7okZ3Jowpmc/eJd3O3KffzjPbt3N1zdIcPHvOXYiPz
        D1pxX18yY4xkhwYux6pw6cdMPX3+NbhUYVOrnrePh0dpC+oQaIYrYszIkYQ+
        a2ub9b6dnf1dtOb6ykqu11d08MDHx4+3tZh8tE5OzrymyhiTrmlbbcixQT6A
        1+R7/i4oGf1VRnK5GS6zUw24nCN41rLs0hmIl9RMnVS/Ql8U88yslphPR8zd
        +Y5/6fw35i1ciZpfJ3G29l+z7qY077rNytm/l3tMf/CH51467HDLppZm+d7v
        v/WyDPnfzz9yneHfW3c+FRJRt/pT3dQHuP/sR7n3zId/W736iqbPP3mjB713
        w+WN9fJnl0wvmrdqUX1xbanuQGG6y6t1bc2vBFQuel7Rw0+6jrj7Rb68by+/
        rxNigxDn3Z2aVf5cfOW0z0c31H01NHPqX+7Va7l3/XU8rrD293Hj6pD7m6hQ
        8Jxp9dPEWtitu3JxyvwpBfPqR/rflzE8+EBhcc4b0UWLv6YcAfASvlbMNER/
        bMxgysnIes2jcuUf7s27/9a13cmDo1uNuBQzO8MV299KXFIvzY0SLqFvyPKS
        8n0oD+9fhUsNbHbkHOjcyxsNMyMdyrgZNhW5aW+fd9jHOxAxsY6+3LNntPK9
        t0znM6ZP6OBD49hSM740j683PNbWT+M+47Yack2QQ1fVesWR9MRRf8q4zEwp
        5YH+EfDLYw5s87IlUyEvRyn41gclj3oxNKvu/cIpGw9PXnwd4ijezNxHcEnm
        qMW3eNXfyN2Usz9nyu6HtqydP+29A7ty/vj53QG477tuXlwqr8Pe3Sty33nv
        y0uCyhbtLqibsT+65vJfXRWZgesbv3T3Dx+8fbAjprlny+Rx8mfXXTYSNtd5
        wUlF1yckJ++vnzn75dFT1uN6PzjZPdz/Ej8DLqMjAl+zs7O7wy2s4NOo+Kgv
        MlMDPo+IT/h2aMly7lSzQcHm5r9r6xoha+KDg0NCx1REIlfD0B9/VEnyxtbG
        gg0j8qKfzk4NPpxVqP/TL6DGjJdknwwenPtTeHBqW0tjzaq8olHP+9Rd/z/3
        ltv+9k2aYYbLqLBCa3D5vsAlfFNUf7lQ8JjkJeUVEC4pb73TPL0zjZ3/J2x2
        +GdddfqxpLeosengMOLt4MAo6mtowGXN2Am8bf7lvGHWWr5nm6mua0pbvTku
        m+p5celIntW6nUc3Xvenrulm7tZ+L0eOamZO008yLnPS6nhoYPLGvLz8OfNm
        t167cOZY2CSjSvVNW50Sqj9wSW16Pyyt5jHF/iLdjXJDDLi8Yftuu4ARcz7y
        gtxsvZ1v3HDpGz++rv/r85eqX8E9b79qSq28Bg/ftnXkqBk7F3tOf3hfWPPm
        fXaZk/a5N95kyIl3n/4Qv3LrM0X03mtX1DfQ3/uefWLwFfOr5s1fcFW4Tezo
        54cMb3zBPX/KI+uuvQG+7MdOdg+jPhV6bFpK8vMRicVvecRXfRkRm/55fGL6
        y7Gxcev808d/4VC1lrvUbfm7ZsY22J/xNjaDgme2jbhJ4HJQ5Yio+5PyR+60
        jWv4ry5uNE/Om2Hk5TBzHWjYsIJ3ggIjEWdBnVyZi4vrmLiKhe83zLrh2dDg
        ka/JuIyNtCrfBzFrqou+VuBykcAl5GXpieDyn45NJvWmcXGrrDPNGidsGucu
        u7iWPBwRlgj/Agh83NmvXz/unzuJO5Yt564j1/GmlmkGHtjb2/PdN8znzc1G
        XVanc+dTmlL45Ysred3cTTy2ds1fflkTf/cat417NO/mHg3b/3Z2zDbg0kOX
        zwuzGnhoQPKGmzbPu/s/61O/PfTY+Lcbx09bG1J91Sdh9dd/5JI89u0JbVOR
        B0o98s1wCR6uWrMt0b5w/p/AZlL7jX+/8kD1z+89P+U23PM7L9w75s/fvzDo
        pke+++L8uW2lGxrnrTyUOPeBV2yzp7y09uqdQdljV93jPn4Xd514F3dvvf2r
        F/e/Oxjvf+aurQb78v4bcjbtXp/4yu7tK0d5ZLTcOyx3xkuDkupfqq6fg7Mf
        Nvq2k93DwwYZcZmUkPCSZ2zpN56pjV86+Gd8GhUdu1V5/bLS8jE3OhYu+sth
        9Ia/XMbd8Mfkmctb2ydNyrx61cLGnduWjrx62eiVBdnxLw6NqPzJNqqO20WP
        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
        yb53yYmxz/linIh5g8yH0PkvjfJMrO+ALVapJKnmKgqr3Fe5n2q8LoUL2TyI
        el3nbuj6TdwvfXCs72U/KVtPhuZppbF1BnUn8flqE+hP8lus3/UmyGXReSw6
        76xS1kbSkpSI/XqVKcxlYz5NufrJtZaUNdjoqh/t755Gr1PWUzHnBXyK65jn
        H+mzTYnL3zYO/biM4+3rS4/DTrtz70bZI545wWMmXPNi52+3mbD/Sp/drfON
        uKdyrjiOel+191Tm9JKr2HyInIjSZf3fxYFaTDiHSucz6vsnHrWOpL1fwF2X
        Y45Wi5hCXUCOpbmtHXdjbFznsDDvbKNaA8RsmZItJshdo2wW4d9wPV2uZIO6
        Pn1Sfh5EbV3mCYdf9mkT7o/B/hMal1pf2v72vB15BsbaH2PX+d0eVp80hWfG
        6D4adv4Le45rvwjjrewfgZ/sHcH+EXquOC/kQJUmvI+SqxC3FOKXYvOfKA6E
        axGH+C7uv6HeBibcn8bmrQV2/kpj6Z8Jm9SdWn9qHUqMkuvqfDO9BtabIPfM
        JXx/nRLmtOlr2vls5NSJzrbZGyLi8uCxzHWx+8KQxzK/UtcInZHYtBrfyHOC
        21vm7pB7w95h9xpifgB0KM9D4nNUmsCnzDyki0yQZ6T7Ltl9Xlz7KvdTzVPm
        40JREsWDuBdzD9Z7Ce9dx4pte/+sxeQi8EkdqvUocarnn2uAcnGE6M/ov+V6
        0rltvD5tXub2Jly9N1FDagr7wugeYvWyZux6hGWv3bPG9NxkdfYNEX7Yerk3
        1tvp+gpbV7Ina8qEe4fq/oTnKeGeqm0Tzo2NVa1XabO4uJDNh2xxcSFeU+cz
        FsvdiMTjcs3PWhD7uU1h7gox6sIp51+vgYXKq5Wcp677L9T/mXNKf7Gf21tT
        k/ljByd8Xq1p1nJB17AHBWvdz3gPCmssz01UZf+bfb94BpPDGe7Nrtnn/bOu
        QtcFV5qgLpj1+npNR82VniN7X7XxavMgzYWi+JDNhWwepM8zni93w+kPX2ms
        rBJ82npU49SFV50bWkxeZQmvEXVdnaPk96ft757e6+KFyD2w1jXsNbtnk10r
        pG3M5cblOTW1mS4n7+6cQf6O5rDsb6PvX+8r7F+u87fts8zPsX7a3MeFVXtP
        dXEhzYcWw4UuUNfXe7Adn3pF68glYFTrU41XLTZ2bXH9zTmO6+o8NmKT9d3l
        yVT2+458tn+r1jXWuLYxWVtbaQpzoJetLsEap3NdMUvUk5gch4XfUveM0Lnm
        mofrcynK1L27fFe/ZaLnKWrvs/dOmwuVyoc0D3Lh0JW7EeNx+TA6H25dUsp1
        dD0MfU9+fXdr89xbnXqoZ9sRtbajfJn6rN1l5bLq/s9pb5lb793TP9j32NyU
        uVfuhRwWel1zWDs2ojks9pSo3jbzjW/UnuraP11caD4+NB8PKqoXlzr2r2Qp
        AVeLFnt9m3AfMb9X9fDAeLPL/9PeOvtpU2ijYX0jFsE4JrjsGenXrdb/uVWp
        7O2uvWOofwz30mJKO7ub923vJ0Xvu4Rxn48HLZYLLWhPXul1/UqQ5Rh7C5fn
        mnBP3Ej/D3JoPZ05ZYpz2SoT9Dpctj4x1hqHv+cnEf4ecFJwU9fZ3VrP6z4a
        2me16HsuAadnkgvFWFzDote4CffEBZfdPNR/eNilj/o6Zx4w4Ro/6p6ovmrL
        0pNC369n/17jrB3pnUbuA/QfOewuE81hGbOcj8MuOSdiiXiNcfgKEjWnLi4L
        HlqBnBl77Sers39r6cxithrP3F2SzlT3GqkrcY6ECfw9uBfoccYsbQ6rzzuM
        2kdWpAYxxtwrW+y1bgIum++72dMxk3XppY7W2c+Y4ryQ/h87lrkonan3kChd
        KbERxiyZp651Je+VvJvndtMPuyZ6jMdy9ouFS9b95/2y0HupdPab89iZjGWy
        n5OukY7q3b0Ym+23onywoiuhpxmz1OfKRfmp7Pukbl/2mGsssZQipjiXRXy9
        HHrIpZ9am+e+aMLxQPbDs+OBi+59aMJ87hxvP7jbqSs7Zo6KnmZN13bj9veQ
        w1KvI+c7Kq4T4zKWFRETzWVDZ5W5dKbkAN2kdBGwoGMPOpd9Uf0P9f3V1mZS
        rntQuhL6r8eE+2jb/h72PaEdDN/xGcuDiCWWxYiljxjvJpflWWUVQ/0T212Y
        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
        SyyxxBJLLLHEEkssscQSSyyxxFIo/x/outWB
        "], {{0, 230}, {230, 0}}, {0, 
                       255},
        ColorFunction->RGBColor],
        BoxForm`ImageTag["Byte", ColorSpace -> "RGB", Interleaving -> True],
        Selectable->False],
        BaseStyle->"ImageGraphics",
        ImageSizeRaw->{230, 230},
        PlotRange->{{0, 230}, {0, 230}}]\)}, 1]]], Alignment -> Center, 2]},
            
            {Insert[
              Grid[Partition[
                Join[{Insert[
                   Grid[{{Style["MetaData:", 15, Underlined, Bold], ""}}], 
                   Alignment -> Left, 2]},
                 {Framed[
                   Insert[
                    Insert[
                     Grid[
                      Transpose[
                       Partition[{
                         Insert[
                          (*these are 4 functions that are output together. 
                          The 4 functions are: mzML file version, 
                          discription, software, and spectrum list.*)
            
                                        Grid[
                          Partition[
                          
                          Join[{MzmlVersion[
                          overall]}, {MzmlFileDescriptionFileContent[
                          overall]}, {MzmlSoftware[
                          overall]}, {MzmlSpectrumList[overall]}], 1], 
                          Dividers -> Center], Alignment -> Left, 2],
                         
                         Insert[
                          (*these are 2 other functions that are output \
        together. The 2 functions are: 
                          mzML file data processing and instrument \
        configuration.*)
                          Grid[
                          Partition[
                          
                          Join[{MzmlDataProcessing[
                          overall]}, {MzmlInstrumentConfiguration[
                          overall]}], 1], Dividers -> Center], 
                          Alignment -> Left, 2]}, 1]],
                      Dividers -> Center], Alignment -> Top, 
                     2], {Background -> {None, {{Lighter[
                          Blend[{Lighter[LightBlue, 0.1], 
                          Lighter[LightGray, 0.9]}], .1]}}}, 
                     Spacings -> {2, {2, {0.7}, 0.9}}}, 2]]}], 1]], 
              Alignment -> Left, 2]}], 1]];
        Return[metaDataOverall]
    ]


(* ::Subsubsection:: *)
(* MZML Back-End *)

(* ::Function:: *)
(* f:MzmlDictionaryAccessionNames *)
(*dictionary takes one "string" input and that is the file name \
constant*)
(***Function***)
MzmlDictionaryAccessionNames[mzmlControlledVocabulary_] :=
    Module[ {
     inputString = mzmlControlledVocabulary,
     msDictionary},
        inputString;
        (*find keys "accessions (id):" and "name:" but extract their \
      corresponding values in a list and map them to their corresponding [
        Term] positions, 
        once the "accessions (id):" and its corresponding "name:" in a \
      list, make an association between them so you can serch this \
      association using the key "accessions (id):" to get the value "name:"*)
        msDictionary = 
        Association[#[[1]] -> #[[
            2]] & /@ (StringTrim[inputString[[#]], 
             "id: " | "name: "] & /@ ({#[[1]] + 1, #[[1]] + 2} & /@ 
             Position[inputString, "[Term]"]))];
        Return[msDictionary]
    ]

(* ::Function:: *)
(* f:MzmlDictionaryMSUOLabeling *)
(***Function***)
MzmlDictionaryMSUOLabeling[mzmlControlledVocabulary_] :=
    Module[ {
      inputString = mzmlControlledVocabulary,
      msDictionary},
     (*Cases will first find all the "has_units". 
     StringCases will find all Accession Numbers with "MS:" or "UO:" and \
   assign them to "x" and then we map to Cases. 
     Union will get ready of all the repeats for each Accession Number. 
     Finaly, make an association between each Accession Number and its \
   definition*)
        msDictionary = 
         Association[
          Flatten[DeleteCases[
            Union[(StringCases[#, 
                 x : ("MS:" | "UO:" ~~ NumberString) ~~ " ! " ~~ 
                   y : __ -> {x -> y}] & /@ (Cases[inputString, 
                 x_ /; StringCases[x, "has_units"] != {}]))], {}], 1]];
        Return[msDictionary]
    ]

(* ::Function:: *)
(* f:MzmlTimeAndMassLabeling *)
(***Function***)
MzmlTimeAndMassLabeling[mzMLFilePath_, mzMLdictionaryUO_] :=
    Module[ {inputFile = mzMLFilePath, emptyList, filePath, allMetaData, 
      testXML, scanElement, unitAccession, unitAccessionValues, 
      labelingUO = mzMLdictionaryUO},
        emptyList = {};
        filePath = inputFile;
        (*open and read your file*)
        inputFile = OpenRead[filePath];
        
        (*Set the data streaming position to 0 to stream from the begining \
      of the file*)
        SetStreamPosition[filePath, 0];
        While[True,
         allMetaData = Read[filePath, String];
         
         (*Using <StringCases> to search "scan" and <..> 
         is to keep keep representing "scan" until it matches and while it \
      doesn't match keep appending <AppendTo> into the emptylist. Using <
         If> and <Break[]> to stop onc "scan" matchs.*)
         If[ StringCases[allMetaData, WordCharacter ..][[1]] == 
           "binaryDataArrayList",
             Break[],
             AppendTo[emptyList, allMetaData]
         ]];
        
        (*</mzML> as closing tags for the "XML" format.*)
        emptyList = 
         AppendTo[emptyList, 
          "</spectrum></spectrumList></run></mzML></indexedmzML>"];
        
        (*Import as "XML"*)
        testXML = ImportString[StringJoin[emptyList], "XML"];
        scanElement = If[ Cases[testXML,
            XMLElement["scan", _, _], Infinity] == {},
                          {"Unknown"},
                          Cases[testXML,
                           XMLElement["scan", _, _], Infinity]
                      ];
        unitAccession = {If[ Flatten[Cases[#, 
                  x_ /; Length[x] > 1 && x[[1]] == "unitAccession" :> 
                   x[[2]], All] & /@ #] == {},
                             {"Unknown"},
                             First[Flatten[
                               Cases[#, 
                                  x_ /; Length[x] > 1 && x[[1]] == "unitAccession" :> 
                                   x[[2]], All] & /@ #]]
                         ] &@scanElement};
        unitAccessionValues = Values[labelingUO[[#]] &@unitAccession][[1]];
        Close[inputFile];
        Return[unitAccessionValues]
    ]

(* ::Function:: *)
(* f:MzmlSpectrumBinaryDecoder *)
(***Function***)
MzmlSpectrumBinaryDecoder[string_, inputAccuracy_, form_, endian_,compression_] :=
    Module[ {
    (*1st input, the user selected spectrum binary data*)   
    inputString = string,
    (*2nd input, spectrum input bit accuracy*)
    
    binaryFormat = inputAccuracy,
    (*3rd input, data enocoding formate (Base64)*)
    inFormat = form,
    (*4th input, 
    interpret the bytes making up the data word as little endian \
    correspond to (-1)*)
    inputEndian = endian,
    (*5th input, binary data is compressed or not*)
    
    compressionStatus = compression,
    result, str, returnListNoCompression, returnListZLIPCompression, 
    characterList, intensities},
        result =
         If[(*check the spectrum compression status*) compressionStatus == "no compression",
          
          (*imports binary data (inputString) in Base64 format from a \
      string, and StringToStream is to open a stream for reading from a \
      string.*)
             str = StringToStream[ImportString[inputString , inFormat]];
             
             (*str: is the imported spectrum binary data. binaryFormat: 
             spectrum input bit accuracy. 
             ByteOrdering: 
             is an option for BinaryReadList that specifies the ordering of \
         bytes should be assumed for your computer system*)
             returnListNoCompression = 
              BinaryReadList[str, binaryFormat, ByteOrdering -> inputEndian],
             returnListZLIPCompression = 
              BinaryReadList[
               StringToStream[
                ImportString[inputString , {"Base64", "String"}]], "Byte", 
               ByteOrdering -> -1];
             characterList = 
              FromCharacterCode[
               Developer`RawUncompress[returnListZLIPCompression]];
             intensities = ImportString[characterList, binaryFormat]
         ];
        Return[result]
    ]

(* ::Function:: *)
(* f:MzmlOffSetListCounter *)
(***Function***)
MzmlOffSetListCounter[mzMLfilePath_] :=
    Module[ {inputFile = mzMLfilePath, filePath, streamingFile, 
      spectrumIndex, msLevel, massInfo, scanStartTime, preCursor, 
      selectedIon, chromatogramIndex, emptyList, 
      foundTagsAndAttributesOfIntrest, getThePosition, processedInfo, 
      masterList, offsetReturned},
     
     (*Make an empty list to ppend to*)
        emptyList = {};
        (*setting the input directory to the file path*)
        filePath = inputFile;
        (*OpenRead will open the file that its path is set to filePath, 
        and returns an input stream object*)
        streamingFile = OpenRead[filePath];
        
        (*read the entire file and look for an '"<spectrum index=", 
        "ms level", "base peak m/z", "scan start time", "scan start time", 
        "<precursor ", "selected ion m/z", 
        "<chromatogram index"' in every scan and "<chromatogram index" in \
      the file.*)
        spectrumIndex = "<spectrum index=";
        msLevel = "ms level";
        massInfo = "base peak m/z";
        scanStartTime = "scan start time";
        preCursor = "<precursor ";
        selectedIon = "selected ion m/z";
        chromatogramIndex = "<chromatogram index";
        
        (*Sets the current point in the open stream to the integer 1 "1 \
      correspond to the beginning of the stream - 1st character in our \
      file."*)
        SetStreamPosition[filePath, 1];
        
         (*While test=True, evaluate the body*)
        While[True,
         (*Find takes two input elements, 1st: input stream, 2nd: tags/
         attributes of intrest in file - givin in a list - 
         and IgnoreCase set to True as an option.*)
         foundTagsAndAttributesOfIntrest = 
          Find[mzMLfilePath, {spectrumIndex, msLevel, massInfo, 
            scanStartTime, preCursor, selectedIon, chromatogramIndex}, 
           IgnoreCase -> True];
         (*StreamPosition[
         mzMLfilePath] to return an integer that specifies the position of \
      the found tags/attributes of intrest in the open stream.*)
         getThePosition = StreamPosition[mzMLfilePath];
         (*If statement used to indicate one end of file is reached, 
         then break, 
         otherwise append to emptyList a list with both found found tags/
         attributes of intrest with thier position in file*)
         If[ foundTagsAndAttributesOfIntrest == EndOfFile,
             Break[],
             Null(**),
             AppendTo[
             emptyList, {foundTagsAndAttributesOfIntrest, getThePosition}]
         ]];
        
        (*the processedInfo function uses pattarn matching to match the \
      fallowing for each spactram in a given file:
        ("Index", 
        "Scan"), ("MS Level"), ("BasePeakMass"), ("Scan Start Time"), \
      ("PrecursorScan"), ("SelectedIonMass") and ("Chromatogram index") \
      into the extracted elements from the while loop that was appended to \
      the variable emptyList.*)
        processedInfo = MapIndexed[Flatten[
            {#1[[2]], (StringCases[#1[[1]],
               {___ ~~ "index=\"" ~~ x : NumberString ~~ ___ ~~ "id=\"" ~~ 
                  b : Shortest[(__ ~~ ___ ~~ "\"")] ~~ ___ -> {"index", x, 
                  "id", "\"" ~~ b},
                ___ ~~ "ms level\" value=\"" ~~ c : NumberString ~~ ___ ~~ 
                  "\"" ~~ ___ -> {"ms level value", c},
                ___ ~~ "base peak m/z\" value=\"" ~~ 
                  d : NumberString ~~ ___ ~~ 
                  "\"" ~~ ___ -> {"base peak m/z value", d}, ___ ~~ 
                  "scan start time\" value=\"" ~~ e : NumberString ~~ ___ ~~
                   "\"" ~~ ___ -> {"scan start time value", e},
                ___ ~~ "precursor " ~~ 
                  f : Longest[(__ ~~ ___ ~~ "\"")] ~~ ___ -> {"precursor", 
                  f}, ___ ~~ "selected ion m/z\" value=\"" ~~ 
                  g : NumberString ~~ ___ ~~ 
                  "\"" -> {"selected ion m/z value", g},
                ___ ~~ "chromatogram index=\"" ~~ h : NumberString ~~ ___ ~~
                   "\"" ~~ ___ -> {"chromatogram index", h}}]), #2}] &, 
          emptyList];
        
        (*
        Position[processedInfo,"index"]: 
        Find each scan index position in every list of lists in the \
      processedInfo,
        {#[[1]]-1,#[[1]]}: 
        subtract 1 to make the scan index position that is one less and put \
      it in a list with its original position,
        processedInfo[[#[[-1]];;]]: 
        to account for the total ion chromatogram, 
        we take it using [[-1]] since TIC is at the end of the list ,
        Partition[#[[2;;-2]],2]: 
        Partition to make a list of lists that each one of those sublists \
      has 2 elements, the 1st one correspond to the scan index position, 
        and the 2nd one is the one element that is befor the next scan index,
        Take[processedInfo,#]: 
        maping all the partitioned list of lists into processedInfo.
        *)
        masterList = (Join[((Take[
                  processedInfo, #]) & /@ (Partition[#[[2 ;; -2]], 
                 2])), {processedInfo[[#[[-1]] ;;]]}] &@
           Flatten[(({#[[1]] - 1, #[[1]]} &) /@ 
              Position[processedInfo, "index"])]);
        offsetReturned = 
         DeleteCases[ masterList, x_ /; MemberQ[{2, 3, 5, 7}, Length[x]]];
        Return[offsetReturned]
    ];

(* ::Function:: *)
(* f:MzmlSpectrumSelector *)
(***Function***)
MzmlSpectrumSelector[fileNameI_, msLevel_, preCursor_, timeMin_, timeMax_, massMin_, massMax_] :=
    Module[ {spectraMSLevel = msLevel,
      spectramPRECursor = preCursor,
      basePeakOrPrecursorMass,
      (*N[expr]gives a machine\[Hyphen]precision number,
      so long as its magnitude is between $MinMachineNumber and \
    $MaxMachineNumber. $MinMachineNumber is the smallest hardware \
    floating-point number (10^-308), 
      and $MaxMachineNumber is the largest hardware floating point \
    number (10^308).*)
      scanTimeMin = N[timeMin],
      scanTimeMax = N[timeMax],
      scanLocator,
      spectramassMin = N[massMin],
      spectramassMax = N[massMax],
      inputFile = fileNameI, offsetLimited, specMin},
     
     (*Searching by Spectrum MS Level*)
     (*"spectraMSLevel" is the \
   spectrum ms level, it can be 0,1,2,3,.... 
     "spectramPRECursor" is the spectrum pre-
     cursor meaning that is the parent spectrum for a selected spectrum. 
     "spectramPRECursor" should always be set to "No precursor-ion" \
   unless spectraMSLevel is set to 2 or higher.
     
     Using Which statment to evalute the the user inputs. 
     Spectrum with ms level -user input - less than 2, 
     "No precursor-ion" is used, 
     since any ms level less than 2 has no parent spectrum. 
     And we search for them by base peak mass - inputed by the user-  
     hint: "base peak m/z value". Spectrum with ms level - user input - 
     more or equal to 2, "Yes precursor-ion" is used, instead. 
     Any ms level greater than 2 will have a parent spectrum that it \
   came from 'fragmentation'. We then search by selected ion mass - 
     user input - correspond to the "selected ion m/z value".*)
        Which[
         (spectraMSLevel == 0 && 
            spectramPRECursor == "No precursor-ion") || (spectraMSLevel == 
             1 && spectramPRECursor == 
             "No precursor-ion") || (spectraMSLevel >= 2 && 
            spectramPRECursor == "No precursor-ion"), 
         basePeakOrPrecursorMass = "base peak m/z value",
         (spectraMSLevel >= 2 && spectramPRECursor == "Yes precursor-ion"), 
         basePeakOrPrecursorMass = "selected ion m/z value"];
        
        
        (*If user sets minimal spectrum mass to 0, 
        replace "spectramassMin" to specMin. specMin=10^-300*)
        specMin = 10^-300;
        (*Search by scan time only given the MS level preference*)
        If[ spectramassMin == 0.,
            spectramassMin = specMin
        ];
        Which[
         (*1st searching by scan time: 
         user input for time range while the mass range is unchanging, 
         search by unchanged mass followed by user time input.*)
         
         (*scan 'minimal' time not equal to 0 and scan 'maximal' time not \
        equal to 1000000)*)(! (scanTimeMin == 0. && 
              scanTimeMax == 
               1000000.)) &&
          (*scan 'minimal' mass equal to specMin - 
          specMin=10^-300 - 
          and scan 'maximal' mass equal to 1000000*)
          (spectramassMin == 
             specMin && spectramassMax == 1000000.),
         
         offsetLimited =
          ((*Mapping the "base peak m/z value" positions into the list of \
        offsets*)
           
           inputFile[[#]] & /@
            (*find the position of all "base peak \
        m/z value". *)((Position[inputFile, 
                x_ /; Length[x] > 1 && x[[2]] == basePeakOrPrecursorMass, 2,
                 Heads -> False])[[All, 1]]));
         scanLocator =
          ((*Mapping the "scan start time value" (these are the selected \
         scans based on user scan time input) positions into the list of \
         offsets*)
           
           offsetLimited[[#]] & /@
            (*find the position of all "scan \
         start time value" given the user scanTimeMin and scanTimeMax input.*)
         \
            ((Position[offsetLimited, 
                x_ /; Length[x] > 1 && 
                  x[[2]] == 
                   "scan start time value" && (scanTimeMin <= 
                    ToExpression[x[[3]]] <= scanTimeMax), 2, 
                Heads -> False])[[All, 1]])),
         
         
         (*2nd searching by scan mass: 
         user input for mass range while the time range is unchanging, 
         search by unchanged time followed by user mass input.*)
         
         (*scan 'minimal' time is set equal to 0 and scan 'maximal' time is \
        set equal to 1000000)*)
         (scanTimeMin == 0. && 
            scanTimeMax == 
             1000000.) &&
          (*scan 'minimal' mass not equal to specMin - 
          specMin=10^-300 - 
          and scan 'maximal' mass not equal to 1000000*)
          (! \
        (spectramassMin == specMin && spectramassMax == 1000000.)),
         
         scanLocator =
          (*Mapping the basePeakOrPrecursorMass (these are the selected \
        scans based on user mass input) positions into the list of offsets*)
        
             inputFile[[#]] & /@
           (*find the position of all \
        basePeakOrPrecursorMass given the user spectramassMin and \
        spectramassMax input.*)(Position[inputFile, 
              x_ /; Length[x] > 1 && x[[2]] == basePeakOrPrecursorMass && 
                spectramassMin <= ToExpression[x[[3]]] <= spectramassMax, 2,
               Heads -> False][[All, 1]]),
         
         (*3rd searching by both scan mass and time: 
         user input mass range and time range.*)
         
         True,
         
         (*searching by a mass range*)
         offsetLimited = 
          inputFile[[#]] & /@ 
           Position[inputFile, 
             x_ /; 
              Length[x] > 1 && 
               x[[2]] == 
                basePeakOrPrecursorMass && (spectramassMin <= 
                 ToExpression[x[[3]]] <= spectramassMax), 2, 
             Heads -> False][[All, 1]];
         (*searching the outcome of the mass range search, 
         searching by a time range*)
         scanLocator = 
          offsetLimited[[#]] & /@ 
           Position[offsetLimited, 
             x_ /; Length[x] > 1 && 
               x[[2]] == 
                "scan start time value" && (scanTimeMin <= 
                 ToExpression[x[[3]]] <= scanTimeMax), 2, Heads -> False][[
            All, 1]]];
        If[ spectraMSLevel == 0,
            Null,
            scanLocator =
             (*Mapping the MS Level (these are the selected scans based on \
            user MS Level input) positions into the scanLocator*)
             
             scanLocator[[#]] & /@
              (*find the position of all \
            spectraMSLevel given the user scan ms level input.*)
              (Position[
                 scanLocator, 
                 x_ /; Length[x] > 1 && x[[2]] == "ms level value" && 
                   ToExpression[x[[3]]] == spectraMSLevel, 2, Heads -> False][[
                All, 1]])
        ];
        Return[scanLocator]
    ];

(* ::Function:: *)
(* f:MzmlGetSpectraMSLevel *)
(***Function***)
MzmlGetSpectraMSLevel[spectramIndexList_] :=
    Module[ {outPutMSLevel,
      msLevel = spectramIndexList},
     (*find "ms level value" in a selected scan using Cases*)
        outPutMSLevel = 
         Cases[msLevel, x_ /; x[[2]] == "ms level value" :> x[[3]], 1];
        (*returing the ms level, e.g. {"1"}, {"2"},...*)
        Return[outPutMSLevel]
    ]

(* ::Function:: *)
(* f:MzmlSpectrumMassTimeAssociation *)
(***Function***)
MzmlSpectrumMassTimeAssociation[listOfselectedSpectrums_, 
  selectedIonMassOrBasePeakMass_] :=
    Module[ {selectInput,
      inputSelectedSpectrums = listOfselectedSpectrums, 
      TrueOrFalseForInput = selectedIonMassOrBasePeakMass},
        selectInput =
         Association @@
          ({(*#[[2]]->#[[
              1]]} is reversed becasue manipulate uses the reversed form*)
        \
             #[[2]] -> 
               Flatten[ToExpression[#] & /@ #[[
                  1]]]} & /@
            (*Transpose the above three together.*)
        
                 Transpose[
             {
              {(*Output either "selected ion m/z value" or "base peak m/z \
        value" masses(x[[3]]) if True or False respectively.*)
                 
                 Cases[#, x_ /; x[[2]] ==
                     (*User will be able to navigate the spectra by \
        searching using "selected ion m/z value" or "base peak m/z value". 
                     If click "selected ion m/z value" in viewer interface, 
                     then is True. If click "base peak m/z value", 
                     then that is false.*)
                     
                     
                     If[ TrueOrFalseForInput == True,
                         "selected ion m/z value",
                         "base peak m/z value"
                     ] :> 
                   x[[3]], 1],(*Output the "scan start time value"*)
             
                     Cases[#, 
                  x_ /; x[[2]] == "scan start time value" :> x[[3]], 1],
                 (*Output the "ms level value"*)
                 
                 Cases[#, x_ /; x[[2]] == "ms level value" :> x[[3]], 
                  1]} & /@ inputSelectedSpectrums, 
              inputSelectedSpectrums}])
    ];
        
(* ::Function:: *)
(* f:MzmlGetChromatogramTIC *)
(***Function***)
MzmlGetChromatogramTIC[offsetListForChromatogram_] :=
    Module[ {outPut, chromatogram = offsetListForChromatogram},
        outPut = {(*take the last element in that last. Indexing [[1,1]]*)
        
             Last[
            ((*find all id's in the off set list*)
             Cases[#, x_ /; Length[x] > 1 && x[[4]] == "id", 2] & /@ 
              chromatogram)][[1, 1]]};
        Return[outPut]
    ]

(* ::Function:: *)
(* f:MzmlTotalIonChromatogram *)
(***Function***)
MzmlTotalIonChromatogram[fileName_, selectedChromatogram_] :=
    Module[ {streamingFile,
      inputFile = fileName,
      ChromatogramPosition = selectedChromatogram[[1]],
      whatIReadIn,
      emptyList,
      joiningListToString,
      binaryDATAArray},
     (*OpenRead will open the file that its path set to filePath, 
     and return an input stream object*)
        streamingFile = OpenRead[inputFile];
        (*Setting the streaming position to the chromatogram position \
      outputed by MzmlGetChromatogramTIC[mzMLoffsetList]*)
        SetStreamPosition[inputFile, ChromatogramPosition];
        (*Clossing tag that must be use in order to parse the chromatogram \
      in mathematica*)
        whatIReadIn = "</chromatogram>";
        (*Opening tag that must be use in order to parse the chromatogram \
      in mathematica*)
        emptyList = {"<chromatogram>"};
        
        (*While test=True, evaluate the body*)
        While[True,
         (*Read takes two inputs, 1st: input stream, 2nd: 
         specified type of intrest, 
         Record represent a sequence of characters delimited by record \
        separators.*)
         whatIReadIn = Read[streamingFile, Record];
         (*Append every line to the emptylist*)
         AppendTo[emptyList, whatIReadIn];
         (*If statement used so when chromatogram is reached, then break, 
         otherwise append to emptyList a list with both found found tags/
         attributes of intrest with thier position in file*)
         If[ StringCases[whatIReadIn, WordCharacter ..][[1]] == 
           "chromatogram",
             Break[]
         ]];
        
        (*trims whitespace from the beginning and end of the string, 
        join and then import the data in the "XML" format from a string.*)
        joiningListToString = 
        ImportString[StringJoin[StringTrim[#] & /@ emptyList], "XML"];
       (*
       Constant binaryDATAArray serve as a storage for element <
       binaryDataArray> in a selected scan. Element <binaryDataArray> 
       is where a binary data arrays are located for a given spectrum \
     scan. Data point arrays for default data arrays (m/z, intensity).
       *)
        binaryDATAArray = 
         Cases[joiningListToString, XMLElement["binaryDataArray", _, _], 
          Infinity];
        Return[binaryDATAArray]
    ]

(* ::Function:: *)
(* f:MzmlGetPrecursorSpectraAtMSLevel *)
(***Function***)
MzmlGetPrecursorSpectraAtMSLevel[mzMLOffSetList_, 
  selectedSpectraOfIntrestII_] :=
    Module[ {(*offSetList is equal to the 1st input to this function and that \
    is "mzMLoffsetList"*)
      offSetList = mzMLOffSetList,
      selectedScanPosition,
      (*findMS1 is equal to the 2nd input to this function and that is \
    the user selected spectrum*)
      findMS1 = selectedSpectraOfIntrestII,
      positionsReturned,
      positionsInListOfLists,
      AllScans},
     
     (*add user selected spectum to a list*)
        positionsReturned = {findMS1};
        While[True,
         selectedScanPosition =
          (*using StringCases and pattern object that matches the longest \
        sequence consistent with the pattern 
          "\""~~__~~
          "\"" and this will be used as a referance to search the of 1st \
        input "mzMLoffsetList" to find the parent spectrum*)
          Flatten[StringCases[#, ___ ~~ "=" ~~ 
                 pre : Longest[
                   "\"" ~~ __ ~~ "\""] ~~ ___ -> {"precursor scan number", 
                 pre}]][[2]] &@
           (*use Cases to find the "precursor" tag in the user selected \
        spectrum.*)(Cases[#, x_ /; x[[2]] == "precursor" :> x[[3]], 1][[1]] &@
             findMS1);
         AppendTo[positionsReturned, selectedScanPosition];
         (*find the position of selectedScanPosition in the offSetList*)
         positionsInListOfLists = 
          Position[offSetList, selectedScanPosition][[1, 1]];
         (*find the parent spectrum for a given positionsInListOfLists*)
         findMS1 = offSetList[[positionsInListOfLists]];
         (*if no "precursor" tag found then break*)
         If[(*using Cases to find the "precursor" tag in the parent \
      spectrum*) Cases[findMS1, x_ /; x[[2]] == "precursor" :> x[[3]], 1] == {},
             Break[]
         ]];
        AllScans = 
         offSetList[[#]] & /@(*get the Position of both the user selected \
        spectrum and the parent spectrum of the user selected spectrum, 
          and you do that by maping positionsReturned into the \
        offSetList*)((Position[offSetList, #][[1, 1]]) & /@ positionsReturned);
        Return[AllScans]
    ]

(* ::Function:: *)
(* f:MzmlGetMSSpectraLocationAndBinary *)
(***Function***)
MzmlGetMSSpectraLocationAndBinary[fileName_, 
  selectedSpectraOfIntrest_] :=
    Module[ {streamingFile,
      (*inputFile is equal to the 1st input to this function and that is \
    "mzMLfilePath"*)
      inputFile = fileName,
      (*findMS is equal to the 2st input to this function and that is \
    the user selected spectrum*)
      findMS = selectedSpectraOfIntrest[[1]],
      whatIReadIn, emptyList, joiningListToString},
     
     (*open an mzML file that is set equal to <mzMLfilePath> 
     and read data from it,and return an input stream object.*)
        streamingFile = OpenRead[inputFile];
        (*set the steaming postion to the location of the selected spectrum \
      index location that was past to this function as a second input*)
        SetStreamPosition[inputFile, findMS];
        whatIReadIn = "</spectrum>";
        (*list with the spectrum closing tag <spectrum> 
        in order to convert the selected spectrum to mathematica's XML \
      format*)
        emptyList = {"<spectrum>"};
        
        (*While test=True, evaluate the body*)
        While[True, whatIReadIn = Read[streamingFile, Record];
                    AppendTo[emptyList, whatIReadIn];
                    (*if statment to break the loop once spectrum is found*)
                    If[ StringCases[whatIReadIn, WordCharacter ..][[1]] == "spectrum",
                        Break[]
                    ]];
        
        (*import the scan to mathematica XML format*)
        joiningListToString = 
         ImportString[StringJoin[StringTrim[#] & /@ emptyList], "XML"];
        Return[joiningListToString]
    ];
  
(* ::Function:: *)
(* f:MzmlSpectrumPlotLabels *)
(***Function***)
MzmlSpectrumPlotLabels[fileName_, selectedSpectraOfIntrest_, 
  dictionaryFunction_] :=
    Module[ {
      (*inputFile is equal to the 1st input to this function and that is \
    "mzMLfilePath"*)
      inputFile = fileName,
      (*findMS is equal to the 2st input to this function and that is \
    the user selected spectrum*)
      findMS = selectedSpectraOfIntrest,
      entireSelectedScan, allCvParamsAccessions, AccessionsOfInterest, 
      matchedAccessions,
      (*mzMLdictionaryFunction is equal to the 3rd input and that is the \
    constant set equal to the dictionary function \
    MzmlDictionaryAccessionNames*)
      
      mzMLdictionaryFunction = dictionaryFunction,
      mappedAccessionsToDictionary, addingTagsToMatchedAccessions, 
      tagingStuffFromselectedScan, finalOutPut},
     
     (*entireSelectedScan is equal to MzmlGetMSSpectraLocationAndBinary \
   function that takes two inputs. 
     The first input is the constant that is set to the file path \
   (inputFile). 
     The second input is the list where the spectrum index location is \
   located (findMS). 
     The function will output the entier scan includ all cvParams and \
   binary information.*)
        entireSelectedScan = 
         MzmlGetMSSpectraLocationAndBinary[inputFile, findMS[[1]]];
        
        (*output each "cvParam" accession number*)
        allCvParamsAccessions =
         #[[2, 2, 
            2]] & /@
          (*output all "cvParam" XMLElement in a spectrum*)
        \
         (Cases[#, XMLElement["cvParam", ___, ___], Infinity] &@
            entireSelectedScan);
        
        (*{{"MS:1000127","centroid spectrum"},{"MS:1000128",
        "profile spectrum"},{"MS:1000129","negative scan"},{"MS:1000130",
        "positive scan"}}*)
        AccessionsOfInterest = {"MS:1000127", "MS:1000128", "MS:1000129", 
          "MS:1000130"};
        
        (*Matched accessions numbers between the selected spectrum and \
      constant AccessionsOfInterest*)
        matchedAccessions = 
         Flatten[DeleteCases[
           Cases[allCvParamsAccessions, #] & /@ AccessionsOfInterest, {}]];
        
        (*mapping matchedAccessions to the dictionary and outputting the \
      corresponding values*)
        mappedAccessionsToDictionary = 
         Values[mzMLdictionaryFunction[[#]] &@matchedAccessions];
        
        (*transposing the matched accessions numbers to the their labels*)
        addingTagsToMatchedAccessions = 
        Transpose[{{"spectrum representation", "ionization mode "}, 
          mappedAccessionsToDictionary}];
       
       (*extracting elements from the selected spectrum*)
        tagingStuffFromselectedScan = Union[findMS[[{2, 3, 4, -1}, {2, 3}]]];
        finalOutPut = 
         Join[tagingStuffFromselectedScan, addingTagsToMatchedAccessions];
        Return[finalOutPut]
    ]

(* ::Function:: *)
(* f:MzmlSpectrumPlotter *)
(***Function***)
MzmlSpectrumPlotter[dictionaryFunction_, binaryArray_, 
  userInputForPeakLabeling_] :=
    Module[ {binaryDataArrayElement,
      mzMLdictionaryFunction = dictionaryFunction,
      binaryDATAArray = binaryArray,
      allCvParamsAccession,
      valuesOfAllAccessions,
      allCvParamsAndBinaryData,
      getTheRealBits,
      decodAndTransposeBinaryData,
      trueOrFalseInput = userInputForPeakLabeling,
      final},
     (*
     Constant binaryDataArrayElement serve as a storage for element \
   <binaryDataArrayList> in a selected scan. Element <
     binaryDataArrayList> 
     is where a list of binary data arrays are located for a given \
   spectrum scan. Data point arrays for default data arrays (m/z, 
     intensity) and meta data arrays for a given spectrum.
     
     Cases here takes 3 element: 
     [1st element in an expiration and that is spectrum1 - 
     selected scan, 
     2nd element is XMLElement pattern matching and that is[
     binaryDataArray], 3rd element is levelspece and that is Infinity - 
     all levels (except 0) - 
     to give a list of all parts of the expiration "1st" on this \
   specified level that match the pattern "2nd"]
     *)
        binaryDataArrayElement = 
         Cases[binaryDATAArray, XMLElement["binaryDataArray", _, _], 
          Infinity];
        (* 
        Constant allCvParamsAccession serve as a storage for all accessions \
      in a given spectrum;
        
        Cases[#,XMLElement["cvParam",___,___],Infinity]&/@
        binaryDataArrayElement;
        Extracting all of XMLElement["cvParam"] that are in XMLElement[
        "binaryDataArray"] for a selected scan for a given spectrum. 
        Element <cvParam> holds additional data or annotation. 
        Only controlled values - in the dictionary - are allowed here;
        
        #[[2,2,2]];
        Mapping into XMLElement["cvParam"] to output all accessions e.g. 
        accession \[Rule] MS:1000576.
        *)
        allCvParamsAccession = #[[2, 2, 2]] & /@ 
            Cases[#, XMLElement["cvParam", ___, ___], Infinity] & /@ 
          binaryDataArrayElement;
        (*
        Constant valuesOfAllAccessions serves as a storage for all \
      accessions VALUES that found in dictionary;
        We mapped all found accessions to our dictionary to find the \
      corresponding VALUES. 
        *)
        valuesOfAllAccessions = 
         Values[mzMLdictionaryFunction[[#]] & /@ allCvParamsAccession];
        (*
        Cases[#,XMLElement["binary",_,_],Infinity]&/@
        binaryDataArrayElement;
        
        To extract the selected spectrum XMLElement["binary"]. 
        There will always be 1 binary element for a given spectrum. 
        Two lists in a given binary element. m/
        z and intensity for 1st and 2nd lists respectively. 
        The binary data is base64 encoded binary data and the byte order is \
      always 'little endian'.;
        
        Transposing each binary list (m/z, 
        intensity) with its corresponding accessions VALUES that found in \
      dictionary that is set to valuesOfAllAccessions.
        *)
        allCvParamsAndBinaryData = 
         Transpose[{valuesOfAllAccessions, 
           Cases[#, XMLElement["binary", _, _], Infinity][[All, 3, 1]] & /@ 
            binaryDataArrayElement}];
        (*
        StringCases[#,x:NumberString~~"-bit "];
        Using StringCases to find the NumberString - 
        represents the characters of a number in StringExpression - 
        that is corresponds to the precision floating-point format. 
        This number is followed by "bit".;
        
        (Switch[ToLowerCase[y],"float","Real","integer","Integer"]<>x);
        Using Switch to find the floating point that follow "bit", 
        if is "float" substitute by "Real" and if "integer"substitute with \
      "Integer" then join - <> - number.
        *)
        getTheRealBits = 
         Flatten[StringCases[#, 
             x : NumberString ~~ "-bit " ~~ 
               y : __ :> (Switch[ToLowerCase[y], "float", "Real", "integer",
                  "Integer"] <> x)]] & /@ valuesOfAllAccessions;
        
        (*We decode each binary data separately, then Transpose. Also, 
        the convention used to interpret the bytes making up the data is \
      little endian correspond to (-1).*)
        decodAndTransposeBinaryData = Transpose[
          {
           MzmlSpectrumBinaryDecoder[
            (*This is the the a selected spectrum m/z binary data.*)
            
            allCvParamsAndBinaryData[[1, 2, 1]],
            (*This is m/
            z binary data type either bit flot (#Real) or bit integer \
        (#Integer).*)
            getTheRealBits[[1, 1]],
            (*m/
            z binary data enocoding formate (Base64) in String format.*)
         \
         {"Base64", "String"},
            -1,
            (*This is m/z binary data copression status, 
            binary data is compressed or not*)
            
            valuesOfAllAccessions[[1, 2]]
            ],
           
           MzmlSpectrumBinaryDecoder[
            (*This is the the a selected spectrum intensity binary data.*)
        
                 allCvParamsAndBinaryData[[2, 2, 1]],
            (*This is intensity binary data type either bit flot (#Real) or \
        bit integer (#Integer).*)
            getTheRealBits[[2, 1]],
            (*intensity binary data enocoding formate (Base64) in String \
        format.*)
            {"Base64", "String"},
            -1,
            valuesOfAllAccessions[[2, 2]]
            ]
           }
          ];
        final =
         (*If user select True, then the plot will have a labeled peaks, 
         peaks are not labeled otherwise*)
         If[ trueOrFalseInput == False,
          (*If False, peaks are not labeled.*)
             Panel[Show[
               (*ListLinePlot will line plot the decodAndTransposeBinaryData \
             values.*)
               ListLinePlot[
                decodAndTransposeBinaryData,
                {(*Using Full for the Plot Range to include the full range of \
             the original data*)
                 PlotRange -> Full,
                 (*valuesOfAllAccessions[[1,3]] = "m/z array" is the x-
                 axis and valuesOfAllAccessions[[2,3]] =
                 "intensity array" is the y-axis.*)
                 
                 AxesLabel -> {valuesOfAllAccessions[[1, 3]], 
                   valuesOfAllAccessions[[2, 3]]},
                 LabelStyle -> Directive[Lighter[Black, 0.1]],
                 AxesStyle -> Lighter[Black, 0.1]
                 },
                ImageSize -> 700],
               
               ListPlot[
                (*Tooltip displays labels while the mouse pointer is in the \
             area where the expression (m/z and intensity) is displayed.*)
                
                Tooltip[decodAndTransposeBinaryData],
                {
                 PlotRange -> Full,
                 PlotMarkers -> {"\[FilledCircle]", 0.01},
                 AxesLabel -> {valuesOfAllAccessions[[1, 3]], 
                   valuesOfAllAccessions[[2, 3]]},
                 LabelStyle -> Directive[Lighter[Black, 0.1]],
                 AxesStyle -> Lighter[Black, 0.1]
                 },
                ImageSize -> 700]
               ]],
             Panel[Show[
             (*ListLinePlot will line plot the decodAndTransposeBinaryData \
             values.*)
             ListLinePlot[
             decodAndTransposeBinaryData,
             {(*Using Full for the Plot Range to include the full range of \
             the original data*)
             PlotRange -> Full,
             (*valuesOfAllAccessions[[1,3]] = "m/z array" is the x-
             axis and valuesOfAllAccessions[[2,3]] =
             "intensity array" is the y-axis.*)
             
             AxesLabel -> {valuesOfAllAccessions[[1, 3]], 
             valuesOfAllAccessions[[2, 3]]},
             LabelStyle -> Directive[Lighter[Black, 0.1]],
             AxesStyle -> Lighter[Black, 0.1]
             },
             ImageSize -> 700]
             ]]
         ];
        Return[final]
    ];


(* ::Subsubsection:: *)
(* MZML Interpretation Function *)

(* ::Function:: *)
(* f:MZMLInterpretationFunction *)
(***Function***)
MZMLInterpretationFunction[filePath_, labelingTimeOrMassRange_] :=
    DynamicModule[{
      mzMLOverallRaw, outputMetaData, mzMLoffsetList,
      (*The first input is set equal to the inputfile and that is the \
    user uploaded file path*)
      inputfile = filePath,
      mzMLdictionaryFunction, mzMLdictionaryMSUO,
      (*The second input is set equal to the controlled vocabulary \
    constant and that is "mzmlControlledVocabulary"*)
      
      timeOrMassLabeling = labelingTimeOrMassRange,
      labeling, interpretationOutcome, metAData, mzMLmsLevel, 
      mzMLpreCursorSwitch, mzMLpreCursor, mzMLtimeMin, mzMLtimeMax, 
      mzMLmassMin, mzMLmassMax, mzMLprecursorSpectrum, mzMLchromatogram, 
      optionForPlotLabeling},
     
     (********** Constants "mzML Back-End" *********)
     
     (*mzMLOverallRaw is equal to function <MzmlRawMetaData> 
     that takes the first input which is the user uploaded file path. 
     The output is the raw metadata.*)
     
     mzMLOverallRaw = MzmlRawMetaData[inputfile];
     (*outputMetaData is equal to function MzmlMetaData that takes the \
   raw metadata and it output a readable metadata for the user*)
     outputMetaData = MzmlMetaData[mzMLOverallRaw];
     
     (*mzMLdictionaryFunction is equal to function \
   MzmlDictionaryAccessionNames that takes the second input the \
   controlled vocabulary constant. 
     Function MzmlDictionaryAccessionNames will make an association \
   between keys "accessions (id)" and values "names" giving the \
   controlled vocabulary constant.*)
     mzMLdictionaryFunction = 
      MzmlDictionaryAccessionNames[timeOrMassLabeling];
     (*mzMLdictionaryMSUO is set equal to function <
     MzmlDictionaryMSUOLabeling> 
     that takes the second input the controlled vocabulary constant. 
     Function MzmlDictionaryMSUOLabeling will make an association \
   between keys "accessions '"MS:" or "UO:
     "' and values "its definition".*)
     mzMLdictionaryMSUO = 
      MzmlDictionaryMSUOLabeling[timeOrMassLabeling];
     
     
     (*labeling is used specifically for retention times lebeling in one \
   of the controls added to allow for user interactive interaction*)
     labeling = MzmlTimeAndMassLabeling[inputfile, mzMLdictionaryMSUO];
     (*mzMLoffsetList is a constant that is set equel to function <
     MzmlOffSetListCounter>that takes mzML file path. 
     MzmlOffSetListCounter output a list of offsets for {spectrum index, 
     ms level, base peak m/z, scan start time, precursor, selected ion m/
     z}. This list is used as an input for downstream functions.*)
     mzMLoffsetList = MzmlOffSetListCounter[inputfile];
     
     
     (*
     1: {Define variables by assignment};
     2: how it appear;
     3: execute
     *)
     
     
     (***Returning mzML Front-End Using Interpretation***)
     
     
     (*Using interpretation to display a version that allow for user \
   interaction, but is interpreted as the unevaluated form*)
     interpretationOutcome = Interpretation[
       
       
       
       (*These are the local defined variables by assignment*)
       {
        (*file metadata will always show by default*)
        
        metAData = True,
        (*"MS Level" variable*)
        mzMLmsLevel = 0,
        (*mzMLpreCursorSwitch is a localy defind variable to allow for \
     user "Search By: " method. 
        Search by "Selected Ion m/z Value" when mzMLpreCursorSwitch = 
        True and by "Base Peak m/z Value" when mzMLpreCursorSwitch = 
        False*)
        mzMLpreCursorSwitch = False,
        mzMLtimeMin = 0,
        mzMLtimeMax = Infinity,
        mzMLmassMin = 0,
        mzMLmassMax = Infinity},
       
       
       
       
       (*This is how it appear*)
       Panel[
        Grid[{{Panel[Style["Input File Path: " <> inputfile, 11]]},
          
          
          {Panel[Grid[{
              
              (*This is one of the defined variables that allow for an \
     interactive way for the user to search for different spectra "MS \
     Level" by using PopupMenu and is set by default to show all MS Levels \
     in the uploaded file*)
              {"MS Level: " PopupMenu[
                 Dynamic[mzMLmsLevel], 
                 Join[{0 -> Style["All", 10]}, 
                  ToExpression[#] & /@ 
                   Union[Cases[mzMLoffsetList, 
                     x_ /; (Length[x] > 1) && (x[[2]] == 
                       "ms level value") :> x[[3]], 2]]]]},
              
              (*This is one of the the defined variables that allow for \
     an interactive way for the user to search the spectra by "Selected \
     Ion m/z Value" or "Base Peak m/z Value" by using RadioButton and is \
     set by default to search by "Base Peak m/z Value"*)
              {Style[
                "Search By: ", 10], 
               RadioButtonBar[
                Dynamic[
                 mzMLpreCursorSwitch], {True -> "Selected Ion m/z Value",
                  False -> "Base Peak m/z Value"}, 
                Appearance -> "Horizontal" , 
                BaselinePosition -> Automatic, Appearance -> Small]},
              
              (*This is one of the defined variables that allow for an \
     interactive way for the user to change/
              input "minimum retention time" viewing MS spectrum by \
     using InputField and is set to search starting at time=
              0*)
              {Style[
                " Input \"minimum [\[GreaterEqual] 0]\" retention time \
\"" <> (# & /@ labeling) <> "\": ", 10, Darker[Green, 0.4]], 
               InputField[Dynamic[mzMLtimeMin], FieldSize -> 5]},
              
              (*This is one of the defined variables that allow for an \
     interactive way for the user to change/
              input "maximum retention time" viewing MS spectra by using \
     InputField and is set to time=
              Infinity*)
              {Style[
                " Input \"maximum\" retention time \"" <> (# & /@ 
                   labeling) <> "\": ", 10, Darker[Green, 0.4]], 
               InputField[Dynamic[mzMLtimeMax], FieldSize -> 5]},
              
              (*This is one of the defined variables that allow for an \
     interactive way for the user to change/
              input "minimum spectra mass" by using InputField and is \
     set to mass=
              0*)
              {Style[
                " Input \"minimum [\[GreaterEqual] 0]\" spectra mass \
\"m/z\": ", 10, Lighter[Blue, 0.2]], 
               InputField[Dynamic[mzMLmassMin], FieldSize -> 5]},
              
              (*This is one of the defined variables that allow for an \
     interactive way for the user to change/
              input "maximum spectra mass" by using InputField and is \
     set to mass=
              Infinity*)
              {Style[
                " Input \"maximum\" spectra mass \"m/z\": ", 10, 
                Lighter[Blue, 0.2]], 
               InputField[Dynamic[mzMLmassMax], FieldSize -> 5]}
              
              }, Alignment -> Left]]},
          
          (*This is the metadata in a panel form*)
          {Panel[
            Grid[{{Style["View File Metadata: ", 10], 
               RadioButtonBar[Dynamic[metAData],
                {True -> "Yes", False -> "No"} , 
                BaselinePosition -> Automatic, 
                Appearance -> Small]}}]]},
          {Panel[
            Grid[{{Button[Style["Evaluate", 12], 
                SelectionMove[EvaluationCell[], All, CellContents];
                SelectionEvaluateCreateCell[EvaluationNotebook[]](*,
                Method\[Rule]"Queued"*)]}}]]}},
         
         Alignment -> Left]],
       
       
       
       
       
       
       
       (*This is how it executed*)
       
       
       
       (***Returning mzML Front-End Using Interpretation***)
       
       Interpretation[
        
        
        (*These are another set of locally defind variables*)
        {
         (*mzMLpreCursor is dummy variable in place to allow for passing \
     of the selected spectrum*)
         mzMLpreCursor,
         mzMLchromatogram,
         mzMLprecursorSpectrum,
         optionForPlotLabeling},
        
        
        
        (*This is how it appear*)
        Panel[Grid[{{
            (*Output file metadate once metAData\[Equal]True*)
            
            If[ metAData == True,
                Panel[outputMetaData]
            ]},
           
           
           {Panel[
             Grid[{
               {Grid[
                 Transpose[
                  
                  (*This is one of the variables that allow for an \
     interactive way for the user to search the spectra by "Selected Ion \
     m/z Value" or "Base Peak m/z Value" by using RadioButton*)
              \
       {{Style[
                     "Please select: {" ~~ 
                      If[ mzMLpreCursorSwitch == True,
                          "Selected Ion m/z Value",
                          "Base Peak m/z Value"
                      ] ~~
                       ", retention time, MS Level } : ", 10]},
                   
                   (*This is one of the varables that allow for an \
     interactive way for the user to select a spectrum using either \
     "Selected Ion m/z Value" or "Base Peak m/z Value", max/
                   min retention time and MS Level of intrest using \
     PopupMenu*)
                   {PopupMenu[
                     Dynamic[mzMLpreCursor],
                     Sort[
                      MzmlSpectrumMassTimeAssociation[
                       Which[
                       mzMLpreCursorSwitch == True,
                       (MzmlSpectrumSelector[mzMLoffsetList, mzMLmsLevel,
                        "Yes precursor-ion",
                       Replace[mzMLtimeMin, Null -> 0],
                       
                       Replace[
                       mzMLtimeMax, {Infinity -> 1000000, 
                       Null -> 1000000}],
                       Replace[mzMLmassMin, Null -> 0],
                       
                       Replace[ 
                       mzMLmassMax, {Infinity -> 1000000, 
                       Null -> 1000000}]]),
                       mzMLpreCursorSwitch == False,
                       (MzmlSpectrumSelector[mzMLoffsetList, mzMLmsLevel,
                        "No precursor-ion",
                       Replace[mzMLtimeMin, Null -> 0],
                       
                       Replace[
                       mzMLtimeMax, {Infinity -> 1000000, 
                       Null -> 1000000}],
                       Replace[mzMLmassMin, Null -> 0],
                       
                       Replace[
                       mzMLmassMax, {Infinity -> 1000000, 
                       Null -> 1000000}]])], mzMLpreCursorSwitch]]]}
                   }], Alignment -> Left]},
               
               {Panel[
                 Grid[
                  {
                   (*This is one of the variables that allow for an \
     interactive way for the user to locate the "precursor spectrum" for \
     any selected "daughter spectrum" by using RadioButtonBar*)
              \
        {Style["Find precursor spectrum for MS Level > 1: ", 10], 
                    RadioButtonBar[
                     
                     Dynamic[mzMLprecursorSpectrum], {True -> "No", 
                      False -> "Yes"} , Appearance -> Tiny]},
                   
                   (*This is one of the variables that allow for an \
     interactive way for the user to "Peaks label" all peaks in selected \
     spectrum by using RadioButton*)
                   {Style[
                     "Peaks Labeling: ", 10], 
                    RadioButtonBar[
                     Dynamic[optionForPlotLabeling], {True -> "No", 
                      False -> "Yes"} , Appearance -> Tiny]}}, 
                  Alignment -> Left]]}
               }, Alignment -> Left]]},
           
           (*This is one of the controls that allow for an interactive \
     way for the user to choose viewing the "Total Ion Current (TIC) \
     Chromatogram" by using RadioButton*)
           {Panel[
             Grid[{{Style["View Chromatogram: ", 10], 
                RadioButtonBar[Dynamic[mzMLchromatogram],
                 {True -> "No", False -> "Yes"} , 
                 BaselinePosition -> Automatic, Appearance -> Small]}}]]},
           
           {Panel[
             Grid[{{Button[Style["Evaluate", 12], 
                 SelectionMove[EvaluationCell[], All, CellContents];
                 SelectionEvaluateCreateCell[EvaluationNotebook[]](*,
                 Method\[Rule]"Queued"*)]}}]]}
           
           
           }, Alignment -> Left]],
        
        
        
        
        
        (*This is how it executes*)
        Grid[{
          
          
          {(*if the user choose to see the total ion chromatogram by \
     clicking on a radio button, then pass True. Once True is passed, 
           plot the "Total Ion Current (TIC) Chromatogram" using \
     function <MzmlSpectrumPlotter>*)
           
           If[ mzMLchromatogram == False,
               Panel[Insert[
                 Grid[{{Style["Chromatogram \"TIC\":", 15, Underlined, 
                     Bold], ""}, {If[ mzMLchromatogram == False,
                                      Framed[MzmlSpectrumPlotter[mzMLdictionaryFunction, 
                                        MzmlTotalIonChromatogram[inputfile, 
                                         MzmlGetChromatogramTIC[mzMLoffsetList]], 
                                        optionForPlotLabeling]]
                                  ]}}], Alignment -> Left, 2]]
           ]},
          {
           If[ Length[mzMLpreCursor] < 1,
               Style["Please Make a Selection....", 10],
               If[
                
                (*using function <MzmlGetSpectraMSLevel> 
                to find the selected spectrum MS level. If the MS Level = 
                1 then plot the spectrum using function <
                MzmlSpectrumPlotter> 
                and creat a spectrum note card label using function <
                MzmlSpectrumPlotLabels>*) MzmlGetSpectraMSLevel[mzMLpreCursor] == {"1"},
                   Panel[Grid[
                     {
                      {Insert[
                        Grid[Transpose[
                          MzmlSpectrumPlotLabels[inputfile, mzMLpreCursor, 
                           mzMLdictionaryFunction]]], {Background -> {None, \
                   {White, {Lighter[
                             Blend[{Lighter[LightBlue, 0.1], 
                             Lighter[LightGray, 0.9]}], .1]}}}, 
                         Dividers -> {Black, {2 -> Black}}, Frame -> True, 
                         Spacings -> {2, {2, {0.7}, 0.9}}}, 2]},
                      
                      {MzmlSpectrumPlotter[mzMLdictionaryFunction, 
                        MzmlGetMSSpectraLocationAndBinary[inputfile, 
                         Dynamic[
                          Cases[mzMLpreCursor, 
                            x_ /; x[[2]] == "index" :> x[[1]], 1][[1]]]], 
                        optionForPlotLabeling]}
                      },
                     Dividers -> {{{True, True}}, {{True, False}}}]],
                   
                   
                   
                   (*If the MS Level \[NotEqual] 1, using a Switch to find/
                   locate the "precursor spectrum" for any selected "daughter \
         spectrum"*)
                   Switch[mzMLprecursorSpectrum,
                    
                    
                    
                    (*If user does want to see the "precursor spectrum" and \
                   True is selected then pass "Yes precursor-ion" to get the "precursor \
                   spectrum" for any selected "daughter spectrum". To do so, 
                    we call MzmlGetMSSpectraLocationAndBinary function for \
                   precursor scans locations.*)
                    True,
                    Panel[Grid[
                      {
                       
                       {Insert[
                         Grid[Transpose[
                           MzmlSpectrumPlotLabels[inputfile, mzMLpreCursor, 
                            mzMLdictionaryFunction]]], {Background -> {None, \
                   {White, {Lighter[
                             Blend[{Lighter[LightBlue, 0.1], 
                             Lighter[LightGray, 0.9]}], .1]}}}, 
                          Dividers -> {Black, {2 -> Black}}, Frame -> True, 
                          Spacings -> {2, {2, {0.7}, 0.9}}}, 2]},
                       
                       {MzmlSpectrumPlotter[mzMLdictionaryFunction, 
                         MzmlGetMSSpectraLocationAndBinary[inputfile, 
                          Dynamic[
                           Cases[mzMLpreCursor, 
                             x_ /; x[[2]] == "index" :> x[[1]], 1][[1]]]], 
                         optionForPlotLabeling]}
                       },
                      Dividers -> {{{True, True}}, {{True, False}}}]],
                    
                    
                    
                    (*If MSLevel \[NotEqual] 1, 
                    use a switch and if False is selected then pass "No \
                   precursor-ion" but if True is selected then pass "Yes precursor-ion" \
                   to get the precursor scan AKA. parant scan. To do so, 
                    we call MzmlGetPrecursorSpectraAtMSLevel function for \
                   precursor scans locations. 
                    We then map into it the binary data and call \
                   mzXMLstreamingAndBinaryAndspectraPlotting to decode and plot the scan \
                   binary data.*)
                    False,
                    
                    Panel[Grid[
                      Partition[
                       Flatten[
                        Transpose[{Insert[
                             Grid[Transpose[#]], {Background -> {None, {White, \
                   {Lighter[Blend[{Lighter[LightBlue, 0.1], 
                             Lighter[LightGray, 0.9]}], .1]}}}, 
                             Dividers -> {Black, {2 -> Black}}, Frame -> True, 
                             Spacings -> {2, {2, {0.7}, 0.9}}}, 
                             2] & /@ (MzmlSpectrumPlotLabels[inputfile, #, 
                             mzMLdictionaryFunction] & /@ \
                   (MzmlGetPrecursorSpectraAtMSLevel[mzMLoffsetList, 
                             mzMLpreCursor])), (MzmlSpectrumPlotter[
                             mzMLdictionaryFunction, #, 
                             optionForPlotLabeling] & /@ \
                   (MzmlGetMSSpectraLocationAndBinary[
                             inputfile, #[[1]]] & /@ \
                   (MzmlGetPrecursorSpectraAtMSLevel[mzMLoffsetList, 
                             mzMLpreCursor])))}]], 1],
                      
                      Dividers -> {{{True, True}}, {{True, False}}}]]]
               ]
           ]}}, 
         Alignment -> Center]]
       ]];
    

(* ::Subsection:: *)
(* MZXML *)

(* ::Subsubsection:: *)
(* MZXML MetaData *)

(* ::Function:: *)
(* f:MzxmlRawMetaData *)
(***Function***)
MzxmlRawMetaData[mzXMLFile_] :=
    Module[ {inputFile = mzXMLFile,
      emptyList, filePath, allMetaData, testXML},
        emptyList = {};
        filePath = inputFile;
        (*open and read mzXML file*)
        inputFile = OpenRead[filePath];
        
        (*Set the the streaming position to 0. This will stream metadata \
      from the begining of the file*)
        SetStreamPosition[filePath, 0];
        (*while loop*)
        While[True,
         allMetaData = Read[filePath, String];
         
         (*Using <StringCases> to search "scan". <..> 
         is to keep repeating the pattern matching by representing "scan" \
      until it matches, if it doesn't match, keep appending <AppendTo> 
         to the emptyList={}. Using <If> and <Break[]> 
         to stop while loop once "scan" matchs.*)
         If[ StringCases[allMetaData, WordCharacter ..][[1]] == "scan",
             Break[],
             AppendTo[emptyList, allMetaData]
         ]];
        
        (*</msRun></mzXML> as closing tags for the "XML" format.*)
        emptyList = AppendTo[emptyList, "</msRun></mzXML>"];
        
        (*Import as "XML"*)
        testXML = ImportString[StringJoin[emptyList], "XML"];
        Close[inputFile];
        Return[testXML]
    ]
  
(* ::Function:: *)
(* f:MzxmlVersion *)
(***Function***)
MzxmlVersion[mzXMLSchemaVersion_] :=
    Module[ {mzXMLSchema = mzXMLSchemaVersion, mzXMLSchemaTag, 
      allTogether},
     
     (*using string pattern matching to look for "mzXML_" and <
     NumberString> to find the mzXML version format. 
     "mzXML_" will always be followed by a number that corresponds to \
   the mzXML schema version. This will be outputed in the <MetaData>.*)
        mzXMLSchemaTag = Grid[StringCases[ToString[mzXMLSchema],
          ___ ~~ "mzXML_" ~~ x : NumberString ~~ ___ :> {"mzXML:", x}]];
       
       (*allTogether is to output a readable mzXML schema version *)
        allTogether = 
         Insert[Insert[
           Grid[Transpose[
             Transpose[{{Style["mzXML Schema:", Bold, 
                 Underlined]}, {mzXMLSchemaTag}}]]], Alignment -> Left, 2], 
          Spacings -> {2, {2, {0.7}, 0.9}}, 2];
        Return[allTogether]
    ];

(* ::Function:: *)
(* f:MzxmlMSRun *)
(***Function***)
MzxmlMSRun[msRunXMLElement_] :=
    Module[ {msRunElement,
      msRUN = msRunXMLElement,
      msRunTag1, msRunTag2, msRunTag3, msRunTag4, allTogether},
     
     (*using cases to look for the <msRun> tag in the data file*)
        msRunElement = Cases[msRUN, XMLElement["msRun", _, _], Infinity];
        
        (*using If statments to look for the attributes "scanCount, \
      startTime, endTime" since they're optional. 
        If any of thema are not found, output {"Unknown"}*)
        msRunTag1 = {If[ Cases[#, 
               x_ /; Length[x] > 1 && x[[1]] == "scanCount" :> x[[2]], 
               All] == {},
                         {"Unknown"},
                         Cases[#, 
                          x_ /; Length[x] > 1 && x[[1]] == "scanCount" :> x[[2]], All]
                     ],
             If[ Cases[#, 
               x_ /; Length[x] > 1 && x[[1]] == "startTime" :> x[[2]], 
               All] == {},
                 {"Unknown"},
                 Cases[#, 
                  x_ /; Length[x] > 1 && x[[1]] == "startTime" :> x[[2]], All]
             ],
             If[ Cases[#, 
               x_ /; Length[x] > 1 && x[[1]] == "endTime" :> x[[2]], 
               All] == {},
                 {"Unknown"},
                 Cases[#, x_ /; Length[x] > 1 && x[[1]] == "endTime" :> x[[2]], 
                  All]
             ]} & /@ msRunElement;
        
        (*StringReplace to get rid of letters found in the attributes \
      "above", only numbers are outputed and stored in msRunTag2*)
        msRunTag2 = 
         StringReplace[#, LetterCharacter -> ""] & /@ msRunTag1[[1]];
        
        (*insert appropiate names for each attribute in msRunTag2*)
        msRunTag3 = {StringInsert[#[[1]], "Scan Count: ", 1], 
            StringInsert[#[[2]], "Start Time: ", 1], 
            StringInsert[#[[3]], "End Time: ", 1]} &@msRunTag2;
        msRunTag4 = Insert[Grid[msRunTag3], Alignment -> Left, 2];
        
        (*allTogether is to output a readable mzXML schema version *)
        allTogether = 
         Insert[Insert[
           Grid[Transpose[
             Transpose[{{Style["MS Run:", Bold, 
                 Underlined]}, {msRunTag4}}]]], Alignment -> Left, 2], 
          Spacings -> {2, {2, {0.7}, 0.9}}, 2];
        Return[allTogether]
    ];

(* ::Function:: *)
(* f:MzxmlParentFile *)
(***Function***)
MzxmlParentFile[parentFileXMLElement_] :=
    Module[ {parentFileElement, parentFILE = parentFileXMLElement, 
      parentFileTag1, parentFileTag2, allTogether},
     
     (*Use Cases to locate XMLElement["parentFile"].*)
        parentFileElement = 
         Cases[parentFILE, XMLElement["parentFile", _, _], Infinity];
        
        (*"parentFile" contains 3 required attributes, If [
        "parentFile"]\[Equal]{},output {"Unknown"}.*)
        parentFileTag1 =
         {Cases[#, x_ /; Length[x] > 1 && x[[1]] == "fileName" :> x[[2]], 
             All],
            Cases[#, x_ /; Length[x] > 1 && x[[1]] == "fileType" :> x[[2]],
              All],
            Cases[#, x_ /; Length[x] > 1 && x[[1]] == "fileSha1" :> x[[2]],
              All]} & /@ parentFileElement;
        
        (*use <StringReplace> to edit the attributes: e.g. 
        "file"\[Rule]"File " and use <Transpose> 
        to match Keys with the corresponding Values*)
        parentFileTag2 = 
         TableForm[
          Transpose[{StringInsert[#[[1]], "File Name: ", 1], 
              StringInsert[#[[2]], "File Type: ", 1], 
              StringInsert[#[[3]], "File Sha1: ", 1]} & /@ 
            parentFileTag1, {1, 3, 2}]];
        
        (*allTogether is to output a readable mzXML schema version *)
        allTogether = 
         Insert[Insert[
           Grid[Transpose[
             Transpose[{{Style["Parent File: ", Bold, 
                 Underlined]}, {parentFileTag2}}]]], Alignment -> Left, 2],
           Spacings -> {2, {2, {0.7}, 0.9}}, 2];
        Return[allTogether]
    ];

(* ::Function:: *)
(* f:MzxmlMSInstrument *)
(***Function***)
MzxmlMSInstrument[msInstrumentXMLElement_] :=
    Module[ {msInstrumentElement,
      msINSTRUMENT = msInstrumentXMLElement,
      msInstrumentTag1, msInstrumentTag2, allTogether},
     
     (*ues Cases to locate XMLElement["msInstrument"] in the MetaData.*)
        msInstrumentElement = 
        Cases[msINSTRUMENT, XMLElement["msInstrument", _, _], Infinity];
        msInstrumentTag1 = Transpose[Flatten[{
              If[ Cases[#, 
                 x_ /; Length[x] > 1 && x[[1]] == "category" :> x[[2]], 
                 All] == {},
                  {"Unknown"},
                  Cases[#, 
                   x_ /; Length[x] > 1 && x[[1]] == "category" :> x[[2]], 
                   All]
              ], If[ Cases[#, x_ /; Length[x] > 1 && x[[1]] == "value" :> x[[2]], 
                 All] == {},
                     {"Unknown"},
                     Cases[#, x_ /; Length[x] > 1 && x[[1]] == "value" :> x[[2]], 
                      All]
                 ]} & /@ msInstrumentElement, 1]];
        
        (* Uses Cases to locate XMLElement["msManufacturer"], ["msModel"], [
        "msIonisation"], ["msMassAnalyzer"] and ["msDetector"] in the <
        msInstrument> Element*)
        msInstrumentTag2 = 
         Insert[Grid[
           StringReplace[#, {"msManufacturer" :> "MS-Manufacturer: ", 
               "msModel" :> "MS-Model: ", 
               "msIonisation" :> "MS-Ionisation: ", 
               "msMassAnalyzer" :> "MS-Analyzer: ", 
               "msDetector" :> "MS-Detector: ", 
               "msResolution" :> "MS-Resolution: "}, IgnoreCase -> True] & /@
             msInstrumentTag1], Alignment -> Left, 2];
        
        (*allTogether is to output a readable mzXML schema version *)
        allTogether = 
         Insert[Insert[
           Grid[Transpose[
             Transpose[{{Style["Instrument:", Bold, 
                 Underlined]}, {msInstrumentTag2}}]]], Alignment -> Left, 
           2], Spacings -> {2, {2, {0.7}, 0.9}}, 2];
        
        (*Using <Return> statment to return all six sub-
        elements in a list that is titled "Instrument"*)
        Return[allTogether]
    ];

(* ::Function:: *)
(* f:MzxmlInstrumentSoftware *)
(***Function***)
MzxmlInstrumentSoftware[instrumentSoftwareXMLElement_] :=
    Module[ {msInstrumentElement,
      instrumentSOFTWARE = instrumentSoftwareXMLElement,
      softwareXMLElement, softwareTag1, 
      softwareTag2, softwareTag3},
     
     (* use Cases to locate XMLElement[
     "msInstrument"] first in the MetaData*)
        msInstrumentElement = 
         Cases[instrumentSOFTWARE, XMLElement["msInstrument", _, _], 
          Infinity];
        
        (*we use Cases to locate XMLElement["software"] in the XMLElement[
        "msInstrument"]*)
        softwareXMLElement = 
         Cases[msInstrumentElement, XMLElement["software", _, _], 
          Infinity];
        
        
        (*the software element contain 3 required attributes, type, 
        name and version of the software.*)
        softwareTag1 = {If[ Cases[#, x_ /; Length[x] > 1 && x[[1]] == "type" :> x[[2]], 
               All] == {},
                            {"Unknown"},
                            Cases[#, x_ /; Length[x] > 1 && x[[1]] == "type" :> x[[2]], 
                             All]
                        ], If[ Cases[#, x_ /; Length[x] > 1 && x[[1]] == "name" :> x[[2]], 
               All] == {},
                               {"Unknown"},
                               Cases[#, x_ /; Length[x] > 1 && x[[1]] == "name" :> x[[2]], 
                                All]
                           ], If[ Cases[#, x_ /; Length[x] > 1 && x[[1]] == "version" :> x[[2]], 
               All] == {},
                                  {"Unknown"},
                                  Cases[#, x_ /; Length[x] > 1 && x[[1]] == "version" :> x[[2]], 
                                   All]
                              ]} & /@ softwareXMLElement;
        
        (*edit all 3 titels {"type"\[Rule]"Type:","name"\[Rule]"Name:",
        "version"\[Rule]"Version:"} using StringInsert*)
        softwareTag2 = 
         Insert[Grid[
           Flatten[{StringInsert[#[[1]], "Type: ", 1], 
               StringInsert[#[[2]], "Name: ", 1], 
               StringInsert[#[[3]], "Version: ", 1]} & /@ softwareTag1, 1]],
           Alignment -> Left, 2];
        
        (*allTogether is to output a readable mzXML schema version *)
        softwareTag3 = 
         Insert[Grid[
           Transpose[
            Transpose[{{Style["Instrument Software:", Bold, 
                Underlined]}, {softwareTag2}}]]], Alignment -> Left, 2];
        Return[softwareTag3]
    ];

(* ::Function:: *)
(* f:MzxmlDataProcessingSoftware *)
(***Function***)
MzxmlDataProcessingSoftware[dataProcessingsoftwareXMLElement_] :=
    Module[ {dataPROCESSINGSOFTWARE = dataProcessingsoftwareXMLElement,
      dataProcessingElement, dataProcessingSoftwareTag1, 
      dataProcessingSoftwareTag2, dataProcessingSoftwareTag3, 
      dataProcessingSoftwareTag4, processingOperationTag1, 
      processingOperationTag2, processingOperationTag3, 
      processingOperationTag4, allTogether},
     
     (*use Cases to locate XMLElement[
     "dataProcessing"] first in the MetaData*)
        dataProcessingElement = 
         Cases[dataPROCESSINGSOFTWARE, XMLElement["dataProcessing", _, _], 
          Infinity];
        
        (*we then use Cases again to locate XMLElement[
        "software"] in the XMLElement["dataProcessing"].*)
        dataProcessingSoftwareTag1 = 
         Cases[dataProcessingElement, XMLElement["software", _, _], 
          Infinity];
        
        
        (*the software element contain 3 required attributes, type, 
        name and version of the software.*)
        dataProcessingSoftwareTag2 = {If[ Cases[#, x_ /; Length[x] > 1 && x[[1]] == "type" :> x[[2]], 
               All] == {},
                                          {"Unknown"},
                                          Cases[#, x_ /; Length[x] > 1 && x[[1]] == "type" :> x[[2]], 
                                           All]
                                      ], If[ Cases[#, x_ /; Length[x] > 1 && x[[1]] == "name" :> x[[2]], 
               All] == {},
                                             {"Unknown"},
                                             Cases[#, x_ /; Length[x] > 1 && x[[1]] == "name" :> x[[2]], 
                                              All]
                                         ], If[ Cases[#, x_ /; Length[x] > 1 && x[[1]] == "version" :> x[[2]], 
               All] == {},
                                                {"Unknown"},
                                                Cases[#, x_ /; Length[x] > 1 && x[[1]] == "version" :> x[[2]], 
                                                 All]
                                            ]} & /@ dataProcessingSoftwareTag1;
        
        
        (*edit all 3 titels {"type"\[Rule]"Type:","name"\[Rule]"Name:",
        "version"\[Rule]"Version:"} using StringInsert*)
        dataProcessingSoftwareTag3 = 
         Insert[Grid[
           Flatten[{StringInsert[#[[1]], "Type: ", 1], 
               StringInsert[#[[2]], "Name: ", 1], 
               StringInsert[#[[3]], "Version: ", 1]} & /@ 
             dataProcessingSoftwareTag2, 1]], Alignment -> Left, 2];
        dataProcessingSoftwareTag4 = 
         Insert[Grid[
           Transpose[
            Transpose[{{Style["Data Processing Software:", Bold, 
                Underlined]}, {dataProcessingSoftwareTag3}}]]], 
          Alignment -> Left, 2];
        
        (*The subelement "processingOperation" is any additional \
      manipulation not included elsewhere in the dataProcessing element \
      located, If ["processingOperation"]\[Equal]{},output {"Unknown"}.*)
        processingOperationTag1 = 
         If[ Cases[dataProcessingElement, 
            XMLElement["processingOperation", _, _], 
            Infinity] == {},
             {"Unknown"},
             Cases[dataProcessingElement, 
              XMLElement["processingOperation", _, _], Infinity]
         ];
        
        (*the processingOperation element contain 3 optional attributes, 
        name, value and type of the processing operation.*)
        processingOperationTag2 = {If[ Cases[#, x_ /; Length[x] > 1 && x[[1]] == "name" :> x[[2]], 
               All] == {},
                                       {"Unknown"},
                                       Cases[#, x_ /; Length[x] > 1 && x[[1]] == "name" :> x[[2]], 
                                        All]
                                   ], If[ Cases[#, x_ /; Length[x] > 1 && x[[1]] == "value" :> x[[2]], 
               All] == {},
                                          {"Unknown"},
                                          Cases[#, x_ /; Length[x] > 1 && x[[1]] == "value" :> x[[2]], 
                                           All]
                                      ], If[ Cases[#, x_ /; Length[x] > 1 && x[[1]] == "type" :> x[[2]], 
               All] == {},
                                             {"Unknown"},
                                             Cases[#, x_ /; Length[x] > 1 && x[[1]] == "type" :> x[[2]], 
                                              All]
                                         ]} & /@ processingOperationTag1;
        processingOperationTag3 = 
         Insert[Grid[
           Flatten[{StringInsert[#[[1]], "Name: ", 1], 
               StringInsert[#[[2]], "Value: ", 1], 
               StringInsert[#[[3]], "Type: ", 1]} & /@ 
             processingOperationTag2, 1]], Alignment -> Left, 2];
        processingOperationTag4 = 
         Insert[Grid[
           Transpose[
            Transpose[{{Style["Data Processing Operation", Bold, 
                Underlined]}, {processingOperationTag3}}]]], 
          Alignment -> Left, 2];
        
        (*allTogether is to output a readable mzXML schema version *)
        allTogether = 
         Insert[Grid[
           Transpose[
            Transpose[{{dataProcessingSoftwareTag4}, \
        {processingOperationTag4}}]]], Alignment -> Left, 2];
        Return[allTogether]
    ];

(* ::Function:: *)
(* f:MzxmlMetaData *)
(***Function***)
MzxmlMetaData[input_] :=
    Module[ {
      (*overall is equal to the only input to this function*)
      
      overall = input,
      metaDataOverall},
        metaDataOverall = Grid[Partition[
           (*this is the viewer title joined with the other functions*)
          
            Join[{Insert[
              Grid[Transpose[
                Partition[{Style["MathIOmica MSViewer", 38, Bold], \!\(\*
        GraphicsBox[
        TagBox[RasterBox[CompressedData["
        1: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
        jpmo8wy6M5VMP9O86qJ/D6nlkJYf1sdHfxUz1VfC36OuGSGfCumHsqwcao2s
        PNq1CuxQ3ATYd/Dw0N9qwS+r1rfh/1HHVmX/T1css4vOGMkyyNdH766lB8ZG
        lhQyoy1mKT+d/D3kT/Fn5jVV/a2RlaprtWRnAvtDkfeg5ZfVuF4t2U7+n65Y
        ZhedMVLvc1dd5VgLfljK7yEd1pK/h+IP6pi9xb4dVlyz2s4E9ofEx+T4Wcgx
        mMK0beEIceZo1Zl0iuWcaZ510dlD0l48YTrO7zXscWfXqjvV+9nLS79VyL+j
        6bCyj1POCz/qHrcSl2qZ2VecAfY698oHO+UP+lXsYZ3jrHIsUz5LkPtgtd59
        pvnbRWeePydJtNcMfhTMFVHv57DgsmZmWYfViglSHgF0Qtnfc8p0QhUuL2Km
        HCBDngFmTKvvQ/RmV+uyScxUm6nlpzoluuwZ4GsXXk8/Fs8/BXS07zXISk9P
        fbqW/idwhniClg6rlUMj+3so5+6UxgOle6C8djlvdnBsVHGUZrwkPKFIde1y
        XEcdy9TqM3IiOsiJnI2ng79dOD11WFTXUtC+s5YuUJH6+y9ydatcot7HbrrK
        55mpnkvOUbdU90w5p+o6jVOaP6NaLzlmAt3TkM8Ou7gTLkPLVjBzWa+uDyVd
        Fn7ZY/qQTxJ7al6cav4eE6tnev+fTXQMfmnxh2qDqT7YGqLv6cZM/lD6G6/3
        UPbxk+p9HOBXcRkzr+eyFLMkXfBYMfpTZqNJayj7fzryf5Bjr74fRb+9n5nO
        Fa3rJz+yg4XrP9a+18KeFl9PN2/VeNW83jONibMMj1pYlOvzqUa/h6CeEvWS
        6BINkl/vqSK5P0d36TcMcQbN+EhUcSnTrn3OZp1zZygOCB32tOe0qdazky6L
        uWLq+3HVVX7OzPPs1fKe/LJyjUkf1tlndSy9RsaeGkcdvPXy0mfI5O2tz4yJ
        LB16FP7KPLbEXxmvMk67MGoZj2ZY9PTShyM2AR1SkVd3ObtWPaElt9QE/4yL
        a+VTBnKrvEfnXrkM5OOjr/L3q8gVexPUV+wr6mcl97Qi3vfBZzr9hmvVD8w8
        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
        Z9Prd26eN1E6/7rNm1i8vHfvPic996fnxb58sHPc733DR/P+UbXcLariyODB
        Q3D//5H29BUFhRU7ggtnvx+Y1fJuc0PN4antDasUW4Ly2fKWzRm9a82K9vUF
        +WU1Hp76gyEhY/j4URE8O69cweStOM+4m6I3eKROoXkEBpszPDQxyt6+4E3q
        gU24DPLLsQaX36hwSfYl/LHIK0BOUrrYe+CXk4TLs65H7CnG5oXMVG9CuWLQ
        Z+XYiSEXiImcrOCA8m1auIyJyO9Y89aWsWY8GKUvNDxW17Xx3Cte5KM3vcoT
        5z+i6LNV9B7Yl7BVYePcDVyGpLRzz8r13CF/Ebf3L+QbLi//6IHbl7zS0twA
        W2TJ0kVToKvibIUPdtz8We3X1je0jPMbf+P3wH7hVft5cuv1PwRkjj30yD23
        AieG3ifPPXFn4qGXH/YS937hkw/dEfLso3fhPs+7etOt8WFzH/0jY+ULPGTO
        o39s3LkvS71+23Y95eM5dsvzXvXbXxqSOeXF5YsmLPvoqeK3vn1l5AP0nrVL
        R86Tcbnu0ukTrNi/HHmv/Xr6fesRkvXd4MR6PiRpPE/IHI38b5KVa7OzCzcF
        5k55P2TEnHe94qte3XF15fMv3FV2+OaNLehbYKh5WzCt4dqysnL4qQ114IF+
        FbcXZE3godFjuV1qG9e138t9ZzzIvafexxOGT0b/invCgjKqhg7News9IAy9
        smwLOvrT+3vnWoNLzOvDObxbhcsZEi7TJFxSb4K+7Byo6TpF2NSyNWVsQr/v
        yPeIDC1bcTRcTmkztx9LCzN5dlY2rx8ZzLNb1/LFj3/Jp931Ea+bt9MgU5kx
        v4xspns8vYt5dN3V3KvqGm6fNpW7uhbyebMaH2WiDkHRXZfNmNYKXkJWNqWm
        pjUtmtsG2ypr9qL1SxRM8XE3vMmbd73NRzcvvf7JB3fWMFE7e/dt230+PPxU
        rODp+fueuWvYa/vug7523r4ntg0vbZzzev32N3njzW/zjPkP//f9j7/vS+v2
        7gdfXRI9Yec9EZNuedG9dtPz6SOapr3zxms9X3p8+Yyfvn97kAmXlUsIl2PG
        jO3eVp8/31pc9u3p901AYPDTdnHV3L9kKnfJnqLIOb2h/0pycup1YflT3w4p
        XvCeV3Ld2yWlo5a0tTVNmdOW+uoNmxYtYKK2u7IkZdvMmdPgTzH0m8vJySvU
        l1Q/Zh8//s+hsY3cIW8O9535EE+59GkeMeVOxQatQTzlT/ThMZChj11+x8yI
        QN8ca+Zx/sJM9uU2gUvUPkBe1gtcIhYTwUzyknBJObH/VFzKtqacC0RxzaFi
        PeA392fSjOnI0NKVlnA5eeI4Mx7UVOUb8NdcX8rrxozkRWte5Ne8/hNv2/0u
        Hzd+Kt7zFjP6/B90cHC8P7h0Bdc3TVdweTUPCK/htkOy+KzpLcCsodY+NDR8
        0YiiYvj2DfURlZVVbevXXrq4X79+adMXrl2X1Hbjr3lXvsSvOfg9b9n66ps3
        bpoPeWnIpdyx7Tr711+6PY7u//GH7un76N7/hGJdGkf6vj662IvnzNrxx1UH
        f+RNN7/D29c8s5TWrWnhra2e0x9+fvTc9QeCR115wxdffnvh3Xt2Dnxq79ps
        eX13bJxbS2fejVs3BlQVxzx3svt3wAAbAy5d7EOPRIXY/ODp5/V3WOUsntm2
        nHvmT//dw8Pz+pjs+pfDy5e+l1TU+F5Bcd12saezJzePgCyl/Nfoouyw+/Pz
        CyjvKbk4P2lrfnHRhz45U990SGj+JSi+iYePWsMzVzzHJ/7nMM8Zt8GgF4G/
        wKRhJoXAZ4BvDmIx1tSSApfwW5E/FvFvxC+RUwl5WSrhEjaus8ClHLv8R+FS
        Q2aqc4HkPD0554ByPQy9VgP9im7QwmVTY62ynwZ0/Lt21AhznE5cxsfvOsxv
        fPcn3rTucT6preELBbePKa89VD9x6Ztpk7bxUaOKeEJhM3d3y+ee7iP4xAlj
        7xZ8g+yBz84gKxWqqa2taW5tba6bNWPyrJr6se/MmTL6S78pd/MqRV9ed+A7
        vmTd+ruZ6Yzt8OGp+fna/sf63rJlSXhN/aTtY7a9wTe9+ROvuvblP57a/5n3
        HXsP+kTMffTpxLl7n6udetm+m3Y+5I7PfPXl5+f/+vO3fSzpInOmNhZVl8ed
        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
        7CkPJMRk8yHCrzN0aJ7okZ3Jowpmc/eJd3O3KffzjPbt3N1zdIcPHvOXYiPz
        D1pxX18yY4xkhwYux6pw6cdMPX3+NbhUYVOrnrePh0dpC+oQaIYrYszIkYQ+
        a2ub9b6dnf1dtOb6ykqu11d08MDHx4+3tZh8tE5OzrymyhiTrmlbbcixQT6A
        1+R7/i4oGf1VRnK5GS6zUw24nCN41rLs0hmIl9RMnVS/Ql8U88yslphPR8zd
        +Y5/6fw35i1ciZpfJ3G29l+z7qY077rNytm/l3tMf/CH51467HDLppZm+d7v
        v/WyDPnfzz9yneHfW3c+FRJRt/pT3dQHuP/sR7n3zId/W736iqbPP3mjB713
        w+WN9fJnl0wvmrdqUX1xbanuQGG6y6t1bc2vBFQuel7Rw0+6jrj7Rb68by+/
        rxNigxDn3Z2aVf5cfOW0z0c31H01NHPqX+7Va7l3/XU8rrD293Hj6pD7m6hQ
        8Jxp9dPEWtitu3JxyvwpBfPqR/rflzE8+EBhcc4b0UWLv6YcAfASvlbMNER/
        bMxgysnIes2jcuUf7s27/9a13cmDo1uNuBQzO8MV299KXFIvzY0SLqFvyPKS
        8n0oD+9fhUsNbHbkHOjcyxsNMyMdyrgZNhW5aW+fd9jHOxAxsY6+3LNntPK9
        t0znM6ZP6OBD49hSM740j683PNbWT+M+47Yack2QQ1fVesWR9MRRf8q4zEwp
        5YH+EfDLYw5s87IlUyEvRyn41gclj3oxNKvu/cIpGw9PXnwd4ijezNxHcEnm
        qMW3eNXfyN2Usz9nyu6HtqydP+29A7ty/vj53QG477tuXlwqr8Pe3Sty33nv
        y0uCyhbtLqibsT+65vJfXRWZgesbv3T3Dx+8fbAjprlny+Rx8mfXXTYSNtd5
        wUlF1yckJ++vnzn75dFT1uN6PzjZPdz/Ej8DLqMjAl+zs7O7wy2s4NOo+Kgv
        MlMDPo+IT/h2aMly7lSzQcHm5r9r6xoha+KDg0NCx1REIlfD0B9/VEnyxtbG
        gg0j8qKfzk4NPpxVqP/TL6DGjJdknwwenPtTeHBqW0tjzaq8olHP+9Rd/z/3
        ltv+9k2aYYbLqLBCa3D5vsAlfFNUf7lQ8JjkJeUVEC4pb73TPL0zjZ3/J2x2
        +GdddfqxpLeosengMOLt4MAo6mtowGXN2Am8bf7lvGHWWr5nm6mua0pbvTku
        m+p5celIntW6nUc3Xvenrulm7tZ+L0eOamZO008yLnPS6nhoYPLGvLz8OfNm
        t167cOZY2CSjSvVNW50Sqj9wSW16Pyyt5jHF/iLdjXJDDLi8Yftuu4ARcz7y
        gtxsvZ1v3HDpGz++rv/r85eqX8E9b79qSq28Bg/ftnXkqBk7F3tOf3hfWPPm
        fXaZk/a5N95kyIl3n/4Qv3LrM0X03mtX1DfQ3/uefWLwFfOr5s1fcFW4Tezo
        54cMb3zBPX/KI+uuvQG+7MdOdg+jPhV6bFpK8vMRicVvecRXfRkRm/55fGL6
        y7Gxcev808d/4VC1lrvUbfm7ZsY22J/xNjaDgme2jbhJ4HJQ5Yio+5PyR+60
        jWv4ry5uNE/Om2Hk5TBzHWjYsIJ3ggIjEWdBnVyZi4vrmLiKhe83zLrh2dDg
        ka/JuIyNtCrfBzFrqou+VuBykcAl5GXpieDyn45NJvWmcXGrrDPNGidsGucu
        u7iWPBwRlgj/Agh83NmvXz/unzuJO5Yt564j1/GmlmkGHtjb2/PdN8znzc1G
        XVanc+dTmlL45Ysred3cTTy2ds1fflkTf/cat417NO/mHg3b/3Z2zDbg0kOX
        zwuzGnhoQPKGmzbPu/s/61O/PfTY+Lcbx09bG1J91Sdh9dd/5JI89u0JbVOR
        B0o98s1wCR6uWrMt0b5w/p/AZlL7jX+/8kD1z+89P+U23PM7L9w75s/fvzDo
        pke+++L8uW2lGxrnrTyUOPeBV2yzp7y09uqdQdljV93jPn4Xd514F3dvvf2r
        F/e/Oxjvf+aurQb78v4bcjbtXp/4yu7tK0d5ZLTcOyx3xkuDkupfqq6fg7Mf
        Nvq2k93DwwYZcZmUkPCSZ2zpN56pjV86+Gd8GhUdu1V5/bLS8jE3OhYu+sth
        9Ia/XMbd8Mfkmctb2ydNyrx61cLGnduWjrx62eiVBdnxLw6NqPzJNqqO20WP
        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
        yb53yYmxz/linIh5g8yH0PkvjfJMrO+ALVapJKnmKgqr3Fe5n2q8LoUL2TyI
        el3nbuj6TdwvfXCs72U/KVtPhuZppbF1BnUn8flqE+hP8lus3/UmyGXReSw6
        76xS1kbSkpSI/XqVKcxlYz5NufrJtZaUNdjoqh/t755Gr1PWUzHnBXyK65jn
        H+mzTYnL3zYO/biM4+3rS4/DTrtz70bZI545wWMmXPNi52+3mbD/Sp/drfON
        uKdyrjiOel+191Tm9JKr2HyInIjSZf3fxYFaTDiHSucz6vsnHrWOpL1fwF2X
        Y45Wi5hCXUCOpbmtHXdjbFznsDDvbKNaA8RsmZItJshdo2wW4d9wPV2uZIO6
        Pn1Sfh5EbV3mCYdf9mkT7o/B/hMal1pf2v72vB15BsbaH2PX+d0eVp80hWfG
        6D4adv4Le45rvwjjrewfgZ/sHcH+EXquOC/kQJUmvI+SqxC3FOKXYvOfKA6E
        axGH+C7uv6HeBibcn8bmrQV2/kpj6Z8Jm9SdWn9qHUqMkuvqfDO9BtabIPfM
        JXx/nRLmtOlr2vls5NSJzrbZGyLi8uCxzHWx+8KQxzK/UtcInZHYtBrfyHOC
        21vm7pB7w95h9xpifgB0KM9D4nNUmsCnzDyki0yQZ6T7Ltl9Xlz7KvdTzVPm
        40JREsWDuBdzD9Z7Ce9dx4pte/+sxeQi8EkdqvUocarnn2uAcnGE6M/ov+V6
        0rltvD5tXub2Jly9N1FDagr7wugeYvWyZux6hGWv3bPG9NxkdfYNEX7Yerk3
        1tvp+gpbV7Ina8qEe4fq/oTnKeGeqm0Tzo2NVa1XabO4uJDNh2xxcSFeU+cz
        FsvdiMTjcs3PWhD7uU1h7gox6sIp51+vgYXKq5Wcp677L9T/mXNKf7Gf21tT
        k/ljByd8Xq1p1nJB17AHBWvdz3gPCmssz01UZf+bfb94BpPDGe7Nrtnn/bOu
        QtcFV5qgLpj1+npNR82VniN7X7XxavMgzYWi+JDNhWwepM8zni93w+kPX2ms
        rBJ82npU49SFV50bWkxeZQmvEXVdnaPk96ft757e6+KFyD2w1jXsNbtnk10r
        pG3M5cblOTW1mS4n7+6cQf6O5rDsb6PvX+8r7F+u87fts8zPsX7a3MeFVXtP
        dXEhzYcWw4UuUNfXe7Adn3pF68glYFTrU41XLTZ2bXH9zTmO6+o8NmKT9d3l
        yVT2+458tn+r1jXWuLYxWVtbaQpzoJetLsEap3NdMUvUk5gch4XfUveM0Lnm
        mofrcynK1L27fFe/ZaLnKWrvs/dOmwuVyoc0D3Lh0JW7EeNx+TA6H25dUsp1
        dD0MfU9+fXdr89xbnXqoZ9sRtbajfJn6rN1l5bLq/s9pb5lb793TP9j32NyU
        uVfuhRwWel1zWDs2ojks9pSo3jbzjW/UnuraP11caD4+NB8PKqoXlzr2r2Qp
        AVeLFnt9m3AfMb9X9fDAeLPL/9PeOvtpU2ijYX0jFsE4JrjsGenXrdb/uVWp
        7O2uvWOofwz30mJKO7ub923vJ0Xvu4Rxn48HLZYLLWhPXul1/UqQ5Rh7C5fn
        mnBP3Ej/D3JoPZ05ZYpz2SoT9Dpctj4x1hqHv+cnEf4ecFJwU9fZ3VrP6z4a
        2me16HsuAadnkgvFWFzDote4CffEBZfdPNR/eNilj/o6Zx4w4Ro/6p6ovmrL
        0pNC369n/17jrB3pnUbuA/QfOewuE81hGbOcj8MuOSdiiXiNcfgKEjWnLi4L
        HlqBnBl77Sers39r6cxithrP3F2SzlT3GqkrcY6ECfw9uBfoccYsbQ6rzzuM
        2kdWpAYxxtwrW+y1bgIum++72dMxk3XppY7W2c+Y4ryQ/h87lrkonan3kChd
        KbERxiyZp651Je+VvJvndtMPuyZ6jMdy9ouFS9b95/2y0HupdPab89iZjGWy
        n5OukY7q3b0Ym+23onywoiuhpxmz1OfKRfmp7Pukbl/2mGsssZQipjiXRXy9
        HHrIpZ9am+e+aMLxQPbDs+OBi+59aMJ87hxvP7jbqSs7Zo6KnmZN13bj9veQ
        w1KvI+c7Kq4T4zKWFRETzWVDZ5W5dKbkAN2kdBGwoGMPOpd9Uf0P9f3V1mZS
        rntQuhL6r8eE+2jb/h72PaEdDN/xGcuDiCWWxYiljxjvJpflWWUVQ/0T212Y
        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
        SyyxxBJLLLHEEkssscQSSyyxxFIo/x/outWB
        "], {{0, 230}, {230, 0}}, {0, 
                       255},
        ColorFunction->RGBColor],
        BoxForm`ImageTag["Byte", ColorSpace -> "RGB", Interleaving -> True],
        Selectable->False],
        BaseStyle->"ImageGraphics",
        ImageSizeRaw->{230, 230},
        PlotRange->{{0, 230}, {0, 230}}]\)}, 1]]], Alignment -> Center, 2]},
            
            {Insert[
              Grid[Partition[
                Join[{Insert[
                   Grid[{{Style["MetaData:", 15, Underlined, Bold], ""}}], 
                   Alignment -> Left, 2]}, {Framed[
                   Insert[Insert[Grid[Transpose[Partition[{
                         
                         (*these are 3 functions that are output together. 
                         The 3 functions are: mzXML file version, MS Run, 
                         and instrument list.*)
                         
                         Insert[
                          Grid[Partition[
                          Join[{MzxmlVersion[overall]}, {MzxmlMSRun[
                          overall]}, {MzxmlMSInstrument[overall]}], 1], 
                          Dividers -> Center], Alignment -> Left, 2],
                         
                         (*these are 2 other functions that are output \
        together. The 2 functions are: 
                         mzXML file instrument software and data processin \
        software.*)
                         
                         Insert[
                          Grid[Partition[
                          Join[{MzxmlInstrumentSoftware[
                          overall]}, {MzxmlDataProcessingSoftware[
                          overall]}], 1], Dividers -> Center], 
                          Alignment -> Left, 2]}, 1]], Dividers -> Center], 
                     Alignment -> Top, 
                     2], {Background -> {None, {{Lighter[
                          Blend[{Lighter[LightBlue, 0.1], 
                          Lighter[LightGray, 0.9]}], .1]}}}, 
                     Spacings -> {2, {2, {0.7}, 0.9}}}, 2]]}], 1]], 
              Alignment -> Left, 2]}], 1]];
        Return[metaDataOverall]
    ]


(* ::Subsubsection:: *)
(* MZXML Back-End *)

(* ::Function:: *)
(* f:MzxmlSpectrumBinaryDecoder *)
(***Function***)
MzxmlSpectrumBinaryDecoder[string_, inputAccuracy_, form_, 
  compression_] :=
    Module[ {
     (*1st input, the user selected spectrum binary data*)
     
     inputString = string,
     (*2nd input, spectrum input bit accuracy*)
     
     binaryFormat = inputAccuracy,
     (*3rd input, data enocoding formate (Base64)*)
     inFormat = form,
     (*4th input, binary data is compressed or not*)
     
     compressionStatus = compression,
     result, str, returnListNoCompression, returnListZLIPCompression, 
     characterList},
        result =
         Switch [
          (*check the spectrum compression status*)
          compressionStatus,
          {},
          (*imports binary data (inputString) in Base64 format from a \
        string, and StringToStream is to open a stream for reading from a \
        string.*)
          
          str = StringToStream[
            ImportString[inputString , {"Base64", "String"}]];
          (*str: is the imported spectrum binary data. binaryFormat: 
          spectrum input bit accuracy. 
          ByteOrdering: 
          is an option for BinaryReadList that specifies the ordering of \
      bytes should be assumed for your computer system*)
          returnListNoCompression = 
           Partition[BinaryReadList[str, binaryFormat, ByteOrdering -> 1], 
            2],
          "none",
          str = 
           StringToStream[
            ImportString[inputString , {"Base64", "String"}]];
          returnListNoCompression = 
           Partition[BinaryReadList[str, binaryFormat, ByteOrdering -> 1], 
            2],
          "zlib",
          returnListZLIPCompression = 
           BinaryReadList[
            StringToStream[ImportString[inputString , inFormat]], "Byte", 
            ByteOrdering -> -1];
          characterList = 
           Partition[
            ImportString[
             FromCharacterCode[
              Developer`RawUncompress[returnListZLIPCompression]], 
             binaryFormat, ByteOrdering -> 1], 2]];
        Return[result]
    ];

(* ::Function:: *)
(* f:MzxmlOffSetListCounter *)
(***Function***)
MzxmlOffSetListCounter[mzXMLFileOfIntrest_] :=
    Module[ {inputFile = mzXMLFileOfIntrest, emptyList, filePath,
      scanNumber, msLevel, retentionTime, basePeakMZ, precursorIntensity,
       centroided, polarity,
      allAttributesOfInterest, getThePosition, mapingIntoEmptyList, 
      masterListOfLists},
     
     (*empty list to append to*)
        emptyList = {};
        filePath = inputFile;
        (*OpenRead will read the data from the file to reaturn an \
      InputStream object*)
        inputFile = OpenRead[filePath];
        
        (*read the entire file and look for an spectrum <scan num> 
        tag for every scan*)
        scanNumber = "<scan num";
        (*read the entire file and look for <mslevel> 
        attribute within each scan*)
        msLevel = "msLevel";
        (*read the entire file and look for <retentionTime> 
        attribute within each scan*)
        retentionTime = "retentionTime";
        (*read the entire file and look for the <basepeakMz> 
        attribute within each scan*)
        basePeakMZ = "basePeakMz";
        (*read the entire file and look for the <precursorIntensity> 
        attribute for every scan*)
        precursorIntensity = "precursorIntensity";
        (*read the entire file and look for the <centroided> 
        attribute for every scan*)
        centroided = "centroided";
        (*read the entire file and look for <polarity> 
        attribute within each scan*)
        polarity = "polarity";
        
        (*set the streaming location to character = 
        1 "at the beginning og the file"*)
        SetStreamPosition[filePath, 1];
        
        (*while loop will take all found location of all given constants \
      above and append them with their locations to an emptyList until the \
      end of the file*)
        While[True,
         allAttributesOfInterest = 
          Find[filePath, {scanNumber, msLevel, retentionTime, basePeakMZ, 
            precursorIntensity, centroided, polarity}];
         getThePosition = StreamPosition[filePath];
         If[ allAttributesOfInterest == EndOfFile,
             Break[],
             Null(**),
             AppendTo[emptyList, {allAttributesOfInterest, getThePosition}]
         ]];
        
        (*MapIndexed applies StringCases to the elements of the emptyList \
      and that is the all detacted attributes {(scanNumber, msLevel, 
        retentionTime, basePeakMZ, precursorIntensity), 
        with each attribute detacted location}.*)
        mapingIntoEmptyList =
         MapIndexed[
          Flatten[{#1[[2]],
             (StringCases[#1[[1]],
               
               {(*detect each spectrum "<scan num=", 
                name it "scan number" and place it in a list with its \
        location*)
                ___ ~~ "<scan num=\"" ~~ 
                  sN : NumberString ~~ ___ ~~ "\"" ~~ ___ -> {"scan number",
                   sN},
                (*detect each spectrum "msLevel=", 
                name it "ms level value" and place it in a list with its \
        location*)
                ___ ~~ "msLevel=\"" ~~ 
                  msL : NumberString ~~ ___ ~~ 
                  "\"" ~~ ___ -> {"ms level value", msL},
                (*detect each spectrum "retentionTime=", 
                name it "scan start time value" and place it in a list with \
        its location*)
                ___ ~~ "retentionTime=\"" ~~ rT : ___ ~~ 
                  "\"" ~~ ___ :> {"scan start time value", 
                  StringReplace[rT, LetterCharacter -> ""]},
                (*detect each spectrum "basePeakMz=", 
                name it "base peak m/z value" and place it in a list with \
        its location*)
                ___ ~~ "basePeakMz=\"" ~~ 
                  bPM : NumberString ~~ ___ ~~ 
                  "\"" ~~ ___ -> {"base peak m/z value", bPM},
                (*if a spectrum contains a precursorMz element, 
                then detect the "precursorIntensity=" as well, 
                name it "selected ion m/z value" and place it in a list \
        with its location*)
                ___ ~~ "precursorIntensity=\"" ~~ ___ ~~
                   ">" ~~ pM : NumberString ~~ 
                  "<" ~~ ___ -> {"selected ion m/z value", pM},
                (*detect each spectrum "polarity=", 
                name it "ionization mode" and place it in a list with its \
        location*)
                ___ ~~ "polarity=\"" ~~ pL : _ ~~ ___ ~~ 
                  "\"" ~~ ___ :> {"ionization mode", 
                  If[ pL == "-",
                      "negative scan",
                      "positive scan"
                  ]},
                (*detect each spectrum "centroided=", 
                
                name it "spectrum representation" and place it in a list \
        with its location*)
                ___ ~~ "centroided=\"" ~~ 
                  cN : NumberString ~~ ___ ~~ 
                  "\"" ~~ ___ :> {"spectrum representation", 
                  If[ cN == "1",
                      "centroid spectrum",
                      "profile spectrum"
                  ]}
                }]), #2}] &, emptyList];
        
        (* 
        Position[mapingIntoEmptyList,"scan number"]: 
        find "scan number" position, 
        this will give you a list with 2 indexes, {#[[1]]-1,#[[1]]}: 
        take the first index - 'scan number' location - 
        and subtract 1 from it, 
        then add the output and the first index in a list,
        Partition[#[[2;;]],2]: 
        partition the Flattened output from the latter from the 2nd index \
      until the end of the list in a nonoverlapping sublist of length 2,
        Take[mapingIntoEmptyList,#]: 
        as you map it into the partitioned list of lists, 
        this will give the partitioned elements of mapingIntoEmptyList
         *)
        masterListOfLists = (Take[
               mapingIntoEmptyList, #]) & /@ (Partition[#[[2 ;;]], 2]) &@
          Flatten[{#[[1]] - 1, #[[1]]} & /@ 
            Position[mapingIntoEmptyList, "scan number"]];
        Close[inputFile];
        Return[masterListOfLists]
    ];

(* ::Function:: *)
(* f:MzxmlSpectrumSelector *)
(***Function***)
MzxmlSpectrumSelector[fileNameI_, mzXMLmsLevel_, mzXMLprecursor_, 
  mzXMLtimeMin_, mzXMLtimeMax_, mzXMLmassMin_, mzXMLmassMax_] :=
    Module[ {mzXMLscanMSLevel = mzXMLmsLevel,
      mzXMLscanPrecursor = mzXMLprecursor,
      (*N[expr]gives a machine\[Hyphen]precision number,
      so long as its magnitude is between $MinMachineNumber and \
    $MaxMachineNumber. $MinMachineNumber is the smallest hardware \
    floating-point number (10^-308), 
      and $MaxMachineNumber is the largest hardware floating point \
    number (10^308).*)
      mzXMLscanTimeMin = N[mzXMLtimeMin],
      mzXMLscanTimeMax = N[mzXMLtimeMax],
      mzXMLscanMassMin = N[mzXMLmassMin],
      mzXMLscanMassMax = N[mzXMLmassMax],
      basePeakOrPrecursorMass, scanLocator, inputFile = fileNameI, 
      offsetLimited, specMin},
     
     (*Searching by Spectrum MS Level*)
     (*"mzXMLscanMSLevel" is the \
   spectrum ms level, it can be 0,1,2,3,.... 
     "mzXMLscanPrecursor" is the spectrum pre-
     cursor meaning that is the parent spectrum for a selected spectrum. 
     "mzXMLscanPrecursor" should always be set to "No precursor-ion" \
   unless mzXMLscanMSLevel is set to 2 or higher.
     
     Using Which statment to evalute the the user inputs. 
     Spectrum with ms level - user input - less than 2, 
     "No precursor-ion" is used, 
     since any ms level less than 2 has no parent spectrum. 
     And we search for them by base peak mass - using user input hint: 
     "base peak m/z value". Spectrum with ms level - user input - 
     greater than or equal to 2, "Yes precursor-ion" is used, instead. 
     Any ms level greater than 2 will have a parent spectrum that it \
   came from 'fragmentation'. We then search by selected ion mass - 
     user input - corresponds to the "selected ion m/z value".*)
        Which[
         (mzXMLscanMSLevel == 0 && 
            mzXMLscanPrecursor == 
             "No precursor-ion") || (mzXMLscanMSLevel == 1 && 
            mzXMLscanPrecursor == 
             "No precursor-ion") || (mzXMLscanMSLevel >= 2 && 
            mzXMLscanPrecursor == "No precursor-ion"),
         basePeakOrPrecursorMass = "base peak m/z value",
         (mzXMLscanMSLevel >= 2 && 
           mzXMLscanPrecursor == "Yes precursor-ion"), 
         basePeakOrPrecursorMass = "selected ion m/z value"];
        
        (*If user sets minimal spectrum mass to 0, 
        replace "mzXMLscanMassMin" to specMin. specMin=10^-300*)
        specMin = 10^-300;
        
        (*Search by scan time only given the MS level preference*)
        If[ mzXMLscanMassMin == 0.,
            mzXMLscanMassMin = specMin
        ];
        Which[
         (*1st searching by scan time: 
         user input for time range while the mass range is unchanging, 
         search by unchanged mass followed by user time input.*)
         
         (*scan 'minimal' time not equal to 0 and scan 'maximal' time not \
        equal to 1000000.)*)
         (! (mzXMLscanTimeMin == 0. && 
              mzXMLscanTimeMax == 
               1000000.)) &&
          (*scan 'minimal' mass equal to 0. and scan \
        'maximal' mass equal to 1000000.*)(mzXMLscanMassMin == 0. && 
            mzXMLscanMassMax == 1000000.),
         
         offsetLimited =
          (*Mapping the "base peak m/z value" positions into the list of \
        offsets*)
          (inputFile[[#]] & /@
            (*find the position of all \
        "base peak m/z value". *)
            ((Position[inputFile, 
                x_ /; Length[x] > 1 && x[[2]] == basePeakOrPrecursorMass, 2,
                 Heads -> False])[[All, 1]]));
         scanLocator =
          (*Mapping the "scan start time value" (these are the selected \
         scans based on user scan time input) positions into the list of \
         offsets*)
          
          offsetLimited[[#]] & /@
           (*find the position of all "scan \
         start time value" given the user mzXMLscanTimeMin and \
         mzXMLscanTimeMax input.*)
           (Position[offsetLimited, 
              x_ /; Length[x] > 1 && 
                x[[2]] == 
                 "scan start time value" && (mzXMLscanTimeMin <= 
                  ToExpression[x[[3]]] <= mzXMLscanTimeMax), 2, 
              Heads -> False][[All, 1]]),
         
         (*2nd searching by scan mass: 
         user input for mass range while the time range is unchanging, 
         search by unchanged time followed by user mass input.*)
         
         (*scan 'minimal' time is set equal to 0 and scan 'maximal' time is \
        set equal to 1000000)*)
         (mzXMLscanTimeMin == 0. && 
            mzXMLscanTimeMax == 
             1000000.) &&
          (*scan 'minimal' mass not equal to 0 and scan \
        'maximal' mass not equal to 1000000*)(! (mzXMLscanMassMin == 0. && 
              mzXMLscanMassMax == 1000000.)),
         
         scanLocator =
          (*Mapping the basePeakOrPrecursorMass (these are the selected \
        scans based on user mass input) positions into the list of offsets*)
        
             inputFile[[#]] & /@
           (*find the position of all \
        basePeakOrPrecursorMass given the user mzXMLscanMassMin and \
        mzXMLscanMassMax input.*)(Position[inputFile, 
              x_ /; Length[x] > 1 && x[[2]] == basePeakOrPrecursorMass && 
                mzXMLscanMassMin <= ToExpression[x[[3]]] <= 
                 mzXMLscanMassMax, 2, Heads -> False][[All, 1]]),
         
         (*3rd searching by both scan mass and time: 
         user input mass range and time range.*)
         
         True,
         
         (*searching by a mass range*)
         offsetLimited = 
          inputFile[[#]] & /@ 
           Position[inputFile, 
             x_ /; Length[x] > 1 && 
               x[[2]] == 
                basePeakOrPrecursorMass && (mzXMLscanMassMin <= 
                 ToExpression[x[[3]]] <= mzXMLscanMassMax), 2, 
             Heads -> False][[All, 1]];
         (*searching the outcome of the mass range search, 
         searching by a time range*)
         scanLocator =
          offsetLimited[[#]] & /@ 
           Position[offsetLimited, 
             x_ /; Length[x] > 1 && 
               x[[2]] == 
                "scan start time value" && (mzXMLscanTimeMin <= 
                 ToExpression[x[[3]]] <= mzXMLscanTimeMax), 2, 
             Heads -> False][[All, 1]]];
        If[ mzXMLscanMSLevel == 0,
            Null,
            scanLocator =
             (*Mapping the MS Level (these are the selected scans based on \
            user MS Level input) positions into the scanLocator*)
             
             scanLocator[[#]] & /@
              (*find the position of all \
            spectraMSLevel given the user scan ms level input.*)(Position[
                 scanLocator, 
                 x_ /; Length[x] > 1 && x[[2]] == "ms level value" && 
                   ToExpression[x[[3]]] == mzXMLscanMSLevel, 2, 
                 Heads -> False][[All, 1]])
        ];
        Return[scanLocator]
    ];

(* ::Function:: *)
(* f:MzxmlGetSpectraMSLevel *)
(***Function***)
MzxmlGetSpectraMSLevel[spectramIndexList_] :=
    Module[ {outPutMSLevel,
      msLevel = spectramIndexList},
     (*find "ms level value" in a selected scan using Cases*)
        outPutMSLevel = 
         Cases[msLevel, x_ /; x[[2]] == "ms level value" :> x[[3]], 1];
        (*returing the ms level, e.g. {"1"}, {"2"},...*)
        Return[outPutMSLevel]
    ]

(* ::Function:: *)
(* f:MzxmlSpectrumMassTimeAssociation *)
(***Function***)
MzxmlSpectrumMassTimeAssociation[listOfselectedSpectrums_, 
  selectedIonMassOrBasePeakMass_] :=
    Module[ {selectInput,
      inputSelectedSpectrums = listOfselectedSpectrums,
      TrueOrFalseForInput = selectedIonMassOrBasePeakMass},
        selectInput =
         Association @@
          ({#[[2]] -> 
               Flatten[ToExpression[#] & /@ #[[
                  1]]]} & /@
            (*Transpose the above three together.*)
        
                 Transpose[
             {
              {(*Output either "selected ion m/z value" or "base peak m/z \
        value" masses(x[[3]]) if True or False respectively.*)
                 
                 Cases[#, x_ /; x[[2]] ==
                     (*User will be able to navigate the spectra by \
        searching using "selected ion m/z value" or "base peak m/z value". 
                     If click "selected ion m/z value" in viewer interface, 
                     then is True. If click "base peak m/z value", 
                     then that is false.*)
                     
                     If[ TrueOrFalseForInput == True,
                         "selected ion m/z value",
                         "base peak m/z value"
                     ] :> 
                   x[[3]], 1],
                 (*Output the "scan start time value"*)
                 Cases[#, x_ /; x[[2]] == "scan start time value" :> x[[3]],
                   1],
                 (*Output the "ms level value"*)
                 
                 Cases[#, x_ /; x[[2]] == "ms level value" :> x[[3]], 
                  1]} & /@ inputSelectedSpectrums, 
              inputSelectedSpectrums}])
    ];

(* ::Function:: *)
(* f:MzxmlGetPrecursorSpectraAtMSLevel *)
(***Function***)
MzxmlGetPrecursorSpectraAtMSLevel[mzXMLOffSetList_, 
  SelectedIonMassOfInterest_] :=
    Module[ {
      (*offSetList is equal to the 1st input to this function and that \
    is "mzXMLoffsetList"*)
      offSetList = mzXMLOffSetList,
      (*findMS1 is equal to the 2nd input to this function and that is \
    the user selected spectrum*)
      findMS1 = SelectedIonMassOfInterest,
      findingNumberOfSelectedIoonMassList, positionInTheOffSetList, 
      currentMSLevel, i,
      (*empty list*)
      positionsReturned = {},
      PrecursorScans},
     
     (*first check to see you have a <selected ion m/z value> 
     in the user selected spectrum. If you do, 
     that mean is a scan that has an MSLevel>
     1. We then extract the 4th element in the list where the <
     selected ion m/z value> 
     location and that is the number of the list.*)
        findingNumberOfSelectedIoonMassList =
         If[
          (*use Cases to find the "selected ion m/z value" tag in the user \
      selected spectrum*) Cases[findMS1, x_ /; x[[2]] == "selected ion m/z value", 1][[1, 
            2]] ==
           "selected ion m/z value",
             Cases[findMS1, x_ /; x[[2]] == "selected ion m/z value" :> x[[4]],
                1][[1]]
         ];
        
        (*take the number of the list where the <selected ion m/z value> 
        is located as an input, and get the location in the offsetList = <
        mzXMLOffsetList>*)
        positionInTheOffSetList = 
         Flatten[Position[offSetList, 
            x_ /; x[[4]] == findingNumberOfSelectedIoonMassList, 2, 
            Heads -> False]][[1]];
        
        (*list positionInTheOffSetList*)
        positionsReturned = {positionInTheOffSetList};
        currentMSLevel =
         (*gives the expression obtained by interpreting strings as Wolfram \
        Language input and subtract 1 from it*)
         ToExpression[
           (*find the msLevel of the selected spectrum*)
           
           Cases[findMS1, x_ /; x[[2]] == "ms level value" :> x[[3]], 1][[
            1]]] - 1;
        
        
        (*For[start,test,incr,body] executes start,
        then repeatedly evaluates body and incr until test fails to give \
      True.*)
        For[
         (*i equal to one less of the position of the selected spectrum in \
        the mzXMLOffsetList*)
         i = positionInTheOffSetList - 1,
         (*test if i is greater than 0, 
         this need to be True and once is not, the For loop will stop*)
         
         i > 0,
         (*keep subtracting 1 from i *)
         i--,
         (*if currentMSLevel which is one less the selected spectrum \
        msLevel is equal to the spectrum that we keep subtracting 1 from the \
        initial spectrum position, then append to the positionsReturned, 
         otherwise break.*)
         If[ currentMSLevel ==
           (*find the msLevel of each spectrum as we keep subtracting 1 \
        from the initial spectrum position*)
           Cases[
             (*indexing each i value*) 
             offSetList[[i]], 
             x_ /; x[[2]] == "ms level value" :> ToExpression[x[[3]]] , 1][[
            1]],
          (*appends i to the value of positionsReturned, 
          and resets positionsReturned to the result.*)
             AppendTo[positionsReturned, i];
             (*keep subtracting 1 from currentMSLevel *)
             currentMSLevel--,
             (*break once MSLevel is equal to 0*)
             If[ currentMSLevel == 0,
                 Break[]
             ]
         ]];
        
        (*output list of found spectra including the selected spectrum*)
        PrecursorScans =
         offSetList[[positionsReturned]];
        Return[PrecursorScans]
    ];

(* ::Function:: *)
(* f:MzxmlGetPrecursorSpectraAtMSLevel *)
(***Function***)
MzxmlGetMSSpectraLocationAndBinary[fileNameII_, 
  selectedSpectraOfIntrestII_] :=
    Module[ {
      streamingFile,
      (*inputFile is equal to the 1st input to this function and that is \
    "mzXMLfilePath"*)
      inputFile = fileNameII,
      (*findMS is equal to the 2st input to this function and that is \
    the user selected spectrum*)
      
      findMS = selectedSpectraOfIntrestII[[1]],
      whatIReadIn, emptyList, joiningListToString},
     
     (*open an mzXML file that is set equal to <mzXMLfilePath> 
     and read data from it,and return an input stream object.*)
        streamingFile = OpenRead[inputFile];
        (*set the steaming position to the location of the selected \
      spectrum index location that was past to this function as a second \
      input*)
        SetStreamPosition[inputFile, findMS];
        
        (*list with the spectrum closing tag </scan> 
        in order to convert the selected spectrum to mathematica's XML \
      format*)
        whatIReadIn = "</scan>";
        emptyList = {"<scan "};
        
        (*While test=True, evaluate the body*)
        While[True, whatIReadIn = Read[streamingFile, Record];
                    AppendTo[emptyList, whatIReadIn];
                    (*if statment to break the loop once spectrum is found*)
                    If[ StringTake[StringTrim[whatIReadIn], -8] == "</peaks>",
                        AppendTo[emptyList, "</scan>"];
                        Break[]
                    ]];
        
        (*import the scan to mathematica XML format*)
        joiningListToString = ImportString[StringJoin[emptyList], "XML"];
        Close[inputFile];
        Return[joiningListToString]
    ];

(* ::Function:: *)
(* f:MzxmlSpectrumPlotter *)
(***Function***)
MzxmlSpectrumPlotter[spectrum_, userInputForPeakLabeling_] :=
    Module[ {peaksElement,
      binaryDATAArray = spectrum,
      getTheRealBits,
      binaryDATAString,
      compresiontype,
      decodeBinaryData,
      trueOrFalseInput = userInputForPeakLabeling,
      final},
     (*
     Constant peaksElement serve as a storage for element <peaks> 
     in a selected scan. Element <peaks> 
     is where the actual data encoded using base64. 
     Byte order must be network. The order of the pairs must be m/z-
     intensity.;
     
     Cases here takes 3 element: 
     [1st element in an expiration is spectrum2 - selected scan, 
     2nd element is XMLElement pattern matching and that is[peaks], 
     3rd element is levelspece and that is Infinity - 
     all levels (except 0) - 
     to give a list of all parts of the expiration "1st" on this \
   specified level that match the pattern "2nd"]
     *)
        peaksElement = 
         Cases[binaryDATAArray, XMLElement["peaks", _, _], Infinity];
        (*
        peaksElement[[1,
        2]] will give you a list of three required attributes e.g. \
      {"precision"\[Rule]"32","byteOrder"\[Rule]"network",
        "pairOrder"\[Rule]"m/z-int"},
        
        precision: number of bits used by each component e.g. 32 bit, 64 bit,
        byteOrder: 
        is FIXED byte order of the encoded binary information (must be \
      network),
        pairOrder: is FIXED order of the m/z intensity pairs (must be m/z-
        int).
          *)
        (*precision: number of bits used by each component e.g. 
        32 bit, 64 bit. Take the number and join it to "Real"*)
        getTheRealBits = 
         "Real" <> 
          Cases[peaksElement[[1, 2]], x_ /; x[[1]] == "precision", 1][[1, 
           2]];
        (*Actual data encoded using base64*)
        binaryDATAString = peaksElement[[1, 3, 1]];
        
        (*constant compresiontype checks if the scan encoded data if \
      compresed or not. 
        If attribute compressionType is not found then just output an \
      empity list {}, but if attribute compressionType is found, 
        then use a Switch to check if equal to "none" or "zlip". 
        Pass compression states to the decoder.*)
        compresiontype = 
         If[ Cases[peaksElement[[1, 2]], 
            x_ /; x[[1]] == "compressionType" :> x[[2]], 1, 
            Heads -> False] == {},
             {},
             Switch[Cases[peaksElement[[1, 2]], 
               x_ /; x[[1]] == "compressionType" :> x[[2]], 1, 
               Heads -> False], {"none"}, "none", {"zlib"}, "zlib"]
         ];
        
        (*Decode the scan binary data by calling the \
      MzxmlSpectrumBinaryDecoder.*)
        decodeBinaryData = 
         MzxmlSpectrumBinaryDecoder[binaryDATAString, 
          getTheRealBits, {"Base64", "String"}, compresiontype];
        final =
         (*If user select True, then the plot will have a labeled peaks, 
         peaks are not labeled otherwise*)
         
         If[ trueOrFalseInput == False,
             Panel[Show[
             (*ListLinePlot will line plot the decodAndTransposeBinaryData \
             values.*)
             ListLinePlot[decodeBinaryData,
             {(*Using Full for the Plot Range to include the full range of \
             the original data*)
             PlotRange -> Full,
             (*m/
             z array" is the x-axis and "intensity array" is the y-axis.*)
             
             AxesLabel -> {"m/z array", "Intensity array"},
             LabelStyle -> Directive[Lighter[Black, 0.1]],
             AxesStyle -> Lighter[Black, 0.1]
             }, ImageSize -> 700],
             
             ListPlot[
             (*Tooltip displays labels while the mouse pointer is in the \
             area where the expression (m/z and intensity) is displayed.*)
             
             Tooltip[decodeBinaryData],
             {
             PlotRange -> Full,
             PlotMarkers -> {"\[FilledCircle]", 0.01},
             AxesLabel -> {"m/z array", "Intensity array"},
             LabelStyle -> Directive[Lighter[Black, 0.1]],
             AxesStyle -> Lighter[Black, 0.1]
             }, ImageSize -> 700]]],
    (*If False, peaks are not labeled.*)
             Panel[Show[
               (*ListLinePlot will line plot the decodAndTransposeBinaryData \
             values.*)
               ListLinePlot[decodeBinaryData,
                {(*Using Full for the Plot Range to include the full range of \
             the original data*)
                 PlotRange -> Full,
                 (*m/
                 z array" is the x-axis and "intensity array" is the y-axis.*)
             
                         AxesLabel -> {"m/z array", "Intensity array"},
                 LabelStyle -> Directive[Lighter[Black, 0.1]],
                 AxesStyle -> Lighter[Black, 0.1]
                 }, ImageSize -> 700]
               ]]
         ];
        Return[final]
    ];


(* ::Subsubsection:: *)
(* MZXML Interpretation *)

(* ::Function:: *)
(* f:MZXMLInterpretationFunction *)
(***Function***)
MZXMLInterpretationFunction[filePath_] :=
    DynamicModule[{
      mzXMLOverallRaw, outputMetaData, mzXMLoffsetList, metaData,
      (*The first and only input is set equal to the inputfile and that \
    is the user uploaded file path*)
      inputfile = filePath,
      interpretationOutcome, mzXMLmsLevel, mzXMLpreCursorSwitch, 
      mzXMLtimeMin, mzXMLtimeMax, mzXMLmassMin, mzXMLmassMax, 
      mzXMLprecursorSpectrum, mzXMLOptionForPlotLabeling, mzXMLpreCursor},
     
     (**********Constants "mzXML Back-End"*********)
     
     (*mzXMLOverallRaw is equal to function <MzxmlRawMetaData> 
     that takes the input which is the user uploaded file path. 
     The output is the raw metadata.*)
     
     mzXMLOverallRaw = MzxmlRawMetaData[inputfile];
     (*outputMetaData is equal to function MzxmlMetaData that takes the \
   raw metadata and it output a readable metadata for the user*)
     outputMetaData = MzxmlMetaData[mzXMLOverallRaw];
     
     
     (*mzXMLoffsetList is a constant that is set equel to function <
     MzxmlOffSetListCounter>that takes mzXML file path. 
     MzxmlOffSetListCounter output a list of offsets for {spectrum \
   index, ms level, base peak m/z, scan start time, precursor, 
     selected ion m/z}. 
     This list is used as an input for downstream functions.*)
     mzXMLoffsetList = MzxmlOffSetListCounter[inputfile];
     
     (*
     1: {Define variables by assignment};
     2: how it appear;
     3: execute
     *)
     
     
     (***Returning mzXML Front-End Using Interpretation***)
     
     
     (*Using interpretation to display a version that allow for user \
   interaction, but is interpreted as the unevaluated form*)
     interpretationOutcome = Interpretation[
       
       
       (*These are the local defined variables by assignment*)
       {
        
        (*file metadata will always show by default*)
        
        metaData = True,
        (*"MS Level" variable*)
        mzXMLmsLevel = 0,
        (*mzXMLpreCursorSwitch is a localy defind variable to allow for \
     user "Search By: " method. 
        Search by "Selected Ion m/z Value" when mzXMLpreCursorSwitch = 
        True and by "Base Peak m/z Value" when mzXMLpreCursorSwitch = 
        False*)
        mzXMLpreCursorSwitch = False,
        mzXMLtimeMin = 0,
        mzXMLtimeMax = Infinity,
        mzXMLmassMin = 0,
        mzXMLmassMax = Infinity},
       
       
       
       
       (*This is how it appears*)
       Panel[
        Grid[{
          
          {Panel[Style["Input File Path: " <> inputfile, 11]]}, {Panel[
            Grid[{
              
              (*This is one of the defined variables that allow for an \
     interactive way for the user to search for different spectra "MS \
     Level" by using PopupMenu and is set by default to show all MS Levels \
     in the uploaded file*)
              {"MS Level:" PopupMenu[
                 Dynamic[mzXMLmsLevel], 
                 Join[{0 -> Style["All", 11]}, 
                  ToExpression[#] & /@ 
                   Union[Cases[mzXMLoffsetList, 
                     x_ /; (Length[x] > 1) && (x[[2]] == 
                       "ms level value") :> x[[3]], 2]]]]},
              
              (*This is one of the the defined variables that allow for \
     an interactive way for the user to search the spectra by "Selected \
     Ion m/z Value" or "Base Peak m/z Value" by using RadioButton and is \
     set by default to search by "Base Peak m/z Value"*)
              {Style[
                "Search By: ", 10], 
               RadioButtonBar[
                Dynamic[
                 mzXMLpreCursorSwitch], {True -> 
                  "Selected Ion m/z Value", 
                 False -> "Base Peak m/z Value"}, 
                Appearance -> "Horizontal" , 
                BaselinePosition -> Automatic, Appearance -> Small]},
              
              (*This is one of the defined variables that allow for an \
     interactive way for the user to change/
              input "minimum retention time" viewing MS spectrum by \
     using InputField and is set to search starting at time=
              0*)
              {Style[
                " Input \"minimum [\[GreaterEqual] 0]\" retention time \
\"second\": ", 11, Darker[Green, 0.4]], 
               InputField[Dynamic[mzXMLtimeMin], FieldSize -> 5]},
              
              (*This is one of the defined variables that allow for an \
     interactive way for the user to change/
              input "maximum retention time" viewing MS spectra by using \
     InputField and is set to time=
              Infinity*)
              {Style[
                " Input \"maximum\" retention time \"second\": ", 11, 
                Darker[Green, 0.4]], 
               InputField[Dynamic[mzXMLtimeMax], FieldSize -> 5]},
              
              (*This is one of the defined variables that allow for an \
     interactive way for the user to change/
              input "minimum spectra mass" by using InputField and is \
     set to mass=
              0*)
              {Style[
                " Input \"minimum [\[GreaterEqual] 0]\" spectra mass \
\"m/z\" : ", 11, Lighter[Blue, 0.2]], 
               InputField[Dynamic[mzXMLmassMin], FieldSize -> 5]},
              
              (*This is one of the defined variables that allow for an \
     interactive way for the user to change/
              input "maximum spectra mass" by using InputField and is \
     set to mass=
              Infinity*)
              {Style[
                " Input \"maximum\" spectra mass \"m/z\" : ", 11, 
                Lighter[Blue, 0.2]], 
               InputField[Dynamic[mzXMLmassMax], FieldSize -> 5]}
              
              }, Alignment -> Left]]},
          
          (*This is the metadata in a panel form*)
          {Panel[
            Grid[{{Style["View File Metadata: ", 10], 
               RadioButtonBar[Dynamic[metaData],
                {True -> "Yes", False -> "No"} , 
                BaselinePosition -> Automatic, 
                Appearance -> Small]}}]]},
          {Panel[
            Grid[{{Button[Style["Evaluate", 12], 
                SelectionMove[EvaluationCell[], All, CellContents];
                SelectionEvaluateCreateCell[EvaluationNotebook[]](*,
                Method\[Rule]"Queued"*)]}}]]}},
         
         Alignment -> Left]],
       
       
       
       
       (*This is how it executes*)
       
       
       (***Returning mzXML Front-End Using Interpretation***)
       
       Interpretation[
        
        (*These are another set of locally defind variables*)
        \
     {mzXMLpreCursor,
         mzXMLprecursorSpectrum,
         mzXMLOptionForPlotLabeling},
        
        
        
        (*This is how it appear*)
        Panel[Grid[{{
            (*Output file metadate once metAData\[Equal]True*)
            
            If[ metaData == True,
                Panel[outputMetaData]
            ]},
           
           
           
           {Panel[
             Grid[{
               {Grid[
                 Transpose[
                  
                  (*This is one of the variables that allow for an \
     interactive way for the user to search the spectra by "Selected Ion \
     m/z Value" or "Base Peak m/z Value" by using RadioButton*)
              \
       {{Style[
                     "Please select: {" ~~ 
                      If[ mzXMLpreCursorSwitch == True,
                          "Selected Ion m/z Value",
                          "Base Peak m/z Value"
                      ] ~~
                       ", retention time, MS Level } : ", 10]},
                   
                   
                   (*This is one of the varables that allow for an \
     interactive way for the user to select a spectrum using either \
     "Selected Ion m/z Value" or "Base Peak m/z Value", max/
                   
                   min retention time and MS Level of intrest using \
     PopupMenu*)
                   {PopupMenu[
                     Dynamic[mzXMLpreCursor],
                     Sort[
                      MzxmlSpectrumMassTimeAssociation[
                       Which[
                       mzXMLpreCursorSwitch == True,
                       (MzxmlSpectrumSelector[mzXMLoffsetList, 
                       mzXMLmsLevel, "Yes precursor-ion",
                       Replace[mzXMLtimeMin, Null -> 0],
                       
                       Replace[
                       mzXMLtimeMax, {Infinity -> 1000000, 
                       Null -> 1000000}],
                       Replace[mzXMLmassMin, Null -> 0],
                       
                       Replace[ 
                       mzXMLmassMax, {Infinity -> 1000000, 
                       Null -> 1000000}]]),
                       mzXMLpreCursorSwitch == False,
                       (MzxmlSpectrumSelector[mzXMLoffsetList, 
                       mzXMLmsLevel, "No precursor-ion",
                       Replace[mzXMLtimeMin, Null -> 0],
                       
                       Replace[
                       mzXMLtimeMax, {Infinity -> 1000000, 
                       Null -> 1000000}],
                       Replace[mzXMLmassMin, Null -> 0],
                       
                       Replace[ 
                       mzXMLmassMax, {Infinity -> 1000000, 
                       Null -> 1000000}]])], mzXMLpreCursorSwitch]]]}}], 
                 Alignment -> Left]},
               
               
               {Panel[
                 Grid[
                  {
                   (*This is one of the variables that allow for an \
     interactive way for the user to locate the "precursor spectrum" for \
     any selected "daughter spectrum" by using RadioButtonBar*)
              \
        {Style["Find precursor spectrum for MS Level > 1: ", 10], 
                    RadioButtonBar[
                     Dynamic[mzXMLprecursorSpectrum], {True -> "No", 
                      False -> "Yes"} , Appearance -> Tiny]},
                   
                   (*This is one of the variables that allow for an \
     interactive way for the user to "Peaks label" all peaks in selected \
     spectrum by using RadioButton*)
                   {Style[
                     "Peaks Labeling: ", 10], 
                    
                    RadioButtonBar[
                     Dynamic[mzXMLOptionForPlotLabeling], {True -> "No", 
                      False -> "Yes"} , Appearance -> Tiny]}}, 
                  Alignment -> Left]]}},
              
              Alignment -> Left]]},
           {Panel[
             Grid[{{Button[Style["Evaluate", 12], 
                 SelectionMove[EvaluationCell[], All, CellContents];
                 SelectionEvaluateCreateCell[EvaluationNotebook[]](*,
                 Method\[Rule]"Queued"*)]}}]]}
           }, Alignment -> Left]],
        (*This is how it executes*)
        Panel[Grid[{
           (*Null since no chromatogram is used*)
           {Null},
           {
            If[ Length[mzXMLpreCursor] < 1,
                Style["Please Make a Selection....", 10],
                If[
                 (*using function <MzxmlGetSpectraMSLevel> 
                 to find the selected spectrum MS level. If the MS Level = 
                 1 then plot the spectrum using function <
                 MzxmlSpectrumPlotter> 
                 and creat a spectrum note card label.*) MzxmlGetSpectraMSLevel[mzXMLpreCursor] == {"1"},
                    Grid[
                     {
                      
                      (*For labeling the scan by framming the scan "ms level \
                    value", "scan start time value", 
                      "base peak m/z value".*)
                      {Insert[
                        Grid[Transpose[
                          If[ Flatten[
                             DeleteCases[(Cases[#, 
                             x_ /; Length[x] > 1 && 
                             x[[1]] == "spectrum representation" :> x[[2]], 
                             All] & /@ 
                             mzXMLpreCursor[[2 ;;, {2, 3}]]), {}]] == {},
                              Insert[mzXMLpreCursor[[2 ;;, {2, 
                                3}]], {"spectrum representation", 
                                "profile spectrum"}, 1],
                              mzXMLpreCursor[[2 ;;, {2, 
                                3}]]
                          ]]], {Background -> {None, {White, {Lighter[
                             Blend[{Lighter[LightBlue, 0.1], 
                             Lighter[LightGray, 0.9]}], .1]}}}, 
                         Dividers -> {Black, {2 -> Black}}, Frame -> True, 
                         Spacings -> {2, {2, {0.7}, 0.9}}}, 2]},
                      
                      (*If MSLevel \[Equal] 1, 
                      get the scan location and pass it to the \
                    MzxmlGetMSSpectraLocationAndBinary function for binary data \
                    collection. 
                      Call MzxmlSpectrumPlotter to decode and plot the scan \
                    binary data.*)
                      {MzxmlSpectrumPlotter[
                        MzxmlGetMSSpectraLocationAndBinary[inputfile, 
                         Dynamic[
                          Cases[mzXMLpreCursor, 
                            x_ /; x[[2]] == "scan number" :> x[[1]], 1][[1]]]],
                         mzXMLOptionForPlotLabeling]}
                      },
                     Dividers -> {{{True, True}}, {{True, False}}}],
                     
                    
                    (*If MSLevel \[NotEqual] 1, 
                    use a switch and if False is selected then pass "No \
         precursor-ion" but if True is selected then pass "Yes precursor-ion" \
         to get the precursor scan AKA. parant scan.*)
                    Switch[mzXMLprecursorSpectrum,
                     
                     
                     (*If user does want to see the "precursor spectrum" and \
                    True is selected then pass "Yes precursor-ion" to get the "precursor \
                    spectrum" for any selected "daughter spectrum". To do so, 
                     we call MzxmlGetMSSpectraLocationAndBinary function for \
                    precursor scans locations.*)
                     True,
                     Grid[
                      {
                       
                       {Insert[
                         Grid[Transpose[
                           If[ Flatten[
                             DeleteCases[(Cases[#, 
                             x_ /; Length[x] > 1 && 
                             x[[1]] == "spectrum representation" :> x[[2]], 
                             All] & /@ 
                             mzXMLpreCursor[[2 ;;, {2, 3}]]), {}]] == {},
                               Insert[
                                mzXMLpreCursor[[2 ;;, {2, 
                                3}]], {"spectrum representation", 
                                "profile spectrum"}, 1],
                               mzXMLpreCursor[[2 ;;, {2, 
                                3}]]
                           ]]], {Background -> {None, {White, {Lighter[
                             Blend[{Lighter[LightBlue, 0.1], 
                             Lighter[LightGray, 0.9]}], .1]}}}, 
                          Dividers -> {Black, {2 -> Black}}, Frame -> True, 
                          Spacings -> {2, {2, {0.7}, 0.9}}}, 2]},
                       
                       {MzxmlSpectrumPlotter[
                         MzxmlGetMSSpectraLocationAndBinary[inputfile, 
                          Dynamic[
                           Cases[mzXMLpreCursor, 
                             x_ /; x[[2]] == "scan number" :> x[[1]], 
                             1][[1]]]], mzXMLOptionForPlotLabeling]}
                       },
                      Dividers -> {{{True, True}}, {{True, False}}}],
                     
                     
                     False,
                     (*If MSLevel \[NotEqual] 1, 
                     
                     use a switch and if False is selected then pass "No \
                    precursor-ion" but if True is selected then pass "Yes precursor-ion" \
                    to get the precursor scan AKA. parant scan. To do so, 
                     we call mzXMLGetPrecursorSpectraAtMSLevel function for \
                    precursor scans locations. 
                     We then map into it the binary data and call \
                    mzXMLstreamingAndBinaryAndspectraPlotting to decode and plot the scan \
                    binary data.*)
                     
                     Grid[Partition[
                       Flatten[
                        Transpose[{Insert[
                             Grid[Transpose[#]], {Background -> {None, {White, \
                    {Lighter[Blend[{Lighter[LightBlue, 0.1], 
                             Lighter[LightGray, 0.9]}], .1]}}}, 
                             Dividers -> {Black, {2 -> Black}}, Frame -> True, 
                             Spacings -> {2, {2, {0.7}, 0.9}}}, 
                             2] & /@ (If[ Flatten[
                             DeleteCases[(Cases[#, 
                             x_ /; Length[x] > 1 && 
                             x[[1]] == "spectrum representation" :> x[[2]], 
                             All] & /@ #[[2 ;;, {2, 3}]]), {}]] == {},
                                          Insert[#[[2 ;;, {2, 
                                          3}]], {"spectrum representation", 
                                          "profile spectrum"}, 
                                          1],
                                          #[[2 ;;, {2, 
                                          3}]]
                                      ] & /@ (MzxmlGetPrecursorSpectraAtMSLevel[
                             mzXMLoffsetList, 
                             mzXMLpreCursor])), (MzxmlSpectrumPlotter[#, 
                             mzXMLOptionForPlotLabeling] & /@ \
                    (MzxmlGetMSSpectraLocationAndBinary[
                             inputfile, #[[1]]] & /@ \
                    (MzxmlGetPrecursorSpectraAtMSLevel[mzXMLoffsetList, 
                             mzXMLpreCursor])))}]], 1],
                      
                      Dividers -> {{{True, True}}, {{True, False}}}]
                     ]
                ]
            ]}}, Alignment -> Center]]]
       ]];

(* ::Subsection:: *)
(* MZXML Or MzML *)

(* ::Function:: *)
(* f:MZXMLOrMZML *)
(***Function***)
MZXMLOrMZML[XMLFile_] :=
    Module[ {
      (*inputFile is equal to the only input to this function*)
      
      inputFile = XMLFile,
      emptyList, filePath, allMetaData, fileMZXMLorMZMLOutPut},
        emptyList = {};
        filePath = inputFile;
        (*open and read file*)
        inputFile = OpenRead[filePath];
        
        (*Set the the streaming position to 0. This will stream metadata \
      from the begining of the file*)
        SetStreamPosition[filePath, 0];
        While[True,
         allMetaData = Read[filePath, String];
         (*Using <StringCases> 
         to search "msRun" for mzXML file or "cvList" for mzML files. <..> 
         is to keep repeating the pattern matching by representing "msRun \
      or cvList" until it matches, if it doesn't match, keep appending <
         AppendTo> to the emptyList={}. Using <If> and <Break[]> 
         to stop while loop once "msRun or cvList" matchs.*)
         If[ StringCases[allMetaData, WordCharacter ..][[1]] == "msRun" || 
           StringCases[allMetaData, WordCharacter ..][[1]] == "cvList",
             Break[],
             AppendTo[emptyList, allMetaData]
         ]];
        
        (*If either mzML or mzXML if found in emptyList, 
        then use Which to output mzML if mzML is found and output mzXML if \
      mzXML is found, instead.*)
        fileMZXMLorMZMLOutPut =
         If[ Union[
             Flatten[
              StringCases[
               emptyList, ___ ~~ "mzML" ~~ ___ :> {"mzML"}]]] == {"mzML"} ||
        
                Union[
             Flatten[
              StringCases[
               emptyList, ___ ~~ "mzXML" ~~ ___ :> "mzXML"]]] == {"mzXML"},
             Which[
              Union[
                Flatten[
                 StringCases[
                  emptyList, ___ ~~ "mzML" ~~ ___ :> {"mzML"}]]] == {"mzML"}, 
              "mzML",
              Union[
                Flatten[
                 StringCases[
                  emptyList, ___ ~~ "mzXML" ~~ ___ :> "mzXML"]]] == {"mzXML"}, 
              "mzXML"]
         ];
        Close[inputFile];
        Return[fileMZXMLorMZMLOutPut]
    ];


(* ::Subsection:: *)
(* MathIOmicaMSViewer *)

(* ::Function:: *)
(* f:MSViewer *)
(***Options***)
Options[MSViewer] = {MathIOmicaDataDirectory -> 
    If[ ! DirectoryQ[
       FileNameJoin[
        Flatten[{FileNameSplit[$UserBaseDirectory], "Applications", 
          "MathIOmica", "MathIOmicaData"}]]],
        CreateDirectory[
         FileNameJoin[
          Flatten[{FileNameSplit[$UserBaseDirectory], "Applications", 
            "MathIOmica", "MathIOmicaData"}]]],
        FileNameJoin[
         Flatten[{FileNameSplit[$UserBaseDirectory], "Applications", 
           "MathIOmica", "MathIOmicaData"}]]
    ],
   PSImsOBOFile -> "psi-ms.obo.txt",
   VocabularyVariable -> None};
(***Function***)
MSViewer[inputfile_:None,OptionsPattern[]] :=
    Module[ {dir = OptionValue[MathIOmicaDataDirectory],
      oboFile = OptionValue[PSImsOBOFile],
      vocabConstant = OptionValue[VocabularyVariable],
      inputFile = inputfile,
      mzMLORmzXMLfilePath,
      finalOutPut,
      fileMSOBO,
      mzmlVocab},
      (**** mzMLORmzXMLfilePath is set to the user uploaded file path ****)
        mzMLORmzXMLfilePath = If[ MatchQ[inputFile,None],
                                  SystemDialogInput["FileOpen"(*, SetDirectory[NotebookDirectory[]]*)],
                                  inputFile
                              ];
        
        (*using Switch to check the file format by calling function \
     MZXMLOrMZML and if "mzML" file, 
        call MZMLInterpretationFunction and if "mzXML" file, 
        call MZXMLInterpretationFunction.*)
        If[ mzMLORmzXMLfilePath =!= $Canceled && 
          mzMLORmzXMLfilePath =!= $Failed,
            finalOutPut = 
             Switch[MZXMLOrMZML[mzMLORmzXMLfilePath], "mzML", 
              If[ MatchQ[vocabConstant, None],
                  fileMSOBO = FileNameJoin[Flatten[{dir, oboFile}]];
                  If[ FileExistsQ[fileMSOBO],
                      mzmlVocab = Import[fileMSOBO, "Lines"],
                      Return["Could not find Control Vocabulary File at: " <> 
                        fileMSOBO <> 
                        ". Please place at MathIOmicaData location or provide a \
VocabularyConstant variable."]
                  ],
                  mzmlVocab = vocabConstant
              ];
              MZMLInterpretationFunction[mzMLORmzXMLfilePath, mzmlVocab] , 
              "mzXML", MZXMLInterpretationFunction[mzMLORmzXMLfilePath]],
            "File Selection was Cancelled"
        ]
    ]

(* ::Section:: *)
(*#####SpectralAnalysis#####*)
(* ::Function:: *)
(* f:LombScargle *)
(*Lomb-Scargle Periodogram*)
(***Options***) 
Options[LombScargle] = {FrequenciesOnly-> False,
	NormalizeIntensities->False,
	OversamplingRate -> 1, 
	PairReturn -> False,
	UpperFrequencyFactor -> 1};
(***Function***)
LombScargle[data_, setTimes_, OptionsPattern[]] :=
    Module[ {dataIn = N[data], inputTimes, inputData, dataAll, 
      oRate = N[OptionValue[OversamplingRate]], 
      upperFact = N[OptionValue[UpperFrequencyFactor]], 
      inputSetTimes = N[setTimes], fOnly = OptionValue[FrequenciesOnly], returnPair = OptionValue[PairReturn],norm = OptionValue[NormalizeIntensities],
       associationToggle = AssociationQ[data], key, inputTimesZeroed, n,
       window, f0, inputDataCentered, varianceInputPoints, freqStep, 
      freq, theta, cosAmpSquared, cosNorm, sinAmpSquared, sinNorm, 
      periodogram, returning},
        If[ associationToggle,
            key = Keys[dataIn][[1]];
            dataIn = Values[dataIn][[1]]
        ];
        If[ Length[dataIn[[1]]] == 
          0(*test the input first element whether it has length 0, i.e. 
         list of values with Missing[]/NonNumeric for missing values}*),
            dataAll = 
              Transpose[
               SortBy[#, First] &@
                DeleteCases[
                 Transpose[{inputSetTimes, dataIn}], {_, 
                  
                  Except[_?NumericQ] | _?
                    MissingQ}]];(*cleaned up the data to remove the missing \
points and sorted in sequence*),
            dataAll = Transpose[SortBy[dataIn, First]]
        ];
        (*above,
        we must take the external input and convert them to internal \
     variables,
        and also define our own internal variables that we are going to \
     need*)
        inputTimes = dataAll[[1]];
        inputData = 
         dataAll[[
          2]];(*adjust inputTimes starting point w.r.t.dataset times*)
        inputTimesZeroed = inputTimes - inputSetTimes[[1]];
        (*calculate the number of timepoints in the overall set*)
        n = N[Length[inputSetTimes]];
        (*calculate the time window of observation*);
        window = (Max[#] - Min[#]) &@inputSetTimes;
        (*invert this window to get the fundamental frequency,with the n/
        n-1 correction to account for the difference between the time \
     window and time period of the signal (WE ARE,FOR NOW,
        NOT INCLUDING THIS!)*)
        f0 = n/(n - 1)*1/window;
   (*subtract the mean from the inputData*)
        inputDataCentered = inputData - Mean[inputData];
        (*calculate a variance for the centered data*)
        varianceInputPoints = Variance[N[inputDataCentered]];
        (*define a frequency step*)
        freqStep = 1/(oRate*(N[Floor[n/2]] - 1))*(n/2*upperFact - 1)*f0;
        (*get the list of frequencies,
        adjusting both the lower frequency (to equal f0 0-
        effectively a lowpass filter) and the upper cutoff Nyquist by the \
     upper factor specified*)
        freq = Range[f0, n/2*upperFact*f0, freqStep];
        If[ fOnly,
            Return[Association@ MapIndexed[("f" <> ToString[Sequence @@ #2] -> #1) &]@freq]
        ];
   (*for a given frequency,
   a phaseshift angle \[Theta](f) is calculated based on the formula \
described above*)
        theta[f_] :=
            1/2 ArcTan[(Plus @@ Cos[4. \[Pi]*f*inputTimesZeroed]) + 10^-20, 
              Plus @@ Sin[4. \[Pi]*f*inputTimesZeroed]];
        (*now we calculate the different frequency components of our \
     spectrum,by projecting the cosine component...*)
        cosAmpSquared[
          f_] :=
            (Plus @@ (inputDataCentered*
            Cos[2. \[Pi]*f*inputTimesZeroed - theta[f]]))^2;
        (*...and then normalizing it...*)
        cosNorm[f_] :=
            Plus @@ (Cos[2 \[Pi] f*inputTimesZeroed - theta[f]]^2);
        (*...and then projecting the sine component...*)
        sinAmpSquared[
          f_] :=
            (Plus @@ (inputDataCentered*
            Sin[2. \[Pi]*f*inputTimesZeroed - theta[f]]))^2;
        (*...and normalizing it as well...*)
        sinNorm[f_] :=
            Plus @@ (Sin[2. \[Pi]*f*inputTimesZeroed - theta[f]]^2);
        (*put it together to get the periodogram*)
        periodogram = 
         1/(2 varianceInputPoints) (Chop@cosAmpSquared[#]/cosNorm[#] + 
              Chop@sinAmpSquared[#]/sinNorm[#]) & /@ freq;
        (*the function finally returns:1) the list of frequencies,
        2) the corresponding list of Lomb-Scargle spectra*)
        returning = If[ norm,
                        N[{freq, Normalize@periodogram}],
                        N[{freq, periodogram}]
                    ];
        If[ returnPair,
            returning = Transpose[returning]
        ];
        If[ associationToggle,
            returning = Association[key -> returning]
        ];
        Return[returning]
    ];
   
(* ::Function:: *)
(* f:InverseAutocovariance *)
(***Options***) 
Options[InverseAutocovariance] = {PairReturn -> True,
  UpperFrequencyFactor -> 1}; 
(***Function***)
InverseAutocovariance[data_, setTimes_, OptionsPattern[]] :=
    Module[ {dataIn = N[data], dataAll, 
      associationToggle = AssociationQ[data], key, inputTimes, inputData,
       upperFact = N[OptionValue[UpperFrequencyFactor]], 
      returnPair = OptionValue[PairReturn], inputSetTimes = N[setTimes], 
      inputTimesNormed, n, window, f0, inputDataCentered, 
      varianceInputPoints, freq, freqStep, theta, cosAmpSquared, cosNorm,
       sinAmpSquared, sinNorm, inverseAuto, 
      returning},(*we must take the external input and convert them to \
   internal variables,
     and also define our own internal variables that we are going to \
   need*)
        If[ associationToggle,
            key = Keys[dataIn][[1]];
            dataIn = Values[dataIn][[1]]
        ];
        If[ Length[dataIn[[1]]] == 
          0(*test the input first element whether it has length 0, i.e. 
         list of values with Missing[]/NonNumeric for missing values}*),
            dataAll = 
              Transpose[
               SortBy[#, First] &@
                DeleteCases[
                 Transpose[{inputSetTimes, dataIn}], {_, 
                  Except[_?NumericQ] | _?
                    MissingQ}]];(*cleaned up the data to remove the missing \
points and sorted in sequence*),
            dataAll = Transpose[SortBy[dataIn, First]]
        ];
        (*above,
        we must take the external input and convert them to internal \
      variables,
        and also define our own internal variables that we are going to \
      need*)
        inputTimes = dataAll[[1]];
        inputData = 
         dataAll[[
          2]];(*adjust inputTimes starting point w.r.t.dataset times,
        AND ZERO-PAD THEM*)
        inputTimesNormed = 
         Join[inputTimes, inputSetTimes + inputSetTimes[[-1]]] - 
          inputSetTimes[[1]];
        (*calculate the number of timepoints in the overall set-
        since we cut freqStep to half f0,
        we should compensate n likewise by multiplication by 2*)
        n = 2*N[Length[inputSetTimes]];
        (*calculate the time window of observation*);
        window = (Max[#] - Min[#]) &@inputSetTimes;
        (*invert this window to get the fundamental frequency,with the n/
        n-1 correction to account for the difference between the time \
      window and time period of the signal (WE REMOVE THIS FOR NOW)*)
        f0 = 1/window;
        (*subtract the mean from the inputData,BEFORE YOU ZERO-PAD IT!*)
        inputDataCentered = 
         Join[inputData - Mean[inputData], 
          ConstantArray[0, Length[inputSetTimes]]];
        (*calculate a variance for the centered data*)
        varianceInputPoints = Variance[inputDataCentered];
        (*define the frequency step as HALF the fundamental frequency in \
      order to zero-pad and get an evenly spaced mesh*)
        freqStep = 1/2*f0;
        (*get the list of frequencies*)
        freq = Range[freqStep, n*upperFact*freqStep, freqStep];
        (*for a given frequency,
        a phaseshift angle \[Theta](f) is calculated based on the formula \
      described above*)
        theta[f_] :=
            1/2 ArcTan[(Plus @@ Cos[4. \[Pi]*f*inputTimesNormed]) + 10^-20, 
              Plus @@ Sin[4. \[Pi]*f*inputTimesNormed]];
        (*now we calculate the different frequency components of our \
      spectrum,by projecting the cosine component...*)
        cosAmpSquared[
          f_] :=
            (Plus @@ (inputDataCentered*
            Cos[2. \[Pi]*f*inputTimesNormed - theta[f]]))^2;
        (*...and then normalizing it...*)
        cosNorm[f_] :=
            Plus @@ (Cos[2 \[Pi] f*inputTimesNormed - theta[f]]^2);
        (*...and then projecting the sine component...*)
        sinAmpSquared[
          f_] :=
            (Plus @@ (inputDataCentered*
            Sin[2. \[Pi]*f*inputTimesNormed - theta[f]]))^2;
        (*...and normalizing it as well...*)
        sinNorm[f_] :=
            Plus @@ (Sin[2. \[Pi]*f*inputTimesNormed - theta[f]]^2);
        (*put it together to get the inverse autocorrelation*)
        inverseAuto = 
         1/(2*varianceInputPoints) (Chop@cosAmpSquared[#]/cosNorm[#] + 
              Chop@sinAmpSquared[#]/sinNorm[#]) & /@ freq;
        (*the function finally returns:1) the list of frequencies,
        2) the correspoinding list of inverse autocovariances*)
        returning = N[{freq, inverseAuto}];
        If[ returnPair,
            returning = Transpose[returning]
        ];
        If[ associationToggle,
            returning = Association[key -> returning]
        ];
        Return[returning]
    ];

(* ::Function:: *)
(* f:Autocorrelation *)
(***Options***) 
Options[Autocorrelation] = {SpectrumFunction -> InverseAutocovariance,
	UpperFrequencyFactor -> 1}; 
(***Function***)
Autocorrelation[data_, setTimes_, OptionsPattern[]] :=
    Module[ {dataIn = N[data], associationToggle = AssociationQ[data], 
      key, dataAll, inputTimes, inputData, 
      upperFact = N[OptionValue[UpperFrequencyFactor]],
      spectrumFn = OptionValue[SpectrumFunction],
      inputSetTimes = N[setTimes], inputInverseAuto, nTDelays, 
      inverseAmplitudes, autoCorrs, norm, 
      returning},(*This layer above defines all external variables as \
   internal variables and defines all constants and functions we are \
   going to need to calculate the autocorrelation*)(*create the inverse \
   autocovariance from the inputs*)
        If[ associationToggle,
            key = Keys[dataIn][[1]];
            dataIn = Values[dataIn][[1]]
        ];
        If[ Length[dataIn[[1]]] == 
          0(*test the input first element whether it has length 0, i.e. 
         list of values with Missing[]/NonNumeric for missing values}*),
            dataAll = 
              Transpose[
               SortBy[#, First] &@
                DeleteCases[
                 Transpose[{inputSetTimes, dataIn}], {_, 
                  Except[_?NumericQ] | _?
                    MissingQ}]];(*cleaned up the data to remove the missing \
points and sorted in sequence*),
            dataAll = Transpose[SortBy[dataIn, First]]
        ];
        (*above,
        we must take the external input and convert them to internal \
      variables,
        and also define our own internal variables that we are going to \
      need*)
        inputTimes = dataAll[[1]];
        inputData = dataAll[[2]];
        inputInverseAuto = 
         spectrumFn[Transpose[{inputTimes, inputData}], inputSetTimes, 
          UpperFrequencyFactor -> upperFact, PairReturn -> False];
        (*we simply create the amplitude spectrum from the input data,
        making sure to add a zero at the first element to make the DFT work,
        and sampling only half the points because we have oversampled by \
      two in the inverseAutocovariance*)
        inverseAmplitudes = 
         Join[{0}, 
          inputInverseAuto[[2]][[1 ;; (Length[inputInverseAuto[[2]]]/2)]]];
        (*and then we go ahead and do the transform*)
        autoCorrs = FourierDCT[inverseAmplitudes, 3];
        (*finally we divide everything by a normalization factor so that \
      the autocorrelation at lag 0=1*)
        norm = autoCorrs[[1]];
  (*make sure we are only returning autocorrelations for the lags we \
can rely on,say,for up to N/2 time points*)
        nTDelays = Floor[1/2 Length[autoCorrs]];
        returning = (autoCorrs[[1 ;; nTDelays]]/norm);
        If[ associationToggle,
            returning = Association[key -> returning]
        ];
        Return[returning]
    ]; 


(* ::Section:: *)
(*#####Visualization####*)   

(* ::Function:: *)
(* f:Heatmapper*)
(***Options***)
Options[Heatmapper] = {ColorBlending -> {CMYKColor[1, 0, 1, 0], CMYKColor[0, 1, 1, 0]},
	ScaleShift -> 1}; 
(***Function***)
Heatmapper[data_, OptionsPattern[]] :=
    Module[ {inputData = data,
        shiftOfScale = OptionValue[ScaleShift],
        color = OptionValue[ColorBlending],
        mapperColorFunctionA},
        (*generate a mapper color function*)
        mapperColorFunctionA = Blend[color, #1*shiftOfScale] &;
        (*and do an array plot with it*)
        ArrayPlot[inputData, Mesh -> False, ColorFunctionScaling -> True, 
            ColorRules -> {Missing[] -> Gray, "Empty" -> Transparent}, 
            ColorFunction -> (mapperColorFunctionA), PlotRangePadding -> None, 
            AspectRatio -> 1]
    ];

(* ::Function:: *)
(* f:TimeSeriesDendrogramHeatmap::*)
(*Dendrogram and Heatmap Displays for Time Series*)
(***Options***)
Options[TimeSeriesDendrogramHeatmap] = {ColorBlending -> {CMYKColor[1, 0, 1, 0], CMYKColor[0, 1, 1, 0]},
	DendrogramColor -> Yellow,
	FrameName -> "Dendrogram and Heatmap", 
	GroupSubSize -> {0.1, 0.1},
	HorizontalAxisName -> "Time (arbitrary units)", 
    HorizontalLabels -> None,
    VerticalLabels -> None,
    IndexColor -> "DeepSeaColors",
    ImageSize -> 200, 
    ScaleShift -> None};
(***Function***)
TimeSeriesDendrogramHeatmap[data_, OptionsPattern[]] :=
    Module[ {
    aggInfo = If[ AssociationQ[data],
                  data["Cluster"],
                  data[[1]]
              ],
    splitClusters = 
    If[ AssociationQ[data],
        data["IntermediateClusters"],
        data[[3]]
    ],
    subsplitClusters = 
    If[ AssociationQ[data],
        data["SubsplitClusters"],
        data[[4]]
    ],
    inputData = 
    If[ AssociationQ[data],
        data["Data"][[All, 1]],
        data[[5]][[All, 1]]
    ],
    groupAssociations = 
    If[ AssociationQ[data],
        data["GroupAssociations"],
        data[[6]]
    ],
    horizAxis = OptionValue[HorizontalAxisName],
    xLabels = OptionValue[HorizontalLabels],
    yLabels = OptionValue[VerticalLabels],
    nameOfFrame = OptionValue[FrameName],
    groupSubSize = OptionValue[GroupSubSize],
    imageSize = OptionValue[ImageSize],
    shiftOfScale = OptionValue[ScaleShift],
    dendroColor = OptionValue[DendrogramColor],
    color = OptionValue[ColorBlending],
    indexColor = OptionValue[IndexColor],
    dataDimensions, splitLengths, subsplitLengths, indexLabel, 
    countsAssociationColoringLabels, countsAssociationClusterLabels, 
    heatmapPlot, dendroPlot, intensityLabel, legendLabels, 
    horizAxisLabel,verticalAxisLabel,returning},(*defined local variables for module*)
   (*define the scale shift so as to have zero coloring standardized - 
  this is problematic if Min is 0, so set to 1 if that is the case*)
        shiftOfScale = 
        If[ MatchQ[shiftOfScale, None],
            If[ N[Min[#]] != 0.,
                ((Max[#] - Min[#]))/(-2*Min[#]),
                1
            ] &@
              inputData /. Missing[] -> 0,
            shiftOfScale
        ];
       (*we will need data dimensions to pad*)
        dataDimensions = Dimensions[inputData];
        (*horizontal labels*)
        If[MatchQ[xLabels,None],xLabels=Table[" ", {dataDimensions[[2]]}]];
        (*we create a vector displaying the lengths of the split clusters*)
        splitLengths = (Plus @@ #[[-2 ;;]]) &@
            If[ MatchQ[Head[#], Cluster],
                #,
                {1, 0}
            ] & /@ splitClusters;
        (*and now one for the lengths of the subsplit clusters*)
        subsplitLengths = 
         Query[All /* Values, Length]@
          groupAssociations;(*(Plus@@#[[-2;;]])&@#&/@Flatten[
        subsplitClusters,1];*)
        (*now,
        we take the split and subsplit lengths,
        and create an index label of bars and sub-
        bars that will be displayed next to Heatmap-NEED TO ARRAYPAD*)
        indexLabel = 
         ArrayPlot[
          ArrayPad[
           Transpose[
            Join[ConstantArray[
              Flatten[Join[(ConstantArray[#[[1]], #[[2]]] & /@ 
                  Transpose[{Range[Length[splitLengths]], splitLengths}])]],
               Ceiling[dataDimensions[[2]]*groupSubSize[[1]]]], 
             ConstantArray[
              Flatten[Join[(ConstantArray[#[[1]], #[[2]]] & /@ 
                  Transpose[{Range[Length[subsplitLengths]], 
                    subsplitLengths}])]], 
              Ceiling[dataDimensions[[2]]*groupSubSize[[2]]]]]], {{0, 
             0}, {0, dataDimensions[[2]]}}, "Empty"], 
          ColorFunction -> indexColor, 
          ColorRules -> {"Empty" -> Transparent}];
        (*create a'count' of how many elements in each cluster based on (G,
        S) ID*)
        countsAssociationClusterLabels = (Counts[
        Transpose[
        Join[ConstantArray[
        Flatten[Join[(ConstantArray[#[[1]], #[[2]]] & /@ 
           Transpose[{Range[Length[splitLengths]], splitLengths}])]],
        Ceiling[dataDimensions[[2]]*groupSubSize[[1]]]], 
        ConstantArray[
        Flatten[Join[(ConstantArray[#[[1]], #[[2]]] & /@ 
           Transpose[{Flatten[
              Join[Range[Length[#]] & /@ subsplitClusters]], 
             subsplitLengths}])]], 
        Ceiling[dataDimensions[[2]]*groupSubSize[[2]]]]]]]);
        countsAssociationColoringLabels = (Counts[
           Transpose[
            Join[ConstantArray[
              Flatten[Join[(ConstantArray[#[[1]], #[[2]]] & /@ 
                  Transpose[{Range[Length[splitLengths]], splitLengths}])]],
               Ceiling[dataDimensions[[2]]*groupSubSize[[1]]]], 
             ConstantArray[
              Flatten[Join[(ConstantArray[#[[1]], #[[2]]] & /@ 
                  Transpose[{Range[Length[subsplitLengths]], 
                    subsplitLengths}])]], 
              Ceiling[dataDimensions[[2]]*groupSubSize[[2]]]]]]]);
        (*we now create dendrogram plot to be generated*)
        dendroPlot = 
        If[ MatchQ[Head[aggInfo], Cluster],
         ImageReflect[
          Quiet@DendrogramPlot[
           aggInfo, 
           LeafLabels -> None, TruncateDendrogram -> All, 
           HighlightLevel -> Length[splitClusters], 
           HighlightStyle -> dendroColor, Orientation -> Left]],Image[{{1}}]];
        (*we now create heatmap plot*)
        heatmapPlot = 
         Heatmapper[
          ArrayPad[
           N[inputData], {{0, 
             0}, {Ceiling[dataDimensions[[2]]*groupSubSize[[1]]] + 
              Ceiling[dataDimensions[[2]]*groupSubSize[[2]]], 0}}, "Empty"],
           ScaleShift -> shiftOfScale, ColorBlending -> color];
        (*Return[heatmapPlot];*)
        (*define Heatmap intensity label*)
        intensityLabel = 
         DensityPlot[y, {x, 0, .05}, {y, 0, 1}, AspectRatio -> Automatic, 
          ColorFunction -> (Blend[color, #1*shiftOfScale] &), 
          PlotRangePadding -> None, Frame -> True, 
          FrameTicks -> {{None, {{0, 
               NumberForm[
                First @(DeleteMissing@
                     If[ MatchQ[Head[#], Min],
                         Apply[List, #],
                         {#}
                     ] &@
                   Min@inputData), 2]}, {1, 
               NumberForm[
                First @(DeleteMissing@
                     If[ MatchQ[Head[#], Max],
                         Apply[List, #],
                         {#}
                     ] &@
                   Max@inputData), 2]}}}, {None, None}}];
        (*define legend labels*)
        legendLabels =
         First@(GraphicsGrid[{{ArrayPlot[#, ColorFunction -> indexColor, 
                 FrameTicks -> {{None, 
                    Transpose[{Range[Length@#], #}] &@
                     Query[All /* Normal, Length]@groupAssociations}, {{{1, 
                      "G"}, {Ceiling[
                        dataDimensions[[2]]*groupSubSize[[1]]] + 
                       Ceiling[dataDimensions[[2]]*groupSubSize[[2]]] - 1, 
                      "S"}}, None}}, AlignmentPoint -> Center]}}, 
              AspectRatio -> 1] & /@ {Keys[
              countsAssociationColoringLabels]});
        horizAxisLabel = 
         GraphicsRow[
          Join[{"G"}, 
           ConstantArray["", 
            Ceiling[dataDimensions[[2]]*groupSubSize[[1]]] - 1], {"S"}, 
           ConstantArray["", 
            Ceiling[dataDimensions[[2]]*groupSubSize[[2]]] - 1], xLabels], 
          AspectRatio -> 0.1, ImageSize -> {imageSize}];
        verticalAxisLabel=If[MatchQ[yLabels,None],Null,GraphicsGrid[Transpose@{yLabels},ImageSize-> {Automatic,imageSize},AspectRatio-> 3]];
        (*return Grid of dendrogram, heatmap, legends*)
        returning =
         Grid[{{"\n"},{Style[nameOfFrame, 
             Bold]}, {Grid[{{"Dendrogram", "HeatMap", Null,Null, 
               "GroupAssociations Index\nGroup#Subgroup#\[Rule]Size"}, \
        {Show[dendroPlot, AspectRatio -> 1, 
                ImageSize -> {imageSize*(31.7/30), Automatic}, 
                ImageMargins -> 0, ImagePadding -> 0], 
               Show[heatmapPlot, indexLabel, 
                ImageSize -> {imageSize, imageSize}, ImageMargins -> 0, 
                AspectRatio -> 1],verticalAxisLabel, intensityLabel, 
               Show[legendLabels, 
                ImageSize -> {0.75 imageSize, 0.75 imageSize}]},
              {Null, horizAxisLabel,Null,Null, Null},
              {Null, horizAxis,Null,Null, Null}}, Frame -> True, 
             Spacings -> 1(*Alignment\[Rule]{{Bottom,Bottom}}*)]}}];
        Return[returning]
    ];
  
(* ::Function:: *)
(* f:TimeSeriesDendrogramsHeatmaps*)
(***Options***)
Options[TimeSeriesDendrogramsHeatmaps] = {FunctionOptions ->{ImageSize -> 200}}; 
(***Function***)
TimeSeriesDendrogramsHeatmaps[data_, OptionsPattern[]] :=
    Module[ {dataIn = data, mapOptions = OptionValue[FunctionOptions]},
        Return[Column@(TimeSeriesDendrogramHeatmap[#[[2]], 
              FrameName -> #[[1]], Sequence @@ mapOptions] & /@ 
            Query[{Keys, Values} /* Transpose]@dataIn)]
    ];

(* ::Function:: *)
(* f:TimeSeriesSingleDendrogramHeatmap*)
(***Options***)    
Options[TimeSeriesSingleDendrogramHeatmap] = {ColorBlending -> {CMYKColor[1, 0, 1, 0], CMYKColor[0, 1, 1, 0]},
	DendrogramColor -> Yellow, 
	FrameName -> "Dendrogram and Heatmap",
	GroupSubSize -> 0.1,
	HorizontalAxisName -> "Time (arbitrary units)", 
    HorizontalLabels -> None,
    VerticalLabels-> None,
    IndexColor -> "DeepSeaColors",
    ImageSize -> 200,  
    ScaleShift -> None};
(***Function***)
TimeSeriesSingleDendrogramHeatmap[data_, OptionsPattern[]] := 
  Module[{inputData = 
     If[AssociationQ@data, 
      If[MatchQ[Head[data["Cluster"]], Cluster], 
       ClusterFlatten[data["Cluster"]][[All, 
         1]], {data["Cluster"][[1]]}], 
      If[MatchQ[Head[data[[1]]], Cluster], 
       ClusterFlatten[data[[1]]][[All, 1]], {data[[1, 1]]}]],
    aggInfo = If[AssociationQ[data], data["Cluster"], data[[1]]],
    splitClusters = 
     If[AssociationQ[data], data["InitialSplitCluster"], data[[2]]],
    groupAssociations = 
     If[AssociationQ[data], data["GroupAssociations"], data[[3]]],
    horizAxis = OptionValue[HorizontalAxisName],
    xLabels = OptionValue[HorizontalLabels],
    yLabels = OptionValue[VerticalLabels],
    nameOfFrame = OptionValue[FrameName],
    groupSubSize = OptionValue[GroupSubSize],
    imageSize = OptionValue[ImageSize],
    shiftOfScale = OptionValue[ScaleShift],
    dendroColor = OptionValue[DendrogramColor],
    color = OptionValue[ColorBlending],
    indexColor = OptionValue[IndexColor],
    dataDimensions, splitLengths, indexLabel, countsAssociation, 
    heatmapPlot, dendroPlot, intensityLabel, legendLabels, horizAxisLabel, verticalAxisLabel,
    returning},(*defined local variables*)
   (*define the scale \
shift so as to have zero coloring standardized*)
   shiftOfScale = 
    If[MatchQ[shiftOfScale, None], 
     If[N[Min[#]] != 0., ((Max[#] - Min[#]))/(-2*Min[#]), 1] &@
       inputData /. Missing[] -> 0, shiftOfScale];
   (*we will need data dimensions to pad*)
   
   dataDimensions = Dimensions[inputData];
   If[MatchQ[xLabels,None],xLabels=Table[" ", {dataDimensions[[2]]}]];
   (*we create a vector displaying the lengths of the split clusters*)
   splitLengths = (Plus @@ #[[-2 ;;]]) &@
       If[MatchQ[Head[#], Cluster], #, {1, 0}] & /@ splitClusters;
   (*now,we take the split lengths,
   and create an index label of bars and sub-
   bars that will be displayed next to Heatmap-NEED TO ARRAYPAD*)
   indexLabel = 
    ArrayPlot[
     ArrayPad[
      Transpose[
       ConstantArray[
        Flatten[Join[(ConstantArray[#[[1]], #[[2]]] & /@ 
            Transpose[{Range[Length[splitLengths]], splitLengths}])]],
         1*Ceiling[dataDimensions[[2]]*groupSubSize]]], {{0, 0}, {0, 
        dataDimensions[[2]]}}, "Empty"], ColorFunction -> indexColor, 
     ColorRules -> {"Empty" -> Transparent}];
   (*create a'count' of how many elements in each cluster based on (G,
   S) ID*)countsAssociation = (Counts[
      Transpose[
       ConstantArray[
        Flatten[
         Join[(ConstantArray[#[[1]], #[[2]]] & /@ 
            Transpose[{Range[Length[splitLengths]], splitLengths}])]],
         1*Ceiling[dataDimensions[[2]]*groupSubSize]]]]);
   (*we now create dendrogram plot to be generated*)
   dendroPlot = 
   If[MatchQ[Head[aggInfo], Cluster],ImageReflect[Quiet@DendrogramPlot[
       aggInfo, 
       LeafLabels -> None, TruncateDendrogram -> All, 
       HighlightLevel -> Length[splitClusters], 
       HighlightStyle -> dendroColor, Orientation -> Left]],Image[{{1}}]];
   (*we now create heatmap plot*)
   heatmapPlot = 
    Heatmapper[
     ArrayPad[
      N[inputData], {{0, 
        0}, {Ceiling[dataDimensions[[2]]*groupSubSize], 0}}, "Empty"],
      ScaleShift -> shiftOfScale, ColorBlending -> color];
   (*define Heatmap intensity label*)
   intensityLabel = 
    DensityPlot[y, {x, 0, .05}, {y, 0, 1}, AspectRatio -> Automatic, 
     ColorFunction -> (Blend[color, #1*shiftOfScale] &), 
     PlotRangePadding -> None, Frame -> True, 
     FrameTicks -> {{None, {{0, 
          
          NumberForm[
           First @(DeleteMissing@
                If[MatchQ[Head[#], Min], Apply[List, #], {#}] &@
              Min@inputData), 2]}, {1, 
          NumberForm[
           First @(DeleteMissing@
                If[MatchQ[Head[#], Max], Apply[List, #], {#}] &@
              Max@inputData), 2]}}}, {None, None}}];
   legendLabels =
    First@(GraphicsGrid[{{ArrayPlot[#, ColorFunction -> indexColor, 
            FrameTicks -> {{None, 
               Transpose[{Range[Length@#], #}] &@
                Query[All /* Normal, Length]@groupAssociations}, {{{1,
                  "G"}}, None}}, AlignmentPoint -> Center]}}, 
         AspectRatio -> 1] & /@ {Keys[countsAssociation]});
   (*define horizontal axis label*)
   
   horizAxisLabel = 
    GraphicsRow[
     Join[{"G"}, 
      ConstantArray["", 
       Ceiling[dataDimensions[[2]]*groupSubSize] - 1], xLabels], 
     AspectRatio -> 0.1, ImageSize -> {imageSize}];
   verticalAxisLabel=If[MatchQ[yLabels,None],Null,GraphicsGrid[Transpose@{yLabels},ImageSize-> {Automatic,imageSize},AspectRatio-> 3]];
   (*return what we need to return*)
   returning =
    Grid[{{"\n"},{Style[nameOfFrame, 
             Bold]}, {Grid[{{"Dendrogram", "HeatMap", Null,Null, 
          "GroupAssociations Index\nGroup#\[Rule]Size"}, {Show[
           dendroPlot, AspectRatio -> 1, 
           ImageSize -> {imageSize*(31/30), Automatic}, 
           ImageMargins -> 0, ImagePadding -> 0], 
          Show[heatmapPlot, indexLabel, 
           ImageSize -> {imageSize, imageSize}, ImageMargins -> 0, 
           AspectRatio -> 1],verticalAxisLabel, intensityLabel, 
          Show[legendLabels, 
           ImageSize -> {0.5 imageSize, 0.5 imageSize}]},
         {Null, horizAxisLabel, Null,Null,Null},
         {Null, horizAxis, Null,Null,Null}}, Frame -> True, 
        Spacings -> 1(*Alignment\[Rule]{{Bottom,Bottom}}*)]}}];
   Return[returning]];
   
(* ::Function:: *)
(* f:TimeSeriesSingleDendrogramsHeatmaps*)
(***Options***)
Options[TimeSeriesSingleDendrogramsHeatmaps] = {FunctionOptions -> \
{ImageSize -> 200}}; 
(***Function***)
TimeSeriesSingleDendrogramsHeatmaps[data_, OptionsPattern[]] :=
    Module[ {dataIn = data, mapOptions = OptionValue[FunctionOptions]},
        Return[Column@(TimeSeriesSingleDendrogramHeatmap[#[[2]], 
              FrameName -> #[[1]], Sequence @@ mapOptions] & /@ 
            Query[{Keys, Values} /* Transpose]@dataIn)]
    ];


(* ::Function:: *)
(* f:MatrixDendrogramsHeatmaps*)
(*Dendrogram and Heatmap Displays for Data Matrix Multiple classifications*)
(***Options***)
Options[MatrixDendrogramsHeatmaps] = {FunctionOptions -> \
{ImageSize -> 200}};
(***Function***) 
MatrixDendrogramsHeatmaps[data_, OptionsPattern[]] :=
    Module[ {dataIn = data, mapOptions = OptionValue[FunctionOptions]},
        Return[Column@(MatrixDendrogramHeatmap[#[[2]], 
              FrameName -> #[[1]], Sequence @@ mapOptions] & /@ 
            Query[{Keys, Values} /* Transpose]@dataIn)]
    ];
    
(* ::Function:: *)
(* f:MatrixDendrogramHeatmap*)
(* Dendrogram and Heatmap Displays for Data Matrix single classification *)
(***Options***)   
Options[MatrixDendrogramHeatmap] = {ColorBlending -> {CMYKColor[1, 0, 1, 0], CMYKColor[0, 1, 1, 0]}, 
	DendrogramColor -> {Yellow, Yellow},
	FrameName -> "Dendrogram and Heatmap",
	HorizontalAxisName -> "Groups", 
	HorizontalLabels -> None,
	VerticalLabels -> None,
	IndexColor -> {"DeepSeaColors", "CandyColors"},
	ImageSize -> 200,
	ScaleShift -> None};
(***Function***)   
MatrixDendrogramHeatmap[data_, OptionsPattern[]] := 
 Module[{inputData = 
   If[AssociationQ@data, 
    If[MatchQ[Head[data["RowCluster"]], Cluster], 
     ClusterFlatten[data["RowCluster"]][[All, 
       1]], {data["RowCluster"][[1]]}], 
    If[MatchQ[Head[data[[1]]], Cluster], 
     ClusterFlatten[data[[1]]][[All, 1]], {data[[1, 1]]}]], 
  aggInfoHoriz = If[AssociationQ@data, data["RowCluster"], data[[1]]],
   aggInfoVert = 
   If[AssociationQ@data, data["ColumnCluster"], data[[2]]], 
  splitClustersHoriz = 
   If[AssociationQ@data, data["RowSplitClusters"], data[[3]]], 
  splitClustersVert = 
   If[AssociationQ@data, data["ColumnSplitClusters"], data[[4]]], 
  associationsH = 
   If[AssociationQ@data, data["GroupAssociationsRows"], data[[5]]], 
  associationsV = 
   If[AssociationQ@data, data["GroupAssociationsColumns"], data[[6]]],
   imageSize = OptionValue[ImageSize], 
  horizAxis = OptionValue[HorizontalAxisName], 
  xLabels = OptionValue[HorizontalLabels],
  yLabels = OptionValue[VerticalLabels], 
  nameOfFrame = OptionValue[FrameName], 
  shiftOfScale = OptionValue[ScaleShift], 
  dendroColor = OptionValue[DendrogramColor], 
  color = OptionValue[ColorBlending], 
  indexColor = OptionValue[IndexColor], dataDimensions, 
  splitLengthsVertical, splitLengthsHorizontal, indexLabelVertical, 
  dendroPlotVertical, indexLabelHorizontal, dendroPlotHorizontal, 
  heatmapPlot, intensityLabel, horizAxisLabel,verticalAxisLabel, countsAssociationH, 
  legendLabelsH, countsAssociationV, legendLabelsV,dendroColorH,dendroColorV, 
  returning},(*the layer above transforms all the external variables \
into internal variables and defines all the other variables we are \
going to need*)(*we will need data dimensions to pad*)
 dataDimensions = Dimensions[inputData];
 If[MatchQ[Head[dendroColor],List],dendroColorH=dendroColor[[1]];dendroColorV=dendroColor[[2]],dendroColorH=dendroColor;dendroColorV=dendroColor];
 If[MatchQ[xLabels,None],xLabels={""}];
 If[!MatchQ[Head[indexColor],List],indexColor={indexColor,indexColor}];
 shiftOfScale = 
  If[MatchQ[shiftOfScale, None], 
   If[N[Min[#]] != 0., ((Max[#] - Min[#]))/(-2*Min[#]), 1] &@
     inputData /. Missing[] -> 0, shiftOfScale];
 (*we create vectors displaying the lengths of the split \
clusters,both vert and horiz*)
 splitLengthsVertical = (Plus @@ #[[-2 ;;]]) &@
     If[MatchQ[Head[#], Cluster], #, {1, 0}] & /@ splitClustersVert;
 splitLengthsHorizontal = (Plus @@ #[[-2 ;;]]) &@
     If[MatchQ[Head[#], Cluster], #, {1, 0}] & /@ splitClustersHoriz;
 (*now we create index labels and dendro plots,both vert and horiz*)
 indexLabelHorizontal = 
  ArrayPlot[
   Transpose[
    List@Flatten@
      Join[(ConstantArray[#[[1]], #[[2]]] & /@ 
         Transpose[{Range[Length[splitLengthsHorizontal]], 
           splitLengthsHorizontal}])]], 
   ColorFunction -> indexColor[[1]], 
   ColorRules -> {"Empty" -> Transparent}, Frame -> False, 
   AspectRatio -> 10];
 countsAssociationH = (Counts[
    Transpose[
     ConstantArray[
      Flatten[Join[(ConstantArray[#[[1]], #[[2]]] & /@ 
          Transpose[{Range[Length[splitLengthsHorizontal]], 
            splitLengthsHorizontal}])]], 1]]]);
 indexLabelVertical = 
  ArrayPlot[
   List@Flatten[
     Join[(ConstantArray[#[[1]], #[[2]]] & /@ 
        Transpose[{Range[Length[splitLengthsVertical]], 
          splitLengthsVertical}])]], ColorFunction -> indexColor[[2]],
    ColorRules -> {"Empty" -> Transparent}, Frame -> False, 
   AspectRatio -> 1/7];
 countsAssociationV = (Counts[
    Transpose[
     ConstantArray[
      Flatten[Join[(ConstantArray[#[[1]], #[[2]]] & /@ 
          Transpose[{Range[Length[splitLengthsVertical]], 
            splitLengthsVertical}])]], 1]]]);
 dendroPlotVertical = 
  If[MatchQ[Head[aggInfoVert], Cluster], 
   Quiet@DendrogramPlot[aggInfoVert, LeafLabels -> None, 
     TruncateDendrogram -> All, 
     HighlightLevel -> Length[splitClustersVert], 
     HighlightStyle -> dendroColorV], Image[{{1}}]];
 dendroPlotHorizontal = 
  If[MatchQ[Head[aggInfoHoriz], Cluster], 
   ImageReflect@
    Quiet@DendrogramPlot[aggInfoHoriz, LeafLabels -> None, 
      TruncateDendrogram -> All, 
      HighlightLevel -> Length[splitClustersHoriz], 
      HighlightStyle -> dendroColorH, Orientation -> Left], 
   Image[{{1}}]];
 (*and now a heatmap plot*)
 heatmapPlot = 
  Heatmapper[
   N[Transpose[Transpose[inputData][[ClusterFlatten@aggInfoVert]]]], 
   ScaleShift -> shiftOfScale, ColorBlending -> color];
 (*define Heatmap intensity label*)
 intensityLabel = 
  Evaluate[DensityPlot[y, {x, 0, .05}, {y, 0, 1}, 
    AspectRatio -> Automatic, 
    ColorFunction -> (Blend[color, #1*shiftOfScale] &), 
    PlotRangePadding -> None, Frame -> True, 
    FrameTicks -> {{None, {{0, 
         NumberForm[Min[#] &@Cases[N@Flatten[inputData], _Real], 
          2]}, {1, NumberForm[Max[#] &@Cases[N@Flatten[inputData], _Real], 
          2]}}}, {None, None}}]];
 (*define horizontal axis label*)
 legendLabelsH = 
  First@(GraphicsGrid[{{ArrayPlot[#, ColorFunction -> indexColor[[1]],
           FrameTicks -> {{None, 
             Transpose[{Range[Length@#], #}] &@
              Query[All /* Normal, Length]@associationsH}, {{{1, 
               "GH"}}, None}}, AlignmentPoint -> Center]}}, 
       AspectRatio -> 1] & /@ {Keys[countsAssociationH]});
 legendLabelsV = 
  First@(GraphicsGrid[{{ArrayPlot[#, ColorFunction -> indexColor[[2]],
           FrameTicks -> {{None, 
             Transpose[{Range[Length@#], #}] &@
              Query[All /* Normal, Length]@associationsV}, {{{1, 
               "GV"}}, None}}, AlignmentPoint -> Center]}}, 
       AspectRatio -> 1] & /@ {Keys[countsAssociationV]});
 horizAxisLabel = 
  GraphicsRow[xLabels, AspectRatio -> 0.1, 
   ImageSize -> {imageSize}];
  verticalAxisLabel = verticalAxisLabel=If[MatchQ[yLabels,None],Null,GraphicsGrid[Transpose@{yLabels},ImageSize-> {Automatic,imageSize},AspectRatio-> 3]];
 (*returning*)
 returning = 
  Grid[{{"\n"}, {Style[nameOfFrame, Bold], 
     SpanFromLeft}, {Grid[{{Grid[{{"GroupAssociations", 
           SpanFromLeft}, {"GroupIndex" -> "Size", 
           SpanFromLeft}, {"Rows", 
           "Columns"}, {Show[legendLabelsH, 
            ImageSize -> {Max[100, imageSize/2], 
              Max[100, imageSize/2]}], 
           Show[legendLabelsV, 
            ImageSize -> {Max[100, imageSize/2], 
              Max[100, imageSize/2]}]}}],Null, 
        Show[dendroPlotVertical, ImageSize -> {imageSize, Automatic}],Null,
         Null}, {Null,Null, 
        Show[indexLabelVertical, 
         ImageSize -> {imageSize*(313/300), Automatic}], 
        Null,Null}, {Show[dendroPlotHorizontal, AspectRatio -> 1, 
         ImageSize -> {imageSize*(317/300), Automatic}, 
         ImageMargins -> 0, ImagePadding -> 0], 
        Show[indexLabelHorizontal, 
         ImageSize -> {Automatic, imageSize*4013/4000}], 
        Show[heatmapPlot, ImageSize -> {imageSize, imageSize}], 
       verticalAxisLabel, Show[intensityLabel, 
         ImageSize -> {Automatic, 0.75*imageSize}]}, {Null, Null, 
        horizAxisLabel, Null,Null}, {Null, Null, horizAxis,Null,Null}}, 
      Frame -> True, Spacings -> 1]}}, Alignment -> Center];
 Return[returning]];

(* ::Function:: *)  
(* f:KEGGPathwayVisual*)
(*Possible inputs: KEGGPathway: "KEGG Pathway Name", "ORA Object"
Optional Data Intensities: Association:  "Gene" -> Values
Output -> one of {"movie", "URL","figure"}*)
(***Options***)
Options[KEGGPathwayVisual] = {AnalysisType -> "Genomic"(*options are "Genomic", "Molecular"*),
	AugmentDictionary -> True,
	BlendColors -> {RGBColor[0, 0, 1], RGBColor[0, 0, 1], 
     RGBColor[0.5, 0.5, 0.5], RGBColor[1, 0, 0], 
     RGBColor[1, 0, 0]} (*note RGB style*),
    ColorSelection -> Association@{"RNA" -> "bg", "Protein" -> "fg"},
    DefaultColors -> {"fg" -> RGBColor[0, 0, 0], "bg" -> RGBColor[0, 1, 0]},
	ExportMovieOptions -> {"VideoEncoding" -> "MPEG-4 Video", "FrameRate" -> 1},
	FileExtend -> ".mov",
	GeneDictionary -> None,
	GetGeneDictionaryOptions -> {},
	InputID -> {"UniProt ID", "Gene Symbol"},
	Intensities -> None,
	IntensityFunction -> (Blend[#1, (#2 + 1)/2] &)(*where#1 is BlendColors and #2 is intensity vector*),
	KEGGAnalysisAssignerOptions -> {},
	KEGGDatabase -> "pathway",
	KEGGMolecular -> "cpd",
	KEGGOrganism -> "hsa",
	MathIOmicaDataDirectory -> 
    If[ ! DirectoryQ[
       FileNameJoin[
        Flatten[{FileNameSplit[$UserBaseDirectory], "Applications", 
          "MathIOmica", "MathIOmicaData"}]]],
        CreateDirectory[
         FileNameJoin[
          Flatten[{FileNameSplit[$UserBaseDirectory], "Applications", 
            "MathIOmica", "MathIOmicaData"}]]],
        FileNameJoin[
         Flatten[{FileNameSplit[$UserBaseDirectory], "Applications", 
           "MathIOmica", "MathIOmicaData"}]]
    ],
	MemberSet -> All (*Options are All, "ORA" output for pathway 
	[in which case ORA\[Rule] True must be set], or a list of **Keys** {Key[identifier]}*),
	MissingValueColor ->  RGBColor[0.4, 0.4, 0.4](*note RGB color*),
	MolecularInputID -> {"cpd"},
	MolecularOutputID -> "cpd",
	MolecularSpecies -> "compound",
	MovieFilePath -> None,
	NonUCSC -> False,
	ORA -> False,
	OutputID -> "KEGG Gene ID",
	ResultsFormat -> "URL"(*"URL","Figure","Movie"*),
	SingleColorPlace -> "bg",
	Species -> "human",
	StandardHighlight -> {"fg" -> RGBColor[1, 0, 0], 
     "bg" -> RGBColor[0.5, 0.70, 1]}};
(***Function***)
KEGGPathwayVisual[KEGGPathway_, OptionsPattern[]] :=
    Module[ {
    pathNm = KEGGPathway,
    type = OptionValue[AnalysisType],
    getGeneDictionaryOptions = OptionValue[GetGeneDictionaryOptions],
    augmentDictionary = OptionValue[AugmentDictionary],
    keggAnalysisAssignOpt = OptionValue[KEGGAnalysisAssignerOptions],
    keggOrg = OptionValue[KEGGOrganism],
    keggMol = OptionValue[KEGGMolecular],
    keggDb = OptionValue[KEGGDatabase],
    inputID = OptionValue[InputID],
    outputID = OptionValue[OutputID],
    molInID = OptionValue[MolecularInputID],
    molOutID = OptionValue[MolecularOutputID],
    members = OptionValue[MemberSet],
    ora = OptionValue[ORA],
    geneDict = OptionValue[GeneDictionary],
    species = OptionValue[Species],
    nonUCSC = OptionValue[NonUCSC],
    data = OptionValue[Intensities],
    stdC = Association@OptionValue[StandardHighlight],
    defaultC = Association@OptionValue[DefaultColors],
    intFn = OptionValue[IntensityFunction],
    blendC = OptionValue[BlendColors],
    missingC = OptionValue[MissingValueColor],
    colorBgFg = OptionValue[ColorSelection],
    singleClrP = OptionValue[SingleColorPlace],
    dataReturn = OptionValue[ResultsFormat],
    movieFile = OptionValue[MovieFilePath],
    movieOpts = OptionValue[ExportMovieOptions],
    ext = OptionValue[FileExtend],
    dir = FileNameSplit[OptionValue[MathIOmicaDataDirectory]],
    molSpecies = OptionValue[MolecularSpecies],
    fileMolDict,
    hexMap,
    pathAssign,
    pathIDs,
    mapIDs,
    group,
    returning,
    urls},
  (*If no individul members of pathway provided in either method return the \
pathway only*)
        If[ MatchQ[members, All] && MatchQ[ora, False] && MatchQ[data, None],
            urls = {"http://www.kegg.jp/kegg-bin/show_pathway?map=" <> (StringSplit[#, 
                    ":"][[-1]] &@pathNm)},
            hexMap = 
             StringJoin[(IntegerString[Round[{#[[1]], #[[2]], #[[3]]}*255], 
                  16] /. {"0" -> "00", "1" -> "01", "2" -> "02", "3" -> "03", 
                  "4" -> "04", "5" -> "05", "6" -> "06", "7" -> "07", "8" -> "08", 
                  "9" -> "09", "a" -> "0a", "b" -> "0b", "c" -> "0c", "d" -> "0d", 
                  "e" -> "0e", "f" -> "0f"})] &;
            If[!MatchQ[data,None],If[!ListQ[data[[1]]],data=Query[All,{#}&]@data]];
            Switch[type,
             "Genomic",
             (*Gene based mapping to pathways*)
             
             pathAssign = If[ MatchQ[keggAnalysisAssignOpt, {}],
               (*If no specific options for function use full and species request*)
                              KEGGAnalysisAssigner[KEGGQuery1 -> keggDb, KEGGQuery2 -> keggOrg],
                              KEGGAnalysisAssigner[Sequence @@ keggAnalysisAssignOpt]
                          ];
             If[ ValueQ[ConstantGeneDictionary],(*variable Exists*)
                 If[ augmentDictionary,(*augment*)
                     ConstantGeneDictionary = 
                     Join[ConstantGeneDictionary, 
                     If[ MatchQ[geneDict, None],
                         GetGeneDictionary[Sequence @@ getGeneDictionaryOptions],
                         geneDict
                     ]],
      (*replace*)
                     ConstantGeneDictionary = 
                     If[ MatchQ[geneDict, None],
                         GetGeneDictionary[Sequence @@ getGeneDictionaryOptions],
                         geneDict
                     ]
                 ],
                 (*create/load UCSC based translation dictionary -
                 NB global variable or use specified variable*)
                 ConstantGeneDictionary = 
                  If[ MatchQ[geneDict, None],
                      GetGeneDictionary[Sequence @@ getGeneDictionaryOptions],
                      geneDict
                  ]
             ];,
             "Molecular",
             (*Molecular based mapping to pathways*)
             inputID = molInID;
             outputID = molOutID;
             species = molSpecies;
             nonUCSC = True;
             keggOrg = keggMol;
             fileMolDict = 
              FileNameJoin[Flatten[{dir, {"MathIOmicaMolecularDictionary"}}]];
             If[ ValueQ[ConstantGeneDictionary],(*variable Exists*)
                 If[ augmentDictionary,(*augment*)
                     ConstantGeneDictionary = 
                     Join[ConstantGeneDictionary, 
                     If[ FileExistsQ[fileMolDict],
                         Get[fileMolDict],
                         Return["Could not find annotation file at " <> fileMolDict <> 
                           " Please either obtain an annotation file from mathiomica.org or \
provide a GeneDictionary option variable."]
                     ]],
      (*replace*)
                     ConstantGeneDictionary = 
                     If[ MatchQ[geneDict, None],
                         If[ FileExistsQ[fileMolDict],
                             Get[fileMolDict],
                             Return["Could not find annotation file at " <> fileMolDict <> 
                               " Please either obtain an annotation file from mathiomica.org or \
provide a GeneDictionary option variable."]
                         ],
                         geneDict
                     ]
                 ],
                 (*create/load UCSC based translation dictionary -
                 NB global variable or use specified variable*)
                 ConstantGeneDictionary = 
                  If[ MatchQ[geneDict, None],
                      If[ FileExistsQ[fileMolDict],
                          Get[fileMolDict],
                          Return["Could not find annotation file at " <> fileMolDict <> 
                            " Please either obtain an annotation file from mathiomica.org or \
provide a GeneDictionary option variable."]
                      ],
                      geneDict
                  ]
             ];(*get the right KEGG terms for the background requested and \
         correct species*)
             pathAssign = If[ MatchQ[keggAnalysisAssignOpt, {}],
               (*If no specific options for function use background, species request, 
               length request*)
                              KEGGAnalysisAssigner[KEGGQuery1 -> keggDb, KEGGQuery2 -> keggOrg ],
                              KEGGAnalysisAssigner[Sequence @@ keggAnalysisAssignOpt]
                          ],
             __,
             (*non match*)
             
             Return["AnalysisType\[Rule]" <> 
               If[ StringQ[type],
                   "\"" <> type <> "\"",
                   ToString@type
               ] <> 
               " is not a valid  choice."]];
            pathIDs = 
             DeleteMissing@
              Flatten@Values@
                Values@GeneTranslation[
                  If[ nonUCSC,
                      #,
                      StringSplit[#, ":"][[-1]]
                  ] & /@ 
                   Query[keggOrg, "PathToID", pathNm]@pathAssign, 
                  inputID, ConstantGeneDictionary, InputID -> {outputID}, 
                  Species -> species];
            If[ MatchQ[data, None],
                members = 
                 If[ MatchQ[ora, False],
                     Select[IntersectingQ[If[ ListQ[#],
                                              #,
                                              {#}
                                          ], pathIDs] &]@members,
                     Flatten[#, 1] &@members[[-1, 2]]
                 ];
                (*No data, will highlight with standard colors *)
                mapIDs = DeleteMissing@(If[ Length[#] > 2 || MissingQ[#[[1]]],
                                            Missing[],
                                            #
                                        ] &@(Join @@ {Union@(If[ MissingQ[#],
                                                                 #,
                                                                 If[ nonUCSC,
                                                                     #,
                                                                     keggOrg <> ":" <> #
                                                                 ]
                                                             ] & /@ (Flatten@
                              
                              Values@GeneTranslation[If[ ListQ[#],
                                                         #[[{1}]],
                                                         {#}
                                                     ], 
                                outputID, ConstantGeneDictionary, InputID -> inputID, 
                                Species -> species])), If[ ListQ[#],
                                                           {#[[-1]]},
                                                           {}
                                                       ]}) & /@ 
                    members);
                group = 
                 Join @@@ Query[(GroupBy[First])(*/*Map[
                     Association](*Make an association*)*)/* (Map[
                       Map[(Query[{#}, List@hexMap]@stdC) &@(If[ Length[#] >= 2,
                                                                 colorBgFg[#[[-1]]],
                                                                 singleClrP
                                                             ]) &]])]@ mapIDs;,
                (*Intensities provided, will highlight as given*)
                members = 
                 If[ MatchQ[ora, False],
                     Key[#] & /@ members,
                     Key[#] & /@ Flatten[#, 1] &@members[[-1, 2]]
                 ];
                mapIDs = 
                 KeyDrop[Missing[]]@
                  KeyMap[If[ Length[#] > 2 || MissingQ[#[[1]]],
                             Missing[],
                             #
                         ] &@(Join @@ {Union@(If[ MissingQ[#],
                                                  #,
                                                  If[ nonUCSC,
                                                      #,
                                                      keggOrg <> ":" <> #
                                                  ]
                                              ] & /@ (Flatten@
                             Values@GeneTranslation[If[ ListQ[#],
                                                        #[[{1}]],
                                                        {#}
                                                    ], 
                               outputID, ConstantGeneDictionary, InputID -> inputID, 
                               Species -> species])), If[ ListQ[#],
                                                          {#[[-1]]},
                                                          {}
                                                      ]}) &, 
                   DeleteMissing@
                    Query[
                      If[ MatchQ[members, All],
                          KeySelect[IntersectingQ[If[ ListQ[#],
                                                      #,
                                                      {#}
                                                  ], pathIDs] &],
                          members
                      ], If[ ListQ[#] || MissingQ[#],
                             #,
                             {#}
                         ] &]@data];
                group = 
                 Query[Normal /* ((*Group By Member Name*)GroupBy[#[[1, 1]] &]) /* 
                    Map[Association](*Make an association*)/* ((Map[
                       KeyMap[If[ Length[#] >= 2,
                                  colorBgFg[#[[-1]]],
                                  singleClrP
                              ] &]])), 
                   All,(*Matching the numbers to the string hex*)
                   hexMap /@ ((*Matching the intensities to a color*)
                     If[ MatchQ[#1, _Missing],
                         missingC,
                         intFn[blendC, #]
                     ] &)]@(KeySelect[
                     MemberQ[Query[keggOrg, "PathToID", pathNm]@pathAssign, #[[1]]] &]@
                    mapIDs)
            ];
            urls = Query[
               All /* KeyValueMap[
                 StringJoin[StringReplace[#1, ":" -> "%3A"] <> "+", #2] &] /* 
                StringJoin /* ("http://www.kegg.jp/kegg-bin/show_pathway?map=" <> \
            (StringSplit[#, ":"][[-1]] &@pathNm) <> "&multi_query=" <> # &), 
               If[ MissingQ[#["bg"]],
                   "%23" <> (hexMap@defaultC["bg"]) <> "%2C",
                   "%23" <> #["bg"] <> "%2C"
               ] <> 
                 If[ MissingQ[#["fg"]],
                     "%23" <> (hexMap@defaultC["fg"]) <> "%0D%0A",
                     "%23" <> #["fg"] <> "%0D%0A"
                 ] &] /@ 
              Query[All /* Transpose(*{Map[{StringReplace[#1,":"\[Rule] "%3A"]<>
                "+%23"&/@Keys[#],Values[#]}&]}*), Transpose]@group
        ];
        returning = Switch[dataReturn,
          "URL", urls,
          "Figure",(*return Figure list*)
          (*get Image Location*)
          
          Flatten@((Import[
                  "http://www.kegg.jp" <> #] &@(DeleteCases[
                    StringCases[
                     Cases[#, 
                      XMLElement["img", {___, "src" -> src_, ___}, _] :> 
                       src, \[Infinity]], ___ ~~ ".png"], {}] &@
                  Import[#, "XMLObject"])) & /@ urls),
          "Movie",(*return Movies *)
          If[FileExistsQ[If[ MatchQ[movieFile, None],
                     StringReplace[pathNm,":"->"_"] <> ext,
                     movieFile
                 ]],Return[Print["File " <> If[ MatchQ[movieFile, None],
                     StringReplace[pathNm,":"->"_"] <> ext,
                     movieFile
                 ]<>" already exists. Please choose a different file name or rename/move existing file and retry."]]];
          Export[If[ MatchQ[movieFile, None],
                     StringReplace[pathNm,":"->"_"] <> ext,
                     movieFile
                 ], 
           Flatten@((Import[
                   "http://www.kegg.jp" <> #] &@(DeleteCases[
                     StringCases[
                      Cases[#, 
                       XMLElement["img", {___, "src" -> src_, ___}, _] :> 
                        src, \[Infinity]], "/tmp" ~~ ___ ~~ ".png"], {}] &@
                   Import[#, "XMLObject"])) & /@ urls), Sequence @@ movieOpts],
          __, Print["Improper Return Selection, Returning URL for online patway."];
              urls];
        Return[<|"Pathway"-> pathNm, "Results"-> returning|>]
    ];             

(* ::Section:: *)
(* End Section *)

End[]

EndPackage[]
